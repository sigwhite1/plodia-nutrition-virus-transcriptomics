# R/09_enrichment_GO_clusters.R
# Purpose: GO enrichment for expression clusters via Drosophila ortholog mapping.
#
# PIPELINE:
#   ProbeName -> maker -> Dmel_protein (DIAMOND BLASTP hit) -> UniProt accession
#   -> FlyBase gene ID (FBgn) -> gost() GO enrichment
#
# INPUTS:
#   cluster_annot_all_raw.rds      — all cluster probes with DIAMOND hits (evalue<=1e-5)
#   cluster_annot_all_filtered.rds — quality-filtered hits (evalue<=1e-5, pident>=35%)
#                                    used as the custom background set
#
# COLUMN NAMING:
#   cluster_annot_all uses lowercase 'cluster' (integer 1-4), from R_08.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/cluster_annot_all_raw.rds
#   - intermediate/cluster_annot_all_filtered.rds
#
# Writes:
#   - intermediate/cluster_fbgn.rds
#   - intermediate/GO_cluster_enrichment.rds
#   - results/tables/GO_cluster_enrichment.csv
#   - results/figures/GO_cluster_enrichment_top_terms.png

# -----------------------------
# 0) Load artifacts
# -----------------------------
config <- readRDS(file.path("intermediate", "config.rds"))
paths  <- config$paths

cluster_annot_all          <- readRDS(file.path(paths$intermediate_dir,
                                                "cluster_annot_all_raw.rds"))
cluster_annot_all_filtered <- readRDS(file.path(paths$intermediate_dir,
                                                "cluster_annot_all_filtered.rds"))

dir.create(paths$tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)

# Minimum number of FBgn IDs required to attempt enrichment
min_genes_for_enrichment <- 10

# -----------------------------
# 1) Helper functions
# -----------------------------

# Extract UniProt accession from DIAMOND FASTA-style header
# e.g. "tr|Q9VFH5|Q9VFH5_DROME" -> "Q9VFH5"
#      "sp|P35220|CTNA_DROME"    -> "P35220"
extract_uniprot_one <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- sub("^[^|]+\\|([^|]+)\\|.*$", "\\1", x)
  out[out == x] <- ""   # pattern didn't match -> drop
  out
}

# Clean a vector of FASTA-style headers -> unique UniProt accessions
clean_fasta_ids <- function(ids) {
  ids <- ids[!is.na(ids) & ids != ""]
  ids <- sub("^[^|]+\\|([^|]+)\\|.*$", "\\1", ids)
  ids <- ids[ids != ""]
  unique(ids)
}

# Convert UniProt accessions -> FlyBase gene IDs (FBgn...)
to_fbgn <- function(uniprot_ids) {
  if (length(uniprot_ids) == 0) return(character(0))
  out <- gprofiler2::gconvert(
    query     = unique(uniprot_ids),
    organism  = "dmelanogaster",
    target    = "FLYBASE_GENE_ID",
    filter_na = TRUE
  )
  if (is.null(out) || nrow(out) == 0) return(character(0))
  unique(out$target)
}

# Extract tidy GO table from gost() and compute BH-FDR within label
extract_go <- function(go_obj, label,
                       keep_sources = c("GO:BP", "GO:MF", "GO:CC")) {
  if (is.null(go_obj) || is.null(go_obj$result) ||
      nrow(go_obj$result) == 0) return(NULL)
  
  go_obj$result %>%
    dplyr::filter(source %in% keep_sources) %>%
    dplyr::transmute(
      label,
      source,
      term_id,
      term_name,
      p_value,
      intersection_size,
      term_size
    ) %>%
    dplyr::group_by(label) %>%
    dplyr::mutate(
      fdr_bh       = p.adjust(p_value, method = "BH"),
      neglog10_fdr = -log10(pmax(fdr_bh, 1e-300)),
      term_name    = stringr::str_trunc(term_name, 60)
    ) %>%
    dplyr::ungroup()
}

# -----------------------------
# 2) Cluster -> UniProt -> FBgn
# -----------------------------
# Build batched UniProt -> FBgn mapping in one gconvert call
# to avoid k+1 separate API calls (one per cluster + one for audit)

# Step 2a: extract unique UniProt accessions across all clusters
uniprot_unique <- cluster_annot_all %>%
  dplyr::filter(!is.na(Dmel_protein)) %>%
  dplyr::mutate(uniprot = extract_uniprot_one(Dmel_protein)) %>%
  dplyr::filter(uniprot != "") %>%
  dplyr::distinct(uniprot) %>%
  dplyr::pull(uniprot)

message("Unique UniProt accessions to convert: ", length(uniprot_unique))

# Step 2b: single batched gconvert call (avoids per-cluster API overhead)
old_timeout <- getOption("timeout")
options(timeout = max(old_timeout, 120))

conv <- tryCatch(
  gprofiler2::gconvert(
    query     = uniprot_unique,
    organism  = "dmelanogaster",
    target    = "FLYBASE_GENE_ID",
    filter_na = TRUE
  ),
  error = function(e) {
    message("gconvert failed: ", conditionMessage(e))
    NULL
  }
)

options(timeout = old_timeout)

# Step 2c: build uniprot -> FBgn lookup table
if (is.null(conv) || nrow(conv) == 0) {
  warning("gconvert returned no results. FBgn mapping will be empty.")
  uniprot_to_fbgn <- tibble::tibble(
    uniprot = uniprot_unique,
    fbgn    = rep(list(character(0)), length(uniprot_unique)),
    n_fbgn  = 0L
  )
} else {
  uniprot_to_fbgn <- conv %>%
    dplyr::transmute(
      uniprot = as.character(input),
      fbgn    = as.character(target)
    ) %>%
    dplyr::group_by(uniprot) %>%
    dplyr::summarise(
      fbgn   = list(unique(fbgn)),
      n_fbgn = dplyr::n_distinct(fbgn),
      .groups = "drop"
    ) %>%
    dplyr::right_join(
      tibble::tibble(uniprot = as.character(uniprot_unique)),
      by = "uniprot"
    ) %>%
    dplyr::mutate(
      fbgn   = purrr::map(fbgn,
                          ~ if (is.null(.x) || all(is.na(.x))) character(0) else .x),
      n_fbgn = dplyr::coalesce(n_fbgn, 0L)
    )
}

message("UniProt -> FBgn mapping summary (n_fbgn per UniProt):")
print(table(uniprot_to_fbgn$n_fbgn))

# Step 2d: build per-cluster FBgn lists by joining through the lookup table
cluster_fbgn <- cluster_annot_all %>%
  dplyr::filter(!is.na(Dmel_protein)) %>%
  dplyr::mutate(uniprot = extract_uniprot_one(Dmel_protein)) %>%
  dplyr::filter(uniprot != "") %>%
  dplyr::left_join(uniprot_to_fbgn %>%
                     dplyr::select(uniprot, fbgn),
                   by = "uniprot") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    uniprot_ids = list(unique(uniprot)),
    fbgn_ids    = list(unique(unlist(fbgn[!sapply(fbgn, is.null)]))),
    .groups     = "drop"
  ) %>%
  dplyr::arrange(cluster)

# QA: UniProt and FBgn counts per cluster
qa <- cluster_fbgn %>%
  dplyr::transmute(
    cluster,
    n_uniprot = purrr::map_int(uniprot_ids, length),
    n_fbgn    = purrr::map_int(fbgn_ids,    length)
  )
message("Cluster FBgn mapping QA:")
print(qa)

saveRDS(cluster_fbgn,
        file = file.path(paths$intermediate_dir, "cluster_fbgn.rds"))

# -----------------------------
# 3) Background set: FBgn from filtered microarray probes
# -----------------------------
# Using the filtered annotation as background is more biologically honest:
# we compare each cluster only against genes that were actually measurable
# in this experiment, not against all ~14k Drosophila genes.
bg_uniprot <- cluster_annot_all_filtered %>%
  dplyr::filter(!is.na(Dmel_protein)) %>%
  dplyr::mutate(uniprot = extract_uniprot_one(Dmel_protein)) %>%
  dplyr::filter(uniprot != "") %>%
  dplyr::pull(uniprot) %>%
  unique()

bg_fbgn <- uniprot_to_fbgn %>%
  dplyr::filter(uniprot %in% bg_uniprot) %>%
  dplyr::pull(fbgn) %>%
  unlist() %>%
  unique()
bg_fbgn <- bg_fbgn[!is.na(bg_fbgn) & bg_fbgn != ""]

message("Background FBgn IDs (filtered microarray): ", length(bg_fbgn))

use_custom_bg <- TRUE   # set FALSE to use full Drosophila genome as background

# -----------------------------
# 4) Run GO enrichment per cluster
# -----------------------------
cluster_gost <- cluster_fbgn %>%
  dplyr::mutate(
    gost_obj = purrr::map(
      fbgn_ids,
      ~ {
        ids <- .x
        if (length(ids) < min_genes_for_enrichment) {
          message("  Cluster skipped: only ", length(ids), " FBgn IDs (< ",
                  min_genes_for_enrichment, ")")
          return(NULL)
        }
        if (use_custom_bg && length(bg_fbgn) >= min_genes_for_enrichment) {
          gprofiler2::gost(ids,
                           organism  = "dmelanogaster",
                           custom_bg = bg_fbgn)
        } else {
          gprofiler2::gost(ids, organism = "dmelanogaster")
        }
      }
    )
  )

# Report term counts per cluster
cluster_gost %>%
  dplyr::transmute(
    cluster,
    n_fbgn  = purrr::map_int(fbgn_ids,  length),
    n_terms = purrr::map_int(gost_obj,
                             ~ if (is.null(.x)) 0L else nrow(.x$result))
  ) %>%
  print()

# Extract tidy GO table
cluster_go <- purrr::pmap_dfr(
  list(cluster_gost$gost_obj, cluster_gost$cluster),
  function(go_obj, cl) extract_go(go_obj, paste0("Cluster_", cl))
)

# Guard against empty results before saving / plotting
if (is.null(cluster_go) || nrow(cluster_go) == 0) {
  warning(
    "No significant GO terms found for any cluster.\n",
    "Consider: (1) lowering the FDR threshold in gost(), ",
    "(2) setting use_custom_bg = FALSE, or ",
    "(3) checking that FBgn IDs are valid Drosophila identifiers."
  )
} else {
  message("GO enrichment complete: ", nrow(cluster_go), " terms across ",
          dplyr::n_distinct(cluster_go$label), " clusters")
  
  saveRDS(cluster_go,
          file = file.path(paths$intermediate_dir, "GO_cluster_enrichment.rds"))
  readr::write_csv(cluster_go,
                   file = file.path(paths$tables_dir, "GO_cluster_enrichment.csv"))
  
  message("Saved GO results to: ",
          file.path(paths$intermediate_dir, "GO_cluster_enrichment.rds"))
  message("Exported CSV to: ",
          file.path(paths$tables_dir, "GO_cluster_enrichment.csv"))
  
  # Top terms summary per cluster
  message("\nTop 3 GO terms per cluster (by -log10 FDR):")
  cluster_go %>%
    dplyr::group_by(label) %>%
    dplyr::slice_min(fdr_bh, n = 3, with_ties = FALSE) %>%
    dplyr::select(label, source, term_name, fdr_bh) %>%
    print(n = Inf)
  
  # -----------------------------
  # 5) Plot: top terms per cluster
  # -----------------------------
  # Layout: Cluster 1 top-left, [empty] top-right,
  #          Cluster 3 bottom-left, Cluster 4 bottom-right.
  # Achieved by inserting a dummy " " level between Cluster_1 and Cluster_3
  # so facet_wrap(ncol=2) leaves the top-right panel empty.
  cluster_go_top <- cluster_go %>%
    dplyr::filter(label %in% c("Cluster_1", "Cluster_3", "Cluster_4")) %>%
    dplyr::mutate(
      label = factor(label,
                     levels = c("Cluster_1", " ",
                                "Cluster_3", "Cluster_4"))
    ) %>%
    dplyr::group_by(label) %>%
    dplyr::slice_max(neglog10_fdr, n = 8, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(label) %>%
    dplyr::mutate(
      term_name = forcats::fct_reorder(term_name, neglog10_fdr)
    ) %>%
    dplyr::ungroup()
  
  p_cluster_go <- ggplot2::ggplot(
    cluster_go_top,
    ggplot2::aes(x    = neglog10_fdr,
                 y    = term_name,
                 size = intersection_size,
                 colour = source)
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~ label, ncol = 2, scales = "free_y",
                        drop = FALSE) +
    ggplot2::labs(
      title  = "GO enrichment across expression clusters",
      x      = expression(-log[10]("FDR (BH)")),
      y      = NULL,
      size   = "Genes in term",
      colour = "GO domain"
    ) +
    ggplot2::theme_bw(base_size = 12)
  
  out_png <- file.path(paths$figures_dir,
                       "GO_cluster_enrichment_top_terms.png")
  ggplot2::ggsave(out_png, plot = p_cluster_go,
                  width = 14, height = 10, dpi = 300)
  message("Saved GO plot to: ", out_png)
  
  print(p_cluster_go)
}
