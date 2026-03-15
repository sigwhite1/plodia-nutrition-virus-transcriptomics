# R/11_immune_gene_evaluation.R
# Purpose: Map immune-annotated genes onto Plodia microarray probes using
#          GO-term based grouping rather than Drosophila canonical pathway
#          membership. Probes are assigned to immune categories based on the
#          GO:BP annotations of their best Drosophila BLAST ortholog.
#
# RATIONALE:
#   Canonical Drosophila immune gene sets (Toll, Imd, RNAi, etc.) rely on
#   protein-level orthology between Plodia and Drosophila. DIAMOND BLASTP
#   analysis confirmed that most canonical immune signalling genes (PGRP-LC,
#   MyD88, Ago2, Dcr-2, hop, Stat92E, bsk, DptA) have no detectable Plodia
#   ortholog above the search thresholds (evalue <= 1e-5), consistent with
#   the high evolutionary divergence of immune genes between Lepidoptera and
#   Diptera. A GO-term based approach is therefore used instead:
#   probes are grouped by the immune GO:BP categories of their Dmel orthologs,
#   giving a Plodia-native view of immune-related gene expression.
#
# IMMUNE GO CATEGORIES (GO:BP):
#   GO:0006955  immune response
#   GO:0045087  innate immune response
#   GO:0006952  defense response
#   GO:0051607  defense response to virus        <- most relevant for PigV
#   GO:0009615  response to virus
#   GO:0050829  defense response to Gram-negative bacterium
#   GO:0050830  defense response to Gram-positive bacterium
#
# PIPELINE:
#   cluster_annot_all_raw.rds
#     -> UniProt accession (via extract_uniprot_one)
#     -> FBgn (gconvert)
#     -> GO:BP annotations (gost, evcodes=TRUE, significant=FALSE)
#     -> filter immune GO terms
#     -> ProbeName -> immune category mapping
#     -> heatmap of logFC across diet x time
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/results_all.rds
#   - intermediate/cluster_annot_all_raw.rds
#
# Writes:
#   - intermediate/immune_go_probe_tbl.rds
#   - intermediate/immune_plot_df.rds
#   - intermediate/immune_map_qc.rds
#   - results/tables/immune_go_probe_summary.csv
#   - results/tables/immune_mapping_qc.csv
#   - results/figures/immune_heatmap_fullrange.png
#   - results/figures/immune_heatmap_contrast_robust.png

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(tibble)
  library(scales)
  library(purrr)
  library(gprofiler2)
  library(readr)
  library(stringr)
})

stopifnot(dir.exists("intermediate"))

# -----------------------------
# 0) Load config
# -----------------------------
config            <- readRDS(file.path("intermediate", "config.rds"))
paths             <- config$paths
time_levels_hours <- config$time_levels_hours

dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,  showWarnings = FALSE, recursive = TRUE)

# Helper: extract UniProt accession from DIAMOND FASTA-style header
extract_uniprot_one <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- sub("^[^|]+\\|([^|]+)\\|.*$", "\\1", x)
  out[out == x] <- ""
  out
}

# -----------------------------
# A) Define immune GO categories
# -----------------------------
# Ordered from most specific (virus defense) to broadest.
# Each probe will be assigned to ALL matching categories, then labelled
# by its most specific one for heatmap faceting.
immune_go_terms <- tibble::tibble(
  term_id   = c("GO:0051607", "GO:0009615",
                "GO:0050829", "GO:0050830",
                "GO:0045087", "GO:0006952",
                "GO:0006955"),
  term_name = c("Defense response to virus",
                "Response to virus",
                "Defense response to Gram-neg. bacteria",
                "Defense response to Gram-pos. bacteria",
                "Innate immune response",
                "Defense response",
                "Immune response"),
  # Specificity rank: 1 = most specific, 7 = broadest
  specificity = c(1L, 2L, 3L, 4L, 5L, 6L, 7L)
)

message("Immune GO categories defined: ", nrow(immune_go_terms))
print(immune_go_terms %>% dplyr::select(term_id, term_name))

# -----------------------------
# B) Load probe -> Dmel FBgn mapping
# -----------------------------
annot_raw <- readRDS(file.path(paths$intermediate_dir,
                               "cluster_annot_all_raw.rds")) %>%
  tibble::as_tibble()

probe_uniprot <- annot_raw %>%
  dplyr::filter(!is.na(Dmel_protein)) %>%
  dplyr::transmute(
    ProbeName = as.character(ProbeName),
    uniprot   = extract_uniprot_one(Dmel_protein)
  ) %>%
  dplyr::filter(uniprot != "", !is.na(uniprot)) %>%
  dplyr::distinct()

message("Probes with Dmel UniProt hit: ",
        dplyr::n_distinct(probe_uniprot$ProbeName))

# UniProt -> FBgn (batched single call)
old_timeout <- getOption("timeout")
options(timeout = max(old_timeout, 120))

fbgn_conv <- tryCatch(
  gprofiler2::gconvert(
    query     = unique(probe_uniprot$uniprot),
    organism  = "dmelanogaster",
    target    = "FLYBASE_GENE_ID",
    filter_na = TRUE
  ),
  error = function(e) {
    message("gconvert failed: ", conditionMessage(e)); NULL
  }
)

options(timeout = old_timeout)

if (is.null(fbgn_conv) || nrow(fbgn_conv) == 0) {
  stop("UniProt -> FBgn conversion failed. Check internet connection.")
}

uniprot_to_fbgn <- fbgn_conv %>%
  dplyr::transmute(uniprot = as.character(input),
                   FBgn    = as.character(target)) %>%
  dplyr::distinct()

# Build probe -> FBgn table
probe_fbgn <- probe_uniprot %>%
  dplyr::left_join(uniprot_to_fbgn, by = "uniprot",
                   relationship = "many-to-many") %>%
  dplyr::filter(!is.na(FBgn), FBgn != "") %>%
  dplyr::distinct(ProbeName, FBgn)

message("Probe -> FBgn pairs: ", nrow(probe_fbgn),
        " (", dplyr::n_distinct(probe_fbgn$ProbeName), " probes, ",
        dplyr::n_distinct(probe_fbgn$FBgn), " FBgn IDs)")

fbgn_ids <- unique(probe_fbgn$FBgn)

# -----------------------------
# C) Get GO:BP annotations for all Dmel orthologs
# -----------------------------
message("\nRunning gost to retrieve GO:BP annotations ",
        "(significant=FALSE, evcodes=TRUE)...")
message("  Query: ", length(fbgn_ids), " FBgn IDs")

go_result <- tryCatch(
  gprofiler2::gost(
    query       = fbgn_ids,
    organism    = "dmelanogaster",
    sources     = "GO:BP",
    evcodes     = TRUE,
    significant = FALSE
  ),
  error = function(e) {
    message("gost failed: ", conditionMessage(e)); NULL
  }
)

if (is.null(go_result) || is.null(go_result$result)) {
  stop("gost returned no results. Check internet connection and FBgn IDs.")
}

message("GO:BP terms retrieved: ", nrow(go_result$result))

# -----------------------------
# D) Filter for immune GO terms -> probe -> category table
# -----------------------------
immune_go_results <- go_result$result %>%
  dplyr::filter(term_id %in% immune_go_terms$term_id) %>%
  # Drop gost's term_name before joining so our curated names don't conflict
  dplyr::select(-term_name) %>%
  dplyr::left_join(immune_go_terms, by = "term_id") %>%
  dplyr::select(term_id, term_name, specificity,
                intersection_size, intersection)

message("\nImmune GO terms recovered:")
print(immune_go_results %>%
        dplyr::select(term_id, term_name, specificity, intersection_size) %>%
        dplyr::arrange(specificity))

# Expand: one row per FBgn x GO term
probe_go_long <- immune_go_results %>%
  dplyr::mutate(
    FBgn = purrr::map(intersection, ~ stringr::str_trim(strsplit(.x, ",")[[1]]))
  ) %>%
  tidyr::unnest(FBgn) %>%
  dplyr::select(FBgn, term_id, term_name, specificity) %>%
  dplyr::distinct()

message("FBgn IDs with at least one immune GO annotation: ",
        dplyr::n_distinct(probe_go_long$FBgn))

# Join to probes
probe_immune <- probe_fbgn %>%
  dplyr::inner_join(probe_go_long, by = "FBgn",
                    relationship = "many-to-many") %>%
  dplyr::distinct(ProbeName, FBgn, term_id, term_name, specificity)

message("Probe-GO term pairs: ", nrow(probe_immune))
message("Unique probes with immune annotation: ",
        dplyr::n_distinct(probe_immune$ProbeName))

# Assign most specific category per probe x FBgn
probe_category <- probe_immune %>%
  dplyr::group_by(ProbeName, FBgn) %>%
  dplyr::slice_min(specificity, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(ProbeName, FBgn,
                category    = term_name,
                category_go = term_id,
                specificity)

message("\nProbes per immune category (most specific assignment):")
print(probe_category %>% dplyr::count(category, sort = TRUE))

saveRDS(probe_immune,
        file.path(paths$intermediate_dir, "immune_go_probe_tbl.rds"))
readr::write_csv(
  probe_category %>%
    dplyr::count(category, category_go, name = "n_probes") %>%
    dplyr::arrange(dplyr::desc(n_probes)),
  file.path(paths$tables_dir, "immune_go_probe_summary.csv")
)

# -----------------------------
# E) Load results_all and build immune_plot_df
# -----------------------------
results_all <- readRDS(file.path(paths$intermediate_dir,
                                 "results_all.rds")) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(ProbeName = as.character(ProbeName))

req_cols <- c("ProbeName", "logFC", "adj.P.Val", "diet", "time_h")
missing_cols <- setdiff(req_cols, colnames(results_all))
if (length(missing_cols) > 0) {
  stop("results_all missing columns: ", paste(missing_cols, collapse = ", "))
}

# Representative probe per FBgn per category: lowest median adj.P.Val
fbgn_rep <- probe_category %>%
  dplyr::inner_join(
    results_all %>% dplyr::select(ProbeName, adj.P.Val),
    by = "ProbeName"
  ) %>%
  dplyr::group_by(category, FBgn, ProbeName) %>%
  dplyr::summarise(med_adjP = median(adj.P.Val, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(category, FBgn, med_adjP) %>%
  dplyr::group_by(category, FBgn) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

# Build plot dataframe
immune_plot_df <- fbgn_rep %>%
  dplyr::select(category, FBgn, ProbeName) %>%
  dplyr::inner_join(
    results_all %>%
      dplyr::select(ProbeName, diet, time_h, logFC, adj.P.Val),
    by = "ProbeName"
  ) %>%
  dplyr::distinct()

message("\nimmune_plot_df: ", nrow(immune_plot_df), " rows | ",
        dplyr::n_distinct(immune_plot_df$ProbeName), " probes | ",
        dplyr::n_distinct(immune_plot_df$FBgn), " FBgn IDs | ",
        dplyr::n_distinct(immune_plot_df$category), " categories")

# Order FBgn within each category by mean |logFC|
gene_order <- immune_plot_df %>%
  dplyr::group_by(category, FBgn) %>%
  dplyr::summarise(mean_abs = mean(abs(logFC), na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(category, dplyr::desc(mean_abs))

# Factor levels: categories ordered by specificity, FBgn by mean |logFC|
category_levels <- immune_go_terms %>%
  dplyr::filter(term_name %in% unique(immune_plot_df$category)) %>%
  dplyr::arrange(specificity) %>%
  dplyr::pull(term_name)

immune_plot_df <- immune_plot_df %>%
  dplyr::left_join(gene_order, by = c("category", "FBgn")) %>%
  dplyr::group_by(category) %>%
  dplyr::mutate(
    FBgn = factor(FBgn,
                  levels = unique(FBgn[order(mean_abs, decreasing = TRUE)]))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    time_h   = factor(time_h, levels = time_levels_hours, ordered = TRUE),
    diet     = factor(diet, levels = c("low", "high"),
                      labels = c("Low resource", "High resource")),
    category = factor(category, levels = category_levels)
  )

saveRDS(immune_plot_df,
        file.path(paths$intermediate_dir, "immune_plot_df.rds"))

# -----------------------------
# F) QC summary
# -----------------------------
qc_summary <- immune_go_terms %>%
  dplyr::left_join(
    probe_immune %>%
      dplyr::group_by(term_id) %>%
      dplyr::summarise(n_probes = dplyr::n_distinct(ProbeName),
                       n_fbgn   = dplyr::n_distinct(FBgn),
                       .groups  = "drop"),
    by = "term_id"
  ) %>%
  dplyr::mutate(
    n_probes = dplyr::coalesce(n_probes, 0L),
    n_fbgn   = dplyr::coalesce(n_fbgn,  0L)
  ) %>%
  dplyr::arrange(specificity)

message("\nImmune probe mapping QC:")
print(qc_summary %>% dplyr::select(term_name, n_probes, n_fbgn))

readr::write_csv(qc_summary,
                 file.path(paths$tables_dir, "immune_mapping_qc.csv"))
saveRDS(list(qc_summary    = qc_summary,
             probe_category = probe_category,
             probe_immune   = probe_immune),
        file.path(paths$intermediate_dir, "immune_map_qc.rds"))

# -----------------------------
# G) Heatmap plots
# -----------------------------
if (nrow(immune_plot_df) == 0) {
  warning("No immune probes to plot — check GO term coverage.")
} else {
  lim_full   <- max(abs(immune_plot_df$logFC), na.rm = TRUE)
  lim_robust <- max(
    as.numeric(quantile(abs(immune_plot_df$logFC), 0.95, na.rm = TRUE)),
    0.5
  )
  
  base_heat <- function(df, lim, title_suffix = "") {
    ggplot2::ggplot(df, ggplot2::aes(x = time_h, y = FBgn, fill = logFC)) +
      ggplot2::geom_tile() +
      ggplot2::facet_grid(category ~ diet, scales = "free_y",
                          space = "free_y") +
      ggplot2::scale_fill_gradient2(
        low      = scales::muted("steelblue"),
        mid      = "white",
        high     = scales::muted("firebrick"),
        midpoint = 0,
        limits   = c(-lim, lim),
        oob      = scales::squish,
        name     = "logFC"
      ) +
      ggplot2::labs(
        title    = paste0("Immune-annotated genes (GO:BP) across diet \u00d7 time",
                          title_suffix),
        subtitle = paste0("Representative probe per Dmel ortholog; ",
                          "grouped by most specific immune GO category"),
        x        = "Time after infection (hours)",
        y        = "FlyBase gene ID (Dmel ortholog)"
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y      = ggplot2::element_text(size = 8),
        strip.background = ggplot2::element_rect(fill = "grey85",
                                                 color = "grey30"),
        strip.text.y     = ggplot2::element_text(face = "bold", size = 9,
                                                 angle = 0),
        strip.text.x     = ggplot2::element_text(face = "bold")
      )
  }
  
  p_full   <- base_heat(immune_plot_df, lim_full,   " (full range)")
  p_robust <- base_heat(immune_plot_df, lim_robust,
                        paste0(" (\u00b1", round(lim_robust, 2), " robust)"))
  
  print(p_robust)
  
  # Auto-scale figure height to number of gene rows + category label space
  # Cap at 48 inches (ggsave default limit); use limitsize=FALSE if exceeded
  n_rows <- dplyr::n_distinct(immune_plot_df$FBgn)
  n_cats <- dplyr::n_distinct(immune_plot_df$category)
  fig_h  <- min(48, max(8, ceiling(n_rows * 0.18 + n_cats * 0.6)))
  
  out_full   <- file.path(paths$figures_dir, "immune_heatmap_fullrange.png")
  out_robust <- file.path(paths$figures_dir, "immune_heatmap_contrast_robust.png")
  
  ggplot2::ggsave(out_full,   p_full,   width = 14, height = fig_h,
                  dpi = 300, limitsize = FALSE)
  ggplot2::ggsave(out_robust, p_robust, width = 14, height = fig_h,
                  dpi = 300, limitsize = FALSE)
  
  message("Saved: ", out_full)
  message("Saved: ", out_robust)
  message("Figure height auto-scaled to: ", fig_h, " inches (",
          n_rows, " gene rows, ", n_cats, " categories)")
}

