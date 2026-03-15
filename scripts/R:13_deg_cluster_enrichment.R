# R/13_deg_cluster_enrichment.R
# Purpose: Test whether DEGs at each diet x time contrast are non-randomly
#          distributed across expression clusters using Fisher's exact test.
#          Produces a log2-enrichment heatmap and a stacked DEG composition plot.
#
# RATIONALE:
#   If clusters capture biologically coherent expression patterns, we expect
#   DEGs at specific timepoints to be over-represented in particular clusters.
#   For example, early-response DEGs (1h, 4h) should enrich in clusters with
#   early-peaking trajectories. This connects the DEG analysis (R_04/R_05) to
#   the clustering analysis (R_08).
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/results_all.rds     (logFC + adj.P.Val per probe x contrast)
#   - intermediate/cluster_df.rds      (ProbeName -> cluster)
#
# Writes:
#   - results/tables/deg_cluster_enrichment_by_contrast.csv
#   - results/figures/deg_cluster_enrichment_heatmap.png
#   - results/figures/deg_cluster_composition_stacked.png

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(readr)
  library(scales)   # squish()
})

# -----------------------------
# 0) Load config + inputs
# -----------------------------
config            <- readRDS(file.path("intermediate", "config.rds"))
paths             <- config$paths
time_levels_hours <- config$time_levels_hours
fdr_cutoff        <- config$thresholds$fdr

dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,  showWarnings = FALSE, recursive = TRUE)

results_all <- readRDS(file.path(paths$intermediate_dir, "results_all.rds")) %>%
  tibble::as_tibble()

cluster_df <- readRDS(file.path(paths$intermediate_dir, "cluster_df.rds")) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(
    ProbeName = as.character(ProbeName),
    cluster   = as.factor(cluster)
  ) %>%
  dplyr::distinct()

# Validate required columns
req_cols <- c("ProbeName", "diet", "time_h", "adj.P.Val")
missing_cols <- setdiff(req_cols, colnames(results_all))
if (length(missing_cols) > 0) {
  stop("results_all missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

# -----------------------------
# 1) Join clusters onto results and flag DEGs
# -----------------------------
res_clust <- results_all %>%
  dplyr::mutate(ProbeName = as.character(ProbeName)) %>%
  dplyr::inner_join(cluster_df, by = "ProbeName") %>%
  dplyr::mutate(
    is_deg = !is.na(adj.P.Val) & adj.P.Val < fdr_cutoff,
    diet   = factor(diet,   levels = c("low", "high")),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

message("Probes with cluster assignments: ",
        dplyr::n_distinct(res_clust$ProbeName))
message("Total DEG calls at FDR < ", fdr_cutoff, ": ",
        sum(res_clust$is_deg, na.rm = TRUE))

# -----------------------------
# 2) Count DEGs per cluster per diet x time
# -----------------------------
deg_counts <- res_clust %>%
  dplyr::filter(!is.na(adj.P.Val)) %>%
  dplyr::group_by(diet, time_h, cluster) %>%
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_deg    = sum(is_deg),
    .groups  = "drop"
  )

totals <- res_clust %>%
  dplyr::filter(!is.na(adj.P.Val)) %>%
  dplyr::group_by(diet, time_h) %>%
  dplyr::summarise(
    N_tested = dplyr::n(),
    N_deg    = sum(is_deg),
    .groups  = "drop"
  )

deg_enrich <- deg_counts %>%
  dplyr::left_join(totals, by = c("diet", "time_h")) %>%
  dplyr::mutate(
    exp_deg      = N_deg * (n_tested / N_tested),
    enrich_ratio = (n_deg    / pmax(n_tested, 1)) /
      (N_deg    / pmax(N_tested, 1))
  )

# -----------------------------
# 3) Fisher's exact test per cluster x diet x time
# -----------------------------
fisher_tbl <- deg_enrich %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p_fisher = {
      a <- n_deg
      b <- n_tested - n_deg
      c <- N_deg    - n_deg
      d <- (N_tested - n_tested) - (N_deg - n_deg)
      if (any(c(a, b, c, d) < 0)) NA_real_
      else fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(diet, time_h) %>%
  dplyr::mutate(fdr_fisher = p.adjust(p_fisher, method = "BH")) %>%
  dplyr::ungroup()

readr::write_csv(fisher_tbl,
                 file.path(paths$tables_dir,
                           "deg_cluster_enrichment_by_contrast.csv"))
message("Saved enrichment table: ",
        file.path(paths$tables_dir, "deg_cluster_enrichment_by_contrast.csv"))

# Quick summary: which clusters are significantly enriched most often?
message("\nSignificant enrichments (FDR < 0.05) per cluster:")
print(fisher_tbl %>%
        dplyr::filter(!is.na(fdr_fisher), fdr_fisher < 0.05) %>%
        dplyr::count(cluster, name = "n_sig_contrasts") %>%
        dplyr::arrange(dplyr::desc(n_sig_contrasts)))

# -----------------------------
# 4) Enrichment heatmap
# -----------------------------
plot_df <- fisher_tbl %>%
  dplyr::mutate(
    log2_enrich = log2(enrich_ratio),
    sig         = dplyr::if_else(!is.na(fdr_fisher) & fdr_fisher < 0.05, "*", "")
  )

# Diet labels for facet
diet_facet_labels <- c("low" = "Low resource", "high" = "High resource")

p_heat <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = time_h, y = cluster, fill = log2_enrich)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sig), size = 4) +
  ggplot2::facet_wrap(~ diet, ncol = 1,
                      labeller = ggplot2::as_labeller(diet_facet_labels)) +
  ggplot2::scale_fill_gradient2(
    midpoint = 0,
    limits   = c(-2, 2),
    oob      = scales::squish,
    low      = "#4575b4",
    mid      = "white",
    high     = "#d73027",
    name     = "log2 enrichment\n(DEG fraction)"
  ) +
  ggplot2::labs(
    title    = paste0("DEG enrichment by cluster within each diet \u00d7 time",
                      " (FDR < ", fdr_cutoff, ")"),
    subtitle = "Tile = log2 enrichment of DEGs vs background; * = Fisher BH-FDR < 0.05",
    x        = "Time after infection (hours)",
    y        = "Cluster"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey85", color = "grey30"),
    strip.text       = ggplot2::element_text(face = "bold")
  )

print(p_heat)
out_heat <- file.path(paths$figures_dir, "deg_cluster_enrichment_heatmap.png")
ggplot2::ggsave(out_heat, p_heat, width = 12, height = 6, dpi = 300)
message("Saved: ", out_heat)

# -----------------------------
# 5) Stacked DEG composition plot
#    (what fraction of DEGs at each timepoint belong to each cluster?)
# -----------------------------
deg_comp <- res_clust %>%
  dplyr::filter(!is.na(adj.P.Val), is_deg) %>%
  dplyr::count(diet, time_h, cluster, name = "n_deg") %>%
  dplyr::group_by(diet, time_h) %>%
  dplyr::mutate(frac_deg = n_deg / sum(n_deg)) %>%
  dplyr::ungroup()

p_comp <- ggplot2::ggplot(
  deg_comp,
  ggplot2::aes(x = time_h, y = frac_deg, fill = cluster)
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ diet, ncol = 1,
                      labeller = ggplot2::as_labeller(diet_facet_labels)) +
  ggplot2::labs(
    title    = paste0("Composition of DEGs across clusters (FDR < ", fdr_cutoff, ")"),
    subtitle = "Fraction of all DEGs at each timepoint attributed to each cluster",
    x        = "Time after infection (hours)",
    y        = "Fraction of DEGs",
    fill     = "Cluster"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid  = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey85", color = "grey30"),
    strip.text       = ggplot2::element_text(face = "bold")
  )

print(p_comp)
out_comp <- file.path(paths$figures_dir, "deg_cluster_composition_stacked.png")
ggplot2::ggsave(out_comp, p_comp, width = 12, height = 6, dpi = 300)
message("Saved: ", out_comp)
