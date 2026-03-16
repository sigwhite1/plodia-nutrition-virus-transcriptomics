# R/10_cluster_summaries.R
# Purpose: Create cluster-wise trajectory summaries and basic statistics/plots.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/clustering_inputs.rds   (logFC_mat, cluster_df)
#   - intermediate/results_all.rds         (raw logFC long table)
#   - intermediate/contrasts.rds           (contrast_key with diet/time_h)
#
# Writes:
#   - intermediate/plot_df_logFC.rds
#   - intermediate/cluster_summary.rds
#   - results/tables/cluster_summary.csv
#   - results/figures/cluster_spaghetti_rawlogFC.png
#   - results/figures/cluster_mean_CI_rawlogFC.png

# -----------------------------
# 0) Load artifacts
# -----------------------------
config <- readRDS(file.path("intermediate", "config.rds"))
paths             <- config$paths
time_levels_hours <- config$time_levels_hours   # numeric: 0.5 1 2 ... 168

clust_in      <- readRDS(file.path(paths$intermediate_dir, "clustering_inputs.rds"))
results_all   <- readRDS(file.path(paths$intermediate_dir, "results_all.rds"))
contrasts_obj <- readRDS(file.path(paths$intermediate_dir, "contrasts.rds"))

dir.create(paths$tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)

# Cluster assignments with factor-encoded cluster
cluster_df <- clust_in$cluster_df %>%
  dplyr::mutate(cluster = factor(cluster))

# Contrast metadata — diet/time_h columns (lowercase, numeric time)
contrast_key <- contrasts_obj$contrast_key %>%
  dplyr::filter(is_estimable) %>%
  dplyr::select(diet, time_h, contrast_name) %>%
  dplyr::mutate(
    diet   = factor(diet,   levels = c("low", "high")),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

# Diet color palette and labels
diet_colors <- c("high" = "#F8766D", "low" = "#00BFC4")
diet_labels <- c("high" = "High resource", "low" = "Low resource")

# -----------------------------
# 1) Build long table of raw logFC for genes used in clustering
# -----------------------------
deg_probes <- rownames(clust_in$logFC_mat)

plot_df <- results_all %>%
  dplyr::filter(ProbeName %in% deg_probes) %>%
  dplyr::select(ProbeName, contrast, logFC) %>%      # 'contrast' is lowercase
  dplyr::left_join(cluster_df,   by = "ProbeName") %>%
  dplyr::left_join(contrast_key, by = c("contrast" = "contrast_name")) %>%
  dplyr::filter(!is.na(diet), !is.na(time_h)) %>%
  dplyr::mutate(
    diet   = factor(diet,   levels = c("low", "high")),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

message("plot_df: ", nrow(plot_df), " rows (",
        dplyr::n_distinct(plot_df$ProbeName), " probes x ",
        dplyr::n_distinct(plot_df$contrast), " contrasts)")

saveRDS(plot_df, file = file.path(paths$intermediate_dir, "plot_df_logFC.rds"))

# -----------------------------
# 2) Cluster summary stats (mean, SE, CI) per diet x time_h
# -----------------------------
cluster_summary <- plot_df %>%
  dplyr::group_by(cluster, diet, time_h) %>%
  dplyr::summarise(
    n_genes    = dplyr::n_distinct(ProbeName),
    mean_logFC = mean(logFC, na.rm = TRUE),
    sd_logFC   = sd(logFC,   na.rm = TRUE),
    se         = sd_logFC / sqrt(n_genes),
    ci95       = 1.96 * se,
    .groups    = "drop"
  )

saveRDS(cluster_summary,
        file = file.path(paths$intermediate_dir, "cluster_summary.rds"))
readr::write_csv(cluster_summary,
                 file = file.path(paths$tables_dir, "cluster_summary.csv"))

message("Cluster sizes:")
print(cluster_summary %>%
        dplyr::distinct(cluster, diet, n_genes) %>%
        tidyr::pivot_wider(names_from = diet, values_from = n_genes))

# -----------------------------
# 3) Plot A: spaghetti (raw logFC) + mean line
# -----------------------------
spaghetti_plot <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = time_h, y = logFC)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = ProbeName, color = diet),
    alpha = 0.08, linewidth = 0.4
  ) +
  ggplot2::geom_line(
    data = cluster_summary,
    ggplot2::aes(y = mean_logFC, group = diet, color = diet),
    linewidth = 1.2
  ) +
  ggplot2::geom_point(
    data = cluster_summary,
    ggplot2::aes(y = mean_logFC, color = diet),
    size = 2
  ) +
  ggplot2::facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::labs(
    title = "Expression trajectories by cluster (raw logFC)",
    x     = "Time after infection (hours)",
    y     = "logFC",
    color = "Diet"
  )

out1 <- file.path(paths$figures_dir, "cluster_spaghetti_rawlogFC.png")
ggplot2::ggsave(out1, plot = spaghetti_plot, width = 14, height = 8, dpi = 300)
message("Saved: ", out1)

# -----------------------------
# 4) Plot B: mean ± 95% CI ribbons
# -----------------------------
ci_plot <- ggplot2::ggplot(
  cluster_summary,
  ggplot2::aes(x = time_h, y = mean_logFC,
               group = diet, color = diet, fill = diet)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_logFC - ci95, ymax = mean_logFC + ci95),
    alpha = 0.20, color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
  ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::scale_fill_manual(values  = diet_colors, labels = diet_labels) +
  ggplot2::theme_classic(base_size = 16) +
  ggplot2::labs(
    title = "Cluster mean trajectories \u00b1 95% CI (raw logFC)",
    x     = "Time after infection (hours)",
    y     = "Mean logFC",
    color = "Diet",
    fill  = "Diet"
  )

out2 <- file.path(paths$figures_dir, "cluster_mean_CI_rawlogFC.png")
ggplot2::ggsave(out2, plot = ci_plot, width = 14, height = 8, dpi = 300)
message("Saved: ", out2)

message("Saved cluster summary table to: ",
        file.path(paths$tables_dir, "cluster_summary.csv"))

# -----------------------------------------------------------------------
# ADDITION FOR R_10_cluster_summaries.R
# Purpose: Cluster-stratified Spearman logFC correlation between diets.
#
# Add this section at the END of your existing R_10 script, after your
# existing cluster trajectory summary code.
#
# Reads (already loaded in R_10, listed here for clarity):
#   - intermediate/config.rds          -> config
#   - intermediate/results_all.rds     -> results_all
#   - intermediate/cluster_df.rds      -> cluster_df
#   - intermediate/clustering_inputs.rds -> clustering_inputs (for logFC_mat)
#
# Writes:
#   - intermediate/spearman_cor_by_cluster_time.rds
#   - results/figures/logFC_spearman_cor_by_cluster.png
#   - results/figures/logFC_spearman_cor_cluster_vs_global.png
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# Load objects (remove any already loaded in your R_10 preamble)
# -----------------------------------------------------------------------
config            <- readRDS(file.path("intermediate", "config.rds"))
results_all       <- readRDS(file.path("intermediate", "results_all.rds"))
cluster_df        <- readRDS(file.path("intermediate", "cluster_df.rds"))
clustering_inputs <- readRDS(file.path("intermediate", "clustering_inputs.rds"))

paths             <- config$paths
time_levels_hours <- config$time_levels_hours
k                 <- config$clustering$k

diet_colors <- c("high" = "#F8766D", "low" = "#00BFC4")
diet_labels <- c("high" = "High resource", "low" = "Low resource")

# Cluster label lookup for plot facets
# Functional identities from GO enrichment (update if your labels differ)
cluster_labels <- c(
  "1" = "Cluster 1\n(Ribosomal / Translation)",
  "2" = "Cluster 2\n(No enrichment)",
  "3" = "Cluster 3\n(Cuticle / ECM)",
  "4" = "Cluster 4\n(Chromatin / Nucleosome)"
)

# -----------------------------------------------------------------------
# 1) Build a long-format logFC table joined to cluster assignments
# -----------------------------------------------------------------------
# We work from results_all rather than logFC_mat so we can use the
# diet and time_h columns that are already parsed there.
# Restrict to probes in the clustering universe.

lfc_clustered <- results_all %>%
  dplyr::filter(
    !is.na(diet),
    !is.na(time_h),
    !is.na(logFC),
    ProbeName %in% cluster_df$ProbeName
  ) %>%
  dplyr::group_by(ProbeName, diet, time_h) %>%
  dplyr::summarise(logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(cluster_df, by = "ProbeName")

# -----------------------------------------------------------------------
# 2) Compute Spearman correlation between diets at each cluster x timepoint
# -----------------------------------------------------------------------
# For each cluster and timepoint, pivot to wide (one column per diet)
# and compute Spearman r across all probes in that cluster.

spearman_by_cluster_time <- lfc_clustered %>%
  tidyr::pivot_wider(names_from = diet, values_from = logFC) %>%
  dplyr::filter(!is.na(high), !is.na(low)) %>%
  dplyr::group_by(cluster, time_h) %>%
  dplyr::summarise(
    n_probes   = dplyr::n(),
    spearman_r = suppressWarnings(
      stats::cor(high, low, method = "spearman", use = "pairwise.complete.obs")
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    time_h  = factor(time_h,  levels = time_levels_hours, ordered = TRUE),
    cluster = factor(cluster, levels = seq_len(k))
  )

message("Spearman logFC correlation by cluster x timepoint:")
print(spearman_by_cluster_time, n = Inf)

saveRDS(
  spearman_by_cluster_time,
  file = file.path(paths$intermediate_dir, "spearman_cor_by_cluster_time.rds")
)

# -----------------------------------------------------------------------
# 3) Also compute the global (all-probe) Spearman correlation per timepoint
#    to use as a reference line in the plot
# -----------------------------------------------------------------------
spearman_global <- results_all %>%
  dplyr::filter(!is.na(diet), !is.na(time_h), !is.na(logFC)) %>%
  dplyr::group_by(ProbeName, diet, time_h) %>%
  dplyr::summarise(logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = diet, values_from = logFC) %>%
  dplyr::filter(!is.na(high), !is.na(low)) %>%
  dplyr::group_by(time_h) %>%
  dplyr::summarise(
    spearman_r_global = suppressWarnings(
      stats::cor(high, low, method = "spearman", use = "pairwise.complete.obs")
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

# -----------------------------------------------------------------------
# 4) Plot: cluster-stratified Spearman correlation faceted by cluster
# -----------------------------------------------------------------------
# Add cluster label column for facet display
spearman_plot_df <- spearman_by_cluster_time %>%
  dplyr::mutate(
    cluster_label = dplyr::recode(as.character(cluster), !!!cluster_labels)
  )

# Join global correlation for reference line
spearman_plot_df <- spearman_plot_df %>%
  dplyr::left_join(spearman_global, by = "time_h")

cluster_spearman_plot <- ggplot2::ggplot(
  spearman_plot_df,
  ggplot2::aes(x = time_h, y = spearman_r, group = 1)
) +
  # Global reference line
  ggplot2::geom_line(
    ggplot2::aes(y = spearman_r_global),
    colour    = "grey60",
    linewidth = 0.7,
    linetype  = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth  = 0.5,
    colour     = "grey40",
    linetype   = "dotted"
  ) +
  ggplot2::geom_line(linewidth = 1.0, colour = "black") +
  ggplot2::geom_point(
    ggplot2::aes(size = n_probes),
    colour = "black"
  ) +
  ggplot2::facet_wrap(~ cluster_label, ncol = 2) +
  ggplot2::scale_size_continuous(
    name   = "N probes",
    range  = c(1.5, 4),
    breaks = c(100, 500, 1000, 2000)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      min(c(spearman_plot_df$spearman_r,
            spearman_plot_df$spearman_r_global), na.rm = TRUE) - 0.05,
      max(c(spearman_plot_df$spearman_r,
            spearman_plot_df$spearman_r_global), na.rm = TRUE) + 0.05
    )
  ) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title    = "Spearman logFC correlation between diets by cluster",
    subtitle = paste0(
      "Solid line = within-cluster correlation; ",
      "dashed grey = global all-probe correlation"
    ),
    x = "Time after infection (hours)",
    y = "Spearman r (high vs low resource logFC)"
  )

out_cluster_spearman <- file.path(
  paths$figures_dir, "logFC_spearman_cor_by_cluster.png"
)
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_cluster_spearman, width = 3200, height = 2400, res = 300)
  print(cluster_spearman_plot); grDevices::dev.off()
} else {
  grDevices::png(out_cluster_spearman, width = 3200, height = 2400, res = 300)
  print(cluster_spearman_plot); grDevices::dev.off()
}
message("Saved cluster-stratified Spearman plot to: ", out_cluster_spearman)

# -----------------------------------------------------------------------
# 5) Summary table: mean Spearman r per cluster averaged across timepoints
#    Split into early (<=24h) and late (>=48h) to highlight temporal patterns
# -----------------------------------------------------------------------
spearman_summary <- spearman_by_cluster_time %>%
  dplyr::mutate(
    phase = dplyr::case_when(
      as.numeric(as.character(time_h)) <= 24  ~ "Early (<=24h)",
      as.numeric(as.character(time_h)) >= 48  ~ "Late (>=48h)",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(phase)) %>%
  dplyr::group_by(cluster, phase) %>%
  dplyr::summarise(
    mean_spearman_r = round(mean(spearman_r, na.rm = TRUE), 3),
    n_timepoints    = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cluster, phase)

message("\nMean Spearman r by cluster and infection phase:")
print(spearman_summary, n = Inf)

# The key prediction to check:
# Cluster 3 (Cuticle/ECM - shared late program) should show higher
# between-diet correlation than Clusters 1 and 4 (diet-specific early programs),
# particularly at late timepoints. If this holds, it supports the argument
# that near-zero Jaccard similarity reflects genuine biological divergence
# in early programs rather than power differences alone.
