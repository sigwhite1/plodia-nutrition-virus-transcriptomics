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

