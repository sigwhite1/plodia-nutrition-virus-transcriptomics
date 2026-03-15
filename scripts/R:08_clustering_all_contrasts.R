# R/08_clustering_all_contrasts.R
# Purpose: Cluster DEGs across ALL diet x time contrasts using scaled logFC
#          patterns. Uses the union of DEGs across all contrasts and cuts
#          the dendrogram at k (from config$clustering$k, default 4).
#
# K SELECTION RATIONALE (from R/07_cluster_k_sweep.R):
#   k=2: highest silhouette (0.189) but biologically too coarse — separates
#        only "strong responders" vs "weak responders"
#   k=3: lowest silhouette (0.143); produces an unstable cluster with an
#        implausible trajectory (deep downregulation reversing to upregulation)
#   k=4: silhouette recovers (0.153); yields four biologically interpretable
#        patterns — early responders, diet-specific suppressors, broadly
#        suppressed genes, and oscillating responders. Chosen.
#   k=5: silhouette flat (0.151); fragments the early-responder cluster from
#        k=4 without revealing new biological patterns
#   config$clustering$k is set to 4 in R/00_config.R.
#
# CLUSTER LABEL STABILITY:
#   Labels are anchored by mean logFC at high.1h.virus_vs_ctrl (the strongest
#   contrast in this dataset) so cluster 1 = most downregulated and cluster k =
#   most upregulated at that timepoint. This ensures labels are consistent
#   across re-runs even if the dendrogram structure is unchanged.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/results_all.rds
#   - intermediate/contrasts.rds
#   - intermediate/gene_sets.rds
#
# Writes:
#   - intermediate/clustering_inputs.rds    (logFC_mat, scaled_mat, hc,
#                                             cluster_df, contrast_key_ordered)
#   - intermediate/cluster_df.rds           (ProbeName -> cluster)
#   - results/figures/cluster_spaghetti_mean_profiles.png
#   - results/figures/cluster_heatmap_ggplot2.png

# -----------------------------
# 0) Load artifacts
# -----------------------------
config        <- readRDS(file.path("intermediate", "config.rds"))
results_all   <- readRDS(file.path("intermediate", "results_all.rds"))
contrasts_obj <- readRDS(file.path("intermediate", "contrasts.rds"))
gene_sets     <- readRDS(file.path("intermediate", "gene_sets.rds"))

paths             <- config$paths
time_levels_hours <- config$time_levels_hours   # numeric: 0.5 1 2 ... 168
fdr_cutoff        <- config$thresholds$fdr
k                 <- config$clustering$k        # from config (default 4)

dir.create(paths$figures_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(paths$intermediate_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Prepare contrast metadata for ordering and plotting
# -----------------------------
# contrast_key columns: diet, time_h, time_label, contrast_name, is_estimable
contrast_key <- contrasts_obj$contrast_key %>%
  dplyr::filter(is_estimable) %>%
  dplyr::select(diet, time_h, contrast_name) %>%
  dplyr::mutate(
    diet   = factor(diet,   levels = c("low", "high")),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

# Diet color palette (keys match results_all$diet)
diet_colors <- c("high" = "#F8766D", "low" = "#00BFC4")
diet_labels <- c("high" = "High resource", "low" = "Low resource")

# -----------------------------
# 2) Choose DEG universe to cluster
# -----------------------------
deg_probes <- gene_sets$all_deg_probes

message("Clustering DEG universe: ", length(deg_probes),
        " probes (union across all contrasts at FDR < ", fdr_cutoff, ")")

if (length(deg_probes) < 2) {
  stop("Too few DEGs to cluster (n < 2). Check fdr_cutoff or upstream results.")
}

# -----------------------------
# 3) Build logFC matrix: rows = ProbeName, cols = contrast
# -----------------------------
lfc_long <- results_all %>%
  dplyr::select(ProbeName, contrast, logFC) %>%   # 'contrast' is lowercase
  dplyr::filter(ProbeName %in% deg_probes)

# Guard against duplicate ProbeName x contrast rows
lfc_long <- lfc_long %>%
  dplyr::group_by(ProbeName, contrast) %>%
  dplyr::summarise(logFC = dplyr::first(logFC), .groups = "drop")

lfc_wide <- lfc_long %>%
  tidyr::pivot_wider(names_from = contrast, values_from = logFC)

logFC_mat <- lfc_wide %>%
  tibble::column_to_rownames("ProbeName") %>%
  as.matrix()

# -----------------------------
# 4) Order columns by time then diet
# -----------------------------
ordered_contrasts <- contrast_key %>%
  dplyr::arrange(time_h, diet) %>%
  dplyr::pull(contrast_name)

ordered_contrasts <- ordered_contrasts[ordered_contrasts %in% colnames(logFC_mat)]

if (length(ordered_contrasts) >= 2) {
  logFC_mat <- logFC_mat[, ordered_contrasts, drop = FALSE]
}

message("logFC matrix: ", nrow(logFC_mat), " probes x ", ncol(logFC_mat), " contrasts")

# -----------------------------
# 5) Handle missing values
# -----------------------------
na_count <- sum(is.na(logFC_mat))
if (na_count > 0) {
  message("Imputing ", na_count, " NA values in logFC matrix to 0.")
  logFC_mat[is.na(logFC_mat)] <- 0
}

# -----------------------------
# 6) Row-wise scaling for pattern clustering
# -----------------------------
# Centre and scale each gene across contrasts so clustering captures
# expression pattern shape rather than magnitude
scaled_mat <- t(scale(t(logFC_mat)))
# Genes with zero variance across all contrasts get NA after scaling -> set to 0
n_zero_var <- sum(is.na(rowSums(scaled_mat)))
if (n_zero_var > 0) {
  message(n_zero_var, " probes have zero variance across contrasts — set to 0.")
  scaled_mat[is.na(scaled_mat)] <- 0
}

# -----------------------------
# 7) Hierarchical clustering
# -----------------------------
d        <- stats::dist(scaled_mat, method = "euclidean")
hc       <- stats::hclust(d, method = "ward.D2")
clusters <- stats::cutree(hc, k = k)

# Stabilise cluster labels across re-runs by anchoring to a biological anchor:
# rank clusters by mean logFC at the strongest contrast (high.1h.virus_vs_ctrl).
# Cluster 1 = most downregulated, cluster k = most upregulated at that timepoint.
# This prevents label permutation if data or pipeline parameters change.
stabilise_cluster_labels <- function(clusters, logFC_mat,
                                     anchor_contrast = "high.1h.virus_vs_ctrl") {
  if (!anchor_contrast %in% colnames(logFC_mat)) {
    warning("Anchor contrast '", anchor_contrast,
            "' not found in logFC_mat — cluster labels left as-is.")
    return(as.integer(clusters))
  }
  anchor_means <- tapply(logFC_mat[, anchor_contrast],
                         clusters, mean, na.rm = TRUE)
  new_order    <- rank(anchor_means, ties.method = "first")
  as.integer(new_order[as.character(clusters)])
}

# Capture probe names before stabilising (stabilise returns unnamed integer vector)
probe_names <- names(clusters)
clusters    <- stabilise_cluster_labels(clusters, logFC_mat)

cluster_df <- tibble::tibble(
  ProbeName = probe_names,
  cluster   = clusters
)

message("Cluster sizes (k = ", k, ") — labels anchored to high.1h.virus_vs_ctrl:")
print(table(cluster_df$cluster))

# -----------------------------
# 8) Save clustering objects
# -----------------------------
saveRDS(cluster_df, file = file.path(paths$intermediate_dir, "cluster_df.rds"))

clustering_inputs <- list(
  logFC_mat            = logFC_mat,
  scaled_mat           = scaled_mat,
  hc                   = hc,
  cluster_df           = cluster_df,
  k                    = k,
  fdr_cutoff           = fdr_cutoff,
  contrast_key_ordered = contrast_key %>%
    dplyr::arrange(time_h, diet)
)

saveRDS(clustering_inputs,
        file = file.path(paths$intermediate_dir, "clustering_inputs.rds"))

message("Saved cluster_df to:        ",
        file.path(paths$intermediate_dir, "cluster_df.rds"))
message("Saved clustering_inputs to: ",
        file.path(paths$intermediate_dir, "clustering_inputs.rds"))

# -----------------------------
# 9) Spaghetti + mean trajectory plot
# -----------------------------
scaled_long <- as.data.frame(scaled_mat) %>%
  tibble::rownames_to_column("ProbeName") %>%
  tidyr::pivot_longer(-ProbeName,
                      names_to  = "contrast",
                      values_to = "scaled_logFC") %>%
  dplyr::left_join(cluster_df, by = "ProbeName") %>%
  dplyr::left_join(
    contrast_key %>% dplyr::select(diet, time_h, contrast_name),
    by = c("contrast" = "contrast_name")
  ) %>%
  dplyr::mutate(
    time_h  = factor(time_h,  levels = time_levels_hours, ordered = TRUE),
    diet    = factor(diet,    levels = c("low", "high")),
    cluster = factor(cluster)
  ) %>%
  dplyr::filter(!is.na(diet), !is.na(time_h))

mean_df <- scaled_long %>%
  dplyr::group_by(cluster, diet, time_h) %>%
  dplyr::summarise(mean_scaled_logFC = mean(scaled_logFC, na.rm = TRUE),
                   .groups = "drop")

traj_plot <- ggplot2::ggplot(
  scaled_long,
  ggplot2::aes(x = time_h, y = scaled_logFC)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = ProbeName, color = diet),
    alpha = 0.08, linewidth = 0.4
  ) +
  ggplot2::geom_line(
    data = mean_df,
    ggplot2::aes(y = mean_scaled_logFC, group = diet, color = diet),
    linewidth = 1.2
  ) +
  ggplot2::geom_point(
    data = mean_df,
    ggplot2::aes(y = mean_scaled_logFC, color = diet),
    size = 2
  ) +
  ggplot2::facet_wrap(~ cluster, ncol = 2) +
  ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title = paste0("Cluster trajectories across diet \u00d7 time (k = ", k, ")"),
    x     = "Time after infection (hours)",
    y     = "Scaled logFC",
    color = "Diet"
  )

profile_file <- file.path(paths$figures_dir,
                          "cluster_spaghetti_mean_profiles.png")
ggplot2::ggsave(profile_file, plot = traj_plot,
                width = 14, height = 8, dpi = 300)
message("Saved spaghetti+mean trajectories to: ", profile_file)

# -----------------------------
# 10) Scaled logFC heatmap
# -----------------------------
# Order rows by cluster then dendrogram order within each cluster
hc_order    <- rownames(scaled_mat)[hc$order]
cluster_vec <- cluster_df$cluster[match(hc_order, cluster_df$ProbeName)]
row_order   <- hc_order[order(cluster_vec, seq_along(cluster_vec))]

heat_long <- scaled_mat[row_order, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ProbeName") %>%
  tidyr::pivot_longer(-ProbeName,
                      names_to  = "contrast",
                      values_to = "scaled_logFC") %>%
  dplyr::left_join(cluster_df, by = "ProbeName") %>%
  dplyr::left_join(
    contrast_key %>% dplyr::select(diet, time_h, contrast_name),
    by = c("contrast" = "contrast_name")
  ) %>%
  dplyr::mutate(
    time_h  = factor(time_h,  levels = time_levels_hours, ordered = TRUE),
    diet    = factor(diet,    levels = c("low", "high")),
    cluster = factor(cluster)
  )

# Stable x-axis ordering: time then diet
contrast_levels <- contrast_key %>%
  dplyr::arrange(time_h, diet) %>%
  dplyr::pull(contrast_name)

heat_long$contrast  <- factor(heat_long$contrast,  levels = contrast_levels)
heat_long$ProbeName <- factor(heat_long$ProbeName, levels = row_order)

heat_plot <- ggplot2::ggplot(
  heat_long,
  ggplot2::aes(x = contrast, y = ProbeName, fill = scaled_logFC)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low  = "#00BFC4", mid = "white", high = "#F8766D", midpoint = 0,
    name = "Scaled\nlogFC"
  ) +
  ggplot2::facet_wrap(~ cluster, scales = "free_y", ncol = 2) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid   = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = paste0("Scaled logFC heatmap by cluster (k = ", k, ")"),
    x     = "Contrast (time \u00d7 diet)",
    y     = "Genes"
  )

heatmap_file <- file.path(paths$figures_dir, "cluster_heatmap_ggplot2.png")
ggplot2::ggsave(heatmap_file, plot = heat_plot,
                width = 14, height = 10, dpi = 300)
message("Saved heatmap to: ", heatmap_file)
