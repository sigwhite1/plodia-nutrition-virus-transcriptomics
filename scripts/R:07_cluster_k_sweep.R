# R/07_cluster_k_sweep.R
# Purpose:
#   1) Build logFC matrix, scale rows, run hierarchical clustering (Ward.D2)
#   2) Sweep k = 2:5 — compute silhouette + within-cluster SSE (elbow)
#   3) Auto-generate trajectory plots for each k (spaghetti+mean AND mean±95% CI)
#   4) Bootstrap stability assessment for the chosen k (config$clustering$k)
#      to confirm cluster assignments are robust to gene resampling
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/results_all.rds
#   - intermediate/contrasts.rds
#   - intermediate/gene_sets.rds
#
# Writes:
#   - intermediate/clustering_inputs.rds     (logFC_mat, scaled_mat, hc — reused by R_08)
#   - results/tables/k_sweep_metrics.csv
#   - intermediate/k_sweep_assignments.rds
#   - results/figures/k_sweep/               (silhouette + elbow + per-k trajectory plots)
#   - intermediate/bootstrap_co_matrix_k4.rds
#   - intermediate/cluster_stability_k4.rds
#   - results/figures/k_sweep/bootstrap_stability_heatmap_k4.png

# -----------------------------
# 0) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(cluster)   # silhouette()
})

# -----------------------------
# 1) Load artifacts
# -----------------------------
config <- readRDS(file.path("intermediate", "config.rds"))
paths             <- config$paths
time_levels_hours <- config$time_levels_hours   # numeric: 0.5 1 2 ... 168
k_range           <- 2:5                        # sweep range — edit here if needed

results_all   <- readRDS(file.path(paths$intermediate_dir, "results_all.rds"))
contrasts_obj <- readRDS(file.path(paths$intermediate_dir, "contrasts.rds"))
gene_sets     <- readRDS(file.path(paths$intermediate_dir, "gene_sets.rds"))

dir.create(paths$figures_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(paths$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- file.path(paths$figures_dir, "k_sweep")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) Build logFC matrix and scaled matrix from scratch
# -----------------------------
# R_07 runs before R_08 so clustering_inputs.rds does not exist yet.
# We build the matrices here from results_all + gene_sets, then save
# clustering_inputs.rds so R_08 can reuse them without rebuilding.

deg_probes <- gene_sets$all_deg_probes
message("DEG universe for sweep: ", length(deg_probes), " probes")

lfc_long <- results_all %>%
  dplyr::select(ProbeName, contrast, logFC) %>%
  dplyr::filter(ProbeName %in% deg_probes) %>%
  dplyr::group_by(ProbeName, contrast) %>%
  dplyr::summarise(logFC = dplyr::first(logFC), .groups = "drop")

lfc_wide <- lfc_long %>%
  tidyr::pivot_wider(names_from = contrast, values_from = logFC)

logFC_mat <- lfc_wide %>%
  tibble::column_to_rownames("ProbeName") %>%
  as.matrix()

# Impute NAs
na_count <- sum(is.na(logFC_mat))
if (na_count > 0) {
  message("Imputing ", na_count, " NA values in logFC matrix to 0.")
  logFC_mat[is.na(logFC_mat)] <- 0
}

# Row-wise scaling
scaled_mat <- t(scale(t(logFC_mat)))
n_zero_var <- sum(is.na(rowSums(scaled_mat)))
if (n_zero_var > 0) {
  message(n_zero_var, " probes have zero variance — set to 0.")
  scaled_mat[is.na(scaled_mat)] <- 0
}

message("logFC matrix: ", nrow(logFC_mat), " probes x ", ncol(logFC_mat), " contrasts")

# -----------------------------
# 3) Contrast metadata for ordering and plotting
# -----------------------------
# contrast_key columns: diet, time_h, time_label, contrast_name
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

# Ensure columns ordered by time then diet
ordered_contrasts <- contrast_key %>%
  dplyr::arrange(time_h, diet) %>%
  dplyr::pull(contrast_name)
ordered_contrasts <- ordered_contrasts[ordered_contrasts %in% colnames(logFC_mat)]

if (length(ordered_contrasts) >= 2) {
  logFC_mat  <- logFC_mat[,  ordered_contrasts, drop = FALSE]
  scaled_mat <- scaled_mat[, ordered_contrasts, drop = FALSE]
}

# -----------------------------
# 4) Build hierarchical clustering (from scratch)
# -----------------------------
message("Building hierarchical clustering (Ward.D2, euclidean)...")
d  <- dist(scaled_mat, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

# Save clustering_inputs.rds so R_08 can reuse matrices + hc
# without rebuilding. R_08 will overwrite this with the final chosen k.
clustering_inputs_pre <- list(
  logFC_mat  = logFC_mat,
  scaled_mat = scaled_mat,
  hc         = hc
)
saveRDS(clustering_inputs_pre,
        file.path(paths$intermediate_dir, "clustering_inputs.rds"))
message("Saved clustering_inputs.rds (pre-sweep) for reuse by R_08.")

# Distance matrix for silhouette (computed once outside loop)
dmat <- dist(scaled_mat, method = "euclidean")

# -----------------------------
# 5) Metric helpers + cluster label stabilisation
# -----------------------------

# Anchor cluster labels to mean logFC at the strongest contrast so that
# label assignments are consistent across k values and re-runs.
# Cluster 1 = most downregulated, cluster k = most upregulated at anchor.
stabilise_cluster_labels <- function(clusters, logFC_mat,
                                     anchor_contrast = "high.1h.virus_vs_ctrl") {
  if (!anchor_contrast %in% colnames(logFC_mat)) {
    warning("Anchor contrast '", anchor_contrast,
            "' not found — cluster labels left as-is.")
    return(as.integer(clusters))
  }
  anchor_means <- tapply(logFC_mat[, anchor_contrast],
                         clusters, mean, na.rm = TRUE)
  new_order    <- rank(anchor_means, ties.method = "first")
  as.integer(new_order[as.character(clusters)])
}

within_sse <- function(mat, clusters) {
  cl  <- as.integer(clusters)
  sse <- 0
  for (ki in sort(unique(cl))) {
    idx <- which(cl == ki)
    if (length(idx) <= 1) next
    sub <- mat[idx, , drop = FALSE]
    cen <- colMeans(sub)
    sse <- sse + sum(rowSums(
      (sub - matrix(cen, nrow = nrow(sub), ncol = ncol(sub), byrow = TRUE))^2
    ))
  }
  sse
}

mean_sil <- function(clusters, d) {
  sil <- cluster::silhouette(as.integer(clusters), d)
  mean(sil[, 3])
}

# -----------------------------
# 6) Prepare base long table for plotting (raw logFC, built once)
# -----------------------------
plot_df_base <- results_all %>%
  dplyr::filter(ProbeName %in% deg_probes) %>%
  dplyr::select(ProbeName, contrast, logFC) %>%       # 'contrast' is lowercase
  dplyr::left_join(
    contrast_key,
    by = c("contrast" = "contrast_name")
  ) %>%
  dplyr::filter(!is.na(diet), !is.na(time_h)) %>%
  dplyr::mutate(
    diet   = factor(diet,   levels = c("low", "high")),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

# -----------------------------
# 7) Sweep k
# -----------------------------
metrics     <- list()
assignments <- list()

for (k in k_range) {
  message("Processing k = ", k)
  
  clusters <- as.integer(cutree(hc, k = k))
  clusters <- stabilise_cluster_labels(clusters, logFC_mat)
  # Note: rownames(scaled_mat) used for ProbeName since stabilise returns
  # an unnamed vector. Row order of scaled_mat is preserved by cutree.
  
  cl_df <- tibble::tibble(
    ProbeName = rownames(scaled_mat),
    cluster   = factor(clusters, levels = sort(unique(clusters)))
  )
  assignments[[paste0("k", k)]] <- cl_df
  
  # Metrics
  sse <- within_sse(scaled_mat, clusters)
  sil <- mean_sil(clusters, dmat)
  
  metrics[[paste0("k", k)]] <- tibble::tibble(
    k               = k,
    n_clusters      = length(unique(clusters)),
    within_sse      = sse,
    mean_silhouette = sil
  )
  
  message("  SSE = ", round(sse, 1), "  |  Mean silhouette = ", round(sil, 4))
  message("  Cluster sizes: ", paste(table(clusters), collapse = " / "))
  
  # Join cluster assignments onto plot base
  plot_df <- plot_df_base %>%
    dplyr::left_join(cl_df, by = "ProbeName") %>%
    dplyr::filter(!is.na(cluster))
  
  # Cluster summary: mean ± CI per cluster x diet x time
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
  
  # ---- Plot 1: spaghetti + mean (raw logFC)
  p_spaghetti <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = time_h, y = logFC)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = ProbeName, color = diet),
      alpha = 0.06, linewidth = 0.35
    ) +
    ggplot2::geom_line(
      data = cluster_summary,
      ggplot2::aes(y = mean_logFC, group = diet, color = diet),
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = cluster_summary,
      ggplot2::aes(y = mean_logFC, color = diet),
      size = 1.6
    ) +
    ggplot2::facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::labs(
      title = paste0("Spaghetti + mean trajectories (raw logFC) \u2014 k = ", k),
      x     = "Time after infection (hours)",
      y     = "logFC",
      color = "Diet"
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("k", k, "_spaghetti_mean_rawlogFC.png")),
    plot     = p_spaghetti,
    width = 14, height = 8, dpi = 300
  )
  
  # ---- Plot 2: mean ± 95% CI (raw logFC)
  p_ci <- ggplot2::ggplot(
    cluster_summary,
    ggplot2::aes(x = time_h, y = mean_logFC,
                 group = diet, color = diet, fill = diet)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_logFC - ci95, ymax = mean_logFC + ci95),
      alpha = 0.18, color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
    ggplot2::scale_fill_manual(values  = diet_colors, labels = diet_labels) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(
      title = paste0("Cluster mean trajectories \u00b1 95% CI (raw logFC) \u2014 k = ", k),
      x     = "Time after infection (hours)",
      y     = "Mean logFC",
      color = "Diet",
      fill  = "Diet"
    )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("k", k, "_mean_CI_rawlogFC.png")),
    plot     = p_ci,
    width = 14, height = 8, dpi = 300
  )
}

# -----------------------------
# 8) Save metrics and assignments
# -----------------------------
metrics_df <- dplyr::bind_rows(metrics)

readr::write_csv(metrics_df,
                 file.path(paths$tables_dir, "k_sweep_metrics.csv"))
saveRDS(assignments,
        file.path(paths$intermediate_dir, "k_sweep_assignments.rds"))

message("Saved metrics to:     ",
        file.path(paths$tables_dir, "k_sweep_metrics.csv"))
message("Saved assignments to: ",
        file.path(paths$intermediate_dir, "k_sweep_assignments.rds"))

# Print metrics table for immediate inspection
message("k-sweep metrics:")
print(metrics_df)

# -----------------------------
# 9) Silhouette and elbow plots
# -----------------------------
p_sil <- ggplot2::ggplot(
  metrics_df,
  ggplot2::aes(x = k, y = mean_silhouette)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::labs(
    title = "Mean silhouette width vs k",
    x     = "k (number of clusters)",
    y     = "Mean silhouette width"
  )

ggplot2::ggsave(
  file.path(out_dir, "k_sweep_mean_silhouette.png"),
  plot = p_sil, width = 8, height = 5, dpi = 300
)

p_elbow <- ggplot2::ggplot(
  metrics_df,
  ggplot2::aes(x = k, y = within_sse)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::labs(
    title = "Elbow plot (within-cluster SSE) vs k",
    x     = "k (number of clusters)",
    y     = "Within-cluster SSE (scaled profiles)"
  )

ggplot2::ggsave(
  file.path(out_dir, "k_sweep_elbow_withinSSE.png"),
  plot = p_elbow, width = 8, height = 5, dpi = 300
)

message("Saved k-sweep figures to: ", out_dir)

# =============================================================================
# 10) Bootstrap stability for chosen k
# =============================================================================
# Assess whether the k=4 cluster assignments are reproducible by resampling
# genes with replacement and re-clustering 100 times. For each bootstrap
# replicate, we record which genes are assigned to the same cluster
# (co-clustering frequency). High co-clustering within the original clusters
# indicates stable, well-separated groups.
#
# NOTE ON RUNTIME: with n_genes=6595 and n_boot=100, the co-clustering matrix
# is 6595x6595 (~350MB). The heatmap step converts this to a long-format data
# frame (43.5M rows) which is memory-intensive. A summary heatmap based on
# cluster-level means is used instead of the full gene x gene matrix.
#
# n_boot = 100 is appropriate for exploratory analysis.

set.seed(123)

k_chosen  <- config$clustering$k          # 4
n_boot    <- 1000
n_genes   <- nrow(scaled_mat)
gene_ids  <- rownames(scaled_mat)

message("\nRunning bootstrap stability (k = ", k_chosen,
        ", n_boot = ", n_boot, ", n_genes = ", n_genes, ")...")

# Retrieve original cluster assignments (from the full-data hclust in section 4)
original_clusters   <- as.integer(cutree(hc, k = k_chosen))
original_clusters   <- stabilise_cluster_labels(original_clusters, logFC_mat)
# gene_ids = rownames(scaled_mat), preserved since stabilise returns unnamed vector
original_cluster_df <- tibble::tibble(
  ProbeName = gene_ids,
  cluster   = original_clusters
)

# Accumulate co-clustering counts
co_matrix <- matrix(0L, n_genes, n_genes,
                    dimnames = list(gene_ids, gene_ids))

for (b in seq_len(n_boot)) {
  if (b %% 25 == 0) message("  Bootstrap replicate ", b, " / ", n_boot)
  
  boot_idx     <- sample(n_genes, replace = TRUE)
  boot_mat     <- scaled_mat[boot_idx, , drop = FALSE]
  boot_hc      <- hclust(dist(boot_mat, method = "euclidean"),
                         method = "ward.D2")
  boot_cl      <- cutree(boot_hc, k = k_chosen)
  boot_gene_ids <- gene_ids[boot_idx]
  
  for (cl in unique(boot_cl)) {
    members <- unique(boot_gene_ids[boot_cl == cl])
    co_matrix[members, members] <- co_matrix[members, members] + 1L
  }
}

# Normalise to co-clustering frequency [0, 1]
co_matrix_freq <- co_matrix / n_boot

saveRDS(co_matrix_freq,
        file.path(paths$intermediate_dir, "bootstrap_co_matrix_k4.rds"))
message("Saved co-clustering matrix to: ",
        file.path(paths$intermediate_dir, "bootstrap_co_matrix_k4.rds"))

# -----------------------------
# Cluster-wise stability scores
# -----------------------------
stability_scores <- purrr::map_dfr(
  sort(unique(original_cluster_df$cluster)),
  function(cl) {
    members <- original_cluster_df$ProbeName[original_cluster_df$cluster == cl]
    submat  <- co_matrix_freq[members, members]
    tibble::tibble(
      cluster         = cl,
      mean_stability  = mean(submat[upper.tri(submat)]),
      size            = length(members)
    )
  }
)

message("\nCluster stability scores (k = ", k_chosen, "):")
print(stability_scores)

saveRDS(stability_scores,
        file.path(paths$intermediate_dir, "cluster_stability_k4.rds"))
readr::write_csv(stability_scores,
                 file.path(paths$tables_dir, "cluster_stability_k4.csv"))

# -----------------------------
# Stability summary heatmap
# (cluster-mean co-clustering, not full gene x gene — avoids memory crash)
# -----------------------------
# Build a k x k matrix of mean co-clustering between cluster pairs
k_labels <- sort(unique(original_cluster_df$cluster))
mean_co  <- matrix(NA_real_, length(k_labels), length(k_labels),
                   dimnames = list(paste0("C", k_labels),
                                   paste0("C", k_labels)))

for (i in seq_along(k_labels)) {
  for (j in seq_along(k_labels)) {
    mi <- original_cluster_df$ProbeName[original_cluster_df$cluster == k_labels[i]]
    mj <- original_cluster_df$ProbeName[original_cluster_df$cluster == k_labels[j]]
    mean_co[i, j] <- mean(co_matrix_freq[mi, mj])
  }
}

mean_co_long <- as.data.frame(as.table(mean_co)) %>%
  dplyr::rename(Cluster_i = Var1, Cluster_j = Var2, mean_co_freq = Freq)

p_stability <- ggplot2::ggplot(
  mean_co_long,
  ggplot2::aes(x = Cluster_i, y = Cluster_j, fill = mean_co_freq)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = round(mean_co_freq, 2)),
                     size = 5, colour = "white", fontface = "bold") +
  ggplot2::scale_fill_viridis_c(
    name   = "Mean co-clustering\nfrequency",
    limits = c(0, 1),
    option = "D"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::labs(
    title    = paste0("Bootstrap co-clustering stability (k = ", k_chosen,
                      ", n_boot = ", n_boot, ")"),
    subtitle = "Diagonal = within-cluster stability; off-diagonal = between-cluster overlap",
    x = NULL, y = NULL
  )

out_stability <- file.path(out_dir, "bootstrap_stability_heatmap_k4.png")
ggplot2::ggsave(out_stability, plot = p_stability,
                width = 6, height = 5, dpi = 300)
message("Saved stability heatmap to: ", out_stability)
