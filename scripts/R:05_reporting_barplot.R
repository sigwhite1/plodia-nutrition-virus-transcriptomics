# R/05_reporting_barplot.R
# Purpose: DEG count barplots, Up/Down stacked barplot, diet overlap analysis,
#          Jaccard similarity plot, and export of full results table.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/results_all.rds
#
# Writes:
#   - results/tables/limma_results_all_loopdesign.xlsx
#   - results/tables/limma_results_all_loopdesign.csv
#   - results/figures/number_of_DEGs_per_contrast.png
#   - results/figures/number_of_DEGs_per_contrast_up_down_stacked.png
#   - results/figures/DEG_overlap_high_vs_low.png
#   - results/figures/DEG_overlap_jaccard.png
#   - intermediate/sig_df2.rds
#   - intermediate/sig_dir_df2.rds
#   - intermediate/deg_overlap_high_low.rds
#   - intermediate/direction_by_time.rds

# -----------------------------
# 0) Load config + results
# -----------------------------
config      <- readRDS(file.path("intermediate", "config.rds"))
results_all <- readRDS(file.path("intermediate", "results_all.rds"))

time_levels_hours <- config$time_levels_hours   # numeric: 0.5 1 2 ... 168
fdr_cutoff        <- config$thresholds$fdr
paths             <- config$paths

# Diet label lookup for plot titles (lowercase keys match results_all$diet)
diet_labels <- c("low" = "Low resource", "high" = "High resource")
diet_colors <- c("low" = "#00BFC4", "high" = "#F8766D")

# Ensure output dirs exist (safe to call even when running standalone)
dir.create(paths$tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Export full results table
# -----------------------------
openxlsx::write.xlsx(
  results_all,
  file = file.path(paths$tables_dir, "limma_results_all_loopdesign.xlsx")
)
readr::write_csv(
  results_all,
  file = file.path(paths$tables_dir, "limma_results_all_loopdesign.csv")
)
message("Exported results table: limma_results_all_loopdesign (.xlsx + .csv)")

# -----------------------------
# 2) Count significant DEGs per diet x time
# -----------------------------
sig_df <- results_all %>%
  dplyr::filter(!is.na(diet), !is.na(time_h)) %>%
  dplyr::group_by(diet, time_h) %>%
  dplyr::summarise(
    DEGs = sum(adj.P.Val < fdr_cutoff, na.rm = TRUE),
    .groups = "drop"
  )

# Apply consistent factor levels — time_h is numeric so factor directly
sig_df <- sig_df %>%
  dplyr::mutate(
    diet   = factor(diet,   levels = c("high", "low")),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

# Fill missing diet x time combinations with 0
sig_df2 <- sig_df %>%
  tidyr::complete(time_h, diet, fill = list(DEGs = 0)) %>%
  dplyr::arrange(time_h, diet)

saveRDS(sig_df2, file = file.path(paths$intermediate_dir, "sig_df2.rds"))

# -----------------------------
# 3) DEG count barplot
# -----------------------------
degs_plot <- ggplot2::ggplot(
  sig_df2,
  ggplot2::aes(x = time_h, y = DEGs, fill = diet)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.9),
    width = 0.85
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = DEGs),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = -0.3,
    size = config$plot_cfg$deg_label_size
  ) +
  ggplot2::scale_fill_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title = paste0("Number of DEGs per Contrast (FDR < ", fdr_cutoff, ")"),
    x     = "Time after infection (hours)",
    y     = "Significant DEGs",
    fill  = "Diet"
  ) +
  ggplot2::ylim(0, max(sig_df2$DEGs, na.rm = TRUE) * 1.15)

print(degs_plot)

out_png <- file.path(paths$figures_dir, "number_of_DEGs_per_contrast.png")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_png, width = 3200, height = 1800, res = 300)
  print(degs_plot); grDevices::dev.off()
} else {
  grDevices::png(out_png, width = 3200, height = 1800, res = 300)
  print(degs_plot); grDevices::dev.off()
}
message("Saved DEG count plot to: ", out_png)

# -----------------------------
# 4) Up vs Down stacked barplot
# -----------------------------
sig_dir_df <- results_all %>%
  dplyr::filter(!is.na(diet), !is.na(time_h)) %>%
  dplyr::mutate(
    is_sig    = !is.na(adj.P.Val) & adj.P.Val < fdr_cutoff,
    Direction = dplyr::case_when(
      is_sig & logFC >  0 ~ "Up",
      is_sig & logFC <  0 ~ "Down",
      is_sig & logFC == 0 ~ "Zero",
      TRUE                ~ NA_character_
    )
  ) %>%
  dplyr::filter(Direction %in% c("Up", "Down")) %>%
  dplyr::group_by(diet, time_h, Direction) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

sig_dir_df <- sig_dir_df %>%
  dplyr::mutate(
    diet      = factor(diet,      levels = c("high", "low")),
    time_h    = factor(time_h,    levels = time_levels_hours, ordered = TRUE),
    Direction = factor(Direction, levels = c("Down", "Up"))
  )

sig_dir_df2 <- sig_dir_df %>%
  tidyr::complete(time_h, diet, Direction, fill = list(n = 0L)) %>%
  dplyr::arrange(time_h, diet, Direction)

saveRDS(sig_dir_df2, file = file.path(paths$intermediate_dir, "sig_dir_df2.rds"))

dir_alpha <- c("Down" = 0.45, "Up" = 1.0)

degs_stack_plot <- ggplot2::ggplot(
  sig_dir_df2,
  ggplot2::aes(x = time_h, y = n, fill = diet, alpha = Direction)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
    width = 0.85
  ) +
  ggplot2::geom_text(
    data = sig_dir_df2 %>%
      dplyr::group_by(diet, time_h) %>%
      dplyr::summarise(total = sum(n), .groups = "drop"),
    ggplot2::aes(x = time_h, y = total, label = total, group = diet),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = -0.3,
    size = config$plot_cfg$deg_label_size,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::scale_alpha_manual(values = dir_alpha) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title = paste0("DEGs per Contrast by Direction (FDR < ", fdr_cutoff, ")"),
    x     = "Time after infection (hours)",
    y     = "Significant DEGs",
    fill  = "Diet",
    alpha = "Direction"
  ) +
  ggplot2::ylim(0, max(
    sig_dir_df2 %>%
      dplyr::group_by(diet, time_h) %>%
      dplyr::summarise(total = sum(n), .groups = "drop") %>%
      dplyr::pull(total),
    na.rm = TRUE
  ) * 1.15)

print(degs_stack_plot)

out_png2 <- file.path(paths$figures_dir,
                      "number_of_DEGs_per_contrast_up_down_stacked.png")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_png2, width = 3200, height = 1800, res = 300)
  print(degs_stack_plot); grDevices::dev.off()
} else {
  grDevices::png(out_png2, width = 3200, height = 1800, res = 300)
  print(degs_stack_plot); grDevices::dev.off()
}
message("Saved Up/Down stacked plot to: ", out_png2)

# -----------------------------
# 5) Diet overlap analysis
# -----------------------------
deg_sets <- results_all %>%
  dplyr::filter(!is.na(diet), !is.na(time_h), adj.P.Val < fdr_cutoff) %>%
  dplyr::group_by(time_h, diet) %>%
  dplyr::summarise(genes = list(unique(ProbeName)), .groups = "drop") %>%
  dplyr::mutate(time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE))

deg_sets_wide <- deg_sets %>%
  tidyr::pivot_wider(names_from = diet, values_from = genes)

overlap_df <- deg_sets_wide %>%
  dplyr::mutate(
    high      = purrr::map(high, ~ if (is.null(.x)) character(0) else .x),
    low       = purrr::map(low,  ~ if (is.null(.x)) character(0) else .x),
    Shared    = purrr::map2(high, low,  intersect),
    High_only = purrr::map2(high, low,  setdiff),
    Low_only  = purrr::map2(low,  high, setdiff),
    n_shared     = purrr::map_int(Shared,    length),
    n_high_only  = purrr::map_int(High_only, length),
    n_low_only   = purrr::map_int(Low_only,  length),
    union_size   = purrr::map2_int(high, low, ~ length(union(.x, .y))),
    jaccard      = ifelse(union_size > 0, n_shared / union_size, 0)
  )

saveRDS(overlap_df,
        file = file.path(paths$intermediate_dir, "deg_overlap_high_low.rds"))

# Overall overlap across all timepoints
high_any <- results_all %>%
  dplyr::filter(diet == "high", adj.P.Val < fdr_cutoff) %>%
  dplyr::pull(ProbeName) %>% unique()

low_any <- results_all %>%
  dplyr::filter(diet == "low", adj.P.Val < fdr_cutoff) %>%
  dplyr::pull(ProbeName) %>% unique()

shared_any <- intersect(high_any, low_any)

message("Overall diet overlap (any timepoint):")
print(c(
  high_only_any = length(setdiff(high_any, low_any)),
  low_only_any  = length(setdiff(low_any,  high_any)),
  shared_any    = length(shared_any),
  jaccard_any   = length(shared_any) / length(union(high_any, low_any))
))

# -----------------------------
# 6) Direction concordance across timepoints
# -----------------------------
# For each probe x timepoint, compare the sign of logFC between diets
direction_by_time <- results_all %>%
  dplyr::group_by(time_h, ProbeName, diet) %>%
  dplyr::summarise(logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = diet, values_from = logFC) %>%
  dplyr::filter(!is.na(high), !is.na(low)) %>%
  dplyr::group_by(time_h) %>%
  dplyr::summarise(
    n_genes       = dplyr::n(),
    frac_same     = mean(sign(high) == sign(low)),
    frac_opposite = mean(sign(high) == -sign(low)),
    cor_logFC     = suppressWarnings(
      stats::cor(high, low, use = "pairwise.complete.obs")
    ),
    .groups = "drop"
  )

message("Direction concordance by timepoint:")
print(direction_by_time)

saveRDS(direction_by_time,
        file = file.path(paths$intermediate_dir, "direction_by_time.rds"))

# Correlation at 120h (5d equivalent) specifically
# 120h was previously referred to as "5d" — using numeric hours throughout
cor_120h <- results_all %>%
  dplyr::filter(time_h == 120, adj.P.Val < fdr_cutoff) %>%
  dplyr::pull(ProbeName) %>% unique()

if (length(cor_120h) > 1) {
  cor_120h_val <- results_all %>%
    dplyr::filter(time_h == 120, ProbeName %in% cor_120h) %>%
    dplyr::group_by(ProbeName, diet) %>%
    dplyr::summarise(logFC = median(logFC), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = diet, values_from = logFC) %>%
    dplyr::filter(!is.na(high), !is.na(low)) %>%
    dplyr::summarise(cor = cor(high, low)) %>%
    dplyr::pull(cor)
  message("logFC correlation between diets at 120h: ", round(cor_120h_val, 3))
} else {
  message("Too few shared DEGs at 120h for correlation.")
}

# -----------------------------
# 7) Stacked overlap plot
# -----------------------------
overlap_long <- overlap_df %>%
  dplyr::select(time_h, n_shared, n_high_only, n_low_only) %>%
  tidyr::pivot_longer(
    cols      = -time_h,
    names_to  = "Category",
    values_to = "Count"
  ) %>%
  dplyr::mutate(
    Category = factor(
      Category,
      levels = c("n_high_only", "n_shared", "n_low_only"),
      labels = c("High resource only", "Shared", "Low resource only")
    ),
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

overlap_colors <- c(
  "High resource only" = "#F8766D",
  "Low resource only"  = "#00BFC4",
  "Shared"             = "grey30"
)

overlap_plot <- ggplot2::ggplot(
  overlap_long,
  ggplot2::aes(x = time_h, y = Count, fill = Category)
) +
  ggplot2::scale_fill_manual(values = overlap_colors) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title = paste0("Overlap of DEGs Between Diets (FDR < ", fdr_cutoff, ")"),
    x     = "Time after infection (hours)",
    y     = "Number of DEGs",
    fill  = "Category"
  )

print(overlap_plot)

out_overlap_png <- file.path(paths$figures_dir, "DEG_overlap_high_vs_low.png")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_overlap_png, width = 3200, height = 1800, res = 300)
  print(overlap_plot); grDevices::dev.off()
} else {
  grDevices::png(out_overlap_png, width = 3200, height = 1800, res = 300)
  print(overlap_plot); grDevices::dev.off()
}
message("Saved overlap plot to: ", out_overlap_png)

# -----------------------------
# 8) Jaccard similarity plot
# -----------------------------
# Labels derived from data rather than hardcoded to avoid stale time references
nonzero_jaccard <- overlap_df %>% dplyr::filter(jaccard > 0)

jaccard_plot <- ggplot2::ggplot(
  overlap_df,
  ggplot2::aes(x = time_h, y = jaccard, group = 1)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, colour = "grey40") +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_point(
    data    = nonzero_jaccard,
    size    = 3,
    colour  = "#F8766D"
  ) +
  ggplot2::geom_text(
    data = nonzero_jaccard,
    ggplot2::aes(
      label = paste0(time_h, "h: ", signif(jaccard, 2))
    ),
    nudge_y = 0.00015,
    nudge_x = 0.2,
    size    = 4
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, max(overlap_df$jaccard, na.rm = TRUE) * 1.2),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title    = "Jaccard Similarity of DEGs Between Diets",
    subtitle = paste0(
      nrow(nonzero_jaccard), " timepoint(s) with shared DEGs: ",
      paste(paste0(nonzero_jaccard$time_h, "h"), collapse = ", ")
    ),
    x = "Time after infection (hours)",
    y = "Jaccard index (Shared / Union)"
  )

print(jaccard_plot)

out_jaccard_png <- file.path(paths$figures_dir, "DEG_overlap_jaccard.png")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_jaccard_png, width = 3200, height = 1800, res = 300)
  print(jaccard_plot); grDevices::dev.off()
} else {
  grDevices::png(out_jaccard_png, width = 3200, height = 1800, res = 300)
  print(jaccard_plot); grDevices::dev.off()
}
message("Saved Jaccard plot to: ", out_jaccard_png)

# ----------------------------- 
# 9) Spearman logFC correlation between diets across all probes
# ----------------------------- 
# Compute Spearman correlation between high and low diet logFC vectors
# at each timepoint, across ALL probes with non-missing values in both diets.
# This complements the Jaccard analysis by asking whether fold-change
# magnitudes and directions are concordant between diets even among
# non-significant genes, providing a direct assessment of whether
# near-zero DEG overlap reflects power differences or genuine divergence.

spearman_cor_by_time <- results_all %>%
  dplyr::filter(!is.na(diet), !is.na(time_h), !is.na(logFC)) %>%
  dplyr::group_by(time_h, ProbeName, diet) %>%
  dplyr::summarise(logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = diet, values_from = logFC) %>%
  dplyr::filter(!is.na(high), !is.na(low)) %>%
  dplyr::group_by(time_h) %>%
  dplyr::summarise(
    n_probes      = dplyr::n(),
    spearman_r    = suppressWarnings(
      stats::cor(high, low, method = "spearman", use = "pairwise.complete.obs")
    ),
    pearson_r     = suppressWarnings(
      stats::cor(high, low, method = "pearson",  use = "pairwise.complete.obs")
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    time_h = factor(time_h, levels = time_levels_hours, ordered = TRUE)
  )

message("Spearman logFC correlation between diets by timepoint:")
print(spearman_cor_by_time)

saveRDS(
  spearman_cor_by_time,
  file = file.path(paths$intermediate_dir, "spearman_cor_by_time.rds")
)

# Plot Spearman correlation across timepoints
spearman_plot <- ggplot2::ggplot(
  spearman_cor_by_time,
  ggplot2::aes(x = time_h, y = spearman_r, group = 1)
) +
  ggplot2::geom_hline(
    yintercept = 0, linewidth = 0.6, colour = "grey40", linetype = "dashed"
  ) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(spearman_r, 2)),
    vjust  = -0.8,
    size   = 3.5,
    color  = "#d6604d",
    fontface = "bold" # blue, visible against black points and white background
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      min(spearman_cor_by_time$spearman_r, na.rm = TRUE) * 1.2,
      max(spearman_cor_by_time$spearman_r, na.rm = TRUE) * 1.2
    ),
    expand = ggplot2::expansion(mult = c(0.05, 0.1))
  ) +
  ggplot2::theme_minimal(base_size = config$plot_cfg$base_size) +
  ggplot2::labs(
    title    = "Spearman Correlation of logFC Between Diets",
    subtitle = "Computed across all probes with non-missing values in both diets",
    x        = "Time after infection (hours)",
    y        = "Spearman r (high vs low resource logFC)"
  )

print(spearman_plot)

out_spearman_png <- file.path(paths$figures_dir, "logFC_spearman_cor_by_time.png")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_spearman_png, width = 3200, height = 1800, res = 300)
  print(spearman_plot); grDevices::dev.off()
} else {
  grDevices::png(out_spearman_png, width = 3200, height = 1800, res = 300)
  print(spearman_plot); grDevices::dev.off()
}
message("Saved Spearman correlation plot to: ", out_spearman_png)
