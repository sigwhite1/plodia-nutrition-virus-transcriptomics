# =========================================================
# 01_import_and_qc.R
#
# Purpose:
# Import Agilent two-color microarray data, perform normalization and
# array-level QC, remove flagged low-quality arrays, and save cleaned
# objects for downstream limma analysis.
#
# Inputs:
# - data/raw/E-MTAB-5868/
# - data/metadata/pigv_infection_timeseries_diet_microarray_slides.txt
#
# Outputs:
# - results/rds/MA_filtered.rds
# - results/rds/targets_filtered.rds
# - results/rds/cor_A_filtered.rds
# - results/rds/cor_M_filtered.rds
# - results/tables/qc_mean_array_correlations.csv
# - results/tables/qc_removed_arrays.csv
# - results/microarray_QC_clean.RData   # legacy output for downstream code
# =========================================================

# ----------------------------
# 1. Load packages
# ----------------------------
library(limma)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ----------------------------
# 2. Define file paths
# ----------------------------
# =========================================================
# Local paths (edit these for your machine)
# =========================================================

raw_dir <- "data/raw/E-MTAB-5868"
targets_file <- "data/metadata/pigv_infection_timeseries_diet_microarray_slides.txt"

# Output directories
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("results/rds", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/qc", showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 3. Read targets / slide design file
# ----------------------------
targets <- read.delim(
  targets_file,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(targets)[1] <- "FileName"

# Adjust file names to match unzipped .txt files if necessary
targets$FileName <- sub("\\.txt\\.gz$", ".txt", targets$FileName)

# Full paths to Agilent files
files_full <- file.path(raw_dir, targets$FileName)

# Check for missing files
missing_files <- targets$FileName[!file.exists(files_full)]
if (length(missing_files) > 0) {
  message("Missing files:")
  print(missing_files)
}

# Keep only files that exist
keep_existing <- file.exists(files_full)
targets <- targets[keep_existing, , drop = FALSE]
files_full <- files_full[keep_existing]

# ----------------------------
# 4. Import and normalize arrays
# ----------------------------
RG <- read.maimages(files_full, source = "agilent", green.only = FALSE)

# Background correction
RG <- backgroundCorrect(RG, method = "normexp", offset = 50)

# Within-array normalization (red/green bias correction)
MA <- normalizeWithinArrays(RG, method = "loess")

# Between-array normalization
MA <- normalizeBetweenArrays(MA, method = "Aquantile")

# Store targets explicitly in MA object for clarity
MA$targets <- targets

# Store probe identifiers as row names for downstream limma + annotation mapping
if (!is.null(MA$genes) && "ProbeName" %in% colnames(MA$genes)) {
  rownames(MA$M) <- MA$genes$ProbeName
  rownames(MA$A) <- MA$genes$ProbeName
}


# ----------------------------
# 5. Quick QC summaries
# ----------------------------

# A-value correlation between arrays
cor_A <- cor(MA$A, use = "pairwise.complete.obs")

# M-value correlation between arrays
cor_M <- cor(MA$M, use = "pairwise.complete.obs")

# Mean A-value correlation per array
mean_corr_A <- rowMeans(cor_A, na.rm = TRUE)

qc_summary <- data.frame(
  FileName = targets$FileName,
  Mean_A_Correlation = mean_corr_A
)

write.csv(
  qc_summary,
  file = "results/tables/qc_mean_array_correlations.csv",
  row.names = FALSE
)

# ----------------------------
# 6. Identify and remove bad arrays
# ----------------------------
# Based on prior QC inspection, remove all slides from batch 310002
# and one problematic slide from batch 310014.

batch002 <- apply(expand.grid(1:2, 1:4), 1, function(x) {
  paste0("US84700254_256408310002_S01_GE2_1100_Jul11_", x[1], "_", x[2])
})

batch014 <- "US84700254_256408310014_S01_GE2_1100_Jul11_2_4"

bad_arrays <- c(batch002, batch014)

# Strip .txt extension from target file names for matching
array_names <- sub("\\.txt$", "", basename(targets$FileName))

keep <- !(array_names %in% bad_arrays)

removed_arrays <- data.frame(
  FileName = targets$FileName[!keep],
  Removed = TRUE
)

write.csv(
  removed_arrays,
  file = "results/tables/qc_removed_arrays.csv",
  row.names = FALSE
)

# Filter MA object manually
MA_filtered <- MA
MA_filtered$M <- MA$M[, keep, drop = FALSE]
MA_filtered$A <- MA$A[, keep, drop = FALSE]
MA_filtered$targets <- targets[keep, , drop = FALSE]

targets_filtered <- targets[keep, , drop = FALSE]

# Preserve probe identifiers after filtering
if (!is.null(MA_filtered$genes) && "ProbeName" %in% colnames(MA_filtered$genes)) {
  rownames(MA_filtered$M) <- MA_filtered$genes$ProbeName
  rownames(MA_filtered$A) <- MA_filtered$genes$ProbeName
}

# ----------------------------
# 7. Recompute correlation matrices after filtering
# ----------------------------
cor_A_filtered <- cor(MA_filtered$A, use = "pairwise.complete.obs")
cor_M_filtered <- cor(MA_filtered$M, use = "pairwise.complete.obs")

message("Post-QC mean A-value correlation: ",
        round(mean(cor_A_filtered[upper.tri(cor_A_filtered)]), 3))
message("Post-QC A-value correlation range: ",
        paste(round(range(cor_A_filtered[upper.tri(cor_A_filtered)]), 3), collapse = " - "))

# ----------------------------
# 8. QC plots
# ----------------------------

# Boxplot of normalized M-values
png("results/figures/qc/boxplot_normalized_M_values.png", width = 1800, height = 1200, res = 200)
boxplot(MA$M, outline = FALSE, main = "Normalized log-ratios (M)")
dev.off()

# MA plot before and after normalization for first array
png("results/figures/qc/MAplot_before_after_normalization_array1.png", width = 1800, height = 900, res = 200)
par(mfrow = c(1, 2))
plotMA(RG, array = 1, main = "Before normalization")
plotMA(MA, array = 1, main = "After normalization")
dev.off()

# Histogram of mean A-value correlations
png("results/figures/qc/hist_mean_A_correlations.png", width = 1600, height = 1200, res = 200)
hist(
  mean_corr_A,
  breaks = 30,
  col = "lightblue",
  main = "Distribution of mean A-value correlations",
  xlab = "Mean Pearson correlation with other arrays"
)
abline(v = c(0.8, 0.7), col = c("orange", "red"), lty = 2)
dev.off()

# A-value correlation heatmap after QC
png("results/figures/qc/heatmap_Avalue_correlations_postQC.png", width = 1800, height = 1600, res = 200)
draw(
  Heatmap(
    cor_A_filtered,
    name = "Pearson r",
    col = colorRamp2(c(0.75, 0.9, 1), c("navy", "white", "firebrick3")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_dend_width = unit(1.2, "cm"),
    column_dend_height = unit(1.2, "cm"),
    use_raster = FALSE,
    heatmap_legend_param = list(
      title = "A-value correlation",
      at = c(0.75, 0.85, 0.95, 1.0),
      labels = c("0.75", "0.85", "0.95", "1.0"),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    ),
    column_title = "Between-array correlation (post-QC)",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    border = TRUE
  )
)
dev.off()

# M-value correlation heatmap after QC
png("results/figures/qc/heatmap_Mvalue_correlations_postQC.png", width = 1800, height = 1600, res = 200)
draw(
  Heatmap(
    cor_M_filtered,
    name = "M-value correlation",
    col = colorRamp2(c(0, 0.5, 1), c("navy", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = FALSE,
    column_title = "Between-array correlation of M-values",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(title = "Pearson r")
  )
)
dev.off()

# ----------------------------
# 9. Save cleaned objects
# ----------------------------

# Preferred modern outputs
saveRDS(MA_filtered, "results/rds/MA_filtered.rds")
saveRDS(targets_filtered, "results/rds/targets_filtered.rds")
saveRDS(cor_A_filtered, "results/rds/cor_A_filtered.rds")
saveRDS(cor_M_filtered, "results/rds/cor_M_filtered.rds")

# Legacy output used by downstream scripts
save(MA_filtered, targets_filtered, file = "results/microarray_QC_clean.RData")

message("QC complete. Saved cleaned objects and legacy RData file.")
