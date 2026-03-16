# R/QC_microarray.R
# Purpose: Load raw Agilent two-color microarray data, normalize, perform QC,
#          remove bad arrays, and save a clean MAList + targets for downstream analysis.
#
# Organism: Plodia interpunctella (Indian meal moth)
# Array:    Agilent custom 2-color (design ID 064083, A-MTAB-618)
#           ~56,737 experimental probes; ~11,295 _X flagged probes removed via ControlType
# Design:   Loop design; factors are diet (high/low), virus (virus/control), time (hours)
#
# Reads:
#   - Raw Agilent Feature Extraction .txt files (one per array)
#   - Targets/metadata file (Paterson_lab_pigv_infection_timeseries_diet_microarray_slides.txt)
#
# Writes:
#   - microarray_QC_clean.RData  (MA_filtered, targets_filtered)
#   - QC figures (optional, displayed inline)
#
# TIME POINT ENCODING:
#   The original targets file uses a mixed encoding (numeric = hours, "Nd" suffix = days).
#   This script rebuilds cy3_time and cy5_time from the SDRF (E-MTAB-5868), which records
#   all times in hours. Canonical values: 0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168
#   After section 8, targets_filtered$cy3_time and cy5_time are numeric hours throughout.
#
# KNOWN BAD ARRAYS (identified by low between-array A-value correlation, mean < 0.7):
#   - All 8 slides from batch 310002 (Jul11)
#   - Slide US84700254_256408310014_S01_GE2_1100_Jul11_2_4
#   These are excluded from MA_filtered and targets_filtered.

# -----------------------------
# 0) Dependencies
# -----------------------------
suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
# -----------------------------
# 1) Paths — edit these to match your local directory structure
# -----------------------------
# file_path:    path to the directory containing raw Agilent .txt files
#               (downloaded from ArrayExpress accession E-MTAB-5868)
# targets_path: path to the targets/metadata file
#               (provided in the repository as data/targets.txt)
# sdrf_path:    path to the SDRF metadata file from EBI
#               (provided in the repository as data/E-MTAB-5868_sdrf.txt)
# output_rdata: name of the output .RData file written to the working directory
 
file_path    <- "path/to/E-MTAB-5868"          # directory of raw Agilent .txt files
targets_path <- "data/targets.txt"             # targets/metadata file
sdrf_path    <- "data/E-MTAB-5868_sdrf.txt"    # SDRF metadata from EBI
output_rdata <- "microarray_QC_clean.RData"    # output file (written to working directory)

# -----------------------------
# 2) Load and clean targets metadata
# -----------------------------
targets <- read.delim(targets_path, sep = "\t", stringsAsFactors = FALSE)
colnames(targets)[1] <- "FileName"

# Strip trailing whitespace from all character columns (e.g. "control " in cy5_virus)
targets <- targets %>% mutate(across(where(is.character), str_trim))

# Match filenames to what read.maimages expects: strip .gz, keep .txt
targets$FileName <- sub("\\.txt\\.gz$", ".txt", targets$FileName)

# Build full paths and check for missing files
files_full    <- file.path(file_path, targets$FileName)
missing_files <- targets$FileName[!file.exists(files_full)]

if (length(missing_files) > 0) {
  message("WARNING: ", length(missing_files), " file(s) not found:")
  print(missing_files)
}

# Keep only arrays with existing files
targets    <- targets[file.exists(files_full), ]
files_full <- files_full[file.exists(files_full)]

message("Arrays to load: ", nrow(targets))

# -----------------------------
# 3) Read raw data
# -----------------------------
# source = "agilent" reads Agilent Feature Extraction files.
# green.only = FALSE retains both Cy3 (green) and Cy5 (red) channels for two-color analysis.
RG <- read.maimages(files_full, source = "agilent", green.only = FALSE)

# Attach targets metadata to the MAList targets slot so they stay in sync
RG$targets <- targets

# -----------------------------
# 4) Normalization
# -----------------------------
# Step 1: Background correction using normexp with offset = 50.
#         normexp is recommended for Agilent data; offset stabilizes variance at low intensities.
RG <- backgroundCorrect(RG, method = "normexp", offset = 50)

# Step 2: Within-array loess normalization corrects red/green dye bias per slide.
MA <- normalizeWithinArrays(RG, method = "loess")

# Step 3: Between-array quantile normalization (A-quantile) makes log-ratios comparable
#         across slides. Aquantile normalizes on A-values, appropriate for loop designs.
MA <- normalizeBetweenArrays(MA, method = "Aquantile")

# -----------------------------
# 5) Pre-QC visualization
# -----------------------------
# Boxplot of M-values across all arrays before filtering
boxplot(MA$M, outline = FALSE, las = 2,
        main = "Normalized log-ratios (M) — all arrays pre-QC",
        ylab = "M-value")

# MA plots for first array: before and after within-array normalization
par(mfrow = c(1, 2))
limma::plotMA(RG, array = 1, main = "Before normalization (array 1)")
limma::plotMA(MA, array = 1, main = "After loess normalization (array 1)")
par(mfrow = c(1, 1))

# -----------------------------
# 6) Between-array QC using A-value correlations
# -----------------------------
# A-values (average log-intensity) reflect overall array quality independent of
# the biological contrast. Low correlation with other arrays flags technical failures.
cor_A <- cor(MA$A, use = "pairwise.complete.obs")

message("A-value correlation summary (all arrays):")
message("  Mean:  ", round(mean(cor_A[upper.tri(cor_A)]), 4))
message("  Range: ", paste(round(range(cor_A[upper.tri(cor_A)]), 4), collapse = " – "))

# Distribution of mean per-array correlations
mean_corr <- rowMeans(cor_A)
hist(mean_corr, breaks = 30, col = "lightblue",
     main = "Distribution of mean A-value correlations (all arrays)",
     xlab = "Mean Pearson correlation with other arrays")
abline(v = c(0.8, 0.7), col = c("orange", "red"), lty = 2)
legend("topleft", legend = c("threshold = 0.8", "threshold = 0.7"),
       col = c("orange", "red"), lty = 2, bty = "n")

# Heatmap of A-value correlations
Heatmap(
  cor_A,
  name = "Pearson r",
  col = colorRamp2(c(0.6, 0.8, 0.9, 1), c("darkblue", "blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = FALSE,
  row_dend_width = unit(1.5, "cm"),
  column_dend_height = unit(1.5, "cm"),
  column_title = "Between-array A-value correlation (pre-QC)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Pearson r")
)

# Report arrays below each threshold for informed decision-making
message("Arrays with mean A-value correlation < 0.8:")
print(names(mean_corr)[mean_corr < 0.8])

message("Arrays with mean A-value correlation < 0.7:")
print(names(mean_corr)[mean_corr < 0.7])

# -----------------------------
# 7) Define and apply array exclusions
# -----------------------------
# Bad arrays identified by low A-value correlation (mean < 0.7):
#   - All 8 slides from batch 310002
#   - Slide 310014_2_4
#
# NOTE: bad_array names must match colnames(MA$A) exactly.
# read.maimages sets colnames to the full file path. We compare using basename
# after stripping the .txt extension from both sides.

bad_arrays <- c(
  # All 8 slides from batch 310002 (blocks 1-2, positions 1-4)
  paste0("US84700254_256408310002_S01_GE2_1100_Jul11_",
         rep(1:2, each = 4), "_", rep(1:4, times = 2)),
  # Single bad slide from batch 310014
  "US84700254_256408310014_S01_GE2_1100_Jul11_2_4"
)

# Strip .txt from array column names so they match bad_arrays (which have no extension)
array_basenames <- gsub("\\.txt$", "", basename(colnames(MA$A)))

# Verify all bad_arrays are actually present — warn if any are missing
not_found <- setdiff(bad_arrays, array_basenames)
if (length(not_found) > 0) {
  warning("The following bad_arrays entries were not found in the data and will be ignored:\n",
          paste(not_found, collapse = "\n"))
}

keep <- !(array_basenames %in% bad_arrays)
message("Arrays before exclusion: ", length(keep))
message("Arrays excluded:         ", sum(!keep))
message("Arrays retained:         ", sum(keep))

# Use the [.MAList method to subset ALL slots consistently:
# $M, $A, $R, $G, $Rb, $Gb, and $targets are all subsetted together.
# This avoids the silent mismatch that occurs when only $M and $A are manually subsetted.
MA_filtered <- MA[, keep]

# Subset the standalone targets data frame to match
targets_filtered <- targets[keep, ]

# Confirm alignment
stopifnot(nrow(targets_filtered) == ncol(MA_filtered$M))
message("Alignment check passed: nrow(targets_filtered) == ncol(MA_filtered$M) == ",
        nrow(targets_filtered))

# -----------------------------
# 8) Rebuild time columns from SDRF (authoritative source)
# -----------------------------
# The original targets file used a mixed encoding: numeric values for hours
# (e.g. "0.5", "1", "16") and "Nd" suffixes for days (e.g. "1d" = 24h, "7d" = 168h).
# The SDRF uses hours throughout and is the ground truth.
# We replace cy3_time and cy5_time in targets_filtered with clean numeric hours
# sourced directly from the SDRF, matched by FileName and channel (Cy3/Cy5).

if (!file.exists(sdrf_path)) {
  stop("SDRF file not found at: ", sdrf_path,
       "\nSet sdrf_path correctly in section 1.")
}

sdrf <- read.delim(sdrf_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Build a lookup: FileName + Label -> time in hours
sdrf_time <- sdrf %>%
  select(
    FileName = `Array Data File`,
    Label,
    time_h = `Factor Value[time]`
  ) %>%
  mutate(time_h = as.numeric(time_h))

# Pivot to wide format: one row per array, columns cy3_time_h and cy5_time_h
sdrf_wide <- sdrf_time %>%
  pivot_wider(names_from = Label, values_from = time_h) %>%
  dplyr::rename(cy3_time_h = "Cy3", cy5_time_h = "Cy5")

# Join onto targets_filtered by FileName
n_before <- nrow(targets_filtered)
targets_filtered <- targets_filtered %>%
  left_join(sdrf_wide, by = "FileName")

# Verify join was 1:1 — no rows should have been lost or duplicated
stopifnot(nrow(targets_filtered) == n_before)

# Check for any NAs introduced by the join (arrays in targets but not in SDRF)
na_cy3 <- sum(is.na(targets_filtered$cy3_time_h))
na_cy5 <- sum(is.na(targets_filtered$cy5_time_h))
if (na_cy3 > 0 || na_cy5 > 0) {
  warning("SDRF join produced NAs: ", na_cy3, " in cy3_time_h, ",
          na_cy5, " in cy5_time_h.\n",
          "These arrays had no matching SDRF entry — check FileName alignment.")
} else {
  message("SDRF time join: all ", nrow(targets_filtered),
          " arrays matched successfully.")
}

# Replace the old mixed-encoding time columns with the clean SDRF hours
# and drop the now-redundant originals
targets_filtered <- targets_filtered %>%
  mutate(
    cy3_time = cy3_time_h,
    cy5_time = cy5_time_h
  ) %>%
  select(-cy3_time_h, -cy5_time_h)

# Confirm canonical values only
canonical_hours <- c(0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168)
unexpected <- setdiff(
  unique(c(targets_filtered$cy3_time, targets_filtered$cy5_time)),
  canonical_hours
)
if (length(unexpected) > 0) {
  warning("Unexpected time values after SDRF join: ",
          paste(unexpected, collapse = ", "))
} else {
  message("Time column check passed: all values are canonical hours (numeric).")
}

# Sync the updated targets back into the MAList targets slot so they stay consistent
MA_filtered$targets <- targets_filtered

# -----------------------------
# 9) Post-QC visualization
# -----------------------------
cor_A_filtered <- cor(MA_filtered$A, use = "pairwise.complete.obs")

message("A-value correlation summary (post-QC arrays):")
message("  Mean:  ", round(mean(cor_A_filtered[upper.tri(cor_A_filtered)]), 4))
message("  Range: ", paste(round(range(cor_A_filtered[upper.tri(cor_A_filtered)]), 4),
                           collapse = " – "))

Heatmap(
  cor_A_filtered,
  name = "Pearson r",
  col = colorRamp2(c(0.75, 0.9, 1), c("navy", "white", "firebrick3")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = FALSE,
  row_dend_width = unit(1.2, "cm"),
  column_dend_height = unit(1.2, "cm"),
  column_title = paste0("Between-array A-value correlation (post-QC, n = ",
                        sum(keep), " arrays)"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "A-value\ncorrelation",
    at = c(0.75, 0.85, 0.95, 1.0),
    labels = c("0.75", "0.85", "0.95", "1.0"),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# M-value correlation heatmap (post-QC)
cor_M_filtered <- cor(MA_filtered$M, use = "pairwise.complete.obs")

Heatmap(
  cor_M_filtered,
  name = "M-value\ncorrelation",
  col = colorRamp2(c(0, 0.5, 1), c("navy", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = FALSE,
  column_title = "Between-array M-value correlation (post-QC)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Pearson r")
)

# Boxplot of M-values post-QC
boxplot(MA_filtered$M, outline = FALSE, las = 2,
        main = paste0("Normalized log-ratios (M) — post-QC (n = ", sum(keep), " arrays)"),
        ylab = "M-value")

# -----------------------------
# 10) Save
# -----------------------------
save(MA_filtered, targets_filtered, file = output_rdata)
message("Saved: ", output_rdata)
message("  MA_filtered:      ", nrow(MA_filtered$M), " probes x ", ncol(MA_filtered$M), " arrays")
message("  targets_filtered: ", nrow(targets_filtered), " rows x ", ncol(targets_filtered), " columns")
