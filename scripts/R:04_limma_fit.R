# R/04_limma_fit.R
# Purpose: Fit limma model using the saved design + contrasts, then create a
#          single canonical results table (results_all) with diet/time metadata.
#
# PROBE FILTERING (applied before lmFit):
#   Step 1: Remove Agilent control probes (ControlType != 0).
#           Keeps ~56,737 experimental probes; removes 56 spike controls.
#   Step 2: Remove _X flagged probes (Agilent low-quality / cross-hybridising).
#           Removes ~11,295 probes; leaves ~45,442 probes for modelling.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/inputs.rds
#   - intermediate/design.rds
#   - intermediate/contrasts.rds
#
# Writes:
#   - intermediate/fit2.rds         (eBayes fit object, ~45k probes x 26 contrasts)
#   - intermediate/results_all.rds  (long-format results table with diet/time columns)

# -----------------------------
# 0) Load config + artifacts
# -----------------------------
config        <- readRDS(file.path("intermediate", "config.rds"))
inputs        <- readRDS(file.path("intermediate", "inputs.rds"))
design_obj    <- readRDS(file.path("intermediate", "design.rds"))
contrasts_obj <- readRDS(file.path("intermediate", "contrasts.rds"))

MA_filtered     <- inputs$MA_filtered
design          <- design_obj$design
contrast_matrix <- contrasts_obj$contrast_matrix
contrast_key    <- contrasts_obj$contrast_key
fdr_cutoff      <- config$thresholds$fdr

# -----------------------------
# 1) Filter probes before fitting
# -----------------------------
# Step 1: Remove Agilent control probes (positive and negative spike controls).
#         ControlType == 0 are experimental probes; != 0 are array controls.
if (is.null(MA_filtered$genes$ControlType)) {
  stop(
    "MA_filtered$genes does not contain a ControlType column.\n",
    "This column is needed to remove Agilent control probes before fitting.\n",
    "Check that read.maimages() loaded the genes annotation correctly."
  )
}

keep_experimental <- MA_filtered$genes$ControlType == 0
MA_fit <- MA_filtered[keep_experimental, ]
message("After removing control probes: ", nrow(MA_fit$M), " probes")

# Step 2: Remove _X flagged probes (Agilent flag for cross-hybridising /
#         low-confidence probes). Standard practice is to exclude these.
if (is.null(MA_fit$genes$ProbeName)) {
  warning(
    "MA_fit$genes does not contain a ProbeName column — cannot filter _X probes.\n",
    "Proceeding with all ControlType==0 probes. Check genes annotation."
  )
} else {
  keep_quality <- !grepl("_X$", MA_fit$genes$ProbeName)
  MA_fit <- MA_fit[keep_quality, ]
  message("After removing _X flagged probes: ", nrow(MA_fit$M), " probes")
}

message("Probes entering lmFit: ", nrow(MA_fit$M),
        " (from ", nrow(MA_filtered$M), " in MA_filtered)")

# -----------------------------
# 2) Fit limma model
# -----------------------------
# lmFit fits a linear model to each probe using the loop design matrix.
# contrasts.fit re-expresses coefficients as the specified contrasts.
# eBayes computes moderated t-statistics using empirical Bayes variance shrinkage.
fit  <- limma::lmFit(MA_fit, design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

message("lmFit complete: ", nrow(fit2$coefficients), " probes x ",
        ncol(fit2$coefficients), " contrasts")

# -----------------------------
# 3) Extract results for every contrast (canonical results_all)
# -----------------------------
# sort.by = "none" keeps all contrasts in the same probe order (MA_fit row order),
# making the long-format table consistent and easier to join downstream.
contrast_names <- colnames(contrast_matrix)

results_list <- lapply(contrast_names, function(cn) {
  tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
  tt$contrast <- cn
  tt
})

results_all <- dplyr::bind_rows(results_list)

message("results_all: ", nrow(results_all), " rows x ", ncol(results_all), " columns")
message("Expected: ", nrow(MA_fit$M), " probes x ", length(contrast_names),
        " contrasts = ", nrow(MA_fit$M) * length(contrast_names), " rows")

# -----------------------------
# 4) Merge in diet/time metadata
# -----------------------------
# contrast_key columns: diet, time_h, time_label, contrast_name, is_estimable
ck <- contrast_key %>%
  dplyr::filter(is_estimable) %>%
  dplyr::select(diet, time_h, time_label, contrast_name)

results_all <- results_all %>%
  dplyr::left_join(ck, by = c("contrast" = "contrast_name"))

# Sanity check: every result row must have a diet/time assignment
n_unmapped <- sum(is.na(results_all$diet) | is.na(results_all$time_h))
if (n_unmapped > 0) {
  warning(
    n_unmapped, " rows in results_all did not map to diet/time via contrast_key.\n",
    "Expected 0. Check that colnames(contrast_matrix) match contrast_key$contrast_name."
  )
} else {
  message("Metadata join check passed: all rows mapped to diet/time.")
}

# Make diet and time ordered factors for consistent downstream plotting
results_all$diet <- factor(results_all$diet,
                           levels = c("low", "high"))
results_all$time_h <- factor(results_all$time_h,
                             levels = config$time_levels_hours,
                             ordered = TRUE)

# -----------------------------
# 5) Save
# -----------------------------
saveRDS(fit2,         file = file.path(config$paths$intermediate_dir, "fit2.rds"))
saveRDS(results_all,  file = file.path(config$paths$intermediate_dir, "results_all.rds"))

message("Saved fit2 to:        ",
        file.path(config$paths$intermediate_dir, "fit2.rds"))
message("Saved results_all to: ",
        file.path(config$paths$intermediate_dir, "results_all.rds"))

# -----------------------------
# 6) Quick DEG count summary
# -----------------------------
sig_by_contrast <- vapply(contrast_names, function(cn) {
  sum(limma::topTable(fit2, coef = cn, number = Inf)$adj.P.Val < fdr_cutoff,
      na.rm = TRUE)
}, numeric(1))

message("DEG counts at FDR < ", fdr_cutoff, " by contrast:")
print(sig_by_contrast)
message("Total DEGs (any contrast): ",
        sum(results_all$adj.P.Val < fdr_cutoff, na.rm = TRUE))
