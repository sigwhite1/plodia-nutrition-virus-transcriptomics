# R/02_design_matrix.R
# Purpose: Build the design matrix for a two-colour loop design using limma's
#          modelMatrix() function, which correctly encodes both channels of each
#          array as +1 (Cy5) and -1 (Cy3) for the relevant condition columns.
#
# DESIGN RATIONALE:
#   This is a proper two-colour loop design (confirmed by authors). Each array's
#   M-value = log2(Cy5/Cy3) represents a contrast between two conditions. The
#   design matrix must encode BOTH channels per array:
#     +1 in the column for the Cy5 condition
#     -1 in the column for the Cy3 condition
#      0 everywhere else
#   limma::modelMatrix() handles this automatically given a targets dataframe
#   with Cy3 and Cy5 condition label columns.
#
#   The previous approach of model.matrix(~ 0 + Group) using only cy5_* metadata
#   was incorrect — it discarded all Cy3 channel information.
#
# CONDITION LABELS:
#   52 unique conditions of the form: diet.virus.timeh
#   e.g. "high.virus.24h", "low.control.96h"
#   diet:  high / low
#   virus: virus / control
#   time:  canonical hours (0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168)
#
# SPARSITY NOTE:
#   high.control.4h has only 1 observation across 102 post-QC arrays.
#   All other conditions have 2-7 observations. Contrasts involving
#   high.control.4h will have reduced power.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/inputs.rds
#
# Writes:
#   - intermediate/design.rds  (list: targets_model, design, coef_names,
#                                      time_levels_hours, diets, virus_levels,
#                                      timepoints, condition_counts)

# -----------------------------
# 0) Load config + inputs
# -----------------------------
config           <- readRDS(file.path("intermediate", "config.rds"))
inputs           <- readRDS(file.path("intermediate", "inputs.rds"))

time_levels_hours <- config$time_levels_hours   # numeric vector: 0.5 1 2 4 8 16 24 48 72 96 120 144 168
targets_filtered  <- inputs$targets_filtered

# -----------------------------
# 1) Build condition labels for both channels
# -----------------------------
# Condition label format: diet.virus.timeh  (e.g. "high.virus.24h", "low.control.0.5h")
# Built from the clean individual columns (cy3_*/cy5_*) rather than the
# inconsistently formatted Cy3/Cy5 label columns in the original targets file.

make_condition <- function(food, virus, time_h) {
  # Convert numeric hours to label: integer hours drop decimal (e.g. 1.0 -> "1h")
  t_str <- ifelse(time_h == as.integer(time_h),
                  paste0(as.integer(time_h), "h"),
                  paste0(time_h, "h"))
  paste(food, virus, t_str, sep = ".")
}

targets_filtered$cy3_condition <- make_condition(
  targets_filtered$cy3_food,
  targets_filtered$cy3_virus,
  targets_filtered$cy3_time
)

targets_filtered$cy5_condition <- make_condition(
  targets_filtered$cy5_food,
  targets_filtered$cy5_virus,
  targets_filtered$cy5_time
)

# Peek at first few rows to confirm labels look right
message("Sample condition pairs (first 5 arrays):")
print(targets_filtered[1:5, c("FileName", "cy3_condition", "cy5_condition")])

# -----------------------------
# 2) Build the modelMatrix targets dataframe
# -----------------------------
# limma::modelMatrix() requires a dataframe with columns named exactly
# "Cy3" and "Cy5" containing the condition label for each channel.
targets_model <- data.frame(
  FileName = targets_filtered$FileName,
  Cy3      = targets_filtered$cy3_condition,
  Cy5      = targets_filtered$cy5_condition,
  stringsAsFactors = FALSE
)

# Confirm all conditions are present in both channels as expected
all_conditions <- sort(unique(c(targets_model$Cy3, targets_model$Cy5)))
message("Total unique conditions: ", length(all_conditions))

# Warn about any conditions appearing in only one channel
# (these reduce connectivity of the loop)
cy3_only <- setdiff(unique(targets_model$Cy3), unique(targets_model$Cy5))
cy5_only <- setdiff(unique(targets_model$Cy5), unique(targets_model$Cy3))
# NOTE: conditions appearing in only one channel is normal in a loop design.
# It means those conditions always play the same role (reference or test) across
# all arrays they appear in. The loop remains connected and estimation is unaffected.
if (length(cy3_only) > 0)
  message("NOTE - Conditions appearing only as Cy3 channel: ",
          paste(cy3_only, collapse = ", "))
if (length(cy5_only) > 0)
  message("NOTE - Conditions appearing only as Cy5 channel: ",
          paste(cy5_only, collapse = ", "))

# -----------------------------
# 3) Build the design matrix using modelMatrix()
# -----------------------------
# modelMatrix() produces a matrix with:
#   nrow = number of arrays (102)
#   ncol = number of unique conditions (52)
#   entries: +1 where condition is Cy5, -1 where condition is Cy3, 0 otherwise
#
# ref argument: the reference condition is subtracted from all others.
# For a loop design with no natural reference, we omit ref (ref = NULL)
# so all 52 conditions get their own column. This is the standard approach
# per the limma User Guide section on loop designs.
# Build design matrix using limma::modelMatrix() with a reference condition.
#
# WHY A REFERENCE IS NEEDED:
#   A connected loop design with n=52 conditions has only n-1=51 estimable
#   degrees of freedom (M-values are ratios, so absolute levels are not
#   identifiable). Without a reference, lmFit drops one coefficient arbitrarily
#   via QR pivoting, which breaks contrasts.fit if that coefficient appears in
#   any contrast.
#
# REFERENCE CHOICE:
#   low.control.144h — chosen as the most observed condition (7 arrays),
#   giving the most stable baseline. The reference cancels out of all
#   pairwise (virus - control) contrasts, so the choice does not affect
#   biological interpretation.
#
# RESULT:
#   102 x 51 design matrix of full rank. Each coefficient represents the
#   difference between that condition and low.control.144h.
#   All 26 virus vs control contrasts remain estimable.

design_ref <- "low.control.144h"

design <- limma::modelMatrix(targets_model, ref = design_ref)

message("Design matrix dimensions: ", nrow(design), " arrays x ", ncol(design), " conditions")
message("Expected: 102 arrays x 51 conditions (52 - 1 reference)")

# Verify row count matches MA_filtered
stopifnot(nrow(design) == ncol(inputs$MA_filtered$M))

# -----------------------------
# 4) Sanity checks
# -----------------------------

# 4a) Check matrix encodes exactly +1 and -1 per row (each array has one Cy5, one Cy3)
row_sums <- rowSums(design)
if (!all(row_sums == 0)) {
  warning("Some design matrix rows do not sum to 0. ",
          "Expected each row to have exactly one +1 and one -1.\n",
          "Non-zero row sums at arrays: ",
          paste(which(row_sums != 0), collapse = ", "))
} else {
  message("Row sum check passed: all rows sum to 0 (one +1 and one -1 per array).")
}

# 4b) Check for all-zero columns (unobserved conditions — should not exist)
zero_cols <- colnames(design)[colSums(abs(design)) == 0]
if (length(zero_cols) > 0) {
  warning("All-zero columns in design (unobserved conditions): ",
          paste(zero_cols, collapse = ", "))
} else {
  message("Column check passed: all ", ncol(design), " conditions observed.")
}

# 4c) Report condition observation counts (how many arrays involve each condition)
condition_counts <- sort(colSums(abs(design)))
message("Condition observation counts (min/median/max): ",
        min(condition_counts), " / ",
        median(condition_counts), " / ",
        max(condition_counts))

sparse_conditions <- names(condition_counts)[condition_counts <= 2]
if (length(sparse_conditions) > 0) {
  message("Sparse conditions (<=2 observations): ",
          paste(sparse_conditions, collapse = ", "))
}

# 4d) Check design matrix rank — must equal ncol(design) for lmFit to work
# With a reference condition, the design matrix should have full rank.
# rank = ncol(design) = 51 is required for lmFit to estimate all coefficients.
design_rank <- qr(design)$rank
if (design_rank == ncol(design)) {
  message("Rank check passed: design matrix has full rank (",
          design_rank, " = ", ncol(design), ").")
} else {
  warning("Design matrix is rank deficient: rank = ", design_rank,
          " but ncol = ", ncol(design), ".\n",
          "This should not happen with a reference condition set. ",
          "Check that the loop is connected and ref = '", design_ref,
          "' appears in both Cy3 and Cy5 columns of targets_model.")
}

# 4e) Confirm time levels in conditions match config
# Extract time values from condition labels (format: diet.virus.Xh)
# Uses a lookbehind for the last dot to correctly handle decimals like "0.5h"
times_in_design <- sort(unique(as.numeric(
  gsub("h$", "",
       regmatches(all_conditions,
                  regexpr("(?<=\\.)[0-9]+(?:\\.[0-9]+)?h$",
                          all_conditions, perl = TRUE)))
)))
unexpected_times <- setdiff(times_in_design, time_levels_hours)
if (length(unexpected_times) > 0) {
  warning("Unexpected time values in condition labels: ",
          paste(unexpected_times, collapse = ", "))
} else {
  message("Time level check passed: all condition times match config$time_levels_hours.")
}

# -----------------------------
# 5) Save exactly what the next block needs
# -----------------------------
# The contrasts script needs:
#   - targets_model: the Cy3/Cy5 condition label dataframe
#   - design:        the 102 x 52 modelMatrix design
#   - coef_names:    condition names (column names of design)
#   - all_conditions: sorted vector of all 52 condition labels
#   - condition_counts: named vector of observation counts per condition
#   - time_levels_hours, diets, virus_levels, timepoints: for contrast construction

design_obj <- list(
  targets_model     = targets_model,
  design            = design,
  coef_names        = colnames(design),
  all_conditions    = all_conditions,
  condition_counts  = condition_counts,
  time_levels_hours = time_levels_hours,
  diets             = c("low", "high"),
  virus_levels      = c("control", "virus"),
  timepoints        = sort(unique(targets_filtered$cy5_time))
)

saveRDS(design_obj, file = file.path(config$paths$intermediate_dir, "design.rds"))
message("Saved design object to: ",
        file.path(config$paths$intermediate_dir, "design.rds"))
