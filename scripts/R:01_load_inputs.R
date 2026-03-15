# R/01_load_inputs.R
# Purpose: Load QC-cleaned microarray objects and save exactly what the next block needs.
#
# Reads:
#   - intermediate/config.rds
#   - config$paths$input_rdata  (e.g., "microarray_QC_clean.RData")
#
# Writes:
#   - intermediate/inputs.rds   (list: MA_filtered, targets_filtered)

# -----------------------------
# 0) Load config
# -----------------------------
config <- readRDS(file.path("intermediate", "config.rds"))
paths  <- config$paths

if (!dir.exists(paths$intermediate_dir)) {
  dir.create(paths$intermediate_dir, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# 1) Load the .RData into a private environment
# -----------------------------
if (!file.exists(paths$input_rdata)) {
  stop(
    "Input .RData not found at: ", paths$input_rdata, "\n",
    "Edit config$paths$input_rdata in R/00_config.R if needed."
  )
}

e <- new.env(parent = emptyenv())
loaded_names <- load(paths$input_rdata, envir = e)

message("Loaded from ", paths$input_rdata, ":")
message("  - ", paste(loaded_names, collapse = "\n  - "))

# -----------------------------
# 2) Pull the exact objects we expect
# -----------------------------
required_objects <- c("MA_filtered", "targets_filtered")
missing_objects  <- setdiff(required_objects, loaded_names)

if (length(missing_objects) > 0) {
  stop(
    "The .RData file did not contain required object(s): ",
    paste(missing_objects, collapse = ", "),
    "\nIt contained: ", paste(loaded_names, collapse = ", ")
  )
}

MA_filtered      <- get("MA_filtered",      envir = e)
targets_filtered <- get("targets_filtered", envir = e)

# -----------------------------
# 3) Sanity checks
# -----------------------------

# 3a) MA_filtered must look like a limma-ready object (MAList has $M; EList has $E)
if (!any(c("M", "E") %in% names(MA_filtered))) {
  stop(
    "MA_filtered does not look like a limma-ready object.\n",
    "Expected a list-like object with component 'M' (MAList) or 'E' (EList).\n",
    "Names present: ", paste(names(MA_filtered), collapse = ", ")
  )
}

# 3b) targets_filtered must contain all required columns from both channels.
#     Both cy3_* and cy5_* columns are needed for loop-design matrix construction.
required_target_cols <- c(
  "cy3_food", "cy3_virus", "cy3_time",
  "cy5_food", "cy5_virus", "cy5_time"
)
missing_cols <- setdiff(required_target_cols, colnames(targets_filtered))

if (length(missing_cols) > 0) {
  stop(
    "targets_filtered is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    "\nColumns present: ", paste(colnames(targets_filtered), collapse = ", ")
  )
}

# 3c) Check that array count and targets row count align.
#     A mismatch here will silently corrupt the design matrix — hard stop is intentional.
expr_mat <- NULL
if ("M" %in% names(MA_filtered) && is.matrix(MA_filtered$M)) expr_mat <- MA_filtered$M
if (is.null(expr_mat) && "E" %in% names(MA_filtered) && is.matrix(MA_filtered$E)) expr_mat <- MA_filtered$E

if (!is.null(expr_mat) && ncol(expr_mat) != nrow(targets_filtered)) {
  stop(
    "Sample count mismatch: ncol(expression) = ", ncol(expr_mat),
    " but nrow(targets_filtered) = ", nrow(targets_filtered), "\n",
    "Design matrix construction requires these to align exactly."
  )
}

# 3d) If MA_filtered carries its own targets slot, verify filenames match targets_filtered.
if (!is.null(MA_filtered$targets)) {
  if (!identical(rownames(targets_filtered), rownames(MA_filtered$targets))) {
    warning(
      "Row names of targets_filtered do not match MA_filtered$targets.\n",
      "Verify that array ordering is consistent between the two objects."
    )
  }
}

# 3e) Confirm time columns contain only known canonical hour values.
#     Per SDRF cross-check, all times should be numeric hours:
#     0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168
#     The "Nd" day-suffix encoding in the original targets file was a labelling
#     convenience — the QC script converts everything to hours. Flag anything unexpected.
canonical_hours <- as.character(c(0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168))
time_vals       <- unique(c(as.character(targets_filtered$cy3_time),
                             as.character(targets_filtered$cy5_time)))
unexpected_times <- setdiff(time_vals, canonical_hours)

if (length(unexpected_times) > 0) {
  warning(
    "Unexpected time values found (not in canonical hours list): ",
    paste(unexpected_times, collapse = ", "), "\n",
    "Check whether time columns have been correctly converted to hours in the QC script."
  )
} else {
  message("Time values check passed: all values are canonical hours.")
}

# -----------------------------
# 4) Save exactly what the next block needs
# -----------------------------
inputs <- list(
  MA_filtered      = MA_filtered,
  targets_filtered = targets_filtered
)

saveRDS(inputs, file = file.path(paths$intermediate_dir, "inputs.rds"))
message("Saved inputs to: ", file.path(paths$intermediate_dir, "inputs.rds"))
