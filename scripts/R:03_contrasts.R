# R/03_contrasts.R
# Purpose: Build the limma contrast matrix for virus vs control within each
#          diet x time combination, using the loop-design modelMatrix coefficients.
#
# CONTRAST STRUCTURE:
#   For each diet (low, high) x time (13 timepoints), one contrast:
#     virus_effect = diet.virus.timeh - diet.control.timeh
#   e.g. "low.virus.24h - low.control.24h"
#   This gives 2 x 13 = 26 contrasts total.
#
# COEFFICIENT NAME FORMAT:
#   Coefficients from modelMatrix(ref = "low.control.144h") are named:
#   diet.virus.timeh, e.g. "high.control.0.5h", "low.virus.168h"
#   The reference condition (low.control.144h) has no column in the design —
#   its coefficient is implicitly 0. All other coefficients represent the
#   difference from low.control.144h.
#
# REFERENCE CONDITION HANDLING:
#   The one contrast involving the reference as its control coefficient is:
#     low.144h.virus_vs_ctrl = low.virus.144h - low.control.144h
#   Since low.control.144h is the reference (coef = 0), this simplifies to:
#     low.144h.virus_vs_ctrl = low.virus.144h  (just the coefficient)
#   makeContrasts handles this correctly when coef_ctrl is absent from coef_names.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/design.rds
#
# Writes:
#   - intermediate/contrasts.rds (list: contrast_matrix, contrast_key,
#                                        missing_contrasts)

# -----------------------------
# 0) Load config + design
# -----------------------------
config     <- readRDS(file.path("intermediate", "config.rds"))
design_obj <- readRDS(file.path("intermediate", "design.rds"))

design            <- design_obj$design
coef_names        <- design_obj$coef_names
diets             <- design_obj$diets             # c("low", "high")
time_levels_hours <- design_obj$time_levels_hours # numeric: 0.5 1 2 ... 168
all_conditions    <- design_obj$all_conditions    # all 52 condition labels

# Derive timepoints from all_conditions rather than cy5_time alone,
# to ensure we don't miss timepoints that only appear in the Cy3 channel
timepoints <- sort(unique(as.numeric(
  gsub("h$", "",
       regmatches(all_conditions,
                  regexpr("(?<=\\.)[0-9]+(?:\\.[0-9]+)?h$",
                          all_conditions, perl = TRUE)))
)))
message("Timepoints derived from all conditions: ",
        paste(timepoints, collapse = ", "))

# -----------------------------
# 1) Time label helper
# -----------------------------
# Converts numeric hours to the label used in coef_names
# e.g. 1.0 -> "1h", 0.5 -> "0.5h", 24.0 -> "24h"
make_time_label <- function(t) {
  ifelse(t == as.integer(t),
         paste0(as.integer(t), "h"),
         paste0(t, "h"))
}

# -----------------------------
# 2) Build contrast expressions + key table
# -----------------------------
contrast_list     <- list()
key_rows          <- list()
missing_contrasts <- character(0)

for (d in diets) {
  for (t in timepoints) {
    
    t_label <- make_time_label(t)
    
    # Coefficient names matching the modelMatrix format: diet.virus.timeh
    coef_ctrl  <- paste0(d, ".control.", t_label)
    coef_virus <- paste0(d, ".virus.",   t_label)
    
    virus_ok <- coef_virus %in% coef_names
    ctrl_ok  <- coef_ctrl  %in% coef_names
    ok       <- virus_ok   # virus coef must always be present
    
    # Contrast name: human-readable, e.g. "low.24h.virus_vs_ctrl"
    contrast_name <- paste0(d, ".", t_label, ".virus_vs_ctrl")
    
    # If coef_ctrl is absent it is the reference condition (implicitly 0).
    # The contrast simplifies to just the virus coefficient.
    if (virus_ok && !ctrl_ok) {
      contrast_expr <- coef_virus          # reference condition: ctrl coef = 0
      is_ref_ctrl   <- TRUE
    } else if (virus_ok && ctrl_ok) {
      contrast_expr <- paste0(coef_virus, " - ", coef_ctrl)
      is_ref_ctrl   <- FALSE
    } else {
      contrast_expr <- NA_character_
      is_ref_ctrl   <- FALSE
    }
    
    # Record key row regardless — useful for debugging and downstream joins
    key_rows[[length(key_rows) + 1]] <- data.frame(
      diet          = d,
      time_h        = t,
      time_label    = t_label,
      coef_ctrl     = coef_ctrl,
      coef_virus    = coef_virus,
      contrast_name = contrast_name,
      contrast_expr = contrast_expr,
      is_estimable  = ok,
      is_ref_ctrl   = is_ref_ctrl,
      stringsAsFactors = FALSE
    )
    
    if (ok) {
      contrast_list[[contrast_name]] <- contrast_expr
    } else {
      missing_contrasts <- c(
        missing_contrasts,
        paste0(
          contrast_name,
          " (missing coef(s): ",
          paste(setdiff(c(coef_virus), coef_names), collapse = ", "),
          ")"
        )
      )
    }
  }
}

contrast_key <- do.call(rbind, key_rows)
rownames(contrast_key) <- NULL

# Report what was built
message("Contrasts built:   ", length(contrast_list))
message("Contrasts skipped: ", length(missing_contrasts))

if (length(contrast_list) == 0) {
  stop(
    "No contrasts were constructible.\n",
    "Coefficient names in design_obj do not match the expected format ",
    "'diet.virus.timeh'.\n",
    "First few coef_names: ",
    paste(head(coef_names), collapse = ", ")
  )
}

if (length(missing_contrasts) > 0) {
  message("Skipped contrasts (missing coefficients):")
  message("  - ", paste(missing_contrasts, collapse = "\n  - "))
}

# -----------------------------
# 3) Build contrast matrix
# -----------------------------
contrast_vec <- unlist(contrast_list)

contrast_matrix <- limma::makeContrasts(
  contrasts = contrast_vec,
  levels    = coef_names
)

# Enforce readable column names (makeContrasts may drop them)
colnames(contrast_matrix) <- names(contrast_vec)

message("Contrast matrix dimensions: ",
        nrow(contrast_matrix), " coefficients x ",
        ncol(contrast_matrix), " contrasts")

# Quick check: each contrast column should sum to 0
# (+1 for virus coef, -1 for control coef)
col_sums <- colSums(contrast_matrix)
bad_cols <- names(col_sums)[col_sums != 0]
if (length(bad_cols) > 0) {
  warning("These contrast columns do not sum to 0 (unexpected): ",
          paste(bad_cols, collapse = ", "))
} else {
  message("Contrast column sum check passed: all contrasts sum to 0.")
}

# -----------------------------
# 4) Save exactly what the next block needs
# -----------------------------
contrasts_obj <- list(
  contrast_matrix   = contrast_matrix,
  contrast_key      = contrast_key,
  missing_contrasts = missing_contrasts
)

saveRDS(contrasts_obj,
        file = file.path(config$paths$intermediate_dir, "contrasts.rds"))
message("Saved contrasts to: ",
        file.path(config$paths$intermediate_dir, "contrasts.rds"))
