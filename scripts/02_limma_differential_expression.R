# =========================================================
# 02_limma_differential_expression.R
#
# Purpose:
# Load QC-filtered two-color microarray data, fit limma models using a
# proper two-color direct-comparison design (Cy5 vs Cy3), estimate
# virus-vs-control contrasts within each diet x time combination, and
# save differential expression results for downstream analyses.
#
# Inputs:
# - results/microarray_QC_clean.RData
#
# Outputs:
# - results/rds/design_matrix.rds
# - results/rds/fit_two_color_model.rds
# - results/rds/contrast_matrix.rds
# - results/rds/fit_contrasts.rds
# - results/rds/results_list.rds
# - results/tables/all_virus_vs_control_contrasts.csv
# - results/tables/sig_gene_counts_by_contrast.csv
# =========================================================

# ----------------------------
# 1. Load packages
# ----------------------------
library(limma)
library(tidyverse)
library(stringr)

# ----------------------------
# 2. Load QC-filtered data
# ----------------------------
load("results/microarray_QC_clean.RData")

stopifnot(exists("MA_filtered"))
stopifnot(exists("targets_filtered"))

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("results/rds", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 3. Remove control probes
# ----------------------------
if (!is.null(MA_filtered$genes) && "ControlType" %in% colnames(MA_filtered$genes)) {
  
  keep_probes <- MA_filtered$genes$ControlType == 0
  
  MA_filtered$M <- MA_filtered$M[keep_probes, , drop = FALSE]
  MA_filtered$A <- MA_filtered$A[keep_probes, , drop = FALSE]
  MA_filtered$genes <- MA_filtered$genes[keep_probes, , drop = FALSE]
  
  if ("ProbeName" %in% colnames(MA_filtered$genes)) {
    MA_filtered$genes$ProbeName <- as.character(MA_filtered$genes$ProbeName)
    rownames(MA_filtered$M) <- MA_filtered$genes$ProbeName
    rownames(MA_filtered$A) <- MA_filtered$genes$ProbeName
  }
  
  message("Removed control probes. Remaining probes: ", nrow(MA_filtered$M))
}

# Ensure probe identifiers exist
if (is.null(rownames(MA_filtered$M))) {
  if (!is.null(MA_filtered$genes) && "ProbeName" %in% colnames(MA_filtered$genes)) {
    MA_filtered$genes$ProbeName <- as.character(MA_filtered$genes$ProbeName)
    rownames(MA_filtered$M) <- MA_filtered$genes$ProbeName
    rownames(MA_filtered$A) <- MA_filtered$genes$ProbeName
  } else {
    stop("No probe identifiers found in MA_filtered.")
  }
}

# ----------------------------
# 4. Clean and standardize metadata
# ----------------------------
time_levels <- c("0.5", "1", "2", "4", "8", "16",
                 "1d", "2d", "3d", "4d", "5d", "6d", "7d")

targets_filtered <- targets_filtered %>%
  mutate(
    cy3_food  = trimws(tolower(cy3_food)),
    cy3_virus = trimws(tolower(cy3_virus)),
    cy3_time  = trimws(as.character(cy3_time)),
    cy5_food  = trimws(tolower(cy5_food)),
    cy5_virus = trimws(tolower(cy5_virus)),
    cy5_time  = trimws(as.character(cy5_time))
  )

targets_filtered$cy3_time <- factor(targets_filtered$cy3_time, levels = time_levels)
targets_filtered$cy5_time <- factor(targets_filtered$cy5_time, levels = time_levels)

complete_idx <- complete.cases(
  targets_filtered[, c("cy3_food", "cy3_virus", "cy3_time",
                       "cy5_food", "cy5_virus", "cy5_time")]
)

targets_filtered <- targets_filtered[complete_idx, , drop = FALSE]
MA_filtered$M <- MA_filtered$M[, complete_idx, drop = FALSE]
MA_filtered$A <- MA_filtered$A[, complete_idx, drop = FALSE]
MA_filtered$targets <- targets_filtered

stopifnot(ncol(MA_filtered$M) == nrow(targets_filtered))

# ----------------------------
# 5. Build sample labels for BOTH channels
# ----------------------------
targets_filtered <- targets_filtered %>%
  mutate(
    Cy3Group = paste(cy3_food, cy3_virus, as.character(cy3_time), sep = "_"),
    Cy5Group = paste(cy5_food, cy5_virus, as.character(cy5_time), sep = "_")
  )

# Keep a clean mapping table for sanity checks
group_counts <- bind_rows(
  targets_filtered %>% count(Cy3Group, name = "n") %>% mutate(Channel = "Cy3") %>% rename(Group = Cy3Group),
  targets_filtered %>% count(Cy5Group, name = "n") %>% mutate(Channel = "Cy5") %>% rename(Group = Cy5Group)
)

message("Example array comparisons (Cy5 - Cy3):")
print(head(targets_filtered[, c("FileName", "Cy3Group", "Cy5Group")]))

# limma::modelMatrix expects columns named Cy3 and Cy5
targets_for_design <- targets_filtered %>%
  mutate(
    Cy3 = factor(make.names(Cy3Group)),
    Cy5 = factor(make.names(Cy5Group))
  )

all_groups <- sort(unique(c(as.character(targets_for_design$Cy3),
                            as.character(targets_for_design$Cy5))))

targets_for_design$Cy3 <- factor(targets_for_design$Cy3, levels = all_groups)
targets_for_design$Cy5 <- factor(targets_for_design$Cy5, levels = all_groups)

# ----------------------------
# 6. Build proper two-color design matrix manually
# ----------------------------

targets_for_design <- targets_filtered %>%
  mutate(
    Cy3 = factor(make.names(Cy3Group)),
    Cy5 = factor(make.names(Cy5Group))
  )

design <- modelMatrix(targets_for_design, ref = "low_control_0.5")

group_levels <- sort(unique(c(targets_filtered$Cy3Group, targets_filtered$Cy5Group)))

design_full <- matrix(
  0,
  nrow = nrow(targets_filtered),
  ncol = length(group_levels),
  dimnames = list(targets_filtered$FileName, group_levels)
)

for (i in seq_len(nrow(targets_filtered))) {
  design_full[i, targets_filtered$Cy5Group[i]] <-  1
  design_full[i, targets_filtered$Cy3Group[i]] <- -1
}

design_full <- as.matrix(design_full)

# Drop one baseline column to make the design full rank
baseline_group <- "low_control_0.5"

if (!baseline_group %in% colnames(design_full)) {
  stop("Baseline group not found in design_full.")
}

design <- design_full[, colnames(design_full) != baseline_group, drop = FALSE]

stopifnot(nrow(design) == ncol(MA_filtered$M))

saveRDS(design, "results/rds/design_matrix.rds")

message("Baseline group: ", baseline_group)
message("Design matrix dimensions: ", nrow(design), " arrays x ", ncol(design), " coefficients")
message("First few array comparisons:")
print(head(targets_filtered[, c("FileName", "Cy3Group", "Cy5Group")]))

# ----------------------------
# 7. Fit limma model
# ----------------------------
fit <- lmFit(MA_filtered, design)
saveRDS(fit, "results/rds/fit_two_color_model.rds")

# ----------------------------
# 8. Build virus-vs-control contrasts
# ----------------------------
coef_names <- colnames(fit$coefficients)

time_levels <- c("0.5", "1", "2", "4", "8", "16",
                 "1d", "2d", "3d", "4d", "5d", "6d", "7d")

contrast_list <- list()

for (diet in c("high", "low")) {
  for (time in time_levels) {
    
    virus_name   <- make.names(paste(diet, "virus", time, sep = "_"))
    control_name <- make.names(paste(diet, "control", time, sep = "_"))
    
    contrast_name <- paste0("virus_vs_control_", diet, "_", time)
    
    if (virus_name %in% coef_names && control_name %in% coef_names) {
      contrast_list[[contrast_name]] <- paste0(virus_name, " - ", control_name)
      
    } else if (virus_name %in% coef_names && control_name == baseline_group) {
      contrast_list[[contrast_name]] <- virus_name
      
    } else if (control_name %in% coef_names && virus_name == baseline_group) {
      contrast_list[[contrast_name]] <- paste0("-", control_name)
    }
  }
}

if (length(contrast_list) == 0) {
  stop("No valid virus-vs-control contrasts could be constructed.")
}

contrast_matrix <- makeContrasts(
  contrasts = unlist(contrast_list),
  levels = coef_names
)

colnames(contrast_matrix) <- names(contrast_list)

saveRDS(contrast_matrix, "results/rds/contrast_matrix.rds")

message("Contrasts fitted:")
print(colnames(contrast_matrix))

# ----------------------------
# 9. Apply contrasts and empirical Bayes moderation
# ----------------------------
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

saveRDS(fit2, "results/rds/fit_contrasts.rds")

# ----------------------------
# 10. Extract full results for each contrast
# ----------------------------
results_list <- lapply(colnames(contrast_matrix), function(cn) {
  
  tt <- topTable(
    fit2,
    coef = cn,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  if ("ProbeName" %in% colnames(tt)) {
    tt$ProbeID <- as.character(tt$ProbeName)
  } else {
    tt$ProbeID <- rownames(tt)
  }
  
  tt$Contrast <- cn
  tt
})

names(results_list) <- colnames(contrast_matrix)

saveRDS(results_list, "results/rds/results_list.rds")

results_all <- bind_rows(results_list)

write.csv(
  results_all,
  file = "results/tables/all_virus_vs_control_contrasts.csv",
  row.names = FALSE
)

# ----------------------------
# 11. Summarize significant probes
# ----------------------------
sig_summary <- data.frame(
  Contrast = names(results_list),
  SigProbes_FDR_0.05 = sapply(results_list, function(x) sum(x$adj.P.Val < 0.05, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

write.csv(
  sig_summary,
  file = "results/tables/sig_gene_counts_by_contrast.csv",
  row.names = FALSE
)

# ----------------------------
# 12. Print quick summary
# ----------------------------
message("Differential expression analysis complete.")
message("Number of contrasts: ", ncol(contrast_matrix))
message("Significant probe counts (BH-FDR < 0.05):")
print(sig_summary)
