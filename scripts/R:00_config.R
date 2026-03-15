# R/00_config.R
# Purpose: shared configuration + dependencies for the limma pipeline.
# Outputs: intermediate/config.rds (objects needed by the next block)
#
# NOTES ON EXPERIMENTAL DESIGN:
#   - Organism: Plodia interpunctella (Indian meal moth)
#   - Array:    Agilent 2-color custom microarray (design ID 064083)
#               384 rows x 160 columns; ~56,737 experimental probes
#   - Design:   Loop design (two-color Cy3/Cy5 co-hybridisation)
#   - Factors:  diet (high/low fat), virus (virus/control), time (see time_levels below)
#
# TIME POINT ENCODING (resolved via SDRF cross-check):
#   All times are in HOURS. The original targets file used a mixed encoding
#   (numeric = hours, "Nd" suffix = days) which was a labelling convenience only.
#   The SDRF (E-MTAB-5868) is the authoritative source and uses hours throughout.
#   Canonical values: 0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168
#
# METADATA NOTES (cross-checked against SDRF — both confirmed as labelling artefacts only):
#   1. Slide 310012_2_1 Cy3: targets label "LFV 1" but cy3_time = 16 — SDRF confirms 16h, label is wrong
#   2. Slide 310004_2_4 Cy3: targets label "HFC 5D" but cy3_time = 3d (72h) — SDRF confirms 72h, label is wrong
#   3. Slide 310005_1_4 cy5_food: trailing whitespace ("control ") — stripped in QC script
#
# BAD ARRAYS (excluded in QC script due to low A-value correlation, mean < 0.7):
#   - All 8 slides from batch 310002 (Jul11)
#   - Slide US84700254_256408310014_S01_GE2_1100_Jul11_2_4
#
# PROBE ANNOTATION NOTES:
#   - Probes with _X suffix (~11,295) are Agilent-flagged low-quality / cross-hybridising probes.
#     Ensure these are filtered during QC (ControlType or ProbeQuality filtering in read.maimages).
#   - Many probes have Uniprot = "unknown" — expected for a non-model organism.
#     Downstream enrichment analysis will rely on SystematicName linking to lepbase.org genome.
#   - Probe gene models derive from two annotation pipelines (maker / augustus) — be aware
#     when collapsing probes to genes.

# -----------------------------
# 0) Housekeeping
# -----------------------------
rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE
)

set.seed(1)

# -----------------------------
# 1) Dependencies
# -----------------------------
# Keep this list centralized so you don't duplicate library() calls across blocks.
required_pkgs <- c(
  "limma",
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "forcats",
  "purrr",
  "readr",
  "openxlsx",
  "patchwork",
  "Biostrings",
  "gprofiler2"
)

missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    "\nInstall them (e.g., install.packages(...) or BiocManager::install(...)) and re-run."
  )
}

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(purrr)
  library(readr)
  library(openxlsx)
  library(patchwork)
  library(Biostrings)
  library(gprofiler2)
})

# -----------------------------
# 2) Paths (edit these once)
# -----------------------------
paths <- list(
  input_rdata      = "microarray_QC_clean.RData",
  intermediate_dir = "intermediate",
  results_dir      = "results",
  figures_dir      = file.path("results", "figures"),
  tables_dir       = file.path("results", "tables")
)

# Create output dirs if needed
dir.create(paths$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$results_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(paths$figures_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,       showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3) Shared analysis constants
# -----------------------------
# All time points in HOURS (confirmed via SDRF cross-check).
# Order is chronological. Use as ordered factor levels in the design matrix script.
time_levels_hours <- c(0.5, 1, 2, 4, 8, 16, 24, 48, 72, 96, 120, 144, 168)

thresholds <- list(
  fdr = 0.05
)

clustering <- list(
  k = 4
)

# -----------------------------
# 4) Minimal plotting defaults
# -----------------------------
plot_cfg <- list(
  base_size      = 30,
  deg_label_size = 8
)

# -----------------------------
# 5) Save exactly what the next block needs
# -----------------------------
config <- list(
  paths             = paths,
  time_levels_hours = time_levels_hours,
  thresholds        = thresholds,
  clustering        = clustering,
  plot_cfg          = plot_cfg
)

saveRDS(config, file = file.path(paths$intermediate_dir, "config.rds"))
message("Saved config to: ", file.path(paths$intermediate_dir, "config.rds"))
