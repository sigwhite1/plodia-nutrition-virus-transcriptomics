# R/06_define_gene_sets.R
# Purpose: Define reusable DEG gene sets across ALL diet x time contrasts.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/results_all.rds
#
# Writes:
#   - intermediate/gene_sets.rds        (all_deg_probes, low_deg_probes,
#                                         high_deg_probes, shared_deg_probes)
#   - results/tables/gene_sets_summary.csv  (counts for quick inspection)

# -----------------------------
# 0) Load artifacts
# -----------------------------
config      <- readRDS(file.path("intermediate", "config.rds"))
results_all <- readRDS(file.path("intermediate", "results_all.rds"))

paths      <- config$paths
fdr_cutoff <- config$thresholds$fdr

dir.create(paths$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,       showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Filter to significant rows
# -----------------------------
sig <- results_all %>%
  dplyr::filter(adj.P.Val < fdr_cutoff) %>%
  dplyr::filter(!is.na(diet), !is.na(time_h))

if (nrow(sig) == 0) {
  stop("No significant genes found at FDR < ", fdr_cutoff, ".")
}

message("Significant probe-contrast rows: ", nrow(sig))

# -----------------------------
# 2) Define gene sets
# -----------------------------
all_deg_probes <- sig %>%
  dplyr::pull(ProbeName) %>%
  unique()

low_deg_probes <- sig %>%
  dplyr::filter(diet == "low") %>%
  dplyr::pull(ProbeName) %>%
  unique()

high_deg_probes <- sig %>%
  dplyr::filter(diet == "high") %>%
  dplyr::pull(ProbeName) %>%
  unique()

shared_deg_probes <- intersect(low_deg_probes, high_deg_probes)

# Sanity check: all + shared must be subsets of the full set
stopifnot(all(low_deg_probes    %in% all_deg_probes))
stopifnot(all(high_deg_probes   %in% all_deg_probes))
stopifnot(all(shared_deg_probes %in% all_deg_probes))

gene_sets <- list(
  fdr_cutoff        = fdr_cutoff,
  all_deg_probes    = all_deg_probes,
  low_deg_probes    = low_deg_probes,
  high_deg_probes   = high_deg_probes,
  shared_deg_probes = shared_deg_probes
)

saveRDS(gene_sets, file = file.path(paths$intermediate_dir, "gene_sets.rds"))

# -----------------------------
# 3) Summary table
# -----------------------------
summary_df <- tibble::tibble(
  set        = c("all_deg_probes", "low_deg_probes",
                 "high_deg_probes", "shared_deg_probes"),
  n          = c(length(all_deg_probes),    length(low_deg_probes),
                 length(high_deg_probes),   length(shared_deg_probes)),
  fdr_cutoff = fdr_cutoff
)

readr::write_csv(summary_df,
                 file = file.path(paths$tables_dir, "gene_sets_summary.csv"))

message("Saved gene sets to:   ",
        file.path(paths$intermediate_dir, "gene_sets.rds"))
message("Saved summary to:     ",
        file.path(paths$tables_dir, "gene_sets_summary.csv"))
print(summary_df)
