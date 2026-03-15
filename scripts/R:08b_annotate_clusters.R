# R/08b_annotate_clusters.R
# Purpose: Join k=4 cluster assignments to Drosophila ortholog annotations via
#          DIAMOND BLASTP results, producing the master annotation tables needed
#          by R_09_enrichment_GO_clusters.R.
#
# PIPELINE OVERVIEW:
#   ProbeName
#     -> maker (Plodia gene model ID, from A-MTAB-618_comments.txt)
#     -> Plodia protein sequence (Plodia_interpunctella_v1_-_proteins.fa)
#     -> DIAMOND BLASTP -> Drosophila melanogaster best hit (plodia_to_dmel.tsv)
#     -> Dmel_protein (tr|Q9VFH5|Q9VFH5_DROME format)
#     -> UniProt accession (Q9VFH5) -> FlyBase ID -> GO enrichment (in R_09)
#
# PREREQUISITE — DIAMOND BLASTP (run once in Terminal before this script):
#   This script is the R side of the pipeline. The DIAMOND alignment step must
#   be run separately in Terminal. Commands (adjust paths as needed):
#
#   # 1. Build Drosophila protein database
#   diamond makedb \
#     --in uniprotkb_proteome_UP000000803_2026_01_21.fasta \
#     -d dmel_proteome
#
#   # 2. Run BLASTP (Plodia proteins vs Drosophila database)
#   diamond blastp \
#     --db dmel_proteome.dmnd \
#     --query Plodia_interpunctella_v1_-_proteins.fa \
#     --out plodia_to_dmel.tsv \
#     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
#     --max-target-seqs 1 \
#     --evalue 1e-5 \
#     --threads 8
#
#   This produces plodia_to_dmel.tsv: one row per Plodia protein,
#   reporting the single best-matching Drosophila protein.
#
# FILTER THRESHOLDS (from PDF documentation):
#   evalue  <= 1e-5   (removes random weak hits)
#   pident  >= 35%    (removes distant low-confidence matches)
#   These are applied to produce cluster_annot_all_filtered.rds.
#   cluster_annot_all_raw.rds retains all hits (evalue <= 1e-5 from DIAMOND,
#   no pident filter) for auditing purposes.
#
# NOTE ON COLUMN NAMING:
#   cluster_df from R_08 uses lowercase 'cluster' (integer 1-4).
#   This script preserves that convention. The old annotation files
#   (generated with k=6) used uppercase 'Cluster' — those are now superseded.
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/cluster_df.rds          (ProbeName + cluster, from R_08)
#   - data/A-MTAB-618_comments.txt         (probe annotation — adjust path)
#   - data/plodia_to_dmel.tsv              (DIAMOND output — adjust path)
#
# Writes:
#   - intermediate/cluster_annot_all_raw.rds       (all DIAMOND hits joined to clusters)
#   - intermediate/cluster_annot_all_filtered.rds  (evalue<=1e-5 AND pident>=35)
#   - results/tables/cluster_annotation_summary.csv

# -----------------------------
# 0) Load config
# -----------------------------
config <- readRDS(file.path("intermediate", "config.rds"))
paths  <- config$paths

dir.create(paths$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,       showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Paths — edit these if needed
# -----------------------------
annot_path   <- "~/Desktop/Microarray/A-MTAB-618_comments.txt"   # probe annotation file
diamond_path <- "~/Desktop/Microarray/plodia_to_dmel.tsv"         # DIAMOND BLASTP output

# -----------------------------
# 2) Load cluster assignments
# -----------------------------
cluster_df <- readRDS(file.path(paths$intermediate_dir, "cluster_df.rds"))

message("Cluster assignments loaded: ", nrow(cluster_df), " probes across ",
        dplyr::n_distinct(cluster_df$cluster), " clusters")
print(table(cluster_df$cluster))

# -----------------------------
# 3) Load probe annotation
# -----------------------------
if (!file.exists(annot_path)) {
  stop("Probe annotation file not found: ", annot_path,
       "\nUpdate annot_path in section 1.")
}

annot_raw <- read.delim(annot_path, header = TRUE, sep = "\t",
                        quote = "", stringsAsFactors = FALSE)

# Keep only experimental probes (ControlType == 0) and the columns we need
annot_clean <- annot_raw %>%
  dplyr::filter(ControlType == 0) %>%
  dplyr::select(ProbeName, maker)

message("Probe annotation loaded: ", nrow(annot_clean),
        " experimental probes with maker IDs")

# Sanity check: all cluster probes should be in annotation
missing_from_annot <- setdiff(cluster_df$ProbeName, annot_clean$ProbeName)
if (length(missing_from_annot) > 0) {
  warning(length(missing_from_annot),
          " cluster probes not found in annotation file:\n",
          paste(head(missing_from_annot, 5), collapse = "\n"))
} else {
  message("Annotation join check passed: all cluster probes found in annotation.")
}

# -----------------------------
# 4) Load DIAMOND BLASTP results
# -----------------------------
if (!file.exists(diamond_path)) {
  stop(
    "DIAMOND output file not found: ", diamond_path, "\n",
    "Run the DIAMOND BLASTP step in Terminal first (see script header for commands).\n",
    "Update diamond_path in section 1 if the file is in a different location."
  )
}

diamond_res <- read.table(diamond_path, sep = "\t", header = FALSE,
                          stringsAsFactors = FALSE)

colnames(diamond_res) <- c(
  "Plodia_protein", "Dmel_protein", "pident", "length",
  "mismatch", "gapopen", "qstart", "qend",
  "sstart", "send", "evalue", "bitscore"
)

message("DIAMOND results loaded: ", nrow(diamond_res), " hits")
message("  evalue range: ", signif(min(diamond_res$evalue), 3),
        " – ", signif(max(diamond_res$evalue), 3))
message("  pident range: ", round(min(diamond_res$pident), 1),
        " – ", round(max(diamond_res$pident), 1), "%")

# DIAMOND was run with --evalue 1e-5 so all hits already pass that threshold,
# but apply it explicitly here as a guard in case parameters were changed
diamond_res <- diamond_res %>%
  dplyr::filter(evalue <= 1e-5)

message("After evalue <= 1e-5 filter: ", nrow(diamond_res), " hits retained")

# -----------------------------
# 5) Build master annotation table (raw)
# -----------------------------
# Join order:
#   cluster_df (ProbeName + cluster)
#   -> annot_clean (ProbeName -> maker)
#   -> diamond_res (maker -> Dmel_protein + BLAST stats)
cluster_annot_all_raw <- cluster_df %>%
  dplyr::left_join(annot_clean, by = "ProbeName") %>%
  dplyr::left_join(diamond_res,
                   by = c("maker" = "Plodia_protein"))

message("\nMaster annotation table (raw):")
message("  Rows:                 ", nrow(cluster_annot_all_raw))
message("  Probes with maker ID: ",
        sum(!is.na(cluster_annot_all_raw$maker)))
message("  Probes with Dmel hit: ",
        sum(!is.na(cluster_annot_all_raw$Dmel_protein)))
message("  Probes without hit:   ",
        sum(is.na(cluster_annot_all_raw$Dmel_protein)))

# Sanity check: 0 DEG probes should have missing sequences
# (confirmed in original analysis — log this for reproducibility)
probes_no_dmel <- cluster_annot_all_raw %>%
  dplyr::filter(is.na(Dmel_protein)) %>%
  dplyr::pull(ProbeName)
message("  (", length(probes_no_dmel),
        " probes have no Dmel hit — expected for genomic interval probes)")

saveRDS(cluster_annot_all_raw,
        file.path(paths$intermediate_dir, "cluster_annot_all_raw.rds"))

# -----------------------------
# 6) Apply quality filters -> filtered table
# -----------------------------
# evalue <= 1e-5 : removes random weak hits (already applied above, retained here)
# pident >= 35%  : removes low-confidence distant matches
# These thresholds follow the original analysis (documented in PDF notes).

cluster_annot_all_filtered <- cluster_annot_all_raw %>%
  dplyr::filter(!is.na(Dmel_protein)) %>%
  dplyr::filter(evalue  <= 1e-5) %>%
  dplyr::filter(pident  >= 35)

message("\nFiltered annotation table (evalue <= 1e-5, pident >= 35%):")
message("  Rows:     ", nrow(cluster_annot_all_filtered))
message("  Clusters: ", dplyr::n_distinct(cluster_annot_all_filtered$cluster))
print(table(cluster_annot_all_filtered$cluster))

saveRDS(cluster_annot_all_filtered,
        file.path(paths$intermediate_dir, "cluster_annot_all_filtered.rds"))

# -----------------------------
# 7) Summary table
# -----------------------------
summary_df <- cluster_annot_all_raw %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_probes        = dplyr::n(),
    n_with_maker    = sum(!is.na(maker)),
    n_with_dmel_hit = sum(!is.na(Dmel_protein)),
    n_filtered      = sum(!is.na(Dmel_protein) &
                            evalue <= 1e-5 & pident >= 35),
    pct_annotated   = round(100 * n_with_dmel_hit / n_probes, 1),
    .groups = "drop"
  )

readr::write_csv(summary_df,
                 file.path(paths$tables_dir, "cluster_annotation_summary.csv"))

message("\nAnnotation summary per cluster:")
print(summary_df)

message("\nSaved:")
message("  ", file.path(paths$intermediate_dir, "cluster_annot_all_raw.rds"))
message("  ", file.path(paths$intermediate_dir, "cluster_annot_all_filtered.rds"))
message("  ", file.path(paths$tables_dir, "cluster_annotation_summary.csv"))
message()
message("Next step: run R_09_enrichment_GO_clusters.R")
