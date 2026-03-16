#!/bin/bash
# =============================================================================
# diamond_blastp.sh
# Purpose: Map Plodia interpunctella proteins to Drosophila melanogaster
#          best-hit orthologs using DIAMOND BLASTP.
#
# Prerequisites:
#   - DIAMOND installed (brew install diamond  OR  conda install -c bioconda diamond)
#   - Plodia_interpunctella_v1_-_proteins.fa  (from LepBase v4 S3 archive)
#   - uniprotkb_proteome_UP000000803_*.fasta  (Dmel proteome from UniProt)
#
# Output:
#   - plodia_to_dmel.tsv  (tab-separated, 12 columns, one row per Plodia protein)
#
# Run from the directory containing both FASTA files:
#   cd ~/Desktop/Microarray
#   bash diamond_blastp.sh
# =============================================================================

set -euo pipefail  # exit on error, undefined variable, or pipe failure

# -----------------------------------------------------------------------------
# Paths — edit these if your files are named or located differently
# -----------------------------------------------------------------------------
PLODIA_FA="Plodia_interpunctella_v1_-_proteins.fa"
DMEL_FA="uniprotkb_proteome_UP000000803_2026_01_21.fasta"
DB_NAME="dmel_proteome"
OUT_TSV="plodia_to_dmel.tsv"
THREADS=8       # adjust to your machine (use nproc to check available cores)

# -----------------------------------------------------------------------------
# Step 1: Build Drosophila protein database
# (only needs to be done once; skip if dmel_proteome.dmnd already exists)
# -----------------------------------------------------------------------------
if [ -f "${DB_NAME}.dmnd" ]; then
    echo "[$(date +%T)] Database ${DB_NAME}.dmnd already exists — skipping makedb."
else
    echo "[$(date +%T)] Building DIAMOND database from ${DMEL_FA}..."
    diamond makedb \
        --in  "${DMEL_FA}" \
        -d    "${DB_NAME}" \
        --threads "${THREADS}"
    echo "[$(date +%T)] Database built: ${DB_NAME}.dmnd"
fi

# -----------------------------------------------------------------------------
# Step 2: Run DIAMOND BLASTP
#
# Key parameters:
#   --max-target-seqs 1   report only the single best Drosophila hit per Plodia protein
#   --evalue 1e-5         discard weak/random hits at source
#   --outfmt 6 ...        tab-separated output with the 12 columns R_08b expects
#   --threads 8           parallelise across 8 CPU cores
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Running DIAMOND BLASTP..."
echo "  Query:    ${PLODIA_FA}  (Plodia proteins)"
echo "  Database: ${DB_NAME}.dmnd  (Drosophila melanogaster)"
echo "  Output:   ${OUT_TSV}"

diamond blastp \
    --db            "${DB_NAME}.dmnd" \
    --query         "${PLODIA_FA}" \
    --out           "${OUT_TSV}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen \
               qstart qend sstart send evalue bitscore \
    --max-target-seqs 1 \
    --evalue        1e-5 \
    --threads       "${THREADS}"

echo "[$(date +%T)] DIAMOND complete."
echo ""
echo "Output file: ${OUT_TSV}"
echo "Row count (= number of Plodia proteins with a Drosophila hit):"
wc -l < "${OUT_TSV}"
echo ""
echo "First 3 rows:"
head -3 "${OUT_TSV}"
echo ""
echo "Next step: run R_08b_annotate_clusters.R in R."
