#!/usr/bin/env bash
# Stage 09: lncRNA classification (FEELnc)
# Three-step pipeline to identify, filter, and classify novel lncRNAs
# 1. FEELnc_filter: remove mono-exonic transcripts, keep protein-coding-class
# 2. FEELnc_codpot: predict coding potential
# 3. FEELnc_classifier: classify by genomic location (intergenic, sense, antisense, etc.)

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
GTF_FILE="${GTF_FILE:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.115.gtf}"
GFFCOMPARE_OUTPUT="${GFFCOMPARE_OUTPUT:-results/07_gffcompare}"
FEELNC_OUTPUT="${FEELNC_OUTPUT:-results/09_feelnc_classify}"
FEELNC_WINDOW="${FEELNC_WINDOW:-1000}"
FEELNC_MAXWINDOW="${FEELNC_MAXWINDOW:-100000}"

# Ensure output directory exists
mkdir -p "$FEELNC_OUTPUT"

echo "[INFO] FEELnc Stage 09: lncRNA Classification"
echo "[INFO] Reference GTF: $GTF_FILE"
echo "[INFO] Merged GTF: $GFFCOMPARE_OUTPUT/merged.annotated.gtf"
echo "[INFO] Output directory: $FEELNC_OUTPUT"
echo "[INFO] FEELnc window: $FEELNC_WINDOW (max: $FEELNC_MAXWINDOW)"
echo ""

# Check if input files exist
if [ ! -f "$GTF_FILE" ]; then
  echo "[ERROR] Reference GTF not found: $GTF_FILE" >&2
  exit 1
fi

MERGED_GTF="$GFFCOMPARE_OUTPUT/merged.annotated.gtf"
if [ ! -f "$MERGED_GTF" ]; then
  echo "[ERROR] Merged GTF not found: $MERGED_GTF" >&2
  exit 1
fi

# Step 1: Filter — remove mono-exonic, keep candidates
echo "[INFO] Step 1: Filtering candidate lncRNAs..."

FEELnc_filter.pl \
  -i "$MERGED_GTF" \
  -a "$GTF_FILE" \
  -b "transcript_biotype=protein_coding" \
  --monoex=1 \
  --outlog "$FEELNC_OUTPUT/feelnc_filter.log" \
  > "$FEELNC_OUTPUT/filtered_candidate_lncRNA.gtf"

echo "[INFO] Filtered candidates: $FEELNC_OUTPUT/filtered_candidate_lncRNA.gtf"

# Step 2: Coding Potential — predict which candidates are truly non-coding
echo "[INFO] Step 2: Computing coding potential..."

# Create directory for FEELnc coding potential (expects a specific structure)
CODING_POT_DIR="$FEELNC_OUTPUT/codpot"
mkdir -p "$CODING_POT_DIR"

# FEELnc_codpot requires chromosomes in separate files
# For simplicity, we'll use the merged GTF directly
# Note: the codon usage table may need to be species-specific

FEELnc_codpot.pl \
  -i "$FEELNC_OUTPUT/filtered_candidate_lncRNA.gtf" \
  -a "$GTF_FILE" \
  -g "$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.toplevel.fa" \
  -b "transcript_biotype=protein_coding" \
  -b "transcript_status=KNOWN" \
  --mode=intergenic \
  --outdir "$CODING_POT_DIR" \
  2>&1 | tee "$FEELNC_OUTPUT/feelnc_codpot.log"

echo "[INFO] Coding potential computed"

# Step 3: Classifier — annotate genomic context (intergenic, sense, antisense, etc.)
echo "[INFO] Step 3: Classifying lncRNA genomic locations..."

CODPOT_OUTPUT=$(ls "$CODING_POT_DIR"/*_POTENTIAL_lncRNA.txt 2>/dev/null | head -1)
if [ -z "$CODPOT_OUTPUT" ]; then
  echo "[WARNING] Coding potential output not found, using filtered candidates directly"
  CLASSIFIER_INPUT="$FEELNC_OUTPUT/filtered_candidate_lncRNA.gtf"
else
  CLASSIFIER_INPUT="$CODPOT_OUTPUT"
fi

FEELnc_classifier.pl \
  -i "$CLASSIFIER_INPUT" \
  -a "$GTF_FILE" \
  --window="$FEELNC_WINDOW" \
  --maxwindow="$FEELNC_MAXWINDOW" \
  2>&1 | tee "$FEELNC_OUTPUT/feelnc_classifier.txt"

echo "[INFO] Classification output: $FEELNC_OUTPUT/feelnc_classifier.txt"

echo "[INFO] FEELnc Stage 09 complete"
echo "[INFO] All FEELnc outputs saved to: $FEELNC_OUTPUT"
