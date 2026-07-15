#!/usr/bin/env bash
# Stage 01: Quality control on raw FASTQ files (FastQC)
# Produces HTML/ZIP reports for all raw FASTQ files

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
FASTQ_DIR="${FASTQ_DIR:-.}"
THREADS="${THREADS:-8}"
QC_OUTPUT="${QC_OUTPUT:-results/01_fastqc_raw}"

# Ensure output directory exists
mkdir -p "$QC_OUTPUT"

echo "[INFO] FastQC Stage 01: Raw FASTQ Quality Control"
echo "[INFO] Input directory: $FASTQ_DIR"
echo "[INFO] Output directory: $QC_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo ""

# Find all FASTQ files (both .fastq and .fastq.gz formats)
FASTQ_FILES=$(find "$FASTQ_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) 2>/dev/null | sort)

if [ -z "$FASTQ_FILES" ]; then
  echo "[WARNING] No FASTQ files found in $FASTQ_DIR"
  echo "[INFO] Skipping FastQC (no files to process)"
  exit 0
fi

# Count files
FILE_COUNT=$(echo "$FASTQ_FILES" | wc -l)
echo "[INFO] Found $FILE_COUNT FASTQ file(s) to process"
echo ""

# Run FastQC on all files
fastqc \
  --threads "$THREADS" \
  --outdir "$QC_OUTPUT" \
  $FASTQ_FILES

echo "[INFO] FastQC Stage 01 complete"
echo "[INFO] Results saved to: $QC_OUTPUT"
