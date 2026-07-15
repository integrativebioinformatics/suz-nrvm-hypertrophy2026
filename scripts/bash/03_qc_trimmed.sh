#!/usr/bin/env bash
# Stage 03: Quality control on trimmed FASTQ files (FastQC)
# Produces HTML/ZIP reports for all trimmed FASTQ files

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
TRIM_OUTPUT="${TRIM_OUTPUT:-results/02_trimmed_fastq}"
THREADS="${THREADS:-8}"
QC_TRIMMED="${QC_TRIMMED:-results/03_fastqc_trimmed}"

# Ensure output directory exists
mkdir -p "$QC_TRIMMED"

echo "[INFO] FastQC Stage 03: Trimmed FASTQ Quality Control"
echo "[INFO] Input directory: $TRIM_OUTPUT"
echo "[INFO] Output directory: $QC_TRIMMED"
echo "[INFO] Threads: $THREADS"
echo ""

# Find all trimmed FASTQ files
FASTQ_FILES=$(find "$TRIM_OUTPUT" -maxdepth 1 \( -name "*_trimmed.fastq*" \) 2>/dev/null | sort)

if [ -z "$FASTQ_FILES" ]; then
  echo "[WARNING] No trimmed FASTQ files found in $TRIM_OUTPUT"
  echo "[INFO] Did stage 02 (fastp) complete successfully?"
  exit 1
fi

# Count files
FILE_COUNT=$(echo "$FASTQ_FILES" | wc -l)
echo "[INFO] Found $FILE_COUNT trimmed FASTQ file(s) to process"
echo ""

# Run FastQC on all trimmed files
fastqc \
  --threads "$THREADS" \
  --outdir "$QC_TRIMMED" \
  $FASTQ_FILES

echo "[INFO] FastQC Stage 03 complete"
echo "[INFO] Results saved to: $QC_TRIMMED"
