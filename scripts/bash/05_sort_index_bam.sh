#!/usr/bin/env bash
# Stage 05: Convert SAM to BAM, sort, and index
# Converts HISAT2 output (.sam) to compressed BAM format,
# sorts by coordinate, and creates index files

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
ALIGN_OUTPUT="${ALIGN_OUTPUT:-results/04_hisat2_sam}"
THREADS="${THREADS:-8}"
BAM_OUTPUT="${BAM_OUTPUT:-results/05_hisat2_bam}"

# Ensure output directory exists
mkdir -p "$BAM_OUTPUT"

echo "[INFO] Samtools Stage 05: Convert SAM to BAM, Sort, Index"
echo "[INFO] Input directory: $ALIGN_OUTPUT"
echo "[INFO] Output directory: $BAM_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo ""

# Find all SAM files
SAM_FILES=$(find "$ALIGN_OUTPUT" -maxdepth 1 -name "*.sam" 2>/dev/null | sort)

if [ -z "$SAM_FILES" ]; then
  echo "[ERROR] No SAM files found in $ALIGN_OUTPUT" >&2
  echo "[ERROR] Did stage 04 (HISAT2) complete successfully?" >&2
  exit 1
fi

# Count files
FILE_COUNT=$(echo "$SAM_FILES" | wc -l)
echo "[INFO] Found $FILE_COUNT SAM file(s) to process"
echo ""

# Process each SAM file
for sam_file in $SAM_FILES; do
  sample_name=$(basename "$sam_file" .sam)
  bam_file="$BAM_OUTPUT/${sample_name}.bam"

  echo "[INFO] Processing $sample_name..."

  # Convert SAM to BAM and sort in one step
  samtools sort \
    -@ "$THREADS" \
    -o "$bam_file" \
    "$sam_file"

  # Index the BAM file
  samtools index \
    -@ "$THREADS" \
    "$bam_file"

  echo "[INFO]   Created: $bam_file and ${bam_file}.bai"
done

echo "[INFO] Samtools Stage 05 complete"
echo "[INFO] BAM files saved to: $BAM_OUTPUT"
echo "[INFO] Optional: SAM files in $ALIGN_OUTPUT can be deleted to save space"
