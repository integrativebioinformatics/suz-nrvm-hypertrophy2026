#!/usr/bin/env bash
# Stage 02: Adapter and quality trimming (fastp)
# Handles both paired-end and single-end reads with parameters:
#   -f 10 -F 10 (trim front)
#   -3 (trim tail, 3' end)
#   --cut_tail_mean_quality 28 (quality threshold)
#   -q 30 (minimum base quality)

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
FASTQ_DIR="${FASTQ_DIR:-.}"
THREADS="${THREADS:-8}"
TRIM_OUTPUT="${TRIM_OUTPUT:-results/02_trimmed_fastq}"
TRIM_FRONT="${TRIM_FRONT:-10}"
TRIM_QUALITY="${TRIM_QUALITY:-28}"
MIN_BASE_QUALITY="${MIN_BASE_QUALITY:-30}"

# Ensure output directory exists
mkdir -p "$TRIM_OUTPUT"

echo "[INFO] fastp Stage 02: Adapter and Quality Trimming"
echo "[INFO] Input directory: $FASTQ_DIR"
echo "[INFO] Output directory: $TRIM_OUTPUT"
echo "[INFO] Parameters: trim_front=$TRIM_FRONT, cut_tail_mean_quality=$TRIM_QUALITY, min_base_quality=$MIN_BASE_QUALITY"
echo ""

# Find paired-end files (assuming naming: *_1.fastq/*_2.fastq or *_R1/*_R2)
declare -a PE_R1_FILES
while IFS= read -r file; do
  if [[ "$file" =~ (_1\.fastq|_R1\.fastq|_1\.fq|_R1\.fq) ]]; then
    PE_R1_FILES+=("$file")
  fi
done < <(find "$FASTQ_DIR" -maxdepth 1 \( -name "*_1.fastq*" -o -name "*_R1.fastq*" \) 2>/dev/null | sort)

# Process paired-end files
if [ ${#PE_R1_FILES[@]} -gt 0 ]; then
  echo "[INFO] Processing ${#PE_R1_FILES[@]} paired-end sample(s)..."
  for r1_file in "${PE_R1_FILES[@]}"; do
    # Infer R2 file
    r2_file="${r1_file/_1/_2}"
    r2_file="${r2_file/_R1/_R2}"

    if [ ! -f "$r2_file" ]; then
      echo "[WARNING] R2 file not found for $r1_file (expected: $r2_file), skipping"
      continue
    fi

    # Extract sample name
    sample_name=$(basename "$r1_file" | sed 's/_[12R].*\.\(fastq\|fq\).*//')

    echo "[INFO] Processing: $sample_name"

    # Run fastp for paired-end
    fastp \
      -i "$r1_file" \
      -I "$r2_file" \
      -o "$TRIM_OUTPUT/${sample_name}_1_trimmed.fastq.gz" \
      -O "$TRIM_OUTPUT/${sample_name}_2_trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      -f "$TRIM_FRONT" \
      -F "$TRIM_FRONT" \
      -3 \
      --cut_tail_mean_quality "$TRIM_QUALITY" \
      -q "$MIN_BASE_QUALITY" \
      -w "$THREADS" \
      -h "$TRIM_OUTPUT/${sample_name}_fastp.html" \
      -j "$TRIM_OUTPUT/${sample_name}_fastp.json"
  done
fi

# Find single-end files (excluding R1/R2 paired-end files)
declare -a SE_FILES
while IFS= read -r file; do
  # Skip if it looks like a paired-end file
  if [[ ! "$file" =~ (_1\.|_R1\.|_2\.|_R2\.) ]]; then
    SE_FILES+=("$file")
  fi
done < <(find "$FASTQ_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) 2>/dev/null | sort)

# Process single-end files
if [ ${#SE_FILES[@]} -gt 0 ]; then
  echo "[INFO] Processing ${#SE_FILES[@]} single-end sample(s)..."
  for se_file in "${SE_FILES[@]}"; do
    sample_name=$(basename "$se_file" | sed 's/\.\(fastq\|fq\).*//')

    echo "[INFO] Processing: $sample_name"

    # Run fastp for single-end
    fastp \
      -i "$se_file" \
      -o "$TRIM_OUTPUT/${sample_name}_trimmed.fastq.gz" \
      -f "$TRIM_FRONT" \
      -3 \
      --cut_tail_mean_quality "$TRIM_QUALITY" \
      -q "$MIN_BASE_QUALITY" \
      -w "$THREADS" \
      -h "$TRIM_OUTPUT/${sample_name}_fastp.html" \
      -j "$TRIM_OUTPUT/${sample_name}_fastp.json"
  done
fi

echo "[INFO] fastp Stage 02 complete"
echo "[INFO] Trimmed files saved to: $TRIM_OUTPUT"
