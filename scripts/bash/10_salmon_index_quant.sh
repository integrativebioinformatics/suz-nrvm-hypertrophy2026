#!/usr/bin/env bash
# Stage 10: Transcript quantification (Salmon)
# Builds a Salmon index from the merged transcriptome and quantifies
# all samples against it. This is the primary quantification method
# for the R DESeq2 analysis (via tximport).

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
FASTA_OUTPUT="${FASTA_OUTPUT:-results/08_gffread_fasta}"
TRIM_OUTPUT="${TRIM_OUTPUT:-results/02_trimmed_fastq}"
THREADS="${THREADS:-8}"
SALMON_OUTPUT="${SALMON_OUTPUT:-results/10_salmon_quant}"

# Ensure output directory exists
mkdir -p "$SALMON_OUTPUT"

echo "[INFO] Salmon Stage 10: Index and Quantify Transcriptome"
echo "[INFO] Transcriptome FASTA: $FASTA_OUTPUT/merged_transcripts.fasta"
echo "[INFO] Trimmed reads: $TRIM_OUTPUT"
echo "[INFO] Output directory: $SALMON_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo ""

# Check if transcriptome FASTA exists
TRANSCRIPTOME_FASTA="$FASTA_OUTPUT/merged_transcripts.fasta"
if [ ! -f "$TRANSCRIPTOME_FASTA" ]; then
  echo "[ERROR] Transcriptome FASTA not found: $TRANSCRIPTOME_FASTA" >&2
  echo "[ERROR] Did stage 08 (gffread) complete successfully?" >&2
  exit 1
fi

# Build Salmon index if it doesn't exist
SALMON_INDEX="$SALMON_OUTPUT/salmon_index"
if [ ! -d "$SALMON_INDEX" ]; then
  echo "[INFO] Building Salmon index..."

  salmon index \
    -t "$TRANSCRIPTOME_FASTA" \
    -i "$SALMON_INDEX" \
    -p "$THREADS" \
    -k 31

  echo "[INFO] Salmon index built: $SALMON_INDEX"
else
  echo "[INFO] Salmon index already exists, skipping build"
fi

echo ""
echo "[INFO] Quantifying samples with Salmon..."
echo ""

# Find and quantify all samples (paired-end and single-end)
declare -a PE_R1_FILES
while IFS= read -r file; do
  if [[ "$file" =~ _1_trimmed\.fastq ]]; then
    PE_R1_FILES+=("$file")
  fi
done < <(find "$TRIM_OUTPUT" -maxdepth 1 -name "*_1_trimmed.fastq*" 2>/dev/null | sort)

# Quantify paired-end samples
if [ ${#PE_R1_FILES[@]} -gt 0 ]; then
  echo "[INFO] Quantifying ${#PE_R1_FILES[@]} paired-end sample(s)..."
  for r1_file in "${PE_R1_FILES[@]}"; do
    r2_file="${r1_file/_1_trimmed/_2_trimmed}"

    if [ ! -f "$r2_file" ]; then
      echo "[WARNING] R2 file not found: $r2_file, skipping"
      continue
    fi

    sample_name=$(basename "$r1_file" | sed 's/_1_trimmed.*//')
    echo "[INFO] Quantifying $sample_name (paired-end)..."

    salmon quant \
      -i "$SALMON_INDEX" \
      -l A \
      -1 "$r1_file" \
      -2 "$r2_file" \
      -p "$THREADS" \
      --validateMappings \
      -o "$SALMON_OUTPUT/${sample_name}"
  done
fi

# Find and quantify single-end samples (not in paired-end naming)
declare -a SE_FILES
while IFS= read -r file; do
  if [[ ! "$file" =~ _[12]_trimmed\.fastq ]]; then
    SE_FILES+=("$file")
  fi
done < <(find "$TRIM_OUTPUT" -maxdepth 1 -name "*_trimmed.fastq*" 2>/dev/null | sort)

if [ ${#SE_FILES[@]} -gt 0 ]; then
  echo "[INFO] Quantifying ${#SE_FILES[@]} single-end sample(s)..."
  for se_file in "${SE_FILES[@]}"; do
    sample_name=$(basename "$se_file" | sed 's/_trimmed.*//')
    echo "[INFO] Quantifying $sample_name (single-end)..."

    salmon quant \
      -i "$SALMON_INDEX" \
      -l A \
      -r "$se_file" \
      -p "$THREADS" \
      --validateMappings \
      -o "$SALMON_OUTPUT/${sample_name}"
  done
fi

echo "[INFO] Salmon Stage 10 complete"
echo "[INFO] Quantification results in: $SALMON_OUTPUT"
echo "[INFO] Next: R analysis scripts will use these Salmon outputs via tximport"
