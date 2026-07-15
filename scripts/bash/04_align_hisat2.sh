#!/usr/bin/env bash
# Stage 04: Genome alignment (HISAT2)
# Aligns trimmed reads to the reference genome with --dta-cufflinks flag
# (suitable for StringTie transcript assembly)

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
GENOME_FASTA="${GENOME_FASTA:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.toplevel.fa}"
TRIM_OUTPUT="${TRIM_OUTPUT:-results/02_trimmed_fastq}"
THREADS="${THREADS:-8}"
ALIGN_OUTPUT="${ALIGN_OUTPUT:-results/04_hisat2_sam}"

# Ensure output directory exists
mkdir -p "$ALIGN_OUTPUT"

echo "[INFO] HISAT2 Stage 04: Genome Alignment"
echo "[INFO] Reference genome: $GENOME_FASTA"
echo "[INFO] Trimmed input: $TRIM_OUTPUT"
echo "[INFO] Output directory: $ALIGN_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo ""

# Check if genome FASTA exists
if [ ! -f "$GENOME_FASTA" ]; then
  echo "[ERROR] Genome FASTA not found: $GENOME_FASTA" >&2
  echo "[ERROR] Run 00_download_references.sh first or provide --genome-fasta flag" >&2
  exit 1
fi

# Build HISAT2 index if it doesn't exist
HISAT2_INDEX="${ALIGN_OUTPUT}/hisat2_index"
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
  echo "[INFO] Building HISAT2 index..."
  hisat2-build \
    -p "$THREADS" \
    "$GENOME_FASTA" \
    "$HISAT2_INDEX"
  echo "[INFO] HISAT2 index built"
fi

echo "[INFO] Aligning reads to genome..."
echo ""

# Find paired-end and single-end trimmed files and align
declare -a PE_R1_FILES
while IFS= read -r file; do
  if [[ "$file" =~ _1_trimmed\.fastq ]]; then
    PE_R1_FILES+=("$file")
  fi
done < <(find "$TRIM_OUTPUT" -maxdepth 1 -name "*_1_trimmed.fastq*" 2>/dev/null | sort)

# Process paired-end files
if [ ${#PE_R1_FILES[@]} -gt 0 ]; then
  echo "[INFO] Processing ${#PE_R1_FILES[@]} paired-end sample(s)..."
  for r1_file in "${PE_R1_FILES[@]}"; do
    r2_file="${r1_file/_1_trimmed/_2_trimmed}"

    if [ ! -f "$r2_file" ]; then
      echo "[WARNING] R2 file not found: $r2_file, skipping"
      continue
    fi

    sample_name=$(basename "$r1_file" | sed 's/_1_trimmed.*//')
    echo "[INFO] Aligning $sample_name (paired-end)..."

    hisat2 \
      --dta-cufflinks \
      -p "$THREADS" \
      -x "$HISAT2_INDEX" \
      -1 "$r1_file" \
      -2 "$r2_file" \
      -S "$ALIGN_OUTPUT/${sample_name}.sam" \
      --summary-file "$ALIGN_OUTPUT/${sample_name}_align_summary.txt"
  done
fi

# Find and process single-end trimmed files (not in paired-end R1/R2 naming)
declare -a SE_FILES
while IFS= read -r file; do
  if [[ ! "$file" =~ _[12]_trimmed\.fastq ]]; then
    SE_FILES+=("$file")
  fi
done < <(find "$TRIM_OUTPUT" -maxdepth 1 -name "*_trimmed.fastq*" 2>/dev/null | sort)

if [ ${#SE_FILES[@]} -gt 0 ]; then
  echo "[INFO] Processing ${#SE_FILES[@]} single-end sample(s)..."
  for se_file in "${SE_FILES[@]}"; do
    sample_name=$(basename "$se_file" | sed 's/_trimmed.*//')
    echo "[INFO] Aligning $sample_name (single-end)..."

    hisat2 \
      --dta-cufflinks \
      -p "$THREADS" \
      -x "$HISAT2_INDEX" \
      -U "$se_file" \
      -S "$ALIGN_OUTPUT/${sample_name}.sam" \
      --summary-file "$ALIGN_OUTPUT/${sample_name}_align_summary.txt"
  done
fi

echo "[INFO] HISAT2 Stage 04 complete"
echo "[INFO] SAM files saved to: $ALIGN_OUTPUT"
