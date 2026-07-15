#!/usr/bin/env bash
# Stage 06: Transcript assembly and merge (StringTie)
# Assembles transcripts per sample using guided mode, then merges all samples

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
GTF_FILE="${GTF_FILE:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.115.gtf}"
BAM_OUTPUT="${BAM_OUTPUT:-results/05_hisat2_bam}"
THREADS="${THREADS:-8}"
STRINGTIE_OUTPUT="${STRINGTIE_OUTPUT:-results/06_stringtie_assembly}"
STRINGTIE_MERGED="${STRINGTIE_MERGED:-results/06_stringtie_assembly/merged.gtf}"

# Ensure output directory exists
mkdir -p "$STRINGTIE_OUTPUT"

echo "[INFO] StringTie Stage 06: Transcript Assembly and Merge"
echo "[INFO] Reference GTF: $GTF_FILE"
echo "[INFO] BAM input: $BAM_OUTPUT"
echo "[INFO] Output directory: $STRINGTIE_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo ""

# Check if GTF exists
if [ ! -f "$GTF_FILE" ]; then
  echo "[ERROR] GTF file not found: $GTF_FILE" >&2
  echo "[ERROR] Run 00_download_references.sh first or provide --gtf flag" >&2
  exit 1
fi

# Find all BAM files
BAM_FILES=$(find "$BAM_OUTPUT" -maxdepth 1 -name "*.bam" 2>/dev/null | sort)

if [ -z "$BAM_FILES" ]; then
  echo "[ERROR] No BAM files found in $BAM_OUTPUT" >&2
  echo "[ERROR] Did stage 05 (samtools) complete successfully?" >&2
  exit 1
fi

FILE_COUNT=$(echo "$BAM_FILES" | wc -l)
echo "[INFO] Found $FILE_COUNT BAM file(s) to process"
echo ""

# Phase 1: Assemble transcripts per sample (guided mode)
echo "[INFO] Phase 1: Per-sample transcript assembly..."
declare -a GTF_LIST
for bam_file in $BAM_FILES; do
  sample_name=$(basename "$bam_file" .bam)
  output_gtf="$STRINGTIE_OUTPUT/${sample_name}.gtf"

  echo "[INFO] Assembling $sample_name..."

  stringtie \
    -p "$THREADS" \
    -G "$GTF_FILE" \
    -o "$output_gtf" \
    "$bam_file"

  GTF_LIST+=("$output_gtf")
done

echo "[INFO] Per-sample assembly complete"
echo ""

# Phase 2: Merge all sample GTFs into a single reference transcriptome
echo "[INFO] Phase 2: Merging all sample GTFs..."

# Create a temporary file with list of GTF files
GTF_MERGE_LIST=$(mktemp)
trap "rm -f $GTF_MERGE_LIST" EXIT

printf '%s\n' "${GTF_LIST[@]}" > "$GTF_MERGE_LIST"

stringtie \
  --merge \
  -G "$GTF_FILE" \
  -p "$THREADS" \
  -o "$STRINGTIE_MERGED" \
  "$GTF_MERGE_LIST"

echo "[INFO] Merged GTF saved to: $STRINGTIE_MERGED"
echo "[INFO] StringTie Stage 06 complete"
