#!/usr/bin/env bash
# Stage 11: Small RNA/miRNA quantification (miRge3.0)
# Identifies and quantifies miRNAs from small RNA-seq libraries
# Uses miRBase database (confirmed final choice)
# Runs in isolated conda environment (env_mirge)

set -euo pipefail

# Source conda activation helper (uses separate env_mirge)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_mirge

# Configuration (from environment or defaults)
MIRGE_LIB="${MIRGE_LIB:-.}"  # path to miRge3.0 library (default: current dir, coauthor to provide)
SMALLRNA_FASTQ_DIR="${SMALLRNA_FASTQ_DIR:-.}"  # directory containing small RNA FASTQ files
MIRGE_OUTPUT="${MIRGE_OUTPUT:-results/11_mirge3_output}"
MIRGE_ORGANISM="${MIRGE_ORGANISM:-rat}"  # organism: rat (default) or custom
MIRGE_DB="${MIRGE_DB:-mirbase}"  # database: mirbase (confirmed final) or mirgenedb
THREADS="${THREADS:-8}"
ADAPTER_3="${ADAPTER_3:-AGATCGGAAGAGCACACGTCT}"  # 3' Illumina TruSeq adapter
ADAPTER_G="${ADAPTER_G:-GTTCAGAGTTCTACAGTCCGACGATC}"  # Additional adapter/UMI sequence
METADATA_CSV="${METADATA_CSV:-.}"  # optional metadata CSV for differential expression

# Ensure output directory exists
mkdir -p "$MIRGE_OUTPUT"

echo "[INFO] miRge3.0 Stage 11: Small RNA / miRNA Quantification"
echo "[INFO] Organism: $MIRGE_ORGANISM"
echo "[INFO] Database: $MIRGE_DB"
echo "[INFO] Output directory: $MIRGE_OUTPUT"
echo "[INFO] Threads: $THREADS"
echo ""

# Check for small RNA FASTQ files
if [ "$SMALLRNA_FASTQ_DIR" = "." ]; then
  echo "[WARNING] SMALLRNA_FASTQ_DIR not set; using current directory"
fi

SMALLRNA_FILES=$(find "$SMALLRNA_FASTQ_DIR" -maxdepth 1 \( -name "*_mi.fastq" -o -name "*_mi.fastq.gz" \) 2>/dev/null | sort)

if [ -z "$SMALLRNA_FILES" ]; then
  echo "[ERROR] No small RNA FASTQ files found (expected naming: *_mi.fastq*)" >&2
  echo "[ERROR] Set SMALLRNA_FASTQ_DIR to the directory containing small RNA libraries" >&2
  exit 1
fi

# Build comma-separated list of sample files for miRge3.0
SAMPLE_LIST=$(echo "$SMALLRNA_FILES" | tr '\n' ',' | sed 's/,$//')

echo "[INFO] Found small RNA samples:"
echo "$SMALLRNA_FILES" | sed 's/^/  /'
echo ""

echo "[INFO] Running miRge3.0..."
echo ""

# Build miRge3.0 command
MIRGE_CMD="miRge3.0 \
  -s $SAMPLE_LIST \
  -lib $MIRGE_LIB \
  -on $MIRGE_ORGANISM \
  -db $MIRGE_DB \
  -o $MIRGE_OUTPUT \
  -a '$ADAPTER_3' \
  -g '$ADAPTER_G' \
  -cpu $THREADS \
  -dex \
  -mEC"

# Add metadata if provided
if [ -n "$METADATA_CSV" ] && [ -f "$METADATA_CSV" ]; then
  MIRGE_CMD="$MIRGE_CMD -mdt $METADATA_CSV"
  echo "[INFO] Using metadata CSV for differential expression: $METADATA_CSV"
fi

echo "[INFO] Command: $MIRGE_CMD"
echo ""

# Execute miRge3.0
eval "$MIRGE_CMD"

echo "[INFO] miRge3.0 Stage 11 complete"
echo "[INFO] Results saved to: $MIRGE_OUTPUT"
echo "[INFO] Key outputs:"
echo "[INFO]   - miR.Counts_* (miRNA count matrices)"
echo "[INFO]   - De_* (differential expression results if -dex used)"
echo "[INFO]   - HTML report"
