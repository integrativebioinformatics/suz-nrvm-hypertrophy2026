#!/usr/bin/env bash
# Stage 08: Extract transcript sequences (gffread)
# Extracts FASTA sequences of all reconstructed transcripts from the genome
# Produces both genomic coordinates (cDNA) for lncRNA/novel transcripts

set -euo pipefail

# Source conda activation helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_qc_align

# Configuration (from environment or defaults)
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
GENOME_FASTA="${GENOME_FASTA:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.toplevel.fa}"
GFFCOMPARE_OUTPUT="${GFFCOMPARE_OUTPUT:-results/07_gffcompare}"
FASTA_OUTPUT="${FASTA_OUTPUT:-results/08_gffread_fasta}"

# Ensure output directory exists
mkdir -p "$FASTA_OUTPUT"

echo "[INFO] gffread Stage 08: Extract Transcript FASTA Sequences"
echo "[INFO] Reference genome: $GENOME_FASTA"
echo "[INFO] Annotated GTF: $GFFCOMPARE_OUTPUT/merged.annotated.gtf"
echo "[INFO] Output directory: $FASTA_OUTPUT"
echo ""

# Check if input files exist
if [ ! -f "$GENOME_FASTA" ]; then
  echo "[ERROR] Genome FASTA not found: $GENOME_FASTA" >&2
  exit 1
fi

ANNOTATED_GTF="$GFFCOMPARE_OUTPUT/merged.annotated.gtf"
if [ ! -f "$ANNOTATED_GTF" ]; then
  echo "[ERROR] Annotated GTF not found: $ANNOTATED_GTF" >&2
  echo "[ERROR] Did stage 07 (gffcompare) complete successfully?" >&2
  exit 1
fi

# Extract all transcript sequences (cDNA) from merged/annotated GTF
echo "[INFO] Extracting cDNA sequences..."

gffread \
  "$ANNOTATED_GTF" \
  -g "$GENOME_FASTA" \
  -w "$FASTA_OUTPUT/merged_transcripts.fasta"

echo "[INFO] cDNA FASTA extracted: $FASTA_OUTPUT/merged_transcripts.fasta"

# Optionally extract only novel/intergenic transcripts (class_code='u') for lncRNA analysis
echo "[INFO] Extracting novel intergenic transcripts only..."

# Create a temporary GTF with only novel transcripts (class_code "u")
NOVEL_GTF=$(mktemp --suffix=.gtf)
trap "rm -f $NOVEL_GTF" EXIT

grep 'class_code "u"' "$ANNOTATED_GTF" > "$NOVEL_GTF" || true

if [ -s "$NOVEL_GTF" ]; then
  gffread \
    "$NOVEL_GTF" \
    -g "$GENOME_FASTA" \
    -w "$FASTA_OUTPUT/merged_novel_lncRNA_candidates.fasta"

  NOVEL_COUNT=$(grep -c "^>" "$FASTA_OUTPUT/merged_novel_lncRNA_candidates.fasta" || echo "0")
  echo "[INFO] Novel lncRNA sequences extracted: $NOVEL_COUNT transcripts"
else
  echo "[WARNING] No novel transcripts found (no class_code 'u' entries)"
fi

echo "[INFO] gffread Stage 08 complete"
echo "[INFO] FASTA files saved to: $FASTA_OUTPUT"
