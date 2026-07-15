#!/usr/bin/env bash
# Stage 12: lncRNA-miRNA interaction prediction (IntaRNA)
# Runs IntaRNA to predict binding interactions between lncRNA candidates
# and miRNA sequences. This is the final bash stage before R-based ceRNA integration.
#
# Note: This stage depends on R script 10 (extract_lncRNA_for_intarna.R)
# being run first to generate the lncRNA FASTA candidates.
#
# Runs in isolated conda environment (env_intarna)

set -euo pipefail

# Source conda activation helper (uses separate env_intarna)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_activate_env.sh" env_intarna

# Configuration (from environment or defaults)
INTARNA_CANDIDATES_DIR="${INTARNA_CANDIDATES_DIR:-results/10_intarna_candidates}"  # output from R script 10
MIRNA_FASTA_DIR="${MIRNA_FASTA_DIR:-results/11_mirge3_output}"  # miRge3.0 output (contains miRNA FASTA)
INTARNA_OUTPUT="${INTARNA_OUTPUT:-results/12_intarna_interactions}"

# Ensure output directory exists
mkdir -p "$INTARNA_OUTPUT"

echo "[INFO] IntaRNA Stage 12: lncRNA-miRNA Interaction Prediction"
echo "[INFO] lncRNA candidates: $INTARNA_CANDIDATES_DIR"
echo "[INFO] miRNA FASTA: $MIRNA_FASTA_DIR"
echo "[INFO] Output directory: $INTARNA_OUTPUT"
echo ""

# Check if lncRNA candidates directory exists
if [ ! -d "$INTARNA_CANDIDATES_DIR" ]; then
  echo "[ERROR] lncRNA candidates directory not found: $INTARNA_CANDIDATES_DIR" >&2
  echo "[ERROR] Ensure R script 10 (extract_lncRNA_for_intarna.R) was run successfully" >&2
  exit 1
fi

# Find all lncRNA FASTA candidate files
# Expected naming pattern from R script: lncRNA_candidates_{timepoint}_{direction}.fa
# Timepoints: 6h, 24h; Directions: UP, DOWN

CANDIDATE_FILES=$(find "$INTARNA_CANDIDATES_DIR" -maxdepth 1 -name "*.fa" 2>/dev/null | sort)

if [ -z "$CANDIDATE_FILES" ]; then
  echo "[ERROR] No lncRNA candidate FASTA files found in $INTARNA_CANDIDATES_DIR" >&2
  exit 1
fi

echo "[INFO] Found lncRNA candidate files:"
echo "$CANDIDATE_FILES" | sed 's/^/  /'
echo ""

# Process each lncRNA candidate file
# For each lncRNA set, find corresponding miRNA FASTA files
# (from miRge3.0 output: miRNAs_UP_*.fa, miRNAs_DOWN_*.fa)

for lncrna_file in $CANDIDATE_FILES; do
  lncrna_basename=$(basename "$lncrna_file" .fa)

  # Extract timepoint and direction from filename
  # Expected: lncRNA_candidates_6h_UP.fa → timepoint=6h, direction=UP
  if [[ "$lncrna_basename" =~ ([0-9]+h)_([A-Z]+) ]]; then
    timepoint="${BASH_REMATCH[1]}"
    lncrna_direction="${BASH_REMATCH[2]}"
  else
    echo "[WARNING] Could not parse timepoint/direction from: $lncrna_basename"
    continue
  fi

  # For ceRNA logic: lncRNA UP pairs with miRNA DOWN, lncRNA DOWN pairs with miRNA UP
  if [ "$lncrna_direction" = "UP" ]; then
    mirna_direction="DOWN"
  else
    mirna_direction="UP"
  fi

  # Find corresponding miRNA FASTA file(s)
  # Try multiple naming patterns that might come from miRge3.0 output
  mirna_files=""
  for pattern in "miRNAs_${mirna_direction}_${timepoint}*.fa" "miRNAs_${mirna_direction}_${timepoint}*.fasta"; do
    matches=$(find "$MIRNA_FASTA_DIR" -maxdepth 1 -name "$pattern" 2>/dev/null | head -1)
    if [ -n "$matches" ]; then
      mirna_files="$matches"
      break
    fi
  done

  if [ -z "$mirna_files" ]; then
    echo "[WARNING] No miRNA FASTA found for timepoint=$timepoint direction=$mirna_direction"
    echo "[WARNING] Skipping IntaRNA for: $lncrna_basename"
    continue
  fi

  mirna_basename=$(basename "$mirna_files" .fa)
  output_name="Interact_${lncrna_basename}_${mirna_basename}"

  echo "[INFO] Running IntaRNA: lncRNA=$lncrna_basename vs miRNA=$mirna_basename"

  # Run IntaRNA with mode M (mfold-based)
  IntaRNA \
    --mode=M \
    -t "$lncrna_file" \
    -q "$mirna_files" \
    --outMode=C \
    --out "$INTARNA_OUTPUT/${output_name}.csv"

  echo "[INFO]   Output: $INTARNA_OUTPUT/${output_name}.csv"
done

echo "[INFO] IntaRNA Stage 12 complete"
echo "[INFO] Interaction predictions saved to: $INTARNA_OUTPUT"
echo "[INFO] Next: R script 11 (ceRNA_coherent_pairs.R) integrates IntaRNA results with DESeq2"
