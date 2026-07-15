#!/usr/bin/env bash
# Download default Ensembl references (rat genome + GTF) if they don't exist
# This is automatically called by suz_pipeline.bash unless --genome-fasta and --gtf flags are provided
#
# Parameters:
#   $1 = genome FASTA path (optional; defaults to references/Rattus_norvegicus.GRCr8.dna.toplevel.fa)
#   $2 = GTF path (optional; defaults to references/Rattus_norvegicus.GRCr8.115.gtf)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

GENOME_FASTA="${1:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.toplevel.fa}"
GTF_FILE="${2:-$REPO_ROOT/references/Rattus_norvegicus.GRCr8.115.gtf}"

# Ensembl FTP URLs (rat genome GRCr8, release 115)
ENSEMBL_FASTA_URL="https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz"
ENSEMBL_GTF_URL="https://ftp.ensembl.org/pub/release-115/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gtf.gz"

echo "[INFO] Reference download script starting..."
echo "[INFO] Target genome FASTA: $GENOME_FASTA"
echo "[INFO] Target GTF: $GTF_FILE"

# Function to download with curl or wget fallback
download_file() {
  local url="$1"
  local output="$2"

  if command -v curl &>/dev/null; then
    echo "[INFO] Downloading with curl: $url"
    curl -L -o "$output" "$url"
  elif command -v wget &>/dev/null; then
    echo "[INFO] Downloading with wget: $url"
    wget -O "$output" "$url"
  else
    echo "[ERROR] Neither curl nor wget found. Install one and try again." >&2
    return 1
  fi
}

# Download genome FASTA if not present
if [ ! -f "$GENOME_FASTA" ]; then
  echo "[INFO] Genome FASTA not found, downloading..."
  TEMP_FASTA="$GENOME_FASTA.gz"
  download_file "$ENSEMBL_FASTA_URL" "$TEMP_FASTA"
  echo "[INFO] Decompressing genome FASTA..."
  gunzip -v "$TEMP_FASTA"
  echo "[INFO] Genome FASTA ready: $GENOME_FASTA"
else
  echo "[INFO] Genome FASTA already exists, skipping download."
fi

# Download GTF if not present
if [ ! -f "$GTF_FILE" ]; then
  echo "[INFO] GTF not found, downloading..."
  TEMP_GTF="$GTF_FILE.gz"
  download_file "$ENSEMBL_GTF_URL" "$TEMP_GTF"
  echo "[INFO] Decompressing GTF..."
  gunzip -v "$TEMP_GTF"
  echo "[INFO] GTF ready: $GTF_FILE"
else
  echo "[INFO] GTF already exists, skipping download."
fi

echo "[INFO] Reference download complete."
