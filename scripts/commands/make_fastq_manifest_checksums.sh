#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash make_fastq_manifest_checksums.sh /path/to/raw_fastq metadata/fastq_manifest_template.tsv
#
# This script creates two files:
#   checksums/fastq_md5_checksums.txt
#   checksums/fastq_file_sizes.tsv
#
# It does not modify your FASTQ files.

RAW_FASTQ_DIR="${1:-/path/to/raw_fastq}"
mkdir -p checksums

find "$RAW_FASTQ_DIR" -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) \
  -print0 | sort -z | xargs -0 md5sum > checksums/fastq_md5_checksums.txt

find "$RAW_FASTQ_DIR" -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) \
  -printf "%p\t%s\n" | sort > checksums/fastq_file_sizes.tsv

echo "Done."
echo "MD5 checksums: checksums/fastq_md5_checksums.txt"
echo "File sizes: checksums/fastq_file_sizes.tsv"
