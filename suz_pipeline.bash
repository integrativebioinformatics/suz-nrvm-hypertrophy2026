#!/usr/bin/env bash
# NE_NRVM Reproducibility Pipeline — Top-Level Orchestrator
#
# This script orchestrates the entire 12-stage pipeline from raw FASTQ to
# quantified transcriptomes (Salmon), miRNA profiles (miRge3.0), and ceRNA
# networks (IntaRNA + R analysis).
#
# It activates appropriate conda environments automatically and runs stages
# sequentially, passing configuration flags to each stage.

set -euo pipefail

# ============================================================================
# Defaults and Configuration
# ============================================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$REPO_ROOT/scripts/bash"
R_SCRIPTS_DIR="$REPO_ROOT/scripts/R"

# Default file paths
LONGRNA_SAMPLES_DEFAULT="$REPO_ROOT/metadata/longRNA_samples.tsv"
SMALLRNA_SAMPLES_DEFAULT="$REPO_ROOT/metadata/smallRNA_samples.tsv"
FASTQ_MANIFEST_DEFAULT="$REPO_ROOT/metadata/fastq_manifest_template.tsv"
GENOME_FASTA_DEFAULT="$REPO_ROOT/references/Rattus_norvegicus.GRCr8.dna.toplevel.fa"
GTF_DEFAULT="$REPO_ROOT/references/Rattus_norvegicus.GRCr8.115.gtf"

# Default runtime parameters
THREADS_DEFAULT=8
OUTPUT_DIR_DEFAULT="$REPO_ROOT/results"
DPI_DEFAULT=600

# Placeholders for user-provided values
LONGRNA_SAMPLES=""
SMALLRNA_SAMPLES=""
FASTQ_MANIFEST=""
GENOME_FASTA="$GENOME_FASTA_DEFAULT"
GTF_FILE="$GTF_DEFAULT"
TFLINK_FILE=""
THREADS="$THREADS_DEFAULT"
OUTPUT_DIR="$OUTPUT_DIR_DEFAULT"
DPI="$DPI_DEFAULT"
STAGES_TO_RUN="all"
SKIP_DOWNLOAD_REFS=false
RUN_R_ANALYSIS=false

# ============================================================================
# Help Function
# ============================================================================

usage() {
  cat <<EOF
usage: $0 [OPTIONS]

NE_NRVM Reproducibility Pipeline

This script orchestrates the complete 12-stage pipeline for NE/Control NRVM
RNA-seq analysis:
  Stages 01–10: Long RNA (QC → quantification via Salmon)
  Stage 11: Small RNA (miRNA via miRge3.0)
  Stage 12: ceRNA predictions (lncRNA-miRNA interactions via IntaRNA)
  Stages R: R-based analysis (DESeq2, visualization, integration)

OPTIONS:
  --longrna-samples FILE       Path to long RNA sample sheet (default: metadata/longRNA_samples.tsv)
  --smallrna-samples FILE      Path to small RNA sample sheet (default: metadata/smallRNA_samples.tsv)
  --fastq-manifest FILE        Path to FASTQ manifest (default: metadata/fastq_manifest_template.tsv)
  --genome-fasta FILE          Path to reference genome FASTA
                               (default: references/Rattus_norvegicus.GRCr8.dna.toplevel.fa)
  --gtf FILE                   Path to reference GTF annotation
                               (default: references/Rattus_norvegicus.GRCr8.115.gtf)
  --tflink FILE                Path to TFLink All interactions TSV (optional)
  --threads N                  Number of threads (default: $THREADS_DEFAULT)
  --out DIR                    Output directory (default: results/)
  --stages SPEC                Run only specific stages (e.g., "1,2,3" or "all"; default: all)
  --skip-download-refs         Skip reference download (use if refs already exist)
  --run-r-analysis             Automatically run all 14 R analysis scripts after bash stages
  --dpi N                      PNG resolution for figure output (default: 600)
  --help                       Show this help message

EXAMPLES:
  # Full pipeline with defaults
  $0

  # Custom paths
  $0 \\
    --longrna-samples my_samples.tsv \\
    --genome-fasta /data/custom_genome.fasta \\
    --threads 16 \\
    --out my_results

  # Run only QC and trimming stages
  $0 --stages "1,2,3"

ENVIRONMENT:
  Four conda environments are required (created via envs/environment_*.yml):
    - env_qc_align: QC, trimming, alignment, assembly, quantification
    - env_mirge: miRNA quantification (miRge3.0)
    - env_intarna: lncRNA-miRNA interactions
    - env_r: R-based analysis

  Environments are automatically activated per stage (no manual activation needed).

EOF
  exit 1
}

# ============================================================================
# Parse Command-Line Arguments
# ============================================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    --longrna-samples)
      LONGRNA_SAMPLES="$2"
      shift 2
      ;;
    --smallrna-samples)
      SMALLRNA_SAMPLES="$2"
      shift 2
      ;;
    --fastq-manifest)
      FASTQ_MANIFEST="$2"
      shift 2
      ;;
    --genome-fasta)
      GENOME_FASTA="$2"
      shift 2
      ;;
    --gtf)
      GTF_FILE="$2"
      shift 2
      ;;
    --tflink)
      TFLINK_FILE="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --out)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --stages)
      STAGES_TO_RUN="$2"
      shift 2
      ;;
    --skip-download-refs)
      SKIP_DOWNLOAD_REFS=true
      shift
      ;;
    --run-r-analysis)
      RUN_R_ANALYSIS=true
      shift
      ;;
    --dpi)
      DPI="$2"
      shift 2
      ;;
    --help)
      usage
      ;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      ;;
  esac
done

# ============================================================================
# Validate and Set Defaults
# ============================================================================

# Use defaults if not provided
LONGRNA_SAMPLES="${LONGRNA_SAMPLES:-$LONGRNA_SAMPLES_DEFAULT}"
SMALLRNA_SAMPLES="${SMALLRNA_SAMPLES:-$SMALLRNA_SAMPLES_DEFAULT}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-$FASTQ_MANIFEST_DEFAULT}"

echo ""
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║ NE_NRVM Reproducibility Pipeline                                      ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "[CONFIG] Repository root: $REPO_ROOT"
echo "[CONFIG] Output directory: $OUTPUT_DIR"
echo "[CONFIG] Threads: $THREADS"
echo "[CONFIG] Genome FASTA: $GENOME_FASTA"
echo "[CONFIG] GTF: $GTF_FILE"
echo "[CONFIG] Long RNA samples: $LONGRNA_SAMPLES"
echo "[CONFIG] Small RNA samples: $SMALLRNA_SAMPLES"
echo "[CONFIG] FASTQ manifest: $FASTQ_MANIFEST"
if [ -n "$TFLINK_FILE" ]; then
  echo "[CONFIG] TFLink: $TFLINK_FILE"
fi
echo ""

# ============================================================================
# Pre-flight Checks
# ============================================================================

echo "[CHECK] Verifying sample sheets and manifests..."

for f in "$LONGRNA_SAMPLES" "$SMALLRNA_SAMPLES" "$FASTQ_MANIFEST"; do
  if [ ! -f "$f" ]; then
    echo "[ERROR] Required file not found: $f" >&2
    exit 1
  fi
done

echo "[CHECK] Sample sheets verified ✓"
echo ""

# ============================================================================
# Stage Execution Function
# ============================================================================

run_stage() {
  local stage_num="$1"
  local stage_script="$SCRIPTS_DIR/${stage_num}*.sh"
  local stage_file=$(ls "$stage_script" 2>/dev/null | head -1)

  if [ -z "$stage_file" ]; then
    echo "[ERROR] Stage script not found: $stage_script" >&2
    return 1
  fi

  echo ""
  echo "════════════════════════════════════════════════════════════════════════"
  echo "Stage $(printf "%02d" "$stage_num"): $(basename "$stage_file" .sh)"
  echo "════════════════════════════════════════════════════════════════════════"
  echo ""

  # Export environment variables for the stage script to use
  export THREADS GENOME_FASTA GTF_FILE OUTPUT_DIR TFLINK_FILE

  # Export standard output paths
  export QC_OUTPUT="${OUTPUT_DIR}/01_fastqc_raw"
  export TRIM_OUTPUT="${OUTPUT_DIR}/02_trimmed_fastq"
  export QC_TRIMMED="${OUTPUT_DIR}/03_fastqc_trimmed"
  export ALIGN_OUTPUT="${OUTPUT_DIR}/04_hisat2_sam"
  export BAM_OUTPUT="${OUTPUT_DIR}/05_hisat2_bam"
  export STRINGTIE_OUTPUT="${OUTPUT_DIR}/06_stringtie_assembly"
  export STRINGTIE_MERGED="${OUTPUT_DIR}/06_stringtie_assembly/merged.gtf"
  export GFFCOMPARE_OUTPUT="${OUTPUT_DIR}/07_gffcompare"
  export FASTA_OUTPUT="${OUTPUT_DIR}/08_gffread_fasta"
  export FEELNC_OUTPUT="${OUTPUT_DIR}/09_feelnc_classify"
  export SALMON_OUTPUT="${OUTPUT_DIR}/10_salmon_quant"
  export MIRGE_OUTPUT="${OUTPUT_DIR}/11_mirge3_output"
  export INTARNA_OUTPUT="${OUTPUT_DIR}/12_intarna_interactions"

  # Run the stage
  bash "$stage_file" || {
    echo "[ERROR] Stage failed: $stage_file" >&2
    return 1
  }
}

# ============================================================================
# Download References (if needed)
# ============================================================================

if [ "$SKIP_DOWNLOAD_REFS" = false ]; then
  echo "[INFO] Checking references..."
  bash "$SCRIPTS_DIR/00_download_references.sh" "$GENOME_FASTA" "$GTF_FILE" || {
    echo "[ERROR] Reference download failed" >&2
    exit 1
  }
  echo ""
fi

# ============================================================================
# Execute Specified Stages
# ============================================================================

if [ "$STAGES_TO_RUN" = "all" ]; then
  STAGE_LIST="1 2 3 4 5 6 7 8 9 10 11 12"
else
  STAGE_LIST="$(echo "$STAGES_TO_RUN" | tr ',' ' ')"
fi

echo "[INFO] Stages to run: $STAGE_LIST"
echo ""

FAILED_STAGES=""

for stage in $STAGE_LIST; do
  if ! run_stage "$stage"; then
    echo "[ERROR] Stage $stage failed" >&2
    FAILED_STAGES="$FAILED_STAGES $stage"
  fi
done

# ============================================================================
# Final Report
# ============================================================================

echo ""
echo "════════════════════════════════════════════════════════════════════════"
echo "Pipeline Execution Complete"
echo "════════════════════════════════════════════════════════════════════════"
echo ""

if [ -z "$FAILED_STAGES" ]; then
  echo "[✓] All bash stages completed successfully!"
  echo ""

  # Auto-run R analysis if requested
  if [ "$RUN_R_ANALYSIS" = true ]; then
    echo "[INFO] Starting R analysis pipeline..."
    echo ""

    # Run R stages 01-02, 10 (pre-IntaRNA)
    bash "$SCRIPTS_DIR/_run_r_analysis.sh" "$OUTPUT_DIR" "$DPI" || {
      echo "[ERROR] R analysis pipeline failed" >&2
      exit 1
    }

    echo ""
    echo "════════════════════════════════════════════════════════════════════════"
    echo "Complete Pipeline Execution Finished Successfully!"
    echo "════════════════════════════════════════════════════════════════════════"
    echo ""
    echo "[✓] All bash stages (00-12) + R analysis (01-14) complete!"
    echo ""
    echo "Generated outputs in: $OUTPUT_DIR"
    echo "  - See MASTER_results_report.txt for comprehensive summary"
    exit 0
  else
    echo "[NEXT STEPS]"
    echo "  1. Review bash pipeline outputs in: $OUTPUT_DIR"
    echo "  2. Run R analysis scripts automatically:"
    echo "     $0 --run-r-analysis --out $OUTPUT_DIR"
    echo "  3. Or run manually:"
    echo "     bash $SCRIPTS_DIR/_run_r_analysis.sh $OUTPUT_DIR $DPI"
    echo "  4. See README.md for detailed usage and output descriptions"
    exit 0
  fi
else
  echo "[✗] Pipeline failed at stages:$FAILED_STAGES"
  echo "[TROUBLESHOOTING]"
  echo "  1. Check stage log files in: $OUTPUT_DIR"
  echo "  2. Verify conda environments are installed"
  echo "  3. Ensure all input files are accessible"
  echo "  4. Re-run failed stages individually for debugging"
  exit 1
fi
