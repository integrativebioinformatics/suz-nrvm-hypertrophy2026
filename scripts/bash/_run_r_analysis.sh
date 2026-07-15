#!/bin/bash
# R Analysis Pipeline Orchestrator
# Runs all 14 R analysis scripts in dependency order with automatic conda activation.
#
# External figure-input tables (GRN edges, KEGG GSEA tables, miRWalk targets,
# lncRNA-mRNA correlation tables) are expected under metadata/ and must be
# provided by the coauthor (see README "figure-input tables still to provide").
# Scripts that consume them degrade gracefully (skip the affected panel) when
# a file is absent.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_SCRIPTS_DIR="${SCRIPT_DIR}/R"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${1:-results}"
DPI="${2:-600}"
META="${REPO_ROOT}/metadata"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

source "$(conda info --base)/etc/profile.d/conda.sh" || { log_error "Failed to initialize conda"; exit 1; }
conda activate r || { log_error "Failed to activate 'r' env (conda env create -f envs/environment_r.yml)"; exit 1; }

log_info "R analysis pipeline started (results=${RESULTS_DIR}, dpi=${DPI})"

run_r_script() {
  local num="$1"; local name="$2"; shift 2
  log_info "==== Stage R${num}: ${name} ===="
  [[ -f "${R_SCRIPTS_DIR}/${num}_${name}.R" ]] || { log_error "Script not found: ${num}_${name}.R"; return 1; }
  Rscript "${R_SCRIPTS_DIR}/${num}_${name}.R" "$@" || { log_error "Failed at stage R${num}"; return 1; }
  log_info "Stage R${num} complete"; echo ""
}

DE01="${RESULTS_DIR}/01_deseq2_analysis"
MI02="${RESULTS_DIR}/02_mirna_deseq2"
MERGED_GTF="${RESULTS_DIR}/06_stringtie_assembly/merged.annotated.gtf"
FEELNC_OUT="${RESULTS_DIR}/09_feelnc_classify/feelnc_classes.txt"

# ============================================================================
# R01: DESeq2 mRNA/lncRNA + biotype split (REQUIRED — most scripts depend on it)
# ============================================================================
# Curated annotation files (original inputs; produced by the metadata generators).
run_r_script "01" "deseq2_mRNA_lncRNA" \
  --gtf "${MERGED_GTF}" \
  --salmon-dir "${RESULTS_DIR}/10_salmon_quant" \
  --sample-sheet "${META}/longRNA_samples.tsv" \
  --ann-annotated "${META}/hisat2_6_24_meta_ann.csv" \
  --ann-novel "${META}/hisat2_6_24_meta_novel_lncRNA.csv" \
  --output-dir "${DE01}" \
  --dpi "${DPI}"

# ============================================================================
# R02: miRNA DESeq2 (add --combat if smallRNA sample sheet has a `batch` column)
# ============================================================================
run_r_script "02" "mirna_deseq2" \
  --counts-dir "${RESULTS_DIR}/11_mirge3_output" \
  --sample-sheet "${META}/smallRNA_samples.tsv" \
  --output-dir "${MI02}" \
  --dpi "${DPI}"

# ============================================================================
# R10: Extract lncRNA candidates for IntaRNA (MUST run before bash stage 12)
# ============================================================================
run_r_script "10" "extract_lncRNA_for_intarna" \
  --lnc-dir "${DE01}" \
  --de-dir "${DE01}" \
  --fasta-file "${RESULTS_DIR}/08_gffread_fasta/merged_transcripts.fasta" \
  --output-dir "${RESULTS_DIR}/10_intarna_candidates"

log_warn "CHECKPOINT: now run bash stages 11-12 (miRge + IntaRNA), then re-run with --continue:"
log_warn "  ./suz_pipeline.bash --stages 11,12"
log_warn "  bash scripts/bash/_run_r_analysis.sh ${RESULTS_DIR} ${DPI} --continue"
if [[ "${3:-}" != "--continue" ]]; then
  log_info "Pausing before figure stages."; conda deactivate; exit 0
fi

# ============================================================================
# R03: mRNA temporal programs (Figure 3)
# ============================================================================
run_r_script "03" "figure3_mRNA_programs" \
  --de-dir "${DE01}" --output-dir "${RESULTS_DIR}/03_figure3" --dpi "${DPI}" --top-per-phase 12

# ============================================================================
# R05: lncRNA volcano plots (Figure 4)
# ============================================================================
run_r_script "05" "figure4_volcanoes" \
  --de-dir "${DE01}" --output-dir "${RESULTS_DIR}/05_figure4_volcanoes" --dpi "${DPI}"

# ============================================================================
# R06: lncRNA structural characterization (Supplementary)
# ============================================================================
run_r_script "06" "lncRNA_characterization" \
  --gtf "${MERGED_GTF}" --output-dir "${RESULTS_DIR}/06_lncRNA_characterization" --dpi "${DPI}"

# ============================================================================
# R07: FEELnc classification panels (Supplementary)
# ============================================================================
run_r_script "07" "feelnc_supplement_panels" \
  --feelnc-file "${FEELNC_OUT}" \
  --mrna-de-dir "${DE01}" \
  --output-dir "${RESULTS_DIR}/07_feelnc_supplement_panels" --dpi "${DPI}"

# ============================================================================
# R08: TF-KEGG integration + lncRNA panels + hubs (Figures 4/5/6)
#      Produces pathways_selected_for_TF.csv + hub/count tables used downstream.
# ============================================================================
TF08="${RESULTS_DIR}/08_tf_kegg"
run_r_script "08" "tf_kegg_integration" \
  --edge-6h "${META}/GRN_6h_edges.csv" \
  --edge-24h "${META}/GRN_24h_edges.csv" \
  --edge-merged "${META}/GRN_merged_edges.csv" \
  --kegg-6h "${META}/KEGG_Only6h.txt" \
  --kegg-24h "${META}/KEGG_Only24h.txt" \
  --kegg-shared "${META}/KEGG_Shared.txt" \
  --de-6h "${DE01}/NE_vs_Ctrl_6h_mRNA_DE.txt" \
  --de-24h "${DE01}/NE_vs_Ctrl_24h_mRNA_DE.txt" \
  --early "${DE01}/early_mRNA.csv" \
  --sustained "${DE01}/sustained_mRNA.csv" \
  --late "${DE01}/late_mRNA.csv" \
  --lnc-dir "${DE01}" \
  --cor-dir "${META}/correlations" \
  --output-dir "${TF08}" --dpi "${DPI}"

# ============================================================================
# R04: lncRNA heatmap + hubs + rewiring (Figure 4).
#      Hub/rewiring panels use co-expression tables produced by R08 / provided
#      externally; copy them into ${DE01} to enable, else those panels skip.
# ============================================================================
cp -f "${TF08}/lncRNA_hubs_by_time.csv" "${DE01}/" 2>/dev/null || true
run_r_script "04" "figure4_lncRNA_heatmap_hubs" \
  --de-dir "${DE01}" --output-dir "${RESULTS_DIR}/04_figure4_lncRNA_heatmap_hubs" --dpi "${DPI}"

# ============================================================================
# R09: miRNA-mRNA reciprocal integration (Figure 6). Uses R08 pathway selection.
# ============================================================================
MIR09="${RESULTS_DIR}/09_mirna_mRNA"
run_r_script "09" "mirna_mRNA_integration" \
  --mirna-de-dir "${MI02}" \
  --mrna-de-dir "${DE01}" \
  --mirwalk-dir "${META}/miRWalk" \
  --kegg-6h "${META}/KEGG_Only6h.txt" \
  --kegg-24h "${META}/KEGG_Only24h.txt" \
  --kegg-shared "${META}/KEGG_Shared.txt" \
  --pathways-selected "${TF08}/pathways_selected_for_TF.csv" \
  --phase-dir "${DE01}" \
  --output-dir "${MIR09}" --dpi "${DPI}"

# ============================================================================
# R11: ceRNA coherent pairs & triplets (Figure 7; lncRNA-only)
# ============================================================================
run_r_script "11" "ceRNA_coherent_pairs" \
  --intarna-dir "${RESULTS_DIR}/12_intarna_interactions" \
  --lnc-annotation "${RESULTS_DIR}/10_intarna_candidates/lncRNA_DE_with_phase_and_transcriptID.csv" \
  --mirna-de-dir "${MI02}" \
  --mirwalk-dir "${META}/miRWalk" \
  --mrna-de-dir "${DE01}" \
  --output-dir "${RESULTS_DIR}/11_ceRNA" --dpi "${DPI}"

# ============================================================================
# R12: Curated main-text panels (Fig4A/4B/6C)
# ============================================================================
run_r_script "12" "curated_fig4a_fig6c" \
  --tf-kegg-file "${TF08}/SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv" \
  --mir24-down "${MIR09}/SUPP_edges_24h_miRNADown_mRNAUp_KEGG_phase.csv" \
  --mir24-up "${MIR09}/SUPP_edges_24h_miRNAUp_mRNADown_KEGG_phase.csv" \
  --output-dir "${RESULTS_DIR}/12_curated" --dpi "${DPI}"

# ============================================================================
# R13: Redraw remaining panels (Fig4A/4B, Fig5D/5E, Fig6C + SuppFig4)
# ============================================================================
run_r_script "13" "redraw_remaining_panels" \
  --tf-dir "${TF08}" --mir-dir "${MIR09}" \
  --output-dir "${RESULTS_DIR}/13_final_panels" --dpi "${DPI}"

# ============================================================================
# R14: Master results report (data-driven)
# ============================================================================
run_r_script "14" "results_report" \
  --results-root "${RESULTS_DIR}" \
  --output-file "${RESULTS_DIR}/MASTER_results_report.txt"

log_info "==== R Analysis Pipeline Complete ===="
log_info "Report: ${RESULTS_DIR}/MASTER_results_report.txt"
conda deactivate
exit 0
