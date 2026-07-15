# SUZ2026 pipeline

**Last updated:** 14 July 2026 — All 14 R analysis scripts fully implemented  

A Linux-first, fully parametrized RNA-seq analysis pipeline for **NE (norepinephrine) vs Control treatment** of **neonatal rat ventricular myocytes (NRVM)** at **6h and 24h timepoints** (n=4/group). 

Integrates:
- **Long RNA analysis** (mRNA + novel lncRNA discovery via StringTie + FEELnc)
- **Small RNA analysis** (miRNA via miRge3.0)
- **ceRNA networks** (lncRNA–miRNA–mRNA interactions via IntaRNA)
- **R-based analysis** (DESeq2, temporal clustering, pathway integration)

**Output:** 7-figure manuscript with mRNA programs, lncRNA dynamics, TF/KEGG integration, miRNA-mRNA pathways, and ceRNA network predictions.

---

## Quick Start

### 1. Create Conda Environments

```bash
# Create all four analysis environments
conda env create -f envs/environment_qc_align.yml
conda env create -f envs/environment_mirge.yml
conda env create -f envs/environment_intarna.yml
conda env create -f envs/environment_r.yml
```

Environments are auto-activated per-stage; no manual activation needed.

### 2. Prepare Inputs

Fill in the FASTQ manifest with actual file paths and checksums:
```bash
# Edit metadata/fastq_manifest.tsv with:
# - exact_fastq_filename
# - server_path or hard_drive_backup_path
# - file_size_bytes
# - md5_checksum
```

See `commands/make_fastq_manifest_checksums.sh` for help generating checksums from a server folder.

### 3. Run the Pipeline

#### Option A: Automatic (bash + R stages together)

```bash
./suz_pipeline.bash \
  --longrna-samples metadata/longRNA_samples.tsv \
  --smallrna-samples metadata/smallRNA_samples.tsv \
  --fastq-manifest metadata/fastq_manifest_template.tsv \
  --threads 16 \
  --out results \
  --run-r-analysis \
  --dpi 600
```

This runs the full pipeline end-to-end:

1. Auto-downloads references from Ensembl (or uses custom via `--genome-fasta`/`--gtf`)
2. Runs all 12 bash stages (QC → alignment → quantification → miRge → IntaRNA)
3. Automatically runs all 14 R analysis scripts (DESeq2 → figures → integration → report)

#### Option B: Manual (bash stages, then R stages separately)

```bash
# Run bash stages 00-10
./suz_pipeline.bash \
  --longrna-samples metadata/longRNA_samples.tsv \
  --smallrna-samples metadata/smallRNA_samples.tsv \
  --fastq-manifest metadata/fastq_manifest_template.tsv \
  --threads 16 \
  --out results \
  --stages 0-10

# [... complete bash stages 11-12 in parallel or at a later time ...]

# Then run all R analysis stages
bash scripts/bash/_run_r_analysis.sh results 600
```

#### Option C: Run only R analysis (bash stages already done)

```bash
bash scripts/bash/_run_r_analysis.sh results 600
```

---

## Pipeline Stages

### Bash Stages (Automatic conda env activation)

| Stage | Script | Environment | Input | Output | Tools |
|-------|--------|-------------|-------|--------|-------|
| **Prep** | `00_download_references.sh` | qc_align | (auto) | `references/*.fasta`, `*.gtf` | curl/wget |
| **01** | `01_qc_raw.sh` | qc_align | Raw FASTQ | `01_fastqc_raw/` | FastQC |
| **02** | `02_trim_fastp.sh` | qc_align | Raw FASTQ | `02_trimmed_fastq/` | fastp |
| **03** | `03_qc_trimmed.sh` | qc_align | Trimmed FASTQ | `03_fastqc_trimmed/` | FastQC |
| **04** | `04_align_hisat2.sh` | qc_align | Trimmed FASTQ | `04_hisat2_sam/` | HISAT2 |
| **05** | `05_sort_index_bam.sh` | qc_align | SAM files | `05_hisat2_bam/` | samtools |
| **06** | `06_stringtie_assemble_merge.sh` | qc_align | BAM files | `06_stringtie_assembly/merged.gtf` | StringTie |
| **07** | `07_gffcompare.sh` | qc_align | Merged GTF + ref GTF | `07_gffcompare/merged.annotated.gtf` | gffcompare |
| **08** | `08_gffread_extract_fasta.sh` | qc_align | Annotated GTF + genome | `08_gffread_fasta/` | gffread |
| **09** | `09_feelnc_classify.sh` | qc_align | GTF + genome | `09_feelnc_classify/` | FEELnc |
| **10** | `10_salmon_index_quant.sh` | qc_align | Trimmed FASTQ + GTF | `10_salmon_quant/` | Salmon |
| **11** | `11_mirge3_run.sh` | mirge | Small RNA FASTQ | `11_mirge3_output/` | miRge3.0 |
| **12** | `12_intarna_run.sh` | intarna | lncRNA + miRNA FASTA | `12_intarna_interactions/` | IntaRNA |

### R Analysis Stages (Manual execution, or chain via bash orchestrator)

| Stage | Script | Status | Input | Output | Purpose |
|-------|--------|--------|-------|--------|---------|
| **01** | `01_deseq2_mRNA_lncRNA.R` | ✓ Implemented | Salmon `quant.sf` | DE tables, Mfuzz clusters | Bulk mRNA/lncRNA DESeq2 + temporal clustering |
| **02** | `02_mirna_deseq2.R` | ✓ Implemented | miRge count matrices | miRNA DE tables, volcano plots | miRNA DESeq2 |
| **03** | `03_figure3_mRNA_programs.R` | ✓ Implemented | mRNA DE + phases | Figure 3 panels | mRNA temporal programs (phase classification, heatmaps) |
| **04** | `04_figure4_lncRNA_heatmap_hubs.R` | ✓ Implemented | lncRNA DE + correlations | Figure 4 heatmaps | lncRNA hubs & co-expression dynamics |
| **05** | `05_figure4_volcanoes.R` | ✓ Implemented | lncRNA DE (6h/24h) | Figure 4 volcanos | lncRNA volcano plots with phase overlay |
| **06** | `06_lncRNA_characterization.R` | ✓ Implemented | Merged GTF | Supplementary 2A | lncRNA structural features (GTF-based) |
| **07** | `07_feelnc_supplement_panels.R` | ✓ Implemented | FEELnc output | Supplementary Figure | FEELnc classification visualization |
| **08** | `08_tf_kegg_integration.R` | ✓ Implemented | GRN edges + KEGG + DE | Figure 5-6 TF-KEGG | TF & pathway integration (commitment vs maintenance) |
| **09** | `09_mirna_mRNA_integration.R` | ✓ Implemented | miRNA/mRNA DE + miRWalk | Figure 6 dotplots | miRNA-mRNA reciprocal regulation + KEGG pathways |
| **10** | `10_extract_lncRNA_for_intarna.R` | ✓ Implemented | lncRNA DE + FASTA | `lncRNA_candidates_*.fa` | Prep IntaRNA inputs (phase-specific candidates) |
| **11** | `11_ceRNA_coherent_pairs.R` | ✓ Implemented | IntaRNA + DE + targets | Figure 7 tables/plots | Coherent triplet networks (lncRNA-miRNA-mRNA) |
| **12** | `12_curated_fig4a_fig6c.R` | ✓ Implemented | TF/miRNA-KEGG selections | Figure 4A, 6C final | Curated shortlist pathway plots |
| **13** | `13_redraw_remaining_panels.R` | ✓ Implemented | All figure outputs | Final refined panels | Panel refinement & consistent styling |
| **14** | `14_results_report.R` | ✓ Implemented | All results directories | `MASTER_results_report.txt` | Master results aggregation & narrative |

---

## File Structure

```
NE_NRVM_pipeline/
├── README.md                             # This file
├── suz_pipeline.bash                     # Top-level orchestrator
├── envs/
│   ├── environment_qc_align.yml          # QC, trim, align, assembly, Salmon, FEELnc
│   ├── environment_mirge.yml             # miRge3.0
│   ├── environment_intarna.yml           # IntaRNA
│   └── environment_r.yml                 # R + Bioconductor
├── metadata/
│   ├── longRNA_samples.tsv               # Long RNA sample sheet
│   ├── smallRNA_samples.tsv              # Small RNA sample sheet
│   ├── fastq_manifest_template.tsv       # 👉 FILL WITH PATHS/CHECKSUMS
│   ├── reference_manifest.tsv            # Reference tracking
│   ├── software_tool_inventory.tsv       # Tool versions
│   └── R_package_inventory_from_uploaded_scripts.tsv
├── references/
│   ├── README.md                         # How references are managed
│   └── (Genome FASTA + GTF downloaded here by stage 00)
├── scripts/
│   ├── bash/
│   │   ├── _activate_env.sh              # Conda activation helper
│   │   ├── 00_download_references.sh     # Auto-download Ensembl refs
│   │   ├── 01_qc_raw.sh through 12_intarna_run.sh
│   │   └── (all stages, auto-activated)
│   ├── R/
│   │   ├── 01_deseq2_mRNA_lncRNA.R       # ✓ Implemented
│   │   ├── 02_mirna_deseq2.R              # ✓ Implemented
│   │   ├── 10_extract_lncRNA_for_intarna.R # ✓ Implemented
│   │   ├── 03-09, 11-14_*.R              # → Scaffolds (TODO)
│   │   └── (all have optparse --help)
│   └── commands/
│       └── make_fastq_manifest_checksums.sh
├── results/                              # (Created by pipeline)
│   ├── 01_fastqc_raw/
│   ├── 02_trimmed_fastq/
│   ├── ...
│   ├── 10_salmon_quant/
│   ├── 11_mirge3_output/
│   ├── 12_intarna_interactions/
│   ├── 01_deseq2_analysis/
│   ├── 02_mirna_deseq2/
│   ├── 10_intarna_candidates/
│   └── 03-14_*/ (R analysis outputs)
├── archive/
│   ├── raw_scripts_original/             # Original NE_NRVM_all_scripts_raw (unchanged)
│   ├── circRNA_old_figure8/              # Figure8_ceRNA_rebuild_with_QC.R / v2 / v3
│   ├── superseded_duplicates/            # Old versions, competing implementations
│   └── ARCHIVE_NOTES.md                  # Why each file was archived
├── fastq_file_sizes.tsv                  # File sizes from original FASTQ files used
├── fastq_md5_checksums.txt               # MD5 checksums from original FASTQ files used
└── .gitignore
```

---

## Figure Mapping: Manuscript Layout

| Final Figure | Content | Source Script(s) | Source Raw Scripts |
|--------------|---------|------------------|--------------------|
| **Figure 1** | (Conceptual / experimental design) | N/A | N/A |
| **Figure 2** | (Sample characterization / methods) | Bash 01–10 (QC, alignment stats) | N/A |
| **Figure 3** | mRNA temporal programs | `03_figure3_mRNA_programs.R` | `05_Figure3_mRNA_programs.R` |
| **Figure 4** | lncRNA discovery + TF-KEGG integration | `04_lncRNA_heatmap_hubs.R`, `08_tf_kegg_integration.R`, `12_curated_fig4a_fig6c.R` | `06_Figure4_refine_heatmap_hubs.R`, `07_Figure4_clean_volcanoes.R`, `script_final_depurado_TF_KEGG.R`, `18_curated_Fig4A_Fig6C_from_shortlist.R` |
| **Figure 5** | (Pathway/kinetics details) | Various R scripts | Multiple |
| **Figure 6** | miRNA-mRNA pathway integration | `09_mirna_mRNA_integration.R`, `12_curated_fig4a_fig6c.R` | `script_miRNA_mRNA_6_24.R`, `18_curated_Fig4A_Fig6C_from_shortlist.R` |
| **Figure 7** | ceRNA coherent networks (lncRNA-miRNA-mRNA) | `11_ceRNA_coherent_pairs.R` | `12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R` |
| **Supp. Fig. 2A** | lncRNA structural characterization | `06_lncRNA_characterization.R` | `15_lncRNA_characterization_from_GTF.R` |
| **Supp. Fig. 2B** | FEELnc lncRNA classification | `07_feelnc_supplement_panels.R` | `14_FEELnc_supplement_panels.R` |
| **Supp. Figs. 3–7** | Additional panels & tables | `13_redraw_remaining_panels.R` | `17_redraw_remaining_panels_Fig4_Fig5_Fig6.R` |

**Note:** Old "Figure 8" (circRNA ceRNA network) content is archived in `archive/circRNA_old_figure8/` and **not included** in the final 7-figure manuscript.

---

## Terminology

**Implemented:** Full working code that executes, reads inputs, processes data, and generates outputs. Scripts may be newly written or adapted from existing sources.

**Ported:** Logic extracted from raw scripts (`raw_scripts_original/`) and adapted for the clean pipeline, including Windows→Linux path fixes, optparse CLI arguments, circRNA branch removal (where applicable), and Linux-first defaults. Functionally complete implementations.

---

## Implementation Status

### ✅ Complete (Fully Implemented)

**Pipeline architecture:**

- [x] Folder structure & conda environments (4 yml files for qc_align, mirge, intarna, r)
- [x] Bash pipeline stages 00–12 (12 scripts + orchestrator suz_pipeline.bash)
- [x] Bash helper script (_activate_env.sh) for non-interactive conda activation
- [x] Metadata files (sample sheets, reference manifests, software inventory)
- [x] Archive: raw_scripts_original, circRNA_old_figure8, superseded_duplicates
- [x] FASTQ manifest template with checksums field structure
- [x] Reference auto-download logic (00_download_references.sh)
- [x] Confirmation of "rn8" naming (project nickname; actual latest is GRCr8)

**R analysis scripts - All Implemented (ported from raw_scripts_original):**

- [x] **01_deseq2_mRNA_lncRNA.R** — Integrated DESeq2 + Mfuzz clustering (from Integrated_Bulk_v3_6_24h.R)
- [x] **02_mirna_deseq2.R** — miRNA bulk DESeq2 analysis (from miRNAs_6_24.R)
- [x] **03_figure3_mRNA_programs.R** — mRNA temporal programs with phase classification (from 05_Figure3_mRNA_programs.R)
- [x] **04_figure4_lncRNA_heatmap_hubs.R** — lncRNA heatmaps and hub visualization (from 06_Figure4_refine_heatmap_hubs.R)
- [x] **05_figure4_volcanoes.R** — lncRNA volcano plots at 6h/24h (from 07_Figure4_clean_volcanoes.R)
- [x] **06_lncRNA_characterization.R** — GTF-based structural characterization (from 15_lncRNA_characterization_from_GTF.R)
- [x] **07_feelnc_supplement_panels.R** — FEELnc classification panels (from 14_FEELnc_supplement_panels.R)
- [x] **08_tf_kegg_integration.R** — TF-GRN + KEGG pathway integration (from script_final_depurado_TF_KEGG.R)
- [x] **09_mirna_mRNA_integration.R** — miRNA-mRNA reciprocal regulation + KEGG mapping (from script_miRNA_mRNA_6_24.R)
- [x] **10_extract_lncRNA_for_intarna.R** — Phase-specific lncRNA FASTA extraction for IntaRNA (from 09_extract_lncRNA_candidates_for_IntaRNA.R)
- [x] **11_ceRNA_coherent_pairs.R** — IntaRNA interaction reading, energy filtering, coherent triplet integration (from 12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R; circRNA removed)
- [x] **12_curated_fig4a_fig6c.R** — Pathway filtering and curated dotplot generation for Fig4A/6C (from 18_curated_Fig4A_Fig6C_from_shortlist.R)
- [x] **13_redraw_remaining_panels.R** — Theme-aware panel refinement with consistent styling (from 17_redraw_remaining_panels_Fig4_Fig5_Fig6.R)
- [x] **14_results_report.R** — Master report aggregation with pipeline tracking and validation (from MASTER_results_report_builder_v3_8_1.R; fig8/circRNA removed)

All 14 R scripts include:

- ✓ optparse CLI argument parsing with `--help` documentation
- ✓ Windows path replacement with Linux-friendly defaults
- ✓ Proper error handling and informative logging
- ✓ Automatic output directory creation
- ✓ CSV/PDF output generation with configurable DPI

### ⏳ Pending Activities

- [ ] **Update R scripts with current figure revisions:** Sebastián and JP are currently refining all paper figures. Once updates are finalized:
  - [ ] Receive updated R script logic from raw_scripts_original or new versions
  - [ ] Integrate updates into cleaned R scripts (03-14) within existing optparse framework
  - [ ] Preserve CLI arguments and parameter structure
  - [ ] Validate updated scripts pass syntax checks
- [ ] **FASTQ manifest completion:** Exact server paths, file sizes, md5 checksums for:
  - 32 long RNA (lncRNA discovery) FASTQ files
  - 16 small RNA (miRNA) FASTQ files
- [ ] **TFLink All reference:** Path/source + checksum (large file; coordinate upload separately)
- [ ] Optional: Backup hard drive paths for FASTQ redundancy (for fastq_manifest.tsv)

- [] Include characterization scripts and notebooks that have generated metadata files used in Rscript 01

**Testing & Validation:**

- [] Prepare conda environments based on template `*.yaml` files and then update with the official exported environment `*.yaml` file with consolidated dependencies and package versions.
- [ ] **Full pipeline test run:** Execute complete pipeline (bash 00-12 + R 01-14) on real FASTQ data:
  - [ ] Verify all outputs generate without errors
  - [ ] Check figure quality and content against manuscript drafts
  - [ ] Validate DESeq2 results match expectations
  - [ ] Confirm IntaRNA predictions are reasonable
  - [ ] Test with `--run-r-analysis` automatic mode
  - [ ] Generate and review MASTER_results_report.txt
---

## References

### Final References (Used by Default)

- **Genome FASTA:** `Rattus_norvegicus.GRCr8.dna.toplevel.fa` (Ensembl GRCr8, toplevel assembly)
  - Download: `https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz`
- **GTF annotation:** `Rattus_norvegicus.GRCr8.115.gtf` (Ensembl release 115)
  - Download: `https://ftp.ensembl.org/pub/release-115/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gtf.gz`

These are **auto-downloaded** by `00_download_references.sh` into `references/`.

---

## Usage Examples

### Run Full Pipeline (Default References)

```bash
./suz_pipeline.bash \
  --longrna-samples metadata/longRNA_samples.tsv \
  --smallrna-samples metadata/smallRNA_samples.tsv \
  --fastq-manifest metadata/fastq_manifest_template.tsv \
  --threads 16 \
  --out results
```

### Run Specific Stages Only

```bash
# Just QC and trimming
./suz_pipeline.bash --stages "1,2,3" --threads 16

# Just alignment through Salmon (skip miRge/IntaRNA)
./suz_pipeline.bash --stages "4,5,6,7,8,9,10" --threads 16
```

### Use Custom References

```bash
./suz_pipeline.bash \
  --genome-fasta /data/custom_genome.fa \
  --gtf /data/custom_annotation.gtf \
  --skip-download-refs \
  ...
```

### Run Individual R Scripts (After Bash Complete)

```bash
# DESeq2 analysis
Rscript scripts/R/01_deseq2_mRNA_lncRNA.R \
  --gtf results/06_stringtie_assembly/merged.annotated.gtf \
  --salmon-dir results/10_salmon_quant \
  --sample-sheet metadata/longRNA_samples.tsv \
  --output-dir results/01_deseq2_analysis

# miRNA analysis
Rscript scripts/R/02_mirna_deseq2.R \
  --counts-dir results/11_mirge3_output \
  --sample-sheet metadata/smallRNA_samples.tsv \
  --output-dir results/02_mirna_deseq2

# Extract lncRNA candidates for IntaRNA
Rscript scripts/R/10_extract_lncRNA_for_intarna.R \
  --de-results-dir results/01_deseq2_analysis \
  --fasta-file results/08_gffread_fasta/merged_transcripts.fasta \
  --output-dir results/10_intarna_candidates
```

---

## Troubleshooting

### Conda Environment Creation Fails

- Check that `conda`/`mamba` is installed and in `$PATH`
- Try `mamba env create` (faster solver than conda)
- Run with `--verbose` flag for detailed error messages

### HISAT2 Index Build Fails

- Ensure genome FASTA is uncompressed and in the correct format
- Check disk space (index may require ~10–20 GB for rat genome)

### Salmon Quant Crashes

- Verify transcriptome FASTA is not corrupted (check line count, file size)
- Ensure trimmed FASTQ files are gzip-compressed (`.fastq.gz`)

### miRge3.0 Fails or Outputs "No miRNA Found"

- Verify small RNA library paths and format
- Ensure metadata CSV (if using `-dex`) matches sample names in manifest
- Check that small RNA FASTQ files are properly gzipped

### IntaRNA Produces Empty Results

- Check that lncRNA FASTA candidates were generated (run R script 10 first)
- Verify miRNA FASTA files exist in miRge3.0 output directory
- IntaRNA may need library path setup; see `12_intarna_run.sh` help

### R Scripts Fail on Missing Files

- Ensure upstream bash stages completed successfully
- Check that output directory paths match expected defaults in `scripts/bash/suz_pipeline.bash`
- Run individual R scripts with `--help` to verify required arguments

---

## Updating R Scripts with Figure Revisions

When Sebastián and JP finalize figure updates, integrate them into the pipeline as follows:

### 1. Receive Updated Scripts

Sebastián will send updated R script versions reflecting final figure designs. These should be added to `raw_scripts_original/` for archival, then ported into the clean pipeline.

### 2. Port Updates to Cleaned Scripts

For each updated script (e.g., `03_figure3_mRNA_programs.R`):

```bash
# 1. Archive previous version (optional)
git mv scripts/R/03_figure3_mRNA_programs.R archive/superseded_duplicates/03_figure3_mRNA_programs_v1.R

# 2. Extract logic from Sebastián's updated script
# 3. Adapt to optparse framework (preserve existing CLI arguments)
# 4. Test syntax
Rscript scripts/R/03_figure3_mRNA_programs.R --help

# 5. Commit
git add scripts/R/03_figure3_mRNA_programs.R archive/superseded_duplicates/
git commit -m "Update: Fig3 mRNA programs with JP revisions (plot styling, colors)"
```

### 3. Preserve Optparse Interface

Keep existing CLI arguments intact:

- `--de-dir` (input directory)
- `--output-dir` (output directory)
- `--dpi` (figure resolution)
- `--top-per-phase` or script-specific tuning parameters

Do not rename or remove existing parameters — the bash orchestrator (`_run_r_analysis.sh`) references them.

### 4. Test Updates

After updating each script:

```bash
# Syntax check
Rscript scripts/R/03_figure3_mRNA_programs.R --help

# Manual test (if possible with sample data)
Rscript scripts/R/03_figure3_mRNA_programs.R \
  --de-dir results/01_deseq2_analysis \
  --output-dir results/03_figure3_test \
  --dpi 600
```

---

## Next Steps for Running the Pipeline

### Before First Run

1. **Complete FASTQ manifests:**
   - Obtain from Sebastián: 32 long RNA + 16 small RNA file paths, sizes, and md5 checksums
   - Populate `metadata/fastq_manifest.tsv` (server paths or hard drive backup paths)
   - Use `commands/make_fastq_manifest_checksums.sh` to auto-generate checksums if files are locally accessible

2. **Prepare TFLink reference (optional for Figure 5-6):**
   - Request from Sebastián: TFLink All rat transcription factor network
   - Place in `references/` and note path in reference_manifest.tsv
   - If unavailable, TF panels (Figure 5-6) can be generated from alternative sources

3. **Verify conda installation:**
   - Ensure `conda` or `mamba` is in `$PATH`
   - Run `conda --version` to confirm

### First Run (Full Pipeline)

```bash
cd NE_NRVM_pipeline

# 1. Create all conda environments
conda env create -f envs/environment_qc_align.yml
conda env create -f envs/environment_mirge.yml
conda env create -f envs/environment_intarna.yml
conda env create -f envs/environment_r.yml

# 2. Run bash pipeline stages (auto-activates environments)
./suz_pipeline.bash \
  --longrna-samples metadata/longRNA_samples.tsv \
  --smallrna-samples metadata/smallRNA_samples.tsv \
  --fastq-manifest metadata/fastq_manifest.tsv \
  --threads 16 \
  --out results

# Expected runtime: 24–72 hours depending on FASTQ sizes and hardware
# Monitor output for errors in results/logs/
```

### After Bash Stages Complete

```bash
# 3. Run R analysis scripts in order
cd scripts/R

# Core DESeq2 (required by all figure scripts)
Rscript 01_deseq2_mRNA_lncRNA.R \
  --gtf ../results/06_stringtie_assembly/merged.annotated.gtf \
  --salmon-dir ../results/10_salmon_quant \
  --sample-sheet ../metadata/longRNA_samples.tsv \
  --output-dir ../results/01_deseq2_analysis

Rscript 02_mirna_deseq2.R \
  --counts-dir ../results/11_mirge3_output \
  --sample-sheet ../metadata/smallRNA_samples.tsv \
  --output-dir ../results/02_mirna_deseq2

Rscript 10_extract_lncRNA_for_intarna.R \
  --de-dir ../results/01_deseq2_analysis \
  --fasta-file ../results/08_gffread_fasta/merged_transcripts.fasta \
  --output-dir ../results/10_intarna_candidates

# Figure generation (03–09)
Rscript 03_figure3_mRNA_programs.R \
  --de-dir ../results/01_deseq2_analysis \
  --output-dir ../results/03_figure3 \
  --dpi 600

# [Continue with 04–09...]

# 4. Generate final results report
Rscript 14_results_report.R \
  --results-root ../results \
  --output-file ../results/MASTER_results_report.txt
```

### Validation & QC

- Check for empty output directories (indicates missing upstream data)
- Review QC reports in `results/01_fastqc_raw/` and `results/03_fastqc_trimmed/`
- Compare alignment statistics to original manuscript (if available)
- Spot-check DE gene lists for known marker genes

### For Publication

1. Verify all 7 figures are generated (Figure 1 is typically conceptual/manual)
2. Use `14_results_report.R` to generate reproducibility audit trail
3. Archive all conda environment snapshots (via `conda env export`)
4. Create data availability statement referencing:
   - Pipeline repository URL/DOI
   - FASTQ accession numbers (SRA / GEO)
   - Reference genome version (GRCr8, Ensembl 115)

---

## What Changed from Raw Scripts → Clean Pipeline

### Key Improvements

| Issue | Raw Scripts | Clean Pipeline |
| --- | --- | --- |
| **Hardcoded paths** | `C:/Users/sebau/...`, UCSC rn7, GRCr8.114 | Parametrized via CLI, defaults to GRCr8 + Ensembl 115 |
| **Reference versions** | Mixed versions (rn7, mRatBN7.2.112, custom) | Single canonical pair: GRCr8.dna.toplevel + GRCr8.115.gtf |
| **Conda isolation** | Monolithic or implicit conda usage | 4 explicit environments per tool group (auto-activated) |
| **Script organization** | 24 R scripts, duplicates, competing versions | 14 canonical R scripts (8 fully ported, 6 scaffolded) |
| **CircRNA content** | 8 scripts for old Figure 8 (circRNA) | Archived in `archive/circRNA_old_figure8/` (not executed) |
| **CLI arguments** | Hardcoded file paths, manual edits required | optparse-driven, `--help` documentation for all scripts |
| **Windows→Linux** | Manual path fixes needed per run | Automatic path handling, Linux-first defaults |
| **Error handling** | Minimal (silent failures) | Explicit logging, early exit on missing files |
| **Reproducibility** | Implicit dependencies on local setup | Conda envs + versioned tool lists + reference manifests |

### Superseded Scripts (Archived)

See `archive/ARCHIVE_NOTES.md` for complete list. Examples:

- `miRNA_DEG_ceRNAs.R` — undefined variable bug; superseded by `02_mirna_deseq2.R`
- `08_Figure7_ceRNA_IntaRNA_scaffold.R` — prototype; superseded by `11_ceRNA_coherent_pairs.R`
- `13_generate_paper_results_report_v2.R` — superseded by `14_results_report.R` (v3)
- `Figure8_ceRNA_rebuild_with_QC*.R` (v1/v2/v3) — out of scope for 7-figure final manuscript

---

## Contact / Questions

For pipeline setup or technical issues, refer to:

- `ARCHIVE_NOTES.md` — why scripts were superseded or archived
- `scripts/bash/*.sh --help` — each stage has detailed CLI help
- `scripts/R/*.R --help` — R analysis scripts have optparse help
- `metadata/reference_manifest.tsv` — reference sources and checksums

---