# Archive Notes: Superseded and Out-of-Scope Scripts

This document explains why certain R scripts were archived and what superseded them in the final pipeline.

## Out-of-Scope: circRNA and Old "Figure 8" Content

The following scripts are archived because they implement analysis for a previous version of the manuscript that included circRNA predictions (Figure 8). **The current manuscript has 7 figures, not 8**, and circRNA content is explicitly out of scope per the coauthor's instructions.

- **`R_scripts/Figure8_ceRNA_rebuild_with_QC.R`** (v1) — Old circRNA+ceRNA integration; archived.
- **`R_scripts/Figure8_ceRNA_rebuild_with_QC_v2.R`** (v2) — Refined version of above; archived.
- **`R_scripts/Figure8_ceRNA_rebuild_with_QC_v3.R`** (v3) — Latest version, still circRNA; archived.

All three are preserved in `archive/circRNA_old_figure8/` for historical reference only.

Note: `R_scripts/12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R` (in the final pipeline) originally contained mixed lncRNA + circRNA branches. The lncRNA-only logic was extracted for the final `scripts/R/11_ceRNA_coherent_pairs.R`, and circRNA branches were removed.

---

## Superseded Duplicates: Competing Analysis Implementations

The following scripts are archived because they were superseded by more complete or newer versions of the same analysis step.

### Core Bulk DESeq2 + Mfuzz Pipeline

- **`R_scripts/Comparación_Cinética_6_24.R`** — Early DESeq2 + kinetics analysis (2 script versions pasted in one file, mixed Spanish/English comments, Windows hardcoded paths). Superseded by `Integrated_Bulk_v3_6_24h.R`.
- **`R_scripts/Versión_Nueva_Nov.R`** — Updated DESeq2/enrichment pipeline, still early draft status. Superseded by `Integrated_Bulk_v3_6_24h.R`.
- **`Integrated_Bulk_v3_6_24h.R`** is the canonical source → ported to final pipeline as `scripts/R/01_deseq2_mRNA_lncRNA.R` (internally deduplicated, Windows paths fixed, optparse-ified).

All three are preserved in `archive/superseded_duplicates/` for reference.

### miRNA DESeq2 Pipeline

- **`R_scripts/miRNA_DEG_ceRNAs.R`** — Contains miRNA + ceRNA integration, but has a confirmed bug (undefined `vsd_mi` variable) and partially duplicates the next script. Superseded by `miRNAs_6_24.R`.
- **`R_scripts/miRNAs_6_24.R`** is the canonical source → ported to final pipeline as `scripts/R/02_mirna_deseq2.R` (optparse-ified).

Archived in `archive/superseded_duplicates/`.

### ceRNA/lncRNA-miRNA-mRNA Integration

- **`R_scripts/08_Figure7_ceRNA_IntaRNA_scaffold.R`** — Simpler, single-timepoint scaffold prototype for ceRNA triplets. Superseded by `12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R`.
- **`R_scripts/12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R`** is the canonical, multi-timepoint version → ported to final pipeline as `scripts/R/11_ceRNA_coherent_pairs.R` (circRNA branches removed).

Archived in `archive/superseded_duplicates/`.

### lncRNA Characterization (Supplementary Figure)

- **`R_scripts/15_lncRNA_characterization_from_GTF.R`** — GTF-based characterization (reproducible, fully derived from pipeline).
- **`R_scripts/16_lncRNA_characterization_curatedCSV_only.R`** — CSV-based characterization (requires manual curation, different source).

Both are marked `useful_final_candidate` in the script inventory with no explicit preference. **The final pipeline uses `scripts/R/06_lncRNA_characterization.R` (ported from script 15, GTF-based)** because it's fully reproducible. Script 16 is archived in `archive/superseded_duplicates/` with a note suggesting the coauthor should confirm this choice.

### Results Report Generator

- **`R_scripts/13_generate_paper_results_report_v2.R`** — Simple, hardcoded-paths report builder. Superseded by `MASTER_results_report_builder_v3_8_1.R`.
- **`R_scripts/MASTER_results_report_builder_v3_8_1.R`** is the canonical version (CLI-driven, configurable) → ported to final pipeline as `scripts/R/14_results_report.R` (fig8/circRNA slots removed).

Archived in `archive/superseded_duplicates/`.

---

## Bash Command History

- **`bash_scripts/bash_scripts.txt`** — A Spanish-language working notebook with command history, notes, tests, and dead-end attempts. This was manually analyzed to extract the canonical pipeline order and rewritten into clean, modular bash scripts (`scripts/bash/01`–`12`). The original file is archived verbatim in `archive/raw_scripts_original/bash_scripts/` for provenance.

---

## Reference

When investigating a choice or confirming a decision, see:
- `scripts_inventory/script_inventory_preliminary.tsv` — official triage of all scripts (status: final/archive/extract-logic, etc.)
- `NE_NRVM_coauthor_reproducibility_starter_v1/README.md` — coauthor's notes on scope and decisions
