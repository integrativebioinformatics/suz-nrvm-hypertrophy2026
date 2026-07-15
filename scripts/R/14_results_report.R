#!/usr/bin/env Rscript
# Stage 14: Master results report — data-driven aggregator
# Scans the pipeline output tree, reads the actual result tables, and reports
# real counts and statistics per figure/analysis (DE genes by biotype and
# direction, phase-set sizes, miRNA DE, lncRNA characterization inventory,
# FEELnc classes, KEGG pathways, TF/miRNA-KEGG edges, ceRNA coherent pairs/
# triplets). circRNA / 8-figure content is intentionally excluded.
#
# Ported from MASTER_results_report_builder_v3_8_1.R with the config-driven
# result extraction preserved (fig8/circRNA slots removed) and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

option_list <- list(
  make_option(c("-r", "--results-root"), type = "character", default = "results",
              help = "Root results directory [default: %default]"),
  make_option(c("-o", "--output-file"), type = "character", default = "results/MASTER_results_report.txt",
              help = "Output report file [default: %default]"),
  make_option(c("--summary-csv"), type = "character", default = NULL,
              help = "Optional machine-readable metrics CSV [default: <report_dir>/MASTER_results_metrics.csv]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
RESULTS_ROOT <- opt$`results-root`
REPORT_FILE <- opt$`output-file`
if (!dir.exists(RESULTS_ROOT)) { cat("[ERROR] Results root does not exist:", RESULTS_ROOT, "\n"); quit(status = 1) }
out_dir <- dirname(REPORT_FILE); if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
SUMMARY_CSV <- if (is.null(opt$`summary-csv`)) file.path(out_dir, "MASTER_results_metrics.csv") else opt$`summary-csv`

cat("[INFO] Results root:", RESULTS_ROOT, "\n[INFO] Report output:", REPORT_FILE, "\n\n")

# ---------------------------------------------------------------------------
# Helpers: locate + read result tables, extract real numbers
# ---------------------------------------------------------------------------

metrics <- tibble(metric = character(), value = character(), source_file = character())
add_metric <- function(metric, value, src = "") {
  metrics <<- bind_rows(metrics, tibble(metric = metric, value = as.character(value), source_file = src))
}

find_one <- function(fname) {
  hits <- list.files(RESULTS_ROOT, pattern = paste0("^", str_replace_all(fname, "([.\\^$|()\\[\\]{}*+?\\\\])", "\\\\\\1"), "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) NA_character_ else hits[1]
}
read_tbl <- function(fname) {
  p <- find_one(fname); if (is.na(p)) return(NULL)
  suppressWarnings(tryCatch(
    if (str_detect(fname, "\\.txt$")) read_tsv(p, show_col_types = FALSE) else read_csv(p, show_col_types = FALSE),
    error = function(e) NULL))
}
count_dir <- function(df, col) {
  if (is.null(df) || !col %in% names(df)) return(NULL)
  as.data.frame(table(str_to_upper(as.character(df[[col]]))), stringsAsFactors = FALSE)
}

report <- c(
  "===============================================================================",
  "NE_NRVM RNA-seq Analysis — Master Results Report (data-driven)",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Results root:", normalizePath(RESULTS_ROOT)),
  "===============================================================================", "")

section <- function(title) report <<- c(report, "", title, strrep("-", nchar(title)))

# ---------------------------------------------------------------------------
# Figure 3 / Stage 01: mRNA & lncRNA DE + temporal programs
# ---------------------------------------------------------------------------

section("Figure 3 — mRNA temporal programs (stage 01)")
for (bt in c("mRNA", "lncRNA")) for (tp in c("6h", "24h")) {
  df <- read_tbl(paste0("NE_vs_Ctrl_", tp, "_", bt, "_DE.txt"))
  if (is.null(df)) { report <- c(report, sprintf("  [ ] %s %s DE table not found", bt, tp)); next }
  dcount <- count_dir(df, "direction")
  up <- if (!is.null(dcount)) sum(dcount$Freq[dcount$Var1 == "UP"]) else NA
  dn <- if (!is.null(dcount)) sum(dcount$Freq[dcount$Var1 == "DOWN"]) else NA
  report <- c(report, sprintf("  [x] %s %s: %d rows (UP=%s, DOWN=%s)", bt, tp, nrow(df), up, dn))
  add_metric(paste0("DE_", bt, "_", tp, "_total"), nrow(df))
  add_metric(paste0("DE_", bt, "_", tp, "_UP"), up); add_metric(paste0("DE_", bt, "_", tp, "_DOWN"), dn)
}
for (ph in c("early", "sustained", "late")) for (bt in c("mRNA", "lncRNA", "novel_lncRNA")) {
  df <- read_tbl(paste0(ph, "_", bt, ".csv"))
  if (!is.null(df)) { report <- c(report, sprintf("  [x] phase %s / %s: %d genes", ph, bt, nrow(df)))
    add_metric(paste0("phase_", ph, "_", bt), nrow(df)) }
}
lrt <- read_tbl("kinetic_genes_lrt.csv")
if (!is.null(lrt)) { report <- c(report, sprintf("  [x] kinetic genes (LRT): %d", nrow(lrt))); add_metric("kinetic_genes_LRT", nrow(lrt)) }
mf <- read_tbl("mfuzz_cluster_assignment.csv")
if (!is.null(mf) && "cluster" %in% names(mf)) {
  report <- c(report, sprintf("  [x] Mfuzz: %d genes in %d clusters", nrow(mf), length(unique(mf$cluster))))
  add_metric("mfuzz_genes", nrow(mf)); add_metric("mfuzz_clusters", length(unique(mf$cluster)))
}

# ---------------------------------------------------------------------------
# Figure 6 / Stage 02: miRNA DE
# ---------------------------------------------------------------------------

section("Figure 6 — miRNA differential expression (stage 02)")
for (tp in c("6h", "24h")) {
  up <- read_tbl(paste0("miRNA_DE_", tp, "_UP.txt")); dn <- read_tbl(paste0("miRNA_DE_", tp, "_DOWN.txt"))
  nu <- if (is.null(up)) NA else nrow(up); nd <- if (is.null(dn)) NA else nrow(dn)
  report <- c(report, sprintf("  [%s] miRNA %s: UP=%s, DOWN=%s", ifelse(is.null(up) && is.null(dn), " ", "x"), tp, nu, nd))
  add_metric(paste0("miRNA_DE_", tp, "_UP"), nu); add_metric(paste0("miRNA_DE_", tp, "_DOWN"), nd)
}
tt <- read_tbl("miRNA_temporal_summary_table.csv")
if (!is.null(tt) && "temporal_pattern" %in% names(tt)) {
  de <- tt %>% filter(temporal_pattern != "Non_DE")
  report <- c(report, sprintf("  [x] miRNA temporal patterns: %d DE miRNAs across %d patterns", nrow(de), length(unique(de$temporal_pattern))))
}

# ---------------------------------------------------------------------------
# Supp / Stages 06-07: lncRNA characterization + FEELnc
# ---------------------------------------------------------------------------

section("Supplementary — lncRNA characterization (stages 06-07)")
inv <- read_tbl("lncRNA_inventory_summary.tsv")
if (!is.null(inv) && all(c("metric", "source", "value") %in% names(inv))) {
  tg <- inv %>% filter(metric == "total_genes")
  for (i in seq_len(nrow(tg))) { report <- c(report, sprintf("  [x] %s lncRNA genes: %s", tg$source[i], tg$value[i]))
    add_metric(paste0("lncRNA_genes_", tg$source[i]), tg$value[i]) }
}
fe <- read_tbl("SuppFig2B_FEELnc_class_counts.tsv")
if (!is.null(fe)) { report <- c(report, sprintf("  [x] FEELnc classified lncRNA-gene rows: %d", nrow(fe))); add_metric("FEELnc_class_rows", nrow(fe)) }

# ---------------------------------------------------------------------------
# Figures 4/5 / Stage 08: KEGG + TF integration
# ---------------------------------------------------------------------------

section("Figures 4/5 — KEGG dynamics + TF integration (stage 08)")
kegg <- read_tbl("SUPP_Fig5_full_KEGG_stats_wide.csv")
if (!is.null(kegg) && "KEGG_group" %in% names(kegg)) {
  for (g in c("Only6h", "Shared", "Only24h")) { n <- sum(kegg$KEGG_group == g)
    report <- c(report, sprintf("  [x] KEGG pathways (%s): %d", g, n)); add_metric(paste0("KEGG_", g), n) }
}
tfc <- read_tbl("SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv")
if (!is.null(tfc)) {
  ntf <- length(unique(if ("TF_key" %in% names(tfc)) tfc$TF_key else tfc$TF))
  report <- c(report, sprintf("  [x] TF-KEGG edges: %d rows across %d TFs", nrow(tfc), ntf))
  add_metric("TF_KEGG_edge_rows", nrow(tfc)); add_metric("TF_KEGG_unique_TFs", ntf)
}

# ---------------------------------------------------------------------------
# Figure 6 / Stage 09: miRNA-mRNA integration
# ---------------------------------------------------------------------------

section("Figure 6 — miRNA-mRNA-KEGG integration (stage 09)")
allc <- read_tbl("SUPP_edges_ALL_COMBOS_KEGG_phase.csv")
if (!is.null(allc)) { report <- c(report, sprintf("  [x] miRNA-mRNA-KEGG edges (all combos): %d", nrow(allc))); add_metric("miRNA_mRNA_edges", nrow(allc)) }

# ---------------------------------------------------------------------------
# Figure 7 / Stages 10-11: IntaRNA candidates + ceRNA
# ---------------------------------------------------------------------------

section("Figure 7 — ceRNA coherent networks (stages 10-11)")
cand <- read_tbl("SUMMARY_lncRNA_candidates_by_time_direction_phase_source.csv")
if (!is.null(cand) && "n_lncRNAs" %in% names(cand)) {
  report <- c(report, sprintf("  [x] lncRNA IntaRNA candidates (sum): %d", sum(cand$n_lncRNAs)))
  add_metric("lncRNA_intarna_candidates", sum(cand$n_lncRNAs))
}
pairs <- read_tbl("COHERENT_lncRNA_miRNA_pairs.csv")
if (!is.null(pairs)) { report <- c(report, sprintf("  [x] coherent lncRNA-miRNA pairs: %d", nrow(pairs))); add_metric("coherent_pairs", nrow(pairs)) }
trip <- read_tbl("COHERENT_lncRNA_miRNA_mRNA_triplets.csv")
if (!is.null(trip)) { report <- c(report, sprintf("  [x] coherent lncRNA-miRNA-mRNA triplets: %d", nrow(trip))); add_metric("coherent_triplets", nrow(trip)) }

# ---------------------------------------------------------------------------
# Stage inventory (file counts)
# ---------------------------------------------------------------------------

section("Pipeline output inventory")
stage_dirs <- list.dirs(RESULTS_ROOT, recursive = FALSE, full.names = TRUE)
for (d in sort(stage_dirs)) {
  n <- length(list.files(d, recursive = TRUE))
  report <- c(report, sprintf("  %-40s %d files", basename(d), n))
}

report <- c(report, "", "",
  "DATA AVAILABILITY", "-----------------",
  "Raw sequencing data: [TO_FILL: SRA/GEO accession]",
  "Reference genome: Rattus_norvegicus.GRCr8.dna.toplevel (Ensembl 115)",
  "Reference annotation: Rattus_norvegicus.GRCr8.115.gtf (Ensembl 115)",
  "Conda environments: see envs/",
  "", "===============================================================================",
  "End of report")

writeLines(report, REPORT_FILE)
if (nrow(metrics) > 0) write_csv(metrics, SUMMARY_CSV)

cat("\n[SUCCESS] Master results report generated (", nrow(metrics), "metrics extracted).\n", sep = "")
cat("[OUTPUT]", REPORT_FILE, "\n")
if (nrow(metrics) > 0) cat("[OUTPUT]", SUMMARY_CSV, "\n")
