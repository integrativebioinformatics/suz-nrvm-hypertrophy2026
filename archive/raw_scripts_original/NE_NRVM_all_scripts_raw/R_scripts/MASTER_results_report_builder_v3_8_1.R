#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

required_pkgs <- c("dplyr", "purrr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Please install them first, e.g. install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))."
    )
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
PROJECT_ROOT <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else getwd()
REPORT_DIR <- if (length(args) >= 2) file.path(PROJECT_ROOT, args[[2]]) else file.path(PROJECT_ROOT, "MASTER_results_report")

TXT_OUT <- file.path(REPORT_DIR, "MASTER_results_report.txt")
MD_OUT <- file.path(REPORT_DIR, "MASTER_results_report.md")
FIG_REPORT_DIR <- file.path(REPORT_DIR, "per_figure_reports")
TABLE_OUT_DIR <- file.path(REPORT_DIR, "candidate_tables")
LOG_FILE <- file.path(REPORT_DIR, "report_log.txt")
INPUT_AUDIT <- file.path(REPORT_DIR, "selected_input_files.tsv")

for (x in c(REPORT_DIR, FIG_REPORT_DIR, TABLE_OUT_DIR)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

writeLines(character(), LOG_FILE)

log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# =========================
# USER-EDITABLE OVERRIDES
# =========================
# Use relative paths from PROJECT_ROOT if auto-discovery picks the wrong files.
FILE_OVERRIDES <- list(
  fig2_mrna = c(
    "DESeq2_Multifactorial_Results_fixed/NE_vs_Ctrl_6h_mRNA_DE.txt",
    "DESeq2_Multifactorial_Results_fixed/NE_vs_Ctrl_24h_mRNA_DE.txt"
  ),
  fig3_mrna_phase = character(0),
  fig4_lnc_de = c(
    "DESeq2_Multifactorial_Results_fixed/NE_vs_Ctrl_6h_lncRNA_DE.txt",
    "DESeq2_Multifactorial_Results_fixed/NE_vs_Ctrl_24h_lncRNA_DE.txt",
    "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT/SUPP_Fig4B2_lncDE_with_phase.csv",
    "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT/SUPP_Fig4B2_lncDE_counts_phase_time_direction.csv",
    "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT/SUPP_Fig4B1_lnc_counts_phase_source.csv",
    "Paper_Fig4_REFINED/SUPP_Fig4D_top_lnc_for_heatmap.csv"
  ),
  fig4_lnc_network = character(0),
  fig5_kegg = c(
    "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT/SUPP_Fig5_full_KEGG_stats_wide.csv",
    "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT/SUPP_Fig5_selected_pathways_dotplot_data.csv"
  ),
  fig6_tf = character(0),
  fig7_mirna_mrna = character(0),
  fig8_cerna = character(0),
  supp_circ = character(0)
)

# Maximum number of auto-selected files per slot.
MAX_FILES_PER_SLOT <- 12

# =========================
# SLOT CONFIGURATION
# =========================
SLOTS <- list(
  fig2_mrna = list(
    title = "Figure 2 - mRNA differential expression at 6 h and 24 h",
    dirs = c("DESeq2_Multifactorial_Results_fixed", "Paper_Fig3_FINAL"),
    include = c("deseq2", "multifactor", "de", "mrna", "protein", "coding", "sig", "significant", "6h", "24h"),
    exclude = c("kegg", "gsea", "tf", "mirna", "lnc", "lncrna", "cerna", "circ", "correlation", "correlations", "rewir", "rewiring", "hub", "network", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("time", "gene_symbol", "gene_id", "log2fc", "padj")
  ),
  fig3_mrna_phase = list(
    title = "Figure 3 - Early / Sustained / Late / Reversal mRNA classes",
    dirs = c("Paper_Fig3_FINAL", "DESeq2_Multifactorial_Results_fixed"),
    include = c("early", "sustained", "late", "reversal", "class", "phase", "kinetic", "shared", "only6", "only24", "mrna", "gene"),
    exclude = c("kegg", "gsea", "tf", "mirna", "lnc", "lncrna", "cerna", "circ", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("class", "time", "gene_symbol", "gene_id")
  ),
  fig4_lnc_de = list(
    title = "Figure 4 - lncRNA differential expression and response classes",
    dirs = c("Paper_Fig4_REFINED", "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT"),
    include = c("lnc", "lncrna", "noncoding", "novel", "known", "heatmap", "de", "sig", "significant", "class", "phase"),
    exclude = c("mirna", "cerna", "kegg", "gsea", "tf", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("time", "class", "gene_id", "log2fc", "padj", "source")
  ),
  fig4_lnc_network = list(
    title = "Figure 4 - lncRNA hubs and rewiring",
    dirs = c("Paper_Fig4_REFINED", "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT"),
    include = c("hub", "rewir", "network", "degree", "centrality", "lnc"),
    exclude = c("mirna", "cerna", "kegg", "gsea", "tf", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("gene_id", "time", "score", "source")
  ),
  fig5_kegg = list(
    title = "Figure 5 - KEGG programs across commitment and maintenance",
    dirs = c("Paper_Fig4_Fig5_FINAL_ONE_SCRIPT", "Paper_Fig3_FINAL"),
    include = c("kegg", "gsea", "pathway", "enrich", "only6", "only24", "shared", "6h", "24h"),
    exclude = c("tf", "mirna", "cerna", "lnc", "lncrna", "circ", "hub", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("pathway", "class", "time", "direction", "nes", "padj", "count")
  ),
  fig6_tf = list(
    title = "Figure 6 - TF layer controlling commitment vs maintenance",
    dirs = c("Fig6_TF_panels_out_FINAL"),
    include = c("tf", "transcription", "regulator", "motif", "network", "panel"),
    exclude = c("mirna", "cerna", "kegg", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("tf", "class", "time", "padj", "score", "count")
  ),
  fig7_mirna_mrna = list(
    title = "Figure 7 - miRNA-mRNA layer",
    dirs = c("Paper_miRNA_mRNA_AllCombos_Dotplots"),
    include = c("mirna", "mrna", "pair", "interaction", "combo", "dotplot", "candidate"),
    exclude = c("cerna", "lnc", "circ", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("mirna", "gene_symbol", "gene_id", "time", "class", "direction")
  ),
  fig8_cerna = list(
    title = "Figure 8 - integrative ceRNA candidates",
    dirs = c("Paper_ceRNA_candidates"),
    include = c("cerna", "triplet", "triple", "priority", "candidate", "coherent", "lnc", "circ", "mirna", "integrat"),
    exclude = c("pair", "pairs", "debug", "annotation_map", "raw", "all", "full", "complete", "edge", "edges", "png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("mirna", "gene_symbol", "gene_id", "circ_id", "class", "time", "energy")
  ),
  supp_circ = list(
    title = "Supplementary circRNA outputs",
    dirs = c("Paper_ceRNA_candidates", "circRNA", "circ", "Supplementary"),
    include = c("circ", "circrna"),
    exclude = c("png", "pdf", "svg", "jpg", "jpeg"),
    wanted_families = c("circ_id", "time", "class", "log2fc", "padj")
  )
)

# =========================
# COLUMN ALIASES
# =========================
ALIASES <- list(
  time = c("^timepoint$", "^time$", "^tp$", "^phase_time$", "^comparison$", "^contrast$", "sample_group", "condition"),
  class = c("^class$", "^response_class$", "^phase$", "^kinetic_class$", "^group$", "^category$", "^logic_class$", "^module$"),
  direction = c("^direction$", "^class_direction$", "^regulation$", "^trend$", "^status$", "^sign$", "^activation_state$"),
  gene_symbol = c("^gene_symbol$", "^gene_symbol_x$", "^gene_symbol_y$", "^symbol$", "^symbol_x$", "^symbol_y$", "^gene_name$", "^gene_name_x$", "^gene_name_y$", "^external_gene_name$", "^gene$", "^mgi_symbol$"),
  gene_id = c("^gene_id$", "^gene_id_x$", "^gene_id_y$", "^feature_id$", "^feature_id_x$", "^feature_id_y$", "^transcript_id$", "^ensembl_gene_id$", "^ensgene$", "^id$"),
  feature_id = c("^feature_id$", "^id$", "^identifier$", "^rowname$"),
  source = c("^source$", "^source_x$", "^source_y$", "^origin$", "^lncrna_source$", "^annotation_source$", "^source_class$"),
  biotype = c("^biotype$", "^gene_biotype$", "^transcript_biotype$", "^type$", "^gene_type$"),
  padj = c("^padj$", "^padj_x$", "^padj_y$", "^adj_p$", "^adj_p_val$", "^fdr$", "^fdr_x$", "^fdr_y$", "^qvalue$", "^q_value$", "^p_adj$", "^fdr_bh$", "^padjust$"),
  pvalue = c("^pvalue$", "^pvalue_x$", "^pvalue_y$", "^p_value$", "^p_val$", "^p$"),
  log2fc = c("^log2fc$", "^log2fc_x$", "^log2fc_y$", "^log2fc_6h$", "^log2fc_24h$", "^log2foldchange$", "^log2foldchange_x$", "^log2foldchange_y$", "^log2foldchange_6h$", "^log2foldchange_24h$", "^log2_fold_change$", "^log2_fold_change_x$", "^log2_fold_change_y$", "^log2_fold_change_6h$", "^log2_fold_change_24h$", "^logfc$", "^lfc$", "^avg_log2fc$", "^fold_change$"),
  pathway = c("^description$", "^term$", "^pathway$", "^pathway_name$", "^gs_description$", "^gs_name$", "^name$"),
  nes = c("^nes$", "^normalized_enrichment_score$", "^enrichment_score$", "^score$", "^combined_score$"),
  count = c("^n$", "^count$", "^count_x$", "^count_y$", "^freq$", "^size$", "^n_genes$", "^n_features$", "^n_lncrnas$", "^n_lncrnas$", "^n_mirnas$", "^n_circrnas$", "^n_targets$", "^total_targets$", "^overlap$", "^set_size$"),
  tf = c("^tf$", "^tf_key$", "^transcription_factor$", "^regulator$", "^motif$", "^tf_name$", "^factor$"),
  mirna = c("^mirna$", "^mirna_x$", "^mirna_y$", "^mirna_id$", "^mir$", "^mir_id$", "^mature_id$", "^mature_mirna$", "^small_rna$"),
  circ_id = c("^circ_id$", "^circ_id_x$", "^circ_id_y$", "^circrna_id$", "^circrna_id2$", "^circRNA_id$", "^stable_circ_id$", "^circular_rna_id$", "^circrna$"),
  energy = c("^e$", "^energy$", "^mfe$", "^inta_energy$", "^interaction_energy$"),
  significant = c("^significant$", "^sig$", "^is_significant$", "^de$", "^selected$", "^keep$"),
  ncRNA = c("^ncrna$", "^ncrna_x$", "^ncrna_y$", "^ncRNA$", "^lncrna$", "^lncrna_id$", "^lncrna_label$", "^circrna$", "^circ_id$", "^rna1$", "^rna_a$"),
  mrna = c("^mrna$", "^mrna_x$", "^mrna_y$", "^target_gene$", "^target_symbol$", "^target_gene_symbol$", "^target$", "^rna2$", "^rna_b$"),
  score = c("^score$", "^score_x$", "^score_y$", "^degree$", "^centrality$", "^weight$", "^rank_score$", "^importance$", "^hub_score$"),
  ncrna_direction = c("^ncrna_direction$", "^lncrna_direction$", "^circrna_direction$", "^circ_direction$", "^rna1_direction$"),
  mirna_direction = c("^mirna_direction$", "^mir_direction$", "^small_rna_direction$"),
  mrna_direction = c("^mrna_direction$", "^target_direction$", "^gene_direction$", "^rna2_direction$"),
  ncrna_log2fc = c("^ncrna_log2fc$", "^lncrna_log2fc$", "^circrna_log2fc$", "^circ_log2fc$", "^rna1_log2fc$"),
  mirna_log2fc = c("^mirna_log2fc$", "^mir_log2fc$", "^small_rna_log2fc$"),
  mrna_log2fc = c("^mrna_log2fc$", "^target_log2fc$", "^gene_log2fc$", "^rna2_log2fc$")
)

normalize_names <- function(x) {
  x <- gsub("\\.+", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

trim_ws <- function(x) {
  x <- as.character(x)
  x <- gsub("^\\s+|\\s+$", "", x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "<NA>")] <- NA_character_
  x
}

as_numeric_safe <- function(x) {
  if (is.null(x)) return(rep(NA_real_, 0))
  x <- trim_ws(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

coalesce_vec <- function(...) {
  xs <- list(...)
  if (length(xs) == 0) return(NULL)
  out <- xs[[1]]
  if (is.null(out)) return(NULL)
  is_blankish <- function(z) {
    if (is.null(z)) return(logical(0))
    if (is.character(z)) return(is.na(z) | trim_ws(z) == "")
    is.na(z)
  }
  for (i in seq_along(xs)[-1]) {
    yi <- xs[[i]]
    if (is.null(yi)) next
    idx <- is_blankish(out)
    out[idx] <- yi[idx]
  }
  out
}

first_nonempty <- function(x) {
  x <- trim_ws(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

truthy <- function(x) {
  y <- tolower(trim_ws(x))
  y %in% c("1", "true", "t", "yes", "y", "significant", "sig", "selected", "keep", "de")
}

find_matching_cols <- function(df, alias_key) {
  if (!alias_key %in% names(ALIASES)) return(character(0))
  nms <- names(df)
  pats <- ALIASES[[alias_key]]
  hits <- nms[vapply(nms, function(nm) any(vapply(pats, function(pt) grepl(pt, nm, ignore.case = TRUE, perl = TRUE), logical(1))), logical(1))]
  unique(hits)
}

family_present <- function(colnames_vec, alias_key) {
  if (!alias_key %in% names(ALIASES)) return(FALSE)
  nms <- normalize_names(colnames_vec)
  pats <- ALIASES[[alias_key]]
  any(vapply(nms, function(nm) any(vapply(pats, function(pt) grepl(pt, nm, ignore.case = TRUE, perl = TRUE), logical(1))), logical(1)))
}

coalesce_family <- function(df, alias_key) {
  nr <- if (is.null(df)) 0L else nrow(df)
  cols <- find_matching_cols(df, alias_key)
  if (length(cols) == 0) return(rep(NA_character_, nr))
  vals <- lapply(cols, function(cc) trim_ws(df[[cc]]))
  out <- vals[[1]]
  if (length(vals) > 1) {
    for (i in 2:length(vals)) {
      idx <- is.na(out) | out == ""
      if (length(vals[[i]]) == length(out)) out[idx] <- vals[[i]][idx]
    }
  }
  if (length(out) == 0 && nr > 0) out <- rep(NA_character_, nr)
  out
}

infer_time_from_text <- function(x) {
  x <- tolower(trim_ws(x))
  out <- rep(NA_character_, length(x))
  pat6 <- "(^|[^0-9])6[[:space:]]*(h|hr|hrs|hour|hours)([^a-z]|$)|6_h|6_hr|6_hrs|6_hour|6_hours|6hr|6hrs|6hour|6hours|ne6|ctrl6|only6|early|t6"
  pat24 <- "(^|[^0-9])24[[:space:]]*(h|hr|hrs|hour|hours)([^a-z]|$)|24_h|24_hr|24_hrs|24_hour|24_hours|24hr|24hrs|24hour|24hours|ne24|ctrl24|only24|late|t24"
  out[grepl(pat6, x, perl = TRUE)] <- "6h"
  out[grepl(pat24, x, perl = TRUE)] <- "24h"
  out
}

normalize_time <- function(x, source_name = NULL) {
  out <- trim_ws(x)
  out <- coalesce_vec(out, infer_time_from_text(out))
  if (!is.null(source_name)) {
    idx <- is.na(out) | out == ""
    out[idx] <- infer_time_from_text(rep(source_name, sum(idx)))
  }
  out <- ifelse(out %in% c("6", "6.0"), "6h", out)
  out <- ifelse(out %in% c("24", "24.0"), "24h", out)
  out
}

infer_class_from_text <- function(x) {
  x <- tolower(trim_ws(x))
  out <- rep(NA_character_, length(x))
  out[grepl("early|only6|6h_only|commit", x)] <- "Early"
  out[grepl("sustained|shared|both|persistent", x)] <- "Sustained/Shared"
  out[grepl("late|only24|24h_only|maint", x)] <- "Late"
  out[grepl("reversal|switch|opposite", x)] <- "Reversal"
  out
}

normalize_class <- function(x, source_name = NULL) {
  out <- trim_ws(x)
  norm <- infer_class_from_text(out)
  idx <- is.na(norm) | norm == ""
  if (!is.null(source_name)) norm[idx] <- infer_class_from_text(rep(source_name, sum(idx)))
  norm
}

infer_direction_from_text <- function(x) {
  x <- tolower(trim_ws(x))
  out <- rep(NA_character_, length(x))
  out[grepl("up|upregulated|activated|induced|positive", x)] <- "UP"
  out[grepl("down|downregulated|suppressed|repressed|negative", x)] <- "DOWN"
  out
}

normalize_direction <- function(x, log2fc = NULL, source_name = NULL, nes = NULL) {
  out <- infer_direction_from_text(x)
  if (!is.null(log2fc)) {
    lfc <- as_numeric_safe(log2fc)
    idx <- is.na(out) & !is.na(lfc)
    out[idx & lfc > 0] <- "UP"
    out[idx & lfc < 0] <- "DOWN"
  }
  if (!is.null(nes)) {
    nscore <- as_numeric_safe(nes)
    idx <- is.na(out) & !is.na(nscore)
    out[idx & nscore > 0] <- "UP"
    out[idx & nscore < 0] <- "DOWN"
  }
  if (!is.null(source_name)) {
    idx <- is.na(out)
    out[idx] <- infer_direction_from_text(rep(source_name, sum(idx)))
  }
  out
}

energy_tier <- function(x) {
  e <- as_numeric_safe(x)
  out <- rep(NA_character_, length(e))
  out[!is.na(e) & e <= -20] <- "High"
  out[!is.na(e) & e > -20 & e <= -15] <- "Moderate"
  out[!is.na(e) & e > -15 & e <= -10] <- "Exploratory"
  out
}

entity_sig_keep <- function(df, entity = c("mrna", "lnc", "mirna", "circ")) {
  entity <- match.arg(entity)
  if (is.null(df) || nrow(df) == 0) return(df)
  lfc_cut <- c(mrna = 1.0, lnc = 0.5, mirna = 0.5, circ = 0.5)[[entity]]
  keep <- df$std_significant
  if ("std_log2fc" %in% names(df)) {
    keep <- keep & (is.na(df$std_log2fc) | abs(df$std_log2fc) >= lfc_cut)
  }
  df[keep %in% TRUE, , drop = FALSE]
}


has_required_raw_de_cols <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  need <- c("gene_id", "log2foldchange", "padj")
  have <- normalize_names(names(df))
  all(need %in% have)
}

extract_raw_de_slot_tables <- function(slot, entity = c("mRNA", "lncRNA")) {
  entity <- match.arg(entity)
  out <- list()
  files <- slot$files %||% character(0)
  if (length(files) == 0) return(out)
  target_pat <- if (entity == "mRNA") "mrna_de" else "lncrna_de"

  for (f in files) {
    base_low <- tolower(basename(f))
    if (!grepl(target_pat, base_low, perl = TRUE)) next
    df <- read_table_safe(f)
    if (is.null(df) || nrow(df) == 0) next
    if (!has_required_raw_de_cols(df)) next
    key <- if (grepl("6h", base_low, perl = TRUE)) "6h" else if (grepl("24h", base_low, perl = TRUE)) "24h" else base_low
    out[[key]] <- standardize_df(df)
  }
  out
}

make_raw_de_candidates <- function(df, label_col = c("std_label_mrna", "std_label_lnc")) {
  label_col <- label_col[label_col %in% names(df)][1]
  if (is.na(label_col) || is.null(label_col)) label_col <- "std_gene_id"
  out <- tibble::tibble(
    timepoint = df$std_time,
    direction = df$std_direction,
    label = df[[label_col]],
    gene_id = df$std_gene_id,
    gene_symbol = df$std_gene_symbol,
    source = df$std_source,
    log2FC = df$std_log2fc,
    padj = df$std_padj,
    source_file = df$.source_name
  ) %>% distinct()
  out
}

assert_directional_summary <- function(df, figure_name = "unknown") {
  if (is.null(df) || nrow(df) == 0) return(invisible(TRUE))
  n_sig <- sum(!is.na(df$std_padj) & df$std_padj < 0.05 & !is.na(df$std_log2fc), na.rm = TRUE)
  n_dir <- sum(!is.na(df$std_direction) & df$std_direction %in% c("UP", "DOWN"), na.rm = TRUE)
  if (n_sig > 0 && n_dir == 0) {
    stop("Sanity check failed in ", figure_name, ": significant rows were found in a raw DE table but no UP/DOWN directions were recovered.")
  }
  invisible(TRUE)
}

looks_like_pathway_label <- function(x) {
  x <- trim_ws(x)
  if (length(x) == 0) return(logical(0))
  good <- !is.na(x) &
    !grepl("^ensrnog", x, ignore.case = TRUE) &
    !grepl("^mstrg\\.", x, ignore.case = TRUE) &
    !grepl("^loc\\d+", x, ignore.case = TRUE) &
    !grepl("^rno\\d{5}$", x, ignore.case = TRUE) &
    nchar(x) >= 4
  good
}

normalize_kegg_group <- function(x, source_name = NULL) {
  raw <- trim_ws(x)
  out <- rep(NA_character_, length(raw))
  low <- tolower(raw)
  out[grepl("only6|early", low)] <- "Early"
  out[grepl("shared|sustained|common", low)] <- "Sustained/Shared"
  out[grepl("only24|late", low)] <- "Late"
  out[grepl("reversal|switch|opposite", low)] <- "Reversal"
  if (!is.null(source_name)) {
    idx <- is.na(out)
    out[idx] <- infer_class_from_text(rep(source_name, sum(idx)))
  }
  out
}

extract_entity_direction <- function(df, which = c("ncrna", "mirna", "mrna")) {
  which <- match.arg(which)
  dir_key <- paste0(which, "_direction")
  lfc_key <- paste0(which, "_log2fc")
  raw_dir <- coalesce_family(df, dir_key)
  raw_lfc <- coalesce_family(df, lfc_key)
  normalize_direction(raw_dir, log2fc = raw_lfc, source_name = df$.source_name)
}

safe_nrow <- function(x) {
  if (is.null(x)) return(0L)
  nrow(x)
}

guess_delim <- function(file) {
  lines <- tryCatch(readLines(file, n = 5, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0) return("\t")
  count_pat <- function(pat) sum(vapply(gregexpr(pat, lines, fixed = TRUE), function(z) sum(z > 0), numeric(1)))
  scores <- c(tab = count_pat("\t"), semicolon = count_pat(";"), comma = count_pat(","))
  if (scores[["tab"]] >= max(scores)) return("\t")
  if (scores[["semicolon"]] >= max(scores)) return(";")
  ","
}

read_table_safe <- function(file) {
  ext <- tolower(tools::file_ext(file))
  out <- NULL

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      log_message("Skipping Excel file (readxl not installed): ", file)
      return(NULL)
    }
    out <- tryCatch(readxl::read_excel(file), error = function(e) NULL)
  } else {
    delims <- unique(c(guess_delim(file), "\t", ";", ","))
    parsed <- vector("list", length(delims))

    for (i in seq_along(delims)) {
      delim <- delims[[i]]
      attempt <- NULL
      if (requireNamespace("readr", quietly = TRUE)) {
        attempt <- tryCatch(
          readr::read_delim(file, delim = delim, show_col_types = FALSE, progress = FALSE, guess_max = 100000),
          error = function(e) NULL
        )
      }
      if (is.null(attempt)) {
        attempt <- tryCatch(
          utils::read.table(file, header = TRUE, sep = delim, quote = "\"", comment.char = "", fill = TRUE, check.names = FALSE),
          error = function(e) NULL
        )
      }
      parsed[[i]] <- attempt
    }

    score_parse <- function(x) {
      if (is.null(x)) return(-Inf)
      nc <- ncol(x)
      nr <- nrow(x)
      if (is.null(nc) || is.null(nr)) return(-Inf)
      bonus <- if (nc > 1) 1000 else 0
      bonus + nc * 10 + min(nr, 50)
    }

    scores <- vapply(parsed, score_parse, numeric(1))
    if (all(is.infinite(scores))) {
      out <- NULL
    } else {
      out <- parsed[[which.max(scores)]]
      if (!is.null(out) && ncol(out) <= 1) {
        log_message("Warning: suspicious one-column parse retained for file: ", file)
      }
    }
  }

  if (is.null(out)) {
    log_message("Could not read table: ", file)
    return(NULL)
  }

  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- make.unique(normalize_names(names(out)), sep = "_dup_")
  nr <- nrow(out)
  out$.source_file <- rep(normalizePath(file, mustWork = FALSE), nr)
  out$.source_name <- rep(basename(file), nr)
  out$.source_dir <- rep(basename(dirname(file)), nr)
  tibble::as_tibble(out)
}

list_table_files <- function(root, dirs) {
  files <- unlist(lapply(dirs, function(d) {
    full <- file.path(root, d)
    if (!dir.exists(full)) return(character(0))
    list.files(full, recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  files <- unique(files)
  files[grepl("\\.(csv|tsv|txt|tab|xls|xlsx)$", files, ignore.case = TRUE)]
}

score_file_for_slot <- function(file, slot_cfg, preview_df = NULL) {
  path_low <- tolower(normalizePath(file, mustWork = FALSE))
  base_low <- tolower(basename(file))
  score <- 0
  if (grepl("correlation|correlations|rewir|rewiring", base_low, perl = TRUE)) score <- score - 8
  if (grepl("ne_vs_ctrl.*(6h|24h).*(mrna|de)|mrna.*(6h|24h).*(de)|deseq2.*(6h|24h)", base_low, perl = TRUE)) score <- score + 6
  if (length(slot_cfg$include) > 0) {
    score <- score + 3 * sum(vapply(slot_cfg$include, function(p) grepl(p, path_low, perl = TRUE), logical(1)))
  }
  if (length(slot_cfg$exclude) > 0) {
    score <- score - 5 * sum(vapply(slot_cfg$exclude, function(p) grepl(p, path_low, perl = TRUE), logical(1)))
  }
  if (identical(slot_cfg$title, "Figure 8 - integrative ceRNA candidates")) {
    if (grepl("pair|pairs|debug|annotation_map|all|full|raw|complete|interaction|edge|edges", base_low, perl = TRUE)) score <- score - 10
    if (grepl("candidate|triplet|triplets|triple|priorit|summary|coherent", base_low, perl = TRUE)) score <- score + 10
  }
  if (!is.null(preview_df)) {
    nms <- names(preview_df)
    score <- score + 4 * sum(vapply(slot_cfg$wanted_families, function(fm) family_present(nms, fm), logical(1)))
    if (nrow(preview_df) > 0) score <- score + 1
  }
  score
}

standardize_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)

  df$std_gene_symbol <- coalesce_family(df, "gene_symbol")
  df$std_gene_id <- coalesce_vec(coalesce_family(df, "gene_id"), coalesce_family(df, "feature_id"))
  df$std_source <- coalesce_family(df, "source")
  df$std_biotype <- coalesce_family(df, "biotype")
  df$std_padj <- as_numeric_safe(coalesce_family(df, "padj"))
  df$std_pvalue <- as_numeric_safe(coalesce_family(df, "pvalue"))
  df$std_log2fc <- as_numeric_safe(coalesce_family(df, "log2fc"))
  df$std_pathway <- coalesce_family(df, "pathway")
  df$std_nes <- as_numeric_safe(coalesce_family(df, "nes"))
  df$std_count <- as_numeric_safe(coalesce_family(df, "count"))
  df$std_tf <- coalesce_family(df, "tf")
  df$std_mirna <- coalesce_family(df, "mirna")
  df$std_circ_id <- coalesce_family(df, "circ_id")
  df$std_energy <- as_numeric_safe(coalesce_family(df, "energy"))
  df$std_energy_tier <- energy_tier(df$std_energy)
  df$std_ncRNA <- coalesce_vec(coalesce_family(df, "ncRNA"), df$std_gene_id, df$std_circ_id)
  df$std_mrna <- coalesce_vec(coalesce_family(df, "mrna"), df$std_gene_symbol, df$std_gene_id)
  df$std_score <- as_numeric_safe(coalesce_family(df, "score"))

  raw_time <- coalesce_family(df, "time")
  raw_class <- coalesce_family(df, "class")
  raw_direction <- coalesce_family(df, "direction")

  df$std_time <- normalize_time(raw_time, source_name = first_nonempty(df$.source_name))
  df$std_class <- normalize_class(raw_class, source_name = first_nonempty(df$.source_name))
  df$std_direction <- normalize_direction(raw_direction, log2fc = df$std_log2fc, source_name = first_nonempty(df$.source_name), nes = df$std_nes)

  sig_cols <- find_matching_cols(df, "significant")
  sig_col <- if (length(sig_cols) > 0) coalesce_family(df, "significant") else NULL
  has_real_sig_col <- !is.null(sig_col) && any(!(is.na(trim_ws(sig_col)) | trim_ws(sig_col) == ""))
  if (has_real_sig_col) {
    df$std_significant <- truthy(sig_col)
  } else if (all(is.na(df$std_padj))) {
    df$std_significant <- TRUE
  } else {
    df$std_significant <- !is.na(df$std_padj) & df$std_padj < 0.05
  }

  df$std_label_mrna <- coalesce_vec(df$std_gene_symbol, df$std_gene_id)
  df$std_label_lnc <- coalesce_vec(df$std_gene_id, df$std_gene_symbol)
  df$std_label_circ <- coalesce_vec(df$std_circ_id, df$std_gene_id)
  df$std_label_tf <- coalesce_vec(df$std_tf, df$std_gene_symbol, df$std_gene_id)

  tibble::as_tibble(df)
}

select_slot_files <- function(slot_name) {
  slot_cfg <- SLOTS[[slot_name]]
  overrides <- FILE_OVERRIDES[[slot_name]] %||% character(0)

  if (length(overrides) > 0) {
    files <- file.path(PROJECT_ROOT, overrides)
    files <- files[file.exists(files)]
    if (length(files) > 0) {
      log_message("Using manual overrides for ", slot_name, ": ", paste(basename(files), collapse = "; "))
      return(unique(files))
    }
  }

  candidates <- list_table_files(PROJECT_ROOT, slot_cfg$dirs)
  if (length(candidates) == 0) return(character(0))

  previews <- lapply(candidates, read_table_safe)
  scores <- mapply(score_file_for_slot, file = candidates, preview_df = previews, MoreArgs = list(slot_cfg = slot_cfg))
  keep <- which(scores > 0)
  if (length(keep) == 0) return(character(0))

  ord <- order(scores[keep], decreasing = TRUE)
  selected <- candidates[keep][ord]

  if (identical(slot_name, "fig8_cerna")) {
    base_low <- tolower(basename(selected))
    preferred <- grepl("triplet|triplets|triple|priority|coherent", base_low, perl = TRUE) &
      !grepl("pair|pairs|debug|annotation_map|raw|all|full|complete|edge|edges", base_low, perl = TRUE)
    if (any(preferred)) {
      selected <- selected[preferred]
    } else {
      selected <- selected[!grepl("debug|annotation_map", base_low, perl = TRUE)]
    }
  }

  max_n <- if (identical(slot_name, "fig8_cerna")) min(3L, MAX_FILES_PER_SLOT) else MAX_FILES_PER_SLOT
  selected <- head(selected, max_n)
  unique(selected)
}

load_slot_data <- function(slot_name) {
  files <- select_slot_files(slot_name)
  if (length(files) == 0) {
    return(list(data = tibble(), files = character(0)))
  }
  dfs <- lapply(files, function(f) {
    x <- read_table_safe(f)
    if (is.null(x)) return(NULL)
    standardize_df(x)
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(list(data = tibble(), files = files))
  out <- bind_rows(dfs)
  list(data = out, files = files)
}

collapse_items <- function(x, n = 12) {
  x <- trim_ws(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0) return("none detected")
  paste(head(x, n), collapse = ", ")
}

prefer_top <- function(df, label_col, n = 12) {
  if (is.null(df) || nrow(df) == 0 || !label_col %in% names(df)) return(character(0))
  ord <- seq_len(nrow(df))
  if ("std_padj" %in% names(df)) {
    padj_vec <- df$std_padj
    if (is.null(padj_vec) || length(padj_vec) != nrow(df)) padj_vec <- rep(Inf, nrow(df))
    fc_vec <- if ("std_log2fc" %in% names(df)) df$std_log2fc else rep(0, nrow(df))
    if (is.null(fc_vec) || length(fc_vec) != nrow(df) || all(is.na(fc_vec))) fc_vec <- rep(0, nrow(df))
    ord <- order(padj_vec, -abs(fc_vec), na.last = TRUE)
  }
  if (length(ord) != nrow(df) || all(is.na(ord))) ord <- seq_len(nrow(df))
  vals <- trim_ws(df[[label_col]][ord])
  vals <- unique(vals[!is.na(vals) & vals != ""])
  vals[seq_len(min(n, length(vals)))]
}

write_table_if_any <- function(df, file) {
  if (is.null(df) || nrow(df) == 0) return(invisible(FALSE))
  utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(TRUE)
}

entity_filter <- function(df, entity = c("mrna", "lnc", "mirna", "circ", "any")) {
  entity <- match.arg(entity)
  if (nrow(df) == 0 || entity == "any") return(df)

  text_blob <- paste(
    trim_ws(df$std_biotype), trim_ws(df$std_source), trim_ws(df$.source_name),
    trim_ws(df$std_gene_id), trim_ws(df$std_gene_symbol), trim_ws(df$std_circ_id), sep = " | "
  )
  text_blob <- tolower(text_blob)

  keep <- rep(TRUE, nrow(df))
  if (entity == "mrna") {
    keep <- !grepl("lnc|noncoding|mirna|miRNA|circ|circrna|ncrna", text_blob, ignore.case = TRUE)
    if (all(!keep)) keep <- rep(TRUE, nrow(df))
  }
  if (entity == "lnc") {
    keep <- grepl("lnc|lncrna|noncoding|non_coding|novel|known", text_blob, ignore.case = TRUE)
    if (all(!keep) && all(!is.na(df$std_gene_id))) {
      keep <- !grepl("mir|mirna|circ|circrna", text_blob, ignore.case = TRUE)
    }
    if (all(!keep)) keep <- rep(TRUE, nrow(df))
  }
  if (entity == "mirna") keep <- grepl("mir|mirna", text_blob, ignore.case = TRUE) | !is.na(df$std_mirna)
  if (entity == "circ") keep <- grepl("circ", text_blob, ignore.case = TRUE) | !is.na(df$std_circ_id)

  df[keep, , drop = FALSE]
}

make_input_audit_row <- function(slot_name, files) {
  if (length(files) == 0) return(tibble(slot = slot_name, file = NA_character_))
  tibble(slot = slot_name, file = normalizePath(files, mustWork = FALSE))
}

audit_slot_structure <- function(slot_name, slot_obj) {
  if (is.null(slot_obj)) {
    log_message("Audit ", slot_name, ": slot object is NULL")
    return(invisible(NULL))
  }
  files_n <- if (!is.null(slot_obj$files)) length(slot_obj$files) else 0L
  data_n  <- if (!is.null(slot_obj$data) && is.data.frame(slot_obj$data)) nrow(slot_obj$data) else 0L
  cols <- if (!is.null(slot_obj$data) && is.data.frame(slot_obj$data)) names(slot_obj$data) else character(0)
  cols_preview <- if (length(cols) == 0) "none" else paste(head(cols, 12), collapse = ", ")
  log_message("Audit ", slot_name, ": files=", files_n, ", rows=", data_n, ", columns=", cols_preview)
  invisible(NULL)
}

log_df_state <- function(tag, df) {
  n <- if (is.null(df)) 0L else nrow(df)
  log_message(tag, ": ", n, " rows")
  invisible(n)
}

# =========================
# REPORT BUILDERS
# =========================
report_fig2 <- function(slot) {
  raw_list <- extract_raw_de_slot_tables(slot, "mRNA")
  raw <- bind_rows(raw_list)
  if (nrow(raw) == 0) {
    lines <- c(
      "## Figure 2 - mRNA differential expression at 6 h and 24 h",
      "No suitable raw mRNA differential-expression tables were found automatically.",
      paste0("Candidate files inspected: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
    )
    return(list(lines = lines, candidates = list(mrna = tibble())))
  }

  raw$std_time <- ifelse(is.na(raw$std_time), infer_time_from_text(raw$.source_name), raw$std_time)
  raw$std_time <- ifelse(is.na(raw$std_time), infer_time_from_text(raw$.source_file), raw$std_time)
  raw$std_direction <- ifelse(is.na(raw$std_direction), ifelse(raw$std_log2fc > 0, "UP", ifelse(raw$std_log2fc < 0, "DOWN", NA)), raw$std_direction)
  raw$std_label_mrna <- coalesce_vec(raw$std_gene_symbol, raw$std_gene_id)

  log_df_state("Figure 2 raw tables before entity filter", raw)
  raw <- entity_filter(raw, "mrna")
  log_df_state("Figure 2 raw tables after entity filter", raw)
  sig <- entity_sig_keep(raw, "mrna")
  log_df_state("Figure 2 rows after entity_sig_keep", sig)
  sig <- sig %>%
    filter(!is.na(std_time), !is.na(std_label_mrna)) %>%
    distinct(std_time, std_label_mrna, std_direction, .keep_all = TRUE)
  log_df_state("Figure 2 rows after time/label filtering", sig)

  assert_directional_summary(sig, "Figure 2")

  counts <- sig %>%
    filter(!is.na(std_direction)) %>%
    count(std_time, std_direction, name = "n")

  get_count <- function(tp, dir) {
    val <- counts %>% filter(std_time == tp, std_direction == dir) %>% pull(n)
    ifelse(length(val) == 0, 0, val[[1]])
  }

  s6_up <- get_count("6h", "UP"); s6_down <- get_count("6h", "DOWN")
  s24_up <- get_count("24h", "UP"); s24_down <- get_count("24h", "DOWN")

  top_6_up <- collapse_items(prefer_top(sig %>% filter(std_time == "6h", std_direction == "UP"), "std_label_mrna", 15))
  top_6_down <- collapse_items(prefer_top(sig %>% filter(std_time == "6h", std_direction == "DOWN"), "std_label_mrna", 15))
  top_24_up <- collapse_items(prefer_top(sig %>% filter(std_time == "24h", std_direction == "UP"), "std_label_mrna", 15))
  top_24_down <- collapse_items(prefer_top(sig %>% filter(std_time == "24h", std_direction == "DOWN"), "std_label_mrna", 15))

  lines <- c(
    "## Figure 2 - mRNA differential expression at 6 h and 24 h",
    paste0("Input files used: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; "))),
    paste0("At 6 h, the report detected ", s6_up + s6_down, " significant mRNAs (", s6_up, " upregulated and ", s6_down, " downregulated)."),
    paste0("Representative 6 h upregulated genes: ", top_6_up, "."),
    paste0("Representative 6 h downregulated genes: ", top_6_down, "."),
    paste0("At 24 h, the report detected ", s24_up + s24_down, " significant mRNAs (", s24_up, " upregulated and ", s24_down, " downregulated)."),
    paste0("Representative 24 h upregulated genes: ", top_24_up, "."),
    paste0("Representative 24 h downregulated genes: ", top_24_down, "."),
    "Interpretation anchor: this section should be written as the direct transcriptional counterpart of the phenotypic commitment at 6 h and the maintenance state at 24 h, without implying a full time-course clustering framework."
  )

  cand <- make_raw_de_candidates(sig, "std_label_mrna")

  list(lines = lines, candidates = list(mrna = cand))
}

report_fig3 <- function(slot) {
  df <- entity_filter(slot$data, "mrna")
  if (nrow(df) == 0) {
    lines <- c(
      "## Figure 3 - Early / Sustained / Late / Reversal mRNA classes",
      "No suitable mRNA phase/class tables were found automatically.",
      paste0("Candidate files inspected: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
    )
    return(list(lines = lines, candidates = list(mrna_phase = tibble())))
  }

  df <- df %>%
    mutate(std_class = ifelse(is.na(std_class), infer_class_from_text(.source_name), std_class),
           std_direction = ifelse(is.na(std_direction), ifelse(std_log2fc > 0, "UP", ifelse(std_log2fc < 0, "DOWN", NA)), std_direction)) %>%
    filter(!is.na(std_class), !is.na(std_label_mrna)) %>%
    distinct(std_class, std_label_mrna, std_direction, .keep_all = TRUE)

  counts <- df %>% count(std_class, std_direction, name = "n")
  get_count <- function(cls, dir) {
    val <- counts %>% filter(std_class == cls, std_direction == dir) %>% pull(n)
    ifelse(length(val) == 0, 0, val[[1]])
  }

  raw_de <- slot$data %>%
    filter(grepl("ne_vs_ctrl_.*_mrna_de", .source_name, ignore.case = TRUE)) %>%
    entity_filter("mrna")
  if (nrow(raw_de) > 0) {
    raw_de <- raw_de %>%
      mutate(std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time),
             std_direction = ifelse(is.na(std_direction), ifelse(std_log2fc > 0, "UP", ifelse(std_log2fc < 0, "DOWN", NA)), std_direction),
             std_label_mrna = coalesce_vec(std_label_mrna, std_gene_symbol, std_gene_id))
    raw_de <- entity_sig_keep(raw_de, "mrna") %>%
      filter(!is.na(std_label_mrna)) %>%
      distinct(std_time, std_label_mrna)
  }
  n_class_union <- nrow(df %>% distinct(std_label_mrna))
  n_de_union <- if (nrow(raw_de) > 0) nrow(raw_de %>% distinct(std_label_mrna)) else NA_integer_
  note_line <- NULL
  if (!is.na(n_de_union) && n_class_union > n_de_union) {
    note_line <- paste0("Audit note: the Figure 3 class table contains ", n_class_union, " unique genes, which exceeds the ", n_de_union, " unique significant mRNAs recovered from the raw Figure 2 DE tables. This suggests that the class counts may reflect a broader phase/program table rather than a strict DE-only partition, and these totals should be verified against the upstream Fig. 3 generation table before being used literally in the manuscript.")
    log_message("Figure 3 audit: class_union=", n_class_union, ", raw_DE_union=", n_de_union, " -> verify upstream class table definition.")
  }

  classes <- c("Early", "Sustained/Shared", "Late", "Reversal")
  lines <- c(
    "## Figure 3 - Early / Sustained / Late / Reversal mRNA classes",
    paste0("Input files used: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
  )

  for (cls in classes) {
    sub <- df %>% filter(std_class == cls)
    lines <- c(
      lines,
      paste0(cls, ": ", nrow(sub), " genes (", get_count(cls, "UP"), " UP; ", get_count(cls, "DOWN"), " DOWN)."),
      paste0("Representative ", cls, " genes: ", collapse_items(prefer_top(sub, "std_label_mrna", 18)), ".")
    )
  }
  if (!is.null(note_line)) lines <- c(lines, note_line)

  sustained <- df %>% filter(std_class == "Sustained/Shared")
  lines <- c(
    lines,
    paste0("Commitment-to-maintenance interpretation: Early genes should be discussed as commitment-associated responses at 6 h; Sustained/Shared genes as the transcriptional bridge between commitment and maintenance; and Late genes as maintenance-associated reinforcement or remodeling outputs."),
    paste0("Representative sustained/shared genes that may connect commitment with maintenance: ", collapse_items(prefer_top(sustained, "std_label_mrna", 20)), ".")
  )

  cand <- df %>%
    transmute(class = std_class, direction = std_direction, gene_symbol = std_label_mrna, gene_id = std_gene_id, timepoint = std_time, log2FC = std_log2fc, padj = std_padj, source_file = .source_name) %>%
    distinct()

  list(lines = lines, candidates = list(mrna_phase = cand))
}

report_fig4 <- function(slot_de, slot_net) {
  raw_list <- extract_raw_de_slot_tables(slot_de, "lncRNA")
  raw <- bind_rows(raw_list)
  net <- entity_filter(slot_net$data, "lnc")

  lines <- c(
    "## Figure 4 - lncRNA differential expression, response classes, and regulatory layer",
    paste0("Input DE files used: ", ifelse(length(slot_de$files) == 0, "none", paste(basename(slot_de$files), collapse = "; "))),
    paste0("Input network files used: ", ifelse(length(slot_net$files) == 0, "none", paste(basename(slot_net$files), collapse = "; ")))
  )

  cand_de <- tibble()
  cand_hubs <- tibble()

  if (nrow(raw) > 0) {
    raw$std_time <- ifelse(is.na(raw$std_time), infer_time_from_text(raw$.source_name), raw$std_time)
    raw$std_time <- ifelse(is.na(raw$std_time), infer_time_from_text(raw$.source_file), raw$std_time)
    raw$std_direction <- ifelse(is.na(raw$std_direction), ifelse(raw$std_log2fc > 0, "UP", ifelse(raw$std_log2fc < 0, "DOWN", NA)), raw$std_direction)
    raw$std_label_lnc <- coalesce_vec(raw$std_gene_id, raw$std_gene_symbol)

    log_df_state("Figure 4 raw tables before entity filter", raw)
    raw <- entity_filter(raw, "lnc")
    log_df_state("Figure 4 raw tables after entity filter", raw)
    sig <- entity_sig_keep(raw, "lnc")
    log_df_state("Figure 4 rows after entity_sig_keep", sig)
    sig <- sig %>%
      filter(!is.na(std_time), !is.na(std_label_lnc)) %>%
      distinct(std_time, std_label_lnc, std_direction, .keep_all = TRUE)
    log_df_state("Figure 4 rows after time/label filtering", sig)

    assert_directional_summary(sig, "Figure 4")

    counts_time <- sig %>% filter(!is.na(std_direction)) %>% count(std_time, std_direction, name = "n")
    get_count <- function(tp, dir) {
      val <- counts_time %>% filter(std_time == tp, std_direction == dir) %>% pull(n)
      ifelse(length(val) == 0, 0, val[[1]])
    }

    src_tab <- sig %>% filter(!is.na(std_source)) %>% count(std_source, sort = TRUE)

    lines <- c(
      lines,
      paste0("At 6 h, the report detected ", get_count("6h", "UP") + get_count("6h", "DOWN"), " significant lncRNAs (", get_count("6h", "UP"), " UP; ", get_count("6h", "DOWN"), " DOWN)."),
      paste0("Representative 6 h lncRNAs (gene_id labels): ", collapse_items(prefer_top(sig %>% filter(std_time == "6h"), "std_label_lnc", 20)), "."),
      paste0("At 24 h, the report detected ", get_count("24h", "UP") + get_count("24h", "DOWN"), " significant lncRNAs (", get_count("24h", "UP"), " UP; ", get_count("24h", "DOWN"), " DOWN)."),
      paste0("Representative 24 h lncRNAs (gene_id labels): ", collapse_items(prefer_top(sig %>% filter(std_time == "24h"), "std_label_lnc", 20)), ".")
    )

    if (nrow(src_tab) > 0) {
      src_text <- paste(paste0(src_tab$std_source, "=", src_tab$n), collapse = ", ")
      lines <- c(lines, paste0("Detected lncRNA source/annotation composition in the raw DE tables: ", src_text, "."))
    }

    cand_de <- sig %>%
      transmute(timepoint = std_time, direction = std_direction, gene_id = std_label_lnc, gene_symbol = std_gene_symbol, source = std_source, log2FC = std_log2fc, padj = std_padj, source_file = .source_name) %>%
      distinct()
  } else {
    lines <- c(lines, "No raw lncRNA differential-expression tables could be summarized automatically.")
  }

  phase_df <- entity_filter(slot_de$data, "lnc")
  log_df_state("Figure 4 phase/class rows after entity filter", phase_df)
  if (nrow(phase_df) > 0) {
    phase_df <- phase_df %>%
      mutate(std_class = ifelse(is.na(std_class), infer_class_from_text(.source_name), std_class),
             std_direction = ifelse(is.na(std_direction), ifelse(std_log2fc > 0, "UP", ifelse(std_log2fc < 0, "DOWN", NA)), std_direction),
             std_label_lnc = coalesce_vec(std_label_lnc, std_gene_id, std_gene_symbol)) %>%
      filter(!is.na(std_class), !is.na(std_label_lnc))

    class_tab <- phase_df %>% count(std_class, std_direction, name = "n")
    if (nrow(class_tab) > 0) {
      cl_text <- paste(paste0(class_tab$std_class, " ", class_tab$std_direction, "=", class_tab$n), collapse = ", ")
      lines <- c(lines, paste0("Response-class composition among lncRNAs: ", cl_text, "."))
    }

    for (cls in c("Early", "Sustained/Shared", "Late", "Reversal")) {
      sub_cls <- phase_df %>% filter(std_class == cls)
      if (nrow(sub_cls) == 0) next
      lines <- c(lines, paste0("Representative ", cls, " lncRNAs (gene_id labels): ", collapse_items(prefer_top(sub_cls, "std_label_lnc", 18)), "."))
    }
  }

  log_df_state("Figure 4 network rows after entity filter", net)
  if (nrow(net) > 0) {
    hub_like <- net %>%
      mutate(std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time),
             std_label_lnc = coalesce_vec(std_label_lnc, std_gene_id, std_gene_symbol),
             hub_metric = coalesce_vec(std_score, std_count)) %>%
      filter(!is.na(std_label_lnc))

    top_6_hubs <- hub_like %>% filter(std_time == "6h") %>% arrange(desc(hub_metric)) %>% slice_head(n = 15)
    top_24_hubs <- hub_like %>% filter(std_time == "24h") %>% arrange(desc(hub_metric)) %>% slice_head(n = 15)

    lines <- c(
      lines,
      paste0("Top 6 h lncRNA hubs or high-connectivity candidates: ", collapse_items(top_6_hubs$std_label_lnc, 15), "."),
      paste0("Top 24 h lncRNA hubs or high-connectivity candidates: ", collapse_items(top_24_hubs$std_label_lnc, 15), ".")
    )

    rewire_hits <- hub_like %>% filter(grepl("rewir", .source_name, ignore.case = TRUE) | grepl("rewir", .source_dir, ignore.case = TRUE))
    if (nrow(rewire_hits) > 0) {
      lines <- c(lines, paste0("Potential rewiring candidates detected from rewiring-oriented tables: ", collapse_items(rewire_hits$std_label_lnc, 15), "."))
    }

    cand_hubs <- hub_like %>%
      transmute(timepoint = std_time, gene_id = std_label_lnc, score = hub_metric, source = std_source, source_file = .source_name) %>%
      distinct()
  } else {
    lines <- c(lines, "No lncRNA hub/rewiring tables could be summarized automatically.")
  }

  lines <- c(lines, "Narrative constraint: keep gene_id as the primary label for lncRNAs in the report and manuscript-facing text.")
  list(lines = lines, candidates = list(lnc_de = cand_de, lnc_hubs = cand_hubs))
}

report_fig5 <- function(slot) {
  df <- slot$data
  lines <- c(
    "## Figure 5 - KEGG pathway programs across commitment and maintenance",
    paste0("Input files used: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
  )

  if (nrow(df) == 0) {
    lines <- c(lines, "No KEGG/GSEA pathway tables were found automatically.")
    return(list(lines = lines, candidates = list(pathways = tibble())))
  }

  kegg_group_raw <- coalesce_family(df, "class")
  if ("kegg_group" %in% names(df)) kegg_group_raw <- coalesce_vec(trim_ws(df$kegg_group), kegg_group_raw)
  if ("kegg_group_x" %in% names(df)) kegg_group_raw <- coalesce_vec(trim_ws(df$kegg_group_x), kegg_group_raw)
  if ("kegg_group_y" %in% names(df)) kegg_group_raw <- coalesce_vec(trim_ws(df$kegg_group_y), kegg_group_raw)

  transition_raw <- rep(NA_character_, nrow(df))
  if ("transition_shared" %in% names(df)) transition_raw <- trim_ws(df$transition_shared)

  df <- df %>%
    mutate(
      std_class = normalize_kegg_group(kegg_group_raw, source_name = .source_name),
      std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time),
      std_direction = ifelse(is.na(std_direction), ifelse(std_nes > 0, "UP", ifelse(std_nes < 0, "DOWN", NA)), std_direction),
      pathway_label = std_pathway,
      transition_shared = transition_raw
    ) %>%
    filter(!is.na(pathway_label), looks_like_pathway_label(pathway_label))

  if (nrow(df) == 0) {
    lines <- c(lines, "KEGG files were detected, but no valid pathway labels could be recovered. Check the selected input files and pathway columns.")
    return(list(lines = lines, candidates = list(pathways = tibble())))
  }

  phases <- c("Early", "Sustained/Shared", "Late", "Reversal")
  for (ph in phases) {
    sub <- df %>% filter(std_class == ph)
    if (nrow(sub) == 0) next
    top_up <- sub %>% filter(std_direction == "UP" | (is.na(std_direction) & !is.na(std_nes) & std_nes > 0)) %>% arrange(std_padj, desc(std_nes)) %>% slice_head(n = 10)
    top_dn <- sub %>% filter(std_direction == "DOWN" | (is.na(std_direction) & !is.na(std_nes) & std_nes < 0)) %>% arrange(std_padj, std_nes) %>% slice_head(n = 10)
    lines <- c(
      lines,
      paste0(ph, " pathway program: ", nrow(sub), " rows detected."),
      paste0("Top activated pathways in ", ph, ": ", collapse_items(top_up$pathway_label, 10), "."),
      paste0("Top suppressed pathways in ", ph, ": ", collapse_items(top_dn$pathway_label, 10), ".")
    )
  }

  shared_like <- df %>% filter(std_class == "Sustained/Shared")
  if (nrow(shared_like) > 0) {
    tr_tab <- shared_like %>% filter(!is.na(transition_shared)) %>% count(transition_shared, sort = TRUE)
    if (nrow(tr_tab) > 0) {
      tr_text <- paste(paste0(tr_tab$transition_shared, "=", tr_tab$n), collapse = ", ")
      lines <- c(lines, paste0("Shared-pathway transition composition: ", tr_text, "."))
    }
    lines <- c(lines, "Shared-pathway scatter interpretation: pathways classified as Sustained/Shared should be discussed as programs that persist from the 6 h commitment stage into the 24 h maintenance state, rather than as independent late-only events.")
  }

  cand <- df %>%
    transmute(class = std_class, timepoint = std_time, direction = std_direction, pathway = pathway_label, NES = std_nes, padj = std_padj, n_genes = std_count, transition_shared = transition_shared, source_file = .source_name) %>%
    distinct()

  list(lines = lines, candidates = list(pathways = cand))
}

report_fig6 <- function(slot) {
  df <- slot$data
  lines <- c(
    "## Figure 6 - TF layer controlling commitment vs maintenance",
    paste0("Input files used: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
  )

  if (nrow(df) == 0) {
    lines <- c(lines, "No TF tables were found automatically.")
    return(list(lines = lines, candidates = list(tfs = tibble())))
  }

  df <- df %>%
    mutate(tf_label = coalesce_vec(std_tf, std_gene_symbol, std_gene_id),
           std_class = ifelse(is.na(std_class), infer_class_from_text(.source_name), std_class),
           std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time)) %>%
    filter(!is.na(tf_label))

  for (ph in c("Early", "Sustained/Shared", "Late", "Reversal")) {
    sub <- df %>% filter(std_class == ph)
    if (nrow(sub) == 0) next
    top_tf <- sub %>% arrange(std_padj, desc(std_score), desc(std_count)) %>% slice_head(n = 15)
    lines <- c(lines, paste0("TFs standing out in ", ph, " modules/programs: ", collapse_items(top_tf$tf_label, 15), "."))
  }

  if (all(is.na(df$std_class)) && any(!is.na(df$std_time))) {
    for (tp in c("6h", "24h")) {
      sub <- df %>% filter(std_time == tp)
      if (nrow(sub) == 0) next
      top_tf <- sub %>% arrange(std_padj, desc(std_score), desc(std_count)) %>% slice_head(n = 15)
      lines <- c(lines, paste0("TFs standing out at ", tp, ": ", collapse_items(top_tf$tf_label, 15), "."))
    }
  }

  lines <- c(lines, "Interpretation anchor: TFs should be described as candidate upstream controllers of commitment-associated versus maintenance-associated gene programs, not as definitive master regulators unless independently validated.")

  cand <- df %>%
    transmute(class = std_class, timepoint = std_time, tf = tf_label, score = std_score, n_targets = std_count, padj = std_padj, source_file = .source_name) %>%
    distinct()

  list(lines = lines, candidates = list(tfs = cand))
}

report_fig7 <- function(slot) {
  df <- slot$data
  lines <- c(
    "## Figure 7 - miRNA-mRNA layer",
    paste0("Input files used: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
  )

  if (nrow(df) == 0) {
    lines <- c(lines, "No miRNA-mRNA interaction tables were found automatically.")
    return(list(lines = lines, candidates = list(mirna_mrna = tibble())))
  }

  log_message("Figure 7 audit: columns seen -> ", paste(names(df), collapse = ", "))

  raw_gene_key <- if ("gene_key" %in% names(df)) trim_ws(df$gene_key) else rep(NA_character_, nrow(df))
  raw_path_id <- if ("id" %in% names(df)) trim_ws(df$id) else rep(NA_character_, nrow(df))
  raw_path_desc <- coalesce_vec(if ("description" %in% names(df)) trim_ws(df$description) else NULL, df$std_pathway)

  raw_mirna_col <- if ("mirna" %in% names(df)) trim_ws(df$mirna) else NULL
  df <- df %>%
    mutate(
      miRNA = coalesce_vec(std_mirna, raw_mirna_col),
      mRNA_gene = coalesce_vec(std_gene_symbol, raw_gene_key),
      pathway_id = raw_path_id,
      pathway_description = raw_path_desc,
      std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time),
      std_class = ifelse(is.na(std_class), infer_class_from_text(.source_name), std_class)
    )

  valid_gene <- is_valid_gene_label(df$mRNA_gene)
  valid_pathway <- !is.na(trim_ws(df$pathway_description)) | is_kegg_id(df$pathway_id)
  valid_mi <- !is.na(trim_ws(df$miRNA))

  log_message("Figure 7 audit: rows total=", nrow(df), ", miRNA_nonmissing=", sum(valid_mi), ", valid_gene_labels=", sum(valid_gene), ", valid_pathway_annotations=", sum(valid_pathway))
  if (sum(valid_gene) == 0 && sum(valid_pathway) > 0) {
    log_message("Figure 7 audit warning: no recoverable gene-level labels found; table appears pathway-expanded only.")
  }

  df_gene <- df %>% filter(valid_mi, valid_gene)
  df_path <- df %>% filter(valid_mi, valid_pathway)

  pair_rank <- df_gene %>% distinct(miRNA, mRNA_gene, std_time, std_class, .keep_all = TRUE) %>% count(miRNA, mRNA_gene, std_time, std_class, sort = TRUE, name = "support_rows")
  mir_rank <- if (nrow(df_gene) > 0) df_gene %>% distinct(miRNA, mRNA_gene) %>% count(miRNA, sort = TRUE, name = "n_targets") else tibble(miRNA = character(), n_targets = integer())
  mrna_rank <- if (nrow(df_gene) > 0) df_gene %>% distinct(miRNA, mRNA_gene) %>% count(mRNA_gene, sort = TRUE, name = "n_miRNAs") else tibble(mRNA_gene = character(), n_miRNAs = integer())
  top_pairs_24 <- pair_rank %>% filter(std_time == "24h" | std_class %in% c("Late", "Sustained/Shared")) %>% slice_head(n = 20)
  top_pairs_6 <- pair_rank %>% filter(std_time == "6h" | std_class == "Early") %>% slice_head(n = 20)

  path_counts <- if (nrow(df_path) > 0) {
    df_path %>%
      mutate(pathway_id = ifelse(is.na(pathway_id), NA_character_, pathway_id),
             pathway_description = ifelse(is.na(pathway_description), pathway_id, pathway_description)) %>%
      distinct(combo, std_time, std_class, pathway_id, pathway_description, .keep_all = TRUE) %>%
      count(combo, std_time, std_class, pathway_id, pathway_description, sort = TRUE, name = "n_edges")
  } else {
    tibble()
  }

  total_pairs <- nrow(df_gene %>% distinct(miRNA, mRNA_gene, std_time, std_class))
  n_mir <- nrow(df_gene %>% distinct(miRNA))
  n_mrna <- nrow(df_gene %>% distinct(mRNA_gene))

  lines <- c(lines, paste0("Audit note: valid gene-level mRNA labels were recovered for ", sum(valid_gene), " of ", nrow(df), " rows; pathway annotations were recovered for ", sum(valid_pathway), " rows."))

  if (nrow(df_gene) > 0) {
    lines <- c(
      lines,
      paste0("The report detected ", total_pairs, " distinct miRNA-mRNA candidate pairs after collapsing repeated rows across files, involving ", n_mir, " miRNAs and ", n_mrna, " mRNA targets."),
      paste0("Most recurrent miRNAs across the interaction layer: ", collapse_items(mir_rank$miRNA, 15), "."),
      paste0("Most recurrent mRNA targets across the interaction layer: ", collapse_items(mrna_rank$mRNA_gene, 15), ".")
    )
    if (nrow(top_pairs_24) > 0) lines <- c(lines, paste0("Representative 24 h / maintenance-relevant miRNA-mRNA pairs: ", collapse_items(paste(top_pairs_24$miRNA, top_pairs_24$mRNA_gene, sep = "->"), 20), "."))
    if (nrow(top_pairs_6) > 0) lines <- c(lines, paste0("Representative 6 h / commitment-relevant miRNA-mRNA pairs: ", collapse_items(paste(top_pairs_6$miRNA, top_pairs_6$mRNA_gene, sep = "->"), 20), "."))
  } else {
    lines <- c(lines, "No trustworthy gene-level miRNA-mRNA labels could be recovered from the selected Figure 7 tables, so gene-pair summaries were not reported.")
  }

  if (nrow(path_counts) > 0) {
    for (cb in unique(path_counts$combo)) {
      txt <- path_counts %>% filter(combo == cb) %>% mutate(label = paste0(pathway_description, " [", pathway_id, "; edges=", n_edges, "]")) %>% pull(label)
      lines <- c(lines, paste0("Pathway annotations represented in ", cb, ": ", paste(head(txt, 20), collapse = "; "), "."))
    }
  }

  lines <- c(lines, "Interpretation anchor: unless 6 h evidence is strong and explicit, this layer should be framed primarily as maintenance reinforcement rather than as the primary initiator of the phenotype. If the selected table is pathway-expanded, pathway associations should not be mislabeled as direct gene identities.")

  cand <- df %>%
    transmute(timepoint = std_time, class = std_class, miRNA = miRNA, mRNA_gene = ifelse(is_valid_gene_label(mRNA_gene), mRNA_gene, NA_character_), pathway_id = pathway_id, pathway_description = pathway_description, source_file = .source_name) %>%
    distinct()

  list(lines = lines, candidates = list(mirna_mrna = cand))
}

report_supp_circ <- function(slot) {
  df <- entity_filter(slot$data, "circ")
  lines <- c(
    "## Supplementary circRNA note",
    paste0("Input files used: ", ifelse(length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
  )

  if (nrow(df) == 0) {
    lines <- c(lines, "No circRNA-oriented supplementary tables were summarized automatically.")
    return(list(lines = lines, candidates = list(circ = tibble())))
  }

  df <- df %>%
    mutate(std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time),
           std_class = ifelse(is.na(std_class), infer_class_from_text(.source_name), std_class),
           std_direction = ifelse(is.na(std_direction), ifelse(std_log2fc > 0, "UP", ifelse(std_log2fc < 0, "DOWN", NA)), std_direction)) %>%
    filter(!is.na(std_label_circ))

  lines <- c(
    lines,
    paste0("Distinct circRNA-labelled candidates detected: ", nrow(df %>% distinct(std_label_circ, std_time, std_class, .keep_all = FALSE)), "."),
    paste0("Representative circRNA IDs: ", collapse_items(prefer_top(df, "std_label_circ", 20)), ".")
  )

  cand <- df %>%
    transmute(timepoint = std_time, class = std_class, direction = std_direction, circ_id = std_label_circ, log2FC = std_log2fc, padj = std_padj, source_file = .source_name) %>%
    distinct()

  list(lines = lines, candidates = list(circ = cand))
}

report_fig8 <- function(slot) {
  lines <- c(
    "## Figure 8 - integrative ceRNA candidates",
    paste0("Input files used: ", ifelse(is.null(slot$files) || length(slot$files) == 0, "none", paste(basename(slot$files), collapse = "; ")))
  )

  if (is.null(slot) || length(slot) == 0 || is.null(slot$data) || nrow(slot$data) == 0) {
    lines <- c(
      lines,
      "No ceRNA candidate tables were found automatically.",
      "Conservative interpretation: this layer should not be emphasized unless coherent triplet-level evidence is clearly available."
    )
    return(list(lines = lines, candidates = list(cerna = tibble())))
  }

  df <- slot$data
  log_message("Figure 8 audit: columns seen -> ", paste(names(df), collapse = ", "))

  source_names <- if (".source_name" %in% names(df)) tolower(trim_ws(df$.source_name)) else rep(NA_character_, nrow(df))
  keep_triplet_sources <- grepl("triplet|triplets|triple|priority|coherent", source_names, perl = TRUE) &
    !grepl("pair|pairs|debug|annotation_map|raw|all|full|complete|edge|edges", source_names, perl = TRUE)

  if (any(keep_triplet_sources, na.rm = TRUE)) {
    df <- df[keep_triplet_sources %in% TRUE, , drop = FALSE]
    lines <- c(lines, "Figure 8 was restricted to triplet/priority candidate files to avoid inflating the report with auxiliary pair/debug tables.")
  } else {
    lines <- c(lines, "No clearly triplet-focused ceRNA file subset was detected; the report falls back to all selected ceRNA tables.")
  }

  if (nrow(df) == 0) {
    lines <- c(lines, "After restricting to triplet/priority candidate files, no ceRNA rows remained.")
    return(list(lines = lines, candidates = list(cerna = tibble())))
  }

  raw_ncrna_label <- if ("ncrna_label" %in% names(df)) trim_ws(df$ncrna_label) else NULL
  raw_mirna_col <- if ("mirna" %in% names(df)) trim_ws(df$mirna) else NULL
  raw_mrna_gene_symbol <- if ("mRNA_gene_symbol" %in% names(df)) trim_ws(df$mRNA_gene_symbol) else NULL
  raw_gene_symbol <- if ("gene_symbol" %in% names(df)) trim_ws(df$gene_symbol) else NULL
  raw_gene_key <- if ("gene_key" %in% names(df)) trim_ws(df$gene_key) else NULL
  df <- df %>%
    mutate(
      ncRNA = coalesce_vec(raw_ncrna_label, std_ncRNA, std_circ_id),
      miRNA = coalesce_vec(std_mirna, raw_mirna_col),
      mRNA_gene = coalesce_vec(raw_mrna_gene_symbol, raw_gene_symbol, raw_gene_key),
      std_time = ifelse(is.na(std_time), infer_time_from_text(.source_name), std_time),
      std_class = ifelse(is.na(std_class), infer_class_from_text(.source_name), std_class)
    )

  ncrna_dir <- extract_entity_direction(df, "ncrna")
  mirna_dir <- extract_entity_direction(df, "mirna")
  mrna_dir <- extract_entity_direction(df, "mrna")
  df$ncRNA_direction <- ncrna_dir
  df$miRNA_direction <- mirna_dir
  df$mRNA_direction <- mrna_dir

  valid_nc <- !is.na(trim_ws(df$ncRNA))
  valid_mi <- !is.na(trim_ws(df$miRNA))
  valid_gene <- is_valid_gene_label(df$mRNA_gene)

  log_message("Figure 8 audit: rows total=", nrow(df), ", ncRNA_nonmissing=", sum(valid_nc), ", miRNA_nonmissing=", sum(valid_mi), ", valid_mRNA_gene_labels=", sum(valid_gene))

  coherent_all <- df %>%
    filter(valid_nc, valid_mi, !is.na(ncRNA_direction), !is.na(miRNA_direction), !is.na(mRNA_direction)) %>%
    mutate(coherent_cerna = ncRNA_direction == mRNA_direction & ncRNA_direction != miRNA_direction) %>%
    filter(coherent_cerna)

  coherent_all$evidence_tier <- energy_tier(coherent_all$std_energy)
  all_tier_tab <- as.data.frame(table(coherent_all$evidence_tier), stringsAsFactors = FALSE)
  if (nrow(all_tier_tab) > 0) names(all_tier_tab) <- c("evidence_tier", "n")
  tier_text <- if (nrow(all_tier_tab) == 0) "no energy tiering available" else paste(paste0(all_tier_tab$evidence_tier, "=", all_tier_tab$n), collapse = ", ")

  lines <- c(
    lines,
    paste0("The report detected ", nrow(coherent_all), " directionally coherent ceRNA triplets after applying the conservative logic ncRNA UP - miRNA DOWN - mRNA UP or ncRNA DOWN - miRNA UP - mRNA DOWN."),
    paste0("Energy-evidence tier composition among coherent triplets: ", tier_text, "."),
    paste0("Audit note: recoverable mRNA gene labels were available for ", sum(valid_gene), " of ", nrow(df), " rows in the exported triplet tables.")
  )

  if (nrow(coherent_all) > 0) {
    class_tab <- as.data.frame(table(coherent_all$std_class), stringsAsFactors = FALSE)
    if (nrow(class_tab) > 0) {
      names(class_tab) <- c("std_class", "n")
      class_tab <- class_tab[class_tab$std_class != "" & !is.na(class_tab$std_class), , drop = FALSE]
      if (nrow(class_tab) > 0) lines <- c(lines, paste0("Coherent triplet distribution by response class: ", paste(paste0(class_tab$std_class, "=", class_tab$n), collapse = ", "), "."))
    }
  }

  coherent_labeled <- coherent_all %>%
    mutate(mRNA_gene = ifelse(is_valid_gene_label(mRNA_gene), mRNA_gene, NA_character_)) %>%
    filter(!is.na(mRNA_gene), ncRNA != mRNA_gene) %>%
    mutate(triplet = paste(ncRNA, miRNA, mRNA_gene, sep = " | "))

  if (nrow(coherent_labeled) > 0) {
    padj_ord <- if ("std_padj" %in% names(coherent_labeled)) coherent_labeled$std_padj else rep(Inf, nrow(coherent_labeled))
    coherent_labeled <- coherent_labeled[order(match(coherent_labeled$evidence_tier, c("High", "Moderate", "Exploratory")), padj_ord, na.last = TRUE), , drop = FALSE]
    coherent_labeled <- coherent_labeled[!duplicated(coherent_labeled$triplet), , drop = FALSE]
    lines <- c(lines, paste0("Representative coherent triplets with recoverable mRNA gene labels: ", collapse_items(coherent_labeled$triplet, 20), "."))
  } else {
    lines <- c(lines, "Representative labeled ceRNA triplets could not be reported because the exported triplet tables do not currently contain recoverable mRNA gene labels (the mRNA-side gene_symbol/gene_key fields are empty in these files). This should be fixed upstream if triplet examples are needed for the manuscript.")
  }

  lines <- c(lines, "Conservative interpretation: emphasize only biologically meaningful triplets that are both directionally coherent and supported by explicit mRNA labels. If the exported triplet files lack recoverable mRNA identities, this layer should be described as provisional and not overstated.")

  cand <- coherent_labeled %>%
    transmute(timepoint = std_time, class = std_class, ncRNA = ncRNA, ncRNA_direction = ncRNA_direction, miRNA = miRNA, miRNA_direction = miRNA_direction, mRNA_gene = mRNA_gene, mRNA_direction = mRNA_direction, energy = std_energy, evidence_tier = evidence_tier, source_file = .source_name) %>%
    distinct()

  list(lines = lines, candidates = list(cerna = tibble::as_tibble(cand)))
}

report_cross_layer <- function(candidate_store) {
  lines <- c("## Cross-layer recurrence and integrative notes")

  collect_gene_vectors <- function(tbls, cols) {
    vec <- character(0)
    for (tb in tbls) {
      if (is.null(tb) || nrow(tb) == 0) next
      for (cc in cols) {
        if (cc %in% names(tb)) vec <- c(vec, trim_ws(tb[[cc]]))
      }
    }
    vec <- vec[!is.na(vec)]
    vec
  }

  gene_vec <- collect_gene_vectors(candidate_store[c("mrna", "mrna_phase", "pathways", "tfs", "mirna_mrna", "cerna")], c("gene_symbol", "mRNA", "mRNA_gene"))
  gene_vec <- gene_vec[is_valid_gene_label(gene_vec)]
  mir_vec <- collect_gene_vectors(candidate_store[c("mirna_mrna", "cerna")], c("miRNA"))
  lnc_vec <- collect_gene_vectors(candidate_store[c("lnc_de", "lnc_hubs", "cerna")], c("gene_id", "ncRNA"))
  circ_vec <- collect_gene_vectors(candidate_store[c("cerna")], c("ncRNA"))
  circ_vec <- circ_vec[grepl("^[0-9XYM]+:[0-9]+[|][0-9]+$", circ_vec, ignore.case = TRUE)]

  make_top_freq <- function(x, n = 15) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return("none detected")
    tb <- sort(table(x), decreasing = TRUE)
    paste(paste0(names(head(tb, n)), " (", as.integer(head(tb, n)), ")"), collapse = ", ")
  }

  lines <- c(
    lines,
    paste0("Recurrent mRNA candidates across layers: ", make_top_freq(gene_vec), "."),
    paste0("Recurrent miRNA candidates across layers: ", make_top_freq(mir_vec), "."),
    paste0("Recurrent lncRNA / ncRNA candidates across layers: ", make_top_freq(lnc_vec), ".")
  )

  if (length(circ_vec) > 0) {
    lines <- c(lines, paste0("circRNA-labelled recurrent candidates: ", make_top_freq(circ_vec), "."))
  } else {
    lines <- c(lines, "circRNA-labelled recurrent candidates: none detected.")
  }

  lines <- c(lines, "Use this section conservatively. Cross-layer recurrence is useful for prioritization, but should not be presented as independent validation unless the same candidates are supported by clearly interpretable gene-level evidence across layers.")

  list(lines = lines)
}

# =========================
# MAIN EXECUTION
# =========================
log_message("Project root: ", PROJECT_ROOT)
log_message("Report directory: ", REPORT_DIR)

slot_order <- c("fig2_mrna", "fig3_mrna_phase", "fig4_lnc_de", "fig4_lnc_network", "fig5_kegg", "fig6_tf", "fig7_mirna_mrna", "fig8_cerna", "supp_circ")

# Defensive helper re-definitions in case the script is run in a partially sourced session
if (!exists("infer_time_from_text", mode = "function")) {
  infer_time_from_text <- function(x) {
    x <- tolower(trim_ws(x))
    out <- rep(NA_character_, length(x))
    pat6 <- "(^|[^0-9])6[[:space:]]*(h|hr|hrs|hour|hours)([^a-z]|$)|6_h|6_hr|6_hrs|6_hour|6_hours|6hr|6hrs|6hour|6hours|ne6|ctrl6|only6|early|t6"
    pat24 <- "(^|[^0-9])24[[:space:]]*(h|hr|hrs|hour|hours)([^a-z]|$)|24_h|24_hr|24_hrs|24_hour|24_hours|24hr|24hrs|24hour|24hours|ne24|ctrl24|only24|late|t24"
    out[grepl(pat6, x, perl = TRUE)] <- "6h"
    out[grepl(pat24, x, perl = TRUE)] <- "24h"
    out
  }
}
if (!exists("infer_class_from_text", mode = "function")) {
  infer_class_from_text <- function(x) {
    x <- tolower(trim_ws(x))
    out <- rep(NA_character_, length(x))
    out[grepl("early|only6|6h_only|commit", x)] <- "Early"
    out[grepl("sustained|shared|both|persistent", x)] <- "Sustained/Shared"
    out[grepl("late|only24|24h_only|maint", x)] <- "Late"
    out[grepl("reversal|switch|opposite", x)] <- "Reversal"
    out
  }
}
if (!exists("infer_direction_from_text", mode = "function")) {
  infer_direction_from_text <- function(x) {
    x <- tolower(trim_ws(x))
    out <- rep(NA_character_, length(x))
    out[grepl("up|upregulated|activated|induced|positive", x)] <- "UP"
    out[grepl("down|downregulated|suppressed|repressed|negative", x)] <- "DOWN"
    out
  }
}

pb <- utils::txtProgressBar(min = 0, max = length(slot_order), style = 3)
slot_cache <- list()
audit_rows <- list()

for (i in seq_along(slot_order)) {
  nm <- slot_order[[i]]
  log_message("Loading slot: ", nm)
  slot_cache[[nm]] <- load_slot_data(nm)
  audit_rows[[nm]] <- make_input_audit_row(nm, slot_cache[[nm]]$files)
  audit_slot_structure(nm, slot_cache[[nm]])
  utils::setTxtProgressBar(pb, i)
}
close(pb)

audit_tbl <- bind_rows(audit_rows)
write_table_if_any(audit_tbl, INPUT_AUDIT)

candidate_store <- list()
all_sections <- list()

fig2 <- report_fig2(slot_cache$fig2_mrna)
fig3 <- report_fig3(slot_cache$fig3_mrna_phase)
fig4 <- report_fig4(slot_cache$fig4_lnc_de, slot_cache$fig4_lnc_network)
fig5 <- report_fig5(slot_cache$fig5_kegg)
fig6 <- report_fig6(slot_cache$fig6_tf)
fig7 <- report_fig7(slot_cache$fig7_mirna_mrna)
fig8 <- report_fig8(slot_cache$fig8_cerna)
supp_circ_rep <- report_supp_circ(slot_cache$supp_circ)

section_list <- list(fig2 = fig2, fig3 = fig3, fig4 = fig4, fig5 = fig5, fig6 = fig6, fig7 = fig7, fig8 = fig8, supp_circ = supp_circ_rep)

for (nm in names(section_list)) {
  sec <- section_list[[nm]]
  all_sections[[nm]] <- sec$lines
  per_fig_file_txt <- file.path(FIG_REPORT_DIR, paste0(nm, "_report.txt"))
  per_fig_file_md <- file.path(FIG_REPORT_DIR, paste0(nm, "_report.md"))
  writeLines(sec$lines, per_fig_file_txt)
  writeLines(sec$lines, per_fig_file_md)
  if (!is.null(sec$candidates)) {
    for (cand_nm in names(sec$candidates)) {
      candidate_store[[cand_nm]] <- sec$candidates[[cand_nm]]
      write_table_if_any(sec$candidates[[cand_nm]], file.path(TABLE_OUT_DIR, paste0(cand_nm, ".tsv")))
    }
  }
}

cross <- report_cross_layer(candidate_store)
all_sections[["cross_layer"]] <- cross$lines
writeLines(cross$lines, file.path(FIG_REPORT_DIR, "cross_layer_report.txt"))
writeLines(cross$lines, file.path(FIG_REPORT_DIR, "cross_layer_report.md"))

master_lines <- c(
  "# MASTER RESULTS REPORT (v3)",
  paste0("Project root: ", PROJECT_ROOT),
  paste0("Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Storyline preserved in this report:",
  "Figure 1 = phenotypic commitment kinetics and irreversible hypertrophy",
  "Figure 2 = mRNA differential expression at 6 h and 24 h",
  "Figure 3 = Early / Sustained / Late / Reversal mRNA classes",
  "Figure 4 = lncRNA differential expression, response classes, and regulatory layer",
  "Figure 5 = KEGG pathway programs across commitment and maintenance",
  "Figure 6 = TF layer controlling commitment vs maintenance",
  "Figure 7 = miRNA-mRNA layer, mainly as maintenance reinforcement",
  "Figure 8 = integrative ceRNA layer only if coherent enough",
  ""
)

for (nm in names(all_sections)) {
  master_lines <- c(master_lines, all_sections[[nm]], "")
}

writeLines(master_lines, TXT_OUT)
writeLines(master_lines, MD_OUT)

log_message("Master report written: ", TXT_OUT)
log_message("Master markdown written: ", MD_OUT)
log_message("Per-figure reports: ", FIG_REPORT_DIR)
log_message("Candidate tables: ", TABLE_OUT_DIR)
log_message("Input audit: ", INPUT_AUDIT)

cat("\nDone. Outputs created in:\n", REPORT_DIR, "\n", sep = "")
