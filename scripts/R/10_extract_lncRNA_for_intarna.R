#!/usr/bin/env Rscript
# Stage 10: Extract DE lncRNA candidate transcripts for IntaRNA
# Builds a gene_id -> transcript_id map from the lncRNA phase/source class
# files, intersects with the lncRNA DE tables, subsets the transcriptome FASTA
# by transcript_id, and exports candidate tables + FASTA files broken down by
# timepoint / direction / phase, plus summary and priority-ranking tables.
#
# Ported from 09_extract_lncRNA_candidates_for_IntaRNA.R with Windows paths
# removed and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

option_list <- list(
  make_option(c("-l", "--lnc-dir"), type = "character", default = NULL,
              help = "Directory with lncRNA class files (early/sustained/late [novel]_lncRNA.csv)"),
  make_option(c("-d", "--de-dir"), type = "character", default = NULL,
              help = "Directory with lncRNA DE tables (NE_vs_Ctrl_{6h,24h}_lncRNA_DE.txt)"),
  make_option(c("-f", "--fasta-file"), type = "character", default = NULL,
              help = "Transcriptome FASTA containing lncRNA transcripts"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/10_intarna_candidates",
              help = "Output directory [default: %default]"),
  make_option(c("--alpha"), type = "numeric", default = 0.05,
              help = "Adjusted p-value threshold [default: %default]"),
  make_option(c("--min-abs-log2fc"), type = "numeric", default = 0,
              help = "abs(log2FC) threshold (original pipeline used 0) [default: %default]"),
  make_option(c("--export-by-phase"), action = "store_true", default = TRUE,
              help = "Also export candidate FASTA split by phase [default: %default]"),
  make_option(c("--priority-top-n"), type = "integer", default = 25,
              help = "Top lncRNAs per group in the priority table [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$`lnc-dir`) || is.null(opt$`de-dir`) || is.null(opt$`fasta-file`)) {
  print_help(parser)
  stop("\n[ERROR] Required arguments: --lnc-dir, --de-dir, --fasta-file")
}

LNC_DIR <- opt$`lnc-dir`
DE_DIR <- opt$`de-dir`
FASTA_FILE <- opt$`fasta-file`
OUTPUT_DIR <- opt$`output-dir`
PADJ_CUTOFF <- opt$alpha
MIN_ABS_LFC <- opt$`min-abs-log2fc`
EXPORT_BY_PHASE <- opt$`export-by-phase`
PRIORITY_N <- opt$`priority-top-n`
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!file.exists(FASTA_FILE)) stop("[ERROR] FASTA file not found: ", FASTA_FILE)

cat("[INFO] lncRNA class dir:", LNC_DIR, "\n")
cat("[INFO] DE results dir:", DE_DIR, "\n")
cat("[INFO] Transcriptome FASTA:", FASTA_FILE, "\n")
cat("[INFO] Output dir:", OUTPUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

read_fasta_simple <- function(fasta_path) {
  x <- readLines(fasta_path)
  hdr_idx <- grep("^>", x)
  if (length(hdr_idx) == 0) stop("[ERROR] No FASTA headers in: ", fasta_path)
  end_idx <- c(hdr_idx[-1] - 1, length(x))
  headers <- x[hdr_idx]
  seqs <- vapply(seq_along(hdr_idx), function(i)
    paste0(x[(hdr_idx[i] + 1):end_idx[i]], collapse = ""), character(1))
  tibble(header = headers, fasta_id = sub("^>(\\S+).*$", "\\1", headers), sequence = seqs)
}

write_fasta_simple <- function(df, out_file, id_col = "transcript_id", header_suffix_cols = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  lines <- character(nrow(df) * 2)
  for (i in seq_len(nrow(df))) {
    extra <- ""
    if (!is.null(header_suffix_cols)) {
      vals <- paste(paste0(header_suffix_cols, "=", as.character(df[i, header_suffix_cols, drop = TRUE])), collapse = " | ")
      extra <- paste0(" | ", vals)
    }
    lines[(2 * i) - 1] <- paste0(">", df[[id_col]][i], extra)
    lines[(2 * i)] <- df$sequence[i]
  }
  writeLines(lines, out_file)
}

# ---------------------------------------------------------------------------
# 1. gene_id -> transcript_id map from lncRNA class files (annotated + novel)
# ---------------------------------------------------------------------------

class_specs <- tribble(
  ~file,                          ~phase,       ~source,
  "early_lncRNA.csv",             "Early",      "annotated",
  "sustained_lncRNA.csv",         "Sustained",  "annotated",
  "late_lncRNA.csv",              "Late",       "annotated",
  "early_novel_lncRNA.csv",       "Early",      "novel",
  "sustained_novel_lncRNA.csv",   "Sustained",  "novel",
  "late_novel_lncRNA.csv",        "Late",       "novel"
)

lnc_map <- bind_rows(lapply(seq_len(nrow(class_specs)), function(i) {
  p <- file.path(LNC_DIR, class_specs$file[i])
  if (!file.exists(p)) return(NULL)
  df <- read_csv(p, show_col_types = FALSE)
  keep <- intersect(c("gene_id", "transcript_id", "gene_name", "gene_biotype", "transcript_biotype"), names(df))
  df %>% select(all_of(keep)) %>%
    mutate(phase = class_specs$phase[i], source = class_specs$source[i])
}))

if (is.null(lnc_map) || nrow(lnc_map) == 0)
  stop("[ERROR] No lncRNA class files found in ", LNC_DIR)

if (!"transcript_id" %in% names(lnc_map)) lnc_map$transcript_id <- NA_character_
if (!"gene_name" %in% names(lnc_map)) lnc_map$gene_name <- lnc_map$gene_id
lnc_map <- lnc_map %>% distinct(gene_id, transcript_id, .keep_all = TRUE)
write_csv(lnc_map, file.path(OUTPUT_DIR, "lncRNA_gene_to_transcript_map_from_classes.csv"))

map_check <- lnc_map %>% count(gene_id, name = "n_tx")
if (any(map_check$n_tx > 1)) cat("[WARN] Some gene_id map to multiple transcript_id; exporting all mapped transcripts.\n")

# ---------------------------------------------------------------------------
# 2. lncRNA DE tables -> annotate with phase/source, filter DE candidates
# ---------------------------------------------------------------------------

read_de <- function(fname, tp) {
  p <- file.path(DE_DIR, fname)
  if (!file.exists(p)) return(NULL)
  read_tsv(p, show_col_types = FALSE) %>% mutate(timepoint = tp)
}
de_all <- bind_rows(read_de("NE_vs_Ctrl_6h_lncRNA_DE.txt", "6h"),
                    read_de("NE_vs_Ctrl_24h_lncRNA_DE.txt", "24h"))
if (is.null(de_all)) stop("[ERROR] lncRNA DE tables not found in ", DE_DIR)

lnc_de <- de_all %>%
  filter(!is.na(padj), padj <= PADJ_CUTOFF, abs(log2FoldChange) >= MIN_ABS_LFC) %>%
  left_join(lnc_map, by = "gene_id", suffix = c("_de", "_map")) %>%
  mutate(
    source = coalesce(if ("source" %in% names(.)) source else NA_character_, "unknown"),
    gene_name = coalesce(if ("gene_name_map" %in% names(.)) gene_name_map else NULL,
                         if ("gene_name_de" %in% names(.)) gene_name_de else NULL, gene_id),
    direction = case_when(
      !is.na(direction) & str_to_upper(direction) %in% c("UP", "DOWN") ~ str_to_upper(direction),
      log2FoldChange > 0 ~ "UP", log2FoldChange < 0 ~ "DOWN", TRUE ~ "NS")
  ) %>%
  filter(direction %in% c("UP", "DOWN")) %>%
  arrange(timepoint, direction, phase, desc(abs(log2FoldChange)))

sel_cols <- intersect(c("gene_id", "transcript_id", "gene_name", "gene_biotype", "transcript_biotype",
                        "source", "phase", "timepoint", "baseMean", "log2FoldChange", "lfcSE",
                        "stat", "pvalue", "padj", "direction"), names(lnc_de))
lnc_de <- lnc_de %>% select(all_of(sel_cols))
write_csv(lnc_de, file.path(OUTPUT_DIR, "lncRNA_DE_with_phase_and_transcriptID.csv"))

n_missing_tx <- sum(is.na(lnc_de$transcript_id))
if (n_missing_tx > 0) cat("[WARN]", n_missing_tx, "DE lncRNAs had no transcript_id mapping.\n")

# ---------------------------------------------------------------------------
# 3. Match transcripts against FASTA (first token, then header fallback)
# ---------------------------------------------------------------------------

fasta_tbl <- read_fasta_simple(FASTA_FILE)
matched <- lnc_de %>% filter(!is.na(transcript_id)) %>%
  left_join(fasta_tbl, by = c("transcript_id" = "fasta_id"))

still_missing <- matched %>% filter(is.na(sequence)) %>% distinct(transcript_id)
if (nrow(still_missing) > 0) {
  fallback <- bind_rows(lapply(still_missing$transcript_id, function(tx) {
    rx <- paste0("(^|[|;: =])", str_replace_all(tx, "([.^$|()\\[\\]{}*+?\\\\-])", "\\\\\\1"), "($|[|;: =])")
    idx <- which(str_detect(fasta_tbl$header, rx))
    if (length(idx) == 0) return(NULL)
    fasta_tbl[idx[1], ] %>% mutate(transcript_id = tx)
  }))
  if (nrow(fallback) > 0) {
    matched <- matched %>% select(-sequence, -header) %>%
      left_join(fallback %>% select(transcript_id, header, sequence), by = "transcript_id")
  }
}
matched_fasta <- matched %>% filter(!is.na(sequence))
write_csv(matched_fasta %>% select(-sequence), file.path(OUTPUT_DIR, "lncRNA_DE_matched_to_fasta.csv"))

# ---------------------------------------------------------------------------
# 4. Summary tables
# ---------------------------------------------------------------------------

write_csv(lnc_de %>% count(timepoint, direction, phase, source, name = "n_lncRNAs") %>%
            arrange(timepoint, direction, phase, source),
          file.path(OUTPUT_DIR, "SUMMARY_lncRNA_candidates_by_time_direction_phase_source.csv"))
write_csv(matched_fasta %>% count(timepoint, direction, phase, source, name = "n_transcripts_with_fasta") %>%
            arrange(timepoint, direction, phase, source),
          file.path(OUTPUT_DIR, "SUMMARY_lncRNA_transcripts_with_fasta.csv"))

# ---------------------------------------------------------------------------
# 5. Candidate FASTA export by time / direction / phase
# ---------------------------------------------------------------------------

suffix_cols <- intersect(c("gene_id", "phase", "source", "direction", "timepoint"), names(matched_fasta))
for (tp in c("6h", "24h")) {
  tp_tbl <- matched_fasta %>% filter(timepoint == tp)
  if (nrow(tp_tbl) == 0) next
  write_csv(tp_tbl %>% select(-sequence, -header),
            file.path(OUTPUT_DIR, paste0("Candidates_lncRNA_", tp, "_ALL_with_tx_phase.csv")))
  write_fasta_simple(tp_tbl %>% distinct(transcript_id, .keep_all = TRUE),
                     file.path(OUTPUT_DIR, paste0("Candidates_lncRNA_", tp, "_ALL.fa")),
                     "transcript_id", suffix_cols)
  for (dirn in c("UP", "DOWN")) {
    sub_tbl <- tp_tbl %>% filter(direction == dirn)
    if (nrow(sub_tbl) == 0) next
    write_csv(sub_tbl %>% select(-sequence, -header),
              file.path(OUTPUT_DIR, paste0("Candidates_lncRNA_", tp, "_", dirn, "_with_tx_phase.csv")))
    write_fasta_simple(sub_tbl %>% distinct(transcript_id, .keep_all = TRUE),
                       file.path(OUTPUT_DIR, paste0("Candidates_lncRNA_", tp, "_", dirn, ".fa")),
                       "transcript_id", suffix_cols)
    if (EXPORT_BY_PHASE) {
      for (ph in c("Early", "Sustained", "Late")) {
        ph_tbl <- sub_tbl %>% filter(phase == ph)
        if (nrow(ph_tbl) == 0) next
        write_fasta_simple(ph_tbl %>% distinct(transcript_id, .keep_all = TRUE),
                           file.path(OUTPUT_DIR, paste0("Candidates_lncRNA_", tp, "_", dirn, "_", ph, ".fa")),
                           "transcript_id", suffix_cols)
      }
    }
  }
}

# ---------------------------------------------------------------------------
# 6. Priority ranking table for manual prioritization
# ---------------------------------------------------------------------------

priority_tbl <- lnc_de %>%
  mutate(abs_log2FC = abs(log2FoldChange)) %>%
  group_by(timepoint, direction, phase, source) %>%
  slice_max(order_by = abs_log2FC, n = PRIORITY_N, with_ties = FALSE) %>%
  ungroup() %>%
  select(any_of(c("timepoint", "direction", "phase", "source", "gene_id", "transcript_id",
                  "log2FoldChange", "padj", "gene_name")))
write_csv(priority_tbl, file.path(OUTPUT_DIR, "PRIORITY_top_lncRNAs_per_group_for_IntaRNA.csv"))

cat("\n[SUCCESS] lncRNA candidate extraction complete!\n")
cat("[OUTPUT] FASTA + tables in:", OUTPUT_DIR, "\n")
