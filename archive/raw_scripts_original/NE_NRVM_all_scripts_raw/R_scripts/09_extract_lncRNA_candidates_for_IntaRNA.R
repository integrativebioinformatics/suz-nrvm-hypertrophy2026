#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Extract DE lncRNA candidate transcripts for IntaRNA
# Purpose:
# 1) Use existing class files (Early / Sustained / Late; annotated / novel)
#    as the intermediate mapping gene_id -> transcript_id.
# 2) Intersect with DE lncRNA tables at 6h and 24h.
# 3) Export candidate tables and FASTA files by timepoint/direction,
#    with phase/source annotations for downstream lncRNA-miRNA analysis.
#
# Notes:
# - This script assumes the class files already exist from the backbone.
# - In your current outputs, gene_id <-> transcript_id is effectively 1:1,
#   so an extra GTF step is not required unless your local files differ.
# - lncRNAs are labeled with gene_id in the outputs; transcript_id is used
#   to subset FASTA sequences.

# User parameters
base_dir   <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"
out_dir  <- file.path(base_dir, "IntaRNA_lncRNA_candidates")

# Edit this path to your transcript FASTA containing the lncRNA transcripts.
# Examples:
# fasta_file <- file.path(base_dir, "merged_lncRNAs.fa")
# fasta_file <- file.path(base_dir, "Transcriptome_Reconstruction", "all_transcripts.fa")
fasta_file <- file.path(base_dir, "merged_annotated_rn8_CNE624_8_H.fa")

padj_cutoff      <- 0.05
min_abs_log2fc   <- 0
export_by_phase  <- TRUE

# Optional: if TRUE, will also write a single combined FASTA per timepoint
# pooling UP and DOWN together.
export_combined_per_time <- TRUE

# Helper functions
find_first_existing <- function(filename, search_dirs) {
  # Direct path check first
  for (d in search_dirs) {
    p <- file.path(d, filename)
    if (file.exists(p)) return(normalizePath(p))
  }
  # Recursive search as fallback
  for (d in search_dirs) {
    hits <- list.files(d, pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
                       recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) return(normalizePath(hits[1]))
  }
  stop("Could not find file: ", filename)
}

read_fasta_simple <- function(fasta_path) {
  x <- readLines(fasta_path)
  hdr_idx <- grep("^>", x)
  if (length(hdr_idx) == 0) stop("No FASTA headers found in: ", fasta_path)
  end_idx <- c(hdr_idx[-1] - 1, length(x))
  seqs <- vector("list", length(hdr_idx))
  headers <- x[hdr_idx]
  for (i in seq_along(hdr_idx)) {
    seqs[[i]] <- paste0(x[(hdr_idx[i] + 1):end_idx[i]], collapse = "")
  }
  ids <- sub("^>(\\S+).*$", "\\1", headers)
  tibble(header = headers, fasta_id = ids, sequence = unlist(seqs))
}

write_fasta_simple <- function(df, out_file, id_col = "transcript_id", header_suffix_cols = NULL) {
  if (nrow(df) == 0) {
    warning("No sequences to write for: ", out_file)
    return(invisible(NULL))
  }
  lines <- character(nrow(df) * 2)
  for (i in seq_len(nrow(df))) {
    extra <- ""
    if (!is.null(header_suffix_cols)) {
      vals <- paste(paste0(header_suffix_cols, "=", as.character(df[i, header_suffix_cols, drop = TRUE])), collapse = " | ")
      extra <- paste0(" | ", vals)
    }
    lines[(2*i)-1] <- paste0(">", df[[id_col]][i], extra)
    lines[(2*i)]   <- df$sequence[i]
  }
  writeLines(lines, out_file)
}

message_if_missing <- function(df, label) {
  miss <- sum(is.na(df$transcript_id))
  if (miss > 0) {
    message(label, ": ", miss, " genes had no transcript_id mapping and were excluded from FASTA export.")
  }
}

# Locate required files
search_dirs <- unique(c(
  base_dir,
  file.path(base_dir, "DESeq2_Multifactorial_Results_fixed"),
  file.path(base_dir, "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT"),
  file.path(base_dir, "Paper_Fig4_REFINED")
))
search_dirs <- search_dirs[dir.exists(search_dirs)]

required_class_files <- c(
  "early_lncRNA.csv",
  "sustained_lncRNA.csv",
  "late_lncRNA.csv",
  "early_novel_lncRNA.csv",
  "sustained_novel_lncRNA.csv",
  "late_novel_lncRNA.csv"
)

class_paths <- setNames(lapply(required_class_files, find_first_existing, search_dirs = search_dirs), required_class_files)
de6_path    <- find_first_existing("NE_vs_Ctrl_6h_lncRNA_DE.txt",  search_dirs)
de24_path   <- find_first_existing("NE_vs_Ctrl_24h_lncRNA_DE.txt", search_dirs)

if (!file.exists(fasta_file)) {
  stop("FASTA file not found. Edit 'fasta_file' at the top of the script. Current value: ", fasta_file)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Build gene_id -> transcript_id mapping
phase_from_file <- c(
  early_lncRNA.csv = "Early",
  sustained_lncRNA.csv = "Sustained",
  late_lncRNA.csv = "Late",
  early_novel_lncRNA.csv = "Early",
  sustained_novel_lncRNA.csv = "Sustained",
  late_novel_lncRNA.csv = "Late"
)

source_from_file <- c(
  early_lncRNA.csv = "annotated",
  sustained_lncRNA.csv = "annotated",
  late_lncRNA.csv = "annotated",
  early_novel_lncRNA.csv = "novel",
  sustained_novel_lncRNA.csv = "novel",
  late_novel_lncRNA.csv = "novel"
)

lnc_map <- bind_rows(lapply(names(class_paths), function(fn) {
  read_csv(class_paths[[fn]], show_col_types = FALSE) %>%
    mutate(phase = phase_from_file[[fn]],
           source = source_from_file[[fn]]) %>%
    select(gene_id, transcript_id, gene_name, gene_biotype, transcript_biotype, source, phase)
})) %>%
  distinct(gene_id, transcript_id, .keep_all = TRUE)

# Sanity check: in your current outputs this should be mostly 1:1
map_check <- lnc_map %>% count(gene_id, name = "n_tx")
if (any(map_check$n_tx > 1)) {
  warning("Some gene_id values map to multiple transcript_id values. This script will export all mapped transcripts.")
}

write_csv(lnc_map, file.path(out_dir, "lncRNA_gene_to_transcript_map_from_classes.csv"))

# Read DE tables and annotate
de6 <- read_tsv(de6_path, show_col_types = FALSE) %>% mutate(timepoint = "6h")
de24 <- read_tsv(de24_path, show_col_types = FALSE) %>% mutate(timepoint = "24h")

required_de_cols <- c("gene_id", "log2FoldChange", "padj", "direction")
for (nm in required_de_cols) {
  if (!(nm %in% colnames(de6)) || !(nm %in% colnames(de24))) {
    stop("Missing required DE column: ", nm)
  }
}

lnc_de <- bind_rows(de6, de24) %>%
  filter(!is.na(padj), padj <= padj_cutoff, abs(log2FoldChange) >= min_abs_log2fc) %>%
  left_join(lnc_map, by = "gene_id", suffix = c("_de", "_map")) %>%
  mutate(
    source = coalesce(source_map, source_de),
    gene_name = coalesce(gene_name_map, gene_name_de),
    gene_biotype = coalesce(gene_biotype_map, gene_biotype_de),
    direction = case_when(
      str_to_upper(direction) %in% c("UP", "DOWN") ~ str_to_upper(direction),
      log2FoldChange > 0 ~ "UP",
      log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    )
  ) %>%
  select(
    gene_id, transcript_id, gene_name, gene_biotype, transcript_biotype,
    source, phase, timepoint, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, direction
  ) %>%
  filter(direction %in% c("UP", "DOWN")) %>%
  arrange(timepoint, direction, phase, desc(abs(log2FoldChange)))

message_if_missing(lnc_de, "DE lncRNA table")

write_csv(lnc_de, file.path(out_dir, "lncRNA_DE_with_phase_and_transcriptID.csv"))

colnames(lnc_de)

# Read transcript FASTA and match transcript IDs
fasta_tbl <- read_fasta_simple(fasta_file)

# First try exact first-token matches
matched <- lnc_de %>%
  filter(!is.na(transcript_id)) %>%
  left_join(fasta_tbl, by = c("transcript_id" = "fasta_id"))

# Fallback: if some transcript IDs are not matched by first token, search full header
still_missing <- matched %>% filter(is.na(sequence)) %>% distinct(transcript_id)
if (nrow(still_missing) > 0) {
  message("Fallback header matching for ", nrow(still_missing), " transcript IDs not found as first FASTA token.")
  fallback_hits <- lapply(still_missing$transcript_id, function(tx) {
    rx <- paste0("(^|[|;: =])", stringr::str_replace_all(tx, "([.^$|()\\[\\]{}*+?\\\\-])", "\\\\\\1"), "($|[|;: =])")
    idx <- which(str_detect(fasta_tbl$header, rx))
    if (length(idx) == 0) return(NULL)
    fasta_tbl[idx[1], ] %>% mutate(transcript_id = tx)
  })
  fallback_hits <- bind_rows(fallback_hits)
  if (nrow(fallback_hits) > 0) {
    matched <- matched %>%
      select(-sequence, -header) %>%
      left_join(fallback_hits %>% select(transcript_id, header, sequence), by = "transcript_id", suffix = c("", ""))
  }
}

# Keep only rows with sequences for FASTA export
matched_fasta <- matched %>% filter(!is.na(sequence))
write_csv(matched_fasta %>% select(-sequence), file.path(out_dir, "lncRNA_DE_with_phase_and_transcriptID_matched_to_fasta.csv"))

# Export summary tables
summary_counts <- lnc_de %>%
  count(timepoint, direction, phase, source, name = "n_lncRNAs") %>%
  arrange(timepoint, direction, phase, source)
write_csv(summary_counts, file.path(out_dir, "SUMMARY_lncRNA_candidates_by_time_direction_phase_source.csv"))

summary_tx_counts <- matched_fasta %>%
  count(timepoint, direction, phase, source, name = "n_transcripts_with_fasta") %>%
  arrange(timepoint, direction, phase, source)
write_csv(summary_tx_counts, file.path(out_dir, "SUMMARY_lncRNA_transcripts_with_fasta_by_time_direction_phase_source.csv"))

# Export candidate tables and FASTA files
for (tp in c("6h", "24h")) {
  tp_tbl <- matched_fasta %>% filter(timepoint == tp)

  if (nrow(tp_tbl) == 0) next

  if (export_combined_per_time) {
    write_csv(tp_tbl %>% select(-sequence, -header),
              file.path(out_dir, paste0("Candidates_lncRNA_", tp, "_ALL_with_tx_phase.csv")))
    write_fasta_simple(
      tp_tbl %>% distinct(transcript_id, .keep_all = TRUE),
      out_file = file.path(out_dir, paste0("Candidates_lncRNA_", tp, "_ALL.fa")),
      id_col = "transcript_id",
      header_suffix_cols = c("gene_id", "phase", "source", "direction", "timepoint")
    )
  }

  for (dirn in c("UP", "DOWN")) {
    sub_tbl <- tp_tbl %>% filter(direction == dirn)
    if (nrow(sub_tbl) == 0) next

    write_csv(sub_tbl %>% select(-sequence, -header),
              file.path(out_dir, paste0("Candidates_lncRNA_", tp, "_", dirn, "_with_tx_phase.csv")))
    write_fasta_simple(
      sub_tbl %>% distinct(transcript_id, .keep_all = TRUE),
      out_file = file.path(out_dir, paste0("Candidates_lncRNA_", tp, "_", dirn, ".fa")),
      id_col = "transcript_id",
      header_suffix_cols = c("gene_id", "phase", "source", "direction", "timepoint")
    )

    if (export_by_phase) {
      for (ph in c("Early", "Sustained", "Late")) {
        ph_tbl <- sub_tbl %>% filter(phase == ph)
        if (nrow(ph_tbl) == 0) next

        write_csv(ph_tbl %>% select(-sequence, -header),
                  file.path(out_dir, paste0("Candidates_lncRNA_", tp, "_", dirn, "_", ph, "_with_tx_phase.csv")))
        write_fasta_simple(
          ph_tbl %>% distinct(transcript_id, .keep_all = TRUE),
          out_file = file.path(out_dir, paste0("Candidates_lncRNA_", tp, "_", dirn, "_", ph, ".fa")),
          id_col = "transcript_id",
          header_suffix_cols = c("gene_id", "phase", "source", "direction", "timepoint")
        )
      }
    }
  }
}

# Export a compact candidate ranking for manual prioritization
priority_tbl <- lnc_de %>%
  mutate(abs_log2FC = abs(log2FoldChange)) %>%
  group_by(timepoint, direction, phase, source) %>%
  slice_max(order_by = abs_log2FC, n = 25, with_ties = FALSE) %>%
  ungroup() %>%
  select(timepoint, direction, phase, source, gene_id, transcript_id, log2FoldChange, padj, gene_name)
write_csv(priority_tbl, file.path(out_dir, "PRIORITY_top25_lncRNAs_per_group_for_IntaRNA.csv"))

message("Done. Outputs written to: ", normalizePath(out_dir))
message("Key files:")
message(" - lncRNA_gene_to_transcript_map_from_classes.csv")
message(" - lncRNA_DE_with_phase_and_transcriptID.csv")
message(" - SUMMARY_lncRNA_candidates_by_time_direction_phase_source.csv")
message(" - Candidates_lncRNA_<time>_<direction>.fa")
message(" - Candidates_lncRNA_<time>_<direction>_<phase>.fa")
