#!/usr/bin/env Rscript
# Stage 11: Figure 7 — ceRNA coherent lncRNA-miRNA pairs and lncRNA-miRNA-mRNA
#           triplets (lncRNA focus; circRNA branches excluded per project scope).
# Integrates IntaRNA lncRNA-miRNA predictions (energy-tiered) with miRNA DE
# direction and miRWalk targets + mRNA DE, then keeps only direction-COHERENT
# ceRNA relationships:
#   * coherent pair    : lncRNA and miRNA regulated in OPPOSITE directions
#   * coherent triplet : lncRNA UP / miRNA DOWN / mRNA UP  (or the mirror)
#
# Ported from 12_ceRNA_coherent_pairs_and_triplets_cleanFinal.R with all
# circRNA logic removed, Windows paths removed, and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--intarna-dir"), type = "character",
              help = "Directory with IntaRNA results (Interact_*.csv)"),
  make_option(c("--lnc-annotation"), type = "character", default = NULL,
              help = "lncRNA_DE_with_phase_and_transcriptID.csv from script 10 (lncRNA direction/phase per timepoint)"),
  make_option(c("--mirna-de-dir"), type = "character", default = NULL,
              help = "Directory with miRNA_DE_{6h,24h}_{UP,DOWN}.txt (from script 02)"),
  make_option(c("--mirwalk-dir"), type = "character", default = NULL,
              help = "Directory with miRWalk_miRNA_Targets_*.csv"),
  make_option(c("--mrna-de-dir"), type = "character", default = NULL,
              help = "Directory with NE_vs_Ctrl_{6h,24h}_mRNA_DE.txt + early/sustained/late_mRNA.csv"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/11_ceRNA",
              help = "Output directory [default: %default]"),
  make_option(c("--energy-high"), type = "double", default = -20, help = "High tier E cutoff [default: %default]"),
  make_option(c("--energy-moderate"), type = "double", default = -15, help = "Moderate tier E cutoff [default: %default]"),
  make_option(c("--energy-exploratory"), type = "double", default = -10, help = "Exploratory tier E cutoff [default: %default]"),
  make_option(c("--moderate-or-better"), action = "store_true", default = FALSE,
              help = "Keep only High/Moderate interaction tiers"),
  make_option(c("--dpi"), type = "integer", default = 600, help = "PNG DPI")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$`intarna-dir`)) {
  cat("[WARN] --intarna-dir not specified.\n[INFO] Exiting with placeholder status.\n"); quit(status = 0)
}
INTARNA_DIR <- opt$`intarna-dir`
OUTPUT_DIR <- opt$`output-dir`
energy_high <- opt$`energy-high`; energy_moderate <- opt$`energy-moderate`; energy_exploratory <- opt$`energy-exploratory`
DPI_PNG <- opt$dpi
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
cat("[INFO] IntaRNA directory:", INTARNA_DIR, "\n[INFO] Output directory:", OUTPUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

normalize_key <- function(x) x %>% as.character() %>% str_trim() %>% str_replace_all("\\s+", "") %>%
  str_replace_all("[\"'`]", "") %>% str_to_lower()
normalize_mirna <- function(x) x %>% as.character() %>% str_trim() %>% str_to_lower()

canonical_direction <- function(x, logfc = NULL) {
  y <- str_to_upper(str_trim(as.character(x))); n <- length(y)
  pos_fc <- rep(FALSE, n); neg_fc <- rep(FALSE, n)
  if (!is.null(logfc)) {
    logfc <- suppressWarnings(as.numeric(logfc))
    if (length(logfc) == 1 && n > 1) logfc <- rep(logfc, n)
    if (length(logfc) == n) { pos_fc <- !is.na(logfc) & logfc > 0; neg_fc <- !is.na(logfc) & logfc < 0 }
  }
  case_when(y %in% c("UP", "UPREGULATED", "POSITIVE") ~ "UP",
            y %in% c("DOWN", "DOWNREGULATED", "NEGATIVE") ~ "DOWN",
            pos_fc ~ "UP", neg_fc ~ "DOWN", TRUE ~ NA_character_)
}
energy_tier <- function(E) case_when(
  is.na(E) ~ "Unscored", E <= energy_high ~ "High",
  E > energy_high & E <= energy_moderate ~ "Moderate",
  E > energy_moderate & E <= energy_exploratory ~ "Exploratory", TRUE ~ "Weak")
rank_tier <- function(x) recode(x, High = 1L, Moderate = 2L, Exploratory = 3L, Weak = 4L, Unscored = 5L, .default = 99L)

read_delim_auto <- function(path) {
  fl <- readLines(path, n = 1, warn = FALSE)
  if (length(fl) == 0) return(read_csv(path, show_col_types = FALSE))
  if (str_detect(fl, ";") && !str_detect(fl, ",")) return(read_delim(path, delim = ";", show_col_types = FALSE, trim_ws = TRUE))
  if (str_detect(fl, "\\t")) return(read_tsv(path, show_col_types = FALSE))
  read_csv(path, show_col_types = FALSE)
}
choose_col <- function(df, cands, fallback = NULL) { hit <- cands[cands %in% names(df)][1]; if (is.na(hit)) fallback else hit }

# Parse a lncRNA IntaRNA id1: leading transcript_id token + optional "key=value"
# suffixes (matching script 10's FASTA header format) or pipe-positional fields.
parse_lnc_id1 <- function(x) {
  x <- as.character(x)
  tx <- str_trim(str_replace(x, "[ |].*$", ""))
  getkv <- function(key) { m <- str_match(x, paste0(key, "\\s*=\\s*([^|]+)")); str_trim(m[, 2]) }
  gene <- getkv("gene_id"); phase <- getkv("phase"); src <- getkv("source")
  dir <- getkv("direction"); tp <- getkv("timepoint")
  # Fallback to pipe-positional: transcript|gene|phase|source|direction|timepoint
  if (all(is.na(gene))) {
    parts <- str_split_fixed(x, "\\|", 6)
    gene <- na_if(str_trim(parts[, 2]), ""); phase <- na_if(str_trim(parts[, 3]), "")
    src <- na_if(str_trim(parts[, 4]), ""); dir <- na_if(str_trim(parts[, 5]), "")
    tp <- na_if(str_trim(parts[, 6]), "")
  }
  tibble(transcript_id_id1 = na_if(tx, ""), gene_id_id1 = gene, phase_id1 = phase,
         source_id1 = src, direction_id1 = str_to_upper(dir), timepoint_id1 = tp,
         tx_key_id1 = normalize_key(tx), gene_key_id1 = normalize_key(gene))
}

# ---------------------------------------------------------------------------
# 1. IntaRNA lncRNA-miRNA interactions
# ---------------------------------------------------------------------------

read_intarna_file <- function(path) {
  raw <- read_delim_auto(path)
  if (ncol(raw) == 1 && str_detect(names(raw)[1], "id1;")) {
    raw <- raw %>% rename(raw_col = 1) %>%
      separate_wider_delim(raw_col, delim = ";",
        names = c("id1", "start1", "end1", "id2", "start2", "end2", "subseqDP", "hybridDP", "E"),
        too_few = "align_start")
  }
  if (!all(c("id1", "id2", "E") %in% names(raw))) { cat("[WARN] Unrecognized IntaRNA file:", basename(path), "\n"); return(NULL) }
  fn <- basename(path)
  # lncRNA focus: skip files explicitly labeled circ
  if (str_detect(str_to_lower(fn), "circ")) return(NULL)
  raw %>% transmute(
    source_file = fn, id1 = as.character(id1),
    miRNA = as.character(id2), miRNA_key = normalize_mirna(id2),
    E = suppressWarnings(as.numeric(E)),
    interaction_tier = energy_tier(E),
    timepoint_file = case_when(str_detect(fn, regex("24h", TRUE)) ~ "24h", str_detect(fn, regex("6h", TRUE)) ~ "6h", TRUE ~ NA_character_),
    expected_ncrna_direction = case_when(str_detect(fn, regex("UP", TRUE)) ~ "UP", str_detect(fn, regex("DOWN", TRUE)) ~ "DOWN", TRUE ~ NA_character_)
  )
}

intarna_files <- list.files(INTARNA_DIR, pattern = "^Interact.*\\.csv$", full.names = TRUE, recursive = TRUE)
intarna_tbl <- bind_rows(lapply(intarna_files, read_intarna_file))
if (nrow(intarna_tbl) == 0) { cat("[WARN] No usable lncRNA IntaRNA interactions found.\n"); quit(status = 0) }
cat("[INFO] Loaded", nrow(intarna_tbl), "lncRNA-miRNA interactions from", length(intarna_files), "files\n")

id1_parsed <- parse_lnc_id1(intarna_tbl$id1)
intarna_tbl <- bind_cols(intarna_tbl, id1_parsed) %>%
  mutate(timepoint = coalesce(timepoint_id1, timepoint_file))

# ---------------------------------------------------------------------------
# 2. lncRNA annotation (direction + phase per timepoint)
# ---------------------------------------------------------------------------

lnc_ann <- NULL
if (!is.null(opt$`lnc-annotation`) && file.exists(opt$`lnc-annotation`)) {
  lnc_ann <- read_csv(opt$`lnc-annotation`, show_col_types = FALSE) %>%
    transmute(
      tx_key = normalize_key(if ("transcript_id" %in% names(.)) transcript_id else gene_id),
      gene_key = normalize_key(gene_id),
      gene_id, transcript_id = if ("transcript_id" %in% names(.)) transcript_id else NA_character_,
      gene_name = if ("gene_name" %in% names(.)) gene_name else gene_id,
      source = if ("source" %in% names(.)) source else "lncRNA",
      phase = if ("phase" %in% names(.)) phase else "Unassigned",
      timepoint, ncRNA_direction = canonical_direction(direction, if ("log2FoldChange" %in% names(.)) log2FoldChange else NULL)
    )
}

annotate_lnc <- function() {
  base <- intarna_tbl
  if (!is.null(lnc_ann)) {
    by_tx <- lnc_ann %>% distinct(tx_key, timepoint, .keep_all = TRUE)
    by_gene <- lnc_ann %>% distinct(gene_key, timepoint, .keep_all = TRUE)
    base <- base %>%
      left_join(by_tx %>% transmute(timepoint, tx_key, gene_id_t = gene_id, gene_name_t = gene_name,
                                    source_t = source, phase_t = phase, dir_t = ncRNA_direction),
                by = c("timepoint", "tx_key_id1" = "tx_key")) %>%
      left_join(by_gene %>% transmute(timepoint, gene_key, gene_id_g = gene_id, gene_name_g = gene_name,
                                      source_g = source, phase_g = phase, dir_g = ncRNA_direction),
                by = c("timepoint", "gene_key_id1" = "gene_key")) %>%
      mutate(gene_id = coalesce(gene_id_t, gene_id_g, gene_id_id1),
             gene_name = coalesce(gene_name_t, gene_name_g, gene_id_id1),
             source = coalesce(source_t, source_g, source_id1, "lncRNA"),
             phase = coalesce(phase_t, phase_g, phase_id1, "Unassigned"),
             ncRNA_direction = coalesce(dir_t, dir_g, canonical_direction(direction_id1), expected_ncrna_direction))
  } else {
    base <- base %>% mutate(gene_id = coalesce(gene_id_id1, transcript_id_id1),
                            gene_name = gene_id, source = coalesce(source_id1, "lncRNA"),
                            phase = coalesce(phase_id1, "Unassigned"),
                            ncRNA_direction = coalesce(canonical_direction(direction_id1), expected_ncrna_direction))
  }
  base %>% mutate(ncRNA_label = coalesce(transcript_id_id1, gene_id),
                  ncrna_key = normalize_key(coalesce(transcript_id_id1, gene_id)))
}
ann <- annotate_lnc()

# ---------------------------------------------------------------------------
# 3. miRNA DE direction
# ---------------------------------------------------------------------------

read_mirna_de <- function(dir) {
  if (is.null(dir)) return(NULL)
  specs <- expand.grid(tp = c("6h", "24h"), dir = c("UP", "DOWN"), stringsAsFactors = FALSE)
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    p <- file.path(dir, paste0("miRNA_DE_", specs$tp[i], "_", specs$dir[i], ".txt"))
    if (!file.exists(p)) return(NULL)
    x <- suppressWarnings(read_tsv(p, show_col_types = FALSE))
    mcol <- choose_col(x, c("miRNA_id", "miRNA", "mirna"), names(x)[1])
    tibble(miRNA_key = normalize_mirna(x[[mcol]]), timepoint = specs$tp[i], miRNA_direction = specs$dir[i],
           miRNA_log2FC = suppressWarnings(as.numeric(if ("log2FoldChange" %in% names(x)) x[["log2FoldChange"]] else NA)),
           miRNA_padj = suppressWarnings(as.numeric(if ("padj" %in% names(x)) x[["padj"]] else NA)))
  })) %>% distinct(miRNA_key, timepoint, .keep_all = TRUE)
}
miRNA_de <- read_mirna_de(opt$`mirna-de-dir`)

annot_pairs <- ann
if (!is.null(miRNA_de)) {
  annot_pairs <- annot_pairs %>% left_join(miRNA_de, by = c("miRNA_key", "timepoint"))
} else {
  annot_pairs <- annot_pairs %>% mutate(miRNA_direction = NA_character_, miRNA_log2FC = NA_real_, miRNA_padj = NA_real_)
}
annot_pairs <- annot_pairs %>%
  mutate(miRNA_direction = canonical_direction(coalesce(miRNA_direction)),
         coherent_pair = (ncRNA_direction == "UP" & miRNA_direction == "DOWN") |
                         (ncRNA_direction == "DOWN" & miRNA_direction == "UP"),
         coherent_pair = replace_na(coherent_pair, FALSE),
         phase = coalesce(phase, "Unassigned")) %>%
  distinct()

if (opt$`moderate-or-better`) annot_pairs <- annot_pairs %>% filter(interaction_tier %in% c("High", "Moderate"))

# Interaction state across timepoints (Early / Late / Sustained)
pair_presence <- annot_pairs %>% filter(coherent_pair) %>%
  distinct(ncrna_key, miRNA_key, timepoint) %>%
  count(ncrna_key, miRNA_key, timepoint) %>%
  pivot_wider(names_from = timepoint, values_from = n, values_fill = 0)
if (!"6h" %in% names(pair_presence)) pair_presence$`6h` <- 0
if (!"24h" %in% names(pair_presence)) pair_presence$`24h` <- 0
pair_presence <- pair_presence %>%
  mutate(interaction_state = case_when(`6h` > 0 & `24h` == 0 ~ "Early", `6h` == 0 & `24h` > 0 ~ "Late",
                                       `6h` > 0 & `24h` > 0 ~ "Sustained", TRUE ~ "Unassigned"))

coherent_pairs <- annot_pairs %>% filter(coherent_pair) %>%
  left_join(pair_presence %>% select(ncrna_key, miRNA_key, interaction_state), by = c("ncrna_key", "miRNA_key")) %>%
  mutate(interaction_state = coalesce(interaction_state, timepoint))
write_csv(coherent_pairs, file.path(OUTPUT_DIR, "COHERENT_lncRNA_miRNA_pairs.csv"))
cat("[INFO] Coherent lncRNA-miRNA pairs:", nrow(coherent_pairs), "\n")

# ---------------------------------------------------------------------------
# 4. miRWalk targets + mRNA DE -> coherent triplets
# ---------------------------------------------------------------------------

read_miRWalk <- function(dir) {
  if (is.null(dir)) return(NULL)
  files <- list.files(dir, pattern = "^miRWalk_miRNA_Targets_.*\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(NULL)
  map_dfr(files, function(p) {
    x <- read_delim_auto(p); fn <- basename(p)
    mcol <- choose_col(x, c("miRNA", "mirna", "miRNA_name"), names(x)[1])
    gcol <- choose_col(x, c("gene_symbol", "gene", "target_gene", "gene_name", "Symbol"), names(x)[min(2, ncol(x))])
    tibble(miRNA_key = normalize_mirna(x[[mcol]]), gene_key = normalize_key(x[[gcol]]),
           region = case_when(str_detect(fn, regex("3UTR", TRUE)) ~ "3UTR", str_detect(fn, regex("CDS", TRUE)) ~ "CDS",
                              str_detect(fn, regex("5UTR", TRUE)) ~ "5UTR", TRUE ~ "Unknown"),
           timepoint = case_when(str_detect(fn, regex("24", TRUE)) ~ "24h", str_detect(fn, regex("6", TRUE)) ~ "6h", TRUE ~ NA_character_))
  }) %>% distinct()
}
read_mRNA_map <- function(dir) {
  if (is.null(dir)) return(NULL)
  phase_map <- bind_rows(
    { p <- file.path(dir, "early_mRNA.csv"); if (file.exists(p)) read_csv(p, show_col_types = FALSE) %>% mutate(phase = "Early") else NULL },
    { p <- file.path(dir, "sustained_mRNA.csv"); if (file.exists(p)) read_csv(p, show_col_types = FALSE) %>% mutate(phase = "Sustained") else NULL },
    { p <- file.path(dir, "late_mRNA.csv"); if (file.exists(p)) read_csv(p, show_col_types = FALSE) %>% mutate(phase = "Late") else NULL }
  )
  phase_map <- if (!is.null(phase_map) && nrow(phase_map) > 0) {
    phase_map %>% mutate(gene_key = normalize_key(if ("gene_name" %in% names(.)) gene_name else gene_id)) %>%
      transmute(gene_key, phase) %>% distinct(gene_key, .keep_all = TRUE)
  } else tibble(gene_key = character(), phase = character())
  read_de <- function(fname, tp) {
    p <- file.path(dir, fname); if (!file.exists(p)) return(NULL)
    x <- read_tsv(p, show_col_types = FALSE)
    scol <- choose_col(x, c("gene_name", "gene_symbol", "gene_id"), names(x)[1])
    tibble(gene_key = normalize_key(x[[scol]]), timepoint = tp,
           mRNA_log2FC = suppressWarnings(as.numeric(if ("log2FoldChange" %in% names(x)) x[["log2FoldChange"]] else NA)),
           mRNA_direction = canonical_direction(if ("direction" %in% names(x)) x[["direction"]] else NA,
                                                if ("log2FoldChange" %in% names(x)) x[["log2FoldChange"]] else NULL),
           mRNA_padj = suppressWarnings(as.numeric(if ("padj" %in% names(x)) x[["padj"]] else NA)))
  }
  bind_rows(read_de("NE_vs_Ctrl_6h_mRNA_DE.txt", "6h"), read_de("NE_vs_Ctrl_24h_mRNA_DE.txt", "24h")) %>%
    left_join(phase_map, by = "gene_key") %>% mutate(phase_mRNA = coalesce(phase, "Unassigned")) %>%
    distinct(gene_key, timepoint, .keep_all = TRUE)
}

miRWalk_targets <- read_miRWalk(opt$`mirwalk-dir`)
mRNA_map <- read_mRNA_map(opt$`mrna-de-dir`)

if (!is.null(miRWalk_targets) && !is.null(mRNA_map) && nrow(coherent_pairs) > 0) {
  triplets <- coherent_pairs %>%
    select(source_file, gene_id, transcript_id = transcript_id_id1, ncRNA_label, phase, timepoint,
           ncRNA_direction, miRNA, miRNA_key, miRNA_direction, E, interaction_tier, interaction_state) %>%
    inner_join(miRWalk_targets %>% select(miRNA_key, gene_key, region, timepoint),
               by = c("miRNA_key", "timepoint"), relationship = "many-to-many") %>%
    inner_join(mRNA_map %>% select(gene_key, timepoint, phase_mRNA, mRNA_direction, mRNA_log2FC, mRNA_padj),
               by = c("gene_key", "timepoint")) %>%
    mutate(coherent_triplet = (ncRNA_direction == "UP" & miRNA_direction == "DOWN" & mRNA_direction == "UP") |
                              (ncRNA_direction == "DOWN" & miRNA_direction == "UP" & mRNA_direction == "DOWN"),
           coherent_triplet = replace_na(coherent_triplet, FALSE),
           phase_relation = paste0(phase, "_lncRNA__", phase_mRNA, "_mRNA")) %>%
    filter(coherent_triplet)
  write_csv(triplets, file.path(OUTPUT_DIR, "COHERENT_lncRNA_miRNA_mRNA_triplets.csv"))
  cat("[INFO] Coherent lncRNA-miRNA-mRNA triplets:", nrow(triplets), "\n")

  priority_triplets <- triplets %>% mutate(tier_rank = rank_tier(interaction_tier)) %>%
    arrange(tier_rank, E, desc(abs(mRNA_log2FC))) %>%
    group_by(ncRNA_label, miRNA, gene_key) %>% slice(1) %>% ungroup() %>% arrange(tier_rank, E)
  write_csv(priority_triplets, file.path(OUTPUT_DIR, "PRIORITY_ceRNA_triplets_best_per_lncRNA_miRNA_mRNA.csv"))

  summary_triplets <- triplets %>%
    count(timepoint, phase, phase_mRNA, phase_relation, interaction_tier, interaction_state, name = "n_triplets") %>%
    arrange(timepoint, phase, phase_mRNA, interaction_tier)
  write_csv(summary_triplets, file.path(OUTPUT_DIR, "SUMMARY_coherent_triplets_by_phase.csv"))

  if (nrow(summary_triplets) > 0) {
    p2 <- ggplot(summary_triplets, aes(x = phase_relation, y = n_triplets, fill = interaction_tier)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) + facet_wrap(~ timepoint, scales = "free_x") +
      labs(x = "lncRNA-mRNA phase relation", y = "Coherent ceRNA triplets", fill = "Tier") +
      theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUTPUT_DIR, "FIG_ceRNA_triplets_by_phase_relation_and_tier.png"), p2, width = 11, height = 6.5, dpi = DPI_PNG)
    ggsave(file.path(OUTPUT_DIR, "FIG_ceRNA_triplets_by_phase_relation_and_tier.pdf"), p2, width = 11, height = 6.5)
  }
} else {
  cat("[WARN] miRWalk targets or mRNA DE not provided; wrote pairs only (no triplets)\n")
}

# ---------------------------------------------------------------------------
# 5. Pair summaries + figure
# ---------------------------------------------------------------------------

summary_pairs <- coherent_pairs %>%
  count(timepoint, phase, ncRNA_direction, miRNA_direction, interaction_tier, interaction_state, name = "n_pairs") %>%
  arrange(timepoint, phase, ncRNA_direction, interaction_tier)
write_csv(summary_pairs, file.path(OUTPUT_DIR, "SUMMARY_coherent_lncRNA_miRNA_pairs.csv"))

if (nrow(summary_pairs) > 0) {
  p1 <- ggplot(summary_pairs, aes(x = phase, y = n_pairs, fill = interaction_tier)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) + facet_wrap(~ timepoint) +
    labs(x = "lncRNA phase", y = "Coherent lncRNA-miRNA pairs", fill = "Tier") + theme_bw(base_size = 11)
  ggsave(file.path(OUTPUT_DIR, "FIG_ceRNA_pairs_by_phase_and_tier.png"), p1, width = 9, height = 5.5, dpi = DPI_PNG)
  ggsave(file.path(OUTPUT_DIR, "FIG_ceRNA_pairs_by_phase_and_tier.pdf"), p1, width = 9, height = 5.5)
}

cat("\n[SUCCESS] ceRNA coherent pairs/triplets analysis complete.\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
