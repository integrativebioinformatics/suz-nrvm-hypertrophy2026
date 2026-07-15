#!/usr/bin/env Rscript
# Stage 09: Figure 6 — miRNA-mRNA reciprocal integration mapped to KEGG
# For each timepoint (6h/24h) and reciprocal combo (miRNA UP -> mRNA DOWN and
# miRNA DOWN -> mRNA UP), intersects DE miRNAs with miRWalk targets (3'UTR + CDS
# + 5'UTR), keeps mRNA targets DE in the opposite direction, annotates with
# temporal phase and KEGG membership/sign, and emits the SUPP_edges_* tables
# and dotplots consumed by scripts 12/13.
#
# Ported from script_miRNA_mRNA_6_24.R with Windows paths removed, inputs
# parameterized as directories, header-aware DE readers, and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--mirna-de-dir"), type = "character",
              help = "Directory with miRNA_DE_{6h,24h}_{UP,DOWN}.txt (from script 02)"),
  make_option(c("--mrna-de-dir"), type = "character",
              help = "Directory with NE_vs_Ctrl_{6h,24h}_mRNA_DE.txt (from script 01)"),
  make_option(c("--mirwalk-dir"), type = "character",
              help = "Directory with miRWalk_miRNA_Targets_{3UTR,CDS,5UTR}_{UP,DOWN}{6,24}.csv"),
  make_option(c("--kegg-6h"), type = "character", help = "KEGG Only6h (long, SYMBOLS)"),
  make_option(c("--kegg-24h"), type = "character", help = "KEGG Only24h (long, SYMBOLS)"),
  make_option(c("--kegg-shared"), type = "character", help = "KEGG Shared (long, SYMBOLS)"),
  make_option(c("--pathways-selected"), type = "character", default = NULL,
              help = "pathways_selected_for_TF.csv from script 08 (for sign_class)"),
  make_option(c("--phase-dir"), type = "character", default = NULL,
              help = "Directory with early/sustained/late_mRNA.csv (phase map)"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/09_mirna_mRNA",
              help = "Output directory [default: %default]"),
  make_option(c("--alpha"), type = "numeric", default = 0.05, help = "DE FDR threshold"),
  make_option(c("--top-mirnas"), type = "integer", default = 25, help = "Top miRNAs per plot"),
  make_option(c("--top-pathways-auto"), type = "integer", default = 12,
              help = "Top KEGG pathways/group when no --pathways-selected"),
  make_option(c("--dpi"), type = "integer", default = 600, help = "PNG DPI")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

OUTPUT_DIR <- opt$`output-dir`
alpha <- opt$alpha
DPI_PNG <- opt$dpi
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

if (is.null(opt$`kegg-shared`) || is.null(opt$`mirna-de-dir`) || is.null(opt$`mirwalk-dir`) || is.null(opt$`mrna-de-dir`)) {
  cat("[WARN] Required inputs missing (--mirna-de-dir, --mrna-de-dir, --mirwalk-dir, --kegg-shared).\n")
  cat("[INFO] Exiting with placeholder status.\n"); quit(status = 0)
}
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

detect_col <- function(nms, patterns) {
  hit <- which(str_detect(nms, regex(paste(patterns, collapse = "|"), ignore_case = TRUE)))
  if (length(hit) == 0) NA_character_ else nms[hit[1]]
}

# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

read_mirna_de <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- suppressWarnings(read_tsv(path, show_col_types = FALSE))
  nms <- names(df)
  mir_col <- if ("miRNA_id" %in% nms) "miRNA_id" else detect_col(nms, c("^mirna", "mirnaid", "miRNA"))
  lfc_col <- if ("log2FoldChange" %in% nms) "log2FoldChange" else detect_col(nms, c("log2fold", "log2FC"))
  padj_col <- if ("padj" %in% nms) "padj" else detect_col(nms, c("padj", "fdr", "adj"))
  if (is.na(mir_col) || is.na(lfc_col) || is.na(padj_col)) return(NULL)
  df %>% transmute(miRNA = as.character(.data[[mir_col]]),
                   log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
                   padj = suppressWarnings(as.numeric(.data[[padj_col]]))) %>%
    filter(!is.na(miRNA), miRNA != "") %>%
    mutate(miRNA_key = str_to_upper(str_trim(miRNA)), de_sig = !is.na(padj) & padj < alpha) %>%
    group_by(miRNA_key) %>% arrange(padj, desc(abs(log2FC))) %>% slice(1) %>% ungroup()
}

read_mrna_de <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- suppressWarnings(read_tsv(path, show_col_types = FALSE))
  nms <- names(df)
  gene_col <- if ("gene_name" %in% nms) "gene_name" else detect_col(nms, c("gene_name", "symbol", "gene"))
  lfc_col  <- if ("log2FoldChange" %in% nms) "log2FoldChange" else detect_col(nms, c("log2fold", "log2FC"))
  padj_col <- if ("padj" %in% nms) "padj" else detect_col(nms, c("padj", "fdr", "adj"))
  dir_col  <- if ("direction" %in% nms) "direction" else detect_col(nms, c("direction"))
  if (is.na(gene_col) || is.na(lfc_col) || is.na(padj_col)) return(NULL)
  df %>% transmute(gene_name = as.character(.data[[gene_col]]),
                   log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
                   padj = suppressWarnings(as.numeric(.data[[padj_col]])),
                   direction = if (!is.na(dir_col)) as.character(.data[[dir_col]]) else NA_character_) %>%
    mutate(gene_key = str_to_upper(str_trim(gene_name)), de_sig = !is.na(padj) & padj < alpha) %>%
    filter(!is.na(gene_key), gene_key != "") %>%
    group_by(gene_key) %>% arrange(padj, desc(abs(log2FC))) %>% slice(1) %>% ungroup()
}

read_mirwalk_region <- function(path, region_label) {
  if (!file.exists(path)) return(NULL)
  df <- read_csv(path, show_col_types = FALSE); nms <- names(df)
  mir_col <- detect_col(nms, c("^mirna", "mirnaid", "miRNA"))
  gene_col <- detect_col(nms, c("genesymbol", "gene_symbol", "target", "symbol", "gene"))
  if (is.na(mir_col) || is.na(gene_col)) return(NULL)
  df %>% transmute(miRNA = as.character(.data[[mir_col]]), gene_symbol = as.character(.data[[gene_col]]), region = region_label) %>%
    filter(!is.na(miRNA), miRNA != "", !is.na(gene_symbol), gene_symbol != "") %>%
    mutate(miRNA_key = str_to_upper(str_trim(miRNA)), gene_key = str_to_upper(str_trim(gene_symbol))) %>%
    distinct(miRNA_key, gene_key, region)
}

read_phase_map <- function(dir) {
  if (is.null(dir)) return(tibble(gene_key = character(), phase = character()))
  safe <- function(fname, lab) {
    p <- file.path(dir, fname); if (!file.exists(p)) return(NULL)
    x <- suppressWarnings(read_csv(p, show_col_types = FALSE)); if (nrow(x) == 0) return(NULL)
    gcol <- if ("gene_name" %in% names(x)) "gene_name" else detect_col(names(x), c("gene", "symbol", "name"))
    if (is.na(gcol)) return(NULL)
    x %>% transmute(gene_key = str_to_upper(str_trim(as.character(.data[[gcol]]))), phase = lab) %>% filter(gene_key != "")
  }
  bind_rows(safe("sustained_mRNA.csv", "Sustained"), safe("early_mRNA.csv", "Early"), safe("late_mRNA.csv", "Late")) %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_kegg_map <- function() {
  bind_rows(
    read_tsv(opt$`kegg-6h`, show_col_types = FALSE) %>% transmute(KEGG_group = "Only6h", ID, Description,
      padj = as.numeric(`p.adjust`), gene_key = str_to_upper(str_trim(as.character(gene_symbol)))),
    read_tsv(opt$`kegg-24h`, show_col_types = FALSE) %>% transmute(KEGG_group = "Only24h", ID, Description,
      padj = as.numeric(`p.adjust`), gene_key = str_to_upper(str_trim(as.character(gene_symbol)))),
    read_tsv(opt$`kegg-shared`, show_col_types = FALSE) %>% transmute(KEGG_group = "Shared", ID, Description,
      padj = pmin(as.numeric(`p.adjust_6h`), as.numeric(`p.adjust_24h`), na.rm = TRUE),
      gene_key = str_to_upper(str_trim(as.character(gene_symbol))))
  ) %>% filter(!is.na(gene_key), gene_key != "") %>% distinct(KEGG_group, ID, Description, gene_key, padj)
}

# ---------------------------------------------------------------------------
# Build KEGG map (+ sign_class from script 08 selection, if provided)
# ---------------------------------------------------------------------------

kegg_map <- read_kegg_map()
if (!is.null(opt$`pathways-selected`) && file.exists(opt$`pathways-selected`)) {
  sel <- read_csv(opt$`pathways-selected`, show_col_types = FALSE) %>%
    select(KEGG_group, ID, Description, sign_class) %>% distinct()
  kegg_map <- kegg_map %>% inner_join(sel, by = c("KEGG_group", "ID", "Description"))
} else {
  top_paths <- kegg_map %>% group_by(KEGG_group, Description) %>%
    summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
    group_by(KEGG_group) %>% arrange(padj_min, .by_group = TRUE) %>%
    slice_head(n = opt$`top-pathways-auto`) %>% ungroup() %>% pull(Description) %>% unique()
  kegg_map <- kegg_map %>% filter(Description %in% top_paths) %>% mutate(sign_class = "Auto")
}

phase_map <- read_phase_map(opt$`phase-dir`)

# ---------------------------------------------------------------------------
# Reciprocal combos
# ---------------------------------------------------------------------------

mirna_de_file <- function(tp, dir) file.path(opt$`mirna-de-dir`, paste0("miRNA_DE_", tp, "_", toupper(dir), ".txt"))
mrna_de_file  <- function(tp) file.path(opt$`mrna-de-dir`, paste0("NE_vs_Ctrl_", tp, "_mRNA_DE.txt"))
mirwalk_file  <- function(region, dir, tpnum) file.path(opt$`mirwalk-dir`,
  paste0("miRWalk_miRNA_Targets_", region, "_", toupper(dir), tpnum, ".csv"))

run_combo <- function(timepoint, mir_dir) {
  mr_dir <- if (mir_dir == "Up") "Down" else "Up"
  tpnum <- if (timepoint == "24h") "24" else "6"
  tag <- paste0(timepoint, "_miRNA", mir_dir, "_mRNA", mr_dir)

  mir_de <- read_mirna_de(mirna_de_file(timepoint, mir_dir))
  mr_de  <- read_mrna_de(mrna_de_file(timepoint))
  t3 <- read_mirwalk_region(mirwalk_file("3UTR", mir_dir, tpnum), "3UTR")
  tc <- read_mirwalk_region(mirwalk_file("CDS", mir_dir, tpnum), "CDS")
  t5 <- read_mirwalk_region(mirwalk_file("5UTR", mir_dir, tpnum), "5UTR")

  if (is.null(mir_de) || is.null(mr_de) || (is.null(t3) && is.null(tc) && is.null(t5))) {
    cat("[SKIP]", tag, "(missing miRNA DE / mRNA DE / miRWalk targets)\n"); return(NULL)
  }

  mir_de <- mir_de %>% filter(de_sig) %>% mutate(time = timepoint, miRNA_dir = mir_dir)
  mr_de <- mr_de %>%
    mutate(mRNA_dir = case_when(
      !is.na(direction) & str_detect(direction, regex("DOWN", TRUE)) ~ "Down",
      !is.na(direction) & str_detect(direction, regex("UP", TRUE)) ~ "Up",
      !is.na(log2FC) & log2FC < 0 ~ "Down",
      !is.na(log2FC) & log2FC > 0 ~ "Up", TRUE ~ NA_character_)) %>%
    filter(de_sig, mRNA_dir == mr_dir)

  targets_collapsed <- bind_rows(t3, tc, t5) %>%
    group_by(miRNA_key, gene_key) %>%
    summarize(region_support = paste(sort(unique(region)), collapse = ";"), .groups = "drop")

  edges <- mir_de %>%
    select(miRNA, miRNA_key, miRNA_log2FC = log2FC, miRNA_padj = padj, time, miRNA_dir) %>%
    inner_join(targets_collapsed, by = "miRNA_key") %>%
    inner_join(mr_de %>% select(gene_key, mRNA_log2FC = log2FC, mRNA_padj = padj, mRNA_dir), by = "gene_key") %>%
    left_join(phase_map, by = "gene_key") %>% mutate(phase = replace_na(phase, "Unclassified"), combo = tag)

  edges_kegg <- edges %>% inner_join(kegg_map, by = "gene_key", relationship = "many-to-many") %>%
    mutate(KEGG_group = factor(KEGG_group, levels = c("Only6h", "Shared", "Only24h")))
  if (nrow(edges_kegg) == 0) { cat("[SKIP]", tag, "(no edges after KEGG join)\n"); return(NULL) }

  write_csv(edges_kegg, file.path(OUTPUT_DIR, paste0("SUPP_edges_", tag, "_KEGG_phase.csv")))

  counts <- edges_kegg %>% group_by(KEGG_group, Description, sign_class, miRNA) %>%
    summarize(n_targets = n_distinct(gene_key), .groups = "drop")
  top_mir <- counts %>% group_by(miRNA) %>% summarize(total = sum(n_targets), .groups = "drop") %>%
    arrange(desc(total)) %>% slice_head(n = opt$`top-mirnas`) %>% pull(miRNA)
  plot_df <- counts %>% filter(miRNA %in% top_mir) %>%
    mutate(miRNA = fct_reorder(miRNA, n_targets, .fun = sum, .desc = TRUE), Description = fct_inorder(Description))

  p <- ggplot(plot_df, aes(x = Description, y = miRNA)) +
    geom_point(aes(size = n_targets), alpha = 0.85) +
    facet_grid(sign_class ~ KEGG_group, scales = "free_x", space = "free_x") +
    labs(x = "KEGG pathway (DE targets in this context)", y = paste0("miRNAs ", mir_dir, " (", timepoint, ")"),
         size = "# targets", title = paste0(tag, ": miRNA", mir_dir, " -> mRNA", mr_dir, " (KEGG)")) +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  ggsave(file.path(OUTPUT_DIR, paste0("Fig_miRNA_", tag, "_dotplot.pdf")), p, width = 13.5, height = 9)
  ggsave(file.path(OUTPUT_DIR, paste0("Fig_miRNA_", tag, "_dotplot.png")), p, width = 13.5, height = 9, dpi = DPI_PNG)

  rank <- counts %>% group_by(miRNA) %>%
    summarize(total_targets = sum(n_targets), n_pathways_hit = sum(n_targets > 0), .groups = "drop") %>%
    arrange(desc(total_targets))
  write_csv(rank, file.path(OUTPUT_DIR, paste0("SUPP_miRNA_rank_", tag, ".csv")))
  cat("[OK]", tag, "->", nrow(edges_kegg), "edges\n")
  invisible(edges_kegg)
}

all_edges <- bind_rows(
  run_combo("24h", "Up"), run_combo("24h", "Down"),
  run_combo("6h", "Up"),  run_combo("6h", "Down")
)
if (!is.null(all_edges) && nrow(all_edges) > 0) {
  write_csv(all_edges, file.path(OUTPUT_DIR, "SUPP_edges_ALL_COMBOS_KEGG_phase.csv"))
}

cat("\n[SUCCESS] miRNA-mRNA integration complete.\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
