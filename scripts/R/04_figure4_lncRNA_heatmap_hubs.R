#!/usr/bin/env Rscript
# Stage 04: Figure 4/5 — lncRNA heatmaps, co-expression hubs, and rewiring
# Builds the lncRNA log2FC heatmap by response class (annotated + novel),
# the top-hub bar plot per timepoint, and the 6h->24h rewiring edge counts.
#
# Ported from 06_Figure4_refine_heatmap_hubs.R with Windows paths removed,
# optparse added, and external co-expression tables read gracefully (skipped
# with a warning if not provided).

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(tibble)
  library(purrr)
})

option_list <- list(
  make_option(c("-d", "--de-dir"), type = "character", default = NULL,
              help = "Directory with lncRNA phase sets, DE results, hubs, correlations"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/04_figure4",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 600,
              help = "PNG DPI [default: %default]"),
  make_option(c("--top-per-phase"), type = "integer", default = 10,
              help = "Top lncRNAs per phase for heatmap [default: %default]"),
  make_option(c("--top-hubs"), type = "integer", default = 15,
              help = "Top hubs per timepoint [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
if (is.null(opt$`de-dir`)) { print_help(parser); stop("\n[ERROR] Required argument: --de-dir") }

DE_DIR <- opt$`de-dir`
OUTPUT_DIR <- opt$`output-dir`
DPI_PNG <- opt$dpi
TOP_LNC_PER_PHASE <- opt$`top-per-phase`
TOP_HUBS_PER_TIME <- opt$`top-hubs`
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[INFO] DE directory:", DE_DIR, "\n")
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

locate_file <- function(fname) {
  p <- file.path(DE_DIR, fname)
  if (!file.exists(p)) return(NA_character_)
  p
}
std_chr <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }

save_gg <- function(plot_obj, base_name, width, height, dpi = DPI_PNG) {
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".pdf")), plot = plot_obj,
         width = width, height = height, units = "in")
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".png")), plot = plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

read_lnc_set <- function(fname, phase_label, source_label) {
  path <- locate_file(fname)
  if (is.na(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      gene_id = std_chr(gene_id),
      gene_key = str_to_upper(gene_id),
      phase = phase_label,
      lnc_source = source_label
    ) %>%
    filter(gene_key != "") %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_lnc_de <- function(fname, time_label) {
  path <- locate_file(fname)
  if (is.na(path)) return(NULL)
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      gene_id = std_chr(gene_id),
      gene_key = str_to_upper(gene_id),
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      time = time_label
    ) %>%
    group_by(gene_key, time) %>% arrange(padj, desc(abs(log2FC))) %>% slice(1) %>% ungroup()
}

# ---------------------------------------------------------------------------
# 1. lncRNA phase sets (annotated + novel) and DE results
# ---------------------------------------------------------------------------

cat("[INFO] Loading lncRNA phase assignments...\n")
lnc_all <- bind_rows(
  read_lnc_set("early_lncRNA.csv", "Early", "annotated"),
  read_lnc_set("sustained_lncRNA.csv", "Sustained", "annotated"),
  read_lnc_set("late_lncRNA.csv", "Late", "annotated"),
  read_lnc_set("early_novel_lncRNA.csv", "Early", "novel"),
  read_lnc_set("sustained_novel_lncRNA.csv", "Sustained", "novel"),
  read_lnc_set("late_novel_lncRNA.csv", "Late", "novel")
)
if (is.null(lnc_all) || nrow(lnc_all) == 0)
  stop("[ERROR] No lncRNA phase sets found in ", DE_DIR, " (expected early_lncRNA.csv, ...)")

lnc_all <- lnc_all %>%
  mutate(phase = factor(phase, levels = c("Early", "Sustained", "Late")),
         lnc_source = factor(lnc_source, levels = c("annotated", "novel")))
lnc_phase_map <- lnc_all %>% distinct(gene_key, .keep_all = TRUE)

lnc_de_all <- bind_rows(
  read_lnc_de("NE_vs_Ctrl_6h_lncRNA_DE.txt", "6h"),
  read_lnc_de("NE_vs_Ctrl_24h_lncRNA_DE.txt", "24h")
)
if (is.null(lnc_de_all)) stop("[ERROR] lncRNA DE tables not found (NE_vs_Ctrl_*_lncRNA_DE.txt)")
lnc_de_all <- lnc_de_all %>% mutate(time = factor(time, levels = c("6h", "24h")))

lnc_wide <- lnc_phase_map %>%
  select(gene_id, gene_key, phase, lnc_source) %>%
  left_join(lnc_de_all, by = c("gene_id", "gene_key")) %>%
  pivot_wider(names_from = time, values_from = c(log2FC, padj), names_sep = "_") %>%
  mutate(padj_min = suppressWarnings(pmin(padj_6h, padj_24h, na.rm = TRUE)),
         max_absFC = pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE))
write_csv(lnc_wide, file.path(OUTPUT_DIR, "Fig4_lnc_wide.csv"))

# ---------------------------------------------------------------------------
# 2. Fig4D: heatmap of top lncRNAs by response class
# ---------------------------------------------------------------------------

if (requireNamespace("pheatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pheatmap))
  heatmap_df <- lnc_wide %>%
    group_by(phase) %>% arrange(padj_min, desc(max_absFC), .by_group = TRUE) %>%
    slice_head(n = TOP_LNC_PER_PHASE) %>% ungroup() %>%
    distinct(gene_id, .keep_all = TRUE) %>% arrange(phase, desc(max_absFC))
  write_csv(heatmap_df, file.path(OUTPUT_DIR, "Fig4D_heatmap_genes.csv"))

  mat <- heatmap_df %>% select(gene_id, log2FC_6h, log2FC_24h) %>%
    column_to_rownames("gene_id") %>% as.matrix()
  ann <- heatmap_df %>% select(gene_id, phase, lnc_source) %>% column_to_rownames("gene_id")

  pdf(file.path(OUTPUT_DIR, "Fig4D_heatmap_lnc.pdf"), width = 7.2, height = 10)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
           main = "Top lncRNAs by response class", fontsize_row = 7, border_color = "grey70")
  dev.off()
  png(file.path(OUTPUT_DIR, "Fig4D_heatmap_lnc.png"), width = 7.2, height = 10, units = "in", res = DPI_PNG)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
           main = "Top lncRNAs by response class", fontsize_row = 7, border_color = "grey70")
  dev.off()
  cat("[INFO] Fig4D heatmap saved\n")
}

# ---------------------------------------------------------------------------
# 3. Fig4E: top lncRNA co-expression hubs per timepoint
#    Reads external co-expression table lncRNA_hubs_by_time.csv (skip if absent)
# ---------------------------------------------------------------------------

hubs_path <- locate_file("lncRNA_hubs_by_time.csv")
if (!is.na(hubs_path)) {
  hubs_by_time <- read_csv(hubs_path, show_col_types = FALSE) %>%
    transmute(
      time = factor(std_chr(time), levels = c("6h", "24h", "Early", "Late")),
      lncRNA = std_chr(lncRNA),
      lnc_source = factor(std_chr(if ("lnc_source" %in% names(.)) lnc_source else "unknown"),
                          levels = c("annotated", "novel", "unknown")),
      degree = dplyr::coalesce(suppressWarnings(as.numeric(if ("n_mRNAs" %in% names(.)) n_mRNAs else NA)),
                               suppressWarnings(as.numeric(if ("n_edges" %in% names(.)) n_edges else NA))),
      mean_rho = suppressWarnings(as.numeric(if ("mean_rho" %in% names(.)) mean_rho else NA))
    ) %>%
    filter(!is.na(degree), !is.na(time))

  hubs_top <- hubs_by_time %>%
    group_by(time) %>% arrange(desc(degree), .by_group = TRUE) %>%
    slice_head(n = TOP_HUBS_PER_TIME) %>% ungroup() %>%
    mutate(label = fct_reorder(lncRNA, degree))
  write_csv(hubs_top, file.path(OUTPUT_DIR, "Fig4E_top_hubs_by_time.csv"))

  p_hubs <- ggplot(hubs_top, aes(x = degree, y = label, fill = lnc_source)) +
    geom_col() + facet_wrap(~ time, scales = "free_y") +
    labs(x = "# mRNA targets (degree)", y = NULL, fill = "Source",
         title = "Top lncRNA hubs across commitment and maintenance") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  save_gg(p_hubs, "Fig4E_lncRNA_hubs_by_time", 11, 6.4)
  cat("[INFO] Fig4E hub plot saved\n")
} else {
  cat("[WARN] lncRNA_hubs_by_time.csv not found; skipping Fig4E hub plot\n")
}

# ---------------------------------------------------------------------------
# 4. Rewiring edge counts between 6h and 24h (external correlation tables)
# ---------------------------------------------------------------------------

edge_specs <- c(
  Lost_after_6h = "correlations_lost_after_6h.csv",
  Gained_at_24h = "correlations_gained_at_24h.csv",
  Late_only     = "correlations_late_only.csv"
)
edge_counts <- imap_dfr(edge_specs, function(fname, nm) {
  path <- locate_file(fname)
  if (is.na(path)) return(NULL)
  tibble(category = nm, n_edges = nrow(read_csv(path, show_col_types = FALSE)))
})

if (nrow(edge_counts) > 0) {
  write_csv(edge_counts, file.path(OUTPUT_DIR, "Fig4_rewiring_edge_counts.csv"))
  p_rewire <- ggplot(edge_counts, aes(x = category, y = n_edges)) +
    geom_col(fill = "grey35") +
    labs(x = NULL, y = "# lncRNA-mRNA correlations", title = "Rewiring categories between 6h and 24h") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
  save_gg(p_rewire, "Fig4_rewiring_edge_counts", 7.8, 4.4)
  cat("[INFO] Rewiring edge-count plot saved\n")
} else {
  cat("[WARN] No correlations_*.csv found; skipping rewiring plot\n")
}

cat("\n[SUCCESS] Figure 4 lncRNA heatmap + hubs complete!\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
