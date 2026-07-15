#!/usr/bin/env Rscript
# Stage 05: Figure 4 — Volcano Plots (lncRNA Differential Expression)
# Clean volcano plots showing lncRNA response at 6h and 24h timepoints.
# Ported from 07_Figure4_clean_volcanoes.R with Windows paths removed and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

option_list <- list(
  make_option(c("-d", "--de-dir"), type = "character", default = NULL,
              help = "Directory with lncRNA DE results"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/05_figure4_volcanoes",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 600,
              help = "PNG DPI [default: %default]"),
  make_option(c("--padj-cutoff"), type = "numeric", default = 0.05,
              help = "Adjusted p-value cutoff [default: %default]"),
  make_option(c("--top-labels"), type = "integer", default = 6,
              help = "Top genes to label per condition [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$`de-dir`)) {
  print_help(parser)
  stop("\n[ERROR] Required argument: --de-dir")
}

DE_DIR <- opt$`de-dir`
OUTPUT_DIR <- opt$`output-dir`
DPI_PNG <- opt$dpi
PADJ_CUTOFF <- opt$`padj-cutoff`
N_LABELS <- opt$`top-labels`

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[INFO] DE directory:", DE_DIR, "\n")
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

locate_file <- function(fname, dirs = c(DE_DIR)) {
  hits <- file.path(dirs, fname)
  hit <- hits[file.exists(hits)][1]
  if (is.na(hit)) stop("[ERROR] File not found: ", fname)
  hit
}

std_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

save_gg <- function(plot_obj, base_name, width, height, dpi = DPI_PNG) {
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".pdf")), plot = plot_obj,
         width = width, height = height, units = "in")
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".png")), plot = plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

make_lnc_volcano <- function(file_name, time_label) {
  df <- read_tsv(locate_file(file_name), show_col_types = FALSE) %>%
    transmute(
      gene_id = std_chr(gene_id),
      gene_name = if ("gene_name" %in% names(.)) std_chr(gene_name) else "",
      label = gene_id,
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      source = if ("source" %in% names(.)) std_chr(source) else "unknown",
      direction = case_when(
        !is.na(padj) & padj < PADJ_CUTOFF & log2FC > 0 ~ "UP",
        !is.na(padj) & padj < PADJ_CUTOFF & log2FC < 0 ~ "DOWN",
        TRUE ~ "NS"
      ),
      neglog10padj = -log10(pmax(padj, 1e-300))
    ) %>%
    mutate(
      source = factor(ifelse(source == "", "unknown", source),
                      levels = c("annotated", "novel", "unknown")),
      direction = factor(direction, levels = c("DOWN", "NS", "UP"))
    )

  # Label the top N per transcript source x direction (matches original)
  label_df <- df %>%
    filter(direction %in% c("UP", "DOWN")) %>%
    group_by(source, direction) %>%
    arrange(padj, desc(abs(log2FC)), .by_group = TRUE) %>%
    slice_head(n = N_LABELS) %>%
    ungroup()

  write_csv(df, file.path(OUTPUT_DIR, paste0("Fig4_", time_label, "_volcano_data.csv")))

  p <- ggplot(df, aes(x = log2FC, y = neglog10padj)) +
    geom_vline(xintercept = c(-1, 1), linetype = 3, color = "grey65") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = 3, color = "grey65") +
    geom_point(aes(color = direction, shape = source), alpha = 0.75, size = 1.8, stroke = 0.15) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 2.7,
      max.overlaps = Inf,
      box.padding = 0.2,
      point.padding = 0.15,
      segment.size = 0.25,
      min.segment.length = 0
    ) +
    scale_color_manual(values = c(DOWN = "#3B82F6", NS = "grey78", UP = "#EF4444")) +
    scale_shape_manual(values = c(annotated = 16, novel = 17, unknown = 15), drop = FALSE) +
    labs(
      x = expression(log[2]*"FC (NE vs Ctrl)"),
      y = expression(-log[10]*" adjusted P"),
      color = "Direction",
      shape = "Source",
      title = paste0("lncRNA differential expression at ", time_label)
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))

  save_gg(p, paste0("Fig4_volcano_", time_label), 6.4, 5.4)
}

make_lnc_volcano("NE_vs_Ctrl_6h_lncRNA_DE.txt", "6h")
make_lnc_volcano("NE_vs_Ctrl_24h_lncRNA_DE.txt", "24h")

cat("\n[SUCCESS] Volcano plots complete!\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
