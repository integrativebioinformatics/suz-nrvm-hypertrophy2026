#!/usr/bin/env Rscript
# Stage 03: Figure 3 — mRNA Temporal Programs (Early/Sustained/Late)
# Visualizes mRNA response phenotypes over time with pathway integration.
# Ported from 05_Figure3_mRNA_programs.R with Windows paths removed and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(ggrepel)
  library(tibble)
})

# ============================================================================
# Command-line Arguments
# ============================================================================

option_list <- list(
  make_option(c("-d", "--de-dir"), type = "character", default = NULL,
              help = "Directory containing DESeq2 results (mRNA DE + phase CSVs)"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/03_figure3",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 600,
              help = "PNG DPI for high-res output [default: %default]"),
  make_option(c("--top-per-phase"), type = "integer", default = 12,
              help = "Number of top mRNAs per phase for heatmap [default: %default]"),
  make_option(c("--top-labels"), type = "integer", default = 18,
              help = "Number of top genes to label in scatter plot [default: %default]")
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
TOP_MRNA_PER_PHASE <- opt$`top-per-phase`
TOP_LABELS_SCATTER <- opt$`top-labels`

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[INFO] DE directory:", DE_DIR, "\n")
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

# ============================================================================
# Helper Functions
# ============================================================================

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

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit)) return(NA_character_)
  hit
}

read_mrna_set <- function(path, phase_label) {
  df <- read_csv(path, show_col_types = FALSE)
  gene_id_col   <- first_existing_col(df, c("gene_id", "GeneID", "ENSEMBL"))
  gene_name_col <- first_existing_col(df, c("gene_name", "gene_symbol", "symbol"))
  df %>%
    transmute(
      gene_id = std_chr(.data[[gene_id_col]]),
      gene_symbol = ifelse(is.na(gene_name_col), gene_id, std_chr(.data[[gene_name_col]])),
      phase = phase_label,
      gene_key = str_to_upper(gene_id)
    ) %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_mrna_de <- function(path, time_label) {
  df <- read_tsv(path, show_col_types = FALSE)
  gene_id_col   <- first_existing_col(df, c("gene_id", "GeneID", "ENSEMBL"))
  gene_name_col <- first_existing_col(df, c("gene_name", "gene_symbol", "symbol"))
  lfc_col       <- first_existing_col(df, c("log2FoldChange", "log2FC"))
  padj_col      <- first_existing_col(df, c("padj", "FDR", "adj.P.Val"))

  df %>%
    transmute(
      gene_id = std_chr(.data[[gene_id_col]]),
      gene_symbol = ifelse(is.na(gene_name_col), gene_id, std_chr(.data[[gene_name_col]])),
      gene_key = str_to_upper(gene_id),
      log2FC = as.numeric(.data[[lfc_col]]),
      padj = as.numeric(.data[[padj_col]]),
      direction = ifelse(log2FC > 0, "UP", ifelse(log2FC < 0, "DOWN", "NS")),
      time = time_label
    ) %>%
    group_by(gene_key, time) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup()
}

save_gg <- function(plot_obj, base_name, width, height, dpi = DPI_PNG) {
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".pdf")),
         plot = plot_obj, width = width, height = height, units = "in")
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".png")),
         plot = plot_obj, width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

# ============================================================================
# Load mRNA Phase Assignments and DE Results
# ============================================================================

cat("[INFO] Loading mRNA phase assignments...\n")
early_set     <- read_mrna_set(locate_file("early_mRNA.csv"), "Early")
sustained_set <- read_mrna_set(locate_file("sustained_mRNA.csv"), "Sustained")
late_set      <- read_mrna_set(locate_file("late_mRNA.csv"), "Late")

mrna_phase_map <- bind_rows(early_set, sustained_set, late_set) %>%
  mutate(phase = factor(phase, levels = c("Early", "Sustained", "Late"))) %>%
  distinct(gene_key, .keep_all = TRUE)

cat("[INFO] Loaded mRNA phases:", nrow(mrna_phase_map), "genes\n")

cat("[INFO] Loading DE results...\n")
mrna_de <- bind_rows(
  read_mrna_de(locate_file("NE_vs_Ctrl_6h_mRNA_DE.txt"), "6h"),
  read_mrna_de(locate_file("NE_vs_Ctrl_24h_mRNA_DE.txt"), "24h")
) %>%
  mutate(time = factor(time, levels = c("6h", "24h")))

# Join phases with DE results
mrna_join <- mrna_phase_map %>%
  select(gene_key, gene_id, gene_symbol, phase) %>%
  full_join(mrna_de %>% select(gene_key, time, log2FC, padj, direction), by = "gene_key") %>%
  group_by(gene_key) %>%
  fill(gene_id, gene_symbol, phase, .direction = "downup") %>%
  ungroup()

# Create wide format with classification
mrna_wide <- mrna_join %>%
  select(gene_key, gene_id, gene_symbol, phase, time, log2FC, padj, direction) %>%
  pivot_wider(names_from = time, values_from = c(log2FC, padj, direction), names_sep = "_") %>%
  mutate(
    phase = factor(phase, levels = c("Early", "Sustained", "Late")),
    class_direction = case_when(
      phase == "Early" ~ direction_6h,
      phase == "Late" ~ direction_24h,
      phase == "Sustained" & !is.na(direction_6h) & !is.na(direction_24h) & direction_6h == direction_24h ~ direction_6h,
      phase == "Sustained" & !is.na(log2FC_6h) & !is.na(log2FC_24h) & sign(log2FC_6h) != sign(log2FC_24h) ~ "REVERSAL",
      phase == "Sustained" & rowMeans(cbind(log2FC_6h, log2FC_24h), na.rm = TRUE) > 0 ~ "UP",
      phase == "Sustained" & rowMeans(cbind(log2FC_6h, log2FC_24h), na.rm = TRUE) < 0 ~ "DOWN",
      TRUE ~ "UNDEF"
    ),
    padj_min = suppressWarnings(pmin(padj_6h, padj_24h, na.rm = TRUE)),
    max_absFC = pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE)
  )

write_csv(mrna_wide, file.path(OUTPUT_DIR, "Fig3_mRNA_programs_wide.csv"))

# ============================================================================
# Figure 3A: Scheme
# ============================================================================

scheme_df <- tibble(
  phase = factor(c("Early", "Sustained", "Late"), levels = c("Early", "Sustained", "Late")),
  definition = c("DE only at 6h", "DE at 6h and 24h", "DE only at 24h"),
  x = c(1, 2, 3),
  y = 1
)

p_Fig3A <- ggplot(scheme_df, aes(x = x, y = y, fill = phase)) +
  geom_tile(width = 0.9, height = 0.5, color = "grey25") +
  geom_text(aes(label = paste0(as.character(phase), "\n", definition)), size = 4) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(title = "mRNA response classes: commitment vs maintenance") +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

save_gg(p_Fig3A, "Fig3A_mRNA_class_scheme", 9, 2.4)

# ============================================================================
# Figure 3B: Counts by Phase and Direction
# ============================================================================

counts_df <- mrna_wide %>%
  filter(class_direction %in% c("UP", "DOWN", "REVERSAL")) %>%
  count(phase, class_direction, name = "n_genes") %>%
  mutate(class_direction = factor(class_direction, levels = c("UP", "DOWN", "REVERSAL")))

write_csv(counts_df, file.path(OUTPUT_DIR, "Fig3B_mRNA_counts_by_phase.csv"))

p_Fig3B <- ggplot(counts_df, aes(x = phase, y = n_genes, fill = class_direction)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "# mRNA genes", fill = "Direction",
       title = "mRNA response across timepoints") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_gg(p_Fig3B, "Fig3B_mRNA_counts", 7.5, 4.6)

# ============================================================================
# Figure 3C: Heatmap of Top mRNAs
# ============================================================================

if (requireNamespace("pheatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pheatmap))

  heatmap_df <- mrna_wide %>%
    group_by(phase) %>%
    arrange(padj_min, desc(max_absFC), .by_group = TRUE) %>%
    slice_head(n = TOP_MRNA_PER_PHASE) %>%
    ungroup() %>%
    mutate(row_label = gene_symbol) %>%
    distinct(row_label, .keep_all = TRUE) %>%
    arrange(phase, desc(max_absFC))

  mat <- heatmap_df %>%
    select(row_label, log2FC_6h, log2FC_24h) %>%
    column_to_rownames("row_label") %>%
    as.matrix()

  ann <- heatmap_df %>%
    select(row_label, phase, class_direction) %>%
    column_to_rownames("row_label")

  write_csv(heatmap_df, file.path(OUTPUT_DIR, "Fig3C_heatmap_genes.csv"))

  pdf(file.path(OUTPUT_DIR, "Fig3C_heatmap_mRNA.pdf"), width = 7.2, height = 10)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
           main = "Top mRNAs by response class (log2FC)", fontsize_row = 7, border_color = "grey70")
  dev.off()

  png(file.path(OUTPUT_DIR, "Fig3C_heatmap_mRNA.png"), width = 7.2, height = 10, units = "in", res = DPI_PNG)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
           main = "Top mRNAs by response class (log2FC)", fontsize_row = 7, border_color = "grey70")
  dev.off()

  cat("[INFO] Heatmap saved\n")
}

# ============================================================================
# Figure 3D: Scatter of Sustained mRNAs (6h vs 24h)
# ============================================================================

scatter_df <- mrna_wide %>%
  filter(phase == "Sustained", !is.na(log2FC_6h), !is.na(log2FC_24h)) %>%
  mutate(
    transition = case_when(
      log2FC_6h > 0 & log2FC_24h > 0 ~ "Maintained_UP",
      log2FC_6h < 0 & log2FC_24h < 0 ~ "Maintained_DOWN",
      TRUE ~ "Reversal"
    ),
    label_rank = rank(padj_min, ties.method = "first"),
    label_me = label_rank <= TOP_LABELS_SCATTER
  )

write_csv(scatter_df, file.path(OUTPUT_DIR, "Fig3D_scatter_data.csv"))

p_Fig3D <- ggplot(scatter_df, aes(x = log2FC_6h, y = log2FC_24h, color = transition)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey35") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey35") +
  geom_point(alpha = 0.75, size = 2.1) +
  ggrepel::geom_text_repel(
    data = subset(scatter_df, label_me),
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(x = "log2FC (6h)", y = "log2FC (24h)", color = NULL,
       title = "Sustained mRNA: commitment vs maintenance") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_gg(p_Fig3D, "Fig3D_scatter", 7.2, 6.2)

cat("\n[SUCCESS] Figure 3 complete!\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
