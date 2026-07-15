suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(ggrepel)
  library(tibble)
})

base_dir   <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"
bulk_dir   <- file.path(base_dir, "DESeq2_Multifactorial_Results_fixed")
out_dir    <- file.path(base_dir, "Paper_Fig3_FINAL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dpi_png <- 600
top_mrna_per_phase_heatmap <- 12
top_labels_scatter <- 18

pdf_device <- function(...) {
  if (capabilities("cairo")) grDevices::cairo_pdf(...) else grDevices::pdf(...)
}

save_gg <- function(plot_obj, base_name, width, height, dpi = 600) {
  ggsave(file.path(out_dir, paste0(base_name, ".pdf")),
         plot = plot_obj, width = width, height = height, units = "in",
         device = pdf_device)
  ggsave(file.path(out_dir, paste0(base_name, "_600dpi.png")),
         plot = plot_obj, width = width, height = height, units = "in",
         dpi = dpi, bg = "white")
}

locate_file <- function(fname, dirs = c(bulk_dir, base_dir)) {
  hits <- file.path(dirs, fname)
  hit <- hits[file.exists(hits)][1]
  if (is.na(hit)) stop("No se encontró el archivo: ", fname)
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
  gene_name_col <- first_existing_col(df, c("gene_name", "gene_symbol", "symbol", "Gene.Symbol"))
  source_col    <- first_existing_col(df, c("source"))
  df %>%
    transmute(
      gene_id = std_chr(.data[[gene_id_col]]),
      gene_symbol_raw = if (!is.na(gene_name_col)) std_chr(.data[[gene_name_col]]) else "",
      gene_symbol = ifelse(gene_symbol_raw == "" | str_to_upper(gene_symbol_raw) == str_to_upper(gene_id), gene_id, gene_symbol_raw),
      source = if (!is.na(source_col)) std_chr(.data[[source_col]]) else NA_character_,
      phase = phase_label,
      gene_key = str_to_upper(gene_id)
    ) %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_mrna_de <- function(path, time_label) {
  df <- read_tsv(path, show_col_types = FALSE)
  gene_id_col   <- first_existing_col(df, c("gene_id", "GeneID", "ENSEMBL"))
  gene_name_col <- first_existing_col(df, c("gene_name", "gene_symbol", "symbol", "Gene.Symbol"))
  lfc_col       <- first_existing_col(df, c("log2FoldChange", "log2FC"))
  padj_col      <- first_existing_col(df, c("padj", "FDR", "adj.P.Val"))
  dir_col       <- first_existing_col(df, c("direction", "Direction"))

  df %>%
    transmute(
      gene_id = std_chr(.data[[gene_id_col]]),
      gene_symbol_raw = if (!is.na(gene_name_col)) std_chr(.data[[gene_name_col]]) else "",
      gene_symbol = ifelse(gene_symbol_raw == "" | str_to_upper(gene_symbol_raw) == str_to_upper(gene_id), gene_id, gene_symbol_raw),
      gene_key = str_to_upper(gene_id),
      log2FC = as.numeric(.data[[lfc_col]]),
      padj = as.numeric(.data[[padj_col]]),
      direction = if (!is.na(dir_col)) std_chr(.data[[dir_col]]) else ifelse(log2FC > 0, "UP", ifelse(log2FC < 0, "DOWN", "NS")),
      time = time_label
    ) %>%
    group_by(gene_key, time) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup()
}

early_set     <- read_mrna_set(locate_file("early_mRNA.csv"), "Early")
sustained_set <- read_mrna_set(locate_file("sustained_mRNA.csv"), "Sustained")
late_set      <- read_mrna_set(locate_file("late_mRNA.csv"), "Late")

mrna_phase_map <- bind_rows(early_set, sustained_set, late_set) %>%
  mutate(phase = factor(phase, levels = c("Early", "Sustained", "Late"))) %>%
  distinct(gene_key, .keep_all = TRUE)

mrna_de <- bind_rows(
  read_mrna_de(locate_file("NE_vs_Ctrl_6h_mRNA_DE.txt"), "6h"),
  read_mrna_de(locate_file("NE_vs_Ctrl_24h_mRNA_DE.txt"), "24h")
) %>%
  mutate(time = factor(time, levels = c("6h", "24h")))

mrna_join <- mrna_phase_map %>%
  select(gene_key, gene_id, gene_symbol, phase) %>%
  full_join(mrna_de %>% select(gene_key, time, log2FC, padj, direction), by = "gene_key") %>%
  group_by(gene_key) %>%
  fill(gene_id, gene_symbol, phase, .direction = "downup") %>%
  ungroup()

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

write_csv(mrna_wide, file.path(out_dir, "SUPP_Fig3_mRNA_programs_wide.csv"))

scheme_df <- tibble(
  phase = factor(c("Early", "Sustained", "Late"), levels = c("Early", "Sustained", "Late")),
  definition = c("DE only at 6h", "DE at 6h and 24h", "DE only at 24h"),
  x = c(1, 2, 3),
  y = 1
)

p_Fig3A <- ggplot(scheme_df, aes(x = x, y = y, fill = phase)) +
  geom_tile(width = 0.9, height = 0.5, color = "grey25") +
  geom_text(aes(label = paste0(as.character(phase), "\n", definition)), size = 4, lineheight = 1.1) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(title = "mRNA response classes linking commitment and maintenance") +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))
save_gg(p_Fig3A, "Fig3A_mRNA_class_scheme", 9, 2.4, dpi_png)

counts_df <- mrna_wide %>%
  filter(class_direction %in% c("UP", "DOWN", "REVERSAL")) %>%
  count(phase, class_direction, name = "n_genes") %>%
  mutate(class_direction = factor(class_direction, levels = c("UP", "DOWN", "REVERSAL")))

write_csv(counts_df, file.path(out_dir, "SUPP_Fig3B_mRNA_counts_by_phase_direction.csv"))

p_Fig3B <- ggplot(counts_df, aes(x = phase, y = n_genes, fill = class_direction)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "# mRNA genes", fill = "Direction",
       title = "mRNA classes across commitment and maintenance") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
save_gg(p_Fig3B, "Fig3B_mRNA_counts_by_phase_direction", 7.5, 4.6, dpi_png)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pheatmap))
  heatmap_df <- mrna_wide %>%
    group_by(phase) %>%
    arrange(padj_min, desc(max_absFC), .by_group = TRUE) %>%
    slice_head(n = top_mrna_per_phase_heatmap) %>%
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

  write_csv(heatmap_df, file.path(out_dir, "SUPP_Fig3C_top_mRNAs_for_heatmap.csv"))

  grDevices::pdf(file.path(out_dir, "Fig3C_heatmap_mRNA_log2FC_top_by_phase.pdf"), width = 7.2, height = 10)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann,
    main = "Top mRNAs by response class (log2FC: NE vs Ctrl)",
    fontsize_row = 7,
    border_color = "grey70"
  )
  dev.off()

  grDevices::png(file.path(out_dir, "Fig3C_heatmap_mRNA_log2FC_top_by_phase_600dpi.png"),
                 width = 7.2, height = 10, units = "in", res = dpi_png)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann,
    main = "Top mRNAs by response class (log2FC: NE vs Ctrl)",
    fontsize_row = 7,
    border_color = "grey70"
  )
  dev.off()
}

scatter_df <- mrna_wide %>%
  filter(phase == "Sustained", !is.na(log2FC_6h), !is.na(log2FC_24h)) %>%
  mutate(
    transition = case_when(
      log2FC_6h > 0 & log2FC_24h > 0 ~ "Maintained_UP",
      log2FC_6h < 0 & log2FC_24h < 0 ~ "Maintained_DOWN",
      TRUE ~ "Reversal"
    ),
    label_rank = rank(padj_min, ties.method = "first"),
    label_me = label_rank <= top_labels_scatter
  )

write_csv(scatter_df, file.path(out_dir, "SUPP_Fig3D_sustained_mRNAs_scatter_data.csv"))

p_Fig3D <- ggplot(scatter_df, aes(x = log2FC_6h, y = log2FC_24h, color = transition)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey35") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey35") +
  geom_point(alpha = 0.75, size = 2.1) +
  ggrepel::geom_text_repel(
    data = subset(scatter_df, label_me),
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  labs(
    x = "log2FC (6h)",
    y = "log2FC (24h)",
    color = NULL,
    title = "Sustained mRNA program: commitment vs maintenance"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
save_gg(p_Fig3D, "Fig3D_scatter_sustained_mRNA_log2FC_6h_vs_24h", 7.2, 6.2, dpi_png)

message("Figure 3 lista en: ", out_dir)
