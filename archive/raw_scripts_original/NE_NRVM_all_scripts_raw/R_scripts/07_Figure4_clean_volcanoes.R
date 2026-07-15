suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

base_dir  <- "."
bulk_dir  <- file.path(base_dir, "DESeq2_Multifactorial_Results_fixed")
out_dir   <- file.path(base_dir, "Paper_Fig4_FINAL_VOLCANOES")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dpi_png <- 600
padj_cutoff <- 0.05
absfc_cutoff <- 0
n_labels_per_source_direction <- 6

pdf_device <- function(...) {
  if (capabilities("cairo")) grDevices::cairo_pdf(...) else grDevices::pdf(...)
}

save_gg <- function(plot_obj, base_name, width, height, dpi = 600) {
  ggsave(file.path(out_dir, paste0(base_name, ".pdf")), plot = plot_obj,
         width = width, height = height, units = "in", device = pdf_device)
  ggsave(file.path(out_dir, paste0(base_name, "_600dpi.png")), plot = plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
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
        !is.na(padj) & padj < padj_cutoff & log2FC > absfc_cutoff ~ "UP",
        !is.na(padj) & padj < padj_cutoff & log2FC < -absfc_cutoff ~ "DOWN",
        TRUE ~ "NS"
      ),
      neglog10padj = -log10(pmax(padj, 1e-300))
    ) %>%
    mutate(
      source = factor(source, levels = c("annotated", "novel", "unknown")),
      direction = factor(direction, levels = c("DOWN", "NS", "UP"))
    )

  label_df <- df %>%
    filter(direction %in% c("UP", "DOWN")) %>%
    group_by(source, direction) %>%
    arrange(padj, desc(abs(log2FC)), .by_group = TRUE) %>%
    slice_head(n = n_labels_per_source_direction) %>%
    ungroup()

  write_csv(df, file.path(out_dir, paste0("SUPP_", time_label, "_lnc_volcano_data.csv")))

  p <- ggplot(df, aes(x = log2FC, y = neglog10padj)) +
    geom_vline(xintercept = c(-1, 1), linetype = 3, color = "grey65") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = 3, color = "grey65") +
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

  save_gg(p, paste0("Fig4", time_label, "_lncRNA_volcano_clean"), 6.4, 5.4, dpi_png)
}

make_lnc_volcano("NE_vs_Ctrl_6h_lncRNA_DE.txt", "6h")
make_lnc_volcano("NE_vs_Ctrl_24h_lncRNA_DE.txt", "24h")

message("Volcanos limpios de Figure 4 listos en: ", out_dir)
