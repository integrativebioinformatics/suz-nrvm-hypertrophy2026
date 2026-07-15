suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(tibble)
  library(purrr)
})

base_dir   <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"
bulk_dir   <- file.path(base_dir, "DESeq2_Multifactorial_Results_fixed")
fig45_dir  <- file.path(base_dir, "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT")
out_dir    <- file.path(base_dir, "Paper_Fig4_REFINED")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dpi_png <- 600
top_lnc_per_phase_heatmap <- 10
top_hubs_per_time <- 15

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

locate_file <- function(fname, dirs = c(fig45_dir, bulk_dir, base_dir)) {
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

read_lnc_set <- function(path, phase_label, source_label) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      gene_id = std_chr(gene_id),
      transcript_id = std_chr(transcript_id),
      gene_key = str_to_upper(gene_id),
      phase = phase_label,
      lnc_source = source_label
    ) %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_lnc_de <- function(path, time_label) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      gene_id = std_chr(gene_id),
      gene_key = str_to_upper(gene_id),
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      direction = std_chr(direction),
      time = time_label
    ) %>%
    group_by(gene_key, time) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup()
}

lnc_all <- bind_rows(
  read_lnc_set(locate_file("early_lncRNA.csv"), "Early", "annotated"),
  read_lnc_set(locate_file("sustained_lncRNA.csv"), "Sustained", "annotated"),
  read_lnc_set(locate_file("late_lncRNA.csv"), "Late", "annotated"),
  read_lnc_set(locate_file("early_novel_lncRNA.csv"), "Early", "novel"),
  read_lnc_set(locate_file("sustained_novel_lncRNA.csv"), "Sustained", "novel"),
  read_lnc_set(locate_file("late_novel_lncRNA.csv"), "Late", "novel")
) %>%
  mutate(phase = factor(phase, levels = c("Early", "Sustained", "Late")),
         lnc_source = factor(lnc_source, levels = c("annotated", "novel")))

lnc_phase_map <- lnc_all %>% distinct(gene_key, .keep_all = TRUE)

lnc_de_all <- bind_rows(
  read_lnc_de(locate_file("NE_vs_Ctrl_6h_lncRNA_DE.txt"), "6h"),
  read_lnc_de(locate_file("NE_vs_Ctrl_24h_lncRNA_DE.txt"), "24h")
) %>% mutate(time = factor(time, levels = c("6h", "24h")))

lnc_wide <- lnc_phase_map %>%
  select(gene_id, gene_key, phase, lnc_source) %>%
  left_join(lnc_de_all, by = c("gene_id", "gene_key")) %>%
  pivot_wider(names_from = time, values_from = c(log2FC, padj, direction), names_sep = "_") %>%
  mutate(
    padj_min = suppressWarnings(pmin(padj_6h, padj_24h, na.rm = TRUE)),
    max_absFC = pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE)
  )

write_csv(lnc_wide, file.path(out_dir, "SUPP_Fig4D_lnc_wide_for_heatmap.csv"))

if (requireNamespace("pheatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pheatmap))

  heatmap_df <- lnc_wide %>%
    group_by(phase) %>%
    arrange(padj_min, desc(max_absFC), .by_group = TRUE) %>%
    slice_head(n = top_lnc_per_phase_heatmap) %>%
    ungroup() %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(phase, desc(max_absFC))

  write_csv(heatmap_df, file.path(out_dir, "SUPP_Fig4D_top_lnc_for_heatmap.csv"))

  mat <- heatmap_df %>%
    select(gene_id, log2FC_6h, log2FC_24h) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()

  ann <- heatmap_df %>%
    select(gene_id, phase, lnc_source) %>%
    column_to_rownames("gene_id")

  grDevices::pdf(file.path(out_dir, "Fig4D_heatmap_lnc_log2FC_top_by_phase_geneID.pdf"), width = 7.2, height = 10)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann,
    main = "Top lncRNAs by response class (gene_id labels)",
    fontsize_row = 7,
    border_color = "grey70"
  )
  dev.off()

  grDevices::png(file.path(out_dir, "Fig4D_heatmap_lnc_log2FC_top_by_phase_geneID_600dpi.png"),
                 width = 7.2, height = 10, units = "in", res = dpi_png)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann,
    main = "Top lncRNAs by response class (gene_id labels)",
    fontsize_row = 7,
    border_color = "grey70"
  )
  dev.off()
}

hubs_by_time <- read_csv(locate_file("lncRNA_hubs_by_time.csv"), show_col_types = FALSE) %>%
  transmute(
    time = factor(std_chr(time), levels = c("6h", "24h", "Early", "Late")),
    lncRNA = std_chr(lncRNA),
    lnc_source = factor(std_chr(lnc_source), levels = c("annotated", "novel", "unknown")),
    degree = dplyr::coalesce(as.numeric(n_mRNAs), as.numeric(n_edges)),
    mean_rho = as.numeric(mean_rho)
  ) %>%
  filter(!is.na(degree), !is.na(time))

hubs_top <- hubs_by_time %>%
  group_by(time) %>%
  arrange(desc(degree), .by_group = TRUE) %>%
  slice_head(n = top_hubs_per_time) %>%
  ungroup() %>%
  mutate(label = fct_reorder(lncRNA, degree))

write_csv(hubs_top, file.path(out_dir, "SUPP_Fig4E_top_hubs_by_time.csv"))

p_Fig4E <- ggplot(hubs_top, aes(x = degree, y = label, fill = lnc_source)) +
  geom_col() +
  facet_wrap(~ time, scales = "free_y") +
  labs(
    x = "# mRNA targets (degree)",
    y = NULL,
    fill = "Source",
    title = "Top lncRNA hubs across commitment and maintenance"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
save_gg(p_Fig4E, "Fig4E_lncRNA_hubs_by_time_geneID", 11, 6.4, dpi_png)

edge_files <- c(
  Lost_after_6h = locate_file("correlations_lost_after_6h.csv"),
  Gained_at_24h = locate_file("correlations_gained_at_24h.csv"),
  Late_only = locate_file("correlations_late_only.csv")
)

edge_counts <- purrr::imap_dfr(edge_files, function(path, nm) {
  df <- read_csv(path, show_col_types = FALSE)
  tibble(category = nm, n_edges = nrow(df))
})

write_csv(edge_counts, file.path(out_dir, "SUPP_Fig4_rewiring_edge_counts_curated.csv"))

p_SUPP_rewire <- ggplot(edge_counts, aes(x = category, y = n_edges)) +
  geom_col(fill = "grey35") +
  labs(x = NULL, y = "# lncRNA-mRNA correlations", title = "Rewiring categories between 6h and 24h") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))
save_gg(p_SUPP_rewire, "SUPP_Fig4_rewiring_edge_counts_curated", 7.8, 4.4, dpi_png)

message("Figure 4 refinada lista en: ", out_dir)
