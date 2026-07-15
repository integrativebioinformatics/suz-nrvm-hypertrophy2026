suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

early_lnc_path <- "early_lncRNA.csv"
late_lnc_path  <- "late_lncRNA.csv"
sust_lnc_path  <- "sustained_lncRNA.csv"
early_novel_lnc_path <- "early_novel_lncRNA.csv"
late_novel_lnc_path  <- "late_novel_lncRNA.csv"
sust_novel_lnc_path  <- "sustained_novel_lncRNA.csv"

lnc_de6_path  <- "NE_vs_Ctrl_6h_lncRNA_DE.txt"
lnc_de24_path <- "NE_vs_Ctrl_24h_lncRNA_DE.txt"

cor_summary_path <- "correlation_summary_counts.csv"
cor_early_path   <- "cor_early_annotated.csv"
cor_late_path    <- "cor_late_annotated.csv"

cor_gained_path    <- "correlations_gained_at_24h.csv"
cor_lost_path      <- "correlations_lost_after_6h.csv"
cor_late_only_path <- "correlations_late_only.csv"

kegg_only6_path  <- "KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.txt"
kegg_only24_path <- "KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt"
kegg_shared_path <- "KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt"

out_dir <- "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05
top_lnc_per_phase_heatmap <- 20
top_hubs_plot <- 20
topN_enrich <- 10
topN_reversal_shared <- 5
dpi_png <- 600

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

read_lnc_set <- function(path, phase_label, source_label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      phase = phase_label,
      source = source_label,
      gene_id = as.character(gene_id),
      transcript_id = as.character(transcript_id),
      gene_key = toupper(gene_id),
      tx_key = toupper(transcript_id)
    )
}

read_lnc_de <- function(path, time_label, alpha = 0.05) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      gene_key = toupper(as.character(gene_id)),
      gene_name_raw = as.character(gene_name),
      gene_name = ifelse(is.na(gene_name_raw) | gene_name_raw == "", gene_key, gene_name_raw),
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      direction = as.character(direction),
      time = time_label
    ) %>%
    group_by(gene_key, time) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(de_sig = !is.na(padj) & padj < alpha)
}

standardize_cor <- function(df, label) {
  nms <- names(df)
  
  lnc_col <- if ("lncRNA" %in% nms) "lncRNA" else nms[str_detect(nms, regex("^lnc", TRUE))][1]
  lncname_col <- if ("lnc_gene_name" %in% nms) "lnc_gene_name" else nms[str_detect(nms, regex("lnc.*gene.*name", TRUE))][1]
  src_col <- if ("lnc_source" %in% nms) "lnc_source" else nms[str_detect(nms, regex("lnc.*source", TRUE))][1]
  mrna_col <- if ("mRNA_gene_name" %in% nms) "mRNA_gene_name" else nms[str_detect(nms, regex("mrna.*gene.*name|^mRNA$|target|symbol", TRUE))][1]
  rho_col <- if ("rho" %in% nms) "rho" else nms[str_detect(nms, regex("^rho$|spearman|pearson|corr", TRUE))][1]
  
  df %>%
    transmute(
      lncRNA = if (!is.na(lnc_col)) as.character(.data[[lnc_col]]) else NA_character_,
      lnc_gene_name = if (!is.na(lncname_col)) as.character(.data[[lncname_col]]) else NA_character_,
      lnc_source = if (!is.na(src_col)) as.character(.data[[src_col]]) else "unknown",
      mRNA_gene_name = if (!is.na(mrna_col)) as.character(.data[[mrna_col]]) else NA_character_,
      rho = if (!is.na(rho_col)) as.numeric(.data[[rho_col]]) else NA_real_,
      set = label
    ) %>%
    mutate(
      lnc_source = as.character(lnc_source),
      lnc_source = ifelse(is.na(lnc_source) | lnc_source == "" | lnc_source == "NA", "unknown", lnc_source),
      lnc_key = toupper(lncRNA),
      mRNA_key = toupper(mRNA_gene_name)
    )
}

hub_from_cor <- function(df) {
  df %>%
    filter(!is.na(lncRNA), !is.na(mRNA_key)) %>%
    group_by(set, lncRNA, lnc_gene_name, lnc_source) %>%
    summarize(
      n_targets = n_distinct(mRNA_key),
      n_pos = n_distinct(mRNA_key[!is.na(rho) & rho > 0]),
      n_neg = n_distinct(mRNA_key[!is.na(rho) & rho < 0]),
      .groups = "drop"
    ) %>%
    arrange(set, desc(n_targets))
}

read_edge_set <- function(path, label) {
  df <- read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0) {
    message("WARN: ", label, " está vacío (0 filas): ", path)
    return(tibble(
      lncRNA = character(), lnc_gene_name = character(), lnc_source = character(),
      mRNA_gene_name = character(), rho = numeric(), set = label,
      lnc_key = character(), mRNA_key = character()
    ))
  }
  standardize_cor(df, label) %>% filter(!is.na(lncRNA), !is.na(mRNA_key))
}

lnc_all <- bind_rows(
  read_lnc_set(early_lnc_path, "Early", "annotated"),
  read_lnc_set(sust_lnc_path,  "Sustained", "annotated"),
  read_lnc_set(late_lnc_path,  "Late", "annotated"),
  read_lnc_set(early_novel_lnc_path, "Early", "novel"),
  read_lnc_set(sust_novel_lnc_path,  "Sustained", "novel"),
  read_lnc_set(late_novel_lnc_path,  "Late", "novel")
) %>%
  mutate(
    phase = factor(phase, levels = c("Early", "Sustained", "Late")),
    source = factor(source, levels = c("annotated", "novel"))
  )

lnc_phase_map <- lnc_all %>%
  transmute(gene_key, phase, lnc_source = source) %>%
  distinct(gene_key, .keep_all = TRUE)

lnc_counts <- lnc_all %>%
  group_by(phase, source) %>%
  summarize(
    n_genes = n_distinct(gene_key),
    n_transcripts = n_distinct(tx_key),
    .groups = "drop"
  )
write_csv(lnc_counts, file.path(out_dir, "SUPP_Fig4B1_lnc_counts_phase_source.csv"))

p_Fig4B1 <- ggplot(lnc_counts, aes(x = phase, y = n_genes, fill = source)) +
  geom_col(position = "dodge") +
  labs(x = NULL, y = "# lncRNA genes", title = "lncRNA temporal programs: annotated vs novel") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_gg(p_Fig4B1, "Fig4B1_lnc_gene_counts", 7.5, 4.5, dpi_png)

lnc_de_all <- bind_rows(
  read_lnc_de(lnc_de6_path,  "6h", alpha),
  read_lnc_de(lnc_de24_path, "24h", alpha)
) %>% mutate(time = factor(time, levels = c("6h", "24h")))

lnc_de_phase <- lnc_de_all %>%
  inner_join(lnc_phase_map, by = "gene_key") %>%
  mutate(
    direction = ifelse(is.na(direction) | direction == "", "NA", direction),
    phase = factor(phase, levels = c("Early","Sustained","Late")),
    lnc_source = factor(lnc_source, levels = c("annotated","novel"))
  )

write_csv(lnc_de_phase, file.path(out_dir, "SUPP_Fig4B2_lncDE_with_phase.csv"))

lnc_de_counts <- lnc_de_phase %>%
  filter(de_sig) %>%
  group_by(time, phase, lnc_source, direction) %>%
  summarize(n = n_distinct(gene_key), .groups = "drop")

write_csv(lnc_de_counts, file.path(out_dir, "SUPP_Fig4B2_lncDE_counts_phase_time_direction.csv"))

p_Fig4B2 <- ggplot(lnc_de_counts, aes(x = phase, y = n, fill = direction)) +
  geom_col(position = "stack") +
  facet_grid(time ~ lnc_source) +
  labs(x = NULL, y = "# DE lncRNA genes", title = "DE lncRNAs across temporal programs (6h vs 24h)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_gg(p_Fig4B2, "Fig4B2_lncDE_counts_phase_time_direction", 10.5, 5.2, dpi_png)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pheatmap))
  
  lnc_de_wide <- lnc_de_phase %>%
    select(gene_key, gene_name, phase, lnc_source, time, log2FC, padj) %>%
    pivot_wider(names_from = time, values_from = c(log2FC, padj), names_sep = "_")
  
  lnc_top <- lnc_de_wide %>%
    mutate(padj_min = pmin(padj_6h, padj_24h, na.rm = TRUE)) %>%
    group_by(phase) %>%
    arrange(padj_min, desc(pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE))) %>%
    slice_head(n = top_lnc_per_phase_heatmap) %>%
    ungroup()
  
  mat_df <- lnc_top %>%
    select(gene_name, phase, lnc_source, log2FC_6h, log2FC_24h) %>%
    mutate(absmax = pmax(abs(log2FC_6h), abs(log2FC_24h), na.rm = TRUE)) %>%
    arrange(phase, desc(absmax))
  
  mat <- mat_df %>%
    select(gene_name, log2FC_6h, log2FC_24h) %>%
    tibble::column_to_rownames("gene_name") %>% as.matrix()
  
  ann <- mat_df %>%
    select(gene_name, phase, lnc_source) %>%
    tibble::column_to_rownames("gene_name")
  
  keep <- apply(mat, 1, function(x) all(is.finite(x)))
  mat <- mat[keep, , drop = FALSE]
  ann <- ann[rownames(mat), , drop = FALSE]
  
  grDevices::pdf(file.path(out_dir, "Fig4B3_heatmap_lnc_log2FC_top_by_phase.pdf"), width = 6.8, height = 10)
  pheatmap::pheatmap(
    mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
    main = "Top lncRNAs by temporal program (log2FC: NE vs Ctrl)", fontsize_row = 7
  )
  dev.off()
  
  grDevices::png(file.path(out_dir, "Fig4B3_heatmap_lnc_log2FC_top_by_phase_600dpi.png"),
                 width = 6.8, height = 10, units = "in", res = dpi_png)
  pheatmap::pheatmap(
    mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann,
    main = "Top lncRNAs by temporal program (log2FC: NE vs Ctrl)", fontsize_row = 7
  )
  dev.off()
}

cor_summary <- read_csv(cor_summary_path, show_col_types = FALSE) %>%
  filter(n > 0) %>%
  filter(!tolower(category) %in% c("stable", "sustained", "stable_6h_24h", "correlations_stable_6h_24h"))

p_Fig4C1 <- ggplot(cor_summary, aes(x = category, y = n)) +
  geom_col() +
  labs(x = NULL, y = "# correlations", title = "lncRNA–mRNA correlation summary (rewiring)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
save_gg(p_Fig4C1, "Fig4C1_correlation_summary_counts_NO_STABLE", 7.5, 4.5, dpi_png)

cor_early <- standardize_cor(read_csv(cor_early_path, show_col_types = FALSE), "Early")
cor_late  <- standardize_cor(read_csv(cor_late_path,  show_col_types = FALSE), "Late")
hubs_early_late <- hub_from_cor(bind_rows(cor_early, cor_late))
write_csv(hubs_early_late, file.path(out_dir, "SUPP_Fig4C2_hubs_early_vs_late.csv"))

hubs_plot_df <- hubs_early_late %>%
  group_by(set) %>% slice_head(n = top_hubs_plot) %>% ungroup() %>%
  mutate(
    label = ifelse(!is.na(lnc_gene_name) & lnc_gene_name != "", lnc_gene_name, lncRNA),
    label = fct_reorder(label, n_targets, .desc = TRUE),
    set = factor(set, levels = c("Early", "Late")),
    lnc_source = factor(lnc_source, levels = c("novel", "annotated", "unknown"))
  )

p_Fig4C2 <- ggplot(hubs_plot_df, aes(x = n_targets, y = label, fill = lnc_source)) +
  geom_col() +
  facet_wrap(~ set, scales = "free_y") +
  labs(x = "# mRNA targets (degree)", y = NULL,
       title = "Top lncRNA hubs: commitment (Early) vs maintenance (Late)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_gg(p_Fig4C2, "Fig4C2_hubs_early_vs_late", 11, 6, dpi_png)

cor_gained   <- read_edge_set(cor_gained_path, "Gained_at_24h")
cor_lost     <- read_edge_set(cor_lost_path,   "Lost_after_6h")
cor_lateonly <- read_edge_set(cor_late_only_path, "Late_only")

edges_list <- list(cor_lost, cor_gained, cor_lateonly)
edges_list <- edges_list[vapply(edges_list, nrow, integer(1)) > 0]
if (length(edges_list) == 0) stop("Rewiring edges vacíos: revisa gained/lost/late-only.")

rewiring_hubs <- hub_from_cor(bind_rows(edges_list)) %>%
  mutate(set = factor(set, levels = c("Lost_after_6h", "Gained_at_24h", "Late_only")))
write_csv(rewiring_hubs, file.path(out_dir, "SUPP_Fig4C3_hubs_by_rewiring_category_NO_STABLE.csv"))

rewiring_top <- rewiring_hubs %>%
  group_by(set) %>% slice_head(n = 15) %>% ungroup() %>%
  mutate(
    label = ifelse(!is.na(lnc_gene_name) & lnc_gene_name != "", lnc_gene_name, lncRNA),
    label = fct_reorder(label, n_targets, .desc = TRUE),
    lnc_source = factor(lnc_source, levels = c("novel", "annotated", "unknown"))
  )

p_Fig4C3 <- ggplot(rewiring_top, aes(x = n_targets, y = label, fill = lnc_source)) +
  geom_col() +
  facet_wrap(~ set, scales = "free_y") +
  labs(x = "# mRNA targets (degree)", y = NULL,
       title = "lncRNA hubs by rewiring category (lost / gained / late-only)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_gg(p_Fig4C3, "Fig4C3_hubs_by_rewiring_category_NO_STABLE", 12.5, 6.5, dpi_png)

kegg_only6 <- read_tsv(kegg_only6_path, show_col_types = FALSE) %>%
  mutate(KEGG_group = "Only6h") %>%
  transmute(
    ID = as.character(ID),
    Description = as.character(Description),
    NES_6h = as.numeric(NES),
    padj_6h = as.numeric(`p.adjust`),
    NES_24h = NA_real_,
    padj_24h = NA_real_,
    gene_symbol = as.character(gene_symbol),
    KEGG_group = KEGG_group
  )

kegg_only24 <- read_tsv(kegg_only24_path, show_col_types = FALSE) %>%
  mutate(KEGG_group = "Only24h") %>%
  transmute(
    ID = as.character(ID),
    Description = as.character(Description),
    NES_6h = NA_real_,
    padj_6h = NA_real_,
    NES_24h = as.numeric(NES),
    padj_24h = as.numeric(`p.adjust`),
    gene_symbol = as.character(gene_symbol),
    KEGG_group = KEGG_group
  )

kegg_shared <- read_tsv(kegg_shared_path, show_col_types = FALSE) %>%
  mutate(KEGG_group = "Shared") %>%
  transmute(
    ID = as.character(ID),
    Description = as.character(Description),
    NES_6h = as.numeric(NES_6h),
    NES_24h = as.numeric(NES_24h),
    padj_6h = as.numeric(`p.adjust_6h`),
    padj_24h = as.numeric(`p.adjust_24h`),
    gene_symbol = as.character(gene_symbol),
    KEGG_group = KEGG_group
  )

kegg_stats <- bind_rows(kegg_only6, kegg_shared, kegg_only24) %>%
  distinct(KEGG_group, ID, Description, NES_6h, NES_24h, padj_6h, padj_24h) %>%
  mutate(
    KEGG_group = factor(KEGG_group, levels = c("Only6h", "Shared", "Only24h")),
    padj_min = pmin(padj_6h, padj_24h, na.rm = TRUE),
    transition_shared = case_when(
      KEGG_group != "Shared" ~ NA_character_,
      NES_6h > 0 & NES_24h > 0 ~ "Maintained_Activated",
      NES_6h < 0 & NES_24h < 0 ~ "Maintained_Suppressed",
      NES_6h * NES_24h < 0     ~ "Reversal",
      TRUE ~ "Other"
    )
  )

write_csv(kegg_stats, file.path(out_dir, "SUPP_Fig5_full_KEGG_stats_wide.csv"))

kegg_long <- kegg_stats %>%
  select(KEGG_group, ID, Description, transition_shared, NES_6h, NES_24h, padj_6h, padj_24h) %>%
  pivot_longer(
    cols = c(NES_6h, NES_24h, padj_6h, padj_24h),
    names_to = c(".value", "time"),
    names_pattern = "(NES|padj)_(6h|24h)"
  ) %>%
  filter(!is.na(NES), !is.na(padj)) %>%
  mutate(
    time = factor(time, levels = c("6h", "24h")),
    sign_class = case_when(
      KEGG_group == "Shared" & transition_shared == "Reversal" ~ "Reversal",
      NES > 0 ~ "Activated (NES>0)",
      NES < 0 ~ "Suppressed (NES<0)",
      TRUE ~ "Other"
    ),
    sign_class = factor(sign_class, levels = c("Activated (NES>0)", "Suppressed (NES<0)", "Reversal", "Other"))
  )

sel_only6_act <- kegg_stats %>%
  filter(KEGG_group == "Only6h", !is.na(NES_6h), !is.na(padj_6h), padj_6h < alpha, NES_6h > 0) %>%
  arrange(padj_6h) %>% slice_head(n = topN_enrich) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Activated (NES>0)", transition_shared = NA_character_)

sel_only6_sup <- kegg_stats %>%
  filter(KEGG_group == "Only6h", !is.na(NES_6h), !is.na(padj_6h), padj_6h < alpha, NES_6h < 0) %>%
  arrange(padj_6h) %>% slice_head(n = topN_enrich) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Suppressed (NES<0)", transition_shared = NA_character_)

sel_only24_act <- kegg_stats %>%
  filter(KEGG_group == "Only24h", !is.na(NES_24h), !is.na(padj_24h), padj_24h < alpha, NES_24h > 0) %>%
  arrange(padj_24h) %>% slice_head(n = topN_enrich) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Activated (NES>0)", transition_shared = NA_character_)

sel_only24_sup <- kegg_stats %>%
  filter(KEGG_group == "Only24h", !is.na(NES_24h), !is.na(padj_24h), padj_24h < alpha, NES_24h < 0) %>%
  arrange(padj_24h) %>% slice_head(n = topN_enrich) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Suppressed (NES<0)", transition_shared = NA_character_)

sel_shared_maint_act <- kegg_stats %>%
  filter(KEGG_group == "Shared", padj_min < alpha, transition_shared == "Maintained_Activated") %>%
  arrange(padj_min) %>% slice_head(n = topN_enrich) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Activated (NES>0)", transition_shared)

sel_shared_maint_sup <- kegg_stats %>%
  filter(KEGG_group == "Shared", padj_min < alpha, transition_shared == "Maintained_Suppressed") %>%
  arrange(padj_min) %>% slice_head(n = topN_enrich) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Suppressed (NES<0)", transition_shared)

sel_shared_rev <- kegg_stats %>%
  filter(KEGG_group == "Shared", padj_min < alpha, transition_shared == "Reversal") %>%
  arrange(padj_min) %>% slice_head(n = topN_reversal_shared) %>%
  transmute(KEGG_group, ID, Description, sign_class = "Reversal", transition_shared)

pathways_selected <- bind_rows(
  sel_only6_act, sel_only6_sup,
  sel_shared_maint_act, sel_shared_maint_sup, sel_shared_rev,
  sel_only24_act, sel_only24_sup
) %>%
  mutate(sign_class = factor(sign_class, levels = c("Activated (NES>0)", "Suppressed (NES<0)", "Reversal"))) %>%
  distinct(KEGG_group, ID, Description, sign_class, .keep_all = TRUE)

write_csv(pathways_selected, file.path(out_dir, "pathways_selected_for_TF.csv"))

dot_df <- kegg_long %>%
  inner_join(pathways_selected, by = c("KEGG_group", "ID", "Description", "sign_class")) %>%
  mutate(Description = fct_reorder(Description, NES, .desc = TRUE))

p_Fig5A <- ggplot(dot_df, aes(x = time, y = Description, size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.9) +
  facet_grid(sign_class ~ KEGG_group, scales = "free_y", space = "free_y") +
  labs(
    x = NULL, y = NULL,
    title = "KEGG GSEA across commitment vs maintenance (Only6h / Shared / Only24h)",
    size = "-log10(FDR)"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9))
save_gg(p_Fig5A, "Fig5A_KEGG_GSEA_dotplot_TopNESpos_neg_by_group", 12.8, 10, dpi_png)

shared_df <- kegg_stats %>%
  filter(KEGG_group == "Shared", padj_min < alpha) %>%
  mutate(transition_shared = factor(
    transition_shared,
    levels = c("Maintained_Activated", "Maintained_Suppressed", "Reversal", "Other")
  ))

p_Fig5B <- ggplot(shared_df, aes(x = NES_6h, y = NES_24h, size = -log10(padj_min), color = transition_shared)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    x = "NES (6h)", y = "NES (24h)",
    title = "Shared KEGG pathways: maintained vs reversal between 6h and 24h",
    size = "-log10(min FDR)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
save_gg(p_Fig5B, "Fig5B_Shared_scatter_NES6h_vs_NES24h", 7.6, 6, dpi_png)

write_csv(dot_df, file.path(out_dir, "SUPP_Fig5_selected_pathways_dotplot_data.csv"))

message("OK. Outputs: ", out_dir)
message("Export TF pathways: ", file.path(out_dir, "pathways_selected_for_TF.csv"))

# Fig6 TF panels with DE + phase filtering:
#  - Keep only target genes that are:
#      (i) DE at the relevant timepoint (6h or 24h; padj < alpha)
#      (ii) classified as Early / Sustained / Late
#      (iii) belong to selected KEGG pathways (Only6h / Shared / Only24h)
#
# Outputs:
#   - Fig6A_TF_KEGG_dotplot_commitment_vs_maintenance.pdf
#   - Fig6B_TF_ranking_commitment_vs_maintenance.pdf
#   - Supplementary CSVs with annotated/filtered edges and summary tables

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

# Paths (uploaded files)
edge_6h_path   <- "GRN_NE06_NORA_only_6h_GO_default_edge.csv"
edge_24h_path  <- "GRN_NE24_NORA_only_24h_go_default_edge.csv"
edge_merged_path <- "Merged_Network_6_24_default_edge.csv"

kegg_only6_path <- "KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.txt"
kegg_only24_path <- "KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt"
kegg_shared_path <- "KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt"

de6_path  <- "NE_vs_Ctrl_6h_mRNA_DE.txt"
de24_path <- "NE_vs_Ctrl_24h_mRNA_DE.txt"

early_path <- "early_mRNA.csv"
late_path  <- "late_mRNA.csv"
sustained_path <- "sustained_mRNA.csv"

out_dir <- "TF_panels_out_DE_phase"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters
alpha <- 0.05                 # DE threshold
top_n_tfs_for_dotplot <- 20   # TFs shown per context
top_n_tfs_ranking <- 15       # TFs shown per context in ranking

# Recommended: focus on a curated set of pathways for main figure clarity
pathways_focus <- c(
  "Focal adhesion",
  "Cytokine-cytokine receptor interaction",
  "Hippo signaling pathway",
  "Apoptosis",
  "Regulation of actin cytoskeleton",
  "PI3K-Akt signaling pathway",
  "Cell cycle",
  "Protein processing in endoplasmic reticulum",
  "TGF-beta signaling pathway",
  "Signaling pathways regulating pluripotency of stem cells",
  "Peroxisome",
  "Fatty acid degradation",
  "Adipocytokine signaling pathway",
  "Calcium signaling pathway",
  "Valine, leucine and isoleucine degradation",
  "PPAR signaling pathway"
)

# If you want automatic selection instead, set:
# pathways_focus <- NULL
top_n_pathways_auto <- 12

# Helpers
read_edges <- function(path, group_label) {
  df <- read_csv(path, show_col_types = FALSE)
  name_col <- if ("shared name" %in% names(df)) "shared name" else if ("name" %in% names(df)) "name" else NA_character_
  if (is.na(name_col)) stop("Edge table missing 'shared name'/'name' column: ", path)
  
  df %>%
    mutate(edge_group = group_label,
           edge_name = .data[[name_col]]) %>%
    mutate(
      TF_raw = str_split_fixed(edge_name, " \\(interacts with\\) ", 2)[, 1],
      Target_raw = str_split_fixed(edge_name, " \\(interacts with\\) ", 2)[, 2]
    ) %>%
    mutate(
      TF = str_trim(TF_raw),
      Target = str_trim(Target_raw),
      TF_key = str_to_upper(TF),
      Target_key = str_to_upper(Target)
    )
}

read_kegg_only <- function(path, group_label) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      kegg_group = group_label,
      ID, Description,
      NES = as.numeric(NES),
      padj = as.numeric(`p.adjust`),
      gene_symbol,
      gene_key = str_to_upper(gene_symbol)
    )
}

read_kegg_shared <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      kegg_group = "Shared",
      ID, Description,
      NES_6h = as.numeric(NES_6h),
      NES_24h = as.numeric(NES_24h),
      padj_6h = as.numeric(`p.adjust_6h`),
      padj_24h = as.numeric(`p.adjust_24h`),
      gene_symbol,
      gene_key = str_to_upper(gene_symbol)
    )
}

read_phase_sets <- function(early_path, late_path, sustained_path) {
  early <- read_csv(early_path, show_col_types = FALSE) %>%
    transmute(gene_key = str_to_upper(gene_name), phase = "Early")
  late <- read_csv(late_path, show_col_types = FALSE) %>%
    transmute(gene_key = str_to_upper(gene_name), phase = "Late")
  sust <- read_csv(sustained_path, show_col_types = FALSE) %>%
    transmute(gene_key = str_to_upper(gene_name), phase = "Sustained")
  
  # If a gene appears in more than one set (shouldn't, but just in case),
  # we enforce priority: Sustained > Early > Late
  bind_rows(sust, early, late) %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_de <- function(path, time_label) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      gene_key = str_to_upper(gene_name),
      gene_name,
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      direction = direction,
      time = time_label,
      de_sig = !is.na(padj) & padj < alpha
    )
}

pick_top_pathways <- function(kegg_tbl, group_label, n = 12) {
  kegg_tbl %>%
    filter(kegg_group == group_label) %>%
    group_by(Description) %>%
    summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(padj_min) %>%
    slice_head(n = n) %>%
    pull(Description)
}

# Load inputs
edges_only6  <- read_edges(edge_6h_path,  "Only6h_edges")
edges_only24 <- read_edges(edge_24h_path, "Only24h_edges")
edges_merged <- read_edges(edge_merged_path, "Merged_edges")
edges_all <- bind_rows(edges_only6, edges_only24, edges_merged)

kegg_only6  <- read_kegg_only(kegg_only6_path,  "Only6h")
kegg_only24 <- read_kegg_only(kegg_only24_path, "Only24h")
kegg_shared <- read_kegg_shared(kegg_shared_path)

kegg_all_genes <- bind_rows(
  kegg_only6 %>% select(kegg_group, ID, Description, gene_key, NES, padj),
  kegg_only24 %>% select(kegg_group, ID, Description, gene_key, NES, padj),
  kegg_shared %>% select(kegg_group, ID, Description, gene_key) %>% mutate(NES = NA_real_, padj = NA_real_)
)

phase_map <- read_phase_sets(early_path, late_path, sustained_path)

de6  <- read_de(de6_path, "6h")
de24 <- read_de(de24_path, "24h")
de_all <- bind_rows(de6, de24)

de_all_unique <- de_all %>%
  group_by(gene_key, time) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice(1) %>%
  ungroup()

de_wide <- de_all_unique %>%
  select(gene_key, time, log2FC, padj, direction, de_sig) %>%
  pivot_wider(
    names_from = time,
    values_from = c(log2FC, padj, direction, de_sig),
    names_sep = "_"
  )

dups <- de_all %>%
  count(gene_key, time) %>%
  filter(n > 1)

write_csv(dups, file.path(out_dir, "QC_DE_duplicates_by_gene_and_time.csv"))

de_wide <- de_wide %>%
  mutate(
    de_sig_6h = ifelse(is.na(de_sig_6h), FALSE, de_sig_6h),
    de_sig_24h = ifelse(is.na(de_sig_24h), FALSE, de_sig_24h)
  )


# Choose KEGG pathways (focus or auto)
if (is.null(pathways_focus)) {
  top6  <- pick_top_pathways(kegg_only6,  "Only6h",  n = top_n_pathways_auto)
  top24 <- pick_top_pathways(kegg_only24, "Only24h", n = top_n_pathways_auto)
  top_shared <- kegg_shared %>%
    group_by(Description) %>%
    summarize(padj_min = min(padj_6h, padj_24h, na.rm = TRUE), .groups = "drop") %>%
    arrange(padj_min) %>%
    slice_head(n = top_n_pathways_auto) %>%
    pull(Description)
  
  pathways_focus <- unique(c(top6, top_shared, top24))
}

# Build Target->(KEGG_group, pathway) mapping for focused pathways
target2kegg <- kegg_all_genes %>%
  filter(Description %in% pathways_focus) %>%
  distinct(kegg_group, Description, gene_key)

# Annotate edges with KEGG membership, phase, and DE stats

edges_kegg <- edges_all %>%
  inner_join(target2kegg,
             by = c("Target_key" = "gene_key"),
             relationship = "many-to-many") %>%
  rename(KEGG_group = kegg_group) %>%
  mutate(KEGG_group = factor(KEGG_group, levels = c("Only6h", "Shared", "Only24h"))) %>%
  left_join(phase_map, by = c("Target_key" = "gene_key")) %>%
  mutate(phase = replace_na(phase, "Unclassified")) %>%
  left_join(de_wide, by = c("Target_key" = "gene_key"))

# Keep only genes in a time category (Early/Sustained/Late)
edges_kegg <- edges_kegg %>%
  filter(phase %in% c("Early", "Sustained", "Late"))

# Export: annotated (pre time-splitting) for supplements
write_csv(edges_kegg, file.path(out_dir, "SUPP_edges_TF_target_annotated_KEGG_phase_DE_raw.csv"))

# Create a time-context variable for commitment vs maintenance
# Rules:
#  - Only6h pathways -> evaluate DE at 6h (Commitment)
#  - Only24h pathways -> evaluate DE at 24h (Maintenance)
#  - Shared pathways -> duplicate into both contexts:
#       * Commitment (6h): keep only targets DE at 6h
#       * Maintenance (24h): keep only targets DE at 24h
edges_context <- edges_kegg %>%
  mutate(
    context = case_when(
      KEGG_group == "Only6h"  ~ "Commitment_6h",
      KEGG_group == "Only24h" ~ "Maintenance_24h",
      KEGG_group == "Shared"  ~ "Shared"
    )
  )

# Split Shared into two contexts
edges_shared_6h <- edges_context %>%
  filter(KEGG_group == "Shared") %>%
  mutate(context = "Commitment_6h")

edges_shared_24h <- edges_context %>%
  filter(KEGG_group == "Shared") %>%
  mutate(context = "Maintenance_24h")

edges_nonshared <- edges_context %>%
  filter(KEGG_group %in% c("Only6h", "Only24h"))

edges_context2 <- bind_rows(edges_nonshared, edges_shared_6h, edges_shared_24h) %>%
  mutate(context = factor(context, levels = c("Commitment_6h", "Maintenance_24h")))

# Apply DE filtering by context timepoint
edges_context2 <- edges_context2 %>%
  filter(
    (context == "Commitment_6h"   & de_sig_6h == TRUE) |
      (context == "Maintenance_24h" & de_sig_24h == TRUE)
  )

# Export the final filtered edge list (this is your key supplement)
write_csv(edges_context2, file.path(out_dir, "SUPP_edges_TF_target_FINAL_filtered_DE_and_phase.csv"))

# Build TF x KEGG counts within each context
tf_kegg_counts <- edges_context2 %>%
  group_by(context, KEGG_group, Description, TF) %>%
  summarize(n_targets = n_distinct(Target_key), .groups = "drop")

write_csv(tf_kegg_counts, file.path(out_dir, "SUPP_TF_by_KEGG_counts_context_DE_phase.csv"))

# Panel A: Dotplot TF x KEGG (DE+phase filtered)
# Strategy:
#   - choose top TFs per context based on total targets across pathways
#   - plot TF (y) vs pathway (x), size = #targets
#   - facet by context; optionally also facet by KEGG_group (comment/uncomment)
top_tfs <- tf_kegg_counts %>%
  group_by(context, TF) %>%
  summarize(total_targets = sum(n_targets), .groups = "drop") %>%
  arrange(context, desc(total_targets)) %>%
  group_by(context) %>%
  slice_head(n = top_n_tfs_for_dotplot) %>%
  ungroup()

tf_kegg_plot <- tf_kegg_counts %>%
  inner_join(top_tfs %>% select(context, TF), by = c("context", "TF")) %>%
  mutate(
    TF = fct_reorder(TF, n_targets, .fun = sum, .desc = TRUE),
    Description = fct_inorder(Description)
  )

# Pathway order: prioritize significance within each KEGG_group
path_order <- bind_rows(
  kegg_only6 %>% filter(Description %in% pathways_focus) %>%
    group_by(Description) %>% summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
    mutate(KEGG_group = "Only6h"),
  kegg_only24 %>% filter(Description %in% pathways_focus) %>%
    group_by(Description) %>% summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
    mutate(KEGG_group = "Only24h"),
  kegg_shared %>% filter(Description %in% pathways_focus) %>%
    group_by(Description) %>% summarize(padj_min = min(padj_6h, padj_24h, na.rm = TRUE), .groups = "drop") %>%
    mutate(KEGG_group = "Shared")
)

tf_kegg_plot <- tf_kegg_plot %>%
  left_join(path_order, by = c("KEGG_group", "Description")) %>%
  mutate(Description = fct_reorder(Description, padj_min, .desc = FALSE))

p_dot <- ggplot(tf_kegg_plot, aes(x = Description, y = TF)) +
  geom_point(aes(size = n_targets), alpha = 0.85) +
  facet_wrap(~ context, scales = "free_x") +
  # If you want more granularity, replace the line above with:
  # facet_grid(context ~ KEGG_group, scales = "free_x") +
  scale_size_continuous(range = c(1.5, 8)) +
  labs(
    x = "KEGG pathway (genes DE + phase-filtered)",
    y = "Transcription factor (TF)",
    size = "# DE targets\nin pathway",
    title = "Commitment vs maintenance: TF coverage of phase-specific KEGG programs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

ggsave(file.path(out_dir, "Fig6A_TF_KEGG_dotplot_commitment_vs_maintenance.pdf"),
       p_dot, width = 13, height = 7)

# Panel B: TF ranking per context
tf_rank <- tf_kegg_counts %>%
  group_by(context, TF) %>%
  summarize(
    total_targets = sum(n_targets),
    n_pathways_hit = sum(n_targets > 0),
    .groups = "drop"
  ) %>%
  arrange(context, desc(total_targets))

write_csv(tf_rank, file.path(out_dir, "SUPP_TF_ranking_context_DE_phase_full.csv"))

tf_rank_plot <- tf_rank %>%
  group_by(context) %>%
  slice_head(n = top_n_tfs_ranking) %>%
  ungroup() %>%
  mutate(TF = fct_reorder(TF, total_targets, .desc = TRUE))

p_rank <- ggplot(tf_rank_plot, aes(x = total_targets, y = TF)) +
  geom_col() +
  facet_wrap(~ context, scales = "free_y") +
  labs(
    x = "Total # DE targets across selected KEGG pathways",
    y = "TF",
    title = "Top TFs governing commitment (6h) vs maintenance (24h) programs"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Fig6B_TF_ranking_commitment_vs_maintenance.pdf"),
       p_rank, width = 11, height = 6)

# Supplement: KEGG enrichment figure suggestion (single integrated plot)
# (Keeps enrichment separate from TF panels)
kegg_pathway_stats <- bind_rows(
  kegg_only6 %>%
    group_by(Description, ID) %>% summarize(NES = first(NES), padj = first(padj), .groups = "drop") %>%
    mutate(KEGG_group = "Only6h", time = "6h"),
  kegg_only24 %>%
    group_by(Description, ID) %>% summarize(NES = first(NES), padj = first(padj), .groups = "drop") %>%
    mutate(KEGG_group = "Only24h", time = "24h"),
  kegg_shared %>%
    group_by(Description, ID) %>%
    summarize(NES_6h = first(NES_6h), NES_24h = first(NES_24h),
              padj_6h = first(padj_6h), padj_24h = first(padj_24h), .groups = "drop") %>%
    pivot_longer(cols = c(NES_6h, NES_24h, padj_6h, padj_24h),
                 names_to = c(".value", "time"),
                 names_pattern = "(NES|padj)_(6h|24h)") %>%
    mutate(KEGG_group = "Shared")
) %>%
  filter(Description %in% pathways_focus) %>%
  mutate(
    time = factor(time, levels = c("6h", "24h")),
    KEGG_group = factor(KEGG_group, levels = c("Only6h", "Shared", "Only24h")),
    neglog10 = -log10(padj)
  )

# Dotplot: size = -log10(FDR); facet = KEGG_group; x = time
p_kegg <- ggplot(kegg_pathway_stats,
                 aes(x = time, y = fct_reorder(Description, neglog10, .fun = max, .desc = TRUE))) +
  geom_point(aes(size = neglog10), alpha = 0.85) +
  facet_wrap(~ KEGG_group, scales = "free_y") +
  scale_size_continuous(range = c(1.5, 8)) +
  labs(
    x = "Time point",
    y = "KEGG pathway",
    size = "-log10(FDR)",
    title = "KEGG enrichment across commitment (6h), shared, and maintenance (24h)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Fig5_KEGG_dotplot_unified.pdf"),
       p_kegg, width = 12, height = 8)

write_csv(kegg_pathway_stats, file.path(out_dir, "SUPP_KEGG_pathway_level_stats_for_plot.csv"))

# Supplement / optional main: KEGG gene-set composition by Early/Sustained/Late
# (very useful to connect processes with your time-classification)
kegg_gene_phase <- target2kegg %>%
  left_join(phase_map, by = c("gene_key")) %>%
  mutate(phase = replace_na(phase, "Unclassified")) %>%
  filter(phase %in% c("Early", "Sustained", "Late")) %>%
  mutate(kegg_group = factor(kegg_group, levels = c("Only6h", "Shared", "Only24h"))) %>%
  filter(Description %in% pathways_focus)

kegg_phase_comp <- kegg_gene_phase %>%
  group_by(kegg_group, Description, phase) %>%
  summarize(n_genes = n_distinct(gene_key), .groups = "drop") %>%
  group_by(kegg_group, Description) %>%
  mutate(total = sum(n_genes), prop = n_genes / total) %>%
  ungroup()

write_csv(kegg_phase_comp, file.path(out_dir, "SUPP_KEGG_phase_composition_counts_props.csv"))

p_comp <- ggplot(kegg_phase_comp,
                 aes(x = fct_reorder(Description, total, .desc = TRUE), y = prop, fill = phase)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ kegg_group, scales = "free_y") +
  labs(
    x = "KEGG pathway",
    y = "Proportion of genes",
    title = "Composition of KEGG gene sets by Early/Sustained/Late classification"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "FigS_KEGG_phase_composition_stackedbars.pdf"),
       p_comp, width = 12, height = 8)

message("Done. Outputs in: ", out_dir)

#####Nueva sección de la figura 6
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

edge_6h_path     <- "GRN_NE06_NORA_only_6h_GO_default_edge.csv"
edge_24h_path    <- "GRN_NE24_NORA_only_24h_go_default_edge.csv"
edge_merged_path <- "Merged_Network_6_24_default_edge.csv"

kegg_only6_path  <- "KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.txt"
kegg_only24_path <- "KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt"
kegg_shared_path <- "KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt"

de6_path  <- "NE_vs_Ctrl_6h_mRNA_DE.txt"
de24_path <- "NE_vs_Ctrl_24h_mRNA_DE.txt"

early_path     <- "early_mRNA.csv"
late_path      <- "late_mRNA.csv"
sustained_path <- "sustained_mRNA.csv"

# IMPORTANTE: este archivo debe venir del script Fig5
pathways_selected_path <- file.path("Paper_Fig4_Fig5_FINAL_ONE_SCRIPT", "pathways_selected_for_TF.csv")

out_dir <- "Fig6_TF_panels_out_FINAL"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05
top_n_tfs_for_dotplot <- 20
top_n_tfs_ranking <- 15
dpi_png <- 600

pdf_device <- function(...) {
  if (capabilities("cairo")) grDevices::cairo_pdf(...) else grDevices::pdf(...)
}

save_gg <- function(p, base_name, w, h, dpi = 600) {
  ggsave(file.path(out_dir, paste0(base_name, ".pdf")),
         plot = p, width = w, height = h, units = "in", device = pdf_device)
  ggsave(file.path(out_dir, paste0(base_name, "_600dpi.png")),
         plot = p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
}

pick_first_existing <- function(nms, candidates_regex) {
  hit <- nms[str_detect(nms, regex(candidates_regex, ignore_case = TRUE))][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

parse_shared_name <- function(x) {
  x <- str_squish(as.character(x))
  m <- str_match(x, "^(.+?)\\s*\\([^)]*\\)\\s*(.+?)$")
  tf <- m[,2]
  tgt <- m[,3]
  
  # Fallback si no matchea (por si viene con separadores raros)
  bad <- is.na(tf) | is.na(tgt) | tf == "" | tgt == ""
  if (any(bad)) {
    # intenta separar por " ) " o " )" si hay paréntesis
    tf2 <- ifelse(bad, str_trim(str_extract(x, "^.+?(?=\\s*\\()")), tf)
    tgt2 <- ifelse(bad, str_trim(str_replace(x, "^.+?\\)\\s*", "")), tgt)
    tf <- tf2; tgt <- tgt2
  }
  
  tibble(TF = tf, Target = tgt)
}

read_edges <- function(path, group_label) {
  df <- read_csv(path, show_col_types = FALSE)
  nms <- names(df)
  
  src_col <- pick_first_existing(nms, "^(source|Source)$")
  tgt_col <- pick_first_existing(nms, "^(target|Target)$")
  
  if (!is.na(src_col) && !is.na(tgt_col)) {
    out <- df %>%
      transmute(
        edge_group = group_label,
        TF = as.character(.data[[src_col]]),
        Target = as.character(.data[[tgt_col]])
      )
  } else {
    name_col <- pick_first_existing(nms, "^shared name$|^name$|shared\\s*name|^interaction name$")
    if (is.na(name_col)) {
      stop("No pude detectar columnas Source/Target ni 'shared name'/'name' en: ", path,
           "\nColumnas disponibles:\n", paste(nms, collapse = ", "))
    }
    parsed <- parse_shared_name(df[[name_col]])
    out <- bind_cols(tibble(edge_group = group_label), parsed)
  }
  
  out %>%
    mutate(
      TF = str_trim(TF),
      Target = str_trim(Target),
      TF_key = str_to_upper(TF),
      Target_key = str_to_upper(Target)
    ) %>%
    filter(!is.na(TF_key), TF_key != "", !is.na(Target_key), Target_key != "")
}

read_kegg_only <- function(path, group_label) {
  df <- read_tsv(path, show_col_types = FALSE)
  nms <- names(df)
  if (!all(c("ID","Description") %in% nms)) stop("KEGG file sin ID/Description: ", path)
  
  pcol <- if ("p.adjust" %in% nms) "p.adjust" else pick_first_existing(nms, "p\\.adjust|padj|FDR")
  gcol <- if ("gene_symbol" %in% nms) "gene_symbol" else pick_first_existing(nms, "gene_symbol|symbol|gene")
  
  df %>%
    transmute(
      kegg_group = group_label,
      ID = as.character(ID),
      Description = as.character(Description),
      NES = as.numeric(NES),
      padj = as.numeric(.data[[pcol]]),
      gene_symbol = as.character(.data[[gcol]]),
      gene_key = str_to_upper(gene_symbol)
    ) %>%
    filter(!is.na(gene_key), gene_key != "")
}

read_kegg_shared <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  nms <- names(df)
  if (!all(c("ID","Description") %in% nms)) stop("KEGG shared file sin ID/Description: ", path)
  
  gcol <- if ("gene_symbol" %in% nms) "gene_symbol" else pick_first_existing(nms, "gene_symbol|symbol|gene")
  
  df %>%
    transmute(
      kegg_group = "Shared",
      ID = as.character(ID),
      Description = as.character(Description),
      NES_6h = as.numeric(NES_6h),
      NES_24h = as.numeric(NES_24h),
      padj_6h = as.numeric(`p.adjust_6h`),
      padj_24h = as.numeric(`p.adjust_24h`),
      gene_symbol = as.character(.data[[gcol]]),
      gene_key = str_to_upper(gene_symbol)
    ) %>%
    filter(!is.na(gene_key), gene_key != "")
}

read_phase_sets <- function(early_path, late_path, sustained_path) {
  read_phase_one <- function(path, phase_label) {
    df <- read_csv(path, show_col_types = FALSE)
    nms <- names(df)
    gcol <- if ("gene_name" %in% nms) "gene_name" else if ("gene" %in% nms) "gene" else nms[1]
    df %>%
      transmute(gene_key = str_to_upper(as.character(.data[[gcol]])), phase = phase_label) %>%
      filter(!is.na(gene_key), gene_key != "")
  }
  
  early <- read_phase_one(early_path, "Early")
  late  <- read_phase_one(late_path, "Late")
  sust  <- read_phase_one(sustained_path, "Sustained")
  
  bind_rows(sust, early, late) %>%
    distinct(gene_key, .keep_all = TRUE)
}

read_de <- function(path, time_label, alpha = 0.05) {
  df <- read_tsv(path, show_col_types = FALSE)
  nms <- names(df)
  gcol <- if ("gene_name" %in% nms) "gene_name" else if ("gene" %in% nms) "gene" else nms[1]
  df %>%
    transmute(
      gene_key = str_to_upper(as.character(.data[[gcol]])),
      log2FC = as.numeric(log2FoldChange),
      padj = as.numeric(padj),
      direction = if ("direction" %in% nms) as.character(direction) else NA_character_,
      time = time_label
    ) %>%
    filter(!is.na(gene_key), gene_key != "") %>%
    group_by(gene_key, time) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(de_sig = !is.na(padj) & padj < alpha)
}

pathways_selected <- read_csv(pathways_selected_path, show_col_types = FALSE) %>%
  transmute(
    KEGG_group = as.character(KEGG_group),
    ID = as.character(ID),
    Description = as.character(Description),
    sign_class = as.character(sign_class)
  ) %>%
  mutate(
    KEGG_group = factor(KEGG_group, levels = c("Only6h","Shared","Only24h")),
    sign_class = factor(sign_class, levels = c("Activated (NES>0)","Suppressed (NES<0)","Reversal"))
  ) %>%
  distinct(KEGG_group, ID, Description, sign_class, .keep_all = TRUE)

write_csv(pathways_selected, file.path(out_dir, "SUPP_pathways_selected_for_TF_used.csv"))

edges_all <- bind_rows(
  read_edges(edge_6h_path, "Only6h_edges"),
  read_edges(edge_24h_path, "Only24h_edges"),
  read_edges(edge_merged_path, "Merged_edges")
) %>%
  distinct(TF_key, Target_key, .keep_all = TRUE)

write_csv(edges_all, file.path(out_dir, "SUPP_edges_all_parsed_dedup.csv"))

kegg_only6  <- read_kegg_only(kegg_only6_path,  "Only6h")
kegg_only24 <- read_kegg_only(kegg_only24_path, "Only24h")
kegg_shared <- read_kegg_shared(kegg_shared_path)

kegg_map_only6  <- kegg_only6  %>% select(kegg_group, ID, Description, gene_key)
kegg_map_only24 <- kegg_only24 %>% select(kegg_group, ID, Description, gene_key)
kegg_map_shared <- kegg_shared %>% select(kegg_group, ID, Description, gene_key)

kegg_all_genes <- bind_rows(kegg_map_only6, kegg_map_shared, kegg_map_only24) %>%
  rename(KEGG_group = kegg_group) %>%
  inner_join(pathways_selected, by = c("KEGG_group","ID","Description")) %>%
  distinct(KEGG_group, ID, Description, sign_class, gene_key)

phase_map <- read_phase_sets(early_path, late_path, sustained_path)

de6  <- read_de(de6_path, "6h", alpha)
de24 <- read_de(de24_path, "24h", alpha)
de_wide <- bind_rows(de6, de24) %>%
  select(gene_key, time, log2FC, padj, direction, de_sig) %>%
  pivot_wider(names_from = time, values_from = c(log2FC, padj, direction, de_sig), names_sep = "_") %>%
  mutate(
    de_sig_6h  = ifelse(is.na(de_sig_6h),  FALSE, de_sig_6h),
    de_sig_24h = ifelse(is.na(de_sig_24h), FALSE, de_sig_24h)
  )

# Join edges -> KEGG membership -> phase -> DE
edges_all <- edges_all %>%
  distinct(TF_key, Target_key, edge_group, .keep_all = TRUE)  # o sin edge_group si no lo quieres

kegg_all_genes <- kegg_all_genes %>%
  distinct(KEGG_group, ID, Description, sign_class, gene_key)

edges_kegg <- edges_all %>%
  inner_join(kegg_all_genes,
             by = c("Target_key" = "gene_key"),
             relationship = "many-to-many") %>%
  left_join(phase_map, by = c("Target_key" = "gene_key")) %>%
  left_join(de_wide,  by = c("Target_key" = "gene_key")) %>%
  mutate(
    phase = ifelse(is.na(phase), "Unclassified", phase),
    phase = factor(phase, levels = c("Early","Sustained","Late","Unclassified")),
    KEGG_group = factor(KEGG_group, levels = c("Only6h","Shared","Only24h")),
    context = case_when(
      KEGG_group == "Only6h"  ~ "Commitment_6h",
      KEGG_group == "Only24h" ~ "Maintenance_24h",
      TRUE                    ~ "Shared"
    )
  ) %>%
  filter(phase %in% c("Early","Sustained","Late"))


write_csv(edges_kegg, file.path(out_dir, "SUPP_edges_TF_target_annotated_KEGGsign_phase_DE_raw.csv"))

# Split Shared into both contexts, then apply DE-by-context filter
edges_shared_6h <- edges_kegg %>% filter(KEGG_group == "Shared") %>% mutate(context = "Commitment_6h")
edges_shared_24h <- edges_kegg %>% filter(KEGG_group == "Shared") %>% mutate(context = "Maintenance_24h")
edges_nonshared <- edges_kegg %>% filter(KEGG_group %in% c("Only6h","Only24h")) %>%
  mutate(context = ifelse(KEGG_group == "Only6h","Commitment_6h","Maintenance_24h"))

edges_context <- bind_rows(edges_nonshared, edges_shared_6h, edges_shared_24h) %>%
  mutate(context = factor(context, levels = c("Commitment_6h","Maintenance_24h"))) %>%
  filter(
    (context == "Commitment_6h"   & de_sig_6h == TRUE) |
      (context == "Maintenance_24h" & de_sig_24h == TRUE)
  )

write_csv(edges_context, file.path(out_dir, "SUPP_edges_TF_target_FINAL_filtered_DE_and_phase.csv"))

# Summaries
tf_kegg_counts <- edges_context %>%
  group_by(context, KEGG_group, sign_class, Description, TF_key) %>%
  summarize(n_targets = n_distinct(Target_key), .groups = "drop")

write_csv(tf_kegg_counts, file.path(out_dir, "SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv"))

top_tfs <- tf_kegg_counts %>%
  group_by(context, TF_key) %>%
  summarize(total_targets = sum(n_targets), .groups = "drop") %>%
  arrange(context, desc(total_targets)) %>%
  group_by(context) %>%
  slice_head(n = top_n_tfs_for_dotplot) %>%
  ungroup()

# Dotplot (Fig6A): facet por contexto; color por sign_class; x = pathway
tf_plot <- tf_kegg_counts %>%
  inner_join(top_tfs %>% select(context, TF_key), by = c("context","TF_key")) %>%
  mutate(
    TF_key = fct_reorder(TF_key, n_targets, .fun = sum, .desc = TRUE),
    Description = fct_inorder(Description)
  )

# Orden de pathways: por mínima FDR disponible dentro de cada grupo (solo para ordenar el eje)
path_order_only6 <- kegg_only6 %>%
  inner_join(pathways_selected %>% filter(KEGG_group == "Only6h"), by = c("ID","Description")) %>%
  group_by(Description) %>% summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
  mutate(KEGG_group = "Only6h")
path_order_only24 <- kegg_only24 %>%
  inner_join(pathways_selected %>% filter(KEGG_group == "Only24h"), by = c("ID","Description")) %>%
  group_by(Description) %>% summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
  mutate(KEGG_group = "Only24h")
path_order_shared <- kegg_shared %>%
  inner_join(pathways_selected %>% filter(KEGG_group == "Shared"), by = c("ID","Description")) %>%
  group_by(Description) %>% summarize(padj_min = min(padj_6h, padj_24h, na.rm = TRUE), .groups = "drop") %>%
  mutate(KEGG_group = "Shared")

path_order <- bind_rows(path_order_only6, path_order_shared, path_order_only24)

tf_plot <- tf_plot %>%
  left_join(path_order, by = c("KEGG_group","Description")) %>%
  mutate(Description = fct_reorder(Description, padj_min, .desc = FALSE))

p_dot <- ggplot(tf_plot, aes(x = Description, y = TF_key)) +
  geom_point(aes(size = n_targets, color = sign_class), alpha = 0.85) +
  facet_grid(context ~ KEGG_group, scales = "free_x", space = "free_x") +
  scale_size_continuous(range = c(1.5, 8)) +
  labs(
    x = "Selected KEGG pathways (targets: DE + phase-filtered)",
    y = "TF",
    size = "# DE targets",
    color = "Sign",
    title = "TF governors of commitment (6h) vs maintenance (24h) programs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

save_gg(p_dot, "Fig6A_TF_KEGG_dotplot_commitment_vs_maintenance", 13.5, 8, dpi_png)

# Ranking (Fig6B)
tf_rank <- tf_kegg_counts %>%
  group_by(context, TF_key) %>%
  summarize(
    total_targets = sum(n_targets),
    n_pathways_hit = sum(n_targets > 0),
    .groups = "drop"
  ) %>%
  arrange(context, desc(total_targets))

write_csv(tf_rank, file.path(out_dir, "SUPP_TF_ranking_context_DE_phase_full.csv"))

tf_rank_plot <- tf_rank %>%
  group_by(context) %>%
  slice_head(n = top_n_tfs_ranking) %>%
  ungroup() %>%
  mutate(TF_key = fct_reorder(TF_key, total_targets, .desc = TRUE))

p_rank <- ggplot(tf_rank_plot, aes(x = total_targets, y = TF_key)) +
  geom_col() +
  facet_wrap(~ context, scales = "free_y") +
  labs(
    x = "Total # DE targets across selected KEGG pathways",
    y = "TF",
    title = "Top TFs governing commitment (6h) vs maintenance (24h)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_gg(p_rank, "Fig6B_TF_ranking_commitment_vs_maintenance", 11, 6, dpi_png)

message("OK Fig6. Outputs en: ", out_dir)
message("Usando pathways seleccionados de: ", pathways_selected_path)
