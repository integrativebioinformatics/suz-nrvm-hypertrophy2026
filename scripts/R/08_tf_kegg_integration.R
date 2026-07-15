#!/usr/bin/env Rscript
# Stage 08: Figures 4/5/6 — lncRNA temporal programs, co-expression hubs,
#            KEGG GSEA dynamics, and TF-by-KEGG integration.
# Combines the lncRNA phase panels, lncRNA-mRNA correlation hubs/rewiring,
# KEGG NES dynamics (Only6h/Shared/Only24h + sign classification), and the
# TF GRN -> KEGG integration filtered by mRNA DE and temporal phase.
#
# Ported from script_final_depurado_TF_KEGG.R with Windows paths removed,
# optparse added, external co-expression tables read gracefully, and the
# downstream count/ranking/hub CSVs exported for scripts 04/12/13.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option(c("--edge-6h"), type = "character", help = "TF GRN edges (Only6h)"),
  make_option(c("--edge-24h"), type = "character", help = "TF GRN edges (Only24h)"),
  make_option(c("--edge-merged"), type = "character", help = "Merged GRN edges"),
  make_option(c("--kegg-6h"), type = "character", help = "KEGG Only6h GSEA (long, SYMBOLS)"),
  make_option(c("--kegg-24h"), type = "character", help = "KEGG Only24h GSEA (long, SYMBOLS)"),
  make_option(c("--kegg-shared"), type = "character", help = "KEGG Shared GSEA (long, SYMBOLS)"),
  make_option(c("--de-6h"), type = "character", help = "mRNA DE at 6h (NE_vs_Ctrl_6h_mRNA_DE.txt)"),
  make_option(c("--de-24h"), type = "character", help = "mRNA DE at 24h"),
  make_option(c("--early"), type = "character", help = "Early mRNA phase set"),
  make_option(c("--sustained"), type = "character", help = "Sustained mRNA phase set"),
  make_option(c("--late"), type = "character", help = "Late mRNA phase set"),
  make_option(c("--lnc-dir"), type = "character", default = NULL,
              help = "Directory with lncRNA phase sets + lncRNA DE tables (for Fig4B panels)"),
  make_option(c("--cor-dir"), type = "character", default = NULL,
              help = "Directory with lncRNA-mRNA correlation tables (for hub/rewiring panels)"),
  make_option(c("--pathways-focus"), type = "character", default = NULL,
              help = "Comma-separated pathways or file; default = manuscript curated shortlist"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/08_tf_kegg",
              help = "Output directory [default: %default]"),
  make_option(c("--alpha"), type = "numeric", default = 0.05, help = "DE/enrichment FDR threshold"),
  make_option(c("--dpi"), type = "integer", default = 600, help = "PNG DPI [default: %default]"),
  make_option(c("--top-tfs"), type = "integer", default = 20, help = "Top TFs per context"),
  make_option(c("--top-hubs"), type = "integer", default = 15, help = "Top hubs per set"),
  make_option(c("--top-per-phase"), type = "integer", default = 12, help = "Top lncRNAs per phase (heatmap)"),
  make_option(c("--topn-enrich"), type = "integer", default = 10, help = "Top pathways per sign for selection")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

OUTPUT_DIR <- opt$`output-dir`
DPI_PNG <- opt$dpi
alpha <- opt$alpha
top_lnc_per_phase_heatmap <- opt$`top-per-phase`
top_hubs_plot <- opt$`top-hubs`
topN_enrich <- opt$`topn-enrich`
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

save_gg <- function(p, name, w, h, dpi = DPI_PNG) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, units = "in")
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
}
exists_opt <- function(x) !is.null(x) && file.exists(x)
in_dir <- function(dir, fname) { if (is.null(dir)) return(NA_character_); p <- file.path(dir, fname); if (file.exists(p)) p else NA_character_ }

default_pathways <- c(
  "Focal adhesion", "Cytokine-cytokine receptor interaction", "Hippo signaling pathway",
  "Apoptosis", "Regulation of actin cytoskeleton", "PI3K-Akt signaling pathway", "Cell cycle",
  "Protein processing in endoplasmic reticulum", "TGF-beta signaling pathway",
  "Signaling pathways regulating pluripotency of stem cells", "Peroxisome", "Fatty acid degradation",
  "Adipocytokine signaling pathway", "Calcium signaling pathway",
  "Valine, leucine and isoleucine degradation", "PPAR signaling pathway"
)
pathways_focus <- default_pathways
if (!is.null(opt$`pathways-focus`)) {
  pathways_focus <- if (file.exists(opt$`pathways-focus`)) read_lines(opt$`pathways-focus`) else str_split(opt$`pathways-focus`, ",")[[1]]
  pathways_focus <- str_trim(pathways_focus); pathways_focus <- pathways_focus[pathways_focus != ""]
}

# ===========================================================================
# A. lncRNA temporal-program panels (Fig4B) — optional (needs --lnc-dir)
# ===========================================================================

read_lnc_set <- function(path, phase_label, source_label) {
  if (is.na(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    transmute(phase = phase_label, source = source_label,
              gene_key = toupper(as.character(gene_id)),
              tx_key = toupper(as.character(if ("transcript_id" %in% names(.)) transcript_id else gene_id)))
}
read_lnc_de <- function(path, time_label) {
  if (is.na(path)) return(NULL)
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(gene_key = toupper(as.character(gene_id)),
              gene_name = ifelse(is.na(gene_name) | gene_name == "", toupper(as.character(gene_id)), gene_name),
              log2FC = as.numeric(log2FoldChange), padj = as.numeric(padj),
              direction = as.character(direction), time = time_label,
              de_sig = !is.na(padj) & padj < alpha) %>%
    group_by(gene_key, time) %>% arrange(padj, desc(abs(log2FC))) %>% slice(1) %>% ungroup()
}

ld <- opt$`lnc-dir`
lnc_all <- bind_rows(
  read_lnc_set(in_dir(ld, "early_lncRNA.csv"), "Early", "annotated"),
  read_lnc_set(in_dir(ld, "sustained_lncRNA.csv"), "Sustained", "annotated"),
  read_lnc_set(in_dir(ld, "late_lncRNA.csv"), "Late", "annotated"),
  read_lnc_set(in_dir(ld, "early_novel_lncRNA.csv"), "Early", "novel"),
  read_lnc_set(in_dir(ld, "sustained_novel_lncRNA.csv"), "Sustained", "novel"),
  read_lnc_set(in_dir(ld, "late_novel_lncRNA.csv"), "Late", "novel")
)

if (!is.null(lnc_all) && nrow(lnc_all) > 0) {
  lnc_all <- lnc_all %>% mutate(phase = factor(phase, levels = c("Early", "Sustained", "Late")),
                                source = factor(source, levels = c("annotated", "novel")))
  lnc_phase_map <- lnc_all %>% transmute(gene_key, phase, lnc_source = source) %>% distinct(gene_key, .keep_all = TRUE)
  lnc_counts <- lnc_all %>% group_by(phase, source) %>%
    summarize(n_genes = n_distinct(gene_key), n_transcripts = n_distinct(tx_key), .groups = "drop")
  write_csv(lnc_counts, file.path(OUTPUT_DIR, "SUPP_Fig4B1_lnc_counts_phase_source.csv"))
  p_b1 <- ggplot(lnc_counts, aes(x = phase, y = n_genes, fill = source)) + geom_col(position = "dodge") +
    labs(x = NULL, y = "# lncRNA genes", title = "lncRNA temporal programs: annotated vs novel") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  save_gg(p_b1, "Fig4B1_lnc_gene_counts", 7.5, 4.5)

  lnc_de_all <- bind_rows(read_lnc_de(in_dir(ld, "NE_vs_Ctrl_6h_lncRNA_DE.txt"), "6h"),
                          read_lnc_de(in_dir(ld, "NE_vs_Ctrl_24h_lncRNA_DE.txt"), "24h"))
  if (!is.null(lnc_de_all)) {
    lnc_de_phase <- lnc_de_all %>% inner_join(lnc_phase_map, by = "gene_key") %>%
      mutate(direction = ifelse(is.na(direction) | direction == "", "NA", direction),
             phase = factor(phase, levels = c("Early", "Sustained", "Late")),
             lnc_source = factor(lnc_source, levels = c("annotated", "novel")))
    write_csv(lnc_de_phase, file.path(OUTPUT_DIR, "SUPP_Fig4B2_lncDE_with_phase.csv"))
    lnc_de_counts <- lnc_de_phase %>% filter(de_sig) %>%
      group_by(time, phase, lnc_source, direction) %>% summarize(n = n_distinct(gene_key), .groups = "drop")
    write_csv(lnc_de_counts, file.path(OUTPUT_DIR, "SUPP_Fig4B2_lncDE_counts_phase_time_direction.csv"))
    p_b2 <- ggplot(lnc_de_counts, aes(x = phase, y = n, fill = direction)) + geom_col(position = "stack") +
      facet_grid(time ~ lnc_source) + labs(x = NULL, y = "# DE lncRNA genes",
      title = "DE lncRNAs across temporal programs (6h vs 24h)") +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
    save_gg(p_b2, "Fig4B2_lncDE_counts_phase_time_direction", 10.5, 5.2)
  }
  cat("[INFO] lncRNA temporal-program panels written\n")
} else {
  cat("[WARN] --lnc-dir not provided or empty; skipping Fig4B lncRNA panels\n")
}

# ===========================================================================
# B. lncRNA-mRNA co-expression hubs + rewiring (Fig4C) — optional (--cor-dir)
# ===========================================================================

standardize_cor <- function(df, label) {
  nms <- names(df)
  pick <- function(exact, rx) if (exact %in% nms) exact else nms[str_detect(nms, regex(rx, TRUE))][1]
  lnc_col <- pick("lncRNA", "^lnc"); lncname_col <- pick("lnc_gene_name", "lnc.*gene.*name")
  src_col <- pick("lnc_source", "lnc.*source"); mrna_col <- pick("mRNA_gene_name", "mrna.*gene.*name|^mRNA$|target|symbol")
  rho_col <- pick("rho", "^rho$|spearman|pearson|corr")
  df %>% transmute(
    lncRNA = if (!is.na(lnc_col)) as.character(.data[[lnc_col]]) else NA_character_,
    lnc_gene_name = if (!is.na(lncname_col)) as.character(.data[[lncname_col]]) else NA_character_,
    lnc_source = if (!is.na(src_col)) as.character(.data[[src_col]]) else "unknown",
    mRNA_gene_name = if (!is.na(mrna_col)) as.character(.data[[mrna_col]]) else NA_character_,
    rho = if (!is.na(rho_col)) as.numeric(.data[[rho_col]]) else NA_real_, set = label) %>%
    mutate(lnc_source = ifelse(is.na(lnc_source) | lnc_source %in% c("", "NA"), "unknown", lnc_source),
           lnc_key = toupper(lncRNA), mRNA_key = toupper(mRNA_gene_name))
}
hub_from_cor <- function(df) df %>% filter(!is.na(lncRNA), !is.na(mRNA_key)) %>%
  group_by(set, lncRNA, lnc_gene_name, lnc_source) %>%
  summarize(n_targets = n_distinct(mRNA_key),
            n_pos = n_distinct(mRNA_key[!is.na(rho) & rho > 0]),
            n_neg = n_distinct(mRNA_key[!is.na(rho) & rho < 0]), .groups = "drop") %>%
  arrange(set, desc(n_targets))
read_cor <- function(path, label) { if (is.na(path)) return(NULL); standardize_cor(read_csv(path, show_col_types = FALSE), label) }

cd <- opt$`cor-dir`
cor_early <- read_cor(in_dir(cd, "correlations_early.csv"), "Early")
cor_late  <- read_cor(in_dir(cd, "correlations_late.csv"),  "Late")
if (!is.null(cor_early) || !is.null(cor_late)) {
  hubs_early_late <- hub_from_cor(bind_rows(cor_early, cor_late))
  write_csv(hubs_early_late, file.path(OUTPUT_DIR, "SUPP_Fig4C2_hubs_early_vs_late.csv"))
  # Also emit the time-based hub table consumed by script 04 (Fig4E)
  hubs_early_late %>%
    transmute(time = set, lncRNA, lnc_gene_name, lnc_source, n_mRNAs = n_targets, mean_rho = NA_real_) %>%
    write_csv(file.path(OUTPUT_DIR, "lncRNA_hubs_by_time.csv"))
  hubs_plot_df <- hubs_early_late %>% group_by(set) %>% slice_head(n = top_hubs_plot) %>% ungroup() %>%
    mutate(label = ifelse(!is.na(lnc_gene_name) & lnc_gene_name != "", lnc_gene_name, lncRNA),
           label = fct_reorder(label, n_targets, .desc = TRUE),
           set = factor(set, levels = c("Early", "Late")),
           lnc_source = factor(lnc_source, levels = c("novel", "annotated", "unknown")))
  p_c2 <- ggplot(hubs_plot_df, aes(x = n_targets, y = label, fill = lnc_source)) + geom_col() +
    facet_wrap(~ set, scales = "free_y") + labs(x = "# mRNA targets (degree)", y = NULL,
    title = "Top lncRNA hubs: commitment (Early) vs maintenance (Late)") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  save_gg(p_c2, "Fig4C2_hubs_early_vs_late", 11, 6)
  cat("[INFO] lncRNA hub panels written\n")
} else {
  cat("[WARN] --cor-dir not provided or no correlation tables; skipping Fig4C hub panels\n")
}

# ===========================================================================
# C. KEGG NES dynamics + pathway sign selection (Fig5)
# ===========================================================================

need_kegg <- exists_opt(opt$`kegg-6h`) && exists_opt(opt$`kegg-24h`) && exists_opt(opt$`kegg-shared`)
kegg_stats <- NULL
pathways_selected <- NULL
if (need_kegg) {
  kegg_only6 <- read_tsv(opt$`kegg-6h`, show_col_types = FALSE) %>%
    transmute(ID = as.character(ID), Description = as.character(Description),
              NES_6h = as.numeric(NES), padj_6h = as.numeric(`p.adjust`),
              NES_24h = NA_real_, padj_24h = NA_real_,
              gene_symbol = as.character(gene_symbol), KEGG_group = "Only6h")
  kegg_only24 <- read_tsv(opt$`kegg-24h`, show_col_types = FALSE) %>%
    transmute(ID = as.character(ID), Description = as.character(Description),
              NES_6h = NA_real_, padj_6h = NA_real_,
              NES_24h = as.numeric(NES), padj_24h = as.numeric(`p.adjust`),
              gene_symbol = as.character(gene_symbol), KEGG_group = "Only24h")
  kegg_shared <- read_tsv(opt$`kegg-shared`, show_col_types = FALSE) %>%
    transmute(ID = as.character(ID), Description = as.character(Description),
              NES_6h = as.numeric(NES_6h), NES_24h = as.numeric(NES_24h),
              padj_6h = as.numeric(`p.adjust_6h`), padj_24h = as.numeric(`p.adjust_24h`),
              gene_symbol = as.character(gene_symbol), KEGG_group = "Shared")

  kegg_stats <- bind_rows(kegg_only6, kegg_shared, kegg_only24) %>%
    distinct(KEGG_group, ID, Description, NES_6h, NES_24h, padj_6h, padj_24h) %>%
    mutate(KEGG_group = factor(KEGG_group, levels = c("Only6h", "Shared", "Only24h")),
           padj_min = pmin(padj_6h, padj_24h, na.rm = TRUE),
           transition_shared = case_when(
             KEGG_group != "Shared" ~ NA_character_,
             NES_6h > 0 & NES_24h > 0 ~ "Maintained_Activated",
             NES_6h < 0 & NES_24h < 0 ~ "Maintained_Suppressed",
             NES_6h * NES_24h < 0 ~ "Reversal", TRUE ~ "Other"))
  write_csv(kegg_stats, file.path(OUTPUT_DIR, "SUPP_Fig5_full_KEGG_stats_wide.csv"))

  sel <- function(grp, nes, padjc, cond, sign, n = topN_enrich, trans = NA_character_) {
    kegg_stats %>% filter(KEGG_group == grp, !is.na(.data[[nes]]), !is.na(.data[[padjc]]),
                          .data[[padjc]] < alpha, cond(.data[[nes]])) %>%
      arrange(.data[[padjc]]) %>% slice_head(n = n) %>%
      transmute(KEGG_group, ID, Description, sign_class = sign, transition_shared = trans)
  }
  sh <- function(trans, sign, n = topN_enrich) kegg_stats %>%
    filter(KEGG_group == "Shared", padj_min < alpha, transition_shared == trans) %>%
    arrange(padj_min) %>% slice_head(n = n) %>%
    transmute(KEGG_group, ID, Description, sign_class = sign, transition_shared)
  pathways_selected <- bind_rows(
    sel("Only6h", "NES_6h", "padj_6h", function(z) z > 0, "Activated (NES>0)"),
    sel("Only6h", "NES_6h", "padj_6h", function(z) z < 0, "Suppressed (NES<0)"),
    sh("Maintained_Activated", "Activated (NES>0)"), sh("Maintained_Suppressed", "Suppressed (NES<0)"),
    sh("Reversal", "Reversal"),
    sel("Only24h", "NES_24h", "padj_24h", function(z) z > 0, "Activated (NES>0)"),
    sel("Only24h", "NES_24h", "padj_24h", function(z) z < 0, "Suppressed (NES<0)")
  ) %>% distinct(KEGG_group, ID, Description, sign_class, .keep_all = TRUE)
  write_csv(pathways_selected, file.path(OUTPUT_DIR, "pathways_selected_for_TF.csv"))

  shared_df <- kegg_stats %>% filter(KEGG_group == "Shared", padj_min < alpha) %>%
    mutate(transition_shared = factor(transition_shared,
      levels = c("Maintained_Activated", "Maintained_Suppressed", "Reversal", "Other")))
  if (nrow(shared_df) > 0) {
    p_5b <- ggplot(shared_df, aes(x = NES_6h, y = NES_24h, size = -log10(padj_min), color = transition_shared)) +
      geom_point(alpha = 0.85) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "NES (6h)", y = "NES (24h)", title = "Shared KEGG pathways: maintained vs reversal", size = "-log10(min FDR)") +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
    save_gg(p_5b, "Fig5B_Shared_scatter_NES6h_vs_NES24h", 7.6, 6)
  }
  cat("[INFO] KEGG stats + pathway selection written\n")
} else {
  cat("[WARN] KEGG inputs incomplete; skipping Fig5 KEGG panels\n")
}

# ===========================================================================
# D. TF GRN -> KEGG integration filtered by DE + phase (Fig6)
# ===========================================================================

if (!(need_kegg && exists_opt(opt$`de-6h`) && exists_opt(opt$early))) {
  cat("[WARN] TF-KEGG integration needs KEGG + --de-6h + --early; skipping Fig6.\n")
  cat("\n[DONE] Stage 08 partial.\n"); quit(status = 0)
}

read_edges <- function(path, grp) {
  if (!exists_opt(path)) return(NULL)
  df <- read_csv(path, show_col_types = FALSE)
  name_col <- if ("shared name" %in% names(df)) "shared name" else if ("name" %in% names(df)) "name" else NA_character_
  if (is.na(name_col)) { cat("[WARN] edge table missing name column:", path, "\n"); return(NULL) }
  parts <- str_split_fixed(df[[name_col]], " \\(interacts with\\) ", 2)
  tibble(edge_group = grp, TF = str_trim(parts[, 1]), Target = str_trim(parts[, 2])) %>%
    mutate(TF_key = str_to_upper(TF), Target_key = str_to_upper(Target)) %>%
    filter(TF_key != "", Target_key != "")
}
read_kegg_genes <- function() {
  bind_rows(
    read_tsv(opt$`kegg-6h`, show_col_types = FALSE) %>% transmute(KEGG_group = "Only6h", ID, Description, gene_key = str_to_upper(gene_symbol)),
    read_tsv(opt$`kegg-24h`, show_col_types = FALSE) %>% transmute(KEGG_group = "Only24h", ID, Description, gene_key = str_to_upper(gene_symbol)),
    read_tsv(opt$`kegg-shared`, show_col_types = FALSE) %>% transmute(KEGG_group = "Shared", ID, Description, gene_key = str_to_upper(gene_symbol))
  ) %>% filter(gene_key != "") %>% distinct()
}
read_de <- function(path, t) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(gene_key = str_to_upper(gene_name), log2FC = as.numeric(log2FoldChange),
              padj = as.numeric(padj), direction = as.character(direction), time = t,
              de_sig = !is.na(padj) & padj < alpha) %>%
    group_by(gene_key, time) %>% arrange(padj, desc(abs(log2FC))) %>% slice(1) %>% ungroup()
}
read_phase <- function(path, lab) {
  if (!exists_opt(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% transmute(gene_key = str_to_upper(gene_name), phase = lab) %>% filter(gene_key != "")
}

edges_all <- bind_rows(read_edges(opt$`edge-6h`, "Only6h"), read_edges(opt$`edge-24h`, "Only24h"),
                       read_edges(opt$`edge-merged`, "Merged")) %>% distinct(TF_key, Target_key, .keep_all = TRUE)
kegg_genes <- read_kegg_genes()
target2kegg <- kegg_genes %>% filter(Description %in% pathways_focus) %>% distinct(KEGG_group, Description, gene_key)

phase_map <- bind_rows(read_phase(opt$sustained, "Sustained"), read_phase(opt$early, "Early"),
                       read_phase(opt$late, "Late")) %>% distinct(gene_key, .keep_all = TRUE)
de_wide <- bind_rows(read_de(opt$`de-6h`, "6h"), read_de(opt$`de-24h`, "24h")) %>%
  select(gene_key, time, log2FC, padj, direction, de_sig) %>%
  pivot_wider(names_from = time, values_from = c(log2FC, padj, direction, de_sig), names_sep = "_") %>%
  mutate(de_sig_6h = replace_na(de_sig_6h, FALSE), de_sig_24h = replace_na(de_sig_24h, FALSE))

edges_kegg <- edges_all %>%
  inner_join(target2kegg, by = c("Target_key" = "gene_key"), relationship = "many-to-many") %>%
  mutate(KEGG_group = factor(KEGG_group, levels = c("Only6h", "Shared", "Only24h"))) %>%
  left_join(phase_map, by = c("Target_key" = "gene_key")) %>% mutate(phase = replace_na(phase, "Unclassified")) %>%
  left_join(de_wide, by = c("Target_key" = "gene_key")) %>%
  filter(phase %in% c("Early", "Sustained", "Late"))

# Split Shared into both contexts, then DE-filter by context timepoint
edges_ctx <- bind_rows(
  edges_kegg %>% filter(KEGG_group %in% c("Only6h", "Only24h")) %>%
    mutate(context = ifelse(KEGG_group == "Only6h", "Commitment_6h", "Maintenance_24h")),
  edges_kegg %>% filter(KEGG_group == "Shared") %>% mutate(context = "Commitment_6h"),
  edges_kegg %>% filter(KEGG_group == "Shared") %>% mutate(context = "Maintenance_24h")
) %>%
  mutate(context = factor(context, levels = c("Commitment_6h", "Maintenance_24h"))) %>%
  filter((context == "Commitment_6h" & de_sig_6h) | (context == "Maintenance_24h" & de_sig_24h))
write_csv(edges_ctx, file.path(OUTPUT_DIR, "SUPP_edges_TF_target_FINAL_filtered_DE_and_phase.csv"))

# TF x KEGG counts, with pathway sign_class joined (consumed by scripts 12/13)
pathway_sign <- if (!is.null(pathways_selected)) {
  pathways_selected %>% distinct(KEGG_group, Description, sign_class)
} else tibble(KEGG_group = character(), Description = character(), sign_class = character())

tf_kegg_counts <- edges_ctx %>%
  group_by(context, KEGG_group, Description, TF) %>%
  summarize(n_targets = n_distinct(Target_key), .groups = "drop") %>%
  mutate(TF_key = str_to_upper(TF)) %>%
  left_join(pathway_sign, by = c("KEGG_group", "Description")) %>%
  mutate(sign_class = replace_na(sign_class, "Other"))
write_csv(tf_kegg_counts, file.path(OUTPUT_DIR, "SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv"))

tf_rank <- tf_kegg_counts %>%
  group_by(context, TF, TF_key) %>%
  summarize(total_targets = sum(n_targets), n_pathways_hit = sum(n_targets > 0), .groups = "drop") %>%
  arrange(context, desc(total_targets))
write_csv(tf_rank, file.path(OUTPUT_DIR, "SUPP_TF_ranking_context_DE_phase_full.csv"))

# Panel A: TF x KEGG dotplot
top_tfs <- tf_rank %>% group_by(context) %>% slice_head(n = opt$`top-tfs`) %>% ungroup()
tf_plot <- tf_kegg_counts %>% semi_join(top_tfs, by = c("context", "TF")) %>%
  mutate(TF = fct_reorder(TF, n_targets, .fun = sum, .desc = TRUE))
p_dot <- ggplot(tf_plot, aes(x = Description, y = TF)) +
  geom_point(aes(size = n_targets, color = sign_class), alpha = 0.85) +
  scale_color_manual(values = c(`Activated (NES>0)` = "#F8766D", `Suppressed (NES<0)` = "#00BFC4",
                                Reversal = "#7C3AED", Other = "grey55")) +
  facet_wrap(~ context, scales = "free_x") + scale_size_continuous(range = c(1.5, 8)) +
  labs(x = "KEGG pathway", y = "TF", size = "# DE targets", color = "Pathway sign",
       title = "TF governance: commitment (6h) vs maintenance (24h)") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
save_gg(p_dot, "Fig6A_TF_KEGG_dotplot", 13, 7)

# Panel B: TF ranking
tf_rank_plot <- tf_rank %>% group_by(context) %>% slice_head(n = opt$`top-tfs`) %>% ungroup() %>%
  mutate(TF = fct_reorder(TF, total_targets, .desc = TRUE))
p_rank <- ggplot(tf_rank_plot, aes(x = total_targets, y = TF)) + geom_col(fill = "grey45") +
  facet_wrap(~ context, scales = "free_y") +
  labs(x = "Total # DE targets across selected KEGG pathways", y = "TF",
       title = "Top TFs: commitment (6h) vs maintenance (24h)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_gg(p_rank, "Fig6B_TF_ranking", 11, 6)

cat("\n[SUCCESS] TF-KEGG integration complete.\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
