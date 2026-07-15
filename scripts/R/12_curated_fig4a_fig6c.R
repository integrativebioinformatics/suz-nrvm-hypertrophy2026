#!/usr/bin/env Rscript
# Stage 12: Figures 4A/4B/6C — curated main-text pathway panels
# Filters the TF-by-KEGG and 24h miRNA-mRNA-KEGG integration tables to the
# manuscript's curated biological pathway shortlist and draws the main-text
# TF bubble plot (4A), TF ranking (4B), and 24h miRNA-KEGG bubble plot (6C),
# with harmonized 6h/24h-state nomenclature.
#
# Ported from 18_curated_Fig4A_Fig6C_from_shortlist.R with Windows paths
# removed, the curated shortlist kept as an overridable default, and optparse.

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
  make_option(c("--tf-kegg-file"), type = "character", default = NULL,
              help = "SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv (from script 08)"),
  make_option(c("--mir24-down"), type = "character", default = NULL,
              help = "SUPP_edges_24h_miRNADown_mRNAUp_KEGG_phase.csv (from script 09)"),
  make_option(c("--mir24-up"), type = "character", default = NULL,
              help = "SUPP_edges_24h_miRNAUp_mRNADown_KEGG_phase.csv (from script 09)"),
  make_option(c("--fig4-shortlist"), type = "character", default = NULL,
              help = "Override Fig4 pathway shortlist (comma list or file)"),
  make_option(c("--fig6-shortlist"), type = "character", default = NULL,
              help = "Override Fig6 pathway shortlist (comma list or file)"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/12_curated",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 600, help = "PNG DPI [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
OUTPUT_DIR <- opt$`output-dir`
DPI_PNG <- opt$dpi
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
cat("[INFO] Output directory:", OUTPUT_DIR, "\n")

save_dual <- function(p, stem, w = 12, h = 7, dpi = DPI_PNG) {
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".pdf")), p, width = w, height = h, bg = "white")
}
theme_set(theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")))
relabel_state <- function(x) case_when(x == "Commitment_6h" ~ "6 h state", x == "Maintenance_24h" ~ "24 h state", TRUE ~ x)
relabel_group <- function(x) case_when(x == "Only6h" ~ "6 h state-only", x == "Only24h" ~ "24 h state-only", x == "Shared" ~ "Shared", TRUE ~ x)
parse_list_opt <- function(val, default) {
  if (is.null(val)) return(default)
  v <- if (file.exists(val)) read_lines(val) else str_split(val, ",")[[1]]
  v <- str_trim(v); v[v != ""]
}

# ---------------------------------------------------------------------------
# Curated biological shortlists (manuscript defaults; overridable)
# ---------------------------------------------------------------------------

fig4_only6h  <- c("Focal adhesion", "Cytokine-cytokine receptor interaction", "TGF-beta signaling pathway")
fig4_shared  <- c("Hippo signaling pathway", "Protein processing in endoplasmic reticulum", "Tight junction", "Leukocyte transendothelial migration")
fig4_only24h <- c("Cell cycle", "Cytoskeleton in muscle cells", "PPAR signaling pathway", "Ether lipid metabolism")
fig4_shortlist <- parse_list_opt(opt$`fig4-shortlist`, c(fig4_only6h, fig4_shared, fig4_only24h))
fig4_order <- fig4_shortlist

fig6_shared  <- c("Hippo signaling pathway", "Protein processing in endoplasmic reticulum", "Tight junction", "Leukocyte transendothelial migration")
fig6_only24h <- c("Cell cycle", "Cytoskeleton in muscle cells", "PPAR signaling pathway", "Peroxisome", "Ether lipid metabolism", "Propanoate metabolism")
fig6_shortlist <- parse_list_opt(opt$`fig6-shortlist`, c(fig6_shared, fig6_only24h))
fig6_order <- fig6_shortlist

# ---------------------------------------------------------------------------
# Figure 4A/4B: curated TF-by-KEGG panels
# ---------------------------------------------------------------------------

if (!is.null(opt$`tf-kegg-file`) && file.exists(opt$`tf-kegg-file`)) {
  tf <- read_csv(opt$`tf-kegg-file`, show_col_types = FALSE) %>%
    mutate(state = relabel_state(context), group = relabel_group(KEGG_group), pathway = Description,
           sign = case_when(str_detect(sign_class, "Activated") ~ "Activated",
                            str_detect(sign_class, "Suppressed") ~ "Suppressed", TRUE ~ sign_class)) %>%
    filter(pathway %in% fig4_shortlist)

  if (nrow(tf) > 0) {
    tf_keep <- tf %>% group_by(state, TF_key) %>% summarise(total_targets = sum(n_targets), .groups = "drop") %>%
      arrange(state, desc(total_targets), TF_key) %>% group_by(state) %>% slice_head(n = 10) %>% ungroup() %>%
      pull(TF_key) %>% unique()

    tf_plot <- tf %>% filter(TF_key %in% tf_keep) %>%
      mutate(state = factor(state, levels = c("6 h state", "24 h state")),
             group = factor(group, levels = c("6 h state-only", "Shared", "24 h state-only")),
             pathway = factor(pathway, levels = fig4_order),
             pathway_label = str_wrap(as.character(pathway), width = 20))
    tf_order <- tf_plot %>% group_by(TF_key) %>% summarise(total = sum(n_targets), .groups = "drop") %>%
      arrange(total) %>% pull(TF_key)
    tf_plot$TF_key <- factor(tf_plot$TF_key, levels = tf_order)

    p_tf4a <- ggplot(tf_plot, aes(x = pathway_label, y = TF_key, size = n_targets, color = sign)) +
      geom_point(alpha = 0.9) +
      facet_grid(state ~ group, scales = "free_x", space = "free_x") +
      scale_color_manual(values = c("Activated" = "#F08A81", "Suppressed" = "#25B9C4"), na.value = "grey60") +
      scale_size_continuous(range = c(2.2, 10), breaks = pretty_breaks()) +
      labs(title = "TF governors of curated 6 h and 24 h state-associated pathway programs",
           x = "Selected KEGG pathway class", y = "Transcription factor",
           color = "Pathway sign", size = "# DE targets") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), legend.position = "right")
    save_dual(p_tf4a, "Figure4A_curated_TF_by_KEGG_6h_24h_states", w = 15.5, h = 8.8)
    write_csv(tf_plot, file.path(OUTPUT_DIR, "Figure4A_curated_plot_table.csv"))

    tf_rank <- tf %>% group_by(state, TF_key) %>% summarise(total_targets = sum(n_targets), .groups = "drop") %>%
      group_by(state) %>% arrange(desc(total_targets), .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup()
    p_tf4b <- ggplot(tf_rank, aes(x = total_targets, y = reorder(TF_key, total_targets), fill = state)) +
      geom_col(show.legend = FALSE) + facet_wrap(~ state, scales = "free_y") +
      scale_fill_manual(values = c("6 h state" = "grey50", "24 h state" = "grey50")) +
      labs(title = "Top TFs across curated 6 h and 24 h state-associated pathway programs",
           x = "Total # DE targets across curated KEGG pathways", y = "Transcription factor")
    save_dual(p_tf4b, "Figure4B_curated_TF_ranking_6h_24h_states", w = 12.5, h = 6.2)
    write_csv(tf_rank, file.path(OUTPUT_DIR, "Figure4B_curated_ranking_table.csv"))
    cat("[OUTPUT] Figure4A + Figure4B curated panels\n")
  } else cat("[WARN] No TF-KEGG rows match the Fig4 shortlist.\n")
} else {
  cat("[WARN] --tf-kegg-file not provided or missing; skipping Fig4A/4B.\n")
}

# ---------------------------------------------------------------------------
# Figure 6C: curated 24h miRNA-mRNA-KEGG bubble plot
# ---------------------------------------------------------------------------

prepare_mir24 <- function(path, direction_block_label) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    mutate(group = relabel_group(KEGG_group), pathway = Description, direction_block = direction_block_label) %>%
    filter(group %in% c("Shared", "24 h state-only"), pathway %in% fig6_shortlist) %>%
    group_by(direction_block, group, miRNA, pathway) %>%
    summarise(n_targets = n_distinct(gene_key), .groups = "drop")
}

mi24_dn <- prepare_mir24(opt$`mir24-down`, "miRNA down -> mRNA up")
mi24_up <- prepare_mir24(opt$`mir24-up`, "miRNA up -> mRNA down")

if (!is.null(mi24_dn) || !is.null(mi24_up)) {
  mi24 <- bind_rows(mi24_dn, mi24_up) %>%
    mutate(group = factor(group, levels = c("Shared", "24 h state-only")),
           direction_block = factor(direction_block, levels = c("miRNA down -> mRNA up", "miRNA up -> mRNA down")),
           pathway = factor(pathway, levels = fig6_order),
           pathway_label = str_wrap(as.character(pathway), width = 18))
  mir_keep <- mi24 %>% group_by(direction_block, miRNA) %>% summarise(total = sum(n_targets), .groups = "drop") %>%
    arrange(direction_block, desc(total), miRNA) %>% group_by(direction_block) %>% slice_head(n = 8) %>% ungroup() %>%
    pull(miRNA) %>% unique()
  mi24 <- mi24 %>% filter(miRNA %in% mir_keep)
  mir_order <- mi24 %>% group_by(direction_block, miRNA) %>% summarise(total = sum(n_targets), .groups = "drop") %>%
    arrange(direction_block, total) %>% pull(miRNA) %>% unique()
  mi24$miRNA <- factor(mi24$miRNA, levels = mir_order)

  if (nrow(mi24) > 0) {
    p_mi24 <- ggplot(mi24, aes(x = pathway_label, y = miRNA, size = n_targets)) +
      geom_point(shape = 21, fill = "black", color = "black", alpha = 0.9) +
      facet_grid(direction_block ~ group, scales = "free_x", space = "free_x") +
      scale_size_continuous(range = c(2.5, 9), breaks = pretty_breaks()) +
      labs(title = "24 h miRNA-associated DE target genes mapped to curated KEGG pathway classes",
           x = "Selected KEGG pathway class", y = "miRNAs at 24 h", size = "# targets") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), legend.position = "right")
    save_dual(p_mi24, "Figure6C_curated_24h_miRNA_mRNA_KEGG_integration", w = 12.8, h = 7.2)
    write_csv(mi24, file.path(OUTPUT_DIR, "Figure6C_curated_plot_table.csv"))
    cat("[OUTPUT] Figure6C curated panel\n")
  }
} else {
  cat("[WARN] --mir24-down/--mir24-up not provided; skipping Fig6C.\n")
}

# ---------------------------------------------------------------------------
# Summary note
# ---------------------------------------------------------------------------

summary_lines <- c(
  "Curated shortlist used for main-text pathway figures", "",
  "Figure 4A curated shortlist",
  paste0("- 6 h state-only: ", paste(fig4_only6h, collapse = "; ")),
  paste0("- Shared: ", paste(fig4_shared, collapse = "; ")),
  paste0("- 24 h state-only: ", paste(fig4_only24h, collapse = "; ")), "",
  "Figure 6C curated shortlist",
  paste0("- Shared: ", paste(fig6_shared, collapse = "; ")),
  paste0("- 24 h state-only: ", paste(fig6_only24h, collapse = "; ")), "",
  "Nomenclature harmonization",
  "- Commitment_6h -> 6 h state", "- Maintenance_24h -> 24 h state",
  "- Only6h -> 6 h state-only", "- Only24h -> 24 h state-only", "- Shared -> Shared"
)
writeLines(summary_lines, file.path(OUTPUT_DIR, "README_curated_shortlist_and_nomenclature.txt"))

cat("\n[SUCCESS] Curated KEGG figures complete.\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
