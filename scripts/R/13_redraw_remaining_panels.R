#!/usr/bin/env Rscript
# Stage 13: Final figure panels — consistent redraw (Fig4A/4B, Fig5D/5E, Fig6C)
# Redraws the TF-by-KEGG bubble/ranking, the lncRNA temporal-program counts and
# hub panels, and the 24h/6h miRNA-mRNA-KEGG integration bubbles with unified
# 6h/24h-state nomenclature and styling, from the SUPP_* tables produced by
# scripts 08 and 09.
#
# Ported from 17_redraw_remaining_panels_Fig4_Fig5_Fig6.R with Windows paths
# removed, inputs parameterized, and optparse added.

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
  make_option(c("--tf-dir"), type = "character", default = NULL,
              help = "Directory with script 08 outputs (SUPP_TF_*, SUPP_Fig4B2_*, SUPP_Fig4C2_*)"),
  make_option(c("--mir-dir"), type = "character", default = NULL,
              help = "Directory with script 09 outputs (SUPP_edges_*_KEGG_phase.csv)"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/13_final_panels",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 600, help = "PNG DPI [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
OUTPUT_DIR <- opt$`output-dir`
DPI_PNG <- opt$dpi
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

save_dual <- function(p, stem, w = 9, h = 6, dpi = DPI_PNG) {
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".pdf")), p, width = w, height = h, bg = "white")
}
wrap_lab <- function(x, width = 24) str_wrap(x, width = width)
theme_set(theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), legend.title = element_text(face = "bold")))

context_map <- c("Commitment_6h" = "6 h state", "Maintenance_24h" = "24 h state")
group_map   <- c("Only6h" = "6 h state-only", "Shared" = "Shared", "Only24h" = "24 h state-only")
sign_map    <- c("Activated (NES>0)" = "Activated", "Suppressed (NES<0)" = "Suppressed")
tf_path  <- function(f) if (is.null(opt$`tf-dir`)) NA_character_ else { p <- file.path(opt$`tf-dir`, f); if (file.exists(p)) p else NA_character_ }
mir_path <- function(f) if (is.null(opt$`mir-dir`)) NA_character_ else { p <- file.path(opt$`mir-dir`, f); if (file.exists(p)) p else NA_character_ }

# ---------------------------------------------------------------------------
# Figure 4A / 4B: TF panels
# ---------------------------------------------------------------------------

tf_counts_path <- tf_path("SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv")
tf_rank_path   <- tf_path("SUPP_TF_ranking_context_DE_phase_full.csv")

if (!is.na(tf_counts_path)) {
  tf_counts <- read_csv(tf_counts_path, show_col_types = FALSE) %>%
    mutate(context_label = factor(recode(context, !!!context_map), levels = c("6 h state", "24 h state")),
           KEGG_group_label = factor(recode(KEGG_group, !!!group_map), levels = c("6 h state-only", "Shared", "24 h state-only")),
           sign_label = factor(recode(sign_class, !!!sign_map), levels = c("Activated", "Suppressed")),
           pathway_label = wrap_lab(Description, 24),
           TF_label = if ("TF_key" %in% names(.)) TF_key else TF)
  p_tf_bubble <- ggplot(tf_counts, aes(x = pathway_label, y = TF_label, size = n_targets, color = sign_label)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Activated" = "#F8766D", "Suppressed" = "#00BFC4"), na.value = "grey60") +
    scale_size_continuous(range = c(1.8, 8)) +
    facet_grid(context_label ~ KEGG_group_label, scales = "free_x", space = "free_x") +
    labs(title = "TF governors of 6 h and 24 h state-associated pathway programs",
         x = "Selected KEGG pathway class", y = "Transcription factor", size = "# DE targets", color = "Pathway sign") +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1), legend.position = "right")
  save_dual(p_tf_bubble, "Figure4A_TF_by_KEGG_6h_24h_states", w = 14, h = 8)
  cat("[OUTPUT] Figure4A\n")

  if (!is.na(tf_rank_path)) {
    tf_rank <- read_csv(tf_rank_path, show_col_types = FALSE) %>%
      mutate(context_label = factor(recode(context, !!!context_map), levels = c("6 h state", "24 h state")),
             TF_label = if ("TF_key" %in% names(.)) TF_key else TF)
    rank_top <- tf_rank %>% group_by(context_label) %>%
      slice_max(order_by = total_targets, n = 10, with_ties = FALSE) %>% ungroup() %>%
      arrange(context_label, total_targets) %>%
      mutate(TF_plot = factor(paste(TF_label, context_label, sep = "___"),
                              levels = paste(TF_label, context_label, sep = "___")))
    p_tf_rank <- ggplot(rank_top, aes(x = total_targets, y = TF_plot)) + geom_col(fill = "grey45") +
      facet_wrap(~ context_label, scales = "free_y", nrow = 1) +
      scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
      labs(title = "Top TFs across 6 h and 24 h state-associated pathway programs",
           x = "Total # DE targets across selected KEGG pathways", y = "Transcription factor")
    save_dual(p_tf_rank, "Figure4B_TF_ranking_6h_24h_states", w = 11, h = 5.5)
    cat("[OUTPUT] Figure4B\n")
  }
} else cat("[WARN] TF counts table not found; skipping Fig4A/4B.\n")

# ---------------------------------------------------------------------------
# Figure 5D / 5E: lncRNA temporal-program counts and hubs
# ---------------------------------------------------------------------------

lnc_phase_path <- tf_path("SUPP_Fig4B2_lncDE_with_phase.csv")
if (!is.na(lnc_phase_path)) {
  lnc_phase <- read_csv(lnc_phase_path, show_col_types = FALSE) %>%
    mutate(source_fixed = if_else(grepl("^MSTRG", gene_key), "novel", as.character(lnc_source)),
           source_label = factor(source_fixed, levels = c("annotated", "novel")),
           time_label = factor(recode(time, "6h" = "6 h state", "24h" = "24 h state"), levels = c("6 h state", "24 h state")),
           phase = factor(phase, levels = c("Early", "Sustained", "Late")),
           direction = factor(direction, levels = c("UP", "DOWN"))) %>%
    filter(!is.na(direction)) %>%
    count(source_label, time_label, phase, direction, name = "n_genes") %>%
    complete(source_label, time_label, phase, direction, fill = list(n_genes = 0))
  p_lnc_counts <- ggplot(lnc_phase, aes(x = phase, y = n_genes, fill = direction)) +
    geom_col(position = "stack", width = 0.74) + facet_grid(source_label ~ time_label) +
    scale_fill_manual(values = c("UP" = "#00BFC4", "DOWN" = "#F8766D")) +
    scale_y_continuous(labels = comma_format()) +
    labs(title = "Differentially expressed lncRNAs across Early, Sustained, and Late response classes",
         x = "Response class", y = "# DE lncRNA genes", fill = "Direction")
  save_dual(p_lnc_counts, "Figure5D_lncRNA_temporal_program_counts", w = 10, h = 6.8)
  write_csv(lnc_phase, file.path(OUTPUT_DIR, "CHECK_lnc_counts_by_source_time_phase_direction.csv"))
  cat("[OUTPUT] Figure5D\n")
} else cat("[WARN] SUPP_Fig4B2_lncDE_with_phase.csv not found; skipping Fig5D.\n")

hubs_path <- tf_path("SUPP_Fig4C2_hubs_early_vs_late.csv")
if (!is.na(hubs_path)) {
  lnc_hubs <- read_csv(hubs_path, show_col_types = FALSE) %>%
    mutate(state_label = factor(recode(set, "Early" = "Early (6 h state)", "Late" = "Late (24 h state)"),
                                levels = c("Early (6 h state)", "Late (24 h state)")),
           source_label = factor(lnc_source, levels = c("annotated", "novel")),
           lnc_label = if_else(is.na(lnc_gene_name) | lnc_gene_name == "", lncRNA, lnc_gene_name)) %>%
    group_by(state_label) %>% slice_max(order_by = n_targets, n = 12, with_ties = FALSE) %>% ungroup() %>%
    arrange(state_label, n_targets) %>%
    mutate(label_plot = factor(paste(lnc_label, state_label, sep = "___"),
                               levels = paste(lnc_label, state_label, sep = "___")))
  p_lnc_hubs <- ggplot(lnc_hubs, aes(x = n_targets, y = label_plot, fill = source_label)) + geom_col() +
    facet_wrap(~ state_label, scales = "free_y", nrow = 1) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = c("annotated" = "#F8766D", "novel" = "#00BFC4"), na.value = "grey60") +
    labs(title = "Top lncRNA hubs across Early and Late state-associated programs",
         x = "# mRNA targets (degree)", y = "lncRNA", fill = "Source")
  save_dual(p_lnc_hubs, "Figure5E_lncRNA_hubs_early_late", w = 12, h = 6)
  cat("[OUTPUT] Figure5E\n")
} else cat("[WARN] SUPP_Fig4C2_hubs_early_vs_late.csv not found; skipping Fig5E.\n")

# ---------------------------------------------------------------------------
# Figure 6C + Supplementary Figure 4: miRNA-mRNA-KEGG integration
# ---------------------------------------------------------------------------

prepare_mir_plot_tbl <- function(path) {
  if (is.na(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    mutate(group_label = factor(recode(KEGG_group, !!!group_map), levels = c("6 h state-only", "Shared", "24 h state-only")),
           direction_block = case_when(grepl("miRNADown_mRNAUp", combo) ~ "miRNA down -> mRNA up",
                                       grepl("miRNAUp_mRNADown", combo) ~ "miRNA up -> mRNA down", TRUE ~ combo),
           pathway_label = wrap_lab(Description, 24)) %>%
    group_by(direction_block, group_label, pathway_label, miRNA) %>%
    summarise(n_targets = n_distinct(gene_key), .groups = "drop")
}
draw_mir <- function(dn_path, up_path, title, stem, h) {
  dn <- prepare_mir_plot_tbl(dn_path); up <- prepare_mir_plot_tbl(up_path)
  if (is.null(dn) && is.null(up)) return(invisible(NULL))
  mi <- bind_rows(dn, up) %>%
    mutate(direction_block = factor(direction_block, levels = c("miRNA down -> mRNA up", "miRNA up -> mRNA down")))
  p <- ggplot(mi, aes(x = pathway_label, y = miRNA, size = n_targets)) +
    geom_point(shape = 21, fill = "black", color = "black", alpha = 0.9) +
    facet_grid(direction_block ~ group_label, scales = "free_x", space = "free_x") +
    scale_size_continuous(range = c(2, 8), breaks = pretty_breaks()) +
    labs(title = title, x = "KEGG pathway (DE target genes in this context)", y = "miRNAs", size = "# targets") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "right")
  save_dual(p, stem, w = 14, h = h)
  write_csv(mi, file.path(OUTPUT_DIR, paste0("CHECK_", stem, "_aggregated.csv")))
  cat("[OUTPUT]", stem, "\n")
}
draw_mir(mir_path("SUPP_edges_24h_miRNADown_mRNAUp_KEGG_phase.csv"),
         mir_path("SUPP_edges_24h_miRNAUp_mRNADown_KEGG_phase.csv"),
         "24 h miRNA-associated DE target genes mapped to KEGG pathway classes",
         "Figure6C_24h_miRNA_mRNA_KEGG_integration", 8)
draw_mir(mir_path("SUPP_edges_6h_miRNADown_mRNAUp_KEGG_phase.csv"),
         mir_path("SUPP_edges_6h_miRNAUp_mRNADown_KEGG_phase.csv"),
         "6 h miRNA-associated DE target genes mapped to KEGG pathway classes",
         "SupplementaryFigure4_6h_miRNA_mRNA_KEGG_integration", 7.2)

cat("\n[SUCCESS] Final panel refinement complete.\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
