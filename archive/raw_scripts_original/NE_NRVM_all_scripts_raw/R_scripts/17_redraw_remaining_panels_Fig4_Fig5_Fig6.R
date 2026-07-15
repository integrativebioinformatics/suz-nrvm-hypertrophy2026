
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

out_dir <- "redrawn_panels_Fig4_Fig5_Fig6"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_dual <- function(p, stem, w = 9, h = 6, dpi = 600){
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h, bg = "white")
}

wrap_lab <- function(x, width = 24) stringr::str_wrap(x, width = width)

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

# -------------------------
# Figure 4A-B: TF panels
# -------------------------
tf_counts <- read.csv("SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv", stringsAsFactors = FALSE)
tf_rank   <- read.csv("SUPP_TF_ranking_context_DE_phase_full.csv", stringsAsFactors = FALSE)

context_map <- c("Commitment_6h" = "6 h state", "Maintenance_24h" = "24 h state")
group_map   <- c("Only6h" = "6 h state-only", "Shared" = "Shared", "Only24h" = "24 h state-only")
sign_map    <- c("Activated (NES>0)" = "Activated", "Suppressed (NES<0)" = "Suppressed")

tf_counts2 <- tf_counts %>%
  mutate(
    context_label = factor(recode(context, !!!context_map), levels = c("6 h state", "24 h state")),
    KEGG_group_label = factor(recode(KEGG_group, !!!group_map), levels = c("6 h state-only", "Shared", "24 h state-only")),
    sign_label = factor(recode(sign_class, !!!sign_map), levels = c("Activated", "Suppressed")),
    pathway_label = wrap_lab(Description, 24),
    TF_label = TF_key
  )

p_tf_bubble <- ggplot(tf_counts2, aes(x = pathway_label, y = TF_label, size = n_targets, color = sign_label)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c("Activated" = "#F8766D", "Suppressed" = "#00BFC4")) +
  scale_size_continuous(range = c(1.8, 8)) +
  facet_grid(context_label ~ KEGG_group_label, scales = "free_x", space = "free_x") +
  labs(
    title = "TF governors of 6 h and 24 h state-associated pathway programs",
    x = "Selected KEGG pathway class",
    y = "Transcription factor",
    size = "# DE targets",
    color = "Pathway sign"
  ) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    legend.position = "right"
  )
save_dual(p_tf_bubble, "Figure4A_TF_by_KEGG_6h_24h_states", w = 14, h = 8)

rank_top <- tf_rank %>%
  mutate(context_label = factor(recode(context, !!!context_map), levels = c("6 h state", "24 h state"))) %>%
  group_by(context_label) %>%
  slice_max(order_by = total_targets, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(context_label, total_targets) %>%
  mutate(TF_plot = factor(paste(TF_key, context_label, sep = "___"), levels = paste(TF_key, context_label, sep = "___")))

p_tf_rank <- ggplot(rank_top, aes(x = total_targets, y = TF_plot)) +
  geom_col(fill = "grey45") +
  facet_wrap(~context_label, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  labs(
    title = "Top TFs across 6 h and 24 h state-associated pathway programs",
    x = "Total # DE targets across selected KEGG pathways",
    y = "Transcription factor"
  )
save_dual(p_tf_rank, "Figure4B_TF_ranking_6h_24h_states", w = 11, h = 5.5)

# -------------------------
# Figure 5D-E: lncRNA panels
# -------------------------
lnc_phase <- read.csv("SUPP_Fig4B2_lncDE_with_phase.csv", stringsAsFactors = FALSE)
lnc_hubs  <- read.csv("SUPP_Fig4C2_hubs_early_vs_late.csv", stringsAsFactors = FALSE)

lnc_phase2 <- lnc_phase %>%
  mutate(
    source_fixed = if_else(grepl("^MSTRG", gene_key), "novel", "annotated"),
    source_label = factor(source_fixed, levels = c("annotated", "novel")),
    time_label = factor(recode(time, "6h" = "6 h state", "24h" = "24 h state"), levels = c("6 h state", "24 h state")),
    phase = factor(phase, levels = c("Early", "Sustained", "Late")),
    direction = factor(direction, levels = c("UP", "DOWN"))
  ) %>%
  count(source_label, time_label, phase, direction, name = "n_genes") %>%
  complete(source_label, time_label, phase, direction, fill = list(n_genes = 0))

p_lnc_counts <- ggplot(lnc_phase2, aes(x = phase, y = n_genes, fill = direction)) +
  geom_col(position = "stack", width = 0.74) +
  facet_grid(source_label ~ time_label) +
  scale_fill_manual(values = c("UP" = "#00BFC4", "DOWN" = "#F8766D")) +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title = "Differentially expressed lncRNAs across Early, Sustained, and Late response classes",
    x = "Response class",
    y = "# DE lncRNA genes",
    fill = "Direction"
  )
save_dual(p_lnc_counts, "Figure5D_lncRNA_temporal_program_counts", w = 10, h = 6.8)

lnc_hubs2 <- lnc_hubs %>%
  mutate(
    state_label = factor(recode(set, "Early" = "Early (6 h state)", "Late" = "Late (24 h state)"),
                         levels = c("Early (6 h state)", "Late (24 h state)")),
    source_label = factor(lnc_source, levels = c("annotated", "novel")),
    lnc_label = if_else(is.na(lnc_gene_name) | lnc_gene_name == "", lncRNA, lnc_gene_name)
  ) %>%
  group_by(state_label) %>%
  slice_max(order_by = n_targets, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(state_label, n_targets) %>%
  mutate(label_plot = factor(paste(lnc_label, state_label, sep = "___"),
                             levels = paste(lnc_label, state_label, sep = "___")))

p_lnc_hubs <- ggplot(lnc_hubs2, aes(x = n_targets, y = label_plot, fill = source_label)) +
  geom_col() +
  facet_wrap(~state_label, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(values = c("annotated" = "#F8766D", "novel" = "#00BFC4")) +
  labs(
    title = "Top lncRNA hubs across Early and Late state-associated programs",
    x = "# mRNA targets (degree)",
    y = "lncRNA",
    fill = "Source"
  )
save_dual(p_lnc_hubs, "Figure5E_lncRNA_hubs_early_late", w = 12, h = 6)

# -------------------------
# Figure 6C and Supplementary Figure 4: miRNA-mRNA KEGG integration
# -------------------------
prepare_mir_plot_tbl <- function(path, main_time = c("24h", "6h")) {
  main_time <- match.arg(main_time)
  df <- read.csv(path, stringsAsFactors = FALSE) %>%
    mutate(
      time_label = if_else(time == "24h", "24 h state", "6 h state"),
      group_label = factor(recode(KEGG_group, !!!group_map), levels = c("6 h state-only", "Shared", "24 h state-only")),
      direction_block = case_when(
        grepl("miRNADown_mRNAUp", combo) ~ "miRNA down → mRNA up",
        grepl("miRNAUp_mRNADown", combo) ~ "miRNA up → mRNA down",
        TRUE ~ combo
      ),
      pathway_label = wrap_lab(Description, 24)
    ) %>%
    group_by(direction_block, group_label, pathway_label, miRNA) %>%
    summarise(
      n_targets = n_distinct(gene_key),
      max_abs_miRNA_fc = max(abs(miRNA_log2FC), na.rm = TRUE),
      .groups = "drop"
    )
  df
}

mi24_dn <- prepare_mir_plot_tbl("SUPP_edges_24h_miRNADown_mRNAUp_KEGG_phase.csv", "24h")
mi24_up <- prepare_mir_plot_tbl("SUPP_edges_24h_miRNAUp_mRNADown_KEGG_phase.csv", "24h")
mi24 <- bind_rows(mi24_dn, mi24_up) %>%
  mutate(
    direction_block = factor(direction_block, levels = c("miRNA down → mRNA up", "miRNA up → mRNA down"))
  )

p_mi24 <- ggplot(mi24, aes(x = pathway_label, y = miRNA, size = n_targets)) +
  geom_point(shape = 21, fill = "black", color = "black", alpha = 0.9) +
  facet_grid(direction_block ~ group_label, scales = "free_x", space = "free_x") +
  scale_size_continuous(range = c(2, 8), breaks = pretty_breaks()) +
  labs(
    title = "24 h miRNA-associated DE target genes mapped to KEGG pathway classes",
    x = "KEGG pathway (DE target genes in this context)",
    y = "miRNAs at 24 h",
    size = "# targets"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )
save_dual(p_mi24, "Figure6C_24h_miRNA_mRNA_KEGG_integration", w = 14, h = 8)

mi6_dn <- prepare_mir_plot_tbl("SUPP_edges_6h_miRNADown_mRNAUp_KEGG_phase.csv", "6h")
mi6_up <- prepare_mir_plot_tbl("SUPP_edges_6h_miRNAUp_mRNADown_KEGG_phase.csv", "6h")
mi6 <- bind_rows(mi6_dn, mi6_up) %>%
  mutate(
    direction_block = factor(direction_block, levels = c("miRNA down → mRNA up", "miRNA up → mRNA down"))
  )

p_mi6 <- ggplot(mi6, aes(x = pathway_label, y = miRNA, size = n_targets)) +
  geom_point(shape = 21, fill = "black", color = "black", alpha = 0.9) +
  facet_grid(direction_block ~ group_label, scales = "free_x", space = "free_x") +
  scale_size_continuous(range = c(2, 8), breaks = pretty_breaks()) +
  labs(
    title = "6 h miRNA-associated DE target genes mapped to KEGG pathway classes",
    x = "KEGG pathway (DE target genes in this context)",
    y = "miRNAs at 6 h",
    size = "# targets"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )
save_dual(p_mi6, "SupplementaryFigure4_6h_miRNA_mRNA_KEGG_integration", w = 14, h = 7.2)

# -------------------------
# Export a brief summary table for checking
# -------------------------
write.csv(lnc_phase2, file.path(out_dir, "CHECK_lnc_counts_by_source_time_phase_direction.csv"), row.names = FALSE)
write.csv(mi24, file.path(out_dir, "CHECK_24h_miRNA_KEGG_aggregated.csv"), row.names = FALSE)
write.csv(mi6, file.path(out_dir, "CHECK_6h_miRNA_KEGG_aggregated.csv"), row.names = FALSE)
