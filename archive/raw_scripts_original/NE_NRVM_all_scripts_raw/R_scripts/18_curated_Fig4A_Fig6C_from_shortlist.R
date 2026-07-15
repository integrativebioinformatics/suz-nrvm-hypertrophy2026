suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(stringr)
})

base_dir <- "C:/Users/sebau/Downloads"
out_dir  <- file.path(base_dir, "Curated_shortlist_Fig4A_Fig6C")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_dual <- function(p, stem, w = 12, h = 7, dpi = 320) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h, bg = "white")
}

theme_set(
  theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
)

relabel_state <- function(x) {
  case_when(
    x == "Commitment_6h" ~ "6 h state",
    x == "Maintenance_24h" ~ "24 h state",
    TRUE ~ x
  )
}

relabel_group <- function(x) {
  case_when(
    x == "Only6h" ~ "6 h state-only",
    x == "Only24h" ~ "24 h state-only",
    x == "Shared" ~ "Shared",
    TRUE ~ x
  )
}

# ----------------------------
# Curated biological shortlist
# ----------------------------

fig4_only6h <- c(
  "Focal adhesion",
  "Cytokine-cytokine receptor interaction",
  "TGF-beta signaling pathway"
)

fig4_shared <- c(
  "Hippo signaling pathway",
  "Protein processing in endoplasmic reticulum",
  "Tight junction",
  "Leukocyte transendothelial migration"
)

fig4_only24h <- c(
  "Cell cycle",
  "Cytoskeleton in muscle cells",
  "PPAR signaling pathway",
  "Ether lipid metabolism"
)

fig4_shortlist <- c(fig4_only6h, fig4_shared, fig4_only24h)

fig6_shared <- c(
  "Hippo signaling pathway",
  "Protein processing in endoplasmic reticulum",
  "Tight junction",
  "Leukocyte transendothelial migration"
)

fig6_only24h <- c(
  "Cell cycle",
  "Cytoskeleton in muscle cells",
  "PPAR signaling pathway",
  "Peroxisome",
  "Ether lipid metabolism",
  "Propanoate metabolism"
)

fig6_shortlist <- c(fig6_shared, fig6_only24h)

# pathway display order
fig4_order <- c(fig4_only6h, fig4_shared, fig4_only24h)
fig6_order <- c(fig6_shared, fig6_only24h)

# ----------------------------
# Figure 4A curated TF bubble
# ----------------------------

tf_file <- file.path(base_dir, "SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv")
if (file.exists(tf_file)) {
  tf <- read_csv(tf_file, show_col_types = FALSE) %>%
    mutate(
      state = relabel_state(context),
      group = relabel_group(KEGG_group),
      pathway = Description,
      sign = case_when(
        sign_class == "Activated (NES>0)" ~ "Activated",
        sign_class == "Suppressed (NES<0)" ~ "Suppressed",
        TRUE ~ sign_class
      )
    ) %>%
    filter(pathway %in% fig4_shortlist)

  # keep the most informative TFs for shortlisted pathways
  tf_keep <- tf %>%
    group_by(state, TF_key) %>%
    summarise(total_targets = sum(n_targets), .groups = "drop") %>%
    arrange(state, desc(total_targets), TF_key) %>%
    group_by(state) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    pull(TF_key) %>%
    unique()

  tf_plot <- tf %>%
    filter(TF_key %in% tf_keep) %>%
    mutate(
      state = factor(state, levels = c("6 h state", "24 h state")),
      group = factor(group, levels = c("6 h state-only", "Shared", "24 h state-only")),
      pathway = factor(pathway, levels = fig4_order),
      pathway_label = str_wrap(as.character(pathway), width = 20),
      TF_key = factor(TF_key, levels = rev(sort(unique(TF_key))))
    )

  # reorder TFs by total counts across shortlist
  tf_order <- tf_plot %>%
    group_by(TF_key) %>%
    summarise(total = sum(n_targets), .groups = "drop") %>%
    arrange(total) %>%
    pull(TF_key)
  tf_plot$TF_key <- factor(tf_plot$TF_key, levels = tf_order)

  p_tf4a <- ggplot(tf_plot, aes(x = pathway_label, y = TF_key, size = n_targets, color = sign)) +
    geom_point(alpha = 0.9) +
    facet_grid(state ~ group, scales = "free_x", space = "free_x") +
    scale_color_manual(values = c("Activated" = "#F08A81", "Suppressed" = "#25B9C4")) +
    scale_size_continuous(range = c(2.2, 10), breaks = pretty_breaks()) +
    labs(
      title = "TF governors of curated 6 h and 24 h state-associated pathway programs",
      x = "Selected KEGG pathway class",
      y = "Transcription factor",
      color = "Pathway sign",
      size = "# DE targets"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      legend.position = "right"
    )

  save_dual(p_tf4a, "Figure4A_curated_TF_by_KEGG_6h_24h_states", w = 15.5, h = 8.8)

  # Optional curated ranking matching the shortlist
  tf_rank <- tf %>%
    group_by(state, TF_key) %>%
    summarise(total_targets = sum(n_targets), .groups = "drop") %>%
    group_by(state) %>%
    arrange(desc(total_targets), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()

  p_tf4b_cur <- ggplot(tf_rank, aes(x = total_targets, y = reorder(TF_key, total_targets), fill = state)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ state, scales = "free_y") +
    scale_fill_manual(values = c("6 h state" = "grey50", "24 h state" = "grey50")) +
    labs(
      title = "Top TFs across curated 6 h and 24 h state-associated pathway programs",
      x = "Total # DE targets across curated KEGG pathways",
      y = "Transcription factor"
    )
  save_dual(p_tf4b_cur, "Figure4B_curated_TF_ranking_6h_24h_states", w = 12.5, h = 6.2)

  write_csv(tf_plot, file.path(out_dir, "Figure4A_curated_plot_table.csv"))
  write_csv(tf_rank, file.path(out_dir, "Figure4B_curated_ranking_table.csv"))
}

# ----------------------------
# Figure 6C curated 24 h miRNA
# ----------------------------

prepare_mir24_tbl <- function(path, direction_block_label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      group = relabel_group(KEGG_group),
      pathway = Description,
      direction_block = direction_block_label
    ) %>%
    filter(group %in% c("Shared", "24 h state-only"), pathway %in% fig6_shortlist) %>%
    group_by(direction_block, group, miRNA, pathway) %>%
    summarise(
      n_targets = n_distinct(gene_key),
      .groups = "drop"
    )
}

mi24_dn_file <- file.path(base_dir, "SUPP_edges_24h_miRNADown_mRNAUp_KEGG_phase.csv")
mi24_up_file <- file.path(base_dir, "SUPP_edges_24h_miRNAUp_mRNADown_KEGG_phase.csv")

if (file.exists(mi24_dn_file) && file.exists(mi24_up_file)) {
  mi24_dn <- prepare_mir24_tbl(mi24_dn_file, "miRNA down -> mRNA up")
  mi24_up <- prepare_mir24_tbl(mi24_up_file, "miRNA up -> mRNA down")

  mi24 <- bind_rows(mi24_dn, mi24_up) %>%
    mutate(
      group = factor(group, levels = c("Shared", "24 h state-only")),
      direction_block = factor(direction_block, levels = c("miRNA down -> mRNA up", "miRNA up -> mRNA down")),
      pathway = factor(pathway, levels = fig6_order),
      pathway_label = str_wrap(as.character(pathway), width = 18)
    )

  # keep the most informative miRNAs per direction block
  mir_keep <- mi24 %>%
    group_by(direction_block, miRNA) %>%
    summarise(total_targets = sum(n_targets), .groups = "drop") %>%
    arrange(direction_block, desc(total_targets), miRNA) %>%
    group_by(direction_block) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    pull(miRNA) %>%
    unique()

  mi24 <- mi24 %>% filter(miRNA %in% mir_keep)

  # order miRNAs within direction by total targets
  mir_order <- mi24 %>%
    group_by(direction_block, miRNA) %>%
    summarise(total_targets = sum(n_targets), .groups = "drop") %>%
    arrange(direction_block, total_targets) %>%
    pull(miRNA) %>%
    unique()

  mi24$miRNA <- factor(mi24$miRNA, levels = mir_order)

  p_mi24 <- ggplot(mi24, aes(x = pathway_label, y = miRNA, size = n_targets)) +
    geom_point(shape = 21, fill = "black", color = "black", alpha = 0.9) +
    facet_grid(direction_block ~ group, scales = "free_x", space = "free_x") +
    scale_size_continuous(range = c(2.5, 9), breaks = pretty_breaks()) +
    labs(
      title = "24 h miRNA-associated DE target genes mapped to curated KEGG pathway classes",
      x = "Selected KEGG pathway class",
      y = "miRNAs at 24 h",
      size = "# targets"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      legend.position = "right"
    )

  save_dual(p_mi24, "Figure6C_curated_24h_miRNA_mRNA_KEGG_integration", w = 12.8, h = 7.2)
  write_csv(mi24, file.path(out_dir, "Figure6C_curated_plot_table.csv"))
}

# ----------------------------
# Summary note
# ----------------------------

summary_lines <- c(
  "Curated shortlist used for main-text pathway figures",
  "",
  "Figure 4A curated shortlist",
  paste0("- 6 h state-only: ", paste(fig4_only6h, collapse = "; ")),
  paste0("- Shared: ", paste(fig4_shared, collapse = "; ")),
  paste0("- 24 h state-only: ", paste(fig4_only24h, collapse = "; ")),
  "",
  "Figure 6C curated shortlist",
  paste0("- Shared: ", paste(fig6_shared, collapse = "; ")),
  paste0("- 24 h state-only: ", paste(fig6_only24h, collapse = "; ")),
  "",
  "Nomenclature harmonization",
  "- Commitment_6h -> 6 h state",
  "- Maintenance_24h -> 24 h state",
  "- Only6h -> 6 h state-only",
  "- Only24h -> 24 h state-only",
  "- Shared -> Shared",
  "- miRNA arrows use ASCII (->) to avoid graphics-device warnings"
)

writeLines(summary_lines, file.path(out_dir, "README_curated_shortlist_and_nomenclature.txt"))
