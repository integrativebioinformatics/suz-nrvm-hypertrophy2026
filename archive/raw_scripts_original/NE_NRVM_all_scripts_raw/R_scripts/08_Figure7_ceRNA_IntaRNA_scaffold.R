suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
})

# =========================================================
# Figure 7 scaffold: lncRNA-miRNA-mRNA ceRNA integration
# Requires:
# 1) IntaRNA output linking lncRNAs and miRNAs
# 2) miRNA -> mRNA target table (e.g., miRWalk merged table)
# 3) DE tables for lncRNAs, miRNAs and mRNAs in the same context (e.g., 24h)
# Optional:
# 4) lncRNA and mRNA phase/class tables for Early/Sustained/Late annotation
# =========================================================

base_dir <- "."
bulk_dir <- file.path(base_dir, "DESeq2_Multifactorial_Results_fixed")
fig45_dir <- file.path(base_dir, "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT")
out_dir <- file.path(base_dir, "Paper_Fig7_ceRNA_FINAL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------- user-editable file names -----------------
context_label <- "24h"

# IntaRNA table. Expected columns after renaming:
# lncRNA, miRNA, energy
intarna_file <- file.path(base_dir, "IntaRNA_lncRNA_miRNA_24h.csv")

# miRNA->mRNA target table. Expected columns after renaming:
# miRNA, target_gene, binding_region, score(optional)
mirna_target_file <- file.path(base_dir, "miRNA_mRNA_targets_24h_merged.csv")

# DE tables
lnc_de_file   <- file.path(base_dir, "NE_vs_Ctrl_24h_lncRNA_DE.txt")
mrna_de_file  <- file.path(bulk_dir, "NE_vs_Ctrl_24h_mRNA_DE.txt")
mirna_de_file <- file.path(base_dir, "DE_miRNA_24h.txt")

# Optional phase tables exported/available from previous scripts
lnc_phase_files <- c(
  Early = file.path(base_dir, "early_lncRNA.csv"),
  Sustained = file.path(base_dir, "sustained_lncRNA.csv"),
  Late = file.path(base_dir, "late_lncRNA.csv"),
  Early_novel = file.path(base_dir, "early_novel_lncRNA.csv"),
  Sustained_novel = file.path(base_dir, "sustained_novel_lncRNA.csv"),
  Late_novel = file.path(base_dir, "late_novel_lncRNA.csv")
)
mrna_phase_files <- c(
  Early = file.path(bulk_dir, "early_mRNA.csv"),
  Sustained = file.path(bulk_dir, "sustained_mRNA.csv"),
  Late = file.path(bulk_dir, "late_mRNA.csv")
)

# Filtering thresholds
padj_cutoff <- 0.05
absfc_cutoff <- 0
energy_cutoff <- -15
max_triplets_for_network <- 60
max_labels_per_type <- 8

dpi_png <- 600
pdf_device <- function(...) {
  if (capabilities("cairo")) grDevices::cairo_pdf(...) else grDevices::pdf(...)
}
save_gg <- function(plot_obj, base_name, width, height, dpi = 600) {
  ggsave(file.path(out_dir, paste0(base_name, ".pdf")), plot = plot_obj,
         width = width, height = height, units = "in", device = pdf_device)
  ggsave(file.path(out_dir, paste0(base_name, "_600dpi.png")), plot = plot_obj,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

std_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

read_phase_map <- function(files_named, kind = c("lnc", "mrna")) {
  kind <- match.arg(kind)
  phase_list <- lapply(names(files_named), function(nm) {
    f <- files_named[[nm]]
    if (!file.exists(f)) return(NULL)
    df <- read_csv(f, show_col_types = FALSE)
    if (kind == "lnc") {
      df %>% transmute(
        id = std_chr(gene_id),
        phase = case_when(
          str_detect(nm, "^Early") ~ "Early",
          str_detect(nm, "^Sustained") ~ "Sustained",
          str_detect(nm, "^Late") ~ "Late",
          TRUE ~ "Unassigned"
        ),
        source = if ("source" %in% names(df)) std_chr(source) else ifelse(str_detect(nm, "novel"), "novel", "annotated")
      )
    } else {
      name_col <- c("gene_name", "gene_symbol", "symbol")[c("gene_name", "gene_symbol", "symbol") %in% names(df)][1]
      df %>% transmute(
        id = std_chr(gene_id),
        gene_symbol = if (!is.na(name_col)) std_chr(.data[[name_col]]) else std_chr(gene_id),
        phase = case_when(
          nm == "Early" ~ "Early",
          nm == "Sustained" ~ "Sustained",
          nm == "Late" ~ "Late",
          TRUE ~ "Unassigned"
        )
      )
    }
  })
  bind_rows(phase_list) %>% distinct(id, .keep_all = TRUE)
}

# ---------- read and standardize DE tables ----------
if (!file.exists(intarna_file)) stop("Falta el archivo de IntaRNA: ", intarna_file)
if (!file.exists(mirna_target_file)) stop("Falta el archivo miRNA->mRNA: ", mirna_target_file)
if (!file.exists(lnc_de_file)) stop("Falta DE lncRNA: ", lnc_de_file)
if (!file.exists(mrna_de_file)) stop("Falta DE mRNA: ", mrna_de_file)
if (!file.exists(mirna_de_file)) stop("Falta DE miRNA: ", mirna_de_file)

lnc_phase_map <- read_phase_map(lnc_phase_files, kind = "lnc")
mrna_phase_map <- read_phase_map(mrna_phase_files, kind = "mrna")

lnc_de <- read_tsv(lnc_de_file, show_col_types = FALSE) %>%
  transmute(
    lncRNA = std_chr(gene_id),
    lnc_log2FC = as.numeric(log2FoldChange),
    lnc_padj = as.numeric(padj),
    lnc_direction = case_when(lnc_log2FC > 0 ~ "UP", lnc_log2FC < 0 ~ "DOWN", TRUE ~ "NS")
  ) %>%
  filter(!is.na(lnc_padj), lnc_padj < padj_cutoff, abs(lnc_log2FC) > absfc_cutoff)

mrna_de_raw <- read_tsv(mrna_de_file, show_col_types = FALSE)
mrna_name_col <- c("gene_name", "gene_symbol", "symbol")[c("gene_name", "gene_symbol", "symbol") %in% names(mrna_de_raw)][1]
mrna_de <- mrna_de_raw %>%
  transmute(
    target_gene = std_chr(gene_id),
    gene_symbol = if (!is.na(mrna_name_col)) std_chr(.data[[mrna_name_col]]) else std_chr(gene_id),
    mrna_log2FC = as.numeric(log2FoldChange),
    mrna_padj = as.numeric(padj),
    mrna_direction = case_when(mrna_log2FC > 0 ~ "UP", mrna_log2FC < 0 ~ "DOWN", TRUE ~ "NS")
  ) %>%
  filter(!is.na(mrna_padj), mrna_padj < padj_cutoff, abs(mrna_log2FC) > absfc_cutoff)

mirna_de_raw <- read_tsv(mirna_de_file, show_col_types = FALSE)
mirna_col <- c("miRNA", "mirna", "miRNA_id", "name")[c("miRNA", "mirna", "miRNA_id", "name") %in% names(mirna_de_raw)][1]
mirna_lfc_col <- c("log2FoldChange", "log2FC")[c("log2FoldChange", "log2FC") %in% names(mirna_de_raw)][1]
mirna_padj_col <- c("padj", "FDR", "adj.P.Val")[c("padj", "FDR", "adj.P.Val") %in% names(mirna_de_raw)][1]
mirna_de <- mirna_de_raw %>%
  transmute(
    miRNA = std_chr(.data[[mirna_col]]),
    miRNA_log2FC = as.numeric(.data[[mirna_lfc_col]]),
    miRNA_padj = as.numeric(.data[[mirna_padj_col]]),
    miRNA_direction = case_when(miRNA_log2FC > 0 ~ "UP", miRNA_log2FC < 0 ~ "DOWN", TRUE ~ "NS")
  ) %>%
  filter(!is.na(miRNA_padj), miRNA_padj < padj_cutoff, abs(miRNA_log2FC) > absfc_cutoff)

# ---------- read interactions ----------
intarna_raw <- read_csv(intarna_file, show_col_types = FALSE)
col_lnc <- c("lncRNA", "gene_id", "query_id", "query")[c("lncRNA", "gene_id", "query_id", "query") %in% names(intarna_raw)][1]
col_mir <- c("miRNA", "target_id", "target", "sRNA")[c("miRNA", "target_id", "target", "sRNA") %in% names(intarna_raw)][1]
col_energy <- c("energy", "E", "interaction_energy", "hybridE")[c("energy", "E", "interaction_energy", "hybridE") %in% names(intarna_raw)][1]

intarna_tbl <- intarna_raw %>%
  transmute(
    lncRNA = std_chr(.data[[col_lnc]]),
    miRNA = std_chr(.data[[col_mir]]),
    energy = as.numeric(.data[[col_energy]])
  ) %>%
  filter(!is.na(energy), energy <= energy_cutoff) %>%
  group_by(lncRNA, miRNA) %>%
  slice_min(order_by = energy, n = 1, with_ties = FALSE) %>%
  ungroup()

mir_target_raw <- read_csv(mirna_target_file, show_col_types = FALSE)
col_mir2 <- c("miRNA", "mirna")[c("miRNA", "mirna") %in% names(mir_target_raw)][1]
col_gene <- c("target_gene", "gene_id", "Target Gene", "mRNA")[c("target_gene", "gene_id", "Target Gene", "mRNA") %in% names(mir_target_raw)][1]
col_region <- c("binding_region", "region", "Site Region")[c("binding_region", "region", "Site Region") %in% names(mir_target_raw)][1]
col_score <- c("score", "Score", "binding_score")[c("score", "Score", "binding_score") %in% names(mir_target_raw)][1]

mir_target_tbl <- mir_target_raw %>%
  transmute(
    miRNA = std_chr(.data[[col_mir2]]),
    target_gene = std_chr(.data[[col_gene]]),
    binding_region = if (!is.na(col_region)) std_chr(.data[[col_region]]) else "unknown",
    target_score = if (!is.na(col_score)) as.numeric(.data[[col_score]]) else NA_real_
  ) %>%
  distinct(miRNA, target_gene, binding_region, .keep_all = TRUE)

# ---------- build coherent triplets ----------
triplets <- intarna_tbl %>%
  inner_join(mirna_de, by = "miRNA") %>%
  inner_join(lnc_de, by = "lncRNA") %>%
  inner_join(mir_target_tbl, by = "miRNA") %>%
  inner_join(mrna_de, by = "target_gene") %>%
  left_join(lnc_phase_map %>% rename(lncRNA = id, lnc_phase = phase, lnc_source = source), by = "lncRNA") %>%
  left_join(mrna_phase_map %>% rename(target_gene = id, mRNA_phase = phase), by = "target_gene") %>%
  mutate(
    ceRNA_pattern = case_when(
      lnc_direction == "UP" & miRNA_direction == "DOWN" & mrna_direction == "UP" ~ "lncUP-miRDOWN-mRNAUP",
      lnc_direction == "DOWN" & miRNA_direction == "UP" & mrna_direction == "DOWN" ~ "lncDOWN-miRUP-mRNADOWN",
      TRUE ~ "incoherent"
    ),
    coherence = ceRNA_pattern != "incoherent",
    triplet_score = (-energy) + abs(lnc_log2FC) + abs(miRNA_log2FC) + abs(mrna_log2FC)
  ) %>%
  filter(coherence)

if (nrow(triplets) == 0) stop("No quedaron tripletas coherentes con los filtros actuales.")

triplets <- triplets %>%
  arrange(desc(triplet_score), miRNA_padj, mrna_padj, lnc_padj)

write_csv(triplets, file.path(out_dir, paste0("Fig7_ceRNA_triplets_", context_label, ".csv")))

# ---------- summary panel ----------
summary_phase <- triplets %>%
  count(lnc_phase, ceRNA_pattern, name = "n_triplets") %>%
  mutate(lnc_phase = factor(lnc_phase, levels = c("Early", "Sustained", "Late", "Unassigned")))

write_csv(summary_phase, file.path(out_dir, paste0("SUPP_Fig7_ceRNA_summary_by_lnc_phase_", context_label, ".csv")))

p_sum <- ggplot(summary_phase, aes(x = lnc_phase, y = n_triplets, fill = ceRNA_pattern)) +
  geom_col(position = "stack") +
  labs(
    x = NULL,
    y = "# coherent ceRNA triplets",
    fill = "Pattern",
    title = paste0("ceRNA-supporting triplets across lncRNA classes (", context_label, ")")
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

save_gg(p_sum, paste0("Fig7A_ceRNA_summary_by_lnc_phase_", context_label), 7.2, 4.8, dpi_png)

# ---------- tripartite network-like plot without extra packages ----------
plot_triplets <- triplets %>%
  slice_head(n = max_triplets_for_network)

lnc_nodes <- plot_triplets %>% distinct(name = lncRNA, group = "lncRNA", phase = lnc_phase)
mir_nodes <- plot_triplets %>% distinct(name = miRNA, group = "miRNA", phase = "miRNA")
mrna_nodes <- plot_triplets %>% distinct(name = target_gene, label = gene_symbol, group = "mRNA", phase = mRNA_phase)

node_tbl <- bind_rows(
  lnc_nodes %>% mutate(label = name),
  mir_nodes %>% mutate(label = name),
  mrna_nodes %>% select(name, label, group, phase)
) %>% distinct(name, .keep_all = TRUE)

node_tbl <- node_tbl %>%
  mutate(group = factor(group, levels = c("lncRNA", "miRNA", "mRNA")))

# simple manual layout: three vertical columns
node_tbl <- node_tbl %>%
  group_by(group) %>%
  arrange(phase, name, .by_group = TRUE) %>%
  mutate(y = rev(seq_len(n())),
         x = c(1, 2, 3)[match(as.character(group), c("lncRNA", "miRNA", "mRNA"))]) %>%
  ungroup()

edges_lnc_mir <- plot_triplets %>%
  distinct(from = lncRNA, to = miRNA, edge_type = "lncRNA-miRNA")
edges_mir_mrna <- plot_triplets %>%
  distinct(from = miRNA, to = target_gene, edge_type = "miRNA-mRNA")
edge_tbl <- bind_rows(edges_lnc_mir, edges_mir_mrna) %>%
  left_join(node_tbl %>% select(name, x_from = x, y_from = y), by = c("from" = "name")) %>%
  left_join(node_tbl %>% select(name, x_to = x, y_to = y), by = c("to" = "name"))

write_csv(node_tbl, file.path(out_dir, paste0("SUPP_Fig7_ceRNA_nodes_", context_label, ".csv")))
write_csv(edge_tbl, file.path(out_dir, paste0("SUPP_Fig7_ceRNA_edges_", context_label, ".csv")))

# label only the highest-degree nodes in each group
node_degree <- bind_rows(
  edge_tbl %>% count(from, name = "deg") %>% rename(name = from),
  edge_tbl %>% count(to, name = "deg") %>% rename(name = to)
) %>% group_by(name) %>% summarise(deg = sum(deg), .groups = "drop")

label_nodes <- node_tbl %>%
  left_join(node_degree, by = "name") %>%
  mutate(deg = dplyr::coalesce(deg, 0L)) %>%
  group_by(group) %>%
  arrange(desc(deg), .by_group = TRUE) %>%
  slice_head(n = max_labels_per_type) %>%
  ungroup()

p_net <- ggplot() +
  geom_segment(data = edge_tbl,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = edge_type),
               alpha = 0.35, linewidth = 0.35) +
  geom_point(data = node_tbl,
             aes(x = x, y = y, fill = group),
             shape = 21, color = "grey20", size = 3.3, stroke = 0.2) +
  ggrepel::geom_text_repel(data = label_nodes,
                           aes(x = x, y = y, label = ifelse(group == "mRNA", label, name)),
                           size = 2.8, max.overlaps = Inf,
                           box.padding = 0.2, point.padding = 0.15,
                           min.segment.length = 0) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("lncRNAs", "miRNAs", "mRNAs")) +
  scale_color_manual(values = c("lncRNA-miRNA" = "#9CA3AF", "miRNA-mRNA" = "#6B7280")) +
  scale_fill_manual(values = c("lncRNA" = "#F59E0B", "miRNA" = "#10B981", "mRNA" = "#3B82F6")) +
  labs(
    x = NULL, y = NULL,
    color = "Edge type", fill = "Node type",
    title = paste0("Candidate ceRNA network linking lncRNAs, miRNAs and mRNAs (", context_label, ")")
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

save_gg(p_net, paste0("Fig7B_ceRNA_network_", context_label), 11, 8.5, dpi_png)

message("Figure 7 ceRNA scaffold listo en: ", out_dir)
