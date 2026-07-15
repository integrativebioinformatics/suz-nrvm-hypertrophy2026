suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

out_dir <- "Paper_miRNA_mRNA_AllCombos_Dotplots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha_mirna <- 0.05
alpha_mrna  <- 0.05
dpi_png <- 600
top_mirnas_for_plot <- 25
top_pathways_for_plot <- 18

# Optional: restrict pathways to those exported from Fig5/6
pathways_selected_for_TF <- NULL  # e.g. "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT/pathways_selected_for_TF.csv"

# Phase sets for mRNAs (optional; if absent, targets tagged as Unclassified)
early_mrna_path <- "early_mRNA.csv"
sust_mrna_path  <- "sustained_mRNA.csv"
late_mrna_path  <- "late_mRNA.csv"

# KEGG split (same used before)
kegg_only6_path  <- "KEGG_6h_vs_24h_Only6h_long_format_SYMBOLS.txt"
kegg_only24_path <- "KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt"
kegg_shared_path <- "KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt"

# mRNA DE files
mrna_de_paths <- list(
  "6h"  = "NE_vs_Ctrl_6h_mRNA_DE.txt",
  "24h" = "NE_vs_Ctrl_24h_mRNA_DE.txt"
)

# miRNA DE files (you already provided 24h; add 6h when you have them)
mirna_de_paths <- list(
  "24h" = list(
    "Up"   = "DE_miRNA_Up_NE24.txt",
    "Down" = "DE_miRNA_Down_NE24.txt"
  ),
  "6h" = list(
    "Up"   = "DE_miRNA_Up_NE6.txt",    # <-- update when available
    "Down" = "DE_miRNA_Down_NE6.txt"   # <-- update when available
  )
)

# miRWalk target files per combo and region (fill missing when you have them)
# Convention: each combo provides 3 paths: 3UTR, CDS, 5UTR
mirwalk_targets <- list(
  "24h" = list(
    "Up" = list(
      "3UTR" = "miRWalk_miRNA_Targets_3UTR_UP24.csv",
      "CDS"  = "miRWalk_miRNA_Targets_CDS_UP24.csv",
      "5UTR" = "miRWalk_miRNA_Targets_5UTR_UP24.csv"
    ),
    "Down" = list(
      "3UTR" = "miRWalk_miRNA_Targets_3UTR_DOWN24.csv",  # <-- update when available
      "CDS"  = "miRWalk_miRNA_Targets_CDS_DOWN24.csv",
      "5UTR" = "miRWalk_miRNA_Targets_5UTR_DOWN24.csv"
    )
  ),
  "6h" = list(
    "Up" = list(
      "3UTR" = "miRWalk_miRNA_Targets_3UTR_UP6.csv",      # <-- update when available
      "CDS"  = "miRWalk_miRNA_Targets_CDS_UP6.csv",
      "5UTR" = "miRWalk_miRNA_Targets_5UTR_UP6.csv"
    ),
    "Down" = list(
      "3UTR" = "miRWalk_miRNA_Targets_3UTR_DOWN6.csv",    # <-- update when available
      "CDS"  = "miRWalk_miRNA_Targets_CDS_DOWN6.csv",
      "5UTR" = "miRWalk_miRNA_Targets_5UTR_DOWN6.csv"
    )
  )
)

detect_col <- function(nms, patterns) {
  hit <- which(str_detect(nms, regex(paste(patterns, collapse="|"), ignore_case = TRUE)))
  if (length(hit) == 0) return(NA_character_)
  nms[hit[1]]
}

read_mirna_de <- function(path, alpha = 0.05) {
  df <- suppressWarnings(read_tsv(path, show_col_types = FALSE, col_names = TRUE))
  nms <- names(df)
  mir_col <- if ("miRNA_id" %in% nms) "miRNA_id" else detect_col(nms, c("^mirna", "mirnaid", "miRNA"))
  lfc_col <- if ("log2FoldChange" %in% nms) "log2FoldChange" else detect_col(nms, c("log2fold", "log2FC"))
  padj_col <- if ("padj" %in% nms) "padj" else detect_col(nms, c("padj", "fdr", "adj"))
  status_col <- if ("status" %in% nms) "status" else detect_col(nms, c("status", "direction"))
  
  out <- df %>%
    transmute(
      miRNA = as.character(.data[[mir_col]]),
      log2FC = suppressWarnings(as.numeric(.data[[lfc_col]])),
      padj = suppressWarnings(as.numeric(.data[[padj_col]])),
      status = if (!is.na(status_col)) as.character(.data[[status_col]]) else NA_character_
    ) %>%
    filter(!is.na(miRNA), miRNA != "") %>%
    mutate(
      miRNA_key = str_to_upper(str_trim(miRNA)),
      de_sig = !is.na(padj) & padj < alpha
    ) %>%
    group_by(miRNA_key) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup()
  
  out
}

read_mirwalk_region <- function(path, region_label) {
  df <- read_csv(path, show_col_types = FALSE)
  nms <- names(df)
  mir_col <- detect_col(nms, c("^mirna", "mirnaid", "miRNA"))
  gene_col <- detect_col(nms, c("genesymbol", "gene_symbol", "target", "symbol", "gene"))
  if (is.na(mir_col) || is.na(gene_col)) {
    stop("No detecté columnas miRNA/genes en miRWalk: ", path, "\nCols: ", paste(nms, collapse = ", "))
  }
  df %>%
    transmute(
      miRNA = as.character(.data[[mir_col]]),
      gene_symbol = as.character(.data[[gene_col]]),
      region = region_label
    ) %>%
    filter(!is.na(miRNA), miRNA != "", !is.na(gene_symbol), gene_symbol != "") %>%
    mutate(
      miRNA_key = str_to_upper(str_trim(miRNA)),
      gene_key  = str_to_upper(str_trim(gene_symbol))
    ) %>%
    distinct(miRNA_key, gene_key, region, .keep_all = TRUE)
}

read_mrna_de <- function(path, alpha = 0.05) {
  df <- suppressWarnings(read_tsv(path, show_col_types = FALSE, col_names = FALSE))
  # heuristic for your format: last col direction, gene_name near col 8
  out <- df %>%
    transmute(
      log2FC = suppressWarnings(as.numeric(X3)),
      padj = suppressWarnings(as.numeric(X7)),
      gene_name = as.character(X8),
      direction = as.character(X11)
    ) %>%
    mutate(
      gene_key = str_to_upper(str_trim(gene_name)),
      de_sig = !is.na(padj) & padj < alpha
    ) %>%
    filter(!is.na(gene_key), gene_key != "") %>%
    group_by(gene_key) %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice(1) %>%
    ungroup()
  out
}

read_phase_map_optional <- function(early_path, sust_path, late_path) {
  safe_read <- function(p, label) {
    if (!file.exists(p)) return(tibble(gene_key = character(), phase = character()))
    x <- suppressWarnings(read_csv(p, show_col_types = FALSE))
    if (nrow(x) == 0) return(tibble(gene_key = character(), phase = character()))
    nms <- names(x)
    gcol <- if ("gene_name" %in% nms) "gene_name" else detect_col(nms, c("gene", "symbol", "name"))
    if (is.na(gcol)) return(tibble(gene_key = character(), phase = character()))
    x %>%
      transmute(gene_key = str_to_upper(str_trim(as.character(.data[[gcol]]))), phase = label) %>%
      filter(!is.na(gene_key), gene_key != "")
  }
  early <- safe_read(early_path, "Early")
  sust  <- safe_read(sust_path,  "Sustained")
  late  <- safe_read(late_path,  "Late")
  bind_rows(sust, early, late) %>% distinct(gene_key, .keep_all = TRUE)
}

read_kegg_map <- function() {
  only6 <- read_tsv(kegg_only6_path, show_col_types = FALSE) %>%
    transmute(KEGG_group = "Only6h", ID, Description,
              padj = as.numeric(`p.adjust`),
              gene_key = str_to_upper(str_trim(as.character(gene_symbol))))
  only24 <- read_tsv(kegg_only24_path, show_col_types = FALSE) %>%
    transmute(KEGG_group = "Only24h", ID, Description,
              padj = as.numeric(`p.adjust`),
              gene_key = str_to_upper(str_trim(as.character(gene_symbol))))
  shared <- read_tsv(kegg_shared_path, show_col_types = FALSE) %>%
    transmute(KEGG_group = "Shared", ID, Description,
              padj = pmin(as.numeric(`p.adjust_6h`), as.numeric(`p.adjust_24h`), na.rm = TRUE),
              gene_key = str_to_upper(str_trim(as.character(gene_symbol))))
  bind_rows(only6, shared, only24) %>%
    filter(!is.na(gene_key), gene_key != "") %>%
    distinct(KEGG_group, ID, Description, gene_key, padj)
}

kegg_map <- read_kegg_map()
if (!is.null(pathways_selected_for_TF) && file.exists(pathways_selected_for_TF)) {
  sel <- read_csv(pathways_selected_for_TF, show_col_types = FALSE) %>%
    select(KEGG_group, ID, Description, sign_class) %>% distinct()
  kegg_map <- kegg_map %>% inner_join(sel, by = c("KEGG_group","ID","Description"))
} else {
  top_paths <- kegg_map %>%
    group_by(KEGG_group, Description) %>%
    summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(padj_min) %>%
    group_by(KEGG_group) %>%
    slice_head(n = top_pathways_for_plot) %>%
    ungroup() %>%
    pull(Description) %>% unique()
  kegg_map <- kegg_map %>% filter(Description %in% top_paths) %>% mutate(sign_class = "Auto")
}

phase_map <- read_phase_map_optional(early_mrna_path, sust_mrna_path, late_mrna_path)

get_expected_mrna_dir <- function(mir_dir) if (mir_dir == "Up") "Down" else "Up"

mrna_dir_flag <- function(direction_str) {
  if (is.na(direction_str)) return(NA_character_)
  if (str_detect(direction_str, regex("DOWN", TRUE))) return("Down")
  if (str_detect(direction_str, regex("UP", TRUE))) return("Up")
  NA_character_
}

run_combo <- function(timepoint, mir_dir) {
  mr_dir <- get_expected_mrna_dir(mir_dir)
  
  mir_de_path <- mirna_de_paths[[timepoint]][[mir_dir]]
  mr_de_path  <- mrna_de_paths[[timepoint]]
  
  targ_paths <- mirwalk_targets[[timepoint]][[mir_dir]]
  needed <- c(mir_de_path, mr_de_path, targ_paths[["3UTR"]], targ_paths[["CDS"]], targ_paths[["5UTR"]])
  missing <- needed[!file.exists(needed)]
  tag <- paste0(timepoint, "_miRNA", mir_dir, "_mRNA", mr_dir)
  
  if (length(missing) > 0) {
    message("SKIP ", tag, " (faltan archivos): ", paste(basename(missing), collapse = ", "))
    return(NULL)
  }
  
  mir_de <- read_mirna_de(mir_de_path, alpha = alpha_mirna) %>%
    filter(de_sig) %>%
    mutate(time = timepoint, miRNA_dir = mir_dir)
  
  mr_de <- read_mrna_de(mr_de_path, alpha = alpha_mrna) %>%
    mutate(mRNA_dir = vapply(direction, mrna_dir_flag, character(1)))
  
  # fallback si direction está ausente: usar signo de log2FC
  mr_de <- mr_de %>%
    mutate(mRNA_dir = ifelse(is.na(mRNA_dir) & !is.na(log2FC) & log2FC < 0, "Down", mRNA_dir),
           mRNA_dir = ifelse(is.na(mRNA_dir) & !is.na(log2FC) & log2FC > 0, "Up",   mRNA_dir))
  
  mr_de <- mr_de %>%
    filter(de_sig, mRNA_dir == mr_dir)
  
  t3 <- read_mirwalk_region(targ_paths[["3UTR"]], "3UTR")
  tc <- read_mirwalk_region(targ_paths[["CDS"]],  "CDS")
  t5 <- read_mirwalk_region(targ_paths[["5UTR"]], "5UTR")
  
  targets_collapsed <- bind_rows(t3, tc, t5) %>%
    group_by(miRNA_key, gene_key) %>%
    summarize(region_support = paste(sort(unique(region)), collapse = ";"), .groups = "drop")
  
  edges <- mir_de %>%
    select(miRNA, miRNA_key, miRNA_log2FC = log2FC, miRNA_padj = padj, time, miRNA_dir) %>%
    inner_join(targets_collapsed, by = "miRNA_key") %>%
    inner_join(mr_de %>% select(gene_key, mRNA_log2FC = log2FC, mRNA_padj = padj, mRNA_dir), by = "gene_key") %>%
    left_join(phase_map, by = "gene_key") %>%
    mutate(phase = replace_na(phase, "Unclassified"),
           combo = tag)
  
  edges_kegg <- edges %>%
    inner_join(kegg_map, by = "gene_key", relationship = "many-to-many") %>%
    mutate(KEGG_group = factor(KEGG_group, levels = c("Only6h","Shared","Only24h")))
  
  write_csv(edges_kegg, file.path(out_dir, paste0("SUPP_edges_", tag, "_KEGG_phase.csv")))
  
  counts <- edges_kegg %>%
    group_by(KEGG_group, Description, sign_class, miRNA) %>%
    summarize(n_targets = n_distinct(gene_key), .groups = "drop")
  
  top_mir <- counts %>%
    group_by(miRNA) %>% summarize(total = sum(n_targets), .groups = "drop") %>%
    arrange(desc(total)) %>% slice_head(n = top_mirnas_for_plot) %>% pull(miRNA)
  
  plot_df <- counts %>%
    filter(miRNA %in% top_mir) %>%
    left_join(kegg_map %>% group_by(KEGG_group, Description) %>%
                summarize(padj_min = min(padj, na.rm = TRUE), .groups = "drop"),
              by = c("KEGG_group","Description")) %>%
    mutate(
      miRNA = fct_reorder(miRNA, n_targets, .fun = sum, .desc = TRUE),
      Description = fct_reorder(Description, padj_min, .desc = FALSE)
    )
  
  p <- ggplot(plot_df, aes(x = Description, y = miRNA)) +
    geom_point(aes(size = n_targets), alpha = 0.85) +
    facet_grid(sign_class ~ KEGG_group, scales = "free_x", space = "free_x") +
    labs(
      x = "KEGG pathway (targets DE en este contexto)",
      y = paste0("miRNAs ", mir_dir, " (", timepoint, ")"),
      size = "# targets",
      title = paste0(tag, ": miRNA", mir_dir, " → mRNA", mr_dir, " mapped to KEGG")
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  
  pdf_path <- file.path(out_dir, paste0("Fig_miRNA_", tag, "_dotplot.pdf"))
  png_path <- file.path(out_dir, paste0("Fig_miRNA_", tag, "_dotplot.png"))
  
  ggsave(pdf_path, p, width = 13.5, height = 9)
  ggsave(png_path, p, width = 13.5, height = 9, dpi = dpi_png)
  
  rank <- counts %>%
    group_by(miRNA) %>%
    summarize(total_targets = sum(n_targets), n_pathways_hit = sum(n_targets > 0), .groups = "drop") %>%
    arrange(desc(total_targets))
  
  write_csv(rank, file.path(out_dir, paste0("SUPP_miRNA_rank_", tag, ".csv")))
  
  message("OK ", tag, " -> ", basename(pdf_path), " + ", basename(png_path))
  invisible(list(edges = edges_kegg, plot = p))
}

# Run all four combos (will skip those missing files)
results <- list(
  run_combo("24h", "Up"),
  run_combo("24h", "Down"),
  run_combo("6h",  "Up"),
  run_combo("6h",  "Down")
)

# Optional: combined master table across combos that actually ran
all_edges <- bind_rows(lapply(results, function(x) if (!is.null(x)) x$edges else NULL))
if (nrow(all_edges) > 0) {
  write_csv(all_edges, file.path(out_dir, "SUPP_edges_ALL_COMBOS_KEGG_phase.csv"))
  message("Wrote combined edges: SUPP_edges_ALL_COMBOS_KEGG_phase.csv")
} else {
  message("No combos produced edges (likely missing target/DE files).")
}

####FASTA 4 miRNAs
library(Biostrings)
library(readr)
library(dplyr)
library(stringr)

# Files

mirbase_fasta <- "mature.fa"
out_dir <- "miRNA_FASTA_by_condition"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

input_files <- c(
  "miRNAs_UP_24hr.txt",
  "miRNAs_DOWN_24hr.txt",
  "miRNAs_UP_6hr.txt",
  "miRNAs_DOWN_6hr.txt"
)

# Read miRBase mature.fa
fa <- readRNAStringSet(mirbase_fasta)

# miRBase headers: first word is the miRNA ID
fa_ids <- word(names(fa), 1)

fa_tbl <- tibble(
  miRNA = fa_ids,
  sequence = as.character(fa)
) %>%
  filter(str_detect(miRNA, "^rno-"))

cat("Rat mature miRNAs found in mature.fa:", nrow(fa_tbl), "\n")

# Function to process one txt file
process_mirna_file <- function(input_txt, fa_tbl, out_dir) {
  
  cat("\nProcessing:", input_txt, "\n")
  
  # Read one miRNA per line
  mirna_list <- read_lines(input_txt)
  
  mirna_tbl <- tibble(miRNA = mirna_list) %>%
    mutate(
      miRNA = str_trim(miRNA),
      miRNA = str_replace(miRNA, "^>", "")
    ) %>%
    filter(!is.na(miRNA), miRNA != "") %>%
    distinct()
  
  # Keep only rat IDs
  mirna_tbl <- mirna_tbl %>%
    filter(str_detect(miRNA, "^rno-"))
  
  # Match against mature.fa rat entries
  matched_tbl <- mirna_tbl %>%
    left_join(fa_tbl, by = "miRNA")
  
  found_tbl <- matched_tbl %>%
    filter(!is.na(sequence))
  
  missing_tbl <- matched_tbl %>%
    filter(is.na(sequence))
  
  # Base name for outputs
  base_name <- tools::file_path_sans_ext(basename(input_txt))
  
  # Save missing table
  write_csv(missing_tbl, file.path(out_dir, paste0(base_name, "_not_found.csv")))
  
  # Save matched table
  write_csv(found_tbl, file.path(out_dir, paste0(base_name, "_matched.csv")))
  
  # Save FASTA
  if (nrow(found_tbl) > 0) {
    out_fa <- RNAStringSet(found_tbl$sequence)
    names(out_fa) <- found_tbl$miRNA
    writeXStringSet(out_fa, file.path(out_dir, paste0(base_name, ".fa")))
    cat("Matched:", nrow(found_tbl), "\n")
    cat("Missing:", nrow(missing_tbl), "\n")
    cat("FASTA written:", file.path(out_dir, paste0(base_name, ".fa")), "\n")
  } else {
    cat("No miRNAs matched for", input_txt, "\n")
  }
}

# Run all 4 files
for (f in input_files) {
  process_mirna_file(f, fa_tbl, out_dir)
}

cat("\nDone.\n")


###FEELnc
library(tidyverse)

feelnc_file <- "lncRNA_classes_CNE624_rn8_8_annotated_1kb_100kb.txt"
mrna6_file  <- "NE_vs_Ctrl_6h_mRNA_DE.txt"
mrna24_file <- "NE_vs_Ctrl_24h_mRNA_DE.txt"
out_dir <- "FEELnc_summary_panels"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

x <- read_tsv(feelnc_file, show_col_types = FALSE) %>%
  filter(isBest == 1) %>%
  distinct() %>%
  mutate(
    source = if_else(str_detect(lncRNA_gene, "^MSTRG"), "novel", "annotated"),
    class_simple = case_when(
      type == "intergenic" & location == "upstream" ~ "intergenic_upstream",
      type == "intergenic" & location == "downstream" ~ "intergenic_downstream",
      type == "intergenic" & subtype %in% c("divergent", "convergent", "same_strand", "unknow strand(s)") ~ paste0("intergenic_", subtype),
      direction == "sense" & location == "exonic" ~ "sense_exonic",
      direction == "sense" & location == "intronic" ~ "sense_intronic",
      direction == "antisense" & location == "exonic" ~ "antisense_exonic",
      direction == "antisense" & location == "intronic" ~ "antisense_intronic",
      TRUE ~ paste(type, direction, location, sep = "_")
    )
  )

class_counts <- x %>%
  count(source, class_simple, sort = TRUE)

write_tsv(class_counts, file.path(out_dir, "FEELnc_class_counts_best_only.tsv"))
write_tsv(x, file.path(out_dir, "FEELnc_best_only_full.tsv"))

p <- class_counts %>%
  mutate(class_simple = forcats::fct_reorder(class_simple, n, .fun = sum, .desc = TRUE)) %>%
  ggplot(aes(x = source, y = n, fill = class_simple)) +
  geom_col(width = 0.72) +
  labs(x = NULL, y = "# lncRNAs", fill = "FEELnc class") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "SuppFig3A_FEELnc_classification.png"), p, width = 7.2, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "SuppFig3A_FEELnc_classification.pdf"), p, width = 7.2, height = 5.5)

mrna_map <- bind_rows(
  read_tsv(mrna6_file, show_col_types = FALSE) %>% select(gene_id, gene_name),
  read_tsv(mrna24_file, show_col_types = FALSE) %>% select(gene_id, gene_name)
) %>% distinct()

candidate_lnc <- c(
  "ENSRNOG00000090514", "MSTRG.8518", "MSTRG.3251", "MSTRG.17518",
  "ENSRNOG00000071598", "ENSRNOG00000085965", "MSTRG.8263",
  "ENSRNOG00000083943", "MSTRG.13291", "ENSRNOG00000090247",
  "ENSRNOG00000073751"
)

cand_tbl <- x %>%
  filter(lncRNA_gene %in% candidate_lnc) %>%
  left_join(mrna_map, by = c("partnerRNA_gene" = "gene_id")) %>%
  mutate(partner_label = coalesce(gene_name, partnerRNA_gene)) %>%
  arrange(lncRNA_gene, distance, partner_label)

write_tsv(cand_tbl, file.path(out_dir, "FEELnc_selected_module_candidates_best_context.tsv"))
