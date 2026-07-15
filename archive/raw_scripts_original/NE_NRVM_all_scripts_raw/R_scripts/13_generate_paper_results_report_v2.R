suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

base_dir <- "."
out_dir  <- file.path(base_dir, "Paper_Results_Report_v2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bulk_dir   <- file.path(base_dir, "DESeq2_Multifactorial_Results_fixed")
fig3_dir   <- file.path(base_dir, "Paper_Fig3_FINAL")
fig45_dir  <- file.path(base_dir, "Paper_Fig4_Fig5_FINAL_ONE_SCRIPT")
fig4r_dir  <- file.path(base_dir, "Paper_Fig4_REFINED")
fig6_dir   <- file.path(base_dir, "Fig6_TF_panels_out_FINAL")
fig7_dir   <- file.path(base_dir, "Paper_miRNA_mRNA_AllCombos_Dotplots")
fig8_dir   <- file.path(base_dir, "Paper_ceRNA_candidates")
circ_dir   <- file.path(base_dir, "Paper_circRNA_FINAL")

alpha_mrna  <- 0.05
lfc_mrna    <- 1.0
alpha_lnc   <- 0.05
lfc_lnc     <- 0.5
alpha_mirna <- 0.05
lfc_mirna   <- 0.5
alpha_circ  <- 0.05
lfc_circ    <- 0.5

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
std_chr <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
fmt_num <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

find_file <- function(filename, dirs = c(base_dir, bulk_dir, fig3_dir, fig45_dir, fig4r_dir, fig6_dir, fig7_dir, fig8_dir, circ_dir)) {
  hits <- unlist(lapply(unique(dirs), function(d) {
    if (!dir.exists(d)) return(character(0))
    allf <- list.files(d, recursive = TRUE, full.names = TRUE)
    allf[basename(allf) == filename]
  }), use.names = FALSE)
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

find_any <- function(filenames, dirs = c(base_dir, bulk_dir, fig3_dir, fig45_dir, fig4r_dir, fig6_dir, fig7_dir, fig8_dir, circ_dir)) {
  for (f in filenames) {
    p <- find_file(f, dirs = dirs)
    if (!is.na(p)) return(p)
  }
  NA_character_
}

list_files_matching <- function(pattern, dirs = c(base_dir, fig7_dir, fig8_dir, circ_dir)) {
  hits <- unlist(lapply(unique(dirs), function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  unique(hits[file.exists(hits)])
}

safe_read <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  out <- tryCatch({
    if (ext == "csv") {
      read_csv(path, show_col_types = FALSE)
    } else {
      read_tsv(path, show_col_types = FALSE)
    }
  }, error = function(e1) {
    tryCatch(read_delim(path, delim = ";", show_col_types = FALSE), error = function(e2) {
      tryCatch(read_delim(path, delim = ",", show_col_types = FALSE), error = function(e3) NULL)
    })
  })
  out
}

safe_write_lines <- function(lines, path) writeLines(enc2utf8(lines), con = path, useBytes = TRUE)

write_df_txt <- function(df, path, title = NULL) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  if (!is.null(title)) {
    writeLines(title, con)
    writeLines("", con)
  }
  if (is.null(df) || nrow(df) == 0) {
    writeLines("No rows available.", con)
  } else {
    utils::write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

append_lines <- function(lines, ...) c(lines, unlist(list(...), use.names = FALSE))

collapse_labels <- function(x, n = 40) {
  x <- unique(std_chr(x))
  x <- x[x != ""]
  if (length(x) == 0) return("none")
  paste(head(x, n), collapse = ", ")
}

fmt_feature <- function(label, log2fc = NA_real_, padj = NA_real_) {
  bits <- c(label)
  if (!is.na(log2fc)) bits <- c(bits, paste0("log2FC=", fmt_num(log2fc, 3)))
  if (!is.na(padj)) bits <- c(bits, paste0("FDR=", fmt_num(padj, 3)))
  paste0(bits[1], if (length(bits) > 1) paste0(" (", paste(bits[-1], collapse = ", "), ")") else "")
}

pick_label <- function(df, id_col = "feature_id", label_cols = c("gene_symbol", "gene_name", "label", "miRNA", "TF_key", "Description", "ncRNA_label", "gene_id", "lncRNA", "circ_id", "gene_key")) {
  hit <- label_cols[label_cols %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) hit <- id_col
  hit
}

format_top_features <- function(df, n = 12, label_col = NULL, fc_col = NULL, padj_col = NULL) {
  if (is.null(df) || nrow(df) == 0) return("none")
  if (is.null(label_col)) label_col <- pick_label(df)
  if (is.null(fc_col)) fc_col <- first_existing_col(df, c("log2FC", "log2FoldChange", "NES", "mean_absFC", "total_targets", "n_targets", "n_edges", "E", "n_priority_triplets"))
  if (is.null(padj_col)) padj_col <- first_existing_col(df, c("padj", "padj_min", "padj_mRNA", "padj_miRNA", "pvalue"))
  df <- df %>% slice_head(n = n)
  vals <- pmap_chr(list(df[[label_col]], if (!is.na(fc_col)) df[[fc_col]] else rep(NA_real_, nrow(df)), if (!is.na(padj_col)) df[[padj_col]] else rep(NA_real_, nrow(df))), fmt_feature)
  paste(vals, collapse = "; ")
}

standardize_de <- function(df, entity = c("mRNA", "lncRNA", "miRNA", "circRNA"), timepoint = NA_character_) {
  entity <- match.arg(entity)
  if (is.null(df) || nrow(df) == 0) return(tibble())
  id_col <- first_existing_col(df, c("gene_id", "GeneID", "ENSEMBL", "miRNA", "mirna", "feature_id", "circ_id", "circRNA_ID", "circRNA_ID2"))
  name_col <- first_existing_col(df, c("gene_symbol", "gene_name", "symbol", "Gene.Symbol", "miRNA", "mirna", "circRNA_ID2", "host_label"))
  lfc_col <- first_existing_col(df, c("log2FoldChange", "log2FC"))
  padj_col <- first_existing_col(df, c("padj", "FDR", "adj.P.Val"))
  p_col <- first_existing_col(df, c("pvalue", "P.Value", "pvalue_raw"))
  dir_col <- first_existing_col(df, c("direction", "Direction", "regulations", "status"))
  source_col <- first_existing_col(df, c("source", "lnc_source"))
  bio_col <- first_existing_col(df, c("gene_biotype", "biotype", "lnc_biotype"))
  out <- df %>%
    transmute(
      feature_id = std_chr(.data[[id_col]]),
      label_raw = if (!is.na(name_col)) std_chr(.data[[name_col]]) else "",
      label = case_when(
        entity == "mRNA" & label_raw != "" & str_to_upper(label_raw) != str_to_upper(feature_id) ~ label_raw,
        entity == "miRNA" & label_raw != "" ~ label_raw,
        TRUE ~ feature_id
      ),
      log2FC = if (!is.na(lfc_col)) suppressWarnings(as.numeric(.data[[lfc_col]])) else NA_real_,
      padj = if (!is.na(padj_col)) suppressWarnings(as.numeric(.data[[padj_col]])) else NA_real_,
      pvalue = if (!is.na(p_col)) suppressWarnings(as.numeric(.data[[p_col]])) else NA_real_,
      direction = if (!is.na(dir_col)) str_to_upper(std_chr(.data[[dir_col]])) else case_when(
        log2FC > 0 ~ "UP",
        log2FC < 0 ~ "DOWN",
        TRUE ~ "NS"
      ),
      source = if (!is.na(source_col)) std_chr(.data[[source_col]]) else NA_character_,
      biotype = if (!is.na(bio_col)) std_chr(.data[[bio_col]]) else NA_character_,
      timepoint = timepoint
    ) %>%
    mutate(direction = recode(direction, "UPREGULATED" = "UP", "DOWNREGULATED" = "DOWN", "UP " = "UP", "DOWN " = "DOWN")) %>%
    group_by(feature_id, timepoint) %>%
    arrange(padj, desc(abs(log2FC)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  out
}

sig_subset <- function(df, alpha = 0.05, lfc_cut = 0, use_padj = TRUE) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  if (use_padj && "padj" %in% names(df)) {
    df %>% filter(!is.na(padj), padj <= alpha, abs(log2FC) >= lfc_cut)
  } else {
    df %>% filter(!is.na(pvalue), pvalue <= alpha, abs(log2FC) >= lfc_cut)
  }
}

read_phase_table <- function(path, timepoint = NA_character_, entity = c("mRNA","lncRNA","circRNA"), phase_name = NULL) {
  entity <- match.arg(entity)
  df <- safe_read(path)
  if (is.null(df) || nrow(df) == 0) return(tibble())
  id_col <- first_existing_col(df, c("gene_id","gene_symbol","gene_name","circ_id","circRNA_ID","circRNA_ID2","lncRNA","feature_id"))
  label_col <- first_existing_col(df, c("gene_symbol","gene_name","circRNA_ID2","gene_id","circ_id"))
  dir_col <- first_existing_col(df, c("class_direction","direction","Direction"))
  log2fc6_col <- first_existing_col(df, c("log2FC_6h"))
  log2fc24_col <- first_existing_col(df, c("log2FC_24h"))
  padj_col <- first_existing_col(df, c("padj_min","padj"))
  out <- df %>% transmute(
    feature_id = std_chr(.data[[id_col]]),
    label = if (!is.na(label_col)) std_chr(.data[[label_col]]) else std_chr(.data[[id_col]]),
    phase = phase_name %||% first_existing_col(df, c("phase")),
    class_direction = if (!is.na(dir_col)) std_chr(.data[[dir_col]]) else NA_character_,
    log2FC_6h = if (!is.na(log2fc6_col)) suppressWarnings(as.numeric(.data[[log2fc6_col]])) else NA_real_,
    log2FC_24h = if (!is.na(log2fc24_col)) suppressWarnings(as.numeric(.data[[log2fc24_col]])) else NA_real_,
    padj_min = if (!is.na(padj_col)) suppressWarnings(as.numeric(.data[[padj_col]])) else NA_real_,
    timepoint = timepoint,
    entity = entity
  )
  if (length(out$phase) == 1 && is.na(out$phase[1])) out$phase <- phase_name
  out
}

master_lines <- c(
  "Paper-wide results report",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "This report summarizes the outputs already generated by the figure workflow.",
  "It is designed to support Results writing, figure legends, Discussion drafting, and candidate prioritization across molecular layers.",
  ""
)

fig2_lines <- c("Figure 2. Differential mRNA expression and overlap", strrep("=", 54), "")
mrna6 <- standardize_de(safe_read(find_file("NE_vs_Ctrl_6h_mRNA_DE.txt")), "mRNA", "6h")
mrna24 <- standardize_de(safe_read(find_file("NE_vs_Ctrl_24h_mRNA_DE.txt")), "mRNA", "24h")
mrna_all <- bind_rows(mrna6, mrna24)
mrna_sig <- sig_subset(mrna_all, alpha_mrna, lfc_mrna, TRUE)
if (nrow(mrna_sig) > 0) {
  cnts <- mrna_sig %>% count(timepoint, direction, name = "n")
  for (tp in c("6h","24h")) {
    tp_df <- mrna_sig %>% filter(timepoint == tp)
    up <- tp_df %>% filter(direction == "UP")
    dn <- tp_df %>% filter(direction == "DOWN")
    fig2_lines <- append_lines(fig2_lines,
      paste0("At ", tp, ", ", nrow(tp_df), " mRNAs passed the significance threshold (", nrow(up), " upregulated and ", nrow(dn), " downregulated)."),
      paste0("Top ", tp, " upregulated mRNAs: ", format_top_features(up %>% arrange(padj, desc(abs(log2FC))), 20, "label", "log2FC", "padj"), "."),
      paste0("Top ", tp, " downregulated mRNAs: ", format_top_features(dn %>% arrange(padj, desc(abs(log2FC))), 20, "label", "log2FC", "padj"), ".")
    )
  }
  only6 <- safe_read(find_any(c("DE_only6_annotated.csv", "DE_only6.csv")))
  shared <- safe_read(find_any(c("DE_common_annotated.csv", "DE_common.csv")))
  only24 <- safe_read(find_any(c("DE_only24_annotated.csv", "DE_only24.csv")))
  if (!is.null(only6) || !is.null(shared) || !is.null(only24)) {
    get_labels <- function(df) if (is.null(df)) character(0) else std_chr(df[[pick_label(df, id_col = first_existing_col(df, c("gene_id","GeneID","ENSEMBL"))) ]])
    fig2_lines <- append_lines(fig2_lines,
      paste0("Overlap classes derived from 6 h versus 24 h differential expression: ",
             length(unique(get_labels(only6))), " mRNAs were 6 h-only, ",
             length(unique(get_labels(shared))), " were shared between 6 h and 24 h, and ",
             length(unique(get_labels(only24))), " were 24 h-only."),
      paste0("Representative 6 h-only mRNAs: ", collapse_labels(get_labels(only6), 40), "."),
      paste0("Representative shared mRNAs: ", collapse_labels(get_labels(shared), 40), "."),
      paste0("Representative 24 h-only mRNAs: ", collapse_labels(get_labels(only24), 40), ".")
    )
  }
  write_df_txt(mrna_sig %>% arrange(timepoint, direction, padj), file.path(out_dir, "Fig2_all_significant_mRNAs.txt"), "Figure 2 - All significant mRNAs")
} else fig2_lines <- append_lines(fig2_lines, "No mRNA DE summary could be generated from the expected files.")
safe_write_lines(fig2_lines, file.path(out_dir, "Fig2_results_report.txt"))
master_lines <- append_lines(master_lines, fig2_lines, "")

fig3_lines <- c("Figure 3. mRNA response classes linking commitment and maintenance", strrep("=", 67), "")
mrna_wide <- safe_read(find_any(c("SUPP_Fig3_mRNA_programs_wide.csv", "Fig3_mRNA_programs_wide.csv"), dirs = c(fig3_dir, base_dir, bulk_dir)))
heat_mrna <- safe_read(find_any(c("SUPP_Fig3C_top_mRNAs_for_heatmap.csv", "Fig3_top_mRNAs_for_heatmap.csv"), dirs = c(fig3_dir, base_dir)))
if (!is.null(mrna_wide) && nrow(mrna_wide) > 0) {
  counts_phase <- mrna_wide %>% count(phase, class_direction, name = "n") %>% arrange(phase, class_direction)
  fig3_lines <- append_lines(fig3_lines,
    "Response classes were defined operationally as Early (6 h-only), Sustained (shared between 6 h and 24 h), and Late (24 h-only).",
    paste0("Counts by class and direction: ", paste(apply(counts_phase, 1, function(r) paste0(r[["phase"]], " ", r[["class_direction"]], "=", r[["n"]])), collapse = "; "), ".")
  )
  for (ph in c("Early","Sustained","Late")) {
    sub <- mrna_wide %>% filter(phase == ph) %>% arrange(padj_min, desc(max_absFC))
    fig3_lines <- append_lines(fig3_lines, paste0("Representative ", ph, " mRNAs: ", format_top_features(sub, 25, "gene_symbol", "max_absFC", "padj_min"), "."))
  }
  sust <- mrna_wide %>% filter(phase == "Sustained", !is.na(log2FC_6h), !is.na(log2FC_24h)) %>% mutate(mean_absFC = rowMeans(cbind(abs(log2FC_6h), abs(log2FC_24h)), na.rm = TRUE)) %>% arrange(padj_min, desc(mean_absFC))
  fig3_lines <- append_lines(fig3_lines, paste0("Strongest sustained mRNAs, ranked by combined effect across 6 h and 24 h: ", format_top_features(sust, 30, "gene_symbol", "mean_absFC", "padj_min"), "."))
  if (!is.null(heat_mrna)) write_df_txt(heat_mrna, file.path(out_dir, "Fig3_top_mRNAs_heatmap.txt"), "Figure 3 - mRNAs used in the heatmap")
  write_df_txt(mrna_wide %>% arrange(phase, class_direction, padj_min), file.path(out_dir, "Fig3_all_mRNA_response_classes.txt"), "Figure 3 - All mRNA response classes")
} else fig3_lines <- append_lines(fig3_lines, "Figure 3 outputs were not detected.")
safe_write_lines(fig3_lines, file.path(out_dir, "Fig3_results_report.txt"))
master_lines <- append_lines(master_lines, fig3_lines, "")

fig4_lines <- c("Figure 4. lncRNA response classes and regulatory rewiring", strrep("=", 59), "")
lnc6 <- standardize_de(safe_read(find_file("NE_vs_Ctrl_6h_lncRNA_DE.txt")), "lncRNA", "6h")
lnc24 <- standardize_de(safe_read(find_file("NE_vs_Ctrl_24h_lncRNA_DE.txt")), "lncRNA", "24h")
lnc_sig <- sig_subset(bind_rows(lnc6, lnc24), alpha_lnc, lfc_lnc, TRUE)
lnc_counts_phase_source <- safe_read(find_any(c("SUPP_Fig4B1_lnc_counts_phase_source.csv"), dirs = c(fig45_dir, base_dir)))
lnc_counts_phase_time <- safe_read(find_any(c("SUPP_Fig4B2_lncDE_counts_phase_time_direction.csv"), dirs = c(fig45_dir, base_dir)))
lnc_heatmap <- safe_read(find_any(c("SUPP_Fig4D_top_lnc_for_heatmap.csv", "SUPP_Fig4B3_top_lnc_for_heatmap.csv"), dirs = c(fig4r_dir, fig45_dir, base_dir)))
lnc_hubs <- safe_read(find_any(c("SUPP_Fig4E_top_hubs_by_time.csv", "lncRNA_hubs_by_time.csv"), dirs = c(fig4r_dir, base_dir)))
lnc_rewire <- safe_read(find_any(c("SUPP_Fig4_rewiring_edge_counts_curated.csv", "correlation_summary_counts.csv"), dirs = c(fig4r_dir, base_dir)))
if (nrow(lnc_sig) > 0) {
  for (tp in c("6h","24h")) {
    sub <- lnc_sig %>% filter(timepoint == tp)
    up <- sub %>% filter(direction == "UP")
    dn <- sub %>% filter(direction == "DOWN")
    fig4_lines <- append_lines(fig4_lines,
      paste0("At ", tp, ", ", nrow(sub), " lncRNAs were significant (", nrow(up), " upregulated and ", nrow(dn), " downregulated)."),
      paste0("Top ", tp, " upregulated lncRNAs (gene_id labels): ", format_top_features(up %>% arrange(padj, desc(abs(log2FC))), 25, "feature_id", "log2FC", "padj"), "."),
      paste0("Top ", tp, " downregulated lncRNAs (gene_id labels): ", format_top_features(dn %>% arrange(padj, desc(abs(log2FC))), 25, "feature_id", "log2FC", "padj"), ".")
    )
  }
  if (!is.null(lnc_counts_phase_source)) {
    fig4_lines <- append_lines(fig4_lines, paste0("lncRNA response-class composition by phase and source: ", paste(lnc_counts_phase_source %>% mutate(txt = paste0(phase, "-", source, "=", n_genes)) %>% pull(txt), collapse = "; "), "."))
  }
  if (!is.null(lnc_counts_phase_time)) {
    fig4_lines <- append_lines(fig4_lines, paste0("Counts across phase, time, and direction: ", paste(lnc_counts_phase_time %>% mutate(txt = paste0(phase, " ", timepoint, " ", direction, "=", n_genes)) %>% pull(txt), collapse = "; "), "."))
  }
  if (!is.null(lnc_heatmap) && nrow(lnc_heatmap) > 0) {
    for (ph in unique(lnc_heatmap$phase)) {
      fig4_lines <- append_lines(fig4_lines, paste0("Representative ", ph, " lncRNAs in the refined heatmap: ", collapse_labels(lnc_heatmap %>% filter(phase == ph) %>% pull(gene_id), 30), "."))
    }
  }
  if (!is.null(lnc_hubs) && nrow(lnc_hubs) > 0) {
    time_col <- first_existing_col(lnc_hubs, c("time","timepoint"))
    hub_label <- first_existing_col(lnc_hubs, c("lncRNA","gene_id","ncRNA_label"))
    edge_col <- first_existing_col(lnc_hubs, c("n_edges","degree","n_targets"))
    lnc_hubs2 <- lnc_hubs %>% mutate(.time = .data[[time_col]], .label = .data[[hub_label]], .score = .data[[edge_col]])
    for (tp in unique(lnc_hubs2$.time)) {
      sub <- lnc_hubs2 %>% filter(.time == tp) %>% arrange(desc(.score))
      fig4_lines <- append_lines(fig4_lines, paste0("Top ", tp, " lncRNA hubs: ", format_top_features(sub, 20, ".label", ".score", NULL), "."))
    }
    write_df_txt(lnc_hubs, file.path(out_dir, "Fig4_lnc_hubs_by_time.txt"), "Figure 4 - lncRNA hubs")
  }
  if (!is.null(lnc_rewire) && nrow(lnc_rewire) > 0) {
    if (all(c("category","n_edges") %in% names(lnc_rewire))) {
      fig4_lines <- append_lines(fig4_lines, paste0("Curated rewiring summary: ", paste(paste0(lnc_rewire$category, "=", lnc_rewire$n_edges), collapse = "; "), "."))
    }
  }
  write_df_txt(lnc_sig %>% arrange(timepoint, direction, padj), file.path(out_dir, "Fig4_all_significant_lncRNAs.txt"), "Figure 4 - All significant lncRNAs")
} else fig4_lines <- append_lines(fig4_lines, "No lncRNA summary could be generated from the expected files.")
safe_write_lines(fig4_lines, file.path(out_dir, "Fig4_results_report.txt"))
master_lines <- append_lines(master_lines, fig4_lines, "")

fig5_lines <- c("Figure 5. KEGG pathway programs across commitment and maintenance", strrep("=", 67), "")
kegg_stats <- safe_read(find_any(c("SUPP_Fig5_full_KEGG_stats_wide.csv"), dirs = c(fig45_dir, base_dir)))
kegg_sel <- safe_read(find_any(c("SUPP_Fig5_selected_pathways_dotplot_data.csv"), dirs = c(fig45_dir, base_dir)))
if (!is.null(kegg_stats) && nrow(kegg_stats) > 0) {
  if ("transition_shared" %in% names(kegg_stats)) {
    tc <- kegg_stats %>% filter(KEGG_group == "Shared") %>% count(transition_shared, name = "n")
    fig5_lines <- append_lines(fig5_lines, paste0("Shared KEGG transition classes: ", paste(paste0(tc$transition_shared, "=", tc$n), collapse = "; "), "."))
  }
  top_kegg <- function(df, group_name, sign = c("activated","suppressed"), n = 15) {
    sign <- match.arg(sign)
    sub <- df %>% filter(KEGG_group == group_name)
    if (group_name == "Only6h") {
      sub <- if (sign == "activated") sub %>% filter(!is.na(NES_6h), NES_6h > 0) else sub %>% filter(!is.na(NES_6h), NES_6h < 0)
      sub <- sub %>% arrange(padj_6h, desc(abs(NES_6h)))
    } else if (group_name == "Only24h") {
      sub <- if (sign == "activated") sub %>% filter(!is.na(NES_24h), NES_24h > 0) else sub %>% filter(!is.na(NES_24h), NES_24h < 0)
      sub <- sub %>% arrange(padj_24h, desc(abs(NES_24h)))
    } else {
      tr <- if (sign == "activated") "Maintained_Activated" else "Maintained_Suppressed"
      sub <- sub %>% filter(transition_shared == tr) %>% arrange(padj_min, desc(abs(NES_6h) + abs(NES_24h)))
    }
    sub %>% slice_head(n = n)
  }
  fig5_lines <- append_lines(fig5_lines,
    paste0("Top Only6h activated pathways: ", format_top_features(top_kegg(kegg_stats, "Only6h", "activated"), 20, "Description", "NES_6h", "padj_6h"), "."),
    paste0("Top Only6h suppressed pathways: ", format_top_features(top_kegg(kegg_stats, "Only6h", "suppressed"), 20, "Description", "NES_6h", "padj_6h"), "."),
    paste0("Top Shared maintained activated pathways: ", format_top_features(top_kegg(kegg_stats, "Shared", "activated"), 20, "Description", "NES_6h", "padj_min"), "."),
    paste0("Top Shared maintained suppressed pathways: ", format_top_features(top_kegg(kegg_stats, "Shared", "suppressed"), 20, "Description", "NES_6h", "padj_min"), "."),
    paste0("Top Only24h activated pathways: ", format_top_features(top_kegg(kegg_stats, "Only24h", "activated"), 20, "Description", "NES_24h", "padj_24h"), "."),
    paste0("Top Only24h suppressed pathways: ", format_top_features(top_kegg(kegg_stats, "Only24h", "suppressed"), 20, "Description", "NES_24h", "padj_24h"), ".")
  )
  write_df_txt(kegg_stats, file.path(out_dir, "Fig5_full_kegg_stats.txt"), "Figure 5 - Full KEGG statistics")
  if (!is.null(kegg_sel)) write_df_txt(kegg_sel, file.path(out_dir, "Fig5_selected_pathways_dotplot.txt"), "Figure 5 - Selected pathways used in the dotplot")
} else fig5_lines <- append_lines(fig5_lines, "Figure 5 KEGG outputs were not detected.")
safe_write_lines(fig5_lines, file.path(out_dir, "Fig5_results_report.txt"))
master_lines <- append_lines(master_lines, fig5_lines, "")

fig6_lines <- c("Figure 6. TF governors of commitment and maintenance", strrep("=", 54), "")
tf_rank <- safe_read(find_any(c("SUPP_Fig6_TF_rankings.csv", "SUPP_Fig6_TF_ranking.csv", "SUPP_TF_ranking_commitment_vs_maintenance.csv"), dirs = c(fig6_dir, base_dir)))
tf_kegg <- safe_read(find_any(c("SUPP_TF_by_KEGG_counts_context_DE_phase_sign.csv", "SUPP_Fig6_TF_KEGG_counts.csv"), dirs = c(fig6_dir, base_dir)))
if (!is.null(tf_rank) && nrow(tf_rank) > 0) {
  context_col <- first_existing_col(tf_rank, c("context","class","group"))
  tf_col <- first_existing_col(tf_rank, c("TF_key","TF","tf"))
  target_col <- first_existing_col(tf_rank, c("total_targets","n_targets","targets"))
  path_col <- first_existing_col(tf_rank, c("n_pathways_hit","n_pathways"))
  tf_rank2 <- tf_rank %>% mutate(context = .data[[context_col]], TF_key = .data[[tf_col]], total_targets = suppressWarnings(as.numeric(.data[[target_col]])), n_pathways_hit = suppressWarnings(as.numeric(if (!is.na(path_col)) .data[[path_col]] else NA_real_)))
  for (ctx in unique(tf_rank2$context)) {
    sub <- tf_rank2 %>% filter(context == ctx) %>% arrange(desc(total_targets), desc(n_pathways_hit))
    fig6_lines <- append_lines(fig6_lines, paste0("Top TFs in ", ctx, ": ", format_top_features(sub, 20, "TF_key", "total_targets", NULL), "."))
  }
  write_df_txt(tf_rank2, file.path(out_dir, "Fig6_TF_ranking.txt"), "Figure 6 - TF ranking")
}
if (!is.null(tf_kegg) && nrow(tf_kegg) > 0) {
  tf_col <- first_existing_col(tf_kegg, c("TF_key","TF"))
  desc_col <- first_existing_col(tf_kegg, c("Description","pathway"))
  ctx_col <- first_existing_col(tf_kegg, c("context","class"))
  n_col <- first_existing_col(tf_kegg, c("n_targets","count"))
  sign_col <- first_existing_col(tf_kegg, c("sign_class","sign"))
  tf_kegg2 <- tf_kegg %>% mutate(TF_key = .data[[tf_col]], Description = .data[[desc_col]], context = .data[[ctx_col]], n_targets = suppressWarnings(as.numeric(.data[[n_col]])), sign_class = if (!is.na(sign_col)) .data[[sign_col]] else NA_character_)
  for (ctx in unique(tf_kegg2$context)) {
    sub <- tf_kegg2 %>% filter(context == ctx) %>% arrange(desc(n_targets))
    txt <- sub %>% mutate(label = paste0(TF_key, " -> ", Description, " (", sign_class, ", n_targets=", n_targets, ")")) %>% pull(label)
    fig6_lines <- append_lines(fig6_lines, paste0("Highest TF-pathway associations in ", ctx, ": ", paste(head(txt, 25), collapse = "; "), "."))
  }
  write_df_txt(tf_kegg2, file.path(out_dir, "Fig6_TF_pathway_counts.txt"), "Figure 6 - TF by pathway counts")
}
if (all(is.null(tf_rank), is.null(tf_kegg))) fig6_lines <- append_lines(fig6_lines, "Figure 6 outputs were not detected.")
safe_write_lines(fig6_lines, file.path(out_dir, "Fig6_results_report.txt"))
master_lines <- append_lines(master_lines, fig6_lines, "")

fig7_lines <- c("Figure 7. miRNA-mRNA regulatory layer", strrep("=", 39), "")
mirna_rank_files <- list_files_matching("^SUPP_miRNA_rank_.*\\.csv$", dirs = c(fig7_dir, base_dir))
mirna_ranks <- if (length(mirna_rank_files) > 0) bind_rows(lapply(mirna_rank_files, function(p) {
  df <- safe_read(p)
  if (is.null(df)) return(NULL)
  df$combo <- str_remove(basename(p), "^SUPP_miRNA_rank_") %>% str_remove("\\.csv$")
  df
})) else NULL
edges_all <- safe_read(find_any(c("SUPP_edges_ALL_COMBOS_KEGG_phase.csv"), dirs = c(fig7_dir, base_dir)))
mirna_de_up6 <- safe_read(find_any(c("DE_miRNA_Up_NE6.txt")))
mirna_de_dn6 <- safe_read(find_any(c("DE_miRNA_Down_NE6.txt")))
mirna_de_up24 <- safe_read(find_any(c("DE_miRNA_Up_NE24.txt")))
mirna_de_dn24 <- safe_read(find_any(c("DE_miRNA_Down_NE24.txt")))
if (!is.null(mirna_ranks) && nrow(mirna_ranks) > 0) {
  for (cb in unique(mirna_ranks$combo)) {
    sub <- mirna_ranks %>% filter(combo == cb) %>% arrange(desc(total_targets), desc(n_pathways_hit))
    fig7_lines <- append_lines(fig7_lines, paste0("Top miRNAs in ", cb, ": ", format_top_features(sub, 15, "miRNA", "total_targets", NULL), "."))
  }
  write_df_txt(mirna_ranks, file.path(out_dir, "Fig7_miRNA_ranking_by_combo.txt"), "Figure 7 - miRNA ranking by combo")
}
if (!is.null(edges_all) && nrow(edges_all) > 0) {
  path_counts <- edges_all %>% count(combo, KEGG_group, Description, sign_class, sort = TRUE, name = "n_edges")
  gene_label_col <- first_existing_col(edges_all, c("gene_symbol", "gene_name", "label", "mRNA", "gene_key"))
  if (!is.na(gene_label_col)) {
    gene_counts <- edges_all %>% count(combo, .data[[gene_label_col]], sort = TRUE, name = "n_edges")
    names(gene_counts)[2] <- "gene_label"
    for (cb in unique(gene_counts$combo)) {
      fig7_lines <- append_lines(fig7_lines, paste0("Representative mRNA targets in ", cb, ": ", collapse_labels(gene_counts %>% filter(combo == cb) %>% pull(gene_label), 30), "."))
    }
  }
  for (cb in unique(path_counts$combo)) {
    txt <- path_counts %>% filter(combo == cb) %>% mutate(label = paste0(Description, " [", KEGG_group, "; ", sign_class, "; edges=", n_edges, "]")) %>% pull(label)
    fig7_lines <- append_lines(fig7_lines, paste0("Pathways represented in ", cb, ": ", paste(head(txt, 20), collapse = "; "), "."))
  }
  write_df_txt(edges_all, file.path(out_dir, "Fig7_all_miRNA_mRNA_edges.txt"), "Figure 7 - miRNA-mRNA edges")
}
mi_counts <- tibble(
  timepoint = c("6h","6h","24h","24h"),
  direction = c("UP","DOWN","UP","DOWN"),
  n = c(nrow(mirna_de_up6 %||% tibble()), nrow(mirna_de_dn6 %||% tibble()), nrow(mirna_de_up24 %||% tibble()), nrow(mirna_de_dn24 %||% tibble()))
)
if (sum(mi_counts$n) > 0) {
  fig7_lines <- append_lines(fig7_lines, paste0("miRNA differential-expression inventory: ", paste(paste0(mi_counts$timepoint, " ", mi_counts$direction, "=", mi_counts$n), collapse = "; "), "."))
}
if (is.null(mirna_ranks) && is.null(edges_all)) fig7_lines <- append_lines(fig7_lines, "Figure 7 outputs were not detected.")
safe_write_lines(fig7_lines, file.path(out_dir, "Fig7_results_report.txt"))
master_lines <- append_lines(master_lines, fig7_lines, "")

fig8_lines <- c("Figure 8. Integrative ceRNA candidate network", strrep("=", 48), "")
coh_lnc_pairs <- safe_read(find_any(c("COHERENT_lncRNA_miRNA_pairs.csv", "COHERENT_ncRNA_miRNA_pairs.csv"), dirs = c(fig8_dir, base_dir)))
coh_circ_pairs <- safe_read(find_any(c("COHERENT_circRNA_miRNA_pairs.csv"), dirs = c(fig8_dir, base_dir)))
coh_lnc_trip <- safe_read(find_any(c("COHERENT_lncRNA_miRNA_mRNA_triplets.csv", "COHERENT_ceRNA_triplets_all.csv"), dirs = c(fig8_dir, base_dir)))
coh_circ_trip <- safe_read(find_any(c("COHERENT_circRNA_miRNA_mRNA_triplets.csv"), dirs = c(fig8_dir, base_dir)))
priority_tri <- safe_read(find_any(c("PRIORITY_ceRNA_triplets_best_per_ncRNA_miRNA_mRNA.csv"), dirs = c(fig8_dir, base_dir)))
pair_sum <- safe_read(find_any(c("SUMMARY_coherent_ncRNA_miRNA_pairs.csv"), dirs = c(fig8_dir, base_dir)))
trip_sum <- safe_read(find_any(c("SUMMARY_coherent_triplets_by_phase.csv"), dirs = c(fig8_dir, base_dir)))
all_pairs <- bind_rows(coh_lnc_pairs, coh_circ_pairs)
all_trip <- bind_rows(coh_lnc_trip, coh_circ_trip)
if (!is.null(all_pairs) && nrow(all_pairs) > 0) {
  ncrna_col <- first_existing_col(all_pairs, c("ncrna_type","ncRNA_type"))
  time_col <- first_existing_col(all_pairs, c("timepoint","pair_time_state"))
  tier_col <- first_existing_col(all_pairs, c("interaction_tier","tier"))
  pair_counts <- all_pairs %>% count(.data[[ncrna_col]], .data[[time_col]], .data[[tier_col]], name = "n_pairs")
  names(pair_counts)[1:3] <- c("ncrna_type","timepoint","interaction_tier")
  fig8_lines <- append_lines(fig8_lines, paste0("Coherent ncRNA-miRNA pairs by layer, time, and evidence tier: ", paste(apply(pair_counts, 1, function(r) paste0(r[["ncrna_type"]], " ", r[["timepoint"]], " ", r[["interaction_tier"]], "=", r[["n_pairs"]])), collapse = "; "), "."))
}
if (!is.null(trip_sum) && nrow(trip_sum) > 0) {
  nm <- intersect(c("ncrna_type","timepoint","phase_nc","phase_mRNA","n_triplets"), names(trip_sum))
  if (length(nm) >= 5) {
    fig8_lines <- append_lines(fig8_lines, paste0("Triplet counts by phase relation: ", paste(apply(trip_sum[, nm], 1, function(r) paste0(r[[1]], " ", r[[2]], " ", r[[3]], "->", r[[4]], "=", r[[5]])), collapse = "; "), "."))
  }
}
if (!is.null(priority_tri) && nrow(priority_tri) > 0) {
  ncrna_col <- first_existing_col(priority_tri, c("ncRNA_label","ncrna_id_raw","gene_id"))
  gene_col <- first_existing_col(priority_tri, c("gene_symbol","gene_name","gene_key"))
  mi_col <- first_existing_col(priority_tri, c("miRNA","mirna"))
  phase_nc_col <- first_existing_col(priority_tri, c("phase","phase_nc"))
  phase_m_col <- first_existing_col(priority_tri, c("phase_mRNA"))
  tp_col <- first_existing_col(priority_tri, c("timepoint"))
  e_col <- first_existing_col(priority_tri, c("E"))
  priority_tri2 <- priority_tri %>% mutate(triplet = paste0(.data[[ncrna_col]], " -- ", .data[[mi_col]], " -- ", .data[[gene_col]], " [", .data[[tp_col]], ", ", .data[[phase_nc_col]], "->", .data[[phase_m_col]], ", E=", fmt_num(.data[[e_col]], 3), "]"))
  fig8_lines <- append_lines(fig8_lines, paste0("Representative priority ceRNA triplets: ", paste(head(priority_tri2$triplet, 50), collapse = "; "), "."))
  rep_mir <- priority_tri %>% count(.data[[mi_col]], sort = TRUE, name = "n_priority_triplets")
  names(rep_mir)[1] <- "miRNA"
  fig8_lines <- append_lines(fig8_lines, paste0("Most recurrent miRNAs across priority triplets: ", format_top_features(rep_mir, 25, "miRNA", "n_priority_triplets", NULL), "."))
  rep_genes <- priority_tri %>% count(.data[[gene_col]], sort = TRUE, name = "n_priority_triplets")
  names(rep_genes)[1] <- "gene_symbol"
  fig8_lines <- append_lines(fig8_lines, paste0("Most recurrent mRNA targets across priority triplets: ", format_top_features(rep_genes, 25, "gene_symbol", "n_priority_triplets", NULL), "."))
  rep_nc <- priority_tri %>% count(.data[[ncrna_col]], sort = TRUE, name = "n_priority_triplets")
  names(rep_nc)[1] <- "ncRNA_label"
  fig8_lines <- append_lines(fig8_lines, paste0("Most recurrent ncRNA candidates across priority triplets: ", format_top_features(rep_nc, 25, "ncRNA_label", "n_priority_triplets", NULL), "."))
  write_df_txt(priority_tri, file.path(out_dir, "Fig8_priority_ceRNA_triplets.txt"), "Figure 8 - Priority ceRNA triplets")
}
if (!is.null(all_pairs) && nrow(all_pairs) > 0) write_df_txt(all_pairs, file.path(out_dir, "Fig8_all_coherent_ncRNA_miRNA_pairs.txt"), "Figure 8 - All coherent ncRNA-miRNA pairs")
if (!is.null(all_trip) && nrow(all_trip) > 0) write_df_txt(all_trip, file.path(out_dir, "Fig8_all_coherent_ceRNA_triplets.txt"), "Figure 8 - All coherent triplets")
if (!is.null(pair_sum)) write_df_txt(pair_sum, file.path(out_dir, "Fig8_pair_summary.txt"), "Figure 8 - Pair summary")
if (!is.null(trip_sum)) write_df_txt(trip_sum, file.path(out_dir, "Fig8_triplet_summary.txt"), "Figure 8 - Triplet summary")
if (all(is.null(all_pairs), is.null(all_trip), is.null(priority_tri))) fig8_lines <- append_lines(fig8_lines, "Figure 8 outputs were not detected yet.")
safe_write_lines(fig8_lines, file.path(out_dir, "Fig8_results_report.txt"))
master_lines <- append_lines(master_lines, fig8_lines, "")

supp_lines <- c("Supplementary / inventory sections", strrep("=", 33), "")
# circRNA supplementary report
circ_lines <- c("Supplementary circRNA layer", strrep("-", 26), "")
circ6 <- standardize_de(safe_read(find_any(c("06_circRNA_sequences_de6_DE_complete_table.txt"))), "circRNA", "6h")
circ24 <- standardize_de(safe_read(find_any(c("07_circRNA_sequences_de24_DE_complete_table.txt"))), "circRNA", "24h")
use_padj_circ <- any(!is.na(c(circ6$padj, circ24$padj)))
circ_sig <- sig_subset(bind_rows(circ6, circ24), alpha_circ, lfc_circ, use_padj = use_padj_circ)
if (nrow(circ_sig) > 0) {
  for (tp in c("6h","24h")) {
    sub <- circ_sig %>% filter(timepoint == tp)
    if (nrow(sub) == 0) next
    up <- sub %>% filter(direction == "UP")
    dn <- sub %>% filter(direction == "DOWN")
    circ_lines <- append_lines(circ_lines,
      paste0("At ", tp, ", ", nrow(sub), " circRNAs passed the configured threshold (", nrow(up), " upregulated and ", nrow(dn), " downregulated)."),
      paste0("Top ", tp, " circRNAs: ", format_top_features(sub %>% arrange(if (use_padj_circ) padj else pvalue, desc(abs(log2FC))), 25, "feature_id", "log2FC", if (use_padj_circ) "padj" else "pvalue"), ".")
    )
  }
  write_df_txt(circ_sig, file.path(out_dir, "Supp_circRNA_all_significant.txt"), "Supplementary circRNA - significant circRNAs")
}
phase_files_circ <- c(early = find_any(c("early_circRNA.csv"), dirs = c(circ_dir, base_dir)), sustained = find_any(c("sustained_circRNA.csv"), dirs = c(circ_dir, base_dir)), late = find_any(c("late_circRNA.csv"), dirs = c(circ_dir, base_dir)), reversal = find_any(c("reversal_circRNA.csv"), dirs = c(circ_dir, base_dir)))
phase_tables_circ <- imap_dfr(phase_files_circ, ~read_phase_table(.x, entity = "circRNA", phase_name = str_to_title(.y)))
if (nrow(phase_tables_circ) > 0) {
  counts <- phase_tables_circ %>% count(phase, class_direction, name = "n")
  circ_lines <- append_lines(circ_lines, paste0("circRNA response classes: ", paste(apply(counts, 1, function(r) paste0(r[["phase"]], " ", r[["class_direction"]], "=", r[["n"]])), collapse = "; "), "."))
  for (ph in unique(phase_tables_circ$phase)) {
    circ_lines <- append_lines(circ_lines, paste0("Representative ", ph, " circRNAs: ", collapse_labels(phase_tables_circ %>% filter(phase == ph) %>% pull(feature_id), 25), "."))
  }
  write_df_txt(phase_tables_circ, file.path(out_dir, "Supp_circRNA_response_classes.txt"), "Supplementary circRNA - response classes")
}
safe_write_lines(circ_lines, file.path(out_dir, "Supplementary_circRNA_results_report.txt"))
supp_lines <- append_lines(supp_lines, circ_lines, "")

# recurring candidates section
rec_lines <- c("Cross-layer recurring candidates", strrep("-", 30), "")
rec_mrna <- bind_rows(
  mrna_sig %>% transmute(label = label, layer = "mRNA_DE"),
  if (!is.null(priority_tri) && nrow(priority_tri) > 0) {
    gene_col <- first_existing_col(priority_tri, c("gene_symbol","gene_name","gene_key"))
    priority_tri %>% transmute(label = std_chr(.data[[gene_col]]), layer = "ceRNA_triplets")
  } else tibble(label = character(), layer = character())
)
if (nrow(rec_mrna) > 0) {
  rec_tbl <- rec_mrna %>% filter(label != "") %>% count(label, sort = TRUE, name = "n_layers")
  rec_lines <- append_lines(rec_lines, paste0("Recurrent mRNA candidates across DE and integrative layers: ", format_top_features(rec_tbl, 30, "label", "n_layers", NULL), "."))
  write_df_txt(rec_tbl, file.path(out_dir, "Recurring_mRNA_candidates.txt"), "Cross-layer recurring mRNA candidates")
}
if (!is.null(priority_tri) && nrow(priority_tri) > 0) {
  ncrna_col <- first_existing_col(priority_tri, c("ncRNA_label","gene_id","ncrna_id_raw"))
  mi_col <- first_existing_col(priority_tri, c("miRNA"))
  rec_nc <- priority_tri %>% count(.data[[ncrna_col]], sort = TRUE, name = "n_triplets")
  names(rec_nc)[1] <- "ncRNA_label"
  rec_mi <- priority_tri %>% count(.data[[mi_col]], sort = TRUE, name = "n_triplets")
  names(rec_mi)[1] <- "miRNA"
  rec_lines <- append_lines(rec_lines,
    paste0("Recurrent ncRNA candidates across priority triplets: ", format_top_features(rec_nc, 25, "ncRNA_label", "n_triplets", NULL), "."),
    paste0("Recurrent miRNA candidates across priority triplets: ", format_top_features(rec_mi, 25, "miRNA", "n_triplets", NULL), ".")
  )
  write_df_txt(rec_nc, file.path(out_dir, "Recurring_ncRNA_candidates.txt"), "Cross-layer recurring ncRNA candidates")
  write_df_txt(rec_mi, file.path(out_dir, "Recurring_miRNA_candidates.txt"), "Cross-layer recurring miRNA candidates")
}
safe_write_lines(rec_lines, file.path(out_dir, "Cross_layer_recurring_candidates_report.txt"))
supp_lines <- append_lines(supp_lines, rec_lines, "")

existing_dirs <- c(bulk_dir, fig3_dir, fig45_dir, fig4r_dir, fig6_dir, fig7_dir, fig8_dir, circ_dir)
existing_dirs <- existing_dirs[dir.exists(existing_dirs)]
all_out_files <- if (length(existing_dirs) > 0) list.files(existing_dirs, recursive = TRUE, full.names = TRUE) else character(0)
if (length(all_out_files) > 0) {
  inv_tbl <- tibble(file = basename(all_out_files), path = all_out_files) %>% arrange(file)
  supp_lines <- append_lines(supp_lines, paste0(nrow(inv_tbl), " workflow output files were detected across the figure directories."))
  write_df_txt(inv_tbl, file.path(out_dir, "Workflow_output_inventory.txt"), "Workflow output inventory")
} else {
  supp_lines <- append_lines(supp_lines, "No workflow output files were detected in the expected directories.")
}
safe_write_lines(supp_lines, file.path(out_dir, "Supplementary_and_inventory_report.txt"))
master_lines <- append_lines(master_lines, supp_lines, "")

safe_write_lines(master_lines, file.path(out_dir, "MASTER_results_report.txt"))
safe_write_lines(c("# Master results report", "", master_lines), file.path(out_dir, "MASTER_results_report.md"))
message("Done. Results reports written to: ", out_dir)
