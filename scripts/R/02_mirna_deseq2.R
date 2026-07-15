#!/usr/bin/env Rscript
# Stage 02: miRNA DESeq2 analysis (+ optional ComBat-seq batch correction)
# Differential expression of miRNA counts from miRge3.0 output using the
# multifactorial design and interaction contrasts, exports miRWalk input
# lists (UP/DOWN per timepoint), volcano plots, and temporal-pattern summary.
#
# Ported from miRNAs_6_24.R with Windows paths removed, sample renaming
# generalized to a sample sheet, batch correction made optional (driven by a
# `batch` column), and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

# ============================================================================
# Command-line Arguments
# ============================================================================

option_list <- list(
  make_option(c("-c", "--counts-dir"), type = "character", default = NULL,
              help = "Directory containing miRge3.0 count matrices (miR.Counts_*.csv)"),
  make_option(c("-m", "--sample-sheet"), type = "character", default = NULL,
              help = "Sample sheet TSV (sample_id, treatment, time[, batch])"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/02_mirna_deseq2",
              help = "Output directory [default: %default]"),
  make_option(c("--alpha"), type = "numeric", default = 0.05,
              help = "Adjusted p-value threshold [default: %default]"),
  make_option(c("--lfc-cut"), type = "numeric", default = 0.5,
              help = "abs(log2FC) cutoff for UP/DOWN direction [default: %default]"),
  make_option(c("--combat"), action = "store_true", default = FALSE,
              help = "Apply ComBat-seq batch correction (requires a `batch` column) [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "PNG DPI for QC/volcano figures [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 123,
              help = "Random seed [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$`counts-dir`) || is.null(opt$`sample-sheet`)) {
  print_help(parser)
  stop("\n[ERROR] Required arguments: --counts-dir, --sample-sheet")
}

COUNTS_DIR   <- opt$`counts-dir`
SAMPLE_SHEET <- opt$`sample-sheet`
OUTPUT_DIR   <- opt$`output-dir`
ALPHA        <- opt$alpha
LFC_CUT      <- opt$`lfc-cut`
USE_COMBAT   <- opt$combat
DPI_PNG      <- opt$dpi
SEED         <- opt$seed

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
set.seed(SEED)

cat("[INFO] miRNA counts dir:", COUNTS_DIR, "\n")
cat("[INFO] Sample sheet:", SAMPLE_SHEET, "\n")
cat("[INFO] Output dir:", OUTPUT_DIR, "\n\n")

# ============================================================================
# 1. Load miRNA counts from miRge3.0 output
# ============================================================================

cat("[INFO] Loading miRNA count matrix...\n")
count_files <- list.files(COUNTS_DIR, pattern = "^miR.Counts", full.names = TRUE)
if (length(count_files) == 0) stop("[ERROR] No miR.Counts_* files found in: ", COUNTS_DIR)

mirna_counts <- NULL
for (f in count_files) {
  df <- as.data.frame(read_csv(f, show_col_types = FALSE))
  rownames(df) <- df[[1]]; df <- df[, -1, drop = FALSE]
  if (is.null(mirna_counts)) {
    mirna_counts <- df
  } else {
    mirna_counts <- merge(mirna_counts, df, by = "row.names", all = TRUE)
    rownames(mirna_counts) <- mirna_counts$Row.names; mirna_counts$Row.names <- NULL
  }
}
mirna_counts[is.na(mirna_counts)] <- 0
mirna_counts <- as.matrix(mirna_counts)
cat("[INFO] Loaded", nrow(mirna_counts), "miRNAs,", ncol(mirna_counts), "samples\n")

# ============================================================================
# 2. Sample metadata
# ============================================================================

sample_data <- read_tsv(SAMPLE_SHEET, show_col_types = FALSE)
required_cols <- c("sample_id", "treatment", "time")
missing_cols <- setdiff(required_cols, names(sample_data))
if (length(missing_cols) > 0) stop("[ERROR] Sample sheet missing columns: ", paste(missing_cols, collapse = ", "))

has_batch <- "batch" %in% names(sample_data)
samples <- sample_data$sample_id
sampleTable <- data.frame(
  sample = samples,
  treatment = factor(sample_data$treatment),
  time = factor(sample_data$time)
)
if (has_batch) sampleTable$batch <- factor(sample_data$batch)
rownames(sampleTable) <- samples
sampleTable$treatment <- relevel(sampleTable$treatment, ref = "Control")
sampleTable$time <- relevel(sampleTable$time, ref = "6h")
sampleTable$group <- paste(sampleTable$treatment, sampleTable$time, sep = "_")

common_samples <- intersect(colnames(mirna_counts), rownames(sampleTable))
if (length(common_samples) == 0) stop("[ERROR] No overlap between count columns and sample_id")
mirna_counts <- round(mirna_counts[, common_samples, drop = FALSE])
sampleTable <- sampleTable[common_samples, ]
cat("[INFO] Matched", nrow(sampleTable), "samples", if (has_batch) "(batch column present)" else "", "\n")

# ============================================================================
# 3. Optional ComBat-seq batch correction
# ============================================================================

count_matrix <- mirna_counts
if (USE_COMBAT) {
  if (!has_batch) {
    cat("[WARN] --combat requested but no `batch` column; skipping correction\n")
  } else if (!requireNamespace("sva", quietly = TRUE)) {
    cat("[WARN] sva not installed; skipping ComBat-seq\n")
  } else {
    cat("[INFO] Applying ComBat-seq batch correction...\n")
    count_matrix <- sva::ComBat_seq(counts = mirna_counts,
                                    batch = sampleTable$batch,
                                    group = sampleTable$group)
  }
}

# ============================================================================
# 4. DESeq2 (multifactorial design)
# ============================================================================

design_formula <- if (has_batch && !USE_COMBAT) {
  ~batch + time + treatment + time:treatment
} else {
  ~time + treatment + time:treatment
}
cat("[INFO] Design:", deparse(design_formula), "\n")

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sampleTable,
                              design = design_formula)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("[INFO] After filtering:", nrow(dds), "miRNAs retained\n")
dds <- DESeq(dds)

# ============================================================================
# 5. QC
# ============================================================================

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) + ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "%")) + ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() + theme_bw() +
  ggtitle(paste0("miRNA PCA: NE vs Control", if (USE_COMBAT) " (ComBat-seq corrected)" else ""))
ggsave(file.path(OUTPUT_DIR, "mirna_pca.pdf"), p_pca, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "mirna_pca.png"), p_pca, width = 8, height = 6, dpi = DPI_PNG)

sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(file.path(OUTPUT_DIR, "mirna_sample_distance.pdf"), width = 8, height = 8)
pheatmap(sampleDistMatrix, col = colors, main = "miRNA Sample Distance")
dev.off()
cat("[INFO] QC plots saved\n")

# ============================================================================
# 6. Differential Expression via interaction contrasts
#    6h  = main treatment effect (reference time)
#    24h = main treatment effect + time24h:treatmentNE interaction
# ============================================================================

rn <- resultsNames(dds)
inter_name <- rn[str_detect(rn, "time.*treatment") | str_detect(rn, "treatment.*time")][1]

get_res <- function(which_time) {
  if (which_time == "6h") {
    results(dds, name = "treatment_NE_vs_Control", alpha = ALPHA)
  } else {
    if (!is.na(inter_name)) {
      results(dds, contrast = list(c("treatment_NE_vs_Control", inter_name)), alpha = ALPHA)
    } else {
      # Fallback: subset + refit if interaction term is unavailable
      d <- dds[, dds$time == "24h"]; d$time <- droplevels(d$time)
      design(d) <- if (has_batch && !USE_COMBAT) ~batch + treatment else ~treatment
      d <- DESeq(d)
      results(d, contrast = c("treatment", "NE", "Control"), alpha = ALPHA)
    }
  }
}

res_6h  <- get_res("6h")
res_24h <- get_res("24h")

as_tbl <- function(res) {
  as.data.frame(res) %>% rownames_to_column("miRNA_id") %>%
    mutate(direction = case_when(
      !is.na(padj) & padj < ALPHA & log2FoldChange >=  LFC_CUT ~ "UP",
      !is.na(padj) & padj < ALPHA & log2FoldChange <= -LFC_CUT ~ "DOWN",
      TRUE ~ "NS"
    ))
}
res_6h_tbl  <- as_tbl(res_6h)
res_24h_tbl <- as_tbl(res_24h)

write_csv(res_6h_tbl,  file.path(OUTPUT_DIR, "miRNA_res_NE_vs_Ctrl_6h.csv"))
write_csv(res_24h_tbl, file.path(OUTPUT_DIR, "miRNA_res_NE_vs_Ctrl_24h.csv"))

# miRWalk input lists (DE only) + UP/DOWN split files consumed by script 09
export_mir_lists <- function(tbl, time_label) {
  de <- tbl %>% filter(direction != "NS")
  write_tsv(de, file.path(OUTPUT_DIR, paste0("miRNA_DE_NE_vs_Ctrl_", time_label, "_for_miRWalk.txt")))
  write_tsv(tbl %>% filter(direction == "UP"),
            file.path(OUTPUT_DIR, paste0("miRNA_DE_", time_label, "_UP.txt")))
  write_tsv(tbl %>% filter(direction == "DOWN"),
            file.path(OUTPUT_DIR, paste0("miRNA_DE_", time_label, "_DOWN.txt")))
  cat("[INFO]", time_label, "DE miRNAs:", nrow(de),
      "(UP:", sum(tbl$direction == "UP"), "DOWN:", sum(tbl$direction == "DOWN"), ")\n")
}
export_mir_lists(res_6h_tbl,  "6h")
export_mir_lists(res_24h_tbl, "24h")

# ============================================================================
# 7. Volcano plots
# ============================================================================

make_volcano <- function(tbl, title, filename, top_labels = 15) {
  df <- tbl %>%
    mutate(negLog10 = -log10(pmax(padj, 1e-300)),
           sig = ifelse(direction != "NS", direction, "NS"))
  lab <- df %>% filter(direction != "NS") %>% arrange(padj) %>% slice_head(n = top_labels)
  p <- ggplot(df, aes(log2FoldChange, negLog10, color = sig)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(ALPHA), linetype = 3, color = "grey60") +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 3, color = "grey60") +
    ggrepel::geom_text_repel(data = lab, aes(label = miRNA_id), size = 2.7, max.overlaps = Inf) +
    scale_color_manual(values = c(UP = "#EF4444", DOWN = "#3B82F6", NS = "grey75")) +
    theme_bw() + ggtitle(title) + xlab("log2(FC) NE vs Ctrl") + ylab("-log10 adj. P")
  ggsave(file.path(OUTPUT_DIR, filename), p, width = 8, height = 6, dpi = DPI_PNG)
}
make_volcano(res_6h_tbl,  "miRNA Volcano (6h)",  "mirna_volcano_6h.pdf")
make_volcano(res_24h_tbl, "miRNA Volcano (24h)", "mirna_volcano_24h.pdf")

# ============================================================================
# 8. Temporal pattern summary (Early / Sustained / Late; switch classes)
# ============================================================================

red6 <- res_6h_tbl %>% transmute(miRNA_id, log2FC_6h = log2FoldChange, padj_6h = padj, dir_6h = direction)
red24 <- res_24h_tbl %>% transmute(miRNA_id, log2FC_24h = log2FoldChange, padj_24h = padj, dir_24h = direction)

summary_df <- full_join(red6, red24, by = "miRNA_id") %>%
  mutate(dir_6h = replace_na(dir_6h, "NS"), dir_24h = replace_na(dir_24h, "NS"),
    temporal_class = case_when(
      dir_6h != "NS" & dir_24h == "NS" ~ "Early_6h_only",
      dir_6h == "NS" & dir_24h != "NS" ~ "Late_24h_only",
      dir_6h != "NS" & dir_24h != "NS" ~ "Sustained_both",
      TRUE ~ "Non_DE"),
    temporal_pattern = case_when(
      dir_6h == "UP"   & dir_24h == "UP"   ~ "Sustained_UP",
      dir_6h == "DOWN" & dir_24h == "DOWN" ~ "Sustained_DOWN",
      dir_6h == "UP"   & dir_24h == "DOWN" ~ "Switch_UP_to_DOWN",
      dir_6h == "DOWN" & dir_24h == "UP"   ~ "Switch_DOWN_to_UP",
      dir_6h != "NS"   & dir_24h == "NS"   ~ paste0("Early_", dir_6h),
      dir_6h == "NS"   & dir_24h != "NS"   ~ paste0("Late_", dir_24h),
      TRUE ~ "Non_DE"))
write_csv(summary_df, file.path(OUTPUT_DIR, "miRNA_temporal_summary_table.csv"))

plot_df <- summary_df %>% filter(temporal_pattern != "Non_DE") %>%
  count(temporal_pattern, name = "n") %>% arrange(desc(n))
if (nrow(plot_df) > 0) {
  p_temp <- ggplot(plot_df, aes(x = reorder(temporal_pattern, -n), y = n)) +
    geom_col(fill = "grey35") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "miRNA temporal patterns", x = "Temporal pattern", y = "# miRNAs")
  ggsave(file.path(OUTPUT_DIR, "miRNA_temporal_patterns_barplot.png"), p_temp, width = 7, height = 5, dpi = DPI_PNG)
}

# Legacy DE lists (gene_id + direction)
write_csv(res_6h_tbl  %>% filter(direction != "NS") %>% select(mirna_id = miRNA_id, direction),
          file.path(OUTPUT_DIR, "DE_mirna_6h.csv"))
write_csv(res_24h_tbl %>% filter(direction != "NS") %>% select(mirna_id = miRNA_id, direction),
          file.path(OUTPUT_DIR, "DE_mirna_24h.csv"))

cat("\n[SUCCESS] miRNA DESeq2 analysis complete!\n")
cat("[OUTPUT] Results saved to:", OUTPUT_DIR, "\n")
