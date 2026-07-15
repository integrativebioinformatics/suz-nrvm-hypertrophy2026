#!/usr/bin/env Rscript
# Stage 01: Bulk RNA-seq DESeq2 analysis + biotype split + Mfuzz clustering
# Performs multifactorial DESeq2 (time + treatment + interaction),
# computes QC (PCA, sample distances), differential expression,
# annotates by biotype (mRNA / annotated lncRNA / novel lncRNA), exports
# phase x biotype gene sets and biotype-split DE tables consumed by the
# figure scripts, and soft clustering (Mfuzz) for temporal dynamics.
#
# Ported from Integrated_Bulk_v3_6_24h.R with Windows paths removed,
# internal duplication removed, biotype annotation made reproducible and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(rtracklayer)
  library(DESeq2)
  library(apeglm)
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(Mfuzz)
  library(Biobase)
})

# ============================================================================
# Command-line Arguments
# ============================================================================

option_list <- list(
  make_option(c("-g", "--gtf"), type = "character", default = NULL,
              help = "Merged/annotated GTF (StringTie+gffcompare) for tx2gene"),
  make_option(c("-s", "--salmon-dir"), type = "character", default = NULL,
              help = "Directory containing per-sample Salmon quantification directories"),
  make_option(c("-m", "--sample-sheet"), type = "character", default = NULL,
              help = "Sample sheet TSV (columns: sample_id, treatment, time, ...)"),
  make_option(c("--ann-annotated"), type = "character", default = "hisat2_6_24_meta_ann.csv",
              help = "Curated annotated-gene metadata (Ensembl biomaRt export; cols ensembl_gene_id, ensembl_transcript_id, external_gene_name, gene_biotype, transcript_biotype) [default: %default]"),
  make_option(c("--ann-novel"), type = "character", default = "hisat2_6_24_meta_novel_lncRNA.csv",
              help = "Curated novel-lncRNA metadata (StringTie/gffcompare; cols type, gene_id, transcript_id, qry_gene_id, pred_biotype) [default: %default]"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/01_deseq2_analysis",
              help = "Output directory [default: %default]"),
  make_option(c("--mfuzz-clusters"), type = "integer", default = 6,
              help = "Number of Mfuzz clusters [default: %default]"),
  make_option(c("--alpha"), type = "numeric", default = 0.05,
              help = "Adjusted p-value threshold for significance [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "PNG DPI for QC/Venn figures [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 123,
              help = "Random seed for reproducibility [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$gtf) || is.null(opt$`salmon-dir`) || is.null(opt$`sample-sheet`)) {
  print_help(parser)
  stop("\n[ERROR] Required arguments missing: --gtf, --salmon-dir, --sample-sheet")
}

GTF_FILE     <- opt$gtf
SALMON_DIR   <- opt$`salmon-dir`
SAMPLE_SHEET <- opt$`sample-sheet`
ANN_ANNOTATED <- opt$`ann-annotated`
ANN_NOVEL    <- opt$`ann-novel`
OUTPUT_DIR   <- opt$`output-dir`
ALPHA        <- opt$alpha
MFUZZ_CLUSTERS <- opt$`mfuzz-clusters`
DPI_PNG      <- opt$dpi
SEED         <- opt$seed

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[INFO] GTF file:", GTF_FILE, "\n")
cat("[INFO] Salmon dir:", SALMON_DIR, "\n")
cat("[INFO] Sample sheet:", SAMPLE_SHEET, "\n")
cat("[INFO] Output dir:", OUTPUT_DIR, "\n\n")

# ============================================================================
# 1. Generate tx2gene mapping + gene-level annotation from merged GTF
# ============================================================================

cat("[INFO] Parsing merged GTF...\n")
gtf_data <- rtracklayer::import(GTF_FILE)
gtf_df <- as.data.frame(gtf_data)

tx2gene <- gtf_df %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id) %>%
  distinct() %>%
  na.omit()

TX2GENE_FILE <- file.path(OUTPUT_DIR, "tx2gene.tsv")
write_tsv(tx2gene, TX2GENE_FILE)
cat("[INFO] tx2gene saved:", TX2GENE_FILE, "(", nrow(tx2gene), "transcripts )\n")

# ============================================================================
# 2. Build biotype annotation (mRNA / annotated lncRNA / novel lncRNA)
#    Merges the curated annotated-gene metadata (--ann-annotated,
#    hisat2_6_24_meta_ann.csv) with the curated novel-lncRNA metadata
#    (--ann-novel, hisat2_6_24_meta_novel_lncRNA.csv), exactly as in the
#    original Integrated_Bulk_v3_6_24h.R.
# ============================================================================

build_annotation <- function() {
  if (is.null(ANN_ANNOTATED) || !file.exists(ANN_ANNOTATED))
    stop("[ERROR] Annotated-gene metadata not found (--ann-annotated): ", ANN_ANNOTATED)
  if (is.null(ANN_NOVEL) || !file.exists(ANN_NOVEL))
    stop("[ERROR] Novel-lncRNA metadata not found (--ann-novel): ", ANN_NOVEL)

  cat("[INFO] Loading curated annotation files:\n      annotated:", ANN_ANNOTATED,
      "\n      novel:", ANN_NOVEL, "\n")
  ann_all <- read_csv(ANN_ANNOTATED, show_col_types = FALSE)
  ann_novel_raw <- read_csv(ANN_NOVEL, show_col_types = FALSE)

  need_ann <- c("ensembl_gene_id", "external_gene_name", "gene_biotype")
  need_nov <- c("gene_id", "pred_biotype")
  if (!all(need_ann %in% names(ann_all)))
    stop("[ERROR] --ann-annotated must contain: ", paste(need_ann, collapse = ", "))
  if (!all(need_nov %in% names(ann_novel_raw)))
    stop("[ERROR] --ann-novel must contain: ", paste(need_nov, collapse = ", "))

  ann_all_clean <- ann_all %>%
    transmute(
      gene_id      = as.character(ensembl_gene_id),
      gene_name    = as.character(external_gene_name),
      gene_biotype = as.character(gene_biotype),
      source       = "annotated"
    )
  ann_novel_clean <- ann_novel_raw %>%
    { if ("type" %in% names(.)) filter(., type == "transcript") else . } %>%
    transmute(
      gene_id      = as.character(gene_id),
      gene_name    = as.character(if ("qry_gene_id" %in% names(.)) qry_gene_id else gene_id),
      gene_biotype = as.character(pred_biotype),
      source       = "novel"
    )

  bind_rows(ann_all_clean, ann_novel_clean) %>%
    distinct(gene_id, .keep_all = TRUE)
}

ann_merged <- build_annotation()
ann_merged$gene_biotype <- tolower(ann_merged$gene_biotype)
cat("[INFO] Annotation genes:", nrow(ann_merged),
    "| annotated:", sum(ann_merged$source == "annotated"),
    "| novel:", sum(ann_merged$source == "novel"), "\n")

is_mrna_df  <- function(df) df %>% filter(gene_biotype == "protein_coding")
is_lnc_df   <- function(df) df %>% filter(grepl("lnc", gene_biotype))
is_novel_df <- function(df) df %>% filter(source == "novel")

# ============================================================================
# 3. Load sample metadata and Salmon quantifications
# ============================================================================

cat("[INFO] Loading sample metadata...\n")
sample_data <- read_tsv(SAMPLE_SHEET, show_col_types = FALSE)

required_cols <- c("sample_id", "treatment", "time")
missing_cols <- setdiff(required_cols, names(sample_data))
if (length(missing_cols) > 0) {
  stop("[ERROR] Sample sheet missing columns: ", paste(missing_cols, collapse = ", "))
}

samples <- sample_data$sample_id
sampleTable <- data.frame(
  sample = samples,
  treatment = factor(sample_data$treatment),
  time = factor(sample_data$time)
)
rownames(sampleTable) <- samples
sampleTable$treatment <- relevel(sampleTable$treatment, ref = "Control")
sampleTable$time <- relevel(sampleTable$time, ref = "6h")
sampleTable$group <- paste(sampleTable$treatment, sampleTable$time, sep = "_")

cat("[INFO] Samples loaded:", nrow(sampleTable), "\n")

quant_files <- file.path(SALMON_DIR, samples, "quant.sf")
names(quant_files) <- samples
if (!all(file.exists(quant_files))) {
  stop("[ERROR] Missing quant.sf files:\n",
       paste(quant_files[!file.exists(quant_files)], collapse = "\n"))
}

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)
cat("[INFO] Transcriptome loaded:", nrow(txi$abundance), "genes\n")

# ============================================================================
# 4. DESeq2 dataset with multifactorial design
# ============================================================================

cat("[INFO] Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromTximport(txi,
  colData = sampleTable,
  design = ~time + treatment + time:treatment
)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("[INFO] After filtering:", nrow(dds), "genes retained\n")

dds <- DESeq(dds)

# ============================================================================
# 5. QC: PCA and sample distance
# ============================================================================

cat("[INFO] Computing QC plots...\n")
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) + ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() + theme_bw() + ggtitle("PCA: NE vs Control by Time")
ggsave(file.path(OUTPUT_DIR, "pca_plot.pdf"), p_pca, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "pca_plot.png"), p_pca, width = 8, height = 6, dpi = DPI_PNG)

sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(file.path(OUTPUT_DIR, "sample_distance_heatmap.pdf"), width = 8, height = 8)
pheatmap(sampleDistMatrix, col = colors, main = "Sample Distance Heatmap")
dev.off()

# ============================================================================
# 6. Differential Expression (per timepoint: NE vs Control)
# ============================================================================

cat("[INFO] Running DE analysis...\n")

run_timepoint_de <- function(tp) {
  d <- dds[, dds$time == tp]
  d$time <- droplevels(d$time)
  design(d) <- ~treatment
  d <- DESeq(d)
  res <- lfcShrink(d, coef = "treatment_NE_vs_Control", type = "apeglm")
  as.data.frame(res) %>% rownames_to_column("gene_id")
}

res6_tbl  <- run_timepoint_de("6h")
res24_tbl <- run_timepoint_de("24h")

write_csv(res6_tbl,  file.path(OUTPUT_DIR, "DE_NE_vs_Ctrl_6h.csv"))
write_csv(res24_tbl, file.path(OUTPUT_DIR, "DE_NE_vs_Ctrl_24h.csv"))

# Annotate + direction
annotate_de <- function(df) {
  df %>%
    left_join(ann_merged, by = "gene_id") %>%
    mutate(
      gene_name = coalesce(gene_name, gene_id),
      direction = case_when(
        !is.na(padj) & padj < ALPHA & log2FoldChange > 0 ~ "UP",
        !is.na(padj) & padj < ALPHA & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
}
res6_ann  <- annotate_de(res6_tbl)
res24_ann <- annotate_de(res24_tbl)

# Biotype-split DE tables consumed by figure scripts 03/04/05/08/09
export_de_table <- function(df, biotype_fun, time_label, biotype_label) {
  out <- df %>%
    filter(!is.na(padj)) %>%
    biotype_fun() %>%
    select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, pvalue, padj,
           direction, gene_biotype, source) %>%
    arrange(padj)
  fname <- file.path(OUTPUT_DIR, paste0("NE_vs_Ctrl_", time_label, "_", biotype_label, "_DE.txt"))
  write_tsv(out, fname)
  cat("[INFO]", time_label, biotype_label, "DE table:", nrow(out), "genes ->", basename(fname), "\n")
}
export_de_table(res6_ann,  is_mrna_df,  "6h",  "mRNA")
export_de_table(res6_ann,  is_lnc_df,   "6h",  "lncRNA")
export_de_table(res24_ann, is_mrna_df,  "24h", "mRNA")
export_de_table(res24_ann, is_lnc_df,   "24h", "lncRNA")

# ============================================================================
# 7. Phase assignment (Early / Sustained / Late) x biotype gene sets
# ============================================================================

cat("[INFO] Building phase x biotype gene sets...\n")

sig6  <- res6_ann  %>% filter(!is.na(padj), padj < ALPHA)
sig24 <- res24_ann %>% filter(!is.na(padj), padj < ALPHA)

genes_6h  <- sig6$gene_id
genes_24h <- sig24$gene_id
intersect_genes <- intersect(genes_6h, genes_24h)
only_6h  <- setdiff(genes_6h, genes_24h)
only_24h <- setdiff(genes_24h, genes_6h)

# annotation lookup for a set of gene_ids; carries transcript_id (from tx2gene)
# so downstream IntaRNA candidate extraction (script 10) can subset the FASTA.
ann_lookup <- function(ids) {
  tibble(gene_id = ids) %>%
    left_join(ann_merged, by = "gene_id") %>%
    left_join(tx2gene, by = "gene_id") %>%
    select(gene_id, transcript_id, gene_name, gene_biotype, source) %>%
    distinct()
}
early_ann     <- ann_lookup(only_6h)
late_ann      <- ann_lookup(only_24h)
sustained_ann <- ann_lookup(intersect_genes)

write_phase_biotype <- function(ann_df, phase_prefix) {
  write_csv(is_mrna_df(ann_df),  file.path(OUTPUT_DIR, paste0(phase_prefix, "_mRNA.csv")))
  write_csv(is_lnc_df(ann_df),   file.path(OUTPUT_DIR, paste0(phase_prefix, "_lncRNA.csv")))
  write_csv(is_novel_df(ann_df), file.path(OUTPUT_DIR, paste0(phase_prefix, "_novel_lncRNA.csv")))
}
write_phase_biotype(early_ann,     "early")
write_phase_biotype(late_ann,      "late")
write_phase_biotype(sustained_ann, "sustained")

# Backward-compatible generic phase sets (gene_id + direction)
write_csv(sig6  %>% select(gene_id, gene_name, direction), file.path(OUTPUT_DIR, "DE_genes_6h.csv"))
write_csv(sig24 %>% select(gene_id, gene_name, direction), file.path(OUTPUT_DIR, "DE_genes_24h.csv"))
write_csv(tibble(gene_id = intersect_genes), file.path(OUTPUT_DIR, "DE_genes_common_sustained.csv"))

cat("[INFO] DE summary  Early:", length(only_6h),
    "| Sustained:", length(intersect_genes),
    "| Late:", length(only_24h), "\n")

# ============================================================================
# 8. Venn diagrams by biotype (mRNA / annotated lncRNA / novel lncRNA)
# ============================================================================

if (requireNamespace("ggvenn", quietly = TRUE)) {
  mk_venn <- function(s6, s24, title, fname) {
    p <- ggvenn::ggvenn(list(`6h` = s6, `24h` = s24),
                        fill_alpha = 0.6, stroke_size = 0.7, set_name_size = 5) +
      ggtitle(title)
    ggsave(file.path(OUTPUT_DIR, fname), p, width = 5, height = 5, dpi = DPI_PNG)
  }
  mk_venn(is_mrna_df(sig6)$gene_id,  is_mrna_df(sig24)$gene_id,  "mRNA DE",  "Venn_mRNA_6h_24h.png")
  mk_venn(is_lnc_df(sig6)$gene_id,   is_lnc_df(sig24)$gene_id,   "lncRNA DE","Venn_lncRNA_6h_24h.png")
  mk_venn(is_novel_df(sig6)$gene_id, is_novel_df(sig24)$gene_id, "novel lncRNA DE", "Venn_novel_lncRNA_6h_24h.png")
  cat("[INFO] Venn diagrams saved\n")
} else {
  cat("[WARN] ggvenn not installed; skipping Venn diagrams\n")
}

# ============================================================================
# 9. Kinetics: LRT for time:treatment interaction
# ============================================================================

cat("[INFO] Computing kinetic genes (LRT)...\n")
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~time + treatment)
res_lrt <- results(dds_lrt)
res_lrt_sig <- res_lrt[which(res_lrt$padj < ALPHA & !is.na(res_lrt$padj)), ]
write_csv(as.data.frame(res_lrt_sig) %>% rownames_to_column("gene_id"),
          file.path(OUTPUT_DIR, "kinetic_genes_lrt.csv"))
cat("[INFO] Kinetic genes (padj <", ALPHA, "):", nrow(res_lrt_sig), "\n")

top_n <- min(50, nrow(res_lrt_sig))
if (top_n > 1) {
  top_kinetic <- rownames(res_lrt_sig)[order(res_lrt_sig$padj)[1:top_n]]
  mat <- assay(vsd)[top_kinetic, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  ann_col <- as.data.frame(colData(dds)[, c("treatment", "time")])
  pdf(file.path(OUTPUT_DIR, "kinetic_genes_heatmap.pdf"), width = 10, height = 12)
  pheatmap(mat, annotation_col = ann_col, main = "Top Kinetic Genes (LRT)",
           fontsize_row = 8, fontsize_col = 10, scale = "row")
  dev.off()
}

# ============================================================================
# 10. Mfuzz: soft clustering for temporal patterns
# ============================================================================

if (nrow(res_lrt_sig) >= MFUZZ_CLUSTERS) {
  cat("[INFO] Computing Mfuzz soft clustering...\n")
  expr_kinetic <- assay(vsd)[rownames(res_lrt_sig), , drop = FALSE]
  expr_kinetic <- expr_kinetic[complete.cases(expr_kinetic), , drop = FALSE]

  pheno <- as.data.frame(colData(dds)[, c("treatment", "time", "group")])
  pheno <- pheno[colnames(expr_kinetic), ]
  stopifnot(all(rownames(pheno) == colnames(expr_kinetic)))

  eset <- ExpressionSet(assayData = expr_kinetic)
  pData(eset) <- pheno
  eset_std <- standardise(eset)
  m_est <- mestimate(eset_std)

  set.seed(SEED)
  cl <- mfuzz(eset_std, c = MFUZZ_CLUSTERS, m = m_est)

  write_csv(as.data.frame(cl$membership) %>% rownames_to_column("gene_id"),
            file.path(OUTPUT_DIR, "mfuzz_membership.csv"))
  write_csv(data.frame(gene_id = names(cl$cluster), cluster = cl$cluster),
            file.path(OUTPUT_DIR, "mfuzz_cluster_assignment.csv"))

  pdf(file.path(OUTPUT_DIR, "mfuzz_clusters.pdf"), width = 12, height = 10)
  mfuzz.plot(eset_std, cl = cl, mfrow = c(3, 2),
             time.labels = as.character(pData(eset_std)$group))
  dev.off()
  cat("[INFO] Mfuzz completed with", MFUZZ_CLUSTERS, "clusters\n")
} else {
  cat("[WARN] Too few kinetic genes for Mfuzz; skipping\n")
}

cat("\n[SUCCESS] DESeq2 + biotype split + Mfuzz analysis complete!\n")
cat("[OUTPUT] Results saved to:", OUTPUT_DIR, "\n")
