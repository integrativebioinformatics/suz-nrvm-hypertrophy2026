library(tximport)
library(rtracklayer)
library(DESeq2)
library(apeglm)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(Mfuzz)
library(Biobase)

# 1. Generar tx2gene
gtf_file <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis/merged_CNE_6_24_HISAT2.annotated.gtf"
tx2gene_file <- "tx2gene.tsv"

if (!file.exists(tx2gene_file)) {
  gtf_data <- rtracklayer::import(gtf_file)
  gtf_df <- as.data.frame(gtf_data)
  tx2gene <- gtf_df %>% filter(type == "transcript") %>% select(transcript_id, gene_id) %>% distinct()
  tx2gene <- na.omit(tx2gene)
  write_tsv(tx2gene, tx2gene_file)
} 

# 2. Metadata
base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/salmon_quants_NE624/"

samples <- c("Ctrl6_rep1","Ctrl6_rep2","Ctrl6_rep4",
             "NE6_rep1","NE6_rep2","NE6_rep4",
             "Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
             "NE24_rep1","NE24_rep2","NE24_rep4")

sampleTable <- data.frame(
  sample = samples,
  treatment = factor(rep(c("Control","NE","Control","NE"), each=3)),
  time = factor(rep(c("6h","24h"), each=6))
)
rownames(sampleTable) <- samples

sampleTable$treatment <- relevel(sampleTable$treatment, "Control")
sampleTable$time <- relevel(sampleTable$time, "6h")

# 👉 Agregar esta línea:
sampleTable$group <- paste(sampleTable$treatment, sampleTable$time, sep="_")


files <- file.path(base_dir, paste0(samples, "_quant"), "quant.sf")
names(files) <- samples
if (!all(file.exists(files))) stop("Faltan archivos quant.sf")

tx2gene <- read_tsv(tx2gene_file, col_names=c("transcript_id","gene_id"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# 3. Crear DESeq2 multifactorial
dds <- DESeqDataSetFromTximport(txi, colData=sampleTable,
                                design = ~ time + treatment + time:treatment)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)

# 4. QC con PCA
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group, label=name)) +
  geom_point(size=3) +
  ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA – NE vs Control por tiempo")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, col=colors, main="Sample distance heatmap")

# 5. DE por tiempo (early y late)
output_dir <- "DESeq2_Multifactorial_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 5A. 6h: NE vs Ctrl
dds_6 <- dds[, dds$time=="6h"]
dds_6$time <- droplevels(dds_6$time)
design(dds_6) <- ~ treatment
dds_6 <- DESeq(dds_6)

res_NE_vs_Ctrl_6h <- results(dds_6, contrast=c("treatment","NE","Control"))
res_NE_vs_Ctrl_6h <- lfcShrink(dds_6, coef="treatment_NE_vs_Control", type="apeglm")
write.csv(as.data.frame(res_NE_vs_Ctrl_6h),
          file.path(output_dir,"res_NE_vs_Ctrl_6h_simple.csv"))

# 5B. 24h: NE vs Ctrl
dds_24 <- dds[, dds$time=="24h"]
dds_24$time <- droplevels(dds_24$time)
design(dds_24) <- ~ treatment
dds_24 <- DESeq(dds_24)

res_NE_vs_Ctrl_24h <- results(dds_24, contrast=c("treatment","NE","Control"))
res_NE_vs_Ctrl_24h <- lfcShrink(dds_24, coef="treatment_NE_vs_Control", type="apeglm")
write.csv(as.data.frame(res_NE_vs_Ctrl_24h),
          file.path(output_dir,"res_NE_vs_Ctrl_24h_simple.csv"))

# 6. Contraste 24h desde el diseño factorial (efecto combinado)
resultsNames(dds)

res_NE_vs_Ctrl_24h_factorial <- results(
  dds,
  contrast = list(c("treatment_NE_vs_Control", "time24h.treatmentNE")),
  alpha = 0.05
)

write.csv(as.data.frame(res_NE_vs_Ctrl_24h_factorial),
          file.path(output_dir, "res_NE_vs_Ctrl_24h_factorial.csv"))

# 7. LRT – CINÉTICA (interacción time:treatment)
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ time + treatment)
res_lrt <- results(dds_lrt)
res_lrt_sig <- res_lrt[which(res_lrt$padj < 0.05 & !is.na(res_lrt$padj)), ]
write.csv(as.data.frame(res_lrt_sig),
          file.path(output_dir,"res_LRT_all_kinetic_genes.csv"))

top_n <- min(50, nrow(res_lrt_sig))
top_kinetic <- rownames(res_lrt_sig)[order(res_lrt_sig$padj)[1:top_n]]

mat <- assay(vsd)[top_kinetic, ]
mat <- mat - rowMeans(mat)
ann_col <- as.data.frame(colData(dds)[,c("treatment","time")])

pdf(file.path(output_dir, "LRT_top_kinetic_heatmap.pdf"), 
    width = 15, height = 20)  # puedes ajustar tamaño

pheatmap(mat,
         annotation_col = ann_col,
         main = "Top kinetic genes (LRT)",
         fontsize_row = 8,        # ajusta tamaño de letra
         fontsize_col = 10,
         scale = "row")

dev.off()

png(file.path(output_dir, "LRT_top_kinetic_heatmap.png"),
    width = 3000, height = 4000, res = 300)  # tamaño gigante

pheatmap(mat,
         annotation_col = ann_col,
         main = "Top kinetic genes (LRT)",
         fontsize_row = 8,
         fontsize_col = 10,
         scale = "row")

dev.off()

# 8. Preparación para MFuzz
expr_vsd <- assay(vsd)

kinetic_genes <- rownames(res_lrt_sig)
expr_kinetic <- expr_vsd[kinetic_genes, ]

# PhenoData: usamos las mismas columnas del dds
pheno <- as.data.frame(colData(dds)[, c("treatment", "time", "group")])
# Asegurarnos que las filas del pheno coinciden con las columnas de la matriz
pheno <- pheno[colnames(expr_kinetic), ]
stopifnot(all(rownames(pheno) == colnames(expr_kinetic)))

eset <- ExpressionSet(assayData = expr_kinetic)
pData(eset) <- pheno

# Estandarizar expresión
eset_std <- standardise(eset)

# Estimar parámetro m
m_est <- mestimate(eset_std)

# Número de clusters
c <- 6
set.seed(123)
cl <- mfuzz(eset_std, c = c, m = m_est)

write.csv(cl$membership, file.path(output_dir, "Mfuzz_membership.csv"))
write.csv(data.frame(cluster = cl$cluster),
          file.path(output_dir, "Mfuzz_clusters.csv"))

# Etiquetas para el eje X: una por muestra (longitud 12)
time_labels <- as.character(pData(eset_std)$group)
# También podrías usar sólo time si quieres: as.character(pData(eset_std)$time)

pdf(file.path(output_dir, "Mfuzz_clusters_plot.pdf"), width = 10, height = 10)
mfuzz.plot(eset_std,
           cl = cl,
           mfrow = c(3, 2),
           time.labels = time_labels)
dev.off()


# Etiquetas de tiempo (una por muestra)
time_labels <- as.character(pData(eset_std)$group)
# o si prefieres sólo tiempo:
# time_labels <- as.character(pData(eset_std)$time)

mfuzz.plot(eset_std,
           cl = cl,
           mfrow = c(3, 2),
           time.labels = time_labels)

png(file.path(output_dir, "Mfuzz_clusters_plot.png"),
    width = 2000, height = 1600, res = 200)

mfuzz.plot(eset_std,
           cl = cl,
           mfrow = c(3, 2),
           time.labels = time_labels)

dev.off()

n_kinetic <- nrow(res_lrt_sig)
n_kinetic

head(res_lrt_sig[order(res_lrt_sig$padj), ])

sig_6h <- res_NE_vs_Ctrl_6h[which(res_NE_vs_Ctrl_6h$padj < 0.05 &
                                    !is.na(res_NE_vs_Ctrl_6h$padj)), ]
nrow(sig_6h)  # número de genes DE a 6h

# Separar up/down
table(sig_6h$log2FoldChange > 0)
# FALSE = down, TRUE = up

sig_24h <- res_NE_vs_Ctrl_24h[which(res_NE_vs_Ctrl_24h$padj < 0.05 &
                                      !is.na(res_NE_vs_Ctrl_24h$padj)), ]
nrow(sig_24h)  # número de genes DE a 24h

table(sig_24h$log2FoldChange > 0)

table(cl$cluster)

# 9. Definir genes DE y listas para Venn/UpSet

# Umbral de significancia
alpha <- 0.05

# Si quisieras agregar un filtro de fold-change, descomenta y ajusta:
# lfc_cutoff <- 0.58  # ~1.5-fold
# filtro_fc <- abs(res_NE_vs_Ctrl_6h$log2FoldChange) >= lfc_cutoff

# Genes DE a 6h (NE vs Ctrl)
sig_6h <- res_NE_vs_Ctrl_6h[
  which(res_NE_vs_Ctrl_6h$padj < alpha & !is.na(res_NE_vs_Ctrl_6h$padj)), ]

# Genes DE a 24h (NE vs Ctrl)
sig_24h <- res_NE_vs_Ctrl_24h[
  which(res_NE_vs_Ctrl_24h$padj < alpha & !is.na(res_NE_vs_Ctrl_24h$padj)), ]

genes_6h  <- rownames(sig_6h)
genes_24h <- rownames(sig_24h)

cat("Genes DE a 6h:",  length(genes_6h),  "\n")
cat("Genes DE a 24h:", length(genes_24h), "\n")

# Intersección y específicos
intersect_genes <- intersect(genes_6h, genes_24h)
only_6h         <- setdiff(genes_6h,  genes_24h)
only_24h        <- setdiff(genes_24h, genes_6h)

cat("DE en ambos tiempos (sostenidos):", length(intersect_genes), "\n")
cat("DE solo en 6h (early):",          length(only_6h),         "\n")
cat("DE solo en 24h (late):",          length(only_24h),        "\n")

# Guardar listas para usar después (anotación, enriquecimiento, etc.)
write.table(genes_6h,
            file = file.path(output_dir, "genes_DE_6h_all.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(genes_24h,
            file = file.path(output_dir, "genes_DE_24h_all.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(intersect_genes,
            file = file.path(output_dir, "genes_DE_common_6h_24h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(only_6h,
            file = file.path(output_dir, "genes_DE_only_6h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(only_24h,
            file = file.path(output_dir, "genes_DE_only_24h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# 10. Venn simple 6h vs 24h (ggvenn)

# install.packages("ggvenn") # si no lo tienes
library(ggvenn)

venn_list <- list(
  DE_6h  = genes_6h,
  DE_24h = genes_24h
)

p_venn <- ggvenn(venn_list,
                 fill_alpha = 0.6,
                 stroke_size = 0.7,
                 set_name_size = 5)

print(p_venn)  # se ve en el panel de Plots

ggsave(file.path(output_dir, "Venn_DE_6h_24h.png"),
       plot = p_venn, width = 5, height = 5, dpi = 300)


# 11. UpSet para DE 6h vs 24h

# install.packages("UpSetR") # si no lo tienes
library(UpSetR)

upset_data <- fromList(venn_list)

upset(upset_data,
      nsets = 2,
      nintersects = 3,
      order.by = "freq",
      sets = c("DE_6h","DE_24h"),
      mainbar.y.label = "Number of genes",
      sets.x.label = "DEGs by time")

# 12. Opcional: UpSet separando UP y DOWN

up_6h   <- rownames(sig_6h)[sig_6h$log2FoldChange > 0]
down_6h <- rownames(sig_6h)[sig_6h$log2FoldChange < 0]

up_24h   <- rownames(sig_24h)[sig_24h$log2FoldChange > 0]
down_24h <- rownames(sig_24h)[sig_24h$log2FoldChange < 0]

cat("6h: UP =", length(up_6h), "DOWN =", length(down_6h), "\n")
cat("24h: UP =", length(up_24h), "DOWN =", length(down_24h), "\n")

updown_list <- list(
  UP_6h    = up_6h,
  DOWN_6h  = down_6h,
  UP_24h   = up_24h,
  DOWN_24h = down_24h
)

updown_data <- fromList(updown_list)

upset(updown_data,
      nsets = 4,
      nintersects = 10,
      order.by = "freq",
      mainbar.y.label = "Number of genes (UP/DOWN)",
      sets.x.label = "UP/DOWN per time")

# Partimos desde gtf_file (el mismo que ya usaste)
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)

# Nos quedamos solo con transcritos
gtf_tx <- gtf_df %>% dplyr::filter(type == "transcript")

colnames(gtf_tx)
# [1] "seqnames"      "start"         "end"           "width"         "strand"       
# [6] "source"        "type"          "score"         "phase"         "transcript_id"
# [11] "gene_id"       "gene_name"     "xloc"          "ref_gene_id"   "cmp_ref"      
# [16] "class_code"    "tss_id"        "exon_number"   "contained_in"  "cmp_ref_gene" 
# [21] "ref_gene_name"

# Construimos el mapa: si hay ref_gene_id → usarlo, si no, dejamos gene_id (MSTRG o ENSRNOG)
gtf_map <- gtf_tx %>%
  transmute(
    transcript_id,
    original_gene_id   = gene_id,
    original_gene_name = gene_name,
    ref_gene_id   = dplyr::na_if(ref_gene_id,   ""),
    ref_gene_id   = ifelse(ref_gene_id %in% c("-", "."), NA, ref_gene_id),
    ref_gene_name = dplyr::na_if(ref_gene_name, ""),
    ref_gene_name = ifelse(ref_gene_name %in% c("-", "."), NA, ref_gene_name),
    final_gene_id   = ifelse(!is.na(ref_gene_id), ref_gene_id, gene_id),
    final_gene_name = ifelse(!is.na(ref_gene_name), ref_gene_name, gene_name)
  )

# Puedes mirar un resumen rápido:
table(is.na(gtf_map$ref_gene_id))
# FALSE = transcritos con match a Ensembl
# TRUE  = transcritos realmente "novel"

# tx2gene original (lo lees de tu archivo o usas el objeto si sigue en memoria)
tx2gene <- readr::read_tsv(tx2gene_file,
                           col_names = c("transcript_id", "gene_id"))

# Corregimos gene_id usando el mapa del GTF
tx2gene_fixed <- tx2gene %>%
  left_join(gtf_map %>% dplyr::select(transcript_id, final_gene_id),
            by = "transcript_id") %>%
  mutate(
    gene_id = ifelse(!is.na(final_gene_id), final_gene_id, gene_id)
  ) %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()

# Guardar el tx2gene corregido
readr::write_tsv(tx2gene_fixed, "tx2gene_fixed.tsv")

tx2gene_fixed <- readr::read_tsv("tx2gene_fixed.tsv",
                                 col_names = c("transcript_id","gene_id"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene_fixed)


###Nuevo
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggvenn)
library(UpSetR)

# 1. Metadata y archivos
base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/salmon_quants_NE624/"
gtf_file <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis/merged_CNE_6_24_HISAT2.annotated.gtf"
output_dir <- "DESeq2_Multifactorial_Results_fixed"
if (!dir.exists(output_dir)) dir.create(output_dir)

samples <- c("Ctrl6_rep1","Ctrl6_rep2","Ctrl6_rep4",
             "NE6_rep1","NE6_rep2","NE6_rep4",
             "Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
             "NE24_rep1","NE24_rep2","NE24_rep4")

sampleTable <- data.frame(
  sample    = samples,
  treatment = factor(rep(c("Control","NE","Control","NE"), each = 3)),
  time      = factor(rep(c("6h","24h"), each = 6))
)
rownames(sampleTable) <- samples

sampleTable$treatment <- relevel(sampleTable$treatment, "Control")
sampleTable$time      <- relevel(sampleTable$time, "6h")
sampleTable$group     <- paste(sampleTable$treatment, sampleTable$time, sep = "_")

files <- file.path(base_dir, paste0(samples, "_quant"), "quant.sf")
names(files) <- samples
if (!all(file.exists(files))) stop("Faltan archivos quant.sf")

# 2. tximport usando tx2gene_fixed
tx2gene_fixed <- read_tsv("tx2gene_fixed.tsv",
                          col_names = c("transcript_id","gene_id"))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene_fixed)

dds <- DESeqDataSetFromTximport(txi,
                                colData = sampleTable,
                                design = ~ time + treatment + time:treatment)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# 3. QC básico: PCA y distancia entre muestras
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA – NE vs Control por tiempo")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, col = colors, main = "Sample distance heatmap")

# 4. DE por tiempo: NE vs Control a 6h y 24h (early / late)
dds_6  <- dds[, dds$time == "6h"]
dds_6$time <- droplevels(dds_6$time)
design(dds_6) <- ~ treatment
dds_6 <- DESeq(dds_6)

dds_24 <- dds[, dds$time == "24h"]
dds_24$time <- droplevels(dds_24$time)
design(dds_24) <- ~ treatment
dds_24 <- DESeq(dds_24)

res_NE_vs_Ctrl_6h <- results(dds_6,  contrast = c("treatment","NE","Control"))
res_NE_vs_Ctrl_24h <- results(dds_24, contrast = c("treatment","NE","Control"))

write.csv(as.data.frame(res_NE_vs_Ctrl_6h),
          file.path(output_dir, "res_NE_vs_Ctrl_6h.csv"))
write.csv(as.data.frame(res_NE_vs_Ctrl_24h),
          file.path(output_dir, "res_NE_vs_Ctrl_24h.csv"))

# 5. Definir listas DE y Venn early/late/sostenidos
alpha <- 0.05

sig_6h <- res_NE_vs_Ctrl_6h[
  which(res_NE_vs_Ctrl_6h$padj < alpha & !is.na(res_NE_vs_Ctrl_6h$padj)), ]
sig_24h <- res_NE_vs_Ctrl_24h[
  which(res_NE_vs_Ctrl_24h$padj < alpha & !is.na(res_NE_vs_Ctrl_24h$padj)), ]

genes_6h  <- rownames(sig_6h)
genes_24h <- rownames(sig_24h)

cat("Genes DE a 6h:",  length(genes_6h),  "\n")
cat("Genes DE a 24h:", length(genes_24h), "\n")

intersect_genes <- intersect(genes_6h, genes_24h)
only_6h         <- setdiff(genes_6h,  genes_24h)
only_24h        <- setdiff(genes_24h, genes_6h)

cat("DE en ambos tiempos (sostenidos):", length(intersect_genes), "\n")
cat("DE solo en 6h (early):",          length(only_6h),         "\n")
cat("DE solo en 24h (late):",          length(only_24h),        "\n")

write.table(genes_6h,
            file = file.path(output_dir, "genes_DE_6h_all.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_24h,
            file = file.path(output_dir, "genes_DE_24h_all.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intersect_genes,
            file = file.path(output_dir, "genes_DE_common_6h_24h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(only_6h,
            file = file.path(output_dir, "genes_DE_only_6h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(only_24h,
            file = file.path(output_dir, "genes_DE_only_24h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

venn_list <- list(
  DE_6h  = genes_6h,
  DE_24h = genes_24h
)

p_venn <- ggvenn(venn_list,
                 fill_alpha = 0.6,
                 stroke_size = 0.7,
                 set_name_size = 5)
print(p_venn)
ggsave(file.path(output_dir, "Venn_DE_6h_24h.png"),
       plot = p_venn, width = 5, height = 5, dpi = 300)

upset_data <- fromList(venn_list)
upset(upset_data,
      nsets = 2,
      nintersects = 3,
      order.by = "freq",
      sets = c("DE_6h","DE_24h"),
      mainbar.y.label = "Number of genes",
      sets.x.label = "DEG per time")

# 12. Opcional: UpSet separando UP y DOWN

up_6h   <- rownames(sig_6h)[sig_6h$log2FoldChange > 0]
down_6h <- rownames(sig_6h)[sig_6h$log2FoldChange < 0]

up_24h   <- rownames(sig_24h)[sig_24h$log2FoldChange > 0]
down_24h <- rownames(sig_24h)[sig_24h$log2FoldChange < 0]

cat("6h: UP =", length(up_6h), "DOWN =", length(down_6h), "\n")
cat("24h: UP =", length(up_24h), "DOWN =", length(down_24h), "\n")

updown_list <- list(
  UP_6h    = up_6h,
  DOWN_6h  = down_6h,
  UP_24h   = up_24h,
  DOWN_24h = down_24h
)

updown_data <- fromList(updown_list)

upset(updown_data,
      nsets = 4,
      nintersects = 10,
      order.by = "freq",
      mainbar.y.label = "Number of genes (UP/DOWN)",
      sets.x.label = "UP/DOWN per time")

# 6. Cargar anotaciones curadas y armonizar columnas
ann_all <- read.csv("hisat2_6_24_meta_ann.csv", stringsAsFactors = FALSE)
ann_novel_raw <- read.csv("hisat2_6_24_meta_novel_lncRNA.csv", stringsAsFactors = FALSE)

ann_all_clean <- ann_all %>%
  transmute(
    gene_id           = ensembl_gene_id,
    transcript_id     = ensembl_transcript_id,
    gene_name         = external_gene_name,
    gene_biotype      = gene_biotype,
    transcript_biotype = transcript_biotype,
    source            = "annotated"
  )

ann_novel_clean <- ann_novel_raw %>%
  filter(type == "transcript") %>%
  transmute(
    gene_id           = gene_id,        # MSTRG.XX
    transcript_id     = transcript_id,  # MSTRG.XX.Y
    gene_name         = qry_gene_id,    # o gene_id
    gene_biotype      = pred_biotype,   # lncRNA
    transcript_biotype = pred_biotype,
    source            = "novel"
  )

ann_merged <- bind_rows(ann_all_clean, ann_novel_clean) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Genes en anotación unificada:", nrow(ann_merged), "\n")

# 7. Anotar listas DE (early, late, sostenidos)
df_6h      <- data.frame(gene_id = genes_6h)
df_24h     <- data.frame(gene_id = genes_24h)
df_common  <- data.frame(gene_id = intersect_genes)
df_only6   <- data.frame(gene_id = only_6h)
df_only24  <- data.frame(gene_id = only_24h)

df_6h_ann      <- df_6h     %>% left_join(ann_merged, by = "gene_id")
df_24h_ann     <- df_24h    %>% left_join(ann_merged, by = "gene_id")
df_common_ann  <- df_common %>% left_join(ann_merged, by = "gene_id")
df_only6_ann   <- df_only6  %>% left_join(ann_merged, by = "gene_id")
df_only24_ann  <- df_only24 %>% left_join(ann_merged, by = "gene_id")

write.csv(df_6h_ann,     file.path(output_dir, "DE_6h_annotated.csv"),     row.names = FALSE)
write.csv(df_24h_ann,    file.path(output_dir, "DE_24h_annotated.csv"),    row.names = FALSE)
write.csv(df_common_ann, file.path(output_dir, "DE_common_annotated.csv"), row.names = FALSE)
write.csv(df_only6_ann,  file.path(output_dir, "DE_only6_annotated.csv"),  row.names = FALSE)
write.csv(df_only24_ann, file.path(output_dir, "DE_only24_annotated.csv"), row.names = FALSE)

cat("Sin anotación en 6h:",  sum(is.na(df_6h_ann$gene_biotype)),  "\n")
cat("Sin anotación en 24h:", sum(is.na(df_24h_ann$gene_biotype)), "\n")

# 8. Separar por biotipo (mRNA, lncRNA conocidas, lncRNA novel)
ann_merged$gene_biotype <- tolower(ann_merged$gene_biotype)

is_mrna <- function(df) {
  df %>% filter(gene_biotype == "protein_coding")
}

is_lncRNA <- function(df) {
  df %>% filter(grepl("lnc", gene_biotype))
}

is_novel <- function(df) {
  df %>% filter(source == "novel")
}

only6_mrna   <- is_mrna(df_only6_ann)
only6_lnc    <- is_lncRNA(df_only6_ann)
only6_novel  <- is_novel(df_only6_ann)

only24_mrna  <- is_mrna(df_only24_ann)
only24_lnc   <- is_lncRNA(df_only24_ann)
only24_novel <- is_novel(df_only24_ann)

common_mrna  <- is_mrna(df_common_ann)
common_lnc   <- is_lncRNA(df_common_ann)
common_novel <- is_novel(df_common_ann)

cat("=== Early (6h-only) ===\n")
cat("mRNAs:", nrow(only6_mrna),
    "lncRNAs:", nrow(only6_lnc),
    "novel lncRNAs:", nrow(only6_novel), "\n")

cat("=== Late (24h-only) ===\n")
cat("mRNAs:", nrow(only24_mrna),
    "lncRNAs:", nrow(only24_lnc),
    "novel lncRNAs:", nrow(only24_novel), "\n")

cat("=== Sostenidos (6h & 24h) ===\n")
cat("mRNAs:", nrow(common_mrna),
    "lncRNAs:", nrow(common_lnc),
    "novel lncRNAs:", nrow(common_novel), "\n")

write.csv(only6_mrna,   file.path(output_dir,"early_mRNA.csv"),             row.names = FALSE)
write.csv(only6_lnc,    file.path(output_dir,"early_lncRNA.csv"),           row.names = FALSE)
write.csv(only6_novel,  file.path(output_dir,"early_novel_lncRNA.csv"),     row.names = FALSE)

write.csv(only24_mrna,  file.path(output_dir,"late_mRNA.csv"),              row.names = FALSE)
write.csv(only24_lnc,   file.path(output_dir,"late_lncRNA.csv"),            row.names = FALSE)
write.csv(only24_novel, file.path(output_dir,"late_novel_lncRNA.csv"),      row.names = FALSE)

write.csv(common_mrna,  file.path(output_dir,"sustained_mRNA.csv"),         row.names = FALSE)
write.csv(common_lnc,   file.path(output_dir,"sustained_lncRNA.csv"),       row.names = FALSE)
write.csv(common_novel, file.path(output_dir,"sustained_novel_lncRNA.csv"), row.names = FALSE)

####venss
alpha <- 0.05

# Tabla DE 6h y 24h con anotación y dirección UP/DOWN
res6_tbl <- as.data.frame(res_NE_vs_Ctrl_6h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_merged, by = "gene_id") %>%
  mutate(direction = case_when(
    padj < alpha & log2FoldChange > 0 ~ "UP",
    padj < alpha & log2FoldChange < 0 ~ "DOWN",
    TRUE ~ "NS"
  ))

res24_tbl <- as.data.frame(res_NE_vs_Ctrl_24h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_merged, by = "gene_id") %>%
  mutate(direction = case_when(
    padj < alpha & log2FoldChange > 0 ~ "UP",
    padj < alpha & log2FoldChange < 0 ~ "DOWN",
    TRUE ~ "NS"
  ))

# Funciones de biotipo (como antes)
ann_merged$gene_biotype <- tolower(ann_merged$gene_biotype)

is_mrna_vec <- function(df) {
  df %>%
    filter(gene_biotype == "protein_coding", padj < alpha, !is.na(padj)) %>%
    pull(gene_id) %>% unique()
}

is_lnc_vec <- function(df) {
  df %>%
    filter(grepl("lnc", gene_biotype), padj < alpha, !is.na(padj)) %>%
    pull(gene_id) %>% unique()
}

is_novel_vec <- function(df) {
  df %>%
    filter(source == "novel", padj < alpha, !is.na(padj)) %>%
    pull(gene_id) %>% unique()
}

# Listas para Venn
mRNA_6  <- is_mrna_vec(res6_tbl)
mRNA_24 <- is_mrna_vec(res24_tbl)

lnc_6   <- is_lnc_vec(res6_tbl)
lnc_24  <- is_lnc_vec(res24_tbl)

novel_6  <- is_novel_vec(res6_tbl)
novel_24 <- is_novel_vec(res24_tbl)

cat("mRNA DE 6h:",  length(mRNA_6),  " | 24h:", length(mRNA_24),  "\n")
cat("lncRNA DE 6h:", length(lnc_6),  " | 24h:", length(lnc_24),   "\n")
cat("novel lncRNA DE 6h:", length(novel_6), " | 24h:", length(novel_24), "\n")

# Venn mRNA
venn_mRNA <- list(
  mRNA_6h  = mRNA_6,
  mRNA_24h = mRNA_24
)
p_venn_mRNA <- ggvenn(venn_mRNA,
                      fill_alpha = 0.6,
                      stroke_size = 0.7,
                      set_name_size = 5)
print(p_venn_mRNA)
ggsave(file.path(output_dir, "Venn_mRNA_6h_24h.png"),
       plot = p_venn_mRNA, width = 5, height = 5, dpi = 300)

# Venn lncRNAs conocidas
venn_lnc <- list(
  lnc_6h  = lnc_6,
  lnc_24h = lnc_24
)
p_venn_lnc <- ggvenn(venn_lnc,
                     fill_alpha = 0.6,
                     stroke_size = 0.7,
                     set_name_size = 5)
print(p_venn_lnc)
ggsave(file.path(output_dir, "Venn_lncRNA_6h_24h.png"),
       plot = p_venn_lnc, width = 5, height = 5, dpi = 300)

# Venn novel lncRNAs
venn_novel <- list(
  novel_6h  = novel_6,
  novel_24h = novel_24
)
p_venn_novel <- ggvenn(venn_novel,
                       fill_alpha = 0.6,
                       stroke_size = 0.7,
                       set_name_size = 5)
print(p_venn_novel)
ggsave(file.path(output_dir, "Venn_novel_lncRNA_6h_24h.png"),
       plot = p_venn_novel, width = 5, height = 5, dpi = 300)

# Helper para exportar tablas limpias
export_DE_table <- function(df, biotype_fun, time_label, biotype_label, outname_prefix) {
  out_df <- df %>%
    filter(padj < alpha, !is.na(padj)) %>%
    biotype_fun() %>%
    select(gene_id, gene_name, log2FoldChange, padj, direction, gene_biotype, source) %>%
    arrange(direction, padj)
  
  message(time_label, " - ", biotype_label, ": ", nrow(out_df), " genes")
  
  write.csv(out_df,
            file.path(output_dir, paste0(outname_prefix, "_", time_label, "_", biotype_label, ".csv")),
            row.names = FALSE)
}

# Biotype helpers que operan sobre df ya filtrado en padj
filter_mrna_tbl <- function(df) df %>% filter(gene_biotype == "protein_coding")
filter_lnc_tbl  <- function(df) df %>% filter(grepl("lnc", gene_biotype))
filter_novel_tbl <- function(df) df %>% filter(source == "novel")

# 6h
export_DE_table(res6_tbl, filter_mrna_tbl,  "6h", "mRNA",       "DE")
export_DE_table(res6_tbl, filter_lnc_tbl,   "6h", "lncRNA",     "DE")
export_DE_table(res6_tbl, filter_novel_tbl, "6h", "novel_lnc",  "DE")

# 24h
export_DE_table(res24_tbl, filter_mrna_tbl,  "24h", "mRNA",      "DE")
export_DE_table(res24_tbl, filter_lnc_tbl,   "24h", "lncRNA",    "DE")
export_DE_table(res24_tbl, filter_novel_tbl, "24h", "novel_lnc", "DE")


##Mfuzz
# === LRT: Genes con cinética significativa ===
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ time + treatment)
res_lrt <- results(dds_lrt)
res_lrt_sig <- res_lrt[which(res_lrt$padj < 0.05 & !is.na(res_lrt$padj)), ]

cat("Genes con respuesta cinética (LRT, padj < 0.05):", nrow(res_lrt_sig), "\n")

write.csv(as.data.frame(res_lrt_sig),
          file.path(output_dir, "res_LRT_all_kinetic_genes.csv"),
          row.names = TRUE)

# === Preparar matriz vst solamente con genes cinéticos ===
kinetic_genes <- rownames(res_lrt_sig)
expr_vsd <- assay(vsd)

expr_kinetic <- expr_vsd[kinetic_genes, ]

# Filtrar posibles NA o filas vacías
expr_kinetic <- expr_kinetic[complete.cases(expr_kinetic), ]

# === Crear ExpressionSet para MFuzz ===
library(Biobase)
library(Mfuzz)

pheno <- as.data.frame(colData(dds)[, c("treatment", "time")])

eset <- ExpressionSet(assayData = expr_kinetic)
pData(eset) <- pheno

# Standardize
eset_std <- standardise(eset)

# Estimate fuzzification parameter m
m_est <- mestimate(eset_std)

# Elegir número de clusters (seis es razonable, pero puedes ajustar)
c <- 6
set.seed(123)

cl <- mfuzz(eset_std, c = c, m = m_est)

# Guardar información
write.csv(cl$membership, file.path(output_dir, "Mfuzz_membership.csv"))
write.csv(data.frame(gene_id = rownames(cl$membership),
                     cluster = cl$cluster),
          file.path(output_dir, "Mfuzz_clusters.csv"),
          row.names = FALSE)

# === Graficar clusters ===
png(file.path(output_dir, "Mfuzz_clusters_plot.png"), width = 1800, height = 1800, res = 200)
mfuzz.plot(eset_std, cl = cl, mfrow = c(3, 2),
           time.labels = colnames(expr_kinetic))  # usa nombres de columnas reales
dev.off()

cat("MFuzz completado. Figuras y tablas exportadas.\n")

##mfuzz por cluster
# Asegurarnos de tener la misma fila en expr_kinetic y en cl
kinetic_genes <- rownames(expr_kinetic)

cluster_df <- data.frame(
  gene_id = kinetic_genes,
  cluster = cl$cluster[kinetic_genes]
)

# Membership máximo (cuán "fuerte" es la asignación al clúster)
max_membership <- apply(cl$membership[kinetic_genes, ], 1, max)
cluster_df$membership_max <- max_membership

# Anotación
cluster_df <- cluster_df %>%
  left_join(ann_merged, by = "gene_id")

# Resumen por clúster y biotipo grueso
cluster_summary <- cluster_df %>%
  mutate(
    biotype_group = case_when(
      gene_biotype == "protein_coding"              ~ "mRNA",
      grepl("lnc", gene_biotype)                    ~ "lncRNA",
      source == "novel"                             ~ "novel_lncRNA",
      TRUE                                          ~ "other"
    )
  ) %>%
  group_by(cluster, biotype_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(cluster, desc(n))

write.csv(cluster_summary,
          file.path(output_dir, "Mfuzz_cluster_biotype_summary.csv"),
          row.names = FALSE)

cluster_summary

##Hagamos la parte de expresión media por cluster
library(dplyr)

expr_vsd <- expr_kinetic  # matriz vst de genes cinéticos (filas = genes, cols = muestras)

# Info de condiciones por muestra
cond_info <- data.frame(
  sample    = colnames(expr_vsd),
  treatment = as.character(colData(dds)$treatment),
  time      = as.character(colData(dds)$time)
)
cond_info$group <- paste(cond_info$treatment, cond_info$time, sep = "_")
group_levels <- unique(cond_info$group)

# Medias por gen y grupo (Control_6h, NE_6h, etc.)
gene_group_means <- sapply(group_levels, function(g) {
  cols <- cond_info$sample[cond_info$group == g]
  rowMeans(expr_vsd[, cols, drop = FALSE])
})
rownames(gene_group_means) <- rownames(expr_vsd)

# Tabla de clusters si aún no la tienes
cluster_df <- data.frame(
  gene_id = rownames(expr_vsd),
  cluster = cl$cluster[rownames(expr_vsd)],
  membership_max = apply(cl$membership[rownames(expr_vsd), , drop = FALSE], 1, max)
) %>%
  left_join(ann_merged, by = "gene_id")

# Medias por cluster (promedio de medias de genes)
cluster_ids <- split(cluster_df$gene_id, cluster_df$cluster)

cluster_means <- lapply(cluster_ids, function(genes) {
  colMeans(gene_group_means[genes, , drop = FALSE], na.rm = TRUE)
})
cluster_means_mat <- do.call(rbind, cluster_means)
rownames(cluster_means_mat) <- names(cluster_ids)
colnames(cluster_means_mat) <- group_levels

# Identificar nombres de grupos
g_Ctrl6  <- grep("Control_6h", colnames(cluster_means_mat),  value = TRUE)
g_NE6    <- grep("NE_6h",      colnames(cluster_means_mat),  value = TRUE)
g_Ctrl24 <- grep("Control_24h",colnames(cluster_means_mat),  value = TRUE)
g_NE24   <- grep("NE_24h",     colnames(cluster_means_mat),  value = TRUE)

early_effect <- cluster_means_mat[, g_NE6]  - cluster_means_mat[, g_Ctrl6]
late_effect  <- cluster_means_mat[, g_NE24] - cluster_means_mat[, g_Ctrl24]

cluster_class <- data.frame(
  cluster       = as.integer(rownames(cluster_means_mat)),
  early_effect  = as.numeric(early_effect),
  late_effect   = as.numeric(late_effect)
)

# Clasificación heurística (ajusta el umbral si quieres)
thr <- 0.5  # umbral de efecto

cluster_class$pattern <- dplyr::case_when(
  abs(early_effect) >= thr & abs(late_effect) < thr/2              ~ "early",
  abs(late_effect)  >= thr & abs(early_effect) < thr/2             ~ "late",
  abs(early_effect) >= thr & abs(late_effect) >= thr               ~ "sustained",
  TRUE                                                             ~ "other"
)

cluster_class <- cluster_class %>% arrange(cluster)
cluster_class

write.csv(cluster_class,
          file.path(output_dir, "Mfuzz_cluster_temporal_pattern.csv"),
          row.names = FALSE)

# IDs de lncRNAs según DE (early / late / sostenidos)
early_lnc_ids <- df_only6_ann  %>%
  filter(grepl("lnc", tolower(gene_biotype))) %>%
  pull(gene_id) %>% unique()

late_lnc_ids <- df_only24_ann %>%
  filter(grepl("lnc", tolower(gene_biotype))) %>%
  pull(gene_id) %>% unique()

sust_lnc_ids <- df_common_ann %>%
  filter(grepl("lnc", tolower(gene_biotype))) %>%
  pull(gene_id) %>% unique()

# Grupo de biotipo (mRNA / lncRNA / novel_lncRNA / other)
cluster_df <- cluster_df %>%
  mutate(
    biotype_group = case_when(
      gene_biotype == "protein_coding"               ~ "mRNA",
      source == "novel"                              ~ "novel_lncRNA",
      grepl("lnc", tolower(gene_biotype))            ~ "lncRNA",
      TRUE                                           ~ "other"
    ),
    time_category = case_when(
      gene_id %in% early_lnc_ids ~ "early_only",
      gene_id %in% late_lnc_ids  ~ "late_only",
      gene_id %in% sust_lnc_ids  ~ "sustained",
      TRUE                       ~ "non_DE_or_other"
    )
  ) %>%
  left_join(cluster_class, by = "cluster")

# Solo lncRNAs (conocidas o novel)
cluster_lncs <- cluster_df %>%
  filter(biotype_group %in% c("lncRNA","novel_lncRNA"))

# Resumen: cuántas lncRNAs por cluster y categoría temporal
lnc_cluster_summary <- cluster_lncs %>%
  group_by(cluster, pattern, biotype_group, time_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(cluster, pattern, biotype_group, desc(n))

write.csv(lnc_cluster_summary,
          file.path(output_dir, "Mfuzz_lncRNA_cluster_temporal_summary.csv"),
          row.names = FALSE)

lnc_cluster_summary

write.csv(cluster_lncs,
          file.path(output_dir, "Mfuzz_lncRNA_cluster_gene_list.csv"),
          row.names = FALSE)

library(pheatmap)

# Usamos expr_vsd (expr_kinetic), cluster_df, cond_info ya definidos antes
expr_mat <- expr_vsd  # alias

# Función para obtener top genes por cluster y biotipo
get_top_genes <- function(cluster_id, biotype, n = 10) {
  cluster_df %>%
    filter(cluster == cluster_id,
           biotype_group == biotype,
           gene_id %in% rownames(expr_mat)) %>%
    arrange(desc(membership_max)) %>%
    slice_head(n = n) %>%
    pull(gene_id)
}

all_clusters <- sort(unique(cluster_df$cluster))
biotypes_to_plot <- c("mRNA", "lncRNA", "novel_lncRNA")

for (cl_id in all_clusters) {
  for (bt in biotypes_to_plot) {
    top_genes <- get_top_genes(cl_id, bt, n = 10)
    if (length(top_genes) == 0) next
    
    mat <- expr_mat[top_genes, , drop = FALSE]
    mat <- mat - rowMeans(mat)  # centrar por gen
    
    ann_col <- cond_info %>%
      select(sample, treatment, time) %>%
      tibble::column_to_rownames("sample")
    
    out_name <- paste0("Heatmap_cluster", cl_id, "_", bt, ".png")
    
    png(file.path(output_dir, out_name),
        width = 2000, height = 1600, res = 200)
    pheatmap(mat,
             annotation_col = ann_col,
             main = paste("Cluster", cl_id, "-", bt),
             scale = "row",
             fontsize_row = 8,
             fontsize_col = 10)
    dev.off()
  }
}

top_list <- list()

for (cl_id in all_clusters) {
  for (bt in biotypes_to_plot) {
    top_genes <- get_top_genes(cl_id, bt, n = 10)
    if (length(top_genes) == 0) next
    
    sub <- cluster_df %>%
      filter(gene_id %in% top_genes,
             cluster == cl_id,
             biotype_group == bt) %>%
      arrange(desc(membership_max))
    
    top_list[[length(top_list) + 1]] <- sub
  }
}

top_genes_all <- bind_rows(top_list)

write.csv(top_genes_all,
          file.path(output_dir, "Mfuzz_top10_genes_per_cluster_biotype.csv"),
          row.names = FALSE)

##Análisis de enriquecimiento
#   - res_NE_vs_Ctrl_6h  : objeto DESeq2 resultados (DESeqResults)
#   - res_NE_vs_Ctrl_24h : idem para 24h
#   - ann_merged         : data.frame con columnas gene_id (ENSEMBL) y gene_name

# 1. Cargar librerías --------------------------------------
library(clusterProfiler)
library(org.Rn.eg.db)
library(dplyr)
library(tibble)
library(ggplot2)

# 2. Directorio de salida ----------------------------------
outdir <- "GSEA_results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 3. Función para preparar el vector rankeado (ENTREZ, único) ----
prep_gsea_list <- function(res_de, ann_merged, OrgDb) {
  # res_de: objeto DESeqResults
  # ann_merged: data.frame con gene_id (ENSEMBL) y, opcionalmente, gene_name
  
  # 3.1 Convertir a tabla y agregar anotación
  res_tbl <- as.data.frame(res_de) %>%
    tibble::rownames_to_column("gene_id") %>%
    filter(!is.na(log2FoldChange)) %>%
    left_join(ann_merged %>% dplyr::select(gene_id, gene_name),
              by = "gene_id")
  
  # 3.2 Crear vector inicial: nombres ENSEMBL, valores log2FC
  geneList <- res_tbl$log2FoldChange
  names(geneList) <- res_tbl$gene_id
  geneList <- sort(geneList, decreasing = TRUE)
  
  # 3.3 Mapear ENSEMBL -> ENTREZ
  gene_df <- bitr(names(geneList),
                  fromType = "ENSEMBL",
                  toType   = "ENTREZID",
                  OrgDb    = OrgDb)
  # Puede haber 1:many mappings
  
  # 3.4 Incorporar log2FC y quedarnos con un ENTREZID único
  gene_df2 <- gene_df %>%
    mutate(log2FC = geneList[ENSEMBL]) %>%
    arrange(desc(abs(log2FC))) %>%        # primero los más extremos
    distinct(ENTREZID, .keep_all = TRUE)  # un único ENTREZID
  
  # 3.5 Construir vector final rankeado por ENTREZID
  geneList_entrez <- gene_df2$log2FC
  names(geneList_entrez) <- gene_df2$ENTREZID
  
  geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)
  return(geneList_entrez)
}

# 4. Preparar listas para 6h y 24h -------------------------
geneList_6h_entrez  <- prep_gsea_list(res_NE_vs_Ctrl_6h,  ann_merged, org.Rn.eg.db)
geneList_24h_entrez <- prep_gsea_list(res_NE_vs_Ctrl_24h, ann_merged, org.Rn.eg.db)

# 5. Correr GSEA para GO BP -------------------------------
gsea_bp_6h <- gseGO(
  geneList     = geneList_6h_entrez,
  ont          = "BP",
  OrgDb        = org.Rn.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 800,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_bp_24h <- gseGO(
  geneList     = geneList_24h_entrez,
  ont          = "BP",
  OrgDb        = org.Rn.eg.db,
  keyType      = "ENTREZID",
  minGSSize    = 10,
  maxGSSize    = 800,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# 6. Correr GSEA para KEGG (organismo: rno) ---------------
gsea_kegg_6h <- gseKEGG(
  geneList     = geneList_6h_entrez,
  organism     = "rno",
  keyType      = "ncbi-geneid",  # ENTREZID
  minGSSize    = 10,
  maxGSSize    = 800,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_kegg_24h <- gseKEGG(
  geneList     = geneList_24h_entrez,
  organism     = "rno",
  keyType      = "ncbi-geneid",
  minGSSize    = 10,
  maxGSSize    = 800,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# 7. Función para extraer top rutas activadas/suprimidas ---
get_top_terms <- function(gsea_res, n = 8, p_cutoff = 0.05) {
  df <- as.data.frame(gsea_res)
  if (!"p.adjust" %in% colnames(df)) {
    stop("El objeto GSEA no tiene columna 'p.adjust'. Revisa gsea_res.")
  }
  df <- df %>% filter(p.adjust <= p_cutoff)
  
  up <- df %>%
    filter(NES > 0) %>%
    arrange(p.adjust) %>%
    dplyr::slice_head(n = n)
  
  down <- df %>%
    filter(NES < 0) %>%
    arrange(p.adjust) %>%
    dplyr::slice_head(n = n)
  
  list(up = up, down = down)
}

# 8. Función para hacer dotplots (GO/KEGG, up/down) -------
plot_gsea_dot <- function(df, title) {
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # Ordenar términos para que el dotplot quede ordenado
  df <- df %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  p <- ggplot(df, aes(x = NES, y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    labs(
      title = title,
      x     = "Normalized Enrichment Score (NES)",
      y     = NULL,
      size  = "Set size",
      color = "adj. p"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y  = element_text(size = 9),
      plot.title   = element_text(hjust = 0.5)
    )
  return(p)
}

# 9. Función para guardar PDF + PNG -----------------------
save_dotplot <- function(p, filename_base, outdir, width = 7, height = 5) {
  if (is.null(p)) return(invisible(NULL))
  
  pdf_file <- file.path(outdir, paste0(filename_base, ".pdf"))
  png_file <- file.path(outdir, paste0(filename_base, ".png"))
  
  ggsave(pdf_file, plot = p, width = width, height = height)
  ggsave(png_file, plot = p, width = width, height = height, dpi = 300)
}

# 10. Generar y exportar dotplots -------------------------
## 10.1 GO BP 6h
bp6_top <- get_top_terms(gsea_bp_6h, n = 8, p_cutoff = 0.05)

p_bp6_up <- plot_gsea_dot(bp6_top$up,
                          "GO BP 6h – Top 8 rutas activadas (NES > 0)")
p_bp6_down <- plot_gsea_dot(bp6_top$down,
                            "GO BP 6h – Top 8 rutas suprimidas (NES < 0)")

save_dotplot(p_bp6_up,   "GSEA_GO_BP_6h_Activated",   outdir)
save_dotplot(p_bp6_down, "GSEA_GO_BP_6h_Suppressed", outdir)

## 10.2 GO BP 24h
bp24_top <- get_top_terms(gsea_bp_24h, n = 8, p_cutoff = 0.05)

p_bp24_up <- plot_gsea_dot(bp24_top$up,
                           "GO BP 24h – Top 8 rutas activadas (NES > 0)")
p_bp24_down <- plot_gsea_dot(bp24_top$down,
                             "GO BP 24h – Top 8 rutas suprimidas (NES < 0)")

save_dotplot(p_bp24_up,   "GSEA_GO_BP_24h_Activated",   outdir)
save_dotplot(p_bp24_down, "GSEA_GO_BP_24h_Suppressed", outdir)

## 10.3 KEGG 6h
kegg6_top <- get_top_terms(gsea_kegg_6h, n = 8, p_cutoff = 0.05)

p_kegg6_up <- plot_gsea_dot(kegg6_top$up,
                            "KEGG 6h – Top 8 rutas activadas (NES > 0)")
p_kegg6_down <- plot_gsea_dot(kegg6_top$down,
                              "KEGG 6h – Top 8 rutas suprimidas (NES < 0)")

save_dotplot(p_kegg6_up,   "GSEA_KEGG_6h_Activated",   outdir)
save_dotplot(p_kegg6_down, "GSEA_KEGG_6h_Suppressed", outdir)

## 10.4 KEGG 24h
kegg24_top <- get_top_terms(gsea_kegg_24h, n = 8, p_cutoff = 0.05)

p_kegg24_up <- plot_gsea_dot(kegg24_top$up,
                             "KEGG 24h – Top 8 rutas activadas (NES > 0)")
p_kegg24_down <- plot_gsea_dot(kegg24_top$down,
                               "KEGG 24h – Top 8 rutas suprimidas (NES < 0)")

save_dotplot(p_kegg24_up,   "GSEA_KEGG_24h_Activated",   outdir)
save_dotplot(p_kegg24_down, "GSEA_KEGG_24h_Suppressed", outdir)

# === 10-BIS. Comparación 6h vs 24h para identificar rutas únicas y compartidas ===

compare_times <- function(gsea_6h, gsea_24h, ontology_label, out_prefix) {
  df6 <- as.data.frame(gsea_6h) %>% filter(p.adjust <= 0.05)
  df24 <- as.data.frame(gsea_24h) %>% filter(p.adjust <= 0.05)
  
  ids6  <- df6$ID
  ids24 <- df24$ID
  
  shared   <- intersect(ids6, ids24)
  only_6h  <- setdiff(ids6, ids24)
  only_24h <- setdiff(ids24, ids6)
  
  shared_df <- df6[df6$ID %in% shared, ] %>%
    rename(NES_6h = NES, p.adjust_6h = p.adjust) %>%
    left_join(df24 %>% dplyr::select(ID, NES, p.adjust) %>% 
                rename(NES_24h = NES, p.adjust_24h = p.adjust),
              by = "ID")
  
  only6_df <- df6[df6$ID %in% only_6h, ]
  only24_df <- df24[df24$ID %in% only_24h, ]
  
  write.table(shared_df,
              file = file.path(outdir, paste0(out_prefix, "_Shared.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(only6_df,
              file = file.path(outdir, paste0(out_prefix, "_Only6h.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(only24_df,
              file = file.path(outdir, paste0(out_prefix, "_Only24h.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  list(shared = shared_df, only6 = only6_df, only24 = only24_df)
}

# === GO BP ===
bp_overlap <- compare_times(
  gsea_6h = gsea_bp_6h,
  gsea_24h = gsea_bp_24h,
  ontology_label = "GO_BP",
  out_prefix = "GO_BP_6h_vs_24h"
)

# === KEGG ===
kegg_overlap <- compare_times(
  gsea_6h = gsea_kegg_6h,
  gsea_24h = gsea_kegg_24h,
  ontology_label = "KEGG",
  out_prefix = "KEGG_6h_vs_24h"
)

##plots nuevos
make_gsea_faceted_plot <- function(gsea_res,
                                   geneList_entrez,
                                   time_label,
                                   ontology_label,
                                   out_basename,
                                   top_n    = 8,
                                   p_cutoff = 0.05) {
  df <- as.data.frame(gsea_res) %>%
    dplyr::filter(p.adjust <= p_cutoff)
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  total_genes <- length(geneList_entrez)
  
  up <- df %>%
    dplyr::filter(NES > 0) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(Direction = "activated")
  
  down <- df %>%
    dplyr::filter(NES < 0) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(Direction = "suppressed")
  
  plot_df <- dplyr::bind_rows(up, down)
  if (nrow(plot_df) == 0) return(invisible(NULL))
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      GeneRatio   = setSize / total_genes,
      Description = factor(Description, levels = rev(unique(Description)))
    )
  
  p <- ggplot(plot_df, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    facet_grid(. ~ Direction) +
    labs(
      title = paste0(ontology_label, " ", time_label,
                     " – top ", top_n, " activated/suppressed"),
      x     = "GeneRatio",
      y     = NULL,
      size  = "Count",
      color = "p.adjust"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y       = element_text(size = 9),
      strip.background  = element_rect(fill = "grey90"),
      strip.text        = element_text(face = "bold"),
      plot.title        = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(outdir, paste0(out_basename, ".pdf")),
         plot = p, width = 7, height = 6)
  ggsave(file.path(outdir, paste0(out_basename, ".png")),
         plot = p, width = 7, height = 6, dpi = 300)
}

# GO BP 6h
make_gsea_faceted_plot(
  gsea_res      = gsea_bp_6h,
  geneList_entrez = geneList_6h_entrez,
  time_label    = "6h",
  ontology_label = "GO BP",
  out_basename  = "GSEA_GO_BP_6h"
)

# GO BP 24h
make_gsea_faceted_plot(
  gsea_res      = gsea_bp_24h,
  geneList_entrez = geneList_24h_entrez,
  time_label    = "24h",
  ontology_label = "GO BP",
  out_basename  = "GSEA_GO_BP_24h"
)

# KEGG 6h
make_gsea_faceted_plot(
  gsea_res      = gsea_kegg_6h,
  geneList_entrez = geneList_6h_entrez,
  time_label    = "6h",
  ontology_label = "KEGG",
  out_basename  = "GSEA_KEGG_6h"
)

# KEGG 24h
make_gsea_faceted_plot(
  gsea_res      = gsea_kegg_24h,
  geneList_entrez = geneList_24h_entrez,
  time_label    = "24h",
  ontology_label = "KEGG",
  out_basename  = "GSEA_KEGG_24h"
)

##Correlaciones
# Para ver si los vectores siguen siendo enormes
length(lnc_DE_ids)
length(mrna_DE_ids)

# Ver si la función vieja sigue viva
ls(pattern = "correlate_sets")
ls(pattern = "get_correlations_rcorr")

fast_cor_lnc_mrna <- function(lnc_ids, mrna_ids, expr_mat,
                              rho_cut = 0.8, p_cut = 0.05) {
  # Sub-matrices
  lnc_ids  <- intersect(lnc_ids,  rownames(expr_mat))
  mrna_ids <- intersect(mrna_ids, rownames(expr_mat))
  
  if (length(lnc_ids) == 0 || length(mrna_ids) == 0) {
    message("No hay genes suficientes para correlación.")
    return(data.frame())
  }
  
  X <- t(expr_mat[lnc_ids, , drop = FALSE])   # muestras x lnc
  Y <- t(expr_mat[mrna_ids, , drop = FALSE])  # muestras x mRNAs
  
  # Matriz de correlación lnc x mRNA
  R <- cor(X, Y, method = "pearson", use = "pairwise.complete.obs")
  # R tiene dimensiones: length(lnc_ids) x length(mrna_ids)
  
  n <- nrow(X)  # número de muestras
  # t-statistics y p-values
  Tvals <- R * sqrt((n - 2) / (1 - R^2))
  Pvals <- 2 * pt(-abs(Tvals), df = n - 2)
  
  # Manejar posibles NAs
  R[is.na(R)] <- 0
  Pvals[is.na(Pvals)] <- 1
  
  # Máscara de correlaciones significativas
  keep <- (abs(R) >= rho_cut) & (Pvals < p_cut)
  idx <- which(keep, arr.ind = TRUE)
  
  if (nrow(idx) == 0) {
    message("No se encontraron correlaciones que cumplan los umbrales.")
    return(data.frame())
  }
  
  lnc_sel  <- lnc_ids[idx[, 1]]
  mrna_sel <- mrna_ids[idx[, 2]]
  rho_sel  <- R[keep]
  p_sel    <- Pvals[keep]
  
  out <- data.frame(
    lncRNA  = lnc_sel,
    mRNA    = mrna_sel,
    rho     = rho_sel,
    pvalue  = p_sel
  )
  
  # ordenar por |rho|
  out <- out[order(-abs(out$rho)), ]
  return(out)
}

expr <- assay(vsd)

# Tablas DE (ya las teníamos, pero repito por claridad)
res6 <- as.data.frame(res_NE_vs_Ctrl_6h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_small, by = "gene_id")

res24 <- as.data.frame(res_NE_vs_Ctrl_24h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_small, by = "gene_id")

# lncRNAs DE fuertes (6h o 24h)
lnc_DE <- bind_rows(
  res6 %>% filter(padj < 0.05,
                  grepl("lnc", gene_biotype, ignore.case = TRUE),
                  abs(log2FoldChange) >= 0.5),
  res24 %>% filter(padj < 0.05,
                   grepl("lnc", gene_biotype, ignore.case = TRUE),
                   abs(log2FoldChange) >= 0.5)
) %>% distinct(gene_id)

lnc_DE_ids <- intersect(lnc_DE$gene_id, rownames(expr))

# mRNAs DE fuertes (6h o 24h)
mrna_DE <- bind_rows(
  res6 %>% filter(padj < 0.05,
                  gene_biotype == "protein_coding",
                  abs(log2FoldChange) >= 1),
  res24 %>% filter(padj < 0.05,
                   gene_biotype == "protein_coding",
                   abs(log2FoldChange) >= 1)
) %>% distinct(gene_id)

mrna_DE_ids <- intersect(mrna_DE$gene_id, rownames(expr))

length(lnc_DE_ids); length(mrna_DE_ids)

cor_DE <- fast_cor_lnc_mrna(
  lnc_ids  = lnc_DE_ids,
  mrna_ids = mrna_DE_ids,
  expr_mat = expr,
  rho_cut  = 0.8,
  p_cut    = 0.05
)

write.csv(cor_DE,
          file.path(output_dir, "cor_lncRNA_mRNA_DE_filtered_fastcor.csv"),
          row.names = FALSE)

early_lnc_filt <- intersect(lnc_DE_ids, df_only6_ann$gene_id)
early_mrna_filt <- intersect(mrna_DE_ids, df_only6_ann$gene_id)

cor_early_fast <- fast_cor_lnc_mrna(
  lnc_ids  = early_lnc_filt,
  mrna_ids = early_mrna_filt,
  expr_mat = expr,
  rho_cut  = 0.8,
  p_cut    = 0.05
)

write.csv(cor_early_fast,
          file.path(output_dir, "cor_early_fast_DE_filtered.csv"),
          row.names = FALSE)

late_lnc_filt <- intersect(lnc_DE_ids, df_only24_ann$gene_id)
late_mrna_filt <- intersect(mrna_DE_ids, df_only24_ann$gene_id)

cor_late_fast <- fast_cor_lnc_mrna(
  lnc_ids  = late_lnc_filt,
  mrna_ids = late_mrna_filt,
  expr_mat = expr,
  rho_cut  = 0.8,
  p_cut    = 0.05
)

write.csv(cor_late_fast,
          file.path(output_dir, "cor_late_fast_DE_filtered.csv"),
          row.names = FALSE)

sust_lnc_filt <- intersect(lnc_DE_ids, df_common_ann$gene_id)
sust_mrna_filt <- intersect(mrna_DE_ids, df_common_ann$gene_id)

cor_sust_fast <- fast_cor_lnc_mrna(
  lnc_ids  = sust_lnc_filt,
  mrna_ids = sust_mrna_filt,
  expr_mat = expr,
  rho_cut  = 0.8,
  p_cut    = 0.05
)

write.csv(cor_sust_fast,
          file.path(output_dir, "cor_sustained_fast_DE_filtered.csv"),
          row.names = FALSE)

annotate_cor <- function(cor_table, ann_table) {
  if (nrow(cor_table) == 0) return(cor_table)
  
  cor_table %>%
    left_join(ann_table, by = c("lncRNA" = "gene_id")) %>%
    rename(
      lnc_gene_name     = gene_name,
      lnc_biotype       = gene_biotype,
      lnc_source        = source
    ) %>%
    left_join(ann_table, by = c("mRNA" = "gene_id")) %>%
    rename(
      mRNA_gene_name    = gene_name,
      mRNA_biotype      = gene_biotype,
      mRNA_source       = source
    )
}

cor_early_ann  <- annotate_cor(cor_early_fast, ann_small)
cor_late_ann   <- annotate_cor(cor_late_fast,  ann_small)
cor_sust_ann   <- annotate_cor(cor_sust_fast,  ann_small)

write.csv(cor_early_ann, file.path(output_dir,"cor_early_annotated.csv"), row.names=FALSE)
write.csv(cor_late_ann,  file.path(output_dir,"cor_late_annotated.csv"),  row.names=FALSE)
write.csv(cor_sust_ann,  file.path(output_dir,"cor_sustained_annotated.csv"), row.names=FALSE)

PAIR <- function(df) paste(df$lncRNA, df$mRNA, sep="__")

pairs_early <- PAIR(cor_early_fast)
pairs_late  <- PAIR(cor_late_fast)
pairs_sust  <- PAIR(cor_sust_fast)

pairs_stable <- intersect(pairs_early, pairs_late)
length(pairs_stable)

cor_stable <- cor_early_ann[PAIR(cor_early_ann) %in% pairs_stable,]
write.csv(cor_stable, file.path(output_dir, "correlations_stable_6h_24h.csv"), row.names=FALSE)

pairs_early_only <- setdiff(pairs_early, pairs_late)

cor_early_only <- cor_early_ann[PAIR(cor_early_ann) %in% pairs_early_only,]

write.csv(cor_early_only,
          file.path(output_dir, "correlations_early_only.csv"),
          row.names=FALSE)

pairs_late_only <- setdiff(pairs_late, pairs_early)

cor_late_only <- cor_late_ann[PAIR(cor_late_ann) %in% pairs_late_only,]

write.csv(cor_late_only,
          file.path(output_dir, "correlations_late_only.csv"),
          row.names=FALSE)

###Lost = fuerte correlación en 6h pero NO en 24h
### Gained = aparece en 24h pero NO estaba en 6h

cor_lost <- cor_early_ann[PAIR(cor_early_ann) %in% pairs_early_only,]
cor_gained <- cor_late_ann[PAIR(cor_late_ann) %in% pairs_late_only,]

write.csv(cor_lost, file.path(output_dir,"correlations_lost_after_6h.csv"), row.names=FALSE)
write.csv(cor_gained, file.path(output_dir,"correlations_gained_at_24h.csv"), row.names=FALSE)

summary_df <- data.frame(
  category = c("early", "late", "sustained", "early_only", "late_only"),
  n = c(
    nrow(cor_early_ann),
    nrow(cor_late_ann),
    length(pairs_stable),
    length(pairs_early_only),
    length(pairs_late_only)
  )
)

write.csv(summary_df, file.path(output_dir,"correlation_summary_counts.csv"), row.names=FALSE)
summary_df

library(dplyr)

# Asegurarse que las tablas tengan columna "time" para saber de dónde vienen
cor_early_ann$time  <- "6h"
cor_late_ann$time   <- "24h"
cor_sust_ann$time   <- "both"   # por ejemplo, correlations en genes sostenidos

cor_all <- bind_rows(
  cor_early_ann,
  cor_late_ann,
  cor_sust_ann
)

# Hubs globales: cuántos mRNAs únicos por lncRNA (en cualquier tiempo)
lnc_hubs_global <- cor_all %>%
  group_by(lncRNA, lnc_gene_name, lnc_biotype, lnc_source) %>%
  summarise(
    n_mRNAs = n_distinct(mRNA),
    mean_rho = mean(rho),
    max_rho = max(abs(rho)),
    n_edges = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_mRNAs), desc(abs(mean_rho)))

# Top 50 hubs globales
lnc_hubs_top50 <- lnc_hubs_global %>% slice_head(n = 50)

write.csv(lnc_hubs_global,
          file.path(output_dir, "lncRNA_hubs_global.csv"),
          row.names = FALSE)
write.csv(lnc_hubs_top50,
          file.path(output_dir, "lncRNA_hubs_global_top50.csv"),
          row.names = FALSE)

# Hubs por tiempo (6h / 24h)
lnc_hubs_by_time <- cor_all %>%
  group_by(time, lncRNA, lnc_gene_name, lnc_biotype, lnc_source) %>%
  summarise(
    n_mRNAs = n_distinct(mRNA),
    mean_rho = mean(rho),
    max_rho = max(abs(rho)),
    n_edges = n(),
    .groups = "drop"
  ) %>%
  arrange(time, desc(n_mRNAs))

write.csv(lnc_hubs_by_time,
          file.path(output_dir, "lncRNA_hubs_by_time.csv"),
          row.names = FALSE)

# Asegúrate que cluster_df tenga biotype_group
# (si ya lo habíamos definido antes, esto sólo lo reitera por si acaso)
cluster_df <- cluster_df %>%
  mutate(
    biotype_group = dplyr::case_when(
      gene_biotype == "protein_coding"               ~ "mRNA",
      source == "novel"                              ~ "novel_lncRNA",
      grepl("lnc", tolower(gene_biotype))            ~ "lncRNA",
      TRUE                                           ~ "other"
    )
  )

# Unir correlaciones con info de cluster según lncRNA
cor_all_cluster <- cor_all %>%
  left_join(
    cluster_df %>% dplyr::select(gene_id, cluster, pattern, biotype_group),
    by = c("lncRNA" = "gene_id")
  ) %>%
  rename(
    lnc_cluster      = cluster,
    lnc_pattern      = pattern,
    lnc_biotype_grp  = biotype_group
  )

# Sólo lncRNAs (conocidas + novel)
cor_all_lnc <- cor_all_cluster %>%
  filter(lnc_biotype_grp %in% c("lncRNA", "novel_lncRNA"))

cluster_lnc_summary <- cor_all_lnc %>%
  group_by(lnc_cluster, lnc_pattern, lnc_biotype_grp) %>%
  summarise(
    n_lncRNAs = n_distinct(lncRNA),
    n_mRNAs   = n_distinct(mRNA),
    n_edges   = n(),
    mean_abs_rho = mean(abs(rho)),
    .groups = "drop"
  ) %>%
  arrange(lnc_cluster, desc(n_edges))

write.csv(cluster_lnc_summary,
          file.path(output_dir, "Mfuzz_clusters_lncRNA_regulation_summary.csv"),
          row.names = FALSE)

res6_full <- as.data.frame(res_NE_vs_Ctrl_6h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_small, by = "gene_id")

res24_full <- as.data.frame(res_NE_vs_Ctrl_24h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_small, by = "gene_id")

add_direction <- function(df, lfc_cut) {
  df %>%
    mutate(
      direction = dplyr::case_when(
        padj < 0.05 & log2FoldChange >=  lfc_cut  ~ "UP",
        padj < 0.05 & log2FoldChange <= -lfc_cut  ~ "DOWN",
        TRUE                                      ~ "NS"
      )
    )
}

# 6h mRNAs
res6_mrna <- res6_full %>%
  filter(gene_biotype == "protein_coding") %>%
  add_direction(lfc_cut = 1) %>%
  filter(direction != "NS")   # solo los significativos

write.table(res6_mrna,
            file = file.path(output_dir, "NE_vs_Ctrl_6h_mRNA_DE.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6h lncRNAs (conocidas + novel; cualquier cosa que tenga "lnc" en el biotipo o source == novel)
res6_lnc <- res6_full %>%
  filter(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel") %>%
  add_direction(lfc_cut = 0.5) %>%
  filter(direction != "NS")

write.table(res6_lnc,
            file = file.path(output_dir, "NE_vs_Ctrl_6h_lncRNA_DE.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 24h mRNAs
res24_mrna <- res24_full %>%
  filter(gene_biotype == "protein_coding") %>%
  add_direction(lfc_cut = 1) %>%
  filter(direction != "NS")

write.table(res24_mrna,
            file = file.path(output_dir, "NE_vs_Ctrl_24h_mRNA_DE.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 24h lncRNAs
res24_lnc <- res24_full %>%
  filter(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel") %>%
  add_direction(lfc_cut = 0.5) %>%
  filter(direction != "NS")

write.table(res24_lnc,
            file = file.path(output_dir, "NE_vs_Ctrl_24h_lncRNA_DE.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

###Volcano plots
library(dplyr)
library(ggplot2)

# Tablas DE con anotación
res6_full <- as.data.frame(res_NE_vs_Ctrl_6h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_small, by = "gene_id")

res24_full <- as.data.frame(res_NE_vs_Ctrl_24h) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_small, by = "gene_id")

add_direction <- function(df, lfc_cut) {
  df %>%
    mutate(
      direction = case_when(
        padj < 0.05 & log2FoldChange >=  lfc_cut  ~ "UP",
        padj < 0.05 & log2FoldChange <= -lfc_cut  ~ "DOWN",
        TRUE                                      ~ "NS"
      )
    )
}

# mRNAs (protein_coding)
all_mrna_lfc <- c(
  res6_full  %>% filter(gene_biotype == "protein_coding") %>% pull(log2FoldChange),
  res24_full %>% filter(gene_biotype == "protein_coding") %>% pull(log2FoldChange)
)
all_mrna_lfc <- all_mrna_lfc[!is.na(all_mrna_lfc)]
L_mrna <- max(abs(all_mrna_lfc))
L_mrna <- ceiling(L_mrna * 10) / 10   # redondear para que quede bonito

# lncRNAs (conocidas + novel)
all_lnc_lfc <- c(
  res6_full  %>% filter(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel") %>% pull(log2FoldChange),
  res24_full %>% filter(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel") %>% pull(log2FoldChange)
)
all_lnc_lfc <- all_lnc_lfc[!is.na(all_lnc_lfc)]
L_lnc <- max(abs(all_lnc_lfc))
L_lnc <- ceiling(L_lnc * 10) / 10

# 6h lncRNAs (conocidas + novel)
p_lnc_6h <- make_volcano(
  df = res6_full,
  biotype_filter = quote(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel"),
  lfc_cut = 0.5,
  xlim_range = L_lnc,
  title = "NE vs Ctrl 6h – lncRNAs (known + novel)",
  outfile = "Volcano_NE_vs_Ctrl_6h_lncRNA.png"
)

# 24h lncRNAs (conocidas + novel)
p_lnc_24h <- make_volcano(
  df = res24_full,
  biotype_filter = quote(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel"),
  lfc_cut = 0.5,
  xlim_range = L_lnc,
  title = "NE vs Ctrl 24h – lncRNAs (known + novel)",
  outfile = "Volcano_NE_vs_Ctrl_24h_lncRNA.png"
)

x_axis_limit <- 15

library(ggrepel)

make_volcano <- function(df,
                         biotype_filter,
                         lfc_cut,
                         title,
                         outfile,
                         x_axis_limit = 15,
                         top_labels = 10   # cuántos genes etiquetar
) {
  
  df_sub <- df %>%
    filter(!!biotype_filter) %>%
    mutate(
      padj_plot = ifelse(is.na(padj) | padj <= 0, NA, padj),
      negLog10Padj = -log10(padj_plot)
    ) %>%
    add_direction(lfc_cut = lfc_cut) %>%
    filter(!is.na(negLog10Padj))
  
  if (nrow(df_sub) == 0) {
    message("No hay genes para plotear en: ", title)
    return(NULL)
  }
  
  # Selección inteligente de genes a etiquetar:
  df_sub <- df_sub %>%
    mutate(absLFC = abs(log2FoldChange))
  
  label_df <- bind_rows(
    df_sub %>% filter(direction == "UP")   %>% arrange(desc(log2FoldChange)) %>% head(top_labels),
    df_sub %>% filter(direction == "DOWN") %>% arrange(log2FoldChange)       %>% head(top_labels),
    df_sub %>% arrange(desc(negLog10Padj)) %>% head(top_labels)
  ) %>% distinct(gene_id, .keep_all = TRUE)
  
  # Volcano plot
  p <- ggplot(df_sub, aes(
    x = log2FoldChange,
    y = negLog10Padj,
    color = direction
  )) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = c(
      "DOWN" = "royalblue",
      "NS"   = "grey70",
      "UP"   = "firebrick"
    )) +
    xlim(-x_axis_limit, x_axis_limit) +
    labs(
      title = title,
      x = "log2 Fold Change (NE vs Ctrl)",
      y = expression(-log[10]("padj"))
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.position = "right",
      text = element_text(size = 12)
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = gene_name),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0.1,
      box.padding = 0.3,
      point.padding = 0.2,
      seed = 123,
      fontface = "bold"
    )
  
  # Guardar PNG en alta resolución
  png(file.path(output_dir, outfile),
      width = 2200, height = 2000, res = 300)
  print(p)
  dev.off()
  
  return(p)
}


# 6h mRNA
p_mrna_6h <- make_volcano(
  df = res6_full,
  biotype_filter = quote(gene_biotype == "protein_coding"),
  lfc_cut = 1,
  title = "NE vs Ctrl 6h – mRNAs",
  outfile = "Volcano_NE_vs_Ctrl_6h_mRNA.png",
  x_axis_limit = 15,
  top_labels = 15
)


# 24h mRNA
p_mrna_24h <- make_volcano(
  df = res24_full,
  biotype_filter = quote(gene_biotype == "protein_coding"),
  lfc_cut = 1,
  title = "NE vs Ctrl 24h – mRNAs",
  outfile = "Volcano_NE_vs_Ctrl_24h_mRNA.png",
  x_axis_limit = 15,
  top_labels = 15
)


# 6h lncRNAs (conocidas + novel)
p_lnc_6h <- make_volcano(
  df = res6_full,
  biotype_filter = quote(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel"),
  lfc_cut = 0.5,
  title = "NE vs Ctrl 6h – lncRNAs (known + novel)",
  outfile = "Volcano_NE_vs_Ctrl_6h_lncRNA.png",
  x_axis_limit = 15,
  top_labels = 15
)


# 24h lncRNAs (conocidas + novel)
p_lnc_24h <- make_volcano(
  df = res24_full,
  biotype_filter = quote(grepl("lnc", gene_biotype, ignore.case = TRUE) | source == "novel"),
  lfc_cut = 0.5,
  title = "NE vs Ctrl 24h – lncRNAs (known + novel)",
  outfile = "Volcano_NE_vs_Ctrl_24h_lncRNA.png",
  x_axis_limit = 15,
  top_labels = 15
)

## ==== 0. Librerías ====
library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(org.Rn.eg.db)   # Para rata; cambia a org.Mm.eg.db u org.Hs.eg.db si corresponde

## ==== 1. Leer archivo de GSEA ====
# Cambia el nombre del archivo por el tuyo
gsea_file <- "KEGG_6h_vs_24h_Shared.txt"

gsea_res <- read.delim(
  gsea_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Verifica que exista la columna core_enrichment
stopifnot("core_enrichment" %in% colnames(gsea_res))

## ==== 2. Obtener todos los ENTREZID únicos del core_enrichment ====
all_entrez <- gsea_res$core_enrichment %>%
  na.omit() %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

length(all_entrez)
head(all_entrez)

## ==== 3. Mapear ENTREZID -> SYMBOL (gene symbol) ====
# Para rata: OrgDb = org.Rn.eg.db
id_map <- bitr(
  all_entrez,
  fromType = "ENTREZID",
  toType   = "SYMBOL",
  OrgDb    = org.Rn.eg.db
)

# Crear un vector nombrado: nombres = ENTREZID, valores = SYMBOL
entrez2symbol <- id_map$SYMBOL
names(entrez2symbol) <- id_map$ENTREZID

# Si hay IDs sin symbol, los dejamos con el mismo ID como fallback
missing_ids <- setdiff(all_entrez, names(entrez2symbol))
if (length(missing_ids) > 0) {
  entrez2symbol[missing_ids] <- missing_ids
}

## ==== 4. ARCHIVO 1: misma tabla GSEA + core_enrichment en SYMBOL separado por "/" ====

gsea_res_symbols <- gsea_res %>%
  mutate(
    core_enrichment_symbol = sapply(
      strsplit(core_enrichment, "/"),
      function(x) paste(entrez2symbol[x], collapse = "/")
    )
  )

# Guardar archivo 1
write.table(
  gsea_res_symbols,
  file = "KEGG_6h_vs_24h_Shared_with_core_symbols.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ==== 5. ARCHIVO 2: formato largo, un gen por fila con proceso ====
# Aquí dejamos por fila: ID (KEGG), Description (proceso),
# ENTREZID, gene_symbol (SYMBOL) y lo que más quieras conservar (ej. NES, pvalue)

library(dplyr)
library(tidyr)
library(readr)

##opción de llamado
gsea_res <- read_tsv(
  "KEGG_6h_vs_24h_Shared_with_core_symbols.txt",
  show_col_types = FALSE
)

# === 2. Generar formato largo ===
core_long <- gsea_res %>%
  dplyr::select(ID, Description, NES_6h, NES_24h, pvalue, p.adjust_6h, p.adjust_24h, core_enrichment_symbol) %>%
  dplyr::mutate(
    core_enrichment_symbol = as.character(.data[["core_enrichment_symbol"]])
  ) %>%
  dplyr::filter(
    !is.na(.data[["core_enrichment_symbol"]]),
    .data[["core_enrichment_symbol"]] != ""
  ) %>%
  tidyr::separate_rows(core_enrichment_symbol, sep = "/") %>%
  dplyr::rename(gene_symbol = core_enrichment_symbol) %>%
  dplyr::arrange(ID, Description, gene_symbol)

# === 3. Guardar archivo ===
write.table(
  core_long,
  file = "KEGG_6h_vs_24h_Shared_long_format_SYMBOLS.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

##Para los no shareds
core_long <- gsea_res %>%
  dplyr::select(ID, Description, NES, pvalue, p.adjust, core_enrichment_symbol) %>%
  # Forzamos a character la columna, usando .data para evitar problemas de búsqueda
  dplyr::mutate(
    core_enrichment_symbol = as.character(.data[["core_enrichment_symbol"]])
  ) %>%
  dplyr::filter(
    !is.na(.data[["core_enrichment_symbol"]]),
    .data[["core_enrichment_symbol"]] != ""
  ) %>%
  # separate_rows sí puede usar el nombre de la columna directamente
  tidyr::separate_rows(core_enrichment_symbol, sep = "/") %>%
  dplyr::rename(gene_symbol = core_enrichment_symbol) %>%
  dplyr::arrange(ID, Description, gene_symbol)

# Guardar el archivo
write.table(
  core_long,
  file = "KEGG_6h_vs_24h_Only24h_long_format_SYMBOLS.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


###versión con las cuentas crudas
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggvenn)
library(UpSetR)

# 1. Paths, metadata y archivos
base_dir   <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/salmon_quants_NE624/"
gtf_file   <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis/merged_CNE_6_24_HISAT2.annotated.gtf"
output_dir <- "DESeq2_Multifactorial_Results_fixed2"
if (!dir.exists(output_dir)) dir.create(output_dir)

samples <- c("Ctrl6_rep1","Ctrl6_rep2","Ctrl6_rep4",
             "NE6_rep1","NE6_rep2","NE6_rep4",
             "Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
             "NE24_rep1","NE24_rep2","NE24_rep4")

sampleTable <- data.frame(
  sample    = samples,
  treatment = factor(rep(c("Control","NE","Control","NE"), each = 3)),
  time      = factor(rep(c("6h","24h"), each = 6))
)
rownames(sampleTable) <- samples
sampleTable$treatment <- relevel(sampleTable$treatment, "Control")
sampleTable$time      <- relevel(sampleTable$time, "6h")
sampleTable$group     <- paste(sampleTable$treatment, sampleTable$time, sep = "_")

files <- file.path(base_dir, paste0(samples, "_quant"), "quant.sf")
names(files) <- samples
if (!all(file.exists(files))) stop("Faltan archivos quant.sf")

# 2. tximport + DESeq multifactorial
tx2gene_fixed <- read_tsv("tx2gene_fixed.tsv",
                          col_names = c("transcript_id","gene_id"))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene_fixed)

dds <- DESeqDataSetFromTximport(txi,
                                colData = sampleTable,
                                design = ~ time + treatment + time:treatment)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# 3. QC: PCA y distancias
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA – NE vs Control por tiempo")

ggsave(file.path(output_dir, "PCA_NE_Control_6h_24h.png"),
       p_pca, width = 6, height = 5, dpi = 300)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, col = colors, main = "Sample distance heatmap",
         filename = file.path(output_dir, "SampleDistanceHeatmap.png"),
         width = 6, height = 5)

# 4. DE por tiempo: NE vs Control a 6h y 24h
dds_6  <- dds[, dds$time == "6h"]
dds_6$time <- droplevels(dds_6$time)
design(dds_6) <- ~ treatment
dds_6 <- DESeq(dds_6)

dds_24 <- dds[, dds$time == "24h"]
dds_24$time <- droplevels(dds_24$time)
design(dds_24) <- ~ treatment
dds_24 <- DESeq(dds_24)

res_NE_vs_Ctrl_6h  <- results(dds_6,  contrast = c("treatment","NE","Control"))
res_NE_vs_Ctrl_24h <- results(dds_24, contrast = c("treatment","NE","Control"))

# (Opcional) tablas básicas de resultados
write.csv(as.data.frame(res_NE_vs_Ctrl_6h),
          file.path(output_dir, "res_NE_vs_Ctrl_6h_basic.csv"))
write.csv(as.data.frame(res_NE_vs_Ctrl_24h),
          file.path(output_dir, "res_NE_vs_Ctrl_24h_basic.csv"))

# 5. Criterio DE y listas (padj < 0.05 & |log2FC| >= 1)
alpha   <- 0.05
lfc_cut <- 1

res_6_df_base  <- as.data.frame(res_NE_vs_Ctrl_6h)
res_24_df_base <- as.data.frame(res_NE_vs_Ctrl_24h)

sig_6h <- res_6_df_base[!is.na(res_6_df_base$padj) &
                          res_6_df_base$padj < alpha &
                          abs(res_6_df_base$log2FoldChange) >= lfc_cut, ]

sig_24h <- res_24_df_base[!is.na(res_24_df_base$padj) &
                            res_24_df_base$padj < alpha &
                            abs(res_24_df_base$log2FoldChange) >= lfc_cut, ]

genes_6h  <- rownames(sig_6h)
genes_24h <- rownames(sig_24h)

cat("Genes DE (padj<0.05 & |log2FC|>=1) a 6h:",  length(genes_6h),  "\n")
cat("Genes DE (padj<0.05 & |log2FC|>=1) a 24h:", length(genes_24h), "\n")

intersect_genes <- intersect(genes_6h, genes_24h)
only_6h         <- setdiff(genes_6h,  genes_24h)
only_24h        <- setdiff(genes_24h, genes_6h)

cat("DE en ambos tiempos (sostenidos):", length(intersect_genes), "\n")
cat("DE solo en 6h (early):",          length(only_6h),         "\n")
cat("DE solo en 24h (late):",          length(only_24h),        "\n")

write.table(genes_6h,
            file = file.path(output_dir, "genes_DE_6h_all.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_24h,
            file = file.path(output_dir, "genes_DE_24h_all.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intersect_genes,
            file = file.path(output_dir, "genes_DE_common_6h_24h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(only_6h,
            file = file.path(output_dir, "genes_DE_only_6h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(only_24h,
            file = file.path(output_dir, "genes_DE_only_24h.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 6. Venn y UpSet (DE global)
venn_list <- list(
  DE_6h  = genes_6h,
  DE_24h = genes_24h
)

p_venn <- ggvenn(venn_list,
                 fill_alpha = 0.6,
                 stroke_size = 0.7,
                 set_name_size = 5)
ggsave(file.path(output_dir, "Venn_DE_6h_24h.png"),
       plot = p_venn, width = 5, height = 5, dpi = 300)

upset_data <- fromList(venn_list)
png(file.path(output_dir, "UpSet_DE_6h_24h.png"),
    width = 2000, height = 1500, res = 300)
upset(upset_data,
      nsets = 2,
      nintersects = 3,
      order.by = "freq",
      sets = c("DE_6h","DE_24h"),
      mainbar.y.label = "Number of genes",
      sets.x.label = "DEG per time")
dev.off()

# 7. UP/DOWN con mismo criterio
up_6h   <- rownames(sig_6h)[sig_6h$log2FoldChange >=  lfc_cut]
down_6h <- rownames(sig_6h)[sig_6h$log2FoldChange <= -lfc_cut]

up_24h   <- rownames(sig_24h)[sig_24h$log2FoldChange >=  lfc_cut]
down_24h <- rownames(sig_24h)[sig_24h$log2FoldChange <= -lfc_cut]

cat("6h: UP =", length(up_6h), "DOWN =", length(down_6h), "\n")
cat("24h: UP =", length(up_24h), "DOWN =", length(down_24h), "\n")

updown_list <- list(
  UP_6h    = up_6h,
  DOWN_6h  = down_6h,
  UP_24h   = up_24h,
  DOWN_24h = down_24h
)

updown_data <- fromList(updown_list)
png(file.path(output_dir, "UpSet_UP_DOWN_6h_24h.png"),
    width = 2500, height = 1800, res = 300)
upset(updown_data,
      nsets = 4,
      nintersects = 10,
      order.by = "freq",
      mainbar.y.label = "Number of genes (UP/DOWN)",
      sets.x.label = "UP/DOWN per time")
dev.off()

# 8. Anotaciones curadas
ann_all       <- read.csv("hisat2_6_24_meta_ann.csv", stringsAsFactors = FALSE)
ann_novel_raw <- read.csv("hisat2_6_24_meta_novel_lncRNA.csv", stringsAsFactors = FALSE)

ann_all_clean <- ann_all %>%
  transmute(
    gene_id            = ensembl_gene_id,
    transcript_id      = ensembl_transcript_id,
    gene_name          = external_gene_name,
    gene_biotype       = gene_biotype,
    transcript_biotype = transcript_biotype,
    source             = "annotated"
  )

ann_novel_clean <- ann_novel_raw %>%
  filter(type == "transcript") %>%
  transmute(
    gene_id            = gene_id,
    transcript_id      = transcript_id,
    gene_name          = qry_gene_id,
    gene_biotype       = pred_biotype,
    transcript_biotype = pred_biotype,
    source             = "novel"
  )

ann_merged <- bind_rows(ann_all_clean, ann_novel_clean) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Genes en anotación unificada:", nrow(ann_merged), "\n")

# 9. Anotar listas DE (early, late, sostenidos)
df_6h     <- data.frame(gene_id = genes_6h)
df_24h    <- data.frame(gene_id = genes_24h)
df_common <- data.frame(gene_id = intersect_genes)
df_only6  <- data.frame(gene_id = only_6h)
df_only24 <- data.frame(gene_id = only_24h)

df_6h_ann     <- df_6h     %>% left_join(ann_merged, by = "gene_id")
df_24h_ann    <- df_24h    %>% left_join(ann_merged, by = "gene_id")
df_common_ann <- df_common %>% left_join(ann_merged, by = "gene_id")
df_only6_ann  <- df_only6  %>% left_join(ann_merged, by = "gene_id")
df_only24_ann <- df_only24 %>% left_join(ann_merged, by = "gene_id")

write.csv(df_6h_ann,     file.path(output_dir, "DE_6h_annotated.csv"),     row.names = FALSE)
write.csv(df_24h_ann,    file.path(output_dir, "DE_24h_annotated.csv"),    row.names = FALSE)
write.csv(df_common_ann, file.path(output_dir, "DE_common_annotated.csv"), row.names = FALSE)
write.csv(df_only6_ann,  file.path(output_dir, "DE_only6_annotated.csv"),  row.names = FALSE)
write.csv(df_only24_ann, file.path(output_dir, "DE_only24_annotated.csv"), row.names = FALSE)

cat("Sin anotación en 6h:",  sum(is.na(df_6h_ann$gene_biotype)),  "\n")
cat("Sin anotación en 24h:", sum(is.na(df_24h_ann$gene_biotype)), "\n")

# 10. Separar por biotipo (mRNA, lncRNA conocidas, lncRNA novel)
ann_merged$gene_biotype <- tolower(ann_merged$gene_biotype)

is_mrna <- function(df) {
  df %>% filter(gene_biotype == "protein_coding")
}
is_lncRNA <- function(df) {
  df %>% filter(grepl("lnc", gene_biotype))
}
is_novel <- function(df) {
  df %>% filter(source == "novel")
}

only6_mrna   <- is_mrna(df_only6_ann)
only6_lnc    <- is_lncRNA(df_only6_ann)
only6_novel  <- is_novel(df_only6_ann)

only24_mrna  <- is_mrna(df_only24_ann)
only24_lnc   <- is_lncRNA(df_only24_ann)
only24_novel <- is_novel(df_only24_ann)

common_mrna  <- is_mrna(df_common_ann)
common_lnc   <- is_lncRNA(df_common_ann)
common_novel <- is_novel(df_common_ann)

cat("=== Early (6h-only) ===\n")
cat("mRNAs:", nrow(only6_mrna),
    "lncRNAs:", nrow(only6_lnc),
    "novel lncRNAs:", nrow(only6_novel), "\n")

cat("=== Late (24h-only) ===\n")
cat("mRNAs:", nrow(only24_mrna),
    "lncRNAs:", nrow(only24_lnc),
    "novel lncRNAs:", nrow(only24_novel), "\n")

cat("=== Sostenidos (6h & 24h) ===\n")
cat("mRNAs:", nrow(common_mrna),
    "lncRNAs:", nrow(common_lnc),
    "novel lncRNAs:", nrow(common_novel), "\n")

write.csv(only6_mrna,   file.path(output_dir,"early_mRNA.csv"),             row.names = FALSE)
write.csv(only6_lnc,    file.path(output_dir,"early_lncRNA.csv"),           row.names = FALSE)
write.csv(only6_novel,  file.path(output_dir,"early_novel_lncRNA.csv"),     row.names = FALSE)

write.csv(only24_mrna,  file.path(output_dir,"late_mRNA.csv"),              row.names = FALSE)
write.csv(only24_lnc,   file.path(output_dir,"late_lncRNA.csv"),            row.names = FALSE)
write.csv(only24_novel, file.path(output_dir,"late_novel_lncRNA.csv"),      row.names = FALSE)

write.csv(common_mrna,  file.path(output_dir,"sustained_mRNA.csv"),         row.names = FALSE)
write.csv(common_lnc,   file.path(output_dir,"sustained_lncRNA.csv"),       row.names = FALSE)
write.csv(common_novel, file.path(output_dir,"sustained_novel_lncRNA.csv"), row.names = FALSE)

# 11. Tablas finales 6h y 24h con anotación + status + counts
# 6h
raw_counts_6  <- counts(dds_6, normalized = FALSE)
norm_counts_6 <- counts(dds_6, normalized = TRUE)

raw_6_df  <- as.data.frame(raw_counts_6)
norm_6_df <- as.data.frame(norm_counts_6)

colnames(raw_6_df)  <- paste0(colnames(raw_6_df), "_raw")
colnames(norm_6_df) <- paste0(colnames(norm_6_df), "_norm")

res_6_df <- as.data.frame(res_NE_vs_Ctrl_6h)
res_6_df$gene_id <- rownames(res_6_df)

res_6_df <- res_6_df %>%
  left_join(ann_merged, by = "gene_id") %>%
  mutate(
    status = case_when(
      !is.na(padj) & padj < alpha & log2FoldChange >=  lfc_cut ~ "UP",
      !is.na(padj) & padj < alpha & log2FoldChange <= -lfc_cut ~ "DOWN",
      TRUE ~ "NS"
    )
  ) %>%
  select(
    gene_id, gene_name, gene_biotype, transcript_biotype, source,
    status,
    baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
    everything()
  )

res_6_full <- cbind(
  res_6_df,
  raw_6_df[res_6_df$gene_id, ],
  norm_6_df[res_6_df$gene_id, ]
)

write.csv(res_6_full,
          file.path(output_dir, "res_NE_vs_Ctrl_6h_annot_counts_status.csv"),
          row.names = FALSE)

write.table(res_6_full,
            file = file.path(output_dir,
                             "res_NE_vs_Ctrl_6h_annot_counts_status.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 24h
raw_counts_24  <- counts(dds_24, normalized = FALSE)
norm_counts_24 <- counts(dds_24, normalized = TRUE)

raw_24_df  <- as.data.frame(raw_counts_24)
norm_24_df <- as.data.frame(norm_counts_24)

colnames(raw_24_df)  <- paste0(colnames(raw_24_df), "_raw")
colnames(norm_24_df) <- paste0(colnames(norm_24_df), "_norm")

res_24_df <- as.data.frame(res_NE_vs_Ctrl_24h)
res_24_df$gene_id <- rownames(res_24_df)

res_24_df <- res_24_df %>%
  left_join(ann_merged, by = "gene_id") %>%
  mutate(
    status = case_when(
      !is.na(padj) & padj < alpha & log2FoldChange >=  lfc_cut ~ "UP",
      !is.na(padj) & padj < alpha & log2FoldChange <= -lfc_cut ~ "DOWN",
      TRUE ~ "NS"
    )
  ) %>%
  select(
    gene_id, gene_name, gene_biotype, transcript_biotype, source,
    status,
    baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
    everything()
  )

res_24_full <- cbind(
  res_24_df,
  raw_24_df[res_24_df$gene_id, ],
  norm_24_df[res_24_df$gene_id, ]
)

write.csv(res_24_full,
          file.path(output_dir, "res_NE_vs_Ctrl_24h_annot_counts_status.csv"),
          row.names = FALSE)

write.table(res_24_full,
            file = file.path(output_dir,
                             "res_NE_vs_Ctrl_24h_annot_counts_status.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
