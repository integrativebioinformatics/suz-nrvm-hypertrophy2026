# ==============================
# Análisis de miRNAs NE 6h / 24h
# ==============================

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Opcional: trabajar en la carpeta donde están los archivos
# setwd("C:/Users/.../TU/CARPETA")

output_dir <- "miRNA_DE_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 1) Cargar matrices de cuentas (6h y 24h) --------------------

file_6h  <- "miR.Counts_C6_NE6_3v3_new.csv"
file_24h <- "miR.Counts_C24_NE24_3v3_new.csv"

counts_6 <- read.csv(file_6h,  header = TRUE, row.names = 1, check.names = FALSE)
counts_24 <- read.csv(file_24h, header = TRUE, row.names = 1, check.names = FALSE)

# Revisar nombres originales de columnas
colnames(counts_6)
colnames(counts_24)

# 2) Renombrar columnas a nombres coherentes con el resto del análisis ----
# 6 h:
# H1_mi  -> Ctrl6_rep1
# H2_mi  -> NE6_rep1
# SQ1_mi -> Ctrl6_rep2
# SQ2_mi -> NE6_rep2
# SQ13_mi -> Ctrl6_rep4
# SQ14_mi -> NE6_rep4

colnames(counts_6) <- c("Ctrl6_rep1", "NE6_rep1",
                        "Ctrl6_rep2", "NE6_rep2",
                        "Ctrl6_rep4", "NE6_rep4")

# 24 h:
# H3_mi  -> Ctrl24_rep1
# H4_mi  -> NE24_rep1
# SQ3_mi -> Ctrl24_rep2
# SQ4_mi -> NE24_rep2
# SQ15_mi -> Ctrl24_rep4
# SQ16_mi -> NE24_rep4

colnames(counts_24) <- c("Ctrl24_rep1", "NE24_rep1",
                         "Ctrl24_rep2", "NE24_rep2",
                         "Ctrl24_rep4", "NE24_rep4")

# 3) Alinear miRNAs entre matrices y fusionar -----------------

common_ids <- intersect(rownames(counts_6), rownames(counts_24))
counts_6   <- counts_6[common_ids, ]
counts_24  <- counts_24[common_ids, ]

mirna_counts <- cbind(counts_6, counts_24)

# Orden de columnas (muestras) consistente
samples_mirna <- c("Ctrl6_rep1", "Ctrl6_rep2", "Ctrl6_rep4",
                   "NE6_rep1",   "NE6_rep2",   "NE6_rep4",
                   "Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
                   "NE24_rep1",  "NE24_rep2",  "NE24_rep4")

mirna_counts <- mirna_counts[, samples_mirna]

# 4) Metadatos de las muestras (colData) ----------------------

sampleTable_mirna <- data.frame(
  sample    = samples_mirna,
  treatment = factor(rep(c("Control", "NE", "Control", "NE"), each = 3)),
  time      = factor(rep(c("6h", "24h"), each = 6))
)
rownames(sampleTable_mirna) <- samples_mirna

sampleTable_mirna$treatment <- relevel(sampleTable_mirna$treatment, ref = "Control")
sampleTable_mirna$time      <- relevel(sampleTable_mirna$time,      ref = "6h")

# 5) Crear DESeqDataSet y correr DESeq2 -----------------------

dds_mirna <- DESeqDataSetFromMatrix(
  countData = mirna_counts,
  colData   = sampleTable_mirna,
  design    = ~ time + treatment + time:treatment
)

# Prefiltrado: eliminar miRNAs con muy pocas cuentas
keep_mi <- rowSums(counts(dds_mirna)) >= 10
dds_mirna <- dds_mirna[keep_mi, ]

dds_mirna <- DESeq(dds_mirna)
resultsNames(dds_mirna)
# Esperado: "Intercept", "time_24h_vs_6h", "treatment_NE_vs_Control", "time24h.treatmentNE"

colData(dds_mirna)$group <- with(colData(dds_mirna),
                                 paste(treatment, time, sep = "_"))

# 6) Transformación VST y QC ----------------------------------

vsd_mirna <- varianceStabilizingTransformation(dds_mirna, blind = FALSE)

pcaData_mi <- plotPCA(vsd_mirna, intgroup = "group", returnData = TRUE)
percentVar_mi <- round(100 * attr(pcaData_mi, "percentVar"))

p_pca_mi <- ggplot(pcaData_mi,
                   aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar_mi[1], "%")) +
  ylab(paste0("PC2: ", percentVar_mi[2], "%")) +
  theme_bw() +
  ggtitle("PCA – miRNAs NE vs Control (6h / 24h)")

ggsave(file.path(output_dir, "PCA_miRNA_6h_24h.png"),
       p_pca_mi, width = 7, height = 6, dpi = 300)

# Heatmap de distancias
sampleDists_mi <- dist(t(assay(vsd_mirna)))
sampleDistMatrix_mi <- as.matrix(sampleDists_mi)
colors_mi <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(file.path(output_dir, "Heatmap_distancias_miRNA.png"),
    width = 2000, height = 1800, res = 250)
pheatmap(sampleDistMatrix_mi,
         clustering_distance_rows = sampleDists_mi,
         clustering_distance_cols = sampleDists_mi,
         col = colors_mi,
         main = "Distancias entre muestras – miRNAs")
dev.off()

# 7) Contrastes: NE vs Control a 6h y a 24h -------------------

# NE vs Ctrl 6h = coeficiente principal de tratamiento
res_mi_6h <- results(dds_mirna, name = "treatment_NE_vs_Control", alpha = 0.05)

# NE vs Ctrl 24h = tratamiento + interacción
res_mi_24h <- results(
  dds_mirna,
  contrast = list(c("treatment_NE_vs_Control", "time24h.treatmentNE")),
  alpha = 0.05
)

res_mi_6h  <- res_mi_6h[order(res_mi_6h$padj), ]
res_mi_24h <- res_mi_24h[order(res_mi_24h$padj), ]

write.csv(as.data.frame(res_mi_6h),
          file.path(output_dir, "miRNA_res_NE_vs_Ctrl_6h.csv"))
write.csv(as.data.frame(res_mi_24h),
          file.path(output_dir, "miRNA_res_NE_vs_Ctrl_24h.csv"))

# 8) Definir dirección (UP/DOWN) y preparar listas para miRWalk -----

add_direction_mi <- function(df, lfc_cut = 0.5) {
  df %>%
    mutate(
      direction = case_when(
        padj < 0.05 & log2FoldChange >=  lfc_cut  ~ "UP",
        padj < 0.05 & log2FoldChange <= -lfc_cut  ~ "DOWN",
        TRUE                                      ~ "NS"
      )
    )
}

res_mi_6h_tbl <- as.data.frame(res_mi_6h) %>%
  rownames_to_column("miRNA_id") %>%
  add_direction_mi(lfc_cut = 0.5)

res_mi_24h_tbl <- as.data.frame(res_mi_24h) %>%
  rownames_to_column("miRNA_id") %>%
  add_direction_mi(lfc_cut = 0.5)

# Filtrar sólo miRNAs DE (UP o DOWN)
mi_6h_DE <- res_mi_6h_tbl %>% filter(direction != "NS")
mi_24h_DE <- res_mi_24h_tbl %>% filter(direction != "NS")

write.table(mi_6h_DE,
            file = file.path(output_dir, "miRNA_DE_NE_vs_Ctrl_6h_for_miRWalk.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mi_24h_DE,
            file = file.path(output_dir, "miRNA_DE_NE_vs_Ctrl_24h_for_miRWalk.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 9) Volcano plots para miRNAs -------------------------------

make_volcano_mi <- function(df,
                            title,
                            outfile,
                            x_axis_limit = 6,
                            top_labels = 15) {
  
  df_sub <- df %>%
    mutate(
      padj_plot = ifelse(is.na(padj) | padj <= 0, NA, padj),
      negLog10Padj = -log10(padj_plot)
    ) %>%
    filter(!is.na(negLog10Padj))
  
  if (nrow(df_sub) == 0) {
    message("No hay miRNAs para plotear en: ", title)
    return(NULL)
  }
  
  df_sub <- df_sub %>% mutate(absLFC = abs(log2FoldChange))
  
  label_df <- bind_rows(
    df_sub %>% filter(direction == "UP")   %>% arrange(desc(log2FoldChange)) %>% head(top_labels),
    df_sub %>% filter(direction == "DOWN") %>% arrange(log2FoldChange)       %>% head(top_labels),
    df_sub %>% arrange(desc(negLog10Padj)) %>% head(top_labels)
  ) %>%
    distinct(miRNA_id, .keep_all = TRUE)
  
  p <- ggplot(df_sub,
              aes(x = log2FoldChange,
                  y = negLog10Padj,
                  color = direction)) +
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
      aes(label = miRNA_id),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0.1,
      box.padding = 0.3,
      point.padding = 0.2,
      seed = 123,
      fontface = "bold"
    )
  
  png(file.path(output_dir, outfile),
      width = 2200, height = 2000, res = 300)
  print(p)
  dev.off()
  
  return(p)
}

# Volcano 6h
p_mi_6h <- make_volcano_mi(
  df = res_mi_6h_tbl,
  title   = "NE vs Ctrl 6h – miRNAs",
  outfile = "Volcano_miRNA_NE_vs_Ctrl_6h.png",
  x_axis_limit = 6,
  top_labels   = 15
)

# Volcano 24h
p_mi_24h <- make_volcano_mi(
  df = res_mi_24h_tbl,
  title   = "NE vs Ctrl 24h – miRNAs",
  outfile = "Volcano_miRNA_NE_vs_Ctrl_24h.png",
  x_axis_limit = 6,
  top_labels   = 15
)

message("Análisis de miRNAs completado (DESeq2, QC, volcanos y listas para miRWalk).")


##Veamos con efecto batch
samples_mirna <- c("Ctrl6_rep1", "Ctrl6_rep2", "Ctrl6_rep4",
                   "NE6_rep1",   "NE6_rep2",   "NE6_rep4",
                   "Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
                   "NE24_rep1",  "NE24_rep2",  "NE24_rep4")

mirna_counts <- mirna_counts[, samples_mirna]

sampleTable_mirna <- data.frame(
  sample    = samples_mirna,
  treatment = factor(rep(c("Control", "NE", "Control", "NE"), each = 3)),
  time      = factor(rep(c("6h", "24h"), each = 6)),
  batch     = factor(c(
    # 6h:   Ctrl6_rep1, Ctrl6_rep2, Ctrl6_rep4, NE6_rep1, NE6_rep2, NE6_rep4
    "H","SQ","SQ", "H","SQ","SQ",
    # 24h:  Ctrl24_rep1, Ctrl24_rep2, Ctrl24_rep4, NE24_rep1, NE24_rep2, NE24_rep4
    "H","SQ","SQ", "H","SQ","SQ"
  ))
)
rownames(sampleTable_mirna) <- samples_mirna

sampleTable_mirna$treatment <- relevel(sampleTable_mirna$treatment, ref = "Control")
sampleTable_mirna$time      <- relevel(sampleTable_mirna$time,      ref = "6h")
sampleTable_mirna$batch     <- relevel(sampleTable_mirna$batch,     ref = "H")

dds_mirna <- DESeqDataSetFromMatrix(
  countData = mirna_counts,
  colData   = sampleTable_mirna,
  design    = ~ batch + time + treatment + time:treatment
)

keep_mi   <- rowSums(counts(dds_mirna)) >= 10
dds_mirna <- dds_mirna[keep_mi, ]

dds_mirna <- DESeq(dds_mirna)
resultsNames(dds_mirna)

vsd_mirna <- varianceStabilizingTransformation(dds_mirna, blind = FALSE)

colData(vsd_mirna)$group <- with(colData(vsd_mirna),
                                 paste(treatment, time, sep = "_"))

pcaData_mi <- plotPCA(vsd_mirna, intgroup = "group", returnData = TRUE)
percentVar_mi <- round(100 * attr(pcaData_mi, "percentVar"))

p_pca_mi <- ggplot(pcaData_mi,
                   aes(PC1, PC2, color = group, shape = colData(vsd_mirna)$batch[name])) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar_mi[1], "%")) +
  ylab(paste0("PC2: ", percentVar_mi[2], "%")) +
  theme_bw() +
  ggtitle("PCA – miRNAs NE vs Control (6h / 24h, corregido por batch)") +
  labs(shape = "batch")

ggsave(file.path(output_dir, "PCA_miRNA_6h_24h_batchCorrected.png"),
       p_pca_mi, width = 7, height = 6, dpi = 300)


####Nuevo
# ==============================
# QC miRNAs NE 6h / 24h (PCA + heatmap)
# ==============================

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# 1) Archivos y carpetas -------------------------------------

file_6h  <- "miR.Counts_C6_NE6_3v3_new.csv"
file_24h <- "miR.Counts_C24_NE24_3v3_new.csv"

output_dir_mi <- "miRNA_QC_multifactorial"
if (!dir.exists(output_dir_mi)) dir.create(output_dir_mi)

# 2) Cargar matrices de cuentas ------------------------------

counts_6  <- read.csv(file_6h,  header = TRUE, row.names = 1, check.names = FALSE)
counts_24 <- read.csv(file_24h, header = TRUE, row.names = 1, check.names = FALSE)

# Revisa por si acaso:
cat("Columnas 6h originales:\n");  print(colnames(counts_6))
cat("Columnas 24h originales:\n"); print(colnames(counts_24))

# IMPORTANTE: los .csv están como:
#  columnas 1-3 = Ctrl,  columnas 4-6 = NE

colnames(counts_6) <- c("Ctrl6_rep1", "Ctrl6_rep2", "Ctrl6_rep4",
                        "NE6_rep1",   "NE6_rep2",   "NE6_rep4")

colnames(counts_24) <- c("Ctrl24_rep1", "Ctrl24_rep2", "Ctrl24_rep4",
                         "NE24_rep1",   "NE24_rep2",   "NE24_rep4")

# 3) Alinear miRNAs entre matrices y fusionar ----------------

common_ids <- intersect(rownames(counts_6), rownames(counts_24))
counts_6   <- counts_6[common_ids, ]
counts_24  <- counts_24[common_ids, ]

mirna_counts <- cbind(counts_6, counts_24)

# Orden coherente de muestras
samples <- c("Ctrl6_rep1","Ctrl6_rep2","Ctrl6_rep4",
             "NE6_rep1","NE6_rep2","NE6_rep4",
             "Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
             "NE24_rep1","NE24_rep2","NE24_rep4")

mirna_counts <- mirna_counts[, samples]

# 4) Metadata (colData) --------------------------------------

sampleTable_mi <- data.frame(
  sample    = samples,
  treatment = c(rep("Control", 3),  # Ctrl6
                rep("NE",      3),  # NE6
                rep("Control", 3),  # Ctrl24
                rep("NE",      3)), # NE24
  time      = c(rep("6h", 6),
                rep("24h", 6))
)

sampleTable_mi$treatment <- factor(sampleTable_mi$treatment,
                                   levels = c("Control","NE"))
sampleTable_mi$time <- factor(sampleTable_mi$time,
                              levels = c("6h","24h"))

rownames(sampleTable_mi) <- sampleTable_mi$sample
sampleTable_mi$group <- with(sampleTable_mi,
                             paste(treatment, time, sep = "_"))

# 5) DESeqDataSet y filtrado básico --------------------------

dds_mi <- DESeqDataSetFromMatrix(
  countData = mirna_counts,
  colData   = sampleTable_mi,
  design    = ~ time + treatment + time:treatment
)

# Prefiltrado un poco más estricto para miRNAs
keep_mi <- rowSums(counts(dds_mi) >= 5) >= 3
dds_mi  <- dds_mi[keep_mi, ]

dds_mi <- DESeq(dds_mi)

# 6) Transformación VST (para QC) ----------------------------

vsd_mi <- varianceStabilizingTransformation(dds_mi, blind = TRUE)

# 7) PCA -----------------------------------------------------

pcaData_mi    <- plotPCA(vsd_mi, intgroup = c("time","treatment"),
                         returnData = TRUE)
percentVar_mi <- round(100 * attr(pcaData_mi, "percentVar"))

p_pca_mi <- ggplot(pcaData_mi,
                   aes(PC1, PC2,
                       color = time,
                       shape = treatment,
                       label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar_mi[1], "%")) +
  ylab(paste0("PC2: ", percentVar_mi[2], "%")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA – miRNAs NE vs Control (6h / 24h)")

ggsave(file.path(output_dir_mi, "PCA_miRNA_multifactorial.png"),
       p_pca_mi, width = 7, height = 6, dpi = 300)

# 8) Heatmap de distancias entre muestras --------------------

sampleDists_mi      <- dist(t(assay(vsd_mi)))
sampleDistMatrix_mi <- as.matrix(sampleDists_mi)
colors_mi <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(file.path(output_dir_mi, "Heatmap_distancias_miRNA_multifactorial.png"),
    width = 2000, height = 1800, res = 250)
pheatmap(sampleDistMatrix_mi,
         clustering_distance_rows = sampleDists_mi,
         clustering_distance_cols = sampleDists_mi,
         col = colors_mi,
         main = "Distancias entre muestras – miRNAs")
dev.off()


library(DESeq2)
library(sva)        # para ComBat_seq
library(pheatmap)
library(RColorBrewer)

# 1) Cargar conteos
countdata <- read.csv("miR.Counts_C6_NE6_3v3_new.csv",
                      header = TRUE, row.names = 1)

# Seleccionar sólo las columnas de interés (ajusta si es necesario)
countdata <- countdata[, 6:13]
countdata <- countdata[, -3]   # según tu caso particular

# Matriz y +1 para evitar ceros
countdata <- as.matrix(countdata)
countdata <- countdata + 1

# Eliminar filas con NA
countdata_clean <- na.omit(countdata)

# 2) Definir condición y batch (para ESTAS 6 muestras, 24 h)
condition <- factor(c("Ctrl_24", "Ctrl_24", "Ctrl_24",
                      "NE_24",  "NE_24",  "NE_24"),
                    levels = c("Ctrl_24", "NE_24"))

condition <- factor(c("Ctrl_6", "Ctrl_6", "Ctrl_6",
                      "NE_6",  "NE_6",  "NE_6"),
                    levels = c("Ctrl_6", "NE_6"))

batch <- factor(c("Batch1", "Batch2", "Batch3",
                  "Batch1", "Batch2", "Batch3"))

coldata <- data.frame(
  row.names  = colnames(countdata_clean),
  condition  = condition,
  batch      = batch
)

# 3) Combat-Seq (opcional)
combat_counts <- ComBat_seq(counts = countdata_clean,
                            batch  = coldata$batch,
                            group  = coldata$condition)

# 4) Crear DESeqDataSet (usando counts corregidos)
dds <- DESeqDataSetFromMatrix(
  countData = combat_counts,
  colData   = coldata,
  design    = ~ condition
)

# 5) Filtro de baja expresión (para miRNA)
# Al menos 5 lecturas en al menos 3 muestras (ajusta a tu gusto)
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep,]

# En este punto, revisa:
nrow(dds)          # cuántos miRNAs quedan
# Si aquí sale 0, el filtro es demasiado estricto para tus datos.

# 6) DESeq2
dds <- DESeq(dds)

# 7) Transformación vst para QC
vsdata <- vst(dds, blind = FALSE)

# Distancias entre muestras
sampleDists <- dist(t(assay(vsdata)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsdata$condition
colnames(sampleDistMatrix) <- vsdata$condition

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         fontsize = 10)

plotPCA(vsdata, intgroup = "condition")

# 8) Resultados y tabla combinada
res   <- results(dds)         # contraste NE_24 vs Ctrl_24 (por defecto)
sigs  <- res                  # si quieres, puedes filtrar por padj mas abajo

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)

resdata <- merge(
  as.data.frame(sigs),
  as.data.frame(norm_counts),
  by = "row.names",
  sort = FALSE
)
names(resdata)[1] <- "miRNA_id"

dim(resdata)  # debería ser (N genes) x (estadísticos + 6 columnas de counts)
head(resdata)

# Crear status basado en padj y log2FC
resdata$status <- ifelse(
  resdata$padj < 0.05 & resdata$log2FoldChange >= 0.5, "Upregulated",
  ifelse(
    resdata$padj < 0.05 & resdata$log2FoldChange <= -0.5, "Downregulated",
    "NotSig"
  )
)

res_sig <- subset(resdata, status != "NotSig")

write.table(resdata,
            file = "DE_miRNA_All_with_status_NE6.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(res_sig,
            file = "DE_miRNA_Significant_only_NE6.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(subset(res_sig, status == "Upregulated"),
            "DE_miRNA_Up_NE6.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(subset(res_sig, status == "Downregulated"),
            "DE_miRNA_Down_NE6.txt", sep="\t", quote=FALSE, row.names=FALSE)

table(resdata$status)

#####tiempo análisis
# miRNA_temporal_summary.R
# Integra resultados DE de 6h y 24h (ya calculados)
# y genera:
#   - Tabla resumen por miRNA
#   - Clase temporal (Early / Late / Sustained / Non_DE)
#   - Patrón detallado (Sustained_UP, Switch_UP_to_DOWN, etc.)
#   - Barplot de conteos por patrón
# =========================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(readr)

# ------------ 0. INPUTS: CAMBIA AQUÍ LAS RUTAS  ------------

de_6_file  <- "DE_miRNA_All_with_status_NE6.txt"   # archivo 6h
de_24_file <- "DE_miRNA_All_with_status_NE24.txt"  # archivo 24h

outdir <- "miRNA_temporal_summary"
if (!dir.exists(outdir)) dir.create(outdir)

# Si tus txt son tabulados:
read_fun <- function(x) read.delim(x, header = TRUE, check.names = FALSE)
# Si son CSV, cambia por:
# read_fun <- function(x) read.csv(x, header = TRUE, check.names = FALSE)

# ------------ 1. Leer resultados de 6h y 24h  ------------

de6  <- read_fun(de_6_file)
de24 <- read_fun(de_24_file)

# Asegurar nombre de ID
if (!"miRNA_id" %in% colnames(de6)) {
  de6 <- de6 %>% rename(miRNA_id = 1)
}
if (!"miRNA_id" %in% colnames(de24)) {
  de24 <- de24 %>% rename(miRNA_id = 1)
}

# Chequeo básico
message("Columnas detectadas en 6h:");  print(colnames(de6))
message("Columnas detectadas en 24h:"); print(colnames(de24))

# ------------ 2. Normalizar la columna 'status' a UP/DOWN/NS ------------

# Esta función convierte el texto de `status`
# en una dirección uniforme: "UP", "DOWN" o "NS"
status_to_dir <- function(x) {
  x <- as.character(x)
  case_when(
    is.na(x)                          ~ "NS",
    grepl("not",  x, ignore.case=TRUE) ~ "NS",
    grepl("ns",   x, ignore.case=TRUE) ~ "NS",
    grepl("up",   x, ignore.case=TRUE) ~ "UP",
    grepl("down", x, ignore.case=TRUE) ~ "DOWN",
    TRUE                               ~ "NS"
  )
}

de6 <- de6 %>%
  mutate(direction = status_to_dir(status))

de24 <- de24 %>%
  mutate(direction = status_to_dir(status))

# Reducimos a columnas clave
de6_red <- de6 %>%
  dplyr::select(miRNA_id, log2FoldChange, padj, direction) %>%
  rename(log2FC_6h = log2FoldChange,
         padj_6h   = padj,
         dir_6h    = direction)

de24_red <- de24 %>%
  dplyr::select(miRNA_id, log2FoldChange, padj, direction) %>%
  rename(log2FC_24h = log2FoldChange,
         padj_24h   = padj,
         dir_24h    = direction)

# ------------ 3. Unir 6h y 24h, crear clases temporales ------------

summary_df <- full_join(de6_red, de24_red, by = "miRNA_id") %>%
  # Asegurar que las direcciones vacías queden como NS
  mutate(
    dir_6h  = ifelse(is.na(dir_6h),  "NS", dir_6h),
    dir_24h = ifelse(is.na(dir_24h), "NS", dir_24h)
  ) %>%
  # Clase temporal gruesa
  mutate(
    temporal_class = case_when(
      dir_6h != "NS" & dir_24h == "NS" ~ "Early_6h_only",
      dir_6h == "NS" & dir_24h != "NS" ~ "Late_24h_only",
      dir_6h != "NS" & dir_24h != "NS" ~ "Sustained_both",
      TRUE                             ~ "Non_DE"
    ),
    # Patrón más detallado usando dirección
    temporal_pattern = case_when(
      dir_6h == "UP"   & dir_24h == "UP"   ~ "Sustained_UP",
      dir_6h == "DOWN" & dir_24h == "DOWN" ~ "Sustained_DOWN",
      dir_6h == "UP"   & dir_24h == "DOWN" ~ "Switch_UP_to_DOWN",
      dir_6h == "DOWN" & dir_24h == "UP"   ~ "Switch_DOWN_to_UP",
      dir_6h != "NS"   & dir_24h == "NS"   ~ paste0("Early_", dir_6h),
      dir_6h == "NS"   & dir_24h != "NS"   ~ paste0("Late_", dir_24h),
      TRUE                                 ~ "Non_DE"
    )
  )

# ------------ 4. Guardar tabla resumen ------------

write.csv(summary_df,
          file.path(outdir, "miRNA_temporal_summary_table.csv"),
          row.names = FALSE)

message("Guardado: miRNA_temporal_summary_table.csv")

# ------------ 5. Gráfico de conteos por patrón temporal ------------

plot_df <- summary_df %>%
  filter(temporal_pattern != "Non_DE") %>%
  group_by(temporal_pattern) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df,
              aes(x = reorder(temporal_pattern, -n), y = n)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Resumen temporal de miRNAs",
         x = "Patrón temporal",
         y = "Número de miRNAs")
  
  ggsave(file.path(outdir, "miRNA_temporal_patterns_barplot.png"),
         p, width = 7, height = 5, dpi = 300)
  
  message("Guardado: miRNA_temporal_patterns_barplot.png")
} else {
  message("No hay miRNAs DE para graficar patrones temporales.")
}

# =========================================================
# miRNA_volcano_upset.R
# - Lee resultados DE de miRNA a 6h y 24h (ya calculados)
# - Hace volcanos con eje X simétrico y homogéneo
# - Genera UpSet de UP/DOWN por tiempo (UP_6h, DOWN_6h, UP_24h, DOWN_24h)
# =========================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(UpSetR)

# ---------- 0. INPUTS: ajusta rutas y tipo de lectura ----------

de_6_file  <- "DE_miRNA_All_with_status_NE6.txt"
de_24_file <- "DE_miRNA_All_with_status_NE24.txt"

outdir <- "miRNA_volcano_upset"
if (!dir.exists(outdir)) dir.create(outdir)

# Si tus archivos son tabulados:
read_fun <- function(x) read.delim(x, header = TRUE, check.names = FALSE)
# Si fueran CSV:
# read_fun <- function(x) read.csv(x, header = TRUE, check.names = FALSE)

alpha_cut <- 0.05
lfc_cut   <- 0.5   # |log2FC| para considerar UP/DOWN

# ---------- 1. Leer resultados de 6h y 24h ----------

de6  <- read_fun(de_6_file)
de24 <- read_fun(de_24_file)

if (!"miRNA_id" %in% colnames(de6))  de6  <- de6  %>% rename(miRNA_id = 1)
if (!"miRNA_id" %in% colnames(de24)) de24 <- de24 %>% rename(miRNA_id = 1)

# Normalizar status -> UP / DOWN / NS
status_to_dir <- function(x, padj, lfc,
                          alpha = alpha_cut, lfc_cut = lfc_cut) {
  x <- as.character(x)
  # Si ya tienes status bien armado, esto lo respeta
  dir_status <- case_when(
    grepl("up",   x, ignore.case = TRUE) ~ "UP",
    grepl("down", x, ignore.case = TRUE) ~ "DOWN",
    grepl("not",  x, ignore.case = TRUE) ~ "NS",
    grepl("ns",   x, ignore.case = TRUE) ~ "NS",
    TRUE                                  ~ NA_character_
  )
  # Si status es NA o raro, usa padj + log2FC
  dir_final <- ifelse(
    is.na(dir_status),
    ifelse(!is.na(padj) & padj < alpha & lfc >=  lfc_cut, "UP",
           ifelse(!is.na(padj) & padj < alpha & lfc <= -lfc_cut, "DOWN", "NS")),
    dir_status
  )
  dir_final
}

de6 <- de6 %>%
  mutate(
    direction = status_to_dir(status, padj, log2FoldChange,
                              alpha = alpha_cut, lfc_cut = lfc_cut)
  )

de24 <- de24 %>%
  mutate(
    direction = status_to_dir(status, padj, log2FoldChange,
                              alpha = alpha_cut, lfc_cut = lfc_cut)
  )

# ---------- 2. Función de volcano CON eje simétrico ----------

make_volcano <- function(df, title, outfile,
                         xlim_global = NULL,
                         top_labels = 15) {
  
  df_sub <- df %>%
    mutate(
      padj_plot    = ifelse(is.na(padj) | padj <= 0, NA, padj),
      negLog10Padj = -log10(padj_plot)
    ) %>%
    filter(!is.na(negLog10Padj))
  
  if (nrow(df_sub) == 0) {
    message("No hay miRNAs para plotear en: ", title)
    return(NULL)
  }
  
  # Calcular límite simétrico si no se pasa explícito
  if (is.null(xlim_global)) {
    max_lfc <- max(abs(df_sub$log2FoldChange), na.rm = TRUE)
    xlim_global <- ceiling(max_lfc)  # redondea hacia arriba
  }
  
  label_df <- bind_rows(
    df_sub %>% filter(direction == "UP")   %>%
      arrange(desc(log2FoldChange)) %>% head(top_labels),
    df_sub %>% filter(direction == "DOWN") %>%
      arrange(log2FoldChange) %>% head(top_labels),
    df_sub %>% arrange(desc(negLog10Padj)) %>% head(top_labels)
  ) %>% distinct(miRNA_id, .keep_all = TRUE)
  
  p <- ggplot(df_sub,
              aes(x = log2FoldChange,
                  y = negLog10Padj,
                  color = direction)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = c(
      "DOWN" = "royalblue",
      "NS"   = "grey70",
      "UP"   = "firebrick"
    )) +
    xlim(-xlim_global, xlim_global) +
    labs(
      title = title,
      x     = "log2 Fold Change (NE vs Ctrl)",
      y     = expression(-log[10]("padj"))
    ) +
    theme_bw() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      text         = element_text(size = 12)
    ) +
    geom_text_repel(
      data            = label_df,
      aes(label = miRNA_id),
      size            = 3,
      max.overlaps    = Inf,
      min.segment.length = 0.1,
      box.padding     = 0.3,
      point.padding   = 0.2,
      seed            = 123,
      fontface        = "bold"
    )
  
  png(file.path(outdir, outfile),
      width = 2200, height = 2000, res = 300)
  print(p)
  dev.off()
  
  invisible(p)
}

# ---------- 3. Volcanos 6 h y 24 h con MISMO xlim ----------

# Límite global usando ambos tiempos
max_lfc_global <- max(
  abs(c(de6$log2FoldChange, de24$log2FoldChange)),
  na.rm = TRUE
)
max_lfc_global <- ceiling(max_lfc_global) # ej: 3.2 -> 4

make_volcano(de6,
             title       = "miRNAs 6 h – NE vs Ctrl",
             outfile     = "Volcano_miRNA_6h.png",
             xlim_global = max_lfc_global)

make_volcano(de24,
             title       = "miRNAs 24 h – NE vs Ctrl",
             outfile     = "Volcano_miRNA_24h.png",
             xlim_global = max_lfc_global)

###nuevo volcano y upsetplot con otros cutoffs
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(UpSetR)

# ---------- 0. INPUTS: ajusta rutas y tipo de lectura ----------

de_6_file  <- "DE_miRNA_All_with_status_NE6.txt"
de_24_file <- "DE_miRNA_All_with_status_NE24.txt"

outdir <- "miRNA_volcano_upset"
if (!dir.exists(outdir)) dir.create(outdir)

# Si tus archivos son tabulados:
read_fun <- function(x) read.delim(x, header = TRUE, check.names = FALSE)
# Si fueran CSV:
# read_fun <- function(x) read.csv(x, header = TRUE, check.names = FALSE)

alpha_cut <- 0.05
lfc_cut   <- 0.5   # |log2FC| para considerar UP/DOWN

# ---------- 1. Leer resultados de 6h y 24h ----------

de6  <- read_fun(de_6_file)
de24 <- read_fun(de_24_file)

if (!"miRNA_id" %in% colnames(de6))  de6  <- de6  %>% rename(miRNA_id = 1)
if (!"miRNA_id" %in% colnames(de24)) de24 <- de24 %>% rename(miRNA_id = 1)

# Normalizar status -> UP / DOWN / NS
status_to_dir <- function(x, padj, lfc,
                          alpha = alpha_cut, lfc_cut = lfc_cut) {
  x <- as.character(x)
  # Si ya tienes status bien armado, esto lo respeta
  dir_status <- dplyr::case_when(
    grepl("up",   x, ignore.case = TRUE) ~ "UP",
    grepl("down", x, ignore.case = TRUE) ~ "DOWN",
    grepl("not",  x, ignore.case = TRUE) ~ "NS",
    grepl("ns",   x, ignore.case = TRUE) ~ "NS",
    TRUE                                  ~ NA_character_
  )
  # Si status es NA o raro, usa padj + log2FC
  dir_final <- ifelse(
    is.na(dir_status),
    ifelse(!is.na(padj) & padj < alpha & lfc >=  lfc_cut, "UP",
           ifelse(!is.na(padj) & padj < alpha & lfc <= -lfc_cut, "DOWN", "NS")),
    dir_status
  )
  dir_final
}

de6 <- de6 %>%
  mutate(
    direction = status_to_dir(status, padj, log2FoldChange,
                              alpha = alpha_cut, lfc_cut = lfc_cut)
  )

de24 <- de24 %>%
  mutate(
    direction = status_to_dir(status, padj, log2FoldChange,
                              alpha = alpha_cut, lfc_cut = lfc_cut)
  )

# ---------- 2. Función de volcano CON eje simétrico
#     y POCAS etiquetas seleccionadas ----------

make_volcano <- function(df, title, outfile,
                         xlim_global = NULL,
                         top_up = 8,
                         top_down = 8,
                         top_most_sig = 5,
                         min_neglog10 = 2) {
  
  df_sub <- df %>%
    mutate(
      padj_plot    = ifelse(is.na(padj) | padj <= 0, NA, padj),
      negLog10Padj = -log10(padj_plot)
    ) %>%
    filter(!is.na(negLog10Padj))
  
  if (nrow(df_sub) == 0) {
    message("No hay miRNAs para plotear en: ", title)
    return(NULL)
  }
  
  # Límite simétrico global
  if (is.null(xlim_global)) {
    max_lfc <- max(abs(df_sub$log2FoldChange), na.rm = TRUE)
    xlim_global <- ceiling(max_lfc)  # redondea hacia arriba
  }
  
  # --- selección de etiquetas ---
  # Sólo miRNAs DE (UP/DOWN) y por encima de cierto -log10(padj)
  de_only <- df_sub %>%
    filter(direction != "NS",
           negLog10Padj >= min_neglog10)
  
  # Top por log2FC (UP y DOWN)
  label_up <- de_only %>%
    filter(direction == "UP") %>%
    arrange(desc(log2FoldChange)) %>%
    head(top_up)
  
  label_down <- de_only %>%
    filter(direction == "DOWN") %>%
    arrange(log2FoldChange) %>%
    head(top_down)
  
  # Top por significancia global
  label_sig <- de_only %>%
    arrange(desc(negLog10Padj)) %>%
    head(top_most_sig)
  
  label_df <- bind_rows(label_up, label_down, label_sig) %>%
    distinct(miRNA_id, .keep_all = TRUE)
  
  # --- plot ---
  p <- ggplot(df_sub,
              aes(x = log2FoldChange,
                  y = negLog10Padj,
                  color = direction)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c(
      "DOWN" = "royalblue",
      "NS"   = "grey80",
      "UP"   = "firebrick"
    )) +
    xlim(-xlim_global, xlim_global) +
    labs(
      title = title,
      x     = "log2 Fold Change (NE vs Ctrl)",
      y     = expression(-log[10]("padj"))
    ) +
    theme_bw() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      text         = element_text(size = 12)
    ) +
    geom_text_repel(
      data            = label_df,
      aes(label = miRNA_id),
      size            = 3,
      max.overlaps    = 50,       # límite para que no explote
      min.segment.length = 0.1,
      box.padding     = 0.3,
      point.padding   = 0.2,
      seed            = 123,
      fontface        = "bold"
    )
  
  png(file.path(outdir, outfile),
      width = 2200, height = 2000, res = 300)
  print(p)
  dev.off()
  
  invisible(p)
}

# ---------- 3. Volcanos 6 h y 24 h con MISMO xlim ----------

# Límite global usando ambos tiempos
max_lfc_global <- max(
  abs(c(de6$log2FoldChange, de24$log2FoldChange)),
  na.rm = TRUE
)
max_lfc_global <- ceiling(max_lfc_global) # ej: 3.2 -> 4

make_volcano(de6,
             title       = "miRNAs 6 h – NE vs Ctrl",
             outfile     = "Volcano_miRNA_6h.png",
             xlim_global = max_lfc_global)

make_volcano(de24,
             title       = "miRNAs 24 h – NE vs Ctrl",
             outfile     = "Volcano_miRNA_24h.png",
             xlim_global = max_lfc_global)

library(UpSetR)
library(dplyr)
library(readr)

# Cargar tus archivos procesados (los que tú generaste)
res6  <- read.table("DE_miRNA_All_with_status_NE6.txt",  header = TRUE, sep = "\t")
res24 <- read.table("DE_miRNA_All_with_status_NE24.txt", header = TRUE, sep = "\t")

# Filtrar solo los Up y Down reales (padj < 0.05 & |log2FC| >= 0.5)
up_6h    <- res6  %>% filter(status == "Upregulated")   %>% pull(miRNA_id)
down_6h  <- res6  %>% filter(status == "Downregulated") %>% pull(miRNA_id)

up_24h   <- res24 %>% filter(status == "Upregulated")   %>% pull(miRNA_id)
down_24h <- res24 %>% filter(status == "Downregulated") %>% pull(miRNA_id)

# Listas para UpSet
updown_list <- list(
  UP_6h    = up_6h,
  DOWN_6h  = down_6h,
  UP_24h   = up_24h,
  DOWN_24h = down_24h
)

updown_data <- fromList(updown_list)

upset(
  updown_data,
  nsets = 4,
  nintersects = 10,
  order.by = "freq",
  mainbar.y.label = "Number of miRNAs (UP/DOWN)",
  sets.x.label    = "miRNAs per group (UP/DOWN)"
)

library(DESeq2)
library(dplyr)
library(tibble)
library(Mfuzz)
library(Biobase)
library(sva)   # <-- para ComBat_seq

# -------- 0. INPUTS --------

file_6h  <- "miR.Counts_C6_NE6_3v3_new.csv"
file_24h <- "miR.Counts_C24_NE24_3v3_new.csv"

outdir_lrt   <- "miRNA_LRT_temporal"
outdir_mfuzz <- "miRNA_Mfuzz_temporal"

if (!dir.exists(outdir_lrt))   dir.create(outdir_lrt)
if (!dir.exists(outdir_mfuzz)) dir.create(outdir_mfuzz)

# Parámetros de filtrado (miRNAs)
min_counts  <- 5    # >=5 lecturas
min_samples <- 3    # en al menos 3 muestras

# -------- 1. Cargar matrices de counts (6h y 24h) --------

counts_6_raw  <- read.csv(file_6h,  header = TRUE,
                          row.names = 1, check.names = FALSE)
counts_24_raw <- read.csv(file_24h, header = TRUE,
                          row.names = 1, check.names = FALSE)

# Selección de columnas de counts (ajusta si cambia el formato)
# Aquí asumo que 6:13 son H1, H2, SQ1, SQ2, SQ13, SQ14 (6h)
# y análogo para 24h: H3, H4, SQ3, SQ4, SQ15, SQ16
counts_6  <- counts_6_raw[, 6:13]
counts_6  <- counts_6[, -3]   # deja 6 columnas (H1,H2,SQ1,SQ2,SQ13,SQ14)

counts_24 <- counts_24_raw[, 6:13]
counts_24 <- counts_24[, -3]  # análogo 24h

# Convierte a matriz (del subset, no del raw completo)
counts_6  <- as.matrix(counts_6_raw)
counts_24 <- as.matrix(counts_24_raw)

# Mantener sólo miRNAs presentes en ambos tiempos
common_ids <- intersect(rownames(counts_6), rownames(counts_24))
counts_6   <- counts_6[common_ids, ]
counts_24  <- counts_24[common_ids, ]

# Renombrar columnas según condición (usando tu mapeo H/SQ)
# 6 h:
# H1 -> Ctrl6_rep1 (N1)
# H2 -> NE6_rep1   (N1)
# SQ1 -> Ctrl6_rep2 (N2)
# SQ2 -> NE6_rep2   (N2)
# SQ13 -> Ctrl6_rep4 (N4)
# SQ14 -> NE6_rep4   (N4)
colnames(counts_6) <- c("Ctrl6_rep1","NE6_rep1",
                        "Ctrl6_rep2","NE6_rep2",
                        "Ctrl6_rep4","NE6_rep4")

# 24 h:
# H3 -> Ctrl24_rep1 (N1)
# H4 -> NE24_rep1   (N1)
# SQ3 -> Ctrl24_rep2 (N2)
# SQ4 -> NE24_rep2   (N2)
# SQ15 -> Ctrl24_rep4 (N4)
# SQ16 -> NE24_rep4   (N4)
colnames(counts_24) <- c("Ctrl24_rep1","NE24_rep1",
                         "Ctrl24_rep2","NE24_rep2",
                         "Ctrl24_rep4","NE24_rep4")

# Reordenar para tener primero todos los Ctrl y luego los NE por tiempo
counts_6  <- counts_6[,  c("Ctrl6_rep1","Ctrl6_rep2","Ctrl6_rep4",
                           "NE6_rep1","NE6_rep2","NE6_rep4")]
counts_24 <- counts_24[, c("Ctrl24_rep1","Ctrl24_rep2","Ctrl24_rep4",
                           "NE24_rep1","NE24_rep2","NE24_rep4")]

# Combinar en una sola matriz (12 muestras)
counts_all <- cbind(counts_6, counts_24)
samples_all <- colnames(counts_all)

# -------- 2. colData multifactorial (time + treatment + batch) --------

sampleTable_all <- data.frame(
  sample    = samples_all,
  treatment = factor(c(rep("Control", 3),  # Ctrl6 (3 reps)
                       rep("NE",      3),  # NE6
                       rep("Control", 3),  # Ctrl24
                       rep("NE",      3)), # NE24
                     levels = c("Control","NE")),
  time      = factor(c(rep("6h", 6),
                       rep("24h", 6)),
                     levels = c("6h","24h"))
)

# Batches según corridas de secuenciación:
# N1 -> H1, H2, H3, H4  (rep1 de cada condición/tiempo)
# N2 -> SQ1, SQ2, SQ3, SQ4  (rep2)
# N4 -> SQ13, SQ14, SQ15, SQ16  (rep4)
sampleTable_all$batch <- factor(c("N1","N2","N4",
                                  "N1","N2","N4",
                                  "N1","N2","N4",
                                  "N1","N2","N4"))

rownames(sampleTable_all) <- sampleTable_all$sample

# -------- 3. ComBat_seq para corregir batch --------

# Grupo biológico: combinación de time y treatment
group_all <- with(sampleTable_all,
                  interaction(time, treatment, drop = TRUE))

combat_counts_all <- ComBat_seq(
  counts = as.matrix(counts_all),
  batch  = sampleTable_all$batch,
  group  = group_all
)

# -------- 4. DESeqDataSet (con counts corregidos por batch) --------

dds_all <- DESeqDataSetFromMatrix(
  countData = combat_counts_all,
  colData   = sampleTable_all,
  design    = ~ time + treatment + time:treatment
)
# (no incluimos batch en el diseño porque ya fue corregido por ComBat_seq)

# Filtro suave para miRNAs
keep_all <- rowSums(counts(dds_all) >= min_counts) >= min_samples
dds_all  <- dds_all[keep_all, ]

# -------- 5. LRT: ¿cambia el efecto de NE entre 6h y 24h? --------

dds_all <- DESeq(dds_all,
                 test    = "LRT",
                 reduced = ~ time + treatment)  # quitamos la interacción

res_LRT <- results(dds_all)
res_LRT_tbl <- as.data.frame(res_LRT) %>%
  rownames_to_column("miRNA_id") %>%
  arrange(padj)

write.csv(res_LRT_tbl,
          file.path(outdir_lrt,
                    "miRNA_LRT_time_treatment_interaction_full.csv"),
          row.names = FALSE)

alpha_lrt <- 0.05
sig_LRT <- res_LRT_tbl %>%
  filter(!is.na(padj) & padj < alpha_lrt)

write.csv(sig_LRT,
          file.path(outdir_lrt,
                    paste0("miRNA_LRT_signif_padj_", alpha_lrt, ".csv")),
          row.names = FALSE)

message("LRT completado. Resultados en carpeta: ", outdir_lrt)

# -------- 6. VST y matriz de promedios por grupo (para Mfuzz) --------

vsd_all <- varianceStabilizingTransformation(dds_all, blind = TRUE)
vst_mat <- assay(vsd_all)

group_all2 <- with(colData(vsd_all),
                   paste(treatment, time, sep = "_"))
groups <- unique(group_all2)

vst_group_means <- sapply(groups, function(g) {
  rowMeans(vst_mat[, group_all2 == g, drop = FALSE])
})
colnames(vst_group_means) <- groups   # Control_6h, NE_6h, Control_24h, NE_24h

write.table(vst_group_means,
            file = file.path(outdir_mfuzz,
                             "miRNA_VST_group_means_for_Mfuzz.txt"),
            sep = "\t", quote = FALSE)

# -------- 7. Mfuzz: clustering de patrones temporales --------

eset <- ExpressionSet(assayData = vst_group_means)
eset <- standardise(eset)

m_val <- mestimate(eset)
c_clusters <- 4   # ajusta si quieres

set.seed(123)
cl <- mfuzz(eset, c = c_clusters, m = m_val)

png(file.path(outdir_mfuzz, "miRNA_Mfuzz_clusters.png"),
    width = 2000, height = 2000, res = 250)
mfuzz.plot(eset, cl = cl,
           mfrow       = c(2, 2),
           time.labels = colnames(vst_group_means),
           xlab = "Condición / tiempo",
           ylab = "Expresión VST estandarizada")
dev.off()

membership_tbl <- as.data.frame(cl$membership) %>%
  rownames_to_column("miRNA_id")

cluster_tbl <- data.frame(
  miRNA_id = names(cl$cluster),
  cluster  = cl$cluster
)

write.csv(membership_tbl,
          file.path(outdir_mfuzz, "miRNA_Mfuzz_membership.csv"),
          row.names = FALSE)

write.csv(cluster_tbl,
          file.path(outdir_mfuzz, "miRNA_Mfuzz_clusters.csv"),
          row.names = FALSE)

message("Mfuzz completado. Resultados en carpeta: ", outdir_mfuzz)
message("Script LRT + Mfuzz para miRNAs finalizado.")

summary(res_LRT)
table(res_LRT$pvalue < 0.05, useNA = "ifany")
sum(res_LRT$padj < 0.1, na.rm = TRUE)

###intento de genración de clusters

# =========================================================
# miRNA_Mfuzz_DEonly.R
# - Usa solo miRNAs DE (UP/DOWN en 6h o 24h)
# - Lee:
#     * DE_miRNA_All_with_status_NE6.txt
#     * DE_miRNA_All_with_status_NE24.txt
#     * miRNA_VST_group_means_for_Mfuzz.txt
# - Corre Mfuzz sólo sobre esos miRNAs
# - Genera:
#     * miRNA_Mfuzz_clusters.csv (miRNA_id, cluster)
#     * miRNA_Mfuzz_membership.csv
# =========================================================

library(dplyr)
library(tibble)
library(Mfuzz)
library(Biobase)

# -------- 0. INPUTS --------

# Archivos DE por tiempo (los que tú generaste con status)
de_6_file  <- "DE_miRNA_All_with_status_NE6.txt"
de_24_file <- "DE_miRNA_All_with_status_NE24.txt"

# Carpeta donde ya guardaste la matriz VST por grupo
outdir_mfuzz <- "miRNA_Mfuzz_temporal"

vst_file <- file.path(outdir_mfuzz,
                      "miRNA_VST_group_means_for_Mfuzz.txt")

if (!dir.exists(outdir_mfuzz)) dir.create(outdir_mfuzz)

# -------- 1. Leer DE 6h y 24h, definir miRNAs DE --------

# tus archivos son tabulados:
read_fun <- function(x) read.delim(x, header = TRUE, check.names = FALSE)

de6  <- read_fun(de_6_file)
de24 <- read_fun(de_24_file)

if (!"miRNA_id" %in% colnames(de6))  de6  <- de6  %>% rename(miRNA_id = 1)
if (!"miRNA_id" %in% colnames(de24)) de24 <- de24 %>% rename(miRNA_id = 1)

# miRNAs DE (status != NotSig)
de6_ids  <- de6  %>% filter(status != "NotSig") %>% pull(miRNA_id)
de24_ids <- de24 %>% filter(status != "NotSig") %>% pull(miRNA_id)

miRNA_DE_union <- union(de6_ids, de24_ids)

cat("miRNAs DE en 6h:",  length(de6_ids),  "\n")
cat("miRNAs DE en 24h:", length(de24_ids), "\n")
cat("Unión de miRNAs DE (6h U 24h):", length(miRNA_DE_union), "\n")

# -------- 2. Leer matriz VST por grupo --------

if (!file.exists(vst_file)) {
  stop("No se encontró el archivo VST: ", vst_file,
       "\nPrimero corre el script que genera miRNA_VST_group_means_for_Mfuzz.txt")
}

vst_group_means <- as.matrix(
  read.delim(vst_file, header = TRUE,
             row.names = 1, check.names = FALSE)
)

cat("Dimensión VST total: ",
    nrow(vst_group_means), "miRNAs x", ncol(vst_group_means), "condiciones\n")

# -------- 3. Filtrar VST sólo a miRNAs DE --------

common_ids <- intersect(rownames(vst_group_means), miRNA_DE_union)

cat("miRNAs DE con VST disponible (intersección):",
    length(common_ids), "\n")

if (length(common_ids) < 5) {
  stop("Muy pocos miRNAs DE tras la intersección con la matriz VST. ",
       "Revisa que los IDs coincidan entre DE y VST.")
}

vst_de <- vst_group_means[common_ids, , drop = FALSE]

# -------- 4. Mfuzz sobre miRNAs DE --------

# ExpressionSet y estandarización
eset_de <- ExpressionSet(assayData = vst_de)
eset_de <- standardise(eset_de)

# Estimar parámetro de fuzziness
m_val <- mestimate(eset_de)

# Número de clusters (ajusta a gusto: 3–6 suele ser razonable)
c_clusters <- 6

set.seed(123)
cl_de <- mfuzz(eset_de, c = c_clusters, m = m_val)

# Plot de los clusters
png(file.path(outdir_mfuzz, "miRNA_Mfuzz_clusters_DEonly.png"),
    width = 2000, height = 2000, res = 250)
mfuzz.plot(eset_de, cl = cl_de,
           mfrow       = c(2, 2),
           time.labels = colnames(vst_de))
dev.off()

# -------- 5. Guardar membresía y cluster (formato compatible con mega tabla) --------

membership_tbl <- as.data.frame(cl_de$membership) %>%
  rownames_to_column("miRNA_id")

cluster_tbl <- data.frame(
  miRNA_id = names(cl_de$cluster),
  cluster  = cl_de$cluster
)

write.csv(membership_tbl,
          file.path(outdir_mfuzz, "miRNA_Mfuzz_membership.csv"),
          row.names = FALSE)

write.csv(cluster_tbl,
          file.path(outdir_mfuzz, "miRNA_Mfuzz_clusters.csv"),
          row.names = FALSE)

cat("Mfuzz DE-only completado. Archivos generados en:", outdir_mfuzz, "\n")

##Ahora si continuar con esto
# =============================================================
# miRNA_temporal_summary_FINAL.R
# Integra:
#   - DESeq2 6h
#   - DESeq2 24h
#   - dirección UP/DOWN/NS
#   - categoría temporal (Early/Late/Sustained/Switch)
#   - Mfuzz cluster (solo miRNAs DE)
# Genera la mega tabla final
# =============================================================

library(dplyr)
library(tibble)
library(readr)

# -------- 0. INPUTS --------

de6_file  <- "DE_miRNA_All_with_status_NE6.txt"
de24_file <- "DE_miRNA_All_with_status_NE24.txt"
mfuzz_file <- "miRNA_Mfuzz_temporal/miRNA_Mfuzz_clusters.csv"

out_file <- "miRNA_temporal_MEGA_TABLE.csv"

# -------- 1. Leer datos DE 6h y 24h --------

read_fun <- function(x) read.delim(x, header = TRUE, check.names = FALSE)

de6  <- read_fun(de6_file)
de24 <- read_fun(de24_file)

if (!"miRNA_id" %in% colnames(de6))
  de6 <- de6 %>% rename(miRNA_id = 1)

if (!"miRNA_id" %in% colnames(de24))
  de24 <- de24 %>% rename(miRNA_id = 1)

# -------- 2. Mantener columnas relevantes --------

de6_clean <- de6 %>%
  dplyr::select(miRNA_id, log2FoldChange, padj, status) %>%
  rename(
    log2FC_6h = log2FoldChange,
    padj_6h   = padj,
    status_6h = status
  )

de24_clean <- de24 %>%
  dplyr::select(miRNA_id, log2FoldChange, padj, status) %>%
  rename(
    log2FC_24h = log2FoldChange,
    padj_24h   = padj,
    status_24h = status
  )

# -------- 3. Merge 6h + 24h --------

mega <- full_join(de6_clean, de24_clean, by = "miRNA_id")

# -------- 4. Categorías temporales --------

mega <- mega %>%
  mutate(
    temporal_category = case_when(
      status_6h != "NotSig" & status_24h == "NotSig" ~ "Early",
      status_6h == "NotSig" & status_24h != "NotSig" ~ "Late",
      status_6h != "NotSig" & status_24h != "NotSig" &
        sign(log2FC_6h) == sign(log2FC_24h) ~ "Sustained",
      status_6h != "NotSig" & status_24h != "NotSig" &
        sign(log2FC_6h) != sign(log2FC_24h) ~ "Switch",
      TRUE ~ "Not_DE"
    )
  )

# -------- 5. Agregar cluster Mfuzz --------

mfuzz <- read.delim(mfuzz_file, header = TRUE)

mega <- mega %>%
  left_join(mfuzz, by = "miRNA_id")

mfuzz <- read.csv("miRNA_Mfuzz_temporal/miRNA_Mfuzz_clusters.csv",
                  header = TRUE, check.names = FALSE)

mega <- mega %>%
  left_join(mfuzz, by = "miRNA_id")


# -------- 6. Guardar tabla final --------

write.csv(mega, out_file, row.names = FALSE)

cat("MEGA TABLA generada en:", out_file, "\n")
