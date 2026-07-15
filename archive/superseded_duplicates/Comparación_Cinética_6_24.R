###     Análisis de Expresión Diferencial Multifactorial con DESeq2          ###
###        Respuesta de Cardiomiocitos a Norepinefrina (6h vs 24h)           ###

# Para importar datos de Salmon
library(tximport)
library(rtracklayer)
library(DESeq2)
library(apeglm) # Para LFC shrinkage
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Rn.eg.db) # Base de datos de anotación para rata

# --- 1. Creación del mapa de Transcrito a Gen (tx2gene) ---
# Este paso es crucial para sumar las cuentas de los transcritos a nivel de gen.
# Se lee un archivo GTF de anotación y se extrae la correspondencia entre
# los IDs de transcritos y los IDs de genes.

gtf_file <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis/merged_CNE_6_24_HISAT2.annotated.gtf"
tx2gene_file <- "tx2gene.tsv"

# Para no repetir este paso, que puede ser lento, comprobamos si el archivo ya existe.
if (!file.exists(tx2gene_file)) {
  message("Creando el archivo tx2gene desde el GTF...")
  gtf_data <- rtracklayer::import(gtf_file)
  gtf_df <- as.data.frame(gtf_data)
  
  tx2gene <- gtf_df %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(transcript_id, gene_id) %>%
    distinct()
  
  # Nos aseguramos de que no haya NAs que puedan causar problemas
  tx2gene <- na.omit(tx2gene)
  
  write_tsv(tx2gene, tx2gene_file)
  message("Archivo tx2gene.tsv creado.")
} else {
  message("Usando el archivo tx2gene.tsv existente.")
}

# --- 2. Preparación de Metadatos y Rutas de Archivos ---
# Aquí definimos TODAS nuestras muestras y creamos una única tabla de metadatos
# que describe el experimento completo. Esta es la clave del análisis multifactorial.

base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/salmon_quants_NE624/"

# Definimos TODAS las muestras en el orden que queramos procesar.
# Es buena práctica agrupar por condición.
samples <- c("Ctrl6_rep1", "Ctrl6_rep2", "Ctrl6_rep4",
             "NE6_rep1", "NE6_rep2", "NE6_rep4",
             "Ctrl24_rep1", "Ctrl24_rep2", "Ctrl24_rep4",
             "NE24_rep1", "NE24_rep2", "NE24_rep4")

# Creamos la tabla de metadatos (colData) para el experimento completo
sampleTable <- data.frame(
  sample = samples,
  treatment = factor(rep(c("Control", "NE", "Control", "NE"), each = 3)),
  time = factor(rep(c("6h", "24h"), each = 6))
)
rownames(sampleTable) <- samples

# ¡Importante! Establecemos los niveles de referencia para que los contrastes
# sean intuitivos (e.g., NE vs Control, 24h vs 6h).
sampleTable$treatment <- relevel(sampleTable$treatment, ref = "Control")
sampleTable$time <- relevel(sampleTable$time, ref = "6h")

print(sampleTable)

# Generamos las rutas a los archivos quant.sf de Salmon
files <- file.path(base_dir, paste0(samples, "_quant"), "quant.sf")
names(files) <- samples

# Verificación de que todos los archivos existen
if (!all(file.exists(files))) {
  stop("¡Faltan algunos archivos quant.sf! Revisa las rutas y los nombres de las muestras.")
}

# --- 3. Importación de datos con tximport ---
# Leemos los datos de cuantificación de Salmon y los colapsamos a nivel de gen.

tx2gene <- read_tsv(tx2gene_file, col_names = c("transcript_id", "gene_id"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


# --- 4. Creación del objeto DESeq2 y Ejecución del Análisis ---
# Aquí configuramos el objeto DESeqDataSet con nuestro diseño de interacción
# y ejecutamos el pipeline de DESeq2.

# El diseño ~ time + treatment + time:treatment modela los efectos principales
# del tiempo y el tratamiento, y lo más importante, la INTERACCIÓN entre ellos.
dds <- DESeqDataSetFromTximport(txi,
                                colData = sampleTable,
                                design = ~ time + treatment + time:treatment)

# Pre-filtrado: mantenemos genes que tienen al menos 10 cuentas en total.
# Esto reduce el ruido y acelera el análisis.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Ejecutamos el análisis de DESeq2 (normalización, dispersión, ajuste de modelo)
# Esto usa el test de Wald por defecto.
dds <- DESeq(dds)
resultsNames(dds) # Muestra los coeficientes que podemos testear

# --- 5. Control de Calidad (QC) Post-Normalización ---
# Es fundamental visualizar cómo se agrupan las muestras para detectar outliers
# o problemas en el experimento.

# Transformación de datos para visualización (estabiliza la varianza)
vsd <- vst(dds, blind = FALSE)

# Gráfico de PCA
# Creamos una variable 'group' para colorear por las 4 condiciones
sampleTable$group <- paste(sampleTable$treatment, sampleTable$time, sep="_")
plotPCA(vsd, intgroup = "group") + 
  geom_text_repel(aes(label=name)) +
  labs(title="Análisis de Componentes Principales (PCA)") +
  theme_bw()

# Heatmap de distancias entre muestras
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Heatmap of Distances between Samples")

# --- 6. Extracción de Resultados Específicos (Contrastes) ---
# Ahora que el modelo está ajustado, podemos hacerle preguntas específicas.

# Directorio para guardar los resultados
output_dir <- "DESeq2_Multifactorial_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# A) Efecto de NE vs Control a las 6h (efecto base del tratamiento)
res_NE_vs_Ctrl_6h <- results(dds, name = "treatment_NE_vs_Control", alpha = 0.05)
res_NE_vs_Ctrl_6h <- lfcShrink(dds, coef = "treatment_NE_vs_Control", type = "apeglm")
summary(res_NE_vs_Ctrl_6h)
write.csv(as.data.frame(res_NE_vs_Ctrl_6h), file.path(output_dir, "res_NE_vs_Ctrl_6h.csv"))

# B) Cambios entre 24h y 6h en las muestras CONTROL (efecto del tiempo)
res_Ctrl_24h_vs_6h <- results(dds, name = "time_24h_vs_6h", alpha = 0.05)
res_Ctrl_24h_vs_6h <- lfcShrink(dds, coef = "time_24h_vs_6h", type = "apeglm")
summary(res_Ctrl_24h_vs_6h)
write.csv(as.data.frame(res_Ctrl_24h_vs_6h), file.path(output_dir, "res_Ctrl_24h_vs_6h.csv"))

# C) Efecto de NE vs Control a las 24h
# Este efecto es la suma del efecto base del tratamiento MÁS el término de interacción
res_NE_vs_Ctrl_24h <- results(dds, contrast = list(c("treatment_NE_vs_Control", "time24h.treatmentNE")), alpha = 0.05)
# Para lfcShrink en un contraste, se necesita un poco más de trabajo o usar un coeficiente específico si es posible.
# Por simplicidad, guardamos el resultado sin shrinkage.
summary(res_NE_vs_Ctrl_24h)
write.csv(as.data.frame(res_NE_vs_Ctrl_24h), file.path(output_dir, "res_NE_vs_Ctrl_24h.csv"))

# D) Término de INTERACCIÓN: genes cuya respuesta a NE es DIFERENTE a las 24h vs 6h (la CINÉTICA)
res_interaction <- results(dds, name = "time24h.treatmentNE", alpha = 0.05)
res_interaction <- lfcShrink(dds, coef = "time24h.treatmentNE", type = "apeglm")
summary(res_interaction)
write.csv(as.data.frame(res_interaction), file.path(output_dir, "res_interaction_kinetics.csv"))

# --- 7. Test de Razón de Verosimilitud (LRT) para Análisis de Cinética Global ---
# Este test identifica CUALQUIER gen que tenga un patrón de expresión a lo largo
# del tiempo que sea dependiente del tratamiento. Es ideal para encontrar perfiles cinéticos.

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ time + treatment)
res_lrt <- results(dds_lrt)
res_lrt_sig <- subset(res_lrt, padj < 0.05)

message(paste(nrow(res_lrt_sig), "genes muestran una respuesta cinética significativa al tratamiento (p.adj < 0.05)."))
write.csv(as.data.frame(res_lrt_sig), file.path(output_dir, "res_LRT_all_kinetic_genes.csv"))

# Visualización de los genes con cinética: Heatmap
# Tomamos los 50 genes más significativos del LRT para visualizarlos
top_kinetic_genes <- rownames(res_lrt_sig)[order(res_lrt_sig$padj)[1:50]]
mat <- assay(vsd)[top_kinetic_genes, ]
mat <- mat - rowMeans(mat) # Centrar por la media de cada gen

pheatmap(mat, 
         annotation_col = as.data.frame(colData(dds)[, c("treatment", "time")]),
         main = "Top 50 Kinetic Response Genes (LRT)",
         scale = "row") # Escalar filas para ver patrones

# --- 8. Visualización de Genes Individuales ---
# Graficar un gen específico nos permite entender su comportamiento a través de las 4 condiciones.

# Escogemos el gen más significativo del test de interacción
top_interaction_gene <- rownames(res_interaction)[which.min(res_interaction$padj)]

plotCounts(dds, gene = top_interaction_gene, intgroup = c("time", "treatment"), 
           returnData = TRUE) %>%
  ggplot(aes(x = time, y = count, color = treatment, group = treatment)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) +
  stat_summary(fun = mean, geom = "line", aes(group = treatment), size = 1) +
  labs(title = paste("Expresión de:", top_interaction_gene),
       subtitle = "Gen con mayor significancia en la interacción",
       x = "Time",
       y = "Normalized Counts") +
  theme_bw() +
  scale_y_log10()

message("Script de análisis completado.")

############################################################
################################################################################
###                                                                          ###
### SCRIPT FINAL DEFINITIVO Y COMPLETO: Análisis de RNA-seq en Cardiomiocitos  ###
###                                                                          ###
###   Pipeline de principio a fin con Anotación Curada y todas las             ###
###          visualizaciones con calidad de publicación en inglés.           ###
###                                                                          ###
################################################################################

# --- Step 0: Load Libraries ---
message("--- Step 0: Loading libraries ---")
suppressPackageStartupMessages({
  library(tximport); library(rtracklayer); library(DESeq2); library(apeglm)
  library(readr); library(dplyr); library(tidyverse); library(ggplot2)
  library(pheatmap); library(RColorBrewer); library(ggrepel); library(AnnotationDbi)
  library(org.Rn.eg.db); library(biomaRt); library(Hmisc)
})


# --- Step 1: Define Paths and Plotting Theme ---
message("--- Step 1: Configuring paths and plotting theme ---")
# --- ADJUST THESE PATHS TO YOUR SYSTEM!!! ---
gtf_file <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis/merged_CNE_6_24_HISAT2.annotated.gtf"
lncrna_ann_all_file <- "hisat2_6_24_meta_ann.csv"
lncrna_ann_novel_file <- "hisat2_6_24_meta_novel_lncRNA.csv"
base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/salmon_quants_NE624/"
output_dir <- "Publication_Ready_Analysis_Final_Definitive"
# --- End of path section ---

# Create output directories
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "1_DE_Tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "2_QC_and_Global_Viz"), showWarnings = FALSE)
dir.create(file.path(output_dir, "3_DE_Visualization"), showWarnings = FALSE)
dir.create(file.path(output_dir, "4_Expression_Patterns"), showWarnings = FALSE)
dir.create(file.path(output_dir, "5_Gene_Lists_for_Enrichr"), showWarnings = FALSE)
dir.create(file.path(output_dir, "6_Coexpression_Analysis"), showWarnings = FALSE)
dir.create(file.path(output_dir, "6_Coexpression_Analysis/Top_Correlation_Plots"), showWarnings = FALSE)

# ggplot2 theme for publication-quality aesthetics
theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      legend.title = element_blank()
    )
}

# --- Steps 2-5: Main Pipeline (Import, DESeq2 & Curated Annotation) ---
message("--- Steps 2-5: Running main pipeline (Import, DESeq2, Annotation) ---")
ann_all <- read_csv(lncrna_ann_all_file, col_types = cols(.default = "c")); ann_novel <- read_csv(lncrna_ann_novel_file, col_types = cols(.default = "c"))
gtf_data <- rtracklayer::import(gtf_file); gtf_df <- as.data.frame(gtf_data)
tx2gene <- gtf_df %>% dplyr::filter(type == "transcript") %>% dplyr::select(transcript_id, gene_id) %>% distinct() %>% na.omit()
samples <- c("Ctrl6_rep1", "Ctrl6_rep2", "Ctrl6_rep4", "NE6_rep1", "NE6_rep2", "NE6_rep4", "Ctrl24_rep1", "Ctrl24_rep2", "Ctrl24_rep4", "NE24_rep1", "NE24_rep2", "NE24_rep4")
sampleTable <- data.frame(treatment = factor(rep(c("Control", "NE"), each = 3, times = 2)), time = factor(rep(c("6h", "24h"), each = 6), levels=c("6h", "24h")))
rownames(sampleTable) <- samples; sampleTable$treatment <- relevel(sampleTable$treatment, ref = "Control")
files <- file.path(base_dir, paste0(samples, "_quant"), "quant.sf"); names(files) <- samples
txi <- tximport(files, type = "salmon", tx2gene = tx2gene); dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ time + treatment + time:treatment)
dds <- dds[rowSums(counts(dds)) >= 10, ]; dds <- DESeq(dds); vsd <- vst(dds, blind = FALSE)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ time + treatment); res_lrt <- results(dds_lrt)
genes_in_analysis <- data.frame(gene_id = rownames(dds))
annotation_from_gtf <- gtf_df %>% dplyr::filter(type == "transcript") %>% dplyr::select(gene_id, gene_name) %>% dplyr::filter(!is.na(gene_name) & gene_name != "") %>% distinct(gene_id, .keep_all = TRUE)
full_annotation <- left_join(genes_in_analysis, annotation_from_gtf, by = "gene_id")
curate_mstrg_annotation_with_biomart <- function(annotation_table) {
  genes_to_curate <- annotation_table %>% dplyr::filter(startsWith(gene_id, "MSTRG") & !is.na(gene_name) & !startsWith(gene_name, "MSTRG"))
  if (nrow(genes_to_curate) == 0) return(annotation_table)
  message(paste("... Curating", nrow(genes_to_curate), "MSTRG genes with known names via biomaRt..."))
  ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl")
  biomart_results <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'gene_biotype'), filters = 'external_gene_name', values = unique(genes_to_curate$gene_name), mart = ensembl)
  if (nrow(biomart_results) > 0) {
    biomart_results_clean <- biomart_results %>% rename(gene_name = external_gene_name, official_ensembl_id = ensembl_gene_id, official_gene_biotype = gene_biotype) %>% distinct(gene_name, .keep_all = TRUE)
    annotation_table <- left_join(annotation_table, biomart_results_clean, by = "gene_name")
  }
  return(annotation_table)
}
full_annotation <- curate_mstrg_annotation_with_biomart(full_annotation)
if (!"official_gene_biotype" %in% names(full_annotation)) { full_annotation$official_gene_biotype <- NA }
if (!"official_ensembl_id" %in% names(full_annotation)) { full_annotation$official_ensembl_id <- NA }
gene_ids_no_version_all <- gsub("\\..*$", "", full_annotation$gene_id); entrez_map <- AnnotationDbi::select(org.Rn.eg.db, keys = gene_ids_no_version_all, columns = c("ENTREZID"), keytype = "ENSEMBL", multiVals = "first") %>% rename(gene_id_no_version = ENSEMBL) %>% distinct(gene_id_no_version, .keep_all = TRUE)
full_annotation$gene_id_no_version <- gene_ids_no_version_all; full_annotation <- left_join(full_annotation, entrez_map, by = "gene_id_no_version")
known_lncRNAs <- ann_all %>% dplyr::filter(gene_biotype == "lncRNA") %>% dplyr::select(ensembl_gene_id) %>% dplyr::distinct(); novel_lncRNAs <- ann_novel %>% dplyr::select(gene_id) %>% dplyr::distinct() %>% dplyr::rename(ensembl_gene_id = gene_id)
all_lncRNAs_ids <- bind_rows(known_lncRNAs, novel_lncRNAs) %>% distinct(); protein_coding_ids <- ann_all %>% dplyr::filter(gene_biotype == "protein_coding") %>% dplyr::select(ensembl_gene_id) %>% distinct()
full_annotation <- full_annotation %>%
  mutate(biotype_clean = case_when(!is.na(official_gene_biotype) ~ official_gene_biotype, gene_id %in% all_lncRNAs_ids$ensembl_gene_id ~ "lncRNA", gene_id %in% protein_coding_ids$ensembl_gene_id ~ "protein_coding", startsWith(gene_id, "ENSRNOG") ~ "protein_coding", TRUE ~ "other"),
         SYMBOL = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)) %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, biotype_clean, official_ensembl_id, official_gene_biotype)
write.csv(full_annotation, file.path(output_dir, "curated_full_annotation.csv"), row.names = FALSE)

# --- Step 6: Process and Save Differential Expression Results ---
message("--- Step 6: Processing and saving DE results ---")
add_status_column <- function(res_df, annotation_df) {
  res_df %>% as.data.frame() %>% rownames_to_column("gene_id") %>% left_join(annotation_df, by = "gene_id") %>%
    mutate(status = case_when(padj < 0.05 & abs(log2FoldChange) >= 1 & biotype_clean == "protein_coding" ~ if_else(log2FoldChange > 0, "Upregulated", "Downregulated"), padj < 0.05 & abs(log2FoldChange) >= 0.5 & biotype_clean == "lncRNA" ~ if_else(log2FoldChange > 0, "Upregulated", "Downregulated"), TRUE ~ "Not Significant"))
}
res_NE6_processed <- add_status_column(results(dds, name = "treatment_NE_vs_Control"), full_annotation)
res_NE24_processed <- add_status_column(results(dds, contrast = list(c("treatment_NE_vs_Control", "time24h.treatmentNE"))), full_annotation)
res_lrt_processed <- as.data.frame(res_lrt) %>% rownames_to_column("gene_id") %>% left_join(full_annotation, by="gene_id")
write.csv(res_NE6_processed, file.path(output_dir, "1_DE_Tables/DE_results_6h_vs_control.csv"), row.names = FALSE)
write.csv(res_NE24_processed, file.path(output_dir, "1_DE_Tables/DE_results_24h_vs_control.csv"), row.names = FALSE)


# --- Step 7: Export Gene Lists for External Functional Analysis ---
message("--- Step 7: Exporting gene lists for external tools (e.g., Enrichr) ---")
export_gene_lists <- function(res_processed, timepoint) {
  significant_genes <- res_processed %>% dplyr::filter(status %in% c("Upregulated", "Downregulated"))
  up_symbols <- significant_genes %>% dplyr::filter(status == "Upregulated") %>% pull(SYMBOL); write.table(up_symbols, file.path(output_dir, "5_Gene_Lists_for_Enrichr", paste0("list_upregulated_symbols_", timepoint, ".txt")), row.names=F, col.names=F, quote=F)
  down_symbols <- significant_genes %>% dplyr::filter(status == "Downregulated") %>% pull(SYMBOL); write.table(down_symbols, file.path(output_dir, "5_Gene_Lists_for_Enrichr", paste0("list_downregulated_symbols_", timepoint, ".txt")), row.names=F, col.names=F, quote=F)
  message(paste("... Gene lists for", timepoint, "saved in:", file.path(output_dir, "5_Gene_Lists_for_Enrichr")))
}
export_gene_lists(res_NE6_processed, "6h"); export_gene_lists(res_NE24_processed, "24h")

# --- Step 8: QC and Global Distribution Visualization ---
message("--- Step 8: Generating QC and global distribution plots ---")
pca_plot <- plotPCA(vsd, intgroup=c("treatment", "time")) + geom_text_repel(aes(label=name)) + labs(title="Principal Component Analysis (PCA)") + theme_publication()
ggsave(file.path(output_dir, "2_QC_and_Global_Viz/PCA_plot.png"), plot = pca_plot, width = 8, height = 6)
vsd_df <- as.data.frame(assay(vsd)) %>% rownames_to_column("gene_id"); vsd_long <- vsd_df %>% pivot_longer(cols = -gene_id, names_to = "sample", values_to = "expression")
plot_data_violin <- vsd_long %>% left_join(full_annotation %>% dplyr::select(gene_id, biotype_clean), by = "gene_id") %>% left_join(sampleTable %>% rownames_to_column("sample"), by = "sample") %>% dplyr::filter(biotype_clean %in% c("protein_coding", "lncRNA")) %>% mutate(group = paste(treatment, time, sep = "_"))
expression_violin_plot <- ggplot(plot_data_violin, aes(x = group, y = expression, fill = group)) + geom_violin(trim = FALSE) + facet_wrap(~ biotype_clean, scales = "free_y") + stat_summary(fun = median, geom = "point", shape = 18, size = 4, color = "white") + labs(title = "Gene Expression Distribution", x = "Experimental Condition", y = "Expression Level (VST)") + theme_publication(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(output_dir, "2_QC_and_Global_Viz/Expression_Distribution_Violin_Plot.png"), plot = expression_violin_plot, width = 12, height = 7)

# --- Step 9: DEGs Summary and Volcano Plot Visualization ---
message("--- Step 9: Creating DEG summary and Volcano plots ---")
# 9.1 Stacked Bar Plot of DEG Counts
counts_6h <- res_NE6_processed %>% dplyr::filter(status != "Not Significant", biotype_clean %in% c("protein_coding", "lncRNA")) %>% group_by(biotype_clean, status) %>% summarise(count = n(), .groups = 'drop') %>% mutate(timepoint = "6h")
counts_24h <- res_NE24_processed %>% dplyr::filter(status != "Not Significant", biotype_clean %in% c("protein_coding", "lncRNA")) %>% group_by(biotype_clean, status) %>% summarise(count = n(), .groups = 'drop') %>% mutate(timepoint = "24h")
summary_counts <- bind_rows(counts_6h, counts_24h) %>% mutate(timepoint = factor(timepoint, levels=c("6h", "24h")), status = factor(status, levels = c("Upregulated", "Downregulated")))
total_counts <- summary_counts %>% group_by(timepoint, biotype_clean) %>% summarise(total = sum(count), .groups = 'drop')
deg_summary_plot <- ggplot(summary_counts, aes(x = timepoint, y = count)) +
  geom_bar(aes(fill = status), stat = "identity", position = "stack") +
  geom_text(data = total_counts, aes(x = timepoint, y = total, label = total), vjust = -0.5, size = 5) +
  facet_wrap(~ biotype_clean, scales = "free_y") +
  scale_fill_manual(values = c("Upregulated" = "#e41a1c", "Downregulated" = "#377eb8")) +
  labs(title = "Number of Differentially Expressed Genes", x = "Timepoint", y = "Number of Genes") +
  theme_publication(base_size = 16) + scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave(file.path(output_dir, "3_DE_Visualization/DEG_Summary_Stacked_Bar_Plot.png"), plot = deg_summary_plot, width = 10, height = 7)

# 9.2 Volcano Plots
create_volcano_plot <- function(res_data, biotype_filter, title) {
  plot_data <- res_data %>% dplyr::filter(biotype_clean == biotype_filter, !is.na(padj)) %>% mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj), label = ifelse(status != "Not Significant" & (abs(log2FoldChange) > 2 | padj < 1e-15), SYMBOL, ""))
  ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = status, label = label)) + geom_point(alpha = 0.7, size = 2) + geom_text_repel(max.overlaps = 15, size = 3.5, box.padding = 0.5, min.segment.length = 0, seed = 42) + scale_color_manual(values = c("Upregulated" = "#e41a1c", "Downregulated" = "#377eb8", "Not Significant" = "grey80")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed") + labs(title = title, x = "Log2 Fold Change", y = "-log10(p-adj)") + theme_publication() + theme(legend.position = "top")
}
ggsave(file.path(output_dir, "3_DE_Visualization/Volcano_6h_protein_coding.png"), create_volcano_plot(res_NE6_processed, "protein_coding", "6h: Protein-Coding Genes"), width=8, height=7)
ggsave(file.path(output_dir, "3_DE_Visualization/Volcano_6h_lncRNA.png"), create_volcano_plot(res_NE6_processed, "lncRNA", "6h: lncRNAs"), width=8, height=7)
ggsave(file.path(output_dir, "3_DE_Visualization/Volcano_24h_protein_coding.png"), create_volcano_plot(res_NE24_processed, "protein_coding", "24h: Protein-Coding Genes"), width=8, height=7)
ggsave(file.path(output_dir, "3_DE_Visualization/Volcano_24h_lncRNA.png"), create_volcano_plot(res_NE24_processed, "lncRNA", "24h: lncRNAs"), width=8, height=7)

# --- Step 10: Expression Pattern Visualization ---
message("--- Step 10: Generating expression pattern plots ---")
# 10.1 Heatmap
top_kinetic_genes <- res_lrt_processed %>% dplyr::filter(padj < 0.05) %>% arrange(padj) %>% head(50)
if (nrow(top_kinetic_genes) > 0) {
  mat <- assay(vsd)[top_kinetic_genes$gene_id, ]; rownames(mat) <- top_kinetic_genes$SYMBOL
  pheatmap(mat, scale = "row", annotation_col = as.data.frame(colData(dds)[,c("treatment", "time")]), main = "Top 50 Genes with Dynamic Patterns (LRT)", border_color = "grey60", filename = file.path(output_dir, "4_Expression_Patterns/heatmap_top50_kinetic_genes.png"), width=8, height=12)
}
# 10.2 Individual Gene Plots
plot_gene <- function(gene_symbol, dds_obj, annotation_df) {
  gene_id <- annotation_df$gene_id[which(annotation_df$SYMBOL == gene_symbol)]; if(length(gene_id)==0) return(NULL)
  p <- plotCounts(dds_obj, gene = gene_id[1], intgroup = c("time", "treatment"), returnData = TRUE) %>%
    ggplot(aes(x = time, y = count, color = treatment, group = treatment)) + geom_point(position = position_jitter(width = 0.1, height = 0), size = 3, alpha=0.7) + stat_summary(fun = mean, geom = "line", linewidth = 1.2) + labs(title = gene_symbol, x = "Timepoint", y = "Normalized Counts") + scale_y_log10() + theme_publication()
  return(p)
}
res_interaction_processed <- add_status_column(results(dds, name="time24h.treatmentNE"), full_annotation)
top_genes_interaction <- res_interaction_processed %>% arrange(padj) %>% na.omit()
top_pc_genes_plot <- top_genes_interaction %>% dplyr::filter(biotype_clean == "protein_coding") %>% head(6)
for (gene in top_pc_genes_plot$SYMBOL) {
  p <- plot_gene(gene, dds, full_annotation); if (!is.null(p)) { ggsave(file.path(output_dir, "4_Expression_Patterns/Genes_Individuales", paste0("PC_", gsub("[/:]", "_", gene), ".png")), plot = p, width = 6, height = 5) }
}

# --- Step 11: Rigorous Co-expression Analysis with p-value ---
message("--- Step 11: Running rigorous co-expression analysis with p-value ---")
run_correlation_analysis <- function(res_processed, vsd_matrix, timepoint, full_annot_table) {
  de_lncRNAs <- res_processed %>% dplyr::filter(biotype_clean == "lncRNA", status != "Not Significant"); de_mRNAs <- res_processed %>% dplyr::filter(biotype_clean == "protein_coding", status != "Not Significant")
  if(nrow(de_lncRNAs) < 2 | nrow(de_mRNAs) < 2) { message(paste("... Skipping", timepoint, "correlation: not enough DE lncRNAs or mRNAs.")); return(NULL) }
  lncRNA_counts <- assay(vsd_matrix)[de_lncRNAs$gene_id, ]; mRNA_counts <- assay(vsd_matrix)[de_mRNAs$gene_id, ]
  cor_results <- Hmisc::rcorr(t(lncRNA_counts), t(mRNA_counts), type = "pearson"); cor_matrix <- cor_results$r; p_matrix <- cor_results$P
  sig_indices <- which(abs(cor_matrix) > 0.8 & p_matrix < 0.05, arr.ind = TRUE); if (nrow(sig_indices) == 0) { message(paste("... No significant correlations found at", timepoint)); return(NULL) }
  sig_cor_df <- data.frame(lncRNA_id=rownames(cor_matrix)[sig_indices[,1]], mRNA_id=colnames(cor_matrix)[sig_indices[,2]], correlation=cor_matrix[sig_indices], p_value=p_matrix[sig_indices])
  sig_cor_df <- sig_cor_df %>%
    left_join(full_annot_table %>% dplyr::select(gene_id, SYMBOL), by = c("lncRNA_id" = "gene_id")) %>% rename(lncRNA_symbol = SYMBOL) %>%
    left_join(full_annot_table %>% dplyr::select(gene_id, SYMBOL), by = c("mRNA_id" = "gene_id")) %>% rename(mRNA_symbol = SYMBOL) %>%
    dplyr::select(lncRNA_id, lncRNA_symbol, mRNA_id, mRNA_symbol, correlation, p_value)
  write.csv(sig_cor_df, file.path(output_dir, "6_Coexpression_Analysis", paste0("significant_correlations_", timepoint, ".csv")), row.names = FALSE)
  return(sig_cor_df)
}
correlations_6h <- run_correlation_analysis(res_NE6_processed, vsd, "6h", full_annotation)
correlations_24h <- run_correlation_analysis(res_NE24_processed, vsd, "24h", full_annotation)

# --- Step 12: Hub Identification and Correlation Visualization ---
message("--- Step 12: Identifying Hubs and visualizing top correlations ---")
all_significant_correlations <- bind_rows(correlations_6h, correlations_24h) %>% distinct(lncRNA_id, mRNA_id, .keep_all = TRUE)
if(!is.null(all_significant_correlations) && nrow(all_significant_correlations) > 0) {
  # 12.1 Generate Scatter Plots for the Top 5 Correlated Pairs
  top_n_pairs <- 5
  top_pairs <- all_significant_correlations %>% arrange(p_value, desc(abs(correlation))) %>% head(top_n_pairs)
  
  message(paste("... Generating correlation scatter plots for the top", top_n_pairs, "lncRNA-mRNA pairs..."))
  for(i in 1:nrow(top_pairs)) {
    pair <- top_pairs[i, ]
    scatter_data <- as.data.frame(assay(vsd)[c(pair$lncRNA_id, pair$mRNA_id), ]) %>% t() %>% as.data.frame() %>% `colnames<-`(c("lncRNA_expr", "mRNA_expr")) %>% rownames_to_column("sample") %>% left_join(sampleTable %>% rownames_to_column("sample"), by="sample")
    correlation_scatter_plot <- ggplot(scatter_data, aes(x=lncRNA_expr, y=mRNA_expr)) + 
      geom_point(aes(color = treatment), size = 4) + 
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) + 
      facet_wrap(~ time) + 
      scale_color_manual(values = c("Control" = "grey50", "NE" = "#e41a1c")) + 
      labs(title = paste("Correlation:", pair$lncRNA_symbol, "&", pair$mRNA_symbol), subtitle = paste("Overall r =", round(pair$correlation, 2), "| p-value =", format.pval(pair$p_value, digits=2)), x = paste("Expression of", pair$lncRNA_symbol), y = paste("Expression of", pair$mRNA_symbol)) + 
      theme_publication(base_size=16)
    ggsave(file.path(output_dir, "6_Coexpression_Analysis/Top_Correlation_Plots", paste0("Correlation_Plot_", i, "_", pair$lncRNA_symbol, "_", pair$mRNA_symbol, ".png")), plot = correlation_scatter_plot, width = 12, height = 7)
  }
  
  # 12.2 Identify and Generate Detailed Table for lncRNA Hubs
  hub_summary <- all_significant_correlations %>% group_by(lncRNA_id, lncRNA_symbol) %>% summarise(n_correlated_mRNAs = n(), .groups = 'drop') %>% arrange(desc(n_correlated_mRNAs))
  hub_threshold <- 10
  lncrna_hubs <- hub_summary %>% dplyr::filter(n_correlated_mRNAs >= hub_threshold)
  if(nrow(lncrna_hubs) > 0) {
    message(paste("... ", nrow(lncrna_hubs), "lncRNA hubs identified. Generating detailed interaction table..."))
    hub_interaction_table <- all_significant_correlations %>%
      dplyr::filter(lncRNA_id %in% lncrna_hubs$lncRNA_id) %>%
      left_join(lncrna_hubs, by=c("lncRNA_id", "lncRNA_symbol")) %>%
      arrange(desc(n_correlated_mRNAs), lncRNA_symbol, desc(abs(correlation))) %>%
      dplyr::select(
        lncRNA_hub_Symbol = lncRNA_symbol,
        Correlated_mRNA_Symbol = mRNA_symbol,
        Correlation = correlation,
        P_value = p_value,
        lncRNA_hub_ID = lncRNA_id,
        Correlated_mRNA_ID = mRNA_id
      )
    write.csv(hub_interaction_table, file.path(output_dir, "6_Coexpression_Analysis/lncRNA_hub_to_mRNA_target_list.csv"), row.names = FALSE)
  } else { message("... No lncRNA hubs were found with the current thresholds.") }
} else { message("No significant correlations were found to build the hub network.") }

message("\n¡ANALYSIS COMPLETE! Check the output folder '", output_dir, "'.")
