# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(circlize)
library(DOSE)
library(enrichplot)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(EnhancedVolcano)
library(sva)
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(RColorBrewer)
library(EnhancedVolcano)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Rn.eg.db) # Annotation database for Rattus norvegicus
library(rtracklayer)
library(ggrepel) # For non-overlapping labels in plots
library(apeglm)             # for lfcShrink
library(tibble)

countdata = read.csv("miR.Counts_C24_NE24_3v3_new.csv", header=TRUE, row.names=1)
countdata = countdata[c(6:13)]  #selecciona columnas a analizar
countdata = countdata[, -3]

# Convert to matrix and add 1 to rows that have 0's
countdata = as.matrix(countdata)
countdata = countdata + 1
# Eliminar filas con valores NA
countdata_clean <- na.omit(countdata)

# Assign condition (affected versus unaffected), you can assign any name, just keep the number of conditions
condition = factor(c("Ctrl_24", "Ctrl_24", "Ctrl_24", "NE_24", "NE_24", "NE_24"),
                   levels=c("NE_24", "Ctrl_24"))
condition = relevel(condition, ref = "Ctrl_24")

condition = factor(c("Ctrl6", "Ctrl6", "Ctrl6", "NE_6", "NE_6", "NE_6"),
                   levels=c("NE_6", "Ctrl6"))
condition = relevel(condition, ref = "Ctrl6")

# Agregar información de batch (4 batches distintos)
batch = factor(c("Batch1", "Batch2", "Batch3", "Batch1", "Batch2", "Batch3"))

# Crear coldata con condición y batch
coldata = data.frame(row.names = colnames(countdata), condition, batch)
coldata = data.frame(row.names = colnames(countdata_clean), condition, batch)
coldata = data.frame(row.names = colnames(countdata), condition)
coldata

combat_counts <- ComBat_seq(counts = countdata_clean, batch = coldata$batch, group = condition)

# Crear DESeq2 dataset usando los datos corregidos por Combat-Seq
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds = DESeqDataSetFromMatrix(countData = combat_counts, colData = coldata, design = ~ condition)

# Filtrar filas con bajo conteo
keep <- rowSums(counts(dds) >= 10)
dds <- dds[keep,]
dds <- DESeq(dds)

vsdata <- vst(dds, blind = FALSE)

#HEATMAP
sampleDists <- dist(t(assay(vsd_mi)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata$condition, vsdata$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         fontsize = 10)

plotPCA(vsdata, intgroup = 'condition')

res <- results(dds)
res
sigs <- na.omit(res)
sigs

resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized=T)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Geneid"
head(resdata)

# Extract normalized counts
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
names(normalized_counts)[-1] <- "Geneid"


# Add Geneid column (rownames are gene IDs)
normalized_counts <- cbind(Geneid = rownames(normalized_counts), normalized_counts)

# Save as TXT (tab-separated)
write.table(normalized_counts, file = "miRNAs_c6_ne6_3v3_normalized_10.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Save as CSV (comma-separated)
write.csv(normalized_counts, file = "normalized_counts_10_NE6_3x3.csv", row.names = FALSE, quote = FALSE)

# Regularized log transformation for clustering/heatmaps
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots below
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists = as.matrix(dist(t(assay(rld))))

library(gplots)

#png("qc-heatmap_baker.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(8, 8), main="Distance Matrix C6-NE6 3vs3 miRNAs")

#packages that  need to be loaded  to make pca plot
library(MASS)
library(genefilter)
library(calibrate)

# Principal components analysis (PCA)
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot C24 vs NE24 miRNAs SM", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

png("C24_NE24_3v3_miRNAs_SM_batch_counts_new.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-30, 30))
dev.off()

# add a column of NAs
resdata$diffexpressed <- "NO"

resdata <- normalized_counts$diffexpressed <- "NO"

# if log2Foldchange >= 0.5 and pvalue < 0.05, set as "UP" 
resdata$diffexpressed[resdata$log2FoldChange >= 1 & resdata$pvalue < 0.05] <- "Upregulated"
resdata$diffexpressed[resdata$log2FoldChange <= -1 & resdata$pvalue < 0.05] <- "Downregulated"

resdata$diffexpressed[resdata$log2FoldChange >= 0.5 & resdata$padj < 0.05] <- "Upregulated"
resdata$diffexpressed[resdata$log2FoldChange <= -0.5 & resdata$padj < 0.05] <- "Downregulated"

normalized_counts$diffexpressed[normalized_counts$log2FoldChange >= 0.5 & normalized_counts$padj < 0.05] <- "Upregulated"
normalized_counts$diffexpressed[normalized_counts$log2FoldChange <= -0.5 & normalized_counts$padj < 0.05] <- "Downregulated"

write.table(resdata, file = "Deseq2_C6_vs_NE6_3v3_miRNAs_norm_new.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

miRNAs_up = subset(resdata, resdata$padj < 0.05 & resdata$log2FoldChange >= 1)
miRNAs_down = subset(resdata, resdata$padj < 0.05 & resdata$log2FoldChange <= -1)

write.table(miRNAs_down, file = "miRNAs_down_c6_ne6_norm_3v3.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

resdata2 <- read.delim("Deseq2_C6_vs_NE6_3v3_miRNAs_norm_new.txt", header = TRUE, sep = "\t")

# add a column of NAs
resdata2$diffexpressed <- "NO"

# if log2Foldchange >= 0.5 and pvalue < 0.05, set as "UP" 
resdata2$diffexpressed[resdata2$log2FoldChange >= 1 & resdata2$pvalue < 0.05] <- "Upregulated"
resdata2$diffexpressed[resdata2$log2FoldChange <= -1 & resdata2$pvalue < 0.05] <- "Downregulated"

#We can select some specific genes to the volcano plot
keyvals <- ifelse(
  resdata$log2FoldChange <= -1 & resdata$padj < 0.05, 'royalblue',
  ifelse(resdata$log2FoldChange >= 1 & resdata$padj < 0.05, 'red',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey'] <- 'NA'
names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'

EnhancedVolcano(resdata, lab = resdata$Geneid,
                x = 'log2FoldChange', y = 'pvalue',
                selectLab = rownames(resdata$Row.names)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                xlab = bquote(~Log[2]~ 'Fold Change'),
                ylab = bquote(~-log[10]~ 'padj'),
                title = '', pCutoff = 0.05, FCcutoff = 1, pointSize = 4,
                labSize = 5, colCustom = keyvals,
                #colAlpha = 1,
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 4.0,
                drawConnectors = F, widthConnectors = 1.0, colConnectors = 'black',
                arrowheads = F, gridlines.major = F, gridlines.minor = F,
                border = 'partial', borderWidth = 1,
                borderColour = 'black',
                hline = NULL, vline = NULL)

# Add a column for labeling only the most significant genes on the plot
final_res <- resdata %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >= 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    # Create a label column: show gene_name only for significant genes
    gene_label = if_else(padj < 0.05 & abs(log2FoldChange) > 1, Geneid, "")
  )

# Create the volcano plot
volcano_plot <- ggplot(final_res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  geom_text_repel(aes(label = gene_label), max.overlaps = 20, size = 3.5) + # Add labels
  scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "grey")) +
  theme_bw(base_size = 14) +
  labs(
    title = "",
    x = "Log2FoldChange",
    y = "-Log10(padj)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

print(volcano_plot)
ggsave("Volcano_NE24_miRNAs.png", plot = volcano_plot, width = 8, height = 7)

library(readxl)
library(dplyr)

#####Nuevo script con miRNAs en columna
# 1. Leer los datos
mirwalk <- read.csv("miRWalk_miRNA_Targets_15up_CDS.csv", header = TRUE)
deseq <- read.table("MBS4_annonated_prot_C24_NE24_new_normalized_3v3.txt", header = TRUE, sep = "\t")

# 2. Filtrar DESeq2: genes upregulated y protein-coding
deseq_up <- deseq[deseq$log2FoldChange <= -1 & deseq$padj < 0.05 & deseq$gene_biotype == "protein_coding", ]

# 3. Extraer genes únicos de miRWalk
mirwalk_targets <- unique(mirwalk$genesymbol)

# 4. Intersección
genes_comunes <- intersect(deseq_up$external_gene_name, mirwalk_targets)

# 5. Subconjunto final
genes_final <- deseq_up[deseq_up$external_gene_name %in% genes_comunes, ]

# 6. Obtener miRNAs correspondientes a cada gen
library(dplyr)

mirna_por_gen <- mirwalk %>%
  filter(genesymbol %in% genes_comunes) %>%
  group_by(genesymbol) %>%
  summarise(miRNAs = paste(unique(mirnaid), collapse = "; "))

# 7. Unir los miRNAs a genes_final
genes_final_con_miRNAs <- merge(genes_final, mirna_por_gen, by.x = "external_gene_name", by.y = "genesymbol", all.x = TRUE)

# 8. Guardar resultado final
write.table(genes_final_con_miRNAs, "Genes_downregulated_targets_miRWalk_CDS_15up_miRNAs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Para integrar los tres sitios donde posiblemente se une
# Leer archivos e identificar región
mirwalk_3utr <- read.csv("miRWalk_miRNA_Targets_15up_3UTR.csv", header = TRUE)
mirwalk_3utr$Region <- "3UTR"

mirwalk_5utr <- read.csv("miRWalk_miRNA_Targets_15up_5UTR.csv", header = TRUE)
mirwalk_5utr$Region <- "5UTR"

mirwalk_cds <- read.csv("miRWalk_miRNA_Targets_15up_CDS.csv", header = TRUE)
mirwalk_cds$Region <- "CDS"

# Unir todos en uno solo
mirwalk_all <- rbind(mirwalk_3utr, mirwalk_5utr, mirwalk_cds)

# Agrupar por miRNA y genesymbol, y colapsar regiones
library(dplyr)

mirwalk_consolidado <- mirwalk_all %>%
  group_by(mirnaid, genesymbol) %>%
  summarise(regiones = paste(unique(Region), collapse = "; "), .groups = "drop")

# Ver ejemplo
head(mirwalk_consolidado)

library(dplyr)
library(tidyr)

# Separar cada miRNA en una fila distinta
genes_final_expandido <- genes_final_con_miRNAs %>%
  separate_rows(miRNAs, sep = ";\\s*")  # separa por "; "

# Ahora hacer el merge correctamente
genes_final_con_regiones <- merge(genes_final_expandido, mirwalk_consolidado, 
                                  by.x = c("external_gene_name", "miRNAs"),
                                  by.y = c("genesymbol", "mirnaid"),
                                  all.x = TRUE)

# Opcional: volver a agrupar los miRNAs y regiones por gen
genes_final_regrouped <- genes_final_con_regiones %>%
  group_by(external_gene_name) %>%
  summarise(across(where(is.character), ~ paste(unique(.), collapse = "; ")), .groups = "drop")

##############################################Nuevo intento, integrando todo
# Cargar paquetes
library(dplyr)
library(tidyr)

# 1. Leer resultados de miRWalk para cada región
mirwalk_3utr <- read.csv("miRWalk_miRNA_Targets_15up_3UTR.csv", header = TRUE)
mirwalk_cds  <- read.csv("miRWalk_miRNA_Targets_15up_CDS.csv", header = TRUE)
mirwalk_5utr <- read.csv("miRWalk_miRNA_Targets_15up_5UTR.csv", header = TRUE)

# 2. Agregar columna de región a cada uno
mirwalk_3utr$region <- "3UTR"
mirwalk_cds$region  <- "CDS"
mirwalk_5utr$region <- "5UTR"

# 3. Unir todos los resultados
mirwalk_combined <- bind_rows(mirwalk_3utr, mirwalk_cds, mirwalk_5utr)

# 4. Filtrar columnas relevantes y eliminar duplicados
mirwalk_filtered <- mirwalk_combined %>%
  select(mirnaid, genesymbol, region) %>%
  distinct()

# 5. Agrupar por miRNA y gene, colapsar regiones múltiples
mirwalk_grouped <- mirwalk_filtered %>%
  group_by(genesymbol, mirnaid) %>%
  summarise(region = paste(sort(unique(region)), collapse = ";"), .groups = "drop")

# 6. Leer archivo DESeq2 anotado
deseq <- read.table("MBS4_annotated_all_C24_NE24_new_normalized_3v3_10.txt", header = TRUE, sep = "\t")

# 7. Filtrar genes upregulated y protein-coding
deseq_up <- deseq %>%
  filter(log2FoldChange <= -1, padj < 0.05, gene_biotype == "protein_coding")

# 8. Cruzar con genes en común
mirwalk_final <- mirwalk_grouped %>%
  filter(genesymbol %in% deseq_up$external_gene_name)

# 9. Unir con datos de DESeq2
final_result <- merge(mirwalk_final, deseq_up, 
                      by.x = "genesymbol", by.y = "external_gene_name")

# 10. Reordenar columnas
final_result <- final_result %>%
  select(genesymbol, mirnaid, region, everything())

# 11. Guardar resultado final
write.table(final_result, "miRNAup_collapsed_regions_DESeq2_filtered_mRNAdown_NE24_10counts.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

library(readxl)
library(dplyr)

# Leer archivo final de miRWalk con DESeq2
mirwalk_final <- read.table("miRNA_gene_collapsed_regions_DESeq2_filtered.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Leer archivo Excel de miRDB
mirDB <- read_excel("miRDB target prediction data_4miRNAs_correlation.xls.xlsx")  # Cambia el nombre si es necesario

colnames(mirDB)

# Elegir solo columnas relevantes y eliminar duplicados
mirDB_clean <- mirDB %>%
  select(`miRNA Name`, `Gene Symbol`) %>%
  distinct()

# Cruzar usando tanto el nombre del gen como del miRNA
validacion_mirdb <- final_result %>%
  inner_join(mirDB_clean, 
             by = c("genesymbol" = "Gene Symbol", "mirnaid" = "miRNA Name"))

write.table(validacion_mirdb, "miRNA_gene_validated_by_miRWalk_and_miRDB.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Nuevo script lncRNA-miRNA-mRNAs
# 1. Leer archivos
mirwalk <- read.table("miRNA_gene_collapsed_regions_DESeq2_filtered_downmiRNAs_NE24.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

correlaciones <- read.table("correlacion_lnc_up_mrna_significativos.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Filtrar correlaciones significativas y positivas
cor_sig <- correlaciones[correlaciones$rho > 0.8 & correlaciones$p.value < 0.05, ]

# 3. Extraer genes comunes
genes_comunes <- intersect(mirwalk$Geneid, cor_sig$mRNA)

# 4. Filtrar ambos datasets
mirwalk_filtrado <- mirwalk[mirwalk$Geneid %in% genes_comunes, ]
cor_filtrado <- cor_sig[cor_sig$mRNA %in% genes_comunes, ]

# 5. Unir datos por mRNA
merged <- merge(mirwalk_filtrado, cor_filtrado, 
                by.x = "Geneid", by.y = "mRNA")

# 7. Guardar archivo final
write.table(merged, "lncRNA_miRNA_mRNA_ceRNA_filtered_final_4downmiRNAs_lncRNAs_up.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

####Acá se debería filtrar con las correlaciones significativas de lncRNA-miRNA significativa
###Para luego ver las que sí son posiblemente validables por interacción
######Nuevo script integrativo
# 1. Leer archivos existentes
#Interacción miRNA-mRNA
mirwalk <- read.table("miRNA_gene_collapsed_regions_DESeq2_filtered_downmiRNAs_NE24.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Correlación lncRNA-mRNA, debe ser lncRNAup-mRNAup
correlaciones_lnc_mrna <- read.table("correlacion_lnc_up_mrna_significativos.txt", 
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Correlación lncRNA-miRNA
correlaciones_lnc_miRNA <- read.table("correlacion_miRNA_vs_lncRNA_all_DE.txt",  # cambia el nombre si es diferente
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Filtrar correlaciones significativas
cor_sig_lnc_mrna <- correlaciones_lnc_mrna[correlaciones_lnc_mrna$rho > 0.8 & correlaciones_lnc_mrna$p.value < 0.05, ]
cor_sig_lnc_miRNA <- correlaciones_lnc_miRNA[correlaciones_lnc_miRNA$rho < -0.8 & correlaciones_lnc_miRNA$p.value < 0.05, ]
cor_sig_lnc_miRNA2 <- correlaciones_lnc_miRNA[correlaciones_lnc_miRNA$rho < -0.8 & correlaciones_lnc_miRNA$p.value < 0.05, ]

# 3. Filtrar mRNAs comunes entre miRWalk y lncRNA-mRNA correlaciones
genes_comunes <- intersect(mirwalk$Geneid, cor_sig_lnc_mrna$mRNA)

# 4. Filtrar datasets por esos genes
mirwalk_filtrado <- mirwalk[mirwalk$Geneid %in% genes_comunes, ]
cor_lnc_mrna_filtrado <- cor_sig_lnc_mrna[cor_sig_lnc_mrna$mRNA %in% genes_comunes, ]

# 5. Merge de lncRNA–mRNA con miRNA–mRNA
merged1 <- merge(mirwalk_filtrado, cor_lnc_mrna_filtrado, by.x = "Geneid", by.y = "mRNA")

colnames(merged1)
colnames(cor_sig_lnc_miRNA2)

colnames(merged1)[colnames(merged1) == "mirnaid"] <- "miRNA"

final <- merge(merged1, cor_sig_lnc_miRNA2, by = c("miRNA", "lncRNA"))

head(final)

# 7. Guardar archivo final
write.table(final, "lncRNA_miRNA_mRNA_ceRNA_filtered_FINAL_4downmiRNAs_lncRNAs_up.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

####crear archivos de nodos con
# 1. Leer archivo final
final_data <- read.table("lncRNA_miRNA_mRNA_ceRNA_filtered_final_4downmiRNAs_lncRNAs_up.txt", 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

###Quedé aquí
# 2. Crear data frames de interacciones
# lncRNA → miRNA
edges_lnc_miR <- data.frame(source = final_data$lncRNA, target = final_data$mirnaid, 
                            interaction = "lncRNA_miRNA")

# miRNA → mRNA
edges_miR_mRNA <- data.frame(source = final_data$mirnaid,
  target = final_data$genesymbol,
  interaction = "miRNA_mRNA")

# 3. Combinar ambas relaciones
edges_final <- rbind(edges_lnc_miR, edges_miR_mRNA)

# 4. Eliminar duplicados por si acaso
edges_final <- unique(edges_final)

# 5. Guardar archivo para Cytoscape
write.table(edges_final, "ceRNA_network_edges_for_cytoscape_lncRNAsup_miRNAdown4_new.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

####Este es para extraer los nodos oficiales
# Leer el archivo con la red final
final_data <- read.table("lncRNA_miRNA_mRNA_ceRNA_filtered_final_4downmiRNAs_lncRNAs_up.txt", 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extraer los identificadores únicos
lncRNAs <- unique(final_data$lncRNA)
miRNAs <- unique(final_data$mirnaid)
mRNAs  <- unique(final_data$genesymbol)

# Crear data frames por separado con la columna 'type'
lncRNA_nodes <- data.frame(id = lncRNAs, type = "lncRNA", stringsAsFactors = FALSE)
miRNA_nodes  <- data.frame(id = miRNAs,  type = "miRNA",  stringsAsFactors = FALSE)
mRNA_nodes   <- data.frame(id = mRNAs,   type = "mRNA",   stringsAsFactors = FALSE)

# Unir todos los nodos
nodes_all <- rbind(lncRNA_nodes, miRNA_nodes, mRNA_nodes)

# Verificar si hay duplicados
nodes_all <- unique(nodes_all)

# Guardar archivo final para Cytoscape
write.table(nodes_all, "ceRNA_network_nodes_for_cytoscape_lncRNAsup_mirnas4down_new.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

###Nuevo script para generar los nodos y edges nuevos a partir del filtro más grande
# 1. Cargar el archivo final
final <- read.table("lncRNA_miRNA_mRNA_ceRNA_filtered_final_4downmiRNAs_lncRNAs_up.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Crear archivo de edges

# lncRNA–miRNA edges (negativa)
edges_lnc_miRNA <- unique(data.frame(
  source = final$lncRNA,
  target = final$miRNA,
  interaction = "lncRNA-miRNA",
  rho = final$rho.y,
  p.value = final$p.value.y
))

# miRNA–mRNA edges (no tienen rho/pvalue aquí)
edges_miRNA_mRNA <- unique(data.frame(
  source = final$miRNA,
  target = final$genesymbol,
  interaction = "miRNA-mRNA",
  rho = NA,
  p.value = NA
))

# lncRNA–mRNA edges (positiva)
edges_lnc_mRNA <- unique(data.frame(
  source = final$lncRNA,
  target = final$genesymbol,
  interaction = "lncRNA-mRNA",
  rho = final$rho.x,
  p.value = final$p.value.x
))

# Unir todos los edges
edges_all <- rbind(edges_lnc_miRNA, edges_miRNA_mRNA, edges_lnc_mRNA)

# Guardar edges
write.table(edges_all, "ceRNA_edges_final_miRNAsdown_lncRNAsup_new.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

#Opcional - las interacciones de miRNA-mRNA y lncRNA-miRNA
# miRNA–lncRNA edges (negativa correlación)
edges_miRNA_lncRNA <- unique(data.frame(
  source = final$miRNA,
  target = final$lncRNA,
  interaction = "miRNA-lncRNA",
  rho = final$rho.y,
  p.value = final$p.value.y
))

# miRNA–mRNA edges (sin rho/pvalue)
edges_miRNA_mRNA <- unique(data.frame(
  source = final$miRNA,
  target = final$genesymbol,  # Usamos genesymbol para legibilidad en Cytoscape
  interaction = "miRNA-mRNA",
  rho = NA,
  p.value = NA
))

# Unir ambos tipos de edges
edges_final <- rbind(edges_miRNA_lncRNA, edges_miRNA_mRNA)

# Guardar archivo
write.table(edges_final, "ceRNA_edges_miRNA_lncRNA_mRNA_4mirnas_lncRNAsup_doblenew.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#####
# 3. Crear archivo de nodes

# Extraer nodos únicos con tipo
lncRNAs <- unique(final$lncRNA)
miRNAs <- unique(final$miRNA)
mRNAs  <- unique(final$genesymbol)

nodes <- data.frame(
  id = c(lncRNAs, miRNAs, mRNAs),
  type = c(
    rep("lncRNA", length(lncRNAs)),
    rep("miRNA", length(miRNAs)),
    rep("mRNA",  length(mRNAs))
  )
)

lncRNA_list = subset(nodes, nodes$type == "lncRNA")
write.table(lncRNA_list$id, "lncRNA_list_NE24_miRNAdown-mRNAup_corr.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Guardar nodes
write.table(nodes, "ceRNA_nodes_final_miRNAsdown_lncRNAsup_new.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

# Carga las librerías
library(circlize)
library(dplyr)
library(readr)

# 📥 Leer tu archivo
# Cambia el nombre del archivo al correcto
interacciones <- read_delim("miRNA_gene_collapsed_regions_DESeq2_filtered_downmiRNAs_NE24.txt", delim = "\t")

# 📑 Seleccionar las columnas necesarias y crear las columnas de regulación
interacciones_clean <- interacciones %>%
  select(miRNA = mirnaid,
         mRNA = genesymbol,
         mRNA_level = diffexpressed) %>%
  mutate(miRNA_level = "Downregulated")  # Todos los miRNA son downregulated según tu filtro

# ✅ Comprobar la estructura
print(interacciones_clean)

# 🎨 Definir colores según regulación
# Rojo = Upregulated, Azul = Downregulated
color_miRNA <- rep("#377EB8", length(unique(interacciones_clean$miRNA)))
names(color_miRNA) <- unique(interacciones_clean$miRNA)

color_mRNA <- ifelse(interacciones_clean$mRNA_level == "Upregulated", "#E41A1C", "#377EB8")
names(color_mRNA) <- interacciones_clean$mRNA

#para miRNAs up
color_miRNA <- rep("#E41A1C", length(unique(interacciones_clean$miRNA)))
names(color_miRNA) <- unique(interacciones_clean$miRNA)

color_mRNA <- ifelse(interacciones_clean$mRNA_level == "Downregulated", "#377EB8", "#D55E00")
names(color_mRNA) <- interacciones_clean$mRNA

# 🔗 Definir colores únicos
grid.col <- c(color_miRNA, color_mRNA)

# 🔗 Crear tabla de enlaces
links <- interacciones_clean %>%
  select(from = miRNA, to = mRNA)

# 🔄 Limpiar gráfico anterior si existe
circos.clear()

# 🧠 Crear el chord diagram
chordDiagram(
  links,
  grid.col = grid.col,
  transparency = 0.3,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

# ✍️ Añadir nombres a los sectores
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.6)
}, bg.border = NA)

png("miRNAdown_mRNAup_C24_NE24_3v3_miRWalk.png", res = 600, width = 600, height = 600)

#####Nuevo script ceRNA integrative
# ----------------------------
# 📦 Librerías necesarias
library(dplyr)
library(readr)

# ----------------------------
# 📂 1. Leer archivos de entrada

# Interacción miRNA–mRNA desde miRWalk
mirwalk <- read.table("miRNAdown_collapsed_regions_DESeq2_filtered_mRNAup_NE24_10counts.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Definir las muestras de interés
muestras <- c("genesymbol", "mirnaid", "diffexpressed")

# Subset de columnas de la matriz de expresión
mirwalk <- mirwalk[, muestras]

# Correlación lncRNA–mRNA (positiva)
corr_lnc_mrna <- read.table("lncRNA_mRNA_significant_correlations_lnc_up_C24_NE24.txt", 
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Correlación lncRNA–miRNA (negativa)
corr_lnc_miRNA <- read.table("correlacion_miRNA_vs_lncRNA_all_DE.txt", 
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Interacción lncRNA–miRNA desde IntaRNA (archivo que generamos antes)
intarna_lnc_miRNA <- read_delim("lncRNAup_miRNAdown_filt_C24NE24_default.csv", delim = ",")

####El archivo de geneide y transcript id puede ser el de feelnc o el gtf filtrado
#Usar la función merge

lncRNA_all <- read.table("lncRNAs_candidates_NE624_merge.txt", 
           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

intarna_lnc_miRNA <- merge(lncRNA_all, intarna_lnc_miRNA, 
                           by.x = "transcript_id", by.y = "lncRNA")

# Definir las muestras de interés
samples <- c("transcript_id", "gene_id", "miRNA",
             "lncRNA_regulation", "miRNA_regulation")

# Subset de columnas de la matriz de expresión
intarna_lnc_miRNA <- intarna_lnc_miRNA[, samples]

#Filtrar correlaciones significativas

# lncRNA–mRNA: correlación positiva (rho > 0.8, p < 0.05)
cor_sig_lnc_mrna <- corr_lnc_mrna %>%
  filter(correlation > 0.8, p_value < 0.05)

# lncRNA–miRNA: correlación negativa (rho < -0.8, p < 0.05)
cor_sig_lnc_miRNA <- corr_lnc_miRNA %>%
  filter(rho < -0.8, p.value < 0.05)

# ----------------------------
# 🔗 3. Filtrar mRNAs comunes entre miRWalk y lncRNA–mRNA
genes_comunes <- intersect(mirwalk$genesymbol, cor_sig_lnc_mrna$mRNA_symbol)

# Filtrar datasets por esos genes
mirwalk_filtrado <- mirwalk %>%
  filter(genesymbol %in% genes_comunes)

cor_lnc_mrna_filtrado <- cor_sig_lnc_mrna %>%
  filter(mRNA_symbol %in% genes_comunes)

# ----------------------------
# 🔀 4. Merge miRNA–mRNA con lncRNA–mRNA (nodo común: mRNA)
merged1 <- merge(mirwalk_filtrado, cor_lnc_mrna_filtrado, 
                 by.x = "genesymbol", by.y = "mRNA_symbol")

# ----------------------------
# 🔗 5. Integrar con lncRNA–miRNA

# Primero unimos usando la correlación (negativa) lncRNA–miRNA
red_cerna_cor <- merge(merged1, cor_sig_lnc_miRNA, 
                       by.x = "mirnaid", by.y = "miRNA")

head(red_cerna_cor)

# Renombrar las columnas lncRNA.y y lncRNA.x para mayor claridad
red_cerna_cor <- red_cerna_cor %>%
  rename(
    lncRNA_mRNA = lncRNA.x,  # LncRNA que está correlacionado positivamente con el mRNA
    lncRNA_miRNA = lncRNA.y  # LncRNA que está correlacionado negativamente con el miRNA
  )

# Hacer el merge con el archivo de IntaRNA

red_cerna_final <- merge(red_cerna_cor, intarna_lnc_miRNA, 
                         by.x = c("lncRNA_miRNA", "mirnaid"), 
                         by.y = c("gene_id", "miRNA"))

#columnas de red_cerna_cor: mirnaid, genesymbol, diffexpressed, lncRNA_mRNA, mRNA, correlation, p_value,
#lncRNA_miRNA, rho, p.value
#columnas de intarna_lnc_miRNA: transcript_id, gene_id, miRNA, lncRNA_regulation, miRNA_regulation

# ----------------------------
# 📝 6. Guardar el resultado
write.csv(red_cerna_final, "ceRNA_network_integrada_lncRNAup_miRNAdown_NE24.csv", row.names = FALSE)

write.table(red_cerna_final, "ceRNA_final_miRNAsdown_lncRNAsup_NE24_interact.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

#Estas son las columnas, lncRNA-miRNA, mirnaid, genesymbol, differexpressed, lncRNA_mRNA,
#mRNA, correlation, p_value, rho, p.value, transcript_id, lncRNA_regulation, miRNA_regulation

# ✔️ Resultado: Archivo con columnas de lncRNA, miRNA, mRNA y toda la info de interacciones y correlaciones


# 🔥 9. Filtro biológico clásico ceRNA
# (Ejemplo: lncRNA upregulated, miRNA downregulated, mRNA upregulated)

red_cerna_final_filtrado <- red_cerna_final %>%
  filter(lncRNA_regulation == "upregulated",
         miRNA_regulation == "downregulated",
         diffexpressed == "Upregulated")

# ------------------------------------------
# 🧹 10. Remover duplicados si hay
red_cerna_final_filtrado <- distinct(red_cerna_final_filtrado)

############################################################Aquí se podría tomar
# Crear tabla de edges (lncRNA–miRNA y miRNA–mRNA)

# Edges lncRNA–miRNA
edges_lnc_miRNA <- red_cerna_final_filtrado %>%
  select(Source = lncRNA_miRNA, Target = mirnaid) %>%
  mutate(Interaction = "lncRNA-miRNA")

# Edges miRNA–mRNA
edges_miRNA_mRNA <- red_cerna_final_filtrado %>%
  select(Source = mirnaid, Target = genesymbol) %>%
  mutate(Interaction = "miRNA-mRNA")

# Combinar ambos
edges_cytoscape <- bind_rows(edges_lnc_miRNA, edges_miRNA_mRNA)

# Exportar
write.csv(edges_cytoscape, "edges_ceRNA_cytoscape.csv", row.names = FALSE)

# Crear lista única de nodos
nodos <- unique(c(
  red_cerna_final_filtrado$lncRNA_miRNA,
  red_cerna_final_filtrado$mirnaid,
  red_cerna_final_filtrado$genesymbol
))

# Clasificar tipo de nodo
nodes_cytoscape <- data.frame(
  ID = nodos,
  Type = case_when(
    nodos %in% red_cerna_final_filtrado$lncRNA_miRNA ~ "lncRNA",
    nodos %in% red_cerna_final_filtrado$mirnaid ~ "miRNA",
    nodos %in% red_cerna_final_filtrado$genesymbol ~ "mRNA",
    TRUE ~ "unknown"
  )
)

# Agregar atributos de regulación
nodes_cytoscape <- nodes_cytoscape %>%
  mutate(Regulation = case_when(
    Type == "lncRNA" ~ red_cerna_final_filtrado$lncRNA_regulation[match(ID, red_cerna_final_filtrado$lncRNA_miRNA)],
    Type == "miRNA" ~ red_cerna_final_filtrado$miRNA_regulation[match(ID, red_cerna_final_filtrado$mirnaid)],
    Type == "mRNA"  ~ red_cerna_final_filtrado$diffexpressed[match(ID, red_cerna_final_filtrado$genesymbol)],
    TRUE ~ NA
  ))

# Exportar
write.csv(nodes_cytoscape, "nodes_ceRNA_cytoscape.csv", row.names = FALSE)

##Otra opción
####Este es para extraer los nodos oficiales
# Leer el archivo con la red final
final_data <- read.table("lncRNA_miRNA_mRNA_ceRNA_filtered_final_4downmiRNAs_lncRNAs_up.txt", 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extraer los identificadores únicos
lncRNAs <- unique(final_data$lncRNA)
miRNAs <- unique(final_data$mirnaid)
mRNAs  <- unique(final_data$genesymbol)

# Crear data frames por separado con la columna 'type'
lncRNA_nodes <- data.frame(id = lncRNAs, type = "lncRNA", stringsAsFactors = FALSE)
miRNA_nodes  <- data.frame(id = miRNAs,  type = "miRNA",  stringsAsFactors = FALSE)
mRNA_nodes   <- data.frame(id = mRNAs,   type = "mRNA",   stringsAsFactors = FALSE)

# Unir todos los nodos
nodes_all <- rbind(lncRNA_nodes, miRNA_nodes, mRNA_nodes)

# Verificar si hay duplicados
nodes_all <- unique(nodes_all)

# Guardar archivo final para Cytoscape
write.table(nodes_all, "ceRNA_network_nodes_for_cytoscape_lncRNAsup_mirnas4down_new.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

###Nuevo script para generar los nodos y edges nuevos a partir del filtro más grande
# 1. Cargar el archivo final
final <- read.table("lncRNA_miRNA_mRNA_ceRNA_filtered_final_4downmiRNAs_lncRNAs_up.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

red_cerna_final

# 🔗 Edges lncRNA–miRNA (correlación negativa)
edges_lnc_miRNA <- unique(data.frame(
  source = red_cerna_final$lncRNA_miRNA,     # Columna con el lncRNA
  target = red_cerna_final$mirnaid,      # Columna con el miRNA
  interaction = "lncRNA-miRNA",
  rho = red_cerna_cor$rho_lnc_miRNA, # Asegúrate de que estas columnas correspondan
  p.value = final$pval_lnc_miRNA
))

# 🔗 Edges miRNA–mRNA (sin rho, solo interacción)
edges_miRNA_mRNA <- unique(data.frame(
  source = final$miRNA,
  target = final$genesymbol,
  interaction = "miRNA-mRNA",
  rho = NA,
  p.value = NA
))

# 🔗 Edges lncRNA–mRNA (correlación positiva)
edges_lnc_mRNA <- unique(data.frame(
  source = final$lncRNA,
  target = final$genesymbol,
  interaction = "lncRNA-mRNA",
  rho = final$rho_lnc_mRNA,
  p.value = final$pval_lnc_mRNA
))

# 🔥 Unir todos los edges
edges_all <- rbind(edges_lnc_miRNA, edges_miRNA_mRNA, edges_lnc_mRNA)

# 💾 Guardar edges
write.table(edges_all, "ceRNA_edges_final.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 🧠 Extraer nodos únicos por tipo
lncRNAs <- unique(final$lncRNA)
miRNAs <- unique(final$miRNA)
mRNAs  <- unique(final$genesymbol)

# 🔥 Crear tabla de nodos
nodes <- data.frame(
  id = c(lncRNAs, miRNAs, mRNAs),
  type = c(
    rep("lncRNA", length(lncRNAs)),
    rep("miRNA", length(miRNAs)),
    rep("mRNA",  length(mRNAs))
  )
)

# 💾 Guardar la lista de nodos completa
write.table(nodes, "ceRNA_nodes_final.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# (Opcional) Exportar solo la lista de lncRNAs si deseas
lncRNA_list <- subset(nodes, nodes$type == "lncRNA")
write.table(lncRNA_list$id, "lncRNA_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)


##Otra opción
library(dplyr)
#######Graficar sólo las interacciones miRNA-lncRNA y miRNA-mRNA
# ✅ Asumiendo que tu dataframe se llama 'final_ceRNA' y tiene las siguientes columnas:
# lncRNA_miRNA, mirnaid, genesymbol, diffexpressed, lncRNA_mRNA, mRNA, correlation, p_value,
# rho, p.value, transcript_id, lncRNA_regulation, miRNA_regulation

# 🔥 1. Crear archivo de edges (sin lncRNA–mRNA)

# ✔️ lncRNA–miRNA (correlación negativa con soporte físico)
edges_lnc_miRNA <- red_cerna_final %>%
  select(Source = lncRNA_miRNA, Target = mirnaid) %>%
  mutate(Interaction = "lncRNA-miRNA")

# ✔️ miRNA–mRNA (target según miRWalk)
edges_miRNA_mRNA <- red_cerna_final %>%
  select(Source = mirnaid, Target = genesymbol) %>%
  mutate(Interaction = "miRNA-mRNA")

# 🔥 Unir interacciones válidas
edges_cytoscape <- bind_rows(edges_lnc_miRNA, edges_miRNA_mRNA) %>%
  distinct()

# 💾 Exportar edges en formato .txt (separado por tabulaciones)
write.table(edges_cytoscape, "edges_ceRNA_cytoscape_lncRNAup_miRNAdown_NE.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ Archivo edges_ceRNA_cytoscape.txt generado correctamente.\n")


# 🔥 2. Crear archivo de nodes

# Generar lista de nodos únicos a partir de las interacciones presentes
nodos <- unique(c(red_cerna_final$lncRNA_miRNA, red_cerna_final$mirnaid, red_cerna_final$genesymbol))

# Asignar tipo de nodo
nodes_cytoscape <- data.frame(
  ID = nodos,
  Type = case_when(
    nodos %in% red_cerna_final$lncRNA_miRNA ~ "lncRNA",
    nodos %in% red_cerna_final$mirnaid ~ "miRNA",
    nodos %in% red_cerna_final$genesymbol ~ "mRNA",
    TRUE ~ "unknown"
  )
)

# ✔️ Agregar estado de regulación (Upregulated / Downregulated)
nodes_cytoscape <- nodes_cytoscape %>%
  mutate(Regulation = case_when(
    Type == "lncRNA" ~ red_cerna_final$lncRNA_regulation[match(ID, red_cerna_final$lncRNA_miRNA)],
    Type == "miRNA"  ~ red_cerna_final$miRNA_regulation[match(ID, red_cerna_final$mirnaid)],
    Type == "mRNA"   ~ red_cerna_final$diffexpressed[match(ID, red_cerna_final$genesymbol)],
    TRUE ~ NA
  ))

# 💾 Exportar nodes como .txt (separado por tabulaciones)
write.table(nodes_cytoscape, "nodes_ceRNA_cytoscape_lncRNAup_miRNAdown_NE24.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ Archivo nodes_ceRNA_cytoscape.txt generado correctamente.\n")

###Volcanoplot
## ============================================
## Volcano plot miRNAs – 0 centrado en el eje X
## ============================================

# Paquetes
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggrepel"))  install.packages("ggrepel")

library(tidyverse)
library(ggrepel)

## 1. Leer archivo ------------------------------------------

mires <- read.table("Deseq2_C24_vs_NE24_3v3_miRNAs_norm_new.txt",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

str(mires)
head(mires)

## 2. Definir columnas para el volcano ----------------------

# Si YA tienes diffexpressed (Upregulated / Downregulated / NO), la usamos.
# Si prefieres recalcular, comenta el bloque de abajo y descomenta el "opcional".

mires <- mires %>%
  mutate(
    negLog10padj = -log10(padj),
    significance = case_when(
      diffexpressed == "Upregulated"   ~ "Upregulated",
      diffexpressed == "Downregulated" ~ "Downregulated",
      TRUE                             ~ "Not Significant"
    )
  )

# OPCIONAL: recalcular diffexpressed aquí (ejemplo con padj < 0.05 y |log2FC| > 0.58)
# mires <- mires %>%
#   mutate(
#     negLog10padj = -log10(padj),
#     significance = case_when(
#       padj < 0.05 & log2FoldChange >  0.58 ~ "Upregulated",
#       padj < 0.05 & log2FoldChange < -0.58 ~ "Downregulated",
#       TRUE                                 ~ "Not Significant"
#     )
#   )

# Aseguramos el orden de la leyenda
mires$significance <- factor(mires$significance,
                             levels = c("Downregulated", "Not Significant", "Upregulated"))

## 3. Definir cortes y límites simétricos del eje X ----------

log2FC_cutoff <- 0.58   # ejemplo FC≈1.5, cámbialo si usas otro
padj_cutoff   <- 0.05

# límite máximo en log2FC para que el 0 quede al centro
xmax <- max(abs(mires$log2FoldChange), na.rm = TRUE)
xmax <- ceiling(xmax)   # redondear hacia arriba (queda más bonito)

## 4. Seleccionar miRNAs a etiquetar (opcional) --------------

mi_label <- mires %>%
  filter(padj < padj_cutoff,
         abs(log2FoldChange) > 1)   # aquí eliges tu criterio de “interesantes”

## 5. Volcano plot con 0 centrado ----------------------------

p_volcano <- ggplot(mires,
                    aes(x = log2FoldChange,
                        y = negLog10padj,
                        color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  # Etiquetas (opcional)
  geom_text_repel(data = mi_label,
                  aes(label = Geneid),
                  size = 3,
                  max.overlaps = 30,
                  min.segment.length = 0) +
  # Líneas de corte
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff),   # <-- AQUÍ EL CAMBIO
             linetype = "dashed") +
  # EJE X SIMÉTRICO ALREDEDOR DE 0
  xlim(-xmax, xmax) +
  labs(
    x = "Log2(Fold Change)",
    y = expression(-log[10]("padj")),
    color = "significance",
    title = "Volcano plot miRNAs C24 vs NE24",
    subtitle = "0 centrado en el eje X"
  ) +
  scale_color_manual(values = c(
    "Downregulated"   = "blue",
    "Not Significant" = "grey70",
    "Upregulated"     = "red"
  )) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_volcano

# Guardar
ggsave("Volcano_miRNAs_C24_vs_NE24_centrado.png",
       plot = p_volcano, width = 7, height = 6, dpi = 300)
ggsave("Volcano_miRNAs_C24_vs_NE24_centrado.pdf",
       plot = p_volcano, width = 7, height = 6)
