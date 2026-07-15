suppressPackageStartupMessages({
  library(tximport); library(rtracklayer); library(DESeq2); library(apeglm); library(ggvenn)
  library(readr); library(dplyr); library(tidyverse); library(ggplot2); library(reshape2); library(scales)
  library(pheatmap); library(RColorBrewer); library(ggrepel); library(AnnotationDbi); library(ggplotify)
  library(org.Rn.eg.db); library(biomaRt); library(Hmisc); library(clusterProfiler); library(ggalluvial)
  library(enrichplot); library(tidyr); library(ggnewscale); library(stringr); library(cowplot); library(UpSetR)
})

### PATHS AND DIRECTORIES
gtf_file <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis/merged_CNE_6_24_HISAT2.annotated.gtf"
lncrna_ann_all_file <- "hisat2_6_24_meta_ann.csv"
lncrna_ann_novel_file <- "hisat2_6_24_meta_novel_lncRNA.csv"
base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/salmon_quants_NE624/"
output_dir <- "Publication_Ready_Analysis_mRNA_lncRNA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "1_DE_Tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "2_QC_and_Global_Viz"), showWarnings = FALSE)
dir.create(file.path(output_dir, "3_DE_Visualization"), showWarnings = FALSE)
dir.create(file.path(output_dir, "4_Expression_Patterns"), showWarnings = FALSE)
dir.create(file.path(output_dir, "5_Gene_Lists_for_Enrichr"), showWarnings = FALSE)
dir.create(file.path(output_dir, "6_Coexpression_Analysis/Top_Correlation_Plots"), showWarnings = FALSE)

### IMPORT GTF AND BUILD ANNOTATION MAP
gtf_df <- as.data.frame(rtracklayer::import(gtf_file))
annotation_map <- gtf_df %>%
  filter(type == "transcript") %>%
  dplyr::select(gene_id, gene_name, ref_gene_id) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  group_by(gene_id) %>%
  summarise(
    gene_name = first(na.omit(gene_name)),
    ref_gene_id = first(na.omit(ref_gene_id))
  ) %>%
  ungroup()

genes_with_ref <- annotation_map %>% filter(!is.na(ref_gene_id))
ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", mirror = "uswest")
biomart_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id"),
                     filters = "ensembl_gene_id",
                     values = unique(genes_with_ref$ref_gene_id),
                     mart = ensembl)
biomart_map <- biomart_map %>% rename(ref_gene_id = ensembl_gene_id, official_gene_name = external_gene_name,
                                      official_gene_biotype = gene_biotype, ENTREZID = entrezgene_id)

biomart_map_clean <- biomart_map %>%
  arrange(ref_gene_id) %>%
  distinct(ref_gene_id, .keep_all = TRUE)

annotation_combined <- annotation_map %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  left_join(biomart_map_clean, by = "ref_gene_id") %>%
  mutate(
    SYMBOL = coalesce(official_gene_name, gene_name, gene_id),
    biotype_clean = case_when(
      !is.na(official_gene_biotype) ~ official_gene_biotype,
      startsWith(gene_id, "MSTRG") ~ "lncRNA",
      TRUE ~ "protein_coding"
    ),
    ENSEMBL_ID = coalesce(ref_gene_id, gene_id)
  )

lnc_ann_all <- read_csv(lncrna_ann_all_file, col_types = cols(.default = "c"))
lnc_ann_novel <- read_csv(lncrna_ann_novel_file, col_types = cols(.default = "c"))
lnc_known <- lnc_ann_all %>% filter(gene_biotype == "lncRNA") %>% pull(ensembl_gene_id)
lnc_novel <- lnc_ann_novel %>% pull(gene_id)
annotation_combined <- annotation_combined %>%
  mutate(biotype_clean = ifelse(ENSEMBL_ID %in% c(lnc_known, lnc_novel), "lncRNA", biotype_clean))
write.csv(annotation_combined, file.path(output_dir, "curated_annotation_hybrid.csv"), row.names = FALSE)

### IMPORT SALMON COUNTS AND RUN DESEQ2
samples <- c("Ctrl6_rep1", "Ctrl6_rep2", "Ctrl6_rep4", "NE6_rep1", "NE6_rep2", "NE6_rep4",
             "Ctrl24_rep1", "Ctrl24_rep2", "Ctrl24_rep4", "NE24_rep1", "NE24_rep2", "NE24_rep4")
sampleTable <- data.frame(treatment = factor(rep(c("Control", "NE"), each = 3, times = 2)),
                          time = factor(rep(c("6h", "24h"), each = 6), levels=c("6h", "24h")))
rownames(sampleTable) <- samples
files <- file.path(base_dir, paste0(samples, "_quant"), "quant.sf"); names(files) <- samples
tx2gene <- gtf_df %>% filter(type == "transcript") %>% select(transcript_id, gene_id) %>% distinct() %>% na.omit()

##nueva opción
tx2gene <- gtf_df %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct() %>%
  na.omit()

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ time + treatment + time:treatment)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ time + treatment)
res_lrt <- results(dds_lrt)
write.csv(as.data.frame(res_lrt), file.path(output_dir, "1_DE_Tables/res_LRT_all.csv"))

### EXTRACT DE RESULTS
res_NE6 <- results(dds, name="treatment_NE_vs_Control")
dds_h24 <- dds; dds_h24$time <- relevel(dds_h24$time, "24h"); dds_h24 <- DESeq(dds_h24)
res_NE24 <- results(dds_h24, name="treatment_NE_vs_Control")
res_inter <- results(dds, name="time24h.treatmentNE")

### ADD ANNOTATION AND STATUS
add_status <- function(res_df, annot){
  res_df %>% as.data.frame() %>% rownames_to_column("gene_id") %>%
    left_join(annot, by="gene_id") %>%
    mutate(status = case_when(
      padj < 0.05 & abs(log2FoldChange) >= 1 & biotype_clean == "protein_coding" ~ if_else(log2FoldChange > 0, "Upregulated", "Downregulated"),
      padj < 0.05 & abs(log2FoldChange) >= 0.5 & biotype_clean == "lncRNA" ~ if_else(log2FoldChange > 0, "Upregulated", "Downregulated"),
      TRUE ~ "Not Significant"))
}
res_NE6_p <- add_status(res_NE6, annotation_combined)
res_NE24_p <- add_status(res_NE24, annotation_combined)
res_INT_p <- add_status(res_inter, annotation_combined)
write.csv(res_NE6_p, file.path(output_dir, "1_DE_Tables/DE_results_NE6_vs_Ctrl.csv"), row.names=FALSE)
write.csv(res_NE24_p, file.path(output_dir, "1_DE_Tables/DE_results_NE24_vs_Ctrl.csv"), row.names=FALSE)
write.csv(res_INT_p, file.path(output_dir, "1_DE_Tables/DE_results_interaction.csv"), row.names=FALSE)

### EXPORT GENE LISTS FOR ENRICHR
export_lists <- function(res_proc, timepoint, outdir){
  sig <- res_proc %>% filter(status %in% c("Upregulated","Downregulated"))
  write.table(sig %>% filter(status=="Upregulated") %>% pull(SYMBOL),
              file.path(outdir, paste0("list_upregulated_", timepoint, ".txt")), row.names=F, col.names=F, quote=F)
  write.table(sig %>% filter(status=="Downregulated") %>% pull(SYMBOL),
              file.path(outdir, paste0("list_downregulated_", timepoint, ".txt")), row.names=F, col.names=F, quote=F)
}
export_lists(res_NE6_p, "6h", file.path(output_dir,"5_Gene_Lists_for_Enrichr"))
export_lists(res_NE24_p, "24h", file.path(output_dir,"5_Gene_Lists_for_Enrichr"))

### QC AND PCA
sampleTable$group <- paste(sampleTable$treatment, sampleTable$time, sep="_")
p <- plotPCA(vsd, intgroup=c("group")) + geom_text_repel(aes(label=name)) + theme_bw()
ggsave(file.path(output_dir,"2_QC_and_Global_Viz/PCA_plot.png"),p,width=8,height=6)

#Segunda opción
colData(vsd)$group <- paste(vsd$treatment, vsd$time, sep="_")
p <- plotPCA(vsd, intgroup="group") +
  geom_text_repel(aes(label=name)) +
  labs(title="Principal Component Analysis (PCA)") +
  theme_bw()

ggsave(file.path(output_dir,"2_QC_and_Global_Viz/PCA_plot.png"), p, width=8, height=6)

### VOLCANO PLOTS
create_volcano <- function(res_data, bio_filter, title){
  df <- res_data %>% filter(biotype_clean==bio_filter,!is.na(padj))
  df <- df %>% mutate(padj = ifelse(padj==0,.Machine$double.xmin,padj),
                      label = ifelse(status!="Not Significant" & (abs(log2FoldChange)>2 | padj<1e-15),SYMBOL,""))
  ggplot(df,aes(x=log2FoldChange,y=-log10(padj),color=status,label=label))+
    geom_point(alpha=0.7,size=2)+
    geom_text_repel(max.overlaps=15,size=3.5)+
    geom_hline(yintercept=-log10(0.05),linetype="dashed")+
    scale_color_manual(values=c("Upregulated"="#e41a1c","Downregulated"="#377eb8","Not Significant"="grey80"))+
    labs(title=title,x="Log2 Fold Change",y="-log10(padj)")+theme_bw()+theme(legend.position="top")
}
ggsave(file.path(output_dir,"3_DE_Visualization/Volcano_NE6_proteinCoding.png"),
       create_volcano(res_NE6_p,"protein_coding","6h: Protein-Coding"),width=8,height=6)
ggsave(file.path(output_dir,"3_DE_Visualization/Volcano_NE6_lncRNA.png"),
       create_volcano(res_NE6_p,"lncRNA","6h: lncRNA"),width=8,height=6)
ggsave(file.path(output_dir,"3_DE_Visualization/Volcano_NE24_proteinCoding.png"),
       create_volcano(res_NE24_p,"protein_coding","24h: Protein-Coding"),width=8,height=6)
ggsave(file.path(output_dir,"3_DE_Visualization/Volcano_NE24_lncRNA.png"),
       create_volcano(res_NE24_p,"lncRNA","24h: lncRNA"),width=8,height=6)

### LRT HEATMAP
lrt_sig <- subset(res_lrt, padj<0.05)
top_kin <- rownames(lrt_sig)[order(lrt_sig$padj)[1:50]]
mat <- assay(vsd)[top_kin,]; mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=as.data.frame(colData(dds)[,c("treatment","time")]),
         main="Top 50 Dynamic Genes (LRT)",
         filename=file.path(output_dir,"4_Expression_Patterns/Heatmap_Top50_LRT.png"),
         scale="row",border_color="grey60",width=8,height=12)

### CO-EXPRESSION (lncRNA–mRNA)
dir.create(file.path(output_dir, "6_Coexpression_Analysis"), showWarnings = FALSE, recursive = TRUE)

run_corr <- function(res_proc, vsd_matrix, timepoint, annot){
  de_lnc <- res_proc %>% dplyr::filter(biotype_clean == "lncRNA", status != "Not Significant")
  de_mr  <- res_proc %>% dplyr::filter(biotype_clean == "protein_coding", status != "Not Significant")
  if(nrow(de_lnc) < 2 | nrow(de_mr) < 2) return(NULL)
  
  lnc <- assay(vsd_matrix)[de_lnc$gene_id, ]
  mr  <- assay(vsd_matrix)[de_mr$gene_id, ]
  cor_res <- Hmisc::rcorr(t(lnc), t(mr), type = "pearson")
  
  cor_mat <- cor_res$r
  p_mat   <- cor_res$P
  idx <- which(abs(cor_mat) > 0.8 & p_mat < 0.05, arr.ind = TRUE)
  if(nrow(idx) == 0) return(NULL)
  
  df <- data.frame(
    lncRNA_id = rownames(cor_mat)[idx[,1]],
    mRNA_id   = colnames(cor_mat)[idx[,2]],
    correlation = cor_mat[idx],
    p_value     = p_mat[idx]
  ) %>%
    dplyr::left_join(annot %>% dplyr::select(gene_id, SYMBOL), by = c("lncRNA_id" = "gene_id")) %>%
    dplyr::rename(lncRNA_symbol = SYMBOL) %>%
    dplyr::left_join(annot %>% dplyr::select(gene_id, SYMBOL), by = c("mRNA_id" = "gene_id")) %>%
    dplyr::rename(mRNA_symbol = SYMBOL) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    dplyr::filter(p_adj < 0.05)
  
  write.csv(df, file.path(output_dir, "6_Coexpression_Analysis",
                          paste0("significant_correlations_", timepoint, ".csv")),
            row.names = FALSE)
  return(df)
}

corr6 <- run_corr(res_NE6_p,vsd,"6h",annotation_combined)
corr24 <- run_corr(res_NE24_p,vsd,"24h",annotation_combined)
all_corr <- bind_rows(corr6,corr24) %>% distinct(lncRNA_id,mRNA_id,.keep_all=TRUE)

### IDENTIFY HUBS (VERSION CORREGIDA)
if(!is.null(all_corr) && nrow(all_corr) > 0) {
  hub_sum <- all_corr %>%
    dplyr::group_by(lncRNA_id, lncRNA_symbol) %>%
    dplyr::summarise(n_correlated_mRNAs = dplyr::n(), .groups = "drop")
  
  hubs <- hub_sum %>%
    dplyr::filter(n_correlated_mRNAs >= 10)
  
  if(nrow(hubs) > 0) {
    hub_table <- all_corr %>%
      dplyr::filter(lncRNA_id %in% hubs$lncRNA_id) %>%
      dplyr::left_join(hubs, by = c("lncRNA_id", "lncRNA_symbol")) %>%
      dplyr::arrange(dplyr::desc(n_correlated_mRNAs),
                     lncRNA_symbol,
                     dplyr::desc(abs(correlation))) %>%
      dplyr::select(
        lncRNA_hub_Symbol = lncRNA_symbol,
        Correlated_mRNA_Symbol = mRNA_symbol,
        Correlation = correlation,
        P_value = p_value,
        Adj_p = p_adj,
        lncRNA_hub_ID = lncRNA_id,
        Correlated_mRNA_ID = mRNA_id
      )
    
    dir.create(file.path(output_dir, "6_Coexpression_Analysis"), 
               showWarnings = FALSE, recursive = TRUE)
    
    write.csv(hub_table,
              file.path(output_dir,
                        "6_Coexpression_Analysis/lncRNA_hubs_to_mRNA.csv"),
              row.names = FALSE)
  }
}

### FUNCTIONAL ENRICHMENT ANALYSIS (GSEA: GO + KEGG)
message("--- Step: Functional Enrichment (GSEA) ---")

# Crear carpetas
dir.create(file.path(output_dir, "7_Functional_Enrichment_GSEA"), showWarnings = FALSE)
dir.create(file.path(output_dir, "7_Functional_Enrichment_GSEA/Plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables/Activated"), showWarnings = FALSE)
dir.create(file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables/Suppressed"), showWarnings = FALSE)

# Función para convertir ENTREZID a SYMBOL dentro de core_enrichment
convert_core_enrichment <- function(df) {
  if (!"core_enrichment" %in% colnames(df)) return(df)
  entrez_to_symbol <- AnnotationDbi::select(
    org.Rn.eg.db,
    keys = unique(unlist(strsplit(paste(df$core_enrichment, collapse = "/"), "/"))),
    columns = "SYMBOL",
    keytype = "ENTREZID"
  ) %>% dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  df$core_enrichment_symbols <- sapply(df$core_enrichment, function(x) {
    ids <- unlist(strsplit(x, "/"))
    syms <- entrez_to_symbol$SYMBOL[match(ids, entrez_to_symbol$ENTREZID)]
    syms <- syms[!is.na(syms)]
    paste(unique(syms), collapse = "/")
  })
  return(df)
}

# --- Función principal GSEA ---
gsea_analysis <- function(res_df, timepoint, output_dir) {
  df <- res_df %>%
    dplyr::filter(!is.na(ENTREZID) & !is.na(log2FoldChange)) %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(log2FC = mean(log2FoldChange)) %>%
    dplyr::ungroup()
  
  geneList <- df$log2FC
  names(geneList) <- df$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  message(paste("... Running GSEA for", timepoint, "with", length(geneList), "genes"))
  
  # --- GSEA GO ---
  gsea_go <- clusterProfiler::gseGO(
    geneList = geneList,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  # --- GSEA KEGG ---
  gsea_kegg <- clusterProfiler::gseKEGG(
    geneList = geneList,
    organism = "rno",
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  # --- Convertir core_enrichment a símbolos ---
  gsea_go_df <- convert_core_enrichment(as.data.frame(gsea_go))
  gsea_kegg_df <- convert_core_enrichment(as.data.frame(gsea_kegg))
  
  # --- Guardar tablas completas ---
  write.csv(gsea_go_df, file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables",
                                  paste0("GSEA_GO_", timepoint, ".csv")), row.names = FALSE)
  write.csv(gsea_kegg_df, file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables",
                                    paste0("GSEA_KEGG_", timepoint, ".csv")), row.names = FALSE)
  
  # --- Guardar activados/suprimidos ---
  if ("NES" %in% colnames(gsea_go_df)) {
    go_up <- gsea_go_df %>% dplyr::filter(NES > 0)
    go_down <- gsea_go_df %>% dplyr::filter(NES < 0)
    write.csv(go_up, file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables/Activated",
                               paste0("GSEA_GO_", timepoint, "_Activated.csv")), row.names = FALSE)
    write.csv(go_down, file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables/Suppressed",
                                 paste0("GSEA_GO_", timepoint, "_Suppressed.csv")), row.names = FALSE)
  }
  if ("NES" %in% colnames(gsea_kegg_df)) {
    kegg_up <- gsea_kegg_df %>% dplyr::filter(NES > 0)
    kegg_down <- gsea_kegg_df %>% dplyr::filter(NES < 0)
    write.csv(kegg_up, file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables/Activated",
                                 paste0("GSEA_KEGG_", timepoint, "_Activated.csv")), row.names = FALSE)
    write.csv(kegg_down, file.path(output_dir, "7_Functional_Enrichment_GSEA/Tables/Suppressed",
                                   paste0("GSEA_KEGG_", timepoint, "_Suppressed.csv")), row.names = FALSE)
  }
  
  # --- Dotplot personalizado: Top 8 activados + 8 suprimidos ---
  if ("NES" %in% colnames(gsea_go_df)) {
    top_up <- gsea_go_df %>% dplyr::filter(NES > 0) %>% dplyr::arrange(desc(NES)) %>% head(8)
    top_down <- gsea_go_df %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES) %>% head(8)
    top_combined <- dplyr::bind_rows(top_up, top_down)
    top_combined$Description <- stringr::str_wrap(top_combined$Description, width = 55)
    top_combined$Description <- factor(top_combined$Description, levels = top_combined$Description[order(top_combined$NES)])
    
    p_dot <- ggplot(top_combined, aes(x = NES, y = Description, size = setSize, color = NES)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c", midpoint = 0) +
      labs(title = paste("Top 8 Activated & Suppressed GO BP -", timepoint),
           x = "Normalized Enrichment Score (NES)",
           y = NULL) +
      theme_bw(base_size = 13) +
      theme(axis.text.y = element_text(size = 10.5),
            axis.text.x = element_text(size = 11),
            plot.title = element_text(face = "bold", hjust = 0.5))
    
    # Guardar PDF y PNG
    plot_base <- file.path(output_dir, "7_Functional_Enrichment_GSEA/Plots",
                           paste0("GSEA_GO_", timepoint, "_Top16_Dotplot"))
    ggsave(paste0(plot_base, ".pdf"), plot = p_dot, width = 8.5, height = 6.5)
    ggsave(paste0(plot_base, ".png"), plot = p_dot, width = 8.5, height = 6.5, dpi = 300)
  }
  
  # --- Ridgeplot (mejorado) ---
  pdf(file.path(output_dir, "7_Functional_Enrichment_GSEA/Plots",
                paste0("GSEA_KEGG_", timepoint, "_Ridgeplot.pdf")), width = 8, height = 6)
  print(ridgeplot(gsea_kegg, showCategory = 20) +
          theme(axis.text.y = element_text(size = 9),
                plot.title = element_text(hjust = 0.5, face = "bold")))
  dev.off()
  ggsave(file.path(output_dir, "7_Functional_Enrichment_GSEA/Plots",
                   paste0("GSEA_KEGG_", timepoint, "_Ridgeplot.png")),
         plot = ridgeplot(gsea_kegg, showCategory = 20) +
           theme(axis.text.y = element_text(size = 9),
                 plot.title = element_text(hjust = 0.5, face = "bold")),
         width = 8, height = 6, dpi = 300)
  
  # --- Emapplots (similitud) ---
  pdf(file.path(output_dir, "7_Functional_Enrichment_GSEA/Plots",
                paste0("GSEA_GO_", timepoint, "_emapplot.pdf")), width = 10, height = 8)
  print(emapplot(pairwise_termsim(gsea_go), showCategory = 20))
  dev.off()
  
  pdf(file.path(output_dir, "7_Functional_Enrichment_GSEA/Plots",
                paste0("GSEA_KEGG_", timepoint, "_emapplot.pdf")), width = 10, height = 8)
  print(emapplot(pairwise_termsim(gsea_kegg), showCategory = 20))
  dev.off()
  
  message(paste("... GSEA completed for", timepoint))
  return(list(gsea_go = gsea_go, gsea_kegg = gsea_kegg))
}

# Ejecutar GSEA para ambos tiempos
gsea6 <- gsea_analysis(res_NE6_p, "6h", output_dir)
gsea24 <- gsea_analysis(res_NE24_p, "24h", output_dir)

message("--- GSEA analyses saved under: 7_Functional_Enrichment_GSEA ---")

### Comparative analysis 6h vs 24h — DEGs + GSEA (GO & KEGG)

# --- INPUT FILES ---
DE6_file <- "Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables/DE_results_NE6_vs_Ctrl.csv"
DE24_file <- "Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables//DE_results_NE24_vs_Ctrl.csv"
GO6_file <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_GO_6h.csv"
GO24_file <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_GO_24h.csv"
KEGG6_file <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_KEGG_6h.csv"
KEGG24_file <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_KEGG_24h.csv"

# --- OUTPUT DIRECTORY ---
outdir <- "Comparative_6h_24h_Analysis"
dir.create(outdir, showWarnings = FALSE)

### 1️⃣  Differential expression comparison
DE6 <- read_csv(DE6_file)
DE24 <- read_csv(DE24_file)

DE6_sig <- DE6 %>% filter(status %in% c("Upregulated","Downregulated"))
DE24_sig <- DE24 %>% filter(status %in% c("Upregulated","Downregulated"))

# Shared and unique genes
shared <- intersect(DE6_sig$gene_id, DE24_sig$gene_id)
unique6 <- setdiff(DE6_sig$gene_id, DE24_sig$gene_id)
unique24 <- setdiff(DE24_sig$gene_id, DE6_sig$gene_id)

# Direction consistency
merged_shared <- DE6_sig %>%
  filter(gene_id %in% shared) %>%
  select(gene_id, SYMBOL_6h = SYMBOL, status_6h = status, log2FoldChange_6h = log2FoldChange) %>%
  left_join(
    DE24_sig %>%
      select(gene_id, SYMBOL_24h = SYMBOL, status_24h = status, log2FoldChange_24h = log2FoldChange),
    by = "gene_id"
  ) %>%
  mutate(direction_class = case_when(
    status_6h == "Upregulated" & status_24h == "Upregulated" ~ "Persistent Up",
    status_6h == "Downregulated" & status_24h == "Downregulated" ~ "Persistent Down",
    TRUE ~ "Switch (Opposite)"
  ))

# Counts summary
summary_DEG <- tibble(
  n_sig_6h = nrow(DE6_sig),
  n_sig_24h = nrow(DE24_sig),
  n_shared = length(shared),
  n_unique_6h = length(unique6),
  n_unique_24h = length(unique24)
)
write_csv(summary_DEG, file.path(outdir, "summary_DEG_counts.csv"))
write_csv(merged_shared, file.path(outdir, "shared_DEGs_direction.csv"))

# Top persistent genes
persistent <- merged_shared %>%
  filter(direction_class %in% c("Persistent Up","Persistent Down")) %>%
  mutate(mean_absLFC = (abs(log2FoldChange_6h) + abs(log2FoldChange_24h))/2) %>%
  arrange(desc(mean_absLFC))
write_csv(persistent, file.path(outdir, "top_persistent_DEGs.csv"))

cat("✓ DEG comparison done\n")

### 2️⃣  Functional enrichment (GSEA GO & KEGG)

prep_gsea <- function(file){
  df <- read_csv(file, show_col_types = FALSE)
  cols <- tolower(names(df))
  desc <- names(df)[which(cols == "description")]
  nes <- names(df)[which(cols == "nes")]
  padj <- names(df)[which(cols == "p.adjust")]
  core <- names(df)[grepl("core_enrichment", cols)]
  keep <- c(desc, nes, padj, core[1])
  df <- df[, keep]
  colnames(df)[1:3] <- c("Description","NES","padj")
  df
}

compare_gsea <- function(file6, file24, label){
  df6 <- prep_gsea(file6)
  df24 <- prep_gsea(file24)
  
  sig6 <- df6 %>% filter(padj < 0.05)
  sig24 <- df24 %>% filter(padj < 0.05)
  
  act6 <- sig6 %>% filter(NES > 0)
  act24 <- sig24 %>% filter(NES > 0)
  sup6 <- sig6 %>% filter(NES < 0)
  sup24 <- sig24 %>% filter(NES < 0)
  
  shared_act <- act6 %>% inner_join(act24, by="Description", suffix=c("_6h","_24h")) %>%
    mutate(mean_absNES=(abs(NES_6h)+abs(NES_24h))/2) %>%
    arrange(desc(mean_absNES))
  shared_sup <- sup6 %>% inner_join(sup24, by="Description", suffix=c("_6h","_24h")) %>%
    mutate(mean_absNES=(abs(NES_6h)+abs(NES_24h))/2) %>%
    arrange(desc(mean_absNES))
  
  uniq6_act <- act6 %>% filter(!Description %in% shared_act$Description)
  uniq24_act <- act24 %>% filter(!Description %in% shared_act$Description)
  uniq6_sup <- sup6 %>% filter(!Description %in% shared_sup$Description)
  uniq24_sup <- sup24 %>% filter(!Description %in% shared_sup$Description)
  
  summary <- tibble(
    set = label,
    n_act_6h = nrow(act6), n_act_24h = nrow(act24), shared_act = nrow(shared_act),
    uniq_act_6h = nrow(uniq6_act), uniq_act_24h = nrow(uniq24_act),
    n_sup_6h = nrow(sup6), n_sup_24h = nrow(sup24), shared_sup = nrow(shared_sup),
    uniq_sup_6h = nrow(uniq6_sup), uniq_sup_24h = nrow(uniq24_sup)
  )
  
  list(summary=summary,
       shared_act=shared_act, uniq6_act=uniq6_act, uniq24_act=uniq24_act,
       shared_sup=shared_sup, uniq6_sup=uniq6_sup, uniq24_sup=uniq24_sup)
}

GO <- compare_gsea(GO6_file, GO24_file, "GO")
KEGG <- compare_gsea(KEGG6_file, KEGG24_file, "KEGG")

# Export
write_csv(bind_rows(GO$summary, KEGG$summary), file.path(outdir, "summary_GSEA_counts.csv"))
write_csv(GO$shared_act, file.path(outdir, "GO_shared_activated.csv"))
write_csv(GO$shared_sup, file.path(outdir, "GO_shared_suppressed.csv"))
write_csv(GO$uniq6_act, file.path(outdir, "GO_unique6_activated.csv"))
write_csv(GO$uniq24_act, file.path(outdir, "GO_unique24_activated.csv"))
write_csv(GO$uniq6_sup, file.path(outdir, "GO_unique6_suppressed.csv"))
write_csv(GO$uniq24_sup, file.path(outdir, "GO_unique24_suppressed.csv"))

write_csv(KEGG$shared_act, file.path(outdir, "KEGG_shared_activated.csv"))
write_csv(KEGG$shared_sup, file.path(outdir, "KEGG_shared_suppressed.csv"))
write_csv(KEGG$uniq6_act, file.path(outdir, "KEGG_unique6_activated.csv"))
write_csv(KEGG$uniq24_act, file.path(outdir, "KEGG_unique24_activated.csv"))
write_csv(KEGG$uniq6_sup, file.path(outdir, "KEGG_unique6_suppressed.csv"))
write_csv(KEGG$uniq24_sup, file.path(outdir, "KEGG_unique24_suppressed.csv"))

cat("✓ Functional comparison done\n")

cat("\nAll comparative tables saved in folder:", outdir, "\n")

### Visualization of Comparative Analysis 6h vs 24h
### (Venn, Scatter, Dotplot GSEA, Sankey)

# --- Paths ---
indir <- "Comparative_6h_24h_Analysis"
outdir <- "Comparative_6h_24h_Figures"
dir.create(outdir, showWarnings = FALSE)

# --- Helper theme ---
theme_pub <- theme_bw(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        axis.text = element_text(color = "black"))

### 1️⃣ Venn Diagrams (mRNAs & lncRNAs)

# Load DE results
DE6 <- read_csv("Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables/DE_results_NE6_vs_Ctrl.csv")
DE24 <- read_csv("Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables/DE_results_NE24_vs_Ctrl.csv")

# Filter by biotype
DE6_mRNA <- DE6 %>% filter(biotype_clean == "protein_coding", padj < 0.05)
DE24_mRNA <- DE24 %>% filter(biotype_clean == "protein_coding", padj < 0.05)
DE6_lnc  <- DE6 %>% filter(biotype_clean == "lncRNA", padj < 0.05)
DE24_lnc <- DE24 %>% filter(biotype_clean == "lncRNA", padj < 0.05)

# Create lists for Venn
venn_mRNA <- list(`6h` = DE6_mRNA$gene_id, `24h` = DE24_mRNA$gene_id)
venn_lnc  <- list(`6h` = DE6_lnc$gene_id, `24h` = DE24_lnc$gene_id)

# Plot Venn for mRNAs
p1 <- ggvenn(venn_mRNA, fill_color = c("#377eb8","#e41a1c"),
             show_percentage = TRUE, stroke_size = 0.5, set_name_size = 5) +
  ggtitle("Differentially Expressed mRNAs (6h vs 24h)") +
  theme_pub
ggsave(file.path(outdir, "Venn_mRNA.png"), p1, width=6, height=5)
ggsave(file.path(outdir, "Venn_mRNA.pdf"), p1, width=6, height=5)

# Plot Venn for lncRNAs
p2 <- ggvenn(venn_lnc, fill_color = c("#4daf4a","#984ea3"),
             show_percentage = TRUE, stroke_size = 0.5, set_name_size = 5) +
  ggtitle("Differentially Expressed lncRNAs (6h vs 24h)") +
  theme_pub
ggsave(file.path(outdir, "Venn_lncRNA.png"), p2, width=6, height=5)
ggsave(file.path(outdir, "Venn_lncRNA.pdf"), p2, width=6, height=5)

cat("✓ Venn diagrams generated\n")

### 2️⃣ Scatter plot of log2FC (6h vs 24h)

shared <- read_csv(file.path(indir, "shared_DEGs_direction.csv"))
scatter <- shared %>%
  mutate(direction = factor(direction_class, levels=c("Persistent Up","Persistent Down","Switch (Opposite)")))

p3 <- ggplot(scatter, aes(x=log2FoldChange_6h, y=log2FoldChange_24h, color=direction)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_point(alpha=0.7, size=2.2) +
  geom_abline(slope=1, intercept=0, linetype="dotted") +
  scale_color_manual(values=c("Persistent Up"="#e41a1c","Persistent Down"="#377eb8","Switch (Opposite)"="grey50")) +
  labs(title="Consistency of Expression Changes (6h vs 24h)",
       x="log2FC (6h)", y="log2FC (24h)") + theme_pub
ggsave(file.path(outdir, "Scatter_log2FC_shared.png"), p3, width=6.5, height=6)
ggsave(file.path(outdir, "Scatter_log2FC_shared.pdf"), p3, width=6.5, height=6)

cat("✓ Scatter plot generated\n")

### 3️⃣ Dot plot of shared and unique GSEA GO terms
### 3️⃣ Optimized GO visualization (Dotplot + NES Barplot)

GO_shared <- read_csv(file.path(indir, "GO_shared_activated.csv"))
GO_uniq6  <- read_csv(file.path(indir, "GO_unique6_activated.csv"))
GO_uniq24 <- read_csv(file.path(indir, "GO_unique24_activated.csv"))

# Combine and clean
GO_shared$Category <- "Shared"
GO_uniq6$Category  <- "Unique_6h"
GO_uniq24$Category <- "Unique_24h"
GO_combined <- bind_rows(GO_shared, GO_uniq6, GO_uniq24)

# --- Limit to top 8 processes per category ---
GO_top <- GO_combined %>%
  mutate(NES_plot = coalesce(NES_6h, NES)) %>%
  group_by(Category) %>%
  slice_max(order_by = abs(NES_plot), n = 8) %>%
  ungroup() %>%
  arrange(desc(NES_plot))

# --- Wrap long labels for better readability ---
GO_top$Description <- str_wrap(GO_top$Description, width = 45)

# --- DOT PLOT (publication-style) ---
p4_dot <- ggplot(GO_top, aes(x=Category, 
                             y=reorder(Description, NES_plot),
                             size=-log10(padj_6h %||% padj), 
                             color=NES_plot)) +
  geom_point(alpha=0.9) +
  scale_color_gradient2(low="#377eb8", mid="grey90", high="#e41a1c") +
  labs(title="Top GO Biological Processes (Activated)",
       x="", y="", size="-log10(padj)", color="NES") +
  theme_pub +
  theme(
    axis.text.y = element_text(size=10, color="black", hjust=1),
    axis.text.x = element_text(size=12),
    legend.position="right",
    plot.title = element_text(face="bold", size=14, hjust=0.5),
    plot.margin = margin(10,180,10,10)
  )

ggsave(file.path(outdir, "GO_dotplot_activated_clean.png"), p4_dot, width=9, height=6)
ggsave(file.path(outdir, "GO_dotplot_activated_clean.pdf"), p4_dot, width=9, height=6)

# --- BARPLOT (alternative style for paper) ---
p4_bar <- ggplot(GO_top, aes(x=reorder(Description, NES_plot), 
                             y=NES_plot, fill=Category)) +
  geom_col(width=0.7, color="black", alpha=0.85) +
  coord_flip() +
  scale_fill_manual(values=c("Shared"="#e41a1c",
                             "Unique_6h"="#377eb8",
                             "Unique_24h"="#4daf4a")) +
  labs(title="Top GO Biological Processes (Activated)",
       x="", y="Normalized Enrichment Score (NES)") +
  theme_pub +
  theme(
    axis.text.y = element_text(size=10, color="black"),
    legend.position="top",
    plot.title = element_text(face="bold", size=14, hjust=0.5)
  )

ggsave(file.path(outdir, "GO_barplot_activated.png"), p4_bar, width=9, height=6)
ggsave(file.path(outdir, "GO_barplot_activated.pdf"), p4_bar, width=9, height=6)

cat("✓ Improved GO dotplot and barplot generated\n")

### 4️⃣ Sankey plot (process transition flow)

# Example with GO shared + unique (Activated)
if(nrow(GO_shared)>0){
  GO_sankey <- bind_rows(
    GO_shared %>% mutate(Source="6h", Target="24h"),
    GO_uniq6 %>% mutate(Source="6h", Target="Unique_6h"),
    GO_uniq24 %>% mutate(Source="24h", Target="Unique_24h")
  ) %>%
    select(Source, Target, Description)
  
  if(nrow(GO_sankey)>0){
    p5 <- ggplot(GO_sankey,
                 aes(axis1 = Source, axis2 = Target, y = 1)) +
      geom_alluvium(aes(fill = Source), width = 1/12) +
      geom_stratum(width = 1/12, fill = "grey80", color = "grey40") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=4) +
      theme_void() +
      ggtitle("Functional Process Flow (Activated GO Terms)") +
      scale_fill_manual(values=c("6h"="#377eb8","24h"="#e41a1c",
                                 "Unique_6h"="#4daf4a","Unique_24h"="#984ea3"))
    ggsave(file.path(outdir, "GO_Sankey_activated.png"), p5, width=7, height=5)
    ggsave(file.path(outdir, "GO_Sankey_activated.pdf"), p5, width=7, height=5)
    cat("✓ Sankey plot generated\n")
  }
}

cat("\nAll visualizations saved in folder:", outdir, "\n")

make_gsea_plot <- function(file6, file24, label, outprefix) {
  df6 <- read_csv(file6, show_col_types = FALSE) %>% mutate(Time = "6h")
  df24 <- read_csv(file24, show_col_types = FALSE) %>% mutate(Time = "24h")
  
  # Detect relevant columns robustly
  padj_col <- intersect(names(df6), c("p.adjust", "p.adj", "padj", "pvalue", "p_val_adj"))[1]
  nes_col  <- intersect(names(df6), c("NES", "nes", "Normalized enrichment score", "normalized_enrichment_score"))[1]
  desc_col <- intersect(names(df6), c("Description", "description", "Term", "term", "Pathway", "pathway"))[1]
  
  if (is.na(padj_col) | is.na(nes_col) | is.na(desc_col)) {
    stop("Missing one of the required columns (Description / NES / p.adjust) in: ", file6)
  }
  
  # Merge and rename
  all_df <- bind_rows(
    df6 %>% select(all_of(c(desc_col, nes_col, padj_col, "Time"))),
    df24 %>% select(all_of(c(desc_col, nes_col, padj_col, "Time")))
  )
  colnames(all_df) <- c("Description", "NES", "padj", "Time")
  
  # Clean numeric columns
  all_df <- all_df %>%
    mutate(
      NES = as.numeric(NES),
      padj = as.numeric(padj),
      Time = factor(Time, levels = c("6h", "24h"))
    ) %>%
    filter(!is.na(NES), !is.na(padj))
  
  # Select top 20 most significant terms per time
  top_df <- all_df %>%
    group_by(Time) %>%
    arrange(padj) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  # Wrap long labels
  top_df$Description <- stringr::str_wrap(top_df$Description, width = 45)
  
  # Comparative GSEA dot plot
  p <- ggplot(top_df, aes(x = -log10(padj),
                          y = reorder(Description, NES),
                          size = -log10(padj),
                          color = NES)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ Time, ncol = 2, scales = "free_y") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "NES") +
    scale_size_continuous(name = "-log10(padj)") +
    labs(title = paste0("Comparative GSEA ", label, " (6h vs 24h)"),
         x = "-log10(p.adjust)",
         y = "Pathway / Process") +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(outdir, paste0(outprefix, "_GSEA_", label, ".png")),
         p, width = 9, height = 8, dpi = 300)
  ggsave(file.path(outdir, paste0(outprefix, "_GSEA_", label, ".pdf")),
         p, width = 9, height = 8)
  
  return(p)
}

# --- Generate comparative GSEA plots ---
p_GO   <- make_gsea_plot(GO6, GO24, "GO Biological Process", "Comparative")
p_KEGG <- make_gsea_plot(KEGG6, KEGG24, "KEGG Pathways", "Comparative")

cat("✓ Comparative GSEA plots (GO + KEGG) generated successfully.\n")

# Comparative lncRNA–mRNA Correlation Analysis (6h vs 24h)
# Extended version with functional summary and core gene lists

# === INPUT FILES ===
corr6_file   <- "Publication_Ready_Analysis_mRNA_lncRNA/6_Coexpression_Analysis/significant_correlations_6h.csv"
corr24_file  <- "Publication_Ready_Analysis_mRNA_lncRNA/6_Coexpression_Analysis/significant_correlations_24h.csv"
annot_file   <- "Publication_Ready_Analysis_mRNA_lncRNA/curated_annotation_hybrid.csv"
GSEA_BP_6h   <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_GO_6h.csv"
GSEA_BP_24h  <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_GO_24h.csv"
GSEA_KEGG_6h <- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_KEGG_6h.csv"
GSEA_KEGG_24h<- "Publication_Ready_Analysis_mRNA_lncRNA/7_Functional_Enrichment_GSEA/Tables/GSEA_KEGG_24h.csv"

outdir <- "Comparative_lncrna_mrna_Analysis"
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "Plots"), showWarnings = FALSE)

# === LOAD DATA ===
corr6  <- read_csv(corr6_file, show_col_types = FALSE)
corr24 <- read_csv(corr24_file, show_col_types = FALSE)
annot  <- read_csv(annot_file, show_col_types = FALSE)

# === IDENTIFY SHARED AND UNIQUE INTERACTIONS ===
shared_pairs <- inner_join(corr6, corr24, by = c("lncRNA_id", "mRNA_id"))
unique_6h    <- anti_join(corr6, corr24, by = c("lncRNA_id", "mRNA_id"))
unique_24h   <- anti_join(corr24, corr6, by = c("lncRNA_id", "mRNA_id"))

write_csv(shared_pairs, file.path(outdir, "shared_pairs.csv"))
write_csv(unique_6h, file.path(outdir, "unique_6h.csv"))
write_csv(unique_24h, file.path(outdir, "unique_24h.csv"))

# === SUMMARY COUNTS ===
summary_df <- tibble(
  Category = c("Shared", "6h only", "24h only"),
  Pairs = c(nrow(shared_pairs), nrow(unique_6h), nrow(unique_24h)),
  Unique_lncRNAs = c(length(unique(shared_pairs$lncRNA_id)),
                     length(unique(unique_6h$lncRNA_id)),
                     length(unique(unique_24h$lncRNA_id))),
  Unique_mRNAs = c(length(unique(shared_pairs$mRNA_id)),
                   length(unique(unique_6h$mRNA_id)),
                   length(unique(unique_24h$mRNA_id)))
)
write_csv(summary_df, file.path(outdir, "comparative_summary_counts.csv"))

# === UPSET PLOT ===
set_list <- list(
  Shared = unique(paste(shared_pairs$lncRNA_id, shared_pairs$mRNA_id, sep = "_")),
  `6h only` = unique(paste(unique_6h$lncRNA_id, unique_6h$mRNA_id, sep = "_")),
  `24h only` = unique(paste(unique_24h$lncRNA_id, unique_24h$mRNA_id, sep = "_"))
)
png(file.path(outdir, "Plots/UpSet_interactions.png"), width = 1200, height = 900, res = 200)
upset(fromList(set_list), order.by = "freq", main.bar.color = "steelblue", sets.bar.color = "grey30")
dev.off()

# === HUBS IDENTIFICATION ===
get_hubs <- function(df, min_targets = 10) {
  df %>% group_by(lncRNA_id, lncRNA_symbol) %>%
    summarise(n_mRNAs = n(), .groups = "drop") %>%
    filter(n_mRNAs >= min_targets) %>%
    arrange(desc(n_mRNAs))
}
hubs_6h <- get_hubs(corr6)
hubs_24h <- get_hubs(corr24)
shared_hubs <- inner_join(hubs_6h, hubs_24h, by = "lncRNA_id") %>% select(lncRNA_id, lncRNA_symbol.x)
write_csv(hubs_6h, file.path(outdir, "lncRNA_hubs_6h.csv"))
write_csv(hubs_24h, file.path(outdir, "lncRNA_hubs_24h.csv"))
write_csv(shared_hubs, file.path(outdir, "lncRNA_hubs_shared.csv"))

# === GSEA CORE ENRICHMENT ===
parse_core <- function(df){
  df %>% separate_rows(core_enrichment, sep = "/") %>%
    rename(Gene = core_enrichment)
}
GSEA_BP6  <- read_csv(GSEA_BP_6h, show_col_types = FALSE)
GSEA_BP24 <- read_csv(GSEA_BP_24h, show_col_types = FALSE)
GSEA_KEGG6  <- read_csv(GSEA_KEGG_6h, show_col_types = FALSE)
GSEA_KEGG24 <- read_csv(GSEA_KEGG_24h, show_col_types = FALSE)

core_BP6  <- parse_core(GSEA_BP6)
core_BP24 <- parse_core(GSEA_BP24)
core_KEGG6  <- parse_core(GSEA_KEGG6)
core_KEGG24 <- parse_core(GSEA_KEGG24)

# === FUNCTIONAL SUMMARY TABLE ===
make_summary_table <- function(hub_table, corr_table, core_BP, core_KEGG, annot, label){
  summary_list <- list()
  for(i in 1:nrow(hub_table)){
    lncrna <- hub_table$lncRNA_symbol[i]
    id <- hub_table$lncRNA_id[i]
    targets <- corr_table %>% filter(lncRNA_id == id) %>% pull(mRNA_symbol)
    biotype <- annot %>% filter(gene_id == id) %>% pull(biotype_clean)
    isNovel <- ifelse(startsWith(id, "MSTRG"), "Novel", "Known")
    
    core_BP_hits <- core_BP %>% filter(Gene %in% targets)
    core_KEGG_hits <- core_KEGG %>% filter(Gene %in% targets)
    
    topBP <- paste(unique(head(core_BP_hits$Description, 3)), collapse = "; ")
    topKEGG <- paste(unique(head(core_KEGG_hits$Description, 3)), collapse = "; ")
    coreGenes <- paste(unique(c(core_BP_hits$Gene, core_KEGG_hits$Gene)), collapse = "/")
    
    summary_list[[i]] <- tibble(
      lncRNA = lncrna,
      Time = label,
      Biotype = biotype,
      Novelty = isNovel,
      n_mRNAs = length(targets),
      Example_mRNAs = paste(head(targets, 4), collapse = ", "),
      GO_BP = topBP,
      KEGG_Pathways = topKEGG,
      Core_genes = coreGenes
    )
  }
  bind_rows(summary_list)
}

summary6  <- make_summary_table(hubs_6h, corr6, core_BP6, core_KEGG6, annot, "6h")
summary24 <- make_summary_table(hubs_24h, corr24, core_BP24, core_KEGG24, annot, "24h")
summary_all <- bind_rows(summary6, summary24)
write_csv(summary_all, file.path(outdir, "lncRNA_hub_functional_summary.csv"))

# === BAR PLOT OF HUB COUNTS ===
hub_summary <- summary_all %>%
  group_by(Time, Novelty) %>%
  summarise(n = n(), .groups = "drop")

ggplot(hub_summary, aes(x = Time, y = n, fill = Novelty)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Novel" = "#e41a1c", "Known" = "#377eb8")) +
  theme_bw(base_size = 14) +
  labs(title = "Number of lncRNA Hubs per Timepoint", x = "", y = "Count")
ggsave(file.path(outdir, "Plots/Hubs_count_comparison.png"), width = 6, height = 5, dpi = 300)

message("✅ Comparative analysis complete! Results saved in: ", outdir)

# === TOP 10 MOST CONNECTED HUBS PER TIMEPOINT (FIXED) ===
top_hubs_6h <- hubs_6h %>%
  arrange(desc(n_mRNAs)) %>%
  slice_head(n = 10) %>%
  left_join(summary6 %>% select(lncRNA, GO_BP, KEGG_Pathways, Core_genes, Novelty), 
            by = c("lncRNA_symbol" = "lncRNA")) %>%
  mutate(Time = "6h") %>%
  select(Time, lncRNA_symbol, n_mRNAs, Novelty, GO_BP, KEGG_Pathways, Core_genes)

top_hubs_24h <- hubs_24h %>%
  arrange(desc(n_mRNAs)) %>%
  slice_head(n = 10) %>%
  left_join(summary24 %>% select(lncRNA, GO_BP, KEGG_Pathways, Core_genes, Novelty), 
            by = c("lncRNA_symbol" = "lncRNA")) %>%
  mutate(Time = "24h") %>%
  select(Time, lncRNA_symbol, n_mRNAs, Novelty, GO_BP, KEGG_Pathways, Core_genes)

top_hubs_combined <- bind_rows(top_hubs_6h, top_hubs_24h)

# === PRINT TO CONSOLE ===
cat("\n🔝 Top 10 lncRNA Hubs per Timepoint:\n")
print(top_hubs_combined %>% select(Time, lncRNA_symbol, n_mRNAs, Novelty) %>% arrange(Time, desc(n_mRNAs)))

# === EXPORT TO FILE ===
write_csv(top_hubs_combined, file.path(outdir, "Top10_lncRNA_Hubs_per_Timepoint.csv"))
cat("\n✅ Top hub summary saved in:", file.path(outdir, "Top10_lncRNA_Hubs_per_Timepoint.csv"), "\n")

####
# ===============================================================
# Separate DEGs by biotype using meta_ann + meta_novel annotation
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# === INPUT FILES ===
DE6_file <- "Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables/DE_results_NE6_vs_Ctrl.csv"
DE24_file <- "Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables/DE_results_NE24_vs_Ctrl.csv"
meta_ann_file <- "hisat2_6_24_meta_ann.csv"
meta_novel_file <- "hisat2_6_24_meta_novel_lncRNA.csv"
output_dir <- "Publication_Ready_Analysis_mRNA_lncRNA/1_DE_Tables/Biotype_Splits"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === LOAD DATA ===
DE6 <- read_csv(DE6_file, show_col_types = FALSE)
DE24 <- read_csv(DE24_file, show_col_types = FALSE)
meta_ann <- read_csv(meta_ann_file, show_col_types = FALSE)
meta_novel <- read_csv(meta_novel_file, show_col_types = FALSE)

# === BUILD ANNOTATION TABLE ===
# Combine known and novel lncRNAs
lnc_known_ids <- meta_ann %>%
  filter(gene_biotype == "lncRNA") %>%
  pull(ensembl_gene_id) %>%
  unique()

lnc_novel_ids <- meta_novel %>%
  pull(gene_id) %>%
  unique()

lnc_all_ids <- unique(c(lnc_known_ids, lnc_novel_ids))

# Extract protein-coding from meta_ann
pc_ids <- meta_ann %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(ensembl_gene_id) %>%
  unique()

# === Function to classify DE table ===
classify_DE <- function(df, lnc_ids, pc_ids) {
  df %>%
    mutate(
      gene_id_clean = gsub("\\..*$", "", gene_id),
      biotype_final = case_when(
        gene_id_clean %in% lnc_ids ~ "lncRNA",
        gene_id_clean %in% pc_ids ~ "protein_coding",
        startsWith(gene_id, "MSTRG") ~ "lncRNA",
        TRUE ~ "other"
      ),
      lncRNA_label = ifelse(biotype_final == "lncRNA" & startsWith(gene_id, "MSTRG"),
                            "lncRNA-MSTRG",
                            ifelse(biotype_final == "lncRNA", "lncRNA", NA))
    )
}

DE6_class <- classify_DE(DE6, lnc_all_ids, pc_ids)
DE24_class <- classify_DE(DE24, lnc_all_ids, pc_ids)

# === FILTER AND EXPORT ===
DE6_lnc <- DE6_class %>% filter(biotype_final == "lncRNA" & status != "Not Significant")
DE6_mrna <- DE6_class %>% filter(biotype_final == "protein_coding" & status != "Not Significant")

DE24_lnc <- DE24_class %>% filter(biotype_final == "lncRNA" & status != "Not Significant")
DE24_mrna <- DE24_class %>% filter(biotype_final == "protein_coding" & status != "Not Significant")

# === Save results ===
write_csv(DE6_lnc, file.path(output_dir, "DEG_lncRNA_6h.csv"))
write_csv(DE6_mrna, file.path(output_dir, "DEG_mRNA_6h.csv"))
write_csv(DE24_lnc, file.path(output_dir, "DEG_lncRNA_24h.csv"))
write_csv(DE24_mrna, file.path(output_dir, "DEG_mRNA_24h.csv"))

# === REPORT ===
cat("\n✅ DEG separation complete:\n")
cat("6h  → lncRNAs:", nrow(DE6_lnc), "| mRNAs:", nrow(DE6_mrna), "\n")
cat("24h → lncRNAs:", nrow(DE24_lnc), "| mRNAs:", nrow(DE24_mrna), "\n")

# Optional: quick summary table
summary_tbl <- data.frame(
  Timepoint = c("6h", "24h"),
  lncRNA = c(nrow(DE6_lnc), nrow(DE24_lnc)),
  mRNA = c(nrow(DE6_mrna), nrow(DE24_mrna))
)
print(summary_tbl)


