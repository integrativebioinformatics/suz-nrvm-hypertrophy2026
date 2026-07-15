suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"
ann_file <- file.path(base_dir, "hisat2_6_24_meta_ann.csv")
nov_file <- file.path(base_dir, "hisat2_6_24_meta_novel_lncRNA.csv")
out_dir  <- file.path(base_dir, "Supp_lncRNA_characterization_curatedCSV")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(ann_file), file.exists(nov_file))

save_gg <- function(plot, name, w = 7.5, h = 5.5, dpi = 300) {
  ggsave(file.path(out_dir, paste0(name, ".png")), plot = plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot = plot, width = w, height = h, bg = "white")
}

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

ann <- read.csv(ann_file, stringsAsFactors = FALSE, check.names = FALSE)
nov <- read.csv(nov_file, stringsAsFactors = FALSE, check.names = FALSE)

ann_lnc <- ann %>%
  mutate(
    is_lnc = grepl("lnc", gene_biotype, ignore.case = TRUE) |
             grepl("lnc", transcript_biotype, ignore.case = TRUE) |
             grepl("lnc", pred_biotype, ignore.case = TRUE)
  ) %>%
  filter(is_lnc) %>%
  transmute(
    source = "annotated",
    gene_id = ensembl_gene_id,
    transcript_id = ensembl_transcript_id,
    gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "", ensembl_gene_id, external_gene_name),
    transcript_length = suppressWarnings(as.numeric(transcript_length)),
    exon_count = suppressWarnings(as.numeric(num_exons)),
    class_code = class_code,
    transcript_class = case_when(
      !is.na(class_code) & class_code == "=" ~ "Reference transcript",
      !is.na(class_code) & class_code != "=" ~ "Alternative isoform",
      TRUE ~ "Annotated transcript"
    )
  ) %>%
  distinct(transcript_id, .keep_all = TRUE)

nov_lnc <- nov %>%
  filter(type == "transcript") %>%
  filter(grepl("lnc", pred_biotype, ignore.case = TRUE)) %>%
  transmute(
    source = "novel",
    gene_id = gene_id,
    transcript_id = transcript_id,
    gene_name = ifelse(is.na(qry_gene_id) | qry_gene_id == "", gene_id, qry_gene_id),
    transcript_length = suppressWarnings(as.numeric(len)),
    exon_count = suppressWarnings(as.numeric(num_exons)),
    class_code = class_code,
    transcript_class = "Novel reconstructed transcript"
  ) %>%
  distinct(transcript_id, .keep_all = TRUE)

tx_tbl <- bind_rows(ann_lnc, nov_lnc) %>%
  filter(!is.na(gene_id), !is.na(transcript_id)) %>%
  mutate(
    transcript_length = ifelse(is.na(transcript_length) | transcript_length <= 0, NA, transcript_length),
    exon_count = ifelse(is.na(exon_count) | exon_count <= 0, NA, exon_count),
    source = factor(source, levels = c("annotated", "novel"))
  )

gene_tbl <- tx_tbl %>%
  group_by(source, gene_id) %>%
  summarise(
    n_transcripts = n_distinct(transcript_id),
    gene_name = first(gene_name),
    .groups = "drop"
  ) %>%
  mutate(
    isoform_class = case_when(
      n_transcripts == 1 ~ "1 transcript",
      n_transcripts == 2 ~ "2 transcripts",
      n_transcripts == 3 ~ "3 transcripts",
      n_transcripts >= 4 ~ "4+ transcripts",
      TRUE ~ "other"
    ),
    isoform_class = factor(isoform_class, levels = c("1 transcript", "2 transcripts", "3 transcripts", "4+ transcripts"))
  )

inventory_tbl <- bind_rows(
  gene_tbl %>% count(source, name = "genes") %>% mutate(metric = "Genes", count = genes) %>% select(source, metric, count),
  tx_tbl %>% count(source, name = "transcripts") %>% mutate(metric = "Transcripts", count = transcripts) %>% select(source, metric, count)
)

write.csv(tx_tbl, file.path(out_dir, "lncRNA_curated_transcript_table.csv"), row.names = FALSE)
write.csv(gene_tbl, file.path(out_dir, "lncRNA_curated_gene_table.csv"), row.names = FALSE)
write.csv(inventory_tbl, file.path(out_dir, "lncRNA_inventory_summary.tsv"), row.names = FALSE)

# Summary text
ann_genes <- gene_tbl %>% filter(source == "annotated") %>% summarise(n = n()) %>% pull(n)
nov_genes <- gene_tbl %>% filter(source == "novel") %>% summarise(n = n()) %>% pull(n)
ann_tx <- tx_tbl %>% filter(source == "annotated") %>% summarise(n = n()) %>% pull(n)
nov_tx <- tx_tbl %>% filter(source == "novel") %>% summarise(n = n()) %>% pull(n)
ann_len_med <- tx_tbl %>% filter(source == "annotated") %>% summarise(m = median(transcript_length, na.rm = TRUE)) %>% pull(m)
nov_len_med <- tx_tbl %>% filter(source == "novel") %>% summarise(m = median(transcript_length, na.rm = TRUE)) %>% pull(m)
ann_ex_med <- tx_tbl %>% filter(source == "annotated") %>% summarise(m = median(exon_count, na.rm = TRUE)) %>% pull(m)
nov_ex_med <- tx_tbl %>% filter(source == "novel") %>% summarise(m = median(exon_count, na.rm = TRUE)) %>% pull(m)

txt <- c(
  "Curated lncRNA transcriptome reconstruction summary",
  "",
  paste0("Annotated lncRNA genes: ", comma(ann_genes)),
  paste0("Novel reconstructed lncRNA loci: ", comma(nov_genes)),
  paste0("Annotated lncRNA transcripts: ", comma(ann_tx)),
  paste0("Novel reconstructed lncRNA transcripts: ", comma(nov_tx)),
  "",
  paste0("Annotated median transcript length (nt): ", comma(round(ann_len_med, 1))),
  paste0("Novel median transcript length (nt): ", comma(round(nov_len_med, 1))),
  paste0("Annotated median exon count: ", comma(round(ann_ex_med, 1))),
  paste0("Novel median exon count: ", comma(round(nov_ex_med, 1))),
  "",
  "Interpretation note:",
  "These counts were derived from manually curated annotated and novel lncRNA metadata tables.",
  "Novel MSTRG entries should be described conservatively as novel reconstructed lncRNA loci/transcripts rather than definitive novel genes."
)
writeLines(txt, file.path(out_dir, "lncRNA_reconstruction_summary.txt"))

# Plot 1: transcript length (log10)
p_len <- tx_tbl %>%
  filter(!is.na(transcript_length), transcript_length > 0) %>%
  ggplot(aes(x = source, y = transcript_length, fill = source)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.9) +
  scale_y_log10(labels = comma_format()) +
  labs(
    title = "Annotated vs novel lncRNA transcript length",
    x = NULL,
    y = "Transcript length (nt, log10 scale)"
  )
save_gg(p_len, "SuppFig2A_lncRNA_transcript_length_log10", 7.5, 5.5)

# Plot 2: exon count binned
tx_exon_binned <- tx_tbl %>%
  filter(!is.na(exon_count), exon_count > 0) %>%
  mutate(
    exon_bin = case_when(
      exon_count == 1 ~ "1",
      exon_count == 2 ~ "2",
      exon_count >= 3 & exon_count <= 5 ~ "3-5",
      exon_count >= 6 & exon_count <= 10 ~ "6-10",
      exon_count > 10 ~ ">10",
      TRUE ~ "other"
    ),
    exon_bin = factor(exon_bin, levels = c("1", "2", "3-5", "6-10", ">10"))
  ) %>%
  count(source, exon_bin)

p_exon <- ggplot(tx_exon_binned, aes(x = exon_bin, y = n, fill = source)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title = "Exon count distribution of annotated vs novel lncRNA transcripts",
    x = "Exons per transcript",
    y = "# transcripts"
  )
save_gg(p_exon, "SuppFig2A_lncRNA_exon_count_binned", 7.5, 5.5)

# Plot 3: inventory
p_inventory <- ggplot(inventory_tbl, aes(x = metric, y = count, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = comma(count)),
    position = position_dodge(width = 0.8),
    vjust = -0.2,
    size = 4.2
  ) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Annotated and novel lncRNA inventory from curated metadata",
    x = NULL,
    y = "Count"
  )
save_gg(p_inventory, "SuppFig2A_lncRNA_gene_transcript_inventory_curated", 7.2, 5.2)

# Plot 4: isoform complexity as percentages
iso_pct <- gene_tbl %>%
  count(source, isoform_class) %>%
  group_by(source) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_iso <- ggplot(iso_pct, aes(x = source, y = pct, fill = isoform_class)) +
  geom_col() +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(
    title = "Isoform complexity per lncRNA gene",
    x = NULL,
    y = "% lncRNA genes",
    fill = "Isoform complexity"
  )
save_gg(p_iso, "SuppFig2A_lncRNA_isoform_complexity_percent", 7.2, 5.2)

# Plot 5: transcript classes
class_tbl <- tx_tbl %>%
  count(source, transcript_class)

p_class <- ggplot(class_tbl, aes(x = source, y = n, fill = transcript_class)) +
  geom_col() +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title = "Curated annotated and novel lncRNA transcript classes",
    x = NULL,
    y = "# transcripts",
    fill = "Transcript class"
  )
save_gg(p_class, "SuppFig2A_lncRNA_transcript_classes", 7.8, 5.2)
