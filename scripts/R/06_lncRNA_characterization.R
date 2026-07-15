#!/usr/bin/env Rscript
# Stage 06: Supplementary — lncRNA structural characterization from GTF
# Parses the merged annotated GTF (StringTie+gffcompare) to characterize
# annotated vs novel lncRNA transcripts: length, exon count, transcript class
# (via gffcompare class_code), isoform complexity, and inventory summaries,
# plus four supplementary panels and a text summary.
#
# Ported from 15_lncRNA_characterization_from_GTF.R with Windows paths removed
# and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option(c("-g", "--gtf"), type = "character", default = NULL,
              help = "Merged annotated GTF file from gffcompare"),
  make_option(c("-o", "--output-dir"), type = "character",
              default = "results/06_lncRNA_characterization",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 400,
              help = "PNG DPI [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
if (is.null(opt$gtf)) { print_help(parser); stop("\n[ERROR] Required argument: --gtf") }

GTF_FILE <- opt$gtf
OUTPUT_DIR <- opt$`output-dir`
DPI_PNG <- opt$dpi
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!file.exists(GTF_FILE)) stop("[ERROR] GTF file not found: ", GTF_FILE)

cat("[INFO] GTF file:", GTF_FILE, "\n")
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# GTF parsing helpers (manual parse preserves gffcompare class_code/cmp_ref)
# ---------------------------------------------------------------------------

read_gtf <- function(path) {
  cols <- c("seqname", "source_col", "feature", "start", "end",
            "score", "strand", "frame", "attribute")
  read_tsv(path, comment = "#", col_names = cols,
           col_types = cols(
             seqname = col_character(), source_col = col_character(),
             feature = col_character(), start = col_integer(), end = col_integer(),
             score = col_character(), strand = col_character(),
             frame = col_character(), attribute = col_character()),
           progress = FALSE)
}
extract_attr <- function(attr, key) str_match(attr, paste0(key, ' "([^"]+)"'))[, 2]
infer_source <- function(gene_id) case_when(
  str_detect(gene_id, "^MSTRG") ~ "novel",
  str_detect(gene_id, "^ENSRNOG") ~ "annotated",
  TRUE ~ "other"
)

cat("[INFO] Loading GTF...\n")
gtf <- read_gtf(GTF_FILE) %>%
  mutate(
    transcript_id = extract_attr(attribute, "transcript_id"),
    gene_id       = extract_attr(attribute, "gene_id"),
    gene_name     = extract_attr(attribute, "gene_name"),
    ref_gene_id   = extract_attr(attribute, "ref_gene_id"),
    cmp_ref       = extract_attr(attribute, "cmp_ref"),
    class_code    = extract_attr(attribute, "class_code"),
    xloc          = extract_attr(attribute, "xloc"),
    exon_number   = suppressWarnings(as.integer(extract_attr(attribute, "exon_number"))),
    source        = infer_source(gene_id)
  )

# ---------------------------------------------------------------------------
# Transcript- and gene-level characterization
# ---------------------------------------------------------------------------

transcripts <- gtf %>%
  filter(feature == "transcript", !is.na(transcript_id), !is.na(gene_id)) %>%
  transmute(seqname, strand, transcript_id, gene_id,
            gene_name = coalesce(gene_name, gene_id),
            ref_gene_id, cmp_ref, class_code, xloc, source,
            transcript_start = start, transcript_end = end,
            transcript_length = end - start + 1L) %>%
  distinct()

exons <- gtf %>%
  filter(feature == "exon", !is.na(transcript_id), !is.na(gene_id)) %>%
  transmute(transcript_id, gene_id, exon_start = start, exon_end = end,
            exon_length = end - start + 1L) %>%
  distinct()

exon_summary <- exons %>%
  group_by(transcript_id, gene_id) %>%
  summarise(exon_count = n(), summed_exon_length = sum(exon_length, na.rm = TRUE),
            .groups = "drop")

tx_tbl <- transcripts %>%
  left_join(exon_summary, by = c("transcript_id", "gene_id")) %>%
  mutate(
    exon_count = replace_na(exon_count, 1L),
    summed_exon_length = coalesce(summed_exon_length, transcript_length),
    transcript_class = case_when(
      source == "annotated" & !is.na(cmp_ref) & !is.na(class_code) & class_code == "=" ~ "known_reference_transcript",
      source == "annotated" & !is.na(cmp_ref) & !is.na(class_code) & class_code != "=" ~ "known_gene_alt_isoform",
      source == "annotated" & is.na(cmp_ref) ~ "annotated_gene_transcript_without_cmp_ref",
      source == "novel" ~ "novel_reconstructed_transcript",
      TRUE ~ "other"
    )
  )

gene_tbl <- tx_tbl %>%
  group_by(gene_id, source) %>%
  summarise(
    gene_name = first(gene_name),
    n_transcripts = n_distinct(transcript_id),
    min_tx_length = min(transcript_length, na.rm = TRUE),
    median_tx_length = median(transcript_length, na.rm = TRUE),
    max_tx_length = max(transcript_length, na.rm = TRUE),
    median_exon_count = median(exon_count, na.rm = TRUE),
    class_code_values = paste(sort(unique(na.omit(class_code))), collapse = ","),
    .groups = "drop"
  )

inventory_summary <- bind_rows(
  tx_tbl %>% count(source, name = "n_transcripts") %>%
    mutate(metric = "total_transcripts") %>% rename(value = n_transcripts),
  gene_tbl %>% count(source, name = "n_genes") %>%
    mutate(metric = "total_genes") %>% rename(value = n_genes),
  tx_tbl %>% count(source, transcript_class, name = "n") %>%
    mutate(metric = paste0("transcript_class__", transcript_class)) %>%
    select(source, metric, value = n)
)

isoform_summary <- gene_tbl %>%
  mutate(isoform_bin = case_when(
    n_transcripts == 1 ~ "1 transcript",
    n_transcripts == 2 ~ "2 transcripts",
    n_transcripts == 3 ~ "3 transcripts",
    n_transcripts >= 4 ~ "4+ transcripts",
    TRUE ~ "other"
  )) %>%
  count(source, isoform_bin, name = "n_genes")

write_tsv(tx_tbl, file.path(OUTPUT_DIR, "lncRNA_transcript_characterization_from_GTF.tsv"))
write_tsv(gene_tbl, file.path(OUTPUT_DIR, "lncRNA_gene_characterization_from_GTF.tsv"))
write_tsv(inventory_summary, file.path(OUTPUT_DIR, "lncRNA_inventory_summary.tsv"))
write_tsv(isoform_summary, file.path(OUTPUT_DIR, "lncRNA_isoform_complexity_summary.tsv"))

# ---------------------------------------------------------------------------
# Supplementary panels
# ---------------------------------------------------------------------------

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.title = element_text(face = "bold"))

p_len <- tx_tbl %>%
  filter(source %in% c("annotated", "novel")) %>%
  mutate(source = factor(source, levels = c("annotated", "novel"))) %>%
  ggplot(aes(x = source, y = transcript_length, fill = source)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.9, color = "grey30") +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Transcript length (nt)", title = "Annotated vs novel lncRNA transcript length") +
  plot_theme + theme(legend.position = "none")
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_transcript_length_violin.png"), p_len, width = 6.4, height = 5.1, dpi = DPI_PNG)
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_transcript_length_violin.pdf"), p_len, width = 6.4, height = 5.1)

p_exon <- tx_tbl %>%
  filter(source %in% c("annotated", "novel")) %>%
  count(source, exon_count, name = "n") %>%
  mutate(source = factor(source, levels = c("annotated", "novel"))) %>%
  ggplot(aes(x = exon_count, y = n, fill = source)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Exons per transcript", y = "# transcripts",
       title = "Exon count distribution of annotated vs novel lncRNA transcripts") +
  plot_theme
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_exon_count_distribution.png"), p_exon, width = 6.8, height = 5.1, dpi = DPI_PNG)
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_exon_count_distribution.pdf"), p_exon, width = 6.8, height = 5.1)

p_inv <- inventory_summary %>%
  filter(metric %in% c("total_transcripts", "total_genes")) %>%
  mutate(metric = recode(metric, total_transcripts = "Transcripts", total_genes = "Genes"),
         source = factor(source, levels = c("annotated", "novel"))) %>%
  ggplot(aes(x = metric, y = value, fill = source)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = scales::comma(value)), position = position_dodge(width = 0.72), vjust = -0.25, size = 3.5) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Count", title = "Annotated and novel lncRNA inventory from transcriptome reconstruction") +
  plot_theme
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_gene_transcript_inventory.png"), p_inv, width = 6.8, height = 5.1, dpi = DPI_PNG)
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_gene_transcript_inventory.pdf"), p_inv, width = 6.8, height = 5.1)

p_iso <- isoform_summary %>%
  mutate(source = factor(source, levels = c("annotated", "novel")),
         isoform_bin = factor(isoform_bin, levels = c("1 transcript", "2 transcripts", "3 transcripts", "4+ transcripts"))) %>%
  ggplot(aes(x = source, y = n_genes, fill = isoform_bin)) +
  geom_col(width = 0.72) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "# lncRNA genes", fill = "Isoform complexity", title = "Isoform complexity per lncRNA gene") +
  plot_theme
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_isoform_complexity.png"), p_iso, width = 6.8, height = 5.1, dpi = DPI_PNG)
ggsave(file.path(OUTPUT_DIR, "SuppFig2A_lncRNA_isoform_complexity.pdf"), p_iso, width = 6.8, height = 5.1)

# ---------------------------------------------------------------------------
# Text summary
# ---------------------------------------------------------------------------

med <- function(v) if (length(v) == 0) NA_real_ else round(median(v, na.rm = TRUE), 1)
text_lines <- c(
  "lncRNA transcriptome reconstruction summary", "",
  paste0("Annotated lncRNA genes: ", scales::comma(sum(gene_tbl$source == 'annotated'))),
  paste0("Novel lncRNA genes: ", scales::comma(sum(gene_tbl$source == 'novel'))),
  paste0("Annotated lncRNA transcripts: ", scales::comma(sum(tx_tbl$source == 'annotated'))),
  paste0("Novel lncRNA transcripts: ", scales::comma(sum(tx_tbl$source == 'novel'))), "",
  paste0("Annotated median transcript length (nt): ", med(tx_tbl$transcript_length[tx_tbl$source == 'annotated'])),
  paste0("Novel median transcript length (nt): ", med(tx_tbl$transcript_length[tx_tbl$source == 'novel'])),
  paste0("Annotated median exon count: ", med(tx_tbl$exon_count[tx_tbl$source == 'annotated'])),
  paste0("Novel median exon count: ", med(tx_tbl$exon_count[tx_tbl$source == 'novel'])), "",
  "Interpretation note:",
  "Because this is a transcriptome reconstruction framework, novel MSTRG entries are best described",
  "conservatively as novel reconstructed lncRNA loci/transcripts rather than definitively novel lncRNA",
  "genes unless orthogonal evidence supports locus-level novelty."
)
writeLines(text_lines, file.path(OUTPUT_DIR, "lncRNA_reconstruction_summary.txt"))

cat("[SUCCESS] lncRNA characterization complete!\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
