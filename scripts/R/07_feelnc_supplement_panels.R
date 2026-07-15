#!/usr/bin/env Rscript
# Stage 07: Supplementary — FEELnc lncRNA classification panels
# Summarizes FEELnc classifier output (best relationship per lncRNA gene),
# plots the class distribution (annotated vs novel), and builds a selected
# module-candidate context table plus a text summary.
#
# Ported from 14_FEELnc_supplement_panels.R with Windows paths removed,
# the candidate list exposed as a parameter, and optparse added.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
})

option_list <- list(
  make_option(c("-f", "--feelnc-file"), type = "character", default = NULL,
              help = "FEELnc classifier output file"),
  make_option(c("--mrna-de-dir"), type = "character", default = NULL,
              help = "Directory with mRNA DE tables (for partner gene names)"),
  make_option(c("--candidates"), type = "character", default = NULL,
              help = "Comma-separated lncRNA gene IDs (or a file) for the candidate context table; defaults to the manuscript shortlist"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/07_feelnc",
              help = "Output directory [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 400,
              help = "PNG DPI [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
if (is.null(opt$`feelnc-file`)) { print_help(parser); stop("\n[ERROR] Required argument: --feelnc-file") }

FEELNC_FILE <- opt$`feelnc-file`
MRNA_DE_DIR <- opt$`mrna-de-dir`
OUTPUT_DIR  <- opt$`output-dir`
DPI_PNG     <- opt$dpi
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!file.exists(FEELNC_FILE)) stop("[ERROR] FEELnc file not found: ", FEELNC_FILE)

cat("[INFO] FEELnc file:", FEELNC_FILE, "\n")
cat("[INFO] Output directory:", OUTPUT_DIR, "\n\n")

# Default manuscript candidate shortlist (overridable via --candidates)
default_candidates <- c(
  "ENSRNOG00000090514", "MSTRG.8518", "MSTRG.3251", "MSTRG.17518",
  "ENSRNOG00000071598", "ENSRNOG00000085965", "MSTRG.8263",
  "ENSRNOG00000083943", "MSTRG.13291", "ENSRNOG00000090247",
  "ENSRNOG00000073751"
)
candidate_lnc <- default_candidates
if (!is.null(opt$candidates)) {
  raw_cand <- if (file.exists(opt$candidates)) read_lines(opt$candidates) else str_split(opt$candidates, ",")[[1]]
  candidate_lnc <- str_trim(raw_cand)
  candidate_lnc <- candidate_lnc[candidate_lnc != ""]
}

# ---------------------------------------------------------------------------
# 1. Load FEELnc classifier output, keep best relationship per lncRNA
# ---------------------------------------------------------------------------

cat("[INFO] Loading FEELnc output...\n")
x <- read_tsv(FEELNC_FILE, show_col_types = FALSE) %>%
  filter(isBest == 1) %>%
  mutate(
    source = if_else(str_detect(lncRNA_gene, "^MSTRG"), "novel", "annotated"),
    class_simple = case_when(
      type == "intergenic" & location == "upstream"   ~ "intergenic_upstream",
      type == "intergenic" & location == "downstream" ~ "intergenic_downstream",
      direction == "sense"     & location == "exonic"   ~ "sense_exonic",
      direction == "sense"     & location == "intronic" ~ "sense_intronic",
      direction == "antisense" & location == "exonic"   ~ "antisense_exonic",
      direction == "antisense" & location == "intronic" ~ "antisense_intronic",
      TRUE ~ "other"
    )
  )

by_gene <- x %>% distinct(lncRNA_gene, source, class_simple, .keep_all = TRUE)

class_counts <- by_gene %>%
  count(source, class_simple, name = "n_lncRNA_genes") %>%
  arrange(source, desc(n_lncRNA_genes))
write_tsv(class_counts, file.path(OUTPUT_DIR, "SuppFig2B_FEELnc_class_counts.tsv"))
write_tsv(by_gene, file.path(OUTPUT_DIR, "FEELnc_best_only_unique_lncRNA_gene.tsv"))

# ---------------------------------------------------------------------------
# 2. Class distribution panel (SuppFig2B)
# ---------------------------------------------------------------------------

class_order <- class_counts %>%
  group_by(class_simple) %>% summarise(total = sum(n_lncRNA_genes), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(class_simple)

class_counts_ordered <- class_counts %>%
  mutate(source = factor(source, levels = c("annotated", "novel")),
         class_simple = factor(class_simple, levels = class_order))

p <- ggplot(class_counts_ordered, aes(x = source, y = n_lncRNA_genes, fill = class_simple)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.2) +
  labs(x = NULL, y = "# unique lncRNA genes", fill = "FEELnc class",
       title = "lncRNA classification by FEELnc") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
ggsave(file.path(OUTPUT_DIR, "SuppFig2B_FEELnc_classification_unique_genes.png"), p, width = 7.2, height = 5.6, dpi = DPI_PNG)
ggsave(file.path(OUTPUT_DIR, "SuppFig2B_FEELnc_classification_unique_genes.pdf"), p, width = 7.2, height = 5.6)

# ---------------------------------------------------------------------------
# 3. Selected module-candidate context table (SuppFig2C)
# ---------------------------------------------------------------------------

mrna_map <- tibble(gene_id = character(), gene_name = character())
if (!is.null(MRNA_DE_DIR) && dir.exists(MRNA_DE_DIR)) {
  mrna_files <- list.files(MRNA_DE_DIR, pattern = "mRNA_DE\\.txt$", full.names = TRUE, ignore.case = TRUE)
  if (length(mrna_files) > 0) {
    mrna_map <- bind_rows(lapply(mrna_files, function(f) {
      df <- read_tsv(f, show_col_types = FALSE)
      if (all(c("gene_id", "gene_name") %in% names(df))) select(df, gene_id, gene_name) else NULL
    })) %>% distinct()
  }
}

has_partner <- "partnerRNA_gene" %in% names(by_gene)
cand_tbl <- by_gene %>% filter(lncRNA_gene %in% candidate_lnc)
if (has_partner) {
  cand_tbl <- cand_tbl %>%
    left_join(mrna_map, by = c("partnerRNA_gene" = "gene_id")) %>%
    mutate(partner_label = coalesce(gene_name, partnerRNA_gene))
}
avail_cols <- intersect(
  c("lncRNA_gene", "lncRNA_transcript", "source", "class_simple", "type", "direction",
    "subtype", "location", "distance", "partnerRNA_gene", "partnerRNA_transcript",
    "gene_name", "partner_label"),
  names(cand_tbl)
)
cand_tbl <- cand_tbl %>% select(all_of(avail_cols)) %>%
  arrange(across(any_of(c("lncRNA_gene", "distance", "partner_label"))))
write_tsv(cand_tbl, file.path(OUTPUT_DIR, "SuppFig2C_FEELnc_selected_module_candidates.tsv"))

# ---------------------------------------------------------------------------
# 4. Text summary
# ---------------------------------------------------------------------------

class_lines <- class_counts %>%
  mutate(line = paste0(source, " | ", class_simple, " = ", n_lncRNA_genes)) %>% pull(line)
cand_lines <- if (nrow(cand_tbl) > 0 && "class_simple" %in% names(cand_tbl)) {
  cand_tbl %>%
    mutate(line = paste0(lncRNA_gene, " [", class_simple, "]")) %>% pull(line)
} else character(0)

summary_lines <- c(
  "FEELnc summary for supplementary lncRNA characterization", "",
  paste0("Best relationships retained: ", nrow(x)),
  paste0("Unique lncRNA genes retained after collapsing isBest=1: ", nrow(by_gene)), "",
  "Counts by simplified class:", class_lines, "",
  "Selected module candidates with FEELnc context:", cand_lines
)
writeLines(summary_lines, file.path(OUTPUT_DIR, "Supp_lncRNA_characterization_summary.txt"))

cat("[SUCCESS] FEELnc supplement panels complete!\n")
cat("[OUTPUT]", OUTPUT_DIR, "\n")
