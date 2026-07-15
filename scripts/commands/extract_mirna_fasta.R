#!/usr/bin/env Rscript
# Helper: extract miRNA query FASTA files for IntaRNA (bash stage 12)
# Subsets miRBase mature sequences to the DE miRNA lists (UP/DOWN per timepoint)
# produced by script 02, writing miRNAs_{UP,DOWN}_{6h,24h}.fa. IntaRNA pairs
# lncRNA UP with miRNA DOWN and lncRNA DOWN with miRNA UP (reciprocal ceRNA).
#
# Ported from the miRNA-FASTA tail of script_miRNA_mRNA_6_24.R.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(Biostrings)
})

option_list <- list(
  make_option(c("--mirna-de-dir"), type = "character", default = NULL,
              help = "Directory with miRNA_DE_{6h,24h}_{UP,DOWN}.txt (from script 02)"),
  make_option(c("--mature-fa"), type = "character", default = NULL,
              help = "miRBase mature.fa (RNA FASTA)"),
  make_option(c("--species-prefix"), type = "character", default = "rno-",
              help = "miRBase species prefix to keep [default: %default]"),
  make_option(c("-o", "--output-dir"), type = "character", default = "results/11_mirge3_output",
              help = "Output directory (where bash stage 12 looks) [default: %default]")
)
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
if (is.null(opt$`mirna-de-dir`) || is.null(opt$`mature-fa`)) {
  print_help(parser); stop("\n[ERROR] Required: --mirna-de-dir, --mature-fa")
}
if (!dir.exists(opt$`output-dir`)) dir.create(opt$`output-dir`, recursive = TRUE)

fa <- readRNAStringSet(opt$`mature-fa`)
fa_tbl <- tibble(miRNA = word(names(fa), 1), sequence = as.character(fa)) %>%
  filter(str_detect(miRNA, paste0("^", opt$`species-prefix`)))
cat("[INFO] miRBase entries for prefix", opt$`species-prefix`, ":", nrow(fa_tbl), "\n")
fa_tbl <- fa_tbl %>% mutate(key = str_to_lower(miRNA))

for (tp in c("6h", "24h")) for (dir in c("UP", "DOWN")) {
  p <- file.path(opt$`mirna-de-dir`, paste0("miRNA_DE_", tp, "_", dir, ".txt"))
  if (!file.exists(p)) { cat("[SKIP] missing", basename(p), "\n"); next }
  df <- suppressWarnings(read_tsv(p, show_col_types = FALSE))
  mcol <- intersect(c("miRNA_id", "miRNA", "mirna"), names(df))[1]
  if (is.na(mcol)) { cat("[SKIP] no miRNA column in", basename(p), "\n"); next }
  ids <- str_to_lower(str_trim(as.character(df[[mcol]])))
  sub <- fa_tbl %>% filter(key %in% ids)
  out <- file.path(opt$`output-dir`, paste0("miRNAs_", dir, "_", tp, ".fa"))
  if (nrow(sub) == 0) { cat("[WARN] no sequence matches for", tp, dir, "\n"); next }
  writeLines(as.vector(rbind(paste0(">", sub$miRNA), sub$sequence)), out)
  cat("[OK]", basename(out), ":", nrow(sub), "miRNAs\n")
}
cat("[SUCCESS] miRNA query FASTA written to", opt$`output-dir`, "\n")
