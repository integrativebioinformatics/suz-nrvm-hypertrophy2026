library(tidyverse)

base_dir <- "C:/Users/sebau/OneDrive/Documentos/Tesis_Doc_SU-Z/Script_Analysis/NE_trancriptome_analysis/Salmon_Quantification_Analysis/NE6_24_analysis"
out_dir <- file.path(base_dir, "Supp_lncRNA_characterization")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

find_first <- function(candidates) {
  hits <- candidates[file.exists(file.path(base_dir, candidates))]
  if (length(hits) == 0) return(NA_character_)
  file.path(base_dir, hits[1])
}

feelnc_file <- find_first(c(
  "lncRNA_classes_CNE624_rn8_8_annotated_1kb_100kb.txt"
))

if (is.na(feelnc_file)) stop("FEELnc file not found in working directory.")

x <- read_tsv(feelnc_file, show_col_types = FALSE) %>%
  filter(isBest == 1) %>%
  mutate(
    source = if_else(str_detect(lncRNA_gene, "^MSTRG"), "novel", "annotated"),
    class_simple = case_when(
      type == "intergenic" & location == "upstream" ~ "intergenic_upstream",
      type == "intergenic" & location == "downstream" ~ "intergenic_downstream",
      direction == "sense" & location == "exonic" ~ "sense_exonic",
      direction == "sense" & location == "intronic" ~ "sense_intronic",
      direction == "antisense" & location == "exonic" ~ "antisense_exonic",
      direction == "antisense" & location == "intronic" ~ "antisense_intronic",
      TRUE ~ "other"
    )
  )

# Collapse at unique lncRNA_gene level for manuscript-level structural summary
by_gene <- x %>%
  distinct(lncRNA_gene, source, class_simple, .keep_all = TRUE)

class_counts <- by_gene %>%
  count(source, class_simple, name = "n_lncRNA_genes") %>%
  arrange(source, desc(n_lncRNA_genes))

write_tsv(class_counts, file.path(out_dir, "SuppFig2B_FEELnc_class_counts_unique_lncRNA_gene.tsv"))
write_tsv(by_gene, file.path(out_dir, "FEELnc_best_only_unique_lncRNA_gene.tsv"))

class_order <- class_counts %>%
  group_by(class_simple) %>%
  summarise(total = sum(n_lncRNA_genes), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(class_simple)

class_counts <- class_counts %>%
  mutate(
    source = factor(source, levels = c("annotated", "novel")),
    class_simple = factor(class_simple, levels = class_order)
  )

p <- ggplot(class_counts, aes(x = source, y = n_lncRNA_genes, fill = class_simple)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.2) +
  labs(
    x = NULL,
    y = "# unique lncRNA genes",
    fill = "FEELnc class"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "SuppFig2B_FEELnc_classification_unique_genes.png"), p, width = 7.2, height = 5.6, dpi = 400)
ggsave(file.path(out_dir, "SuppFig2B_FEELnc_classification_unique_genes.pdf"), p, width = 7.2, height = 5.6)

# Selected module candidates for discussion / cis-context summary
candidate_lnc <- c(
  "ENSRNOG00000090514", "MSTRG.8518", "MSTRG.3251", "MSTRG.17518",
  "ENSRNOG00000071598", "ENSRNOG00000085965", "MSTRG.8263",
  "ENSRNOG00000083943", "MSTRG.13291", "ENSRNOG00000090247",
  "ENSRNOG00000073751"
)

mrna_files <- c(find_first(c("NE_vs_Ctrl_6h_mRNA_DE.txt", "NE_vs_Ctrl_6h_mRNA_DE.txt")),
                find_first(c("NE_vs_Ctrl_24h_mRNA_DE.txt", "NE_vs_Ctrl_24h_mRNA_DE.txt")))
mrna_files <- mrna_files[!is.na(mrna_files)]

mrna_map <- if (length(mrna_files) > 0) {
  bind_rows(lapply(mrna_files, function(f) read_tsv(f, show_col_types = FALSE))) %>%
    select(gene_id, gene_name) %>%
    distinct()
} else {
  tibble(gene_id = character(), gene_name = character())
}

cand_tbl <- by_gene %>%
  filter(lncRNA_gene %in% candidate_lnc) %>%
  left_join(mrna_map, by = c("partnerRNA_gene" = "gene_id")) %>%
  mutate(partner_label = coalesce(gene_name, partnerRNA_gene)) %>%
  select(lncRNA_gene, lncRNA_transcript, source, class_simple, type, direction, subtype, location, distance,
         partnerRNA_gene, partnerRNA_transcript, gene_name, partner_label) %>%
  arrange(lncRNA_gene, distance, partner_label)

write_tsv(cand_tbl, file.path(out_dir, "SuppFig2C_FEELnc_selected_module_candidates_best_context.tsv"))

# Optional: create a very simple text summary for the discussion
summary_lines <- c(
  "FEELnc summary for supplementary lncRNA characterization",
  "",
  paste0("Best relationships retained: ", nrow(x)),
  paste0("Unique lncRNA genes retained after collapsing isBest=1: ", nrow(by_gene)),
  "",
  "Counts by simplified class:",
  paste(apply(class_counts %>% mutate(line = paste0(source, " | ", class_simple, " = ", n_lncRNA_genes)) %>% select(line), 1, identity), collapse = "\n"),
  "",
  "Selected module candidates with FEELnc context:",
  paste(apply(cand_tbl %>% mutate(line = paste0(lncRNA_gene, " -> ", partner_label, " [", class_simple, "; ", subtype, "; ", location, "; distance=", distance, "]")) %>% select(line), 1, identity), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "Supp_lncRNA_characterization_summary.txt"))

