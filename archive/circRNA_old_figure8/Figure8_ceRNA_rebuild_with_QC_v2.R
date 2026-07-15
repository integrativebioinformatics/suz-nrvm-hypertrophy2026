
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
})

base_dir <- "."
out_dir <- file.path(base_dir, "Paper_ceRNA_candidates")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

energy_high <- -20
energy_moderate <- -15
energy_exploratory <- -10
keep_only_moderate_or_better <- FALSE
interaction_files_manual <- NULL

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

find_file <- function(filename, base = ".") {
  hits <- list.files(base, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == filename]
  if (length(hits) == 0) stop("File not found: ", filename)
  hits[1]
}

find_optional_file <- function(filename, base = ".") {
  hits <- list.files(base, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == filename]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

find_all_interaction_files <- function(base = ".") {
  list.files(base, recursive = TRUE, full.names = TRUE, pattern = "^Interact_.*\\.csv$")
}

read_delim_auto <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) == 0) return(read_csv(path, show_col_types = FALSE))
  if (str_detect(first_line, ";") && !str_detect(first_line, ",")) {
    return(read_delim(path, delim = ";", show_col_types = FALSE, trim_ws = TRUE, escape_double = FALSE))
  }
  if (str_detect(first_line, "\\t")) return(read_tsv(path, show_col_types = FALSE))
  read_csv(path, show_col_types = FALSE)
}

normalize_key <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("[\"'`]", "") %>%
    str_to_lower()
}

normalize_mirna <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower()
}

canonicalize_circ_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_

  # remove optional host / suffix pieces after a pipe-with-spaces pattern
  x <- ifelse(!is.na(x), str_replace(x, "\\s*\\|\\s*[^|]+$", function(m) m), x)
  out <- x

  # Case 1: already chr:start|end, optionally followed by host label after ' | '
  has_coord_pipe <- !is.na(x) & str_detect(x, "^[^:|_]+:[0-9]+\\|[0-9]+")
  out[has_coord_pipe] <- str_extract(x[has_coord_pipe], "^[^:|_]+:[0-9]+\\|[0-9]+")

  # Case 2: compact chr_start_end, optionally with extra suffixes
  compact_pat <- "^([^:|_\\s]+)_([0-9]+)_([0-9]+)(?:$|_.+)"
  has_compact <- !is.na(x) & !has_coord_pipe & str_detect(x, compact_pat)
  out[has_compact] <- str_replace(x[has_compact], compact_pat, "\\1:\\2|\\3")

  out
}

choose_col <- function(df, candidates, fallback = NULL) {
  hit <- candidates[candidates %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) return(fallback)
  hit
}

choose_vec <- function(df, candidates, fallback = NA_character_) {
  hit <- choose_col(df, candidates, NULL)
  if (is.null(hit)) return(rep(fallback, nrow(df)))
  as.character(df[[hit]])
}

coalesce_chr <- function(...) {
  vals <- list(...)
  if (length(vals) == 0) return(character())
  out <- vals[[1]]
  if (length(vals) > 1) {
    for (i in 2:length(vals)) out <- dplyr::coalesce(out, vals[[i]])
  }
  out
}

first_nonempty_col <- function(df, candidates) {
  for (nm in candidates) {
    if (!nm %in% names(df)) next
    v <- as.character(df[[nm]])
    if (any(!is.na(v) & str_trim(v) != "")) return(nm)
  }
  NA_character_
}

safe_write_lines <- function(lines, path) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

write_df_txt <- function(df, path, title = NULL) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  if (!is.null(title)) {
    writeLines(title, con)
    writeLines("", con)
  }
  if (is.null(df) || nrow(df) == 0) {
    writeLines("No rows available.", con)
  } else {
    utils::write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

append_lines <- function(lines, ...) c(lines, unlist(list(...), use.names = FALSE))

collapse_labels <- function(x, n = 20) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & str_trim(x) != ""]
  if (length(x) == 0) return("none")
  paste(head(x, n), collapse = "; ")
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

format_top_features <- function(df, n = 20, label_col, score_col = NULL, aux_col = NULL) {
  if (is.null(df) || nrow(df) == 0 || !label_col %in% names(df)) return("none")
  sub <- df %>% filter(!is.na(.data[[label_col]]), str_trim(as.character(.data[[label_col]])) != "")
  if (!is.null(score_col) && score_col %in% names(sub)) {
    sub <- sub %>% arrange(desc(abs(suppressWarnings(as.numeric(.data[[score_col]])))))
  }
  sub <- head(sub, n)
  out <- as.character(sub[[label_col]])
  if (!is.null(score_col) && score_col %in% names(sub)) {
    out <- paste0(out, " (", score_col, "=", fmt_num(suppressWarnings(as.numeric(sub[[score_col]]))), ")")
  }
  if (!is.null(aux_col) && aux_col %in% names(sub)) {
    out <- paste0(out, "; ", aux_col, "=", fmt_num(suppressWarnings(as.numeric(sub[[aux_col]]))))
  }
  paste(out, collapse = "; ")
}

canonical_direction <- function(x, logfc = NULL) {
  y <- str_to_upper(str_trim(as.character(x)))
  n <- length(y)
  pos_fc <- rep(FALSE, n)
  neg_fc <- rep(FALSE, n)
  if (!is.null(logfc)) {
    logfc <- suppressWarnings(as.numeric(logfc))
    if (length(logfc) == 1 && n > 1) logfc <- rep(logfc, n)
    if (length(logfc) == n) {
      pos_fc <- !is.na(logfc) & logfc > 0
      neg_fc <- !is.na(logfc) & logfc < 0
    }
  }
  case_when(
    y %in% c("UP", "UPREGULATED", "POSITIVE") ~ "UP",
    y %in% c("DOWN", "DOWNREGULATED", "NEGATIVE") ~ "DOWN",
    pos_fc ~ "UP",
    neg_fc ~ "DOWN",
    TRUE ~ NA_character_
  )
}

infer_phase_two_timepoints <- function(logfc6, padj6, logfc24, padj24, alpha = 0.05, lfc_cut = 0) {
  sig6 <- !is.na(padj6) & padj6 <= alpha & abs(logfc6) > lfc_cut
  sig24 <- !is.na(padj24) & padj24 <= alpha & abs(logfc24) > lfc_cut
  sign6 <- sign(logfc6)
  sign24 <- sign(logfc24)
  case_when(
    sig6 & !sig24 ~ "Early",
    !sig6 & sig24 ~ "Late",
    sig6 & sig24 & sign6 == sign24 ~ "Sustained",
    sig6 & sig24 & sign6 != sign24 ~ "Reversal",
    TRUE ~ "Unassigned"
  )
}

energy_tier <- function(E) {
  case_when(
    is.na(E) ~ "Unscored",
    E <= energy_high ~ "High",
    E > energy_high & E <= energy_moderate ~ "Moderate",
    E > energy_moderate & E <= energy_exploratory ~ "Exploratory",
    TRUE ~ "Weak"
  )
}

rank_tier <- function(x) recode(x, High = 1L, Moderate = 2L, Exploratory = 3L, Weak = 4L, Unscored = 5L, .default = 99L)

read_intarna_file <- function(path) {
  raw <- read_delim_auto(path)
  if (ncol(raw) == 1 && str_detect(names(raw)[1], "id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E")) {
    raw <- raw %>%
      rename(raw_col = 1) %>%
      separate_wider_delim(
        raw_col,
        delim = ";",
        names = c("id1", "start1", "end1", "id2", "start2", "end2", "subseqDP", "hybridDP", "E"),
        too_few = "align_start"
      )
  }
  if (!all(c("id1", "id2", "E") %in% names(raw))) stop("IntaRNA file not recognized: ", basename(path))
  fn <- basename(path)
  raw %>%
    mutate(
      source_file = fn,
      ncrna_type = case_when(
        str_detect(str_to_lower(fn), "circ") ~ "circRNA",
        str_detect(str_to_lower(fn), "lnc") ~ "lncRNA",
        TRUE ~ "unknown"
      ),
      timepoint_file = case_when(
        str_detect(fn, regex("24h", ignore_case = TRUE)) ~ "24h",
        str_detect(fn, regex("6h", ignore_case = TRUE)) ~ "6h",
        TRUE ~ NA_character_
      ),
      expected_ncrna_direction = case_when(
        str_detect(fn, regex("UP", ignore_case = TRUE)) ~ "UP",
        str_detect(fn, regex("DOWN", ignore_case = TRUE)) ~ "DOWN",
        TRUE ~ NA_character_
      ),
      expected_miRNA_direction = case_when(
        str_detect(fn, regex("miUP", ignore_case = TRUE)) ~ "UP",
        str_detect(fn, regex("miDOWN", ignore_case = TRUE)) ~ "DOWN",
        TRUE ~ NA_character_
      ),
      id1 = as.character(id1),
      miRNA = as.character(id2),
      miRNA_key = normalize_mirna(id2),
      E = suppressWarnings(as.numeric(E)),
      interaction_tier = energy_tier(E)
    )
}

parse_lnc_intarna_id1 <- function(x) {
  x <- as.character(x)
  parts <- str_split_fixed(x, "\\|", 6)
  tibble(
    id1_raw = x,
    transcript_id_from_id1 = na_if(parts[, 1], ""),
    gene_id_from_id1 = na_if(parts[, 2], ""),
    phase_from_id1 = na_if(parts[, 3], ""),
    source_from_id1 = na_if(parts[, 4], ""),
    direction_from_id1 = str_to_upper(na_if(parts[, 5], "")),
    timepoint_from_id1 = na_if(parts[, 6], ""),
    tx_key_from_id1 = normalize_key(na_if(parts[, 1], "")),
    gene_key_from_id1 = normalize_key(na_if(parts[, 2], ""))
  )
}

parse_circ_intarna_id1 <- function(x) {
  x <- as.character(x)
  parts <- str_split_fixed(x, "\\|", 5)
  tibble(
    id1_raw = x,
    timepoint_from_id1 = na_if(parts[, 1], ""),
    expected_direction_from_id1 = str_to_upper(na_if(parts[, 2], "")),
    circ_compact = na_if(parts[, 3], ""),
    circ_compact_host = na_if(parts[, 4], ""),
    circ_len_from_id1 = suppressWarnings(as.integer(na_if(parts[, 5], "")))
  ) %>%
    mutate(
      circ_id_from_id1 = if_else(
        str_detect(circ_compact, "^[^_]+_[0-9]+_[0-9]+$"),
        str_replace(circ_compact, "^([^_]+)_([0-9]+)_([0-9]+)$", "\\1:\\2|\\3"),
        NA_character_
      ),
      host_label_from_id1 = if_else(
        !is.na(circ_compact_host) & str_detect(circ_compact_host, "^[^_]+_[0-9]+_[0-9]+_.+"),
        str_replace(circ_compact_host, "^[^_]+_[0-9]+_[0-9]+_", ""),
        NA_character_
      ),
      circRNA_ID2_from_id1 = if_else(
        !is.na(circ_id_from_id1) & !is.na(host_label_from_id1),
        paste0(circ_id_from_id1, " | ", host_label_from_id1),
        NA_character_
      ),
      circRNA_id_AS_from_id1 = if_else(
        !is.na(circ_id_from_id1) & !is.na(circ_len_from_id1),
        paste0(circ_id_from_id1, "_", circ_len_from_id1),
        NA_character_
      ),
      circRNA_DE_id_from_id1 = if_else(
        !is.na(circRNA_ID2_from_id1) & !is.na(expected_direction_from_id1) & !is.na(circ_len_from_id1),
        paste0(circRNA_ID2_from_id1, "_", str_to_title(str_to_lower(expected_direction_from_id1)), "_", circ_len_from_id1),
        NA_character_
      )
    )
}

read_lnc_annotation_map <- function(base = ".") {
  mapped_file <- find_optional_file("lncRNA_DE_with_phase_and_transcriptID.csv", base)
  if (!is.na(mapped_file)) {
    x <- read_csv(mapped_file, show_col_types = FALSE)
    out <- x %>%
      mutate(
        gene_id = coalesce_chr(choose_vec(x, c("gene_id", "gene_name")), choose_vec(x, c("gene_name"))),
        transcript_id = choose_vec(x, c("transcript_id")),
        gene_name = choose_vec(x, c("gene_name")),
        phase = str_to_title(choose_vec(x, c("phase"))),
        timepoint = choose_vec(x, c("timepoint")),
        direction = canonical_direction(choose_vec(x, c("direction")), choose_vec(x, c("log2FoldChange", "log2FC"))),
        source = coalesce_chr(choose_vec(x, c("source")), rep("annotated", nrow(x))),
        ncrna_type = "lncRNA",
        ncrna_id_raw = coalesce_chr(transcript_id, gene_id),
        ncrna_id_key = normalize_key(coalesce_chr(transcript_id, gene_id)),
        gene_key = normalize_key(gene_id),
        tx_key = normalize_key(transcript_id),
        ncRNA_label = gene_id
      ) %>%
      select(ncrna_type, ncrna_id_raw, ncrna_id_key, gene_id, transcript_id, gene_name, source, phase, timepoint, direction, ncRNA_label, gene_key, tx_key) %>%
      distinct(tx_key, timepoint, .keep_all = TRUE)
    return(out)
  }

  phase_files <- tribble(
    ~phase, ~source, ~file,
    "Early", "annotated", find_optional_file("early_lncRNA.csv", base),
    "Sustained", "annotated", find_optional_file("sustained_lncRNA.csv", base),
    "Late", "annotated", find_optional_file("late_lncRNA.csv", base),
    "Early", "novel", find_optional_file("early_novel_lncRNA.csv", base),
    "Sustained", "novel", find_optional_file("sustained_novel_lncRNA.csv", base),
    "Late", "novel", find_optional_file("late_novel_lncRNA.csv", base)
  )

  map_tbl <- pmap_dfr(phase_files, function(phase, source, file) {
    if (is.na(file)) return(tibble())
    z <- read_csv(file, show_col_types = FALSE)
    tibble(
      gene_id = choose_vec(z, c("gene_id", "gene_name")),
      transcript_id = choose_vec(z, c("transcript_id")),
      gene_name = choose_vec(z, c("gene_name")),
      phase = phase,
      source = source
    )
  }) %>% distinct()

  de6_file <- find_optional_file("NE_vs_Ctrl_6h_lncRNA_DE.txt", base)
  de24_file <- find_optional_file("NE_vs_Ctrl_24h_lncRNA_DE.txt", base)
  de_tbl <- bind_rows(
    if (!is.na(de6_file)) read_tsv(de6_file, show_col_types = FALSE) %>% mutate(timepoint = "6h") else tibble(),
    if (!is.na(de24_file)) read_tsv(de24_file, show_col_types = FALSE) %>% mutate(timepoint = "24h") else tibble()
  )

  if (nrow(map_tbl) == 0 && nrow(de_tbl) == 0) {
    return(tibble(
      ncrna_type = character(), ncrna_id_raw = character(), ncrna_id_key = character(),
      gene_id = character(), transcript_id = character(), gene_name = character(),
      source = character(), phase = character(), timepoint = character(),
      direction = character(), ncRNA_label = character(), gene_key = character(),
      tx_key = character()
    ))
  }

  de_std <- de_tbl %>%
    mutate(
      gene_id = choose_vec(de_tbl, c("gene_id", "gene_name")),
      transcript_id = choose_vec(de_tbl, c("transcript_id")),
      gene_name = choose_vec(de_tbl, c("gene_name")),
      direction = canonical_direction(choose_vec(de_tbl, c("direction")), choose_vec(de_tbl, c("log2FoldChange", "log2FC"))),
      source = coalesce_chr(choose_vec(de_tbl, c("source")), rep("annotated", nrow(de_tbl)))
    ) %>%
    select(gene_id, transcript_id, gene_name, direction, source, timepoint)

  map_tbl %>%
    left_join(de_std, by = c("gene_id", "transcript_id", "source")) %>%
    mutate(
      ncrna_type = "lncRNA",
      ncrna_id_raw = coalesce_chr(transcript_id, gene_id),
      ncrna_id_key = normalize_key(coalesce_chr(transcript_id, gene_id)),
      gene_key = normalize_key(gene_id),
      tx_key = normalize_key(transcript_id),
      ncRNA_label = gene_id,
      timepoint = coalesce_chr(timepoint, rep(NA_character_, n()))
    ) %>%
    select(ncrna_type, ncrna_id_raw, ncrna_id_key, gene_id, transcript_id, gene_name, source, phase, timepoint, direction, ncRNA_label, gene_key, tx_key) %>%
    distinct(tx_key, timepoint, .keep_all = TRUE)
}

read_circ_annotation_map <- function(base = ".") {
  de6_file <- find_optional_file("06_circRNA_sequences_de6_DE_complete_table.txt", base)
  de24_file <- find_optional_file("07_circRNA_sequences_de24_DE_complete_table.txt", base)
  if (is.na(de6_file) && is.na(de24_file)) stop("Need at least one circRNA DE table.")

  read_circ_de <- function(path, tp) {
    if (is.na(path)) return(tibble())
    x <- read_delim_auto(path)

    circ_col <- choose_col(x, c("circ_id", "circRNA_ID", "circRNA_ID2"), names(x)[1])
    id2_col  <- choose_col(x, c("circRNA_ID2"))
    len_col  <- choose_col(x, c("length"))
    deid_col <- choose_col(x, c("circRNA_DE_id"))
    as_col   <- choose_col(x, c("circRNA_id_AS"))
    host_col <- choose_col(x, c("host_label", "host_gene", "gene_symbol", "genesymbol"))
    reg_col  <- choose_col(x, c("regulations", "direction"))
    lfc_col  <- choose_col(x, c("log2FC", "log2FoldChange"))
    p_col    <- choose_col(x, c("padj", "pvalue"))

    circ_raw <- if (!is.null(circ_col)) as.character(x[[circ_col]]) else rep(NA_character_, nrow(x))
    id2_raw  <- if (!is.null(id2_col)) as.character(x[[id2_col]]) else rep(NA_character_, nrow(x))
    as_raw   <- if (!is.null(as_col)) as.character(x[[as_col]]) else rep(NA_character_, nrow(x))
    deid_raw <- if (!is.null(deid_col)) as.character(x[[deid_col]]) else rep(NA_character_, nrow(x))
    host_raw <- if (!is.null(host_col)) as.character(x[[host_col]]) else rep(NA_character_, nrow(x))

    tibble(
      timepoint = tp,
      circ_id = coalesce_chr(
        canonicalize_circ_id(circ_raw),
        canonicalize_circ_id(id2_raw),
        canonicalize_circ_id(as_raw),
        canonicalize_circ_id(deid_raw)
      ),
      circRNA_ID2 = if (!is.null(id2_col)) as.character(x[[id2_col]]) else NA_character_,
      circRNA_id_AS = if (!is.null(as_col)) as.character(x[[as_col]]) else NA_character_,
      circRNA_DE_id = if (!is.null(deid_col)) as.character(x[[deid_col]]) else NA_character_,
      length = suppressWarnings(as.numeric(if (!is.null(len_col)) x[[len_col]] else NA_real_)),
      direction = canonical_direction(if (!is.null(reg_col)) x[[reg_col]] else NA_character_, if (!is.null(lfc_col)) x[[lfc_col]] else NA_real_),
      log2FC = suppressWarnings(as.numeric(if (!is.null(lfc_col)) x[[lfc_col]] else NA_real_)),
      pval_like = suppressWarnings(as.numeric(if (!is.null(p_col)) x[[p_col]] else NA_real_)),
      host_label_raw = host_raw
    ) %>%
      mutate(
        host_label = coalesce_chr(
          host_label_raw,
          if_else(!is.na(circRNA_ID2) & str_detect(circRNA_ID2, "\|"), str_trim(str_replace(circRNA_ID2, "^.*\|", "")), NA_character_)
        ),
        circRNA_ID2 = coalesce_chr(
          circRNA_ID2,
          if_else(!is.na(circ_id) & !is.na(host_label), paste0(circ_id, " | ", host_label), NA_character_)
        ),
        circRNA_id_AS = coalesce_chr(
          circRNA_id_AS,
          if_else(!is.na(circ_id) & !is.na(length), paste0(circ_id, "_", length), NA_character_)
        ),
        circRNA_DE_id = coalesce_chr(
          circRNA_DE_id,
          if_else(!is.na(circRNA_ID2) & !is.na(direction) & !is.na(length),
                  paste0(circRNA_ID2, "_", str_to_title(str_to_lower(direction)), "_", length),
                  NA_character_)
        ),
        host_label = na_if(str_trim(host_label), ""),
        circ_key = normalize_key(circ_id),
        key_id2 = normalize_key(circRNA_ID2),
        key_as = normalize_key(circRNA_id_AS),
        key_deid = normalize_key(circRNA_DE_id)
      ) %>%
      distinct()
  }

  de_tbl <- bind_rows(read_circ_de(de6_file, "6h"), read_circ_de(de24_file, "24h"))

  phase_files <- tribble(
    ~phase, ~file,
    "Early", find_optional_file("early_circRNA.csv", base),
    "Sustained", find_optional_file("sustained_circRNA.csv", base),
    "Late", find_optional_file("late_circRNA.csv", base),
    "Reversal", find_optional_file("reversal_circRNA.csv", base)
  )

  phase_map <- pmap_dfr(phase_files, function(phase, file) {
    if (is.na(file)) return(tibble())
    z <- read_delim_auto(file)
    circ_col <- choose_col(z, c("circ_id", "circRNA_ID", "circRNA_ID2"), names(z)[1])
    circ_raw <- if (!is.null(circ_col)) as.character(z[[circ_col]]) else rep(NA_character_, nrow(z))
    tibble(
      circ_key = normalize_key(canonicalize_circ_id(circ_raw)),
      phase = phase
    )
  }) %>%
    filter(!is.na(circ_key), circ_key != "") %>%
    distinct()

  if (nrow(phase_map) == 0) {
    phase_map <- de_tbl %>%
      select(circ_key, timepoint, log2FC, pval_like) %>%
      distinct() %>%
      pivot_wider(names_from = timepoint, values_from = c(log2FC, pval_like), names_sep = "_") %>%
      mutate(phase = infer_phase_two_timepoints(log2FC_6h, pval_like_6h, log2FC_24h, pval_like_24h)) %>%
      select(circ_key, phase)
  }

  de_tbl %>%
    filter(!is.na(circ_id), circ_id != "") %>%
    left_join(phase_map, by = "circ_key") %>%
    mutate(
      phase = coalesce(phase, "Unassigned"),
      ncrna_type = "circRNA",
      ncrna_id_raw = circ_id,
      ncrna_id_key = circ_key,
      gene_id = circ_id,
      transcript_id = NA_character_,
      source = "circRNA",
      ncRNA_label = coalesce_chr(circRNA_ID2, circ_id)
    ) %>%
    select(ncrna_type, ncrna_id_raw, ncrna_id_key, gene_id, transcript_id, source, phase, timepoint, direction, ncRNA_label, host_label, circ_id, circRNA_ID2, circRNA_id_AS, circRNA_DE_id, pval_like, circ_key, key_id2, key_as, key_deid) %>%
    distinct()
}

read_miRNA_de <- function(base = ".") {
  files <- tribble(
    ~file, ~timepoint, ~direction,
    "DE_miRNA_Up_NE6.txt", "6h", "UP",
    "DE_miRNA_Down_NE6.txt", "6h", "DOWN",
    "DE_miRNA_Up_NE24.txt", "24h", "UP",
    "DE_miRNA_Down_NE24.txt", "24h", "DOWN"
  ) %>% mutate(path = map_chr(file, ~ find_optional_file(.x, base)))

  pmap_dfr(files, function(file, timepoint, direction, path) {
    if (is.na(path)) return(tibble())
    x <- read_delim_auto(path)
    mir_col <- choose_col(x, c("miRNA", "mirna", "miRNA_ID", "mirnaid", "id", "gene", "Gene"), names(x)[1])
    lfc_col <- choose_col(x, c("log2FoldChange", "log2FC"))
    padj_col <- choose_col(x, c("padj", "FDR", "adj.P.Val", "qvalue"))
    p_col <- choose_col(x, c("pvalue", "PValue"))

    tibble(
      miRNA = as.character(x[[mir_col]]),
      miRNA_key = normalize_mirna(x[[mir_col]]),
      timepoint = timepoint,
      direction = direction,
      log2FC = if (!is.null(lfc_col)) suppressWarnings(as.numeric(x[[lfc_col]])) else NA_real_,
      padj = if (!is.null(padj_col)) suppressWarnings(as.numeric(x[[padj_col]])) else NA_real_,
      pvalue = if (!is.null(p_col)) suppressWarnings(as.numeric(x[[p_col]])) else NA_real_
    )
  }) %>% distinct()
}

read_mRNA_map <- function(base = ".") {
  de6_file <- find_optional_file("NE_vs_Ctrl_6h_mRNA_DE.txt", base)
  de24_file <- find_optional_file("NE_vs_Ctrl_24h_mRNA_DE.txt", base)
  if (is.na(de6_file) && is.na(de24_file)) stop("Need at least one mRNA DE file.")

  phase_tbls <- list(
    early = find_optional_file("early_mRNA.csv", base),
    sustained = find_optional_file("sustained_mRNA.csv", base),
    late = find_optional_file("late_mRNA.csv", base)
  )

  phase_map <- bind_rows(
    if (!is.na(phase_tbls$early)) read_csv(phase_tbls$early, show_col_types = FALSE) %>% mutate(phase = "Early") else tibble(),
    if (!is.na(phase_tbls$sustained)) read_csv(phase_tbls$sustained, show_col_types = FALSE) %>% mutate(phase = "Sustained") else tibble(),
    if (!is.na(phase_tbls$late)) read_csv(phase_tbls$late, show_col_types = FALSE) %>% mutate(phase = "Late") else tibble()
  )

  if (nrow(phase_map) > 0) {
    phase_map <- phase_map %>%
      mutate(
        gene_symbol_tmp = coalesce_chr(choose_vec(., c("gene_name", "gene_symbol", "symbol", "gene_id")), choose_vec(., c("gene_id"))),
        gene_key = normalize_key(gene_symbol_tmp)
      ) %>%
      transmute(gene_key, phase) %>%
      distinct(gene_key, .keep_all = TRUE)
  } else {
    phase_map <- tibble(gene_key = character(), phase = character())
  }

  read_mrna_de <- function(path, tp) {
    if (is.na(path)) return(tibble())
    x <- read_delim_auto(path)
    sym_col <- choose_col(x, c("gene_name", "gene_symbol", "genesymbol", "symbol", "gene_id"), names(x)[1])
    lfc_col <- choose_col(x, c("log2FoldChange", "log2FC"))
    dir_col <- choose_col(x, c("direction", "regulations"))
    padj_col <- choose_col(x, c("padj", "FDR", "adj.P.Val"))

    tibble(
      gene_symbol = as.character(x[[sym_col]]),
      gene_key = normalize_key(x[[sym_col]]),
      timepoint = tp,
      log2FC = suppressWarnings(as.numeric(if (!is.null(lfc_col)) x[[lfc_col]] else NA_real_)),
      direction = canonical_direction(if (!is.null(dir_col)) x[[dir_col]] else NA_character_, if (!is.null(lfc_col)) x[[lfc_col]] else NA_real_),
      padj = suppressWarnings(as.numeric(if (!is.null(padj_col)) x[[padj_col]] else NA_real_))
    )
  }

  bind_rows(read_mrna_de(de6_file, "6h"), read_mrna_de(de24_file, "24h")) %>%
    filter(!is.na(gene_key), gene_key != "") %>%
    left_join(phase_map, by = "gene_key") %>%
    mutate(phase = coalesce(phase, "Unassigned")) %>%
    select(gene_symbol, gene_key, timepoint, log2FC, direction, padj, phase) %>%
    distinct(gene_key, timepoint, .keep_all = TRUE)
}

read_miRWalk_targets <- function(base = ".") {
  files <- list.files(base, recursive = TRUE, full.names = TRUE, pattern = "^miRWalk_miRNA_Targets_.*\\.csv$")
  if (length(files) == 0) {
    return(list(targets = tibble(), audit = tibble()))
  }

  out <- map(files, function(path) {
    x <- read_delim_auto(path)
    fn <- basename(path)

    mir_col <- choose_col(
      x,
      c("miRNA", "mirna", "mirnaid", "miRNA_name", "mature_miRNA", "mature_mirna"),
      names(x)[1]
    )

    gene_col <- choose_col(
      x,
      c("gene_symbol", "genesymbol", "gene", "target_gene", "gene_name", "Symbol", "symbol", "GeneSymbol"),
      if (ncol(x) >= 3) names(x)[3] else names(x)[min(2, ncol(x))]
    )

    refseq_col <- choose_col(x, c("refseqid", "refseq_id", "RefSeq", "transcript_id"))
    region_col <- choose_col(x, c("position", "region", "binding_region"))
    score_col <- choose_col(x, c("energy", "score", "binding_score", "probability", "bindingp", "TargetScan", "miRDB"))

    region_from_file <- case_when(
      str_detect(fn, regex("3UTR", ignore_case = TRUE)) ~ "3UTR",
      str_detect(fn, regex("CDS", ignore_case = TRUE)) ~ "CDS",
      str_detect(fn, regex("5UTR", ignore_case = TRUE)) ~ "5UTR",
      TRUE ~ "Unknown"
    )

    timepoint_from_file <- case_when(
      str_detect(fn, regex("24", ignore_case = TRUE)) ~ "24h",
      str_detect(fn, regex("6", ignore_case = TRUE)) ~ "6h",
      TRUE ~ NA_character_
    )

    mir_dir_from_file <- case_when(
      str_detect(fn, regex("UP", ignore_case = TRUE)) ~ "UP",
      str_detect(fn, regex("DOWN", ignore_case = TRUE)) ~ "DOWN",
      TRUE ~ NA_character_
    )

    targets <- tibble(
      miRNA = as.character(x[[mir_col]]),
      miRNA_key = normalize_mirna(x[[mir_col]]),
      gene_symbol_miRWalk = as.character(x[[gene_col]]),
      gene_key = normalize_key(x[[gene_col]]),
      refseqid = if (!is.null(refseq_col)) as.character(x[[refseq_col]]) else NA_character_,
      region = if (!is.null(region_col)) as.character(x[[region_col]]) else region_from_file,
      timepoint = timepoint_from_file,
      miRNA_direction_file = mir_dir_from_file,
      support_file = fn,
      target_score = suppressWarnings(as.numeric(if (!is.null(score_col)) x[[score_col]] else NA_real_))
    ) %>%
      mutate(region = coalesce(region, region_from_file)) %>%
      filter(!is.na(miRNA_key), miRNA_key != "", !is.na(gene_key), gene_key != "")

    audit <- tibble(
      file = fn,
      n_rows = nrow(x),
      n_rows_kept = nrow(targets),
      miRNA_col = mir_col,
      gene_col = gene_col,
      refseq_col = refseq_col %||% NA_character_,
      region_col = region_col %||% NA_character_,
      score_col = score_col %||% NA_character_,
      timepoint = timepoint_from_file,
      miRNA_direction_file = mir_dir_from_file,
      n_unique_miRNAs = n_distinct(targets$miRNA_key),
      n_unique_gene_keys = n_distinct(targets$gene_key),
      sample_columns = paste(names(x), collapse = " | ")
    )

    list(targets = targets, audit = audit)
  })

  list(
    targets = bind_rows(map(out, "targets")) %>% distinct(),
    audit = bind_rows(map(out, "audit"))
  )
}

annotate_intarna_pairs <- function(intarna_tbl, lnc_map, circ_map) {
  if (nrow(intarna_tbl) == 0) return(tibble())

  lnc_lookup_tx <- lnc_map %>% filter(!is.na(tx_key), tx_key != "") %>% distinct(tx_key, timepoint, .keep_all = TRUE)
  lnc_lookup_gene <- lnc_map %>% filter(!is.na(gene_key), gene_key != "") %>% distinct(gene_key, timepoint, .keep_all = TRUE)

  circ_by_deid <- circ_map %>% filter(!is.na(circRNA_DE_id), circRNA_DE_id != "") %>% distinct(timepoint, circRNA_DE_id, .keep_all = TRUE)
  circ_by_as <- circ_map %>% filter(!is.na(circRNA_id_AS), circRNA_id_AS != "") %>% distinct(timepoint, circRNA_id_AS, .keep_all = TRUE)
  circ_by_id2 <- circ_map %>% filter(!is.na(circRNA_ID2), circRNA_ID2 != "") %>% group_by(timepoint, circRNA_ID2) %>% slice_min(order_by = pval_like, n = 1, with_ties = FALSE) %>% ungroup()
  circ_by_circ <- circ_map %>% filter(!is.na(circ_id), circ_id != "") %>% group_by(timepoint, circ_id) %>% slice_min(order_by = pval_like, n = 1, with_ties = FALSE) %>% ungroup()

  lnc_pairs <- intarna_tbl %>%
    filter(ncrna_type == "lncRNA") %>%
    bind_cols(parse_lnc_intarna_id1(.$id1)) %>%
    mutate(join_timepoint = coalesce(timepoint_from_id1, timepoint_file)) %>%
    left_join(
      lnc_lookup_tx %>%
        transmute(join_timepoint = timepoint, tx_key, gene_id_tx = gene_id, transcript_id_tx = transcript_id, gene_name_tx = gene_name, source_tx = source, phase_tx = phase, direction_tx = direction, ncRNA_label_tx = ncRNA_label, ncrna_id_key_tx = ncrna_id_key),
      by = c("join_timepoint", "tx_key_from_id1" = "tx_key")
    ) %>%
    left_join(
      lnc_lookup_gene %>%
        transmute(join_timepoint = timepoint, gene_key, gene_id_gene = gene_id, transcript_id_gene = transcript_id, gene_name_gene = gene_name, source_gene = source, phase_gene = phase, direction_gene = direction, ncRNA_label_gene = ncRNA_label, ncrna_id_key_gene = ncrna_id_key),
      by = c("join_timepoint", "gene_key_from_id1" = "gene_key")
    ) %>%
    mutate(
      gene_id = coalesce(gene_id_tx, gene_id_gene, gene_id_from_id1),
      transcript_id = coalesce(transcript_id_tx, transcript_id_gene, transcript_id_from_id1),
      gene_name = coalesce(gene_name_tx, gene_name_gene),
      source = coalesce(source_tx, source_gene, source_from_id1),
      phase = coalesce(phase_tx, phase_gene, phase_from_id1, "Unassigned"),
      ncRNA_label = coalesce(ncRNA_label_tx, ncRNA_label_gene, gene_id, transcript_id),
      ncRNA_direction = canonical_direction(coalesce(direction_tx, direction_gene, direction_from_id1, expected_ncrna_direction)),
      timepoint = coalesce(join_timepoint, timepoint_from_id1, timepoint_file),
      ncrna_match_method = case_when(
        !is.na(gene_id_tx) ~ "transcript_id_plus_timepoint",
        !is.na(gene_id_gene) ~ "gene_id_plus_timepoint",
        !is.na(gene_id_from_id1) ~ "header_gene_id_only",
        TRUE ~ "raw_lnc_id"
      ),
      ncrna_id_raw = coalesce(transcript_id, gene_id),
      ncrna_id_key = normalize_key(coalesce(transcript_id, gene_id)),
      ncrna_type = "lncRNA"
    ) %>%
    select(source_file, ncrna_type, id1, miRNA, miRNA_key, E, interaction_tier, expected_ncrna_direction, expected_miRNA_direction, timepoint_file, gene_id, transcript_id, gene_name, source, phase, timepoint, ncRNA_direction, ncRNA_label, ncrna_match_method, ncrna_id_raw, ncrna_id_key) %>%
    distinct()

  circ_pairs <- intarna_tbl %>%
    filter(ncrna_type == "circRNA") %>%
    bind_cols(parse_circ_intarna_id1(.$id1)) %>%
    mutate(join_timepoint = coalesce(timepoint_from_id1, timepoint_file)) %>%
    left_join(
      circ_by_deid %>% transmute(join_timepoint = timepoint, circRNA_DE_id, circ_id_deid = circ_id, phase_deid = phase, direction_deid = direction, ncRNA_label_deid = ncRNA_label, host_label_deid = host_label),
      by = c("join_timepoint", "circRNA_DE_id_from_id1" = "circRNA_DE_id")
    ) %>%
    left_join(
      circ_by_as %>% transmute(join_timepoint = timepoint, circRNA_id_AS, circ_id_as = circ_id, phase_as = phase, direction_as = direction, ncRNA_label_as = ncRNA_label, host_label_as = host_label),
      by = c("join_timepoint", "circRNA_id_AS_from_id1" = "circRNA_id_AS")
    ) %>%
    left_join(
      circ_by_id2 %>% transmute(join_timepoint = timepoint, circRNA_ID2, circ_id_id2 = circ_id, phase_id2 = phase, direction_id2 = direction, ncRNA_label_id2 = ncRNA_label, host_label_id2 = host_label),
      by = c("join_timepoint", "circRNA_ID2_from_id1" = "circRNA_ID2")
    ) %>%
    left_join(
      circ_by_circ %>% transmute(join_timepoint = timepoint, circ_id, phase_circ = phase, direction_circ = direction, ncRNA_label_circ = ncRNA_label, host_label_circ = host_label),
      by = c("join_timepoint", "circ_id_from_id1" = "circ_id")
    ) %>%
    mutate(
      gene_id = coalesce(circ_id_deid, circ_id_as, circ_id_id2, circ_id_from_id1),
      transcript_id = NA_character_,
      source = "circRNA",
      phase = coalesce(phase_deid, phase_as, phase_id2, phase_circ, "Unassigned"),
      ncRNA_direction = canonical_direction(coalesce(direction_deid, direction_as, direction_id2, direction_circ, expected_direction_from_id1, expected_ncrna_direction)),
      ncRNA_label = coalesce(ncRNA_label_deid, ncRNA_label_as, ncRNA_label_id2, ncRNA_label_circ, circRNA_ID2_from_id1, circ_id_from_id1),
      host_label = coalesce(host_label_deid, host_label_as, host_label_id2, host_label_circ, host_label_from_id1),
      timepoint = join_timepoint,
      ncrna_match_method = case_when(
        !is.na(circ_id_deid) ~ "circRNA_DE_id_plus_timepoint",
        !is.na(circ_id_as) ~ "circRNA_id_AS_plus_timepoint",
        !is.na(circ_id_id2) ~ "circRNA_ID2_plus_timepoint",
        !is.na(circ_id_from_id1) ~ "circ_id_from_header",
        TRUE ~ "unmatched_circRNA"
      ),
      ncrna_id_raw = coalesce(circRNA_DE_id_from_id1, circRNA_id_AS_from_id1, circRNA_ID2_from_id1, circ_id_from_id1),
      ncrna_id_key = normalize_key(coalesce(gene_id, circ_id_from_id1)),
      ncrna_type = "circRNA"
    ) %>%
    select(source_file, ncrna_type, id1, miRNA, miRNA_key, E, interaction_tier, expected_ncrna_direction, expected_miRNA_direction, timepoint_file, gene_id, transcript_id, source, phase, timepoint, ncRNA_direction, ncRNA_label, ncrna_match_method, ncrna_id_raw, ncrna_id_key, host_label) %>%
    distinct()

  bind_rows(lnc_pairs, circ_pairs) %>% distinct()
}

# -----------------------------
# main
# -----------------------------
interaction_files <- if (is.null(interaction_files_manual)) find_all_interaction_files(base_dir) else file.path(base_dir, interaction_files_manual)
if (length(interaction_files) == 0) stop("No IntaRNA files found. Expected files matching Interact_*.csv")

input_inventory <- tibble(
  category = c(
    rep("IntaRNA", length(interaction_files)),
    "lncRNA_phase_map", "lncRNA_DE_6h", "lncRNA_DE_24h",
    "circRNA_DE_6h", "circRNA_DE_24h",
    "mRNA_DE_6h", "mRNA_DE_24h",
    "mRNA_phase_early", "mRNA_phase_sustained", "mRNA_phase_late",
    "miRNA_DE_up_6h", "miRNA_DE_down_6h", "miRNA_DE_up_24h", "miRNA_DE_down_24h"
  ),
  file = c(
    interaction_files,
    find_optional_file("lncRNA_DE_with_phase_and_transcriptID.csv", base_dir),
    find_optional_file("NE_vs_Ctrl_6h_lncRNA_DE.txt", base_dir),
    find_optional_file("NE_vs_Ctrl_24h_lncRNA_DE.txt", base_dir),
    find_optional_file("06_circRNA_sequences_de6_DE_complete_table.txt", base_dir),
    find_optional_file("07_circRNA_sequences_de24_DE_complete_table.txt", base_dir),
    find_optional_file("NE_vs_Ctrl_6h_mRNA_DE.txt", base_dir),
    find_optional_file("NE_vs_Ctrl_24h_mRNA_DE.txt", base_dir),
    find_optional_file("early_mRNA.csv", base_dir),
    find_optional_file("sustained_mRNA.csv", base_dir),
    find_optional_file("late_mRNA.csv", base_dir),
    find_optional_file("DE_miRNA_Up_NE6.txt", base_dir),
    find_optional_file("DE_miRNA_Down_NE6.txt", base_dir),
    find_optional_file("DE_miRNA_Up_NE24.txt", base_dir),
    find_optional_file("DE_miRNA_Down_NE24.txt", base_dir)
  )
) %>%
  mutate(found = !is.na(file), file = if_else(found, file, NA_character_))
write_csv(input_inventory, file.path(out_dir, "QC_input_inventory.csv"))

intarna_tbl <- map_dfr(interaction_files, read_intarna_file)
lnc_map <- read_lnc_annotation_map(base_dir)
circ_map <- read_circ_annotation_map(base_dir)
miRNA_de <- read_miRNA_de(base_dir)
mRNA_map <- read_mRNA_map(base_dir)
miRWalk_info <- read_miRWalk_targets(base_dir)
miRWalk_targets <- miRWalk_info$targets
miRWalk_audit <- miRWalk_info$audit

write_csv(miRWalk_audit, file.path(out_dir, "QC_miRWalk_file_audit.csv"))
write_csv(miRWalk_targets, file.path(out_dir, "DEBUG_miRWalk_targets_standardized.csv"))
write_csv(bind_rows(lnc_map %>% mutate(map_type = "lncRNA"), circ_map %>% mutate(map_type = "circRNA")), file.path(out_dir, "DEBUG_ncRNA_annotation_map.csv"))
write_csv(mRNA_map, file.path(out_dir, "DEBUG_mRNA_map.csv"))
write_csv(miRNA_de, file.path(out_dir, "DEBUG_miRNA_DE_map.csv"))

annot_pairs <- annotate_intarna_pairs(intarna_tbl, lnc_map, circ_map) %>%
  left_join(
    miRNA_de %>%
      select(miRNA_key, timepoint, miRNA_direction = direction, miRNA_log2FC = log2FC, miRNA_padj = padj) %>%
      distinct(),
    by = c("miRNA_key", "timepoint")
  ) %>%
  mutate(
    miRNA_direction = coalesce(miRNA_direction, expected_miRNA_direction),
    coherent_pair = case_when(
      ncRNA_direction == "UP" & miRNA_direction == "DOWN" ~ TRUE,
      ncRNA_direction == "DOWN" & miRNA_direction == "UP" ~ TRUE,
      TRUE ~ FALSE
    ),
    pair_time_state = timepoint,
    phase = coalesce(phase, "Unassigned")
  ) %>%
  distinct()

if (keep_only_moderate_or_better) {
  annot_pairs <- annot_pairs %>% filter(interaction_tier %in% c("High", "Moderate"))
}

pair_presence <- annot_pairs %>%
  filter(coherent_pair) %>%
  mutate(ncrna_key = normalize_key(ncRNA_label)) %>%
  distinct(ncrna_type, ncrna_key, miRNA_key, timepoint) %>%
  count(ncrna_type, ncrna_key, miRNA_key, timepoint) %>%
  pivot_wider(names_from = timepoint, values_from = n, values_fill = 0) %>%
  mutate(
    interaction_state = case_when(
      `6h` > 0 & `24h` == 0 ~ "Early",
      `6h` == 0 & `24h` > 0 ~ "Late",
      `6h` > 0 & `24h` > 0 ~ "Sustained",
      TRUE ~ "Unassigned"
    )
  )

coherent_pairs <- annot_pairs %>%
  filter(coherent_pair) %>%
  mutate(ncrna_key = normalize_key(ncRNA_label)) %>%
  left_join(pair_presence %>% select(ncrna_type, ncrna_key, miRNA_key, interaction_state), by = c("ncrna_type", "ncrna_key", "miRNA_key")) %>%
  mutate(interaction_state = coalesce(interaction_state, timepoint))

write_csv(coherent_pairs %>% filter(ncrna_type == "lncRNA"), file.path(out_dir, "COHERENT_lncRNA_miRNA_pairs.csv"))
write_csv(coherent_pairs %>% filter(ncrna_type == "circRNA"), file.path(out_dir, "COHERENT_circRNA_miRNA_pairs.csv"))

# strict miRNA-target join
miRWalk_targets2 <- miRWalk_targets %>%
  filter(!is.na(miRNA_key), miRNA_key != "", !is.na(gene_key), gene_key != "", !is.na(timepoint), timepoint != "")

mRNA_map2 <- mRNA_map %>%
  filter(!is.na(gene_key), gene_key != "", !is.na(timepoint), timepoint != "") %>%
  distinct(gene_key, timepoint, .keep_all = TRUE)

triplet_seed <- coherent_pairs %>%
  select(
    source_file, ncrna_type, gene_id, transcript_id, ncRNA_label, phase,
    timepoint, ncRNA_direction, miRNA, miRNA_key, miRNA_direction,
    E, interaction_tier, interaction_state, ncrna_match_method
  ) %>%
  filter(!is.na(miRNA_key), miRNA_key != "", !is.na(timepoint), timepoint != "")

triplets_after_target_join <- triplet_seed %>%
  left_join(
    miRWalk_targets2 %>%
      select(miRNA_key, gene_key, gene_symbol_miRWalk, refseqid, region, timepoint, support_file, target_score),
    by = c("miRNA_key", "timepoint")
  )

triplets_full <- triplets_after_target_join %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  left_join(
    mRNA_map2 %>%
      select(
        gene_symbol, gene_key, timepoint,
        phase_mRNA = phase,
        mRNA_direction = direction,
        mRNA_log2FC = log2FC,
        mRNA_padj = padj
      ),
    by = c("gene_key", "timepoint")
  ) %>%
  mutate(
    mRNA_label = coalesce(
      na_if(gene_symbol, ""),
      na_if(gene_symbol_miRWalk, ""),
      na_if(gene_key, "")
    ),
    miRNA_direction = canonical_direction(miRNA_direction),
    mRNA_direction = canonical_direction(mRNA_direction, mRNA_log2FC),
    coherent_triplet = case_when(
      ncRNA_direction == "UP" & miRNA_direction == "DOWN" & mRNA_direction == "UP" ~ TRUE,
      ncRNA_direction == "DOWN" & miRNA_direction == "UP" & mRNA_direction == "DOWN" ~ TRUE,
      TRUE ~ FALSE
    ),
    phase_mRNA = coalesce(phase_mRNA, "Unassigned"),
    phase_relation = paste0(phase, "_ncRNA__", phase_mRNA, "_mRNA")
  )

triplets <- triplets_full %>%
  filter(coherent_triplet, !is.na(mRNA_label), str_trim(mRNA_label) != "")

write_csv(triplets %>% filter(ncrna_type == "lncRNA"), file.path(out_dir, "COHERENT_lncRNA_miRNA_mRNA_triplets.csv"))
write_csv(triplets %>% filter(ncrna_type == "circRNA"), file.path(out_dir, "COHERENT_circRNA_miRNA_mRNA_triplets.csv"))

priority_triplets <- triplets %>%
  mutate(tier_rank = rank_tier(interaction_tier)) %>%
  arrange(tier_rank, E, desc(abs(mRNA_log2FC))) %>%
  group_by(ncrna_type, ncRNA_label, miRNA_key, gene_key) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(ncrna_type, tier_rank, E)

write_csv(priority_triplets, file.path(out_dir, "PRIORITY_ceRNA_triplets_best_per_ncRNA_miRNA_mRNA.csv"))

pair_summary <- coherent_pairs %>%
  count(ncrna_type, timepoint, phase, ncRNA_direction, miRNA_direction, interaction_tier, interaction_state, name = "n_pairs") %>%
  arrange(ncrna_type, timepoint, phase, ncRNA_direction, miRNA_direction, interaction_tier, interaction_state)

triplet_summary <- triplets %>%
  count(ncrna_type, timepoint, phase, phase_mRNA, ncRNA_direction, miRNA_direction, mRNA_direction, interaction_tier, interaction_state, name = "n_triplets") %>%
  arrange(ncrna_type, timepoint, phase, phase_mRNA, ncRNA_direction, miRNA_direction, mRNA_direction, interaction_tier, interaction_state)

write_csv(pair_summary, file.path(out_dir, "SUMMARY_coherent_ncRNA_miRNA_pairs.csv"))
write_csv(triplet_summary, file.path(out_dir, "SUMMARY_coherent_triplets_by_phase.csv"))

# QC
missing_mirna_targets <- triplet_seed %>%
  distinct(miRNA_key, timepoint) %>%
  anti_join(miRWalk_targets2 %>% distinct(miRNA_key, timepoint), by = c("miRNA_key", "timepoint")) %>%
  arrange(timepoint, miRNA_key)

missing_gene_in_mrna_map <- triplets_after_target_join %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, timepoint, gene_symbol_miRWalk, refseqid) %>%
  anti_join(mRNA_map2 %>% distinct(gene_key, timepoint), by = c("gene_key", "timepoint")) %>%
  arrange(timepoint, gene_key)

join_qc <- tibble(
  metric = c(
    "n_intarna_rows",
    "n_annot_pairs",
    "n_coherent_pairs",
    "n_distinct_pair_miRNA_time",
    "n_miRWalk_targets",
    "n_triplet_seed_rows",
    "n_rows_after_miRWalk_join",
    "n_rows_with_gene_key_after_miRWalk_join",
    "n_rows_after_mRNA_join",
    "n_rows_with_mRNA_direction",
    "n_coherent_triplets",
    "n_priority_triplets",
    "n_missing_miRNA_time_in_miRWalk",
    "n_missing_gene_time_in_mRNA_map"
  ),
  value = c(
    nrow(intarna_tbl),
    nrow(annot_pairs),
    nrow(coherent_pairs),
    coherent_pairs %>% distinct(miRNA_key, timepoint) %>% nrow(),
    nrow(miRWalk_targets2),
    nrow(triplet_seed),
    nrow(triplets_after_target_join),
    sum(!is.na(triplets_after_target_join$gene_key) & triplets_after_target_join$gene_key != ""),
    nrow(triplets_full),
    sum(!is.na(triplets_full$mRNA_direction) & triplets_full$mRNA_direction != ""),
    nrow(triplets),
    nrow(priority_triplets),
    nrow(missing_mirna_targets),
    nrow(missing_gene_in_mrna_map)
  )
)

write_csv(join_qc, file.path(out_dir, "QC_join_summary.csv"))
write_csv(missing_mirna_targets, file.path(out_dir, "QC_missing_miRNA_targets_by_time.csv"))
write_csv(missing_gene_in_mrna_map, file.path(out_dir, "QC_missing_mRNA_map_by_gene_time.csv"))

# simple plots
p1 <- pair_summary %>%
  ggplot(aes(x = interaction_tier, y = n_pairs, fill = interaction_state)) +
  geom_col(position = "stack") +
  facet_grid(ncrna_type ~ timepoint, scales = "free_y") +
  labs(title = "Figure 8 - coherent ncRNA-miRNA pairs", x = "Interaction tier", y = "Number of coherent pairs") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "FIG_ceRNA_pairs_by_tier_and_state.png"), p1, width = 11, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "FIG_ceRNA_pairs_by_tier_and_state.pdf"), p1, width = 11, height = 6.5)

p2 <- triplet_summary %>%
  mutate(phase_relation2 = paste0(phase, " -> ", phase_mRNA)) %>%
  ggplot(aes(x = interaction_tier, y = n_triplets, fill = interaction_state)) +
  geom_col(position = "stack") +
  facet_grid(ncrna_type ~ phase_relation2, scales = "free_y") +
  labs(title = "Figure 8 - coherent ceRNA triplets", x = "Interaction tier", y = "Number of coherent triplets") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "FIG_ceRNA_triplets_by_phase_relation_and_tier.png"), p2, width = 13, height = 7, dpi = 300)
ggsave(file.path(out_dir, "FIG_ceRNA_triplets_by_phase_relation_and_tier.pdf"), p2, width = 13, height = 7)

# readable text exports for report-building
write_df_txt(pair_summary, file.path(out_dir, "Fig8_pair_summary.txt"), "Figure 8 - Pair summary")
write_df_txt(priority_triplets, file.path(out_dir, "Fig8_priority_ceRNA_triplets.txt"), "Figure 8 - Priority ceRNA triplets")
write_df_txt(triplet_summary, file.path(out_dir, "Fig8_triplet_summary.txt"), "Figure 8 - Triplet summary")

fig8_lines <- c("Figure 8. Integrative ceRNA candidate network", strrep("=", 48), "")

if (nrow(pair_summary) > 0) {
  txt <- pair_summary %>%
    mutate(tag = paste0(ncrna_type, " ", timepoint, " ", interaction_tier, "=", format(n_pairs, justify = "right", width = 4))) %>%
    pull(tag)
  fig8_lines <- append_lines(fig8_lines, paste0("Coherent ncRNA-miRNA pairs by layer, time, and evidence tier: ", paste(txt, collapse = "; "), "."))
}

if (nrow(priority_triplets) > 0) {
  fig8_lines <- append_lines(
    fig8_lines,
    paste0(
      "Representative priority ceRNA triplets: ",
      priority_triplets %>%
        mutate(
          label = paste0(
            ncRNA_label, " -- ", miRNA, " -- ", mRNA_label,
            " [", timepoint, ", ", phase, "->", phase_mRNA, ", E=", fmt_num(E, 3), "]"
          )
        ) %>%
        slice_head(n = 40) %>%
        pull(label) %>%
        paste(collapse = "; "),
      "."
    )
  )

  rep_mi <- priority_triplets %>% count(miRNA, sort = TRUE, name = "n_priority_triplets")
  rep_mr <- priority_triplets %>% count(mRNA_label, sort = TRUE, name = "n_priority_triplets")
  rep_nc <- priority_triplets %>% count(ncRNA_label, sort = TRUE, name = "n_priority_triplets")

  fig8_lines <- append_lines(
    fig8_lines,
    paste0("Most recurrent miRNAs across priority triplets: ", format_top_features(rep_mi, 20, "miRNA", "n_priority_triplets", NULL), "."),
    paste0("Most recurrent mRNA targets across priority triplets: ", format_top_features(rep_mr, 20, "mRNA_label", "n_priority_triplets", NULL), "."),
    paste0("Most recurrent ncRNA candidates across priority triplets: ", format_top_features(rep_nc, 20, "ncRNA_label", "n_priority_triplets", NULL), ".")
  )
}

if (nrow(missing_mirna_targets) > 0) {
  fig8_lines <- append_lines(
    fig8_lines,
    paste0(
      "QC: Some coherent miRNA-timepoint combinations were not found in the miRWalk target tables and were excluded from ceRNA triplets: ",
      collapse_labels(paste0(missing_mirna_targets$timepoint, ":", missing_mirna_targets$miRNA_key), 30),
      "."
    )
  )
}

if (nrow(missing_gene_in_mrna_map) > 0) {
  fig8_lines <- append_lines(
    fig8_lines,
    paste0(
      "QC: Some miRWalk target genes were not found in the mRNA DE map for the same timepoint and were excluded from coherent triplets: ",
      collapse_labels(paste0(missing_gene_in_mrna_map$timepoint, ":", missing_gene_in_mrna_map$gene_symbol_miRWalk), 30),
      "."
    )
  )
}

safe_write_lines(fig8_lines, file.path(out_dir, "Fig8_results_report.txt"))

sink(file.path(out_dir, "README_ceRNA_run_summary.txt"))
cat("ceRNA integration run summary\n===========================\n\n")
cat("Input IntaRNA files detected:", length(interaction_files), "\n")
cat(paste0(" - ", basename(interaction_files), collapse = "\n"), "\n\n")
cat("Interaction tiers:\n")
cat(" - High        : E <= ", energy_high, "\n", sep = "")
cat(" - Moderate    : ", energy_high, " < E <= ", energy_moderate, "\n", sep = "")
cat(" - Exploratory : ", energy_moderate, " < E <= ", energy_exploratory, "\n", sep = "")
cat(" - Weak        : E > ", energy_exploratory, "\n\n", sep = "")
cat("Core outputs:\n")
cat(" - COHERENT_lncRNA_miRNA_pairs.csv\n")
cat(" - COHERENT_circRNA_miRNA_pairs.csv\n")
cat(" - COHERENT_lncRNA_miRNA_mRNA_triplets.csv\n")
cat(" - COHERENT_circRNA_miRNA_mRNA_triplets.csv\n")
cat(" - PRIORITY_ceRNA_triplets_best_per_ncRNA_miRNA_mRNA.csv\n")
cat(" - SUMMARY_coherent_ncRNA_miRNA_pairs.csv\n")
cat(" - SUMMARY_coherent_triplets_by_phase.csv\n")
cat(" - Fig8_results_report.txt\n\n")
cat("QC outputs:\n")
cat(" - QC_input_inventory.csv\n")
cat(" - QC_miRWalk_file_audit.csv\n")
cat(" - QC_join_summary.csv\n")
cat(" - QC_missing_miRNA_targets_by_time.csv\n")
cat(" - QC_missing_mRNA_map_by_gene_time.csv\n")
cat(" - DEBUG_miRWalk_targets_standardized.csv\n")
cat(" - DEBUG_mRNA_map.csv\n")
cat(" - DEBUG_miRNA_DE_map.csv\n")
cat(" - DEBUG_ncRNA_annotation_map.csv\n\n")
cat("Notes:\n")
cat(" - miRWalk columns are auto-detected with explicit support for mirnaid and genesymbol.\n")
cat(" - The miRNA-target join is now strict and the mRNA join only runs after a real gene_key was recovered.\n")
cat(" - Priority triplets now use mRNA_label derived from gene_symbol, miRWalk gene symbol, or gene_key.\n")
sink()

message("Done. Outputs written to: ", normalizePath(out_dir))
message("Recommended first QC files:")
message(" - ", file.path(out_dir, "QC_miRWalk_file_audit.csv"))
message(" - ", file.path(out_dir, "QC_join_summary.csv"))
message(" - ", file.path(out_dir, "QC_missing_miRNA_targets_by_time.csv"))
message(" - ", file.path(out_dir, "Fig8_results_report.txt"))
