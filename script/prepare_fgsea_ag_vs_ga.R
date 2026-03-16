#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(ggplot2)
})

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl", repos = "https://cloud.r-project.org")
}

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[grep("--file=", script_args)])
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path[1], mustWork = TRUE))
} else {
  getwd()
}

project_dir <- if (length(args) >= 1 && !startsWith(args[1], "--")) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}

setwd(project_dir)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Optional CLI: --go-file /path/to/annotation.csv
get_arg_value <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  if (idx[length(idx)] == length(args)) return(NULL)
  args[idx[length(idx)] + 1]
}

go_file <- get_arg_value("--go-file")
if (is.null(go_file)) {
  candidates <- c(
    list.files("data", pattern = "GCF_.*background\\.(tab|tsv|txt)$", full.names = TRUE),
    list.files("data", pattern = "background\\.(tab|tsv|txt)$", full.names = TRUE),
    list.files("data", pattern = "funcannot.*\\.(csv|tsv|txt)$", full.names = TRUE),
    list.files("data", pattern = "GO.*\\.(csv|tsv|txt)$", full.names = TRUE)
  )
  candidates <- unique(candidates)
  if (length(candidates) > 0) go_file <- candidates[1]
}
if (is.null(go_file)) {
  stop(
    "No GO annotation file found in data/. Provide one with --go-file <path>. "
  )
}
go_file <- normalizePath(go_file, mustWork = TRUE)

# 2) Load ortholog map (ORTHO_ID -> CG_gene)
ortho_file <- "data/ortholog_1to1.with_mito.tsv"
if (!file.exists(ortho_file)) {
  stop("Missing data/ortholog_1to1.with_mito.tsv")
}
ortho <- fread(ortho_file)
if (!all(c("ORTHO_ID", "CG_gene") %in% names(ortho))) {
  stop("ortholog table must contain ORTHO_ID and CG_gene columns")
}

# Normalize CG IDs (CGmt|... means mitochondrial)
clean_gene_id <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y <- gsub("^gene-", "", y, ignore.case = TRUE)
  y <- gsub("^CGmt\\|", "", y, ignore.case = TRUE)
  y <- toupper(y)
  y
}

get_contrast_table <- function(label, numerator, denominator) {
  out_file <- file.path("output", paste0("Table_S1_DESeq2_", label, ".tsv"))

  if (file.exists(out_file)) {
    res_df <- fread(out_file)
  } else {
    if (!file.exists("output/step1.hybrids.ortho.Rda")) {
      stop("Missing output/step1.hybrids.ortho.Rda needed to compute ", label)
    }
    load("output/step1.hybrids.ortho.Rda")
    if (!exists("dds")) {
      stop("`dds` not found in output/step1.hybrids.ortho.Rda")
    }

    samples_keep <- colnames(dds)[colData(dds)$Cross %in% c(numerator, denominator)]
    dds_sub <- dds[, samples_keep]
    dds_sub$Cross <- droplevels(dds_sub$Cross)
    design(dds_sub) <- ~ Cross
    dds_sub <- DESeq(dds_sub)
    res_obj <- results(dds_sub, contrast = c("Cross", numerator, denominator))
    res_df <- as.data.frame(res_obj)
    res_df$gene <- rownames(res_df)

    write.table(
      res_df,
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (!all(c("gene", "log2FoldChange") %in% names(res_df))) {
    stop(out_file, " must contain columns: gene, log2FoldChange")
  }

  res_df
}

build_ranks <- function(res_df) {
  res_mapped <- res_df %>%
    select(ORTHO_ID = gene, log2FoldChange) %>%
    left_join(ortho %>% select(ORTHO_ID, CG_gene), by = "ORTHO_ID") %>%
    filter(!is.na(CG_gene), !is.na(log2FoldChange)) %>%
    mutate(gene_clean = clean_gene_id(CG_gene))

  ranks_tbl <- res_mapped %>%
    group_by(gene_clean) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(log2FoldChange))

  ranks <- ranks_tbl$log2FoldChange
  names(ranks) <- ranks_tbl$gene_clean
  ranks <- ranks[is.finite(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  sort(ranks, decreasing = TRUE)
}

# 3) Load GO annotations and build term2gene
read_annotation <- function(path) {
  ext <- tolower(tools::file_ext(path))
  tab <- if (ext == "csv") {
    fread(path)
  } else {
    fread(path, sep = "\t")
  }

  # Case 1: simple 2-column background file (with or without header)
  # Normalize to columns: gene, GO
  if (ncol(tab) == 2) {
    setnames(tab, c("gene", "GO"))
    if (nrow(tab) > 0 &&
        tolower(as.character(tab$gene[1])) == "gene" &&
        tolower(as.character(tab$GO[1])) == "go") {
      tab <- tab[-1, ]
    }
    return(tab)
  }

  # Some files are exported without true headers, with the first row containing
  # names such as: gene, proteinID, Name, GO, ...
  if (all(grepl("^V\\d+$", names(tab))) && nrow(tab) > 0) {
    first_row <- as.character(unlist(tab[1, ]))
    if (any(tolower(first_row) %in% c("gene", "go"))) {
      names(tab) <- make.names(first_row, unique = TRUE)
      tab <- tab[-1, ]
    }
  }

  tab
}

go <- read_annotation(go_file)

pick_col <- function(nms, candidates) {
  m <- match(tolower(candidates), tolower(nms))
  m <- m[!is.na(m)]
  if (length(m) == 0) return(NULL)
  nms[m[1]]
}

gene_col <- pick_col(names(go), c("gene", "name", "gene_id", "id", "locus", "cg_gene"))
go_col <- pick_col(names(go), c("GO", "go", "go_terms", "go_ids", "gos"))

if (is.null(gene_col) || is.null(go_col)) {
  stop(
    "Could not detect gene/GO columns in GO file. Expected columns like gene and GO."
  )
}

term2gene <- go %>%
  transmute(gene = .data[[gene_col]], GO = .data[[go_col]]) %>%
  mutate(
    gene = clean_gene_id(gene),
    GO = as.character(GO)
  ) %>%
  separate_rows(GO, sep = ";|,|\\|") %>%
  mutate(GO = trimws(GO)) %>%
  filter(gene != "", GO != "", grepl("^GO:", GO)) %>%
  distinct(GO, gene)

fwrite(as.data.table(term2gene), file = "output/fgsea_term2gene_from_go.tsv", sep = "\t")

run_fgsea_for_comparison <- function(label, numerator, denominator) {
  res_df <- get_contrast_table(label, numerator, denominator)
  ranks <- build_ranks(res_df)

  rank_file <- file.path("output", paste0("fgsea_ranks_", label, "_log2FC.tsv"))
  fgsea_file <- file.path("output", paste0("fgsea_", label, "_log2FC.tsv"))
  plot_file <- file.path("output", paste0("fgsea_top_enrichment_", label, ".png"))

  fwrite(
    data.table(gene = names(ranks), rank_log2FC = as.numeric(ranks)),
    file = rank_file,
    sep = "\t"
  )

  pathways <- split(term2gene$gene, term2gene$GO)
  pathways <- lapply(pathways, intersect, names(ranks))
  pathways <- pathways[sapply(pathways, length) >= 5 & sapply(pathways, length) <= 500]

  if (length(pathways) == 0) {
    fwrite(data.table(), file = fgsea_file, sep = "\t")
    return(list(
      label = label,
      fgsea = data.frame(),
      ranks = data.frame(gene = names(ranks), rank_log2FC = as.numeric(ranks))
    ))
  }

  fg <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    eps = 0
  ) %>%
    arrange(padj)

  fg_xlsx <- as.data.frame(fg)
  if ("leadingEdge" %in% colnames(fg_xlsx)) {
    fg_xlsx$leadingEdge <- vapply(
      fg_xlsx$leadingEdge,
      function(x) paste(as.character(x), collapse = ";"),
      character(1)
    )
  }

  fwrite(as.data.table(fg), file = fgsea_file, sep = "\t")

  if (nrow(fg) > 0) {
    top_term <- fg$pathway[1]
    p <- plotEnrichment(pathways[[top_term]], ranks) +
      labs(title = paste("Top FGSEA term:", top_term, "[", label, "]")) +
      theme_classic()
    ggsave(plot_file, p, width = 7, height = 5, dpi = 300)
  }

  list(
    label = label,
    fgsea = fg_xlsx,
    ranks = data.frame(gene = names(ranks), rank_log2FC = as.numeric(ranks))
  )
}

comparison_specs <- list(
  list(label = "AG_vs_GA", numerator = "AG", denominator = "GA"),
  list(label = "AA_vs_GG", numerator = "AA", denominator = "GG")
)

results <- lapply(comparison_specs, function(spec) {
  run_fgsea_for_comparison(spec$label, spec$numerator, spec$denominator)
})

workbook_sheets <- list()
for (res in results) {
  suffix <- tolower(res$label)
  workbook_sheets[[paste0("fgsea_results_", suffix)]] <- res$fgsea
  workbook_sheets[[paste0("ranks_", suffix)]] <- res$ranks
}

writexl::write_xlsx(workbook_sheets, path = "output/Table_S3.gsea.xlsx")

message("Created:")
message(" - output/fgsea_term2gene_from_go.tsv")
for (res in results) {
  message(" - output/fgsea_ranks_", res$label, "_log2FC.tsv")
  message(" - output/fgsea_", res$label, "_log2FC.tsv")
  message(" - output/fgsea_top_enrichment_", res$label, ".png (if at least one pathway)")
}
message(" - output/Table_S3.gsea.xlsx")
