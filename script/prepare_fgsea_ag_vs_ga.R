#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(ggplot2)
})

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

# 1) Load AG vs GA differential expression table
res_file <- "output/Table_S1_DESeq2_AG_vs_GA.tsv"
if (!file.exists(res_file)) {
  stop("Missing output/Table_S1_DESeq2_AG_vs_GA.tsv. Run DE step first.")
}
res_df <- fread(res_file)
if (!all(c("gene", "log2FoldChange") %in% names(res_df))) {
  stop("Table_S1_DESeq2_AG_vs_GA.tsv must contain columns: gene, log2FoldChange")
}

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

res_mapped <- res_df %>%
  select(ORTHO_ID = gene, log2FoldChange) %>%
  left_join(ortho %>% select(ORTHO_ID, CG_gene), by = "ORTHO_ID") %>%
  filter(!is.na(CG_gene), !is.na(log2FoldChange)) %>%
  mutate(gene_clean = clean_gene_id(CG_gene))

# Deduplicate gene ranks if multiple ORTHO_ID map to same CG gene:
# keep value with largest absolute effect size.
ranks_tbl <- res_mapped %>%
  group_by(gene_clean) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(log2FoldChange))

ranks <- ranks_tbl$log2FoldChange
names(ranks) <- ranks_tbl$gene_clean
ranks <- ranks[is.finite(ranks)]
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

fwrite(
  data.table(gene = names(ranks), rank_log2FC = as.numeric(ranks)),
  file = "output/fgsea_ranks_AG_vs_GA_log2FC.tsv",
  sep = "\t"
)

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

# 4) Restrict pathways to rank universe and run fgsea
pathways <- split(term2gene$gene, term2gene$GO)
pathways <- lapply(pathways, intersect, names(ranks))
pathways <- pathways[sapply(pathways, length) >= 5 & sapply(pathways, length) <= 500]

if (length(pathways) == 0) {
  fwrite(
    data.table(),
    file = "output/fgsea_AG_vs_GA_log2FC.tsv",
    sep = "\t"
  )
  message("No pathways left after intersection and size filtering (5-500).")
  message("Created rank and term2gene files; FGSEA results file is empty.")
  quit(save = "no", status = 0)
}

fg <- fgseaMultilevel(
  pathways = pathways,
  stats = ranks,
  minSize = 5,
  maxSize = 500,
  eps = 0
) %>%
  arrange(padj)

fwrite(as.data.table(fg), file = "output/fgsea_AG_vs_GA_log2FC.tsv", sep = "\t")

# 5) Top enrichment plot
if (nrow(fg) > 0) {
  top_term <- fg$pathway[1]
  p <- plotEnrichment(pathways[[top_term]], ranks) +
    labs(title = paste("Top FGSEA term:", top_term)) +
    theme_classic()
  ggsave("output/fgsea_top_enrichment_AG_vs_GA.png", p, width = 7, height = 5, dpi = 300)
}

message("Created:")
message(" - output/fgsea_ranks_AG_vs_GA_log2FC.tsv")
message(" - output/fgsea_term2gene_from_go.tsv")
message(" - output/fgsea_AG_vs_GA_log2FC.tsv")
message(" - output/fgsea_top_enrichment_AG_vs_GA.png (if at least one pathway)")
