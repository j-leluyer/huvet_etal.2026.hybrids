#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  cran_pkgs <- c("data.table", "ggplot2", "dplyr", "ComplexUpset")
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }

  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(ComplexUpset)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript script/ase/make_upset_mbased_deseq2.R \\",
  "    --mbased-file=output/ase/MBASED_all_samples.tsv \\",
  "    --deseq-file=output/Table_S1_DESeq2_AG_vs_GA.tsv \\",
  "    --out-prefix=output/ase/mbased_deseq2_upset \\",
  "    [--mbased-input=table|list] [--deseq-input=table|list] \\",
  "    [--mbased-gene-col=gene_id] [--mbased-padj-col=qASE] \\",
  "    [--deseq-gene-col=gene] [--deseq-padj-col=padj] \\",
  "    [--mbased-padj-threshold=0.05] [--deseq-padj-threshold=0.05]",
  sep = "\n"
)

if ("--help" %in% args || length(args) == 0) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

assert_file <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
}

detect_gene_col <- function(dt, preferred = NULL) {
  candidates <- unique(c(preferred, "gene_id", "gene", "Gene", "name", "Name", "symbol", "SYMBOL"))
  for (candidate in candidates) {
    if (!is.null(candidate) && candidate %in% names(dt)) return(candidate)
  }
  names(dt)[1]
}

read_gene_set <- function(file, input_type, gene_col = NULL, padj_col = NULL, padj_threshold = 0.05, set_name) {
  dt <- fread(file)
  if (nrow(dt) == 0) return(character())

  if (input_type == "list") {
    gene_col <- detect_gene_col(dt, gene_col)
    genes <- unique(as.character(dt[[gene_col]]))
  } else if (input_type == "table") {
    gene_col <- detect_gene_col(dt, gene_col)
    if (is.null(padj_col) || !(padj_col %in% names(dt))) {
      stop("For ", set_name, ", padj column not found: ", padj_col)
    }
    keep <- !is.na(dt[[padj_col]]) & dt[[padj_col]] < padj_threshold
    genes <- unique(as.character(dt[[gene_col]][keep]))
  } else {
    stop("Unknown input type for ", set_name, ": ", input_type)
  }

  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  sort(unique(genes))
}

mbased_file <- get_arg("--mbased-file", NULL)
deseq_file <- get_arg("--deseq-file", NULL)
out_prefix <- get_arg("--out-prefix", "output/ase/mbased_deseq2_upset")

mbased_input <- tolower(get_arg("--mbased-input", "table"))
deseq_input <- tolower(get_arg("--deseq-input", "table"))

mbased_gene_col <- get_arg("--mbased-gene-col", "gene_id")
mbased_padj_col <- get_arg("--mbased-padj-col", "qASE")
mbased_padj_threshold <- as.numeric(get_arg("--mbased-padj-threshold", "0.05"))

deseq_gene_col <- get_arg("--deseq-gene-col", "gene")
deseq_padj_col <- get_arg("--deseq-padj-col", "padj")
deseq_padj_threshold <- as.numeric(get_arg("--deseq-padj-threshold", "0.05"))

mbased_label <- get_arg("--mbased-label", "MBASED_sig")
deseq_label <- get_arg("--deseq-label", "DESeq2_sig")
plot_title <- get_arg("--title", "Overlap of significant MBASED-biased genes and DESeq2 DEGs")

assert_file(mbased_file, "MBASED input file")
assert_file(deseq_file, "DESeq2 input file")

if (!(mbased_input %in% c("table", "list"))) stop("--mbased-input must be table or list")
if (!(deseq_input %in% c("table", "list"))) stop("--deseq-input must be table or list")
if (is.na(mbased_padj_threshold) || is.na(deseq_padj_threshold)) stop("Padj thresholds must be numeric")

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

mbased_genes <- read_gene_set(
  file = mbased_file,
  input_type = mbased_input,
  gene_col = mbased_gene_col,
  padj_col = mbased_padj_col,
  padj_threshold = mbased_padj_threshold,
  set_name = "MBASED"
)

deseq_genes <- read_gene_set(
  file = deseq_file,
  input_type = deseq_input,
  gene_col = deseq_gene_col,
  padj_col = deseq_padj_col,
  padj_threshold = deseq_padj_threshold,
  set_name = "DESeq2"
)

if (length(mbased_genes) == 0 || length(deseq_genes) == 0) {
  stop("At least one set is empty after filtering. MBASED=", length(mbased_genes), ", DESeq2=", length(deseq_genes))
}

all_genes <- sort(unique(c(mbased_genes, deseq_genes)))
membership <- data.table(
  gene = all_genes,
  MBASED = all_genes %in% mbased_genes,
  DESeq2 = all_genes %in% deseq_genes
)
setnames(membership, c("MBASED", "DESeq2"), c(mbased_label, deseq_label))

membership[[mbased_label]] <- as.logical(membership[[mbased_label]])
membership[[deseq_label]] <- as.logical(membership[[deseq_label]])

intersection_label <- paste(mbased_label, deseq_label, sep = "_AND_")
summary_dt <- data.table(
  set = c(mbased_label, deseq_label, intersection_label),
  n_genes = c(
    sum(membership[[mbased_label]]),
    sum(membership[[deseq_label]]),
    sum(membership[[mbased_label]] & membership[[deseq_label]])
  )
)

intersections_dt <- rbindlist(list(
  data.table(set = mbased_label, gene = membership[get(mbased_label), gene]),
  data.table(set = deseq_label, gene = membership[get(deseq_label), gene]),
  data.table(set = intersection_label, gene = membership[get(mbased_label) & get(deseq_label), gene])
), use.names = TRUE)

plot_df <- as.data.frame(membership)
intersect_cols <- c(mbased_label, deseq_label)
plot_df <- plot_df[rowSums(plot_df[, intersect_cols, drop = FALSE]) > 0, , drop = FALSE]

p_upset <- ComplexUpset::upset(
  plot_df,
  intersect = intersect_cols,
  min_size = 1,
  width_ratio = 0.2,
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(text = list(size = 3.2))
  ),
  set_sizes = ComplexUpset::upset_set_size()
) +
  ggtitle(plot_title)

fwrite(membership, paste0(out_prefix, ".membership.tsv"), sep = "\t")
fwrite(summary_dt, paste0(out_prefix, ".summary.tsv"), sep = "\t")
fwrite(intersections_dt, paste0(out_prefix, ".sets.tsv"), sep = "\t")

ggsave(paste0(out_prefix, ".png"), p_upset, width = 10, height = 7, dpi = 300)
ggsave(paste0(out_prefix, ".pdf"), p_upset, width = 10, height = 7)

message("MBASED significant genes: ", length(mbased_genes))
message("DESeq2 significant genes: ", length(deseq_genes))
message("Shared genes: ", summary_dt[set == intersection_label, n_genes])
message("Wrote: ", paste0(out_prefix, ".png"))
message("Wrote: ", paste0(out_prefix, ".pdf"))
message("Wrote: ", paste0(out_prefix, ".membership.tsv"))
message("Wrote: ", paste0(out_prefix, ".summary.tsv"))
message("Wrote: ", paste0(out_prefix, ".sets.tsv"))
