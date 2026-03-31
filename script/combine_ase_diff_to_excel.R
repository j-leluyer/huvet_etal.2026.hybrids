#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[grep("--file=", script_args)])
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path[1], mustWork = TRUE))
} else {
  getwd()
}

project_dir <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}

setwd(project_dir)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl", repos = "https://cloud.r-project.org")
}

pick_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

diff_candidates <- c(
  "output/ase/ga_ag_diffase.gene_diffASE_GA_vs_AG.tsv",
  "data/aser/mbased.ga_vs_ag.gene_diffASE_GA_vs_AG.tsv"
)

counts_candidates <- c(
  "output/ase/ga_ag_diffase.gene_sample_counts.tsv",
  "data/aser/mbased.ga_vs_ag.gene_sample_counts.tsv"
)

diff_file <- pick_existing(diff_candidates)
counts_file <- pick_existing(counts_candidates)

if (is.na(diff_file)) {
  stop(
    "No GA-vs-AG differential ASE TSV found. Checked: ",
    paste(diff_candidates, collapse = ", ")
  )
}

sheets <- list(
  ase_diffAG_GA = readr::read_tsv(diff_file, show_col_types = FALSE, progress = FALSE)
)

if (!is.na(counts_file)) {
  sheets[["ase_gene_sample_counts"]] <- readr::read_tsv(counts_file, show_col_types = FALSE, progress = FALSE)
}

out_xlsx <- file.path("output", "Table_S6.ASE_diffAG_GA.xlsx")
writexl::write_xlsx(sheets, path = out_xlsx)

message("Created: ", out_xlsx)
message("Source diff TSV: ", diff_file)
if (!is.na(counts_file)) message("Source counts TSV: ", counts_file)
message("Included sheets: ", paste(names(sheets), collapse = ", "))
