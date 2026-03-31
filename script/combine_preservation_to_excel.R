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

run_if_missing <- function(output_file, script_file) {
  if (!file.exists(output_file)) {
    cmd <- c(script_file, project_dir)
    status <- system2("Rscript", cmd)
    if (!identical(status, 0L)) {
      stop("Failed while running ", script_file)
    }
  }
}

run_if_missing(
  "output/module_preservation_summary_GG_to_AA.tsv",
  "script/explore_module_preservation_pure.R"
)
run_if_missing(
  "output/module_preservation_summary_AG_to_GA.tsv",
  "script/explore_module_preservation_hybrids.R"
)

read_tsv_checked <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

pure_summary <- read_tsv_checked("output/module_preservation_summary_GG_to_AA.tsv")
hybrid_summary <- read_tsv_checked("output/module_preservation_summary_AG_to_GA.tsv")
pure_me <- read_tsv_checked("output/module_eigengene_diff_GG_to_AA.tsv")
hybrid_me <- read_tsv_checked("output/module_eigengene_diff_AG_to_GA.tsv")

if (!is.null(pure_summary)) {
  pure_summary$comparison_label <- "AA_vs_GG"
  pure_summary$reference_group <- "GG"
  pure_summary$test_group <- "AA"
}
if (!is.null(hybrid_summary)) {
  hybrid_summary$comparison_label <- "AG_vs_GA"
  hybrid_summary$reference_group <- "AG"
  hybrid_summary$test_group <- "GA"
}
if (!is.null(pure_me)) {
  pure_me$comparison_label <- "AA_vs_GG"
  pure_me$reference_group <- "GG"
  pure_me$test_group <- "AA"
}
if (!is.null(hybrid_me)) {
  hybrid_me$comparison_label <- "AG_vs_GA"
  hybrid_me$reference_group <- "AG"
  hybrid_me$test_group <- "GA"
}

combined_summary <- do.call(
  rbind,
  Filter(Negate(is.null), list(pure_summary, hybrid_summary))
)
combined_me <- do.call(
  rbind,
  Filter(Negate(is.null), list(pure_me, hybrid_me))
)

sheets <- list()
if (!is.null(pure_summary)) sheets[["preservation_summary_aa_vs_gg"]] <- pure_summary
if (!is.null(hybrid_summary)) sheets[["preservation_summary_ag_vs_ga"]] <- hybrid_summary
if (!is.null(combined_summary)) sheets[["preservation_summary_all"]] <- combined_summary
if (!is.null(pure_me)) sheets[["eigengene_diff_aa_vs_gg"]] <- pure_me
if (!is.null(hybrid_me)) sheets[["eigengene_diff_ag_vs_ga"]] <- hybrid_me
if (!is.null(combined_me) && nrow(combined_me) > 0) sheets[["eigengene_diff_all"]] <- combined_me

if (length(sheets) == 0) {
  stop("No preservation summary files were available to write into Table S4.")
}

out_xlsx <- file.path("output", "Table_S5.preservation_summary.xlsx")
writexl::write_xlsx(sheets, path = out_xlsx)

message("Created: ", out_xlsx)
message("Included sheets: ", paste(names(sheets), collapse = ", "))
