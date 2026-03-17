suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

ortho_file <- "data/ortholog_1to1.with_mito.tsv"
if (!file.exists(ortho_file)) {
  stop("Missing ortholog mapping file: ", ortho_file)
}

ortholog_map <- readr::read_tsv(ortho_file, show_col_types = FALSE, progress = FALSE) %>%
  dplyr::select(ORTHO_ID, `gene-AA` = CA_gene, `gene-GG` = CG_gene)

append_gene_names <- function(df) {
  key_col <- NULL
  if ("gene" %in% names(df)) key_col <- "gene"
  if (is.null(key_col) && "ORTHO_ID" %in% names(df)) key_col <- "ORTHO_ID"
  if (is.null(key_col)) return(df)

  df %>%
    dplyr::select(-dplyr::any_of(c("gene-AA", "gene-GG"))) %>%
    dplyr::left_join(ortholog_map, by = setNames("ORTHO_ID", key_col))
}

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl", repos = "https://cloud.r-project.org")
}

# Collect all ePST TSV outputs
all_tsv <- list.files("output", pattern = "^epst_.*\\.tsv$", full.names = TRUE)

if (length(all_tsv) == 0) {
  stop("No ePST TSV files found in output/ (expected files named epst_*.tsv).")
}

# Prioritize common expected tables first
preferred <- c(
  "epst_summary_stats.tsv",
  "epst_gene_values.tsv",
  "epst_with_deseq_dominance.tsv",
  "epst_set_strong.tsv",
  "epst_set_consistent.tsv",
  "epst_set_noisy.tsv"
)

preferred_paths <- file.path("output", preferred)
ordered <- c(preferred_paths[file.exists(preferred_paths)], setdiff(all_tsv, preferred_paths))

sanitize_sheet <- function(x) {
  x <- sub("\\.tsv$", "", basename(x))
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- substr(x, 1, 31)
  x
}

sheet_names <- vapply(ordered, sanitize_sheet, character(1))
# Ensure unique sheet names
if (any(duplicated(sheet_names))) {
  idx <- ave(seq_along(sheet_names), sheet_names, FUN = seq_along)
  sheet_names <- ifelse(idx == 1, sheet_names, substr(paste0(sheet_names, "_", idx), 1, 31))
}

sheets <- list()
for (i in seq_along(ordered)) {
  raw_sheet <- readr::read_tsv(ordered[i], show_col_types = FALSE, progress = FALSE)
  sheets[[sheet_names[i]]] <- append_gene_names(raw_sheet)
}

out_xlsx <- file.path("output", "Table_S2.ePST_summary.xlsx")
writexl::write_xlsx(sheets, path = out_xlsx)

message("Created: ", out_xlsx)
message("Included sheets: ", paste(names(sheets), collapse = ", "))
