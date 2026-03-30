#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")

  bioc_pkgs <- c("MBASED", "GenomicRanges", "IRanges", "SummarizedExperiment", "rtracklayer", "BiocParallel")
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }

  library(data.table)
  library(MBASED)
  library(GenomicRanges)
  library(IRanges)
  library(SummarizedExperiment)
  library(rtracklayer)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

ase_dir      <- get_arg("--ase-dir", "asre")
ase_pattern  <- get_arg("--pattern", "*.ASE.table")
gff_file     <- get_arg("--gff", NULL)
out_dir      <- get_arg("--out-dir", "output/mbased")
feature_type <- get_arg("--feature-type", "gene")
min_depth    <- as.integer(get_arg("--min-depth", "10"))
num_sim      <- as.integer(get_arg("--num-sim", "100000"))
workers      <- as.integer(get_arg("--workers", "1"))
bp_type      <- tolower(get_arg("--bp-type", "auto"))

if (is.null(gff_file)) {
  if (file.exists("data/genomic.mito.gff")) {
    gff_file <- "data/genomic.mito.gff"
  } else if (file.exists("asre/genomic.mito.gff")) {
    gff_file <- "asre/genomic.mito.gff"
  } else {
    gff_file <- "data/genomic.mito.gff"
  }
}

if (!dir.exists(ase_dir)) stop("Missing ASE directory: ", ase_dir)
if (!file.exists(gff_file)) stop("Missing GFF file: ", gff_file)
if (is.na(min_depth) || min_depth < 0) stop("--min-depth must be >= 0")
if (is.na(num_sim) || num_sim < 0) stop("--num-sim must be >= 0")
if (is.na(workers) || workers < 1) stop("--workers must be >= 1")
if (!(bp_type %in% c("auto", "serial", "multicore", "snow"))) {
  stop("--bp-type must be one of: auto, serial, multicore, snow")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading GFF: ", gff_file)
# Some GFF files use '.' for strand (valid in GFF, but GRanges expects '*').
# Sanitize strand column before import.
sanitize_gff_for_import <- function(in_gff) {
  lines <- readLines(in_gff, warn = FALSE)
  if (length(lines) == 0) return(in_gff)

  fix_line <- function(x) {
    if (startsWith(x, "#") || !nzchar(x)) return(x)
    parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 8) return(x)
    st <- parts[7]
    if (!(st %in% c("+", "-", "*"))) parts[7] <- "*"
    paste(parts, collapse = "\t")
  }

  fixed <- vapply(lines, fix_line, character(1))
  tf <- tempfile(pattern = "gff_sanitized_", fileext = ".gff")
  writeLines(fixed, tf)
  tf
}

gff_for_import <- sanitize_gff_for_import(gff_file)
g <- rtracklayer::import(gff_for_import)
genes <- g[g$type == feature_type]
if (length(genes) == 0) stop("No features of type '", feature_type, "' in GFF")

get_first_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  NULL
}

gid <- get_first_col(mcols(genes), c("gene_id", "ID", "Name", "gene", "locus_tag"))
if (is.null(gid)) {
  gid <- paste0(as.character(seqnames(genes)), ":", start(genes), "-", end(genes))
}
gid <- as.character(gid)
fallback_gid <- paste0(as.character(seqnames(genes)), ":", start(genes), "-", end(genes))
gid[is.na(gid) | gid == ""] <- fallback_gid[is.na(gid) | gid == ""]
mcols(genes)$aseID <- gid

# Build gene annotation lookup: aseID -> gene_name, gene_product
gene_name_vec <- if ("Name" %in% names(mcols(genes))) as.character(mcols(genes)$Name) else rep("", length(genes))
gene_name_vec[is.na(gene_name_vec)] <- ""
gene_product_vec <- if ("description" %in% names(mcols(genes))) as.character(mcols(genes)$description) else rep("", length(genes))
gene_product_vec[is.na(gene_product_vec)] <- ""
gene_annot <- data.table(gene_id = gid, gene_name = gene_name_vec, gene_product = gene_product_vec)
gene_annot <- unique(gene_annot, by = "gene_id")

ase_files <- Sys.glob(file.path(ase_dir, ase_pattern))
if (length(ase_files) == 0) stop("No ASE tables found with pattern: ", file.path(ase_dir, ase_pattern))

message("Found ", length(ase_files), " ASE tables")

make_bpparam <- function(workers, bp_type) {
  if (bp_type == "serial" || workers <= 1) return(SerialParam())

  if (bp_type == "multicore") return(MulticoreParam(workers = workers))
  if (bp_type == "snow") return(SnowParam(workers = workers, type = "SOCK"))

  # auto
  if (.Platform$OS.type == "windows") {
    SnowParam(workers = workers, type = "SOCK")
  } else {
    MulticoreParam(workers = workers)
  }
}

bpparam_obj <- make_bpparam(workers = workers, bp_type = bp_type)
message("BiocParallel backend: ", class(bpparam_obj)[1], " (workers=", workers, ")")

run_one_sample <- function(f) {
  sample_id <- sub("\\.ASE\\.table$", "", basename(f))
  dt <- data.table::fread(f)

  required_cols <- c("contig", "position", "refAllele", "altAllele", "refCount", "altCount", "totalCount")
  miss <- setdiff(required_cols, names(dt))
  if (length(miss) > 0) {
    stop("File ", f, " is missing columns: ", paste(miss, collapse = ", "))
  }

  dt <- dt[totalCount >= min_depth & (refCount + altCount) > 0]
  if (nrow(dt) == 0) {
    message("No SNPs kept for ", sample_id)
    return(NULL)
  }

  snv <- GRanges(
    seqnames = dt$contig,
    ranges = IRanges(start = dt$position, width = 1)
  )

  hits <- findOverlaps(snv, genes, ignore.strand = TRUE)
  if (length(hits) == 0) {
    message("No SNP-gene overlaps for ", sample_id)
    return(NULL)
  }

  use <- dt[queryHits(hits)]
  use[, aseID := mcols(genes)$aseID[subjectHits(hits)]]
  use <- use[!is.na(aseID) & aseID != ""]
  if (nrow(use) == 0) {
    message("No valid gene IDs after overlap for ", sample_id)
    return(NULL)
  }
  use[, snv_id := paste(contig, position, refAllele, altAllele, sep = ":")]
  use <- unique(use, by = c("snv_id", "aseID"))

  rr <- GRanges(
    seqnames = use$contig,
    ranges = IRanges(start = use$position, width = 1),
    aseID = use$aseID,
    allele1 = use$refAllele,
    allele2 = use$altAllele
  )
  names(rr) <- use$snv_id

  se <- SummarizedExperiment(
    assays = list(
      lociAllele1Counts = matrix(use$refCount, ncol = 1, dimnames = list(use$snv_id, sample_id)),
      lociAllele2Counts = matrix(use$altCount, ncol = 1, dimnames = list(use$snv_id, sample_id))
    ),
    rowRanges = rr
  )

  fit <- runMBASED(
    ASESummarizedExperiment = se,
    isPhased = FALSE,
    numSim = num_sim,
    BPPARAM = bpparam_obj
  )

  out <- data.table(
    sample_id = sample_id,
    gene_id = rownames(fit),
    majorAF = as.numeric(assays(fit)$majorAlleleFrequency[, 1]),
    pASE = as.numeric(assays(fit)$pValueASE[, 1]),
    pHet = as.numeric(assays(fit)$pValueHeterogeneity[, 1])
  )
  out[, qASE := p.adjust(pASE, method = "BH")]
  out <- merge(out, gene_annot, by = "gene_id", all.x = TRUE)

  data.table::fwrite(out, file.path(out_dir, paste0(sample_id, ".MBASED.tsv")), sep = "\t")
  out
}

all_res <- rbindlist(lapply(ase_files, run_one_sample), fill = TRUE)
if (nrow(all_res) > 0) {
  data.table::fwrite(all_res, file.path(out_dir, "MBASED_all_samples.tsv"), sep = "\t")
  message("Wrote: ", file.path(out_dir, "MBASED_all_samples.tsv"))
} else {
  message("No MBASED output generated.")
}

message("Done.")
