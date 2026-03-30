#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges", ask = FALSE, update = FALSE)
  if (!requireNamespace("IRanges", quietly = TRUE)) BiocManager::install("IRanges", ask = FALSE, update = FALSE)
  if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer", ask = FALSE, update = FALSE)

  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

ase_dir <- get_arg("--ase-dir", "data/aser")
pattern <- get_arg("--pattern", "*.ASE.table")
gff_file <- get_arg("--gff", "data/genomic.mito.gff")
out_prefix <- get_arg("--out-prefix", "output/ase/ga_ag_diffase")
feature_type <- get_arg("--feature-type", "gene")
min_depth <- as.integer(get_arg("--min-depth", "10"))
min_samples <- as.integer(get_arg("--min-samples-per-group", "3"))

if (!dir.exists(ase_dir)) stop("Missing ASE directory: ", ase_dir)
if (!file.exists(gff_file)) stop("Missing GFF: ", gff_file)

out_dir <- dirname(out_prefix)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sanitize_gff_for_import <- function(in_gff) {
  lines <- readLines(in_gff, warn = FALSE)
  if (length(lines) == 0) return(in_gff)
  fix_line <- function(x) {
    if (startsWith(x, "#") || !nzchar(x)) return(x)
    parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 8) return(x)
    if (!(parts[7] %in% c("+", "-", "*"))) parts[7] <- "*"
    paste(parts, collapse = "\t")
  }
  fixed <- vapply(lines, fix_line, character(1))
  tf <- tempfile(pattern = "gff_sanitized_", fileext = ".gff")
  writeLines(fixed, tf)
  tf
}

gff_for_import <- sanitize_gff_for_import(gff_file)
g <- import(gff_for_import)
genes <- g[g$type == feature_type]
if (length(genes) == 0) stop("No features with type=", feature_type)

gid <- NULL
for (nm in c("gene_id", "ID", "Name", "gene", "locus_tag")) {
  if (nm %in% names(mcols(genes))) {
    gid <- as.character(mcols(genes)[[nm]])
    break
  }
}
if (is.null(gid)) gid <- rep(NA_character_, length(genes))
fallback_gid <- paste0(as.character(seqnames(genes)), ":", start(genes), "-", end(genes))
gid[is.na(gid) | gid == ""] <- fallback_gid[is.na(gid) | gid == ""]
mcols(genes)$gene_id <- gid

# Build gene annotation lookup: gene_id -> gene_name, gene_product
gene_name_vec <- if ("Name" %in% names(mcols(genes))) as.character(mcols(genes)$Name) else rep("", length(genes))
gene_name_vec[is.na(gene_name_vec)] <- ""
gene_product_vec <- if ("description" %in% names(mcols(genes))) as.character(mcols(genes)$description) else rep("", length(genes))
gene_product_vec[is.na(gene_product_vec)] <- ""
gene_annot <- data.table(gene_id = gid, gene_name = gene_name_vec, gene_product = gene_product_vec)
gene_annot <- unique(gene_annot, by = "gene_id")

files <- Sys.glob(file.path(ase_dir, pattern))
if (length(files) == 0) stop("No files matched: ", file.path(ase_dir, pattern))

infer_cross <- function(sample_id) {
  if (grepl("^AG", sample_id)) return("AG")
  if (grepl("^GA", sample_id)) return("GA")
  return(NA_character_)
}

agg_list <- list()

for (f in files) {
  sample_id <- sub("\\.ASE\\.table$", "", basename(f))
  cross <- infer_cross(sample_id)
  if (is.na(cross)) next

  dt <- fread(f)
  req <- c("contig", "position", "refCount", "altCount", "totalCount")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    warning("Skipping ", f, " missing columns: ", paste(miss, collapse = ", "))
    next
  }

  dt <- dt[totalCount >= min_depth & (refCount + altCount) > 0]
  if (nrow(dt) == 0) next

  snv <- GRanges(seqnames = dt$contig, ranges = IRanges(start = dt$position, width = 1))
  hits <- findOverlaps(snv, genes, ignore.strand = TRUE)
  if (length(hits) == 0) next

  use <- dt[queryHits(hits)]
  use[, gene_id := mcols(genes)$gene_id[subjectHits(hits)]]
  use <- use[!is.na(gene_id) & gene_id != ""]
  if (nrow(use) == 0) next

  sid <- sample_id
  cr <- cross
  use[, `:=`(sample_id = sid, cross = cr)]

  # aggregate SNPs -> gene within sample
  gene_dt <- use[, .(
    ref_count = sum(refCount),
    alt_count = sum(altCount),
    total_count = sum(totalCount),
    n_snps = .N
  ), by = .(sample_id, cross, gene_id)]

  agg_list[[length(agg_list) + 1]] <- gene_dt
}

if (length(agg_list) == 0) stop("No AG/GA sample data found.")
gene_sample <- rbindlist(agg_list, fill = TRUE)

# write intermediate per-sample per-gene counts
fwrite(gene_sample, paste0(out_prefix, ".gene_sample_counts.tsv"), sep = "\t")

# per-gene differential ASE: alt/(alt+ref) differs between AG and GA
results <- gene_sample[, {
  d <- .SD
  nAG <- sum(d$cross == "AG")
  nGA <- sum(d$cross == "GA")

  if (nAG < min_samples || nGA < min_samples) {
    list(
      n_AG = nAG,
      n_GA = nGA,
      mean_alt_ratio_AG = NA_real_,
      mean_alt_ratio_GA = NA_real_,
      beta_GA_vs_AG = NA_real_,
      p_value = NA_real_
    )
  } else {
    d$cross <- factor(d$cross, levels = c("AG", "GA"))
    fit <- tryCatch(
      glm(cbind(alt_count, ref_count) ~ cross, family = quasibinomial(), data = d),
      error = function(e) NULL
    )

    mean_AG <- with(d[cross == "AG", ], sum(alt_count) / sum(alt_count + ref_count))
    mean_GA <- with(d[cross == "GA", ], sum(alt_count) / sum(alt_count + ref_count))

    if (is.null(fit)) {
      list(
        n_AG = nAG,
        n_GA = nGA,
        mean_alt_ratio_AG = mean_AG,
        mean_alt_ratio_GA = mean_GA,
        beta_GA_vs_AG = NA_real_,
        p_value = NA_real_
      )
    } else {
      co <- summary(fit)$coefficients
      beta <- if ("crossGA" %in% rownames(co)) unname(co["crossGA", "Estimate"]) else NA_real_
      pval <- if ("crossGA" %in% rownames(co)) unname(co["crossGA", "Pr(>|t|)"]) else NA_real_
      list(
        n_AG = nAG,
        n_GA = nGA,
        mean_alt_ratio_AG = mean_AG,
        mean_alt_ratio_GA = mean_GA,
        beta_GA_vs_AG = beta,
        p_value = pval
      )
    }
  }
}, by = gene_id]

results[, q_value := p.adjust(p_value, method = "BH")]
results <- merge(results, gene_annot, by = "gene_id", all.x = TRUE)
fwrite(results, paste0(out_prefix, ".gene_diffASE_GA_vs_AG.tsv"), sep = "\t")

message("Wrote: ", paste0(out_prefix, ".gene_sample_counts.tsv"))
message("Wrote: ", paste0(out_prefix, ".gene_diffASE_GA_vs_AG.tsv"))
message("Done.")
