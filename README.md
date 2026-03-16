# huvet_etal.2026.hybrids

RNA-seq ortholog analysis workflow used for hybrid oyster transcriptome exploration.

## Repository structure

- `data/`: input metadata, count matrix, and annotation files
- `script/rnaseq_ortho.R`: main analysis pipeline
- `output/`: generated files (ignored by git)

## Requirements

R (>= 4.2 recommended) with required packages used by the script:

`vegan, ape, mclust, reshape2, DESeq2, ggplot2, pheatmap, RColorBrewer, cluster, dplyr, PCAtools, KOGMWU, assertthat, scales, WGCNA, tidyverse, ggpubr, rstatix, MCMCglmm, flashClust, adegenet, variancePartition, BiocParallel, limma, car`

## Run

From the project root:

```bash
Rscript script/rnaseq_ortho.R
```

If you launch it from another folder, pass a relative project directory:

```bash
Rscript script/rnaseq_ortho.R .
```

Generate Figure 1 only (combined PCA):

```bash
Rscript script/make_fig1_combined_pca.R
```

Explore module preservation with two distinct workflows:

```bash
Rscript script/explore_module_preservation_pure.R
Rscript script/explore_module_preservation_hybrids.R
```

- `explore_module_preservation_pure.R`: builds a new pure-strain network (AA+GG) and tests preservation GG → AA.
- `explore_module_preservation_hybrids.R`: uses existing hybrid modules from `output/networkConstruction_hybrids_sft10.Rda` and tests preservation AG → GA.

Both scripts write summary tables and preservation plots to `output/`.

Compute dominance effects (nuclear + mitochondrial):

```bash
Rscript script/compute_dominance_effects.R
```

This script writes class summaries, gene-level tables, asymmetry summaries, and dominance plots to `output/`.

Build the ortholog count matrix used by the downstream RNA-seq analysis:

```bash
cd script/ortho
./01_combine_reference.sh
./03_build_salmon_index.sh
./02_trim_reads.sh SAMPLE_BASE
./04_run_salmon_quant.sh SAMPLE_BASE
./05_build_ortholog_table.sh /path/to/project --matrix-script /path/to/make_ortholog_matrix.py
./06_build_ortho_counts_matrix.sh
```

See [script/ortho/README.md](script/ortho/README.md) for the standalone ortho pipeline details. Legacy `jobs/` cluster submission files are excluded from the repository.

## Notes

- The pipeline now auto-detects its project root and avoids hardcoded absolute paths.
- The script checks for missing R packages before running.
- Generated files are written to `output/`.
