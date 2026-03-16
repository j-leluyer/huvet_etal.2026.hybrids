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

## Notes

- The pipeline now auto-detects its project root and avoids hardcoded absolute paths.
- The script checks for missing R packages before running.
- Generated files are written to `output/`.
