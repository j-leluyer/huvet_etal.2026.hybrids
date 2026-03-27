# huvet_etal.2026.hybrids

RNA-seq ortholog analysis workflow used for hybrid oyster transcriptome exploration.

## Repository structure

- `data/`: input metadata, count matrix, and annotation files
- `script/rnaseq_ortho.R`: main analysis pipeline
- `output/`: generated files (ignored by git)

## Requirements

R (>= 4.2 recommended) with required packages used by the script:

`vegan, ape, mclust, reshape2, DESeq2, ggplot2, pheatmap, RColorBrewer, cluster, dplyr, PCAtools, KOGMWU, assertthat, scales, WGCNA, tidyverse, ggpubr, rstatix, MCMCglmm, flashClust, adegenet, variancePartition, BiocParallel, limma, car`

External programs used for ortholog quantification (original workflow versions):

- `bedtools` 2.30.0
- `trimmomatic` 0.36
- `salmon` 1.10.0
- `seqkit` 2.9.0
- `TransDecoder` 5.7.1
- `diamond` 2.1.10
- `python` 3.7

See [script/ortho/TOOL_VERSIONS.md](script/ortho/TOOL_VERSIONS.md) for details.

External programs used for ASE/WASP workflow (original workflow versions):

- `trimmomatic` 0.36
- `STAR` 2.7.11a (mapping) and 2.7.10b (genome indexing)
- `samtools` 1.22.1
- `GATK` 4.4.0.0
- `bcftools` 1.23
- `bedtools` 2.30.0
- `tabix` 0.2.6
- `python3` (helper scripts)

See [script/ase/TOOL_VERSIONS.md](script/ase/TOOL_VERSIONS.md) for details.

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

Run ASE marker filtering (PBS, resumable from selected step):

```bash
qsub script/ase/7-filter_AF_gvcf.pbs
```

Important runtime parameters for `script/ase/7-filter_AF_gvcf.pbs`:

- `START_STEP`: resume point (default `7`) to avoid re-running upstream completed steps (`2,3,5,6`).
- `ENABLE_AB_MASKING`: `0` (default) disables AB masking at step 7; `1` enables AB masking.
- `MIN_DAF`: parent-informative threshold `|AF_AA - AF_GG|` (default `0.40`).
- `MAX_MISS_GROUP`: max missingness in AA and GG group VCFs before deltaAF join (default `0.50`).
- `HYB_MIN_DP`: hybrid genotype DP masking threshold in step 12 (default `8`).
- `HYB_MAX_MISSING`: max missingness allowed across hybrids in step 12 (default `0.50`).
- `AB_MAX_MISSING`: post-AB site missingness threshold (default `1.00`, effectively disabled).

Example: rerun only downstream steps with relaxed thresholds and no AB masking:

```bash
PROFILE=balanced START_STEP=7 ENABLE_AB_MASKING=0 \
MIN_DAF=0.40 MAX_MISS_GROUP=0.50 HYB_MIN_DP=8 HYB_MAX_MISSING=0.50 AB_MAX_MISSING=1.00 \
qsub script/ase/7-filter_AF_gvcf.pbs
```

Example: full rerun from raw biallelic SNP extraction:

```bash
PROFILE=balanced START_STEP=2 ENABLE_AB_MASKING=1 qsub script/ase/7-filter_AF_gvcf.pbs
```

## Notes

- The pipeline now auto-detects its project root and avoids hardcoded absolute paths.
- The script checks for missing R packages before running.
- Generated files are written to `output/`.
