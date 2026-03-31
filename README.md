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

## ASE / allele-specific expression pipeline

The ASE pipeline quantifies hybrid reciprocal allele imbalance (AG vs GA crosses) at the gene level using parent-informative SNPs.
It is composed of four stages: (1) SNP calling, (2) SNP filtration and marker selection, (3) WASP-corrected mapping and GATK ASEReadCounter counts, and (4) GA-vs-AG statistical test.

### Stage 1 – Genome indexing and read trimming

| Script | Description |
|---|---|
| `script/ase/0-genome_index.pbs` | Build STAR and GATK genome index (combined *C. gigas* nuclear + mitochondrial reference). |
| `script/ase/1-trimmomatic.pbs` | Quality-trim raw reads (Trimmomatic PE, ILLUMINACLIP, LEADING/TRAILING 20, SLIDINGWINDOW 30:30, MINLEN 60). |

### Stage 2 – Mapping and variant calling (pure parental lines AA / GG)

Pure-line reads are first mapped without WASP to build a genotype resource used later for marker selection.

| Script | Description |
|---|---|
| `script/ase/3-star_mapping_pures.pbs` | STAR mapping of pure lines (AA, GG) — single-mapping only (`--outFilterMultimapNmax 1`), MAPQ ≥ 10, sorted BAM. |
| `script/ase/4-prepare_bam-cgigas.pbs` | GATK BAM preparation: `CleanSam` → `MarkDuplicates` (remove duplicates) → `BuildBamIndex` → `SplitNCigarReads`. |
| `script/ase/5-variant_calling.gatk.pbs` | Per-sample GATK `HaplotypeCaller` in GVCF mode (`-ERC GVCF`). |
| `script/ase/6-combine_gvcf.pbs` | `CombineGVCFs` across all pure samples, then `GenotypeGVCFs` to produce a joint raw VCF. |

Batch submission:

```bash
# substitute __BASE__ with each sample name
cat script/ase/3-star_mapping_pures.pbs | sed 's/__BASE__/SAMPLE_ID/g' > jobs/map_SAMPLE_ID.pbs
qsub jobs/map_SAMPLE_ID.pbs
```

### Stage 3 – SNP filtration and marker selection

The main filtration script `script/ase/7-filter_AF_gvcf.pbs` is fully resumable.
It applies a series of ordered filters and produces the WASP-compatible marker VCF used in the hybrid mapping step.

**Steps inside the script:**

| Step | Description |
|---|---|
| 2 | Extract biallelic SNPs (`bcftools view -v snps -m2 -M2`). |
| 3 | Diagnose DP vs GQ masking impact; set low-quality genotypes to missing (`FMT/DP < MIN_DP` or `FMT/GQ < MIN_GQ` → `.`). |
| 5 | Fill tags: `AC, AN, AF, MAF, F_MISSING`. |
| 6 | Site-level QC: `QUAL ≥ MIN_QUAL`, `INFO/DP ≥ MIN_INFO_DP`, `F_MISSING ≤ MAX_MISSING`, `AC ≥ 2 && AC ≤ AN-2`. |
| 7 | Optional allelic-balance (AB) masking in parents: mask heterozygous genotypes with `DP ≥ MIN_AB_DP` and `AB ∉ [AB_LO, AB_HI]`. |
| 8 | Keep only simple A/C/G/T alleles (required for WASP). |
| 9 | Auto-detect AA / GG sample names from VCF header. |
| 10 | Split to AA and GG subsets; recompute AF per group. |
| 11 | Select parent-informative loci by `|AF_AA − AF_GG| ≥ MIN_DAF` (default 0.40). |
| 12 | Re-filter selected markers on hybrid samples only: mask `FMT/DP < HYB_MIN_DP`, remove sites with `F_MISSING > HYB_MAX_MISSING`, require `AC ≥ 1`. |
| 13 | Build final STAR/WASP VCF (FORMAT fields stripped) and SNP list files. |

**Threshold profiles:**

| Profile | QUAL | INFO/DP | F_MISSING | DP/GT | GQ/GT | AB DP | AB range |
|---|---|---|---|---|---|---|---|
| `permissive` | ≥ 15 | ≥ 8 | ≤ 0.50 | ≥ 8 | ≥ 15 | ≥ 12 | [0.10–0.90] |
| `balanced` (default) | ≥ 20 | ≥ 10 | ≤ 0.30 | ≥ 10 | ≥ 20 | ≥ 8 | [0.15–0.85] |
| `strict` | ≥ 20 | ≥ 12 | ≤ 0.20 | ≥ 12 | ≥ 20 | ≥ 20 | [0.20–0.80] |

**Run (resumable from step 7, balanced profile, no AB masking):**

```bash
PROFILE=balanced START_STEP=7 ENABLE_AB_MASKING=0 \
MIN_DAF=0.40 MAX_MISS_GROUP=0.50 HYB_MIN_DP=8 HYB_MAX_MISSING=0.50 \
qsub script/ase/7-filter_AF_gvcf.pbs
```

**Full rerun from raw biallelic SNP extraction:**

```bash
PROFILE=balanced START_STEP=2 ENABLE_AB_MASKING=1 qsub script/ase/7-filter_AF_gvcf.pbs
```

**Key runtime parameters:**

- `START_STEP` – resume point (default `7`); set to `2` for full run.
- `ENABLE_AB_MASKING` – `0` (default) skip AB masking; `1` apply it.
- `PROFILE` – `permissive` / `balanced` / `strict`.
- `MIN_DAF` – `|AF_AA − AF_GG|` threshold (default `0.40`).
- `MAX_MISS_GROUP` – max missingness per pure group before deltaAF join (default `0.50`).
- `HYB_MIN_DP` – hybrid genotype DP threshold in step 12 (default `8`).
- `HYB_MAX_MISSING` – max site missingness across hybrids in step 12 (default `0.50`).
- `AB_MAX_MISSING` – post-AB masking missingness threshold (default `1.00`, effectively off).

### Stage 4 – WASP-corrected mapping and GATK ASEReadCounter (hybrids)

Hybrid reads (AG, GA) are re-mapped with STAR's WASP mode to eliminate reference-mapping bias at marker SNPs, then allele counts are extracted per sample with GATK `ASEReadCounter`.

| Script | Description |
|---|---|
| `script/ase/9-star_wasp_mapping.pbs` | STAR mapping with `--varVCFfile` (hybrid marker VCF) and `--waspOutputMode SAMtag`; keep primary aligned reads, MAPQ ≥ 10; coordinate-sort and index. |
| `script/ase/11-asre_counts.pbs` | GATK BAM preparation on WASP BAMs: `CleanSam` → `MarkDuplicates` (remove duplicates) → `BuildBamIndex` → `SplitNCigarReads`. Then per-sample bcftools filter to heterozygous SNPs, followed by GATK `ASEReadCounter` (`--min-depth 10`, `--min-mapping-quality 20`, `--min-base-quality 20`). Produces per-sample `.ASE.table` files. |

Batch submission (using job submission helpers in `script/ase/jobs/`):

```bash
qsub script/ase/jobs/9-starW-mapping.pbs    # generates and submits per-sample STAR-WASP jobs
qsub script/ase/jobs/13-run_mbased_hpc.pbs  # HPC wrapper for run_mbased_from_ase_gff.R
qsub script/ase/jobs/14-run_diffase_hpc.pbs # HPC wrapper for run_diff_ase_ga_ag.R
```

### Stage 5 – MBASED per-sample ASE and GA-vs-AG differential test

Two R scripts perform downstream ASE quantification:

**`script/ase/run_mbased_from_ase_gff.R`** – run MBASED on each sample independently.

```bash
Rscript script/ase/run_mbased_from_ase_gff.R \
  --ase-dir=data/aser \
  --gff=data/genomic.mito.gff \
  --out-dir=output/ase \
  --min-depth=10 \
  --num-sim=100000 \
  --workers=8
```

Outputs: one `<SAMPLE>.MBASED.tsv` per sample + `output/ase/MBASED_all_samples.tsv`.
Columns: `sample_id`, `gene_id`, `majorAF`, `pASE`, `pHet`, `qASE`, `gene_name`, `gene_product`.

**`script/ase/run_diff_ase_ga_ag.R`** – per-gene quasibinomial GLM comparing allele ratios between AG and GA crosses.

```bash
Rscript script/ase/run_diff_ase_ga_ag.R \
  --ase-dir=data/aser \
  --gff=data/genomic.mito.gff \
  --out-prefix=output/ase/ga_ag_diffase \
  --min-depth=10 \
  --min-samples-per-group=3
```

Model: `glm(cbind(alt_count, ref_count) ~ cross, family = quasibinomial())` with `cross` as AG (reference) vs GA.
The coefficient `beta_GA_vs_AG` is the log-odds-ratio of the ALT allele ratio in GA relative to AG.
p-values are adjusted by Benjamini–Hochberg (`q_value`).

Outputs:
- `output/ase/ga_ag_diffase.gene_sample_counts.tsv` – per-sample per-gene aggregated allele counts.
- `output/ase/ga_ag_diffase.gene_diffASE_GA_vs_AG.tsv` – gene-level test results.
- `output/Table_S5.ASE_diffAG_GA.xlsx` – supplementary Excel workbook (regenerated by `script/combine_ase_diff_to_excel.R`).

**Column descriptions for the differential ASE output:**

| Column | Description |
|---|---|
| `gene_id` | Gene identifier (from GFF `ID` / `Name` / `locus_tag`). |
| `n_AG` / `n_GA` | Number of AG / GA samples with data for the gene. |
| `mean_alt_ratio_AG` | Pooled ALT / total ratio across AG samples. |
| `mean_alt_ratio_GA` | Pooled ALT / total ratio across GA samples. |
| `beta_GA_vs_AG` | Log-odds-ratio from quasibinomial GLM (GA vs AG). |
| `p_value` | Wald p-value for `beta_GA_vs_AG`. |
| `q_value` | BH-adjusted p-value. |
| `gene_name` | Gene name from GFF `Name` attribute. |
| `gene_product` | Gene product description from GFF `description` attribute. |

To regenerate the supplementary Excel file only (without re-running the model):

```bash
Rscript script/combine_ase_diff_to_excel.R
```

## Notes

- The pipeline now auto-detects its project root and avoids hardcoded absolute paths.
- The script checks for missing R packages before running.
- Generated files are written to `output/`.
