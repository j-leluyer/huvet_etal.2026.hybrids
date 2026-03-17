# Ortholog count-building pipeline

These standalone scripts build the ortholog count matrix used by downstream analysis in [script/rnaseq_ortho.R](../rnaseq_ortho.R).

They intentionally avoid:
- absolute filesystem paths
- cluster-specific `module` / env-loader paths
- generated `jobs/` submission files

## Required tools

See [TOOL_VERSIONS.md](TOOL_VERSIONS.md) for the tool versions used in the original workflow.
To snapshot exact versions in your current environment, run:

`./script/ortho/00_log_tool_versions.sh output/tool_versions.snapshot.txt`

Commands are expected to be available in `PATH`:
- `bedtools`
- `trimmomatic`
- `salmon`
- `seqkit`
- `TransDecoder.LongOrfs`
- `TransDecoder.Predict`
- `diamond`
- `python3`
- standard Unix tools: `awk`, `sort`, `join`, `comm`, `grep`, `sed`, `gawk`

## Scripts

- `00_log_tool_versions.sh` — records exact versions of required tools in the current environment
- `01_combine_reference.sh` — prepares combined transcript reference and `tx2gene.with_mito.tsv`
- `02_trim_reads.sh` — trims one paired-end sample with Trimmomatic
- `03_build_salmon_index.sh` — builds the Salmon index
- `04_run_salmon_quant.sh` — quantifies one sample with Salmon
- `05_build_ortholog_table.sh` — infers orthologs and appends mitochondrial genes
- `06_build_ortho_counts_matrix.sh` — aggregates Salmon quantifications into `ortho_counts_matrix.tsv`

## Notes

- `05_build_ortholog_table.sh` requires the helper `make_ortholog_matrix.py`, which is not present in this repository. Pass it explicitly with `--matrix-script /path/to/make_ortholog_matrix.py`.
- Scripts write outputs relative to the project root, matching the current repository layout.
- `04_run_salmon_quant.sh` writes the Salmon version used for each sample to `4-mapped-salmon/<sample>/salmon.version.txt`.
