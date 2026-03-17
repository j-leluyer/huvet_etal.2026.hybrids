# Tool versions

Original workflow versions retained for traceability:

- `bedtools` 2.30.0
- `trimmomatic` 0.36
- `salmon` 1.10.0
- `seqkit` 2.9.0
- `TransDecoder` 5.7.1
- `diamond` 2.1.10
- `python` 3.7

These standalone scripts do not load tools via absolute module paths; they assume equivalent commands are available in `PATH`.

For exact run-time provenance in a new environment, generate a snapshot with:

- `./script/ortho/00_log_tool_versions.sh output/tool_versions.snapshot.txt`

In addition, each Salmon quantification now records the exact Salmon version used in:

- `4-mapped-salmon/<sample>/salmon.version.txt`
