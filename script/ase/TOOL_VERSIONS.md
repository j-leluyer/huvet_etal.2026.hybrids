# ASE/WASP tool versions

Original workflow versions retained for traceability:

- `trimmomatic` 0.36
- `STAR` 2.7.11a (mapping scripts) and 2.7.10b (genome-index script)
- `samtools` 1.22.1
- `GATK` 4.4.0.0
- `bcftools` 1.23
- `bedtools` 2.30.0
- `tabix` 0.2.6
- `python3` (for helper scripts; version depends on environment)

Notes:
- Versions come from the original PBS scripts under `script/ase/`.
- Current standalone usage assumes tools are available in `PATH`.
