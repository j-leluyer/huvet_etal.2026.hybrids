#!/usr/bin/env python3
"""
Combine SNP-level ASE counts into gene-level ASE estimates across samples.

Features:
- Maps ASE SNP positions to genes using GFF/GTF coordinates
- Aggregates multiple SNPs per gene (default) or keeps one SNP per gene
- Computes pooled ASE ratio (ALT / total)
- Computes a simple empirical-Bayes shrunk ASE ratio using a Beta prior
- Outputs per-gene/per-sample table suitable for downstream modeling

Example:
  python3 script/ase/combine_gene_ase.py \
    --ase-dir /path/to/starW \
    --pattern "*.ASE.table" \
    --gff /path/to/genomic.mito.gtf \
    --output-prefix output/gene_ase \
    --min-depth 10
"""

from __future__ import annotations

import argparse
import glob
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs = {}
    parts = [x.strip() for x in attr_text.strip().split(";") if x.strip()]
    for part in parts:
        if "=" in part:  # GFF3
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip().strip('"')
        else:  # GTF key "value"
            toks = part.split(None, 1)
            if len(toks) == 2:
                k, v = toks
                attrs[k.strip()] = v.strip().strip('"')
    return attrs


def load_genes(
    gff_path: str,
    feature_type: str = "gene",
    id_attr: str = "ID",
    fallback_attrs: Tuple[str, ...] = ("gene_id", "Name", "gene", "locus_tag"),
):
    genes = []
    with open(gff_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _phase, attrs = cols
            if feature != feature_type:
                continue
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            a = parse_attrs(attrs)
            gid = a.get(id_attr)
            if gid is None:
                for f in fallback_attrs:
                    if f in a:
                        gid = a[f]
                        break
            if gid is None:
                gid = f"{chrom}:{s}-{e}"
            genes.append({"gene_id": gid, "chrom": chrom, "start": s, "end": e, "strand": strand})
    return genes


def build_bin_index(genes: List[dict], bin_size: int = 100_000):
    idx = defaultdict(lambda: defaultdict(list))
    for i, g in enumerate(genes):
        b0 = g["start"] // bin_size
        b1 = g["end"] // bin_size
        for b in range(b0, b1 + 1):
            idx[g["chrom"]][b].append(i)
    return idx


def detect_col(cols: List[str], candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def read_ase_table(path: str):
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")

        chrom_col = detect_col(header, ["contig", "chrom", "chr"])
        pos_col = detect_col(header, ["position", "pos", "start"])
        ref_col = detect_col(header, ["refCount", "ref_count", "refAlleleCount"])
        alt_col = detect_col(header, ["altCount", "alt_count", "altAlleleCount"])
        total_col = detect_col(header, ["totalCount", "total_count", "depth", "DP"])

        if chrom_col is None or pos_col is None or ref_col is None or alt_col is None:
            raise ValueError(
                f"Could not detect required columns in {path}. "
                f"Need chrom/position/refCount/altCount-like columns."
            )

        ix = {name: i for i, name in enumerate(header)}
        records = []
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                chrom = parts[ix[chrom_col]]
                pos = int(parts[ix[pos_col]])
                ref = int(parts[ix[ref_col]])
                alt = int(parts[ix[alt_col]])
                if total_col is not None:
                    tot = int(parts[ix[total_col]])
                else:
                    tot = ref + alt
            except (ValueError, IndexError):
                continue
            records.append((chrom, pos, ref, alt, tot))
    return records


def estimate_beta_prior(proportions: List[float], default_strength: float = 20.0):
    """Method-of-moments Beta prior; fallback to weak centered prior if unstable."""
    if len(proportions) < 10:
        mu = 0.5
        k = default_strength
        return mu * k, (1 - mu) * k

    mu = statistics.mean(proportions)
    var = statistics.pvariance(proportions)

    # guard rails
    if var <= 1e-8 or mu <= 1e-4 or mu >= 1 - 1e-4:
        k = default_strength
        mu = min(max(mu, 0.05), 0.95)
        return mu * k, (1 - mu) * k

    k = mu * (1 - mu) / var - 1.0
    if not math.isfinite(k) or k <= 0:
        k = default_strength

    alpha = mu * k
    beta = (1 - mu) * k
    return alpha, beta


def main():
    ap = argparse.ArgumentParser(description="Combine SNP ASE into gene-level ASE")
    ap.add_argument("--ase-dir", required=True, help="Directory containing ASE tables")
    ap.add_argument("--pattern", default="*.ASE.table", help="Glob pattern for ASE tables")
    ap.add_argument("--gff", required=True, help="Gene annotation GFF/GTF")
    ap.add_argument("--feature-type", default="gene", help="Annotation feature type (default: gene)")
    ap.add_argument("--id-attr", default="ID", help="Primary ID attribute in GFF/GTF")
    ap.add_argument("--output-prefix", default="output/gene_ase", help="Output prefix")
    ap.add_argument("--min-depth", type=int, default=10, help="Minimum SNP depth to keep")
    ap.add_argument("--one-snp-per-gene", action="store_true", help="Use one SNP per gene/sample (max depth)")
    ap.add_argument("--bin-size", type=int, default=100000, help="Bin size for interval indexing")
    args = ap.parse_args()

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ase_files = sorted(glob.glob(str(Path(args.ase_dir) / args.pattern)))
    if not ase_files:
        raise SystemExit(f"No ASE files found: {Path(args.ase_dir) / args.pattern}")

    genes = load_genes(args.gff, feature_type=args.feature_type, id_attr=args.id_attr)
    if not genes:
        raise SystemExit("No genes loaded from annotation.")
    gene_idx = build_bin_index(genes, bin_size=args.bin_size)

    # aggregated counts: (sample, gene_id) -> dict
    agg = defaultdict(lambda: {"ref": 0, "alt": 0, "total": 0, "n_snps": 0})

    for f in ase_files:
        sample = Path(f).name
        # strip known suffix patterns
        for sfx in [".ASE.table", ".ase.table", ".tsv", ".txt"]:
            if sample.endswith(sfx):
                sample = sample[: -len(sfx)]
                break

        records = read_ase_table(f)

        # optional one SNP per gene: keep deepest SNP for each gene in sample
        best_for_gene = {}

        for chrom, pos, ref, alt, total in records:
            if total < args.min_depth:
                continue

            b = pos // args.bin_size
            candidates = gene_idx.get(chrom, {}).get(b, [])
            gene_hits = []
            for gi in candidates:
                g = genes[gi]
                if g["start"] <= pos <= g["end"]:
                    gene_hits.append(g)

            if not gene_hits:
                continue

            for g in gene_hits:
                gid = g["gene_id"]
                if args.one_snp_per_gene:
                    cur = best_for_gene.get(gid)
                    if cur is None or total > cur["total"]:
                        best_for_gene[gid] = {"ref": ref, "alt": alt, "total": total}
                else:
                    key = (sample, gid)
                    agg[key]["ref"] += ref
                    agg[key]["alt"] += alt
                    agg[key]["total"] += total
                    agg[key]["n_snps"] += 1

        if args.one_snp_per_gene:
            for gid, d in best_for_gene.items():
                key = (sample, gid)
                agg[key]["ref"] = d["ref"]
                agg[key]["alt"] = d["alt"]
                agg[key]["total"] = d["total"]
                agg[key]["n_snps"] = 1

    # Estimate empirical prior from observed pooled proportions
    props = []
    for d in agg.values():
        if d["total"] > 0:
            props.append(d["alt"] / d["total"])

    a0, b0 = estimate_beta_prior(props)

    out_gene = str(out_prefix) + ".gene_ase.tsv"
    with open(out_gene, "w") as out:
        out.write(
            "sample\tgene_id\tref_count\talt_count\ttotal_count\tn_snps\t"
            "ase_ratio_alt\tase_ratio_ref\tpost_mean_alt\tpost_sd_alt\t"
            "post_ci95_low\tpost_ci95_high\n"
        )

        for (sample, gid), d in sorted(agg.items()):
            ref = d["ref"]
            alt = d["alt"]
            tot = d["total"]
            if tot <= 0:
                continue

            p_alt = alt / tot
            p_ref = ref / tot

            # Empirical-Bayes posterior under Beta(a0,b0) prior
            a_post = a0 + alt
            b_post = b0 + ref
            post_mean = a_post / (a_post + b_post)
            post_var = (a_post * b_post) / (((a_post + b_post) ** 2) * (a_post + b_post + 1))
            post_sd = math.sqrt(max(post_var, 0.0))
            ci_lo = max(0.0, post_mean - 1.96 * post_sd)
            ci_hi = min(1.0, post_mean + 1.96 * post_sd)

            out.write(
                f"{sample}\t{gid}\t{ref}\t{alt}\t{tot}\t{d['n_snps']}\t"
                f"{p_alt:.6f}\t{p_ref:.6f}\t{post_mean:.6f}\t{post_sd:.6f}\t"
                f"{ci_lo:.6f}\t{ci_hi:.6f}\n"
            )

    out_summary = str(out_prefix) + ".summary.tsv"
    genes_tested = len({gid for (_s, gid) in agg.keys()})
    samples_tested = len({s for (s, _gid) in agg.keys()})
    with open(out_summary, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"n_ase_files\t{len(ase_files)}\n")
        out.write(f"n_samples_with_gene_ase\t{samples_tested}\n")
        out.write(f"n_genes_with_ase\t{genes_tested}\n")
        out.write(f"beta_prior_alpha\t{a0:.6f}\n")
        out.write(f"beta_prior_beta\t{b0:.6f}\n")
        out.write(f"mode\t{'one_snp_per_gene' if args.one_snp_per_gene else 'combined_snps_per_gene'}\n")
        out.write(f"min_depth\t{args.min_depth}\n")

    print(f"Wrote: {out_gene}")
    print(f"Wrote: {out_summary}")


if __name__ == "__main__":
    main()
