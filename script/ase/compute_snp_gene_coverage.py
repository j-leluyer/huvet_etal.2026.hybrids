#!/usr/bin/env python3
"""
Compute SNP-per-gene coverage from a VCF and a GFF/GTF gene annotation,
and optionally keep one representative SNP per gene.

Usage example:
  python3 script/ase/compute_snp_gene_coverage.py \
    --vcf /path/to/pb.AF_AB_markers.ACGT.vcf.gz \
    --gff /path/to/genomic.mito.gtf \
    --output-prefix output/ase_snp_gene \
    --select-one-snp-per-gene
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def open_text(path: str):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")


def parse_attrs(attr_text: str) -> Dict[str, str]:
    """Parse GFF3 or GTF attributes into a dict."""
    attrs = {}

    # GTF-style: key "value";
    parts = [x.strip() for x in attr_text.strip().split(";") if x.strip()]
    for part in parts:
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip().strip('"')
        else:
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
    with open_text(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start, end, _score, _strand, _phase, attrs = cols
            if feature != feature_type:
                continue
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            a = parse_attrs(attrs)
            gid = a.get(id_attr)
            if gid is None:
                for fa in fallback_attrs:
                    if fa in a:
                        gid = a[fa]
                        break
            if gid is None:
                gid = f"{chrom}:{s}-{e}"

            genes.append(
                {
                    "gene_id": gid,
                    "chrom": chrom,
                    "start": s,
                    "end": e,
                }
            )
    return genes


def build_bin_index(genes: List[dict], bin_size: int = 100_000):
    idx = defaultdict(lambda: defaultdict(list))
    for i, g in enumerate(genes):
        b0 = g["start"] // bin_size
        b1 = g["end"] // bin_size
        for b in range(b0, b1 + 1):
            idx[g["chrom"]][b].append(i)
    return idx


def parse_vcf_record(line: str):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
        return None
    chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
    try:
        pos_i = int(pos)
    except ValueError:
        return None
    return {
        "chrom": chrom,
        "pos": pos_i,
        "id": vid,
        "ref": ref,
        "alt": alt,
        "qual": qual,
        "line": line,
    }


def qual_value(q: str) -> float:
    if q in (".", ""):
        return float("-inf")
    try:
        return float(q)
    except ValueError:
        return float("-inf")


def main():
    ap = argparse.ArgumentParser(description="SNP vs gene coverage from VCF + GFF/GTF")
    ap.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ")
    ap.add_argument("--gff", required=True, help="Input GFF/GTF (gene annotation)")
    ap.add_argument("--feature-type", default="gene", help="Feature type to use (default: gene)")
    ap.add_argument("--id-attr", default="ID", help="Primary gene ID attribute (default: ID)")
    ap.add_argument("--bin-size", type=int, default=100000, help="Bin size for interval lookup")
    ap.add_argument("--output-prefix", default="output/ase_snp_gene", help="Output prefix")
    ap.add_argument(
        "--select-one-snp-per-gene",
        action="store_true",
        help="Select one representative SNP per gene (best QUAL)",
    )
    args = ap.parse_args()

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    genes = load_genes(
        gff_path=args.gff,
        feature_type=args.feature_type,
        id_attr=args.id_attr,
    )
    if not genes:
        raise SystemExit("No genes loaded from annotation. Check --feature-type / --id-attr.")

    idx = build_bin_index(genes, bin_size=args.bin_size)

    gene_counts = [0] * len(genes)
    total_snps_seen = 0
    total_snps_in_gene = 0

    # one SNP per gene (best QUAL)
    best_snp_for_gene: Dict[int, dict] = {}

    with open_text(args.vcf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            rec = parse_vcf_record(line)
            if rec is None:
                continue

            total_snps_seen += 1
            chrom = rec["chrom"]
            pos = rec["pos"]
            b = pos // args.bin_size

            candidate_genes = idx.get(chrom, {}).get(b, [])
            overlapped = []
            for gi in candidate_genes:
                g = genes[gi]
                if g["start"] <= pos <= g["end"]:
                    overlapped.append(gi)

            if not overlapped:
                continue

            total_snps_in_gene += 1
            for gi in overlapped:
                gene_counts[gi] += 1
                if args.select_one_snp_per_gene:
                    cur = best_snp_for_gene.get(gi)
                    if cur is None or qual_value(rec["qual"]) > qual_value(cur["qual"]):
                        best_snp_for_gene[gi] = rec

    # Per-gene table
    per_gene_path = str(out_prefix) + ".per_gene.tsv"
    with open(per_gene_path, "w") as out:
        out.write("gene_id\tchrom\tstart\tend\tsnp_count\n")
        for g, c in zip(genes, gene_counts):
            out.write(f"{g['gene_id']}\t{g['chrom']}\t{g['start']}\t{g['end']}\t{c}\n")

    genes_with_snp = sum(1 for c in gene_counts if c > 0)
    ratio_genes_with_snp = genes_with_snp / len(genes)
    snps_per_gene_mean = statistics.mean(gene_counts)
    snps_per_gene_median = statistics.median(gene_counts)

    summary_path = str(out_prefix) + ".summary.tsv"
    with open(summary_path, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_genes\t{len(genes)}\n")
        out.write(f"total_snps_in_vcf\t{total_snps_seen}\n")
        out.write(f"total_snps_overlapping_genes\t{total_snps_in_gene}\n")
        out.write(f"genes_with_>=1_snp\t{genes_with_snp}\n")
        out.write(f"fraction_genes_with_>=1_snp\t{ratio_genes_with_snp:.6f}\n")
        out.write(f"mean_snps_per_gene\t{snps_per_gene_mean:.6f}\n")
        out.write(f"median_snps_per_gene\t{snps_per_gene_median:.6f}\n")

    if args.select_one_snp_per_gene:
        # Gene -> selected SNP table
        sel_table = str(out_prefix) + ".selected_one_snp_per_gene.tsv"
        with open(sel_table, "w") as out:
            out.write("gene_id\tchrom\tpos\tid\tref\talt\tqual\n")
            for gi, rec in best_snp_for_gene.items():
                g = genes[gi]
                out.write(
                    f"{g['gene_id']}\t{rec['chrom']}\t{rec['pos']}\t{rec['id']}\t"
                    f"{rec['ref']}\t{rec['alt']}\t{rec['qual']}\n"
                )

        # Deduplicated VCF of selected SNPs (some SNPs can map to multiple genes)
        selected_keys = set()
        for rec in best_snp_for_gene.values():
            selected_keys.add((rec["chrom"], rec["pos"], rec["ref"], rec["alt"]))

        selected_vcf = str(out_prefix) + ".selected_one_snp_per_gene.vcf"
        with open_text(args.vcf) as fh, open(selected_vcf, "w") as out:
            for line in fh:
                if line.startswith("#"):
                    out.write(line)
                    continue
                rec = parse_vcf_record(line)
                if rec is None:
                    continue
                key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
                if key in selected_keys:
                    out.write(line)

    print(f"Wrote: {per_gene_path}")
    print(f"Wrote: {summary_path}")
    if args.select_one_snp_per_gene:
        print(f"Wrote: {str(out_prefix)}.selected_one_snp_per_gene.tsv")
        print(f"Wrote: {str(out_prefix)}.selected_one_snp_per_gene.vcf")


if __name__ == "__main__":
    main()
