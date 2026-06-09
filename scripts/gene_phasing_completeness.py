#!/usr/bin/env python3
"""Compute % of genes fully phased from a phased VCF + gene BED.

A gene is "fully phased" if all its heterozygous variants share a single
PS (phase set) value.  Genes with 0 or 1 het variant are excluded (they
are trivially phased and uninformative).

Also counts total phased het variants (variants with a PS tag).

Input:
  - Phased VCF (.vcf or .vcf.gz) with GT and PS fields
  - Gene BED (at least 4 columns: chrom, start, end, gene_name)

Output (JSON to stdout):
  {
    "total_het_variants": N,
    "phased_het_variants": N,
    "total_genes": N,
    "genes_with_het": N,
    "genes_fully_phased": N,
    "pct_genes_fully_phased": X.XX
  }

Usage:
  python3 gene_phasing_completeness.py phased.vcf.gz genes.bed
"""
import sys, gzip, collections, json, bisect

def main():
    if len(sys.argv) < 3:
        print("Usage: gene_phasing_completeness.py <phased.vcf[.gz]> <genes.bed>",
              file=sys.stderr)
        sys.exit(1)

    vcf_path = sys.argv[1]
    bed_path = sys.argv[2]

    # ── load gene intervals ──────────────────────────────────────────
    # BED: chrom, start, end, name[, ...]
    genes = []  # list of (chrom, start, end, name)
    gene_by_chrom = collections.defaultdict(list)  # chrom -> [(start, end, idx)]
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            idx = len(genes)
            genes.append((chrom, start, end, name))
            gene_by_chrom[chrom].append((start, end, idx))

    # Sort per-chrom intervals by start for binary search
    for chrom in gene_by_chrom:
        gene_by_chrom[chrom].sort()

    total_genes = len(genes)
    print(f"  Loaded {total_genes} genes from {bed_path}", file=sys.stderr)

    # ── parse VCF for het variants with PS ───────────────────────────
    # For each gene, collect the set of PS values of its het variants.
    gene_ps_sets = [set() for _ in genes]
    gene_het_counts = [0] * len(genes)

    total_het = 0
    phased_het = 0

    _open = gzip.open if vcf_path.endswith(".gz") else open
    with _open(vcf_path, "rt") as fh:
        gt_idx = None
        ps_idx = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # Header line — we'll parse FORMAT per-record
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = int(parts[1])  # 1-based
            fmt_fields = parts[8].split(":")
            sample_fields = parts[9].split(":")

            # Find GT and PS indices in FORMAT
            gt_i = ps_i = -1
            for i, f in enumerate(fmt_fields):
                if f == "GT":
                    gt_i = i
                elif f == "PS":
                    ps_i = i

            if gt_i < 0 or gt_i >= len(sample_fields):
                continue

            gt = sample_fields[gt_i]
            # Check if heterozygous
            sep = "/" if "/" in gt else "|"
            alleles = gt.split(sep)
            if len(alleles) != 2:
                continue
            if alleles[0] == "." or alleles[1] == ".":
                continue
            if alleles[0] == alleles[1]:
                continue  # homozygous

            total_het += 1

            # Get PS value
            ps_val = None
            if ps_i >= 0 and ps_i < len(sample_fields):
                ps_str = sample_fields[ps_i]
                if ps_str != "." and ps_str != "":
                    ps_val = ps_str
                    phased_het += 1
            elif sep == "|":
                # Phased genotype without explicit PS — treat as phased
                # Use a synthetic PS based on position (each such variant
                # is its own block unless grouped by a real PS)
                ps_val = str(pos)
                phased_het += 1

            # Find which genes this variant overlaps (pos is 1-based,
            # BED is 0-based half-open)
            pos0 = pos - 1  # convert to 0-based
            intervals = gene_by_chrom.get(chrom)
            if not intervals:
                continue

            # Binary search: find genes where start <= pos0 < end
            # intervals are sorted by start
            hi = bisect.bisect_right(intervals, (pos0, float('inf'), float('inf')))
            for j in range(max(0, hi - 1), -1, -1):
                s, e, idx = intervals[j]
                if s > pos0:
                    continue
                if pos0 >= e:
                    break
                # pos0 is in [s, e)
                gene_het_counts[idx] += 1
                if ps_val is not None:
                    gene_ps_sets[idx].add(ps_val)

            # Also check forward (overlapping genes with start == pos0)
            for j in range(hi, len(intervals)):
                s, e, idx = intervals[j]
                if s > pos0:
                    break
                if pos0 < e:
                    gene_het_counts[idx] += 1
                    if ps_val is not None:
                        gene_ps_sets[idx].add(ps_val)

    # ── compute gene-level stats ─────────────────────────────────────
    genes_with_het = 0
    genes_fully_phased = 0

    for i in range(len(genes)):
        n_het = gene_het_counts[i]
        if n_het < 2:
            continue  # need at least 2 het variants to evaluate phasing
        genes_with_het += 1
        # Fully phased = all het variants share exactly one PS value
        if len(gene_ps_sets[i]) == 1 and len(gene_ps_sets[i]) <= n_het:
            # All het variants have the same PS
            genes_fully_phased += 1

    pct = 100.0 * genes_fully_phased / genes_with_het if genes_with_het else 0.0

    result = {
        "total_het_variants": total_het,
        "phased_het_variants": phased_het,
        "total_genes": total_genes,
        "genes_with_het": genes_with_het,
        "genes_fully_phased": genes_fully_phased,
        "pct_genes_fully_phased": round(pct, 2),
    }

    json.dump(result, sys.stdout, indent=2)
    print()  # trailing newline

    print(f"  Het variants: {total_het:,} total, {phased_het:,} phased",
          file=sys.stderr)
    print(f"  Genes with >=2 het: {genes_with_het:,}", file=sys.stderr)
    print(f"  Fully phased genes: {genes_fully_phased:,} "
          f"({pct:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
