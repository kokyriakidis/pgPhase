#!/usr/bin/env python3
"""Analyze overlap between graph snarl sites and BAM-called het variants.

Determines how many snarl sites coincide with BAM het sites at the same
GRCh38 position, and how many BAM sites have no snarl counterpart.
This measures the feasibility of augmenting BAM phasing with graph read
observations.

Usage:
    python3 analyze_site_overlap.py \
        --graph-sites sites.vcf.gz \
        --bam-vcf bam_phased.vcf \
        [--bam-tsv variants.tsv] \
        [--contig chr20] \
        [--output report.txt]

Inputs:
    --graph-sites   Graph sites VCF from vg deconstruct (tabix-indexed).
                    CHROM is pangenome path name (e.g. GRCh38#0#chr20).
    --bam-vcf       BAM pipeline phased VCF output (het sites with GT:PS).
    --bam-tsv       BAM pipeline TSV output (alternative to --bam-vcf).
                    Columns: CHROM POS TYPE REF ALT ... CATEGORY ...
    --contig        Restrict analysis to this contig (e.g. chr20).
"""

import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom(chrom):
    """Extract GRCh38 contig from pangenome path name.

    GRCh38#0#chr20 -> chr20
    CHM13#0#chr20  -> chr20
    chr20          -> chr20
    """
    parts = chrom.split("#")
    return parts[-1]


def parse_graph_sites(path, contig_filter=None):
    """Parse graph sites VCF. Returns dict of (chrom, pos) -> list of (ref, alts)."""
    sites = defaultdict(list)
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom_raw = fields[0]
            chrom = normalize_chrom(chrom_raw)
            if contig_filter and chrom != contig_filter:
                continue
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")
            sites[(chrom, pos)].append((ref, alts))
    return sites


def parse_bam_vcf(path, contig_filter=None):
    """Parse BAM phased VCF. Returns dict of (chrom, pos) -> (ref, alt, is_het, ps)."""
    sites = {}
    with open_maybe_gz(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = normalize_chrom(fields[0])
            if contig_filter and chrom != contig_filter:
                continue
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            # Parse GT from FORMAT/SAMPLE
            is_het = False
            ps = 0
            if len(fields) >= 10:
                fmt_keys = fields[8].split(":")
                fmt_vals = fields[9].split(":")
                fmt = dict(zip(fmt_keys, fmt_vals))
                gt = fmt.get("GT", "0/0")
                gt_alleles = gt.replace("|", "/").split("/")
                if len(gt_alleles) == 2 and gt_alleles[0] != gt_alleles[1]:
                    is_het = True
                if "PS" in fmt:
                    try:
                        ps = int(fmt["PS"])
                    except ValueError:
                        pass

            sites[(chrom, pos)] = (ref, alt, is_het, ps)
    return sites


def parse_bam_tsv(path, contig_filter=None):
    """Parse BAM TSV output. Returns dict of (chrom, pos) -> (ref, alt, category)."""
    sites = {}
    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom = normalize_chrom(fields[col["CHROM"]])
            if contig_filter and chrom != contig_filter:
                continue
            pos = int(fields[col["POS"]])
            ref = fields[col["REF"]]
            alt = fields[col["ALT"]]
            cat = fields[col["CATEGORY"]]
            sites[(chrom, pos)] = (ref, alt, cat)
    return sites


def allele_match(graph_ref, graph_alts, bam_ref, bam_alt):
    """Check if any graph allele matches the BAM allele (exact match after normalization)."""
    # Exact REF+ALT match
    if graph_ref == bam_ref and bam_alt in graph_alts:
        return "exact"

    # Position match only (same POS, different alleles)
    return "position_only"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-sites", required=True, help="Graph sites VCF (from vg deconstruct)")
    parser.add_argument("--bam-vcf", help="BAM pipeline phased VCF")
    parser.add_argument("--bam-tsv", help="BAM pipeline TSV output")
    parser.add_argument("--contig", help="Restrict to this contig")
    parser.add_argument("--output", "-o", help="Output report file (default: stdout)")
    args = parser.parse_args()

    if not args.bam_vcf and not args.bam_tsv:
        parser.error("Provide --bam-vcf or --bam-tsv")

    out = open(args.output, "w") if args.output else sys.stdout

    # Load graph sites
    graph_sites = parse_graph_sites(args.graph_sites, args.contig)
    graph_positions = set(graph_sites.keys())

    # Load BAM sites
    if args.bam_vcf:
        bam_raw = parse_bam_vcf(args.bam_vcf, args.contig)
        # Filter to het sites only
        bam_sites = {}
        for key, (ref, alt, is_het, ps) in bam_raw.items():
            if is_het:
                bam_sites[key] = (ref, alt)
        bam_het_count = len(bam_sites)
        bam_total_count = len(bam_raw)
    else:
        bam_raw = parse_bam_tsv(args.bam_tsv, args.contig)
        het_categories = {"CLEAN_HET_SNP", "CLEAN_HET_INDEL", "CleanHetSnp", "CleanHetIndel",
                          "REP_HET_INDEL", "RepeatHetIndel",
                          "NOISY_CAND_HET", "NoisyCandHet"}
        bam_sites = {}
        for key, (ref, alt, cat) in bam_raw.items():
            if cat in het_categories:
                bam_sites[key] = (ref, alt)
        bam_het_count = len(bam_sites)
        bam_total_count = len(bam_raw)

    bam_positions = set(bam_sites.keys())

    # Overlap analysis
    shared_positions = graph_positions & bam_positions
    graph_only = graph_positions - bam_positions
    bam_only = bam_positions - graph_positions

    # Allele-level matching for shared positions
    exact_match = 0
    position_only = 0
    for pos_key in shared_positions:
        bam_ref, bam_alt = bam_sites[pos_key]
        matched = False
        for graph_ref, graph_alts in graph_sites[pos_key]:
            result = allele_match(graph_ref, graph_alts, bam_ref, bam_alt)
            if result == "exact":
                exact_match += 1
                matched = True
                break
        if not matched:
            position_only += 1

    # Variant type breakdown for BAM-only sites
    bam_only_snp = 0
    bam_only_indel = 0
    for key in bam_only:
        ref, alt = bam_sites[key]
        if len(ref) == 1 and len(alt) == 1:
            bam_only_snp += 1
        else:
            bam_only_indel += 1

    # Graph-only site count (sites in graph VCF with no BAM het at same position)
    # These are sites where graph has a variant but BAM either didn't call it
    # or called it as hom/filtered.
    n_graph = len(graph_positions)
    n_bam_het = bam_het_count
    n_shared = len(shared_positions)
    n_graph_only = len(graph_only)
    n_bam_only = len(bam_only)

    # Report
    out.write("=" * 60 + "\n")
    out.write("Site Overlap Analysis: Graph Sites vs BAM Het Variants\n")
    out.write("=" * 60 + "\n\n")

    if args.contig:
        out.write(f"Contig: {args.contig}\n\n")

    out.write(f"Graph sites (total positions):     {n_graph:>8,}\n")
    out.write(f"BAM het variants:                  {n_bam_het:>8,}\n")
    out.write(f"BAM total variants (all cats):     {bam_total_count:>8,}\n\n")

    out.write("--- Position overlap ---\n")
    def pct(num, denom):
        return f"{100*num/denom:.1f}" if denom > 0 else "N/A"

    out.write(f"Shared positions:                  {n_shared:>8,}  "
              f"({pct(n_shared, n_graph)}% of graph, "
              f"{pct(n_shared, n_bam_het)}% of BAM het)\n")
    out.write(f"  Exact allele match:              {exact_match:>8,}  "
              f"({pct(exact_match, n_shared)}% of shared)\n")
    out.write(f"  Position only (diff alleles):    {position_only:>8,}  "
              f"({pct(position_only, n_shared)}% of shared)\n")
    out.write(f"Graph-only (no BAM het):           {n_graph_only:>8,}  "
              f"({pct(n_graph_only, n_graph)}% of graph)\n")
    out.write(f"BAM-only (no graph site):          {n_bam_only:>8,}  "
              f"({pct(n_bam_only, n_bam_het)}% of BAM het)\n")
    out.write(f"  BAM-only SNPs:                   {bam_only_snp:>8,}\n")
    out.write(f"  BAM-only indels:                 {bam_only_indel:>8,}\n\n")

    out.write("--- Interpretation ---\n")
    bridge_pct = 100 * n_shared / n_bam_het if n_bam_het > 0 else 0
    out.write(f"Bridgeable sites (shared):         {n_shared:>8,}  "
              f"({pct(n_shared, n_bam_het)}% of BAM het sites have graph counterparts)\n")
    out.write(f"  -> Graph reads at these sites can augment BAM phasing.\n")
    out.write(f"BAM sites with no graph coverage:  {n_bam_only:>8,}  "
              f"({pct(n_bam_only, n_bam_het)}%)\n")
    out.write(f"  -> These sites rely on BAM reads only.\n")
    out.write(f"Graph sites with no BAM het:       {n_graph_only:>8,}  "
              f"({pct(n_graph_only, n_graph)}%)\n")
    out.write(f"  -> Potential novel variants or graph artifacts.\n\n")

    # Feasibility assessment
    out.write("--- Feasibility ---\n")
    if bridge_pct >= 50:
        out.write("HIGH: >50% of BAM het sites have graph counterparts.\n")
        out.write("Graph reads can provide substantial additional phasing signal.\n")
    elif bridge_pct >= 20:
        out.write("MODERATE: 20-50% of BAM het sites have graph counterparts.\n")
        out.write("Graph reads can help with stitching but won't transform phasing.\n")
    else:
        out.write("LOW: <20% of BAM het sites have graph counterparts.\n")
        out.write("Limited bridging potential — most BAM sites lack graph coverage.\n")

    if args.output:
        out.close()
        print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
