#!/usr/bin/env python3
"""Phase 0: Merge feasibility analysis for hybrid BAM+graph phasing.

Compares two phased BAMs (one from the BAM pipeline, one from the graph
pipeline) by read name.  Reports how many reads fall into each category
and whether haplotype labels agree for doubly-phased reads.

Categories:
  A) Phased by both pipelines
  B) BAM-only (phased by BAM, not by graph)
  C) Graph-only (phased by graph, not by BAM)  <-- rescue candidates
  D) Neither (unphased in both)

For category A, reports haplotype agreement rate per overlapping phase
block pair.  This determines whether labels need flipping and whether
the graph assignments are trustworthy.

Usage:
    python3 merge_feasibility.py \
        --bam-phased bam_pipeline_phased.bam \
        --graph-phased graph_pipeline_phased.bam \
        [--output report.txt]

Both BAMs must have HP:i and PS:i tags on phased reads.  The graph BAM
may be unaligned (flag 0x4); reads are matched by query name.
"""

import argparse
import sys
from collections import defaultdict

import pysam


def load_phase_tags(bam_path):
    """Load HP and PS tags for every read, keyed by query name.

    Returns dict: qname -> (hp, ps)
      hp: 1 or 2 (phased) or 0 (unphased)
      ps: phase set value or -1
    """
    tags = {}
    with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            # Skip secondary/supplementary — one record per read.
            if read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            hp = read.get_tag("HP") if read.has_tag("HP") else 0
            ps = read.get_tag("PS") if read.has_tag("PS") else -1
            tags[qname] = (int(hp), int(ps))
    return tags


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam-phased", required=True,
                        help="BAM pipeline phased BAM (aligned, HP/PS tags)")
    parser.add_argument("--graph-phased", required=True,
                        help="Graph pipeline phased BAM (may be unaligned)")
    parser.add_argument("--output", "-o",
                        help="Output report file (default: stdout)")
    args = parser.parse_args()

    out = open(args.output, "w") if args.output else sys.stdout

    # Load phase tags from both BAMs.
    out.write("Loading BAM pipeline phased BAM...\n")
    out.flush()
    bam_tags = load_phase_tags(args.bam_phased)
    out.write(f"  {len(bam_tags):,} reads loaded\n")

    out.write("Loading graph pipeline phased BAM...\n")
    out.flush()
    graph_tags = load_phase_tags(args.graph_phased)
    out.write(f"  {len(graph_tags):,} reads loaded\n\n")

    # Classify reads into categories.
    all_reads = set(bam_tags.keys()) | set(graph_tags.keys())

    cat_a = []  # Phased by both
    cat_b = []  # BAM-only
    cat_c = []  # Graph-only (rescue candidates)
    cat_d = []  # Neither

    for qname in all_reads:
        bam_hp, bam_ps = bam_tags.get(qname, (0, -1))
        graph_hp, graph_ps = graph_tags.get(qname, (0, -1))
        bam_phased = bam_hp in (1, 2)
        graph_phased = graph_hp in (1, 2)

        if bam_phased and graph_phased:
            cat_a.append(qname)
        elif bam_phased and not graph_phased:
            cat_b.append(qname)
        elif not bam_phased and graph_phased:
            cat_c.append(qname)
        else:
            cat_d.append(qname)

    n_total = len(all_reads)
    n_a = len(cat_a)
    n_b = len(cat_b)
    n_c = len(cat_c)
    n_d = len(cat_d)

    def pct(n, d):
        return f"{100*n/d:.1f}" if d > 0 else "N/A"

    out.write("=" * 60 + "\n")
    out.write("Merge Feasibility: BAM + Graph Phased BAMs\n")
    out.write("=" * 60 + "\n\n")

    out.write(f"Total unique reads:                {n_total:>10,}\n")
    out.write(f"A) Phased by both:                 {n_a:>10,}  ({pct(n_a, n_total)}%)\n")
    out.write(f"B) BAM-only:                       {n_b:>10,}  ({pct(n_b, n_total)}%)\n")
    out.write(f"C) Graph-only (rescue candidates): {n_c:>10,}  ({pct(n_c, n_total)}%)\n")
    out.write(f"D) Neither:                        {n_d:>10,}  ({pct(n_d, n_total)}%)\n\n")

    # For category A: haplotype agreement analysis.
    # Group by (bam_ps, graph_ps) pair to find overlapping phase blocks.
    out.write("--- Haplotype agreement (category A) ---\n\n")

    if n_a == 0:
        out.write("No doubly-phased reads to compare.\n\n")
    else:
        # For each (bam_ps, graph_ps) pair, count agree/disagree.
        block_pairs = defaultdict(lambda: {"agree": 0, "disagree": 0})
        total_agree = 0
        total_disagree = 0

        for qname in cat_a:
            bam_hp, bam_ps = bam_tags[qname]
            graph_hp, graph_ps = graph_tags[qname]
            key = (bam_ps, graph_ps)
            if bam_hp == graph_hp:
                block_pairs[key]["agree"] += 1
                total_agree += 1
            else:
                block_pairs[key]["disagree"] += 1
                total_disagree += 1

        out.write(f"Total doubly-phased reads:         {n_a:>10,}\n")
        out.write(f"  Same HP label:                   {total_agree:>10,}  ({pct(total_agree, n_a)}%)\n")
        out.write(f"  Different HP label:              {total_disagree:>10,}  ({pct(total_disagree, n_a)}%)\n")
        out.write(f"Overlapping block pairs:           {len(block_pairs):>10,}\n\n")

        # Per-block-pair orientation.
        # A block pair is "concordant" if one orientation (agree or disagree)
        # dominates — meaning the labels are either aligned or consistently
        # flipped.  A block pair is "discordant" if neither dominates.
        concordant = 0
        flipped = 0
        discordant = 0
        concordant_reads = 0
        flipped_reads = 0
        discordant_reads = 0

        # Show top block pairs by read count.
        sorted_pairs = sorted(block_pairs.items(),
                              key=lambda x: x[1]["agree"] + x[1]["disagree"],
                              reverse=True)

        for (bam_ps, graph_ps), counts in sorted_pairs:
            n_agree = counts["agree"]
            n_disagree = counts["disagree"]
            n_pair = n_agree + n_disagree
            if n_pair < 2:
                # Single-read pairs are uninformative.
                continue
            majority = max(n_agree, n_disagree)
            minority = min(n_agree, n_disagree)
            concordance = majority / n_pair

            if concordance >= 0.9:
                if n_agree >= n_disagree:
                    concordant += 1
                    concordant_reads += n_pair
                else:
                    flipped += 1
                    flipped_reads += n_pair
            else:
                discordant += 1
                discordant_reads += n_pair

        out.write("Block pair orientation (pairs with ≥2 reads):\n")
        out.write(f"  Concordant (same labels, ≥90%):  {concordant:>10,}  "
                  f"({concordant_reads:,} reads)\n")
        out.write(f"  Flipped (opposite labels, ≥90%): {flipped:>10,}  "
                  f"({flipped_reads:,} reads)\n")
        out.write(f"  Discordant (<90% majority):      {discordant:>10,}  "
                  f"({discordant_reads:,} reads)\n\n")

        # Show top 10 block pairs.
        out.write("Top 10 block pairs by read count:\n")
        out.write(f"  {'BAM_PS':>12}  {'GRAPH_PS':>12}  {'AGREE':>7}  {'DISAGREE':>9}  "
                  f"{'TOTAL':>7}  {'ORIENTATION':>12}\n")
        for (bam_ps, graph_ps), counts in sorted_pairs[:10]:
            n_agree = counts["agree"]
            n_disagree = counts["disagree"]
            n_pair = n_agree + n_disagree
            if n_pair < 2:
                orient = "too_few"
            elif n_agree >= n_disagree:
                orient = f"same ({pct(n_agree, n_pair)}%)"
            else:
                orient = f"flip ({pct(n_disagree, n_pair)}%)"
            out.write(f"  {bam_ps:>12}  {graph_ps:>12}  {n_agree:>7,}  {n_disagree:>9,}  "
                      f"{n_pair:>7,}  {orient:>12}\n")
        out.write("\n")

    # Rescue potential summary.
    out.write("--- Rescue potential ---\n\n")
    bam_phased_total = n_a + n_b
    potential_total = bam_phased_total + n_c
    out.write(f"BAM pipeline phases:               {bam_phased_total:>10,}  "
              f"({pct(bam_phased_total, n_total)}%)\n")
    out.write(f"After graph rescue (potential):     {potential_total:>10,}  "
              f"({pct(potential_total, n_total)}%)\n")
    out.write(f"Reads rescued:                     {n_c:>10,}  "
              f"(+{pct(n_c, n_total)}pp)\n")
    out.write(f"Remaining unphased:                {n_d:>10,}  "
              f"({pct(n_d, n_total)}%)\n\n")

    # Category C breakdown: reads in graph BAM but not phased by BAM.
    # Check if they exist in the BAM at all (mapped but unphased vs absent).
    n_c_in_bam = sum(1 for qname in cat_c if qname in bam_tags)
    n_c_not_in_bam = n_c - n_c_in_bam
    out.write("Graph-only reads breakdown:\n")
    out.write(f"  In BAM (mapped, unphased):       {n_c_in_bam:>10,}\n")
    out.write(f"  Not in BAM (unmapped/absent):     {n_c_not_in_bam:>10,}\n\n")

    out.write("--- Feasibility assessment ---\n\n")
    if n_c > 0 and n_a > 0:
        # Check if block pairs are mostly concordant/flipped (orientable).
        orientable = concordant + flipped
        total_pairs = concordant + flipped + discordant
        if total_pairs > 0 and orientable / total_pairs >= 0.8:
            out.write("FEASIBLE: Block pairs are mostly orientable "
                      f"({orientable}/{total_pairs} = "
                      f"{pct(orientable, total_pairs)}%).\n")
            out.write(f"Graph rescue can add {n_c:,} reads "
                      f"(+{pct(n_c, n_total)}pp coverage).\n")
            if n_c_not_in_bam > 0:
                out.write(f"  {n_c_not_in_bam:,} reads are graph-only "
                          "(not in BAM) — need special handling.\n")
        else:
            out.write("RISKY: Many block pairs are discordant "
                      f"({discordant}/{total_pairs}).\n")
            out.write("Haplotype orientation is unreliable — "
                      "graph rescue may inject errors.\n")
    elif n_c > 0 and n_a == 0:
        out.write("UNKNOWN: No doubly-phased reads to calibrate orientation.\n")
        out.write("Cannot determine if graph HP labels match BAM labels.\n")
    else:
        out.write("NO RESCUE CANDIDATES: Graph pipeline didn't phase any "
                  "reads that BAM missed.\n")

    if args.output:
        out.close()
        print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
