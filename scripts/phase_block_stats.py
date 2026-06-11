#!/usr/bin/env python3
"""Compute VCF-level phasing metrics from a pgphase phased VCF.

Reports the metrics that can be derived from the phased VCF alone (no truth
set required): het/hom counts, phased fraction, phase-block count, and block
N50. This is the lightweight comparison used by the hybrid-model Phase 0
harness to A/B the BAM, graph, and hybrid pipelines on the same input.

For switch/accuracy metrics against a diploid truth assembly, use
evaluate_phase_accuracy.sh instead (requires diplinator).

Usage:
    phase_block_stats.py LABEL phased.vcf [phased.vcf ...]
"""

import sys
from collections import defaultdict


def parse_format(fmt: str, sample: str) -> dict:
    keys = fmt.split(":")
    vals = sample.split(":")
    return dict(zip(keys, vals))


def n50(lengths: list) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def analyze(path: str) -> dict:
    n_records = 0
    n_het = 0
    n_hom = 0
    n_phased = 0
    # phase_set_id -> [positions] to size blocks by span
    block_positions: dict = defaultdict(list)

    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            n_records += 1
            pos = int(fields[1])
            info = parse_format(fields[8], fields[9])
            gt = info.get("GT", "")
            sep = "|" if "|" in gt else "/"
            alleles = gt.split(sep)
            is_het = len(alleles) == 2 and alleles[0] != alleles[1]
            if is_het:
                n_het += 1
            else:
                n_hom += 1
            if "|" in gt and "PS" in info and info["PS"] not in (".", ""):
                n_phased += 1
                block_positions[info["PS"]].append(pos)

    block_spans = [
        max(p) - min(p) + 1 for p in block_positions.values() if len(p) >= 2
    ]

    return {
        "records": n_records,
        "het": n_het,
        "hom": n_hom,
        "phased_het": n_phased,
        "phased_frac": (n_phased / n_het) if n_het else 0.0,
        "blocks": len(block_positions),
        "blocks_ge2": len(block_spans),
        "n50": n50(block_spans),
        "largest_block": max(block_spans) if block_spans else 0,
    }


def main(argv: list) -> int:
    if len(argv) < 3:
        sys.stderr.write(__doc__)
        return 1

    label = argv[1]
    print(f"== {label} ==")
    header = (
        f"{'file':<28} {'het':>7} {'hom':>6} {'phased':>7} "
        f"{'phased%':>8} {'blocks':>7} {'N50':>9} {'largest':>9}"
    )
    print(header)
    for path in argv[2:]:
        try:
            s = analyze(path)
        except OSError as exc:
            sys.stderr.write(f"skip {path}: {exc}\n")
            continue
        name = path.split("/")[-1] if "/" in path else path
        print(
            f"{name:<28} {s['het']:>7} {s['hom']:>6} {s['phased_het']:>7} "
            f"{s['phased_frac'] * 100:>7.1f}% {s['blocks']:>7} "
            f"{s['n50']:>9} {s['largest_block']:>9}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
