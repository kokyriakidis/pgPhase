#!/usr/bin/env python3
"""Compare two hap.py runs to show whether pgphase phasing helps DeepVariant.

Reads the ``*.summary.csv`` files hap.py writes for two arms (baseline =
DeepVariant's own phasing, treatment = DeepVariant fed pgphase HP tags) and
prints SNP and INDEL recall / precision / F1 side by side with the delta.

Usage:
    summarize_happy.py --baseline A.summary.csv --treatment B.summary.csv \
        [--baseline-label STR] [--treatment-label STR]

A positive F1 delta means pgphase phasing improved DeepVariant variant calling.
"""

import argparse
import csv
import sys

# hap.py reports each Type twice (Filter=ALL and Filter=PASS). With --pass-only
# the PASS row is the one to compare.
FILTER = "PASS"
TYPES = ("SNP", "INDEL")
METRICS = (
    ("METRIC.Recall", "Recall"),
    ("METRIC.Precision", "Precision"),
    ("METRIC.F1_Score", "F1"),
)


def load_summary(path: str) -> dict:
    """Return {type: {metric_col: float}} for the PASS rows of a summary.csv."""
    rows: dict = {}
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("Filter") != FILTER:
                    continue
                vtype = row.get("Type")
                if vtype not in TYPES:
                    continue
                rows[vtype] = row
    except FileNotFoundError:
        sys.exit(f"error: summary not found: {path}")
    missing = [t for t in TYPES if t not in rows]
    if missing:
        sys.exit(f"error: {path} missing {FILTER} rows for: {', '.join(missing)}")
    return rows


def as_float(row: dict, col: str) -> float:
    val = row.get(col, "")
    try:
        return float(val)
    except (TypeError, ValueError):
        return float("nan")


def fmt(x: float) -> str:
    return "n/a" if x != x else f"{x:.6f}"


def fmt_delta(x: float) -> str:
    if x != x:
        return "n/a"
    sign = "+" if x >= 0 else ""
    return f"{sign}{x:.6f}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True, help="baseline hap.py summary.csv")
    ap.add_argument("--treatment", required=True, help="treatment hap.py summary.csv")
    ap.add_argument("--baseline-label", default="baseline")
    ap.add_argument("--treatment-label", default="treatment")
    args = ap.parse_args()

    base = load_summary(args.baseline)
    treat = load_summary(args.treatment)

    blabel = args.baseline_label
    tlabel = args.treatment_label
    width = max(len(blabel), len(tlabel), 12)

    header = (f"{'Type':<6} {'Metric':<10} "
              f"{blabel:>{width}} {tlabel:>{width}} {'Delta':>12}")
    print(header)
    print("-" * len(header))

    improved = []
    for vtype in TYPES:
        for col, name in METRICS:
            b = as_float(base[vtype], col)
            t = as_float(treat[vtype], col)
            delta = t - b
            print(f"{vtype:<6} {name:<10} "
                  f"{fmt(b):>{width}} {fmt(t):>{width}} {fmt_delta(delta):>12}")
            if name == "F1" and delta == delta:
                improved.append((vtype, delta))
        print("-" * len(header))

    print()
    for vtype, delta in improved:
        verdict = "improves" if delta > 0 else ("no change" if delta == 0 else "regresses")
        print(f"{vtype}: pgphase HP {verdict} DeepVariant F1 ({fmt_delta(delta)})")


if __name__ == "__main__":
    main()
