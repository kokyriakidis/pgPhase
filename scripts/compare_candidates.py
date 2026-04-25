#!/usr/bin/env python3
"""Compare pgPhase collect-bam-variation candidates with longcallD verbose output."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


LONGCALLD_SITE_RE = re.compile(r"^CandVarSite:\s+([^:]+):(\d+)\s+(\d+)-([XID])-([0-9]+)")
LONGCALLD_CATE_RE = re.compile(
    r"^CandVarCate-([A-Za-z0-9]):\s+([^:]+):(\d+)\s+(\d+)-([XID])-([0-9]+)\s+\d+\tLow-Depth:\s+\d+\t(.*):\s+\d+"
)

PGPHASE_TO_LONGCALLD_CATE = {
    "LOW_COV": "L",
    "STRAND_BIAS": "B",
    "CLEAN_HET_SNP": "N",
    "CLEAN_HET_INDEL": "I",
    "REP_HET_INDEL": "R",
    "CLEAN_HOM": "H",
    "LOW_AF": "l",
    "NON_VAR": "0",
}


@dataclass(frozen=True, order=True)
class Candidate:
    chrom: str
    pos: int
    kind: str
    ref_len: int
    alt_len: int


@dataclass(frozen=True, order=True)
class CategorizedCandidate:
    candidate: Candidate
    alt: str


def normalize_kind(kind: str) -> str:
    kind = kind.upper()
    if kind in {"X", "SNV", "SNP"}:
        return "SNP"
    if kind in {"I", "INS", "INSERTION"}:
        return "INS"
    if kind in {"D", "DEL", "DELETION"}:
        return "DEL"
    raise ValueError(f"unknown variant type: {kind}")


def pgphase_candidate(row: dict[str, str]) -> Candidate:
    chrom = row.get("chrom") or row.get("CHROM")
    raw_pos = row.get("key_pos") or row.get("POS")
    raw_kind = row.get("type") or row.get("TYPE")
    if chrom is None or raw_pos is None or raw_kind is None:
        raise ValueError("pgPhase TSV is missing chrom/POS/type columns")

    kind = normalize_kind(raw_kind)
    if "ref_len" in row and "alt_len" in row:
        ref_len = int(row["ref_len"])
        alt_len = int(row["alt_len"])
        # Some pgPhase debug TSVs are VCF-padded (INS A->AG, DEL CA->C),
        # while longcallD CandVarSite uses its internal unpadded lengths.
        if kind == "INS" and ref_len == 1 and alt_len > 1:
            ref_len = 0
            alt_len -= 1
        elif kind == "DEL" and alt_len == 1 and ref_len > 1:
            ref_len -= 1
            alt_len = 0
    elif kind == "SNP":
        ref_len = 1
        alt_len = 1
    elif kind == "INS":
        ref_len = 0
        alt_len = len(row.get("ALT", "").replace(".", ""))
    else:
        ref_len = int(row.get("REF", "0")) if row.get("REF", "").isdigit() else len(row.get("REF", ""))
        alt_len = 0

    return Candidate(chrom=chrom, pos=int(raw_pos), kind=kind, ref_len=ref_len, alt_len=alt_len)


def load_pgphase(path: Path, chrom: str | None) -> set[Candidate]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        candidates = {pgphase_candidate(row) for row in reader}
    if chrom:
        candidates = {candidate for candidate in candidates if candidate.chrom == chrom}
    return candidates


def pgphase_category_key(row: dict[str, str]) -> CategorizedCandidate:
    alt = row.get("ALT", "")
    if alt == ".":
        alt = ""
    if normalize_kind(row.get("TYPE", row.get("type", ""))) == "DEL":
        alt = ""
    return CategorizedCandidate(candidate=pgphase_candidate(row), alt=alt)


def load_pgphase_categories(path: Path, chrom: str | None) -> dict[CategorizedCandidate, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        categories: dict[CategorizedCandidate, str] = {}
        for row in reader:
            candidate = pgphase_candidate(row)
            if chrom and candidate.chrom != chrom:
                continue
            raw_category = row.get("CATEGORY")
            if not raw_category:
                continue
            mapped = PGPHASE_TO_LONGCALLD_CATE.get(raw_category)
            if mapped is None:
                # Categories such as NOISY_CANDIDATE/NOISY_RESOLVED have no direct
                # one-letter longcallD CandVarCate code in this stage.
                continue
            categories[pgphase_category_key(row)] = mapped
    return categories


def load_longcalld(path: Path, chrom: str | None) -> set[Candidate]:
    candidates: set[Candidate] = set()
    with path.open() as handle:
        for line in handle:
            match = LONGCALLD_SITE_RE.match(line)
            if match is None:
                continue
            site_chrom, pos, ref_len, raw_kind, alt_len = match.groups()
            if chrom and site_chrom != chrom:
                continue
            candidates.add(
                Candidate(
                    chrom=site_chrom,
                    pos=int(pos),
                    kind=normalize_kind(raw_kind),
                    ref_len=int(ref_len),
                    alt_len=int(alt_len),
                )
            )
    return candidates


def load_longcalld_categories(path: Path, chrom: str | None, stage: str) -> dict[CategorizedCandidate, str]:
    categories: dict[CategorizedCandidate, str] = {}
    in_final_block = stage != "final"
    with path.open() as handle:
        for line in handle:
            if stage == "final" and "After classify var:" in line and "candidate vars" in line:
                in_final_block = True
                continue
            if stage == "final" and in_final_block and not line.startswith("CandVarCate-"):
                in_final_block = False
            if not in_final_block:
                continue
            match = LONGCALLD_CATE_RE.match(line)
            if match is None:
                continue
            cate, site_chrom, pos, ref_len, raw_kind, alt_len, alt = match.groups()
            if chrom and site_chrom != chrom:
                continue
            candidate = Candidate(
                chrom=site_chrom,
                pos=int(pos),
                kind=normalize_kind(raw_kind),
                ref_len=int(ref_len),
                alt_len=int(alt_len),
            )
            if normalize_kind(raw_kind) == "DEL":
                alt = ""
            categories[CategorizedCandidate(candidate=candidate, alt=alt)] = cate
    return categories


def format_candidate(candidate: Candidate) -> str:
    return f"{candidate.chrom}:{candidate.pos}:{candidate.kind}:{candidate.ref_len}:{candidate.alt_len}"


def format_categorized(candidate: CategorizedCandidate) -> str:
    alt = candidate.alt if candidate.alt else "."
    return f"{format_candidate(candidate.candidate)}:{alt}"


def print_examples(title: str, candidates: set[Candidate], limit: int) -> None:
    print(f"{title}: {len(candidates)}")
    for candidate in sorted(candidates)[:limit]:
        print(f"  {format_candidate(candidate)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pgphase", required=True, type=Path, help="pgPhase TSV from collect-bam-variation")
    parser.add_argument("--longcalld", required=True, type=Path, help="longcallD verbose log")
    parser.add_argument("--chrom", help="Restrict comparison to one chromosome/contig")
    parser.add_argument("--show", type=int, default=10, help="Number of mismatch examples to print [10]")
    parser.add_argument(
        "--compare-categories",
        action="store_true",
        help="Also compare per-site CandVarCate labels (requires longcallD log with -V 2 and pgPhase CATEGORY column)",
    )
    parser.add_argument(
        "--category-stage",
        choices=("final", "all"),
        default="final",
        help="Which longcallD CandVarCate lines to compare [final]",
    )
    args = parser.parse_args()

    pgphase = load_pgphase(args.pgphase, args.chrom)
    longcalld = load_longcalld(args.longcalld, args.chrom)
    shared = pgphase & longcalld
    only_pgphase = pgphase - longcalld
    only_longcalld = longcalld - pgphase

    print(f"pgPhase candidates: {len(pgphase)}")
    print(f"longcallD candidates: {len(longcalld)}")
    print(f"shared candidates: {len(shared)}")
    print_examples("only in pgPhase", only_pgphase, args.show)
    print_examples("only in longcallD", only_longcalld, args.show)

    category_ok = True
    if args.compare_categories:
        pg_categories = load_pgphase_categories(args.pgphase, args.chrom)
        longcalld_categories = load_longcalld_categories(args.longcalld, args.chrom, args.category_stage)
        shared_keys = sorted(set(pg_categories) & set(longcalld_categories))
        category_mismatches = [
            (candidate, pg_categories[candidate], longcalld_categories[candidate])
            for candidate in shared_keys
            if pg_categories[candidate] != longcalld_categories[candidate]
        ]

        print(f"pgPhase categorized candidates: {len(pg_categories)}")
        print(f"longcallD categorized candidates: {len(longcalld_categories)}")
        print(f"shared categorized candidates: {len(shared_keys)}")
        print(f"category mismatches: {len(category_mismatches)}")
        for candidate, pg_cate, lcd_cate in category_mismatches[: args.show]:
            print(f"  {format_categorized(candidate)} pgPhase={pg_cate} longcallD={lcd_cate}")
        if len(longcalld_categories) == 0:
            print("warning: no CandVarCate lines found in longcallD log (run longcallD with -V 2)")
        category_ok = len(category_mismatches) == 0

    return 0 if (not only_pgphase and not only_longcalld and category_ok) else 2


if __name__ == "__main__":
    raise SystemExit(main())
