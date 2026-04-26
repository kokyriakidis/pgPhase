#!/usr/bin/env python3
"""Compare pgPhase collect-bam-variation candidates with longcallD verbose output.

By default this checks both (1) candidate site identity and (2) per-site category codes
parsed from longcallD CandVarCate lines (requires longcallD built with verbose logging, -V 2).

pgPhase TSV must include **INIT_CAT** (first ``classify_var_cate``) for ``--category-stage initial``;
otherwise the script falls back to CATEGORY for that column.

Category stage ``initial`` (default) uses CandVarCate lines before the first ``After classify var:``
(longcallD's first loop over ``classify_var_cate``). **longcallD does not print those lines for LOW_COV (and skips STRAND_BIAS on ONT before the print)**,
and pgPhase may label sites ``NON_VAR`` while longcallD never printed an initial line — missing LCD
keys are ignored for those cases when using ``initial``. Stage ``final`` uses only the second block
(``After classify var: … candidate vars``),
i.e. survivors after noisy compaction. Stage ``all`` reads every CandVarCate line (duplicates
last-wins, not recommended).

Use --no-compare-categories to only compare candidate sets.

longcallD category letters (collect_var.h LONGCALLD_VAR_CATE_STR / log2 bit index):
  L LOW_COV, B STRAND_BIAS, N CLEAN_HET_SNP, I CLEAN_HET_INDEL, R REP_HET, X (unused),
  S somatic, H CLEAN_HOM, e NOISY_CAND_HET, h NOISY_CAND_HOM, l LOW_AF, 0 NON_VAR
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

LONGCALLD_SITE_RE = re.compile(r"^CandVarSite:\s+([^:]+):(\d+)\s+(\d+)-([XID])-([0-9]+)")
LONGCALLD_CATE_RE = re.compile(
    r"^CandVarCate-([A-Za-z0-9]):\s+([^:]+):(\d+)\s+(\d+)-([XID])-([0-9]+)\s+\d+\tLow-Depth:\s+\d+\t(.*):\s+\d+"
)

# Map pgPhase CATEGORY column to longcallD one-letter code(s). Use frozenset when pgPhase
# collapses multiple longcallD states (e.g. noisy het vs hom).
PGPHASE_TO_LONGCALLD_CATE: dict[str, str | frozenset[str]] = {
    "LOW_COV": "L",
    "STRAND_BIAS": "B",
    "CLEAN_HET_SNP": "N",
    "CLEAN_HET_INDEL": "I",
    "REP_HET_INDEL": "R",
    "CLEAN_HOM": "H",
    "LOW_AF": "l",
    "NON_VAR": "0",
    "NOISY_CAND_HET": "e",
    "NOISY_CAND_HOM": "h",
    # pgPhase promotes large repeat/noisy MSA candidates after classify; longcallD may still print R/e/h.
    "NOISY_RESOLVED": frozenset({"R", "e", "h"}),
}


def lcd_category_matches(
    pg_expect: str | frozenset[str],
    lcd: str,
    *,
    lowcov_lowaf_equivalent: bool = False,
) -> bool:
    if isinstance(pg_expect, frozenset):
        return lcd in pg_expect
    if lowcov_lowaf_equivalent:
        # longcallD: L=LOW_COV, l=LOW_AF; pgPhase sometimes labels LOW_COV where lcd first-pass is l.
        if (pg_expect == "L" and lcd == "l") or (pg_expect == "l" and lcd == "L"):
            return True
    return pg_expect == lcd


def lcd_missing_ok_for_initial_stage(pg_expect: str | frozenset[str]) -> bool:
    """Whether a pgPhase key may legitimately lack a longcallD *initial* CandVarCate line.

    longcallD skips the fprintf for LOW_COV and (on ONT) STRAND_BIAS before that loop's print.
    pgPhase may later classify the same site as NON_VAR (e.g. containment), while longcallD
    never emitted a line for that site in the first pass — so missing LCD is expected.
    """
    if isinstance(pg_expect, frozenset):
        return False
    return pg_expect in ("L", "B", "0")


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


def pgphase_category_token(row: dict[str, str], category_stage: str) -> str:
    """TSV column for longcallD comparison: INIT_CAT (first classify_var_cate) or final CATEGORY."""
    if category_stage == "initial":
        return (row.get("INIT_CAT") or row.get("CATEGORY") or "").strip()
    return (row.get("CATEGORY") or "").strip()


def load_pgphase_categories(
    path: Path, chrom: str | None, category_stage: str
) -> tuple[dict[CategorizedCandidate, str | frozenset[str]], list[tuple[CategorizedCandidate, str]]]:
    """Return (mapped categories, unmapped rows: (key, raw CATEGORY string))."""
    mapped: dict[CategorizedCandidate, str | frozenset[str]] = {}
    unmapped: list[tuple[CategorizedCandidate, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            candidate = pgphase_candidate(row)
            if chrom and candidate.chrom != chrom:
                continue
            raw_category = pgphase_category_token(row, category_stage)
            if not raw_category:
                unmapped.append((pgphase_category_key(row), ""))
                continue
            lcd = PGPHASE_TO_LONGCALLD_CATE.get(raw_category)
            if lcd is None:
                unmapped.append((pgphase_category_key(row), raw_category))
                continue
            mapped[pgphase_category_key(row)] = lcd
    return mapped, unmapped


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
    """Parse CandVarCate lines from longcallD stderr (-V 2).

    initial: lines before the first ``After classify var:`` (classify_var_cate for every site).
    final: only lines after ``After classify var: … candidate vars`` (compacted survivors).
    all: every CandVarCate line in file order (later lines overwrite keys).
    """

    def parse_cate_line(line: str) -> None:
        match = LONGCALLD_CATE_RE.match(line)
        if match is None:
            return
        cate, site_chrom, pos, ref_len, raw_kind, alt_len, alt = match.groups()
        if chrom and site_chrom != chrom:
            return
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

    categories: dict[CategorizedCandidate, str] = {}

    if stage == "initial":
        with path.open() as handle:
            for line in handle:
                if "After classify var:" in line:
                    break
                if line.startswith("CandVarCate-"):
                    parse_cate_line(line)
        return categories

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
            if line.startswith("CandVarCate-"):
                parse_cate_line(line)
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
        "--no-compare-categories",
        action="store_true",
        help="Only compare candidate site sets; skip CandVarCate / CATEGORY checks",
    )
    parser.add_argument(
        "--category-stage",
        choices=("initial", "final", "all"),
        default="initial",
        help="initial=first block (all sites); final=after compaction; all=entire log [initial]",
    )
    parser.add_argument(
        "--lowcov-lowaf-equivalent",
        action="store_true",
        help="Treat longcallD L vs l as matching (LOW_COV vs LOW_AF coarse parity)",
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

    sites_ok = not only_pgphase and not only_longcalld
    category_ok = True

    if args.no_compare_categories:
        return 0 if sites_ok else 2

    pg_categories, unmapped_pg = load_pgphase_categories(args.pgphase, args.chrom, args.category_stage)
    longcalld_categories = load_longcalld_categories(args.longcalld, args.chrom, args.category_stage)

    if not longcalld_categories:
        print(
            "error: no CandVarCate lines found in longcallD log "
            "(run longcallD with verbose enabled, e.g. -V 2, and capture stderr)",
            file=sys.stderr,
        )
        return 2

    if unmapped_pg:
        print(f"error: {len(unmapped_pg)} pgPhase row(s) have missing or unknown CATEGORY for longcallD mapping")
        for key, raw in unmapped_pg[: args.show]:
            print(f"  {format_categorized(key)} CATEGORY={raw!r}", file=sys.stderr)
        return 2

    shared_keys = sorted(set(pg_categories) & set(longcalld_categories))
    only_pg_c = set(pg_categories) - set(longcalld_categories)
    only_lcd_c = set(longcalld_categories) - set(pg_categories)
    if args.category_stage == "initial":
        only_pg_c_expected_missing = {k for k in only_pg_c if lcd_missing_ok_for_initial_stage(pg_categories[k])}
        only_pg_c = {k for k in only_pg_c if k not in only_pg_c_expected_missing}
    else:
        only_pg_c_expected_missing = set()

    category_mismatches = [
        (candidate, pg_categories[candidate], longcalld_categories[candidate])
        for candidate in shared_keys
        if not lcd_category_matches(
            pg_categories[candidate],
            longcalld_categories[candidate],
            lowcov_lowaf_equivalent=args.lowcov_lowaf_equivalent,
        )
    ]

    print(f"pgPhase categorized candidates: {len(pg_categories)}")
    print(f"longcallD categorized candidates: {len(longcalld_categories)}")
    print(f"shared categorized keys: {len(shared_keys)}")
    if only_pg_c_expected_missing:
        print(
            f"only in pgPhase (expected missing initial CandVarCate; see doc — L/B/0): "
            f"{len(only_pg_c_expected_missing)}"
        )
    if only_pg_c:
        print(f"only in pgPhase (categorized keys): {len(only_pg_c)}")
        for c in sorted(only_pg_c)[: args.show]:
            print(f"  {format_categorized(c)}")
    if only_lcd_c:
        print(f"only in longcallD (categorized keys): {len(only_lcd_c)}")
        for c in sorted(only_lcd_c)[: args.show]:
            print(f"  {format_categorized(c)}")
    print(f"category mismatches: {len(category_mismatches)}")
    for candidate, pg_cate, lcd_cate in category_mismatches[: args.show]:
        exp = pg_cate if isinstance(pg_cate, str) else sorted(pg_cate)
        print(f"  {format_categorized(candidate)} pgPhase expects {exp!r} longcallD={lcd_cate!r}")

    category_ok = len(category_mismatches) == 0
    return 0 if (sites_ok and category_ok) else 2


if __name__ == "__main__":
    raise SystemExit(main())
