#!/usr/bin/env python3
"""
Annotate chr20 truth-vs-pgphase missing SNP rows with PRIMARY_REASON and DETAIL.

Reads an align-stats TSV (comment line + header) and the pgphase candidate TSV used for the diff.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def load_candidates_by_pos(cand_path: Path, chrom: str) -> dict[int, list[dict[str, str]]]:
    by_pos: dict[int, list[dict[str, str]]] = defaultdict(list)
    with open(cand_path) as f:
        hdr = f.readline().strip().split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ix.values()):
                continue
            if p[ix["CHROM"]] != chrom:
                continue
            row = {h: p[ix[h]] for h in hdr if ix[h] < len(p)}
            by_pos[int(p[ix["POS"]])].append(row)
    return by_pos


def classify(
    pos: int,
    ref: str,
    alt: str,
    cand_rows: list[dict[str, str]],
    s: dict[str, str],
) -> tuple[str, str]:
    """Return (PRIMARY_REASON, DETAIL)."""
    ref = ref.upper()
    alt = alt.upper()

    snps = [r for r in cand_rows if r.get("TYPE") == "SNP"]
    exact = [r for r in snps if r.get("REF") == ref and r.get("ALT") == alt]
    if exact:
        return (
            "inconsistent_input_should_have_matched",
            "candidate_tsv_contains_exact_truth_allele",
        )

    if snps:
        bits = []
        for r in snps[:5]:
            bits.append(
                f"{r.get('REF')}->{r.get('ALT')} DP={r.get('DP')} AF={r.get('AF')} {r.get('CATEGORY')}"
            )
        return (
            "pgphase_called_different_allele_at_site",
            "; ".join(bits),
        )

    non_snp = [r for r in cand_rows if r.get("TYPE") != "SNP"]
    if non_snp:
        bits = [f"{r.get('TYPE')} {r.get('REF')}->{r.get('ALT')} {r.get('CATEGORY', '')}" for r in non_snp[:4]]
        return (
            "variant_at_position_but_not_this_snp",
            "; ".join(bits),
        )

    def i_(k: str) -> int:
        try:
            return int(float(s[k]))
        except (KeyError, ValueError):
            return 0

    def f_(k: str) -> float:
        try:
            v = s[k].strip()
            if not v:
                return float("nan")
            return float(v)
        except (KeyError, ValueError):
            return float("nan")

    overlap = i_("OVERLAP_PRIMARY_READS_MQ_GE5")
    obs_r = i_("OBS_REF_MQ_GE5_BQ_GE10")
    obs_a = i_("OBS_ALT_MQ_GE5_BQ_GE10")
    obs_o = i_("OBS_OTHER_MQ_GE5_BQ_GE10")
    obs_t = i_("OBS_TOTAL_MQ_GE5_BQ_GE10")
    dels = i_("READS_DEL_AT_POS")
    lowbq = i_("BASES_SKIP_LOWBQ_AT_POS")
    st_dep = i_("SAMTOOLS_DEPTH")
    mp5 = i_("MPILEUP_DP_MQ5_BQ10")

    pieces: list[str] = []

    if overlap == 0:
        return (
            "no_primary_alignments_mapq_ge_5",
            "pgphase_read_loader_excludes_unmapped_secondary_supplementary_and_low_mapq",
        )

    if obs_t == 0:
        if dels > 0:
            return (
                "no_hiq_base_at_ref_coord",
                f"{dels}_primary_reads_span_REF_POS_via_deletion_skipping_coord;_{lowbq}_bases_below_min_bq10_at_coord",
            )
        return (
            "no_hiq_base_at_ref_coord",
            f"overlap_mq5={overlap}_but_no_ACGT_obs_ge_bq10;_lowbq_skips={lowbq}",
        )

    af = obs_a / obs_t if obs_t else 0.0
    pieces.append(f"pileup_af_alt={af:.4f}_alt={obs_a}_ref={obs_r}_other={obs_o}_total={obs_t}")

    if obs_o > obs_a and obs_o > obs_r:
        return (
            "dominant_third_allele_or_mismatch_in_hiq_bases",
            "_".join(pieces),
        )

    if obs_a == 0:
        if obs_r > 0:
            return (
                "hiq_bams_support_reference_only",
                "_".join(pieces)
                + f";_samtools_depth={st_dep}_mpileup_mq5bq10={mp5}_truth_alt_not_seen",
            )
        return (
            "no_truth_alt_in_hiq_primary_reads",
            "_".join(pieces),
        )

    if obs_a < 2:
        return (
            "alt_depth_below_pgphase_min_alt_depth_2",
            "_".join(pieces),
        )

    # classify_cand_vars_pgphase folds LowAlleleFraction -> LowCoverage then prune removes LOW_COV
    if af < 0.12 and obs_a < 3:
        return (
            "low_af_pruned_via_low_cov_fold",
            "_".join(pieces)
            + ";_initial_LOW_AF_rewritten_to_LOW_COV_then_pruned_alt_lt_3_no_isolated_rescue",
        )

    if af < 0.12 and obs_a == 3:
        return (
            "low_af_alt_eq3_neighbor_or_prune",
            "_".join(pieces)
            + ";_isolated_rescue_requires_alt_ge3_and_no_variant_within_4bp_else_stays_LOW_AF_then_folded",
        )

    if af < 0.12:
        return (
            "low_af_unexpected_with_alt_ge4",
            "_".join(pieces)
            + ";_check_collect_allele_counts_vs_simple_pileup_or_merge_across_chunks",
        )

    if st_dep >= 5 and mp5 < max(3, st_dep // 4):
        pieces.append(f"mpileup_mq5_low_vs_depth_{mp5}_vs_{st_dep}")

    if obs_a >= 2 and af >= 0.12 and obs_t < 5:
        return (
            "alt_reads_high_af_but_few_mq5_primary_observations",
            "_".join(pieces)
            + ";_most_primary_reads_low_mapq_or_coord_not_covered_at_mq5;_compare_SAMTOOLS_DEPTH_vs_MPILEUP_MQ5",
        )

    if obs_a >= 2 and af >= 0.12 and obs_t >= 5:
        return (
            "pileup_passes_basic_gates_still_no_row",
            "_".join(pieces)
            + ";_digars_allele_counter_or_merge_differs_from_simple_cigar_walk_see_collect_var.cpp",
        )

    return (
        "moderate_support_check_classification_pipeline",
        "_".join(pieces),
    )


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stats",
        type=Path,
        default=repo / "test_data/deepVariantTest/chr20_truth_missing_snps_align_stats.tsv",
    )
    ap.add_argument(
        "--candidates",
        type=Path,
        default=repo / "test_data/deepVariantTest/pgphase_out/chr20_shortread_relaxed.tsv",
    )
    ap.add_argument("--chrom", default="chr20")
    ap.add_argument("-o", "--output", type=Path, default=None, help="default: overwrite --stats")
    args = ap.parse_args()

    out_path = args.output or args.stats
    cand_by_pos = load_candidates_by_pos(args.candidates, args.chrom)

    lines_in = args.stats.read_text().splitlines()
    header_idx = next(i for i, ln in enumerate(lines_in) if ln.startswith("CHROM\t"))
    preamble = lines_in[:header_idx]
    header = lines_in[header_idx].split("\t")
    if "PRIMARY_REASON" in header:
        print("Already annotated; aborting.", file=sys.stderr)
        return 1

    new_header = header + ["PRIMARY_REASON", "DETAIL"]
    out_lines = preamble + ["\t".join(new_header)]

    key_cols = {"CHROM", "POS", "REF", "ALT"}
    col_i = {h: header.index(h) for h in key_cols}

    for ln in lines_in[header_idx + 1 :]:
        if not ln.strip():
            continue
        cells = ln.split("\t")
        row_dict = {header[i]: cells[i] if i < len(cells) else "" for i in range(len(header))}
        pos = int(row_dict["POS"])
        ref = row_dict["REF"]
        alt = row_dict["ALT"]
        primary, detail = classify(pos, ref, alt, cand_by_pos.get(pos, []), row_dict)
        esc_detail = detail.replace("\t", " ")
        out_lines.append("\t".join(cells + [primary, esc_detail]))

    out_path.write_text("\n".join(out_lines) + "\n")
    print(f"Wrote {out_path} ({len(out_lines) - len(preamble) - 1} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
