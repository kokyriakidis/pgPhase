#!/usr/bin/env python3
"""Rebuild chr20 truth missing SNP alignment stats TSV (no PRIMARY_REASON columns)."""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path

MIN_MAPQ = 5
MIN_BQ = 10
FLAG_BAD = 4 | 256 | 2048


def load_truth(path: Path):
    keys = set()
    meta = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 8:
                continue
            chrom, pos_s, ref, alt_field = t[0], t[1], t[3], t[4]
            qual, filt = t[5], t[6]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            if len(ref) != 1:
                continue
            for al in alt_field.split(","):
                if len(al) != 1:
                    continue
                if ref.upper() not in "ACGT" or al.upper() not in "ACGT":
                    continue
                key = (chrom, pos, ref.upper(), al.upper())
                keys.add(key)
                meta[key] = (qual, filt)
    return keys, meta


def load_cand_snps(path: Path):
    keys = set()
    with open(path) as f:
        hdr = f.readline().strip().split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ix["TYPE"]] != "SNP":
                continue
            keys.add((p[ix["CHROM"]], int(p[ix["POS"]]), p[ix["REF"]], p[ix["ALT"]]))
    return keys


def samtools_depth_batch(bam: Path, positions: list[int]) -> dict[int, int]:
    bed = "\n".join(f"chr20\t{p-1}\t{p}" for p in sorted(positions)) + "\n"
    r = subprocess.run(
        ["samtools", "depth", "-b", "/dev/stdin", "-a", str(bam)],
        input=bed.encode(),
        capture_output=True,
        check=True,
    )
    d: dict[int, int] = {}
    for ln in r.stdout.decode().strip().split("\n"):
        if not ln:
            continue
        _, p, dep = ln.split("\t")
        d[int(p)] = int(dep)
    return d


def mpileup_dp(bam: Path, pos: int, mq: int, bq: int) -> int:
    reg = f"chr20:{pos}-{pos}"
    r = subprocess.run(
        ["samtools", "mpileup", "-r", reg, f"-q{mq}", f"-Q{bq}", str(bam)],
        capture_output=True,
        text=True,
    )
    if r.returncode != 0 or not r.stdout.strip():
        return 0
    return int(r.stdout.split("\t")[3])


def alignment_obs(bam: Path, pos: int, ref_base: str, alt_base: str) -> dict:
    reg = f"chr20:{pos}-{pos}"
    r = subprocess.run(
        ["samtools", "view", str(bam), reg],
        capture_output=True,
        text=True,
        check=True,
    )
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    obs_ref = obs_alt = obs_other = 0
    n_del = n_skip_lowbq = 0
    sum_mq_obs = 0
    n_overlap_primary = 0
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 11:
            continue
        flag = int(p[1])
        if flag & FLAG_BAD:
            continue
        mapq = int(p[4])
        if mapq < MIN_MAPQ:
            continue
        pos_read = int(p[3])
        cigar, seq, qual = p[5], p[9], p[10]
        n_overlap_primary += 1
        ref_coord = pos_read
        qi = 0
        saw = None
        mq_here = None
        for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            ln = int(ln)
            if op in "M=X":
                for _ in range(ln):
                    if ref_coord == pos:
                        if qi >= len(seq):
                            saw = "oob"
                            break
                        qv = ord(qual[qi]) - 33 if qi < len(qual) else -1
                        if qv < MIN_BQ:
                            n_skip_lowbq += 1
                            saw = "lowbq"
                            break
                        qb = seq[qi].upper()
                        mq_here = mapq
                        if qb == ref_base:
                            saw = "ref"
                        elif qb == alt_base:
                            saw = "alt"
                        else:
                            saw = "other"
                        break
                    ref_coord += 1
                    qi += 1
                if saw:
                    break
            elif op == "I":
                qi += ln
            elif op in "DN":
                if op == "D" and ref_coord <= pos < ref_coord + ln:
                    n_del += 1
                    saw = "del"
                    break
                ref_coord += ln
            elif op == "S":
                qi += ln
            elif op == "H":
                pass
            if saw == "oob":
                break
        if saw == "ref":
            obs_ref += 1
            sum_mq_obs += mq_here or 0
        elif saw == "alt":
            obs_alt += 1
            sum_mq_obs += mq_here or 0
        elif saw == "other":
            obs_other += 1
            sum_mq_obs += mq_here or 0

    obs_total = obs_ref + obs_alt + obs_other
    mean_mq = (sum_mq_obs / obs_total) if obs_total else ""
    return {
        "overlap_reads_mq_ge": n_overlap_primary,
        "obs_ref": obs_ref,
        "obs_alt": obs_alt,
        "obs_other": obs_other,
        "obs_total_bq_ge": obs_total,
        "reads_del_at_pos": n_del,
        "bases_skip_lowbq_at_pos": n_skip_lowbq,
        "mean_mapq_over_obs": mean_mq,
    }


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--truth",
        type=Path,
        default=repo / "test_data/deepVariantTest/benchmark/HG003_chr20_truth.vcf.gz",
    )
    ap.add_argument(
        "--candidates",
        type=Path,
        default=repo / "test_data/deepVariantTest/pgphase_out/chr20_shortread_relaxed.tsv",
    )
    ap.add_argument(
        "--bam",
        type=Path,
        default=repo
        / "test_data/deepVariantTest/HG003.novaseq.pcr-free.35x.vg-1.55.0.chr20.bam",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=repo / "test_data/deepVariantTest/chr20_truth_missing_snps_align_stats.tsv",
    )
    args = ap.parse_args()

    truth_set, truth_meta = load_truth(args.truth)
    cand_set = load_cand_snps(args.candidates)
    missing = sorted(truth_set - cand_set)
    positions = sorted({p for _, p, _, _ in missing})
    depth_map = samtools_depth_batch(args.bam, positions)

    hdr = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "TRUTH_QUAL",
        "TRUTH_FILTER",
        "SAMTOOLS_DEPTH",
        "MPILEUP_DP_MQ5_BQ10",
        "MPILEUP_DP_MQ0_BQ10",
        "OBS_REF_MQ_GE5_BQ_GE10",
        "OBS_ALT_MQ_GE5_BQ_GE10",
        "OBS_OTHER_MQ_GE5_BQ_GE10",
        "OBS_TOTAL_MQ_GE5_BQ_GE10",
        "READS_DEL_AT_POS",
        "BASES_SKIP_LOWBQ_AT_POS",
        "OVERLAP_PRIMARY_READS_MQ_GE5",
        "MEAN_MAPQ_OVER_OBS_BASES",
    ]
    note1 = (
        "# OBS_* counts primary alignments (exclude unmap/secondary/supplementary) MAPQ>="
        + str(MIN_MAPQ)
        + " with base QUAL>="
        + str(MIN_BQ)
        + " at POS; dup/QCFAIL reads kept (matches pgphase --include-filtered). "
        "Pgphase short-read loader now defaults min_mapq=0; these columns still use MQ>="
        + str(MIN_MAPQ)
        + " for pileup telemetry."
    )
    note2 = "# Run scripts/annotate_chr20_truth_missing_reasons.py to append PRIMARY_REASON / DETAIL."

    lines_out = [note1, note2, "\t".join(hdr)]
    for chrom, pos, refb, altb in missing:
        qual, filt = truth_meta[(chrom, pos, refb, altb)]
        st_dep = depth_map.get(pos, 0)
        mp5 = mpileup_dp(args.bam, pos, 5, MIN_BQ)
        mp0 = mpileup_dp(args.bam, pos, 0, MIN_BQ)
        o = alignment_obs(args.bam, pos, refb, altb)
        row = [
            chrom,
            pos,
            refb,
            altb,
            qual,
            filt,
            st_dep,
            mp5,
            mp0,
            o["obs_ref"],
            o["obs_alt"],
            o["obs_other"],
            o["obs_total_bq_ge"],
            o["reads_del_at_pos"],
            o["bases_skip_lowbq_at_pos"],
            o["overlap_reads_mq_ge"],
            o["mean_mapq_over_obs"],
        ]
        lines_out.append("\t".join(str(x) for x in row))

    args.output.write_text("\n".join(lines_out) + "\n")
    print(f"Truth alleles {len(truth_set)}, missing {len(missing)}, wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
