#!/usr/bin/env python3
"""Reconstruct read sequences from a coordinate-projected GAF using the cs tag.

For Phase 0 smoke testing only: builds a FASTQ whose read names match the GAF
so the hybrid pipeline can bridge doubly-mapped reads.
"""
import gzip
import re
import sys

import pysam

CS_TOKEN = re.compile(r"(:\d+|\*[acgtn][acgtn]|\+[acgtn]+|-[acgtn]+)", re.I)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def reconstruct(ref, contig, ref_start, cs):
    read = []
    pos = ref_start
    for tok in CS_TOKEN.findall(cs):
        op = tok[0]
        if op == ":":
            n = int(tok[1:])
            read.append(ref.fetch(contig, pos, pos + n).upper())
            pos += n
        elif op == "*":
            read.append(tok[2].upper())
            pos += 1
        elif op == "+":
            read.append(tok[1:].upper())
        elif op == "-":
            pos += len(tok) - 1
    return "".join(read)


def main(argv):
    gaf_path, ref_path, out_path = argv[1], argv[2], argv[3]
    ref = pysam.FastaFile(ref_path)
    # ref contigs are like CHM13#0#chr19; GAF col1 is chr19 -> suffix match
    suffix_map = {}
    for name in ref.references:
        suffix_map[name.split("#")[-1]] = name

    opener = gzip.open if gaf_path.endswith(".gz") else open
    seen = set()
    written = 0
    with opener(gaf_path, "rt") as fh, open(out_path, "w") as out:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            gaf_contig, ref_start, qname, strand = f[0], int(f[1]), f[3], f[7]
            contig = suffix_map.get(gaf_contig)
            if contig is None:
                continue
            cs = None
            bq = None
            for tag in f[11:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                elif tag.startswith("bq:Z:"):
                    bq = tag[5:]
            if cs is None:
                continue
            seq = reconstruct(ref, contig, ref_start, cs)
            if not seq:
                continue
            if strand == "-":
                seq = revcomp(seq)
            # Deduplicate read names (GAF may list paired/multi rows).
            name = qname
            if name in seen:
                continue
            seen.add(name)
            if bq and len(bq) == len(seq):
                qual = bq if strand != "-" else bq[::-1]
            else:
                qual = "F" * len(seq)
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")
            written += 1
    sys.stderr.write(f"wrote {written} reads\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
