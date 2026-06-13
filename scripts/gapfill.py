#!/usr/bin/env python3
"""
Additive gap-fill: phased core + reads only a secondary pipeline phased.

The hybrid pipeline's `skip_noisy_kmeans` default deliberately leaves a set of
hard reads (segdup / low-complexity regions) unphased rather than phasing them
at ~36% error. This recovers them additively WITHOUT degrading the clean core:
take the hybrid phased BAM as-is, and for reads a secondary pipeline phased but
hybrid did not, copy the secondary HP/PS onto the hybrid record. The filled
reads' phase sets are shifted into a disjoint namespace (PS += PS_OFFSET) so no
hybrid phase set is renumbered, merged, or re-oriented — purely additive.

NOTE ON PROVENANCE: the bam / graph / hybrid pipeline outputs are all phasings
of the SAME vg-giraffe alignment, not different aligners. The BAM is the
linear surjection of the graph alignment; the GAF is its graph-space form. So
"fill reads" are not reads from another aligner — they are reads that a
different pgphase PIPELINE managed to phase. The fill BAM is only read for its
HP/PS labels by read name; the core already carries every read as a record, so
no realignment is needed even when the fill BAM is unaligned (e.g. the graph
pipeline's phased BAM has no @SQ header — samtools view still yields its tags).

Chaining: run twice with different ps_offset to add two sources to one core,
e.g. core=hybrid, fill1=bam (offset 1e9), fill2=graph (offset 2e9).

On HG002 HiFi chr20 this Pareto-beats the standalone BAM pipeline: more reads
phased, lower Hamming error, fewer switch errors, higher auN, more genome
covered. See CHECKPOINT.md "hybrid-core + gap-fill" for measured results.

Usage:
  python3 gapfill.py CORE_PHASED_BAM FILL_PHASED_BAM OUT_BAM [samtools] [ps_offset]

  CORE_PHASED_BAM   phased BAM whose reads are kept verbatim (e.g. hybrid);
                    must carry all reads as records (the hybrid pipeline does)
  FILL_PHASED_BAM   phased BAM read only for HP/PS of core-unphased reads
                    (e.g. the bam or graph pipeline output; may be unaligned)
  OUT_BAM           output BAM (core reads unchanged + fill reads stamped)
  samtools          samtools binary [default: samtools]
  ps_offset         disjoint PS namespace offset [default: 1000000000]
"""
import subprocess
import sys

DEFAULT_PS_OFFSET = 1_000_000_000


def load_hp_ps(samtools, path):
    """Return {read_name: (hp, ps)} for primary reads carrying an HP tag."""
    phased = {}
    proc = subprocess.Popen([samtools, "view", "-F", "0x900", path],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        hp = ps = None
        for tag in fields[11:]:
            if tag.startswith("HP:i:"):
                hp = int(tag[5:])
            elif tag.startswith("PS:i:"):
                ps = int(tag[5:])
        if hp is not None:
            phased[fields[0]] = (hp, ps)
    proc.wait()
    return phased


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    core_bam = sys.argv[1]
    fill_bam = sys.argv[2]
    out_bam = sys.argv[3]
    samtools = sys.argv[4] if len(sys.argv) > 4 else "samtools"
    ps_offset = int(sys.argv[5]) if len(sys.argv) > 5 else DEFAULT_PS_OFFSET

    print("loading core HP/PS ...", file=sys.stderr, flush=True)
    core = load_hp_ps(samtools, core_bam)
    print(f"  core phased reads: {len(core)}", file=sys.stderr, flush=True)

    print("loading fill HP/PS ...", file=sys.stderr, flush=True)
    fill = load_hp_ps(samtools, fill_bam)
    print(f"  fill phased reads: {len(fill)}", file=sys.stderr, flush=True)

    # Reads to stamp: phased in fill, not phased in core.
    gap = {name: hp_ps for name, hp_ps in fill.items() if name not in core}
    print(f"  gap-fill reads (fill-only): {len(gap)}", file=sys.stderr, flush=True)

    read_proc = subprocess.Popen([samtools, "view", "-h", core_bam],
                                 stdout=subprocess.PIPE, text=True)
    write_proc = subprocess.Popen([samtools, "view", "-b", "-o", out_bam],
                                  stdin=subprocess.PIPE, text=True)

    stamped = 0
    for line in read_proc.stdout:
        if line.startswith("@"):
            write_proc.stdin.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        qname = fields[0]
        # Only stamp primary records that the core left unphased.
        flag = int(fields[1])
        is_secondary_or_supp = flag & 0x900
        has_hp = any(t.startswith("HP:i:") for t in fields[11:])
        if qname in gap and not is_secondary_or_supp and not has_hp:
            hp, ps = gap[qname]
            fields.append(f"HP:i:{hp}")
            if ps is not None:
                fields.append(f"PS:i:{ps + ps_offset}")
            stamped += 1
        write_proc.stdin.write("\t".join(fields) + "\n")

    write_proc.stdin.close()
    write_proc.wait()
    read_proc.wait()

    print(f"  records stamped: {stamped}", file=sys.stderr, flush=True)
    print(f"  wrote {out_bam}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
