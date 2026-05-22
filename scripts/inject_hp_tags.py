#!/usr/bin/env python3
"""
Inject HP, PS, and YC (color) tags from a phased uBAM into a truth BAM.

After remapping phased reads to the truth assembly (via diplinator),
the HP/PS tags are lost.  This script reads them from the original
phased uBAM and writes them onto the truth BAM records so IGV can
natively color and group by haplotype.

Existing tags (e.g. HO:Z: and hq:i: from diplinator) are preserved.

Tags added:
  HP:i:  haplotype (1 or 2)
  PS:i:  phase set
  YC:Z:  IGV color — blue for HP=1, red for HP=2, grey for unphased

Usage:
  python3 inject_hp_tags.py phased.bam truth.bam output.bam [samtools]
"""
import sys
import subprocess
import os

phased_bam = sys.argv[1]
truth_bam  = sys.argv[2]
out_bam    = sys.argv[3]
samtools   = sys.argv[4] if len(sys.argv) > 4 else "samtools"

# ── load HP/PS from phased uBAM ─────────────────────────────────────
read_phase = {}
proc = subprocess.Popen([samtools, "view", phased_bam],
                        stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    hp = ps = 0
    for tag in f[11:]:
        if tag.startswith("HP:i:"): hp = int(tag[5:])
        elif tag.startswith("PS:i:"): ps = int(tag[5:])
    if hp > 0:
        read_phase[f[0]] = (hp, ps)
proc.wait()
print(f"  Loaded {len(read_phase)} HP/PS entries.", file=sys.stderr)

# ── color map ────────────────────────────────────────────────────────
HP_COLORS = {1: "0,0,255", 2: "255,0,0"}  # blue for HP1, red for HP2
UNPHASED_COLOR = "200,200,200"

# ── stream truth BAM, inject tags, write output ─────────────────────
# Use samtools view -h to get SAM, add tags, pipe back to samtools view -b.
read_proc = subprocess.Popen(
    [samtools, "view", "-h", truth_bam],
    stdout=subprocess.PIPE, text=True)

write_proc = subprocess.Popen(
    [samtools, "view", "-b", "-o", out_bam],
    stdin=subprocess.PIPE, text=True)

n_tagged = 0
n_total = 0
for line in read_proc.stdout:
    if line.startswith("@"):
        write_proc.stdin.write(line)
        continue

    n_total += 1
    fields = line.rstrip("\n").split("\t")
    qname = fields[0]

    # Remove any existing HP/PS/YC tags
    cleaned = []
    for f in fields[11:]:
        if not (f.startswith("HP:") or f.startswith("PS:") or f.startswith("YC:")):
            cleaned.append(f)
    fields = fields[:11] + cleaned

    phase = read_phase.get(qname)
    if phase:
        hp, ps = phase
        fields.append(f"HP:i:{hp}")
        fields.append(f"PS:i:{ps}")
        fields.append(f"YC:Z:{HP_COLORS.get(hp, UNPHASED_COLOR)}")
        n_tagged += 1
    else:
        fields.append(f"YC:Z:{UNPHASED_COLOR}")

    write_proc.stdin.write("\t".join(fields) + "\n")

write_proc.stdin.close()
write_proc.wait()
read_proc.wait()

print(f"  Tagged {n_tagged}/{n_total} reads. Wrote {out_bam}", file=sys.stderr)
