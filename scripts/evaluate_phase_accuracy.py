#!/usr/bin/env python3
"""Evaluate pgphase phase sets against a diploid truth assembly.

Reads a diplinator-merged BAM where each read carries:
  HO:Z:MAT or HO:Z:PAT  — haplotype of origin (which haplotype won)
  hq:i:<N>               — diplinator HapQ confidence score
  HP:i:<N>               — pgphase haplotype assignment (injected)
  PS:i:<N>               — pgphase phase set (injected)
"""
import sys, subprocess, collections, csv, json, os, gzip

phased_bam  = sys.argv[1]
truth_bam   = sys.argv[2]
min_mapq    = int(sys.argv[3])
min_hapq    = int(sys.argv[4])
min_reads   = int(sys.argv[5])
chroms_arg  = sys.argv[6]
outdir      = sys.argv[7]
samtools    = sys.argv[8] if len(sys.argv) > 8 else "samtools"
exclude_bed = sys.argv[9] if len(sys.argv) > 9 else ""
sites_vcf   = sys.argv[10] if len(sys.argv) > 10 else ""

chrom_filter = set()
if chroms_arg:
    chrom_filter = {c.strip() for c in chroms_arg.split(",") if c.strip()}

# ── load exclude regions ─────────────────────────────────────────────
# BED file with difficult regions (censat, segdups, etc.)
# Stored as dict: contig -> sorted list of (start, end)
exclude_regions = collections.defaultdict(list)
n_excluded = 0
if exclude_bed:
    with open(exclude_bed) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            exclude_regions[parts[0]].append((int(parts[1]), int(parts[2])))
    # Sort intervals for binary search
    for contig in exclude_regions:
        exclude_regions[contig].sort()
    total_intervals = sum(len(v) for v in exclude_regions.values())
    print(f"  Loaded {total_intervals} exclude regions from {exclude_bed}",
          file=sys.stderr)

import bisect

def in_excluded_region(contig, pos):
    """Check if a position falls within any excluded region."""
    intervals = exclude_regions.get(contig)
    if not intervals:
        return False
    # Binary search: find rightmost interval with start <= pos
    idx = bisect.bisect_right(intervals, (pos, float('inf'))) - 1
    if idx < 0:
        return False
    return intervals[idx][0] <= pos < intervals[idx][1]

# ── load sites VCF positions ──────────────────────────────────────────
# Sites VCF from graph decomposition — used to detect phase sets with
# no het sites (unphaseable).  Positions are in graph/reference coords.
sites_by_chrom = collections.defaultdict(list)
if sites_vcf:
    import gzip as _gzip
    _open = _gzip.open if sites_vcf.endswith(".gz") else open
    with _open(sites_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 3)
            chrom = parts[0]
            # Strip graph-coordinate prefixes (e.g. CHM13#0#chr20 → chr20)
            if "#" in chrom:
                chrom = chrom.rsplit("#", 1)[-1]
            sites_by_chrom[chrom].append(int(parts[1]))
    for ch in sites_by_chrom:
        sites_by_chrom[ch].sort()
    total_sites = sum(len(v) for v in sites_by_chrom.values())
    print(f"  Loaded {total_sites} sites from {sites_vcf}", file=sys.stderr)

def count_sites_in_range(chrom, start, end):
    """Count VCF sites in [start, end) using binary search."""
    positions = sites_by_chrom.get(chrom)
    if not positions:
        return 0
    lo = bisect.bisect_left(positions, start)
    hi = bisect.bisect_left(positions, end)
    return hi - lo

# ── load HP/PS per read from phased uBAM ─────────────────────────────
read_phase = {}
proc = subprocess.Popen([samtools, "view", phased_bam], stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    hp = ps = 0
    for tag in f[11:]:
        if tag.startswith("HP:i:"): hp = int(tag[5:])
        elif tag.startswith("PS:i:"): ps = int(tag[5:])
    if hp > 0:
        read_phase[f[0]] = (hp, ps)
proc.wait()
print(f"  Loaded {len(read_phase)} phased reads.", file=sys.stderr)

# ── parse truth alignments (diplinator-merged BAM) ───────────────────
# Haplotype origin is determined from:
#   1. Contig name suffix (_MATERNAL/_PATERNAL) if present
#   2. HO:Z: tag (added during merge step) as fallback
ReadTruth = collections.namedtuple("ReadTruth", "hap chrom mapq pos hapq")
read_truth = {}
n_mapped = n_low_mapq = n_low_hapq = 0

def parse_contig(rname):
    """Extract (haplotype, bare_chrom) from contig name.

    Handles both suffixed (chr20_MATERNAL) and bare (chr20) names.
    Returns (hap, chrom) where hap is 'MAT'/'PAT' or None.
    """
    if rname.endswith("_MATERNAL"):
        return "MAT", rname[:-9]
    elif rname.endswith("_PATERNAL"):
        return "PAT", rname[:-9]
    return None, rname

proc = subprocess.Popen([samtools, "view", "-F", "0x904", truth_bam],
                        stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    qname, mapq, rname = f[0], int(f[4]), f[2]
    pos = int(f[3])
    if mapq < min_mapq:
        n_low_mapq += 1; continue

    # Determine haplotype from contig name or HO tag
    contig_hap, chrom = parse_contig(rname)
    ho = None
    hapq = 60  # default if tag missing
    for tag in f[11:]:
        if tag.startswith("HO:Z:"): ho = tag[5:]
        elif tag.startswith("hq:i:"): hapq = int(tag[5:])

    hap = contig_hap or ho
    if hap is None:
        continue  # no haplotype information
    if hapq < min_hapq:
        n_low_hapq += 1; continue

    # Chromosome filter uses bare name (e.g. chr20, not chr20_MATERNAL)
    if chrom_filter and chrom not in chrom_filter:
        continue

    # Exclude difficult regions (BED uses full contig name e.g. chr20_MATERNAL)
    if exclude_regions and in_excluded_region(rname, pos):
        n_excluded += 1
        continue

    if qname in read_phase and qname not in read_truth:
        read_truth[qname] = ReadTruth(hap, chrom, mapq, pos, hapq)
        n_mapped += 1
proc.wait()
excl_msg = f", {n_excluded} in excluded regions" if exclude_bed else ""
print(f"  {n_mapped} mapped to truth (MAPQ>={min_mapq}, HapQ>={min_hapq}), "
      f"{n_low_mapq} MAPQ-filtered, {n_low_hapq} HapQ-filtered{excl_msg}.",
      file=sys.stderr)

# ── group by phase set ──────────────────────────────────────────────
ps_counts = collections.defaultdict(lambda: collections.Counter())
ps_chroms = collections.defaultdict(set)
ps_contig_positions = collections.defaultdict(lambda: collections.defaultdict(list))

for qname, t in read_truth.items():
    hp, ps = read_phase[qname]
    ps_counts[ps][(hp, t.hap)] += 1
    ps_chroms[ps].add(t.chrom)
    ps_contig_positions[ps][t.chrom].append(t.pos)

# ── evaluate each phase set ─────────────────────────────────────────
results = []
total_conc = total_disc = total_reads = total_ps = total_perfect = 0

for ps in sorted(ps_counts):
    c = ps_counts[ps]
    hp1_mat, hp1_pat = c.get((1,"MAT"),0), c.get((1,"PAT"),0)
    hp2_mat, hp2_pat = c.get((2,"MAT"),0), c.get((2,"PAT"),0)
    tot = hp1_mat + hp1_pat + hp2_mat + hp2_pat
    if tot < min_reads: continue

    conc_a = hp1_mat + hp2_pat  # HP1=MAT, HP2=PAT
    conc_b = hp1_pat + hp2_mat  # HP1=PAT, HP2=MAT
    if conc_a >= conc_b:
        conc, disc, orient = conc_a, conc_b, "HP1=MAT,HP2=PAT"
    else:
        conc, disc, orient = conc_b, conc_a, "HP1=PAT,HP2=MAT"

    # Span is the max per-contig span (positions from different contigs
    # are not comparable since reads may map to MAT or PAT).
    span = 0
    for contig_pos in ps_contig_positions[ps].values():
        if len(contig_pos) > 1:
            s = max(contig_pos) - min(contig_pos)
            if s > span:
                span = s
    acc = conc / tot

    total_conc += conc; total_disc += disc; total_reads += tot
    total_ps += 1
    if disc == 0: total_perfect += 1

    # Per-contig min/max for BED output
    contig_ranges = {}
    for ch, positions in ps_contig_positions[ps].items():
        contig_ranges[ch] = (min(positions), max(positions))

    # Count sites from graph VCF in this phase set's range.
    # PS value is the graph-coordinate start; use span as approximate end.
    # Compute site density (sites/kb) to distinguish well-covered from sparse.
    n_sites = 0
    site_density = 0.0
    if sites_vcf:
        for ch in ps_chroms[ps]:
            n_sites += count_sites_in_range(ch, ps, ps + max(span, 1))
        if span > 0:
            site_density = round(n_sites / (span / 1000), 2)

    results.append({
        "chrom": ",".join(sorted(ps_chroms[ps])),
        "phase_set": ps,
        "n_reads": tot,
        "span_bp": span,
        "hp1_mat": hp1_mat, "hp1_pat": hp1_pat,
        "hp2_mat": hp2_mat, "hp2_pat": hp2_pat,
        "concordant": conc, "discordant": disc,
        "accuracy": acc, "orientation": orient,
        "contig_ranges": contig_ranges,
        "n_sites": n_sites,
        "site_density": site_density,
    })

# ── build orientation map for per-read annotation ───────────────────
ps_orient = {}
for r in results:
    ps_orient[r["phase_set"]] = r["orientation"]

def is_concordant(hp, truth_hap, orientation):
    """Check if a read agrees with the phase set's best orientation."""
    if orientation == "HP1=MAT,HP2=PAT":
        return (hp == 1 and truth_hap == "MAT") or (hp == 2 and truth_hap == "PAT")
    else:  # HP1=PAT,HP2=MAT
        return (hp == 1 and truth_hap == "PAT") or (hp == 2 and truth_hap == "MAT")

# ── compute switch errors per phase set ──────────────────────────────
# A switch error is a transition between concordant and discordant
# in consecutive reads sorted by truth position within a phase set.
# This counts actual phase switches, not individual discordant reads.
ps_switches = {}
ps_read_list = collections.defaultdict(list)  # ps -> [(pos, is_concordant)]

for qname in read_truth:
    hp, ps = read_phase[qname]
    t = read_truth[qname]
    orient = ps_orient.get(ps)
    if orient is None:
        continue
    conc = is_concordant(hp, t.hap, orient)
    ps_read_list[ps].append((t.pos, conc))

total_switches = 0
total_switch_opportunities = 0
for ps in ps_read_list:
    reads = sorted(ps_read_list[ps], key=lambda x: x[0])
    switches = 0
    for i in range(1, len(reads)):
        if reads[i][1] != reads[i-1][1]:
            switches += 1
    ps_switches[ps] = switches
    total_switches += switches
    total_switch_opportunities += max(len(reads) - 1, 0)

# Attach switch count to results
for r in results:
    r["switches"] = ps_switches.get(r["phase_set"], 0)

# ── per-read output with concordance ────────────────────────────────
per_read_file = os.path.join(outdir, "per_read.tsv.gz")
with gzip.open(per_read_file, "wt") as fh:
    fh.write("read_name\tHP\tPS\ttruth_hap\ttruth_contig\ttruth_mapq\t"
             "hapq\ttruth_pos\tstatus\torientation\n")
    for qname in sorted(read_truth):
        hp, ps = read_phase[qname]
        t = read_truth[qname]
        orient = ps_orient.get(ps)
        if orient is None:
            status = "skipped"
        elif is_concordant(hp, t.hap, orient):
            status = "concordant"
        else:
            status = "DISCORDANT"
        orient_str = orient if orient else "."
        fh.write(f"{qname}\t{hp}\t{ps}\t{t.hap}\t{t.chrom}\t{t.mapq}\t"
                 f"{t.hapq}\t{t.pos}\t{status}\t{orient_str}\n")

# ── per-phase-set TSV ────────────────────────────────────────────────
ps_fields = ["chrom","phase_set","n_reads","span_bp",
             "hp1_mat","hp1_pat","hp2_mat","hp2_pat",
             "concordant","discordant","accuracy","switches","orientation"]
if sites_vcf:
    ps_fields.extend(["n_sites", "site_density"])
per_ps_file = os.path.join(outdir, "per_phase_set.tsv")
with open(per_ps_file, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=ps_fields,
                       extrasaction="ignore")
    w.writeheader()
    for r in results:
        out = dict(r); out["accuracy"] = f"{r['accuracy']:.4f}"
        w.writerow(out)

# ── per-chromosome TSV ──────────────────────────────────────────────
cs = collections.defaultdict(lambda: dict(n_ps=0, perfect=0, conc=0, disc=0, tot=0, spans=[]))
for r in results:
    for ch in r["chrom"].split(","):
        s = cs[ch]; s["n_ps"] += 1
        s["perfect"] += 1 if r["discordant"] == 0 else 0
        s["conc"] += r["concordant"]; s["disc"] += r["discordant"]
        s["tot"] += r["n_reads"]; s["spans"].append(r["span_bp"])

per_chrom_file = os.path.join(outdir, "per_chromosome.tsv")
with open(per_chrom_file, "w") as fh:
    fh.write("chrom\tphase_sets\tperfect_ps\treads\tconcordant\tdiscordant\t"
             "accuracy\tmedian_span_bp\tmax_span_bp\n")
    for ch in sorted(cs):
        s = cs[ch]; acc = s["conc"]/s["tot"] if s["tot"] else 0
        spans = sorted(s["spans"])
        med = spans[len(spans)//2] if spans else 0
        mx = spans[-1] if spans else 0
        fh.write(f"{ch}\t{s['n_ps']}\t{s['perfect']}\t{s['tot']}\t"
                 f"{s['conc']}\t{s['disc']}\t{acc:.4f}\t{med}\t{mx}\n")

# ── worst phase sets ────────────────────────────────────────────────
worst_file = os.path.join(outdir, "worst_phase_sets.tsv")
worst = sorted(results, key=lambda r: r["accuracy"])[:50]
with open(worst_file, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=ps_fields,
                       extrasaction="ignore")
    w.writeheader()
    for r in worst:
        out = dict(r); out["accuracy"] = f"{r['accuracy']:.4f}"
        w.writerow(out)

# ── bad regions BED ─────────────────────────────────────────────────
# Phase sets below threshold accuracy, output as BED9 for IGV.
# One row per contig span of each bad phase set.
BAD_ACCURACY_THRESHOLD = 0.60
bad_bed_file = os.path.join(outdir, "bad_phase_regions.bed")
n_bad_ps = 0
n_bad_rows = 0
with open(bad_bed_file, "w") as fh:
    for r in sorted(results, key=lambda x: x["accuracy"]):
        if r["accuracy"] >= BAD_ACCURACY_THRESHOLD:
            break
        n_bad_ps += 1
        for chrom, (start, end) in r["contig_ranges"].items():
            name = (f"PS={r['phase_set']};acc={r['accuracy']:.2f};"
                    f"n={r['n_reads']};disc={r['discordant']}")
            score = int(r["accuracy"] * 1000)
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\t"
                     f"{start}\t{end}\t255,0,0\n")
            n_bad_rows += 1

print(f"\n  {n_bad_ps} bad phase sets ({n_bad_rows} regions, "
      f"accuracy < {BAD_ACCURACY_THRESHOLD:.0%}) in {bad_bed_file}",
      file=sys.stderr)

# ── N50 of phase block spans ────────────────────────────────────────
all_spans = sorted([r["span_bp"] for r in results], reverse=True)
genome_span = sum(all_spans)
n50 = 0
running = 0
for sp in all_spans:
    running += sp
    if running >= genome_span / 2:
        n50 = sp; break

# ── JSON summary (machine-readable) ─────────────────────────────────
overall_acc = total_conc / total_reads if total_reads else 0
hamming_rate = total_disc / total_reads if total_reads else 0
switch_rate = total_switches / total_switch_opportunities if total_switch_opportunities else 0
pct_perfect = 100 * total_perfect / total_ps if total_ps else 0

# Accuracy split by het site availability
# Detect unphaseable blocks: accuracy not significantly different from 50%.
# Use a simple threshold: a block is "phaseable" if accuracy > 60%.
# This separates blocks where phasing worked from those that are random.
phaseable_conc = phaseable_disc = phaseable_ps = 0
unphaseable_conc = unphaseable_disc = unphaseable_ps = 0
PHASEABLE_THRESHOLD = 0.60
for r in results:
    if r["accuracy"] > PHASEABLE_THRESHOLD:
        phaseable_conc += r["concordant"]
        phaseable_disc += r["discordant"]
        phaseable_ps += 1
    else:
        unphaseable_conc += r["concordant"]
        unphaseable_disc += r["discordant"]
        unphaseable_ps += 1
phaseable_acc = phaseable_conc / (phaseable_conc + phaseable_disc) if (phaseable_conc + phaseable_disc) else 0
unphaseable_acc = unphaseable_conc / (unphaseable_conc + unphaseable_disc) if (unphaseable_conc + unphaseable_disc) else 0

summary = {
    "chroms": chroms_arg if chroms_arg else "all",
    "min_mapq": min_mapq,
    "min_hapq": min_hapq,
    "min_reads_per_ps": min_reads,
    "total_phased_reads": len(read_phase),
    "reads_mapped_to_truth": n_mapped,
    "reads_filtered_mapq": n_low_mapq,
    "reads_filtered_hapq": n_low_hapq,
    "reads_excluded_regions": n_excluded,
    "phase_sets_evaluated": total_ps,
    "phase_sets_perfect": total_perfect,
    "phase_sets_perfect_pct": round(pct_perfect, 2),
    "total_reads_evaluated": total_reads,
    "concordant_reads": total_conc,
    "discordant_reads": total_disc,
    "overall_accuracy": round(overall_acc, 6),
    "hamming_error_rate": round(hamming_rate, 6),
    "switch_errors": total_switches,
    "switch_opportunities": total_switch_opportunities,
    "switch_error_rate": round(switch_rate, 6),
    "phase_block_n50_bp": n50,
    "phase_block_median_span_bp": all_spans[len(all_spans)//2] if all_spans else 0,
    "phase_block_max_span_bp": all_spans[0] if all_spans else 0,
    "per_chrom": {},
}
summary["phaseable"] = {
    "threshold": PHASEABLE_THRESHOLD,
    "phase_sets": phaseable_ps,
    "reads": phaseable_conc + phaseable_disc,
    "concordant": phaseable_conc,
    "discordant": phaseable_disc,
    "accuracy": round(phaseable_acc, 6),
}
summary["unphaseable"] = {
    "phase_sets": unphaseable_ps,
    "reads": unphaseable_conc + unphaseable_disc,
    "concordant": unphaseable_conc,
    "discordant": unphaseable_disc,
    "accuracy": round(unphaseable_acc, 6),
}
for ch in sorted(cs):
    s = cs[ch]; acc = s["conc"]/s["tot"] if s["tot"] else 0
    summary["per_chrom"][ch] = {
        "phase_sets": s["n_ps"], "perfect": s["perfect"],
        "reads": s["tot"], "accuracy": round(acc, 6),
    }

json_file = os.path.join(outdir, "summary.json")
with open(json_file, "w") as fh:
    json.dump(summary, fh, indent=2)

# ── human-readable summary ──────────────────────────────────────────
summary_file = os.path.join(outdir, "summary.txt")
with open(summary_file, "w") as fh:
    lines = [
        "=" * 64,
        "  pgphase — Phase Accuracy Evaluation (diplinator)",
        "=" * 64, "",
        f"  Chromosomes:              {chroms_arg if chroms_arg else 'all'}",
        f"  Min MAPQ:                 {min_mapq}",
        f"  Min HapQ:                 {min_hapq}",
        f"  Min reads/phase set:      {min_reads}", "",
        f"  Phase sets evaluated:     {total_ps}",
        f"  Phase sets perfect:       {total_perfect} ({pct_perfect:.1f}%)",
        f"  Phase block N50:          {n50:,} bp", "",
        f"  Total reads evaluated:    {total_reads:,}",
        f"  Concordant:               {total_conc:,}",
        f"  Discordant:               {total_disc:,}",
        f"  MAPQ-filtered:            {n_low_mapq:,}",
        f"  HapQ-filtered:            {n_low_hapq:,}",
        f"  Excluded (difficult):     {n_excluded:,}", "",
        f"  Overall accuracy:         {100*overall_acc:.2f}%",
        f"  Hamming error rate:       {100*hamming_rate:.2f}%  "
        f"({total_disc:,} discordant reads / {total_reads:,})",
        f"  Switch error rate:        {100*switch_rate:.2f}%  "
        f"({total_switches:,} switches / {total_switch_opportunities:,} pairs)",
    ]
    lines += [
        "",
        f"  Phaseable (>{100*PHASEABLE_THRESHOLD:.0f}%):     {phaseable_ps} PS, "
        f"{phaseable_conc+phaseable_disc:,} reads, {100*phaseable_acc:.2f}%",
        f"  Unphaseable (≤{100*PHASEABLE_THRESHOLD:.0f}%):    {unphaseable_ps} PS, "
        f"{unphaseable_conc+unphaseable_disc:,} reads, {100*unphaseable_acc:.2f}%",
    ]
    lines += [
        "",
        "  Per-chromosome:",
    ]
    for ch in sorted(cs):
        s = cs[ch]; acc = s["conc"]/s["tot"] if s["tot"] else 0
        lines.append(f"    {ch:20s} {s['n_ps']:5d} PS  "
                     f"{s['tot']:7d} reads  {100*acc:.2f}%")
    lines += ["", "=" * 64]
    for l in lines:
        fh.write(l + "\n"); print(l, file=sys.stderr)

# ── IGV-ready files for worst phase sets ─────────────────────────────
N_IGV_BLOCKS = 20
igv_dir = os.path.join(outdir, "igv")
os.makedirs(igv_dir, exist_ok=True)

# Group reads by PS for quick lookup
ps_reads = collections.defaultdict(list)
for qname in read_truth:
    hp, ps = read_phase[qname]
    t = read_truth[qname]
    orient = ps_orient.get(ps)
    if orient is None:
        continue
    conc = is_concordant(hp, t.hap, orient)
    ps_reads[ps].append((qname, hp, t.hap, t.chrom, t.pos,
                          "concordant" if conc else "DISCORDANT"))

igv_blocks = sorted(results, key=lambda r: r["accuracy"])[:N_IGV_BLOCKS]

block_beds = []
for i, blk in enumerate(igv_blocks):
    ps = blk["phase_set"]
    reads = ps_reads.get(ps, [])
    if not reads:
        continue

    tag = f"ps{ps}_acc{blk['accuracy']:.2f}"
    names_file = os.path.join(igv_dir, f"{tag}.names.txt")
    with open(names_file, "w") as fh:
        for r in reads:
            fh.write(r[0] + "\n")

    bed_file = os.path.join(igv_dir, f"{tag}.bed")
    # Group positions by contig to avoid mixing coordinates
    contig_pos = collections.defaultdict(list)
    chrom_set = set()
    with open(bed_file, "w") as fh:
        for qname, hp, thap, tchrom, tpos, status in reads:
            color = "255,0,0" if status == "DISCORDANT" else "0,0,255"
            label = f"HP{hp}_{thap}_{status}"
            fh.write(f"{tchrom}\t{tpos}\t{tpos+1}\t{label}:{qname}\t0\t.\t"
                     f"{tpos}\t{tpos+1}\t{color}\n")
            chrom_set.add(tchrom)
            contig_pos[tchrom].append(tpos)

    # Pick the contig with the most reads for the IGV region
    best_contig = max(contig_pos, key=lambda c: len(contig_pos[c]))
    min_pos = min(contig_pos[best_contig])
    max_pos = max(contig_pos[best_contig])

    block_beds.append({
        "tag": tag, "ps": ps, "names_file": names_file, "bed_file": bed_file,
        "chroms": sorted(chrom_set), "best_contig": best_contig,
        "min_pos": min_pos, "max_pos": max_pos,
        "accuracy": blk["accuracy"], "n_reads": blk["n_reads"],
        "discordant": blk["discordant"],
    })

# IGV batch script
igv_batch = os.path.join(igv_dir, "load_worst_blocks.bat")
with open(igv_batch, "w") as fh:
    fh.write("new\n")
    fh.write(f"load {os.path.abspath(truth_bam)}\n")
    fh.write("colorBy YC\n")
    fh.write("group TAG HP\n")
    fh.write("sort BASE\n")
    fh.write("squish\n")
    fh.write(f"snapshotDirectory {os.path.abspath(igv_dir)}\n\n")
    for b in block_beds:
        pad = max(500, (b["max_pos"] - b["min_pos"]) // 10)
        region = f"{b['best_contig']}:{max(0,b['min_pos']-pad)}-{b['max_pos']+pad}"
        fh.write(f"goto {region}\n")
        fh.write(f"snapshot {b['tag']}.png\n\n")

# Extraction script (uses tagged BAM so HP/PS/YC tags are preserved)
extract_script = os.path.join(igv_dir, "extract_block_bams.sh")
with open(extract_script, "w") as fh:
    fh.write("#!/usr/bin/env bash\n")
    fh.write("# Extract per-block BAMs from the tagged truth alignment.\n")
    fh.write(f"# Source: {os.path.abspath(truth_bam)}  (has HP/PS/YC/HO/hq tags)\n\n")
    fh.write("set -euo pipefail\n\n")
    fh.write(f'SAMTOOLS="${{SAMTOOLS:-{samtools}}}"\n\n')
    for b in block_beds:
        out_bam_path = os.path.join(igv_dir, f"{b['tag']}.bam")
        fh.write(f"# PS={b['ps']}  accuracy={b['accuracy']:.2f}  "
                 f"reads={b['n_reads']}  discordant={b['discordant']}\n")
        fh.write(f'"$SAMTOOLS" view -b -N {os.path.abspath(b["names_file"])} '
                 f'{os.path.abspath(truth_bam)} > {os.path.abspath(out_bam_path)}\n')
        fh.write(f'"$SAMTOOLS" index {os.path.abspath(out_bam_path)}\n\n')
    fh.write('echo "Done. Load BAMs in IGV or run: igv -b load_worst_blocks.bat"\n')
os.chmod(extract_script, 0o755)

# Block summary README
block_summary = os.path.join(igv_dir, "README.txt")
with open(block_summary, "w") as fh:
    fh.write("IGV inspection files for worst phase sets\n")
    fh.write("=" * 50 + "\n\n")
    fh.write("The tagged truth BAM (diplinator_merged.tagged.bam) has HP/PS/YC\n")
    fh.write("tags injected from the phased uBAM, plus HO/hq from diplinator.\n")
    fh.write("IGV will automatically:\n")
    fh.write("  - Color reads blue (HP=1) / red (HP=2) via the YC tag\n")
    fh.write("  - Group reads by haplotype via the HP tag\n\n")
    fh.write("Files per block:\n")
    fh.write("  <tag>.names.txt   — read names for samtools view -N\n")
    fh.write("  <tag>.bed         — read positions colored red (discordant) / blue (concordant)\n")
    fh.write("  <tag>.bam         — extracted reads with HP/PS/YC/HO/hq tags (after extract)\n\n")
    fh.write("Workflow:\n")
    fh.write("  1. bash extract_block_bams.sh\n")
    fh.write("  2. Open IGV, load the maternal or paternal FASTA\n")
    fh.write("  3. Load a <tag>.bam — reads auto-colored by haplotype\n")
    fh.write("  4. Right-click track > Group by > tag > HP\n")
    fh.write("  5. Or batch: igv -b load_worst_blocks.bat\n\n")
    fh.write(f"{'PS':>12s}  {'Acc':>6s}  {'Reads':>6s}  {'Disc':>5s}  Region\n")
    fh.write("-" * 60 + "\n")
    for b in block_beds:
        region = f"{b['best_contig']}:{b['min_pos']}-{b['max_pos']}"
        fh.write(f"{b['ps']:>12d}  {b['accuracy']:>6.2f}  {b['n_reads']:>6d}  "
                 f"{b['discordant']:>5d}  {region}\n")

n_igv = len(block_beds)
print(f"\n  IGV files for {n_igv} worst blocks in {igv_dir}/", file=sys.stderr)

print(f"\n  Outputs: {summary_file}", file=sys.stderr)
print(f"           {json_file}", file=sys.stderr)
print(f"           {per_ps_file}", file=sys.stderr)
print(f"           {per_chrom_file}", file=sys.stderr)
print(f"           {per_read_file}", file=sys.stderr)
print(f"           {worst_file}", file=sys.stderr)
print(f"           {bad_bed_file}", file=sys.stderr)
print(f"           {igv_dir}/  ({n_igv} blocks)", file=sys.stderr)
