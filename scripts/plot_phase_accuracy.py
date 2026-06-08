#!/usr/bin/env python3
"""Generate multi-page plots from evaluate_phase_accuracy.py outputs.

Page 1: Summary table (the "money slide")
Page 2: Contiguity (Nx curve + block span distribution)
Page 3: Accuracy deep dive (histogram, switch vs hamming, per-category, accuracy vs site density)
Page 4: Genomic distribution (per-chromosome accuracy + block count, accuracy vs span)
"""
import sys, os, csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

outdir = sys.argv[1]

# ── load data ────────────────────────────────────────────────────────
with open(os.path.join(outdir, "summary.json")) as f:
    summary = json.load(f)

ps_rows = []
with open(os.path.join(outdir, "per_phase_set.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        for k in ("accuracy", "hamming_rate", "switch_rate", "site_density"):
            if k in r and r[k]:
                r[k] = float(r[k])
            else:
                r[k] = 0.0
        for k in ("n_reads", "span_bp", "concordant", "discordant", "switches",
                   "n_sites"):
            if k in r and r[k]:
                r[k] = int(r[k])
            else:
                r[k] = 0
        ps_rows.append(r)

chrom_rows = []
with open(os.path.join(outdir, "per_chromosome.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        chrom_rows.append(r)

# Load Nx curve if available
nx_data = []
nx_file = os.path.join(outdir, "nx_curve.csv")
if os.path.exists(nx_file):
    with open(nx_file) as f:
        next(f)  # skip header
        for line in f:
            pct, sp = line.strip().split(",")
            nx_data.append((float(pct), int(sp)))

# Load per-category if available
cat_rows = []
cat_file = os.path.join(outdir, "per_category.tsv")
if os.path.exists(cat_file):
    with open(cat_file) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cat_rows.append(r)

# ── color palette ────────────────────────────────────────────────────
C_BLUE = "#4C72B0"
C_GREEN = "#55A868"
C_RED = "#C44E52"
C_ORANGE = "#DD8452"
C_PURPLE = "#8172B3"
C_GREY = "#999999"

# ════════════════════════════════════════════════════════════════════
# Page 1: Summary table
# ════════════════════════════════════════════════════════════════════
fig1, ax = plt.subplots(figsize=(10, 7))
ax.axis("off")
ax.set_title("pgphase — Phase Accuracy Summary", fontsize=16, fontweight="bold",
             pad=20)

def fmt_bp(bp):
    if bp >= 1e9: return f"{bp/1e9:.2f} Gb"
    if bp >= 1e6: return f"{bp/1e6:.2f} Mb"
    if bp >= 1e3: return f"{bp/1e3:.1f} kb"
    return f"{bp:,} bp"

def fmt_pct(v): return f"{100*v:.2f}%"

rows = [
    ["Metric", "Value"],
    ["", ""],
    ["COMPLETENESS", ""],
    ["Total input reads", f"{summary.get('total_input_reads', 0):,}"],
    ["Phased reads", f"{summary.get('total_phased_reads', 0):,} ({fmt_pct(summary.get('fraction_phased', 0))})"],
    ["Total phase sets", f"{summary.get('total_phase_sets', 0):,}"],
    ["", ""],
    ["CONTIGUITY", ""],
    ["Phase block N50", fmt_bp(summary.get("phase_block_n50_bp", 0))],
    ["Phase block NG50", fmt_bp(summary.get("phase_block_ng50_bp", 0))],
    ["Phase block auN", fmt_bp(summary.get("phase_block_aun_bp", 0))],
    ["Largest block", f"{fmt_bp(summary.get('phase_block_max_span_bp', 0))} ({summary.get('phase_block_max_span_chrom', '')})"],
    ["Genome covered", f"{fmt_bp(summary.get('genome_covered_bp', 0))} ({fmt_pct(summary.get('fraction_genome_covered', 0))})"],
    ["Mean blocks/chrom", f"{summary.get('mean_blocks_per_chrom', 0):.1f}"],
    ["", ""],
    ["ACCURACY", ""],
    ["Phase sets evaluated", f"{summary.get('phase_sets_evaluated', 0):,}"],
    ["Perfect phase sets", f"{summary.get('phase_sets_perfect', 0):,} ({summary.get('phase_sets_perfect_pct', 0):.1f}%)"],
    ["Overall accuracy", fmt_pct(summary.get("overall_accuracy", 0))],
    ["Hamming error rate", fmt_pct(summary.get("hamming_error_rate", 0))],
    ["Switch error rate", fmt_pct(summary.get("switch_error_rate", 0))],
]

# Render as table
cell_text = [r for r in rows]
table = ax.table(cellText=cell_text, cellLoc="left", loc="center",
                 colWidths=[0.45, 0.45])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 1.4)

# Style: bold section headers, hide grid for blank rows
for i, row in enumerate(rows):
    for j in range(2):
        cell = table[i, j]
        cell.set_edgecolor("#DDDDDD")
        cell.set_linewidth(0.5)
        if row[0] in ("Metric", ""):
            cell.set_facecolor("white")
            cell.set_edgecolor("white")
        elif row[0] in ("COMPLETENESS", "CONTIGUITY", "ACCURACY"):
            cell.set_facecolor("#E8EDF2")
            cell.get_text().set_fontweight("bold")
            cell.get_text().set_fontsize(10)
        elif row[0] == "Metric":
            cell.set_facecolor("#4C72B0")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")

fig1.tight_layout()
fig1.savefig(os.path.join(outdir, "page1_summary.png"), dpi=150, bbox_inches="tight")
fig1.savefig(os.path.join(outdir, "page1_summary.pdf"), bbox_inches="tight")
plt.close(fig1)

# ════════════════════════════════════════════════════════════════════
# Page 2: Contiguity
# ════════════════════════════════════════════════════════════════════
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
fig2.suptitle("Phase Block Contiguity", fontsize=14, fontweight="bold")

# Plot 2a: Nx curve
ax = axes2[0]
if nx_data:
    pcts = [d[0] for d in nx_data]
    spans_mb = [d[1] / 1e6 for d in nx_data]
    ax.plot(pcts, spans_mb, color=C_BLUE, linewidth=2)
    ax.fill_between(pcts, spans_mb, alpha=0.15, color=C_BLUE)
    # Mark N50
    n50_mb = summary.get("phase_block_n50_bp", 0) / 1e6
    ax.axhline(y=n50_mb, color=C_RED, linestyle="--", linewidth=1,
               label=f"N50 = {fmt_bp(summary.get('phase_block_n50_bp', 0))}")
    ax.axvline(x=50, color=C_GREY, linestyle=":", linewidth=0.8)
    ax.legend(fontsize=9)
ax.set_xlabel("Cumulative % of phased bases")
ax.set_ylabel("Block span (Mb)")
ax.set_title("Nx Curve")
ax.set_xlim(0, 100)
ax.set_ylim(bottom=0)

# Plot 2b: Block span distribution
ax = axes2[1]
spans_kb = [r["span_bp"] / 1000 for r in ps_rows if r["span_bp"] > 0]
if spans_kb:
    ax.hist(spans_kb, bins=50, color=C_GREEN, edgecolor="white", linewidth=0.5)
    n50_kb = summary.get("phase_block_n50_bp", 0) / 1000
    ax.axvline(x=n50_kb, color=C_RED, linestyle="--", linewidth=1,
               label=f"N50 = {n50_kb:.0f} kb")
    ax.legend(fontsize=9)
ax.set_xlabel("Phase block span (kb)")
ax.set_ylabel("Count")
ax.set_title("Block Span Distribution")
ax.set_xscale("symlog", linthresh=10)

fig2.tight_layout()
fig2.savefig(os.path.join(outdir, "page2_contiguity.png"), dpi=150, bbox_inches="tight")
fig2.savefig(os.path.join(outdir, "page2_contiguity.pdf"), bbox_inches="tight")
plt.close(fig2)

# ════════════════════════════════════════════════════════════════════
# Page 3: Accuracy deep dive (2x2)
# ════════════════════════════════════════════════════════════════════
fig3, axes3 = plt.subplots(2, 2, figsize=(14, 10))
fig3.suptitle("Phase Accuracy Analysis", fontsize=14, fontweight="bold")

# Plot 3a: Accuracy histogram
ax = axes3[0, 0]
accs = [r["accuracy"] for r in ps_rows]
ax.hist(accs, bins=50, range=(0, 1), color=C_BLUE, edgecolor="white", linewidth=0.5)
overall = summary.get("overall_accuracy", 0)
ax.axvline(x=overall, color=C_RED, linestyle="--", linewidth=1,
           label=f"overall = {100*overall:.1f}%")
ax.set_xlabel("Phase set accuracy")
ax.set_ylabel("Count")
ax.set_title("Accuracy Distribution")
ax.legend(fontsize=9)

# Plot 3b: Switch error rate vs Hamming error rate (per phase set)
ax = axes3[0, 1]
hamming = [r["hamming_rate"] for r in ps_rows]
switch = [r["switch_rate"] for r in ps_rows]
sizes = [max(3, min(50, r["n_reads"] / 5)) for r in ps_rows]
ax.scatter(hamming, switch, s=sizes, alpha=0.4, c=C_BLUE, edgecolors="none")
ax.plot([0, 1], [0, 1], color=C_GREY, linestyle=":", linewidth=0.8, alpha=0.5)
ax.set_xlabel("Hamming error rate")
ax.set_ylabel("Switch error rate")
ax.set_title("Switch vs Hamming Error (per phase set)")
ax.set_xlim(-0.02, max(0.5, max(hamming) * 1.1) if hamming else 0.5)
ax.set_ylim(-0.02, max(0.5, max(switch) * 1.1) if switch else 0.5)

# Plot 3c: Per-category accuracy (grouped bar chart)
ax = axes3[1, 0]
if cat_rows:
    cats = [r["category"] for r in cat_rows]
    cat_acc = [float(r["accuracy"]) * 100 for r in cat_rows]
    cat_serr = [float(r["switch_error"]) * 100 for r in cat_rows]
    x = np.arange(len(cats))
    w = 0.35
    bars1 = ax.bar(x - w/2, cat_acc, w, label="Accuracy", color=C_BLUE)
    bars2 = ax.bar(x + w/2, cat_serr, w, label="Switch error", color=C_ORANGE)
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=9)
    ax.set_ylabel("Rate (%)")
    ax.set_title("Per-Category Accuracy & Switch Error")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 105)
    # Add read count labels on top of accuracy bars
    for bar, row in zip(bars1, cat_rows):
        n = int(row["reads"])
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"n={n:,}", ha="center", va="bottom", fontsize=7, color=C_GREY)
else:
    ax.text(0.5, 0.5, "No category data\n(provide --censat-bed / --segdup-bed)",
            ha="center", va="center", transform=ax.transAxes, fontsize=10, color=C_GREY)
    ax.set_title("Per-Category Accuracy")

# Plot 3d: Accuracy vs site density
ax = axes3[1, 1]
has_density = any(r.get("site_density", 0) > 0 for r in ps_rows)
if has_density:
    densities = [r["site_density"] for r in ps_rows]
    accs_d = [r["accuracy"] for r in ps_rows]
    sizes_d = [max(3, min(50, r["n_reads"] / 5)) for r in ps_rows]
    ax.scatter(densities, accs_d, s=sizes_d, alpha=0.4, c=C_GREEN, edgecolors="none")
    ax.set_xlabel("Het site density (sites/kb)")
    ax.set_ylabel("Accuracy")
    ax.set_title("Accuracy vs Site Density")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xscale("symlog", linthresh=1)
else:
    ax.text(0.5, 0.5, "No site density data\n(provide --sites-vcf)",
            ha="center", va="center", transform=ax.transAxes, fontsize=10, color=C_GREY)
    ax.set_title("Accuracy vs Site Density")

fig3.tight_layout()
fig3.savefig(os.path.join(outdir, "page3_accuracy.png"), dpi=150, bbox_inches="tight")
fig3.savefig(os.path.join(outdir, "page3_accuracy.pdf"), bbox_inches="tight")
plt.close(fig3)

# ════════════════════════════════════════════════════════════════════
# Page 4: Genomic distribution
# ════════════════════════════════════════════════════════════════════
fig4, axes4 = plt.subplots(2, 1, figsize=(14, 9))
fig4.suptitle("Genomic Distribution", fontsize=14, fontweight="bold")

# Plot 4a: Per-chromosome accuracy + block count (dual axis)
ax = axes4[0]
chroms = [r["chrom"] for r in chrom_rows]
chr_accs = [float(r["accuracy"]) * 100 for r in chrom_rows]
chr_ps = [int(r["phase_sets"]) for r in chrom_rows]
x = np.arange(len(chroms))
bars = ax.bar(x, chr_accs, color=C_BLUE, edgecolor="white", alpha=0.8, label="Accuracy")
ax.set_xticks(x)
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Accuracy (%)", color=C_BLUE)
ax.set_ylim(max(0, min(chr_accs) - 5) if chr_accs else 0, 105)
ax.axhline(y=overall * 100, color=C_RED, linestyle="--", linewidth=1, alpha=0.7)

ax2 = ax.twinx()
ax2.plot(x, chr_ps, color=C_ORANGE, marker="o", markersize=4, linewidth=1.5,
         label="# blocks")
ax2.set_ylabel("Phase blocks", color=C_ORANGE)
ax2.tick_params(axis="y", labelcolor=C_ORANGE)

# Combined legend
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc="lower left")
ax.set_title("Per-Chromosome Accuracy & Block Count")

# Plot 4b: Accuracy vs block span with marginal histograms effect
ax = axes4[1]
spans = [r["span_bp"] for r in ps_rows]
accs = [r["accuracy"] for r in ps_rows]
n_reads = [r["n_reads"] for r in ps_rows]
sizes = [max(3, min(80, n / 3)) for n in n_reads]
colors = [C_RED if a < 0.6 else C_ORANGE if a < 0.9 else C_BLUE for a in accs]
ax.scatter(spans, accs, s=sizes, alpha=0.4, c=colors, edgecolors="none")
ax.set_xlabel("Phase block span (bp)")
ax.set_ylabel("Accuracy")
ax.set_title("Accuracy vs Block Span (size = read count, red = <60%, orange = <90%)")
ax.set_xscale("symlog", linthresh=1000)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(
    lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6
    else f"{x/1e3:.0f}k" if x >= 1e3 else f"{x:.0f}"))
ax.set_ylim(-0.05, 1.05)
ax.axhline(y=0.6, color=C_RED, linestyle=":", linewidth=0.8, alpha=0.5)
ax.axhline(y=0.9, color=C_ORANGE, linestyle=":", linewidth=0.8, alpha=0.5)

fig4.tight_layout()
fig4.savefig(os.path.join(outdir, "page4_genomic.png"), dpi=150, bbox_inches="tight")
fig4.savefig(os.path.join(outdir, "page4_genomic.pdf"), bbox_inches="tight")
plt.close(fig4)

# Also keep the legacy single-page plot for backward compatibility
fig_legacy, axes_l = plt.subplots(2, 2, figsize=(14, 10))
fig_legacy.suptitle(f"Phase Accuracy — {overall*100:.2f}% overall, "
                    f"N50={summary.get('phase_block_n50_bp',0):,} bp",
                    fontsize=13, fontweight="bold")

ax = axes_l[0, 0]
ax.hist([r["accuracy"] for r in ps_rows], bins=50, range=(0, 1),
        color=C_BLUE, edgecolor="white", linewidth=0.5)
ax.axvline(x=overall, color="red", linestyle="--", linewidth=1,
           label=f"mean={overall*100:.1f}%")
ax.set_xlabel("Phase set accuracy"); ax.set_ylabel("Count")
ax.set_title("Accuracy distribution"); ax.legend(fontsize=9)

ax = axes_l[0, 1]
ax.scatter([r["span_bp"] for r in ps_rows], [r["accuracy"] for r in ps_rows],
           s=8, alpha=0.4, c=C_BLUE, edgecolors="none")
ax.set_xlabel("Phase block span (bp)"); ax.set_ylabel("Accuracy")
ax.set_title("Accuracy vs block span")
ax.set_xscale("symlog", linthresh=1000)
ax.set_ylim(-0.05, 1.05)

ax = axes_l[1, 0]
ax.bar(range(len(chroms)), [float(r["accuracy"])*100 for r in chrom_rows],
       color=C_BLUE, edgecolor="white")
ax.set_xticks(range(len(chroms)))
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Accuracy (%)"); ax.set_title("Per-chromosome accuracy")
ax.axhline(y=overall*100, color="red", linestyle="--", linewidth=1)

ax = axes_l[1, 1]
spans_kb_l = [r["span_bp"]/1000 for r in ps_rows if r["span_bp"] > 0]
if spans_kb_l:
    ax.hist(spans_kb_l, bins=50, color=C_GREEN, edgecolor="white", linewidth=0.5)
    ax.axvline(x=summary.get("phase_block_n50_bp",0)/1000, color="red",
               linestyle="--", linewidth=1)
ax.set_xlabel("Phase block span (kb)"); ax.set_ylabel("Count")
ax.set_title("Block span distribution"); ax.set_xscale("symlog", linthresh=10)

fig_legacy.tight_layout()
fig_legacy.savefig(os.path.join(outdir, "phase_accuracy.png"), dpi=150, bbox_inches="tight")
fig_legacy.savefig(os.path.join(outdir, "phase_accuracy.pdf"), bbox_inches="tight")
plt.close(fig_legacy)

print(f"  Saved plots:", file=sys.stderr)
for f in ["page1_summary.png", "page2_contiguity.png",
          "page3_accuracy.png", "page4_genomic.png",
          "phase_accuracy.png"]:
    p = os.path.join(outdir, f)
    if os.path.exists(p):
        print(f"    {p}", file=sys.stderr)
