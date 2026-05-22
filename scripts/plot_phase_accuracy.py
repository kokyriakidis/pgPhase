#!/usr/bin/env python3
"""Generate plots from evaluate_phase_accuracy.py outputs."""
import sys, os, csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

outdir = sys.argv[1]

# ── load data ────────────────────────────────────────────────────────
with open(os.path.join(outdir, "summary.json")) as f:
    summary = json.load(f)

ps_rows = []
with open(os.path.join(outdir, "per_phase_set.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        r["accuracy"] = float(r["accuracy"])
        r["n_reads"] = int(r["n_reads"])
        r["span_bp"] = int(r["span_bp"])
        r["concordant"] = int(r["concordant"])
        r["discordant"] = int(r["discordant"])
        ps_rows.append(r)

chrom_rows = []
with open(os.path.join(outdir, "per_chromosome.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        chrom_rows.append(r)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(f"Phase Accuracy — {summary['overall_accuracy']*100:.2f}% overall, "
             f"N50={summary['phase_block_n50_bp']:,} bp",
             fontsize=13, fontweight="bold")

# ── plot 1: accuracy distribution ────────────────────────────────────
ax = axes[0, 0]
accs = [r["accuracy"] for r in ps_rows]
ax.hist(accs, bins=50, range=(0, 1), color="#4C72B0", edgecolor="white", linewidth=0.5)
ax.set_xlabel("Phase set accuracy")
ax.set_ylabel("Count")
ax.set_title("Accuracy distribution")
ax.axvline(x=summary["overall_accuracy"], color="red", linestyle="--", linewidth=1,
           label=f"mean={summary['overall_accuracy']*100:.1f}%")
ax.legend(fontsize=9)

# ── plot 2: accuracy vs block span ───────────────────────────────────
ax = axes[0, 1]
spans = [r["span_bp"] for r in ps_rows]
ax.scatter(spans, accs, s=8, alpha=0.4, c="#4C72B0", edgecolors="none")
ax.set_xlabel("Phase block span (bp)")
ax.set_ylabel("Accuracy")
ax.set_title("Accuracy vs block span")
ax.set_xscale("symlog", linthresh=1000)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6
                                                   else f"{x/1e3:.0f}k" if x >= 1e3
                                                   else f"{x:.0f}"))
ax.set_ylim(-0.05, 1.05)

# ── plot 3: per-chromosome accuracy ──────────────────────────────────
ax = axes[1, 0]
chroms = [r["chrom"] for r in chrom_rows]
chr_accs = [float(r["accuracy"]) for r in chrom_rows]
chr_reads = [int(r["reads"]) for r in chrom_rows]
bars = ax.bar(range(len(chroms)), [a * 100 for a in chr_accs], color="#4C72B0", edgecolor="white")
ax.set_xticks(range(len(chroms)))
ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Accuracy (%)")
ax.set_title("Per-chromosome accuracy")
ax.set_ylim(max(0, min(a * 100 for a in chr_accs) - 5) if chr_accs else 0, 105)
ax.axhline(y=summary["overall_accuracy"] * 100, color="red", linestyle="--", linewidth=1)

# ── plot 4: phase block span distribution ────────────────────────────
ax = axes[1, 1]
spans_kb = [s / 1000 for s in spans if s > 0]
if spans_kb:
    ax.hist(spans_kb, bins=50, color="#55A868", edgecolor="white", linewidth=0.5)
    ax.axvline(x=summary["phase_block_n50_bp"] / 1000, color="red", linestyle="--",
               linewidth=1, label=f"N50={summary['phase_block_n50_bp']/1e3:.0f} kb")
    ax.legend(fontsize=9)
ax.set_xlabel("Phase block span (kb)")
ax.set_ylabel("Count")
ax.set_title("Block span distribution")
ax.set_xscale("symlog", linthresh=10)

plt.tight_layout()
plt.savefig(os.path.join(outdir, "phase_accuracy.png"), dpi=150, bbox_inches="tight")
plt.savefig(os.path.join(outdir, "phase_accuracy.pdf"), bbox_inches="tight")
print(f"  Saved phase_accuracy.png and .pdf", file=sys.stderr)
