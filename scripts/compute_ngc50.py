#!/usr/bin/env python3
"""Compute NGC50 (corrected NG50) by subtracting switch error regions
from phase blocks.

Replicates the HiPhase paper methodology:
  1. Convert whatshap stats --block-list TSV to BED
  2. bedtools subtract switch error BED from phase block BED
  3. Compute NG50 on the corrected blocks (genome size = 3.1 Gb)

Usage:
  python3 compute_ngc50.py <phase_blocks.tsv> <switch_errors.bed> <output.json>
"""
import csv, json, subprocess, sys, tempfile, os

def main():
    if len(sys.argv) < 4:
        print("Usage: compute_ngc50.py <blocks.tsv> <errors.bed> <output.json>",
              file=sys.stderr)
        sys.exit(1)

    blocks_tsv = sys.argv[1]
    error_bed  = sys.argv[2]
    output_json = sys.argv[3]

    # Convert whatshap phase blocks TSV to BED
    # whatshap stats --block-list columns: sample, chromosome, phase_set, from, to, variants
    block_bed_lines = []
    with open(blocks_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("chromosome", row.get("#chromosome", ""))
            start = row.get("from", row.get("phase_set", ""))
            end   = row.get("to", "")
            if chrom and start and end:
                block_bed_lines.append(f"{chrom}\t{start}\t{end}")

    if not block_bed_lines:
        json.dump({"ngc50_bp": 0, "corrected_blocks": 0}, open(output_json, "w"))
        return

    # Write phase block BED
    block_bed = tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False)
    block_bed.write("\n".join(block_bed_lines) + "\n")
    block_bed.close()

    corrected_bed = tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False)
    corrected_bed.close()

    try:
        # bedtools subtract: remove switch error regions from phase blocks
        with open(corrected_bed.name, "w") as out:
            subprocess.run(
                ["bedtools", "subtract", "-a", block_bed.name, "-b", error_bed],
                stdout=out, check=True
            )

        # Compute NG50 from corrected blocks (genome size = 3.1 Gb)
        spans = []
        with open(corrected_bed.name) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    spans.append(int(parts[2]) - int(parts[1]))

        spans.sort(reverse=True)
        genome_size = 3_100_000_000
        half = genome_size / 2
        running = 0
        ngc50 = 0
        for s in spans:
            running += s
            if running >= half:
                ngc50 = s
                break

        result = {
            "ngc50_bp": ngc50,
            "corrected_blocks": len(spans),
            "total_corrected_span_bp": sum(spans),
        }
        with open(output_json, "w") as f:
            json.dump(result, f, indent=2)
            f.write("\n")

        print(f"  NGC50: {ngc50:,} bp ({len(spans)} corrected blocks)",
              file=sys.stderr)

    finally:
        os.unlink(block_bed.name)
        os.unlink(corrected_bed.name)


if __name__ == "__main__":
    main()
