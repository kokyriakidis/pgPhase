#!/usr/bin/env bash
#
# Build a merged BED of difficult regions for HG002 v1.1 T2T assembly.
#
# Combines centromere/satellite annotations and segmental duplications
# from the T2T consortium into a single sorted, merged BED file.
#
# Sources:
#   - Centromeric satellites: hg002v1.1.cenSatv2.0.bed
#   - Segmental duplications: hg002v1.1.SDs.013025.bed
#
# These regions are expected to have unreliable phasing due to
# repetitive sequence, low heterozygosity, or ambiguous read mapping.
#
# Usage:
#   ./build_difficult_regions_bed.sh [--outdir DIR]
#
# Requires: bedtools, curl (for download)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${1:-$SCRIPT_DIR/../data/annotations}"
mkdir -p "$OUTDIR"

CENSAT_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.bed"
SEGDUP_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/segdups/hg002v1.1.SDs.013025.bed"

CENSAT="$OUTDIR/hg002v1.1.cenSat.bed"
SEGDUP="$OUTDIR/hg002v1.1.segdups.bed"
OUTPUT="$OUTDIR/hg002v1.1.difficult_regions.bed"

# Download if not present
if [[ ! -f "$CENSAT" ]]; then
    echo "Downloading centromere/satellite annotations ..."
    curl -sL "$CENSAT_URL" -o "$CENSAT"
fi

if [[ ! -f "$SEGDUP" ]]; then
    echo "Downloading segmental duplications ..."
    curl -sL "$SEGDUP_URL" -o "$SEGDUP"
fi

echo "Building merged difficult regions BED ..."

# Extract chrom/start/end, skip track headers, add region type, sort, merge
{
    # Centromere/satellite regions
    grep -v '^track' "$CENSAT" | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "censat"}'

    # Segmental duplications
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "segdup"}' "$SEGDUP"
} | sort -k1,1 -k2,2n | awk '
BEGIN { OFS="\t" }
{
    if ($1 == chr && $2 <= end) {
        # Overlapping or adjacent — extend
        if ($3 > end) end = $3
        if (index(types, $4) == 0) types = types "," $4
    } else {
        if (chr != "") print chr, start, end, types
        chr = $1; start = $2; end = $3; types = $4
    }
}
END { if (chr != "") print chr, start, end, types }
' > "$OUTPUT"

N_REGIONS=$(wc -l < "$OUTPUT")
TOTAL_BP=$(awk '{sum += $3 - $2} END {print sum}' "$OUTPUT")

echo "  $N_REGIONS merged regions, ${TOTAL_BP} bp total"
echo "  Written to: $OUTPUT"

# Also create a per-chromosome summary
echo ""
echo "  Per-chromosome coverage:"
awk '{
    chr = $1
    sub(/_MATERNAL$/, "", chr)
    sub(/_PATERNAL$/, "", chr)
    bp[chr] += $3 - $2
    n[chr]++
}
END {
    for (c in bp) printf "    %-25s %6d regions  %12d bp\n", c, n[c], bp[c]
}' "$OUTPUT" | sort -k1,1V
