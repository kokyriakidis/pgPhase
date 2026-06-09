#!/usr/bin/env bash
#
# prepare_refseq_gene_bed.sh — Download RefSeq GRCh38 GFF3 and convert
# to a gene BED file for phasing completeness evaluation.
#
# Replicates the HiPhase paper (Holt et al. 2024) gene annotation:
#   - NCBI Homo sapiens Annotation Release 110
#   - Feature types: gene, pseudogene
#   - Sources: BestRefSeq, RefSeq, Gnomon, Curated Genomic, BestRefSeq%2CGnomon
#   - Primary chromosomes only (chr1-22, chrX, chrY, chrM)
#
# Output: 4-column BED (chrom, start, end, gene_name)
#
# Usage:
#   ./prepare_refseq_gene_bed.sh [output.bed]
#   Default output: refseq_genes_grch38_ar110.bed

set -euo pipefail

OUTBED="${1:-refseq_genes_grch38_ar110.bed}"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
GFF_GZ="$(mktemp --suffix=.gff.gz)"

trap 'rm -f "$GFF_GZ"' EXIT

echo "Downloading RefSeq GRCh38 Annotation Release 110 GFF3 ..."
curl -sL "$GFF_URL" -o "$GFF_GZ"

echo "Converting to BED ..."
python3 - "$GFF_GZ" "$OUTBED" <<'PYEOF'
import gzip, sys, re

gff_path = sys.argv[1]
out_path = sys.argv[2]

# NCBI RefSeq accession -> UCSC chromosome name
ACCESSION_TO_CHROM = {}
for i in range(1, 23):
    ACCESSION_TO_CHROM[f"NC_000{i:03d}"] = f"chr{i}"
ACCESSION_TO_CHROM["NC_000023"] = "chrX"
ACCESSION_TO_CHROM["NC_000024"] = "chrY"
ACCESSION_TO_CHROM["NC_012920"] = "chrM"

ALLOWED_TYPES = {"gene", "pseudogene"}
ALLOWED_SOURCES = {"BestRefSeq", "RefSeq", "Gnomon", "Curated Genomic",
                   "BestRefSeq%2CGnomon", "BestRefSeq,Gnomon"}

def parse_name(attributes):
    """Extract gene Name from GFF3 attributes."""
    for attr in attributes.split(";"):
        if attr.startswith("Name="):
            return attr[5:]
    # Fallback to gene ID
    for attr in attributes.split(";"):
        if attr.startswith("ID="):
            return attr[3:]
    return "unknown"

def accession_to_chrom(seqid):
    """Convert NC_000001.11 -> chr1 (strip version suffix)."""
    base = seqid.split(".")[0]
    return ACCESSION_TO_CHROM.get(base)

count = 0
with gzip.open(gff_path, "rt") as fh, open(out_path, "w") as out:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        seqid, source, ftype = parts[0], parts[1], parts[2]
        start, end = int(parts[3]), int(parts[4])
        attributes = parts[8]

        if ftype not in ALLOWED_TYPES:
            continue
        if source not in ALLOWED_SOURCES:
            continue

        chrom = accession_to_chrom(seqid)
        if chrom is None:
            continue

        name = parse_name(attributes)
        # GFF3 is 1-based inclusive; BED is 0-based half-open
        out.write(f"{chrom}\t{start - 1}\t{end}\t{name}\n")
        count += 1

print(f"  Wrote {count} gene/pseudogene entries to {out_path}", file=sys.stderr)
PYEOF

echo "Done: $OUTBED"
wc -l "$OUTBED"
