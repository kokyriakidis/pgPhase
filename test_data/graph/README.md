# Graph pipeline test data

Small pangenome graph test case for `collect-graph-variation`, covering the
MICB (chr6) and KIR3DL1 (chr19) loci from the HPRC pangenome.  Based on the
haplotype sampling test data from
[gbz-base](https://github.com/jltsiren/gbz-base/tree/main/test-data).

## Running the graph pipeline

Minimal run (no `--region` needed — the FFI clamps queries to the GBZ path range):

```bash
./pgphase collect-graph-variation \
  --ref test_data/graph/ref.fa.gz \
  --sites test_data/graph/micb-kir3dl1.sites.vcf.gz \
  --gbz-db test_data/graph/micb-kir3dl1.gbz.db \
  --gaf-db test_data/graph/micb-kir3dl1_HG003.gaf.db \
  --sample CHM13 \
  -o output.tsv
```

Expected output: **292 candidate variant sites** (53 on chr6, 239 on chr19).

To restrict to a single locus:

```bash
./pgphase collect-graph-variation \
  --ref test_data/graph/ref.fa.gz \
  --sites test_data/graph/micb-kir3dl1.sites.vcf.gz \
  --gbz-db test_data/graph/micb-kir3dl1.gbz.db \
  --gaf-db test_data/graph/micb-kir3dl1_HG003.gaf.db \
  --sample CHM13 \
  --region "CHM13#0#chr6:31350873-31363898" \
  -o output.tsv
```

## Files

### Source files (checked in)

| File | Description |
|------|-------------|
| `micb-kir3dl1.gbz` | GBZ pangenome graph (CHM13 + GRCh38 refs, 46 samples) |
| `micb-kir3dl1_HG003.gaf` | GAF alignments for sample HG003 (12,439 reads) |

### Derived files (checked in, can be regenerated)

| File | Description |
|------|-------------|
| `micb-kir3dl1.gbz.db` | GBZ database (built by `gbz2db`) |
| `micb-kir3dl1_HG003.gaf.db` | GAF database (built by `gaf2db`) |
| `micb-kir3dl1.sites.vcf.gz` | Sites VCF from `build-snarl-catalog` (917 snarl sites) |
| `micb-kir3dl1.sites.vcf.gz.tbi` | Tabix index for sites VCF |
| `ref.fa.gz` | Reference FASTA (bgzipped, N-padded to genome coordinates) |
| `ref.fa.gz.fai` | FASTA index |
| `ref.fa.gz.gzi` | bgzip index |

### Notes on the reference FASTA

The GBZ subgraph has offset paths (e.g., CHM13 chr6 starts at genome position
31,350,872).  The reference FASTA is N-padded so that FASTA coordinates match
genome coordinates.  This is required because pgphase tiles the FASTA into
chunks and uses those coordinates for GBZ interval queries.

## Regenerating derived files

```bash
# GBZ database
third_party/gbz-base/target/release/gbz2db \
  test_data/graph/micb-kir3dl1.gbz \
  test_data/graph/micb-kir3dl1.gbz.db

# GAF database
third_party/gbz-base/target/release/gaf2db \
  -r test_data/graph/micb-kir3dl1.gbz \
  -o test_data/graph/micb-kir3dl1_HG003.gaf.db \
  --overwrite \
  test_data/graph/micb-kir3dl1_HG003.gaf

# Sites VCF (requires vg)
./pgphase build-snarl-catalog \
  --ref-sample CHM13 \
  -o test_data/graph/micb-kir3dl1.sites.vcf.gz \
  test_data/graph/micb-kir3dl1.gbz

# Reference FASTA (requires vg + samtools)
# 1. Extract CHM13 paths from the GBZ
vg paths -x test_data/graph/micb-kir3dl1.gbz -S CHM13 -F > /tmp/chm13_paths.fa

# 2. N-pad each contig so FASTA coordinates match genome coordinates
python3 -c "
import re, sys
seqs = {}
with open('/tmp/chm13_paths.fa') as f:
    name = None
    for line in f:
        if line.startswith('>'):
            m = re.match(r'>(\S+)\[(\d+)\]', line)
            name, offset = m.group(1), int(m.group(2))
            seqs[name] = (offset, [])
        else:
            seqs[name][1].append(line.strip())
for name, (offset, parts) in seqs.items():
    seq = ''.join(parts)
    print(f'>{name}')
    padded = 'N' * offset + seq
    for i in range(0, len(padded), 80):
        print(padded[i:i+80])
" > /tmp/ref_padded.fa

# 3. Compress and index
bgzip -c /tmp/ref_padded.fa > test_data/graph/ref.fa.gz
samtools faidx test_data/graph/ref.fa.gz
```
