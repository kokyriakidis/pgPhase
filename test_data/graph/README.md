# Graph pipeline test data

Small pangenome graph test case for `collect-graph-variation`, based on the
MICB/KIR3DL1 haplotype sampling test from
[gbz-base](https://github.com/jltsiren/gbz-base/tree/main/test-data).

## Source files (checked in)

| File | Description |
|------|-------------|
| `micb-kir3dl1.gbz` | GBZ pangenome graph (CHM13 + GRCh38 refs, 46 samples) |
| `micb-kir3dl1_HG003.gaf` | GAF alignments for sample HG003 (12,439 reads) |

## Derived files (checked in, can be regenerated)

| File | Description |
|------|-------------|
| `micb-kir3dl1.gbz.db` | GBZ database (built by `gbz2db`) |
| `micb-kir3dl1_HG003.gaf.db` | GAF database (built by `gaf2db`) |
| `micb-kir3dl1.sites.vcf.gz` | Sites VCF from `build-snarl-catalog` (917 snarl sites) |
| `micb-kir3dl1.sites.vcf.gz.tbi` | Tabix index for sites VCF |
| `ref.fa.gz` | Reference FASTA (bgzipped, N-padded to genome coordinates) |
| `ref.fa.gz.fai` | FASTA index |
| `ref.fa.gz.gzi` | bgzip index |

## Running the test

```bash
./pgphase collect-graph-variation \
  --ref test_data/graph/ref.fa.gz \
  --sites test_data/graph/micb-kir3dl1.sites.vcf.gz \
  --gbz-db test_data/graph/micb-kir3dl1.gbz.db \
  --gaf-db test_data/graph/micb-kir3dl1_HG003.gaf.db \
  --sample CHM13 \
  --region "CHM13#0#chr6:31350873-31363898" \
  --region "CHM13#0#chr19:18926906-18941251" \
  --phased-vcf-out phased.vcf \
  -o output.tsv
```

The `--region` flags restrict queries to the subgraph's coordinate range.
Without them, chunks outside the GBZ path range return no reads (the FFI
query silently skips intervals not covered by the reference path).

Expected output: 292 candidate variant sites.

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
vg paths -x test_data/graph/micb-kir3dl1.gbz -S CHM13 -F | \
  python3 -c "
import sys, re
for line in sys.stdin:
    if line.startswith('>'):
        m = re.match(r'>(\S+)\[(\d+)\]', line)
        name, offset = m.group(1), int(m.group(2))
        print(f'>{name}')
        seq = []
    else:
        seq.append(line.strip())
        if not line.strip(): continue
# flush
" > /dev/null  # see the python script in the repo history for the full version
```
