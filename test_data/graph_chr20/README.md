# Graph pipeline test data — chr20:25,000,001–25,500,000

Small chr20 test case for `collect-graph-variation` using a 500 kb slice of
the CHM13 chr20 locus (HG002 HiFi reads aligned to the HPRC pangenome).

## Running the graph pipeline

```bash
./pgphase collect-graph-variation \
  --ref  test_data/graph_chr20/ref.fa.gz \
  --sites test_data/graph_chr20/chr20_25M.sites.vcf.gz \
  --gaf  test_data/graph_chr20/HG002_chr20_25M.coord.gaf.gz \
  --hifi \
  -o output.tsv
```

Expected output: **509 candidate variant sites** from 2,048 reads.

## Files

| File | Size | Description |
|------|------|-------------|
| `ref.fa.gz` | 261 KB | bgzipped reference FASTA; chr20 N-padded to genome coords |
| `ref.fa.gz.fai` | 23 B | FASTA index |
| `ref.fa.gz.gzi` | 6.2 KB | bgzip index |
| `chr20_25M.sites.vcf.gz` | 213 KB | Graph snarl sites VCF (tabix-indexed) |
| `chr20_25M.sites.vcf.gz.tbi` | 488 B | Tabix index |
| `HG002_chr20_25M.coord.gaf.gz` | 5.6 MB | Coordinate-sorted GAF, HG002 HiFi (tabix-indexed) |
| `HG002_chr20_25M.coord.gaf.gz.tbi` | 2.1 KB | Tabix index |

## Notes

**Reference FASTA**: The slice covers positions 25,000,001–25,500,000 on chr20.
The contig is N-padded from position 1 so that FASTA coordinates match genome
coordinates. This is required because pgphase tiles the reference into chunks
and uses absolute positions for tabix queries against the sites VCF and GAF.

**Sites VCF**: Subset from the full CHM13 chr20 snarl catalog with
`tabix ... CHM13#0#chr20:25000001-25500000`. The pipeline resolves
`CHM13#0#chr20` → `chr20` via suffix matching.

**GAF**: Coordinate-sorted and bgzipped, indexed with `tabix -p bed`.
The format has reference contig/start/end in columns 1–3.

## Regenerating

```bash
# Reference
samtools faidx ref_full.fa "CHM13#0#chr20:25000001-25500000" \
    | grep -v "^>" | tr -d '\n' > seq.txt
python3 -c "
seq = open('seq.txt').read().strip()
print('>chr20')
padded = 'N'*25000000 + seq
[print(padded[i:i+80]) for i in range(0, len(padded), 80)]
" | bgzip > ref.fa.gz
samtools faidx ref.fa.gz

# Sites VCF
tabix sites_full.vcf.gz "CHM13#0#chr20:25000001-25500000" \
    | bgzip > chr20_25M.sites.vcf.gz
tabix -p vcf chr20_25M.sites.vcf.gz

# GAF
tabix reads_full.coord.gaf.gz "chr20:25000001-25500000" \
    | bgzip > HG002_chr20_25M.coord.gaf.gz
tabix -p bed HG002_chr20_25M.coord.gaf.gz
```
