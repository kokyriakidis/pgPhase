# pgPhase

## Portable binary bundle (Linux/macOS)

Build a self-contained runtime bundle (binary + copied shared libraries + launcher):

```bash
make portable-bundle
```

Create a release artifact (build + checks + portable bundle + tarball):

```bash
make release
```

Strict variant (fails if validation checks fail):

```bash
make release-strict
```

`release-strict` uses fixture golden outputs:

```text
test_data/expected/hifi_collect_expected.tsv
test_data/expected/hifi_collect_expected.phased.vcf
test_data/expected/ont_collect_expected.tsv
test_data/expected/ont_collect_expected.phased.vcf
```

Regenerate goldens (for intentional baseline updates):

```bash
./pgphase collect-bam-variation -t 1 --hifi --include-filtered test_data/chr11_2M.fa test_data/HG002_chr11_hifi_test.bam -o test_data/expected/hifi_collect_expected.tsv
./pgphase collect-bam-variation -t 1 --ont  --include-filtered test_data/chr11_2M.fa test_data/HG002_chr11_ont_test.bam  -o test_data/expected/ont_collect_expected.tsv
./pgphase collect-bam-variation -t 1 --hifi --include-filtered test_data/chr11_2M.fa test_data/HG002_chr11_hifi_test.bam -o test_data/expected/hifi_collect_expected.tsv --phased-vcf-output test_data/expected/hifi_collect_expected.phased.vcf
./pgphase collect-bam-variation -t 1 --ont  --include-filtered test_data/chr11_2M.fa test_data/HG002_chr11_ont_test.bam  -o test_data/expected/ont_collect_expected.tsv  --phased-vcf-output test_data/expected/ont_collect_expected.phased.vcf
```

Output directory:

```text
dist/pgphase-<os>-<arch>/
  bin/pgphase
  bin/pgphase-portable
  lib/*.(so|dylib)*
  README.portable.txt
dist/pgphase-<os>-<arch>.tar.gz
```

Run:

```bash
./dist/pgphase-<os>-<arch>/bin/pgphase-portable --help
```

Notes:

- Linux bundle portability depends on compatible CPU/kernel/glibc ABI.
- macOS bundle portability depends on compatible CPU architecture/runtime ABI.
- Build separate bundles for `arm64` and `x86_64` macOS targets when needed.