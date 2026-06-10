# AGENTS.md

Project-specific guidelines for AI agents working on pgphase.

## Project Overview

pgphase is a C++17/Rust tool for variant calling and haplotype phasing from long reads. It has three commands:

- `collect-bam-variation` — SNP/indel calling from BAM/CRAM alignments
- `collect-graph-variation` — variant calling from pangenome graphs (GBZ + GAF)
- `build-snarl-catalog` — preprocess GBZ graphs into phasing site catalogs

## Build System

- C++17 compiled with GCC, linked against htslib, zlib, pthreads.
- Third-party libraries: WFA2-lib, abPOA, edlib (git submodules under `third_party/`).
- Rust FFI component: `src/gbz_ffi/` (static library linked into the C++ binary).
- Rust tools: `third_party/gbz-base/` (query, gaf2db, gbz2db binaries).
- Build order: `make third-party-libs` → `make gbz-base` → `make -j$(nproc)`.

## Testing

- `make check` — validation gates (requires `test_data/`).
- `make unit-tests` — builds and runs unit test binaries.

## Key Files

| Path | Purpose |
|---|---|
| `CHECKPOINT.md` | Evaluation results, design decisions, noise filtering rationale |
| `src/noise_filter.hpp/cpp` | Shared noise detection (BAM + graph pipelines) |
| `src/graph_bam_adapter.hpp/cpp` | Graph-to-phasing bridge, noise filter wiring |
| `src/graph_collect.cpp` | Graph pipeline worker loops |
| `src/collect_var.cpp` | BAM pipeline candidate classification |
| `src/collect_phase.cpp` | K-means phasing, category bitmasks |
| `src/phasing_types.hpp` | Core types: VariantKey, Interval, constants |
| `benchmark/` | Snakemake benchmark pipeline |
| `scripts/` | Evaluation and benchmark helper scripts |

## Architecture Notes

### Two pipelines, shared phasing core

Both pipelines produce a `PhasingChunk` with `CandidateVariant` entries, then share the same k-means phasing and chunk stitching code. The graph pipeline converts snarl site observations into candidates via `build_graph_chunk` in `graph_bam_adapter.cpp`.

### Noise filtering

Reference-based noise detection (homopolymer, tandem repeat, SDUST low-complexity) is shared between pipelines via `noise_filter.hpp/cpp`. The BAM pipeline has additional read-level noise filtering (XID clustering, noisy-reads-ratio, MSA recall) that depends on CIGAR data unavailable in the graph pipeline. See `CHECKPOINT.md` "Noise Filtering: BAM vs Graph Pipeline" for the full rationale.

### Category bitmask system

Candidates are classified with both a `VariantCategory` enum and a `lcd_var_i_to_cate` bitmask. The bitmask controls which candidates enter k-means phasing:

- `kCandGermlineClean` (SNP | Indel | Hom) — used for graph pipeline k-means.
- `kLongcalldRepHetVar` (0x010) — repeat/homopolymer indels, excluded from k-means.
- `kCandGermlineVarCate` — includes noisy-recalled candidates, used in BAM pipeline second k-means.

### Thread safety

`faidx_t` is not thread-safe. Each worker thread opens its own handle via `load_reference_index(opts.ref_fasta)`. The graph pipeline does this in both `process_graph_chunk_batch` and `process_graph_chunk_batch_indexed_gaf`.

---

## Behavioral Guidelines

Imported from [andrej-karpathy-skills](https://github.com/multica-ai/andrej-karpathy-skills/blob/main/CLAUDE.md).

### 1. Think Before Coding

Don't assume. Don't hide confusion. Surface tradeoffs.

- State assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them — don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First

Minimum code that solves the problem. Nothing speculative.

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical Changes

Touch only what you must. Clean up only your own mess.

When editing existing code:

- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it — don't delete it.

When your changes create orphans:

- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

Define success criteria. Loop until verified.

Transform tasks into verifiable goals:

- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:

```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

---

## C++ Conventions

Adapted from [XOOS C++ rules](https://github.com/Roche-DIA-RDS-CSI/XOOS).

### Naming

| Element | Convention | Example |
|---|---|---|
| Files | snake_case | `noise_filter.cpp`, `graph_sites.hpp` |
| Variables | snake_case | `ref_beg`, `chunk_rows` |
| Private members | trailing underscore | `loaded_tid_`, `seq_` |
| Functions | snake_case | `find_low_complexity_intervals()` |
| Classes/structs | PascalCase | `GraphSiteMeta`, `PhasingChunk` |
| Constants | kCamelCase | `kSdustThreshold`, `kDefaultNoisyRegMaxXgaps` |
| Namespaces | snake_case | `pgphase_collect` |
| Bitmask flags | kCamelCase | `kCandCleanHetSnp`, `kLongcalldRepHetVar` |

### Header rules

- Use `#ifndef` / `#define` include guards (not `#pragma once`).
- Use `"file.h"` for project-local files.
- All code under `pgphase_collect` namespace.
- Include order: (1) corresponding header, (2) project headers, (3) third-party headers (`sdust.h`, `cgranges.h`), (4) system/htslib headers, (5) standard library.

### Implementation rules

- Implement functions in `.cpp`, not `.hpp` (except templates).
- Const correctness: const parameters, const references, const locals.
- Use RAII for resource management (`std::unique_ptr` with custom deleters for htslib handles).
- No magic numbers — use named `constexpr` constants.
- Prefer `static` file-scope functions over anonymous namespaces for small helpers.
- Match existing style in the file you're editing.

### Error handling

- Use `std::runtime_error` for fatal errors.
- Use `std::cerr` for warnings and verbose output (gated by `opts.verbose`).
- Make error messages actionable — include the file path or parameter that failed.

### Comments

- Explain "why", not "what".
- Use `///` or `/** */` for public API documentation.
- Use `//` for inline implementation notes.
- Do not comment obvious code.

### Testing

- Test binaries: `test_phase_block_stitch`, `test_graph_sites`, `test_graph_bam_adapter`.
- Tests are standalone `.cpp` files in `src/` (e.g., `src/test_phase_block_stitch.cpp`).
- Run all: `make unit-tests`.
- When adding a new `.o` dependency, update both the main `pgphase` target and any test targets that link the dependent object.

---

## Shell Script Conventions

- Start every script with `set -euo pipefail`.
- Use `#!/usr/bin/env bash`.
- Quote all variable expansions: `"${var}"`.
- Use `readonly` for constants.
- Provide `--help` / usage output.
- Use `die()` or similar for fatal errors with exit code 1.

---

## Git Conventions

- Commit messages: imperative mood, first line ≤72 chars, blank line before body.
- Add co-author: `Co-authored-by: Ona <no-reply@ona.com>`.
- Only stage files relevant to the current task.
- Do not commit files modified before the task began unless directly related.

---

## Design Decisions Log

All significant design decisions, evaluation results, and filtering rationale are recorded in `CHECKPOINT.md`. Update it when:

- Adding or removing a filtering stage.
- Changing how a pipeline classifies or handles candidates.
- Running evaluation benchmarks with new results.
- Discovering pitfalls or gotchas worth preserving.

---

## Code Review Checklist

Before submitting changes:

- [ ] Every changed line traces to the user's request.
- [ ] `make -j$(nproc)` compiles with zero errors.
- [ ] `make unit-tests` passes.
- [ ] No new warnings introduced.
- [ ] Makefile updated if new `.cpp` files added (both main target and test targets).
- [ ] `CHECKPOINT.md` updated if design decisions were made.
- [ ] No magic numbers — named constants used.
- [ ] Const correctness maintained.
- [ ] Thread safety considered (no shared mutable state without synchronization).
- [ ] htslib handles wrapped in RAII (`std::unique_ptr` with custom deleter).
