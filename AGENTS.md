# AGENTS.md

Project-specific guidelines for AI agents working on pgphase.

## Project Overview

pgphase is a C++17/Rust tool for variant calling and haplotype phasing from long reads. It has two pipelines:
- `collect-bam-variation` — SNP/indel calling from BAM/CRAM alignments
- `collect-graph-variation` — variant calling from pangenome graphs (GBZ + GAF)
- `build-snarl-catalog` — preprocess GBZ graphs into phasing site catalogs

## Build System

- C++17 compiled with GCC, linked against htslib, zlib, pthreads
- Third-party libraries: WFA2-lib, abPOA, edlib (git submodules under `third_party/`)
- Rust FFI component: `src/gbz_ffi/` (static library linked into the C++ binary)
- Rust tools: `third_party/gbz-base/` (query, gaf2db, gbz2db binaries)
- Build order: `make third-party-libs` → `make gbz-base` → `make -j$(nproc)`

## Testing

- `make check` — validation gates (requires `test_data/`)
- `make unit-tests` — builds and runs unit test binaries

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
