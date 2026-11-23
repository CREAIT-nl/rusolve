# rusolve-io

File formats: readers and writers that turn a model file into a `Model`, plus the format-agnostic entry points that detect compression and format.

## What lives here
- `src/lib.rs`: `read_model_file` and friends: decompression by magic bytes (gzip, bzip2, xz, zstd), then reader choice by extension or content sniff.
- `src/mps.rs` (MPS reader and writer), `src/lp.rs` and `src/lp_write.rs` (CPLEX LP), `src/nl.rs` (AMPL), `src/mof.rs` (MathOptFormat JSON).
- `src/gc.rs`, `src/gc_lin.rs`, `src/conic.rs`: general constraints (indicator, SOS, piecewise) and conic sections.
- `src/model.rs`: the `Model` type that all readers produce.
- `examples/model_json.rs`; `tests/` with `dialects/`, `compressed/`, `cplex_lp/`, `mof/`, `nl/` corpora.

## Invariants
- Behavior must match what HiGHS reads: `tools/diff/lp_read_diff.py` and `tools/diff/corpus_diff.py` compare against highspy.
- All decompressors are pure Rust so the crate builds anywhere (WASM included). Do not add a C dependency.
- A reader bug can move the default search silently, so identity 20/20 (`tools/bench/identity.sh`) must stay green for any parser change.
- `crates/rusolve-ekk/src/io` is a different thing: HiGHS logging, not file formats.

## How to test
- `cargo test --release -p rusolve-io`; round trips in `tests/lp_write_roundtrip.rs` and `tests/gc_roundtrip.rs`.
- Differential: `python tools/diff/lp_read_diff.py` on the build host (needs highspy).

## Change recipes
1. Support a new MPS dialect keyword: add the parse arm in `src/mps.rs`, a sample under `tests/dialects/`, a case in `tests/dialects.rs`.
2. Add a CPLEX LP syntax case: edit `src/lp.rs`, add a sample to `tools/diff/data/lp_samples/` and run `tools/diff/lp_read_diff.py`.
3. Add a compression format: detect its magic in `src/lib.rs`, add a pure-Rust decoder to `Cargo.toml`, a file under `tests/compressed/` and a case in `tests/compressed.rs`.

## Pitfalls
- File order matters: the port reproduces HiGHS only when nonzeros keep their file order and duplicates resolve like HiGHS's MPS reader (`crates/rusolve-cli/src/ekk.rs` header lists the rules).
- Format detection never trusts the extension for compressed files; test a mislabelled file.
- Further reading: `docs/design/diagnostics.md` (parse diagnostics); the module headers of `src/mps.rs`, `src/lp.rs` and `src/lp_write.rs` list the format rules each follows.
