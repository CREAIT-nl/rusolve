# rusolve-model

The append-only model builder: the Rust core of the Python `rusolve.Model` API. Design in `docs/design/python-api.md`; measurements in `crates/rusolve-py/bench/`.

## What lives here
- `src/lib.rs`: `ModelBuilder` (columns and rows are appended, values checked on arrival, duplicates merged and zeros dropped when a row is finished) and the `take` step that moves the model into a `HighsLp` in one O(nnz) pass.
- `src/truncate.rs`: truncating a builder back to an earlier size.
- `src/tests.rs`: unit tests; `benches/builder.rs`: criterion benchmark (`harness = false`).

## Invariants
- Nothing is O(num_col) per row; the builder never stores the model twice after `take`.
- It produces `rusolve_ekk::lp_data::highs_lp::HighsLp`, so it must keep building exactly what the port's own `pass_model` expects.
- Input validation (NaN, infinite coefficients, index out of range) happens at the edge, with a typed error the Python layer maps to exceptions.
- No search code: identity 20/20 is unaffected, but run the tests.

## How to test
- `cargo test --release -p rusolve-model`.
- Performance: `cargo bench -p rusolve-model --bench builder`; a hot-path change needs before and after numbers.
- End to end through Python: `crates/rusolve-py/tests/test_model_core.py` (see `crates/rusolve-py/AGENTS.md`).

## Change recipes
1. Add a model element (say a new bound kind): add the buffer and the validating method in `src/lib.rs`, extend `take`, add a case in `src/tests.rs`, then expose it in `crates/rusolve-py/src/model.rs`.
2. Add a check at the edge: put it where the value arrives, return the existing error type, test the bad input.
3. Speed the builder: change the hot loop, run `benches/builder.rs`, keep the O(row length) merge.

## Pitfalls
- Indices are global across `take` calls: after a `take` the builder keeps only counts, so new columns and rows continue the numbering.
- Exact zeros are dropped and duplicates merged only when a row is finished; do not rely on a half-built row.
- Further reading: `docs/design/python-api.md`.
- Run `cargo fmt --all` before committing; CI checks it.
