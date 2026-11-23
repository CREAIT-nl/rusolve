# rusolve-lp

The legacy native LP engine: dense and sparse LU, product-form and Forrest-Tomlin updates, primal and dual simplex, interior point, presolve. It is NOT the default LP path: `rusolve solve` runs the HiGHS-port simplex (`crates/rusolve-cli/src/ekk.rs`); this crate runs with `--legacy-lp` or `--ipm`.

## What lives here
- `src/lu.rs`, `src/lu_ft.rs`, `src/pfi.rs`, `src/dense.rs`, `src/hvec.rs`: factorization layer.
- `src/simplex.rs`, `src/dual.rs`: primal and dual simplex; `src/presolve.rs`, `src/scaling.rs`, `src/rng.rs`.
- `src/ipm/`: interior point.
- `benches/numerics.rs`, `benches/ipm_chol.rs`: criterion benchmarks (`harness = false` in `Cargo.toml`).

## Invariants
- Reference and fallback only. New LP speed work goes to the port; this crate is the benchmark floor and the engine behind `rusolve-mip`'s native search (`RUSOLVE_MIP_HIGHS=0`).
- Hot-path changes need a criterion run (AGENTS.md at the repo root: performance is measured, never claimed).
- `unsafe` needs a written justification in the commit message.
- The default MIP search does not use this crate, so identity 20/20 is unaffected, but still run the check.

## How to test
- `cargo test --release -p rusolve-lp --lib`.
- `cargo bench -p rusolve-lp --bench numerics` for the numerics.
- Differential: `tools/diff/simplex_diff.py`, `tools/diff/dual_diff.py`, `tools/diff/lu_diff.py` (objective versus highspy, rel 1e-6; they run `target/release/rusolve-cli`).

## Change recipes
1. Change a pivot rule: edit `src/simplex.rs` or `src/dual.rs`, run `tools/diff/simplex_diff.py` for presolve on and off.
2. Add a presolve reduction: add it to `src/presolve.rs` with its postsolve step and a unit test; check `tools/diff/lu_diff.py` on Netlib.
3. Speed a factorization kernel: change `src/lu.rs`, compare with `benches/numerics.rs` before and after.

## Pitfalls
- Both `--legacy-lp` and `--ipm` are native-only flags; they are marked `Native` in the registry of `crates/rusolve-cli/src/commands.rs`.
- A silent wrong answer once hid in a presolve-off LP (the MIP search runs presolve-off LPs), so differential gates run presolve on and off.
- Further reading: the provenance headers of `src/dual.rs`, `src/lu.rs` and `src/ipm/mod.rs` cite the papers and HiGHS sources each part follows.
