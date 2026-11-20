# rusolve-core

Data structures shared by every other crate. No solver logic, no dependencies (the sparse linear algebra is hand-rolled, see the crate docs in `src/lib.rs`).

## What lives here
- `src/lib.rs`: `CscMatrix` and `CsrMatrix` with conversions both ways.
- `src/csc_builder.rs`: incremental CSC construction.
- `src/sparse_lp.rs`: `SparseLp`, the LP front door the legacy path takes (`Model::to_sparse_lp` in rusolve-io).
- `src/qp_convexity.rs`: convexity test for quadratic objectives and rows (read by `crates/rusolve-cli/src/qp.rs` and `crates/rusolve-py/src/quad.rs`).
- `tests/mps_round_trip.rs`: MPS round trip; it uses rusolve-io as a dev-dependency only.

## Invariants
- Nothing here mirrors a HiGHS file; the faithful HiGHS port lives in `crates/rusolve-ekk` (see its `AGENTS.md`).
- Every other crate depends on this one, so keep the public types stable and the crate free of dependencies. Do not add one without a measured reason.
- Does not move the default search: pure data structures, no hooks, no levers.

## How to test
- `cargo test --release -p rusolve-core`.
- A change to a matrix type must also pass `cargo test --release -p rusolve-io -p rusolve-lp -p rusolve-ekk --lib`, because they build on it.
- Identity (`tools/bench/identity.sh`) only matters if the port's input changes; run it locally after `cargo build --release -p rusolve-cli`.

## Change recipes
1. Add a matrix helper: put it next to the type in `src/lib.rs`, add a unit test in the same file, check that `rusolve-ekk` (which has its own `HighsSparseMatrix`) does not need it.
2. Add a field to `SparseLp`: update the constructors in `src/sparse_lp.rs`, then every reader in `crates/rusolve-io/src/model.rs` and every consumer; the compiler lists them.
3. Add a convexity rule: extend `src/qp_convexity.rs` with a test of a convex and a nonconvex case; the caller uses `crates/rusolve-cli/src/qp.rs`.

## Pitfalls
- `CscMatrix` and `CsrMatrix` conversions are used on hot paths in `crates/rusolve-lp`; a change in ordering inside a column breaks pivot-for-pivot comparisons downstream.
- The dev-dependency on rusolve-io is a cycle on purpose (`Cargo.toml` comment); do not turn it into a normal dependency.
- Further reading: the module docs of `src/sparse_lp.rs` and `src/csc_builder.rs`.
