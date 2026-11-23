# rusolve-mip

The legacy native branch-and-bound: node LPs through `rusolve-lp`, cuts, heuristics, restarts, racing. Since commit 4843fcc9 the default MIP engine is the HiGHS port in `crates/rusolve-ekk/src/mip`; this crate runs only with `RUSOLVE_MIP_HIGHS=0` and supplies instance classification (`classify`).

## What lives here
- `src/lib.rs`, `src/problem.rs`: `MipProblem` and the crate overview.
- `src/bnb.rs`, `src/branching.rs`, `src/nodequeue.rs`, `src/pseudocost.rs`: tree search.
- `src/presolve.rs`, `src/probing.rs`, `src/propagate.rs`, `src/conflict.rs`, `src/redcost.rs`, `src/restart.rs`, `src/submip.rs`, `src/alns.rs`, `src/race.rs`.
- `src/cuts/`, `src/heuristics/`, `src/engine/`, `src/engine.rs`, `src/engine_select.rs`: separators, primal heuristics, LP engine adapters.
- `src/classify.rs`: structure classification behind `rusolve classify`.
- `examples/node_lp_repro.rs`.

## Invariants
- New MIP capability goes to the port (`crates/rusolve-ekk/src/mip/AGENTS.md`). A technique that exists only here counts as not in the default solver.
- Nothing here may change the default (port) search: identity 20/20 stays.
- Cite a paper or other public source for every technique (provenance rule in the root `AGENTS.md`).

## How to test
- `cargo test --release -p rusolve-mip --lib`.
- Correctness gate against HiGHS and SCIP: `tools/diff/mip_diff.py`.
- Run the native engine from the CLI: `RUSOLVE_MIP_HIGHS=0 target/release/rusolve-cli mip <file> --time-limit 60`.

## Change recipes
1. Add a cut family: new file in `src/cuts/`, register it where the separators are listed, test validity on a crafted model, then run `tools/diff/mip_diff.py`.
2. Change branching: edit `src/branching.rs`; compare node counts under a node limit, not wall time.
3. Fix a classification: edit `src/classify.rs`; the CLI router reads it from `crates/rusolve-cli/src/route.rs`.

## Pitfalls
- Wall time hides effects under about 12 percent; compare node counts or instructions retired (`tools/bench/abtest.py`).
- Several ideas first built here were later ported; check the port first (`docs/agent/codemap.md`) before extending this crate.
- Further reading: the module docs in `src/lib.rs`.
