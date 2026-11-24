# mip

The HiGHS MIP solver: `HighsMipSolver`, LP relaxation, domain propagation, cut generation and separators, conflict analysis, primal heuristics (RINS, RENS, feasibility pump, feasibility jump, local MIP), symmetry, racing, plus the rusolve levers that sit beside them.

## What lives here
- `mip_solver.rs`, `mip_worker.rs`, `mip_solver_data/`, `search/`, `lp_relaxation/`, `node_queue.rs`, `pseudocost.rs`, `mip_stats.rs`.
- `domain/`, `implications.rs`, `clique_table/`, `conflict_pool.rs`, `conflict_resolution/`, `cdcl_uip.rs`, `cdcl_dive.rs`.
- Cuts: `cut_generation.rs`, `cut_pool.rs`, `cut_selection/`, `separation.rs`, `separator.rs` and one `*_separator.rs` per family.
- Heuristics: `primal_heuristics/`, `fpr_local_mip/`, `feasibility_jump.rs`, `fj_native.rs`, `scip_heuristics/`, `ls_int_lp_repair.rs`, `node_rounding.rs` (LP-free rounding at dive nodes, lever `mip_node_rounding`, called from `search/part2.rs`).
- `hooks.rs`: the `MipHooks` extension trait (`docs/design/hooks.md`); `race.rs`, `redcost_fixing.rs`, `debug_sol.rs`.
- `oa.rs`, `lazy.rs`: convex MIQP outer approximation (separator, exactness oracle, lazy-constraint loop; `docs/design/miqp-oa.md`), lever `miqp_oa` default off. They also hold the cone block (`OaCone`, `Aff`, optional disaggregation, `Polisher` hook, `LazyConstraints::before_resolve`) of the mixed-integer SOCP route, `docs/design/misocp-oa.md`, levers `misocp_oa*` default off; unit tests: `cargo test -p rusolve-ekk --lib mip::oa`.

## HiGHS files mirrored
`mip/HighsMipSolver*.cpp`, `HighsMipSolverData.cpp`, `HighsSearch.cpp`, `HighsLpRelaxation.cpp`, `HighsDomain.cpp`, `HighsCutGeneration.cpp`, `Highs*Separator.cpp`, `HighsPrimalHeuristics.cpp`, `HighsCliqueTable.cpp`, `HighsConflictPool.cpp`. Table: `docs/agent/codemap.md`.

## Invariants
- With one thread and a fixed seed the search is deterministic; node counts and LP iterations are the identity witness (`tools/bench/identity.sh`, 20 instances, must stay 20 of 20).
- Levers are options with default off, or default on only after a measured win with a held-out check; their sources are cited in `src/lp_data/highs_option_records.rs`.
- Hooks are empty by default; with none installed every call site is one `None` check.
- A time-limited instance is no identity witness; compare under a node limit.

## How to test
- Trace: `crates/rusolve-cli/tests/mip_trace.rs` replays `tests/fixtures/mip/*.trace` (needs `benches/data/miplib`); hooks: `tests/hooks_mip.rs`; per-feature: `conflict_resolution.rs`, `cut_selection.rs`, `fpr_local_mip.rs` under `crates/rusolve-cli/tests`.
- Unit tests sit in the files (`cargo test --release -p rusolve-ekk --lib mip`).
- A/B with `python tools/bench/abtest.py` on `tools/bench/pools/t2.txt`; see `docs/skills/contrib-identity-abtest/SKILL.md`.

## Change recipes
1. Add a separator: write `<name>_separator.rs` modeled on `modk_separator.rs`, add a `HighsSeparatorBox` variant in `src/mip/mod.rs`, push it in `separation.rs` behind a default-off option.
2. Add a heuristic: put it in `primal_heuristics/` and wire it in where `primal_heuristics/mod.rs` dispatches heuristics, behind a default-off option, and respect `mip_heuristic_effort`.
3. Add a hook: extend `MipHooks` in `hooks.rs` with a default body and a call site; keep `crates/rusolve-wasm` in step.

## Pitfalls
- `RUSOLVE_MIP_AC_FIXTURE` is not set in normal runs; the `.ac` fixtures are only a debugging reference for the analytic centre (`ipm_replay.rs`).
- Debug a wrong answer with the debug solution (`debug_sol.rs`): each step that cuts it off is reported on stderr as `HighsDebugSol: ...`, and `RUSOLVE_DEBUG_SOL_ABORT=1` panics at the first one.
- Further reading: `docs/design/hooks.md` (extension hooks), `docs/agent/codemap.md` (HiGHS file to Rust file map).
