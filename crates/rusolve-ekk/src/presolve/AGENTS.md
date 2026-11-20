# presolve

The HiGHS LP/MIP presolve (`HPresolve`), postsolve stack, symmetry detection, and the ICrash crash start.

## What lives here
- `hpresolve/`: one `HPresolve` struct, its `impl` split by C++ line range. Rules: `row_presolve.rs`, `col_presolve.rs`, `singletons.rs`, `dual_fixing.rs`, `eq_subst.rs`, `aggregator.rs`, `col_stuffing.rs`, `fourier_motzkin.rs`, `sparsify.rs`, `mip_dominated.rs`, `mip_probing.rs`, `clique_varbound.rs`; driver `presolve_loop.rs` and `run.rs`; HiGHS 1.16 additions in `*_h16.rs`.
- `hpresolve_analysis.rs`: rule logging and the `allow_rule_` switches; `hpresolve_clique_table.rs`; `presolve_component.rs`.
- `highs_postsolve_stack/`: every presolve reduction and its `undo`.
- `highs_symmetry/`: `HighsSymmetry` detection. `icrash.rs`, `icrash_x.rs`, `icrash_util.rs`: the crash start.

## HiGHS files mirrored
`presolve/HPresolve*.cpp`, `HighsPostsolveStack.cpp`, `HighsSymmetry.cpp`, `ICrash*.cpp`, `PresolveComponent.cpp`. See `docs/agent/codemap.md` for the per-file mapping.

## Invariants
- Reductions apply in the HiGHS order; a different order changes the presolved model and moves identity.
- Every reduction pushes a postsolve record, and `undo` restores primal, dual and basis values.
- `mipsolver` is `None` for LP presolve; MIP-only branches translate the C++ conditions unchanged.
- Rules not in HiGHS are options that default off; the 1.16 rules sit behind `h16_presolve_rules` and the `presolve_rule_off` bits (`src/lp_data/hconst.rs`).

## How to test
- `cargo test --release -p rusolve-ekk --lib presolve` (unit tests live in the rule files), and `--test t_g1_presolve_api --test t_g1_clique_normalise --test t_w1_icrash`.
- Profile: `RUSOLVE_PRESOLVE_PROF` (see `hpresolve/presolve_prof.rs`).
- Identity and A/B: `docs/skills/contrib-identity-abtest/SKILL.md`.

## Change recipes
1. Add a presolve rule: follow `docs/skills/contrib-add-presolve-rule/SKILL.md`.
2. A rule differs from HiGHS: compare the presolved dimensions on one model with `RUSOLVE_PRESOLVE_PROF`, then diff the rule body against its C++ lines.
3. Add a postsolve undo: add the record in `highs_postsolve_stack/`, push it where the reduction happens, test a model where the reduction fires and the solution must be restored.

## Pitfalls
- A change that only matters for MIP (`mipsolver` set) is invisible to LP tests; check it on a MIP instance with identity.
- 1.16 rules append extra `allow_rule_` slots under `h16_presolve_rules`; the 1.16 setup no longer decodes `presolve_rule_off` the 1.15 way (`hpresolve_analysis.rs`).
- The 1.16 presolve changes (HiGHS PRs #3167 #3210 #3270 and others) live in `presolve_h16.rs`, `initial_sweep_h16.rs` and `row_bounds_h16.rs`, each citing its PRs.
