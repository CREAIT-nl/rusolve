# lp_data

The `Highs` driver and everything around it: model and solution types, the option registry, info, callbacks, run orchestration, the modification API, rays, ranging, IIS, multi-objective, elasticity, and the glue to presolve and solvers.

## What lives here
- `highs_run.rs` (`Highs::run`), `highs_solve.rs` (`solve_lp`), `highs_api*.rs` (the public surface, including `highs_api_mip.rs` for option set and get).
- Types: `hconst.rs`, `hstruct.rs`, `highs_lp.rs`, `highs_info.rs`, `highs_run_data.rs`, `highs_status.rs`, `highs_callback.rs`.
- Options: `highs_options.rs` (the `HighsOptions` struct and defaults) and `highs_option_records.rs` (the registry: name, range, default, source, `affects`, `tunable`).
- Interface and modification: `highs_interface_*.rs`, `highs_lp_modify.rs`, `highs_modify_api.rs`, `highs_pass_model_raw.rs`, `highs_set_solution.rs`.
- Analysis: `highs_iis/`, `highs_ranging/`, `highs_rays.rs`, `highs_ill_conditioning/`, `highs_elasticity.rs`, `highs_multi_objective.rs`, `highs_qp.rs`.
- Glue: `highs_presolve.rs`, `highs_postsolve.rs`, `highs_decompose.rs`, `highs_lp_concurrent.rs`, `highs_lp_fold.rs`.

## HiGHS files mirrored
`lp_data/Highs.cpp`, `HighsInterface.cpp`, `HighsLpUtils.cpp`, `HighsOptions.*`, `HighsInfo.*`, `HighsSolution.cpp`, `HighsIis.cpp`, `HighsRanging.cpp`, `HConst.h`, `HStruct.h`. Table: `docs/agent/codemap.md`.

## Invariants
- Option names and defaults match HiGHS; a rusolve lever carries a `source` naming its origin (a HiGHS backport, another solver, a paper, or `rusolve`), cites it in its description, and defaults to off (`Affects::Search` marks options that change the tree).
- `OPTION_RECORDS.len() + INTERNAL_FIELDS.len()` is asserted equal to the count in the test at the end of `highs_option_records.rs`; adding an option means bumping it.
- `docs/agent/reference.md` lists every option and is regenerated, never hand-edited.

## How to test
- `cargo test --release -p rusolve-ekk --test t_highs_run --test t_g1_modify --test t_g1_rays --test t_w1_iis --test t_w1_ranging` and the other `tests/t_g1_*.rs`.
- Option registry tests are in `highs_option_records.rs`.

## Change recipes
1. Add an option: `docs/skills/contrib-add-option/SKILL.md` (field, record, mapping, destructure, count, oracle pins, regenerate docs).
2. Add an API call: translate it in the matching `highs_api*.rs` with its HiGHS doc line and a `tests/t_g1_*.rs` golden.
3. Add an analysis that reads the solution: put it in its own file, call it from `highs_run.rs` behind an option, default off.

## Pitfalls
- Reading an option anywhere in the crate is `<expr>.<option>`; the codemap's option table (`docs/agent/codemap.md`) shows which files read it before you change a default.
- `highs_option_records.rs` is over 5000 lines; search by option name, and edit the record, the mapping line and the destructure together.
