# pdlp

First-order LP solvers: cuPDLP-C and HiPDLP (PDHG with restarts), with the wrappers `Highs::run` uses.

## What lives here
- `cupdlp/`: the C solver translated file by file (`solver/`, `linalg/`, `step.rs`, `restart.rs`, `scaling.rs`, `proj.rs`, `cs.rs`, `utils/`); `cupdlp_wrapper/`: `CupdlpWrapper.cpp` (`formulate.rs` builds the PDLP problem).
- `hipdlp/`: HiPDLP (`pdhg/`, `restart.rs`, `scaling.rs`, `logger.rs`, `solver_results.rs`); `hipdlp_wrapper.rs`.
- `oracle_trace.rs`: iteration trace emitter (`RUSOLVE_PDLP_TRACE`).

## HiGHS files mirrored
`pdlp/CupdlpWrapper.cpp`, `pdlp/HiPdlpWrapper.cpp`, `pdlp/cupdlp/cupdlp_*.c`, `pdlp/hipdlp/*`. Per-file table in `docs/agent/codemap.md`.

## Invariants
- Same iterates as the oracle: restart criteria, step sizes and scaling must match bit for bit; no FMA, no reassociated reductions.
- C structs become Rust structs with the C names (`CUPDLPdata`, `CUPDLPiterates`); keep the names so the trace diff points at the right place.
- Not part of the default MIP search; identity is unaffected, but run it for any shared-code change.

## How to test
- `cargo test --release -p rusolve-ekk --test t_w3_cupdlp --test t_w3_hipdlp`.
- Unit tests live next to the code (`cargo test --release -p rusolve-ekk --lib pdlp`).
- Options: the `pdlp_*` entries of `docs/agent/reference.md` (registry names, set with `--opt name=value`).

## Change recipes
1. Trace diverges: use `RUSOLVE_PDLP_TRACE` on the failing instance and diff against the fixture; open the C function named in the doc line.
2. Add a PDLP option: `docs/skills/contrib-add-option/SKILL.md`; the `pdlp_*` fields are read in `cupdlp_wrapper/` (for example `formulate.rs`) and `hipdlp_wrapper.rs`.
3. Translate a missing helper: `docs/skills/contrib-port-function/SKILL.md`; oracle patch `tools/highs-oracle/highs-oracle-w3-hooks.patch`.

## Pitfalls
- `pdlp/mod.rs` also reads `RUSOLVE_PDLP_TRACE`; enable it only on the failing instance.
- The C macros and the dropped GPU branches are listed in `cupdlp/defs.rs` and the module doc of `pdlp/mod.rs`.
