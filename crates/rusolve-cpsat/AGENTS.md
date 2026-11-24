# rusolve-cpsat

Binding to Google OR-Tools CP-SAT (Apache-2.0) through its C API: a model translator, a `solve` entry and a run-time loader of the OR-Tools shared library. A default workspace member: it needs no OR-Tools to build.

## What lives here
- `src/lib.rs`: `solve(lp, cp, params)`, `solve_proto`, `CpParams`, `CpOutcome`, `CpStatus`, `CpError`; the `unsafe` call through the loaded `SolveCpModelWithParameters` and the `free` of its response buffer.
- `src/runtime.rs`: opens the library once (`libloading`), looks the symbol up, finds its version, reports `availability()`. The search order is documented in that module and nowhere else: change it there.
- `src/translate.rs`: `(&HighsLp, Option<&CpModel>) -> CpModelProto`, with typed refusals (`Refusal`) that name the first offending column or row.
- `proto/`: vendored `cp_model.proto` and `sat_parameters.proto` of OR-Tools 9.15, with the license; `build.rs` generates the Rust types (`protox`, `prost-build`); nothing is linked.
- `tests/engine.rs`: translator and solve tests, one per `CpModel` object type.

## Invariants
- No OR-Tools at build time. The library is loaded on first use (`RUSOLVE_ORTOOLS_LIB`, install prefix, `ortools` Python package, system loader) and never unloaded; a missing library is `CpError::Unavailable`, which the CLI maps to RS0012.
- Translation errors (`Refused`) are reported before the library is touched.
- Integer columns with finite bounds, integral coefficients and right-hand sides only; anything else is a refusal, never a rounding.
- Never edit `crates/rusolve-io/src/cp.rs` from here: it is the contract with the CP model layer.
- Every `unsafe` block carries a SAFETY comment; the response buffer is C-allocated and freed with `free`.
- CP-SAT's interval does not enforce `start + size = end`; the translator adds it.

## How to test
- A translator test builds a small `HighsLp` or `CpModel` and compares the `CpModelProto` or the solved objective.
- `RUSOLVE_ORTOOLS_LIB=<dir>/libortools.so.9 cargo test --release -p rusolve-cpsat`; CLI level: `... cargo test --release -p rusolve-cli --test cpsat`. Without the library the solve tests skip and say so.

## Change recipes
1. Support a new CP object: add the arm in `src/translate.rs` and a test in `tests/engine.rs`, after the type exists in `rusolve-io/src/cp.rs`.
2. Pass a new parameter: add a field to `CpParams`, set it on `SatParameters` in `src/lib.rs`, read it in `crates/rusolve-cli/src/cpsat_solve.rs`.

## Pitfalls
- The OR-Tools library finds its own dependencies through `$ORIGIN`, taken from the path it was opened by: copy the whole library directory of the wheel (and, on Linux, the wheel's vendored Abseil directory two levels up), never a symlink to the one file.
- Upgrading OR-Tools means re-vendoring both `.proto` files from the same version as the library.
- `rusolve_ekk::prelude::*` exports lowercase constants (`c`, ...) that clash with patterns; import names explicitly.
- Further reading: `docs/design/cpsat-engine.md`.
