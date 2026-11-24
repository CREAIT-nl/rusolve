# CP-SAT as an opt-in MIP engine

`rusolve mip --opt mip_engine=cpsat` (and the MCP `solve` tool with the same option) solves a
pure integer program with Google OR-Tools CP-SAT (Apache-2.0, https://github.com/google/or-tools)
instead of the HiGHS port. It is opt-in: the registered default is `mip_engine=highs`, and the
default tests and the parity check never need OR-Tools. This is an
integration only: CP-SAT runs as a library and nothing of it is ported.

## What it takes

CP-SAT works on integers. The engine therefore refuses, with a typed error (RS0012, exit 3) that
names the first offender, any model with

* a continuous column ("column y (index 1) is continuous"),
* an infinite column bound,
* a non-integral constraint coefficient or row bound,
* an empty column domain, or a value outside the 64-bit range CP-SAT takes.

Fractional costs are allowed: the objective is then passed as a floating-point objective and the
`cpsat` object says `objective_kind: "floating_point"` (the integer objective is exact and is
preferred when every cost is integral). The refused models stay on the port.

## Finding the library

The engine is always compiled in and needs no OR-Tools at build time. The OR-Tools shared library
is opened on first use (`libloading`, ISC licence, a thin wrapper over dlopen and LoadLibrary),
once per process and never unloaded; the symbol `SolveCpModelWithParameters` is looked up
right after. Nothing is loaded for a model that does not select `mip_engine=cpsat`. `libloading`
was chosen over raw dlopen because it covers the three platforms with one reviewed `unsafe`
surface (and the Windows altered search path, which makes the library find its own neighbours).

The search order below is the single source of truth (code: `crates/rusolve-cpsat/src/runtime.rs`).
The first candidate that loads wins. The library name is `libortools.so.9` (Linux),
`libortools.9.dylib` (macOS) or `ortools.dll` (Windows).

1. `RUSOLVE_ORTOOLS_LIB`, the path of the library file. When it is set, only this path is tried.
2. `<exe dir>/../lib/rusolve/ortools/<libname>`: where `install.sh --with-cpsat` puts it.
3. `<exe dir>/../../ortools/.libs/<libname>`: the `ortools` PyPI wheel, next to
   `site-packages/rusolve/_bin/rusolve-cli`. In-process Python (the exe is the interpreter) passes
   the `ortools` package directory instead.
4. The system loader default (`libortools.so.9` by name).

The library finds its own dependencies through `$ORIGIN`, so the directory has to be kept whole.
When none loads, the error is `RS0012` (exit 3) naming every path searched and the two install
routes, `install.sh --with-cpsat` and `pip install 'rusolve[cpsat]'`. A library without the
symbol is a typed error too. `rusolve-cli version --format json` (and `--version`) reports
availability without solving: `cpsat: {available, library, version, searched, ...}`.

The OR-Tools version in the `cpsat` envelope object comes from the best runtime source: a `VERSION`
file next to the library, else the `ortools-<v>.dist-info` of the wheel, else the library's file
name, else `"unknown"`; `ortools_library` is the resolved path. The protobuf types for the model
and the parameters are generated at build time (`protox` and `prost-build`, about 15 CPU seconds
of a 108 CPU-second cold release build) from the vendored 9.15 `cp_model.proto` and
`sat_parameters.proto` (`crates/rusolve-cpsat/proto/`, with the OR-Tools license). They are not
pre-generated: a checked-in generated file would drift from the protos and need review noise, and
the cost is small. The protos are wire-compatible within OR-Tools major 9, which is why the
Python extra is `ortools>=9.15,<10`.

The binding is the C API (`cp_solver_c.h`): a serialized `CpModelProto` and `SatParameters`
go in, a serialized `CpSolverResponse` comes out. The `unsafe` code is the call through the loaded
symbol, the dlopen and the `free` of the response buffer, each with a SAFETY comment.

An unknown engine name is a usage error (RS0004, exit 2).

## Answers

* Time limit, threads (CP-SAT workers), seed, relative and absolute gap map to the matching
  `SatParameters`. The relative gap default is the port's (`mip_rel_gap=1e-4`); pass
  `--opt mip_rel_gap=0` for a proof of optimality. `--threads` defaults to 1 here as everywhere.
* Statuses: `OPTIMAL` is `Optimal`, `FEASIBLE` is `TimeLimit` with the incumbent, `INFEASIBLE`
  is `Infeasible`, `UNKNOWN` is `TimeLimit` without a point; a model CP-SAT rejects is an error.
* Every point goes through `check.rs` on the original model before it is reported, and its
  objective must equal CP-SAT's. A failure is an internal error (RS0090, exit 4), never an
  answer. A bound that passes the checked point downgrades `Optimal`.
* The envelope has `engine: "cp-sat"` (the class stays `mip`) and a `cpsat` object: CP-SAT's
  `status`, `bound`, `objective`, `wall_time_s`, `workers`, `ortools_version`, `objective_kind`,
  `conflicts`, `branches`, `aux_variables`, the limits passed, and `checked` (the verdict of our
  checker). Schema: `docs/schema/rusolve-result-1.schema.json`.
* Refused with the engine: `--start`, `--basis`, `--add-rows`, `--wasm-hook`, `--watchdog`,
  `--log-format jsonl`, `miqp_oa`, general constraints (RS0003, exit 2). MCP cancellation does not
  interrupt a running CP-SAT solve.

## Calling it from the constraint-programming model layer

```rust
rusolve_cpsat::solve(
    lp: &rusolve_ekk::prelude::HighsLp,
    cp: Option<&rusolve_io::cp::CpModel>,
    params: &rusolve_cpsat::CpParams,
) -> Result<rusolve_cpsat::CpOutcome, rusolve_cpsat::CpError>
```

`cp` carries the interval, no-overlap, cumulative, all-different, element and table objects of
`rusolve_io::cp`; they are translated next to the linear rows (`crates/rusolve-cpsat/src/translate.rs`,
one unit test per object type in `tests/engine.rs`). The MIP CLI glue passes `None`; the
CP router (`crates/rusolve-cli/src/cp_solve.rs`, Python and MCP) calls `solve` with `Some(..)` when
`mip_engine=cpsat`, and checks the point with `CpModel::check` on the user's objects as well as the
rows; a failed check is RS0090, never an answer. The default for CP models stays the linearization,
and its RS4013 refusal hints `mip_engine=cpsat`. Semantics fixed in the translator: CP-SAT's interval does not enforce
`start + size = end`, so the translator adds it; a zero-size interval in `no_overlap` uses an
auxiliary interval to keep the contract of `cp.rs`; a zero-size interval in `cumulative`
consumes nothing.

## Not done

* Bundling OR-Tools in our wheel: it is a separate install (`pip install 'rusolve[cpsat]'`).
* Interrupting a running solve from MCP cancellation and progress callbacks.
* Porting CP-SAT layers.
