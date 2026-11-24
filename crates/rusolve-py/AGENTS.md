# rusolve-py

Python bindings, built with maturin. It is not in the workspace `default-members` because it needs a Python interpreter: only `-p rusolve-py` or maturin builds it.

## What lives here
- `src/lib.rs`: the PyO3 module; `src/model.rs`, `src/expr.rs`, `src/ingest.rs`, `src/quad.rs`, `src/solution.rs`, `src/diag.rs` (diagnostic codes and Python source-line locations), `src/callback.rs` (the `Model.solve(callback=...)` bridge to the HiGHS callback system, `MipHooks` and `run_with_lazy_constraints`; failure policy in `docs/design/hooks.md`, Python).
- `python/rusolve/__init__.py`: `rusolve.Model`, the `Highs` class (highspy-compatible); `python/rusolve/scipy.py`: `milp` and `linprog`; `python/rusolve/errors.py` (exception classes with `.code`/`.diagnostics`); `.pyi` stubs and `py.typed` (checked by `mypy --strict` in `tests/test_typing.py`).
- `tests/`: `test_basic.py`, `test_callbacks.py` (hooks, lazy TSP, `diagnose`), `test_diagnostics.py`, `test_errors.py`, `test_typing.py`, `test_model_api.py`, `test_model_core.py`, `test_qp.py`, `test_scipy.py`, `test_semi.py`, `test_gencons.py`.
- `examples/`: `tsp_lazy.py`, `user_cuts.py`, `prototype_branching.py` (shared `hooks_common.py`); the tests import it.
- `bench/`: `bench_api.py`, `compare_floors.py`, `compare_libs.py`; `README.md` has the build steps.

## Invariants
- The solver underneath is `rusolve-ekk`; this crate converts and forwards, it does not solve.
- `rusolve-cli` is a dependency: diagnostics (`check`, `explain`/`iis`, `ranging`, `feasrelax`, and `lint(with_solve=True)`) call `rusolve_cli::api`, which works on the port's LP with no Hessian, so each must call `refuse_quad` first (RS4001) on a model with a quadratic objective. Plain `lint()` is static and is allowed on a QP.
- General constraints (`Model.add_indicator`, `add_sos1/2`, `add_pwl`, `add_abs/max/min/and/or/not`) are stored as the typed `rusolve_io::GcModel` in `Model.gc` (`src/model/gencon.rs`); a solve with any of them linearizes a copy through `rusolve_cli::gc_solve` on a separate solver and never touches `self.port`. A call that would work on the plain rows and ignore them must call `refuse_gc` (RS4010).
- Scheduling and combinatorial objects (`Model.add_interval`, `add_no_overlap`, `add_cumulative`, `add_all_different`, `add_element`, `add_table`) are stored as the typed `rusolve_io::CpModel` in `Model.cp` (`src/model/cp.rs`); a solve goes through `rusolve_cli::cp_solve::solve_with` on a separate solver. A call that would ignore them must call `refuse_cp` (RS4012). Each `add_*` validates and undoes itself on RS4011.
- Model building goes through `crates/rusolve-model`; do not add a second builder.
- No search change here, so identity 20/20 is unaffected; options reach the port by name through the registry.
- Do not publish a wheel from a development branch (no PyPI).

## How to test
- Build and test locally: `maturin develop --release -m crates/rusolve-py/Cargo.toml` then `python -m pytest crates/rusolve-py/tests` (`README.md` has the venv setup).
- Rust side only: `cargo check -p rusolve-py`.

## Change recipes
1. Expose a model method: implement it in `src/model.rs` (reusing `rusolve-model`), add the Python wrapper in `python/rusolve/__init__.py`, a test in `tests/test_model_api.py`.
2. Expose a new option: nothing to do if it is in the registry (`rusolve options --json`); add a test that sets it by name.
3. Add a SciPy-compatible keyword: edit `python/rusolve/scipy.py` and `tests/test_scipy.py`.

## Pitfalls
- A Python hook must never unwind into the solver: every call goes through `Bridge::call` (`callback.rs`), which records the first exception and raises `CallbackError` (RS4008) after the solve. Hooks of the `separate`, `branching_score`, `select_node` kind force presolve off; keep new hooks in that table (`Has`).
- A callback test needs a tree: use `hooks_common.TREE_OPTIONS`/`CUT_OPTIONS`; under the default profile small models are solved at the root and no hook fires.
- The crate has its own `build.rs` and `pyproject.toml`; abi3 wheels target Python 3.11 and newer (`Cargo.toml`).
- Release wheels are built by `tools/build_wheel.sh`.
- Further reading: `docs/design/python-api.md`.
- Run `cargo fmt --all` before committing; CI checks it.
