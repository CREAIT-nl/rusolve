# rusolve-py: Python bindings

```
maturin develop --release -m crates/rusolve-py/Cargo.toml     # see "Build and test" below
```

`import rusolve` gives three things: `rusolve.Model`, a fast modelling API built for rusolve
(front door, below); `rusolve.scipy`, drop-in `milp` and `linprog`; and `rusolve.Highs`, the
highspy-compatible class. The design and its measurements are in
`docs/design/python-api.md`; `bench/compare_libs.py` and `bench/compare_floors.py` measure the
build paths against other open-source Python modelling layers.

## Quickstart

```python
import rusolve as rs

m = rs.Model()                                  # minimise by default; Model(sense="max")
x = m.add_vars(3, binary=True)                  # a block of 3 binary variables
m.add_constr(x[0] + x[1] + x[2] <= 2)           # expressions, or plain data (below)
m.maximize(5 * x[0] + 4 * x[1] + 3 * x[2])
sol = m.solve(time_limit=60)

sol.status        # "Optimal"
sol.objective     # 9.0
sol.x             # numpy array, one value per variable
sol[x[0]]         # value of one variable or expression
```

Plain data goes in without per-term Python objects, and is the fast path for big models:

```python
y = m.add_vars(n, lb=0, ub=ub_array, cost=cost_array)        # scalars, lists or numpy arrays
m.add_constr({y[0]: 1.0, y[1]: 2.0}, "<=", 10.0)             # dict row ({var or index: coef})
m.add_constr((indices, values), ">=", 1.0)                   # index and value arrays
m.add_constrs(rows, "<=", rhs_vector)                        # many rows at once
m.add_matrix(A, "<=", b)                                     # scipy.sparse csr/csc/coo or dense
m.minimize(cost_array)                                       # dense objective
```

Adding only appends to a buffer; checking (NaN, infinite coefficients, foreign variables) and
merging of duplicate entries happen once, at `solve()`, and errors name the row and column.

After a solve the model stays alive. Changes are applied in place and the next `solve()` restarts
from the previous basis:

```python
m.set_bounds(x, 0, 5); m.set_cost(x, new_cost); m.set_row_bounds(rows, lo, hi); m.set_coef(r, v, 2.0)
sol = m.solve()
```

Helpers for reading solutions: `sol.values({"a": va, "b": vb}, above=1e-6)` gives a dict,
`sol.selected(mapping)` the keys of binaries that are on. `sol.duals`, `sol.reduced_costs`
and `sol.row_activity` are available for LPs. `sol.status` is one of `"Optimal"`,
`"Infeasible"`, `"Unbounded"`, `"InfeasibleOrUnbounded"`, `"TimeLimit"`, ... and `sol.ok` says
whether a primal solution exists.

Measured, not claimed: `python crates/rusolve-py/bench/bench_api.py` runs the acceptance
targets P1 to P6 (tiny-model round trip, dict-row and matrix ingest, the facility-location
model through the expression API against highspy and pyoptinterface, scaling to 3M nonzeros,
`sol.x` for 1M variables) and prints PASS/FAIL; `--quick` runs it in seconds.

## Diagnostics: what an agent asks after a solve

Each CLI diagnostic is a `Model` method. They return the same `rusolve.<command>/1` object that
`--format json` prints, parsed into a dict, so a Python agent and a shell agent read identical
fields (schemas in `docs/schema/`):

```python
import rusolve as rs

m = rs.Model(sense="max", name="knapsack")
x = m.add_vars(4, binary=True, names="x")
m.add_constr(3 * x[0] + 4 * x[1] + 5 * x[2] + 6 * x[3] <= 10, name="weight")
m.maximize(4 * x[0] + 5 * x[1] + 6 * x[2] + 8 * x[3])

report = m.lint()                       # rusolve lint: big-M, coefficient range, free columns
sol = m.solve(time_limit=10)
envelope = sol.to_dict()                # rusolve.result/1, what `rusolve mip --format json` prints
text = sol.to_json(indent=2)            # the same as a JSON string (include_solution=True)
verdict = m.check(sol)["verdict"]       # rusolve check: bounds, rows, integrality, objective
ranges = m.ranging()                    # rusolve ranging: cost / bound / rhs ranges (solves a copy)
m.write("knapsack.mps")                 # .mps, .lp, either with .gz
```

For an infeasible model, `explain()` (alias `iis()`) names the conflicting rows and bounds and
`feasrelax()` finds the cheapest relaxation:

```python
import rusolve as rs

bad = rs.Model(name="bad")
v = bad.add_var(name="v")
bad.add_constr(v >= 3, name="lo")
bad.add_constr(v <= 1, name="hi")
iis = bad.explain()                     # rusolve explain: the irreducible infeasible subsystem
fix = bad.feasrelax(mode="min-count")   # rusolve feasrelax: fewest relaxed rows or bounds
```

`check(solution)` takes a `Solution`, a sequence or array of values, or a dict keyed by `Var`
or column name. `check`, `ranging`, `feasrelax`, `explain` (and `iis`) and `lint(with_solve=True)` refuse a model with a quadratic objective, and `write`
needs `.lp` for one (the MPS writer has no quadratic section). All of these read the model as
the solver holds it and never change it.

The package ships `py.typed` and `.pyi` stubs for the whole public API; `mypy --strict` accepts
the examples in this file (`tests/test_typing.py` checks it) and every `>>>` example in the
docstrings runs as a doctest.

## Typed errors and source lines

Every error `rusolve` raises on purpose is a `rusolve.RusolveError` carrying `.diagnostics`, the
`diagnostics` array of the `rusolve.diagnostics/1` document (`docs/design/diagnostics.md`). Each
entry has a registered `code` (`rusolve --explain RS4001` or `docs/agent/diagnostics.md` says what it
means), `message`, `help` lines, and a `location`: the line of your Python that made the call.

```python
import rusolve as rs

m = rs.Model()
x = m.add_var(name="x")
m.add_constr(x >= 3, name="lo")
m.add_constr(x <= 1, name="hi")
try:
    m.solve(presolve=False).raise_for_status()   # raises for anything but Optimal
except rs.InfeasibleError as e:
    code = e.code                                # "RS4004"
    where = e.location                           # {"kind": "python", "path": ..., "line": ..., "source_line": ...}
    fixes = e.help                               # ["`Model.explain()` returns an irreducible ..."]
```

| class | codes | raised when |
|---|---|---|
| `ModelError` (also `ValueError`) | RS4001, RS4003, RS4009, RS4010, RS4011, RS4012, RS4013 | a model-building call gets a bad argument, or mixes two Models; a general constraint is malformed (RS4009) or a call cannot honour one (RS4010); a scheduling object is malformed (RS4011), a call cannot honour one (RS4012) or the MIP fallback cannot take it (RS4013) |
| `ModelIndexError` (also `IndexError`) | RS4002 | a row or column index is out of range |
| `InfeasibleError`, `UnboundedError`, `LimitError` | RS4004, RS4005, RS4006 | `Solution.raise_for_status()` on that status |
| `NoSolutionError` (also `ValueError`) | RS4007 | a value is read from a Solution without one |
| `CallbackError` | RS4008 | a hook of `Model.solve(callback=...)` raised; `__cause__` is the original, `.solution` the partial result |

`Model.solve()` still returns the status instead of raising; `raise_for_status()` is the opt-in.

`Model.capture_source = True` records the Python line that adds each constraint (off by default).
`Model.source_of(row)` returns `{"path", "line", "source_line"}`, and `explain()`/`iis()` put it on
every row of the infeasible subsystem as `member["source"]`, so an agent sees which line of its own
code to fix. Cost, measured by `tests/test_errors.py::test_capture_cost_is_measured_and_bounded`
(Python 3.14, `add_constr` of a two-term row, min of 3 runs of 20000): about 0.4 us per call off,
1.35 us on (about 0.9 us for the frame walk).

## Quadratic objective (convex QP)

```python
m.set_quadratic_objective(Q)                       # objective = c'x + (1/2) x'Qx
m.set_quadratic_objective({(0, 0): 2.0, (0, 1): 0.5, (1, 1): 1.0})
m.set_quadratic_objective(rows=i, cols=j, vals=v)  # triplets
m.set_quadratic_objective(None)                    # back to the LP
```

`Q` is a scipy.sparse matrix, a 2-D array or a list of lists; duplicates are summed and an
asymmetric `Q` becomes `(Q + Q')/2`, which has the same objective. Quadratic *constraints* are
not supported. Models without a Hessian take the unchanged LP/MIP path.

## General constraints

Indicator, SOS1/SOS2, piecewise-linear, abs/min/max and and/or/not are methods of `Model`; each
returns the object's name. `solve()` linearizes a copy of the model the way `rusolve mip` does for
the same objects in a file (big-M from activity bounds, incremental piecewise-linear, an LP
epigraph when the objective pushes the right way), and `check()`, `write()` (LP; MPS holds
indicator and SOS only) and the result envelope report on your objects, not on generated rows.

```python
import rusolve as rs
m = rs.Model()
z = m.add_var(binary=True, cost=7); x = m.add_var(ub=10, cost=2); y = m.add_var(ub=10, cost=1)
m.add_constr(x + y >= 3)
m.add_indicator(z, 0, x >= 4)                    # z = 0 -> x >= 4
# also add_sos1, add_sos2, add_pwl, add_abs, add_max, add_min, add_and, add_or, add_not
sol = m.solve()
gc = sol.general_constraints                     # None for a model without general constraints
assert gc is not None
gc["verified"]                                   # the checker's verdict on the objects
```

A malformed object raises `ModelError` with code `RS4009` (non-binary indicator variable, unsorted
or short breakpoints, mismatched lengths); an unknown variable raises `ModelIndexError`. The
diagnostics that would ignore the objects (`explain`, `iis`, `ranging`, `feasrelax`, `diagnose`), a
`solve` callback and a quadratic objective refuse such a model with `RS4010`.

## Scheduling and combinatorial objects

Intervals, no-overlap, cumulative, all-different, element and table are methods of `Model`
(spec: [cp-model.md](../../docs/design/cp-model.md)). They range over integer variables with finite
bounds. `add_interval(start, size, end=None, presence=None)` returns an `Interval`; `size` is an
int or an integer Var, and `end` is made for you when you leave it out. `solve()` linearizes a
copy of the model into a MIP (time-indexed for cumulative, so the horizon is limited), checks the
answer on your objects and cuts it back to your variables; `sol.cp` is the report (`None` for a
model without these objects) and `check()` evaluates the objects too. There is no MPS or LP syntax
for them, so `write()` refuses such a model.

```python
import rusolve as rs
m = rs.Model()
cmax = m.add_var(ub=20, integer=True, cost=1)
jobs = []
for p in (3, 2, 4):
    s = m.add_var(ub=20, integer=True)
    job = m.add_interval(s, p)                       # start s, size p, end made for you
    m.add_constr(cmax >= job.end)
    jobs.append(job)
m.add_no_overlap(jobs)                               # one machine
sol = m.solve()
report = sol.cp
assert report is not None and report["verified"]
assert sol.objective == 9.0 and m.check(sol)["verdict"] == "feasible"
# also add_cumulative(jobs, demands, capacity), add_all_different(vars),
# add_element(index, values, target), add_table(vars, tuples, negated=False)
```

A malformed object raises `ModelError` with code `RS4011` (a variable that is not integer or has
no finite bounds, a negative size, a tuple of the wrong length, an interval listed twice).
`solve()` raises `RS4013` for a valid object the MIP fallback cannot take (a cumulative horizon
above 1000 time points, a variable size or demand in a cumulative, a derived big-M above 1e6).
`explain`, `iis`, `ranging`, `feasrelax`, `diagnose`, `write`, a `solve` callback and a quadratic
objective refuse such a model with `RS4012`.

## SciPy drop-in

```python
from rusolve.scipy import milp, linprog        # was: from scipy.optimize import milp, linprog
```

Same arguments, same `OptimizeResult` fields (`x`, `fun`, `status`, `success`, `message`,
`mip_gap`, `mip_node_count`, `mip_dual_bound`; `linprog` also `slack`, `con`, `ineqlin`,
`eqlin`, `lower`, `upper`) and the same status codes (0 optimal, 1 limit, 2 infeasible,
3 unbounded, 4 other). `Bounds` and `LinearConstraint` objects are read by attribute, so scipy
is not needed to import `rusolve.scipy`; sparse matrices are passed through, never densified.
`integrality` 2 and 3 (semi-continuous, semi-integer) are supported (`Model.set_semi`). Not supported: the non-HiGHS `linprog` methods.
`linprog` ignores `callback` and `x0`, as SciPy's HiGHS methods do.

## highspy-compatible class

`rusolve.Highs` drives the ported HiGHS (`crates/rusolve-ekk`) with a `highspy`-shaped API
(HiGHS 1.16.0 `highs_bindings.cpp` / `highspy/highs.py`): `Highs()`, `readModel`, `passModel`
(a `HighsLp` with the `col_cost_`, `a_matrix_.start_` ... fields, or highspy's flat argument
form), `setOptionValue`/`getOptionValue`, `run`, `getModelStatus`, `getInfo`, `getSolution`,
`getObjectiveValue`, `changeColBounds`/`changeRowBounds`/`changeColCost`/`changeCoeff`,
`addRow`/`addCol`, `getLp`, `writeSolution(path)` (raw style only), and the enums
`HighsStatus`, `HighsModelStatus`, `ObjSense`, `MatrixFormat`, `HighsVarType`.
Inputs are plain lists (no numpy needed). Not covered: callbacks, `HighsModel`/Hessians (use `Model.set_quadratic_objective`),
basis and IIS access.

Defaults equal the CLI's: `Highs()` starts from `HighsOptions::recommended()`;
`Highs(profile="highs-1.15.1")` from the HiGHS 1.15.1 defaults. The CLI's `rusolve mip`
additionally sets `output_flag=false`, `threads=1`, `random_seed=0`; set those through
`setOptionValue` to reproduce its numbers (the tests do).

## Build and test

Run from the repository root:

```
python3 -m venv .venv-py && .venv-py/bin/pip install maturin pytest numpy scipy
.venv-py/bin/maturin develop --release -m crates/rusolve-py/Cargo.toml   # editable install
cargo build --release -p rusolve-cli                                         # for the CLI comparison
.venv-py/bin/pytest crates/rusolve-py/tests
.venv-py/bin/python crates/rusolve-py/bench/bench_api.py --quick   # speed targets P1-P6 of docs/design/python-api.md (add highspy pyoptinterface for the comparison rows)
.venv-py/bin/maturin build --release -m crates/rusolve-py/Cargo.toml     # abi3 wheel, Python >= 3.11
```

`maturin develop` needs the venv active (`source .venv-py/bin/activate`) or `VIRTUAL_ENV`
set. `test_read_mip_matches_cli` skips when `target/release/rusolve-cli` is absent.

## Workspace builds

The crate is a workspace member. `build.rs` adds the extension-module link flags
(`-undefined dynamic_lookup` on macOS), so `cargo build`/`clippy -p rusolve-py` need only a
Python interpreter for pyo3's build configuration (abi3: no linking to libpython). A
clippy run over `-p rusolve-ekk -p rusolve-cli` never builds this crate;
add `-p rusolve-py` to lint it. There are no Rust tests; the tests are the
pytest suite above.
