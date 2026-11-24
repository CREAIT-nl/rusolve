# Python API: `rusolve.Model`

Implemented; measured results at the end. The design rules come from measurements of highspy,
SciPy, pyoptinterface, linopy and PuLP building the same models (table below).

## Goal

One Python interface that is **fast** and **simple**, built for rusolve, not a copy of
another library. Users define models with plain Python data (ints, dicts, lists, tuples),
numpy arrays, scipy.sparse matrices, or small expression objects. Rust reads all of them
directly. Compatibility with other libraries (SciPy `milp`/`linprog`, highspy `Highs`) is
provided as thin layers on top. They never shape the core.

What the measurements showed, and the rule each finding gives:

| Finding | Rule |
|---|---|
| highspy's expression API costs 1.5 µs per term: every term is a Python object, and every row gets an argsort and a solver call. | No Python object per term unless the user creates one. Rows go into Rust buffers. |
| HiGHS `addRow` costs ~1.1 µs per row because it validates on every call. pyoptinterface is capped there. | Adding a row only appends to a buffer. Checking and merging duplicates happen once, at `solve()`. |
| SciPy `milp` has ~0.2 ms fixed cost per call; highspy 0.06 ms. On a 1 ms production request that is 20 %. | Fixed cost for a tiny model, round trip, ≤ 20 µs. |
| PuLP and linopy keep their own model and copy it at solve time (44 % of PuLP's time). | One model, in Rust, from the first `add`. |
| A wrapper built a dense `rows × cols` matrix from dict rows. | No dense path anywhere; every path is O(nnz). |
| SciPy rebuilds the solver on every call: no warm start. | The model persists; changes are applied in place and re-solved from the previous basis. |
| The current `rusolve-py` `addRow` loop is quadratic: `assess_matrix` allocates one entry per column on every call. | Fixed in the port, and avoided by design. |

## Performance targets (acceptance; measured by `crates/rusolve-py/bench/bench_api.py`)

| ID | What | Target | For scale |
|---|---|---|---|
| P1 | Tiny model (1 var, 1 row): `Model()`, add, `solve()`, read objective; median | ≤ 20 µs | SciPy `milp` 210 µs, highspy 60 µs |
| P2 | `add_constrs` with dict rows: library time per term, excluding the user's dict construction | ≤ 30 ns | pure-Python dict→CSR 45 ns |
| P3 | `add_matrix(scipy csr)`, from the call to a built `HighsLp` | ≤ 5 ns/nnz | highspy `addRows` 8 ns |
| P4 | Expression API on the 600k-nnz facility model (`bench_api.py fl_expr`) | ≤ 300 ns/nnz | highspy 1500, pyoptinterface 370 |
| P5 | Scaling: ns/nnz at 3M nnz within 1.5× of the value at 600k, for every path | no superlinear path | old `addRow` loop 10 µs/nnz at 3M |
| P6 | `sol.x` for 1M variables | ≤ 2 ms | — |

## Package layout

```
crates/rusolve-model/            pure Rust, no pyo3: ModelBuilder (buffers → HighsLp), Rust tests
crates/rusolve-py/
  src/lib.rs                     module init; the legacy highspy-shaped `Highs` class (unchanged)
  src/ingest.rs                  reading Python values fast: scalars/sequences/buffers → f64/i32,
                                 rows (dict, (idx, val) pair, LinExpr, Var, Constraint) → builder
  src/expr.rs                    Var, Vars, LinExpr, Constraint, quicksum
  src/model.rs                   Model: variables, constraints, objective, solve, changes
  src/solution.rs                Solution
  python/rusolve/__init__.py     re-exports from rusolve._rusolve, plus inf
  python/rusolve/scipy.py        milp() / linprog() with SciPy's signatures (compat layer)
  tests/                         pytest
  bench/bench_api.py             P1–P6 and the comparison table
```

The extension module is `rusolve._rusolve` (maturin mixed layout, `python-source = "python"`).
`import rusolve` re-exports everything, so `rusolve.Highs` keeps working. abi3 minimum becomes
Python 3.11, because `PyBuffer` (zero-copy reads of numpy/array buffers) is in the limited API
only from 3.11. Python 3.10 reached end of life in Oct 2026. No Rust dependency on numpy: inputs are
read through the buffer protocol, and outputs become numpy arrays through
`numpy.frombuffer(bytes)` when numpy is importable, or lists otherwise.

## The API

```python
import rusolve as rs

m = rs.Model(sense="min", name=None)

# --- variables -------------------------------------------------------------------------------
x = m.add_vars(n, lb=0.0, ub=rs.inf, cost=0.0, integer=False, binary=False, names=None)  # -> Vars
v = m.add_var(lb=0.0, ub=rs.inf, cost=0.0, integer=False, binary=False, name=None)       # -> Var
#   lb, ub, cost: scalar, or a sequence / buffer of length n (numpy arrays read without a copy)
#   integer: bool, or a sequence of bool of length n
#   binary=True: integer, with default bounds [0, 1] (an explicit lb/ub still applies, e.g. ub=1 or lb=1)
#   names: None, a list of n str, or a str prefix p giving "p[0]" ... "p[n-1]"

# --- constraints -----------------------------------------------------------------------------
# A "row" is any of:
#   dict {var: coef}            var = Var or int (column index); the fastest scalar form
#   (indices, values) pair      each a sequence or buffer (int32/int64 and float64 buffers zero-copy)
#   LinExpr or Var              from the expression layer
r  = m.add_constr(row, "<=", 3.0, name=None)        # -> int row index; sense "<=" ">=" "==" (or "L" "G" "E")
r  = m.add_constr(x[0] + 2 * x[1] <= 3)             # a Constraint object; then sense/rhs not given
r  = m.add_constr(row, lo=1.0, hi=4.0)              # ranged row
rr = m.add_constrs(rows, "<=", rhs, names=None)     # rows: iterable of rows or Constraints -> range
#   sense: one str or a sequence of str; rhs / lo / hi: scalar or sequence / buffer
rr = m.add_matrix(A, "<=", b, vars=None, names=None)  # -> range
rr = m.add_matrix(A, lo=l, hi=u)                     # ranged rows, as add_constr
#   A: scipy.sparse (csr/csc/coo, read by attributes, no scipy import), or a 2-D float64 buffer
#      (numpy dense; nonzeros only are kept); vars: Vars / sequence mapping A's columns to
#      variables (default: column j is variable j)

# --- objective -------------------------------------------------------------------------------
m.minimize(obj); m.maximize(obj); m.set_objective(obj, sense=None)
m.set_quadratic_objective(Q)     # objective = c'x + (1/2) x'Qx; Q: matrix, {(i, j): v} or
                                 # rows=/cols=/vals= triplets; None clears it. PSD (NSD for max)
                                 # is checked at solve(); a Hessian-free model never enters this path.
#   obj: dict {var: coef}, LinExpr, Var, (indices, values), or a dense sequence / buffer of
#        length num_vars. Replaces the whole objective (other costs become 0); a LinExpr's
#        constant becomes the objective offset. add_vars(cost=...) also sets costs.

# --- general constraints (model/gencon.rs): each returns the object's name -------------------
m.add_indicator(z, 0, x + y <= 5, name=None)   # binary z; value 0/1; one-sided or equality row
m.add_sos1(vars, weights=None); m.add_sos2(vars, weights=None)
m.add_pwl(x, y, xs, ys)                        # y = f(x), xs non-decreasing, x within [xs[0], xs[-1]]
m.add_abs(res, arg); m.add_max(res, args, constant=None); m.add_min(res, args, constant=None)
m.add_and(res, args); m.add_or(res, args); m.add_not(res, arg)     # binaries
#   Stored as the typed rusolve-io objects; solve() linearizes a copy (rusolve-cli gc_solve) and
#   reports `sol.general_constraints`; check() and write() use the objects. Malformed: RS4009;
#   explain/iis/ranging/feasrelax/diagnose, solve(callback=) and a Hessian refuse: RS4010.

# --- scheduling and combinatorial objects (model/cp.rs; docs/design/cp-model.md) -------------
iv = m.add_interval(start, size, end=None, presence=None, name=None)   # -> Interval handle
#   start/end: integer Vars with finite bounds (end=None makes `<name>_end`); size: int or Var >= 0;
#   presence: binary Var, 0 = absent (takes part in nothing).
m.add_no_overlap(intervals, name=None); m.add_cumulative(intervals, demands, capacity, name=None)
m.add_all_different(vars, name=None); m.add_element(index, values, target, name=None)
m.add_table(vars, tuples, negated=False, name=None)       # each returns the object's name
m.num_cp_objects
#   Stored as the typed rusolve-io objects (`CpModel`); solve() goes through
#   `rusolve_cli::cp_solve::solve_with` (engine hook; today the MIP linearization), reports
#   `sol.cp` (None for a plain model) and the envelope `cp` object; check() uses the objects.
#   Malformed: RS4011; the fallback cannot take a valid object: RS4013; explain/iis/ranging/
#   feasrelax/diagnose/write, solve(callback=) and a Hessian refuse: RS4012.

# --- solve -----------------------------------------------------------------------------------
sol = m.solve(time_limit=None, mip_gap=None, threads=None, presolve=None, seed=None,
              verbose=False, options=None)          # options: {raw_option_name: value}
#   presolve: bool (True -> "on", False -> "off") or a raw str ("choose"). Raw HiGHS option
#   names pass through `options` (e.g. mip_max_nodes, simplex_iteration_limit, solver).
#   The GIL is released during the solve. Defaults = HighsOptions::recommended() with
#   output off; `verbose=True` turns the HiGHS log on.

# --- solution --------------------------------------------------------------------------------
sol.status          # str, the closed status enum of rusolve.result/1: "Optimal", "Infeasible", "Unbounded",
                    # "InfeasibleOrUnbounded", "TimeLimit", "IterationLimit", "NodeLimit",
                    # "SolutionLimit", "Interrupted", "Error", ... (same names as the envelope)
sol.ok              # True if a primal solution is available (x is not None)
sol.objective, sol.bound, sol.gap        # float or None
sol.x               # numpy float64 array (list if no numpy), None if no solution
sol.row_activity, sol.duals, sol.reduced_costs   # arrays or None (duals: LP only)
sol.nodes, sol.iterations, sol.time
sol[v]              # value of a Var / int index / LinExpr (including its constant) / dict row
sol.values(mapping, above=None)          # {key: value} for a {key: var} mapping; drop value <= above
sol.selected(mapping, threshold=0.5)     # [key ...] with value > threshold (binary picks)

# --- changes (persistent model; re-solve reuses the basis) -----------------------------------
m.set_bounds(vars, lb, ub); m.set_cost(vars, cost); m.set_integer(vars, flag)
m.set_row_bounds(rows, lo, hi); m.set_coef(row, var, value)
m.num_vars, m.num_constrs, m.num_nonzeros
m.var_name(v), m.constr_name(r)
m.write(path)       # .mps (stretch goal)
```

### Expression layer

- `Var` stores only its column index (and a model id, to catch mixing two models). It supports
  `__index__`, `__int__`, `__hash__` (= hash of the index), `+ - *` with numbers, Vars and
  LinExprs, `/` by a number, unary `-`, and `<= >= ==`, which return a `Constraint`.
  `Var == Var` returns a Constraint whose `__bool__` is index equality, so Vars still work as
  dict keys and in `in` tests. Every other Constraint raises in `__bool__`, with a message
  pointing to `between`.
- `Vars` is the block `add_vars` returns: `start`, `stop`, `len`, `indices` (a `range`), and
  iteration, int and slice indexing (an int gives a `Var`, a slice gives a `Vars`), `sum()`, and
  `dot(coefs)` (sequence or buffer). It is a compact range: no per-element objects until indexed.
- `LinExpr`: a term list plus a constant. `+=`/`-=` work in place. Duplicate terms are kept
  and merged in Rust when the row is built. `expr.between(lo, hi)` gives a ranged Constraint.
- `rs.quicksum(iterable)`: sums Vars, LinExprs, numbers and `(coef, var)` pairs into one
  LinExpr without intermediate objects.
- `Constraint`: holds a LinExpr (the constant moved to the bounds), `lo` and `hi`.

### Semantics

- Indices are 0-based column numbers. A Var from another Model is an error (ValueError).
- Duplicate entries in a row are summed, and exact-zero coefficients are dropped, when the row is
  built. NaN anywhere is a ValueError naming the row (by name if it has one) and the column.
  An infinite coefficient is a ValueError. `±inf` bounds mean unbounded.
- Names: if no name is ever given, the HighsLp gets no names. Otherwise unnamed columns and rows
  get HiGHS's defaults `c<j>` / `r<i>` (HConst.h `kHighsMinimalColNamePrefix`) at build time.
- Errors are Python exceptions (`ValueError`, `TypeError`, `IndexError`) with the offending
  row and column, raised by the `add_*` / `set_*` call that brings the bad value (a failed
  call adds nothing). Nothing fails silently or is deferred to `solve()`.
- `m._sync()` (internal, for benchmarks and tests): hand everything pending to the solver
  without solving, so build cost can be timed apart from solve cost.
- Persistence: until the first `solve()` the model lives in the `ModelBuilder`. The first
  `solve()` passes it to the ported `Highs` (`pass_model`). After that, additions are buffered and
  sent in one `add_cols` / `add_rows` per kind at the next `solve()`. Changes to existing
  entities go straight to the `Highs` change API, which keeps the basis for a warm start.
  The objective is pushed with `change_cols_cost_interval`, the sense with
  `change_objective_sense`, and the offset with `change_objective_offset`.

### Compatibility layers (on top, optional)

- `rusolve.scipy.milp(c, *, integrality=None, bounds=None, constraints=None, options=None)`
  and `rusolve.scipy.linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None,
  method="highs", options=None)`: SciPy's signatures and result fields (`x`, `fun`, `status`,
  `success`, `message`, `mip_gap`, `mip_node_count`, `mip_dual_bound`; `linprog` also
  `ineqlin`/`eqlin` marginals for LPs), with SciPy's status codes (0 optimal, 1 limit,
  2 infeasible, 3 unbounded, 4 other). They duck-type `Bounds` and `LinearConstraint`
  (`.lb`/`.ub`/`.A`), so scipy is never imported by rusolve itself. Switching from SciPy is a change of
  import.
- `rusolve.Highs`: the existing highspy-shaped class, unchanged.

## Tests

`cargo fmt --all --check`, `cargo clippy -p rusolve-model -p rusolve-py -- -D warnings`,
`cargo test -p rusolve-model`, `pytest crates/rusolve-py/tests`, and `bench_api.py` for P1–P6.
The `addRow` fix in the port (`crates/rusolve-ekk/src/util/highs_matrix_utils.rs`) also passes the
parity check `tools/bench/identity.sh` (`IDENTICAL 20 of 20`).

## Results (`bench_api.py` full run)

Apple M-series, Python 3.12, release build. Other benchmark jobs were sharing the machine (load ~6).

| ID | Measured | Target | For scale |
|---|---|---|---|
| P1 | 8.8 µs (`rusolve.scipy.milp` 14.3 µs) | ≤ 20 µs | scipy.optimize.milp 151 µs, highspy 47 µs |
| P2 | 6.9 ns/term after the raw `PyDict_Next` reader (17.0 before; (idx, val) array pairs 5.4) | ≤ 30 | building the dicts themselves: 51 ns/term |
| P3 | 9.1 ns/nnz (add_matrix 2.9 + hand-over to the solver 7.0) | ≤ 5 | **miss**, see below |
| P4 | 61 ns/nnz | ≤ 300 | highspy expressions 1574, pyoptinterface 471 |
| P5 | ratio 1.07 (expr), 1.12 (dicts), 1.08 (CSR) | ≤ 1.5 | — |
| P6 | 0.12 ms | ≤ 2 ms | — |

P3 is measured from the `add_matrix` call to the model inside the solver (`_sync`), which is stricter
than the target as written ("to a built HighsLp"). The builder's own part is 4.7 ns/nnz in Rust
(`crates/rusolve-model/benches/builder.rs`: CSR append 2.0 + take/into_lp 2.6). The rest is
`Highs::pass_model` re-checking the matrix the builder already checked. Assembly included, highspy
`addRows` costs 8.2 ns/nnz on the same model, and rusolve's numpy-CSR path 10.3.
TODO(perf): give the port a pre-validated pass-model entry point (skip `assess_matrix` for a
builder-made LP) to reach ≤ 5. The win is ~4 ns/nnz, i.e. ~2.5 ms on a 600k-nnz model whose LP
solve takes 4.2 s, so this is not urgent.

Deviations of the implementation from the spec above:
- `Vars` takes strided slices too (`x[c::C]`).
- `Var != Var` is a plain bool; strict `<` / `>` raise TypeError.
- `Var == Var` across two models is False, not an error, so mixed dict lookups don't raise.
- `quicksum` takes `(var, coef)` pairs too.
- `sol[Vars]` is an array.
- `add_matrix` takes any scipy format (via `.tocsr()`) and lists of lists.
- `Model.var(j)` and `Model.sense` exist.
- `var_name`/`constr_name` give the `c<j>`/`r<i>` defaults only after the first solve or `_sync`.
- Unknown option names raise ValueError. The SciPy layer warns first, then raises: stricter than
  SciPy, which ignores them.
- `write(path)` is not implemented yet.

## Callbacks and `diagnose`

```python
class SubtourElimination:                      # any object; each method is one hook
    def lazy_constraints(self, ctx: rs.CallbackContext) -> list[tuple[list[int], list[float], float]]: ...

sol = m.solve(callback=SubtourElimination(n), threads=1)     # raises rs.CallbackError (RS4008) if a hook raises
m.callback_stats        # calls per hook, errors, cuts accepted, the lazy report
m.diagnose()            # == `rusolve diagnose <model.mps> --format json`
```

Hooks, what each returns and the decisions (presolve off for the steering hooks, lazy rows per
solve, failure policy, interrupt) are in `docs/design/hooks.md`, section Python. `Model.diagnose`
(`rusolve_cli::api::diagnose`, a thin adapter over the same `diagnose_core` the CLI runs) returns
the `rusolve.diagnose/1` object; a MIP is diagnosed through its LP relaxation, a quadratic
objective raises RS4001 like the other diagnostics; with it every CLI diagnostic has a Python call with
a test (`tests/test_callbacks.py::test_diagnose_equals_the_cli` compares the whole object with
the CLI's JSON). Stubs: `CallbackContext`, `Model.solve(callback=, lazy_rounds=)`,
`Model.callback_stats`, `Model.diagnose` in `_rusolve.pyi`; `CallbackError` in `errors.py`; the
docstring examples type-check under `mypy --strict` (`test_typing.py`).
