# Scheduling and combinatorial objects (CP model)

A model can state scheduling and combinatorial structure natively: intervals, no-overlap,
cumulative resources, all-different, element and table constraints. The default build solves them
through a MIP linearization and checks the answer on the objects you wrote. With
`mip_engine=cpsat` the same objects go to OR-Tools CP-SAT instead, through one hook (below;
`docs/design/cpsat-engine.md`); the answer is checked on your objects either way.

## The objects

All of them range over integer variables with finite bounds. Indices are columns of the model.

| object | meaning |
|---|---|
| interval | `start + size = end` while the interval is present. `size` is a constant or a column. An optional `presence` binary: when it is 0 the interval is absent, takes part in no constraint, and its columns are free within their bounds. |
| no-overlap | no two present intervals of the set overlap. An interval of size 0 overlaps nothing. |
| cumulative | at every time `t`, the demands of the present intervals with `start <= t < end` sum to at most the capacity. |
| all-different | the columns take pairwise different values. |
| element | `target = values[index]`, `index` 0-based, `values` constants or columns. |
| table | the tuple of values of the columns is one of the rows of `tuples` (none of them when `negated`). |

The semantics follow the Global Constraint Catalogue (Beldiceanu, Carlsson and Rampon, SICS
T2005:08) and the public OR-Tools CP-SAT documentation. The objects live in
`crates/rusolve-io/src/cp.rs` (`CpModel`) and stay in the model as you wrote them; nothing derived
from them is ever reported as yours.

## Surfaces

- Python: `Model.add_interval(start, size, end=None, presence=None)` returns an `Interval`;
  `add_no_overlap`, `add_cumulative`, `add_all_different`, `add_element`, `add_table` (see
  [python-api.md](python-api.md)). `Solution.cp` holds the report (`None` for a plain model).
  Python needs the `rusolve` wheel; the `install.sh` route has no Python.
- MCP `build_model`: a `cp` section with `intervals`, `no_overlap`, `cumulative`,
  `all_different`, `element`, `table`. This is the surface of an agent that installed only the
  binary. Variables and intervals are named by string. A bad entry is
  an error that names its JSON path, for example `cp.no_overlap[0].intervals[2]`. `solve`,
  `get_solution`, `check_solution` and `free_model` work on such a handle; the other tools refuse
  it (`RS4012`).
- No file syntax. MPS and LP have no grammar for these objects, so `rusolve solve model.mps`
  never sees them, `Model.write()` refuses a model that has them, and `rusolve check` on a file
  has nothing to check.

### The `cp` section of `build_model`

Entries by section (`mcp_build.rs`, `parse_cp`):

| section | entry |
|---|---|
| `intervals` | `{name, start, size, end, presence?}`: `start` and `end` are names of variables declared in `variables`; `size` an integer or a variable name; `presence` an optional binary variable |
| `no_overlap` | `{name, intervals}`: interval names |
| `cumulative` | `{name, intervals, demands, capacity}`: one demand per interval, each and the capacity an integer or a variable name |
| `all_different` | `{name, variables}` |
| `element` | `{name, index, values, target}`: `index` and `target` variables, `values` integers or variable names |
| `table` | `{name, variables, tuples, negated?}` |

Rules, each checked when the model is built:

- Every variable a `cp` entry names must be declared in `variables` with `type` `integer` (or
  `binary` for a presence), and with both `lower` and `upper` finite. A variable's bounds default
  to `0` and infinity, so write `upper`. The interval does not create its `end`: declare the
  variable and name it.
- `start`, `size`, `end` and `name` of an interval are required; `presence` is the only optional
  field. A missing field is `RS3001` naming `cp.intervals[i].end` and so on; a variable that is
  continuous or has an infinite bound is `RS4011` (the same text the Python call raises).
- A number for `size`, a demand, a value or the capacity must be an integer: `3`, not `3.0`
  (`must be an integer or a variable name`). A size or demand, as a number or through a
  variable's lower bound, must not be negative.
- The reply has the handle's `model_id` (`m1`, `m2`, ... per server process) in
  `structuredContent` and in the JSON text of `content[0].text`.

The same model, two ways (two jobs on one machine, optimum 5):

```python
m = rs.Model()
s1, e1 = m.add_var(ub=9, integer=True), m.add_var(ub=9, integer=True)
s2, e2 = m.add_var(ub=9, integer=True), m.add_var(ub=9, integer=True)
mk = m.add_var(ub=9, integer=True, cost=1)
a, b = m.add_interval(s1, 2, end=e1), m.add_interval(s2, 3, end=e2)
m.add_constr(mk >= a.end)
m.add_constr(mk >= b.end)
m.add_no_overlap([a, b])
```

```json
{"variables": [{"name": "s1", "type": "integer", "lower": 0, "upper": 9},
               {"name": "e1", "type": "integer", "lower": 0, "upper": 9},
               {"name": "s2", "type": "integer", "lower": 0, "upper": 9},
               {"name": "e2", "type": "integer", "lower": 0, "upper": 9},
               {"name": "mk", "type": "integer", "lower": 0, "upper": 9, "objective": 1}],
 "rows": [{"name": "mk1", "terms": {"mk": 1, "e1": -1}, "sense": ">=", "rhs": 0},
          {"name": "mk2", "terms": {"mk": 1, "e2": -1}, "sense": ">=", "rhs": 0}],
 "cp": {"intervals": [{"name": "a", "start": "s1", "size": 2, "end": "e1"},
                      {"name": "b", "start": "s2", "size": 3, "end": "e2"}],
        "no_overlap": [{"name": "one", "intervals": ["a", "b"]}]}}
```

A job-shop recipe with both forms is in [`docs/skills/scheduling`](../skills/scheduling/SKILL.md).

## The checker

`CpModel::check` evaluates every object on a point and reports violations per object (kind,
name, amount, detail). It reads the user's columns, rounds an integer column that is within the
integrality tolerance and counts a column that is farther from an integer as a violation of the
object that uses it. `rusolve check`, the Python `Model.check`, the MCP `check_solution` and the
result envelope all use it (`summary.cp`, findings of kind `interval`, `no_overlap`, `cumulative`,
`all_different`, `element`, `table`; the envelope's `cp.verified`).

A solve never reports an answer as feasible when the checker fails. The linearization can leave a
binary a hair off its value and let a big-M row hold some slack; the solve fixes every integer
column at its rounded value and solves once more. An answer that still fails becomes "no
solution" with status `Unknown`.

## The MIP linearization

`crates/rusolve-io/src/cp_lin.rs`. It appends binaries after your columns and rows after your
rows, so cutting the answer back to your columns is exact. The route of each object is in the
envelope (`cp.routes`).

| object | formulation | source |
|---|---|---|
| interval | `start + size = end` as two big-M rows relaxed by the presence column | Williams, *Model Building in Mathematical Programming*, 5th ed., ch. 9; Bonami, Lodi, Tramontani and Wiese, Math. Programming 151 (2015) |
| no-overlap | one disjunctive binary per pair of intervals, big-M from the bounds; a pair the bounds already order is dropped | Manne, Oper. Res. 8 (1960); Applegate and Cook, ORSA J. Comput. 3 (1991) |
| cumulative | time-indexed: one binary per (interval, start time), one capacity row per time point | Pritsker, Watters and Wolfe, Management Science 16 (1969) |
| all-different | direct encoding: one binary per (variable, value) | Walsh, "SAT v CSP", CP 2000 |
| element | one selection binary per admissible index | Hooker, *Integrated Methods for Optimization*, 2nd ed. (2012) |
| table | one binary per allowed tuple; a negated table forbids each tuple with a disjunction over the columns | the same |

Limits (`CpLinOptions`): a derived big-M above `1e6`, a cumulative horizon above 1000 time points
(the horizon is read from the start and end bounds), more than 100000 enumerated binaries for one
object. The fallback also refuses a variable size or demand in a cumulative and a no-overlap
interval whose size can be 0. Each is a typed `RS4013` error that names the object and says the
CP-SAT engine takes the model, so tighten the bounds, use constants, or select that engine
(`mip_engine=cpsat`, installing OR-Tools if `rusolve-cli version` says it is missing). The error's
help line says which of the two applies.

## Diagnostic codes

| code | meaning |
|---|---|
| `RS4011` | an invalid object: a variable that is not integer or has an infinite bound, a presence that is not binary, a negative size, demand or capacity, an interval index out of range or listed twice, a tuple of the wrong length |
| `RS4012` | a call that would ignore the objects: `explain`, `iis`, `ranging`, `feasrelax`, `diagnose`, `lint(with_solve=True)`, a `solve` callback, a quadratic objective, `write`, and the MCP tools `modify`, `explain`, `lint`, `ranging`, `feasrelax` |
| `RS4013` | the fallback cannot take a valid object (the limits above) |

`rusolve --explain RS4011` (and `RS4012`, `RS4013`) prints the cases.

## Router and the engine hook

The router class is `cp` and the engine is `cp-linearization`. A model that has both CP objects
and general constraints is class `cp_mixed` and is refused. The one place that picks the engine is
`crates/rusolve-cli/src/cp_solve.rs`:

```rust
pub fn select_engine(cp: &CpModel) -> CpEngine
pub fn solve_with(lp: HighsLp, cp: CpModel, options: HighsOptions,
                  run: impl FnOnce(&mut Highs) -> HighsStatus) -> Result<CpSolved, CpError>
```

Python and MCP both call `solve_with`. An optional engine adds a `CpEngine` variant, an arm in
`select_engine` and an arm in `solve_with` that returns a `CpSolved` whose solution is in your
columns. The report comes from `cp_solve::report`, so the answer is checked on your objects
whatever engine produced it.

## Modelling advice

Use the native objects when a model has resources or sequencing and you want an answer you can
check without trusting a big-M. The fallback is a MIP, so a big-M you would write by hand is no
slower; the gain is correctness of the statement and a checker that speaks your objects. See
[`docs/skills/scheduling`](../skills/scheduling/SKILL.md) for when big-M and indicators are still
the right tool.
