---
name: rusolve-scheduling
description: Model machine scheduling, job shops and sequencing with disjunctive constraints in rusolve.Model, choose between big-M, tightened big-M and indicator constraints, and find a weak big-M with lint. Use when jobs share a machine or resource, order matters, and a makespan, tardiness or completion-time objective is minimised. Also covers native intervals (Python add_interval, add_no_overlap, add_cumulative, or the cp section of the MCP build_model tool) instead of hand-written big-M, with a job-shop recipe and when to choose mip_engine=cpsat.
---

# Scheduling: big-M, tight big-M or indicators

Jobs on one machine, with release times. For every pair `(i, j)` one of two things holds: `i`
finishes before `j` starts, or the other way round. The disjunction is written with a binary
`y_ij` (1 means `i` goes first):

    s_i + p_i <= s_j + M_ij * (1 - y_ij)        s_j + p_j <= s_i + M_ji * y_ij

The only thing that matters for speed is the size of `M`. It must be at least as large as the
biggest value the left side can exceed the right side when the row is meant to be inactive, and it
should be no larger. A huge `M` makes the relaxation allow fractional `y` at almost no cost and
turns integrality tolerances into real violations: with `M = 1e6` and an integrality tolerance of
1e-6, `y = 1e-6` buys one unit of overlap (Bonami, Lodi, Tramontani and Wiese, *On mathematical
programming with indicator constraints*, Math. Programming 151, 2015; Queyranne and Schulz,
*Polyhedral approaches to machine scheduling*, 1994).

## Formulation with a computed M

```python
import itertools
import rusolve as rs

p = [4, 3, 5, 2, 6]                  # processing times
r = [0, 2, 1, 8, 3]                  # release times
n = len(p)
H = max(r) + sum(p)                  # a valid horizon: some optimal schedule ends by H

def build(big_m=None):
    """big_m None: M_ij from the bounds of the start times (tight). A number: that M everywhere."""
    m = rs.Model()
    s = [m.add_var(lb=r[i], ub=H - p[i], name=f"s{i}") for i in range(n)]
    cmax = m.add_var(lb=0, ub=H, cost=1.0, name="cmax")
    for i in range(n):
        m.add_constr(cmax - s[i] >= p[i], name=f"mk_{i}")
    y = {}
    for i, j in itertools.combinations(range(n), 2):
        y[i, j] = m.add_var(lb=0, ub=1, binary=True, name=f"y_{i}_{j}")   # 1: i before j
        mij = big_m if big_m else (H - p[i]) + p[i] - r[j]      # max(s_i + p_i) - min(s_j)
        mji = big_m if big_m else (H - p[j]) + p[j] - r[i]
        m.add_constr(s[i] + p[i] - s[j] <= mij * (1 - y[i, j]), name=f"ij_{i}_{j}")
        m.add_constr(s[j] + p[j] - s[i] <= mji * y[i, j], name=f"ji_{i}_{j}")
    return m, s, cmax, y

def brute():                         # all n! job orders, started as early as possible
    best = None
    for order in itertools.permutations(range(n)):
        t = 0
        for i in order:
            t = max(t, r[i]) + p[i]
        best = t if best is None else min(best, t)
    return best

for label, bm in (("computed M", None), ("M = 1000", 1000)):
    m, s, cmax, y = build(bm)
    sol = m.solve()
    assert sol.status == "Optimal" and round(sol.objective) == brute() == 20, (label, sol.objective, brute())
    print(label, "optimum", sol.objective, "nodes", sol.nodes)
```

Both models have the same optimum (the brute force over all 120 orders agrees) and, on five jobs,
the same search (14 nodes each). The weak `M` costs numerical safety, and bound strength as the
instance grows; `lint` flags it before you pay for it (below). Always compute
`M_ij` from the bounds of the variables it multiplies (the `max` of the left side minus the `min`
of the right side).
Bounds on the start times are therefore part of the formulation: without a finite horizon there is
no finite `M`.

## Indicator constraints

When a good `M` is hard to derive, state the logic and let the solver derive `M` from the
presolved bounds. The Python `Model` does not take
indicators; write the LP file, which does (`name: y = 1 -> row`):

```python
import itertools
import json
import os
import subprocess
import tempfile

CLI = os.environ.get("RSV_CLI", "rusolve-cli")
p, r = [4, 3, 5, 2], [0, 2, 1, 8]
n, H = 4, 8 + 14
pairs = list(itertools.combinations(range(n), 2))

def lp_text(kind, big_m=None):
    """kind: "bigm" (rows with M) or "indicator" (logic rows, M derived by the solver)."""
    lines = ["Minimize", " obj: cmax", "Subject To"]
    lines += [f" mk_{i}: cmax - s{i} >= {p[i]}" for i in range(n)]
    for i, j in pairs:
        if kind == "indicator":
            lines.append(f" ij_{i}_{j}: y_{i}_{j} = 1 -> s{i} - s{j} <= {-p[i]}")
            lines.append(f" ji_{i}_{j}: y_{i}_{j} = 0 -> s{j} - s{i} <= {-p[j]}")
        else:
            mij = big_m or (H - p[i]) + p[i] - r[j]
            mji = big_m or (H - p[j]) + p[j] - r[i]
            lines.append(f" ij_{i}_{j}: s{i} - s{j} + {mij} y_{i}_{j} <= {mij - p[i]}")
            lines.append(f" ji_{i}_{j}: s{j} - s{i} - {mji} y_{i}_{j} <= {-p[j]}")
    lines += ["Bounds"] + [f" {r[i]} <= s{i} <= {H - p[i]}" for i in range(n)] + [f" 0 <= cmax <= {H}"]
    lines += ["Binaries", " " + " ".join(f"y_{i}_{j}" for i, j in pairs), "End"]
    return "\n".join(lines) + "\n"

def cli(*args):
    return json.loads(subprocess.run([CLI, *args], capture_output=True, text=True, check=True).stdout)

def brute():
    best = None
    for order in itertools.permutations(range(n)):
        t = 0
        for i in order:
            t = max(t, r[i]) + p[i]
        best = t if best is None else min(best, t)
    return best

with tempfile.TemporaryDirectory() as d:
    files = {}
    for kind, bm in (("tight", None), ("weak", 1000), ("indicator", None)):
        files[kind] = os.path.join(d, kind + ".lp")
        open(files[kind], "w").write(lp_text("indicator" if kind == "indicator" else "bigm", bm))
    # lint: the weak model is flagged row by row, the computed-M model is clean
    flagged = lambda k: sum(g["count"] for g in cli("lint", files[k], "--format", "json")["findings"] if g["rule"] == "big-m")
    assert flagged("tight") == 0 and flagged("weak") == 2 * len(pairs), (flagged("tight"), flagged("weak"))
    print("big-m findings: computed M", flagged("tight"), "M = 1000:", flagged("weak"), "of", 2 * len(pairs), "rows")
    # all three formulations have the optimum of the brute force
    for kind in ("tight", "weak", "indicator"):
        res = cli("mip", files[kind], "--format", "json")
        assert res["status"] == "Optimal" and round(res["objective"]) == brute() == 14, (kind, res["objective"])
        print(kind, res["class"], res["engine"], res["objective"])
    print(res["general_constraints"]["routes"][0])
```

The envelope's `general_constraints.routes` says how each indicator was handled (a big-M row from
the presolved bounds, or a native route when no finite `M` exists), and `verified` says the logic
was re-checked on the rounded binaries. An indicator whose variable has no finite bounds is a typed
error naming the constraint; add the bound. Do not assume indicators are faster than a tight
hand-written `M`; measure both with `rusolve-cli tune` on your own instances.

## Native intervals instead of hand-written big-M

When the model is jobs on machines or tasks on a shared resource, state it with intervals and
let rusolve build the disjunctions (`docs/design/cp-model.md`). An interval is `start + size =
end`; a no-overlap replaces every pair of big-M rows of one machine; a cumulative states a
resource with a capacity. The solve linearizes a copy of the model (big-M per pair from the
bounds, time-indexed for cumulative) and checks the answer on your intervals, so the report's
`verified` is a check of the schedule, not of generated rows.

Two surfaces state the same objects. Python (`pip install rusolve`) has `add_interval`,
`add_no_overlap` and `add_cumulative`. An agent without Python, for example with only
`install.sh` (the binary and the MCP server), uses the `cp` section of the MCP `build_model`
tool. Rules of the `cp` section, all checked when the model is built (a violation is an error
naming its JSON path, for example `cp.intervals[0].end`):

- `start` and `end` are the names of variables you declared in `variables`: type `integer`, and
  both `lower` and `upper` given (an omitted `upper` is infinite, which is refused). The interval
  does not create them. `start`, `size`, `end` and `name` are required; only `presence` is
  optional (a binary variable).
- `size` is an integer written without a fraction (`3`, not `3.0`) or a variable name; it must
  not be negative.
- `no_overlap` and `cumulative` list intervals by their interval `name`, not by variable.
- The reply to `build_model` has the handle `model_id` (`m1`, `m2`, ...) in `structuredContent`
  and in the JSON text of `content[0].text`. Pass it to `solve` and `get_solution`.
- `solve` returns the envelope; `status`, `objective` and `cp.verified` are the three things to
  read. `get_solution` pages the column values by name.

A two-machine flow (every job on machine 1, then machine 2), three jobs:

```python
import itertools
import rusolve as rs

jobs = [(2, 3), (3, 1), (2, 2)]                       # (time on machine 1, time on machine 2)
H = sum(a + b for a, b in jobs)
m = rs.Model()
cmax = m.add_var(ub=H, integer=True, cost=1)
first, second = [], []
for a, b in jobs:
    i1 = m.add_interval(m.add_var(ub=H, integer=True), a)
    i2 = m.add_interval(m.add_var(ub=H, integer=True), b)
    m.add_constr(i2.start >= i1.end)                  # machine 1 before machine 2
    m.add_constr(cmax >= i2.end)
    first.append(i1)
    second.append(i2)
m.add_no_overlap(first)
m.add_no_overlap(second)
sol = m.solve()

def brute():
    best = None
    for p1, p2 in itertools.product(itertools.permutations(range(3)), repeat=2):
        t, end1 = 0, {}
        for j in p1:
            t += jobs[j][0]
            end1[j] = t
        t = 0
        for j in p2:
            t = max(t, end1[j]) + jobs[j][1]
        best = t if best is None else min(best, t)
    return best

assert sol.cp["verified"] and round(sol.objective) == brute()
print(sol.objective, sol.cp["engine"])
```

The same model as the `arguments` of an MCP `build_model` call (optimum 8; `solve` reports
`objective: 8.0` and `cp.verified: true`):

```json
{
  "sense": "minimize",
  "variables": [
    {"name": "cmax", "type": "integer", "lower": 0, "upper": 13, "objective": 1},
    {"name": "s1_0", "type": "integer", "lower": 0, "upper": 13},
    {"name": "e1_0", "type": "integer", "lower": 0, "upper": 13},
    {"name": "s2_0", "type": "integer", "lower": 0, "upper": 13},
    {"name": "e2_0", "type": "integer", "lower": 0, "upper": 13},
    {"name": "s1_1", "type": "integer", "lower": 0, "upper": 13},
    {"name": "e1_1", "type": "integer", "lower": 0, "upper": 13},
    {"name": "s2_1", "type": "integer", "lower": 0, "upper": 13},
    {"name": "e2_1", "type": "integer", "lower": 0, "upper": 13},
    {"name": "s1_2", "type": "integer", "lower": 0, "upper": 13},
    {"name": "e1_2", "type": "integer", "lower": 0, "upper": 13},
    {"name": "s2_2", "type": "integer", "lower": 0, "upper": 13},
    {"name": "e2_2", "type": "integer", "lower": 0, "upper": 13}
  ],
  "rows": [
    {"name": "flow_0", "terms": {"s2_0": 1, "e1_0": -1}, "sense": ">=", "rhs": 0},
    {"name": "mk_0", "terms": {"cmax": 1, "e2_0": -1}, "sense": ">=", "rhs": 0},
    {"name": "flow_1", "terms": {"s2_1": 1, "e1_1": -1}, "sense": ">=", "rhs": 0},
    {"name": "mk_1", "terms": {"cmax": 1, "e2_1": -1}, "sense": ">=", "rhs": 0},
    {"name": "flow_2", "terms": {"s2_2": 1, "e1_2": -1}, "sense": ">=", "rhs": 0},
    {"name": "mk_2", "terms": {"cmax": 1, "e2_2": -1}, "sense": ">=", "rhs": 0}
  ],
  "cp": {
    "intervals": [
      {"name": "m1_0", "start": "s1_0", "size": 2, "end": "e1_0"},
      {"name": "m2_0", "start": "s2_0", "size": 3, "end": "e2_0"},
      {"name": "m1_1", "start": "s1_1", "size": 3, "end": "e1_1"},
      {"name": "m2_1", "start": "s2_1", "size": 1, "end": "e2_1"},
      {"name": "m1_2", "start": "s1_2", "size": 2, "end": "e1_2"},
      {"name": "m2_2", "start": "s2_2", "size": 2, "end": "e2_2"}
    ],
    "no_overlap": [
      {"name": "machine1", "intervals": ["m1_0", "m1_1", "m1_2"]},
      {"name": "machine2", "intervals": ["m2_0", "m2_1", "m2_2"]}
    ]
  }
}
```

## Job shop (precedence within a job, one machine at a time, makespan)

Each job is a fixed sequence of operations; each operation needs one machine for its duration;
a machine does one operation at a time; minimise the time the last operation ends. Recipe:

1. One interval per operation, with integer start and end variables in `[0, H]`, where `H` is
   a valid horizon (the sum of all durations always is).
2. Precedence inside a job: `start(next) >= end(previous)`, one row per consecutive pair.
3. One `no_overlap` per machine, listing every operation that runs on it.
4. Makespan: an integer variable `makespan` with objective 1 and `makespan >= end(last
   operation of the job)` for every job. Forgetting a job gives a too-short schedule with
   status `Optimal`.

Three jobs, two machines (`A`: machine 0 for 3 then machine 1 for 2; `B`: machine 1 for 2 then
machine 0 for 4; `C`: machine 0 for 2 then machine 1 for 3). In Python:

```python
import rusolve as rs

jobs = {"A": [(0, 3), (1, 2)], "B": [(1, 2), (0, 4)], "C": [(0, 2), (1, 3)]}   # [(machine, duration)]
H = sum(d for ops in jobs.values() for _, d in ops)
m = rs.Model()
makespan = m.add_var(ub=H, integer=True, cost=1)
on_machine = {}
for ops in jobs.values():
    prev = None
    for machine, d in ops:
        iv = m.add_interval(m.add_var(ub=H, integer=True), d)
        if prev is not None:
            m.add_constr(iv.start >= prev.end)                       # precedence inside the job
        on_machine.setdefault(machine, []).append(iv)
        prev = iv
    m.add_constr(makespan >= prev.end)                               # the job's last operation
for ivs in on_machine.values():
    m.add_no_overlap(ivs)                                            # one machine, one at a time
sol = m.solve()                       # CP-SAT: m.solve(options={"mip_engine": "cpsat"})
assert sol.cp["verified"] and round(sol.objective) == 9
```

The same model as the `arguments` of an MCP `build_model` call, then `solve {"model_id": "m1"}`
(with the id from the reply). Optimum 9, `cp.verified: true`:

```json
{
  "sense": "minimize",
  "variables": [
    {"name": "makespan", "type": "integer", "lower": 0, "upper": 16, "objective": 1},
    {"name": "s_A0", "type": "integer", "lower": 0, "upper": 16},
    {"name": "e_A0", "type": "integer", "lower": 0, "upper": 16},
    {"name": "s_A1", "type": "integer", "lower": 0, "upper": 16},
    {"name": "e_A1", "type": "integer", "lower": 0, "upper": 16},
    {"name": "s_B0", "type": "integer", "lower": 0, "upper": 16},
    {"name": "e_B0", "type": "integer", "lower": 0, "upper": 16},
    {"name": "s_B1", "type": "integer", "lower": 0, "upper": 16},
    {"name": "e_B1", "type": "integer", "lower": 0, "upper": 16},
    {"name": "s_C0", "type": "integer", "lower": 0, "upper": 16},
    {"name": "e_C0", "type": "integer", "lower": 0, "upper": 16},
    {"name": "s_C1", "type": "integer", "lower": 0, "upper": 16},
    {"name": "e_C1", "type": "integer", "lower": 0, "upper": 16}
  ],
  "rows": [
    {"name": "prec_A1", "terms": {"s_A1": 1, "e_A0": -1}, "sense": ">=", "rhs": 0},
    {"name": "mk_A", "terms": {"makespan": 1, "e_A1": -1}, "sense": ">=", "rhs": 0},
    {"name": "prec_B1", "terms": {"s_B1": 1, "e_B0": -1}, "sense": ">=", "rhs": 0},
    {"name": "mk_B", "terms": {"makespan": 1, "e_B1": -1}, "sense": ">=", "rhs": 0},
    {"name": "prec_C1", "terms": {"s_C1": 1, "e_C0": -1}, "sense": ">=", "rhs": 0},
    {"name": "mk_C", "terms": {"makespan": 1, "e_C1": -1}, "sense": ">=", "rhs": 0}
  ],
  "cp": {
    "intervals": [
      {"name": "A0", "start": "s_A0", "size": 3, "end": "e_A0"},
      {"name": "A1", "start": "s_A1", "size": 2, "end": "e_A1"},
      {"name": "B0", "start": "s_B0", "size": 2, "end": "e_B0"},
      {"name": "B1", "start": "s_B1", "size": 4, "end": "e_B1"},
      {"name": "C0", "start": "s_C0", "size": 2, "end": "e_C0"},
      {"name": "C1", "start": "s_C1", "size": 3, "end": "e_C1"}
    ],
    "no_overlap": [
      {"name": "machine0", "intervals": ["A0", "B1", "C0"]},
      {"name": "machine1", "intervals": ["A1", "B0", "C1"]}
    ]
  }
}
```

For the CP-SAT engine pass `"options": {"mip_engine": "cpsat"}` in the `solve` arguments; the
envelope's top-level `engine` is then `cp-sat` (instead of `cp-linearization`) and it has a
`cpsat` object. Both engines gave 9 on the instance above.

**When to choose `mip_engine=cpsat`.** The default turns every pair of operations on a machine
into a binary and a pair of big-M rows, so it grows with the square of the operations per
machine and its bound weakens with the horizon. Choose CP-SAT when the default stops with a gap at
its time limit, when it fails with `RS4013` (a variable size or demand in a cumulative, a
horizon above 1000 time points, a size that can be 0 in a no-overlap), or when the problem has
many operations per machine. It needs the OR-Tools library (`rusolve-cli version --format json`,
`cpsat.available`; `install.sh --with-cpsat`, or `pip install 'rusolve[cpsat]'`) and every column
of the model integer with finite bounds; otherwise it fails with `RS0012`. Measure both engines
on your own instances before settling on one.

Limits of the default build: a cumulative needs constant sizes and demands and a horizon of at
most 1000 time points, otherwise `solve` raises `RS4013` naming the object. Pass
`options={"mip_engine": "cpsat"}` (MCP `solve`: `options: {"mip_engine": "cpsat"}`) to solve the
model with OR-Tools CP-SAT instead, when it is installed; the answer is checked on your objects
the same way. Bound every start and size: the big-M comes from the bounds. `explain`, `iis`,
`ranging` and `write` refuse a model with these objects (`RS4012`). Keep the hand-written big-M
or indicator form above when the same model must be written to an LP or MPS file.

## Classic mistakes

- **One global `M = 1e6` ("big enough").** The LP relaxation is then useless, and a value near
  `1e-6` in `y` buys a real violation (the trickle-flow failure that the indicator checker
  guards against). Compute `M_ij` per row.
- **No horizon.** Unbounded start times force a huge `M`. Put a valid upper bound `H` on every
  start and the makespan.
- **Symmetry from identical jobs.** Interchangeable jobs produce equivalent `y` assignments;
  fix an order between identical jobs (`y_ij = 1` for `i < j`).
- **Minimising the wrong thing.** Makespan needs `cmax >= s_i + p_i` for every job (the `mk_`
  rows); forgetting one job gives a too-short schedule with status `Optimal`.
- **Time-indexed models for long horizons.** One binary per job and time step grows with
  `sum(p)`. Use it only for small integral horizons.

## What lint says

`rusolve-cli lint model.lp` reports `big-m` per row: the coefficient, the tightest valid value
and the exact fix text ("set the coefficient of y_1_3 in ij_1_3 to 36 and the upper bound to
32"). On the weak model (`M = 1000`) it flags every disjunction row (12 of 12 in the example); on the computed-`M` model it
is silent. In the JSON, `count` is the number of flagged rows and `items` lists only the first
`--max-items` (default 10), so sum `count`, not `len(items)`. Apply the suggested value only after checking it against the start-time bounds you
intend to keep. `coef-range` and `condition-estimate` also grow with a huge `M`.

## Options

- `--opt mip_rel_gap=...`, `--time-limit`: set them; report `gap` from the result envelope.
- `mip_ti_precedence_sepa` and `mip_ti_list_sched` (both default true) act only on a
  time-indexed model (exactly-one start-time sets; time-indexed formulations after Dyer and
  Wolsey 1990, list scheduling after Kolisch 1996). They do nothing for the
  big-M model above.
- Everything else: `rusolve-cli tune` on a family of instances before changing a default.

## Verify

Brute force on a small instance, as in `brute()`, or recompute the schedule from the solution: sort
the jobs by `s`, check that consecutive jobs do not overlap, that every `s_i >= r_i`, and that
`cmax` is the last finish. `rusolve-cli check model.lp solution.json` re-evaluates every row and,
for LP files with indicators, the logic on rounded binaries.
