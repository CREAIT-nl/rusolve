---
name: rusolve-scheduling
description: Model machine scheduling and sequencing with disjunctive constraints in rusolve.Model, choose between big-M, tightened big-M and indicator constraints, and find a weak big-M with lint. Use when jobs share a machine or resource, order matters, and a makespan, tardiness or completion-time objective is minimised.
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
