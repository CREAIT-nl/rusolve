---
name: rusolve-piecewise-semi
description: Model piecewise-linear costs (volume discounts, tiered prices, convex cost curves), semi-continuous variables (minimum lot size, minimum power when on) and indicator-style logic with rusolve, choosing the formulation that needs the fewest integers. Use when a cost or a bound changes shape with the value, or a quantity is either zero or at least some minimum.
---

# Piecewise-linear costs and semi-continuous variables

Three shapes, three formulations, in order of cost:

| shape | formulation | integers added |
|---|---|---|
| **convex** cost, minimised (tiered cost that rises with volume) | epigraph: `c >= a_k x + b_k` for every piece `k` | none: pure LP |
| **semi-continuous** `x = 0` or `lb <= x <= ub` | `Model.set_semi(x)` (no binary written by hand) | handled by the solver |
| **non-convex** cost (volume discount, fixed charge plus slope) | `PWL` general constraint in an LP file, or a binary per piece | solver adds them (incremental or log formulations; Vielma, Ahmed and Nemhauser, *Operations Research* 58, 2010) |

The router (`rusolve mip`) classifies a model with these objects as `gc` / `gc_mip` and reports
how each was handled in the envelope's `general_constraints.routes` (`epigraph_lp`,
`incremental`, ...).

## 1. Convex piecewise-linear cost: an LP

Buy `D = 100` units. The first 40 cost 2 each, the next 40 cost 3, beyond that 5 (slopes rise:
convex). A second source sells at a flat 4.2 up to 70 units. Minimise the total.

```python
import rusolve as rs

pieces = [(2, 0), (3, -40), (5, -200)]       # (slope, intercept): cost >= 2x, 3x - 40, 5x - 200
D = 100
m = rs.Model()
x = m.add_var(lb=0, ub=D, name="x")          # units from source 1
cx = m.add_var(lb=0, cost=1.0, name="cost1") # epigraph variable: its value is source 1's cost
y = m.add_var(lb=0, ub=70, cost=4.2, name="y")
for k, (a, b) in enumerate(pieces):
    m.add_constr(cx - a * x >= b, name=f"piece_{k}")
m.add_constr(x + y >= D, name="demand")
sol = m.solve()

def tiered(v):
    return max(a * v + b for a, b in pieces)

assert sol.status == "Optimal"
xv, yv = sol[x], sol[y]
assert abs(sol[cx] - tiered(xv)) < 1e-6                              # the epigraph variable sits on the curve
best = min(tiered(v) + 4.2 * max(0, D - v) for v in range(0, D + 1) if D - v <= 70)   # grid check
assert abs(sol.objective - best) < 1e-6, (sol.objective, best)
print(sol.objective, "x =", xv, "y =", yv)
```

This only works because the objective pushes `cx` down onto the curve. If the cost were
*maximised*, or a convex curve appeared on the wrong side of a constraint (`cost <= ...`), the
epigraph is wrong and a binary formulation is required.

## 2. Semi-continuous: minimum lot size

Machine `a` is off or produces between 15 and 40 units at cost 1 per unit; `b` makes up the rest at
`cost_b` per unit (up to 50). Demand is 12. Whether `a` should run depends on `cost_b`: running
it wastes `15 - 12 = 3` units.

```python
import rusolve as rs

def solve(cost_b, semi):
    m = rs.Model()
    a = m.add_var(lb=15, ub=40, cost=1.0, name="a")
    b = m.add_var(lb=0, ub=50, cost=cost_b, name="b")
    m.add_constr(a + b >= 12, name="demand")
    if semi:
        m.set_semi(a)        # a = 0 or 15 <= a <= 40;  m.set_semi(a, integer=True): 0 or an integer in range
    sol = m.solve()
    assert sol.status == "Optimal"
    return sol.objective, sol[a], sol[b]

obj, a, b = solve(cost_b=1.6, semi=True)       # b is dear: switch a on, waste 3 units
assert abs(obj - 15.0) < 1e-9 and abs(a - 15.0) < 1e-9
obj, a, b = solve(cost_b=0.9, semi=True)       # b is cheap: leave a off
assert abs(obj - 10.8) < 1e-9 and abs(a) < 1e-9 and abs(b - 12.0) < 1e-9
obj, a, b = solve(cost_b=0.9, semi=False)      # the same data with a plain lower bound: a >= 15 is forced
assert abs(obj - 15.0) < 1e-9                  # 15.0 instead of 10.8: the plain bound is the wrong model
print("semi-continuous and plain bound differ as expected")
```

## 3. Non-convex piecewise-linear cost and semi-continuous together, from an LP file

A cost curve with a cheap first block, a surcharge block, then a cheap block again is not convex.
The Python `Model` has no general constraints, so write the LP file (`General Constraints`
section, `y = PWL ( x ) : (x0, y0) (x1, y1) ...`, and a `Semi-continuous` section) and solve it with
`rusolve mip`. The model: plant `a` is semi-continuous (0 or 15..40) with the non-convex cost
`PWL`, plant `b` has a flat cost of 3.1 per unit (up to 50), demand is 60.

```python
import json
import os
import subprocess
import tempfile

LP = """Minimize
 obj: ca + 3.1 b
Subject To
 demand: a + b >= 60
Bounds
 15 <= a <= 40
 0 <= b <= 50
 0 <= ca <= 400
Semi-continuous
 a
General Constraints
 pa: ca = PWL ( a ) : (0, 0) (15, 30) (25, 90) (40, 110)
End
"""
cost = lambda a: 30 * a / 15 if a <= 15 else 30 + 6 * (a - 15) if a <= 25 else 90 + (a - 25) * 20 / 15   # the curve

# independent check: the objective is piecewise linear in a, so its minimum is at a breakpoint or a bound
cands = [a for a in (0, 15, 25, 40) if 60 - a <= 50]          # a = 0 would need b = 60 > 50
best = min(cost(a) + 3.1 * (60 - a) for a in cands)

with tempfile.TemporaryDirectory() as d:
    path = os.path.join(d, "plant.lp")
    open(path, "w").write(LP)
    r = subprocess.run([os.environ.get("RSV_CLI", "rusolve-cli"), "mip", path, "--format", "json", "--solution"],
                       capture_output=True, text=True, check=True)
res = json.loads(r.stdout)
assert res["status"] == "Optimal" and res["class"] == "gc", res
assert abs(res["objective"] - best) < 1e-6 and abs(best - 169.5) < 1e-9
print(res["class"], res["engine"], res["objective"], res["general_constraints"]["routes"])
```

The envelope's `general_constraints.verified` is `true` when the solver re-checked the logic
(the PWL value and the semi-continuous condition) on the final solution, and `violations` lists any
object that failed. Treat `verified: false` as not solved.

Other general constraints in the LP format: `y = ABS ( x )`, `y = MAX ( a , b , 1 )`, `MIN`,
`AND`, `OR`, indicator rows `z = 1 -> row`, SOS rows (`name: S1 :: x : 1 y : 2`). Convex directions
(minimising `|x|`, or `max(...)` pushed down) become pure LPs; the other direction adds a binary.

## Classic mistakes

- **A hand-written fixed-charge row with a huge `M`** (`x <= 1e6 z`) for a minimum lot size.
  `set_semi` (or the LP `Semi-continuous` section) states the intent and lets the solver derive the
  bound from `ub`; give the column a finite upper bound.
- **Treating `lb` as always binding.** After `set_semi` the lower bound applies only when the
  variable is on; the contrast in block 2 shows the two models give different answers.
- **The epigraph for a non-convex curve.** The rows `c >= a_k x + b_k` describe the convex function
  `max_k(a_k x + b_k)`, which is not your curve once the slopes stop rising. Check the slopes are
  non-decreasing before using block 1; if not, use the `PWL` general constraint of block 3.
- **Breakpoints out of order.** List them in increasing `x`. A cost that jumps (a fixed charge) is
  not continuous piecewise-linear: use a semi-continuous column or a binary.
- **Reading a gc model with `classify`.** `rusolve classify` refuses general constraints and says
  so; use `rusolve mip --format json` (the `class` field) instead.

## What lint says

The big-M lint rules act on written rows; a model that uses `set_semi` or `PWL` has none to flag,
which is the point. If you do see `big-m`, a fixed-charge or logic row was written by hand and
can usually be replaced by a semi-continuous column or an indicator.

## Options

- Nothing special: `--opt mip_rel_gap=...` and `--time-limit` on the `mip` call.
- The incremental PWL route adds one binary per piece; for curves with many breakpoints the
  logarithmic formulation uses fewer. Look at `general_constraints.added_columns`
  in the envelope to see what the model grew by.

## Verify

Evaluate the true cost curve at the solution (the `cost` function above) and compare with the
reported objective; check semi-continuous columns are 0 or inside their bounds; and read
`general_constraints.verified`. `rusolve-cli check model.lp solution.json` also evaluates the PWL,
indicator and semi-continuous logic on a solution file.
