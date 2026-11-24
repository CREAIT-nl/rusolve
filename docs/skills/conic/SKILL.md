---
name: rusolve-conic
description: Model and solve second-order cone programs (SOCP), rotated cones and convex quadratically constrained programs (QCQP) with rusolve, continuous or with integer variables (the experimental misocp_oa lever), and check the answer and the certificates. Use when a constraint bounds a norm or a quadratic form (risk budget, ellipsoid, distance, chance constraint), not only the objective.
---

# Conic programs: SOCP, rotated cones, convex QCQP, and the mixed-integer lever

Three router classes belong here. Read `class` (and `engine`) in the envelope of
`rusolve-cli mip model --format json` (see "What classify says" below for why not `classify`).

| `class` | what the model has | who solves it |
|---|---|---|
| `conic` | quadratic rows and/or cones, no integer column (a Q objective is allowed) | Clarabel interior point (`engine: clarabel`), always on |
| `mi_conic` | the same objects plus integer columns | refused by default; with `--opt misocp_oa=true`, outer approximation over the MIP search (`engine: oa`) |
| `nonconvex_qcqp` | a quadratic row that is not convex on its finite side, a quadratic equality, an indefinite Q | nothing: exit 3 with a typed error naming the row, the variable and the side |

A convex quadratic *objective* with linear rows is not this skill: see `rusolve-convex-qp`.
A quadratic objective with integers (MIQP) is there too, behind its own lever (`miqp_oa`).

## 1. Write the model

**Quadratic row, LP file.** Square brackets, monomial coefficients, **no factor one half**:
`q1: [ 0.04 a ^ 2 + 0.02 a * b + 0.09 b ^ 2 ] <= 1` means `0.04 a^2 + 0.02 ab + 0.09 b^2 <= 1`.
The objective bracket is the other way round: `[ ... ] / 2` is the `1/2 x'Qx` convention. A row
with a quadratic term has one side only (`<=` or `>=`, not a range), and linear terms may sit
beside the bracket. A row is convex on its `<=` side when its Q is positive semidefinite, and on
its `>=` side when `-Q` is (so a `>=` row with a PSD Q is a reverse-convex constraint: refused).

**Cone, MPS file only.** The LP format has no cone syntax. After `BOUNDS`, one section per cone
lists columns in order:

    CSECTION k1 0.0 QUAD     x0 >= ||(x1, x2, ...)||
    CSECTION k2 0.0 RQUAD    2 x0 x1 >= ||(x2, ...)||^2, with x0, x1 >= 0

A `QCMATRIX <row>` section holds the quadratic part of a row (`<col> <col> <coef>` lines, monomial
coefficients, duplicates summed). Exponential and power cones (`PEXP`, `PPOW`) are
reserved and rejected. The columns of a cone keep their own bounds: the MPS default lower bound
0 makes a cone column nonnegative, so declare free columns `FR`. A norm bounded by a *variable*
(`||(x, y)|| <= t`) is a cone; as a quadratic row it would be `x^2 + y^2 - t^2 <= 0`, which is
nonconvex and refused. Squared norms bounded by a constant, or by a linear function through a
rotated cone, are rows.

## 2. A quadratic constraint, checked on a grid

Maximise `3x + 4y` with `x^2 + y^2 <= 25` and `x + y <= 6`. The disc alone is optimal at
`(3, 4)` (value 25), which breaks the linear row, so the answer is the corner where the line
meets the circle: `x = 3 - sqrt(14)/2`, value `21 + sqrt(14)/2`.

```python
import json, math, os, subprocess, tempfile

CLI = os.environ.get("RSV_CLI", "rusolve-cli")
work = tempfile.mkdtemp()

def cli(*args):
    return subprocess.run([CLI, *args], capture_output=True, text=True)

model = os.path.join(work, "disc.lp")
open(model, "w").write(
    "Maximize\n obj: 3 x + 4 y\nSubject To\n disc: [ x ^ 2 + y ^ 2 ] <= 25\n lin: x + y <= 6\n"
    "Bounds\n x free\n y free\nEnd\n")
r = cli("mip", model, "--format", "json", "--include-solution")
assert r.returncode == 0, r.stderr
env = json.loads(r.stdout)
assert (env["class"], env["engine"], env["status"]) == ("conic", "clarabel", "Optimal")
assert env["conic"]["checker"]["verdict"] == "feasible"       # the user's rows, not the conic form
want = 21 + math.sqrt(14) / 2
assert abs(env["objective"] - want) < 1e-6, (env["objective"], want)

# independent check: scan the circle and the line, both active at the optimum
best = max(3 * x + 4 * y
           for k in range(200001)
           for x, y in [(5 * math.cos(2 * math.pi * k / 200000), 5 * math.sin(2 * math.pi * k / 200000))]
           if x + y <= 6 + 1e-12)
assert best <= env["objective"] + 1e-6 and env["objective"] - best < 1e-3, (best, env["objective"])

# the checker the solver ships with, on the saved envelope
sol = os.path.join(work, "sol.json")
open(sol, "w").write(r.stdout)
chk = json.loads(cli("check", model, sol, "--format", "json").stdout)
assert chk["schema_version"] == "rusolve.check-conic/1" and chk["verdict"] == "feasible"
duals = {d["name"]: d["value"] for d in env["conic"]["row_duals"]}
assert duals["lin"] > 0 and duals["disc"] > 0                  # both rows are binding: positive duals
print("objective", env["objective"], "row duals", duals)
```

## 3. A second-order cone and a rotated cone in MPS

`min t` with `x + y = 2` and `t >= ||(x, y)||` has `t = sqrt(2)` at `x = y = 1`. The rotated cone
`2 u v >= w^2` with `w` fixed at 2 gives `u v >= 2`, so `min u + v` is `2 sqrt(2)`.

```python
import json, math, os, subprocess, tempfile

CLI = os.environ.get("RSV_CLI", "rusolve-cli")
work = tempfile.mkdtemp()

def solve(name, text):
    path = os.path.join(work, name)
    open(path, "w").write(text)
    r = subprocess.run([CLI, "mip", path, "--format", "json"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    env = json.loads(r.stdout)
    assert env["class"] == "conic" and env["status"] == "Optimal", env
    return env

soc = solve("soc.mps", """NAME SOC1
ROWS
 N  obj
 E  e1
COLUMNS
    t  obj  1
    x  e1  1
    y  e1  1
RHS
    rhs  e1  2
BOUNDS
 FR bnd t
 FR bnd x
 FR bnd y
CSECTION k1 0.0 QUAD
    t
    x
    y
ENDATA
""")
assert abs(soc["objective"] - math.sqrt(2)) < 1e-6
cone_duals = soc["conic"]["cone_duals"]                        # one dual vector per user cone
assert len(cone_duals) == 1

rot = solve("rot.mps", """NAME RSOC1
ROWS
 N  obj
COLUMNS
    u  obj  1
    v  obj  1
    w  obj  0
BOUNDS
 PL bnd u
 PL bnd v
 FX bnd w 2
CSECTION k1 0.0 RQUAD
    u
    v
    w
ENDATA
""")
assert abs(rot["objective"] - 2 * math.sqrt(2)) < 1e-6
print("soc", soc["objective"], "rotated", rot["objective"])
```

## 4. Reading the result

The usual `rusolve.result/1` envelope, plus `convexity: "convex"` and a `conic` object:

- `clarabel_status` and `reduced_accuracy`: `true` when Clarabel stopped on an `Almost...` status
  (a lower-accuracy answer). Say so when you report it; do not tighten tolerances hoping it goes
  away.
- `checker`: the same check as `rusolve-cli check`, on *your* rows and cones (not the reformulation),
  with a violation scaled by the size of each row. Read `verdict`.
- `row_duals` (by row name; derivative of the objective with respect to the right-hand side, in the
  model's sense) and `cone_duals` (one vector per cone).
- Infeasible and unbounded answers carry a `certificate` that rusolve verified (`verified: true`):
  a Farkas vector for infeasibility, a ray for unboundedness. A certificate that does not verify
  turns the status into `Unknown`, never into a claim.

The conic route takes `--format`, `--time-limit` and `--include-solution`; `--seed` and
`--threads` are accepted and ignored; there is no `--solution-file` (save the envelope instead).
`check model solution.json` accepts that envelope as the solution.

## 5. Integer columns: the `misocp_oa` lever (experimental)

A model with integer columns and cones or convex quadratic rows is class `mi_conic`. By default
`rusolve-cli mip` and `solve` refuse it with exit 3 and a message that names the lever. The
experimental engine is **off by default because it is slower** and its dual bound can stall:
turn it on only when you need a checked answer for such a model, and always give `--time-limit`.

    rusolve-cli mip model.mps --opt misocp_oa=true --time-limit 60 --format json

(`solve` takes the same `--opt`, or set `RUSOLVE_EKK_OPTS=misocp_oa=true`. Through MCP, `load_model` accepts
the model and `solve` takes `"options": {"misocp_oa": true}`.)

The envelope then has `class: "mi_conic"`, `engine: "oa"` and a `misocp_oa` object (rounds, cuts,
polish counters, `accepted`, `bound_inconsistent`). `Optimal` means the point passed the same
checker as `rusolve-cli check` on the original cones and rows, and the bound closed. A point the
checker rejects is never reported; at the time limit the status is `TimeLimit` with the best checked
point and the bound; a bound that contradicts a checked point gives `Unknown`
(`bound_inconsistent: true`). Related options, all off by default: `misocp_oa_disaggregate=true`
(extended formulation for cones of three or more components), `misocp_oa_lazy_only=true` (keeps
presolve), `misocp_oa_no_polish=true`. `--start`, `--basis` and `--add-rows` are refused with the
lever. It stays off by default because on the pool it was measured on it proved 10 of 18 instances
where a general MISOCP solver proves all 18, with no wrong answers.

Whole batches `a`, `b` under an energy row, checked by enumeration:

```python
import json, math, os, subprocess, tempfile

CLI = os.environ.get("RSV_CLI", "rusolve-cli")
work = tempfile.mkdtemp()
model = os.path.join(work, "batches.lp")
open(model, "w").write(
    "Maximize\n obj: 4 a + 3 b + 2 c\nSubject To\n energy: [ 2 a ^ 2 + b ^ 2 + c ^ 2 ] <= 30\n"
    " weight: a + b + c <= 7\nBounds\n 0 <= a <= 5\n 0 <= b <= 5\n 0 <= c <= 8\nGeneral\n a b\nEnd\n")

off = subprocess.run([CLI, "mip", model, "--format", "json"], capture_output=True, text=True)
assert off.returncode == 3                                        # refused by default ...
assert "misocp_oa=true" in off.stderr                             # ... and the message names the lever
err = json.loads(off.stdout)["error"]
assert err["class"] == "mi_conic" and "misocp_oa=true" in err["hint"]

on = subprocess.run([CLI, "mip", model, "--opt", "misocp_oa=true", "--time-limit", "60", "--format", "json"],
                    capture_output=True, text=True)
assert on.returncode == 0, on.stderr
env = json.loads(on.stdout)
assert (env["class"], env["engine"], env["status"]) == ("mi_conic", "oa", "Optimal"), env
assert env["misocp_oa"]["accepted"] and not env["misocp_oa"]["bound_inconsistent"]

best = 0.0                                                         # enumerate a, b; c is then closed-form
for a in range(6):
    for b in range(6):
        room = min(30 - 2 * a * a - b * b, 64)
        c = min(8, 7 - a - b, math.sqrt(room) if room >= 0 else -1)
        if c >= 0:
            best = max(best, 4 * a + 3 * b + 2 * c)
assert abs(env["objective"] - best) < 1e-4 * (1 + best), (env["objective"], best)
print("objective", env["objective"], "enumeration", best)
```

## Reformulations when the model is refused

- **`nonconvex_qcqp`, a `>=` row with PSD Q** (stay outside a disc): not convex. Linearise
  around a fixed point, or enumerate the sign of the choice with a binary and a big-M (then it is
  a MIP with convex rows on each branch), and say it is a restriction.
- **A quadratic equality.** If the objective pushes the row to one side, replace `=` by that
  inequality and verify the equality holds in the answer; otherwise it is nonconvex.
- **A norm on a variable's right side.** Write a cone (`CSECTION`), not a row (section 1).
- **An indefinite Q** from estimates: symmetrise and add a small ridge (`Q + 1e-8 I`) as in
  `rusolve-convex-qp`, or factor `Q = L L'` and use a cone on `L'x`.
- **Nonconvex piecewise or general constraints beside quadratic terms** (`gc_quadratic`): no engine
  takes both. Separate them or drop one.
- **Exponential and power cones**: approximate by a second-order cone model or a piecewise-linear
  one (`rusolve-piecewise-semi`); say which.

## Classic mistakes

- **The factor of one half** differs between a row (none) and the objective (`/ 2`). A disc of
  radius 1 is `[ x ^ 2 + y ^ 2 ] <= 1`.
- **Cone column bounds.** An MPS cone column without a bound line is `>= 0`; a free variable in a
  cone needs `FR`. Cones cannot be written to an LP file (`convert` errors).
- **`rusolve.Model` has no quadratic constraints or cones**: the Python API covers the quadratic
  *objective* only. Write the LP or MPS text and call the CLI or the MCP `load_model`.
- **Calling the answer exact.** An interior-point answer is accurate to a relative 1e-6 or so;
  compare with a tolerance, never `==`.
- **Equality on integers in the oa engine.** The cuts leave cones feasible to about
  `2 * mip_feasibility_tolerance`; the engine polishes the final integer assignment, and the
  reported point is the polished one.

## What classify says

`rusolve-cli classify` reads LP and MIP models and a quadratic objective; it exits 3 on a model with
quadratic rows or cones (and on general constraints), because its reader does not hold them. Get
the class from the router instead: `rusolve-cli mip model --format json` has `class` and `engine`,
and a refusal prints `error.class` and `error.hint`.

## Verify

`rusolve-cli check model.mps solution.json` checks bounds, linear rows, quadratic rows, cone membership
and integrality on the model as written, and the objective (schema `rusolve.check-conic/1`). Pair it
with an independent calculation on a small instance: a closed form, a scan of the boundary, or
enumeration of the integers as above. A conic problem is convex, so a feasible point whose dual
(`row_duals`, `cone_duals`) matches the primal objective is globally optimal; the envelope's
`objective_primal` and `objective_dual` give both.
