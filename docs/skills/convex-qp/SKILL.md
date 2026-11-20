---
name: rusolve-convex-qp
description: Model and solve convex quadratic programs (mean-variance portfolios, least squares with constraints, quadratic costs) with rusolve.Model.set_quadratic_objective, and avoid the classic factor-of-two and non-convexity mistakes. Use when the objective has squared terms or a covariance matrix and the constraints are linear and continuous.
---

# Convex QP and portfolios

`rusolve.Model.set_quadratic_objective(Q)` sets the objective to **`c'x + (1/2) x'Qx`**: the linear
costs stay on the variables (`add_vars(cost=...)`), `Q` is the Hessian. `Q` must be positive
semidefinite for a minimisation (negative semidefinite for a maximisation); `solve()` checks it and
raises a `ValueError` that names the smallest eigenvalue and the variables involved. A QP with
integer variables (a cardinality limit, minimum buy-in) is a MIQP, which the solver refuses today
with a typed error naming the class; for those, use the
epigraph or a heuristic outside the solver and say so.

`rusolve mip` and `rusolve solve` route a convex QP to the QP engine (`class: "convex_qp"`).

## 1. Minimum-variance portfolio, checked on a grid

Four assets, full investment, expected return at least 8 %, no short sales.

```python
import itertools
import rusolve as rs

mu = [0.05, 0.08, 0.11, 0.14]
cov = [[0.04, 0.006, 0.002, 0.0],
       [0.006, 0.09, 0.012, 0.01],
       [0.002, 0.012, 0.16, 0.03],
       [0.0, 0.01, 0.03, 0.25]]
n, target = 4, 0.08

m = rs.Model()
w = m.add_vars(n, lb=0, ub=1, names="w")                    # weights
m.add_constr(rs.quicksum(w[i] for i in range(n)), "==", 1, name="budget")
ret = m.add_constr(rs.quicksum(mu[i] * w[i] for i in range(n)), ">=", target, name="return")
m.set_quadratic_objective(cov)                              # objective = 1/2 w' cov w  (half the variance)
sol = m.solve()
assert sol.status == "Optimal"

def half_var(x):
    return 0.5 * sum(x[i] * cov[i][j] * x[j] for i in range(n) for j in range(n))

x = list(sol.x)
assert abs(sum(x) - 1) < 1e-6 and sum(mu[i] * x[i] for i in range(n)) >= target - 1e-7 and min(x) >= -1e-9
assert abs(sol.objective - half_var(x)) < 1e-9              # the reported objective is 1/2 x'Qx, recomputed

# independent check: no feasible grid portfolio (step 0.025) is better, and the best one is close
step, best = 40, None
for a, b, c in itertools.product(range(step + 1), repeat=3):
    if a + b + c > step:
        continue
    g = [a / step, b / step, c / step, 1 - (a + b + c) / step]
    if sum(mu[i] * g[i] for i in range(n)) < target - 1e-12:
        continue
    z = half_var(g)
    best = z if best is None else min(best, z)
assert sol.objective <= best + 1e-9 and best - sol.objective < 2e-4, (sol.objective, best)
print("weights", [round(v, 4) for v in x], "volatility", round((2 * sol.objective) ** 0.5, 4))
```

## 2. The efficient frontier by re-solving

The model stays alive, so changing the return target is a row-bound edit and the next `solve()`
starts from the previous solution.

```python
import rusolve as rs

mu = [0.05, 0.08, 0.11, 0.14]
cov = [[0.04, 0.006, 0.002, 0.0], [0.006, 0.09, 0.012, 0.01], [0.002, 0.012, 0.16, 0.03], [0.0, 0.01, 0.03, 0.25]]
n = 4
m = rs.Model()
w = m.add_vars(n, lb=0, ub=1, names="w")
m.add_constr(rs.quicksum(w[i] for i in range(n)), "==", 1, name="budget")
r = m.add_constr(rs.quicksum(mu[i] * w[i] for i in range(n)), ">=", 0.06, name="return")
m.set_quadratic_objective(cov)

risk = []
for target in (0.06, 0.08, 0.10, 0.12, 0.14):
    m.set_row_bounds(r, target, rs.inf)
    sol = m.solve()
    assert sol.status == "Optimal"
    risk.append(2 * sol.objective)
assert all(b >= a - 1e-9 for a, b in zip(risk, risk[1:]))   # more return never costs less variance
assert abs(2 * sol.objective - cov[3][3]) < 1e-6            # 14 % is only reachable with all weight on asset 3
print([round(v, 5) for v in risk])
```

## 3. Risk-aversion form and least squares

Same model, different objective: `min -mu'w + (lambda/2) w'Sigma w` puts the return in the linear
costs (`add_vars(cost=-mu)`) and `Q = lambda * Sigma`. Constrained least squares
`min ||Ax - b||^2` is `Q = 2 A'A`, `c = -2 A'b` (drop the constant `b'b`). In both cases the
Hessian is a PSD matrix by construction; if you assemble it from estimates, symmetrise it and
add a small ridge (`Q + 1e-8 I`) when `solve()` reports a tiny negative eigenvalue from rounding.

## Classic mistakes

- **The factor of one half.** The objective is `1/2 x'Qx`. Passing `Q = Sigma` gives half the
  variance (as above); to minimise the variance itself pass `Q = 2 Sigma`. The minimiser is the
  same, the reported objective is not.
- **A Hessian that is not PSD** (a sample covariance from fewer observations than assets, or a
  hand-typed matrix): `ValueError: nonconvex quadratic objective ...`. Repair the matrix; do not
  look for an option that disables the check.
- **Maximising a convex quadratic.** `Model(sense="max")` needs a negative semidefinite `Q`.
- **Only the upper triangle given by mistake.** A `{(i, j): v}` dict with `i <= j` is symmetrised
  (see `crates/rusolve-py/tests/test_qp.py`); a dense matrix is taken as written. Give a full
  symmetric matrix when in doubt.
- **Mixed scales.** Returns of order 0.1 and variances of order 1e-4 are fine; coefficients
  spanning more than 1e6 are not; `rusolve-cli lint` reports `coef-range` and `obj-range`.
- **Integer weights ("lots") with a quadratic objective.** MIQP is refused with a typed error; see
  the top of this page.

## What lint and classify say

`rusolve-cli classify model.lp --format json` gives `class: "convex_qp"`. Lint is about the linear
part (`coef-range`, `obj-range`, `fixed-col`). The Hessian's PSD check happens at `solve()`.

## Options

Defaults first. The QP engine's option names are in `rusolve-cli options --json` (`qp_*`:
`qp_iteration_limit`, `qp_regularization_value`, `qp_allow_hot_start`); the router reports the
engine used in the envelope (`engine`). Set `--time-limit`. If a solve hits an iteration limit, the
status says so; raise the limit before suspecting the model.

## Verify

Recompute `1/2 x'Qx + c'x` from the returned weights (as in `half_var`), check every linear row
and bound, and compare with a coarse brute-force or a second formulation (grid above). Because
the problem is convex, a point that satisfies the first-order (KKT) conditions is globally
optimal; `sol.duals` and `sol.reduced_costs` are not available for QPs, so use the grid or
perturbation check rather than a dual certificate.
