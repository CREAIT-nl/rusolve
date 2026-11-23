---
name: rusolve-facility-location
description: Model capacitated and uncapacitated facility location (which sites to open, which customers each serves) in rusolve.Model with strong linking constraints. Use for plant, depot, warehouse or hub placement, and for any fixed-charge network design where opening an arc or site costs a fixed amount.
---

# Facility location

Binary `open_j` for each site, continuous `x_ji` in `[0, 1]` for the share of customer `i`'s demand
served from site `j`. Rows: every customer fully served, site capacity only when open, and the
**strong linking** `x_ji <= open_j`. The linking is the whole game: it is redundant for integer
points and decisive for the LP bound.

## Formulation, checked by enumeration

```python
import itertools
import rusolve as rs

fixed = [60, 70, 50, 65]                     # cost of opening a site
cap = [45, 50, 40, 50]                       # capacity of a site
dem = [10, 12, 8, 15, 9]                     # customer demand
c = [[4, 6, 9, 5, 8], [7, 3, 5, 6, 4], [6, 8, 3, 7, 5], [5, 7, 6, 3, 9]]    # unit cost site -> customer
J, I = range(len(fixed)), range(len(dem))

def build(linking="strong", relax=False, open_set=None):
    m = rs.Model()
    op = {j: m.add_var(lb=0, ub=1, binary=not relax and open_set is None, cost=fixed[j], name=f"open_{j}") for j in J}
    x = {(j, i): m.add_var(lb=0, ub=1, cost=dem[i] * c[j][i], name=f"x_{j}_{i}") for j in J for i in I}
    for i in I:
        m.add_constr(rs.quicksum(x[j, i] for j in J), "==", 1, name=f"serve_{i}")
    for j in J:
        m.add_constr(rs.quicksum(dem[i] * x[j, i] for i in I) - cap[j] * op[j] <= 0, name=f"cap_{j}")
    if linking == "strong":
        for j in J:
            for i in I:
                m.add_constr(x[j, i] - op[j] <= 0, name=f"link_{j}_{i}")
    if open_set is not None:                       # the enumeration check: sites fixed open or closed
        for j in J:
            m.set_bounds(op[j], 1 if j in open_set else 0, 1 if j in open_set else 0)
    return m, op, x

def enumerate_best():
    best = None
    for k in range(1, len(fixed) + 1):
        for s in itertools.combinations(J, k):
            if sum(cap[j] for j in s) < sum(dem):
                continue
            sol = build(open_set=set(s))[0].solve()
            if sol.status == "Optimal" and (best is None or sol.objective < best[0]):
                best = (sol.objective, s)
    return best

m, op, x = build()
sol = m.solve()
assert sol.status == "Optimal"
best = enumerate_best()
assert abs(sol.objective - best[0]) < 1e-6, (sol.objective, best)
opened = tuple(sol.selected(op))
print(round(sol.objective, 3), "open", opened, "enumeration:", best)

# strong linking versus the aggregated form: same optimum, weaker root bound
weak = build(linking="aggregated")[0].solve()
assert abs(weak.objective - sol.objective) < 1e-6
lp_strong = build(relax=True)[0].solve().objective
lp_weak = build(linking="aggregated", relax=True)[0].solve().objective
assert lp_weak < lp_strong <= sol.objective + 1e-9, (lp_weak, lp_strong)
print("root LP bound: strong", round(lp_strong, 2), "capacity row only", round(lp_weak, 2))
```

`build(open_set=...)` fixes the open sites and leaves a transportation LP; trying every subset is
a slow but independent check on a handful of sites (the solver's MIP answer must equal the best
subset). Never ship the enumeration; it is the test.

## Uncapacitated and single-source variants

- **Uncapacitated:** drop the capacity rows. The strong linking `x_ji <= open_j` is then the only
  thing connecting the sites, and the LP bound is often near-integral; a MIP solve is quick.
- **Single sourcing** (each customer from exactly one site): declare `x_ji` binary
  (`binary=True`), keep `serve_i` as `== 1`, keep the strong links. The model becomes a harder
  MIP; the capacity rows turn into knapsacks (cover cuts apply). `rusolve-cli probe model.lp`
  reports the root behaviour.
- **Fixed-charge network design** (open an arc `a` with capacity `u_a` to carry commodity flow
  `f_a`): the same shape, `f_a <= min(demand, u_a) * y_a`. `mip_netdesign_sepa` (default true)
  adds strong linking cuts `x <= min(d, u) y` on capacity rows and checks each commodity before
  using it (strong linking inequalities; Gendron, Crainic and Frangioni 1999).

## Classic mistakes

- **Only the aggregated link `sum_i x_ji <= n * open_j`** (or only the capacity row). Same integer
  solutions, a much weaker relaxation (the assertion `lp_weak < lp_strong`). Add the
  disaggregated rows `x_ji <= open_j`; they number `|J||I|` but presolve and the cut pool
  handle them.
- **Capacity row without the `open_j` term** (`sum d_i x_ji <= cap_j`): a closed site still serves
  customers for free. Subtract `cap_j * open_j`.
- **Symmetric identical sites** make many equivalent optima. Order them with
  `open_j >= open_{j+1}` only when the sites are truly identical (same cost, capacity, and costs
  to customers).
- **Demands that exceed total capacity** are infeasible; check `sum(cap) >= sum(dem)` first.
- **Distances in different units** make `lint` report `coef-range`; scale demand-times-cost.

## What lint and classify say

`classify` shows many `variable_bound` rows (`detail.row_types.variable_bound`, the links) and a
few knapsack-like capacity rows. `lint` is silent on the strong model; on a model with a huge
constant in the capacity or link rows it reports `big-m` with the tightest value.

## Options

- `mip_rel_gap` and `--time-limit` as usual; report gap and `dual_bound`.
- Heuristic effort (`mip_heuristic_effort`, default 0.05) helps when the optimality proof is slow
  but good solutions are wanted early: raise it and compare with `rusolve-cli tune`.
- Keep `mip_netdesign_sepa` on (default) for network-design models.

## Verify

Enumerate (above) for small cases. Otherwise recompute the cost from the opened sites and shares:
`sum(fixed[j])` over open sites plus `sum(dem[i] * c[j][i] * x[j][i])`, and check each customer's
shares add to 1 and each site's load fits its capacity. `rusolve-cli check model.lp solution.json`
does the row checks.
