---
name: rusolve-lot-sizing
description: Model single-item and multi-item lot sizing (production planning with setup costs and inventory) in rusolve.Model, with a tight setup-linking bound instead of a weak big-M. Use when production happens in batches, a setup is paid when anything is produced, and demand is met from stock or production over several periods.
---

# Lot sizing

Plan production `x_t`, end-of-period stock `I_t` and a setup binary `y_t` over `T` periods with
demand `d_t`. Stock balance is a flow conservation; the setup is the linking row
`x_t <= M_t * y_t`. Everything about solve speed is in `M_t`.

## 1. Uncapacitated: checked against Wagner-Whitin

```python
import rusolve as rs

d = [20, 40, 0, 35, 50, 15]
T, K, h, c = len(d), 90, 1.0, 2.0          # setup cost, holding cost per unit-period, unit production cost

def build(cap=None, big_m=None, relax=False):
    m = rs.Model()
    x = m.add_vars(T, lb=0, ub=cap if cap else rs.inf, cost=c, names="x")
    I = m.add_vars(T, lb=0, cost=h, names="I")
    y = m.add_vars(T, lb=0, ub=1, cost=K, binary=not relax, names="y")
    for t in range(T):
        prev = I[t - 1] if t else 0
        m.add_constr(prev + x[t] - I[t] == d[t], name=f"balance_{t}")
        # tightest valid M: you never produce more than the remaining demand, nor more than capacity
        remaining = sum(d[t:])
        M = big_m if big_m else min(remaining, cap or remaining)
        m.add_constr(x[t] - M * y[t] <= 0, name=f"setup_{t}")
    return m, x, I, y

def wagner_whitin():
    f = [0.0] + [None] * T
    for t in range(1, T + 1):
        f[t] = min(f[j - 1] + K + sum(h * d[k] * (k - (j - 1)) for k in range(j - 1, t)) for j in range(1, t + 1))
    return f[T] + c * sum(d)

m, x, I, y = build()
sol = m.solve()
assert sol.status == "Optimal"
assert abs(sol.objective - wagner_whitin()) < 1e-6, (sol.objective, wagner_whitin())
assert round(sol.objective) == 620
print(sol.objective, "setups in", sol.selected({t: y[t] for t in range(T)}))
```

The independent check is the Wagner-Whitin dynamic program (Wagner and Whitin, *Management
Science* 5, 1958): the optimum of the uncapacitated problem is a shortest path over "produce in
period `j` for demands `j..t`". When a formulation disagrees with that dynamic program, the
formulation is wrong, not the solver.

## 2. Capacitated: tight versus weak M

```python
import rusolve as rs

d = [20, 40, 0, 35, 50, 15]
T, K, h, c, cap = len(d), 90, 1.0, 2.0, 60

def build(big_m=None, relax=False):
    m = rs.Model()
    x = m.add_vars(T, lb=0, ub=cap, cost=c, names="x")
    I = m.add_vars(T, lb=0, cost=h, names="I")
    y = m.add_vars(T, lb=0, ub=1, cost=K, binary=not relax, names="y")
    for t in range(T):
        prev = I[t - 1] if t else 0
        m.add_constr(prev + x[t] - I[t] == d[t], name=f"balance_{t}")
        M = big_m if big_m else min(cap, sum(d[t:]))
        m.add_constr(x[t] - M * y[t] <= 0, name=f"setup_{t}")
    return m, x, I, y

tight = build()[0].solve()
weak = build(big_m=5000)[0].solve()
assert tight.status == weak.status == "Optimal"
assert abs(tight.objective - weak.objective) < 1e-4                    # same optimum
lp_tight = build(relax=True)[0].solve().objective
lp_weak = build(big_m=5000, relax=True)[0].solve().objective
assert lp_weak < lp_tight <= tight.objective + 1e-9, (lp_weak, lp_tight)
print("optimum", round(tight.objective, 3), "root LP bound: tight M", round(lp_tight, 2), "M=5000", round(lp_weak, 2))

# feasibility recomputed from the solution, not from the solver's say-so
xs, Is = tight.x[:T], tight.x[T:2 * T]
stock = 0
for t in range(T):
    stock = stock + xs[t] - d[t]
    assert abs(stock - Is[t]) < 1e-6 and xs[t] <= cap + 1e-9 and stock >= -1e-9
```

The root bound with the tight `M` is the number that matters: it decides how many nodes the search
needs. With `M = 5000` the relaxation sets `y_t = x_t / 5000` and the setup cost nearly vanishes
(the assertion `lp_weak < lp_tight` checks that on this instance).

## Classic mistakes

- **`M` = a large constant for every period.** Use `min(capacity, remaining demand)`; the
  remaining demand shrinks over time, so later periods get a smaller `M`.
- **Capacity on `x` but not in `M`.** If `cap < remaining demand`, `M = remaining demand` is
  valid but weaker than `cap`.
- **Backlog by accident.** A free (negative) stock variable lets demand be met "later" without
  cost. Keep `I_t >= 0` (`lb=0`) unless backlogging is part of the problem and priced.
- **No initial or final stock.** The balance row for `t = 0` has no `I_{-1}`; set an initial stock
  as a constant on the right side.
- **Equality instead of flow conservation after a shift in indices.** A stock row that mixes
  `I_{t}` and `I_{t-1}` wrongly is feasible but gives a plan of the wrong shape; the dynamic-program
  check above catches it on the uncapacitated case.

## Stronger formulations (when the root gap is still large)

The big-M linking row is the weakest valid description. The standard strengthenings are the
`(l, S)` inequalities (Barany, Van Roy and Wolsey, *Mathematical Programming Study* 22, 1984) and
the facility-location (flow) reformulation `w_{ts}` = amount produced in `t` for demand in `s`,
with `w_{ts} <= d_s * y_t` (Krarup and Bilde 1977; Pochet and Wolsey, *Production Planning by
Mixed Integer Programming*, 2006, ch. 7). The reformulation has `T^2/2` columns but needs no
cut separation; try it when `lint` is clean and the gap after 10 000 nodes is still large.

## What lint says

`rusolve-cli lint model.lp` flags the weak linking rows by name with the exact replacement
("set the coefficient of y3 in link_3 to -80 ...") when `M` is larger than the activity bound
implied by the bounds of `x`. That is why `x` should carry its capacity as an upper bound. Also
`rhs-range` when demand is in different units from costs.

## Options

- `--opt mip_rel_gap=...` and `--time-limit`.
- Multi-item problems with a shared capacity row: leave the defaults, use `rusolve-cli tune` on
  a family of instances (`docs/skills/tune/SKILL.md`).
- Symmetry (`mip_detect_symmetry`, true by default) matters when items or machines are identical.

## Verify

Uncapacitated single item: Wagner-Whitin (above). Otherwise recompute the stock balance and
capacity from the solution (the last lines of block 2), and check the solution file with
`rusolve-cli check model.lp solution.json`.
