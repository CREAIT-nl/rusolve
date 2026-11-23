---
name: rusolve-assignment-transport-flow
description: Model and solve assignment, transportation and min-cost network flow problems with rusolve.Model as plain LPs, without integer variables. Use when workers go to jobs, goods ship from supplies to demands, or units flow through a graph with arc capacities.
---

# Assignment, transportation and network flow

Read `docs/skills/start-here/SKILL.md` first if you have not classified the model. These problems
are the easy class: the constraint matrix is a network matrix (each column has at most one `+1`
and one `-1` in the node-balance form), so the **LP relaxation already has an integral optimal
vertex** when supplies, demands and capacities are integers (Ahuja, Magnanti and Orlin, *Network
Flows*, 1993, ch. 9, the integrality property; Schrijver, *Theory of Linear and Integer
Programming*, 1986, ch. 19, total unimodularity). Declare the variables continuous and solve an
LP: no branching, and `rusolve classify` reports `class: "lp"`.

## 1. Transportation

Plants with supply, markets with demand, a unit cost per lane. Total supply equals total demand,
so both sides are equalities.

```python
import rusolve as rs

supply, demand = [30, 40], [20, 25, 25]
cost = [[8, 6, 10], [9, 12, 13]]

m = rs.Model()
ship = {(i, j): m.add_var(lb=0, cost=cost[i][j], name=f"ship_{i}_{j}")
        for i in range(2) for j in range(3)}
for i in range(2):
    m.add_constr(rs.quicksum(ship[i, j] for j in range(3)), "==", supply[i], name=f"supply_{i}")
for j in range(3):
    m.add_constr(rs.quicksum(ship[i, j] for i in range(2)), "==", demand[j], name=f"demand_{j}")
sol = m.solve()

assert sol.status == "Optimal"
flows = sol.values(ship, above=1e-9)
assert all(abs(v - round(v)) < 1e-9 for v in flows.values())     # integral without integer columns
assert sum(cost[i][j] * v for (i, j), v in flows.items()) == round(sol.objective)   # recomputed cost
assert round(sol.objective) == 640      # brute force over all 411 integral plans agrees
print(sol.status, sol.objective, flows)
```

Verify the answer without trusting the solver's objective: recompute the cost from the flows
(second assertion), check every supply and demand row, and read `sol.duals`: for a minimisation
LP the dual of `demand_j` is the price of one more unit at market `j`.

## 2. Assignment

`n` workers, `n` jobs, one each. Use continuous variables in `[0, 1]`; the vertex is 0/1.

```python
import itertools
import rusolve as rs

cost = [[9, 2, 7], [6, 4, 3], [5, 8, 1]]
n = 3
m = rs.Model()
x = {(i, j): m.add_var(lb=0, ub=1, cost=cost[i][j], name=f"x_{i}_{j}") for i in range(n) for j in range(n)}
for i in range(n):
    m.add_constr(rs.quicksum(x[i, j] for j in range(n)), "==", 1, name=f"worker_{i}")
for j in range(n):
    m.add_constr(rs.quicksum(x[i, j] for i in range(n)), "==", 1, name=f"job_{j}")
sol = m.solve()

best = min(sum(cost[i][p[i]] for i in range(n)) for p in itertools.permutations(range(n)))   # brute force, n! cases
assert sol.status == "Optimal" and round(sol.objective) == best == 9
assert sorted(sol.selected(x, threshold=0.5)) == [(0, 1), (1, 0), (2, 2)]
print(sol.objective, sol.selected(x))
```

## 3. Min-cost flow with arc capacities

Node balance: outflow minus inflow equals the node's net supply (positive at sources, negative at
sinks, zero in between). Capacities are upper bounds on the arc variables, not extra rows.

```python
import rusolve as rs

net = {"s": 10, "a": 0, "b": 0, "t": -10}
arcs = {("s", "a"): (2, 6), ("s", "b"): (5, 8), ("a", "b"): (1, 4), ("a", "t"): (6, 5), ("b", "t"): (2, 9)}  # (cost, capacity)
m = rs.Model()
f = {arc: m.add_var(lb=0, ub=cap, cost=c, name=f"f_{arc[0]}_{arc[1]}") for arc, (c, cap) in arcs.items()}
for v, b in net.items():
    out_ = [f[a] for a in arcs if a[0] == v]
    in_ = [f[a] for a in arcs if a[1] == v]
    m.add_constr(rs.quicksum(out_) - rs.quicksum(in_), "==", b, name=f"balance_{v}")
sol = m.solve()
assert sol.status == "Optimal"
flow = sol.values(f, above=1e-9)
assert flow == {("s", "a"): 5.0, ("s", "b"): 5.0, ("a", "b"): 4.0, ("a", "t"): 1.0, ("b", "t"): 9.0}, flow
assert round(sol.objective) == 63       # unique optimum, checked by enumerating all integral flows
print(sol.objective, flow)
```

If total supply exceeds total demand, use `<=` on supply rows (or a dummy sink node with cost 0).
If the capacities cannot carry the demand the status is `Infeasible`: run
`rusolve-cli explain model.lp` (see `docs/skills/debug-infeasible/SKILL.md`); the IIS names the
cut of saturated arcs.

## Classic mistakes

- **Integer or binary columns on a pure network.** They add branching for nothing. Leave them
  continuous; the optimum is integral anyway. (Exception: side constraints that break the network
  structure, such as "at most 3 arcs used"; then it is a MIP and `facility-location` or
  `scheduling` applies.)
- **Unbalanced totals with equality rows.** Supply 70 and demand 75 with `==` everywhere is
  infeasible by arithmetic. Check `sum(supply) == sum(demand)` before solving.
- **Missing or loose capacities.** An arc without `ub` is uncapacitated; the model silently
  routes everything through the cheapest arc.
- **Mixed units.** Costs in cents on some arcs and euros on others make `lint` report
  `coef-range` or `obj-range`; fix the units, do not tune the solver.
- **Reading a degenerate dual as a price.** Several optimal duals exist when a supply row is
  tight with a zero-flow arc; report ranging (`rusolve-cli ranging model.lp`) instead.

## What lint and classify say

`rusolve-cli classify model.lp --format json` gives `class: "lp"` and `detail.structure.network_rows`
equal to the number of rows. `rusolve-cli lint model.lp` is normally empty; `coef-range`,
`rhs-range` or `obj-range` mean mixed units.

## Options

Defaults are right. For a very large LP try `--opt solver=ipm` or `--opt solver=simplex`
(`options --json` lists `solver`, `simplex_strategy`, `presolve`) and compare on the same file; a
network LP is also a candidate for the primal simplex (`simplex_strategy=4`). Measure before you
change anything; the work counters in the result envelope (`stats`) make the comparison exact.

## Verify

`rusolve-cli check model.lp solution.json` re-evaluates every row and bound against a solution
file; for integral data, also confirm integrality of the flows (assertion above).
