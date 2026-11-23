---
name: rusolve-routing
description: Model the traveling salesman problem and capacitated vehicle routing in rusolve.Model with subtour elimination (Miller-Tucker-Zemlin rows or lazily added cut rows). Use for tours, delivery routes, drilling or picking sequences, and any model whose binaries choose arcs of a graph that must form cycles.
---

# Routing: TSP and capacitated VRP

Binary `x_ij` = 1 when arc `i -> j` is used. Degree rows (each node entered once and left once)
describe an assignment problem; its solutions are a set of cycles. The model is only a tour when
**subtours are eliminated**, and that is the one decision in this class.

| formulation | rows | use when |
|---|---|---|
| MTZ (Miller, Tucker and Zemlin, *J. ACM* 7, 1960) | `n^2` big-M rows with an order variable `u_i` | up to a few dozen nodes, simplest to write |
| DFJ cuts (Dantzig, Fulkerson and Johnson, *Operations Research* 2, 1954), added lazily | `sum_{i,j in S} x_ij <= |S| - 1`, few rows | larger instances; solve, find subtours, add them, re-solve |
| load propagation (CVRP) | MTZ-style `u_j >= u_i + q_j - Q (1 - x_ij)` with `q_i <= u_i <= Q` | vehicle capacity |

## 1. TSP with MTZ, checked by brute force

```python
import itertools
import math
import rusolve as rs

pts = [(0, 0), (4, 1), (7, 5), (3, 8), (-2, 5), (-3, 1)]            # node 0 is the depot
n = len(pts)
dist = [[round(math.dist(pts[i], pts[j]) * 10) for j in range(n)] for i in range(n)]

m = rs.Model()
x = {(i, j): m.add_var(binary=True, cost=dist[i][j], name=f"x_{i}_{j}") for i in range(n) for j in range(n) if i != j}
u = {i: m.add_var(lb=1, ub=n - 1, name=f"u_{i}") for i in range(1, n)}      # visit order of node i
for i in range(n):
    m.add_constr(rs.quicksum(x[i, j] for j in range(n) if j != i), "==", 1, name=f"out_{i}")
    m.add_constr(rs.quicksum(x[j, i] for j in range(n) if j != i), "==", 1, name=f"in_{i}")
for i in range(1, n):
    for j in range(1, n):
        if i != j:
            m.add_constr(u[i] - u[j] + (n - 1) * x[i, j] <= n - 2, name=f"mtz_{i}_{j}")
sol = m.solve()

def tour(arcs):
    nxt = dict(arcs)
    t, k = [0], nxt[0]
    while k != 0:
        t.append(k)
        k = nxt[k]
    return t

best = min(sum(dist[a][b] for a, b in zip((0,) + p, p + (0,))) for p in itertools.permutations(range(1, n)))
assert sol.status == "Optimal" and round(sol.objective) == best
t = tour(sol.selected(x))
assert sorted(t) == list(range(n))                                # one tour through every node
print(sol.objective, t)
```

## 2. The classic failure and the lazy-cut fix

Without subtour rows the assignment relaxation reports `Optimal` with several small cycles. The loop
below keeps the model alive, adds the violated subtour rows, and re-solves (the previous basis is
reused); on a symmetric TSP it needs a handful of rounds.

```python
import itertools
import math
import rusolve as rs

pts = [(0, 0), (4, 1), (7, 5), (3, 8), (-2, 5), (-3, 1), (1, 3), (6, 8)]
n = len(pts)
dist = [[round(math.dist(pts[i], pts[j]) * 10) for j in range(n)] for i in range(n)]

m = rs.Model()
x = {(i, j): m.add_var(binary=True, cost=dist[i][j], name=f"x_{i}_{j}") for i in range(n) for j in range(n) if i != j}
for i in range(n):
    m.add_constr(rs.quicksum(x[i, j] for j in range(n) if j != i), "==", 1, name=f"out_{i}")
    m.add_constr(rs.quicksum(x[j, i] for j in range(n) if j != i), "==", 1, name=f"in_{i}")

def cycles(arcs):
    nxt, seen, out = dict(arcs), set(), []
    for s in range(n):
        if s in seen:
            continue
        c, k = [], s
        while k not in seen:
            seen.add(k)
            c.append(k)
            k = nxt[k]
        out.append(c)
    return out

rounds, first = 0, None
while True:
    sol = m.solve()
    assert sol.status == "Optimal"
    cs = cycles(sol.selected(x))
    first = first if first is not None else len(cs)
    if len(cs) == 1:
        break
    for c in cs:                                    # one subtour row per cycle that is not the whole tour
        m.add_constr(rs.quicksum(x[i, j] for i in c for j in c if i != j), "<=", len(c) - 1, name=f"sub_{rounds}_{c[0]}")
    rounds += 1

best = min(sum(dist[a][b] for a, b in zip((0,) + p, p + (0,))) for p in itertools.permutations(range(1, n)))
assert round(sol.objective) == best and rounds >= 1
print("cycles in the first solution:", first, "rounds:", rounds, "tour cost:", sol.objective)
```

The first solution of the degree-only model has more than one cycle (`first` in the output) and
passes every row you wrote: nothing but the subtour rows makes it a tour. Adding rows to a
persistent `Model` and calling `solve()` again re-optimises from the current model; each round is a
full MIP solve, so keep the loop for instances where MTZ is too weak or too large.

## 3. Capacitated VRP with load variables

```python
import itertools
import math
import rusolve as rs

pts = [(0, 0), (4, 1), (7, 5), (3, 8), (-2, 5), (-3, 1)]
q = [0, 4, 6, 5, 3, 4]                                  # demands, node 0 is the depot
Q, K = 12, 2                                            # vehicle capacity, at most K vehicles
n = len(pts)
dist = [[round(math.dist(pts[i], pts[j]) * 10) for j in range(n)] for i in range(n)]

m = rs.Model()
x = {(i, j): m.add_var(binary=True, cost=dist[i][j], name=f"x_{i}_{j}") for i in range(n) for j in range(n) if i != j}
load = {i: m.add_var(lb=q[i], ub=Q, name=f"load_{i}") for i in range(1, n)}   # load on board after serving i
for i in range(1, n):
    m.add_constr(rs.quicksum(x[i, j] for j in range(n) if j != i), "==", 1, name=f"out_{i}")
    m.add_constr(rs.quicksum(x[j, i] for j in range(n) if j != i), "==", 1, name=f"in_{i}")
m.add_constr(rs.quicksum(x[0, j] for j in range(1, n)), "<=", K, name="fleet")
m.add_constr(rs.quicksum(x[0, j] for j in range(1, n)) - rs.quicksum(x[j, 0] for j in range(1, n)), "==", 0, name="depot")
for i in range(1, n):
    for j in range(1, n):
        if i != j:       # load_j >= load_i + q_j when i -> j; M = Q is the tightest valid value
            m.add_constr(load[i] - load[j] + Q * x[i, j] <= Q - q[j], name=f"load_{i}_{j}")
sol = m.solve()

def brute():                                          # best split of every customer order into <= K routes
    best = None
    for p in itertools.permutations(range(1, n)):
        for cut in range(0, n):                       # route 1 = p[:cut], route 2 = p[cut:]; K = 2
            routes = [r for r in (p[:cut], p[cut:]) if r]
            if any(sum(q[i] for i in r) > Q for r in routes):
                continue
            z = sum(dist[a][b] for r in routes for a, b in zip((0,) + r, r + (0,)))
            best = z if best is None else min(best, z)
    return best

assert sol.status == "Optimal" and round(sol.objective) == brute(), (sol.objective, brute())
print(sol.objective, "routes start with", sorted(j for (i, j) in sol.selected(x) if i == 0))
```

The brute force above is written for `K = 2`; it splits each customer order once.

## Classic mistakes

- **Degree rows only.** The model is an assignment problem and returns cycles with status
  `Optimal` (block 2 shows it). Always add subtour rows and check the solution is one cycle.
- **MTZ with `u` unbounded or `n` instead of `n - 1`.** A weak `M` weakens the LP; keep
  `1 <= u_i <= n - 1` and the coefficient `n - 1` (`lint` flags a larger one as `big-m`).
- **Symmetric distances modelled as two directed arcs with no cut between them.** Fine for
  correctness; the directed model has a weaker bound than an undirected one with edge variables.
- **Float distances.** Round to integers (as above) or accept that the objective has solver
  tolerance noise; compare with a tolerance, never `==`.
- **VRP without a fleet bound.** `sum_j x_0j <= K` is needed, or the solver adds vehicles freely.
- **Time windows with a big-M of the whole horizon.** Same weak-`M` mistake as in `scheduling`.

## What lint and classify say

`classify` shows `set_partitioning` rows (the degree rows, an assignment structure) plus
`mixed_binary` rows (MTZ or load rows). `lint` reports `big-m` on MTZ rows whose coefficient is
larger than the activity bound implied by `u`.

## Options

- `--opt mip_rel_gap=...`, `--time-limit`; routing proofs can be slow, so report `gap`.
- Heuristics find good tours early; `mip_heuristic_effort` (default 0.05) can be raised. Compare
  with `rusolve-cli tune` on your own instances.

## Verify

Check the output structure, not only the objective: follow the arcs from the depot (function
`tour`), confirm every node appears once, per-route load is at most `Q`, and recompute the length
from the distance matrix. Brute force (all permutations) is the check up to about 9 nodes.
