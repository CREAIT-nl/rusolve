---
name: rusolve-set-cover-packing
description: Model and solve set covering, set packing and set partitioning MIPs with rusolve.Model (crew and shift cover, facility coverage, column selection, conflict-free choice). Use when binary picks must cover every element at least once, at most once or exactly once.
---

# Set cover, set packing, set partitioning

Binary `x_s` says whether set `s` is picked. For each element `e` the row sums the picks of the
sets that contain `e`:

| class | row | objective | when |
|---|---|---|---|
| cover | `sum >= 1` | min cost | every element must be served, extra service is allowed |
| packing | `sum <= 1` | max value | picks must not overlap, elements may stay unused |
| partition | `sum == 1` | min cost or max value | served exactly once (crew pairings, exact tilings) |

`rusolve classify` names these rows: `detail.row_types.set_covering`, `set_packing`,
`set_partitioning`. When one of them counts most rows, this playbook applies. The LP relaxation of
a cover is usually strong, so the work is in the formulation and the gap tolerance, not in
exotic options.

## Formulation

One function builds all three; `kind` only changes the row sense and the objective sense.

```python
import itertools
import rusolve as rs

universe = range(8)
sets = {"A": ({0, 1, 2}, 4), "B": ({2, 3, 4}, 3), "C": ({4, 5, 6}, 3),
        "D": ({6, 7, 0}, 3), "E": ({1, 3, 5, 7}, 5), "F": ({0, 4}, 2),
        "G": ({1, 2, 3}, 4), "H": ({5, 6, 7}, 4)}      # name: (elements, cost)

def build(kind):
    m = rs.Model(sense="max" if kind == "packing" else "min")
    x = {s: m.add_var(binary=True, cost=c, name=f"pick_{s}") for s, (_, c) in sets.items()}
    sense = {"cover": ">=", "packing": "<=", "partition": "=="}[kind]
    for e in universe:
        row = [x[s] for s, (els, _) in sets.items() if e in els]
        m.add_constr(rs.quicksum(row), sense, 1, name=f"elem_{e}")
    return m, x

def brute(kind):                                   # all 2^8 subsets, the independent check
    best = None
    for pick in itertools.product([0, 1], repeat=len(sets)):
        chosen = [s for s, p in zip(sets, pick) if p]
        hits = [sum(e in sets[s][0] for s in chosen) for e in universe]
        ok = {"cover": all(h >= 1 for h in hits), "packing": all(h <= 1 for h in hits),
              "partition": all(h == 1 for h in hits)}[kind]
        if ok:
            z = sum(sets[s][1] for s in chosen)
            if best is None or (z > best if kind == "packing" else z < best):
                best = z
    return best

for kind in ("cover", "packing", "partition"):
    m, x = build(kind)
    sol = m.solve()
    assert sol.status == "Optimal", (kind, sol.status)
    assert round(sol.objective) == brute(kind), (kind, sol.objective, brute(kind))
    print(kind, sol.objective, sol.selected(x))
```

For a set-partitioning model that has no solution the status is `Infeasible`; relax one family of
rows to `>=` to see whether it is a cover (the test in `brute` shows the three cases differ).

## Larger instances

Never build rows with a Python double loop over all elements and sets. Build the incidence once
(`element -> [set indices]`) and add the rows with `m.add_constrs(rows, ">=", 1)` or one sparse
matrix with `m.add_matrix(A, ">=", 1, vars=x)` (`A` a scipy.sparse 0/1 matrix, elements by sets).

## Classic mistakes

- **`==` where the business wants a cover.** Partition rows can make a feasible cover infeasible
  (no exact tiling exists). Ask whether over-coverage is really forbidden.
- **`<=` rows with a minimising objective.** Packing needs a maximising one; minimising a packing
  returns the empty pick with objective 0 and the status `Optimal`. Check that the objective
  sense matches the class.
- **Duplicate sets or dominated sets** (same elements, higher cost). They add symmetry and no
  value; drop the dominated ones before the solve, and the near-parallel rows that `lint` reports
  as `near-parallel-rows` usually point at duplicated elements.
- **Integer instead of binary** columns with no upper bound: a set picked twice is never better
  for a cover, but the missing bound weakens presolve. Use `binary=True`.
- **An element in no set.** The model is infeasible for a cover and partition; `explain` names
  the row (`elem_k`).

## What lint says

`rusolve-cli lint model.lp`: `near-parallel-rows` (two elements covered by the same sets, merge
them), `empty-row` (an element no set contains), `fixed-col`. A `big-m` finding means another
class was modelled in; see `scheduling` or `facility-location`.

## Options

- `--opt mip_rel_gap=0.01` when a 1 % gap is enough (default 1e-4). Report the gap with the status.
- `--time-limit S` on big instances; the envelope's `dual_bound` and `gap` say how far the best
  cover is from proven optimal.
- Leave symmetry detection (`mip_detect_symmetry`, default true) on for covers with many
  interchangeable sets (identical rows and columns). Measure with `rusolve-cli tune` before
  turning an option off.

## Verify

The brute-force function in the example is the pattern for small data: enumerate or use a second
formulation. For any size, `rusolve-cli check model.lp solution.json` re-evaluates every cover
row and integrality of the picks; recompute the cost yourself from the picked names.
