---
name: rusolve-start-here
description: Entry point for any optimization request with rusolve. Classify the model (or the problem statement) first, then open the matching playbook skill for assignment, transportation, network flow, set cover, scheduling, lot-sizing, facility location, routing, piecewise-linear or semi-continuous costs, or convex QP. Use before formulating or solving when the problem class is not already obvious.
---

# Start here: classify, then pick a playbook

Two inputs, one routing step.

**A model file exists.** Run the classifier once and read three fields:

    rusolve-cli classify model.lp --format json

- `class`: `lp`, `mip`, `convex_qp`, `nonconvex_qp`, `miqp`, `conic`, `nonconvex_qcqp`,
  `mi_conic`. It exits with code 3 and a typed message for general constraints (indicators,
  PWL, SOS, MIN/MAX/ABS); for those, ask the router instead:
  `rusolve-cli mip model.lp --format json` has `class` (`gc`, `gc_mip`, `gc_quadratic`),
  `engine` and `general_constraints.routes`.
- `detail.row_types`: counts of `set_partitioning`, `set_packing`, `set_covering`, `knapsack`,
  `variable_bound`, `mixed_binary`, `precedence`, ... (MIPLIB 2017 row types).
- `detail.structure.network_rows` against `detail.rows`: all rows are network rows means a
  network flow (the LP is integral at a vertex).

**Only a description exists.** Match the wording in the table, write the model with
`rusolve.Model` (the playbook has runnable code), then classify the file you wrote.

## Routing table

| signal | playbook |
|---|---|
| `class` `lp`, `network_rows == rows`; or "ship", "supply and demand", "assign workers", "min-cost flow" | `rusolve-assignment-transport-flow` |
| `row_types.set_covering`, `set_packing` or `set_partitioning` dominate; "cover every item", "pick at most one per" | `rusolve-set-cover-packing` |
| binaries ordering jobs, rows with a big coefficient on a binary (`variable_bound`, `precedence`, lint `big-m`); "machine", "before or after", "makespan" | `rusolve-scheduling` |
| inventory balance rows plus setup binaries linked to production (`variable_bound`); "order quantity", "setup cost", "demand per period" | `rusolve-lot-sizing` |
| "open a site", "which depots", fixed cost per opened facility, assignment shares | `rusolve-facility-location` |
| arc binaries with degree rows (`set_partitioning`), "tour", "route", "visit every" | `rusolve-routing` |
| `class` `gc` / `gc_mip`, a `PWL`, indicator or `Semi-continuous` section; "minimum lot", "volume discount", "tiered price" | `rusolve-piecewise-semi` |
| `class` `convex_qp`; "variance", "least squares", "portfolio", a Hessian | `rusolve-convex-qp` |
| none of the above | `rusolve-solve`, then `rusolve-tune` if it is slow |

Refusals are information, not failures. `nonconvex_qp`, `miqp`, `mi_conic`, `nonconvex_qcqp`,
`gc_quadratic` come back as one typed error naming the class:
tell the user the class is not solvable yet and offer a reformulation (linearise, drop the
integer, or convexify) instead of searching for an option that bypasses it.

## Always, whatever the playbook

1. Write the model with `rusolve.Model` and keep an independent check (brute force, a DP, a
   closed form) for the small instance; see each playbook for the one that fits.
2. Run `rusolve-cli lint model.lp --format json`. `big-m` findings carry the tightest valid
   coefficient: apply it. `coef-range` and `obj-range` mean rescale.
3. Solve with the defaults first (`rusolve-solve`). Report status, objective, `dual_bound` and
   `gap`; a non-`Optimal` status is not a result.
4. Only then tune (`rusolve-tune`): `rusolve-cli probe` and `rusolve-cli tune` measure options on
   your instance family. No option changes the answer, only the time and the tree.
5. When the status is `Infeasible`, switch to `rusolve-debug-infeasible`.

## Verify the router

`rusolve-cli classify` on the model you wrote must report the class the playbook expects (for
example `convex_qp` for the portfolio, `mip` with `set_partitioning` rows for the TSP degree
constraints). If it does not, the model is not the one you meant to write.
