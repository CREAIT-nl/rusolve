---
name: rusolve-debug-infeasible
description: Find why a model is infeasible and repair it with rusolve, by row and column name. Use when a solve returns Infeasible or InfeasibleOrUnbounded, or when asked to debug, explain or fix an infeasible model.
---

# Debug an infeasible model

Work from the solver's own certificates, not from reading the MPS by hand.

## Steps (MCP)

1. `mcp__rusolve__load_model`, then `mcp__rusolve__solve`. Confirm the status really is
   Infeasible (or InfeasibleOrUnbounded: `explain` also settles which).
2. `mcp__rusolve__explain` with the `model_id`. For a MIP it explains the LP relaxation; if the LP
   is feasible but the MIP is not, use the CLI `rusolve-cli explain model.mps --mip-iis` (an
   irreducible infeasible subsystem found with integrality). The result names a small set of rows
   and column bounds that cannot hold together (an IIS) and a verified Farkas certificate
   (`explanation.farkas.proof.statement` says in words why the combination is contradictory).
3. `mcp__rusolve__feasrelax` with `mode: "min-count"` for the fewest changes, or `"weighted"` for
   the cheapest total change (`penalties` per row/column, a negative weight forbids relaxing it).
   It returns named edits ("change lower bound of row dem1 from 30 to 25").
4. Decide which edit is the modelling fix. The solver finds feasible repairs, not the right one:
   a demand row relaxed by 5 and a capacity raised by 5 are both valid, but only the user knows
   which is true. When the task fixes what may change (for example bounds only), choose among
   the edits that respect it.
5. Apply with `mcp__rusolve__modify` (edit kinds: `col_bounds`, `row_bounds`, `cost`, `coef`,
   `add_row`, `add_col`, `del_row`, `del_col`, all by name; a bound left out is unchanged, `null`
   is infinite), then `mcp__rusolve__solve` again. Confirm the status is Optimal; the basis is kept
   so the re-solve is cheap.
6. Optionally `mcp__rusolve__lint` for the cause class (big-M, a bound that contradicts a row);
   its findings carry a fix hint.

## CLI equivalents

    rusolve-cli explain model.mps --format json
    rusolve-cli feasrelax model.mps --mode min-count --format json
    rusolve-cli diagnose model.mps --format json     # status + lint + explanation in one call
    rusolve-cli lint model.mps --format json

## Report

Name the IIS rows and columns, say in one sentence why they conflict (use the certificate, not a
guess), give the edit(s) in `modify` format, and state the status of the repaired model. Do not
claim a fix without re-solving. `rusolve-cli check model.mps solution.json` verifies a solution you
are handed.
