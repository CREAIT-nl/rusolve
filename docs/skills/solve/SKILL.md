---
name: rusolve-solve
description: Solve an LP or MIP with rusolve and report status, objective, bound, gap and the solution by name. Use when asked to solve, optimise or run a model file (MPS, LP, .gz).
---

# Solve a model with rusolve and report

Prefer the MCP tools; the CLI gives the same JSON. Full parameter lists:
`docs/agent/reference.md`.

## With MCP

1. `mcp__rusolve__load_model` with `path` (or `mps_text`). Keep the returned `model_id`; the
   result lists rows, columns, integer columns and names.
2. `mcp__rusolve__solve` with that `model_id`. Optional: `time_limit` (seconds), `options` (a flat
   map such as `{"mip_rel_gap": 0.01}`; names come from `rusolve-cli options --json`),
   `include_solution`. The result is a `rusolve.result/1` envelope.
3. Read `status` first. It is a closed enum (Optimal, Infeasible, Unbounded,
   InfeasibleOrUnbounded, TimeLimit, NodeLimit, SolutionLimit, Interrupted, Error, ...; the full
   list is in `docs/schema/rusolve-result-1.schema.json`). Do not report an objective as optimal unless the
   status is `Optimal`; on a limit status report the incumbent, the dual bound and the gap.
4. `mcp__rusolve__get_solution` for values by column name (`names`, or `offset`/`limit` to page).
5. Optional safety net: `mcp__rusolve__check_solution` re-checks a solution against bounds, rows
   and integrality and lists any violation by name.
6. `mcp__rusolve__free_model` when done.

## With the CLI

    rusolve-cli mip model.mps --format json --solution

`rusolve mip` runs the HiGHS-port path that produces the envelope by default (`RUSOLVE_MIP_HIGHS=0` selects the native search). `--time-limit S`,
`--node-limit N`, `--seed N`, `--opt key=value` (repeatable) and `--options-file FILE` set limits
and options; `--solution-file out.sol` writes the HiGHS solution format (`.json` writes JSON).
stdout is exactly one JSON object; the solver log is on stderr.

## Reporting checklist

- Status string exactly as the envelope gives it; objective with its sense (`min`/`max`).
- For a MIP: `dual_bound` and `gap`, `stats.nodes`. For an LP: iteration count.
- Solution by column name, nonzeros only unless asked otherwise.
- If status is Infeasible or InfeasibleOrUnbounded, do not guess: follow
  `docs/skills/debug-infeasible/SKILL.md`.
- Exit code 0 means the solver ran, not that the model is optimal. 2 is a usage error (read the
  message: it names the nearest valid option), 3 an unreadable model, 4 an internal error.
