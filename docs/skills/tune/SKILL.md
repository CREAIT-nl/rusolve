---
name: rusolve-tune
description: Tune rusolve's options for a family of similar MIPs (fewer nodes or less time, same answer). Use when asked to speed up, tune or choose settings for a set of related models.
---

# Tune the solver for a problem family

The option registry says which options are worth touching and which can change the answer.

## Steps

1. Baseline. Solve every instance of the family with defaults and record, per instance, status,
   objective, dual bound, nodes, LP iterations and time (`rusolve.result/1` fields `status`,
   `objective`, `stats.nodes`, `stats.lp_iterations`, `stats.time_s`). The default search is
   deterministic for a fixed seed, so node counts are a stable metric; use time only as a
   secondary one.
   - MCP: `mcp__rusolve__load_model` and `mcp__rusolve__solve` per instance.
   - CLI: `rusolve-cli mip inst.mps --format json --seed 0`.
2. Look at the family. `rusolve-cli probe inst.mps --format json` runs a short bounded solve and
   reports root LP, presolve and tree features; `rusolve-cli classify inst.mps --json` gives the
   structure. Use them to pick which options to try instead of sweeping all of them.
3. Pick candidates from `rusolve-cli options --json` (or `docs/agent/reference.md#options`):
   - `affects: "search"` and `tunable: true` options change the tree: heuristic effort, cut and
     presolve switches, branching and age limits, symmetry detection.
   - `affects: "tolerance"` (for example `mip_rel_gap`, `mip_abs_gap`) trade proof strength for
     speed and can change the reported optimum: allowed only if the task accepts a gap.
   - `affects: "runtime"` / `"output"` never help the search.
   - `rusolve-cli options --json --profile recommended` shows the project's recommended profile.
4. Try one change at a time, on every instance, same seed:
   `mcp__rusolve__solve` with `options: {"name": value}` or `--opt name=value` (repeatable;
   `--options-file FILE` for a set). A misspelled option exits with code 2 and names the nearest
   valid options. Keep a change only if it helps on the whole family and every instance keeps its
   status and objective (within the agreed tolerance).
5. Guard against overfitting: tune on part of the family and check the rest, and re-run with a
   second `--seed`.
6. Report the final options as a map, the per-instance and total node/time numbers against the
   baseline, and which instances got worse.

## One command for the whole loop

`rusolve-cli tune DIR --seeds 3 --format json -o profile.json` does steps 1 to 5 in one run: it
splits the family into training and held-out instances, races one-exchange candidates (and your
own, `--candidates FILE|-`) over several seeds with a statistical test that drops losers early,
charges every run in deterministic work units (`--budget-work`), checks every answer against
`--known` optima (a wrong answer disqualifies the candidate), and writes a profile only when the
held-out ratio beats the baseline with a confidence interval below 1.0. Among candidates within the
noise it keeps the one that changes the fewest options.

- `rusolve-cli tune --space` lists the options it mutates.
- To propose candidates yourself: run `rusolve-cli tune DIR --format json`, read
  `profile.evidence.train_instances` (static features and the cost of every candidate per instance)
  and `evidence.candidates`, write `{"candidates": [{"name": "...", "options": {...}}]}` and run
  again with `--candidates cand.json --candidates-only`. A bad proposal costs budget, never
  correctness.
- `--iterations 3` also races pairs and triples of options. Read `evidence.verdict` first:
  `baseline_kept` is a valid answer.

## Rules

- Never set `time_limit` or `node_limit` as a "speed-up": that is a different, weaker answer.
- Report the speed-up as a ratio of totals, not of per-instance best cases.
- Options that change the answer (gap, objective cutoff) must be called out explicitly.
