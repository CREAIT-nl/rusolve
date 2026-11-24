# Agent-task evals

How well can an AI agent drive rusolve? Each task gives a real agent (Claude Code headless) a
prompt and a model file, lets it work through the MCP server or the CLI, and grades its final
answer with a deterministic grader whose ground truth is computed by rusolve itself at grading
time (no expected numbers are stored). The runner records success, tool calls, tool errors,
turns, tokens and cost.

| task | kind | agent must | grader checks (ground truth from rusolve) |
|---|---|---|---|
| `solve_report` | solve | solve a 12-item MIP, report status, objective, selected columns | status and objective equal `rusolve mip`; the selection is feasible with that objective (`check`) |
| `debug_infeasible` | debug | find why an LP is infeasible and repair it with bound edits only | at most `feasrelax --mode min-count` edits, accepted by `modify`, model re-solves Optimal |
| `verify_formulation` | verify | check a claimed solution of a MIP | feasible flag and the exact set of violated rows/columns equal `rusolve check` |
| `tune_family` | tune | find options cutting total nodes to 80% on a 3-instance family | options accepted, all Optimal within 1% of the optimum, nodes at most 0.8 x the solve without it |
| `add_cuts` | extend | derive valid cuts that close at least half the root gap of a knapsack | each cut valid (optimised over the original MIP by rusolve), LP bound after the cuts closes at least 50% of the gap |
| `extend_wasm` | extend | write a Wasm (`.wat`) node-selection rule, save it, load it through the MCP `solve` tool's `wasm_hook`, and beat the tree of the HiGHS child rule (`mip_child_selection=-1` on both solves) | the component loads and the solve is Optimal with the hook-less objective; `stats.wasm` calls > 0 and ok > 0 with traps/out_of_fuel/invalid 0; nodes at most 0.8 x the solve without it |
| `class_router` | class-router | route four models (flow, cover, QP, PWL) to class and playbook | class equals the `rusolve mip` envelope `class`; playbook is the skill for that class |
| `class_flow`, `class_cover`, `class_pwl_semi`, `class_qp` | class-* | solve the model, report class, objective, playbook (and weights / engine) | status, objective and class from `rusolve mip`; QP weights within 2e-3 |
| `class_conic` | class-conic | solve a quadratically constrained model, report class, engine, objective, values, playbook | status, objective, class, engine from `rusolve mip`; values within 2e-3 |
| `class_miqp`, `class_misocp` | class-qp, class-conic | solve a convex MIQP / mixed-integer quadratic-row model that is refused by default; find the lever from the refusal message and the skill | truth from `rusolve mip` with `miqp_oa=true` / `misocp_oa=true` (never stored): class `miqp` / `mi_conic`, engine `oa`, status, objective, values (and the held set for the MIQP) |
| `class_scheduling`, `class_lotsizing` | class-* | solve, report the number of rows the `big-m` lint flags | objective from `rusolve mip`; `big_m_rows` equals the sum of `count` of `lint` big-m findings |
| `class_facility`, `class_routing` | class-* | solve, report the open sites / the tour | open set equals rusolve's; the tour is one cycle with the optimal length |

### Word problems: the agent writes the model

Five tasks give no model file, only a word problem; the agent must formulate it and solve it with
`rusolve.Model` (`--mode python`). Kind `model`, `"modes": ["python", ...]` in `task.json`.

| task | class | what a model must get right |
|---|---|---|
| `diet_plan` | LP | four foods, three nutrient minimums, a per-food cap and a share rule (beans at most 40% of calories, a ratio row) |
| `project_selection` | 0-1 | budget, two prerequisites (E needs A, F needs B and C), an exclusion, a minimum count, maximise |
| `fixed_charge_supply` | MIP | open/ship linking (`ship <= cap * open`), exact demands, 15 columns |
| `shift_cover` | integer | five-day shifts wrapping over the weekend, weekend premium, whole workers |
| `supplier_sourcing` | general constraints | an indicator (minimum lot of 25 t and a fixed fee if supplier B delivers), an SOS1 (A and C not both), a max (freight on the largest delivery); the reference is `selfcheck/reference.lp` |

`class_cp_scheduling` is a word problem for the scheduling objects (`--mode python` only): a
three-job, two-machine shop the agent states with `add_interval` and `add_no_overlap`. Those
objects have no file syntax, so the reference is built in Python by `evals/cp_reference.py` (run it
with `RSV_PY`, an interpreter that has the `rusolve` package), solved at grading time, and the
reported plan is judged with `Model.check` on the intervals and no-overlap objects. The selfcheck
feeds it the reference answer, a wrong objective, no plan, and a plan with the same makespan that
overlaps two operations on M1 (resp. M2), which must be rejected naming `machine1` (`machine2`).

Ground truth comes from rusolve itself: `tasks/<id>/selfcheck/reference.mps` is the same problem
written by hand with the variable names the prompt dictates, solved at grading time. The grader
requires the reported status and objective to match, then judges the agent's `plan` against that
reference with `rusolve check` (so a plan that breaks a stated rule fails whatever objective it
reports) and re-derives the plan's own objective. `graders.py --selfcheck` feeds each grader the
reference answer (must pass) and three wrong ones (objective off, plan breaks a rule, no plan; all
must fail). The prompt never mentions the reference files.

`supplier_sourcing` runs in `--mode python` only: the highs-mcp tool has no general constraints.

`mcp_build_model` is a sixth word problem for the MCP server alone (`"modes": ["mcp"]`): a bakery
integer program (two resource rows, a cap on cakes, `rolls <= 1.5 loaves`) that the agent must write
with `build_model` JSON or `load_model` LP text, with no MPS file given. Same grader and reference
layout (`selfcheck/reference.mps`, integer columns with explicit upper bounds, since an MPS integer
column without one defaults to 1). Its rule is binding: a model without the `rolls_cap` row reaches
profit 630, the right model 501.

The selfcheck also feeds each word-problem grader one **rule-breaking** answer per structural row:
the reported objective and the plan's own objective both equal the optimum, but the plan breaks
`beans_share` (diet), `E_needs_A` and `C_or_D` (projects), `cap_2` (supply), `cover_mon` (shifts) or
`rolls_cap` (bakery), or `b_off`, `a_or_c`, `b_min` (suppliers: the indicator, SOS1 and minimum-lot objects, which `check` verifies on the reference model's own objects). The grader must reject it and name that row; a rejection for another reason
counts as BROKEN. The runner checks every task's mode before the first run, so a task that does not
fit `--mode` fails at once instead of after the paid runs before it.

`add_cuts` is the MCP-level proxy of "add a cut family": the agent derives the cuts, rusolve
proves validity and measures the bound. `extend_wasm` is the task where code the agent wrote
really runs inside the solver, as a fuel-metered WebAssembly component on the Wasm hook surface
(`docs/design/wasm-hooks.md`). Writing a Rust `MipHooks::separate` is still not graded here.

## extend_wasm

`model.mps` is a 30-item, 4-constraint 0-1 knapsack (`gen_models.py`: `multi_knapsack(..., 30, 4,
seed 7)`), whose search needs a real tree. Both solves use the HiGHS child rule
(`mip_child_selection=-1`): the default hybrid rule dives so much that the open-node queue is
consulted only about a dozen times, and no simple node-selection rule then beats it by 20 %.
Measured with `RUSOLVE_MIP_HIGHS=1 rusolve mip model.mps --format json --seed 0 --opt
mip_child_selection=-1` (release build):

| run | nodes | `stats.wasm` |
|---|---|---|
| no component (the best-bound pop) | 912 | absent |
| `selfcheck/reference_depth_first.wat` (dive into the deepest open node) | 534 | calls 13, ok 13, traps 0 |
| `selfcheck/no_opinion.wat` (always answers -1) | 912 | calls 22, ok 0, declined 22 |
| `selfcheck/traps.wat` (traps every call) | 912 | calls 22, traps 22 |

The bar is 0.8 x the node count without a component: the reference passes at 0.59, a component with no
opinion misses it at 1.0 (and fails on `ok == 0` as well). `tasks/extend_wasm/selfcheck/` holds
the grader's fixtures — the reference rule and the three components that must fail, one of which
cannot even be loaded; the prompt never points at them.

## Run

    cargo build --release -p rusolve-cli
    python3 evals/graders.py --selfcheck                       # graders: reference passes, wrong answer fails
    python3 evals/run.py --tasks solve_report --model haiku    # one task, one run
    python3 evals/run.py --all --model haiku --mode cli        # Bash on rusolve-cli instead of MCP
    RSV_PY=/path/to/venv/bin/python python3 evals/run.py --mode python --model haiku   # word problems
    python3 evals/run.py --mode python --model sonnet          # second model: any `claude --model` value

`RSV_CLI=/path/to/rusolve-cli` overrides the binary. The agent is started as
`claude -p <prompt> --model <m> --mcp-config <cfg> --strict-mcp-config --allowedTools "mcp__rusolve__*"
--permission-mode dontAsk --output-format stream-json --verbose --no-session-persistence
--setting-sources ""` in an empty directory, so only the rusolve tools are available and no user or
project instructions leak in. Results go to `evals/results/<stamp>-<model>-<mode>.json`; a run is
one sample, so repeat it (and pin the model) before drawing conclusions. A `--max-budget-usd`
cap applies per task.

The answer protocol: the agent ends its message with `ANSWER: {json}`; the shape is in each
`tasks/<id>/task.json` prompt. A missing or unparsable answer is a failure, and a grader bug is
reported as a failure, never a pass.

`python3 evals/gen_models.py` regenerates the model files (seeded; the committed files are its
output).

`--mode python` starts the agent with `Bash(<RSV_PY> *)`, `Write`, `Edit` and `Read` only (no rusolve tools, no
network): it writes code with the package, runs it, and answers. `RSV_PY` is an interpreter with
`rusolve` installed (`maturin develop --release -m crates/rusolve-py/Cargo.toml` in a venv). The
second model is the same command with another `--model`; results files carry the model name.

## highs-mcp comparison

`--mode highs-mcp` runs the same four word problems with the agent connected to the `highs-mcp`
server (github.com/wspringer/highs-mcp, npm package `highs-mcp` 0.3.2, Node >= 16) instead of
rusolve, started by `HIGHS_MCP_CMD` (`HIGHS_MCP_CMD="npx highs-mcp" python3 evals/run.py --mode
highs-mcp`). The grading is identical, since the answer protocol and the ground truth do not
depend on the solver the agent used. `highs-mcp` 0.3.2 exposes one tool, `optimize-mip-lp-tool`.
The word-problem `task.json` files list `highs-mcp` in `modes`. A run is one sample per task, so
compare tool calls and cost across modes, and repeat runs before ranking the two.

## CI

`.github/workflows/ci.yml` runs fmt, clippy `-D warnings`, the ekk and CLI tests, the
`rusolve-py` pytest suite and `tools/agent_ci.sh` (docs match the binary, every grader accepts its
reference and rejects its wrong answers; no LLM). It does not run the parity check
(`tools/bench/identity.sh`): that needs MIPLIB data, which is not in the repository, so contributors
run it locally (`docs/skills/contrib-identity-abtest/SKILL.md`). For the same reason the CLI tests
run with `--skip sp150x300d` (six tests need that gitignored MIPLIB instance). The evals themselves
are never run in CI (they call an LLM); run them by hand (`results/<stamp>-<model>-<mode>.json`
records the commit it ran on).

## Playbook skills

Each `docs/skills/<class>/SKILL.md` has runnable `rusolve.Model` code checked against an independent
answer (brute force, a DP, a grid). `python3 evals/skill_code.py [skill ...]` executes those blocks
(skips a block fenced as ```python no-run) and exits 1 on a failure; it needs the Python package and
`RSV_CLI` for blocks that call the CLI. The `class_*` tasks above are the agent-facing eval for the same skills.
