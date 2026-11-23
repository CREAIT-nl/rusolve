# rusolve documentation

## Using rusolve

| Document | What it covers |
|---|---|
| [`../README.md`](../README.md) | Install, quick start in Python and on the command line |
| [`agent/reference.md`](agent/reference.md) | Every command, MCP tool, option and environment variable (generated) |
| [`agent/diagnostics.md`](agent/diagnostics.md) | Every error and warning code, with fixes (generated; also `rusolve-cli --explain RSnnnn`) |
| [`design/python-api.md`](design/python-api.md) | The Python `rusolve.Model` API |
| [`design/diagnostics.md`](design/diagnostics.md) | The diagnostics contract: codes, JSON shape, exit codes |
| [`design/mcp.md`](design/mcp.md) | The MCP server (`rusolve-cli mcp`) and its tools |
| [`design/hooks.md`](design/hooks.md), [`design/wasm-hooks.md`](design/wasm-hooks.md) | Custom cut separators, branching rules and node selectors, in Python or sandboxed WebAssembly |
| [`schema/`](schema) | JSON Schemas for every `--format json` output |

## For AI agents

| Document | What it covers |
|---|---|
| [`../llms.txt`](../llms.txt) | Short entry point for agents using rusolve |
| [`skills/start-here`](skills/start-here/SKILL.md) | Classifies a problem and routes to a modelling playbook |
| `skills/` playbooks | [solve](skills/solve/SKILL.md), [debug-infeasible](skills/debug-infeasible/SKILL.md), [tune](skills/tune/SKILL.md), [extend](skills/extend/SKILL.md), [assignment-transport-flow](skills/assignment-transport-flow/SKILL.md), [set-cover-packing](skills/set-cover-packing/SKILL.md), [facility-location](skills/facility-location/SKILL.md), [lot-sizing](skills/lot-sizing/SKILL.md), [scheduling](skills/scheduling/SKILL.md), [routing](skills/routing/SKILL.md), [piecewise-semi](skills/piecewise-semi/SKILL.md), [convex-qp](skills/convex-qp/SKILL.md) |
| [`../evals/`](../evals) | Tasks and graders that measure how well agents use rusolve |

## Contributing

| Document | What it covers |
|---|---|
| [`../CONTRIBUTING.md`](../CONTRIBUTING.md) | Setup, checks, the parity check, rules for solver code |
| [`../LICENSING.md`](../LICENSING.md) | What the license allows, in plain questions and answers |
| [`../AGENTS.md`](../AGENTS.md) | The same rules for AI coding agents, plus where to start |
| [`agent/codemap.md`](agent/codemap.md) | HiGHS source file → Rust file, option → files that read it (generated) |
| `skills/contrib-*` | Recipes: [port a HiGHS function](skills/contrib-port-function/SKILL.md), [add an option](skills/contrib-add-option/SKILL.md), [add a presolve rule](skills/contrib-add-presolve-rule/SKILL.md), [parity check and A/B](skills/contrib-identity-abtest/SKILL.md) |
| [`literature.md`](literature.md) | The papers and open-source solvers rusolve draws on |
| [`../THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md) | Licenses of ported and bundled code |
| `crates/*/AGENTS.md` | Per-crate and per-subsystem notes: contents, invariants, how to test |

Generated files are never edited by hand: `python3 tools/gen_agent_docs.py` writes
`agent/reference.md` and `agent/diagnostics.md`, and `python3 tools/gen_codemap.py` writes
`agent/codemap.md`. CI fails when they are stale.
