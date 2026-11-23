# rusolve-cli

The `rusolve-cli` binary: command registry, JSON envelope, MCP server, problem-class router, diagnostics, tuning and the glue that calls the HiGHS port.

## What lives here
- `src/main.rs`: dispatch; `src/commands.rs`: the command registry that drives `--help --json`, flags and `docs/agent/reference.md`.
- `src/envelope.rs`: the `--format json` output; schemas in `docs/schema/`.
- `src/mcp.rs`, `src/mcp_tools.rs`: the MCP server (`rusolve mcp`); `src/env_registry.rs`: `RUSOLVE_*` variables.
- `src/route.rs`, `src/quad_convexity.rs`, `src/qp.rs`, `src/gc_solve.rs`, `src/conic_solve.rs`: problem-class routing and engines.
- `src/mip_highs.rs`: the MIP path into the port, including the oracle-pinned option sets; `src/ekk.rs`, `src/ekk_seq.rs`: LP path.
- `src/check.rs`, `src/explain.rs`, `src/diagnose.rs`, `src/lint.rs`, `src/tune.rs`, `src/tune_race.rs`, `src/options_cmd.rs`, `src/profile.rs`.
- `tests/`: one file per feature (`mip_trace.rs`, `agent_docs.rs`, `envelope.rs`, `mcp.rs`, ...).

## Invariants
- Shared files (`src/commands.rs`, `src/main.rs`, `src/envelope.rs`, `docs/schema/*`) get additive edits only.
- Logs go to stderr; `--format json` prints exactly one object on stdout with a `schema_version`.
- Exit codes: 0 ran, 2 usage, 3 input, 4 internal.
- Default search unchanged: identity 20/20 (`tools/bench/identity.sh`); new levers are default off and pinned off in the oracle option sets of `src/mip_highs.rs`.
- `docs/agent/reference.md` and `docs/agent/codemap.md` are generated; never hand-edit.

## How to test
- `cargo test --release -p rusolve-cli --no-fail-fast`, then `python3 tools/gen_agent_docs.py --check` (needs the built binary, `RSV_CLI`).
- Before a merge, run locally: `cargo fmt --all --check`, `cargo clippy -p rusolve-ekk -p rusolve-cli -p rusolve-wasm --release --tests --benches -- -D warnings`, `cargo test --release -p rusolve-ekk --lib`, the crate tests above, and `tools/bench/identity.sh`.

## Change recipes
1. Add a command: add its entry to the registry in `src/commands.rs`, the handler arm in `src/main.rs`, a schema in `docs/schema/`, a test in `tests/`, regenerate the reference.
2. Add an MCP tool: define it in `src/mcp_tools.rs`, cover it in `tests/mcp.rs`.
3. Add an option: follow `docs/skills/contrib-add-option/SKILL.md`; the CLI part is the oracle pin in `src/mip_highs.rs`.

## Pitfalls
- A new flag must be in the registry or the option parser rejects it with exit 2, and the generated reference goes stale.
- The byte scan in `src/route.rs` keeps plain LP and MIP files on the unchanged legacy path; a router change needs identity 20/20.
- Further reading: `docs/design/mcp.md`, `docs/design/diagnostics.md`.
