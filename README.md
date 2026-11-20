# rusolve

An LP and MIP solver in Rust, built to be used by people and by AI agents.

rusolve starts from a faithful Rust port of [HiGHS](https://highs.dev): its dual simplex,
interior point, presolve and branch-and-cut are translated function by function and checked
against HiGHS itself. On top of that base it adds techniques from the optimization literature
and from other permissively licensed solvers (SCIP, Clarabel, and others, credited in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)). It is also designed for AI agents: every
command can return JSON, every error has a code and an explanation, and an MCP server is built
in.

**Status:** early (0.1.0), not yet on PyPI. LP and MIP are the core. Convex quadratic
objectives and semi-continuous variables are supported. Interfaces may still change.

## Install

From source (works today; needs a Rust toolchain):

    git clone https://github.com/CREAIT-nl/rusolve && cd rusolve
    cargo build --release -p rusolve-cli          # the solver binary: target/release/rusolve-cli

Python package (3.11 or newer; the wheel bundles the CLI and the MCP server):

    pip install maturin && tools/build_wheel.sh    # writes the wheel to dist-wheels/
    pip install dist-wheels/rusolve-*.whl

After the first release this becomes `pip install rusolve`, or `uvx --from rusolve rusolve mcp`
to run the MCP server without installing.

The CLI is called `rusolve-cli` when built with cargo and `rusolve` when installed from the
wheel. The examples below use `rusolve`.

## Quick start

Python:

```python
import rusolve as rs

m = rs.Model()
x = m.add_vars(3, binary=True)
m.add_constr(x[0] + x[1] + x[2] <= 2)
m.maximize(5 * x[0] + 4 * x[1] + 3 * x[2])
sol = m.solve(time_limit=60)
print(sol.status, sol.objective)    # Optimal 9.0
```

More in [`crates/rusolve-py/README.md`](crates/rusolve-py/README.md): bulk data and sparse
matrices, warm re-solves after changes, callbacks, and SciPy- and highspy-compatible layers.

Command line (MPS and LP files, optionally compressed):

    rusolve mip model.mps --format json --solution       # solve, result as one JSON object
    rusolve explain model.mps                            # why is it infeasible? (an IIS, by name)
    rusolve feasrelax model.mps --mode min-count          # the cheapest repair
    rusolve check model.mps solution.json                # which constraints a solution violates
    rusolve lint model.mps                               # big-M values, scaling, modelling smells
    rusolve options                                      # every option, its range and default

`rusolve --help` lists all commands.

## Diagnostics

Problems with a model or a command come back as coded diagnostics, in the style of `rustc`:

    warning[RS2001]: big-m: row 'c1': coefficient -1e6 on binary 'b' is far larger in magnitude than the tightest valid value -5
     --> row 'c1'
      = row: c1: x - 1e6 b <= 0
      = help: set the coefficient of b in c1 to -5 and the upper bound to 0 (was -1000000 and 0)

`rusolve --explain RS2001` explains any code. With `--format json`, the same diagnostics are
part of the JSON output. The list of codes is in
[`docs/agent/diagnostics.md`](docs/agent/diagnostics.md).

Conventions: `--format json` prints exactly one JSON object on stdout (logs go to stderr);
schemas are in [`docs/schema/`](docs/schema). Exit codes: 0 ran (the solve status is in the
output), 2 usage error, 3 input error, 4 internal error. A misspelled option or `RUSOLVE_*`
variable is an error that suggests the nearest valid name.

## For AI agents

- **MCP server** (stdio): `claude mcp add rusolve -- /path/to/rusolve-cli mcp`, or in any MCP
  client `{"mcpServers": {"rusolve": {"command": "/path/to/rusolve-cli", "args": ["mcp"]}}}`.
  Tools: `load_model`, `build_model`, `solve`, `get_solution`, `modify`, `check_solution`,
  `explain`, `lint`, `ranging`, `feasrelax`, `free_model`. Design:
  [`docs/design/mcp.md`](docs/design/mcp.md).
- **Entry point:** [`llms.txt`](llms.txt). Full generated reference of commands, MCP tools,
  options and environment variables: [`docs/agent/reference.md`](docs/agent/reference.md).
- **Skills:** task recipes in [`docs/skills/`](docs/skills) (start with `start-here`, then
  `solve`, `debug-infeasible`, `tune`, `extend`, and playbooks for common model classes such as
  scheduling, routing, lot sizing and facility location). The folder is also a Claude Code plugin
  (`.claude-plugin/`).
- **Evaluations:** [`evals/`](evals) measures how well agents use all of this.

## Documentation

[`docs/README.md`](docs/README.md) is the index.

## Contributing

Contributions are welcome, including AI-assisted ones. See [`CONTRIBUTING.md`](CONTRIBUTING.md).
Agents start at [`AGENTS.md`](AGENTS.md).

## License

Copyright 2026 CREAIT. rusolve is source-available under the [Business Source License 1.1](LICENSE).
In short:

- **Free:** reading, modifying, forking and redistributing the code; research, teaching,
  testing and evaluation; and production use to solve your own (or your company's) optimization
  problems.
- **Needs a commercial license:** using rusolve to provide optimization to others, as a
  product, a hosted service or consulting. Contact [CREAIT](https://creait.nl).
- **Becomes open source:** each release converts to the GNU AGPL v3 four years after it is
  published.

[`LICENSING.md`](LICENSING.md) answers the common questions (can my company use it in
production? can we sell a product built on it?). The [license text](LICENSE) is what applies.

rusolve contains code translated from HiGHS (MIT) and other permissively licensed projects; their notices
are in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
