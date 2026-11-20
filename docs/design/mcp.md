# rusolve MCP server and model handles

`rusolve mcp` is a Model Context Protocol server on stdio: newline-delimited
JSON-RPC 2.0 on stdin/stdout, nothing but protocol on stdout (solver logging is
forced off), diagnostics on stderr. There is no HTTP transport.

## Client configuration

    claude mcp add rusolve -- /path/to/rusolve-cli mcp

or in `.mcp.json`:

    {"mcpServers": {"rusolve": {"command": "/path/to/rusolve-cli", "args": ["mcp"]}}}

## Protocol

Hand-rolled rather than `rmcp`: the server is five methods over a blocking
stdin loop, and `rmcp` brings an async runtime (tokio) the solver does not
otherwise have, for no behavior we need. Both spec eras are served:

- legacy `initialize` / `notifications/initialized` (the requested version is
  echoed if supported, otherwise the latest supported is returned), `ping`,
  `tools/list`, `tools/call`
  (https://modelcontextprotocol.io/specification/2025-11-25/basic/lifecycle#version-negotiation);
- the stateless 2026-07-28 era: `server/discover`, a per-request `_meta`
  protocol version, `resultType: "complete"`, and error -32022 for an
  unsupported version
  (https://modelcontextprotocol.io/specification/2026-07-28/basic/versioning,
  https://modelcontextprotocol.io/specification/2026-07-28/server/discover).

Tool results carry `structuredContent` and the same JSON as a text item. A
tool failure is `isError: true` with a message naming the bad field
(https://modelcontextprotocol.io/specification/2026-07-28/server/tools).
Unknown method is -32601, unknown tool or bad params -32602, bad JSON -32700.

## Handles

`load_model {path | mps_text | lp_text | mof_text}` returns a `model_id`. Each handle owns one live
ported `Highs` (its `HighsLp` and last basis). There is no global current
model. `free_model` releases a handle.

| tool | does |
|---|---|
| load_model | read MPS from a path or text; size and name summary |
| solve | `rusolve.result/1`; `options` (flat map through the HiGHS option setters), `time_limit`, `warm_start`, `include_solution` |
| get_solution | column values by name, paged |
| modify | `col_bounds`, `row_bounds`, `cost`, `coef`, `add_row`, `add_col`, `del_row`, `del_col`, by name |
| check_solution, explain, lint, ranging, feasrelax | the subcommand cores; the same JSON as `--format json` |
| free_model | release the handle |

After `modify`, an LP re-solve starts from the kept basis (HiGHS skips presolve
when a valid basis exists); `solve {warm_start: false}` forces a cold solve for
comparison. A solve of an edited model reports `input.path = "<id> (edited)"`,
`input.sha256` of a canonical dump of the model, and `input.edits_applied`.
For an LP, `stats.nodes` is clamped from HiGHS's -1 to 0 so the result
validates against its schema (the CLI still prints the raw value).

## Models an agent can write

`load_model` takes exactly one of `path`, `mps_text`, `lp_text` (CPLEX LP) or `mof_text`
(MathOptFormat JSON). The text goes to a temp file of the matching extension and through the
same reader as a path, so the format logic exists once (`rusolve-io`).

`build_model` takes JSON and returns the same handle summary:

    {"sense": "maximize",
     "variables": [{"name": "x", "type": "binary", "objective": 5},
                   {"name": "f", "lower": null, "upper": 3}],
     "rows": [{"name": "c1", "terms": {"x": 1, "f": 1}, "sense": "<=", "rhs": 4},
              {"name": "c2", "terms": {"x": 1}, "lower": 1, "upper": 6}],
     "objective": {"terms": {"f": 2}, "offset": 5}}

`mcp_build.rs` writes it as free-format MPS (OBJSENSE, MARKER lines around integers, RANGES,
explicit bounds for every column; <https://lpsolve.sourceforge.net/5.5/mps-format.htm>) and loads
that. A bad entry is a tool error naming its JSON path (`rows[2].terms.z`). Quadratic terms and
general constraints are written as LP text instead.

### Routing

A model that may hold a quadratic term, a cone or a general constraint is classified once by
`route::classify_model`. The summary carries `class` and `engine`. A plain LP or MIP is the
unchanged handle (live `Highs`, warm re-solves, every tool). A convex QP, a conic model or a
general-constraint model gets a routed handle: `solve` runs `rusolve-cli solve <file> --format
json --include-solution` as a child process, so the answer is the CLI's and a cancel can kill it.
`solve` (with `time_limit`, `include_solution`), `get_solution` and `free_model` work on it; the
tools that need the live LP (`modify`, `check_solution`, `explain`, `lint`, `ranging`,
`feasrelax`) answer with the class and the tools that do work, and `options`, `warm_start` and
the wasm arguments are refused rather than ignored. A class nothing solves (nonconvex, MIQP,
mixed-integer conic) is refused at load with the class and the reason.

## Prompts and resources

The skills in `docs/skills/*/SKILL.md` are `prompts/list` and `prompts/get`
(<https://modelcontextprotocol.io/specification/2026-07-28/server/prompts>): name = directory,
`title`/`description` from the front matter, one user message with the body. The schemas in
`docs/schema/*.json` (`rusolve://schema/<name>`) and `docs/agent/reference.md`
(`rusolve://docs/agent/reference`) are `resources/list` and `resources/read`
(<https://modelcontextprotocol.io/specification/2026-07-28/server/resources>); an unknown URI is
-32002. `build.rs` lists the directories and embeds the files (`$OUT_DIR/mcp_docs.rs`), so a
skill or schema added to the tree appears at the next build without touching Rust. The
capabilities advertise `prompts` and `resources`.

## Progress and cancellation

A `tools/call` with `_meta.progressToken` gets `notifications/progress` (at most five a second,
`progress` strictly increasing, the message naming nodes, LP iterations and the bound)
<https://modelcontextprotocol.io/specification/2026-07-28/basic/utilities/progress>.
`notifications/cancelled {requestId}` stops a running or queued call and, as the spec asks, gets no
response <https://modelcontextprotocol.io/specification/2026-07-28/basic/utilities/cancellation>.

A solve blocks its thread, so a reader thread owns stdin: it handles the cancel notification
itself and queues every other line. A handle's `solve` installs the port's own interrupt
callbacks for that solve only (`kCallbackMipInterrupt`, `kCallbackSimplexInterrupt`, HiGHS
`HighsCallbackInput::user_interrupt`) and removes them after; they read the cancel flag and emit
progress, never steer the search. A routed solve polls the flag and kills the child.

Tests (`tests/mcp.rs`): LP, MOF and JSON models load and solve; a convex QP and a
general-constraint model solve and page through `get_solution`; every skill is a prompt and every
schema a resource; progress is monotonic; a cancelled solve ends within seconds with no reply; a cancelled routed solve (a general-
constraint MIP with a 60 s limit) is killed, the ping behind it answers within five seconds and
nothing arrives for the solve's id.

## Tool errors

A failed tool call is `isError: true` whose text is the `rusolve.diagnostics/1` document of the diagnostic
type (`docs/design/diagnostics.md`), the same JSON the CLI writes under `--format json`. The code
follows the failure (`mcp.rs` `error_diagnostic`):

| code | when |
|---|---|
| RS0010 | `load_model` `path` names a file that cannot be opened |
| RS0011 | a model that does not parse (LP, MPS or MathOptFormat text or file), or a `build_model` entry that is invalid |
| RS0012 | a model class nothing solves (nonconvex, MINLP) |
| RS3001 | a missing or invalid tool argument, naming the field |
| RS3002 | a `model_id` that is not loaded |
| RS0090 | a panic inside a tool, caught so the server and the other handles survive |

RS0011 carries a `location` of kind `file`, taken from the reader's own error text: LP and MPS
readers say `at line N` (`lp.rs`, `mps.rs` in `rusolve-io`), MathOptFormat says `at line N column M`
for invalid JSON. `path` is the file for `path`, or `<lp_text>`, `<mps_text>`, `<mof_text>` for
inline text; `source_line` is read back from there. A reader that gives no column gets the first
non-blank one with a one-character caret. A reader error that names no line (a bad MOF field, a
`build_model` entry) has no location; the message names the JSON path instead. This reads only the
message in `mcp.rs`; when `rusolve_io` grows structured reader diagnostics the same
`Location` is filled from them.

There is no code for a cancelled call: a cancelled request gets no reply (MCP cancellation spec),
so there is no such code. Test:
`mcp_load_model_errors_carry_the_reader_line` (bad LP text, a bad MPS file and the same text
inline, a missing file, bad MOF JSON).

## Open

HTTP transport; option schema generated from the option registry; tasks extension, status and
cancel tools; `resource_link` for large solutions; concise/detailed switch;
list-options and features tools.
