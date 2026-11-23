# Diagnostics: the API contract

Every crate, the CLI, the MCP server and the Python package report errors and warnings through the one type below.
Source of the layout: the rustc diagnostic guide (<https://rustc-dev-guide.rust-lang.org/diagnostics.html>) and
error-code index (<https://doc.rust-lang.org/error_codes/>); did-you-mean follows rustc's
`find_best_match_for_name` (distance at most `max(len, 3) / 3`) over optimal string alignment distance.

## Where it lives

`crates/rusolve-core/src/diagnostic.rs` (`rusolve_core::diagnostic`). `rusolve-core` has no dependencies and
every other crate already depends on it, so `rusolve-io`, `rusolve-model`, `rusolve-cli` and `rusolve-py` use it
without a cycle. The CLI funnel is `crates/rusolve-cli/src/diag.rs`.

## The type

```rust
use rusolve_core::diagnostic::{Diagnostic, Location, Severity};

let d = Diagnostic::error("RS1004", "row 'R7' is not declared")
    .with_location(Location::File { path, line, column, len, source_line })
    .with_note("ROWS declared 3 rows")
    .with_help("did you mean 'R1'?");
```

| field | meaning |
|---|---|
| `code` | `&'static str`, `RSnnnn`, must be in `REGISTRY` |
| `severity` | `Error` (fails the command), `Warning`, `Note`, `Help` |
| `message` | one line, lower case start, no trailing period |
| `location` | optional: `File {path, line, column, len, source_line}` (1-based), `Row {name, rendered}`, `Column {name, rendered}`, `Python {path, line, source_line}` |
| `notes` | facts (`= note:`) |
| `help` | concrete fixes (`= help:`), including did-you-mean |

Code ranges: `RS0xxx` command line, files and the process; `RS1xxx` model readers;
`RS2xxx` and up lint and solve warnings. To add a code: add a `CodeEntry` to
`diagnostic::REGISTRY` (title plus an explanation of at least a few sentences: what it means, why, the fix), run
`python3 tools/gen_agent_docs.py` and commit `docs/agent/diagnostics.md`.

Two tests keep the registry honest (`rusolve-core`, `diagnostic::tests`): every `"RSdddd"` string literal anywhere
in `crates/**/*.rs` must be registered, and every registered code needs a title, a long explanation and a
`### RSnnnn` section in `docs/agent/diagnostics.md` (which `tools/gen_agent_docs.py --check` also keeps current).
Write a code as a plain string literal at the emit site so the scan sees it.

## Renderers (one value, three outputs)

- `d.render_text()` rustc-style text for stderr (snapshot test `text_renderer_snapshot`).
- `d.to_json()` / `diagnostics_json(&[d], exit_code)` the `rusolve.diagnostics/1` document
  (`docs/schema/rusolve-diagnostics-1.schema.json`).
- `d.mcp_error_content(exit_code)` the same document, as the text content of an MCP result with `isError: true`.
  Python wraps the same JSON in `.diagnostics` on `rusolve.ModelError` and friends.

## CLI funnel (`crates/rusolve-cli/src/diag.rs`)

- `diag::wants_json(args)` reads `--format json` before the command parser runs.
- `diag::fail(json, diagnostic, exit)` prints the text to stderr and, under `--format json`, the document to
  stdout, then returns the exit code (2 usage, 3 input, 4 internal; the exit-code table). Use it instead of
  `eprintln!` + `ExitCode::from(n)` on every error path; `fail_msg`, `fail_read` (RS0010/RS0011),
  `fail_usage` (a `commands::UsageError`, which carries its code and help), `unknown_command` (RS0001) and
  `unknown_option` (RS0004) are the common shapes.
- `commands::parse_flags` errors carry `code` and `help`; `UsageError::diagnostic()` converts.
- `rusolve --explain RSnnnn` prints the registry entry (exit 0); `--explain all` prints the registry as JSON
  (the source of `docs/agent/diagnostics.md`); an unknown code is RS0003 with the closest codes as help.

## Did you mean

`diagnostic::did_you_mean(name, candidates)` and `suggestions(name, candidates, n)` serve the registries:
command names (RS0001), a command's flags (RS0002), option names (RS0004, via the option registry's
`unknown_option_hint`), `RUSOLVE_*` variables (RS0005, `env_registry::suggest`). Row and column names
(RS1004) use the same functions over the declared names in the readers.

## Coverage

Every error path of `solve`, `mip` (both stacks, including the native-stack flag parser under
`RUSOLVE_MIP_HIGHS=0`), `convert`, `check`, `lint`, `explain` (and `explain --mip-iis`), `feasrelax`, `ranging`,
`classify`, `probe`, `diagnose`, `tune`, `profile`, `options` and the router writes one `rusolve.diagnostics/1`
document under `--format json` (the test `every_error_path_prints_one_diagnostics_document_under_json` in
`crates/rusolve-cli/tests/diagnostics.rs` runs a missing file, a bad flag value and a wrong positional count of each
and validates the document, the exit code and that every code is registered). Codes: RS0010 missing file, RS0011 a
file that does not parse (a solution file in `check`, a profile that is not JSON), RS0002 unknown flag, RS0003 bad
flag value, RS0006 wrong positional count, RS0012 an input the command cannot honour, RS0090 internal. The shared
helpers are `diag::fail_read_str`, `fail_usage_text`, `fail_core` (the `CoreError` of explain, ranging and
feasrelax) and `fail_coded` (the `(exit, message)` pairs of profile and tune); a command derives `json` from its
parsed arguments, or `wants_json(args)` when the parse itself failed. A refused class (`route::refuse`) prints its
`rusolve.error/1` object with an additive `diagnostics` array holding one RS0012. The conic safety net in `main.rs`
re-runs a `solve`/`mip` the legacy reader refused for conic content; the first attempt's JSON read error is held back
(`diag::defer_conic_json`) so the command still ends with one document. Remaining plain-text sites: the
`mip-trace`/`presolve-dump`/`ekk-seq`/`dump` developer commands, the QP front's own `rusolve.error/1` refusals
(reached only with the router bypassed), and progress or warning lines of `tune` (stderr by design). To convert a
site, replace `eprintln!(..); return ExitCode::from(n)` with `return diag::fail_msg(json, code, n, msg)`.

A registered code must be emitted: the tests `every_registered_code_is_emitted` (cli) and
`every_registered_code_is_emitted_somewhere` (core) scan the sources for `"RSdddd"` outside the registry's
`code:` lines and test modules.

## Post-solve diagnostics budget

`adopt::budget_for(solve_time, user_time_limit)` gives the relaxation re-solve and the IIS one shared budget:
`max(0.25 s, half the solve time)` (a 1.2 s solve costs at most 0.6 s more), at most 10 s, at most what is left of the user's `--time-limit`. The status is claimed
(`Tracker::claim_final`) before this work, so a watchdog hard stop cannot turn a proved status into TimeLimit.
"The LP relaxation is feasible" is claimed only for a relaxation status of Optimal or Unbounded; any other status
adds a note that the budget ran out and points to `rusolve explain`.
