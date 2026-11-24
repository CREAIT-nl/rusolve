# AGENTS.md

Instructions for AI coding agents working in this repository. Humans: the same rules are in
[`CONTRIBUTING.md`](CONTRIBUTING.md), which also applies to you, the agent.

## What this is

rusolve is an LP and MIP solver in Rust. Its default engine (`crates/rusolve-ekk`) is a
faithful port of HiGHS, improved with techniques from published papers and other permissively
licensed solvers. Everything is also exposed to agents: JSON output, coded diagnostics, an MCP
server.

If you only want to *use* rusolve to solve a model, install it with [`INSTALL.md`](INSTALL.md), then read [`llms.txt`](llms.txt) and
[`docs/skills/start-here`](docs/skills/start-here/SKILL.md) instead.

## Before you change anything

1. Read the `AGENTS.md` in the directory you are changing. Every crate and every subsystem
   directory of `crates/rusolve-ekk/src` has one: what lives there, the HiGHS files it mirrors,
   its invariants, how to test it, and recipes for common changes.
2. [`docs/agent/codemap.md`](docs/agent/codemap.md) maps every HiGHS source file to its Rust
   file and every option to the files that read it.
3. Recipes for the common contributions are in `docs/skills/`:
   [`contrib-port-function`](docs/skills/contrib-port-function/SKILL.md),
   [`contrib-add-option`](docs/skills/contrib-add-option/SKILL.md),
   [`contrib-add-presolve-rule`](docs/skills/contrib-add-presolve-rule/SKILL.md),
   [`contrib-identity-abtest`](docs/skills/contrib-identity-abtest/SKILL.md).

## Layout

| Path | What it is |
|---|---|
| `crates/rusolve-ekk` | The HiGHS port: simplex, presolve, MIP, interior point, QP, the option registry. Default engine. |
| `crates/rusolve-cli` | The `rusolve-cli` binary: commands, JSON output, MCP server, diagnostics, lint, tuning. |
| `crates/rusolve-py` | Python bindings (`rusolve.Model`), built with maturin. |
| `crates/rusolve-model` | The append-only model builder behind the Python API. |
| `crates/rusolve-core` | Shared data structures and the diagnostic code registry. |
| `crates/rusolve-io` | MPS and LP readers and writers. |
| `crates/rusolve-lp`, `crates/rusolve-mip` | The older native LP and MIP engines (not the default). |
| `crates/rusolve-wasm` | User-written MIP hooks as sandboxed WebAssembly. |
| `docs/` | User and contributor docs; index in `docs/README.md`. Generated files are listed there. |
| `tools/` | Doc generators, the HiGHS oracle (`tools/highs-oracle`), differential tests against HiGHS (`tools/diff`), benchmark and parity scripts (`tools/bench`). |
| `evals/` | Tasks and graders that measure how well agents use rusolve. |

## Rules

- **Ported code stays faithful.** Files translated from HiGHS keep the same algorithm, constants
  and control flow, and a header naming the HiGHS source file and its MIT notice. A behavioural
  difference from HiGHS is a bug, not an improvement. Improvements go behind an option.
- **The default search does not move by accident.** A change that is not meant to alter results
  (refactor, performance, cleanup, a new option that is off by default) must pass the parity
  check in [`CONTRIBUTING.md`](CONTRIBUTING.md#parity-check). A change that *is* meant to alter
  them needs benchmark evidence in the pull request.
- **Cite the source** of every algorithmic idea taken from elsewhere: a paper, a textbook, public
  documentation, or the source of an open-source solver whose license allows it (HiGHS, SCIP,
  Clarabel, …). Put the citation in a comment next to the code. An original idea needs no
  citation: explain how it works in the comment, and give the benchmark evidence in the pull
  request.
- **Only permissively licensed code** may be translated or copied in (MIT, BSD, Apache-2.0, or
  similar), with its notice in the file header and in `THIRD_PARTY_NOTICES.md`. Never copy from
  GPL code or from proprietary software.
- **Do not edit generated files by hand.** `docs/agent/reference.md` and `docs/agent/diagnostics.md` come from
  `python3 tools/gen_agent_docs.py`, and `docs/agent/codemap.md` comes from
  `python3 tools/gen_codemap.py`. Regenerate them, and commit the result together with the change
  that caused it.
- **Performance is measured, never claimed.** "Should be faster" is not evidence.
- **Keep the change small and in scope.** One pull request does one thing. Do not reformat or
  refactor code you are not otherwise changing.
- **Match the existing patterns.** Look at neighbouring code before introducing a new structure,
  dependency or convention.
- `unsafe` needs a written justification in the pull request.

## Verify before you open a pull request

Run the checks from [`CONTRIBUTING.md`](CONTRIBUTING.md#checks) and say in the pull request which
ones you ran and what they printed. Do not claim a check passed if you did not run it.

## The pull request

Fill in the template. State what changed and why, the source you followed, how you verified
it, and that the change was made with AI assistance. A maintainer reviews every pull request;
you or the person running you must be able to answer review questions about it.
