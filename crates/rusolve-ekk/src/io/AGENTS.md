# io

HiGHS logging: the minimal subset of `HighsIO` that the port needs, not file formats (readers and writers are in `crates/rusolve-io`).

## What lives here
- `highs_io.rs`: `HighsLogOptions`, `HighsLogStream`, the log-type levels and the printing helpers.
- `mod.rs`: the module doc and re-exports.
- Related: `src/util/c_printf.rs` formats numbers the way C `printf` does, which log lines and traces depend on.

## HiGHS files mirrored
`io/HighsIO.h` and `io/HighsIO.cpp` (subset). Table: `docs/agent/codemap.md`.

## Invariants
- Log text is part of the contract: traces, the progress stream and several tests compare lines, so message wording and number formatting must match HiGHS byte for byte.
- Output flags (`output_flag`, `log_to_console`) are options in `src/lp_data/highs_options.rs`; callbacks can capture the log (the MIP trace tests rely on that).
- Logging never changes the search; no lever belongs here. Logs go to stderr in the CLI.

## How to test
- `cargo test --release -p rusolve-ekk --lib io`; anything that compares log lines (`crates/rusolve-cli/tests/mip_trace.rs`, `crates/rusolve-cli/tests/progress.rs`).
- Entry points are `highs_log_user` and `highs_log_dev` in `highs_io.rs`; number formatting helpers there are `highs_double_to_string` and `highs_bool_to_string`.

## Change recipes
1. Add a log message: translate the HiGHS call with its log type, check the exact format string and spacing against the C++, and run the tests that compare logs.
2. Add a log destination: add it beside `HighsLogStream`, default to the current behavior, test that default output is unchanged.
3. Port a missing `HighsIO` function: `docs/skills/contrib-port-function/SKILL.md`.

## Pitfalls
- Dev-level messages are gated by `log_dev_level`; do not print from a hot loop unless `highs_log_dev_enabled` says so.
- `HighsLogOptions` carries `output_flag`, which defaults to true; `HighsLogOptions::quiet()` builds a silent one for internal calls.
- Further reading: the module doc in `mod.rs` and the HiGHS file map in `docs/agent/codemap.md`.
