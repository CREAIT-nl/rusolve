---
name: rusolve-contrib-add-option
audience: contributor
description: Add a solver option to rusolve-ekk and wire it everywhere it must appear (field, registry record, setter mapping, coverage test, oracle pins, generated docs). Use when asked to add, expose or gate a new option, switch or tuning parameter.
---

# Add an option

An option is a field of `HighsOptions` plus a registry record. An option for a technique HiGHS does not
have must default to off, name its source, and leave the default search untouched (the parity check
`tools/bench/identity.sh` prints `IDENTICAL 20 of 20`). Worked example to
copy: `mip_dff_separator`, which touches every file below.

## Wiring list

1. `crates/rusolve-ekk/src/lp_data/highs_options.rs`: add `pub <name>: <type>` to `HighsOptions` and its
   default in `impl Default for HighsOptions`. `recommended()` is the CLI default profile; leave a new
   option off there unless a measured, held-out win says otherwise.
2. `crates/rusolve-ekk/src/lp_data/highs_option_records.rs`, four places:
   - an `OptionRecord` in `OPTION_RECORDS` (`name`, `description`, `advanced`, `kind` with its default
     and range, `source` naming a paper or "HiGHS 1.16.0", `affects`, `tunable`);
   - the mapping line `"<name>" => <field> : OptionValue::<Bool|Int|Double|Str>;` in the
     `option_field_accessors!` invocation;
   - `<field>: _,` in the destructure of the coverage test;
   - the count in `assert_eq!(OPTION_RECORDS.len() + INTERNAL_FIELDS.len(), N)` (a field that is not a
     public option goes into `INTERNAL_FIELDS` instead).
3. `crates/rusolve-ekk/src/lp_data/highs_api_mip.rs`: add set and get arms only if the typed setters need a
   special case; plain fields work through the mapping line.
4. Read it where it acts, as `options.<name>` and guarded so the old path is byte-identical when it is off
   (for a separator: `crates/rusolve-ekk/src/mip/separation.rs`).
5. If the default profile turns it on, pin it off in the HiGHS-oracle option sets of
   `crates/rusolve-cli/src/mip_highs.rs` (`trace_options` and `run_mip_trace`), which must reproduce the
   1.15.1 oracle trace.
6. Regenerate the generated docs after `cargo build --release -p rusolve-cli`:
   `python3 tools/gen_agent_docs.py` and `python3 tools/gen_codemap.py`.

## Check

- `python3 tools/gen_agent_docs.py --check` and `python3 tools/gen_codemap.py --check` pass.
- `cargo test --release -p rusolve-ekk --lib highs_option_records` passes (the coverage test).
- Prove the default search did not move: `docs/skills/contrib-identity-abtest/SKILL.md`.
- Run the CI subset:
  `cargo fmt --all --check`,
  `cargo clippy -p rusolve-ekk -p rusolve-cli -p rusolve-wasm --release --tests --benches -- -D warnings`,
  `cargo test --release -p rusolve-ekk --lib`, `cargo test --release -p rusolve-cli`.

## Common mistakes

- Forgetting the count in the coverage test, or the destructure line: the test names the missing field.
- A new option that defaults on without a held-out measurement. Keep it off.
- Editing `docs/agent/reference.md` by hand: it is generated.
