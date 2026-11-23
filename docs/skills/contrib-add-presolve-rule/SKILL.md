---
name: rusolve-contrib-add-presolve-rule
audience: contributor
description: Add or change a presolve reduction in rusolve-ekk (rule body, postsolve record, rule switch, tests, identity proof). Use when asked to add a presolve rule, a reduction, or a postsolve step.
---

# Add a presolve rule

Presolve lives in `crates/rusolve-ekk/src/presolve/`. Read `crates/rusolve-ekk/src/presolve/AGENTS.md` first.
The order of reductions is part of the faithful port: a new rule that fires in the default path changes the
presolved model and moves the default search, so a rule that HiGHS does not have is a default-off option.

## Steps

1. Decide the kind. A HiGHS rule you have not ported: translate it (`docs/skills/contrib-port-function/SKILL.md`).
   A rule HiGHS does not have: cite the paper it comes from, and gate it with a default-off option
   (`docs/skills/contrib-add-option/SKILL.md`).
2. Find the neighbors. Row rules are in `hpresolve/row_presolve.rs`, column rules in
   `hpresolve/col_presolve.rs`, singleton and doubleton rules in `hpresolve/singletons.rs`, dual rules in
   `hpresolve/dual_fixing.rs`, substitutions in `hpresolve/eq_subst.rs` and `hpresolve/aggregator.rs`.
   The loop that calls them is `hpresolve/presolve_loop.rs`.
3. Write the rule as a method of `HPresolve` in the file of its kind (one `impl` block per contiguous
   range). Return the `Result` the neighbors return, and mark changed rows and columns the way they do.
4. Record the postsolve. Every change that removes or rewrites a row or column pushes a record onto the
   postsolve stack (`crates/rusolve-ekk/src/presolve/highs_postsolve_stack/`) with an `undo` that restores
   primal, dual and basis values. No undo, no rule.
5. Rule switch: HiGHS rule slots are `kPresolveRule*` in `crates/rusolve-ekk/src/lp_data/hconst.rs`, honored
   through `allow_rule_` in `crates/rusolve-ekk/src/presolve/hpresolve_analysis.rs`. A rusolve rule gets its
   own option instead of a new bit.
6. Tests: a unit test next to the rule on a tiny model where it fires (check the presolved model and that
   postsolve restores the solution), plus one where it must not fire.

## Check

- `cargo test --release -p rusolve-ekk --lib presolve` and `--test t_g1_presolve_api`.
- Profile the effect with `RUSOLVE_PRESOLVE_PROF` (`hpresolve/presolve_prof.rs`).
- The parity check (`tools/bench/identity.sh`) prints `IDENTICAL 20 of 20` with the rule off, and an A/B on the affected set with it on:
  `docs/skills/contrib-identity-abtest/SKILL.md`.
- Run the CI subset:
  `cargo fmt --all --check`,
  `cargo clippy -p rusolve-ekk -p rusolve-cli -p rusolve-wasm --release --tests --benches -- -D warnings`,
  `cargo test --release -p rusolve-ekk --lib`, `cargo test --release -p rusolve-cli`.
