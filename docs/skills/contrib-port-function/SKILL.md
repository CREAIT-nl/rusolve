---
name: rusolve-contrib-port-function
audience: contributor
description: Port a HiGHS function into rusolve-ekk test-first, with a golden fixture dumped from the real HiGHS internals. Use when asked to translate, add or fix a HiGHS routine so the port matches HiGHS bit for bit.
---

# Port a HiGHS function, test first

The port is faithful: same algorithm, constants, control flow and data layout; every divergence is a bug
found by trace diff. Rules: the crate root doc in `crates/rusolve-ekk/src/lib.rs` (naming, type mapping,
no FMA) and `crates/rusolve-ekk/AGENTS.md`.

## Steps

1. Locate the C++ in `reference/HiGHS` (read-only; gitignored). Get it with
   `git clone https://github.com/ERGO-Code/HiGHS reference/HiGHS` and
   `git -C reference/HiGHS checkout 73cac48c5340d775a477087198611862559be250` (HiGHS 1.15.1, the pinned
   oracle version), or point `HIGHS_REF` at an existing checkout. Note the file and the line range.
   Find where the Rust belongs with `docs/agent/codemap.md` (HiGHS file to Rust file); a long `.cpp` is a
   directory of Rust files, one per contiguous line range.
2. Fixture first. Add a hook to the oracle in `tools/highs-oracle/src/` (one `.cpp` per fixture set:
   `unit.cpp`, `factor.cpp`, `simplex.cpp`, `g1.cpp`, `w1.cpp`, ...) that calls the real HiGHS function and
   writes JSONL (`tools/highs-oracle/FORMAT.md`). Generate with `tools/highs-oracle/gen_fixtures.sh` (setup: `tools/highs-oracle/README.md`);
   commit small fixtures under `crates/rusolve-ekk/tests/fixtures/` (they are listed in `MANIFEST.tsv`).
3. Write the test before the body. Use `tests/common/` (fixture loading, `Checker`, bit comparison) and
   `tests/port/` (builders from fixture records). Compare every f64 by bit pattern, exact counts and order,
   no tolerances. Name tests `<cpp_class>_<cpp_method>__<what>` and name the HiGHS function in the context.
   Model the file on `crates/rusolve-ekk/tests/t_hvector.rs`.
4. Translate statement by statement. Header: `//! HiGHS source: <path in highs/> (lines a-b)` and the MIT
   notice line used by neighbors. Doc line on each function: `/// HiGHS: <File>::<Class>::<method> (lines a-b)`
   and `/// C++ default arguments: ...` where any. Types keep HiGHS names; fields and locals are snake_case
   (`baseLower_` becomes `base_lower_`); overloads get a suffix.
5. Never use `mul_add` or reorder float operations. Never edit the test to make it pass.
6. Regenerate the map: `python3 tools/gen_codemap.py`.

## Check

- The new test passes: `cargo test --release -p rusolve-ekk --test <name>`; a failure prints the first
  diverging record and the HiGHS function.
- After generating the large fixtures under `target/`, run with `RUSOLVE_REQUIRE_TARGET_FIXTURES=1`, so a
  missing target fixture fails instead of being skipped.
- The parity check prints `IDENTICAL 20 of 20` (`docs/skills/contrib-identity-abtest/SKILL.md`).
- Run the CI subset:
  `cargo fmt --all --check`,
  `cargo clippy -p rusolve-ekk -p rusolve-cli -p rusolve-wasm --release --tests --benches -- -D warnings`,
  `cargo test --release -p rusolve-ekk --lib`, `cargo test --release -p rusolve-cli`.
