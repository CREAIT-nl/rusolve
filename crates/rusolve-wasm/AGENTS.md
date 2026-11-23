# rusolve-wasm

Agent-written MIP hooks as fuel-metered WebAssembly with no host access (ABI in `docs/design/wasm-hooks.md`).

## What lives here
- `src/lib.rs`: `WasmHooks`, which implements `rusolve_ekk::mip::hooks::MipHooks` by calling exports of a module run by the `wasmi` interpreter, plus `WasmStats`.
- `components/`: example Wasm hook components; `examples/overhead.rs`: per-call overhead measurement.
- `tests/wasm_hooks.rs`: the contract tests.

## Invariants
- No imports for the module (no WASI, no clock, no randomness), a hard memory cap and a fixed fuel budget per call, so a call costs the same on every run.
- A trap, fuel exhaustion, a missing export or an invalid answer makes that single call fall back to the HiGHS behavior and is counted in `WasmStats`; the solve is never interrupted.
- Hooks are empty by default and the default search does not change when none is installed, so identity stays 20/20.
- Pinned dependency `wasmi`; keep the version in `Cargo.toml`.

## How to test
- `cargo test --release -p rusolve-wasm`.
- CLI side: `crates/rusolve-cli/tests/wasm_cli.rs`.
- Hook trait and call sites: `crates/rusolve-ekk/src/mip/hooks.rs`, tests in `crates/rusolve-ekk/tests/hooks_mip.rs`.

## Change recipes
1. Add a hook to the ABI: extend the `MipHooks` trait in `crates/rusolve-ekk/src/mip/hooks.rs` with a default body, add the export call in `src/lib.rs`, a component under `components/`, a test in `tests/wasm_hooks.rs`, and a section in `docs/design/wasm-hooks.md`.
2. Change a resource limit: edit the constants in `src/lib.rs`, show the overhead with `examples/overhead.rs`.
3. Count a new failure kind: add a counter to `WasmStats` and assert it in a test that provokes it.

## Pitfalls
- Hooks run in the presolved column space; develop with `--opt presolve=off` to see original indices.
- An invalid cut from a `separate` hook loses solutions; the score hooks only cost speed (`docs/skills/extend/SKILL.md`).
- Further reading: `docs/design/hooks.md`, `docs/design/wasm-hooks.md`.
- Run `cargo fmt --all` before committing; CI checks it.
