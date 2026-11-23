---
name: rusolve-contrib-identity-abtest
audience: contributor
description: Prove a rusolve-ekk change leaves the default search unchanged (parity check), and measure a default-off option with the A/B harness. Use when asked to verify that a refactor or a default-off option does not move the search, or to compare two settings or binaries.
---

# Parity check and A/B

Two questions, two tools. Run both arms on the same machine.

## Parity: did the default search move?

`tools/bench/identity.sh <reference cli> <new cli> [pool] [--opt k=v ...]` solves 20 fast MIPLIB instances
(`tools/bench/pools/t2fast.txt`, one thread, seed 0) with both binaries and compares status, primal and dual
bound, nodes and LP iterations. It must print `IDENTICAL 20 of 20`.

- Build the reference from `main` and the new binary from your branch, e.g. with two worktrees:
  `git worktree add ../rusolve-main main && (cd ../rusolve-main && cargo build --release -p rusolve-cli)`,
  then `cargo build --release -p rusolve-cli` in your checkout and
  `tools/bench/identity.sh ../rusolve-main/target/release/rusolve-cli ./target/release/rusolve-cli`.
- The instances are read from `$RUSOLVE_MIPLIB_DIR` (default `~/.cache/rusolve-bench/miplib-benchmark`);
  the script header shows how to download them.
- A reference that lacks your new option: pass `IDENTITY_B_OPTS="--opt <name>=false"` so only the new
  binary gets it.
- A time-limited instance is no witness; use a node limit for it.
- On a DIFF line: find the first instance, then read the first diverging record of the trace
  (`rusolve-cli mip-trace`, `crates/rusolve-cli/tests/mip_trace.rs`).

## A/B: does a default-off option help?

`tools/bench/abtest.py` (its docstring has the full usage). Its metric is instructions retired (macOS
`/usr/bin/time -l`), which barely moves with machine load, so runs can go in parallel.

1. Find the instances the option changes: `python tools/bench/abtest.py probe --pool tools/bench/pools/t2.txt --arm base=port --arm opt=port+<opt>=true --out probe.csv`, then `python tools/bench/abtest.py affected probe.csv --base base`.
2. Run: `python tools/bench/abtest.py run --pool tools/bench/pools/t2.txt --seeds 0,1 --arm base=port --arm opt=port+<opt>=true --out ab.csv`.
3. Report: `python tools/bench/abtest.py report ab.csv --base base`.

Arms are `name=spec` with spec `kind[+key=value...]`; kinds are `port`, `port115`, `highspy` and
`hbin:<path to a highs executable>`.

## Rules

- Report held-out instances separately from the ones used to tune a threshold.
- Use the smallest instance set that decides the question.
- An option that does not win stays default-off; cite the paper it implements in the option's `source`.
