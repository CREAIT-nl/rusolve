# Sandboxed Wasm hooks

Code: `crates/rusolve-wasm` (`WasmHooks`, an implementation of `MipHooks`, see `hooks.md`).
Tests: `crates/rusolve-wasm/tests/wasm_hooks.rs`, `wasm_rust_component.rs`. Example components:
`crates/rusolve-wasm/components/most_fractional.wat` (branching + cut score) and
`components/deepest_first.wat` (node selection). Overhead run:
`cargo run --release -p rusolve-wasm --example overhead`.

Loading: `RUSOLVE_MIP_HIGHS=1 rusolve mip <model> --wasm-hook PATH [--wasm-fuel N]
[--wasm-memory BYTES]` (`mip_highs.rs`; a component that fails to load — unreadable, invalid
module, an import, a bad export signature — is an input error, exit 3, naming the file) and the
MCP `solve` tool's `wasm_hook` / `wasm_fuel` / `wasm_memory` inputs (`mcp_tools.rs`, applied to that
solve only; a load failure is a tool error naming the file, and a budget given without `wasm_hook`
is a tool error saying so). With a component loaded, the `--format json`
envelope gains `stats.wasm` (the component name and the counters below); it is absent otherwise
and the schema change is additive (`rusolve.result/1` unchanged).

The core (`rusolve-ekk`) does not depend on any Wasm code: the sandbox lives in its own crate and
reaches the solver only through `Highs::set_mip_hooks(Box<dyn MipHooks>)`. With no hooks
installed nothing changes (the parity check `tools/bench/identity.sh` passes).

## Runtime choice: wasmi

`wasmi` 2.0 (MIT/Apache-2.0, audited twice per its README, a pure-Rust interpreter) over
`wasmtime`.

- **Fuel metering is built in and deterministic.** wasmi docs, `Config::consume_fuel`: "make Wasmi
  instrument its internal bytecode so that it consumes fuel as it executes. Once an execution runs
  out of fuel a `TrapCode::OutOfFuel` trap is raised. This way users can deterministically halt
  or yield the execution of WebAssembly code." Fuel is set per call with `Store::set_fuel`.
  Wasmtime has the same feature (`Config::consume_fuel`), also deterministic.
- **No extra toolchain and no JIT.** wasmi is plain Rust (a few small crates), so the crate builds
  with `cargo build` everywhere the workspace builds and is in `default-members`. wasmtime pulls
  Cranelift, a much larger dependency tree, and generates machine code at run time, which is
  unneeded for hooks that do a few hundred instructions per call (see overhead below).
- **Fast instantiation and low per-call latency** matter more here than peak throughput.
- wasmi's README: "Simple, correct and deterministic execution of WebAssembly" and "Compiler/JIT
  bomb resisting translation" (a hostile module cannot blow up compile time).

Switching to wasmtime later only touches `crates/rusolve-wasm/src/lib.rs`.

## Sandbox

- **No imports.** A module that imports anything (functions, WASI, memory, table) is rejected when
  built (`WasmHooks::new` returns an error naming the import). The linker is empty; there is no
  clock, randomness, file or network access.
- **Fuel per call** (`WasmLimits::fuel_per_call`, default 1,000,000): reset before every call, so a
  call's cost does not depend on earlier calls or on the machine. The start function runs under the
  same budget.
- **Memory cap** (`WasmLimits::max_memory_bytes`, default 1 MiB, one memory, one table, 1024 table
  elements, one instance) via a wasmi `StoreLimits`; a module declaring more fails to instantiate.
- Deep recursion ends in a stack-overflow trap, also handled as a trap.

## Failure behaviour

Per call, never per solve: a trap, `OutOfFuel`, a missing export, an invalid answer (NaN score, a
non-finite cut score, an out-of-range `select_node` index) returns `None` to the `MipHooks`
caller, which then uses the HiGHS rule for that call exactly as if no hook existed. Counters in
`WasmStats` (shared handle from `WasmHooks::stats()`): `calls`, `ok`, `declined`, `traps`,
`out_of_fuel`, `invalid`, `fuel_used`. The tests show that a component that always traps, one
that loops forever and one that answers NaN all leave the solve at the same optimum and the same
node count as the default search — for the branching-score and the node-selection hook alike.

## ABI v1

Little-endian, linear memory exported as `memory`. All exports are optional.

| export | signature | meaning |
|---|---|---|
| `branch_score` | `(depth i32, n i32, in_ptr i32, out_ptr i32) -> i32` | The host writes `n` records of four `f64` at `in_ptr` (`col`, `value`, `local_lower`, `local_upper`; 32 bytes each; `in_ptr = 0`) and expects `n` `f64` scores at `out_ptr = 32*n`. Highest score wins. Return 0 for "scores written", anything else for "no opinion". Needs an exported `memory`. |
| `cut_score` | `(violation f64, num_nonzeros i32, num_active_nonzeros i32, active_norm f64, default_score f64) -> f64` | Replacement score of a violated pool cut; NaN or infinite means "no opinion". |
| `select_node` | `(n i32, in_ptr i32) -> i32` | The host writes `n` records of four `f64` at `in_ptr = 0` (`lower_bound`, `estimate`, `depth`, `num_domain_changes`; 32 bytes each) and the component returns the index of the open node to dive into next, in `[0, n)`, or -1 for "no opinion" (the HiGHS best-bound pop is kept). An index outside `[0, n)` is invalid. If the `n` records do not fit the component's memory the call is declined without running the export. Needs an exported `memory`. |

| `separate` | `(n i32, in_ptr i32, out_ptr i32, out_cap i32) -> i32` | The user separator, with a variable-size cut result. The host writes the node LP solution, `n` `f64`, at `in_ptr = 0` and offers all remaining memory as `out_cap` `f64` slots at `out_ptr = 8*n`. The component writes cut records `nnz, rhs, (col, val)*nnz` (cut `sum val*x[col] <= rhs`, assumed globally valid), one after the other, and returns the number of `f64` slots written; 0 = no cuts, -1 = no opinion. A count above `out_cap`, a negative count other than -1, a truncated record or a non-integer `nnz`/`col` is `invalid` and drops every cut of that call (nothing partial reaches the solver). If the solution plus one record do not fit the memory the call is declined. Needs an exported `memory`. |

Like every `MipHooks` hook they run in the presolved column space. `cut_score`, `branch_score`
and `select_node` can only reorder, so a Wasm component cannot make the solver wrong. `separate`
can add rows, so its cuts go through the validity gates exactly like a native separator's
(`hooks.md`): a well-formed record with an out-of-range column, a non-finite value, a duplicate
column or a cut that removes the debug solution is rejected by the host and counted in
`HookStats::user_cuts_rejected_*`, while the ABI-level failures above are counted in `WasmStats`.
A valid-looking but wrong cut without a debug solution cannot be detected (same as native): validate
a new separator with `mip_debug_solution_file` on the validity pool before trusting it.

Example component `components/cover_cut.wat`: the knapsack cover inequality `x0+x1+x2 <= 2`,
emitted when the node LP violates it. Tests `crates/rusolve-wasm/tests/wasm_separator.rs` (model
max 10x0+11x1+12x2+3x3, 5x0+5x1+5x2+4x3 <= 14, optimum 26): the working component reaches the
optimum and its cut reaches the gate (HiGHS' own cover separator had already found the same cut, so
it is counted as a duplicate, not accepted twice); a trapping component, an infinite loop
(out of fuel), an oversized answer and a truncated record each fall back on every call (traps,
out_of_fuel or invalid equal calls, no cut offered, optimum unchanged); a record with column 77 is
rejected by the gate as malformed; a two-record answer offers two cuts per call; fuel use is the
same on two runs; a wrong `separate` signature is refused at load.

The example components are hand-written WAT; the same functions written in Rust are in the next
section.
`WasmHooks::new` accepts binary Wasm or WAT text.

## Measured overhead (release build, Apple silicon, one thread)

`cargo run --release -p rusolve-wasm --example overhead`, call cost through the `MipHooks` trait
(200,000 calls, includes memory writes/reads for `branch_score`), and end-to-end on the 28-item,
4-row knapsack of `hooks_mip.rs` (presolve and heuristics off; median of 9 runs):

| | native Rust hook | Wasm hook (wasmi) |
|---|---|---|
| `cut_score` per call | below 1 ns (inlined) | 48 ns |
| `branch_score`, 10 candidates, per call | 12 ns | 187 ns |
| average fuel per call (mix) | n/a | 189 |

| run | nodes | median time |
|---|---|---|
| no hooks | 85 | 0.0242 s |
| native hooks (same two rules) | 155 | 0.0216 s |
| Wasm hooks (same two rules) | 155 | 0.0212 s |

Reading: a Wasm hook call costs tens to a few hundred nanoseconds for small hooks, about 15x a
native call, which is negligible against an LP solve (microseconds to milliseconds). The end-to-end
difference between native and Wasm is below the run-to-run noise on this model (the two rules
change the tree: 155 nodes instead of 85, and the same tree in both, which also confirms the
Wasm component computes the native rule). The ratio will matter only for hooks called per cut
candidate on huge pools; the fuel budget bounds that cost by construction.

### User separator overhead

`cargo run --release -p rusolve-wasm --example separator_overhead` (20-thread x86-64 Linux host, one
solver thread, 200,000 calls for the per-call row, median of 201 solves for the end-to-end rows,
two runs shown; model of `wasm_separator.rs`, one root separation call):

| | native Rust separator | Wasm separator (wasmi) |
|---|---|---|
| `separate` per call, 14 columns, 1 cut returned | 19 ns | 142-145 ns |
| fuel per call | n/a | 45 |
| end to end, whole solve (no hooks: 198 / 205 us) | 189 / 195 us | 240 / 255 us |

Reading: a call costs about 7x the native call, but both are far below one LP solve. End to end the
Wasm separator is about 50 us slower per solve than the native one on this 200 us model, more than
the per-call cost explains (one call per solve here); the residual is a fixed cost of running a
hooked solve through a component (memory write, `Mutex`, parsing), measured, not yet broken down.
On this 4-column model it is a 25% slowdown; it does not grow with the node count because the
hook is called once per separation round, and the fuel budget bounds a call by construction.

## A component written in Rust

An agent that writes Rust rather than WAT compiles a `no_std` cdylib to `wasm32-unknown-unknown`.
Source: `crates/rusolve-wasm/components/rust/src/lib.rs` (standalone crate with its own
`[workspace]`, so the workspace build never touches it). It exports `branch_score` + `cut_score`
(= `most_fractional.wat`), `select_node` (= `deepest_first.wat`) and `separate` (= `cover_cut.wat`)
in one 761-byte module, `components/hooks_rust.wasm`, which is committed so the tests need no wasm
target.

Build (`tools/build_wasm_components.sh`; `--check` rebuilds and compares byte for byte):

```
rustup target add wasm32-unknown-unknown        # user-local, no sudo; toolchain pinned to 1.99.1
tools/build_wasm_components.sh                  # writes crates/rusolve-wasm/components/hooks_rust.wasm
```

Flags that matter, all set by the script (Rust: `-C link-arg`, wasm-ld docs
<https://lld.llvm.org/WebAssembly.html>): `panic = "abort"` and a `#[panic_handler]` that traps
(no unwinding, no `core::fmt` code), `lto`, `strip`, `--remap-path-prefix` (no host paths in the
file), and a memory layout `--global-base=262144 --no-stack-first -zstack-size=65536
--initial-memory=1048576 --max-memory=1048576`. The layout is the one trap of this path: the host
writes its input at address 0, where wasm-ld by default puts the stack, so the module's own data
and stack are moved above the first 256 KiB (inputs of up to 8192 32-byte records fit below it),
and the whole memory is exactly the 1 MiB cap. Reads of address 0 use `read_volatile` because a
plain null-pointer read is undefined behaviour in Rust. There are no imports, so the sandbox
accepts the module as is; `f64::floor` is not in `core`, so the example floors by hand.

Checks (`tests/wasm_rust_component.rs`): on the 28-item knapsack the Rust module reaches the
default optimum with the same tree as the WAT module (3799 nodes against 85 for the default: the
two rules are in use) and every call `ok`; with `fuel_per_call = 3` every call is `out_of_fuel`,
the tree is the default one and nothing crashes; on the cover model its `separate` offers the same
cuts and tree as `cover_cut.wat`; fuel use is identical across runs; `tools/build_wasm_components.sh
--check` rebuilds the committed file byte for byte after deleting the build directory.
CLI: `RUSOLVE_MIP_HIGHS=1 rusolve mip e2e_knap.mps --wasm-hook hooks_rust.wasm --format json`
reports `stats.wasm` = `{calls: 9, ok: 9, fuel_used: 1555, ...}`.

Per-call overhead against the WAT equivalent (Intel Xeon E5-2690 v4, release, 200,000 calls
through `MipHooks`, `cargo test --release -p rusolve-wasm --test wasm_rust_component --
--ignored --nocapture`; three runs agree within about 15%, middle run shown; fuel per call in
brackets):

| per call | WAT (ns, fuel) | Rust (ns, fuel) |
|---|---|---|
| `cut_score` | 176 (2) | 169 (2) |
| `branch_score`, 10 candidates | 532 (376) | 481 (370) |
| `select_node`, 32 open nodes | 795 (1245) | 584 (1158) |
| `separate`, 14 columns, 1 cut | 407 (45) | 414 (69) |

Reading: the compiler's output costs the same as hand-written WAT to within noise; where it
differs it is the Rust build that is faster (`select_node`: tighter loop) or a few fuel more
(`separate`: the 8-slot record is written by a loop, 69 fuel against 45 for 8 straight stores).
A Rust author therefore pays nothing at run time for not writing WAT, and the fuel budget
(`WasmLimits::fuel_per_call`, default 1,000,000) is unaffected. The absolute ns here are higher
than in the tables above because the Xeon E5-2690 v4 is an older, slower core.
