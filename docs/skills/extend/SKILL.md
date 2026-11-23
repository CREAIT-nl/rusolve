---
name: rusolve-extend
description: Add a problem-specific cut separator, branching rule or node selector to rusolve, and validate it. Use when asked to extend the solver for a problem class, add a cut family, or write a custom heuristic or branching hook.
---

# Extend rusolve for a problem class

Extension points are the `MipHooks` trait (`crates/rusolve-ekk/src/mip/hooks.rs`, design and
call sites in `docs/design/hooks.md`) and the ported HiGHS callbacks
(`crates/rusolve-ekk/src/lp_data/highs_callback.rs`). Hooks are empty by default and the default
search is unchanged when none is installed.

## Hooks

| hook | you return | can be wrong? |
|---|---|---|
| `branching_score` | one score per fractional candidate, highest wins | only slower |
| `select_node` | which open node to process next | only slower |
| `cut_score` | the efficacy score of a pool cut | only slower |
| `separate` | cuts `a x <= b` (as `UserCut` rows) in the presolved space | yes: an invalid cut loses solutions |

Hooks run in the presolved column space; use `--opt presolve=off` to get original indices while
developing. Install with `Highs::set_mip_hooks(Box<dyn MipHooks>)`, read counters with
`mip_hook_stats()`. The three score hooks are also reachable without compiling any Rust, as a
sandboxed Wasm component — see "No Rust needed: the Wasm path" below.

## Adding a cut family: workflow

1. Derive the inequality from a source you can cite (a paper or textbook; AGENTS.md requires it
   in a code comment). Try it first through the solver you already have, without writing Rust:
   with MCP, `mcp__rusolve__modify` `add_row` adds the cut as a constraint to a loaded model, and
   `mcp__rusolve__solve` shows whether the root bound moves.
2. Check validity before anything else: a cut is valid only if no integer-feasible point of the
   original model violates it. Test it by optimising the cut's left-hand side over the original
   model (put its coefficients in the objective with `modify` `cost` edits, solve, compare with
   the cut's right-hand side), and, for a real separator, by running with a known optimal
   solution: `--opt mip_debug_solution_file=<sol file>` makes the solver reject and report
   (`HOOK VALIDITY GATE`) any user cut that cuts off that solution
   (`docs/design/hooks.md`, "Validity gates").
3. Implement `separate` in a new type implementing `MipHooks`, return only cuts violated by the
   current LP solution (`SeparationInfo`), keep rows sparse, and cap the number per round.
4. Measure on the family it was built for: root dual bound after cuts, nodes, time, against the
   same run without the hook (`stats.nodes`, `dual_bound`, `stats.time_s` in the `rusolve.result/1` envelope).
   Keep it only if it wins.
5. Run the repository gates (`cargo fmt`, `cargo clippy -- -D warnings`, the test suites) and
   `tools/bench/identity.sh`: it must still report identical results with no hook installed.

## No Rust needed: the Wasm path

For a branching rule, a node selector or a cut-ranking rule you do not touch the repository at
all: write it as a sandboxed WebAssembly component and hand it to the solver. Design and ABI:
`docs/design/wasm-hooks.md`; examples: `crates/rusolve-wasm/components/most_fractional.wat`
(branch and cut score) and `components/deepest_first.wat` (node selection).

- **Write it** as `.wat` (WebAssembly Text) or `.wasm`; both are accepted. No imports of any kind
  (a module that imports anything is rejected before the solve), one linear memory exported as
  `memory` and capped at 1 MiB, and a fuel budget per call.
- **ABI v1 exports**, all optional: `branch_score` (one score per fractional candidate, highest
  wins), `cut_score` (the score of a violated pool cut), `select_node` (which open node to explore
  next). "No opinion" is `-1` for `select_node`, a non-zero return for `branch_score`, NaN for
  `cut_score`; that one call then uses the HiGHS rule.
- **Load it**: `rusolve mip <model> --wasm-hook rule.wat --format json`, or
  the MCP `solve` tool with `wasm_hook` (plus optional `wasm_fuel` and `wasm_memory`, the byte
  cap on the component's memory). A component that will not load is exit 3 naming the file, or a
  tool error naming the file.
- **Read `stats.wasm`** in the `rusolve.result/1` envelope: `{component, calls, ok, declined,
  traps, out_of_fuel, invalid, fuel_used}`. `calls > 0` proves the solver reached your rule;
  `ok` counts the calls where it actually took the decision, `declined` the ones where it said no
  opinion, and `traps` / `out_of_fuel` / `invalid` the ones silently answered by HiGHS. A rule
  that never made a decision has `ok: 0` and changes nothing. `stats.wasm` is absent with no hook.
- **The three ordering hooks can only be slower, never wrong**: they just reorder, so the objective must
  not move. If it does, that is a solver bug, not something to design around.
- `separate` (user cuts) is in the ABI: export `separate(n, in_ptr, out_ptr, out_cap) -> i32`. It
  writes records `nnz, rhs, (col, val)*nnz` at `out_ptr` and returns the number of records; -1
  means no opinion, 0 means no cuts. A trap, out of fuel, or a malformed or oversized answer falls
  back to no cuts. Unlike the other hooks it can be wrong: the cuts go through the validity gate, so
  validate them with `mip_debug_solution_file` (a known optimum must stay feasible). Example:
  `crates/rusolve-wasm/components/cover_cut.wat`.

Minimal component: dive into the deepest open node.

```wasm
(module
  (memory (export "memory") 1)
  ;; select_node(n, in_ptr) -> index of the node to explore next in [0, n), or -1 for "no opinion".
  ;; The host wrote n records of 32 bytes at in_ptr: lower_bound +0, estimate +8, depth +16,
  ;; num_domain_changes +24 (all f64, little-endian).
  (func (export "select_node") (param $n i32) (param $in i32) (result i32)
    (local $i i32) (local $best i32) (local $best-depth f64) (local $d f64)
    (local.set $best (i32.const -1))
    (local.set $best-depth (f64.const -inf))
    (block $done
      (loop $next
        (br_if $done (i32.ge_u (local.get $i) (local.get $n)))
        (local.set $d
          (f64.load offset=16 (i32.add (local.get $in) (i32.shl (local.get $i) (i32.const 5)))))
        (if (f64.gt (local.get $d) (local.get $best-depth))
          (then
            (local.set $best (local.get $i))
            (local.set $best-depth (local.get $d))))
        (local.set $i (i32.add (local.get $i) (i32.const 1)))
        (br $next)))
    (local.get $best)))
```

### The same component in Rust

If you write Rust rather than WAT, compile a `no_std` cdylib to `wasm32-unknown-unknown`
(`rustup target add wasm32-unknown-unknown`, user-local). Full example with all four exports:
`crates/rusolve-wasm/components/rust/src/lib.rs`; build with `tools/build_wasm_components.sh`
(the built `components/hooks_rust.wasm` is committed). The compiled module costs the same per call
as the WAT (`docs/design/wasm-hooks.md`, "A component written in Rust"). Three rules:

- no imports: `#![no_std]`, `panic = "abort"` and a `#[panic_handler]` that traps; no allocation;
- the host writes its input at address 0, so move the module's own data and stack above it:
  `-C link-arg=--global-base=262144 -C link-arg=--no-stack-first -C link-arg=-zstack-size=65536
  -C link-arg=--initial-memory=1048576 -C link-arg=--max-memory=1048576` (the memory cap is 1 MiB),
  and read the input with `read_volatile` (a plain read of address 0 is undefined in Rust);
- `f64::floor` is not in `core`: floor by hand.

```rust
#![no_std]
use core::ptr::read_volatile;

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! { core::arch::wasm32::unreachable() }

/// select_node(n, in_ptr): the host wrote n records of 4 f64 (lower_bound, estimate, depth,
/// num_domain_changes) at in_ptr; return the deepest node, ties by the lowest bound.
#[no_mangle]
pub extern "C" fn select_node(n: i32, inp: i32) -> i32 {
    let (mut best, mut best_depth, mut best_lb) = (-1, f64::NEG_INFINITY, f64::INFINITY);
    for i in 0..n.max(0) as usize {
        let rec = inp as usize + 32 * i;
        // SAFETY: the record is host-written memory inside the component's linear memory.
        let (lb, depth) = unsafe {
            (read_volatile(rec as *const f64), read_volatile((rec + 16) as *const f64))
        };
        if depth > best_depth || (depth == best_depth && lb < best_lb) {
            (best, best_depth, best_lb) = (i as i32, depth, lb);
        }
    }
    best
}
```

The next block checks the committed Rust-built module (no imports, exports `memory` and the four
hooks) and, when a built `rusolve-cli` is at hand, that the solver really calls it
(`stats.wasm.ok > 0`):

```python
import json, os, subprocess

def find(rel):
    d = os.getcwd()
    while True:
        if os.path.exists(os.path.join(d, rel)):
            return os.path.join(d, rel)
        if os.path.dirname(d) == d:
            raise SystemExit("run from inside the rusolve repository: " + rel)
        d = os.path.dirname(d)

wasm = find("crates/rusolve-wasm/components/hooks_rust.wasm")
data = open(wasm, "rb").read()
assert data[:8] == b"\0asm\x01\0\0\0"

def leb(b, i):
    r = s = 0
    while True:
        c = b[i]; i += 1
        r |= (c & 0x7F) << s; s += 7
        if not c & 0x80:
            return r, i

sections, i = {}, 8
while i < len(data):
    sid = data[i]; size, j = leb(data, i + 1)
    sections[sid] = data[j:j + size]; i = j + size
assert 2 not in sections, "an import section means the sandbox rejects the module"
body = sections[7]; count, k = leb(body, 0); exports = []
for _ in range(count):
    n, k = leb(body, k); exports.append(body[k:k + n].decode()); k += n + 1; _, k = leb(body, k)
assert {"memory", "branch_score", "cut_score", "select_node", "separate"} <= set(exports), exports

cli = os.path.join(os.path.dirname(find("Cargo.toml")), "target", "release", "rusolve-cli")
if os.path.exists(cli):
    model = find("crates/rusolve-ekk/tests/fixtures/s5/e2e_knap.mps")
    out = subprocess.run([cli, "mip", model, "--wasm-hook", wasm, "--format", "json"],
                         capture_output=True, text=True, env={**os.environ, "RUSOLVE_MIP_HIGHS": "1"})
    st = json.loads(out.stdout)
    assert abs(st["objective"] + 417.0) < 1e-6 and st["stats"]["wasm"]["ok"] > 0, st
print("rust component ok:", sorted(exports))
```

## Boundaries

- Base a new technique on published literature and cite the paper in the code.
- A hook that can cut off solutions needs the validity gate above before it is trusted.
- The hook trait is Rust-only in-process. The Wasm sandbox above is how you extend the solver
  without building the repository; `crates/rusolve-py` exposes neither yet.
