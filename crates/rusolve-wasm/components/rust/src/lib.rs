//! Wasm hook component in Rust (ABI v1, `crates/rusolve-wasm/src/lib.rs`, docs/design/wasm-hooks.md).
//! Same four exports as the hand-written WAT examples, so the two can be compared call for call:
//! * `branch_score` + `cut_score` = `components/most_fractional.wat` (branch on the most fractional
//!   candidate; cut score = raw violation);
//! * `select_node` = `components/deepest_first.wat` (deepest node, ties by the lowest bound);
//! * `separate` = `components/cover_cut.wat` (cover cut x0 + x1 + x2 <= 2).
//!
//! `no_std`, no imports (the sandbox rejects any), no allocation. Memory layout: the host writes its
//! input at address 0, so `build_wasm_components.sh` links with `--global-base=262144` and
//! `--no-stack-first`, which moves this module's own data and stack above the first 256 KiB; the
//! host-written region (inputs of up to 8192 records, then the `separate` output) cannot clobber
//! them. Reads from address 0 are volatile because a plain null-pointer read is undefined in Rust.
#![no_std]

use core::ptr::{read_volatile, write_volatile};

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    core::arch::wasm32::unreachable()
}

#[inline(always)]
fn rd(base: i32, idx: usize) -> f64 {
    // SAFETY: `base + 8 * idx` is a host-written slot inside the module's linear memory.
    unsafe { read_volatile((base as usize + 8 * idx) as *const f64) }
}

/// `f64::floor` is not in `core` yet; candidates are integer-variable values, far below 2^63.
#[inline(always)]
fn floor(x: f64) -> f64 {
    let t = x as i64 as f64;
    if t > x {
        t - 1.0
    } else {
        t
    }
}

#[inline(always)]
fn wr(base: i32, idx: usize, v: f64) {
    // SAFETY: `base + 8 * idx` is inside the memory region the host offered for output.
    unsafe { write_volatile((base as usize + 8 * idx) as *mut f64, v) }
}

/// Score of each fractional candidate: distance to the nearest integer, 0.5 - |frac - 0.5|.
/// Records at `inp`: col, value, local_lower, local_upper (4 f64); scores go to `out`.
#[no_mangle]
pub extern "C" fn branch_score(_depth: i32, n: i32, inp: i32, out: i32) -> i32 {
    for i in 0..n.max(0) as usize {
        let v = rd(inp, 4 * i + 1);
        let frac = v - floor(v);
        wr(out, i, 0.5 - (frac - 0.5).abs());
    }
    0
}

/// Pool cut score: the raw violation.
#[no_mangle]
pub extern "C" fn cut_score(violation: f64, _nnz: i32, _active: i32, _norm: f64, _default: f64) -> f64 {
    violation
}

/// Deepest open node first, ties by the lowest lower bound. Records: lower_bound, estimate,
/// depth, num_domain_changes (4 f64).
#[no_mangle]
pub extern "C" fn select_node(n: i32, inp: i32) -> i32 {
    let mut best = -1;
    let mut best_depth = f64::NEG_INFINITY;
    let mut best_lb = f64::INFINITY;
    for i in 0..n.max(0) as usize {
        let lb = rd(inp, 4 * i);
        let depth = rd(inp, 4 * i + 2);
        if depth > best_depth || (depth == best_depth && lb < best_lb) {
            best = i as i32;
            best_depth = depth;
            best_lb = lb;
        }
    }
    best
}

/// Cover cut x0 + x1 + x2 <= 2 when the node LP solution violates it; one record
/// `nnz, rhs, (col, val)*nnz` = 3, 2, 0,1, 1,1, 2,1 (8 f64). -1 = no opinion, 0 = no cut.
#[no_mangle]
pub extern "C" fn separate(n: i32, inp: i32, out: i32, cap: i32) -> i32 {
    if n < 3 {
        return -1;
    }
    if cap < 8 {
        return 0;
    }
    if rd(inp, 0) + rd(inp, 1) + rd(inp, 2) <= 2.0001 {
        return 0;
    }
    let rec = [3.0, 2.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0];
    for (i, v) in rec.iter().enumerate() {
        wr(out, i, *v);
    }
    8
}
