//! Golden test for manifest row 39, `src/simplex/dual_row.rs`
//! (HEkkDualRow.cpp: chooseMakepack, choosePossible, chooseJoinpack,
//! chooseFinal, chooseFinalWorkGroupQuad / Heap, chooseFinalLargeAlpha,
//! updateFlip, updateDual, createFreelist, computeDevexWeight, ...).
//!
//! HEkkDualRow has no isolated fixture (the packed row, workData,
//! workGroup are not dumped), so it is gated through the dual traces,
//! iteration by iteration (iteration-limit method, `tests/port/sim.rs`):
//! after iteration k the port's basic set must gain exactly the fixture's
//! `variable_in` (CHUZC + BFRT: `chooseFinal`), and
//! `updated_dual_objective_value` (`updateFlip` / `updateDual`),
//! `row_ap_density` and `h_basicIndex` must match.
//!
//! Like `t_dual_rhs`, this target scans EVERY iteration of every committed
//! dual trace; a divergence attributed to another component still fails
//! here, because HEkkDualRow is unverified after it.

mod common;
mod port;

use common::*;
use port::*;

#[test]
fn dual_row_iterations_match_highs_committed() {
    run_trace_files(
        committed_trace_files(TraceKind::Dual),
        "t_dual_row (committed dual traces, every iteration)",
        |fx| iteration_gate(fx, true, "HEkkDualRow (CHUZC + BFRT)"),
    );
}

#[test]
#[ignore = "target-only fixtures; slow: cargo test -p rusolve-ekk --release --test t_dual_row -- --ignored"]
fn dual_row_iterations_match_highs_all() {
    run_trace_files(
        trace_files(TraceKind::Dual),
        "t_dual_row (all dual traces, bisection)",
        |fx| iteration_gate(fx, false, "HEkkDualRow (CHUZC + BFRT)"),
    );
}
