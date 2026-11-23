//! Golden test for manifest row 38, `src/simplex/dual_rhs.rs`
//! (HEkkDualRHS.cpp: setup, chooseNormal, updatePrimal, updatePivots,
//! updateInfeasList, createArrayOfPrimalInfeasibilities, createInfeasList,
//! assessOptimality), with the DSE weight path it feeds
//! (HEkkDual::updateFtranDSE, HEkk::updateDualSteepestEdgeWeights).
//!
//! HEkkDualRHS has no isolated fixture (`work_infeasibility`,
//! `workCount`, `workCutoff` and the CHUZR random draw are not dumped), so
//! it is gated through the dual traces, iteration by iteration
//! (iteration-limit method, `tests/port/sim.rs`): after iteration k the
//! port's basic set must lose exactly the fixture's `variable_out`
//! (CHUZR: `chooseNormal`), and `dual_edge_weight_[row_out]`,
//! `h_dual_edge_weight`, `row_DSE_density` and `h_basicIndex` must match.
//!
//! This target scans EVERY iteration of every committed dual trace, even
//! when the exit record matches. A failure names the first diverging
//! iteration and which side (CHUZR vs CHUZC vs weights vs control flow) the
//! mismatching fields point at; a divergence attributed to another
//! component still fails here, because HEkkDualRHS is unverified after it.

mod common;
mod port;

use common::*;
use port::*;

#[test]
fn dual_rhs_iterations_match_highs_committed() {
    run_trace_files(
        committed_trace_files(TraceKind::Dual),
        "t_dual_rhs (committed dual traces, every iteration)",
        |fx| iteration_gate(fx, true, "HEkkDualRHS (CHUZR + DSE)"),
    );
}

#[test]
#[ignore = "target-only fixtures; slow: cargo test -p rusolve-ekk --release --test t_dual_rhs -- --ignored"]
fn dual_rhs_iterations_match_highs_all() {
    run_trace_files(
        trace_files(TraceKind::Dual),
        "t_dual_rhs (all dual traces, bisection)",
        |fx| iteration_gate(fx, false, "HEkkDualRHS (CHUZR + DSE)"),
    );
}
