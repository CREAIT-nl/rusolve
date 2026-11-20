//! Golden test for manifest rows 40-44, `src/simplex/dual/{solve,phases,
//! iterate,update,assess}.rs` (HEkkDual.cpp), end to end through
//! `HEkk::solve` with `simplex_strategy = 1`; it also gates everything below
//! (stages 1-3).
//!
//! Fixtures: `simplex/<lp>.dual.jsonl` (adlittle, afiro, blend, kb2, sc105,
//! sc50a, share2b, afirocut, afiromax committed; the rest, including the
//! other `<lp>{cut,max}` variants, target-only, `#[ignore]` test).
//!
//! Per trace:
//! 1. Replay HApp's whole `solveLpSimplex` (`run_happ`): the scaled
//!    `HEkk::solve`, then, when HApp refines, the second, unscaled one
//!    (`solve_call` 1). Compare every `ekk_solve_enter`,
//!    `ekk_solve_initialised` and `ekk_solve_exit` record: model_status,
//!    exit_algorithm, bailout, scalars, status, final basis and every work
//!    array (bit-exact).
//! 2. Only if (1) fails: find the FIRST iteration whose observable state
//!    differs, by rerunning with `simplex_iteration_limit = k` (see
//!    `tests/port/sim.rs`):
//!    iteration count, model status, basic set (with the CHUZR/CHUZC side
//!    that differs), `h_basicIndex`, `update_count`,
//!    `updated_dual_objective_value`, `total_synthetic_tick_`, the four
//!    densities, `dual_edge_weight_[row_out]` and `h_dual_edge_weight`.
//!    Committed traces scan k = 1, 2, ...; target-only ones bisect.
//!
//! The record-by-record stream (rebuild, iter_none, phase_end, cleanup)
//! cannot be compared: the port has no hooks (signature note S1).

mod common;
mod port;

use common::*;
use port::*;

#[test]
fn dual_traces_match_highs_committed() {
    run_trace_files(
        committed_trace_files(TraceKind::Dual),
        "t_dual_trace (committed dual traces)",
        |fx| trace_gate(fx, true),
    );
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh); slow: cargo test -p rusolve-ekk --release --test t_dual_trace -- --ignored"]
fn dual_traces_match_highs_all() {
    run_trace_files(
        trace_files(TraceKind::Dual),
        "t_dual_trace (all dual traces)",
        |fx| trace_gate(fx, false),
    );
}
