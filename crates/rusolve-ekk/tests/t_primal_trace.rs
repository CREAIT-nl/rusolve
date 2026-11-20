//! Golden test for manifest rows 45-50, `src/simplex/primal/*.rs`
//! (HEkkPrimal.cpp), end to end through `HEkk::solve` with
//! `simplex_strategy = 4`.
//!
//! Fixtures: `simplex/<lp>.primal.jsonl` (adlittle, afiro, blend, kb2,
//! sc105, sc50a, share2b, afirocut, afiromax committed; the rest target-only,
//! `#[ignore]`).
//!
//! Per trace, as `t_dual_trace`: the enter/initialised/exit records of every
//! `HEkk::solve` HApp runs (the `<lp>cut` primal traces re-solve unscaled
//! with dual); on a mismatch, the first iteration whose observable state differs
//! (iteration-limit method): iteration count, model status, basic set
//! (CHUZC `chuzc`/`chooseColumn` vs CHUZR `chooseRow`), `h_basicIndex`,
//! `update_count`, `updated_primal_objective_value`, `total_synthetic_tick_`.
//! Not observable: `HEkkPrimal::edge_weight_` (`h_edge_weight`,
//! `edge_weight_in`), which lives in the local solver object.

mod common;
mod port;

use common::*;
use port::*;

#[test]
fn primal_traces_match_highs_committed() {
    run_trace_files(
        committed_trace_files(TraceKind::Primal),
        "t_primal_trace (committed primal traces)",
        |fx| trace_gate(fx, true),
    );
}

#[test]
#[ignore = "target-only fixtures; slow: cargo test -p rusolve-ekk --release --test t_primal_trace -- --ignored"]
fn primal_traces_match_highs_all() {
    run_trace_files(
        trace_files(TraceKind::Primal),
        "t_primal_trace (all primal traces)",
        |fx| trace_gate(fx, false),
    );
}
