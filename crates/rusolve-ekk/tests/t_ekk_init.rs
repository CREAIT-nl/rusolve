//! Golden test for manifest rows 29-37, `src/simplex/ekk/*.rs` on the
//! HEkk initialisation path (HEkk.cpp: moveLp, initialiseEkk,
//! setSimplexOptions, initialiseSimplexLpRandomVectors (both draws),
//! initialiseForSolve, initialiseSimplexLpBasisAndFactor, setBasis,
//! computeFactor, initialisePartitionedRowwiseMatrix,
//! allocateWorkAndBaseArrays, initialiseCost, initialiseBound,
//! initialiseNonbasicValueAndMove, computePrimal, computeDual,
//! computeSimplexInfeasible, compute{Dual,Primal}ObjectiveValue;
//! HEkkControl.cpp: initialiseControl).
//!
//! Fixtures: every `simplex/<lp>.<dual|primal>.jsonl` (7 LPs x 2 committed;
//! the rest target-only, in the `#[ignore]` test). For each trace:
//!
//! 1. Input: `ekk_solve_enter` (scaled LP + every HighsOptions value).
//! 2. HApp::solveLpSimplex prefix (`iteration_count_ = 0;
//!    initialiseSimplexStats(); moveLp(solver_object)`), then the
//!    statements of `HEkk::solve` before the enter hook (`debugInitialise`,
//!    `initialiseAnalysis`, `initialiseControl`, cleanup levels, cycling).
//! 3. Compare `ekk_solve_enter` (status, cost_scale, LP after moveLp).
//! 4. `initialise_for_solve()`, then compare the first `invert` record and
//!    `ekk_solve_initialised` (random vectors, logical basis, work arrays,
//!    perturbation-free costs and bounds, `ar_matrix_`, scalars, strategy
//!    echo).
//!
//! Records are compared with the trace reporter, so a failure prints the
//! first diverging record with every mismatching field (bits + decimals)
//! and the HiGHS function likely responsible.

mod common;
mod port;

use common::*;
use port::*;

fn check_trace(fx: &Fixture) -> Result<(), String> {
    let enter = enter_of(fx);
    let ctx = fx.file.clone();
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut out = Vec::new();
        at("HApp prefix: initialiseSimplexStats + HEkk::moveLp".into());
        let (mut ekk, _parts) = ekk_after_move_lp(&enter, &ctx, None);
        at("HEkk::solve prefix (debugInitialise, initialiseAnalysis, initialiseControl)".into());
        ekk_solve_prefix(&mut ekk);
        out.push(enter_obj(&ekk, &enter));
        at("HEkk::initialiseForSolve".into());
        ekk.initialise_for_solve();
        out.push(invert_obj(&ekk));
        out.push(initialised_obj(&ekk));
        out
    }));
    let actual = r.map_err(|p| {
        let msg = p
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
            .unwrap_or_else(|| "<non-string panic>".into());
        format!("{}: PANIC inside the port: {msg}", fx.file)
    })?;
    compare_traces(
        fx,
        &actual,
        // HEkk::build_synthetic_tick_ is set by resetSyntheticClock()
        // (HEkk.cpp:1533) AFTER the computeFactor hook, so the state after
        // initialise_for_solve cannot show its value at the hook (0). The
        // reset value is checked in ekk_solve_initialised.scalars.
        &TraceCompare::default()
            .skip("invert/build_synthetic_tick")
            .only(&["ekk_solve_enter", "invert", "ekk_solve_initialised"])
            .max_records(3),
    )
    .map_err(|d| d.to_string())
}

fn run(files: Vec<String>, what: &str) {
    let mut fails = Failures::default();
    for file in files {
        let Some(fx) = load(&file) else { continue };
        match check_trace(&fx) {
            Ok(()) => fails.cases += 1,
            Err(m) => fails.fail(m),
        }
    }
    fails.finish(what);
}

#[test]
fn ekk_initialisation_matches_highs_committed() {
    let mut files = committed_trace_files(TraceKind::Dual);
    files.extend(committed_trace_files(TraceKind::Primal));
    run(files, "t_ekk_init (committed traces)");
}

#[test]
#[ignore = "target-only fixtures (tools/highs-oracle/gen_fixtures.sh): cargo test -p rusolve-ekk --test t_ekk_init -- --ignored"]
fn ekk_initialisation_matches_highs_all() {
    let mut files = trace_files(TraceKind::Dual);
    files.extend(trace_files(TraceKind::Primal));
    run(files, "t_ekk_init (all traces)");
}
