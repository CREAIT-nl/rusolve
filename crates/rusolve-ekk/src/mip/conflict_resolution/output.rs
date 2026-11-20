//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! conflict_resolution.c::addConflictRows (2908-2958), SCIPconflictAddConflictCon (1958-2098),
//! scip_conflict.c (735-764), HighsCutGeneration.cpp (1262-1291), HighsDebugSol.cpp (27-30, 70-74, 125, 167-177).

use super::*;
use crate::mip::conflict_resolution::reason::S2DomainCtx;
use crate::mip::conflict_resolution::reason::s2_initial_conflict_row;
use crate::mip::conflict_resolution::resolve::analyze_resolution;
use crate::mip::mip_solver::HighsMipSolver;
use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;

/// SCIP 11.0.0 (32dd425): scip_conflict.c::SCIPanalyzeConflictCons (lines 735-764)
/// and conflict_resolution.c::addConflictRows (lines 2908-2958)
pub fn s2_conflict_analysis(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    globaldom: usize,
    conflict_pool: usize,
    proof: Option<(&[HighsInt], &[f64], HighsInt, f64)>,
) {
    // a. Return if proof.is_some() && !S2_RESOLVE_LP_PROOFS
    if proof.is_some() && !S2_RESOLVE_LP_PROOFS {
        return;
    }

    // b. localdom.replay_pool_events(ms)
    localdom.replay_pool_events(ms);

    // c. Take workspace, ensure_dim, stats.calls += 1
    let mut ws = ms.mipdata_.conflict_resolution.take(conflict_pool);
    // Rows are indexed by the domain's columns (model_, presolved when presolve is on), as S2DomainCtx::num_col.
    let num_col = localdom.col_lower_.len();
    ws.ensure_dim(num_col);
    ws.stats.calls += 1;

    // d. Env flags: stats_env, and the debug solution taken from HighsDebugSol (inert unless
    //    the mip_debug_solution_file option names a file; see mip::debug_sol). Refreshed every
    //    call: HighsDebugSol deactivates itself once an incumbent beats the debug objective,
    //    and maps the solution into the presolved space.
    if ws.stats_env.is_none() {
        ws.stats_env = Some(std::env::var("RUSOLVE_S2_STATS").is_ok_and(|v| v == "1"));
    }
    s2_trace_dump_enable(std::env::var("RUSOLVE_S2_TRACE_DUMP").is_ok_and(|v| v == "1"));
    {
        let ds = &ms.mipdata_.debug_solution;
        if ds.debug_sol_active {
            if ws.debug_sol.x.as_deref() != Some(ds.debug_solution.as_slice()) {
                ws.debug_sol.x = Some(ds.debug_solution.clone());
                ws.debug_sol.announced = false;
            }
        } else {
            ws.debug_sol.x = None;
        }
    }

    // e. Build initial conflict row
    let mut initial = std::mem::take(&mut ws.conflictrow);
    initial.ensure_dim(num_col);
    let upto = {
        let mut ctx = S2DomainCtx {
            ms: &*ms,
            localdom,
            globaldom: &ms.mipdata_.domains[globaldom],
        };
        s2_initial_conflict_row(&mut ctx, proof, &mut initial)
    };

    if upto.is_none() {
        ws.stats.initial_fail += 1;
        ws.conflictrow = initial;
        ms.mipdata_.conflict_resolution.put(conflict_pool, ws);
        return;
    }

    // f. Analysis in its own block with a fresh S2DomainCtx
    {
        let mut ctx = S2DomainCtx {
            ms: &*ms,
            localdom,
            globaldom: &ms.mipdata_.domains[globaldom],
        };
        analyze_resolution(&mut ctx, &mut ws, &initial, upto.unwrap());
    }
    ws.conflictrow = initial;

    // g. add_conflict_rows (S2DomainCtx must already be dropped)
    let rows_added_before = ws.stats.rows_added;
    add_conflict_rows(ms, &mut ws, globaldom, conflict_pool);

    // h. Stats line if env on and rows_added changed
    if ws.stats_env == Some(true) && ws.stats.rows_added != rows_added_before {
        eprintln!(
            "s2-stats calls={} rows_found={} rows_added={} long={} steps={}",
            ws.stats.calls,
            ws.stats.rows_found,
            ws.stats.rows_added,
            ws.stats.nreslongconfs,
            ws.stats.nressteps
        );
    }

    // i. Restore workspace
    ms.mipdata_.conflict_resolution.put(conflict_pool, ws);
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::addConflictRows (lines 2908-2958)
/// and SCIPconflictAddConflictCon (lines 1958-2098).
fn add_conflict_rows(
    ms: &mut HighsMipSolver,
    ws: &mut HighsConflictResolution,
    globaldom: usize,
    conflict_pool: usize,
) {
    let nint = ms.mipdata_.integral_cols.len() as HighsInt;
    let maxsize = 2 * cpp_max((S2_MAXVARSFRACRES * nint as f64) as HighsInt, S2_MINMAXVARS);

    let nconstoadd = std::cmp::min(S2_RESFUIPLEVELS as usize, ws.conflictrows.len());

    for i in 0..nconstoadd {
        let row = &ws.conflictrows[i];

        // SCIP addConflictRows line 2936: IsLT(coefquotient, maxcoefquot)
        if !s2_is_lt(row.coefquotient, S2_MAXCOEFQUOT) {
            continue;
        }

        let nnz = row.nnz();

        // SCIPconflictAddConflictCon line 2006: nnz > maxsize
        if nnz > maxsize as usize {
            ws.stats.nreslongconfs += 1;
            continue;
        }

        // SCIPconflictAddConflictCon line 2014: nnz == 0 && IsGT(lhs, 0.0)
        if nnz == 0 {
            if s2_is_gt(row.lhs, 0.0) {
                ws.stats.empty_infeasible += 1;
            }
            continue;
        }

        // Diagnostic: the row as resolution produced it, before the <= conversion and before
        // tighten_coefficients. Reports the first invalid row even when tightening is innocent.
        s2_trace_row(&format!("final[{i}] pre-tighten"), row, &ws.debug_sol);

        // Otherwise: convert to <=, tighten, check integral, debug check, add cut
        let (mut inds, mut vals, mut rhs) = to_leq(row);
        let len = inds.len() as HighsInt;

        ms.mipdata_.domains[globaldom]
            .tighten_coefficients(&*ms, &mut inds, &mut vals, len, &mut rhs);

        let integral = cut_is_integral(ms, &inds, &vals);

        debug_sol_check(ms, ws, &inds, &vals, rhs);

        let cutindex = detach!(ms.mipdata_.cutpools[conflict_pool], |p| {
            p.add_cut(
                ms, &mut inds, &mut vals, len, rhs, integral, true, true, true,
            )
        });

        if cutindex != -1 {
            ws.stats.rows_added += 1;
        }
    }
}

/// Convert a SCIP conflict row (>= form) to a HiGHS cut (<= form).
/// inds copied, vals = -v, rhs = -lhs.
fn to_leq(row: &S2ConflictRow) -> (Vec<HighsInt>, Vec<f64>, f64) {
    let inds: Vec<HighsInt> = row.inds.clone();
    let vals: Vec<f64> = row.inds.iter().map(|&j| -row.vals[j as usize]).collect();
    let rhs = -row.lhs;
    (inds, vals, rhs)
}

/// All columns integral and every |v - v.round()| <= epsilon.
fn cut_is_integral(ms: &HighsMipSolver, inds: &[HighsInt], vals: &[f64]) -> bool {
    cut_is_integral_by(
        |col: HighsInt| ms.is_col_integral(col),
        ms.mipdata_.epsilon,
        inds,
        vals,
    )
}

/// All columns integral (per is_int) and every |v - v.round()| <= epsilon.
fn cut_is_integral_by(
    is_int: impl Fn(HighsInt) -> bool,
    epsilon: f64,
    inds: &[HighsInt],
    vals: &[f64],
) -> bool {
    for k in 0..inds.len() {
        let col = inds[k];
        if !is_int(col) {
            return false;
        }
        let v = vals[k];
        if (v - v.round()).abs() > epsilon {
            return false;
        }
    }
    true
}

/// Debug-solution check on the row about to enter the conflict pool, in the same style as
/// `HighsCutPool::add_cut`'s `check_cut` (HighsDebugSol.cpp:167-177) and behind the same
/// inert-unless-active switch: `mip_debug_solution_file` (`mip::debug_sol`). Unlike `check_cut`
/// it names the row index and reports the S2-side violation.
fn debug_sol_check(
    ms: &HighsMipSolver,
    ws: &mut HighsConflictResolution,
    inds: &[HighsInt],
    vals: &[f64],
    rhs: f64,
) {
    if ms.submip {
        return;
    }
    let ds = &ms.mipdata_.debug_solution;
    if !ds.debug_sol_active {
        return;
    }
    let x = match &ws.debug_sol.x {
        Some(x) => x,
        None => return,
    };
    // rows are indexed by the local domain's columns, which is the space HighsDebugSol
    // maps the solution into (debug_sol.rs run_setup_block).
    if x.len() != ds.debug_solution.len() {
        return;
    }

    if !ws.debug_sol.announced {
        ws.debug_sol.announced = true;
        eprintln!("s2-debugsol active obj={}", ds.debug_sol_objective);
    }

    ws.stats.debugsol_checked += 1;

    // viol = -rhs + sum vals[k] * x[inds[k]] in HighsCDouble
    let viol = s2_debugsol_viol(x, inds, vals, rhs);

    if viol > ms.mipdata_.feastol {
        ws.stats.debugsol_violations += 1;
        eprintln!(
            "s2-debugsol VIOLATION viol={viol} nnz={} rhs={rhs}",
            inds.len()
        );
    }
}

/// Pure violation computation: viol = -rhs + sum vals[k] * x[inds[k]] in HighsCDouble.
/// HighsDebugSol.cpp:167-177.
pub(crate) fn s2_debugsol_viol(x: &[f64], inds: &[HighsInt], vals: &[f64], rhs: f64) -> f64 {
    let mut violation = HighsCDouble::from(-rhs);
    for k in 0..inds.len() {
        // C++ `debugSolution[Rindex[i]] * Rvalue[i]` is a plain double product.
        violation += x[inds[k] as usize] * vals[k];
    }
    f64::from(violation)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;
    use crate::prelude::*;

    const N_ROWS: usize = 3;

    /// A strongly correlated multi-row binary knapsack (value = mean weight + 10), so the
    /// root LP bound does not close and the tree branches: min -value.x,
    /// sum_j w_ij x_j <= half the row weight, x binary. Weights come from a fixed LCG.
    fn knapsack_lp(n: usize) -> HighsLp {
        let mut state: u64 = 12345;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (10 + (state >> 33) % 91) as f64
        };
        let weights: Vec<[f64; N_ROWS]> = (0..n).map(|_| [next(), next(), next()]).collect();
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = N_ROWS as HighsInt;
        lp.col_cost_ = weights
            .iter()
            .map(|w| -((w[0] + w[1] + w[2]) / 3.0).floor() - 10.0)
            .collect();
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf; N_ROWS];
        lp.row_upper_ = (0..N_ROWS)
            .map(|i| (weights.iter().map(|w| w[i]).sum::<f64>() / 2.0).floor())
            .collect();
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = N_ROWS as HighsInt;
        lp.a_matrix_.start_ = (0..=n).map(|j| (N_ROWS * j) as HighsInt).collect();
        lp.a_matrix_.index_ = (0..n).flat_map(|_| [0, 1, 2]).collect();
        lp.a_matrix_.value_ = weights.iter().flat_map(|w| *w).collect();
        lp
    }

    /// The knapsack optimum by enumeration (for the small instance).
    fn brute_force_optimum(lp: &HighsLp) -> f64 {
        let n = lp.num_col_ as usize;
        let mut best = 0.0f64;
        for mask in 0u32..(1 << n) {
            let mut act = [0.0f64; N_ROWS];
            let mut v = 0.0;
            for j in (0..n).filter(|j| mask & (1 << j) != 0) {
                for (i, a) in act.iter_mut().enumerate() {
                    *a += lp.a_matrix_.value_[N_ROWS * j + i];
                }
                v -= lp.col_cost_[j];
            }
            if act.iter().zip(&lp.row_upper_).all(|(a, u)| a <= u) && v > best {
                best = v;
            }
        }
        -best
    }

    fn base_options() -> HighsOptions {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.presolve = kHighsOffString.to_string();
        // No primal heuristics, so the root does not close the gap and the tree branches.
        options.mip_heuristic_effort = 0.0;
        options.mip_heuristic_run_feasibility_jump = false;
        options.mip_heuristic_run_rens = false;
        options.mip_heuristic_run_rins = false;
        options
    }

    #[test]
    fn s2_pool_to_leq() {
        // x0 + 2x1 >= 3 gives vals [-1,-2], rhs -3.0
        let mut row = S2ConflictRow::default();
        row.vals = vec![1.0, 2.0, 0.0];
        row.inds = vec![0, 1];
        row.lhs = 3.0;

        let (inds, vals, rhs) = to_leq(&row);
        assert_eq!(inds, vec![0, 1]);
        assert_eq!(vals, vec![-1.0, -2.0]);
        assert_eq!(rhs, -3.0);
    }

    #[test]
    fn s2_pool_integral() {
        // integer columns with vals [-1,-2] give true
        let inds = vec![0, 1];
        let vals = vec![-1.0, -2.0];
        assert!(cut_is_integral_by(|_| true, 1e-6, &inds, &vals));

        // vals [-1, -1.5] give false
        let vals2 = vec![-1.0, -1.5];
        assert!(!cut_is_integral_by(|_| true, 1e-6, &inds, &vals2));
    }

    #[test]
    fn s2_pool_debugsol_violation() {
        // The cut x0 <= 0 with x = [1.0] gives viol 1.0
        // In <= form: 1*x0 <= 0, so inds=[0], vals=[1.0], rhs=0.0
        // viol = -rhs + vals[0]*x[inds[0]] = 0 + 1.0*1.0 = 1.0
        let x = vec![1.0];
        let inds = vec![0];
        let vals = vec![1.0];
        let rhs = 0.0;
        let viol = s2_debugsol_viol(&x, &inds, &vals, rhs);
        assert!((viol - 1.0).abs() < 1e-12);
        // counted as a violation (viol > feastol 1e-6)
        assert!(viol > 1e-6);
    }

    #[test]
    fn s2_pool_on_knapsack() {
        let handle = std::thread::spawn(|| {
            initialize_scheduler(1);
            let lp = knapsack_lp(20);
            let expected = brute_force_optimum(&lp);

            let mut options = base_options();
            options.mip_conflict_resolution = true;

            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();

            let obj = ms.solution_objective_;
            assert!(
                (obj - expected).abs() < 1e-6,
                "objective {} != brute force {}",
                obj,
                expected
            );
            assert!(
                ms.mipdata_.conflict_resolution.stats().calls > 0,
                "s2_conflict_analysis was never called"
            );
            println!(
                "rows_found={} rows_added={}",
                ms.mipdata_.conflict_resolution.stats().rows_found,
                ms.mipdata_.conflict_resolution.stats().rows_added
            );
            // Verified: n=20, LCG seed 12345 gives rows_added = 52 (flag on, presolve off).
            assert!(
                ms.mipdata_.conflict_resolution.stats().rows_added > 0,
                "expected rows_added > 0 with n=20 seed=12345"
            );

            drop(ms);
            HighsTaskExecutor::shutdown(true);
        });
        handle.join().unwrap();
    }

    #[test]
    fn s2_pool_off_knapsack() {
        let handle = std::thread::spawn(|| {
            initialize_scheduler(1);

            // Default options (flag off by default)
            let lp1 = knapsack_lp(20);
            let options1 = base_options();
            let solution1 = HighsSolution::default();
            let mut callback1 = HighsCallback::default();
            let mut ms1 =
                HighsMipSolver::new(&mut callback1, &options1, &lp1, &solution1, false, 0);
            ms1.run();

            // Explicit flag = false
            let lp2 = knapsack_lp(20);
            let mut options2 = base_options();
            options2.mip_conflict_resolution = false;
            let solution2 = HighsSolution::default();
            let mut callback2 = HighsCallback::default();
            let mut ms2 =
                HighsMipSolver::new(&mut callback2, &options2, &lp2, &solution2, false, 0);
            ms2.run();

            // Bit-identical objective, num_nodes, total_lp_iterations
            assert_eq!(
                ms1.solution_objective_.to_bits(),
                ms2.solution_objective_.to_bits(),
                "objective bits differ"
            );
            assert_eq!(
                ms1.mipdata_.num_nodes, ms2.mipdata_.num_nodes,
                "num_nodes differ"
            );
            assert_eq!(
                ms1.mipdata_.total_lp_iterations, ms2.mipdata_.total_lp_iterations,
                "total_lp_iterations differ"
            );

            drop(ms1);
            drop(ms2);
            HighsTaskExecutor::shutdown(true);
        });
        handle.join().unwrap();
    }
}
