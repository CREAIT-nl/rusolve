//! Tests for `ls_int_lp_repair`.

use super::ls_int_lp_repair::*;
use crate::mip::Slot;
use crate::mip::mip_worker::HighsMipWorker;
use crate::prelude::*;

/// Capacitated facility location: `nf` facilities (binary y_i, fixed cost f_i, capacity
/// u_i), `nc` customers (demand d_j), continuous shipments x_ij >= 0 at cost c_ij.
///   min sum f_i y_i + sum c_ij x_ij
///   s.t. sum_i x_ij = d_j            (demand rows 0..nc)
///        sum_j x_ij - u_i y_i <= 0   (capacity rows nc..nc+nf)
///        sum_i y_i >= 1              (last row, integer-only)
/// Columns: x_ij at i * nc + j, then y_i at nf * nc + i.
fn facility_location(nf: usize, nc: usize, seed: u64) -> HighsLp {
    let mut s = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let mut next = |lo: u64, hi: u64| {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        lo + (s >> 33) % (hi - lo + 1)
    };
    let d: Vec<f64> = (0..nc).map(|_| next(5, 20) as f64).collect();
    let total: f64 = d.iter().sum();
    let u: Vec<f64> = (0..nf)
        .map(|_| (total * next(30, 60) as f64 / 100.0).ceil())
        .collect();
    let f: Vec<f64> = (0..nf).map(|_| next(50, 150) as f64).collect();
    let cost: Vec<f64> = (0..nf * nc).map(|_| next(1, 20) as f64).collect();

    let nx = nf * nc;
    let ncol = nx + nf;
    let nrow = nc + nf + 1;
    let mut lp = HighsLp::default();
    lp.num_col_ = ncol as HighsInt;
    lp.num_row_ = nrow as HighsInt;
    lp.col_cost_ = cost.iter().copied().chain(f.iter().copied()).collect();
    lp.col_lower_ = vec![0.0; ncol];
    lp.col_upper_ = (0..ncol)
        .map(|k| if k < nx { kHighsInf } else { 1.0 })
        .collect();
    lp.row_lower_ = d
        .iter()
        .copied()
        .chain(vec![-kHighsInf; nf])
        .chain([1.0])
        .collect();
    lp.row_upper_ = d
        .iter()
        .copied()
        .chain(vec![0.0; nf])
        .chain([kHighsInf])
        .collect();
    lp.integrality_ = (0..ncol)
        .map(|k| {
            if k < nx {
                HighsVarType::kContinuous
            } else {
                HighsVarType::kInteger
            }
        })
        .collect();
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = ncol as HighsInt;
    lp.a_matrix_.num_row_ = nrow as HighsInt;
    let mut start = vec![0 as HighsInt];
    let mut index = Vec::new();
    let mut value = Vec::new();
    for i in 0..nf {
        for j in 0..nc {
            index.push(j as HighsInt);
            value.push(1.0);
            index.push((nc + i) as HighsInt);
            value.push(1.0);
            start.push(index.len() as HighsInt);
        }
    }
    for i in 0..nf {
        index.push((nc + i) as HighsInt);
        value.push(-u[i]);
        index.push((nc + nf) as HighsInt);
        value.push(1.0);
        start.push(index.len() as HighsInt);
    }
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp
}

/// Solver data for `lp` without presolve (model_ == lp), as the oneopt tests build it,
/// plus domain propagation, the master worker and `integer_cols`, which run_setup / run()
/// provide.
fn make_ms(lp: &HighsLp, on: bool) -> HighsMipSolver {
    let options = HighsOptions {
        output_flag: false,
        mip_ls_int_lp_repair: on,
        ..Default::default()
    };
    let mut ms = HighsMipSolver::new(
        &mut HighsCallback::default(),
        &options,
        lp,
        &HighsSolution::default(),
        false,
        0,
    );
    ms.new_data();
    ms.init();
    // Row-wise matrix and row activities of the global domain, as run_setup builds them.
    ms.setup_domain_propagation();
    crate::detach!(ms.mipdata_.domains[0], |d| d.compute_row_activities(&ms));
    ms.mipdata_.integer_cols = (0..lp.num_col_)
        .filter(|&j| lp.integrality_[j as usize] == HighsVarType::kInteger)
        .collect();
    let master = HighsMipWorker::new(&mut ms, 0, 0, 0, 0, 0);
    ms.mipdata_.workers.push(Slot::new(master));
    ms
}

fn objective(lp: &HighsLp, x: &[f64]) -> f64 {
    lp.col_cost_.iter().zip(x).map(|(cj, v)| cj * v).sum()
}

/// A feasible point with every facility open and each customer served from the most
/// expensive facility that still has room: right integers, poor continuous values.
fn expensive_feasible_point(lp: &HighsLp, nf: usize, nc: usize) -> Vec<f64> {
    let mut x = vec![0.0; nf * nc + nf];
    for i in 0..nf {
        x[nf * nc + i] = 1.0;
    }
    let cap: Vec<f64> = (0..nf)
        .map(|i| -lp.a_matrix_.value_[2 * nf * nc + 2 * i])
        .collect();
    let mut room = cap.clone();
    for j in 0..nc {
        let mut left = lp.row_lower_[j];
        let mut order: Vec<usize> = (0..nf).collect();
        order.sort_by(|&a, &b| lp.col_cost_[b * nc + j].total_cmp(&lp.col_cost_[a * nc + j]));
        for i in order {
            let q = left.min(room[i]);
            x[i * nc + j] += q;
            room[i] -= q;
            left -= q;
        }
        assert!(left <= 0.0, "instance too tight for the crafted point");
    }
    x
}

#[test]
fn polish_improves_feasible_ls_point() {
    let (nf, nc) = (4, 10);
    let lp = facility_location(nf, nc, 1);
    let mut ms = make_ms(&lp, true);
    let x = expensive_feasible_point(&lp, nf, nc);
    assert!(ms.try_solution(&x, kSolutionSourceFeasibilityJump));
    let ub_ls = ms.mipdata_.upper_bound;
    assert!((ub_ls - objective(&lp, &x)).abs() < 1e-6);

    let out = repair(&mut ms, &x, kSolutionSourceFeasibilityJump, true);
    assert_eq!(out, RepairOutcome::Improved);
    let ub = ms.mipdata_.upper_bound;
    assert!(ub < ub_ls - 1.0, "polish should lower {ub_ls} (got {ub})");
    // The integers are those of the LS point.
    let inc = &ms.mipdata_.incumbent;
    for i in 0..nf {
        assert_eq!(inc[nf * nc + i], 1.0);
    }
    let st = ms.mipdata_.ls_int_lp_repair.as_ref().unwrap();
    assert_eq!((st.calls, st.lp_solves, st.improved), (1, 1, 1));
    assert!(st.lp_iterations > 0);
    assert_eq!(ms.mipdata_.heuristic_lp_iterations, st.lp_iterations);
}

#[test]
fn repair_turns_infeasible_end_state_into_incumbent() {
    let (nf, nc) = (4, 10);
    let lp = facility_location(nf, nc, 2);
    let mut ms = make_ms(&lp, true);
    // Right integers (all open), continuous part at zero: every demand row is violated.
    let mut x = vec![0.0; nf * nc + nf];
    for i in 0..nf {
        x[nf * nc + i] = 1.0;
    }
    assert!(!ms.try_solution(&x, kSolutionSourceFeasibilityJump));
    assert!(ms.mipdata_.incumbent.is_empty());
    let out = repair(&mut ms, &x, kSolutionSourceFeasibilityJump, false);
    assert_eq!(out, RepairOutcome::Improved);
    assert!(!ms.mipdata_.incumbent.is_empty());
    assert!(ms.mipdata_.upper_bound < kHighsInf);
}

#[test]
fn row_prefilter_rejects_integer_only_violation() {
    let (nf, nc) = (3, 6);
    let lp = facility_location(nf, nc, 3);
    let mut ms = make_ms(&lp, true);
    // All facilities closed violates sum y >= 1, which has no LP column.
    let x = vec![0.0; nf * nc + nf];
    assert_eq!(
        repair(&mut ms, &x, kSolutionSourceFeasibilityJump, false),
        RepairOutcome::RowInfeasible
    );
    let st = ms.mipdata_.ls_int_lp_repair.as_ref().unwrap();
    assert_eq!(st.lp_solves, 0);
    assert!(st.prefilter_effort > 0);
}

#[test]
fn propagation_rejects_infeasible_assignment() {
    let (nf, nc) = (3, 6);
    let lp = facility_location(nf, nc, 4);
    let mut ms = make_ms(&lp, true);
    // One facility open: its capacity (30-60 % of total demand) cannot serve everyone.
    // Every single row is satisfiable on its own, so only propagation/LP can see it.
    let mut x = vec![0.0; nf * nc + nf];
    x[nf * nc] = 1.0;
    let out = repair(&mut ms, &x, kSolutionSourceFeasibilityJump, false);
    assert!(
        matches!(
            out,
            RepairOutcome::PropInfeasible | RepairOutcome::LpInfeasible
        ),
        "{out:?}"
    );
    assert!(ms.mipdata_.incumbent.is_empty());
}

#[test]
fn duplicate_assignment_is_skipped() {
    let (nf, nc) = (4, 10);
    let lp = facility_location(nf, nc, 1);
    let mut ms = make_ms(&lp, true);
    let x = expensive_feasible_point(&lp, nf, nc);
    assert!(ms.try_solution(&x, kSolutionSourceFeasibilityJump));
    assert_eq!(
        repair(&mut ms, &x, kSolutionSourceFeasibilityJump, true),
        RepairOutcome::Improved
    );
    // Same integers, other continuous values (the polished incumbent itself).
    let inc = ms.mipdata_.incumbent.clone();
    assert_eq!(
        repair(&mut ms, &inc, kSolutionSourceLocalMIP, true),
        RepairOutcome::Duplicate
    );
    assert_eq!(ms.mipdata_.ls_int_lp_repair.as_ref().unwrap().calls, 1);
}

#[test]
fn no_cost_on_lp_columns_skips_feasible_point() {
    let (nf, nc) = (3, 6);
    let mut lp = facility_location(nf, nc, 5);
    for k in 0..nf * nc {
        lp.col_cost_[k] = 0.0;
    }
    let mut ms = make_ms(&lp, true);
    let x = expensive_feasible_point(&lp, nf, nc);
    assert_eq!(
        repair(&mut ms, &x, kSolutionSourceFeasibilityJump, true),
        RepairOutcome::NoCost
    );
}

#[test]
fn gate_flag_off_and_low_continuous_share() {
    let (nf, nc) = (3, 6);
    let lp = facility_location(nf, nc, 6);
    let mut ms = make_ms(&lp, false);
    let x = expensive_feasible_point(&lp, nf, nc);
    assert_eq!(
        repair(&mut ms, &x, kSolutionSourceFeasibilityJump, true),
        RepairOutcome::Gated
    );
    assert!(ms.mipdata_.ls_int_lp_repair.is_none());

    // 1 continuous column among 11: below the 10 % share.
    let mut lp2 = HighsLp::default();
    lp2.num_col_ = 11;
    lp2.num_row_ = 1;
    lp2.col_cost_ = vec![1.0; 11];
    lp2.col_lower_ = vec![0.0; 11];
    lp2.col_upper_ = vec![1.0; 11];
    lp2.row_lower_ = vec![1.0];
    lp2.row_upper_ = vec![kHighsInf];
    lp2.integrality_ = (0..11)
        .map(|k| {
            if k == 10 {
                HighsVarType::kContinuous
            } else {
                HighsVarType::kInteger
            }
        })
        .collect();
    lp2.a_matrix_.format_ = MatrixFormat::kColwise;
    lp2.a_matrix_.num_col_ = 11;
    lp2.a_matrix_.num_row_ = 1;
    lp2.a_matrix_.start_ = (0..=11).collect();
    lp2.a_matrix_.index_ = vec![0; 11];
    lp2.a_matrix_.value_ = vec![1.0; 11];
    let mut ms2 = make_ms(&lp2, true);
    assert!(!gate(&ms2));
    let mut x2 = vec![0.0; 11];
    x2[0] = 1.0;
    assert_eq!(
        repair(&mut ms2, &x2, kSolutionSourceFeasibilityJump, true),
        RepairOutcome::Gated
    );
}

#[test]
fn call_budget_is_enforced() {
    let (nf, nc) = (4, 10);
    let lp = facility_location(nf, nc, 7);
    let mut ms = make_ms(&lp, true);
    // Distinct infeasible-but-repairable assignments: subsets with >= 3 of 4 open, and
    // then the remaining 2-open subsets; more distinct points than the call budget.
    let mut outcomes = Vec::new();
    for mask in (1u32..16).rev() {
        let mut x = vec![0.0; nf * nc + nf];
        for i in 0..nf {
            if mask & (1 << i) != 0 {
                x[nf * nc + i] = 1.0;
            }
        }
        outcomes.push(repair(&mut ms, &x, kSolutionSourceFeasibilityJump, false));
    }
    let st = ms.mipdata_.ls_int_lp_repair.as_ref().unwrap();
    assert_eq!(st.calls, 8);
    assert!(outcomes[8..].iter().all(|o| *o == RepairOutcome::Budget));
}

// ---- end to end ----

struct E2e {
    obj: f64,
    stats: Option<LsIntLpRepair>,
}

fn e2e_solve(lp: &HighsLp, on: bool, local_mip: bool) -> E2e {
    let lp = lp.clone();
    std::thread::spawn(move || {
        crate::parallel::highs_parallel::initialize_scheduler(1);
        let options = HighsOptions {
            output_flag: false,
            mip_ls_int_lp_repair: on,
            mip_heuristic_run_local_mip: local_mip,
            ..Default::default()
        };
        let mut callback = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
        ms.run();
        let r = E2e {
            obj: ms.solution_objective_,
            stats: ms.mipdata_.ls_int_lp_repair.as_deref().cloned(),
        };
        drop(ms);
        crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
        r
    })
    .join()
    .unwrap()
}

#[test]
fn e2e_facility_location_same_optimum_and_repair_runs() {
    let lp = facility_location(8, 30, 11);
    let off = e2e_solve(&lp, false, false);
    let on = e2e_solve(&lp, true, false);
    assert!(off.stats.is_none());
    assert!(
        (off.obj - on.obj).abs() <= 1e-6 * off.obj.abs().max(1.0),
        "{} vs {}",
        off.obj,
        on.obj
    );
    let st = on
        .stats
        .expect("FJ runs on this model, so repair is called");
    assert!(st.calls >= 1, "{st:?}");
    assert!(
        st.outcomes
            .iter()
            .all(|(s, _)| *s == kSolutionSourceFeasibilityJump)
    );

    // With Local-MIP on as well, the hook also reaches repair.
    let lm_off = e2e_solve(&lp, false, true);
    let lm_on = e2e_solve(&lp, true, true);
    assert!((lm_off.obj - lm_on.obj).abs() <= 1e-6 * lm_off.obj.abs().max(1.0));
    assert!(lm_on.stats.is_some());
}
