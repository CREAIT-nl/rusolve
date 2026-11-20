//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! End-to-end tests for the SCIP oneopt and crossover heuristics.

use super::*;
use crate::lp_data::hconst::HighsVarType;
use crate::mip::mip_solver::HighsMipSolver;

/// Build a 4-variable knapsack LP (maximize, all integer).
/// max 3x0 + 2x1 + 4x2 + 1x3
/// s.t. 2x0 + x1 + 3x2 + x3 <= 4
///      x0 + 2x1 + x2 + 2x3 <= 5
///      0 <= xi <= 1, xi integer
fn e2e_knapsack_lp() -> HighsLp {
    let mut lp = HighsLp::default();
    lp.num_col_ = 4;
    lp.num_row_ = 2;
    lp.sense_ = ObjSense::kMaximize;
    lp.offset_ = 0.0;
    lp.col_cost_ = vec![3.0, 2.0, 4.0, 1.0];
    lp.col_lower_ = vec![0.0; 4];
    lp.col_upper_ = vec![1.0; 4];
    lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
    lp.row_upper_ = vec![4.0, 5.0];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.start_ = vec![0, 2, 4, 6, 8];
    lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1, 0, 1];
    lp.a_matrix_.value_ = vec![2.0, 1.0, 1.0, 2.0, 3.0, 1.0, 1.0, 2.0];
    lp.integrality_ = vec![HighsVarType::kInteger; 4];
    lp
}

/// Build a 6-variable set-cover LP (minimize, all integer).
/// min 2x0 + 3x1 + 1x2 + 4x3 + 2x4 + 1x5
/// s.t. x0 + x1 + x3 >= 1
///      x1 + x2 + x4 >= 1
///      x0 + x2 + x5 >= 1
///      x3 + x4 + x5 >= 1
///      0 <= xi <= 1, xi integer
fn e2e_set_cover_lp() -> HighsLp {
    let mut lp = HighsLp::default();
    lp.num_col_ = 6;
    lp.num_row_ = 4;
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0.0;
    lp.col_cost_ = vec![2.0, 3.0, 1.0, 4.0, 2.0, 1.0];
    lp.col_lower_ = vec![0.0; 6];
    lp.col_upper_ = vec![1.0; 6];
    lp.row_lower_ = vec![1.0, 1.0, 1.0, 1.0];
    lp.row_upper_ = vec![kHighsInf; 4];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    // col 0: rows 0,2
    // col 1: rows 0,1
    // col 2: rows 1,2
    // col 3: rows 0,3
    // col 4: rows 1,3
    // col 5: rows 2,3
    lp.a_matrix_.start_ = vec![0, 2, 4, 6, 8, 10, 12];
    lp.a_matrix_.index_ = vec![0, 2, 0, 1, 1, 2, 0, 3, 1, 3, 2, 3];
    lp.a_matrix_.value_ = vec![1.0; 12];
    lp.integrality_ = vec![HighsVarType::kInteger; 6];
    lp
}

/// Build a mixed facility-location LP with continuous flows.
/// 3 facilities (integer open/close), 4 customers (continuous flow).
/// min 5y0 + 3y1 + 4y2 + 1x0 + 2x1 + 1x2 + 3x3
/// s.t. x0 + x1 + x2 + x3 = 10  (demand)
///      x0 <= 5*y0
///      x1 + x2 <= 8*y1
///      x3 <= 6*y2
///      y0 + y1 + y2 <= 2  (at most 2 facilities open)
///      0 <= y_i <= 1, y integer
///      0 <= x_j <= 10, x continuous
fn e2e_facility_location_lp() -> HighsLp {
    let mut lp = HighsLp::default();
    lp.num_col_ = 7;
    lp.num_row_ = 5;
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0.0;
    lp.col_cost_ = vec![5.0, 3.0, 4.0, 1.0, 2.0, 1.0, 3.0];
    lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    lp.col_upper_ = vec![1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 10.0];
    lp.row_lower_ = vec![10.0, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf];
    lp.row_upper_ = vec![10.0, 0.0, 0.0, 0.0, 2.0];
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    // col 0 (y0): row 1 (x0-5y0<=0 -> coeff -5), row 4 (y0+y1+y2<=2 -> coeff 1)
    // col 1 (y1): row 2 (x1+x2-8y1<=0 -> coeff -8), row 4 (coeff 1)
    // col 2 (y2): row 3 (x3-6y2<=0 -> coeff -6), row 4 (coeff 1)
    // col 3 (x0): row 0 (coeff 1), row 1 (coeff 1)
    // col 4 (x1): row 0 (coeff 1), row 2 (coeff 1)
    // col 5 (x2): row 0 (coeff 1), row 2 (coeff 1)
    // col 6 (x3): row 0 (coeff 1), row 3 (coeff 1)
    lp.a_matrix_.start_ = vec![0, 2, 4, 6, 8, 10, 12, 14];
    lp.a_matrix_.index_ = vec![1, 4, 2, 4, 3, 4, 0, 1, 0, 2, 0, 2, 0, 3];
    lp.a_matrix_.value_ = vec![
        -5.0, 1.0, -8.0, 1.0, -6.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ];
    lp.integrality_ = vec![
        HighsVarType::kInteger,
        HighsVarType::kInteger,
        HighsVarType::kInteger,
        HighsVarType::kContinuous,
        HighsVarType::kContinuous,
        HighsVarType::kContinuous,
        HighsVarType::kContinuous,
    ];
    lp
}

/// Build an n-item, m-constraint 0/1 knapsack (maximize) from a fixed LCG seed.
/// Weights are in [10, 99], each capacity is half its row's weight sum, and the profit of
/// item j is its mean weight plus a value in [0, 49].
fn e2e_mkp_lp(n: usize, m: usize, seed: u64) -> HighsLp {
    let mut s: u64 = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let mut rnd = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 33) % 1000
    };
    let mut lp = HighsLp::default();
    lp.num_col_ = n as HighsInt;
    lp.num_row_ = m as HighsInt;
    lp.sense_ = ObjSense::kMaximize;
    lp.offset_ = 0.0;
    let mut a = vec![vec![0.0; n]; m];
    let mut rowsum = vec![0.0; m];
    for j in 0..n {
        for r in 0..m {
            a[r][j] = (10 + rnd() % 90) as f64;
            rowsum[r] += a[r][j];
        }
    }
    lp.col_cost_ = (0..n)
        .map(|j| (0..m).map(|r| a[r][j]).sum::<f64>() / m as f64 + (rnd() % 50) as f64)
        .collect();
    lp.col_lower_ = vec![0.0; n];
    lp.col_upper_ = vec![1.0; n];
    lp.row_lower_ = vec![-kHighsInf; m];
    lp.row_upper_ = rowsum.iter().map(|v| (v * 0.5).floor()).collect();
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    let mut start = vec![0 as HighsInt];
    for j in 0..n {
        for r in 0..m {
            lp.a_matrix_.index_.push(r as HighsInt);
            lp.a_matrix_.value_.push(a[r][j]);
        }
        start.push(lp.a_matrix_.index_.len() as HighsInt);
    }
    lp.a_matrix_.start_ = start;
    lp.integrality_ = vec![HighsVarType::kInteger; n];
    lp
}

/// The 40-item, 10-constraint knapsack used for crossover: it takes about 1500 nodes and
/// its best solution is found more than 200 nodes (nwaitingnodes) before the end.
fn e2e_large_knapsack_lp() -> HighsLp {
    e2e_mkp_lp(40, 10, 1)
}

/// The 30-item, 5-constraint knapsack used for oneopt. With the HiGHS primal heuristics
/// off (see `e2e_heuristics_off`), the incumbent at the end of the root node leaves slack
/// that 1-opt shifts can use.
fn e2e_oneopt_lp() -> HighsLp {
    e2e_mkp_lp(30, 5, 0)
}

/// A pool entry as compared by the determinism test: index, objective bits, x bits, origin
/// and node number.
type E2ePoolEntry = (i64, u64, Vec<u64>, S6Origin, i64);

struct E2eResult {
    /// mipdata_.upper_bound (minimization sense of the solved model, offset 0 here).
    obj: f64,
    nodes: i64,
    lp_iters: i64,
    x: Vec<f64>,
    xover_ncalls: i64,
    oneopt_ncalls: i64,
    pool: Vec<E2ePoolEntry>,
}

/// Turns off the HiGHS primal heuristics so that the root incumbent is not 1-opt.
fn e2e_heuristics_off(options: &mut HighsOptions) {
    options.mip_heuristic_effort = 0.0;
    options.mip_heuristic_run_feasibility_jump = false;
    options.mip_heuristic_run_rins = false;
    options.mip_heuristic_run_rens = false;
    options.mip_heuristic_run_root_reduced_cost = false;
    options.mip_heuristic_run_zi_round = false;
    options.mip_heuristic_run_shifting = false;
}

fn e2e_solve(lp: &HighsLp, oneopt: bool, crossover: bool, seed: HighsInt) -> E2eResult {
    e2e_solve_with(lp, oneopt, crossover, seed, |_| {})
}

/// Solves `lp` with one thread and presolve off, so that `mipdata_.incumbent` and the pool
/// live in the space of `lp` (same columns) and can be checked against it.
fn e2e_solve_with(
    lp: &HighsLp,
    oneopt: bool,
    crossover: bool,
    seed: HighsInt,
    tweak: fn(&mut HighsOptions),
) -> E2eResult {
    let mut lp = lp.clone();
    // The builders fill only the CSC vectors; the matrix dimensions must match the LP, or
    // `a_matrix_.num_nz()` is 0 and presolve's entry assert (HPresolve.cpp:4570) fires.
    lp.a_matrix_.num_col_ = lp.num_col_;
    lp.a_matrix_.num_row_ = lp.num_row_;
    let handle = std::thread::spawn(move || {
        crate::parallel::highs_parallel::initialize_scheduler(1);
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.presolve = "off".to_string();
        options.mip_heuristic_run_oneopt = oneopt;
        options.mip_heuristic_run_crossover = crossover;
        options.random_seed = seed;
        tweak(&mut options);
        let mut callback = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
        ms.run();
        let s6 = ms.mipdata_.s6.as_ref();
        assert_eq!(s6.is_some(), oneopt || crossover, "S6 state presence");
        let result = E2eResult {
            obj: ms.mipdata_.upper_bound,
            nodes: ms.mipdata_.num_nodes,
            lp_iters: ms.mipdata_.total_lp_iterations,
            x: ms.mipdata_.incumbent.clone(),
            xover_ncalls: s6.map_or(0, |s| s.xover_ncalls),
            oneopt_ncalls: s6.map_or(0, |s| s.oneopt_ncalls),
            pool: s6.map_or(vec![], |s| {
                s.pool
                    .sols
                    .iter()
                    .map(|sol| {
                        (
                            sol.index,
                            sol.obj.to_bits(),
                            sol.x.iter().map(|v| v.to_bits()).collect(),
                            sol.origin,
                            sol.nodenum,
                        )
                    })
                    .collect()
            }),
        };
        drop(ms);
        crate::parallel::highs_task_executor::HighsTaskExecutor::shutdown(true);
        result
    });
    handle.join().unwrap()
}

/// Asserts that x is feasible for `lp` and that its objective (in the sense of `lp`) is
/// `sense * obj`.
fn check_feasible(lp: &HighsLp, x: &[f64], obj: f64) {
    let nc = lp.num_col_ as usize;
    let nr = lp.num_row_ as usize;
    assert_eq!(x.len(), nc, "solution length mismatch");
    let mut activity = vec![0.0; nr];
    let mut cx = lp.offset_;
    for j in 0..nc {
        assert!(
            x[j] >= lp.col_lower_[j] - 1e-6,
            "col {} below lower bound",
            j
        );
        assert!(
            x[j] <= lp.col_upper_[j] + 1e-6,
            "col {} above upper bound",
            j
        );
        if lp.integrality_[j] == HighsVarType::kInteger {
            assert!(
                (x[j] - x[j].round()).abs() < 1e-6,
                "col {} not integer: {}",
                j,
                x[j]
            );
        }
        cx += lp.col_cost_[j] * x[j];
        for k in lp.a_matrix_.start_[j] as usize..lp.a_matrix_.start_[j + 1] as usize {
            activity[lp.a_matrix_.index_[k] as usize] += lp.a_matrix_.value_[k] * x[j];
        }
    }
    for i in 0..nr {
        assert!(
            activity[i] >= lp.row_lower_[i] - 1e-6 && activity[i] <= lp.row_upper_[i] + 1e-6,
            "row {} activity {} outside [{}, {}]",
            i,
            activity[i],
            lp.row_lower_[i],
            lp.row_upper_[i]
        );
    }
    let user_obj = (lp.sense_ as i32) as f64 * obj;
    assert!(
        (cx - user_obj).abs() < 1e-6,
        "objective of x {} differs from the solver objective {}",
        cx,
        user_obj
    );
}

/// Solves `lp` with each flag combination (off/off, on/off, off/on, on/on) and checks that
/// every run reaches `optimum` (in the sense of `lp`) with a feasible x.
fn e2e_all_flags(name: &str, lp: &HighsLp, optimum: f64) {
    for &(oneopt, crossover) in &[(false, false), (true, false), (false, true), (true, true)] {
        let r = e2e_solve(lp, oneopt, crossover, 0);
        let user_obj = (lp.sense_ as i32) as f64 * r.obj;
        assert!(
            (user_obj - optimum).abs() < 1e-6,
            "{} oneopt={} crossover={}: objective {} != {}",
            name,
            oneopt,
            crossover,
            user_obj,
            optimum
        );
        check_feasible(lp, &r.x, r.obj);
    }
}

#[test]
fn s6_e2e_knapsack_all_flags() {
    // x1 + x2 (or x0 + x1 + x3) gives 6.
    e2e_all_flags("knapsack", &e2e_knapsack_lp(), 6.0);
}

#[test]
fn s6_e2e_set_cover_all_flags() {
    // x0 + x2 + x5 gives 2 + 1 + 1.
    e2e_all_flags("set cover", &e2e_set_cover_lp(), 4.0);
}

#[test]
fn s6_e2e_facility_location_all_flags() {
    // y0 and y1 open, x0 = 2, x2 = 8: 5 + 3 + 2 + 8.
    e2e_all_flags("facility location", &e2e_facility_location_lp(), 18.0);
}

#[test]
fn s6_e2e_determinism_flags_on() {
    for lp in [e2e_knapsack_lp(), e2e_large_knapsack_lp()] {
        let r1 = e2e_solve(&lp, true, true, 0);
        let r2 = e2e_solve(&lp, true, true, 0);
        assert_eq!(r1.nodes, r2.nodes, "node count differs");
        assert_eq!(r1.lp_iters, r2.lp_iters, "LP iterations differ");
        assert_eq!(r1.pool, r2.pool, "pool contents differ");
        assert_eq!(r1.obj.to_bits(), r2.obj.to_bits(), "objective differs");
        assert_eq!(r1.xover_ncalls, r2.xover_ncalls, "crossover calls differ");
        assert_eq!(r1.oneopt_ncalls, r2.oneopt_ncalls, "oneopt calls differ");
        assert!(!r1.pool.is_empty(), "pool is empty");
    }
}

#[test]
fn s6_e2e_random_seed_0_vs_1_optimal() {
    for lp in [e2e_knapsack_lp(), e2e_large_knapsack_lp()] {
        let off = e2e_solve(&lp, false, false, 0);
        for seed in [0, 1] {
            let r = e2e_solve(&lp, true, true, seed);
            assert!(
                (r.obj - off.obj).abs() < 1e-6,
                "seed {} objective {} differs from the flags-off optimum {}",
                seed,
                r.obj,
                off.obj
            );
            check_feasible(&lp, &r.x, r.obj);
        }
    }
}

#[test]
fn s6_e2e_crossover_runs_large_knapsack() {
    let lp = e2e_large_knapsack_lp();
    let off = e2e_solve(&lp, false, false, 0);
    let r = e2e_solve(&lp, false, true, 0);
    assert!(
        r.nodes >= 200,
        "too few nodes for crossover: {} (nwaitingnodes is 200)",
        r.nodes
    );
    assert!(
        r.xover_ncalls > 0,
        "crossover did not run: xover_ncalls = {}",
        r.xover_ncalls
    );
    assert!(
        (r.obj - off.obj).abs() < 1e-6,
        "crossover objective {} differs from the flags-off optimum {}",
        r.obj,
        off.obj
    );
    check_feasible(&lp, &r.x, r.obj);
}

#[test]
fn s6_e2e_oneopt_finds_improvement() {
    let lp = e2e_oneopt_lp();
    let off = e2e_solve_with(&lp, false, false, 0, e2e_heuristics_off);
    let on = e2e_solve_with(&lp, true, false, 0, e2e_heuristics_off);
    assert!(
        (off.obj - on.obj).abs() < 1e-6,
        "oneopt objective {} differs from the flags-off optimum {}",
        on.obj,
        off.obj
    );
    check_feasible(&lp, &on.x, on.obj);
    assert!(on.oneopt_ncalls > 0, "oneopt was not called");
    // A oneopt solution in the pool is strictly better than every pool solution found
    // before it, so it improved on the incumbent it started from.
    let oneopt_sols: Vec<&E2ePoolEntry> =
        on.pool.iter().filter(|e| e.3 == S6Origin::Oneopt).collect();
    assert!(!oneopt_sols.is_empty(), "oneopt found no solution");
    for e in &oneopt_sols {
        let obj = f64::from_bits(e.1);
        let x: Vec<f64> = e.2.iter().map(|&b| f64::from_bits(b)).collect();
        check_feasible(&lp, &x, obj);
        for earlier in on.pool.iter().filter(|p| p.0 < e.0) {
            assert!(
                obj < f64::from_bits(earlier.1),
                "oneopt solution {} (obj {}) does not improve on solution {} (obj {})",
                e.0,
                obj,
                earlier.0,
                f64::from_bits(earlier.1)
            );
        }
    }
}
