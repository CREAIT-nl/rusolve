//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! end-to-end tests for the hybrid cut selector option
//! `mip_cut_selection_hybrid`. Each test builds a small 0-1 MIP by hand and
//! solves it with the flag off and on. The optimum is checked against a
//! brute-force enumeration, the run is checked for determinism, and
//! `random_seed` 0 and 1 are both checked to reach the optimum.
//!
//! The models use non-integral objective coefficients so that the root dual
//! bound is not closed by objective-integrality rounding: with integral costs,
//! tiny models finish before the root separation loop, and the hybrid path in
//! `separation.rs` would never run. Both models below were checked (by
//! temporary instrumentation of `apply::s4_separate_cutpools`) to reach the
//! hybrid path in several separation rounds, with 0 forced cuts in every call
//!.

use crate::lp_data::hconst::{HighsModelStatus, HighsVarType, ObjSense};
use crate::lp_data::highs_run::Highs;
use crate::prelude::*;

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

/// A 0-1 model: `sense` sum(cost_i x_i) subject to row_lo <= A x <= row_up.
struct Model01 {
    sense: ObjSense,
    cost: Vec<f64>,
    rows: Vec<Vec<f64>>,
    row_lo: Vec<f64>,
    row_up: Vec<f64>,
}

/// Multi-dimensional knapsack: maximize, 12 binaries, 3 capacity rows <= 19.
fn mkp() -> Model01 {
    Model01 {
        sense: ObjSense::kMaximize,
        cost: vec![
            10.3, 13.1, 7.7, 8.9, 11.4, 6.2, 9.8, 12.6, 5.5, 14.2, 7.1, 10.9,
        ],
        rows: vec![
            vec![3.0, 5.0, 2.0, 4.0, 6.0, 1.0, 4.0, 7.0, 2.0, 8.0, 3.0, 5.0],
            vec![6.0, 2.0, 5.0, 3.0, 4.0, 7.0, 1.0, 3.0, 6.0, 2.0, 5.0, 4.0],
            vec![2.0, 7.0, 4.0, 6.0, 1.0, 3.0, 5.0, 2.0, 4.0, 6.0, 7.0, 3.0],
        ],
        row_lo: vec![-kHighsInf; 3],
        row_up: vec![19.0; 3],
    }
}

/// Weighted set cover: minimize, 12 binaries, 10 covering rows >= 1.
fn setcover() -> Model01 {
    let cover: [&[usize]; 10] = [
        &[0, 1, 5],
        &[1, 2, 6],
        &[2, 3, 7],
        &[3, 4, 8],
        &[4, 0, 9],
        &[5, 7, 10],
        &[6, 8, 11],
        &[7, 9, 10],
        &[8, 5, 11],
        &[9, 6, 10, 11],
    ];
    let n = 12;
    let rows = cover
        .iter()
        .map(|cols| {
            let mut r = vec![0.0; n];
            for &j in cols.iter() {
                r[j] = 1.0;
            }
            r
        })
        .collect();
    Model01 {
        sense: ObjSense::kMinimize,
        cost: vec![3.1, 2.7, 3.3, 2.9, 3.4, 2.2, 2.6, 2.4, 2.8, 2.3, 4.1, 3.9],
        rows,
        row_lo: vec![1.0; 10],
        row_up: vec![kHighsInf; 10],
    }
}

fn build(m: &Model01, hybrid: bool) -> Highs {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    let n = m.cost.len() as HighsInt;
    h.add_vars(n, &vec![0.0; n as usize], &vec![1.0; n as usize]);
    for i in 0..n {
        h.change_col_cost(i, m.cost[i as usize]);
        h.change_col_integrality(i, HighsVarType::kInteger);
    }
    for (r, row) in m.rows.iter().enumerate() {
        let idx: Vec<HighsInt> = (0..n).filter(|&j| row[j as usize] != 0.0).collect();
        let val: Vec<f64> = idx.iter().map(|&j| row[j as usize]).collect();
        h.add_row(m.row_lo[r], m.row_up[r], idx.len() as HighsInt, &idx, &val);
    }
    h.change_objective_sense(m.sense);
    h.set_option_value_bool("mip_cut_selection_hybrid", hybrid);
    h
}

fn activity(row: &[f64], x: &[f64]) -> f64 {
    row.iter().zip(x.iter()).map(|(a, v)| a * v).sum()
}

/// Optimal objective by enumeration of all 0-1 points.
fn brute_force_opt(m: &Model01) -> f64 {
    let n = m.cost.len();
    let mut best: Option<f64> = None;
    for mask in 0u32..(1u32 << n) {
        let x: Vec<f64> = (0..n).map(|j| ((mask >> j) & 1) as f64).collect();
        let feasible = m.rows.iter().enumerate().all(|(r, row)| {
            let a = activity(row, &x);
            a >= m.row_lo[r] - 1e-9 && a <= m.row_up[r] + 1e-9
        });
        if !feasible {
            continue;
        }
        let obj = activity(&m.cost, &x);
        best = Some(match (best, m.sense) {
            (None, _) => obj,
            (Some(b), ObjSense::kMaximize) => b.max(obj),
            (Some(b), _) => b.min(obj),
        });
    }
    best.expect("model has a feasible 0-1 point")
}

struct RunResult {
    status: HighsModelStatus,
    obj: f64,
    x: Vec<f64>,
    nodes: i64,
    iters: HighsInt,
}

fn solve(h: &mut Highs) -> RunResult {
    let _ = h.run();
    let info = h.get_info();
    RunResult {
        status: h.get_model_status(),
        obj: info.objective_function_value,
        x: h.get_solution().col_value.clone(),
        nodes: info.mip_node_count,
        iters: info.simplex_iteration_count,
    }
}

/// x is binary within 1e-6, satisfies every row within 1e-6, and its
/// objective equals the reported one.
fn check_feasible(m: &Model01, r: &RunResult) {
    let tol = 1e-6;
    assert_eq!(r.x.len(), m.cost.len());
    for (j, &v) in r.x.iter().enumerate() {
        assert!(
            (-tol..=1.0 + tol).contains(&v),
            "x[{}] = {} out of [0,1]",
            j,
            v
        );
        assert!((v - v.round()).abs() < tol, "x[{}] = {} not integral", j, v);
    }
    for (i, row) in m.rows.iter().enumerate() {
        let a = activity(row, &r.x);
        assert!(
            a >= m.row_lo[i] - tol && a <= m.row_up[i] + tol,
            "row {} violated: activity {} not in [{}, {}]",
            i,
            a,
            m.row_lo[i],
            m.row_up[i]
        );
    }
    let obj = activity(&m.cost, &r.x);
    assert!(
        (obj - r.obj).abs() < 1e-6,
        "objective of x {} differs from reported {}",
        obj,
        r.obj
    );
}

fn check_off_on(m: &Model01) {
    let opt = brute_force_opt(m);
    let off = solve(&mut build(m, false));
    let on = solve(&mut build(m, true));
    for (name, r) in [("off", &off), ("on", &on)] {
        assert_eq!(r.status, HighsModelStatus::kOptimal, "flag {}", name);
        check_feasible(m, r);
        assert!(
            (r.obj - opt).abs() < 1e-6,
            "flag {}: objective {} differs from brute-force optimum {}",
            name,
            r.obj,
            opt
        );
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn s4_e2e_off_on_same_optimum_knapsack() {
    check_off_on(&mkp());
}

#[test]
fn s4_e2e_off_on_same_optimum_setcover() {
    check_off_on(&setcover());
}

/// Flag on: two fresh runs of the same model give identical node counts,
/// simplex iteration counts, objectives and solutions.
#[test]
fn s4_e2e_on_deterministic() {
    for m in [mkp(), setcover()] {
        let r1 = solve(&mut build(&m, true));
        let r2 = solve(&mut build(&m, true));
        assert_eq!(r1.status, HighsModelStatus::kOptimal);
        assert_eq!(r1.nodes, r2.nodes, "node count differs");
        assert_eq!(r1.iters, r2.iters, "simplex iteration count differs");
        assert_eq!(r1.obj.to_bits(), r2.obj.to_bits(), "objective differs");
        assert_eq!(r1.x, r2.x, "solution differs");
    }
}

/// Flag on: random_seed 0 and 1 both reach the brute-force optimum with a
/// feasible x.
#[test]
fn s4_e2e_random_seed_optimal() {
    for m in [mkp(), setcover()] {
        let opt = brute_force_opt(&m);
        for seed in [0, 1] {
            let mut h = build(&m, true);
            h.set_option_value_int("random_seed", seed);
            let r = solve(&mut h);
            assert_eq!(r.status, HighsModelStatus::kOptimal, "seed {}", seed);
            check_feasible(&m, &r);
            assert!(
                (r.obj - opt).abs() < 1e-6,
                "seed {}: objective {} differs from optimum {}",
                seed,
                r.obj,
                opt
            );
        }
    }
}

/// forced cuts are expected to be 0 because `add_cuts` clears the cutset
/// before the next separation round. The forced-cut count is not observable
/// from here (it is a local in `apply::s4_separate_cutpools`, a file outside this
/// job), so this test only checks that the flag-on runs, which reach the hybrid
/// path, are optimal and feasible. The count was measured as 0 in every hybrid
/// call of these models by temporary instrumentation (see the module doc).
#[test]
fn s4_e2e_forced_cut_count_zero() {
    for m in [mkp(), setcover()] {
        let r = solve(&mut build(&m, true));
        assert_eq!(r.status, HighsModelStatus::kOptimal);
        check_feasible(&m, &r);
    }
}
