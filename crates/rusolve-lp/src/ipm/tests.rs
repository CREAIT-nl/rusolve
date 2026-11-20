//! Tiny-LP tests of the IPM pipeline: bounded, free columns, ranged rows,
//! maximize, infeasible, unbounded; random LPs against the dual simplex.

use super::*;
use crate::rng::XorShift64;
use crate::simplex::{SimplexConfig, SimplexSolver, SimplexStatus};
use rusolve_core::CscMatrix;

const INF: f64 = f64::INFINITY;

#[allow(clippy::too_many_arguments)]
fn lp(
    m: usize,
    n: usize,
    trip: &[(usize, usize, f64)],
    rl: Vec<f64>,
    ru: Vec<f64>,
    cl: Vec<f64>,
    cu: Vec<f64>,
    c: Vec<f64>,
    maximize: bool,
) -> SparseLp {
    SparseLp::new(
        CscMatrix::from_triplets(m, n, trip),
        rl,
        ru,
        cl,
        cu,
        c,
        maximize,
    )
}

fn both(sp: &SparseLp) -> (IpmSolution, IpmSolution) {
    let cfg = IpmConfig::default();
    let a = solve_sparse(sp, &cfg).unwrap();
    let b = solve_presolved(sp, &cfg, &SimplexConfig::default()).unwrap();
    (a, b)
}

fn assert_opt(sp: &SparseLp, want: f64) {
    let (a, b) = both(sp);
    for (tag, s) in [("no-presolve", &a), ("presolve", &b)] {
        assert_eq!(s.status, IpmStatus::Optimal, "{tag}: {s:?}");
        assert!(
            (s.objective - want).abs() <= 1e-6 * (1.0 + want.abs()),
            "{tag}: objective {} want {want}",
            s.objective
        );
    }
    assert!(a.duals_valid);
    let v = a.verify.unwrap();
    assert!(v.ok(1e-6), "{v:?}");
}

#[test]
fn bounded_two_rows() {
    // min -x - y : x + 2y <= 4, 3x + y <= 6, 0 <= x, y <= 10  →  -2.8
    let sp = lp(
        2,
        2,
        &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 1.0)],
        vec![-INF, -INF],
        vec![4.0, 6.0],
        vec![0.0, 0.0],
        vec![10.0, 10.0],
        vec![-1.0, -1.0],
        false,
    );
    assert_opt(&sp, -2.8);
}

#[test]
fn free_columns() {
    // min y : y - x >= -3, y + x >= 1, x, y free  →  -1 at x = 2
    let sp = lp(
        2,
        2,
        &[(0, 0, -1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)],
        vec![-3.0, 1.0],
        vec![INF, INF],
        vec![-INF, -INF],
        vec![INF, INF],
        vec![0.0, 1.0],
        false,
    );
    assert_opt(&sp, -1.0);
    let s = solve_sparse(&sp, &IpmConfig::default()).unwrap();
    assert!((s.x[0] - 2.0).abs() < 1e-6, "{:?}", s.x);
}

#[test]
fn free_column_far_optimum() {
    // min x : x >= 3e6 (free column; optimum beyond the artificial box)
    let sp = lp(
        1,
        1,
        &[(0, 0, 1.0)],
        vec![3e6],
        vec![INF],
        vec![-INF],
        vec![INF],
        vec![1.0],
        false,
    );
    let s = solve_sparse(&sp, &IpmConfig::default()).unwrap();
    assert_eq!(s.status, IpmStatus::Optimal, "{s:?}");
    assert!((s.objective - 3e6).abs() <= 1e-6 * 3e6);
}

#[test]
fn ranged_row_maximize() {
    // max 3x + y : 1 <= x + y <= 2, 0 <= x <= 1.5, y >= 0  →  5.0
    let sp = lp(
        1,
        2,
        &[(0, 0, 1.0), (0, 1, 1.0)],
        vec![1.0],
        vec![2.0],
        vec![0.0, 0.0],
        vec![1.5, INF],
        vec![3.0, 1.0],
        true,
    );
    assert_opt(&sp, 5.0);
}

#[test]
fn equality_fixed_and_upper_only() {
    // max x + y : x + y + z = 5, z in [1, inf), x <= 3 (x has no lower
    // bound), y in [0, 2], w fixed at 4 with x - w <= 0  →  x + y = 4
    let sp = lp(
        2,
        4,
        &[
            (0, 0, 1.0),
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 0, 1.0),
            (1, 3, -1.0),
        ],
        vec![5.0, -INF],
        vec![5.0, 0.0],
        vec![-INF, 0.0, 1.0, 4.0],
        vec![3.0, 2.0, INF, 4.0],
        vec![1.0, 1.0, 0.0, 0.0],
        true,
    );
    assert_opt(&sp, 4.0);
}

#[test]
fn infeasible_is_certified() {
    // x + y <= 1, x + y >= 3, x, y >= 0
    let sp = lp(
        2,
        2,
        &[(0, 0, 1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)],
        vec![-INF, 3.0],
        vec![1.0, INF],
        vec![0.0, 0.0],
        vec![INF, INF],
        vec![1.0, 1.0],
        false,
    );
    let (a, b) = both(&sp);
    assert_eq!(a.status, IpmStatus::PrimalInfeasible, "{a:?}");
    assert_eq!(b.status, IpmStatus::PrimalInfeasible, "{b:?}");
}

#[test]
fn infeasible_three_rows_is_certified() {
    // x + y + z = 1, x - y >= 2, z - x >= 0, all >= 0  → x >= 2 contradicts row 0
    let sp = lp(
        3,
        3,
        &[
            (0, 0, 1.0),
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 0, 1.0),
            (1, 1, -1.0),
            (2, 2, 1.0),
            (2, 0, -1.0),
        ],
        vec![1.0, 2.0, 0.0],
        vec![1.0, INF, INF],
        vec![0.0; 3],
        vec![INF; 3],
        vec![1.0, 2.0, 3.0],
        false,
    );
    let a = solve_sparse(&sp, &IpmConfig::default()).unwrap();
    assert_eq!(a.status, IpmStatus::PrimalInfeasible, "{a:?}");
}

#[test]
fn unbounded_is_certified() {
    // min -x : x - y <= 1, x, y >= 0
    let sp = lp(
        1,
        2,
        &[(0, 0, 1.0), (0, 1, -1.0)],
        vec![-INF],
        vec![1.0],
        vec![0.0, 0.0],
        vec![INF, INF],
        vec![-1.0, 0.0],
        false,
    );
    let (a, b) = both(&sp);
    assert_eq!(a.status, IpmStatus::DualInfeasible, "{a:?}");
    assert_eq!(b.status, IpmStatus::DualInfeasible, "{b:?}");
}

#[test]
fn unbounded_free_column_is_certified() {
    // min x + y : x - y = 0... no: min x : x + y >= 1, x free, y >= 0
    let sp = lp(
        1,
        2,
        &[(0, 0, 1.0), (0, 1, 1.0)],
        vec![1.0],
        vec![INF],
        vec![-INF, 0.0],
        vec![INF, INF],
        vec![1.0, 0.0],
        false,
    );
    let a = solve_sparse(&sp, &IpmConfig::default()).unwrap();
    assert_eq!(a.status, IpmStatus::DualInfeasible, "{a:?}");
}

#[test]
fn random_lps_match_dual_simplex() {
    let mut rng = XorShift64::new(4242);
    let mut checked = 0;
    let mut unb = (0, 0);
    for case in 0..120 {
        let m = 3 + (rng.next_u64() % 25) as usize;
        let n = m + 1 + (rng.next_u64() % 30) as usize;
        let mut trip = Vec::new();
        for j in 0..n {
            for i in 0..m {
                if rng.next_f64() < 0.3 {
                    trip.push((i, j, (rng.uniform(-3.0, 3.0) * 4.0).round() / 4.0));
                }
            }
        }
        // a feasible point inside the box
        let x0: Vec<f64> = (0..n).map(|_| rng.uniform(0.0, 2.0)).collect();
        let a = CscMatrix::from_triplets(m, n, &trip);
        let mut act = vec![0.0; m];
        a.matvec(&x0, &mut act);
        let (mut rl, mut ru) = (vec![0.0; m], vec![0.0; m]);
        for i in 0..m {
            match rng.next_u64() % 4 {
                0 => (rl[i], ru[i]) = (act[i], act[i]),
                1 => (rl[i], ru[i]) = (-INF, act[i] + rng.uniform(0.0, 1.0)),
                2 => (rl[i], ru[i]) = (act[i] - rng.uniform(0.0, 1.0), INF),
                _ => (rl[i], ru[i]) = (act[i] - 1.0, act[i] + 1.0),
            }
        }
        let (mut cl, mut cu) = (vec![0.0; n], vec![0.0; n]);
        for j in 0..n {
            match rng.next_u64() % 5 {
                0 => (cl[j], cu[j]) = (-INF, INF),
                1 => (cl[j], cu[j]) = (0.0, INF),
                2 => (cl[j], cu[j]) = (-INF, 3.0),
                _ => (cl[j], cu[j]) = (0.0, 3.0),
            }
        }
        let c: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let maximize = case % 3 == 0;
        let sp = SparseLp::new(a, rl, ru, cl, cu, c, maximize);
        let sx = SimplexSolver::new(SimplexConfig::default())
            .solve_sparse(&sp)
            .unwrap();
        let ip = solve_sparse(&sp, &IpmConfig::default()).unwrap();
        let ipp = solve_presolved(&sp, &IpmConfig::default(), &SimplexConfig::default()).unwrap();
        match sx.status {
            SimplexStatus::Optimal => {
                for (tag, s) in [("raw", &ip), ("presolve", &ipp)] {
                    assert_eq!(
                        s.status,
                        IpmStatus::Optimal,
                        "case {case} {tag}: {:?}",
                        s.stats
                    );
                    let tol = 1e-6 * (1.0 + sx.objective.abs());
                    assert!(
                        (s.objective - sx.objective).abs() <= tol,
                        "case {case} {tag}: ipm {} simplex {}",
                        s.objective,
                        sx.objective
                    );
                }
                checked += 1;
            }
            SimplexStatus::Unbounded => {
                assert!(
                    matches!(ip.status, IpmStatus::DualInfeasible | IpmStatus::NoProgress),
                    "case {case}: {:?}",
                    ip.status
                );
                unb.0 += 1;
                if ip.status == IpmStatus::DualInfeasible {
                    unb.1 += 1;
                }
            }
            _ => assert_ne!(ip.status, IpmStatus::Optimal),
        }
    }
    eprintln!(
        "optimal {checked}, unbounded {} (certified {})",
        unb.0, unb.1
    );
    assert!(checked >= 30, "only {checked} optimal cases");
}
