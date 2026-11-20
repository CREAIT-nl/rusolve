//! Cut validity: every generated cut is checked against every integer
//! assignment of small random mixed-integer problems (the continuous part is
//! maximized by an LP with the integers fixed), so no feasible point can be
//! cut off; plus end-to-end agreement between cuts on and off.
#![allow(clippy::needless_range_loop, clippy::type_complexity)]

use super::*;
use crate::bnb::{MipStatus, solve};
use crate::engine::RusolveEngine;
use rusolve_core::CscMatrix;
use rusolve_lp::rng::XorShift64;

fn pick(rng: &mut XorShift64, xs: &[f64]) -> f64 {
    xs[(rng.next_u64() % xs.len() as u64) as usize]
}

/// Random MIP with binaries, a few small general integers and bounded
/// continuous columns; linking rows `x ≤ U·y` / `Σx ≤ K·y` and knapsack rows
/// make implied bounds, covers and c-MIR fire.
pub(super) fn random_mixed(
    rng: &mut XorShift64,
    n_int: usize,
    n_cont: usize,
    m: usize,
) -> MipProblem {
    let n = n_int + n_cont;
    let mut trip = Vec::new();
    let inf = f64::INFINITY;
    let (mut row_lower, mut row_upper) = (Vec::new(), Vec::new());
    let mut col_lower = vec![0.0; n];
    let mut col_upper = vec![1.0; n];
    for j in 0..n_int {
        if rng.next_u64().is_multiple_of(4) {
            col_upper[j] = pick(rng, &[2.0, 3.0]);
        }
        if rng.next_u64().is_multiple_of(8) {
            col_lower[j] = -1.0;
        }
    }
    for j in n_int..n {
        col_upper[j] = pick(rng, &[1.0, 2.5, 4.0, 7.0]);
        if rng.next_u64().is_multiple_of(6) {
            col_lower[j] = 0.5;
        }
    }
    let mut row = 0;
    for _ in 0..m {
        match rng.next_u64() % 4 {
            0 if n_cont > 0 => {
                // linking: Σ_{c∈S} x_c − K·y ≤ 0
                let y = (rng.next_u64() % n_int as u64) as usize;
                let mut k = 0.0;
                for c in n_int..n {
                    if rng.next_u64().is_multiple_of(2) {
                        trip.push((row, c, 1.0));
                        k += col_upper[c];
                    }
                }
                trip.push((row, y, -pick(rng, &[k, (k * 0.6).ceil(), k + 1.0])));
                row_lower.push(-inf);
                row_upper.push(0.0);
            }
            1 => {
                // knapsack over integers (and maybe a continuous)
                let mut tot = 0.0;
                for j in 0..n {
                    if !rng.next_u64().is_multiple_of(3) {
                        let a = pick(rng, &[1.0, 2.0, 3.0, 5.0, 4.5, 7.0]);
                        trip.push((row, j, a));
                        tot += a * col_upper[j];
                    }
                }
                row_lower.push(-inf);
                row_upper
                    .push((tot * pick(rng, &[0.3, 0.5, 0.7])).floor() + pick(rng, &[0.0, 0.5]));
            }
            _ => {
                for j in 0..n {
                    let a = pick(rng, &[0.0, 0.0, 1.0, 2.0, -1.0, 3.0, -2.0, 0.5, 1.5]);
                    if a != 0.0 {
                        trip.push((row, j, a));
                    }
                }
                let b = pick(rng, &[-1.0, 0.0, 1.0, 2.5, 4.0, 6.0]);
                let (lo, hi) = match rng.next_u64() % 4 {
                    0 => (-inf, b),
                    1 => (b, inf),
                    2 => (b, b),
                    _ => (b - 2.0, b + 1.5),
                };
                row_lower.push(lo);
                row_upper.push(hi);
            }
        }
        row += 1;
    }
    MipProblem {
        num_col: n,
        num_row: row,
        matrix: CscMatrix::from_triplets(row, n, &trip),
        row_lower,
        row_upper,
        cost: (0..n)
            .map(|_| pick(rng, &[-5.0, -3.0, -1.0, 0.0, 1.0, 2.0, 1.5, 4.0]))
            .collect(),
        maximize: false,
        col_lower,
        col_upper,
        integer: (0..n).map(|j| j < n_int).collect(),
    }
}

/// Every integer assignment of the box (integers first in column order).
fn int_assignments(p: &MipProblem) -> Vec<Vec<f64>> {
    let ints: Vec<usize> = (0..p.num_col).filter(|&j| p.integer[j]).collect();
    let mut out = Vec::new();
    let mut cur: Vec<f64> = ints.iter().map(|&j| p.col_lower[j]).collect();
    loop {
        out.push(cur.clone());
        let mut k = 0;
        loop {
            if k == ints.len() {
                return out;
            }
            if cur[k] < p.col_upper[ints[k]] {
                cur[k] += 1.0;
                break;
            }
            cur[k] = p.col_lower[ints[k]];
            k += 1;
        }
    }
}

/// Max of `Σ cut` over the continuous part with integers fixed to `a`
/// (`None` when that integer assignment admits no feasible point).
fn max_over_fiber(p: &MipProblem, a: &[f64], inds: &[usize], vals: &[f64]) -> Option<f64> {
    let mut q = p.clone();
    q.cost = vec![0.0; p.num_col];
    for (&j, &v) in inds.iter().zip(vals) {
        q.cost[j] = -v;
    }
    let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
    let mut k = 0;
    for j in 0..p.num_col {
        if p.integer[j] {
            lo[j] = a[k];
            hi[j] = a[k];
            k += 1;
        }
    }
    let mut eng = RusolveEngine::new(&q, MipConfig::default().lp);
    let out = eng.solve(&lo, &hi).expect("fiber LP");
    match out.status {
        LpStatus::Optimal => Some(-out.objective),
        LpStatus::Infeasible => None,
        LpStatus::Unbounded => panic!("bounded continuous part"),
        LpStatus::Cutoff => unreachable!("no cutoff given"),
        LpStatus::IterationLimit => unreachable!("not a probe"),
    }
}

/// Feasible integer assignments of `p` (the fiber LP has a point).
fn feasible_assignments(p: &MipProblem) -> Vec<Vec<f64>> {
    int_assignments(p)
        .into_iter()
        .filter(|a| max_over_fiber(p, a, &[], &[]).is_some())
        .collect()
}

fn run_root(
    p: &MipProblem,
    config: &MipConfig,
) -> Option<(LpOutcome, CutStats, Vec<(CutClass, CutRow)>)> {
    let mut eng = RusolveEngine::new(p, config.lp);
    let lp = eng.solve(&p.col_lower, &p.col_upper).unwrap();
    if lp.status != LpStatus::Optimal {
        return None;
    }
    let mut stats = CutStats::default();
    let mut log = Vec::new();
    let out = root_separation_logged(
        p,
        &mut eng,
        &p.col_lower,
        &p.col_upper,
        lp,
        config,
        &|| false,
        &mut stats,
        Some(&mut log),
        f64::INFINITY,
        &mut |_, _, _| (Vec::new(), 0),
    )
    .unwrap()
    .0;
    Some((out, stats, log))
}

#[test]
fn every_cut_is_valid_on_random_mixed_problems() {
    let mut rng = XorShift64::new(5_20261022);
    let mut config = MipConfig::default();
    config.cuts.modk = true;
    config.cuts.gomory = true;
    let mut total = [0u64; CutClass::ALL.len()];
    for case in 0..400 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let Some((lp, stats, log)) = run_root(&p, &config) else {
            continue;
        };
        let feas = feasible_assignments(&p);
        for (class, cut) in &log {
            total[*class as usize] += 1;
            for a in &feas {
                let mx = max_over_fiber(&p, a, &cut.inds, &cut.vals).unwrap();
                assert!(
                    mx <= cut.rhs + 1e-6 * cut.rhs.abs().max(1.0),
                    "case {case}: {class:?} cut {cut:?} cuts off assignment {a:?} (lhs {mx})"
                );
            }
        }
        // The cut LP bound never exceeds the MIP optimum.
        if lp.status == LpStatus::Optimal && stats.rounds > 0 {
            let best = feas
                .iter()
                .filter_map(|a| {
                    max_over_fiber(
                        &p,
                        a,
                        &(0..p.num_col).collect::<Vec<_>>(),
                        &p.cost.iter().map(|c| -c).collect::<Vec<_>>(),
                    )
                })
                .map(|v| -v)
                .fold(f64::INFINITY, f64::min);
            assert!(
                lp.objective <= best + 1e-6,
                "case {case}: bound {} > opt {best}",
                lp.objective
            );
        }
    }
    eprintln!("cuts per class: {total:?}");
    // Every family must be exercised for the check to mean anything.
    for class in CutClass::ALL {
        assert!(total[class as usize] > 10, "{class:?}: {total:?}");
    }
}

/// Node separation runs on a node LP (tightened bounds) but must only ever
/// produce cuts valid for the *global* problem: every pool cut is checked
/// against every feasible assignment of the root box.
#[test]
fn node_cuts_are_globally_valid() {
    let mut rng = XorShift64::new(9_20261022);
    let mut config = MipConfig::default();
    // One root round leaves work for the node rounds.
    config.cuts.max_rounds = Some(1);
    config.cuts.tree = true;
    config.cuts.modk = true;
    let mut tree_cuts = 0;
    for case in 0..400 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        if lp.status != LpStatus::Optimal {
            continue;
        }
        let mut stats = CutStats::default();
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp,
            &config,
            &|| false,
            &mut stats,
            f64::INFINITY,
            &mut |_, _, _| (Vec::new(), 0),
        )
        .unwrap();
        let Some(mut sep) = sep else { continue };
        if root.status != LpStatus::Optimal {
            continue;
        }
        let frac: Vec<usize> = (0..p.num_col)
            .filter(|&j| p.integer[j] && (root.x[j] - root.x[j].round()).abs() > 1e-6)
            .collect();
        for (j, up) in frac.into_iter().flat_map(|j| [(j, false), (j, true)]) {
            let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
            if up {
                lo[j] = root.x[j].ceil();
            } else {
                hi[j] = root.x[j].floor();
            }
            let node = eng.solve(&lo, &hi).unwrap();
            if node.status != LpStatus::Optimal {
                continue;
            }
            sep.separate_node(
                &p,
                &mut eng,
                &lo,
                &hi,
                node,
                &config,
                f64::INFINITY,
                0,
                &mut stats,
            )
            .unwrap();
        }
        tree_cuts += stats.tree_applied;
        let feas = feasible_assignments(&p);
        for cut in sep.pool.cuts.iter().filter(|c| !c.deleted) {
            for a in &feas {
                let mx = max_over_fiber(&p, a, &cut.inds, &cut.vals).unwrap();
                assert!(
                    mx <= cut.rhs + 1e-6 * cut.rhs.abs().max(1.0),
                    "case {case}: {:?} pool cut {:?} ≤ {} cuts off {a:?} (lhs {mx})",
                    cut.class,
                    cut.vals,
                    cut.rhs
                );
            }
        }
    }
    eprintln!("tree cuts applied: {tree_cuts}");
    assert!(tree_cuts > 50, "node separation never fired: {tree_cuts}");
}

/// Node-LP aging: every optimal node LP ticks the row ages
/// ([`Separator::age_node`], HiGHS HighsSearch.cpp:958 + performAging), so a
/// cut slack at the node point leaves the LP once its age passes the
/// (periodically clamped) `lp_age_limit`, while a cut tight at every point
/// it is aged against never ages and stays.
#[test]
fn node_lp_aging_drops_slack_cuts_and_keeps_tight_ones() {
    let mut rng = XorShift64::new(12_20261023);
    let mut config = MipConfig::default();
    config.cuts.tree = true;
    config.cuts.max_rounds = Some(1);
    config.cuts.modk = true;
    config.cuts.gomory = false;
    let feastol = config.feastol;
    let mut exercised = 0;
    for _case in 0..400 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp = match eng.solve(&p.col_lower, &p.col_upper) {
            Ok(lp) if lp.status == LpStatus::Optimal => lp,
            _ => continue,
        };
        let mut stats = CutStats::default();
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp,
            &config,
            &|| false,
            &mut stats,
            f64::INFINITY,
            &mut |_, _, _| (Vec::new(), 0),
        )
        .unwrap();
        let Some(mut sep) = sep else { continue };
        if root.status != LpStatus::Optimal || sep.lp_cuts.is_empty() {
            continue;
        }
        // A point slack for every LP cut: `drop_slack_at_end` kept exactly
        // the rows binding at `root.x`, so aging against `root.x` must keep
        // all of them; aging 20 epochs against the zero point must drop
        // exactly the rows slack there (tracked by pool index, since drops
        // shift positions).
        let zero = vec![0.0; p.num_col];
        let is_slack = |sep: &Separator, c: usize| -> bool {
            let cut = &sep.pool.cuts[c];
            let act: f64 = cut
                .inds
                .iter()
                .zip(&cut.vals)
                .map(|(&j, &v)| v * zero[j])
                .sum();
            act < cut.rhs - feastol
        };
        let slack: Vec<usize> = sep
            .lp_cuts
            .iter()
            .map(|e| e.pool.unwrap())
            .filter(|&c| is_slack(&sep, c))
            .collect();
        let tight: Vec<usize> = sep
            .lp_cuts
            .iter()
            .map(|e| e.pool.unwrap())
            .filter(|&c| !is_slack(&sep, c))
            .collect();
        if slack.is_empty() {
            continue;
        }
        exercised += 1;
        for _ in 0..20 {
            sep.age_node(&mut eng, &root.x, feastol, false, &mut CutStats::default());
        }
        assert_eq!(
            sep.lp_cuts.len(),
            slack.len() + tight.len(),
            "tight rows aged out"
        );
        for _ in 0..20 {
            sep.age_node(&mut eng, &zero, feastol, false, &mut CutStats::default());
        }
        let left: Vec<usize> = sep.lp_cuts.iter().map(|e| e.pool.unwrap()).collect();
        assert!(
            left.iter().all(|c| !slack.contains(c)),
            "rows slack for 20 node LPs survived aging: {left:?} vs slack {slack:?}"
        );
        assert!(
            tight.iter().all(|c| left.contains(c)),
            "rows tight at every aged point were dropped: {left:?} vs tight {tight:?}"
        );
    }
    assert!(exercised > 10, "aging never exercised: {exercised}");
}

/// The effort-capped pass-through in [`Separator::separate_node`] must not
/// age: the caller's [`Separator::age_node`] tick already counted that node
/// LP, and aging again there double-counted the epoch (each node LP moved
/// rows two epochs toward deletion).
#[test]
fn separate_node_without_budget_ages_nothing() {
    let mut rng = XorShift64::new(13_20261023);
    let mut config = MipConfig::default();
    config.cuts.tree = true;
    config.cuts.max_rounds = Some(1);
    config.cuts.modk = true;
    config.cuts.gomory = false;
    for _case in 0..400 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp = match eng.solve(&p.col_lower, &p.col_upper) {
            Ok(lp) if lp.status == LpStatus::Optimal => lp,
            _ => continue,
        };
        let mut stats = CutStats::default();
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp,
            &config,
            &|| false,
            &mut stats,
            f64::INFINITY,
            &mut |_, _, _| (Vec::new(), 0),
        )
        .unwrap();
        let Some(mut sep) = sep else { continue };
        if root.status != LpStatus::Optimal {
            continue;
        }
        // A negative effort with 20,000 total iterations puts the cap at
        // tree_effort·total + 10,000 < 0: every node passes through.
        config.cuts.tree_effort = -1.0;
        let node = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        let epochs = sep.epochs;
        let out = sep
            .separate_node(
                &p,
                &mut eng,
                &p.col_lower,
                &p.col_upper,
                node,
                &config,
                f64::INFINITY,
                20_000,
                &mut stats,
            )
            .unwrap();
        assert_eq!(sep.epochs, epochs, "capped node aged the pool");
        assert_eq!(stats.tree_calls, 0);
        assert_eq!(out.status, LpStatus::Optimal);
        return;
    }
    panic!("no case with an optimal root LP found");
}

#[test]
fn cuts_on_and_off_agree_on_random_mixed_problems() {
    let mut rng = XorShift64::new(77);
    let on = MipConfig::default();
    let mut off = MipConfig::default();
    off.cuts.enabled = false;
    for case in 0..150 {
        let n_int = 2 + (rng.next_u64() % 5) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 5) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let a = solve(&p, &on).unwrap();
        let b = solve(&p, &off).unwrap();
        assert_eq!(a.status, b.status, "case {case}");
        if a.status == MipStatus::Optimal {
            assert!(
                (a.objective - b.objective).abs() <= 1e-6 * b.objective.abs().max(1.0),
                "case {case}: {} vs {}",
                a.objective,
                b.objective
            );
            let v = p.violation(&a.x);
            assert!(v.row <= 1e-6 && v.bound <= 1e-6 && v.integrality <= 1e-6);
        }
    }
}

/// Weak facility location: `Σ_i x_ij ≤ n·y_j` implies `x_ij ≤ y_j`, which
/// closes the root gap.
#[test]
fn implied_bounds_close_weak_facility_location_gap() {
    // 2 facilities (y0,y1 binary, open cost 3), 3 customers, x_fc ∈ [0,1]
    // continuous, cost 1 to facility c mod 2 and 5 otherwise; demand
    // Σ_f x_fc = 1, linking Σ_c x_fc ≤ 3 y_f. Weak LP: y = (2/3, 1/3) ⇒ 6;
    // optimum opens both ⇒ 9, which the x ≤ y relaxation already attains.
    let (nf, nc) = (2, 3);
    let n = nf + nf * nc;
    let mut trip = Vec::new();
    let mut row_lower = Vec::new();
    let mut row_upper = Vec::new();
    let x = |f: usize, c: usize| nf + f * nc + c;
    for c in 0..nc {
        for f in 0..nf {
            trip.push((c, x(f, c), 1.0));
        }
        row_lower.push(1.0);
        row_upper.push(1.0);
    }
    for f in 0..nf {
        for c in 0..nc {
            trip.push((nc + f, x(f, c), 1.0));
        }
        trip.push((nc + f, f, -(nc as f64)));
        row_lower.push(f64::NEG_INFINITY);
        row_upper.push(0.0);
    }
    let m = nc + nf;
    let mut cost = vec![0.0; n];
    cost[0] = 3.0;
    cost[1] = 3.0;
    for f in 0..nf {
        for c in 0..nc {
            cost[x(f, c)] = if c % 2 == f { 1.0 } else { 5.0 };
        }
    }
    let p = MipProblem {
        num_col: n,
        num_row: m,
        matrix: CscMatrix::from_triplets(m, n, &trip),
        row_lower,
        row_upper,
        cost,
        maximize: false,
        col_lower: vec![0.0; n],
        col_upper: vec![1.0; n],
        integer: (0..n).map(|j| j < nf).collect(),
    };
    let mut config = MipConfig::default();
    config.cuts.aggregation = false;
    config.cuts.gomory = false;
    let (lp, stats, log) = run_root(&p, &config).unwrap();
    assert!((stats.root_bound_before - 6.0).abs() < 1e-6, "{stats:?}");
    assert!(
        (lp.objective - 9.0).abs() < 1e-6,
        "{} {stats:?}",
        lp.objective
    );
    assert!(log.iter().all(|(c, _)| *c == CutClass::ImpliedBound));
    assert!(!log.is_empty());
}

/// Single binary knapsack `5x0 + 5x1 + 5x2 + 5x3 + 3x4 ≤ 11`, maximize
/// `Σ x + 0.5·x4`: the LP takes 2.2 of the 5s plus x4; the cover
/// `x0+x1+x2+x3 ≤ 2` lifted on x4 cuts it off (Gu, Nemhauser &
/// Savelsbergh 1998). The row's integral slack is a general integer in the
/// base row, so (as in HiGHS) the mixed-integer cover lifting produces it.
#[test]
fn lifted_cover_on_single_knapsack_row() {
    let n = 5;
    let w = [5.0, 5.0, 5.0, 5.0, 3.0];
    let trip: Vec<(usize, usize, f64)> = (0..n).map(|j| (0, j, w[j])).collect();
    let p = MipProblem {
        num_col: n,
        num_row: 1,
        matrix: CscMatrix::from_triplets(1, n, &trip),
        row_lower: vec![f64::NEG_INFINITY],
        row_upper: vec![11.0],
        cost: vec![-1.0, -1.0, -1.0, -1.0, -1.5],
        maximize: false,
        col_lower: vec![0.0; n],
        col_upper: vec![1.0; n],
        integer: vec![true; n],
    };
    let (lp, stats, log) = run_root(&p, &MipConfig::default()).unwrap();
    assert!(
        log.iter().any(|(c, cut)| {
            matches!(c, CutClass::KnapsackCover | CutClass::MixedIntegerCover)
                && cut.vals == [1.0; 5]
                && cut.rhs == 2.0
        }),
        "{log:?}"
    );
    // Any three items exceed 11, so the optimum is one 5 plus x4 = −2.5;
    // the lifted cover closes the whole gap.
    assert!(stats.root_bound_before < -2.6, "{stats:?}");
    assert!((lp.objective + 2.5).abs() < 1e-6, "{}", lp.objective);
}

/// Edge formulation of a triangle (`x_i + x_j ≤ 1`), maximize `Σ x`: the LP
/// sits at ½ everywhere (1.5); the clique `x0 + x1 + x2 ≤ 1` assembled from
/// the three edges closes the gap (Atamtürk, Nemhauser & Savelsbergh 2000).
#[test]
fn clique_cut_closes_triangle() {
    let edges = [(0, 1), (1, 2), (0, 2)];
    let mut trip = Vec::new();
    for (r, &(a, b)) in edges.iter().enumerate() {
        trip.push((r, a, 1.0));
        trip.push((r, b, 1.0));
    }
    let p = MipProblem {
        num_col: 3,
        num_row: 3,
        matrix: CscMatrix::from_triplets(3, 3, &trip),
        row_lower: vec![f64::NEG_INFINITY; 3],
        row_upper: vec![1.0; 3],
        cost: vec![-1.0; 3],
        maximize: false,
        col_lower: vec![0.0; 3],
        col_upper: vec![1.0; 3],
        integer: vec![true; 3],
    };
    let (lp, stats, log) = run_root(&p, &MipConfig::default()).unwrap();
    assert!((stats.root_bound_before + 1.5).abs() < 1e-6, "{stats:?}");
    assert!((lp.objective + 1.0).abs() < 1e-6, "{}", lp.objective);
    assert!(
        log.iter()
            .any(|(c, cut)| *c == CutClass::Clique && cut.inds.len() == 3 && cut.rhs == 1.0)
    );
}

/// The weak-round test is scale-free: on `gap` the objective (≈ 810) dwarfs
/// the root gap (≈ 40), and rounds adding 1–3% of the gain must keep going.
#[test]
fn weak_round_is_relative_to_gain_not_objective() {
    // gap_s1 rounds 3 and 4: +0.199 and +0.114 on a gain of ≈ 7.3.
    assert!(!weak_round(804.99, 812.08, 812.27, 1e-3));
    assert!(!weak_round(804.99, 812.27, 812.39, 1e-3));
    // A flat round is weak, also before any gain (bin_packing).
    assert!(weak_round(10.24, 10.24, 10.24, 1e-3));
    assert!(weak_round(10.24, 10.24, 10.24 + 1e-12, 1e-3));
    // A tiny step on a large gain is weak.
    assert!(weak_round(0.0, 100.0, 100.05, 1e-3));
}

/// 5-cycle edge formulation (`x_i + x_{i+1} ≤ 1`), maximize `Σ x`: the LP
/// sits at ½ everywhere (2.5) and no clique has more than two nodes; the
/// sum of the five edges halved and rounded, `Σ x ≤ 2`, is the zero-half
/// (mod-2) cut that closes the gap (Caprara & Fischetti 1996).
#[test]
fn modk_closes_odd_cycle() {
    let mut trip = Vec::new();
    for i in 0..5 {
        trip.push((i, i, 1.0));
        trip.push((i, (i + 1) % 5, 1.0));
    }
    let p = MipProblem {
        num_col: 5,
        num_row: 5,
        matrix: CscMatrix::from_triplets(5, 5, &trip),
        row_lower: vec![f64::NEG_INFINITY; 5],
        row_upper: vec![1.0; 5],
        cost: vec![-1.0; 5],
        maximize: false,
        col_lower: vec![0.0; 5],
        col_upper: vec![1.0; 5],
        integer: vec![true; 5],
    };
    let mut on = MipConfig::default();
    on.cuts.modk = true;
    let (lp, stats, log) = run_root(&p, &on).unwrap();
    assert!((stats.root_bound_before + 2.5).abs() < 1e-6, "{stats:?}");
    assert!((lp.objective + 2.0).abs() < 1e-6, "{}", lp.objective);
    assert!(log.iter().any(|(c, _)| *c == CutClass::ModK), "{log:?}");
    let (lp, _, _) = run_root(&p, &MipConfig::default()).unwrap();
    assert!(lp.objective < -2.4, "{}", lp.objective);
}

/// Pure-integer rows with coefficients that are not multiples of 2, 3, 5 or
/// 7, so every modulus gets a chance; mod-k is the only separator and every
/// cut is checked against every feasible point.
#[test]
fn modk_cuts_are_valid_on_random_integer_problems() {
    let mut rng = XorShift64::new(11_20261022);
    let mut config = MipConfig::default();
    config.cuts.implied_bounds = false;
    config.cuts.cliques = false;
    config.cuts.aggregation = false;
    config.cuts.gomory = false;
    config.cuts.modk = true;
    let mut total = 0;
    for case in 0..400 {
        let n = 3 + (rng.next_u64() % 4) as usize;
        let m = 2 + (rng.next_u64() % 4) as usize;
        let mut trip = Vec::new();
        let (mut lo, mut hi) = (Vec::new(), Vec::new());
        for i in 0..m {
            let mut tot = 0.0;
            for j in 0..n {
                let a = pick(&mut rng, &[0.0, 1.0, 1.0, 2.0, 3.0, 5.0, -1.0, -3.0]);
                if a != 0.0 {
                    trip.push((i, j, a));
                    tot += a.max(0.0) * 2.0;
                }
            }
            let b = (tot * pick(&mut rng, &[0.3, 0.5, 0.7])).floor() + pick(&mut rng, &[0.0, 1.0]);
            match rng.next_u64() % 3 {
                0 => {
                    lo.push(f64::NEG_INFINITY);
                    hi.push(b);
                }
                1 => {
                    lo.push(b - 3.0);
                    hi.push(f64::INFINITY);
                }
                _ => {
                    lo.push(b - 2.0);
                    hi.push(b);
                }
            }
        }
        let p = MipProblem {
            num_col: n,
            num_row: m,
            matrix: CscMatrix::from_triplets(m, n, &trip),
            row_lower: lo,
            row_upper: hi,
            cost: (0..n)
                .map(|_| pick(&mut rng, &[-5.0, -3.0, -2.0, -1.0, 1.0]))
                .collect(),
            maximize: false,
            col_lower: vec![0.0; n],
            col_upper: (0..n)
                .map(|_| pick(&mut rng, &[1.0, 1.0, 2.0, 3.0]))
                .collect(),
            integer: vec![true; n],
        };
        let Some((_, _, log)) = run_root(&p, &config) else {
            continue;
        };
        let feas = feasible_assignments(&p);
        for (class, cut) in &log {
            assert_eq!(*class, CutClass::ModK);
            total += 1;
            for a in &feas {
                let lhs: f64 = cut
                    .inds
                    .iter()
                    .zip(&cut.vals)
                    .map(|(&j, &v)| v * a[j])
                    .sum();
                assert!(
                    lhs <= cut.rhs + 1e-6 * cut.rhs.abs().max(1.0),
                    "case {case}: mod-k cut {cut:?} cuts off {a:?} (lhs {lhs})"
                );
            }
        }
    }
    assert!(total > 50, "only {total} mod-k cuts");
}

fn gomory_only() -> MipConfig {
    let mut config = MipConfig::default();
    config.cuts.implied_bounds = false;
    config.cuts.cliques = false;
    config.cuts.aggregation = false;
    config.cuts.modk = false;
    config.cuts.gomory = true;
    config
}

/// Wolsey's textbook integer program, `max 5x + 4y` over `6x + 4y <= 24`,
/// `x + 2y <= 6`: the LP optimum is (3, 1.5) and the integer optimum 20. The
/// tableau rows give Gomory cuts that cut off the LP point and no integer
/// point, and a single round already lifts the bound (Gomory 1960; Wolsey,
/// Integer Programming, 1998, ch. 8).
#[test]
fn gomory_cuts_off_the_fractional_lp_point() {
    let trip = [(0, 0, 6.0), (0, 1, 4.0), (1, 0, 1.0), (1, 1, 2.0)];
    let p = MipProblem {
        num_col: 2,
        num_row: 2,
        matrix: CscMatrix::from_triplets(2, 2, &trip),
        row_lower: vec![f64::NEG_INFINITY; 2],
        row_upper: vec![24.0, 6.0],
        cost: vec![-5.0, -4.0],
        maximize: false,
        col_lower: vec![0.0; 2],
        col_upper: vec![10.0; 2],
        integer: vec![true; 2],
    };
    let mut config = gomory_only();
    config.cuts.max_rounds = Some(1);
    let mut eng = RusolveEngine::new(&p, config.lp);
    let lp0 = eng.solve(&p.col_lower, &p.col_upper).unwrap();
    assert!((lp0.x[0] - 3.0).abs() < 1e-6 && (lp0.x[1] - 1.5).abs() < 1e-6);
    let (lp, stats, log) = run_root(&p, &config).unwrap();
    assert!(stats.generated[CutClass::Gomory as usize] > 0, "{stats:?}");
    assert!(
        log.iter().all(|(c, _)| *c == CutClass::Gomory),
        "only Gomory is on: {log:?}"
    );
    let feas = feasible_assignments(&p);
    for (_, cut) in &log {
        let lhs = |x: &[f64]| -> f64 {
            cut.inds
                .iter()
                .zip(&cut.vals)
                .map(|(&j, &v)| v * x[j])
                .sum()
        };
        assert!(
            lhs(&lp0.x) > cut.rhs + 1e-6,
            "cut {cut:?} does not cut off the LP point"
        );
        for a in &feas {
            assert!(lhs(a) <= cut.rhs + 1e-6, "cut {cut:?} cuts off {a:?}");
        }
    }
    // -22.5 before, strictly better after one round, never past -20.
    assert!(lp.objective > -22.5 + 1e-6, "{}", lp.objective);
    assert!(lp.objective <= -20.0 + 1e-6, "{}", lp.objective);
}

/// Round-level Gomory gate on a set-partitioning system (`Ax = 1`, `x`
/// binary — the set_partition perf family, Balas & Padberg 1976, "Set
/// partitioning: a survey", SIAM Rev. 18): round 1 generates Gomory cuts
/// but the pool selects none of them, so the separator is parked for the
/// rest of the root loop (`TableauState::parked`) instead of re-paying the
/// basis scan in every later round — while the stored cuts stay available
/// to the pool, which still selects some in round 2. Measured motivation:
/// set_partition_s1 runs 19 root rounds with a cut loop going 0.004 s →
/// 0.015 s from the scan alone.
#[test]
fn gomory_parks_after_a_round_without_selected_cuts() {
    let mut rng = XorShift64::new(42_2026_1023);
    let m = 12;
    // planted exact cover: disjoint random subsets over a permutation
    let mut perm: Vec<usize> = (0..m).collect();
    for k in (1..m).rev() {
        let j = (rng.next_u64() % (k as u64 + 1)) as usize;
        perm.swap(k, j);
    }
    let mut cols: Vec<Vec<usize>> = Vec::new();
    let mut k = 0;
    while k < m {
        let size = 1 + (rng.next_u64() % 4) as usize;
        cols.push(perm[k..(k + size).min(m)].to_vec());
        k += size;
    }
    while cols.len() < 48 {
        let size = 1 + (rng.next_u64() % 5) as usize;
        let mut idx: Vec<usize> = (0..m).collect();
        for s in (1..m).rev() {
            let j = (rng.next_u64() % (s as u64 + 1)) as usize;
            idx.swap(s, j);
        }
        cols.push(idx[..size].to_vec());
    }
    let n = cols.len();
    let mut trip = Vec::new();
    for (j, c) in cols.iter().enumerate() {
        for &i in c {
            trip.push((i, j, 1.0));
        }
    }
    let p = MipProblem {
        num_col: n,
        num_row: m,
        matrix: CscMatrix::from_triplets(m, n, &trip),
        row_lower: vec![1.0; m],
        row_upper: vec![1.0; m],
        cost: (0..n)
            .map(|j| {
                (1.0 + (rng.next_u64() % 10) as f64) * cols[j].len() as f64
                    + (rng.next_u64() % 6) as f64
            })
            .collect(),
        maximize: false,
        col_lower: vec![0.0; n],
        col_upper: vec![1.0; n],
        integer: vec![true; n],
    };
    // The default density (0.6·avg row length) rejects every cut this
    // instance's long set-partitioning rows produce; the test targets the
    // parking gate, so it pins the looser density where the separator
    // generates.
    let mut cfg = MipConfig::default();
    cfg.cuts.gomory_park = true; // default on (fcn[fixed] needs it)
    cfg.cuts.gomory_density = 2.0;
    // Extended clique cuts close this set-partitioning root in two
    // rounds before any Gomory cut is selected; the bare cliques leave room
    // for the parking gate to act.
    cfg.cuts.clique_extend = false;
    let (_, stats, _) = run_root(&p, &cfg).unwrap();
    assert!(
        stats.rounds >= 2,
        "cut loop too short to test parking: {stats:?}"
    );
    assert!(
        stats.generated[CutClass::Gomory as usize] > 0,
        "instance should make the separator run: {stats:?}"
    );
    assert!(
        stats.gomory_calls < stats.rounds as u64,
        "gate never parked the separator: {stats:?}"
    );
    assert!(
        stats.applied[CutClass::Gomory as usize] > 0,
        "stored cuts must stay selectable after parking: {stats:?}"
    );
}

/// Gomory is the only separator on random mixed problems (continuous
/// columns, ranged and equality rows, negative lower bounds, general
/// integers): every cut keeps every feasible integer assignment, and each
/// first-round cut is violated at the root LP point.
#[test]
fn gomory_cuts_are_valid_on_random_mixed_problems() {
    let mut rng = XorShift64::new(13_20261022);
    let mut config = gomory_only();
    config.cuts.max_rounds = Some(1);
    let mut total = 0;
    for case in 0..500 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp0 = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        if lp0.status != LpStatus::Optimal {
            continue;
        }
        let Some((_, _, log)) = run_root(&p, &config) else {
            continue;
        };
        let feas = feasible_assignments(&p);
        for (class, cut) in &log {
            assert_eq!(*class, CutClass::Gomory);
            total += 1;
            let at_lp: f64 = cut
                .inds
                .iter()
                .zip(&cut.vals)
                .map(|(&j, &v)| v * lp0.x[j])
                .sum();
            assert!(
                at_lp > cut.rhs,
                "case {case}: Gomory cut {cut:?} does not cut off the LP point (lhs {at_lp})"
            );
            for a in &feas {
                let mx = max_over_fiber(&p, a, &cut.inds, &cut.vals).unwrap();
                assert!(
                    mx <= cut.rhs + 1e-6 * cut.rhs.abs().max(1.0),
                    "case {case}: Gomory cut {cut:?} cuts off {a:?} (lhs {mx})"
                );
            }
        }
    }
    assert!(total > 50, "only {total} Gomory cuts");
}

/// Obsolete-row removal (port of HiGHS `removeObsoleteRows`): dropping
/// the LP rows whose slack is basic (or, without basis status, whose dual is
/// zero) leaves the LP optimum unchanged (their duals contribute nothing),
/// the dropped cuts re-compete from the pool, and the bookkeeping
/// (`lp_cuts`, relax rows, engine rows) stays aligned.
#[test]
fn drop_obsolete_rows_keeps_the_optimum_and_releases_the_cuts() {
    let mut rng = XorShift64::new(13_20261024);
    let mut config = MipConfig::default();
    // The loop must LEAVE obsolete rows in the LP for the drop to act on.
    config.cuts.efficacy_end = false;
    let feastol = config.feastol;
    let mut exercised = 0;
    for _case in 0..1500 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp = match eng.solve(&p.col_lower, &p.col_upper) {
            Ok(lp) if lp.status == LpStatus::Optimal => lp,
            _ => continue,
        };
        let mut stats = CutStats::default();
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp,
            &config,
            &|| false,
            &mut stats,
            f64::INFINITY,
            &mut |_, _, _| (Vec::new(), 0),
        )
        .unwrap();
        let Some(mut sep) = sep else { continue };
        if root.status != LpStatus::Optimal || sep.lp_cuts.is_empty() {
            continue;
        }
        // The final drop_slack_at_end can leave `root`'s basis stale for
        // `obsolete_row_mask`: re-solve the same row set and mask off the
        // fresh solve.
        let root = match eng.solve(&p.col_lower, &p.col_upper) {
            Ok(r) if r.status == LpStatus::Optimal => r,
            _ => continue,
        };
        let Some(mask) = obsolete_row_mask(&mut eng, &root, sep.lp_cuts.len(), feastol) else {
            continue;
        };
        let doomed: Vec<usize> = sep
            .lp_cuts
            .iter()
            .enumerate()
            .filter(|(k, _)| mask[*k])
            .map(|(_, e)| e.pool.unwrap())
            .collect();
        if doomed.is_empty() || doomed.len() == sep.lp_cuts.len() {
            continue; // want both classes present: the invariant is the mix
        }
        exercised += 1;
        let obj0 = root.objective;
        let total = sep.lp_cuts.len();
        assert!(drop_obsolete_rows(&mut sep, &mut eng, &mask));
        assert_eq!(sep.lp_cuts.len(), total - doomed.len());
        for &c in &doomed {
            assert!(!sep.pool.cuts[c].deleted, "released cut deleted from pool");
            assert!(sep.pool.cuts[c].age >= 0, "released cut still marked in LP");
        }
        // Removing obsolete rows cannot move the optimum: solve the shrunken
        // LP cold and compare.
        let re = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        assert_eq!(re.status, LpStatus::Optimal);
        assert!(
            (re.objective - obj0).abs() <= 1e-6 * obj0.abs().max(1.0),
            "dropping obsolete rows moved the LP bound {obj0} -> {}",
            re.objective
        );
    }
    assert!(exercised >= 10, "only {exercised} exercised cases");
}

/// With `efficacy_rounds` + `efficacy_end` the root loop ends with no
/// basic-slack row in the LP (HiGHS `removeObsoleteRows` after every round,
/// HighsMipSolverData.cpp:2098, and before branching, :2481), and the row
/// bookkeeping never changes the LP's optimum.
#[test]
fn efficacy_loop_leaves_only_nonbasic_rows_in_the_lp() {
    let mut rng = XorShift64::new(14_20261024);
    let mut config = MipConfig::default();
    config.cuts.efficacy_rounds = true;
    config.cuts.efficacy_end = true;
    config.cuts.drop_slack_at_end = false;
    let mut exercised = 0;
    for _case in 0..400 {
        let n_int = 2 + (rng.next_u64() % 4) as usize;
        let n_cont = (rng.next_u64() % 4) as usize;
        let m = 1 + (rng.next_u64() % 4) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp = match eng.solve(&p.col_lower, &p.col_upper) {
            Ok(lp) if lp.status == LpStatus::Optimal => lp,
            _ => continue,
        };
        let mut stats = CutStats::default();
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp,
            &config,
            &|| false,
            &mut stats,
            f64::INFINITY,
            &mut |_, _, _| (Vec::new(), 0),
        )
        .unwrap();
        let Some(sep) = sep else { continue };
        if root.status != LpStatus::Optimal || sep.lp_cuts.is_empty() {
            continue;
        }
        exercised += 1;
        // The last drop may have left `root`'s basis stale for `duals`
        // (rows changed after the solve): re-solve the final row set —
        // the same LP, so the objective must be unchanged — and read the
        // duals off the fresh solve.
        let re = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        assert_eq!(re.status, LpStatus::Optimal);
        assert!(
            (re.objective - root.objective).abs() <= 1e-6 * root.objective.abs().max(1.0),
            "row bookkeeping changed the LP: {} -> {}",
            root.objective,
            re.objective
        );
        let basic = eng
            .cut_rows_basic(&re)
            .expect("RusolveEngine provides basis status");
        let n_lp = sep.lp_cuts.len();
        assert_eq!(basic.len(), n_lp);
        for (k, &b) in basic.iter().enumerate() {
            assert!(!b, "obsolete row (basic slack) left in the LP at {k}");
        }
    }
    assert!(exercised >= 10, "only {exercised} exercised cases");
}

/// All fiber optima of `p` (integers fixed per feasible assignment, the
/// continuous part minimized by an LP), sorted by internal objective.
fn feasible_points(p: &MipProblem) -> Vec<(Vec<f64>, f64)> {
    let mut out = Vec::new();
    for a in feasible_assignments(p) {
        let (mut lo, mut hi) = (p.col_lower.clone(), p.col_upper.clone());
        let mut k = 0;
        for j in 0..p.num_col {
            if p.integer[j] {
                lo[j] = a[k];
                hi[j] = a[k];
                k += 1;
            }
        }
        let mut eng = RusolveEngine::new(p, MipConfig::default().lp);
        let lp = eng.solve(&lo, &hi).unwrap();
        if lp.status == LpStatus::Optimal {
            out.push((lp.x.clone(), lp.objective));
        }
    }
    out.sort_by(|a, b| a.1.total_cmp(&b.1));
    out
}

/// The cutoff row of the entry incumbent: it lands before round 1, and is
/// stripped before branching (inert in the tree, pure matrix tax),
/// while its rhs survives on the separator as the guard against worse
/// re-lands. Valid cuts alone keep the bound between the root LP and the
/// incumbent.
#[test]
fn cutoff_row_from_entry_incumbent_lands_then_strips_before_tree() {
    let mut rng = XorShift64::new(16_202_610_241);
    let config = MipConfig::default();
    assert!(config.cuts.root_cutoff_row, "opt-in default is on");
    let mut cases = 0;
    for _ in 0..250 {
        let n_int = 2 + (rng.next_u64() % 3) as usize;
        let n_cont = (rng.next_u64() % 3) as usize;
        let m = 1 + (rng.next_u64() % 3) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        if p.cost.iter().all(|c| *c == 0.0) {
            continue; // a cutoff row needs an objective to cut with
        }
        let Some((_, best_obj)) = feasible_points(&p).first().cloned() else {
            continue;
        };
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp0 = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        if lp0.status != LpStatus::Optimal {
            continue;
        }
        let mut stats = CutStats::default();
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp0.clone(),
            &config,
            &|| false,
            &mut stats,
            best_obj,
            &mut |_, _, _| (Vec::new(), 0),
        )
        .unwrap();
        let Some(sep) = sep else { continue };
        cases += 1;
        assert_eq!(stats.cutoff_rows, 1, "row lands from the entry cutoff");
        assert!(
            !sep.cutoff_row && sep.cutoff_ub == best_obj,
            "the row is stripped before branching; its rhs survives as the re-land guard"
        );
        assert!(
            sep.lp_cuts.iter().all(|e| e.pool.is_some()),
            "no pool-less row reaches the tree"
        );
        if root.status == LpStatus::Optimal {
            assert!(
                root.objective <= best_obj + 1e-6 * best_obj.abs().max(1.0),
                "valid cuts keep the stripped bound below the incumbent: {} > {best_obj}",
                root.objective
            );
            assert!(
                root.objective >= lp0.objective - 1e-6 * lp0.objective.abs().max(1.0),
                "a valid cut cannot lower the bound"
            );
        }
    }
    assert!(cases > 20, "only {cases} usable cases");
}

/// A verified hook point lands the row mid-loop (cold re-solve), a better
/// point tightens it (the hook sees the tightened cutoff on its next call),
/// and a non-improving point is refused; every landed row is stripped
/// before branching.
#[test]
fn cutoff_row_lands_from_hook_and_tightens_on_improvement() {
    let mut rng = XorShift64::new(16_202_610_242);
    let mut config = MipConfig::default();
    config.cuts.root_cutoff_tighten = true;
    let mut exercised = 0;
    for _ in 0..250 {
        let n_int = 2 + (rng.next_u64() % 3) as usize;
        let n_cont = (rng.next_u64() % 3) as usize;
        let m = 1 + (rng.next_u64() % 3) as usize;
        let p = random_mixed(&mut rng, n_int, n_cont, m);
        if p.cost.iter().all(|c| *c == 0.0) {
            continue;
        }
        let points = feasible_points(&p);
        let Some((_, best_obj)) = points.first().cloned() else {
            continue;
        };
        // A strictly worse feasible point to land the first row.
        let Some(worse) = points
            .iter()
            .find(|(_, o)| *o > best_obj + 1e-6 * best_obj.abs().max(1.0))
        else {
            continue;
        };
        let (worse_x, worse_obj) = (worse.0.clone(), worse.1);
        let best_x = points[0].0.clone();
        let mut eng = RusolveEngine::new(&p, config.lp);
        let lp0 = eng.solve(&p.col_lower, &p.col_upper).unwrap();
        if lp0.status != LpStatus::Optimal {
            continue;
        }
        let mut stats = CutStats::default();
        let mut seen_cutoffs = Vec::new();
        let mut calls = 0u32;
        let (root, sep) = root_separation(
            &p,
            &mut eng,
            &p.col_lower,
            &p.col_upper,
            lp0,
            &config,
            &|| false,
            &mut stats,
            f64::INFINITY,
            &mut |_, _, cutoff| {
                seen_cutoffs.push(cutoff);
                calls += 1;
                let x = if calls == 1 { &worse_x } else { &best_x };
                (
                    vec![(x.clone(), crate::heuristics::Source::RandomizedRounding)],
                    0,
                )
            },
        )
        .unwrap();
        let Some(sep) = sep else { continue };
        if stats.rounds < 2 {
            continue; // needs a second hook call to see the tightening
        }
        exercised += 1;
        assert_eq!(seen_cutoffs[0], f64::INFINITY, "first call is blind");
        assert!(
            stats.cutoff_rows >= 2,
            "the hook point lands the row and the better point tightens it"
        );
        assert!(
            seen_cutoffs[1..].iter().all(|c| *c <= worse_obj + 1e-9),
            "after landing, the hook runs under the incumbent's cutoff"
        );
        assert!(
            sep.lp_cuts.iter().all(|e| e.pool.is_some()),
            "every landed row is stripped before branching"
        );
        assert!(
            !sep.cutoff_row && sep.cutoff_ub <= worse_obj + 1e-9,
            "the last landed rhs survives as the re-land guard"
        );
        if root.status == LpStatus::Optimal {
            assert!(
                root.objective <= best_obj + 1e-6 * best_obj.abs().max(1.0),
                "valid cuts keep the stripped bound below the optimum"
            );
        }
        if exercised >= 10 {
            break;
        }
    }
    assert!(exercised >= 5, "only {exercised} exercised cases");
}
