//! Shared by `examples/tsp_lazy.rs`, `examples/user_cuts.rs` and `tests/lazy_cuts.rs`.
#![allow(dead_code)]

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::hooks::UserCut;
use rusolve_ekk::mip::lazy::LazyConstraints;
use rusolve_ekk::prelude::*;

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

/// One row `lo <= sum val*x[col] <= hi`.
pub struct Row {
    pub lo: f64,
    pub hi: f64,
    pub terms: Vec<(usize, f64)>,
}

/// A minimisation MIP over binaries (`maximize` flips the sense).
#[allow(clippy::field_reassign_with_default)]
pub fn build(cost: Vec<f64>, rows: &[Row], maximize: bool) -> HighsLp {
    let n = cost.len();
    let mut lp = HighsLp::default();
    lp.num_col_ = n as HighsInt;
    lp.num_row_ = rows.len() as HighsInt;
    lp.sense_ = if maximize {
        ObjSense::kMaximize
    } else {
        ObjSense::kMinimize
    };
    lp.col_cost_ = cost;
    lp.col_lower_ = vec![0.0; n];
    lp.col_upper_ = vec![1.0; n];
    lp.row_lower_ = rows.iter().map(|r| r.lo).collect();
    lp.row_upper_ = rows.iter().map(|r| r.hi).collect();
    lp.integrality_ = vec![HighsVarType::kInteger; n];
    let mut cols: Vec<Vec<(HighsInt, f64)>> = vec![Vec::new(); n];
    for (i, r) in rows.iter().enumerate() {
        for &(cv, v) in &r.terms {
            cols[cv].push((i as HighsInt, v));
        }
    }
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = n as HighsInt;
    lp.a_matrix_.num_row_ = rows.len() as HighsInt;
    let mut start = vec![0 as HighsInt];
    let (mut index, mut value) = (Vec::new(), Vec::new());
    for col in &cols {
        for &(i, v) in col {
            index.push(i);
            value.push(v);
        }
        start.push(index.len() as HighsInt);
    }
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp
}

/// A silent, single-thread, exact-gap solver holding `lp`. `presolve` off keeps hook column
/// indices equal to the model's.
pub fn quiet_with(lp: HighsLp, presolve: bool) -> Highs {
    let mut h = Highs::default();
    h.set_option_value_bool("output_flag", false);
    h.set_option_value_str("presolve", if presolve { "choose" } else { "off" });
    h.set_option_value_double("mip_rel_gap", 0.0);
    h.set_option_value_double("mip_abs_gap", 0.0);
    h.set_option_value_int("threads", 1);
    h.pass_model(lp);
    h
}

pub fn quiet(lp: HighsLp) -> Highs {
    quiet_with(lp, true)
}

// ---------------------------------------------------------------- TSP (lazy constraints)

pub fn cities(n: usize, seed: u64) -> Vec<(f64, f64)> {
    let mut s = seed * 7919 + 17;
    (0..n)
        .map(|_| ((lcg(&mut s) % 100) as f64, (lcg(&mut s) % 100) as f64))
        .collect()
}

fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

/// Column of the undirected edge {i, j}, i < j, in the TSP model.
pub fn edge_col(n: usize, i: usize, j: usize) -> usize {
    let (i, j) = if i < j { (i, j) } else { (j, i) };
    // edges ordered (0,1) (0,2) ... (0,n-1) (1,2) ...
    i * n - i * (i + 1) / 2 + (j - i - 1)
}

/// Symmetric TSP with degree-2 rows. `full` also lists every subtour row for 3 <= |S| <= n/2.
pub fn tsp_model(pts: &[(f64, f64)], full: bool) -> HighsLp {
    let n = pts.len();
    let mut cost = vec![0.0; n * (n - 1) / 2];
    for i in 0..n {
        for j in i + 1..n {
            cost[edge_col(n, i, j)] = dist(pts[i], pts[j]);
        }
    }
    let mut rows = Vec::new();
    for i in 0..n {
        let terms = (0..n)
            .filter(|&j| j != i)
            .map(|j| (edge_col(n, i, j), 1.0))
            .collect();
        rows.push(Row {
            lo: 2.0,
            hi: 2.0,
            terms,
        });
    }
    if full {
        for mask in 0u32..(1 << n) {
            let s: Vec<usize> = (0..n).filter(|&i| mask >> i & 1 == 1).collect();
            if s.len() >= 3 && s.len() <= n / 2 {
                rows.push(subtour_row(n, &s));
            }
        }
    }
    build(cost, &rows, false)
}

fn subtour_row(n: usize, s: &[usize]) -> Row {
    let mut terms = Vec::new();
    for (a, &i) in s.iter().enumerate() {
        for &j in &s[a + 1..] {
            terms.push((edge_col(n, i, j), 1.0));
        }
    }
    Row {
        lo: -kHighsInf,
        hi: s.len() as f64 - 1.0,
        terms,
    }
}

/// Connected components of the chosen edges (x > 0.5).
pub fn tour_components(x: &[f64], n: usize) -> Vec<Vec<usize>> {
    let mut comp = vec![usize::MAX; n];
    let mut out: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if comp[start] != usize::MAX {
            continue;
        }
        let id = out.len();
        let mut stack = vec![start];
        comp[start] = id;
        let mut members = Vec::new();
        while let Some(u) = stack.pop() {
            members.push(u);
            for v in 0..n {
                if v != u && comp[v] == usize::MAX && x[edge_col(n, u, v)] > 0.5 {
                    comp[v] = id;
                    stack.push(v);
                }
            }
        }
        members.sort();
        out.push(members);
    }
    out
}

/// The lazy-constraint oracle: one subtour-elimination row per subtour of the candidate.
pub struct SubtourOracle {
    n: usize,
}
impl SubtourOracle {
    pub fn new(n: usize) -> Self {
        SubtourOracle { n }
    }
}
impl LazyConstraints for SubtourOracle {
    fn lazy_constraints(&mut self, sol: &[f64]) -> Vec<UserCut> {
        let comps = tour_components(sol, self.n);
        if comps.len() == 1 {
            return Vec::new();
        }
        comps
            .iter()
            .filter(|s| s.len() >= 3)
            .map(|s| {
                let r = subtour_row(self.n, s);
                UserCut {
                    indices: r.terms.iter().map(|t| t.0 as HighsInt).collect(),
                    values: vec![1.0; r.terms.len()],
                    rhs: r.hi,
                }
            })
            .collect()
    }
}

/// Optimal tour length by trying every permutation (city 0 fixed).
pub fn brute_force_tour(pts: &[(f64, f64)]) -> f64 {
    fn rec(
        pts: &[(f64, f64)],
        perm: &mut Vec<usize>,
        used: &mut Vec<bool>,
        len: f64,
        best: &mut f64,
    ) {
        if len >= *best {
            return;
        }
        if perm.len() == pts.len() {
            let total = len + dist(pts[*perm.last().unwrap()], pts[0]);
            if total < *best {
                *best = total;
            }
            return;
        }
        for j in 1..pts.len() {
            if !used[j] {
                used[j] = true;
                let last = *perm.last().unwrap();
                perm.push(j);
                rec(pts, perm, used, len + dist(pts[last], pts[j]), best);
                perm.pop();
                used[j] = false;
            }
        }
    }
    let mut best = f64::INFINITY;
    let mut used = vec![false; pts.len()];
    used[0] = true;
    rec(pts, &mut vec![0], &mut used, 0.0, &mut best);
    best
}

// ---------------------------------------------------------------- knapsack covers (user cuts)

pub fn knapsack_data() -> (Vec<f64>, Vec<f64>, f64) {
    let mut s = 4242u64;
    let n = 14;
    let weights: Vec<f64> = (0..n).map(|_| (10 + lcg(&mut s) % 30) as f64).collect();
    let values: Vec<f64> = weights
        .iter()
        .map(|w| w + (lcg(&mut s) % 7) as f64)
        .collect();
    let cap = (weights.iter().sum::<f64>() * 0.4).floor();
    (values, weights, cap)
}

/// Every minimal cover C of the knapsack row (sum_C w > cap, and dropping any item fits).
pub fn minimal_covers(weights: &[f64], cap: f64) -> Vec<Vec<usize>> {
    let n = weights.len();
    let mut out = Vec::new();
    for mask in 1u32..(1 << n) {
        let items: Vec<usize> = (0..n).filter(|&i| mask >> i & 1 == 1).collect();
        let total: f64 = items.iter().map(|&i| weights[i]).sum();
        if total > cap && items.iter().all(|&i| total - weights[i] <= cap) {
            out.push(items);
        }
    }
    out
}

/// Cover inequality `sum_C x <= |C| - 1` as a user cut.
pub fn cover_cuts(covers: &[Vec<usize>]) -> Vec<UserCut> {
    covers
        .iter()
        .map(|cv| UserCut {
            indices: cv.iter().map(|&i| i as HighsInt).collect(),
            values: vec![1.0; cv.len()],
            rhs: cv.len() as f64 - 1.0,
        })
        .collect()
}

/// The knapsack row plus (optionally) the cover cuts as model rows. Presolve is off so the
/// separator hook sees model indices.
pub fn knapsack_model(values: &[f64], weights: &[f64], cap: f64, covers: &[Vec<usize>]) -> HighsLp {
    let mut rows = vec![Row {
        lo: -kHighsInf,
        hi: cap,
        terms: weights.iter().copied().enumerate().collect(),
    }];
    for cv in cover_cuts(covers) {
        rows.push(Row {
            lo: -kHighsInf,
            hi: cv.rhs,
            terms: cv.indices.iter().map(|&i| (i as usize, 1.0)).collect(),
        });
    }
    build(values.to_vec(), &rows, true)
}
