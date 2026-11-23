//! rusolve `mip_activation_dff`: dual-feasible-function rows for
//! activation knapsacks.
//!
//! Not in HiGHS. A row `sum_i a_i x_i - b y <= 0` (x, y binary, a_i > 0, b > 0) says: items
//! of total size at most `b` may be packed once the bin `y` is open. The LP relaxation only
//! sees the area `sum a_i x_i / b <= y`. For any dual feasible function `f` of the bin
//! packing problem (a set of sizes fitting into `b` weighs at most 1), `sum f(a_i) x_i <= y`
//! is valid. The function used here is the elementary one `f(a) = 1` for `2a > b`, `1/2`
//! for `2a = b`, `0` otherwise: two items above `b/2` never share a bin, and at most two
//! items of size exactly `b/2` do. Sources: Fekete and Schepers, "New classes of fast lower
//! bounds for bin packing problems" (Math. Programming 91, 2001), dual feasible functions;
//! Achterberg et al., "Presolve reductions in mixed integer programming" (INFORMS J.
//! Computing 32, 2020) on the same coefficient/conflict strengthening of knapsack rows with
//! binary activation. Rows of the form `sum a_i x_i <= b` (no activation variable) give the
//! same inequality with right-hand side 1.
//!
//! The rows are implied by the model, so they are added to the presolved model only (no
//! postsolve change). Default off.

use crate::prelude::*;

/// Append the rows to `lp` (column-wise or row-wise matrix); returns the number added.
pub fn add_activation_dff_rows(lp: &mut HighsLp) -> usize {
    let (n, m) = (lp.num_col_ as usize, lp.num_row_ as usize);
    if n == 0 || m == 0 || lp.integrality_.len() != n {
        return 0;
    }
    // a feasibility problem has no bound to strengthen (ns1116954: the rows only perturbed the search)
    if lp.col_cost_.iter().all(|&cost| cost == 0.0) {
        return 0;
    }
    let is_bin = |j: usize| {
        lp.integrality_[j] != HighsVarType::kContinuous
            && lp.col_lower_[j] == 0.0
            && lp.col_upper_[j] == 1.0
    };
    // candidate rows: only upper side finite, rhs 0 (activation) or positive (plain)
    let mut ent: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
    let a = &lp.a_matrix_;
    match a.format_ {
        MatrixFormat::kColwise => {
            for j in 0..n {
                for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                    ent[a.index_[p] as usize].push((j, a.value_[p]));
                }
            }
        }
        _ => return 0,
    }
    let mut starts: Vec<HighsInt> = vec![0];
    let mut index: Vec<HighsInt> = Vec::new();
    let mut value: Vec<f64> = Vec::new();
    let mut ub: Vec<f64> = Vec::new();
    let mut row_idx: Vec<HighsInt> = Vec::new();
    let mut row_val: Vec<f64> = Vec::new();
    for r in 0..m {
        if ent[r].len() < 3 {
            continue;
        }
        // `a x <= rhs`, or `a x >= lo` read as `-a x <= -lo`
        let (sign, rhs) = if lp.row_lower_[r] == -kHighsInf {
            (1.0, lp.row_upper_[r])
        } else if lp.row_upper_[r] == kHighsInf {
            (-1.0, -lp.row_lower_[r])
        } else {
            continue;
        };
        if !(rhs.is_finite()) || rhs != 0.0 || !ent[r].iter().all(|&(j, _)| is_bin(j)) {
            continue;
        }
        let neg: Vec<(usize, f64)> = ent[r]
            .iter()
            .map(|e| (e.0, sign * e.1))
            .filter(|e| e.1 < 0.0)
            .collect();
        // activation rows only: plain knapsacks gained nothing in the sample (neos-957323 lost 33%)
        let (cap, y) = if rhs == 0.0 && neg.len() == 1 {
            (-neg[0].1, Some(neg[0].0))
        } else {
            continue;
        };
        // item sizes relative to the capacity; items above the capacity are fixed at 0 by the row
        let items: Vec<(usize, f64)> = ent[r]
            .iter()
            .map(|e| (e.0, sign * e.1))
            .filter(|e| e.1 > 0.0)
            .map(|(j, v)| (j, v / cap))
            .collect();
        if items.len() < 2 {
            continue;
        }
        let area: f64 = items.iter().map(|e| e.1).sum();
        let mode = dff_mode();
        // candidate functions, scored by the total weight of the row against its area
        let mut scored: Vec<(f64, Dff)> = Vec::new();
        for f in dff_candidates(&items, mode) {
            let w: Vec<f64> = items.iter().map(|e| f.eval(e.1)).collect();
            let total: f64 = w.iter().sum();
            let nonzero = w.iter().filter(|&&v| v > 1e-12).count();
            let fractional = w.iter().any(|&v| v > 1e-9 && v < 1.0 - 1e-9);
            // plain rows without a fractional weight are cliques HiGHS finds itself
            if total <= area * (1.0 + 1e-6) + 1e-9 || nonzero < 2 || (y.is_none() && !fractional) {
                continue;
            }
            scored.push((total / area, f));
        }
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let mut kept: Vec<Dff> = Vec::new();
        let mut kept_w: Vec<Vec<f64>> = Vec::new();
        for (_, f) in scored {
            if kept.len() >= max_rows_per_capacity_row(mode) {
                break;
            }
            let w: Vec<f64> = items.iter().map(|e| f.eval(e.1)).collect();
            // skip a function whose weights equal one already kept
            if kept_w
                .iter()
                .any(|k| k.iter().zip(&w).all(|(a, b)| (a - b).abs() < 1e-9))
            {
                continue;
            }
            kept_w.push(w);
            kept.push(f);
        }
        for w in kept_w {
            for (k, &(j, _)) in items.iter().enumerate() {
                if w[k] > 1e-12 {
                    row_idx.push(j as HighsInt);
                    row_val.push(w[k]);
                }
            }
            if let Some(y) = y {
                row_idx.push(y as HighsInt);
                row_val.push(-1.0);
            }
            index.append(&mut row_idx);
            value.append(&mut row_val);
            starts.push(index.len() as HighsInt);
            ub.push(if y.is_some() { 0.0 } else { 1.0 });
        }
        continue;
    }
    let k = ub.len();
    if k == 0 {
        return 0;
    }
    let new_rows = HighsSparseMatrix {
        format_: MatrixFormat::kRowwise,
        num_col_: n as HighsInt,
        num_row_: k as HighsInt,
        start_: starts,
        p_end_: Vec::new(),
        index_: index,
        value_: value,
    };
    lp.a_matrix_.add_rows(&new_rows, None);
    lp.row_lower_.extend(std::iter::repeat_n(-kHighsInf, k));
    lp.row_upper_.extend(ub);
    if !lp.row_names_.is_empty() {
        lp.row_names_.resize(m + k, String::new());
    }
    lp.num_row_ += k as HighsInt;
    k
}

/// Rows kept per activation row: the two best distinct functions of the family.
fn dff_mode() -> u32 {
    1
}

fn max_rows_per_capacity_row(mode: u32) -> usize {
    if mode == 1 || mode == 0 { 1 } else { 2 }
}

/// A dual feasible function on item sizes `x = a / capacity` in (0, 1]: every set of items with
/// sizes summing to at most 1 has weights summing to at most 1.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Dff {
    /// Fekete and Schepers 2001, `u^(eps)`: 1 above `1 - eps`, 0 below `eps`, else `x`.
    Fs1(f64),
    /// Fekete and Schepers 2001, `f_k`: `x` if `(k+1) x` is an integer, else `floor((k+1) x) / k`.
    FsK(u32),
    /// Carlier and Neron 2007 / Alves et al. 2016 floor family: `floor(lambda x) / lambda`
    /// (the floors of a set of sizes summing to at most 1 sum to at most `lambda`).
    Floor(u32),
}

impl Dff {
    pub(crate) fn eval(self, x: f64) -> f64 {
        const T: f64 = 1e-9;
        match self {
            Dff::Fs1(eps) => {
                if x > 1.0 - eps + T {
                    1.0
                } else if x < eps - T {
                    0.0
                } else {
                    x
                }
            }
            Dff::FsK(k) => {
                let v = (k as f64 + 1.0) * x;
                if (v - v.round()).abs() <= T {
                    x
                } else {
                    v.floor() / k as f64
                }
            }
            Dff::Floor(l) => (l as f64 * x + T).floor() / l as f64,
        }
    }
}

pub(crate) fn dff_candidates(items: &[(usize, f64)], mode: u32) -> Vec<Dff> {
    // mode 0 is `f_1` (Fekete and Schepers 2001)
    if mode == 0 {
        return vec![Dff::Fs1(0.5)];
    }
    let mut out = vec![Dff::Fs1(0.5)];
    // u^(eps) breakpoints at the distinct item sizes up to 1/2
    let mut eps: Vec<f64> = items
        .iter()
        .map(|e| e.1)
        .filter(|&x| x <= 0.5 + 1e-9)
        .collect();
    eps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eps.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let step = (eps.len() / 6).max(1);
    for e in eps.iter().step_by(step) {
        out.push(Dff::Fs1(*e));
    }
    for k in 2..=4 {
        out.push(Dff::FsK(k));
    }
    for l in 2..=8 {
        out.push(Dff::Floor(l));
    }
    out
}
