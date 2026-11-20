//! Mod-k cuts (zero-half for k = 2): combinations of tight LP rows whose
//! integer part vanishes modulo k while the right-hand side does not, so the
//! combined row divided by k has a fractional right-hand side that c-MIR
//! rounds into a cut violated by up to (k−1)/k at the LP point.
//!
//! Follows HiGHS `HighsModkSeparator::separateLpSolution`
//! (reference/HiGHS/highs/mip/HighsModkSeparator.cpp:42–266) and the GF(k)
//! elimination `HighsGFkSolve::solve` (mip/HighsGFkSolve.h:203–434). HiGHS
//! is MIT licensed (Copyright (c) 2016–2024 the HiGHS authors); it was read
//! as a reference text and reimplemented. The underlying separation scheme
//! (maximally violated mod-k cuts from a congruence system over the tight
//! rows, restricted to columns off their bounds) is Caprara, Fischetti &
//! Letchford, "On the separation of maximally violated mod-k cuts", Math.
//! Program. 87 (2000).
//!
//! Deviation: HiGHS aggregates a `≥`-tight row with the weight it found for
//! the negated row (HighsModkSeparator.cpp:209–214 ignore the side), so the
//! row's slack enters the base row with the wrong sign and stays in the cut
//! (for k > 2 even the integer part is a different combination than the one
//! solved for). Here the weight carries the row's sign. Measured on our
//! medium benchmark instances (seeds 1–3): set_cover root gap
//! closed 78.1% signed vs 61.8% unsigned, every other family equal.

use super::generate::{CutGen, GeneratedCut, integral_scale, nearest_integer};
use super::relax::{Aggregator, Relaxation};
use super::transform::TransformedLp;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};

/// Moduli tried in order; the first one that adds a cut ends the round
/// (HighsModkSeparator.cpp:244–265).
const MODULI: [u32; 4] = [2, 3, 5, 7];

/// Rows whose scaled integer part needs a larger factor are skipped
/// (HighsModkSeparator.cpp:166).
const MAX_INT_SCALE: f64 = 1e6;

/// `a^(k−2) mod k`, the inverse of `a ≠ 0` in GF(k) for prime `k` (Fermat;
/// HighsGFkSolve.h:23–57).
fn inverse(a: u32, k: u32) -> u32 {
    let mut r = 1;
    for _ in 0..k - 2 {
        r = r * a % k;
    }
    r
}

/// Solutions `w` of `Σ_u w_u·A[e][u] ≡ rhs[e] (mod k)` over a sparse system
/// whose unknowns are columns `u` and equations rows `e`, enumerated by
/// LU-style elimination and then one basis swap per nonbasic unknown
/// (HighsGFkSolve.h:203–434).
pub(crate) struct GfkSystem {
    k: u32,
    num_unknowns: usize,
    /// Equations, entries sorted by unknown with values in `1..k`.
    rows: Vec<Vec<(usize, u32)>>,
    rhs: Vec<u32>,
    /// Equations that held a nonzero for each unknown at some point
    /// (superset; may repeat).
    col_rows: Vec<Vec<usize>>,
    /// Nonzeros of each unknown among the not-yet-pivoted equations.
    colsize: Vec<usize>,
}

impl GfkSystem {
    /// `columns[u]` lists `(equation, integer coefficient)` of unknown `u`
    /// (HighsGFkSolve::fromCSC, HighsGFkSolve.h:146–186).
    pub(crate) fn new(k: u32, num_eqs: usize, columns: &[Vec<(usize, i64)>]) -> Self {
        let mut rows: Vec<Vec<(usize, u32)>> = vec![Vec::new(); num_eqs];
        let mut col_rows = vec![Vec::new(); columns.len()];
        let mut colsize = vec![0; columns.len()];
        for (u, col) in columns.iter().enumerate() {
            for &(e, a) in col {
                let v = a.rem_euclid(i64::from(k)) as u32;
                if v != 0 {
                    rows[e].push((u, v));
                    col_rows[u].push(e);
                    colsize[u] += 1;
                }
            }
        }
        Self {
            k,
            num_unknowns: columns.len(),
            rows,
            rhs: vec![0; num_eqs],
            col_rows,
            colsize,
        }
    }

    pub(crate) fn set_rhs(&mut self, e: usize, v: i64) {
        self.rhs[e] = (v.unsigned_abs() % u64::from(self.k)) as u32;
    }

    fn value(&self, e: usize, u: usize) -> Option<u32> {
        let row = &self.rows[e];
        row.binary_search_by_key(&u, |&(c, _)| c)
            .ok()
            .map(|p| row[p].1)
    }

    /// Adds `scale·rows[p]` to `rows[e]` (both mod k), keeping `colsize` and
    /// `col_rows` current for the unpivoted equation `e`.
    fn axpy(&mut self, e: usize, p: usize, scale: u32) {
        let k = self.k;
        let (a, b) = (&self.rows[e], &self.rows[p]);
        let mut out = Vec::with_capacity(a.len() + b.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() || j < b.len() {
            let ca = a.get(i).map_or(usize::MAX, |x| x.0);
            let cb = b.get(j).map_or(usize::MAX, |x| x.0);
            if ca < cb {
                out.push(a[i]);
                i += 1;
            } else if cb < ca {
                let v = scale * b[j].1 % k;
                if v != 0 {
                    out.push((cb, v));
                    self.colsize[cb] += 1;
                    self.col_rows[cb].push(e);
                }
                j += 1;
            } else {
                let v = (a[i].1 + scale * b[j].1) % k;
                if v != 0 {
                    out.push((ca, v));
                } else {
                    self.colsize[ca] -= 1;
                }
                i += 1;
                j += 1;
            }
        }
        self.rows[e] = out;
        self.rhs[e] = (self.rhs[e] + scale * self.rhs[p]) % k;
    }

    /// Calls `report` with every solution found (sparse `(unknown, value)`,
    /// nonzero values only). Nothing is reported when the system is
    /// inconsistent.
    pub(crate) fn solve(mut self, mut report: impl FnMut(&mut Vec<(usize, u32)>)) {
        let k = self.k;
        let num_eqs = self.rows.len();
        let n = self.num_unknowns;
        let max_pivot = num_eqs.min(n);
        let mut used = vec![false; num_eqs];
        // 0 nonbasic, 1 basic, -1 queued to enter (HighsGFkSolve.h:345–368).
        let mut status = vec![0i8; n];
        let mut factor: Vec<(usize, usize)> = Vec::with_capacity(max_pivot);
        let mut heap: BinaryHeap<Reverse<(usize, usize)>> =
            (0..n).map(|u| Reverse((self.colsize[u], u))).collect();
        let mut cand = Vec::new();

        // Pivot on the sparsest unknown, in its shortest unpivoted equation
        // (HighsGFkSolve.h:214–313).
        while let Some(Reverse((size, pc))) = heap.pop() {
            if self.colsize[pc] == 0 {
                continue;
            }
            if self.colsize[pc] != size {
                heap.push(Reverse((self.colsize[pc], pc)));
                continue;
            }
            cand.clear();
            cand.extend(
                self.col_rows[pc]
                    .iter()
                    .copied()
                    .filter(|&e| !used[e] && self.value(e, pc).is_some()),
            );
            cand.sort_unstable();
            cand.dedup();
            let pr = *cand
                .iter()
                .min_by_key(|&&e| (self.rows[e].len(), e))
                .expect("colsize counts unpivoted equations");
            let pinv = inverse(self.value(pr, pc).expect("pivot present"), k);
            let old_sizes: Vec<(usize, usize)> = self.rows[pr]
                .iter()
                .map(|&(c, _)| (c, self.colsize[c]))
                .collect();
            for &e in &cand {
                if e == pr {
                    continue;
                }
                let a = self.value(e, pc).expect("candidate holds the pivot column");
                self.axpy(e, pr, pinv * (k - a) % k);
            }
            status[pc] = 1;
            used[pr] = true;
            factor.push((pr, pc));
            if factor.len() == max_pivot {
                break;
            }
            for (c, old) in old_sizes {
                self.colsize[c] -= 1;
                if self.colsize[c] > 0 && self.colsize[c] < old {
                    heap.push(Reverse((self.colsize[c], c)));
                }
            }
        }

        // Dependent equations must have a zero right-hand side
        // (HighsGFkSolve.h:316–336).
        if (0..num_eqs).any(|e| !used[e] && self.rhs[e] != 0) {
            return;
        }

        // Each nonbasic unknown enters once, in the latest factor equation
        // that holds it (HighsGFkSolve.h:345–368).
        let mut swaps = Vec::new();
        for (i, &(e, _)) in factor.iter().enumerate().rev() {
            for &(u, _) in &self.rows[e] {
                if status[u] == 0 {
                    status[u] = -1;
                    swaps.push((i, u));
                }
            }
        }

        let mut x = vec![0u32; n];
        let mut solution = Vec::new();
        let mut next_swap = 0;
        loop {
            for &(u, _) in &solution {
                x[u] = 0;
            }
            solution.clear();
            // Back substitution, nonbasic unknowns at zero
            // (HighsGFkSolve.h:375–418).
            for &(e, pc) in factor.iter().rev() {
                let mut s = 0u64;
                for &(u, a) in &self.rows[e] {
                    s += u64::from(a) * u64::from(x[u]);
                }
                let s = (s % u64::from(k)) as u32;
                let a = self.value(e, pc).expect("factor equation holds its pivot");
                let v = (self.rhs[e] + k - s) * inverse(a, k) % k;
                if v != 0 {
                    x[pc] = v;
                    solution.push((pc, v));
                }
            }
            report(&mut solution);
            // `report` may reorder the entries; `x` is reset from them.
            let Some(&(i, entering)) = swaps.get(next_swap) else {
                break;
            };
            factor[i].1 = entering;
            next_swap += 1;
        }
    }
}

/// One tight LP row in the congruence system.
struct SystemRow {
    row: usize,
    /// Multiplier making the transformed integer part integral.
    scale: f64,
    /// −1 when the row entered as its `≥` side negated.
    sign: f64,
}

/// Separates mod-k cuts at `relax.x`; `add` receives each generated cut and
/// returns whether the pool accepted it.
#[allow(clippy::too_many_arguments)]
pub(crate) fn separate(
    relax: &Relaxation,
    trans: &mut TransformedLp,
    agg: &mut Aggregator,
    cutgen: &mut CutGen,
    feastol: f64,
    epsilon: f64,
    mut add: impl FnMut(GeneratedCut) -> bool,
) {
    let n = relax.num_col;
    let m = relax.rows.len();
    agg.ensure_dim(n + m);
    let max_int_row_len = 1000 + n / 10;

    // Rows with a continuous column off its (variable) bound cannot be in
    // the system (HighsModkSeparator.cpp:49–59).
    let off_bound = |j: usize| !relax.integer[j] && trans.bound_distance(j) != 0.0;
    let skip: Vec<bool> = relax
        .rows
        .iter()
        .map(|r| r.inds.iter().any(|&j| off_bound(j)))
        .collect();

    let mut system: Vec<SystemRow> = Vec::new();
    let mut columns: Vec<Vec<(usize, i64)>> = Vec::new();
    let mut num_nonzero_rhs = 0;
    let mut scale_vals = Vec::new();
    for (i, r) in relax.rows.iter().enumerate() {
        if skip[i] {
            continue;
        }
        let act = relax.activity[i];
        let (sign, rhs) = if r.row_upper - act <= feastol {
            (1.0, r.row_upper)
        } else if act - r.row_lower <= feastol {
            (-1.0, -r.row_lower)
        } else {
            continue;
        };
        let vals: Vec<f64> = r.vals.iter().map(|v| sign * v).collect();
        let Some(t) = trans.transform(relax, r.inds.clone(), vals, rhs, false, true) else {
            continue;
        };
        if t.inds.iter().any(|&j| j >= n) {
            continue;
        }
        let len = t.inds.len();
        if len > max_int_row_len {
            let int_len = (0..len)
                .filter(|&p| t.solval[p] > feastol && relax.integer[t.inds[p]])
                .count();
            if int_len > max_int_row_len || (int_len == 0 && t.rhs.abs() <= epsilon) {
                continue;
            }
        }

        // Integer part over the columns off zero (HighsModkSeparator.cpp:
        // 154–196); columns at zero do not affect the violation.
        let mut col = Vec::new();
        let (scale, intrhs) = if r.integral {
            for p in 0..len {
                if t.solval[p] > feastol {
                    col.push((t.inds[p], nearest_integer(t.vals[p]) as i64));
                }
            }
            (1.0, nearest_integer(t.rhs) as i64)
        } else {
            scale_vals.clear();
            for p in 0..len {
                if relax.integer[t.inds[p]] && t.solval[p] > feastol {
                    scale_vals.push(t.vals[p]);
                }
            }
            if t.rhs.abs() > epsilon {
                scale_vals.push(-t.rhs);
            }
            if scale_vals.is_empty() {
                continue;
            }
            let s = integral_scale(&scale_vals, feastol, epsilon);
            if s == 0.0 || s > MAX_INT_SCALE {
                continue;
            }
            for p in 0..len {
                if relax.integer[t.inds[p]] && t.solval[p] > feastol {
                    col.push((t.inds[p], nearest_integer(s * t.vals[p]) as i64));
                }
            }
            (s, nearest_integer(s * t.rhs) as i64)
        };
        num_nonzero_rhs += usize::from(intrhs != 0);
        col.push((n, intrhs));
        columns.push(col);
        system.push(SystemRow {
            row: i,
            scale,
            sign,
        });
    }
    if system.is_empty() || num_nonzero_rhs == 0 {
        return;
    }

    let mut used_weights: HashSet<Vec<(usize, u32)>> = HashSet::new();
    for k in MODULI {
        used_weights.clear();
        let mut added = false;
        let mut gfk = GfkSystem::new(k, n + 1, &columns);
        gfk.set_rhs(n, 1);
        gfk.solve(|weights| {
            // Single-row cuts are the aggregation heuristic's business.
            if weights.is_empty() {
                return;
            }
            weights.sort_unstable();
            if !used_weights.insert(weights.clone()) {
                return;
            }
            let kf = f64::from(k);
            // Weights (k−1)·w mod k make the rhs fraction (k−1)/k
            // (HighsModkSeparator.cpp:205–222).
            agg.clear();
            for &(u, w) in weights.iter() {
                let s = &system[u];
                let weight = s.scale * f64::from(w * (k - 1) % k) / kf;
                agg.add_row(relax, s.row, s.sign * weight);
            }
            let (inds, vals) = agg.current(n, epsilon, false);
            if let Some(cut) = cutgen.generate(relax, trans, inds, vals, 0.0, true) {
                added |= add(cut);
            }
            // Weights w/k, negated (HighsModkSeparator.cpp:224–237).
            if k != 2 {
                agg.clear();
                for &(u, w) in weights.iter() {
                    let s = &system[u];
                    agg.add_row(relax, s.row, s.sign * s.scale * f64::from(w) / kf);
                }
            }
            let (inds, vals) = agg.current(n, epsilon, true);
            if let Some(cut) = cutgen.generate(relax, trans, inds, vals, 0.0, true) {
                added |= add(cut);
            }
            agg.clear();
        });
        if added {
            return;
        }
    }
}

#[cfg(test)]
mod gfk_tests {
    use super::GfkSystem;
    use rusolve_lp::rng::XorShift64;

    /// Every reported solution satisfies every equation mod k, and systems
    /// built around a planted solution report at least one.
    #[test]
    fn gfk_solutions_satisfy_the_system() {
        let mut rng = XorShift64::new(7_20261022);
        for case in 0..600 {
            let k = [2u32, 3, 5, 7][case % 4];
            let num_eqs = 1 + (rng.next_u64() % 6) as usize;
            let num_unk = 1 + (rng.next_u64() % 8) as usize;
            let mut cols: Vec<Vec<(usize, i64)>> = (0..num_unk)
                .map(|_| {
                    let mut col = Vec::new();
                    for e in 0..num_eqs {
                        if rng.next_u64().is_multiple_of(2) {
                            col.push((e, (rng.next_u64() % 15) as i64 - 7));
                        }
                    }
                    col
                })
                .collect();
            let dense = |cols: &[Vec<(usize, i64)>], e: usize, u: usize| {
                cols[u]
                    .iter()
                    .find(|x| x.0 == e)
                    .map_or(0, |x| x.1.rem_euclid(i64::from(k)))
            };
            // Half the cases plant a solution in the right-hand side.
            let planted = case % 2 == 0;
            let rhs: Vec<i64> = if planted {
                let w: Vec<i64> = (0..num_unk)
                    .map(|_| (rng.next_u64() % u64::from(k)) as i64)
                    .collect();
                (0..num_eqs)
                    .map(|e| {
                        (0..num_unk).map(|u| dense(&cols, e, u) * w[u]).sum::<i64>() % i64::from(k)
                    })
                    .collect()
            } else {
                (0..num_eqs)
                    .map(|_| (rng.next_u64() % u64::from(k)) as i64)
                    .collect()
            };
            // Unknown coefficients that are multiples of k must vanish.
            cols.push(vec![(0, i64::from(k))]);
            let mut sys = GfkSystem::new(k, num_eqs, &cols);
            for (e, &b) in rhs.iter().enumerate() {
                sys.set_rhs(e, b);
            }
            let mut found = 0;
            sys.solve(|sol| {
                found += 1;
                let mut w = vec![0i64; cols.len()];
                for &(u, v) in sol.iter() {
                    assert!(v > 0 && v < k);
                    w[u] = i64::from(v);
                }
                for (e, &b) in rhs.iter().enumerate() {
                    let s: i64 = (0..cols.len()).map(|u| dense(&cols, e, u) * w[u]).sum();
                    assert_eq!(s.rem_euclid(i64::from(k)), b, "case {case} k {k} eq {e}");
                }
            });
            if planted {
                assert!(found > 0, "case {case}: planted system reported nothing");
            }
        }
    }
}
