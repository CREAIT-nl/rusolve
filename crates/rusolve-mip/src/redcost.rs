//! Reduced-cost fixing.
//!
//! With an incumbent, every improving point has objective at most the upper
//! limit `U`. If a relaxation over the box `[l, u]` has duals `y` and reduced
//! costs `d`, its Lagrangian bound `L` ([`LpDuals::bound`]) gives
//! `cᵀx ≥ L + d_j (x_j − l_j)` for `d_j > 0`, so improving points satisfy
//! `x_j ≤ l_j + (U − L)/d_j`, and symmetrically `x_j ≥ u_j + (U − L)/d_j` for
//! `d_j < 0` (Nemhauser & Wolsey, *Integer and Combinatorial Optimization*,
//! 1988, §II.5.1; Achterberg, *Constraint Integer Programming*, 2007, §7.7).
//!
//! The rules, constants and integer rounding follow HiGHS (MIT, read as a
//! reference text and reimplemented):
//!
//! * node-local fixing — `HighsRedcostFixing::propagateRedCost`
//!   (reference/HiGHS/highs/mip/HighsRedcostFixing.cpp:79–137): skip columns
//!   with `|d_j| ≤ max(10·feastol, ε·gap)`, `ε = 1e-9` (`small_matrix_value`,
//!   HighsMipSolverData.cpp:756), tighten only when `d_j (u_j − l_j)`
//!   exceeds the gap, round `⌊gap/d_j + l_j + feastol⌋` /
//!   `⌈gap/d_j + u_j − feastol⌉`; called after every node LP with an upper
//!   limit (HighsSearch.cpp:1000–1029, which re-evaluates the node when
//!   bounds moved);
//! * root fixing — `addRootRedcost` / `propagateRootRedcost` (:244–323,
//!   :38–77): the root LP's reduced costs (`|d_j| > feastol`) are kept and
//!   re-applied to the global domain whenever the upper limit drops
//!   (HighsMipSolverData.cpp:1576–1578, :1946–1953; HighsSeparation.cpp:80,
//!   :169). HiGHS precomputes "lurking bounds" (cutoff → bound tables); we
//!   store the reduced costs and the root bound and evaluate the same
//!   formula on demand, which gives the same bounds without the tables.
//!
//! Unlike HiGHS we use the Lagrangian bound `min(L, z_LP)` in place of the
//! LP objective, so the fixing stays valid for any sign-feasible `y` (the
//! duals of warm solves are recomputed from the basis, see
//! [`crate::engine::LpEngine::duals`]).
// TODO(perf): HiGHS also turns binary fixings into conflicts through a dual
// proof (HighsRedcostFixing.cpp:139–193) and computes basic-degenerate duals
// first (`computeBasicDegenerateDuals`, HighsSearch.cpp:1002); neither exists
// here yet (no conflict pool).

use crate::engine::LpDuals;

/// No improving point exists in the box (the node, or the whole problem at
/// the root, can be pruned).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Infeasible;

/// HiGHS `small_matrix_value` (`mipdata.epsilon`).
const EPSILON: f64 = 1e-9;

/// Root snapshots kept (the first root LP's and the latest ones).
// TODO(perf): HiGHS keeps one lurking-bound table fed by every root LP
// (each cut round); we keep the first root LP and the post-cut LP. Revisit
// with per-round snapshots if root fixing proves weak on the miplib tier.
const MAX_ROOT_SNAPSHOTS: usize = 8;

/// Reduced-cost fixing counters.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RedcostStats {
    /// Node LPs whose duals were available while an incumbent existed.
    pub nodes_with_duals: u64,
    /// Node-local bound tightenings and the nodes that had at least one.
    pub local_tightenings: u64,
    pub nodes_tightened: u64,
    /// Global bound tightenings from the stored root reduced costs.
    pub root_tightenings: u64,
    /// Node LPs re-solved because the tightened box no longer held the LP
    /// point, and nodes pruned by the fixing (with its propagation).
    pub resolves: u64,
    pub pruned: u64,
}

/// The tightening of one integer column with reduced cost `d` against the
/// objective gap `gap = U − L`, or `None`. `Some((true, v))` = new upper
/// bound `v`, `Some((false, v))` = new lower bound.
fn tightening(d: f64, gap: f64, lo: f64, up: f64, feastol: f64) -> Option<(bool, f64)> {
    let max_increase = d * (up - lo);
    if d > 0.0 && max_increase > gap && lo.is_finite() {
        let ub = (gap / d + lo + feastol).floor();
        (ub < up).then_some((true, ub))
    } else if d < 0.0 && max_increase < -gap && up.is_finite() {
        let lb = (gap / d + up - feastol).ceil();
        (lb > lo).then_some((false, lb))
    } else {
        None
    }
}

/// Applies a tightening; `false` if it empties the domain.
fn apply(j: usize, t: (bool, f64), lower: &mut [f64], upper: &mut [f64]) -> bool {
    match t {
        (true, v) => upper[j] = v,
        (false, v) => lower[j] = v,
    }
    lower[j] <= upper[j]
}

/// Node-local reduced-cost fixing on the box the relaxation was solved over:
/// `z` is its objective, `cutoff` the upper limit. Tightens integer columns
/// in place and lists them in `changed`. `Err(Infeasible)` = no improving point in
/// the box.
#[allow(clippy::too_many_arguments)]
pub fn fix_local(
    duals: &LpDuals,
    z: f64,
    cutoff: f64,
    integer: &[bool],
    lower: &mut [f64],
    upper: &mut [f64],
    feastol: f64,
    changed: &mut Vec<usize>,
) -> Result<(), Infeasible> {
    let Some(bound) = duals.bound(lower, upper) else {
        return Ok(());
    };
    let gap = cutoff - bound.min(z);
    if gap < 0.0 {
        return Err(Infeasible);
    }
    let tol = (10.0 * feastol).max(EPSILON * gap);
    for (j, &d) in duals.reduced_costs.iter().enumerate() {
        if !integer[j] || lower[j] == upper[j] || d.abs() <= tol {
            continue;
        }
        if let Some(t) = tightening(d, gap, lower[j], upper[j], feastol) {
            changed.push(j);
            if !apply(j, t, lower, upper) {
                return Err(Infeasible);
            }
        }
    }
    Ok(())
}

struct RootSnapshot {
    /// Lagrangian bound of the root relaxation over the root box.
    bound: f64,
    /// `(column, d_j, anchor)`: integer columns with `|d_j| > feastol`; the
    /// anchor is the root lower bound for `d_j > 0`, the upper for `d_j < 0`.
    cols: Vec<(usize, f64, f64)>,
}

/// Root reduced costs, re-applied to the global domain at every new upper
/// limit (HiGHS `HighsRedcostFixing` lurking bounds).
#[derive(Default)]
pub struct RootRedcost {
    snapshots: Vec<RootSnapshot>,
}

impl RootRedcost {
    /// Stores the duals of a root relaxation with objective `z`, solved over
    /// the global box `[lower, upper]` (HiGHS `addRootRedcost`).
    pub fn add(
        &mut self,
        duals: &LpDuals,
        z: f64,
        lower: &[f64],
        upper: &[f64],
        integer: &[bool],
        feastol: f64,
    ) {
        let Some(bound) = duals.bound(lower, upper) else {
            return;
        };
        let cols: Vec<(usize, f64, f64)> = duals
            .reduced_costs
            .iter()
            .enumerate()
            .filter_map(|(j, &d)| {
                if !integer[j] || lower[j] == upper[j] {
                    return None;
                }
                let anchor = if d > feastol {
                    lower[j]
                } else if d < -feastol {
                    upper[j]
                } else {
                    return None;
                };
                anchor.is_finite().then_some((j, d, anchor))
            })
            .collect();
        if cols.is_empty() {
            return;
        }
        if self.snapshots.len() == MAX_ROOT_SNAPSHOTS {
            self.snapshots.remove(1);
        }
        self.snapshots.push(RootSnapshot {
            bound: bound.min(z),
            cols,
        });
    }

    /// Tightens the global box for the upper limit `cutoff` (HiGHS
    /// `propagateRootRedcost`); returns the number of tightenings, listing
    /// the columns in `changed`. `Err(Infeasible)` = no improving point exists.
    pub fn apply(
        &self,
        cutoff: f64,
        lower: &mut [f64],
        upper: &mut [f64],
        feastol: f64,
        changed: &mut Vec<usize>,
    ) -> Result<u64, Infeasible> {
        let mut n = 0;
        if !cutoff.is_finite() {
            return Ok(0);
        }
        for s in &self.snapshots {
            let gap = cutoff - s.bound;
            if gap < 0.0 {
                return Err(Infeasible);
            }
            for &(j, d, anchor) in &s.cols {
                // The formula against the anchor (root bound), then only
                // kept if it beats the current global bound.
                let (lo, up) = if d > 0.0 {
                    (anchor, f64::INFINITY)
                } else {
                    (f64::NEG_INFINITY, anchor)
                };
                let Some(t) = tightening(d, gap, lo, up, feastol) else {
                    continue;
                };
                let better = match t {
                    (true, v) => v < upper[j],
                    (false, v) => v > lower[j],
                };
                if better {
                    n += 1;
                    changed.push(j);
                    if !apply(j, t, lower, upper) {
                        return Err(Infeasible);
                    }
                }
            }
        }
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounding_matches_the_textbook_rule_with_feastol_slack() {
        // x at lower 0 with d = 2, gap 3: x ≤ 1.5 ⇒ ub 1.
        assert_eq!(tightening(2.0, 3.0, 0.0, 5.0, 1e-6), Some((true, 1.0)));
        // gap/d hits an integer up to round-off: keep it (feastol slack).
        assert_eq!(
            tightening(3.0, 6.0 - 1e-9, 0.0, 5.0, 1e-6),
            Some((true, 2.0))
        );
        // x at upper 5 with d = −2, gap 3: x ≥ 3.5 ⇒ lb 4.
        assert_eq!(tightening(-2.0, 3.0, 0.0, 5.0, 1e-6), Some((false, 4.0)));
        // The whole range fits in the gap: nothing to do.
        assert_eq!(tightening(0.5, 3.0, 0.0, 5.0, 1e-6), None);
        // Zero gap fixes the column at its bound.
        assert_eq!(tightening(1.0, 0.0, 2.0, 5.0, 1e-6), Some((true, 2.0)));
    }
}
