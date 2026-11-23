//! Root probing: bound strengthening on binary implications.
//!
//! The clique table alone is correct but seedless — a clique fixes only
//! when a member sits at its literal in the global box, and nothing at the
//! root puts a binary there. HiGHS's seed for the whole incumbent →
//! fixings → restart chain is probing: fix one binary at a time, propagate,
//! and read off what the fix implies ([Savelsbergh 1994] §"probing";
//! HiGHS `HighsImplications::computeImplications` /
//! `HPresolve::runProbing`, reference/HiGHS/highs/mip/HighsImplications.cpp
//! :30–389 and presolve/HPresolve.cpp:1601–1840, MIT; read as a reference
//! text and reimplemented).
//!
//! What one probe does (HiGHS `computeImplications`, simplified to the
//! root-only scope):
//!
//! * fix binary `j` to `val` in a scratch copy of the global box and
//!   propagate; an infeasible propagation proves `j = 1 − val` in *every*
//!   feasible point — a global fixing (HiGHS backtracks and calls
//!   `cliquetable.vertexInfeasible`, which fixes the opposite literal).
//! * otherwise the scratch box's diffs are `j = val`'s implications:
//!   binaries fixed under the probe become 2-cliques `{j = val, k = ¬b}`
//!   in the clique table (HiGHS `cliquetable.addClique`,
//!   HighsImplications.cpp:101–112), and
//! * a column tightened under **both** directions is tightened globally to
//!   the weaker of the two bounds: `x_j` is 0 or 1, so `x_k ≥
//!   min(lb↓, lb↑)` holds unconditionally (HiGHS `runProbing`'s up/down
//!   merge, :316–383; the substitution branch is not ported — no
//!   postsolve-level column removal at root scope).
//!
//! Effort is gated as in HiGHS (`HPresolve::runProbing`): a probe
//! contingent that starts at 1000 and grows by the reductions found
//! (:155, :1798–1814), an early abort when probing has already fixed too
//! many columns (:1769–1777), a cap on new cliques (:1782–1787), and the
//! solver's time limit (:1749–1759). Not ported: lifting opportunities and
//! variable-bound (VLB/VUB) bookkeeping (in-tree separation scope),
//! sub-MIP probing, and HiGHS's clique-implication-count candidate score
//! (our table has no per-literal implication index before probing;
//! candidates order by lock count then objective weight, after
//! Savelsbergh's probing-candidate criteria).
//!
//! No `unsafe`.

use crate::cuts::CliqueTable;
use crate::heuristics::Locks;
use crate::problem::MipProblem;
use crate::propagate::Propagator;
use std::time::Instant;

/// HiGHS's initial probe contingent (`probingContingent`, HPresolve.cpp:155).
pub const CONTINGENT: i64 = 1000;

/// Root-probing tunables (the constants a test or A/B moves).
#[derive(Clone, Copy)]
pub struct Params {
    pub feastol: f64,
    /// Solver time limit: probing stops when it passes (HiGHS checks
    /// `options->time_limit` per probe, HPresolve.cpp:1749–1759).
    pub deadline: Option<Instant>,
    /// Starting probe contingent (HPresolve.cpp:155); grows by the
    /// reductions each productive probe finds.
    pub contingent: i64,
    /// Objective cutoff for the probe propagations (`cᵀx ≤ cutoff`):
    /// `∞` without an incumbent; with one, a probe whose fix pushes the
    /// box's minimum objective past the cutoff proves the fixed value
    /// cannot belong to an improving point — the fixings HiGHS's root
    /// probing contributes through `separateImpliedBounds`
    /// (HighsSeparation.cpp:80–95 + HighsImplications.cpp:593–641), where
    /// probing runs per separation round with the domain's cutoff active.
    pub cutoff: f64,
}

/// Counters for the root-fixings measurement.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ProbingStats {
    /// Free binaries when probing started.
    pub candidates: u32,
    /// Probes actually run (both directions).
    pub probed: u32,
    /// Directions whose propagation was infeasible (each fixes the binary
    /// to the opposite value globally).
    pub infeasible_probes: u32,
    /// Columns probing fixed to a point (global box collapsed), cascade
    /// (clique + row propagation) included.
    pub fixings: u64,
    /// Non-fixing global bound tightenings.
    pub tightenings: u64,
    /// 2-cliques stored into the clique table from binary implications.
    pub cliques: u32,
    /// Probing emptied the global box: the problem is infeasible.
    pub infeasible: bool,
    /// Wall seconds spent probing (the time-box check).
    pub seconds: f64,
}

impl ProbingStats {
    /// Two-pass accumulation (pre-LP pass with cutoff ∞ + cut-loop-end pass
    /// with the incumbent cutoff): counters sum, the infeasible flag is
    /// sticky, the candidate census keeps its widest value (pass 2 sees only
    /// the binaries still free after pass 1's fixings).
    pub fn merge(self, later: Self) -> Self {
        Self {
            candidates: self.candidates.max(later.candidates),
            probed: self.probed + later.probed,
            infeasible_probes: self.infeasible_probes + later.infeasible_probes,
            fixings: self.fixings + later.fixings,
            tightenings: self.tightenings + later.tightenings,
            cliques: self.cliques + later.cliques,
            infeasible: self.infeasible || later.infeasible,
            seconds: self.seconds + later.seconds,
        }
    }
}

/// One direction's implications: how the propagated scratch box differs
/// from the global one (`j` excluded — it is the probe fix itself).
pub(crate) struct Direction {
    pub(crate) infeasible: bool,
    /// `(col, value)` tightened lower bounds under this direction.
    pub(crate) lower: Vec<(usize, f64)>,
    /// `(col, value)` tightened upper bounds under this direction.
    pub(crate) upper: Vec<(usize, f64)>,
    /// Free binaries the probe fixed to a value `(col, 0.0 | 1.0)`.
    pub(crate) fixed_binary: Vec<(usize, f64)>,
}

/// Runs root probing over the global box, tightening it in place. `cliques`
/// receives the binary implications as 2-cliques (the clique table when the
/// incumbent→fixings link is on). Every tightening holds for all
/// integer-feasible points, so the search may simply continue on the
/// tightened box; [`ProbingStats::infeasible`] means the box is empty.
pub(crate) fn run(
    p: &MipProblem,
    prop: &mut Propagator,
    locks: &Locks,
    mut cliques: Option<&mut CliqueTable>,
    lower: &mut [f64],
    upper: &mut [f64],
    params: Params,
) -> ProbingStats {
    let start = Instant::now();
    let ft = params.feastol;
    let mut stats = ProbingStats::default();
    // Integer bounds rounded inward (root_node's rule; HiGHS's domain keeps
    // integral bounds on integer columns), so binary fixings compare on
    // exact 0/1.
    for j in 0..p.num_col {
        if p.integer[j] {
            lower[j] = (lower[j] - ft).ceil();
            upper[j] = (upper[j] + ft).floor();
            if lower[j] > upper[j] {
                stats.infeasible = true;
                stats.seconds = start.elapsed().as_secs_f64();
                return stats;
            }
        }
    }
    // HiGHS `prepareProbing` propagates before the first probe
    // (HPresolve.cpp:1490).
    if !prop.propagate(lower, upper, None, params.cutoff) {
        stats.infeasible = true;
        stats.seconds = start.elapsed().as_secs_f64();
        return stats;
    }
    // Fixings this pass has produced (for the stats; a column counts once,
    // whatever channel fixed it — probe, merge, propagation or cascade).
    let mut fixed_seen: Vec<bool> = (0..p.num_col).map(|j| lower[j] == upper[j]).collect();
    // Candidates: free binaries, most locks first, then heaviest objective
    // coefficient (HiGHS scores by clique-table
    // implication counts, which need the probing this is seeding).
    let mut candidates: Vec<usize> = (0..p.num_col)
        .filter(|&j| p.integer[j] && lower[j] == 0.0 && upper[j] == 1.0)
        .collect();
    stats.candidates = candidates.len() as u32;
    candidates.sort_unstable_by(|&a, &b| {
        let key = |j: usize| (locks.up[j] + locks.down[j], p.cost[j].abs());
        key(b)
            .partial_cmp(&key(a))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });
    // HiGHS's abort bars (HPresolve.cpp:1769–1787, 1798).
    let max_fixings = 1000.max((p.num_row + p.num_col) / 20) as u64;
    let max_cliques = 1_000_000u64.max(2 * p.matrix.nnz() as u64);
    let mut contingent = params.contingent;
    let mut probed: i64 = 0;
    let (mut lo_s, mut hi_s) = (lower.to_vec(), upper.to_vec());
    for &j in &candidates {
        // Later probes may have fixed j; only free binaries are probed
        // (HiGHS re-checks `domain.isBinary`, HPresolve.cpp:1743).
        if lower[j] != 0.0 || upper[j] != 1.0 {
            continue;
        }
        if let Some(deadline) = params.deadline
            && Instant::now() >= deadline
        {
            break;
        }
        if contingent - probed < 0 {
            break;
        }
        if stats.fixings > max_fixings || stats.cliques as u64 > max_cliques {
            break;
        }
        let up = probe(
            j,
            true,
            prop,
            lower,
            upper,
            &mut lo_s,
            &mut hi_s,
            params.cutoff,
        );
        let down = probe(
            j,
            false,
            prop,
            lower,
            upper,
            &mut lo_s,
            &mut hi_s,
            params.cutoff,
        );
        probed += 1;
        stats.probed += 1;
        // An infeasible direction is the strongest outcome: `j = val` has
        // no feasible completion, so `j = 1 − val` in every feasible point
        // (HiGHS `vertexInfeasible` from `computeImplications` :54–60).
        let mut changed: Vec<usize> = Vec::new();
        for (dir, fix_val) in [(&up, 0.0), (&down, 1.0)] {
            if dir.infeasible {
                stats.infeasible_probes += 1;
                if fix_val == 0.0 {
                    upper[j] = 0.0;
                } else {
                    lower[j] = 1.0;
                }
                changed.push(j);
            }
        }
        if lower[j] > upper[j] {
            stats.infeasible = true;
            return stats;
        }
        // Binary implications → 2-cliques: `j = val ⇒ k = b` means the
        // literals `j = val` and `k = 1 − b` never coexist
        // (HighsImplications.cpp:101–112). Only sound when the direction
        // was feasible (an infeasible direction fixes more than cliques
        // express, and its diffs are the infeasibility trail).
        if let Some(table) = cliques.as_deref_mut() {
            let mut added = 0;
            for (dir, val) in [(&up, true), (&down, false)] {
                if dir.infeasible {
                    continue;
                }
                for &(k, b) in &dir.fixed_binary {
                    added += u32::from(table.add_probing_clique((j, val), (k, b <= 0.5)));
                }
            }
            stats.cliques += added;
        }
        if up.infeasible || down.infeasible {
            // The global fixing still needs propagation; the up/down merge
            // below is skipped (its inputs are one-sided).
            contingent += 1;
            if !apply(
                prop,
                cliques.as_deref_mut(),
                lower,
                upper,
                &changed,
                ft,
                params.cutoff,
                &mut stats,
                &mut fixed_seen,
            ) {
                stats.infeasible = true;
                stats.seconds = start.elapsed().as_secs_f64();
                return stats;
            }
            continue;
        }
        // The up/down merge: columns tightened under both directions take
        // the weaker bound globally (HighsImplications.cpp:316–383).
        let merged = merge(&up, &down, lower, upper, &mut stats);
        if merged.is_empty() {
            continue;
        }
        contingent += merged.len() as i64;
        if !apply(
            prop,
            cliques.as_deref_mut(),
            lower,
            upper,
            &merged,
            ft,
            params.cutoff,
            &mut stats,
            &mut fixed_seen,
        ) {
            stats.infeasible = true;
            stats.seconds = start.elapsed().as_secs_f64();
            return stats;
        }
    }
    stats.seconds = start.elapsed().as_secs_f64();
    stats
}

/// Probes `j = val` on a scratch copy of the global box: fix, propagate,
/// and return the box diffs as this direction's implications. `lo_s`/`hi_s`
/// are scratch buffers (reused across probes).
#[allow(clippy::too_many_arguments)] // one probe is its inputs plus the cutoff row
pub(crate) fn probe(
    j: usize,
    val: bool,
    prop: &mut Propagator,
    lower: &[f64],
    upper: &[f64],
    lo_s: &mut [f64],
    hi_s: &mut [f64],
    cutoff: f64,
) -> Direction {
    let mut dir = Direction {
        infeasible: false,
        lower: Vec::new(),
        upper: Vec::new(),
        fixed_binary: Vec::new(),
    };
    lo_s.copy_from_slice(lower);
    hi_s.copy_from_slice(upper);
    if val {
        lo_s[j] = 1.0;
    } else {
        hi_s[j] = 0.0;
    }
    if !prop.propagate(lo_s, hi_s, Some(j), cutoff) {
        dir.infeasible = true;
        return dir;
    }
    for k in 0..lower.len() {
        if k == j {
            continue;
        }
        let (lo, hi) = (lo_s[k], hi_s[k]);
        if lo > lower[k] {
            dir.lower.push((k, lo));
        }
        if hi < upper[k] {
            dir.upper.push((k, hi));
        }
        // A free binary the probe fixed (exact 0/1 by inward rounding).
        if lo == hi && (lo == 0.0 || lo == 1.0) && lower[k] == 0.0 && upper[k] == 1.0 {
            dir.fixed_binary.push((k, lo));
        }
    }
    dir
}

/// The up/down merge (HiGHS `runProbing` analysis :316–383): for each
/// column tightened under both directions, the weaker bound holds
/// unconditionally (`x_j` is 0 or 1). Applies the tightenings to the
/// global box and returns the columns that moved. The substitution branch
/// (a column fixed to *different* values under the two directions,
/// :356–365) is deliberately not ported; its effect is approximated by the
/// valid weaker-bound tightening.
fn merge(
    up: &Direction,
    down: &Direction,
    lower: &mut [f64],
    upper: &mut [f64],
    stats: &mut ProbingStats,
) -> Vec<usize> {
    let mut changed = Vec::new();
    let tighten = |side: &[(usize, f64)],
                   other: &[(usize, f64)],
                   bounds: &mut [f64],
                   weaker: fn(f64, f64) -> f64,
                   improves: fn(f64, f64) -> bool,
                   changed: &mut Vec<usize>| {
        for &(k, v) in side {
            let Some(&(_, v_other)) = other.iter().find(|&&(k2, _)| k2 == k) else {
                continue;
            };
            let w = weaker(v, v_other);
            if improves(w, bounds[k]) {
                bounds[k] = w;
                changed.push(k);
            }
        }
    };
    tighten(
        &up.lower,
        &down.lower,
        lower,
        f64::min,
        |w, b| w > b,
        &mut changed,
    );
    tighten(
        &up.upper,
        &down.upper,
        upper,
        f64::max,
        |w, b| w < b,
        &mut changed,
    );
    stats.tightenings += changed.iter().filter(|&&k| lower[k] != upper[k]).count() as u64;
    changed
}

/// Propagates a probe's global tightenings to a fixpoint with the row
/// propagator and the clique table's fixings cascade (the incumbent → fixings link —
/// probing's fixings put clique members at their literals, the cascade
/// fixes the rest; HiGHS gets the same interleaving from its event-driven
/// domain). Newly fixed columns are counted once into the stats via
/// `fixed_seen`. `false` = the global box is empty.
#[allow(clippy::too_many_arguments)] // one fixpoint step is its inputs plus its counters
fn apply(
    prop: &mut Propagator,
    mut cliques: Option<&mut CliqueTable>,
    lower: &mut [f64],
    upper: &mut [f64],
    changed: &[usize],
    feastol: f64,
    cutoff: f64,
    stats: &mut ProbingStats,
    fixed_seen: &mut [bool],
) -> bool {
    fn count_fixings(
        lower: &[f64],
        upper: &[f64],
        fixed_seen: &mut [bool],
        stats: &mut ProbingStats,
    ) {
        for j in 0..fixed_seen.len() {
            if !fixed_seen[j] && lower[j] == upper[j] {
                fixed_seen[j] = true;
                stats.fixings += 1;
            }
        }
    }
    let mut changed = changed.to_vec();
    loop {
        if !prop.propagate_cols(lower, upper, &changed, cutoff) {
            return false;
        }
        let Some(table) = cliques.as_deref_mut() else {
            count_fixings(lower, upper, fixed_seen, stats);
            return true;
        };
        let mut round = Vec::new();
        if !table.propagate_fixings(lower, upper, feastol, &mut round) {
            return false;
        }
        if round.is_empty() {
            count_fixings(lower, upper, fixed_seen, stats);
            return true;
        }
        changed = round;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_core::CscMatrix;

    fn problem(
        rows: &[(&[f64], f64, f64)],
        lower: &[f64],
        upper: &[f64],
        integer: &[bool],
    ) -> MipProblem {
        let n = lower.len();
        let mut trip = Vec::new();
        for (i, (coefs, _, _)) in rows.iter().enumerate() {
            for (j, &a) in coefs.iter().enumerate() {
                if a != 0.0 {
                    trip.push((i, j, a));
                }
            }
        }
        MipProblem {
            num_col: n,
            num_row: rows.len(),
            matrix: CscMatrix::from_triplets(rows.len(), n, &trip),
            row_lower: rows.iter().map(|r| r.1).collect(),
            row_upper: rows.iter().map(|r| r.2).collect(),
            cost: vec![0.0; n],
            maximize: false,
            col_lower: lower.to_vec(),
            col_upper: upper.to_vec(),
            integer: integer.to_vec(),
        }
    }

    fn probe_all(p: &MipProblem) -> (Vec<f64>, Vec<f64>, ProbingStats) {
        let mut prop = Propagator::new(p, 1e-6);
        let locks = Locks::new(p);
        let mut lo = p.col_lower.clone();
        let mut hi = p.col_upper.clone();
        let stats = run(
            p,
            &mut prop,
            &locks,
            None,
            &mut lo,
            &mut hi,
            Params {
                feastol: 1e-6,
                deadline: None,
                contingent: CONTINGENT,
                cutoff: f64::INFINITY,
            },
        );
        (lo, hi, stats)
    }

    #[test]
    fn infeasible_probe_fixes_opposite_value() {
        // The combinatorial case probing exists for: x0 = 1 forces
        // x1 = x2 = 0 (r0), which breaks r1 — a conflict no per-row bound
        // strengthening can see (x0 appears in no other row, so nothing
        // bounds it short of the hypothetical fix). r2 is slack and only
        // gives x0 the top lock score so it is probed first.
        let inf = f64::INFINITY;
        let p = problem(
            &[
                (&[1.0, 1.0, 1.0, 0.0], -inf, 1.0), // x0 + x1 + x2 ≤ 1
                (&[0.0, 1.0, 1.0, 0.0], 1.0, inf),  // x1 + x2 ≥ 1
                (&[1.0, 0.0, 0.0, 1.0], -5.0, inf), // x0 + t ≥ −5 (slack)
            ],
            &[0.0, 0.0, 0.0, 0.0],
            &[1.0, 1.0, 1.0, 10.0],
            &[true, true, true, false],
        );
        let (lo, hi, s) = probe_all(&p);
        assert!(!s.infeasible);
        assert_eq!((lo[0], hi[0]), (0.0, 0.0), "x0 = 0 via the up probe");
        assert_eq!(s.infeasible_probes, 1);
        assert_eq!(s.fixings, 1);
        // Both completions stay feasible: x1, x2 keep both values.
        assert_eq!((lo[1], hi[1]), (0.0, 1.0));
        assert_eq!((lo[2], hi[2]), (0.0, 1.0));
    }

    #[test]
    fn up_down_merge_tightens_weaker_bound() {
        // w ≥ x0 (F) and w ≥ 1 − x0 (D): per-row propagation leaves
        // w ≥ 0 (the max of two rows is not a row), but x0 = 1 forces
        // w ≥ 1 through F and x0 = 0 through D — both directions fix
        // w = 1, and E lifts it to k ≥ 2. The merge is the only rule
        // that gets here.
        let inf = f64::INFINITY;
        let p = problem(
            &[
                (&[-1.0, 1.0, 0.0, 0.0], 0.0, inf), // w − x0 ≥ 0
                (&[1.0, 1.0, 0.0, 0.0], 1.0, inf),  // w + x0 ≥ 1
                (&[0.0, -1.0, 1.0, 0.0], 1.0, inf), // k − w ≥ 1
                (&[1.0, 0.0, 0.0, 1.0], -inf, 5.0), // x0 − u ≤ 5 (slack)
            ],
            &[0.0, 0.0, 0.0, 0.0],
            &[1.0, 1.0, 5.0, 10.0],
            &[true, true, true, false],
        );
        let (lo, hi, s) = probe_all(&p);
        assert!(!s.infeasible);
        assert_eq!(
            s.probed, 1,
            "x0 probes first (top lock score), w is fixed by it"
        );
        assert_eq!((lo[1], hi[1]), (1.0, 1.0), "w = 1 via the merge");
        assert_eq!((lo[2], hi[2]), (2.0, 5.0), "k ≥ 2 via the merge");
        assert_eq!(s.fixings, 1);
        assert_eq!(s.tightenings, 1);
        assert_eq!(s.infeasible_probes, 0);
    }

    #[test]
    fn initial_propagation_empties_box_detects_infeasible() {
        // 5·x0 ≥ 3 (x0 ≥ 0.6 ⇒ 1) and 5·x0 ≤ 2 (x0 ≤ 0.4 ⇒ 0).
        let inf = f64::INFINITY;
        let p = problem(
            &[(&[5.0], 3.0, inf), (&[5.0], -inf, 2.0)],
            &[0.0],
            &[1.0],
            &[true],
        );
        let (_, _, s) = probe_all(&p);
        assert!(s.infeasible);
    }

    #[test]
    fn contingent_stops_unproductive_probes() {
        // Slack rows: neither direction of any probe tightens anything
        // (continuous acceptance rejects the 1-of-10 move), so the
        // contingent never grows and the pass stops after it is spent.
        let inf = f64::INFINITY;
        let n = 10;
        let rows: Vec<(Vec<f64>, f64, f64)> = (0..n)
            .map(|i| {
                let mut coefs = vec![0.0; 2 * n];
                coefs[i] = 1.0; // x_i
                coefs[n + i] = 1.0; // s_i ∈ [0, 10]
                (coefs, -inf, 10.0) // x_i + s_i ≤ 10
            })
            .collect();
        let lo = vec![0.0; 2 * n];
        let hi: Vec<f64> = (0..2 * n).map(|i| if i < n { 1.0 } else { 10.0 }).collect();
        let integer: Vec<bool> = (0..2 * n).map(|i| i < n).collect();
        let p = problem(
            &rows
                .iter()
                .map(|(c, l, u)| (c.as_slice(), *l, *u))
                .collect::<Vec<_>>(),
            &lo,
            &hi,
            &integer,
        );
        let mut prop = Propagator::new(&p, 1e-6);
        let locks = Locks::new(&p);
        let mut lo = p.col_lower.clone();
        let mut hi = p.col_upper.clone();
        let s = run(
            &p,
            &mut prop,
            &locks,
            None,
            &mut lo,
            &mut hi,
            Params {
                feastol: 1e-6,
                deadline: None,
                contingent: 3,
                cutoff: f64::INFINITY,
            },
        );
        assert!(!s.infeasible);
        assert_eq!(s.probed, 4, "the check is < 0, so contingent+1 probes");
        assert_eq!(s.fixings, 0);
        assert_eq!(s.tightenings, 0);
        // Untouched binaries stay free.
        assert_eq!((lo[5], hi[5]), (0.0, 1.0));
    }

    #[test]
    fn binary_implications_feed_the_clique_table() {
        // Same shape as `infeasible_probe_fixes_opposite_value`, run with
        // the clique table: probing's binary implications must land in
        // the table as 2-cliques (the seed measured as missing), and
        // the apply fixpoint must survive the cascade channel being live.
        let inf = f64::INFINITY;
        let p = problem(
            &[
                (&[1.0, 1.0, 1.0, 0.0], -inf, 1.0), // x0 + x1 + x2 ≤ 1
                (&[0.0, 1.0, 1.0, 0.0], 1.0, inf),  // x1 + x2 ≥ 1
                (&[1.0, 0.0, 0.0, 1.0], -5.0, inf), // x0 + t ≥ −5 (slack)
            ],
            &[0.0, 0.0, 0.0, 0.0],
            &[1.0, 1.0, 1.0, 10.0],
            &[true, true, true, false],
        );
        let mut prop = Propagator::new(&p, 1e-6);
        let locks = Locks::new(&p);
        // The model rows already contribute the triple clique; the feed
        // under test is probing's own 2-cliques on top.
        let mut table = crate::cuts::root_clique_table(&p, &p.col_lower, &p.col_upper, 1e-6);
        let before = table.num_cliques();
        let mut lo = p.col_lower.clone();
        let mut hi = p.col_upper.clone();
        let s = run(
            &p,
            &mut prop,
            &locks,
            Some(&mut table),
            &mut lo,
            &mut hi,
            Params {
                feastol: 1e-6,
                deadline: None,
                contingent: CONTINGENT,
                cutoff: f64::INFINITY,
            },
        );
        assert!(!s.infeasible);
        assert_eq!((lo[0], hi[0]), (0.0, 0.0));
        assert!(s.cliques > 0, "probing must add 2-cliques");
        assert!(table.num_cliques() > before);
    }

    #[test]
    fn cutoff_probing_fixes_improving_infeasible_directions() {
        // min 3·x0 + 5·x1 + 5·x2 over x0 + x1 ≥ 1 and x0 + x2 ≥ 1, cutoff
        // 9.5: only the down-probe conflicts — fixing x0 = 0 forces
        // x1 = x2 = 1 at objective 10 > 9.5, so no improving point has
        // x0 = 0 and the probe fixes x0 = 1. The cutoff window is load-
        // bearing: below 8 the objective row (5·x1 ≤ cutoff − 3) rounds x1
        // and x2 to 0 in the *initial* propagation and everything is fixed
        // before any probe runs (the 2-var version of this test died that
        // way); at ≥ 10 the down-probe stays feasible. This is the pass the
        // root actually needs (HiGHS `separateImpliedBounds` probes with the
        // domain cutoff active, HighsImplications.cpp:583–641); with cutoff
        // ∞ the same model probes to nothing.
        let inf = f64::INFINITY;
        let mut p = problem(
            &[
                (&[1.0, 1.0, 0.0], 1.0, inf), // x0 + x1 ≥ 1
                (&[1.0, 0.0, 1.0], 1.0, inf), // x0 + x2 ≥ 1
            ],
            &[0.0, 0.0, 0.0],
            &[1.0, 1.0, 1.0],
            &[true, true, true],
        );
        p.cost = vec![3.0, 5.0, 5.0];
        let mut prop = Propagator::new(&p, 1e-6);
        let locks = Locks::new(&p);
        let mut lo = p.col_lower.clone();
        let mut hi = p.col_upper.clone();
        let s = run(
            &p,
            &mut prop,
            &locks,
            None,
            &mut lo,
            &mut hi,
            Params {
                feastol: 1e-6,
                deadline: None,
                contingent: CONTINGENT,
                cutoff: 9.5,
            },
        );
        assert!(!s.infeasible);
        assert_eq!((lo[0], hi[0]), (1.0, 1.0));
        // x1/x2 stay free: every other direction (x0 = 1, or probing x1/x2)
        // leaves objective ≤ 9.5 reachable.
        assert_eq!((lo[1], hi[1]), (0.0, 1.0));
        assert_eq!((lo[2], hi[2]), (0.0, 1.0));
        assert_eq!(s.infeasible_probes, 1);
        assert_eq!(s.fixings, 1);
    }
}
