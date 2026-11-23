//! Gomory mixed-integer cuts from simplex tableau rows.
//!
//! Follows HiGHS `HighsTableauSeparator::separateLpSolution`
//! (reference/HiGHS/highs/mip/HighsTableauSeparator.cpp:39–244): pick the
//! basic integer columns with a fractional value (:62–82), score them by
//! `f(1−f)/‖scaled tableau row weights‖²` (:143–175), and for the best ones
//! aggregate the LP rows with the tableau row's weights `y = e_posᵀB⁻¹`
//! (:200–203, `HighsLpAggregator::addRow`, HighsLpAggregator.cpp:19–28) and
//! hand the aggregation, and its negation, to the c-MIR generator
//! (:231–238). The aggregation `Σ y_r (a_rᵀx − s_r) = 0` over columns and row
//! slacks IS the tableau row with the basic variables eliminated; the
//! generator substitutes the nonbasic variables at their bounds and rounds,
//! which is the Gomory mixed-integer rounding of Gomory (1960) in the
//! mixed-integer-rounding form of Marchand & Wolsey, "Aggregation and mixed
//! integer rounding to solve MIPs", Oper. Res. 49 (2001), and is the
//! combination Balas, Ceria, Cornuéjols & Natraj, "Gomory cuts revisited",
//! Oper. Res. Lett. 19 (1996) argue for. HiGHS is MIT licensed (Copyright (c)
//! 2016–2024 the HiGHS authors); it was read as a reference text and
//! reimplemented.
//!
//! Validity does not depend on the weights: any combination `Σ w_r (a_rᵀx −
//! s_r) = 0` is an identity over the row activities, and the generator
//! derives its cut from the global bounds of columns and slacks, so a
//! numerically poor tableau row can only give a weak cut, never a wrong one.
//!
//! Deviations: HiGHS also separates on basic integral row slacks
//! (:66–71); this does not, because those positions are mostly one-weight
//! rows another separator covers. HiGHS orders candidates beyond its try
//! budget by dual edge weights when it has them (:99–118) and by
//! fractionality otherwise (:119–133); there are no edge weights here.

use super::generate::{CutGen, GeneratedCut};
use super::relax::{Aggregator, Relaxation};
use super::transform::TransformedLp;
use crate::engine::{LpEngine, LpOutcome};

/// Candidates with fractionality below `FRACTIONALITY_FACTOR · feastol` are
/// skipped (HighsTableauSeparator.cpp:79).
const FRACTIONALITY_FACTOR: f64 = 1000.0;

/// Tableau rows whose scaled weights span more than this are numerically
/// unreliable (HighsTableauSeparator.cpp:171).
const MAX_WEIGHT_RATIO: f64 = 1e4;

/// Aggregations whose column coefficients span more than this are skipped
/// (HighsTableauSeparator.cpp:221).
const MAX_COEFFICIENT_RATIO: f64 = 1e6;

/// Score thresholds relative to the best cut-producing score, before and
/// after 50 cuts (HighsTableauSeparator.cpp:190, :195–196).
const SCORE_FACTOR: [f64; 2] = [0.0025, 0.01];

/// Cuts per call (HighsTableauSeparator.cpp:193).
const MAX_CUTS: usize = 1000;

/// How many candidate rows the separator has tried over the run, which
/// budgets later calls (HighsTableauSeparator.cpp:85–97).
#[derive(Debug, Default)]
pub(crate) struct TableauState {
    num_tries: u64,
    num_calls: u64,
    /// Parked for the rest of the root loop: a round ran the separator and
    /// the pool selected none of its cuts, so re-running only pays the
    /// basis scan. Our own rule (HiGHS
    /// re-separates every round but its `getBasisInverseRowSparse` is one
    /// B-solve per candidate, while our weights come from a full
    /// basis-position scan): measured 2026-10-23, medium tier,
    /// set_partition_s1 runs 19 root rounds where Gomory generates nothing
    /// (0 cuts) and the cut loop still quadruples (0.004 s → 0.015 s);
    /// graph_coloring 9 generated 0 applied; gap_s3 69 generated 0 applied.
    pub(crate) parked: bool,
}

struct Candidate {
    col: usize,
    score: f64,
    weights: Vec<(usize, f64)>,
}

/// Longest cut the separator lets into the pool: `max(density · average model
/// row length, 0.1 · numCol)`. This is the rule `HighsCutPool::addCut` applies
/// to decide whether a cut may be used for propagation (alpha = 2,
/// HighsCutPool.cpp:547–563, `minDensityLim = 0.1 * ncols` at
/// HighsCutPool.h:97): `Rlen >= max(alpha * avgModelNzs, minDensityLim)` means
/// dense. HiGHS only withholds dense cuts from propagation and lets them into
/// the LP; here a dense Gomory cut stays out of the pool, because it costs
/// every node LP. `density <= 0` disables the cap.
fn max_cut_len(relax: &Relaxation, density: f64) -> usize {
    if density <= 0.0 {
        return usize::MAX;
    }
    let rows = &relax.rows[..relax.num_model_rows];
    let nnz: usize = rows.iter().map(|r| r.vals.len()).sum();
    let avg = nnz as f64 / rows.len().max(1) as f64;
    (density * avg).max(0.1 * relax.num_col as f64).ceil() as usize
}

/// Distance of `v` to the nearest integer.
fn fractionality(v: f64) -> f64 {
    (v - (v + 0.5).floor()).abs()
}

/// Separates Gomory cuts at `lp.x` (= `relax.x`); `add` receives each
/// generated cut and returns whether the pool accepted it.
#[allow(clippy::too_many_arguments)]
pub(crate) fn separate<E: LpEngine>(
    engine: &mut E,
    lower: &[f64],
    upper: &[f64],
    lp: &LpOutcome,
    relax: &Relaxation,
    trans: &mut TransformedLp,
    agg: &mut Aggregator,
    cutgen: &mut CutGen,
    state: &mut TableauState,
    total_lp_iterations: u64,
    feastol: f64,
    epsilon: f64,
    density: f64,
    mut add: impl FnMut(GeneratedCut) -> bool,
) {
    if state.parked {
        return;
    }
    let n = relax.num_col;
    let max_len = max_cut_len(relax, density);
    agg.ensure_dim(n + relax.rows.len());

    let mut frac: Vec<(usize, f64)> = (0..n)
        .filter(|&j| relax.integer[j])
        .map(|j| (j, fractionality(lp.x[j])))
        .filter(|&(_, f)| f >= FRACTIONALITY_FACTOR * feastol)
        .collect();
    if frac.is_empty() {
        return;
    }
    let calls = state.num_calls;
    let max_tries = 5000 + calls * 50 + total_lp_iterations / 10;
    if state.num_tries >= max_tries {
        return;
    }
    let num_int = relax.integer.iter().filter(|&&b| b).count();
    let per_call = 200 + (0.1 * relax.rows.len().min(num_int) as f64) as u64;
    let budget = (max_tries - state.num_tries).min(per_call) as usize;
    if frac.len() > budget {
        // Most fractional first, index as the tie-break (HiGHS breaks ties
        // by hash; any fixed order will do).
        frac.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        frac.truncate(budget);
    }
    state.num_calls += 1;
    state.num_tries += frac.len() as u64;

    let cols: Vec<usize> = frac.iter().map(|&(j, _)| j).collect();
    let Some(rows) = engine.tableau_weights(lp, lower, upper, &cols) else {
        return;
    };

    // Score by fractionality and the scaled weights
    // (HighsTableauSeparator.cpp:149–175).
    let max_abs: Vec<f64> = relax
        .rows
        .iter()
        .map(|r| r.vals.iter().fold(0.0f64, |m, v| m.max(v.abs())))
        .collect();
    let mut cands: Vec<Candidate> = Vec::new();
    for row in rows {
        // A single weight is a row of the LP itself, which the other
        // separators handle (:149).
        if row.weights.len() <= 1 {
            continue;
        }
        let f = fractionality(lp.x[row.col]);
        let (mut lo, mut hi, mut norm2) = (f64::INFINITY, 0.0f64, 0.0);
        let mut weights = Vec::with_capacity(row.weights.len());
        for &(r, w) in &row.weights {
            let scaled = max_abs[r] * w.abs();
            if scaled <= feastol {
                continue;
            }
            lo = lo.min(scaled);
            hi = hi.max(scaled);
            norm2 += scaled * scaled;
            weights.push((r, w));
        }
        if weights.len() <= 1 || hi / lo > MAX_WEIGHT_RATIO {
            continue;
        }
        let score = f * (1.0 - f) / norm2;
        if score > feastol {
            cands.push(Candidate {
                col: row.col,
                score,
                weights,
            });
        }
    }
    if cands.is_empty() {
        return;
    }
    cands.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.col.cmp(&b.col)));

    let mut best_score = -1.0;
    let mut cuts = 0usize;
    for cand in &cands {
        if cuts >= MAX_CUTS || cand.score < SCORE_FACTOR[usize::from(cuts >= 50)] * best_score {
            break;
        }
        agg.clear();
        for &(r, w) in &cand.weights {
            agg.add_row(relax, r, w);
        }
        let (inds, vals) = agg.current(n, epsilon, false);
        let len = inds.len();
        if 10 * (len - cand.weights.len()) > 10_000 + n {
            agg.clear();
            continue;
        }
        if len > cand.weights.len() {
            let (mut mx, mut mn) = (0.0f64, f64::INFINITY);
            for (&k, &v) in inds.iter().zip(&vals) {
                if k < n {
                    mx = mx.max(v.abs());
                    mn = mn.min(v.abs());
                }
            }
            if mx / mn > MAX_COEFFICIENT_RATIO {
                agg.clear();
                continue;
            }
        }
        let mut found = false;
        if let Some(cut) = cutgen.generate(relax, trans, inds, vals, 0.0, false) {
            let ok = cut.inds.len() <= max_len && add(cut);
            cuts += usize::from(ok);
            found |= ok;
        }
        let (inds, vals) = agg.current(n, epsilon, true);
        if let Some(cut) = cutgen.generate(relax, trans, inds, vals, 0.0, false) {
            let ok = cut.inds.len() <= max_len && add(cut);
            cuts += usize::from(ok);
            found |= ok;
        }
        agg.clear();
        if best_score == -1.0 && found {
            best_score = cand.score;
        }
    }
}
