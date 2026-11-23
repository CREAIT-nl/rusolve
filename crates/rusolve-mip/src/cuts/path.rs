//! Path aggregation: starting from each tight LP row, continuous columns that
//! sit strictly inside their bounds are eliminated by adding a row they
//! appear in with the opposite sign, up to 6 rows; after every step the
//! aggregation (and its negation) is handed to [`super::generate::CutGen`]
//! (lifted covers / c-MIR). With variable-bound substitution in
//! [`super::transform`] this produces flow-cover-type cuts, so there is no
//! separate flow-cover separator (as in HiGHS).
//!
//! Source: Marchand & Wolsey, "Aggregation and mixed integer rounding to
//! solve MIPs", Operations Research 49(3):363–371 (2001), §4 (aggregation
//! heuristic: eliminate the continuous variable farthest from its bounds,
//! path length ≤ 6). Arc bookkeeping, row typing, single-continuous equation
//! substitution, the 'try the negated start' rule and the weight window
//! `[feastol, 1/feastol]` follow HiGHS `HighsPathSeparator::separateLpSolution`
//! (reference/HiGHS/highs/mip/HighsPathSeparator.cpp:26–550, MIT; read as a
//! reference text and reimplemented) and `HighsLpAggregator`
//! (HighsLpAggregator.cpp:19–52).
//!
//! Not ported: path-mixing cuts (HighsPathSeparator.cpp:390–540) and the
//! dual-sign choice of the start direction for equations (no LP duals yet:
//! equations always start with weight +1, then −1).

use super::generate::{CutGen, GeneratedCut};
use super::relax::{Aggregator, Relaxation};
use super::transform::{SMALL, TransformedLp};

const MAX_PATH_LEN: usize = 6;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RowType {
    Unusable,
    Geq,
    Eq,
    Leq,
}

fn next_rand(state: &mut u64) -> u64 {
    // xorshift64*
    let mut x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    x.wrapping_mul(0x2545_f491_4f6c_dd1d)
}

/// Runs the path separator on the current point of `relax`. `add` receives
/// every generated cut with its path length and returns whether the pool
/// accepted it.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub(crate) fn separate(
    relax: &Relaxation,
    trans: &mut TransformedLp,
    agg: &mut Aggregator,
    cutgen: &mut CutGen,
    rng: &mut u64,
    feastol: f64,
    mut add: impl FnMut(GeneratedCut, usize) -> bool,
) {
    let n = relax.num_col;
    let m = relax.rows.len();
    agg.ensure_dim(n + m);
    let mut rowtype = vec![RowType::Unusable; m];
    for (i, r) in relax.rows.iter().enumerate() {
        if r.row_lower == r.row_upper {
            rowtype[i] = RowType::Eq;
            continue;
        }
        let a = relax.activity[i];
        let lowerslack = if r.row_lower.is_finite() {
            a - r.row_lower
        } else {
            f64::INFINITY
        };
        let upperslack = if r.row_upper.is_finite() {
            r.row_upper - a
        } else {
            f64::INFINITY
        };
        rowtype[i] = if lowerslack > feastol && upperslack > feastol {
            RowType::Unusable
        } else if lowerslack < upperslack {
            RowType::Geq
        } else {
            RowType::Leq
        };
    }

    // Column view of the rows, restricted to continuous columns off their
    // bounds (the only ones a path eliminates).
    let bound_dist: Vec<f64> = (0..n).map(|j| trans.bound_distance(j)).collect();
    let active = |j: usize| !relax.integer[j] && bound_dist[j] > 0.0;
    let mut col_rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut num_cont = vec![0usize; m];
    for (i, r) in relax.rows.iter().enumerate() {
        for (&j, &v) in r.inds.iter().zip(&r.vals) {
            if active(j) {
                col_rows[j].push((i, v));
                num_cont[i] += 1;
            }
        }
    }

    // Equations with a single active continuous column substitute it away.
    let mut subst: Vec<Option<(usize, f64)>> = vec![None; n];
    for i in 0..m {
        if rowtype[i] != RowType::Eq || num_cont[i] != 1 {
            continue;
        }
        let r = &relax.rows[i];
        let Some((j, v)) = r
            .inds
            .iter()
            .zip(&r.vals)
            .find(|(j, _)| active(**j))
            .map(|(&j, &v)| (j, v))
        else {
            continue;
        };
        if subst[j].is_some() {
            continue;
        }
        subst[j] = Some((i, v));
        rowtype[i] = RowType::Unusable;
    }

    // In-arcs / out-arcs: ≤ rows with a negative (≥ rows: positive)
    // coefficient are in-arcs, the opposite sign out-arcs; equations both.
    let mut in_arcs: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut out_arcs: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for j in 0..n {
        if subst[j].is_some() {
            continue;
        }
        for &(i, v) in &col_rows[j] {
            match rowtype[i] {
                RowType::Unusable => {}
                RowType::Leq => {
                    if v < 0.0 {
                        in_arcs[j].push((i, v));
                    } else {
                        out_arcs[j].push((i, v));
                    }
                }
                RowType::Geq => {
                    if v > 0.0 {
                        in_arcs[j].push((i, v));
                    } else {
                        out_arcs[j].push((i, v));
                    }
                }
                RowType::Eq => {
                    in_arcs[j].push((i, v));
                    out_arcs[j].push((i, v));
                }
            }
        }
    }

    let max_weight = 1.0 / feastol;
    let min_weight = feastol;
    let check_weight = |w: f64| {
        let w = w.abs();
        w <= max_weight && w >= min_weight
    };

    for start in 0..m {
        let scales = match rowtype[start] {
            RowType::Unusable => continue,
            RowType::Eq | RowType::Leq => [1.0, -1.0],
            RowType::Geq => [-1.0, 1.0],
        };
        for &scale in &scales {
            agg.clear();
            agg.add_row(relax, start, scale);
            let mut path: Vec<usize> = vec![start];
            let mut try_negated = false;
            let mut subst_rounds = 0;

            while path.len() != MAX_PATH_LEN {
                let (inds, vals) = agg.current(n, SMALL, false);
                let mut added_subst = false;
                let mut best_out: Option<(usize, f64, f64)> = None;
                let mut best_in: Option<(usize, f64, f64)> = None;
                for (&j, &v) in inds.iter().zip(&vals) {
                    if j >= n || !active(j) {
                        continue;
                    }
                    if let Some((row, sv)) = subst[j] {
                        added_subst = true;
                        agg.add_row(relax, row, -v / sv);
                        continue;
                    }
                    if added_subst {
                        continue;
                    }
                    // HighsPathSeparator skipCol: note whether the start row
                    // could also be extended in the negated direction, and
                    // skip columns with no arc leaving the current path.
                    let (arcs, other) = if v < 0.0 {
                        (&out_arcs[j], &in_arcs[j])
                    } else {
                        (&in_arcs[j], &out_arcs[j])
                    };
                    if path.len() == 1
                        && !try_negated
                        && (arcs.len() > path.len() || arcs.iter().any(|a| a.0 != start))
                    {
                        try_negated = true;
                    }
                    let skip = other.is_empty()
                        || (other.len() <= path.len() && other.iter().all(|a| path.contains(&a.0)));
                    if skip {
                        continue;
                    }
                    let bd = bound_dist[j];
                    let best = if v < 0.0 { &mut best_out } else { &mut best_in };
                    if best.is_none_or(|b| bd > b.2) {
                        *best = Some((j, v, bd));
                    }
                }
                if added_subst {
                    subst_rounds += 1;
                    if subst_rounds > MAX_PATH_LEN * 4 {
                        break;
                    }
                    continue;
                }

                let mut success = false;
                let len = path.len();
                if let Some(cut) =
                    cutgen.generate(relax, trans, inds.clone(), vals.clone(), 0.0, false)
                {
                    success |= add(cut, len);
                }
                let neg: Vec<f64> = vals.iter().map(|v| -v).collect();
                if let Some(cut) = cutgen.generate(relax, trans, inds, neg, 0.0, false) {
                    success |= add(cut, len);
                }
                if success || (best_out.is_none() && best_in.is_none()) {
                    break;
                }

                // Prefer an out-arc column unless the in-arc one is farther
                // from its bounds by more than feastol.
                let find_row = |val: f64, arcs: &[(usize, f64)], rng: &mut u64| {
                    if arcs.is_empty() {
                        return None;
                    }
                    let s = (next_rand(rng) % arcs.len() as u64) as usize;
                    (s..arcs.len()).chain(0..s).find_map(|k| {
                        let (r, a) = arcs[k];
                        let w = -val / a;
                        (!path.contains(&r) && check_weight(w)).then_some((r, w))
                    })
                };
                let pick = match (best_out, best_in) {
                    (Some(o), i) if i.is_none_or(|i| o.2 >= i.2 - feastol) => {
                        find_row(o.1, &in_arcs[o.0], rng)
                            .or_else(|| i.and_then(|i| find_row(i.1, &out_arcs[i.0], rng)))
                    }
                    (_, Some(i)) => find_row(i.1, &out_arcs[i.0], rng),
                    _ => None,
                };
                let Some((row, weight)) = pick else {
                    break;
                };
                path.push(row);
                agg.add_row(relax, row, weight);
            }
            agg.clear();
            if !try_negated {
                break;
            }
        }
    }
}
