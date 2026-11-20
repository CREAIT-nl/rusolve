//! SCIP linear-row knapsack separation (rusolve lever `mip_knapsack_cover_sepa`, default off).
//!
//! Translated from SCIP 10 `src/scip/cons_knapsack.c` (Apache License 2.0, Copyright (c)
//! 2002-2025 Zuse Institute Berlin; notice in THIRD_PARTY_NOTICES.md): `SCIPseparateRelaxedKnapsack`
//! (relaxation of a general row to a 0-1 knapsack, variable-bound substitution),
//! `SCIPseparateKnapsackCuts` (LMCI1 and LEWI; GUBs and LMCI2 are off by default in SCIP and
//! not ported), `getCover`, `makeCoverMinimal`, `getFeasibleSet`, `getPartition*`,
//! `getLiftingSequence`, `sequentialUpAndDownLifting`, `separateSequLiftedMinimalCoverInequality`
//! (cut name `_mcseq`) and `separateSequLiftedExtendedWeightInequality` (`_ewseq`);
//! `SCIPsolveKnapsackApproximately` (Dantzig rounding) from the same file. Row selection follows
//! `cons_linear.c::separateCons`: only rows with a nonzero LP dual, on the side the dual
//! indicates. Rust port, not a line-by-line copy:
//! SCIP's `SCIPselectWeightedDownRealLongRealInt` partial sort is a full descending sort here.

use crate::prelude::*;

/// cons_knapsack.c `MAXABSVBCOEF`
const MAXABSVBCOEF: f64 = 1e5;
/// cons_knapsack.c `KNAPSACKRELAX_MAXSCALE`
const KNAPSACKRELAX_MAXSCALE: i64 = 1000;
/// cons_knapsack.c `MAXCOVERSIZEITERLEWI`
const MAXCOVERSIZEITERLEWI: usize = 1000;
/// cons_linear.c `DEFAULT_MAXSEPACUTSROOT`
const MAXSEPACUTS: usize = 200;
/// SCIP `SCIP_DEFAULT_SEPA_MINEFFICACY`
const MIN_EFFICACY: f64 = 1e-4;
/// Rows longer than this are skipped (cost guard, not in SCIP).
const MAX_ROW_LEN: usize = 5_000;

#[derive(Clone, Default, Debug)]
pub struct HighsKnapsackCoverSeparator {
    pub base: HighsSeparator,
}

/// Relaxed knapsack `sum w_i y_i <= cap`, `y_i = x_col` or `1 - x_col` (`neg`), `val` = LP value of `y`.
struct Knap {
    col: Vec<usize>,
    neg: Vec<bool>,
    w: Vec<i64>,
    val: Vec<f64>,
    cap: i64,
}

type KCut = (Vec<(usize, i64)>, i64);

impl HighsKnapsackCoverSeparator {
    pub fn new(ms: &HighsMipSolver) -> Self {
        HighsKnapsackCoverSeparator {
            base: HighsSeparator::new(ms, "Separation: Knapsack cover"),
        }
    }

    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        _lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        let g = trans_lp.globaldom_;
        let nrow = ms.num_row();
        let sol: Vec<f64> = lp_relaxation.get_solution().col_value.clone();
        let dual: Vec<f64> = lp_relaxation.get_solution().row_dual.clone();
        let mut cuts: Vec<(Vec<HighsInt>, Vec<f64>, f64)> = Vec::new();
        {
            let lower = &ms.mipdata_.domains[g].col_lower_;
            let upper = &ms.mipdata_.domains[g].col_upper_;
            let lp = lp_relaxation.get_lp();
            'rows: for r in 0..nrow {
                if cuts.len() >= MAXSEPACUTS {
                    break 'rows;
                }
                let d = dual.get(r as usize).copied().unwrap_or(0.0);
                if d.abs() <= 1e-6 {
                    continue;
                }
                let mut len: HighsInt = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                lp_relaxation.get_row(ms, r, &mut len, &mut inds, &mut vals);
                if len < 2 || len as usize > MAX_ROW_LEN {
                    continue;
                }
                let (sign, rhs) = if d < 0.0 {
                    (1.0, lp.row_upper_[r as usize])
                } else {
                    (-1.0, -lp.row_lower_[r as usize])
                };
                if !rhs.is_finite() || rhs.abs() >= kHighsInf {
                    continue;
                }
                let Some(kn) = relax_row(ms, inds, vals, sign, rhs, lower, upper, &sol) else {
                    continue;
                };
                for (terms, rhs_k) in separate_knapsack(&kn) {
                    // y-space `sum a y <= rhs_k` to x-space
                    let mut oi: Vec<HighsInt> = Vec::new();
                    let mut ov: Vec<f64> = Vec::new();
                    let mut orhs = rhs_k as f64;
                    for (i, a) in terms {
                        oi.push(kn.col[i] as HighsInt);
                        if kn.neg[i] {
                            ov.push(-(a as f64));
                            orhs -= a as f64;
                        } else {
                            ov.push(a as f64);
                        }
                    }
                    let (mut act, mut nrm) = (0.0, 0.0);
                    for (j, &cj) in oi.iter().enumerate() {
                        act += ov[j] * sol[cj as usize];
                        nrm += ov[j] * ov[j];
                    }
                    if nrm > 0.0 && (act - orhs) / nrm.sqrt() > MIN_EFFICACY {
                        cuts.push((oi, ov, orhs));
                    }
                }
            }
        }
        for (mut inds, mut vals, rhs) in cuts {
            let len = inds.len() as HighsInt;
            detach!(ms.mipdata_.cutpools[cutpool], |p| p.add_cut(
                ms,
                &mut inds[..],
                &mut vals[..],
                len,
                rhs,
                true,
                true,
                false,
                false,
            ));
        }
    }
}

/// cons_knapsack.c `SCIPseparateRelaxedKnapsack` (relaxation part): `sign * a x <= rhs`.
#[allow(clippy::too_many_arguments)]
fn relax_row(
    ms: &HighsMipSolver,
    inds: &[HighsInt],
    vals: &[f64],
    sign: f64,
    rhs: f64,
    lower: &[f64],
    upper: &[f64],
    sol: &[f64],
) -> Option<Knap> {
    let is_bin = |cc: usize| {
        !matches!(
            ms.variable_type_2(cc as HighsInt),
            HighsVarType::kContinuous
        ) && lower[cc] == 0.0
            && upper[cc] == 1.0
    };
    let imps = &ms.mipdata_.implications;
    let mut rhs = rhs;
    let mut binvals: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
    for (k, &cj) in inds.iter().enumerate() {
        let cc = cj as usize;
        let a = sign * vals[k];
        if is_bin(cc) {
            if sol[cc] < -1e-6 || sol[cc] > 1.0 + 1e-6 {
                return None;
            }
            *binvals.entry(cc).or_insert(0.0) += a;
        } else if a > 0.0 {
            // a > 0: substitute the lower bound or the variable lower bound with largest LP value
            let mut best = lower[cc];
            let mut bt: Option<(usize, f64, f64)> = None;
            let mut lst: Vec<(HighsInt, f64, f64)> = Vec::new();
            imps.vlbs[cc].for_each(&mut |&z, vb| lst.push((z, vb.coef, vb.constant)));
            for (z, b, d) in lst {
                let zc = z as usize;
                if !is_bin(zc) || b.abs() > MAXABSVBCOEF {
                    continue;
                }
                let v = b * sol[zc] + d;
                if v >= best - 1e-9 {
                    best = v;
                    bt = Some((zc, b, d));
                }
            }
            if best <= -kHighsInf {
                return None;
            }
            match bt {
                None => rhs -= a * best,
                Some((zc, b, d)) => {
                    rhs -= a * d;
                    *binvals.entry(zc).or_insert(0.0) += a * b;
                }
            }
        } else {
            let mut best = upper[cc];
            let mut bt: Option<(usize, f64, f64)> = None;
            let mut lst: Vec<(HighsInt, f64, f64)> = Vec::new();
            imps.vubs[cc].for_each(&mut |&z, vb| lst.push((z, vb.coef, vb.constant)));
            for (z, b, d) in lst {
                let zc = z as usize;
                if !is_bin(zc) || b.abs() > MAXABSVBCOEF {
                    continue;
                }
                let v = b * sol[zc] + d;
                if v <= best + 1e-9 {
                    best = v;
                    bt = Some((zc, b, d));
                }
            }
            if best >= kHighsInf {
                return None;
            }
            match bt {
                None => rhs -= a * best,
                Some((zc, b, d)) => {
                    rhs -= a * d;
                    *binvals.entry(zc).or_insert(0.0) += a * b;
                }
            }
        }
    }
    // integral scalar (SCIPcalcIntegralScalar approximated by a search over 1..MAXSCALE)
    let bv: Vec<f64> = binvals
        .values()
        .copied()
        .filter(|v| v.abs() > 1e-9)
        .collect();
    let mut scalar = 1.0;
    for s in 1..=KNAPSACKRELAX_MAXSCALE {
        let sf = s as f64;
        if bv.iter().all(|v| (v * sf - (v * sf).round()).abs() <= 1e-6) {
            scalar = sf;
            break;
        }
    }
    rhs *= scalar;
    let mut kn = Knap {
        col: Vec::new(),
        neg: Vec::new(),
        w: Vec::new(),
        val: Vec::new(),
        cap: 0,
    };
    let (mut minact, mut maxact) = (0i64, 0i64);
    for (&cc, &bvv) in &binvals {
        let mut v = (bvv * scalar + 1e-9).floor() as i64;
        if v == 0 {
            continue;
        }
        let neg = v < 0;
        if neg {
            v = -v;
            rhs += v as f64;
        }
        let y = if neg { 1.0 - sol[cc] } else { sol[cc] };
        // global bounds [0,1]: y can be 0 or 1
        maxact += v;
        kn.col.push(cc);
        kn.neg.push(neg);
        kn.w.push(v);
        kn.val.push(y);
    }
    let _ = &mut minact;
    if kn.w.len() < 2 || !rhs.is_finite() || rhs.abs() > 1e12 {
        return None;
    }
    kn.cap = (rhs + 1e-6).floor() as i64;
    if kn.cap < 0 || maxact <= kn.cap {
        return None; // cutoff (minact > cap) or the knapsack is redundant
    }
    Some(kn)
}

/// cons_knapsack.c `SCIPseparateKnapsackCuts`: LMCI1 then LEWI.
fn separate_knapsack(kn: &Knap) -> Vec<KCut> {
    let n = kn.w.len();
    let mut out: Vec<KCut> = Vec::new();
    // LMCI1: modified transformed separation problem
    if let Some(mut cv) = get_cover(kn, true) {
        if !cv.fractional {
            return out;
        }
        make_cover_minimal(kn, &mut cv);
        if let Some(kc) = sequ_lifted_mincover(kn, &cv.cover, &cv.noncover) {
            out.push(kc);
        }
    }
    // LEWI: transformed separation problem
    if let Some(mut cv) = get_cover(kn, false) {
        // sort by (x-1)/a ascending (SCIPsortRealInt)
        let key = |j: usize| (kn.val[j] - 1.0) / kn.w[j] as f64;
        cv.cover
            .sort_by(|&p, &q| key(p).partial_cmp(&key(q)).unwrap().then(p.cmp(&q)));
        let mut cover = cv.cover.clone();
        let mut non = cv.noncover.clone();
        let mut cw = cv.weight;
        while cover.len() >= 2 {
            let first = cover.remove(0);
            non.push(first);
            cw -= kn.w[first];
            if cw <= kn.cap {
                if let Some(kc) = sequ_lifted_ewi(kn, &cover, &non) {
                    out.push(kc);
                }
            }
            if cover.len() >= MAXCOVERSIZEITERLEWI {
                break;
            }
        }
    }
    let _ = n;
    out
}

struct Cover {
    cover: Vec<usize>,
    noncover: Vec<usize>,
    weight: i64,
    fractional: bool,
    #[allow(dead_code)]
    ntightened: usize,
}

fn feas_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-6
}

/// cons_knapsack.c `getCover`. `None` = no cover found (or not fractional: `fractional=false`
/// is returned as `Some` with empty sets so the caller can stop LMCI1).
fn get_cover(kn: &Knap, modtrans: bool) -> Option<Cover> {
    let n = kn.w.len();
    let mut ones = Vec::new();
    let mut zeros = Vec::new();
    let mut items = Vec::new();
    let (mut onesw, mut itemsw) = (0i64, 0i64);
    let mut nt = 0usize;
    for j in 0..n {
        if kn.w[j] > kn.cap {
            nt += 1;
            continue;
        }
        if feas_eq(kn.val[j], 1.0) {
            ones.push(j);
            onesw += kn.w[j];
        } else if feas_eq(kn.val[j], 0.0) {
            zeros.push(j);
        } else {
            items.push(j);
            itemsw += kn.w[j];
        }
    }
    if items.is_empty() {
        return Some(Cover {
            cover: vec![],
            noncover: vec![],
            weight: 0,
            fractional: false,
            ntightened: nt,
        });
    }
    let tcap = onesw + itemsw - kn.cap - 1;
    if tcap < 0 {
        return None;
    }
    // Dantzig: sort by profit/weight descending, take the prefix that fits (z=1 -> non-cover)
    let profit = |j: usize| {
        let p = 1.0 - kn.val[j];
        if modtrans { p * kn.w[j] as f64 } else { p }
    };
    let mut ord = items.clone();
    ord.sort_by(|&p, &q| {
        let rp = profit(p) / kn.w[p] as f64;
        let rq = profit(q) / kn.w[q] as f64;
        rq.partial_cmp(&rp).unwrap().then(p.cmp(&q))
    });
    let mut sw = 0i64;
    let mut split = 0usize;
    while split < ord.len() && sw + kn.w[ord[split]] <= tcap {
        sw += kn.w[ord[split]];
        split += 1;
    }
    let mut noncover: Vec<usize> = ord[..split].to_vec();
    let mut cover: Vec<usize> = ord[split..].to_vec();
    let mut weight: i64 = cover.iter().map(|&j| kn.w[j]).sum();
    for &j in &ones {
        cover.push(j);
        weight += kn.w[j];
    }
    noncover.extend(zeros.iter().copied());
    if weight <= kn.cap {
        return None;
    }
    Some(Cover {
        cover,
        noncover,
        weight,
        fractional: true,
        ntightened: nt,
    })
}

/// cons_knapsack.c `makeCoverMinimal` (modified transformed variant).
fn make_cover_minimal(kn: &Knap, cv: &mut Cover) {
    // ascending (x, a); ties stable
    cv.cover.sort_by(|&p, &q| {
        kn.val[p]
            .partial_cmp(&kn.val[q])
            .unwrap()
            .then((kn.w[p] as f64).partial_cmp(&(kn.w[q] as f64)).unwrap())
    });
    let minw = |cs: &Vec<usize>| cs.iter().map(|&j| kn.w[j]).min().unwrap_or(0);
    let mut j = 0;
    let mut mw = minw(&cv.cover);
    while j < cv.cover.len() && cv.weight - mw > kn.cap {
        let cj = cv.cover[j];
        if cv.weight - kn.w[cj] <= kn.cap {
            j += 1;
            continue;
        }
        cv.noncover.push(cj);
        cv.weight -= kn.w[cj];
        cv.cover.remove(j);
        if cv.cover.is_empty() {
            break;
        }
        mw = minw(&cv.cover);
    }
}

fn ge1(v: f64) -> bool {
    v >= 1.0 - 1e-9
}

/// Partition of the non-cover set into F (x>0) and R (x=0), `getPartitionNoncovervars`.
fn partition_non(kn: &Knap, non: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let (mut f, mut r) = (Vec::new(), Vec::new());
    for &j in non {
        if feas_eq(kn.val[j], 0.0) {
            r.push(j);
        } else {
            f.push(j);
        }
    }
    (f, r)
}

/// `getLiftingSequence`: F by (x, a) descending, C2 and R by a descending.
fn lifting_sequence(kn: &Knap, f: &mut [usize], c2: &mut [usize], r: &mut [usize]) {
    f.sort_by(|&p, &q| {
        kn.val[q]
            .partial_cmp(&kn.val[p])
            .unwrap()
            .then(kn.w[q].cmp(&kn.w[p]))
    });
    c2.sort_by(|&p, &q| kn.w[q].cmp(&kn.w[p]));
    r.sort_by(|&p, &q| kn.w[q].cmp(&kn.w[p]));
}

/// `separateSequLiftedMinimalCoverInequality` (no GUBs).
fn sequ_lifted_mincover(kn: &Knap, cover: &[usize], non: &[usize]) -> Option<KCut> {
    let mut c1 = Vec::new();
    let mut c2 = Vec::new();
    for &j in cover {
        if ge1(kn.val[j]) {
            c2.push(j);
        } else {
            c1.push(j);
        }
    }
    if c1.len() < 2 && !c2.is_empty() {
        c2.sort_by(|&p, &q| kn.w[q].cmp(&kn.w[p]));
        while c1.len() < 2 && !c2.is_empty() {
            c1.push(c2.pop().unwrap());
        }
    }
    if c1.is_empty() {
        return None;
    }
    let (mut f, mut r) = partition_non(kn, non);
    lifting_sequence(kn, &mut f, &mut c2, &mut r);
    let alpha0 = c1.len() as i64 - 1;
    lift_and_emit(kn, c1, c2, f, r, alpha0)
}

/// `separateSequLiftedExtendedWeightInequality`.
fn sequ_lifted_ewi(kn: &Knap, feas: &[usize], non: &[usize]) -> Option<KCut> {
    let mut t1 = Vec::new();
    let mut t2 = Vec::new();
    for &j in feas {
        if ge1(kn.val[j]) {
            t2.push(j);
        } else {
            t1.push(j);
        }
    }
    if t1.is_empty() && !t2.is_empty() {
        t2.sort_by(|&p, &q| kn.w[q].cmp(&kn.w[p]));
        t1.push(t2.pop().unwrap());
    }
    let (mut f, mut r) = partition_non(kn, non);
    lifting_sequence(kn, &mut f, &mut t2, &mut r);
    let alpha0 = t1.len() as i64;
    lift_and_emit(kn, t1, t2, f, r, alpha0)
}

fn lift_and_emit(
    kn: &Knap,
    m1: Vec<usize>,
    m2: Vec<usize>,
    f: Vec<usize>,
    r: Vec<usize>,
    alpha0: i64,
) -> Option<KCut> {
    let (coefs, cutact, liftrhs) = seq_up_down_lifting(kn, &m1, &m2, &f, &r, alpha0);
    let eff = (cutact - liftrhs as f64) / (liftrhs.max(1) as f64).sqrt();
    if eff <= 1e-6 {
        return None;
    }
    let mut terms: Vec<(usize, i64)> = Vec::new();
    for &j in &m1 {
        terms.push((j, 1));
    }
    for &j in m2.iter().chain(f.iter()).chain(r.iter()) {
        if coefs[j] > 0 {
            terms.push((j, coefs[j]));
        }
    }
    Some((terms, liftrhs))
}

fn enlarge(mw: &mut Vec<i64>, newlen: usize) {
    while mw.len() < newlen {
        mw.push(i64::MAX);
    }
}

/// minweight table update after lifting with `(coef, weight)`, over `lo..=hi` descending.
fn upd_table(mw: &mut [i64], top: usize, coef: usize, weight: i64) {
    for w in (0..=top).rev() {
        let m = if w < coef {
            mw[w].min(weight)
        } else {
            mw[w].min(mw[w - coef].saturating_add(weight))
        };
        mw[w] = m;
    }
}

/// Largest `z` in `0..right` with `mw[z] <= bound` (binary search as in SCIP).
fn bsearch(mw: &[i64], right0: usize, bound: i64) -> usize {
    let (mut left, mut right) = (0usize, right0);
    while left + 1 < right {
        let mid = (left + right) / 2;
        if mw[mid] <= bound {
            left = mid;
        } else {
            right = mid;
        }
    }
    left
}

/// cons_knapsack.c `sequentialUpAndDownLifting`; returns (coefficients, activity, rhs).
fn seq_up_down_lifting(
    kn: &Knap,
    m1: &[usize],
    m2: &[usize],
    f: &[usize],
    r: &[usize],
    alpha0: i64,
) -> (Vec<i64>, f64, i64) {
    let cap = kn.cap;
    let mut coefs = vec![0i64; kn.w.len()];
    let mut cutact = 0.0;
    let mut m1s: Vec<usize> = m1.to_vec();
    for &j in &m1s {
        coefs[j] = 1;
        cutact += kn.val[j];
    }
    m1s.sort_by(|&p, &q| kn.w[p].cmp(&kn.w[q]));
    let mut mw: Vec<i64> = vec![0; m1s.len() + 1];
    for w in 1..=m1s.len() {
        mw[w] = mw[w - 1] + kn.w[m1s[w - 1]];
    }
    let mut fixedones: i64 = m2.iter().map(|&j| kn.w[j]).sum();
    let mut liftrhs = alpha0;
    // up-lift F
    for &lv in f {
        let wt = kn.w[lv];
        let z = if cap - fixedones - wt < 0 {
            0
        } else if mw[liftrhs as usize] <= cap - fixedones - wt {
            liftrhs
        } else {
            bsearch(&mw, liftrhs as usize + 1, cap - fixedones - wt) as i64
        };
        let lc = liftrhs - z;
        coefs[lv] = lc;
        if lc == 0 {
            continue;
        }
        cutact += lc as f64 * kn.val[lv];
        let newlen = mw.len() + lc as usize;
        enlarge(&mut mw, newlen);
        let top = mw.len() - 1;
        upd_table(&mut mw, top, lc as usize, wt);
    }
    // down-lift M2
    for &lv in m2 {
        let wt = kn.w[lv];
        let z = bsearch(&mw, mw.len(), cap - fixedones + wt) as i64;
        let lc = (z - liftrhs).max(0);
        coefs[lv] = lc;
        fixedones -= wt;
        liftrhs += lc;
        if lc == 0 {
            continue;
        }
        cutact += lc as f64 * kn.val[lv];
        let newlen = mw.len() + lc as usize;
        enlarge(&mut mw, newlen);
        let top = mw.len() - 1;
        upd_table(&mut mw, top, lc as usize, wt);
    }
    // up-lift R
    for &lv in r {
        let wt = kn.w[lv];
        if cap - wt < 0 {
            continue;
        }
        let lr = liftrhs as usize;
        let z = if lr < mw.len() && mw[lr] <= cap - wt {
            liftrhs
        } else {
            bsearch(&mw, (lr + 1).min(mw.len()), cap - wt) as i64
        };
        let lc = liftrhs - z;
        coefs[lv] = lc;
        if lc == 0 {
            continue;
        }
        cutact += lc as f64 * kn.val[lv];
        let top = lr.min(mw.len() - 1);
        upd_table(&mut mw, top, lc as usize, wt);
    }
    (coefs, cutact, liftrhs)
}
