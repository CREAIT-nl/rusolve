//! rusolve `mip_dff_separator`: dual-feasible-function cuts for
//! activation knapsacks, separated at the LP point.
//!
//! Not in HiGHS. A row `sum_i a_i x_i - b y <= 0` (x, y binary, a_i > 0, b > 0) implies
//! `sum_i f(a_i / b) x_i <= y` for every dual feasible function `f` of bin packing (Fekete and
//! Schepers, "New classes of fast lower bounds for bin packing problems", Math. Programming 91,
//! 2001; Clautiaux, Alves and Valerio de Carvalho, "A survey of dual-feasible and superadditive
//! functions", Annals of Operations Research 179, 2010). The family is the one of
//! `activation_dff.rs` (`mip_activation_dff` adds one row per capacity row up front); here each
//! round evaluates every member at the current LP point and cuts the row of the function with the
//! largest violation, if it exceeds the efficacy threshold. Candidate rows are detected once on
//! the global model rows.

use crate::mip::activation_dff::{Dff, dff_candidates};
use crate::prelude::*;

/// Same efficacy threshold as the knapsack separator.
const MIN_EFFICACY: f64 = 1e-4;
/// Cuts handed to the pool per round.
const MAX_CUTS: usize = 2000;

struct Cand {
    /// None: the activation variable is fixed at 1 in the global domain (plain knapsack, rhs 1)
    y: Option<usize>,
    items: Vec<(usize, f64)>,
    funcs: Vec<Dff>,
}

#[derive(Clone, Default)]
pub struct HighsDffSeparator {
    pub base: HighsSeparator,
    cands: Option<std::sync::Arc<Vec<Cand>>>,
}

impl std::fmt::Debug for HighsDffSeparator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HighsDffSeparator")
    }
}

impl HighsDffSeparator {
    pub fn new(ms: &HighsMipSolver) -> Self {
        HighsDffSeparator {
            base: HighsSeparator::new(ms, "Separation: DFF"),
            cands: None,
        }
    }

    fn detect(ms: &HighsMipSolver, lp_relaxation: &mut HighsLpRelaxation, g: usize) -> Vec<Cand> {
        let nrow = ms.num_row();
        let lower = &ms.mipdata_.domains[g].col_lower_;
        let upper = &ms.mipdata_.domains[g].col_upper_;
        let is_bin = |cc: usize| {
            !matches!(
                ms.variable_type_2(cc as HighsInt),
                HighsVarType::kContinuous
            ) && lower[cc] == 0.0
                && upper[cc] == 1.0
        };
        let mut out = Vec::new();
        for r in 0..nrow {
            let (lo, up) = {
                let lp = lp_relaxation.get_lp();
                (lp.row_lower_[r as usize], lp.row_upper_[r as usize])
            };
            let sign = if lo <= -kHighsInf && up == 0.0 {
                1.0
            } else if up >= kHighsInf && lo == 0.0 {
                -1.0
            } else {
                continue;
            };
            let mut len: HighsInt = 0;
            let mut inds: &[HighsInt] = &[];
            let mut vals: &[f64] = &[];
            lp_relaxation.get_row(ms, r, &mut len, &mut inds, &mut vals);
            if len < 3 {
                continue;
            }
            let mut neg: Vec<(usize, f64)> = Vec::new();
            let mut items: Vec<(usize, f64)> = Vec::new();
            let mut ok = true;
            let mut fixed_y = false;
            for k in 0..len as usize {
                let cc = inds[k] as usize;
                // a column fixed at 0 in the global domain drops out of the row
                if upper[cc] == 0.0 && lower[cc] == 0.0 {
                    continue;
                }
                if vals[k] * sign < 0.0 && lower[cc] == 1.0 && upper[cc] == 1.0 {
                    fixed_y = true;
                    neg.push((cc, -sign * vals[k]));
                    continue;
                }
                if !is_bin(cc) {
                    ok = false;
                    break;
                }
                let a = sign * vals[k];
                if a < 0.0 {
                    neg.push((cc, -a));
                } else {
                    items.push((cc, a));
                }
            }
            if !ok || neg.len() != 1 {
                continue;
            }
            let (yc, b) = neg[0];
            let y = if fixed_y { None } else { Some(yc) };
            let items: Vec<(usize, f64)> = items.into_iter().map(|(cc, a)| (cc, a / b)).collect();
            if items.len() < 2 {
                continue;
            }
            let funcs = dff_candidates(&items, 1);
            out.push(Cand { y, items, funcs });
        }
        out
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
        if self.cands.is_none() {
            let det = Self::detect(ms, lp_relaxation, g);
            self.cands = Some(std::sync::Arc::new(det));
        }
        let cands = self.cands.clone().unwrap();
        let sol: Vec<f64> = lp_relaxation.get_solution().col_value.clone();
        let mut found: Vec<(f64, Vec<HighsInt>, Vec<f64>, bool)> = Vec::new();
        for cd in cands.iter() {
            let ysol = cd.y.map_or(1.0, |y| sol[y]);
            let mut best: Option<(f64, usize)> = None;
            for (fi, f) in cd.funcs.iter().enumerate() {
                let (mut act, mut nrm) = (0.0f64, if cd.y.is_some() { 1.0f64 } else { 0.0 });
                for &(j, x) in &cd.items {
                    let w = f.eval(x);
                    if w > 1e-12 {
                        act += w * sol[j];
                        nrm += w * w;
                    }
                }
                let eff = (act - ysol) / nrm.sqrt();
                if eff > MIN_EFFICACY && best.is_none_or(|b| eff > b.0) {
                    best = Some((eff, fi));
                }
            }
            if let Some((eff, fi)) = best {
                let f = cd.funcs[fi];
                let mut ind: Vec<HighsInt> = Vec::new();
                let mut val: Vec<f64> = Vec::new();
                for &(j, x) in &cd.items {
                    let w = f.eval(x);
                    if w > 1e-12 {
                        ind.push(j as HighsInt);
                        val.push(w);
                    }
                }
                if let Some(y) = cd.y {
                    ind.push(y as HighsInt);
                    val.push(-1.0);
                }
                found.push((eff, ind, val, cd.y.is_none()));
            }
        }
        found.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        found.truncate(MAX_CUTS);
        for (_, mut inds, mut vals, rhs1) in found {
            let len = inds.len() as HighsInt;
            detach!(ms.mipdata_.cutpools[cutpool], |p| p.add_cut(
                ms,
                &mut inds[..],
                &mut vals[..],
                len,
                if rhs1 { 1.0 } else { 0.0 },
                false,
                true,
                false,
                false,
            ));
        }
    }
}
