//! rusolve ls-gub-swap: GUB swap and ejection-chain compound moves for LocalMIP.
//!
//! Not a translation: original rusolve code, behind
//! `mip_ls_gub_swap` (default off). It plugs into the mip-heuristics LocalMIP port
//! as phase 3b of `lm_search.rs::infeasible_step` and scores with the same rules as
//! `compute_candidate_scores`.
//!
//! Provenance:
//! * Compound moves that keep a "constrained sum" satisfied, k-paths and ejection
//!   chains: Benoist, Estellon, Gardi, Megel, Nouioua, "LocalSolver 1.x: a black-box
//!   local-search solver for 0-1 programming", 4OR 9(3):299-316, 2011.
//! * Ejection chains: Glover, Discrete Appl. Math. 65:223-253, 1996; shift/swap chains
//!   over GAP's assignment rows: Yagiura, Ibaraki, Glover, INFORMS J. Computing
//!   16(2):133-151, 2004.
//! * Compound multi-variable moves in an FJ-style LS: Davies, Didier, Perron,
//!   "ViolationLS", CPAIOR 2024.
//! * The score (Defs. 5-9): Lin, Zou, Cai, "An Efficient Local Search Solver for Mixed
//!   Integer Programming", CP 2024, as ported in `lm_search.rs`.

use super::rng::S5Rng;
use super::*;
use crate::prelude::*;

/// Most flips in one compound move (a depth-bounded ejection chain).
pub(crate) const GUB_MAX_FLIPS: usize = 6;
/// Most extension steps of a chain after its seed flip.
pub(crate) const GUB_MAX_DEPTH: usize = 3;
/// Refill candidates tried for an emptied partition row.
pub(crate) const GUB_BRANCH: usize = 8;
/// Cells (row and column entries) one call may scan.
// TODO(perf): unmeasured; about 2x phase 1's kBmsBudget (2250 row entries).
pub(crate) const GUB_EVAL_BUDGET: usize = 4096;
/// Violated (infeasible point) or GUB (feasible point) rows sampled per call;
/// the phase 1 count of `lm_search.rs`.
pub(crate) const GUB_SEED_ROWS: HighsInt = kBmsConstraints;

const GUB_COEFF_RTOL: f64 = 1e-12;
const GUB_RHS_TOL: f64 = 1e-9;

const KIND_NONE: u8 = 0;
const KIND_PACKING: u8 = 1;
const KIND_PARTITION: u8 = 2;

/// Immutable GUB structure of one problem.
#[derive(Clone, Debug, Default)]
pub(crate) struct GubIndex {
    /// Per row: KIND_NONE, KIND_PACKING (`Σ x ≤ 1`) or KIND_PARTITION (`Σ x = 1`).
    pub(crate) kind: Vec<u8>,
    /// The GUB rows, ascending.
    pub(crate) rows: Vec<HighsInt>,
    /// Column -> its GUB rows (CSR).
    pub(crate) col_start: Vec<HighsInt>,
    pub(crate) col_row: Vec<HighsInt>,
}

/// Counters, for tests and the eval note.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct GubStats {
    pub(crate) calls: u64,
    pub(crate) evaluated: u64,
    pub(crate) applied: u64,
    pub(crate) applied_flips: u64,
}

/// The operator state a worker owns (`WorkerCtx::gub`).
#[derive(Clone, Debug)]
pub(crate) struct GubSwap {
    pub(crate) index: GubIndex,
    row_delta: Vec<f64>,
    row_mark: Vec<bool>,
    touched: Vec<HighsInt>,
    chain: Vec<(HighsInt, f64)>,
    best: Vec<(HighsInt, f64)>,
    best_score: f64,
    best_bonus: f64,
    spent: usize,
    pub(crate) stats: GubStats,
}

/// Row kind under the GUB detection rule: at least 2 binary entries, one common coefficient
/// `v`, scaled upper bound 1; a partition row when the scaled lower bound is positive, else
/// a packing row.
fn classify_row(ctx: &WorkerCtx, i: HighsInt) -> u8 {
    let start = ctx.ar_start[i as usize] as usize;
    let end = ctx.ar_start[(i + 1) as usize] as usize;
    if end - start < 2 {
        return KIND_NONE;
    }
    let v = ctx.ar_value[start];
    if v == 0.0 || !v.is_finite() {
        return KIND_NONE;
    }
    for k in start..end {
        if !ctx.is_binary(ctx.ar_index[k]) {
            return KIND_NONE;
        }
        if (ctx.ar_value[k] - v).abs() > GUB_COEFF_RTOL * v.abs() {
            return KIND_NONE;
        }
    }
    let (mut lo, mut hi) = (ctx.row_lo[i as usize] / v, ctx.row_hi[i as usize] / v);
    if v < 0.0 {
        std::mem::swap(&mut lo, &mut hi);
    }
    if (hi - 1.0).abs() > GUB_RHS_TOL {
        return KIND_NONE;
    }
    if lo > GUB_RHS_TOL {
        // Σ x ≥ lo > 0 over binaries is Σ x ≥ 1.
        if lo > 1.0 + GUB_RHS_TOL {
            return KIND_NONE;
        }
        KIND_PARTITION
    } else {
        KIND_PACKING
    }
}

impl GubIndex {
    pub(crate) fn build(ctx: &WorkerCtx) -> GubIndex {
        let nrow = ctx.nrow;
        let ncol = ctx.ncol as usize;
        let mut kind = vec![KIND_NONE; nrow as usize];
        let mut rows = Vec::new();
        let mut count = vec![0 as HighsInt; ncol + 1];
        for i in 0..nrow {
            let kd = classify_row(ctx, i);
            if kd == KIND_NONE {
                continue;
            }
            kind[i as usize] = kd;
            rows.push(i);
            for k in ctx.ar_start[i as usize]..ctx.ar_start[(i + 1) as usize] {
                count[ctx.ar_index[k as usize] as usize + 1] += 1;
            }
        }
        for j in 0..ncol {
            count[j + 1] += count[j];
        }
        let col_start = count.clone();
        let mut fill = count;
        let mut col_row = vec![0 as HighsInt; col_start[ncol] as usize];
        for &i in &rows {
            for k in ctx.ar_start[i as usize]..ctx.ar_start[(i + 1) as usize] {
                let j = ctx.ar_index[k as usize] as usize;
                col_row[fill[j] as usize] = i;
                fill[j] += 1;
            }
        }
        GubIndex {
            kind,
            rows,
            col_start,
            col_row,
        }
    }

    fn col_rows(&self, j: HighsInt) -> std::ops::Range<usize> {
        self.col_start[j as usize] as usize..self.col_start[(j + 1) as usize] as usize
    }
}

impl GubSwap {
    /// The operator for `ctx`'s problem, or None when it has no GUB row.
    // TODO(perf): rebuilt with every worker (each LocalMIP restart); O(nnz), the
    // order of `WorkerCtx::new` + `rebuild_state`.
    pub(crate) fn build(ctx: &WorkerCtx) -> Option<Box<GubSwap>> {
        let index = GubIndex::build(ctx);
        if index.rows.is_empty() {
            return None;
        }
        Some(Box::new(GubSwap {
            index,
            row_delta: vec![0.0; ctx.nrow as usize],
            row_mark: vec![false; ctx.nrow as usize],
            touched: Vec::new(),
            chain: Vec::with_capacity(GUB_MAX_FLIPS),
            best: Vec::with_capacity(GUB_MAX_FLIPS),
            best_score: f64::NEG_INFINITY,
            best_bonus: 0.0,
            spent: 0,
            stats: GubStats::default(),
        }))
    }

    fn charge(&mut self, ctx: &mut WorkerCtx, cells: usize) {
        self.spent += cells;
        ctx.effort += cells;
    }

    fn exhausted(&self) -> bool {
        self.spent >= GUB_EVAL_BUDGET
    }

    /// Value of `j` in the current point with the chain applied.
    fn virtual_val(&self, ctx: &WorkerCtx, j: HighsInt) -> f64 {
        for &(v, val) in &self.chain {
            if v == j {
                return val;
            }
        }
        ctx.solution[j as usize]
    }

    fn in_chain(&self, j: HighsInt) -> bool {
        self.chain.iter().any(|&(v, _)| v == j)
    }

    /// Joint score of the chain: `compute_candidate_scores` summed over its flips,
    /// with each shared row scored once on its summed delta.
    fn evaluate(
        &mut self,
        ctx: &mut WorkerCtx,
        step: HighsInt,
        best_feasible: bool,
        best_obj: f64,
    ) {
        self.stats.evaluated += 1;
        let mut obj_delta = 0.0;
        let mut tabu = false;
        let mut cells = 0usize;
        for &(j, val) in &self.chain {
            let delta = val - ctx.solution[j as usize];
            obj_delta += ctx.col_cost[j as usize] * delta;
            if ctx.is_tabu(j, delta, step) {
                tabu = true;
            }
            let (c0, c1) = (
                ctx.csc.col_start[j as usize] as usize,
                ctx.csc.col_start[(j + 1) as usize] as usize,
            );
            cells += c1 - c0;
            for p in c0..c1 {
                let i = ctx.csc.col_row[p];
                if !self.row_mark[i as usize] {
                    self.row_mark[i as usize] = true;
                    self.touched.push(i);
                }
                self.row_delta[i as usize] += ctx.csc.col_val[p] * delta;
            }
        }
        self.charge(ctx, cells);

        let eps = ctx.epsilon;
        let new_obj = ctx.current_obj + obj_delta;
        let aspiration = best_feasible && new_obj < best_obj - eps;
        let skip = tabu && !aspiration;

        let mut progress = 0.0;
        let mut bonus = 0.0;
        if !skip {
            // Def 5 / Def 8 (lm_search.rs::compute_candidate_scores)
            if new_obj < ctx.current_obj - eps {
                progress += ctx.obj_weight as f64;
            } else if new_obj > ctx.current_obj + eps {
                progress -= ctx.obj_weight as f64;
            }
            if aspiration {
                bonus += ctx.obj_weight as f64;
            }
        }
        for t in 0..self.touched.len() {
            let i = self.touched[t];
            let d = self.row_delta[i as usize];
            self.row_delta[i as usize] = 0.0;
            self.row_mark[i as usize] = false;
            if skip {
                continue;
            }
            let (p, b) = score_row(ctx, i, d);
            progress += p;
            bonus += b;
        }
        self.touched.clear();
        if skip {
            return;
        }

        // select_best_from_batch's rule
        if progress > self.best_score + kScoreTol
            || (progress > self.best_score - kScoreTol && bonus > self.best_bonus)
        {
            self.best_score = progress;
            self.best_bonus = bonus;
            self.best.clear();
            self.best.extend_from_slice(&self.chain);
        }
    }

    /// Scores the chain (two or more flips) and extends it by one ejection step.
    fn extend(
        &mut self,
        ctx: &mut WorkerCtx,
        rng: &mut S5Rng,
        step: HighsInt,
        best_feasible: bool,
        best_obj: f64,
        depth: usize,
        front: usize,
    ) {
        if self.chain.len() >= 2 {
            self.evaluate(ctx, step, best_feasible, best_obj);
        }
        if depth >= GUB_MAX_DEPTH || self.chain.len() >= GUB_MAX_FLIPS || self.exhausted() {
            return;
        }
        let len = self.chain.len();
        let raising = self.chain[front].1 > 0.5;
        if raising {
            // Forced: eject every other 1 in the GUB rows of the raised variables.
            let mut cells = 0usize;
            for f in front..len {
                let v = self.chain[f].0;
                for r in self.index.col_rows(v) {
                    let g = self.index.col_row[r];
                    let (k0, k1) = (
                        ctx.ar_start[g as usize] as usize,
                        ctx.ar_start[(g + 1) as usize] as usize,
                    );
                    cells += k1 - k0;
                    for k in k0..k1 {
                        let m = ctx.ar_index[k];
                        if m != v && !self.in_chain(m) && self.virtual_val(ctx, m) > 0.5 {
                            if self.chain.len() >= GUB_MAX_FLIPS {
                                break;
                            }
                            self.chain.push((m, 0.0));
                        }
                    }
                }
            }
            self.charge(ctx, cells);
            if self.chain.len() > len {
                self.extend(ctx, rng, step, best_feasible, best_obj, depth + 1, len);
            }
            self.chain.truncate(len);
            return;
        }

        // Lowered variables: refill the first partition row left with no 1.
        let mut cells = 0usize;
        let mut refill: HighsInt = -1;
        'find: for f in front..len {
            let v = self.chain[f].0;
            for r in self.index.col_rows(v) {
                let g = self.index.col_row[r];
                if self.index.kind[g as usize] != KIND_PARTITION {
                    continue;
                }
                let (k0, k1) = (
                    ctx.ar_start[g as usize] as usize,
                    ctx.ar_start[(g + 1) as usize] as usize,
                );
                cells += k1 - k0;
                if (k0..k1).all(|k| self.virtual_val(ctx, ctx.ar_index[k]) < 0.5) {
                    refill = g;
                    break 'find;
                }
            }
        }
        self.charge(ctx, cells);
        if refill < 0 {
            return;
        }
        let (k0, k1) = (
            ctx.ar_start[refill as usize] as usize,
            ctx.ar_start[(refill + 1) as usize] as usize,
        );
        let n = k1 - k0;
        let offset = (rng.next_u64() % (n as u64)) as usize;
        let mut tried = 0usize;
        for s in 0..n {
            if tried >= GUB_BRANCH || self.exhausted() {
                break;
            }
            let m = ctx.ar_index[k0 + (offset + s) % n];
            if self.in_chain(m) {
                continue;
            }
            tried += 1;
            self.chain.push((m, 1.0));
            self.extend(ctx, rng, step, best_feasible, best_obj, depth + 1, len);
            self.chain.truncate(len);
        }
        self.charge(ctx, tried);
    }

    fn seed(
        &mut self,
        ctx: &mut WorkerCtx,
        rng: &mut S5Rng,
        step: HighsInt,
        best_feasible: bool,
        best_obj: f64,
        j: HighsInt,
        val: f64,
    ) {
        self.chain.clear();
        self.chain.push((j, val));
        self.extend(ctx, rng, step, best_feasible, best_obj, 0, 0);
        self.chain.clear();
    }

    fn run(
        &mut self,
        ctx: &mut WorkerCtx,
        rng: &mut S5Rng,
        step: HighsInt,
        best_feasible: bool,
        best_obj: f64,
    ) -> Option<Candidate> {
        self.stats.calls += 1;
        self.best.clear();
        self.best_score = f64::NEG_INFINITY;
        self.best_bonus = 0.0;
        self.spent = 0;

        if !ctx.violated.empty() {
            // Seeds: GUB variables whose flip reduces a sampled violated row.
            let nviol = ctx.violated.size();
            let nsample = GUB_SEED_ROWS.min(nviol);
            for s in 0..nsample {
                if self.exhausted() {
                    break;
                }
                let ci = if nsample == nviol {
                    ctx.violated.elements[s as usize]
                } else {
                    ctx.violated.elements[(rng.next_u64() % (nviol as u64)) as usize]
                };
                let (k0, k1) = (
                    ctx.ar_start[ci as usize] as usize,
                    ctx.ar_start[(ci + 1) as usize] as usize,
                );
                self.charge(ctx, k1 - k0);
                let old_lhs = ctx.lhs[ci as usize];
                let old_viol = ctx.compute_violation(ci, old_lhs);
                for k in k0..k1 {
                    if self.exhausted() {
                        break;
                    }
                    let j = ctx.ar_index[k];
                    if !ctx.is_binary(j) || self.index.col_rows(j).is_empty() {
                        continue;
                    }
                    let val = if ctx.solution[j as usize] < 0.5 {
                        1.0
                    } else {
                        0.0
                    };
                    let new_lhs = old_lhs + ctx.ar_value[k] * (val - ctx.solution[j as usize]);
                    if ctx.compute_violation(ci, new_lhs) < old_viol - kViolDeltaTol {
                        self.seed(ctx, rng, step, best_feasible, best_obj, j, val);
                    }
                }
            }
        } else {
            // Seeds: each 0 of a sampled GUB row that has a 1 (the swap).
            let ngub = self.index.rows.len() as HighsInt;
            let nsample = GUB_SEED_ROWS.min(ngub);
            for s in 0..nsample {
                if self.exhausted() {
                    break;
                }
                let g = if nsample == ngub {
                    self.index.rows[s as usize]
                } else {
                    self.index.rows[(rng.next_u64() % (ngub as u64)) as usize]
                };
                let (k0, k1) = (
                    ctx.ar_start[g as usize] as usize,
                    ctx.ar_start[(g + 1) as usize] as usize,
                );
                self.charge(ctx, k1 - k0);
                if !(k0..k1).any(|k| ctx.solution[ctx.ar_index[k] as usize] > 0.5) {
                    continue;
                }
                for k in k0..k1 {
                    if self.exhausted() {
                        break;
                    }
                    let j = ctx.ar_index[k];
                    if ctx.solution[j as usize] < 0.5 {
                        self.seed(ctx, rng, step, best_feasible, best_obj, j, 1.0);
                    }
                }
            }
        }

        if self.best.is_empty() || self.best_score <= kScoreTol {
            return None;
        }
        // The worker applies the head; the tail is applied here.
        for t in 1..self.best.len() {
            let (j, val) = self.best[t];
            ctx.apply_move_with_tabu(j, val, step, rng);
        }
        self.stats.applied += 1;
        self.stats.applied_flips += self.best.len() as u64;
        Some(Candidate {
            var_idx: self.best[0].0,
            new_val: self.best[0].1,
            score: self.best_score,
            bonus: self.best_bonus,
        })
    }
}

/// One row's share of the score when its activity moves by `d`: the row part of
/// `lm_search.rs::compute_candidate_scores` (Defs. 6 and 9, the authors' variants).
fn score_row(ctx: &WorkerCtx, i: HighsInt, d: f64) -> (f64, f64) {
    let lo = ctx.row_lo[i as usize];
    let hi = ctx.row_hi[i as usize];
    let old_lhs = ctx.lhs[i as usize];
    let new_lhs = old_lhs + d;
    let old_viol = ctx.compute_violation(i, old_lhs);
    let new_viol = ctx.compute_violation(i, new_lhs);
    let w = ctx.weight[i as usize] as f64;
    let was_viol = old_viol > ctx.feastol;
    let now_viol = new_viol > ctx.feastol;
    let mut progress = 0.0;
    if was_viol && !now_viol {
        progress += w;
    } else if !was_viol && now_viol {
        progress -= w;
    } else if was_viol && now_viol {
        if new_viol < old_viol - kViolDeltaTol {
            progress += w;
        } else if new_viol > old_viol + kViolDeltaTol {
            progress -= w;
        }
    }
    let mut bonus = 0.0;
    if !now_viol {
        let strict = |l: f64| {
            (hi >= kHighsInf || l < hi - ctx.feastol) && (lo <= -kHighsInf || l > lo + ctx.feastol)
        };
        let old_strict = !was_viol && strict(old_lhs);
        if !old_strict && strict(new_lhs) {
            bonus += w;
        }
    }
    (progress, bonus)
}

/// Phase 3b of `infeasible_step`: the best compound GUB move when it scores
/// above `kScoreTol`, with its tail flips already applied; the caller applies the
/// returned head flip. None (and no state change) otherwise.
pub(crate) fn gub_phase(
    ctx: &mut WorkerCtx,
    rng: &mut S5Rng,
    step: HighsInt,
    best_feasible: bool,
    best_obj: f64,
) -> Option<Candidate> {
    let mut gub = ctx.gub.take()?;
    let out = gub.run(ctx, rng, step, best_feasible, best_obj);
    ctx.gub = Some(gub);
    out
}

#[cfg(test)]
mod tests {
    use super::super::lm_search::infeasible_step;
    use super::super::local_mip;
    use super::super::rng::{heuristic_effort_budget, make_budget, patience_threshold};
    use super::super::sink::S5Sink;
    use super::*;

    type Row = (Vec<(HighsInt, f64)>, f64, f64);

    /// Binary problem from sparse rows (all columns binary in [0, 1]).
    fn problem(ncol: HighsInt, rows: &[Row], cost: &[f64]) -> S5Problem {
        let nrow = rows.len() as HighsInt;
        let mut p = S5Problem {
            ncol,
            nrow,
            col_lower: vec![0.0; ncol as usize],
            col_upper: vec![1.0; ncol as usize],
            col_cost: cost.to_vec(),
            integrality: vec![HighsVarType::kInteger; ncol as usize],
            binary: vec![1; ncol as usize],
            ar_start: vec![0],
            uplocks: vec![0; ncol as usize],
            downlocks: vec![0; ncol as usize],
            feastol: 1e-6,
            epsilon: 1e-9,
            ..S5Problem::default()
        };
        for (ents, lo, hi) in rows {
            for &(j, a) in ents {
                p.ar_index.push(j);
                p.ar_value.push(a);
            }
            p.ar_start.push(p.ar_index.len() as HighsInt);
            p.row_lower.push(*lo);
            p.row_upper.push(*hi);
        }
        p.nnz = p.ar_index.len();
        let mut col_start = vec![0 as HighsInt; ncol as usize + 1];
        for &j in &p.ar_index {
            col_start[j as usize + 1] += 1;
        }
        for j in 0..ncol as usize {
            col_start[j + 1] += col_start[j];
        }
        let mut fill = col_start.clone();
        let mut col_row = vec![0; p.nnz];
        let mut col_val = vec![0.0; p.nnz];
        for i in 0..nrow as usize {
            for k in p.ar_start[i] as usize..p.ar_start[i + 1] as usize {
                let j = p.ar_index[k] as usize;
                col_row[fill[j] as usize] = i as HighsInt;
                col_val[fill[j] as usize] = p.ar_value[k];
                fill[j] += 1;
            }
        }
        p.csc = CscMatrix {
            col_start,
            col_row,
            col_val,
        };
        p
    }

    fn unit(cols: &[HighsInt]) -> Vec<(HighsInt, f64)> {
        cols.iter().map(|&j| (j, 1.0)).collect()
    }

    fn ctx_at<'a>(p: &'a S5Problem, x: &[f64], gub: bool) -> WorkerCtx<'a> {
        let mut ctx = WorkerCtx::new(p);
        ctx.solution = x.to_vec();
        if gub {
            ctx.gub = GubSwap::build(&ctx);
        }
        ctx.rebuild_state();
        ctx
    }

    /// One host step at `x`: (feasible afterwards, objective afterwards, stats).
    fn host_step(p: &S5Problem, x: &[f64], gub: bool, best: Option<f64>) -> (bool, f64, GubStats) {
        let mut ctx = ctx_at(p, x, gub);
        let costed: Vec<HighsInt> = (0..p.ncol)
            .filter(|&j| p.col_cost[j as usize] != 0.0)
            .collect();
        let binary: Vec<HighsInt> = (0..p.ncol).collect();
        let mut rng = S5Rng::new(7);
        let cand = infeasible_step(
            &mut ctx,
            &mut rng,
            1,
            best.is_some(),
            best.unwrap_or(f64::INFINITY),
            &costed,
            &binary,
        );
        if cand.var_idx != -1 {
            ctx.apply_move_with_tabu(cand.var_idx, cand.new_val, 1, &mut rng);
        }
        let stats = ctx.gub.as_ref().map(|g| g.stats).unwrap_or_default();
        (ctx.full_recheck(false, false), ctx.current_obj, stats)
    }

    #[test]
    fn ls_gub_detects_partition_packing_and_rejects_others() {
        let inf = kHighsInf;
        let rows: Vec<Row> = vec![
            (unit(&[0, 1, 2]), 1.0, 1.0),            // partition
            (unit(&[0, 3]), -inf, 1.0),              // packing
            (vec![(1, 2.0), (3, 2.0)], -inf, 2.0),   // scaled packing
            (vec![(2, -1.0), (3, -1.0)], -1.0, inf), // negated packing
            (vec![(0, 1.0), (1, 2.0)], -inf, 1.0),   // not unit
            (unit(&[0, 1]), 1.0, inf),               // covering
            (unit(&[0, 1, 2, 3]), -inf, 2.0),        // cardinality 2
            (unit(&[2]), -inf, 1.0),                 // singleton
            (vec![(1, 3.0), (2, 3.0)], 1.5, 3.0),    // Σ x ≥ 0.5: partition
        ];
        let p = problem(4, &rows, &[0.0; 4]);
        let ctx = WorkerCtx::new(&p);
        let idx = GubIndex::build(&ctx);
        assert_eq!(idx.kind, vec![2, 1, 1, 1, 0, 0, 0, 0, 2]);
        assert_eq!(idx.rows, vec![0, 1, 2, 3, 8]);
        let r: Vec<HighsInt> = idx.col_rows(3).map(|k| idx.col_row[k]).collect();
        assert_eq!(r, vec![1, 2, 3]);

        // A continuous column disqualifies its row.
        let mut q = problem(2, &[(unit(&[0, 1]), 1.0, 1.0)], &[0.0; 2]);
        q.binary[1] = 0;
        q.integrality[1] = HighsVarType::kContinuous;
        assert!(GubSwap::build(&WorkerCtx::new(&q)).is_none());
        // So does a problem with no GUB row at all.
        let q = problem(2, &[(vec![(0, 1.0), (1, 2.0)], -inf, 2.0)], &[1.0, 1.0]);
        assert!(GubSwap::build(&WorkerCtx::new(&q)).is_none());
    }

    #[test]
    fn ls_gub_feasible_swap_where_single_flips_stall() {
        // min 5 x0 + 3 x1 + x2, x0 + x1 + x2 = 1, at x = (1, 0, 0). Every single
        // flip breaks the row, so no single move scores > 0; the swap does.
        let p = problem(3, &[(unit(&[0, 1, 2]), 1.0, 1.0)], &[5.0, 3.0, 1.0]);
        let x = [1.0, 0.0, 0.0];
        let (feas_off, _, _) = host_step(&p, &x, false, Some(5.0));
        assert!(!feas_off, "a single flip cannot stay feasible");
        let (feas_on, obj_on, st) = host_step(&p, &x, true, Some(5.0));
        assert!(feas_on);
        assert!(obj_on < 5.0, "obj {obj_on}");
        assert_eq!((st.calls, st.applied, st.applied_flips), (1, 1, 2));
    }

    #[test]
    fn ls_gub_assignment_repair_where_single_flips_stall() {
        // 2 jobs x 2 agents, cols (0: job0/agent0, 1: job0/agent1, 2: job1/agent0,
        // 3: job1/agent1). Both jobs on agent 0 overload it (6 > 4); moving one
        // job is a swap.
        let inf = kHighsInf;
        let rows: Vec<Row> = vec![
            (unit(&[0, 1]), 1.0, 1.0),
            (unit(&[2, 3]), 1.0, 1.0),
            (vec![(0, 3.0), (2, 3.0)], -inf, 4.0),
            (vec![(1, 3.0), (3, 3.0)], -inf, 4.0),
        ];
        let p = problem(4, &rows, &[0.0; 4]);
        let x = [1.0, 0.0, 1.0, 0.0];
        let (feas_off, _, _) = host_step(&p, &x, false, None);
        assert!(!feas_off);
        let (feas_on, _, st) = host_step(&p, &x, true, None);
        assert!(feas_on);
        assert_eq!((st.applied, st.applied_flips), (1, 2));
    }

    #[test]
    fn ls_gub_ejection_chain_refills_a_second_row() {
        // x0 + x1 = 1, x0 + x3 = 1, x2 + x3 <= 1. From x0 = 1 the swap x1 <- x0
        // empties row 1, so the chain refills it with x3: three flips.
        let rows: Vec<Row> = vec![
            (unit(&[0, 1]), 1.0, 1.0),
            (unit(&[0, 3]), 1.0, 1.0),
            (unit(&[2, 3]), -kHighsInf, 1.0),
        ];
        let p = problem(4, &rows, &[10.0, 1.0, 0.0, 1.0]);
        let x = [1.0, 0.0, 0.0, 0.0];
        let (feas_on, obj_on, st) = host_step(&p, &x, true, Some(10.0));
        assert!(feas_on);
        assert_eq!(obj_on, 2.0);
        assert_eq!(st.applied_flips, 3);
    }

    #[test]
    fn ls_gub_tabu_blocks_a_compound_without_aspiration() {
        let p = problem(3, &[(unit(&[0, 1, 2]), 1.0, 1.0)], &[5.0, 3.0, 1.0]);
        let mut ctx = ctx_at(&p, &[1.0, 0.0, 0.0], true);
        // x0 may not decrease until step 100, and every swap lowers x0.
        ctx.tabu_dec_until[0] = 100;
        let mut rng = S5Rng::new(1);
        // No feasible best known: no aspiration, so nothing is admissible.
        assert!(gub_phase(&mut ctx, &mut rng, 1, false, f64::INFINITY).is_none());
        assert_eq!(ctx.solution, vec![1.0, 0.0, 0.0]);
        // With a best of 5 the swaps beat it, so aspiration admits them.
        assert!(gub_phase(&mut ctx, &mut rng, 1, true, 5.0).is_some());
    }

    /// Generalized assignment in the shape of Yagiura, Ibaraki, Glover (2004):
    /// one partition row per job, one capacity row per agent.
    fn gap_problem(agents: HighsInt, jobs: HighsInt, seed: u64, slack: f64) -> S5Problem {
        let mut rng = S5Rng::new(seed);
        let col = |a: HighsInt, j: HighsInt| a * jobs + j;
        let ncol = agents * jobs;
        let mut cost = vec![0.0; ncol as usize];
        let mut w = vec![0.0; ncol as usize];
        for jc in 0..ncol as usize {
            cost[jc] = (10 + rng.next_u64() % 40) as f64;
            w[jc] = (5 + rng.next_u64() % 20) as f64;
        }
        let mut rows: Vec<Row> = Vec::new();
        for j in 0..jobs {
            let cols: Vec<HighsInt> = (0..agents).map(|a| col(a, j)).collect();
            rows.push((unit(&cols), 1.0, 1.0));
        }
        for a in 0..agents {
            let ents: Vec<(HighsInt, f64)> = (0..jobs)
                .map(|j| (col(a, j), w[col(a, j) as usize]))
                .collect();
            let total: f64 = ents.iter().map(|e| e.1).sum();
            rows.push((ents, -kHighsInf, (slack * total / agents as f64).floor()));
        }
        problem(ncol, &rows, &cost)
    }

    /// `local_mip::run` under the budget: (accepted offers, best objective, effort).
    fn run_lm(p: &S5Problem, gub: bool, seed: u32) -> (usize, f64, usize) {
        let total = heuristic_effort_budget(p.nnz, S5_LOCAL_MIP_EFFORT);
        let budget = make_budget(
            total,
            1,
            patience_threshold(p.nnz, S5_LOCAL_MIP_PATIENCE, total),
        );
        let exec = ExecutionContext {
            base_seed: seed,
            gub_swap: gub,
            ..ExecutionContext::default()
        };
        let mut sink = S5Sink::new(&p.integrality, 0);
        let out = local_mip::run(p, &budget, &exec, &mut sink);
        let mut best = Vec::new();
        let obj = if sink.copy_best(&mut best) {
            assert!(local_mip::is_solution_feasible(p, &best));
            best.iter().zip(&p.col_cost).map(|(x, cj)| x * cj).sum()
        } else {
            f64::INFINITY
        };
        (sink.accepted(), obj, out.effort)
    }

    #[test]
    fn ls_gub_run_improves_gap_primal() {
        // 5 agents x 30 jobs. Loose capacity (slack 1.0): both arms find points,
        // the swap finds better ones. Tight (0.8): the swap finds points more often.
        for seed in 0..6u64 {
            let p = gap_problem(5, 30, seed, 1.0);
            let off = run_lm(&p, false, seed as u32);
            let on = run_lm(&p, true, seed as u32);
            assert!(on.1 < off.1, "seed {seed}: off {off:?} on {on:?}");
        }
        let found = |gub: bool| {
            (0..6u64)
                .filter(|&s| run_lm(&gap_problem(5, 30, s, 0.8), gub, s as u32).0 > 0)
                .count()
        };
        assert!(found(true) > found(false));
    }

    #[test]
    fn ls_gub_run_without_gub_rows_is_identical() {
        // Multi-knapsack, no unit-coefficient rows: the flag builds no operator.
        let mut rng = S5Rng::new(3);
        let n = 40;
        let cost: Vec<f64> = (0..n)
            .map(|_| -((1 + rng.next_u64() % 30) as f64))
            .collect();
        let rows: Vec<Row> = (0..4)
            .map(|_| {
                let ents: Vec<(HighsInt, f64)> = (0..n)
                    .map(|j| (j, (2 + rng.next_u64() % 20) as f64))
                    .collect();
                let cap = (0.3 * ents.iter().map(|e| e.1).sum::<f64>()).floor();
                (ents, -kHighsInf, cap)
            })
            .collect();
        let p = problem(n, &rows, &cost);
        for seed in 0..3u32 {
            assert_eq!(run_lm(&p, false, seed), run_lm(&p, true, seed));
        }
    }

    fn to_lp(p: &S5Problem) -> crate::lp_data::highs_lp::HighsLp {
        let mut lp = crate::lp_data::highs_lp::HighsLp::default();
        lp.num_col_ = p.ncol;
        lp.num_row_ = p.nrow;
        lp.col_cost_ = p.col_cost.clone();
        lp.col_lower_ = p.col_lower.clone();
        lp.col_upper_ = p.col_upper.clone();
        lp.row_lower_ = p.row_lower.clone();
        lp.row_upper_ = p.row_upper.clone();
        lp.integrality_ = p.integrality.clone();
        lp.a_matrix_.num_col_ = p.ncol;
        lp.a_matrix_.num_row_ = p.nrow;
        lp.a_matrix_.start_ = p.csc.col_start.clone();
        lp.a_matrix_.index_ = p.csc.col_row.clone();
        lp.a_matrix_.value_ = p.csc.col_val.clone();
        lp
    }

    /// Full solve: (objective, nodes, LP iterations, offer objectives).
    fn solve(
        lp: &crate::lp_data::highs_lp::HighsLp,
        lm: bool,
        gub: bool,
    ) -> (f64, i64, i64, Vec<f64>) {
        let mut options = crate::lp_data::highs_options::HighsOptions::default();
        options.output_flag = false;
        options.log_options.output_flag = false;
        options.threads = 1;
        options.mip_heuristic_run_local_mip = lm;
        options.mip_ls_gub_swap = gub;
        let solution = crate::lp_data::hstruct::HighsSolution::default();
        let mut callback = crate::lp_data::highs_callback::HighsCallback::default();
        let mut ms = crate::mip::mip_solver::HighsMipSolver::new(
            &mut callback,
            &options,
            lp,
            &solution,
            false,
            0,
        );
        ms.run();
        let offers = ms
            .mipdata_
            .s5_report
            .as_ref()
            .map(|r| r.offers.iter().map(|o| o.obj).collect())
            .unwrap_or_default();
        (
            ms.solution_objective_,
            ms.mipdata_.num_nodes,
            ms.mipdata_.total_lp_iterations,
            offers,
        )
    }

    #[test]
    fn ls_gub_e2e_flag_inert_without_local_mip_and_keeps_optimum() {
        let lp = to_lp(&gap_problem(3, 12, 1, 0.9));
        let base = solve(&lp, false, false);
        // The flag acts only through Local-MIP: bit-identical solve without it.
        assert_eq!(solve(&lp, false, true), base);
        let lm = solve(&lp, true, false);
        let lm_gub = solve(&lp, true, true);
        assert!((lm.0 - base.0).abs() < 1e-6, "lm {} base {}", lm.0, base.0);
        assert!(
            (lm_gub.0 - base.0).abs() < 1e-6,
            "lm+gub {} base {}",
            lm_gub.0,
            base.0
        );
    }
}
