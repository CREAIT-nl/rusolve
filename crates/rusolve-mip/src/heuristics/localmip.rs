//! Local-MIP: a weighted local search that improves a feasible MIP point.
//!
//! Every constraint is a `aᵀx ≤ b` row (one per finite side) and the
//! objective is treated as one more row, `cᵀx ≤ best − opt_tol`. While the
//! current point is feasible the search applies *lift moves*: single-variable
//! moves to the end of the variable's local feasible domain that improve the
//! objective without breaking any row. When no lift move is left it gives up
//! feasibility on purpose: *breakthrough moves* push the objective below the
//! incumbent and *mixed tight moves* make a violated (or a satisfied) row
//! tight; a weighted, two-level score (progress, then breakthrough bonus and
//! age) picks among a random sample of them, tabu forbids undoing a move for
//! a few steps, and constraint weights grow at local optima so the search
//! drifts back to feasibility. Every time it gets there with a better
//! objective the point is recorded and lifting resumes.
//!
//! Provenance:
//!
//! * P. Lin, M. Zou & S. Cai, "An Efficient Local Search Solver for Mixed
//!   Integer Programming", CP 2024, LIPIcs 307, 19:1–19:19,
//!   doi:10.4230/LIPIcs.CP.2024.19 — breakthrough move (Def. 4), mixed tight
//!   move (Def. 5), the two-level progress / bonus score (Defs. 7–10),
//!   weighting, the candidate-operation cascade (Alg. 2), tabu tenure.
//! * P. Lin, S. Cai, M. Zou & J. Lin, "Local Search for Integer Linear
//!   Programming", arXiv:2305.00188 (2023) — local domain reduction and the
//!   lift move (§5, Defs. 4–6), BMS sampling (§6), restart from the best point.
//! * Local-MIP reference code, github.com/shaowei-cai-group/Local-MIP at
//!   3bd2d611 (MIT, © 2025 Local-MIP), read as a reference text and
//!   reimplemented, not translated:
//!   - lift domain and lift score with age tie-break
//!     (`lift_move.cpp:98`, `lift_scoring.cpp:59`);
//!   - breakthrough and mixed-tight deltas with the out-of-bound fallbacks
//!     (`neighbor.cpp:126`, `:161`);
//!   - the progress score with ±2w / ±w terms and the bonus / age tie-breaks
//!     (`neighbor_scoring.cpp:82`);
//!   - the neighbourhood list unsat-mtm+bm → sat-mtm → flip → easy →
//!     random-unsat-mtm+bm, stopping at the first positive score, with the
//!     weight update right before the last (`Local_Search.cpp:121`, `:1034`);
//!   - `tabu` / `tabu_latest` (`neighbor.h:215`, `:226`);
//!   - monotone and smoothed weights (`weight.cpp:83`) and restart to the
//!     best point with weights reset (`restart.cpp:113`);
//!   - defaults from `default.set`: tabu 4 + rand(7), BMS sizes 10 / 2250 /
//!     1 / 80 / 0 / 5 / 250, restart after 10⁶ steps without improvement,
//!     opt tolerance 1e-4, equalities split into two inequalities.
//!
//! Deviations (deliberate):
//!
//! * No presolve or bound strengthening here (the reference has its own);
//!   the problem is searched as given.
//! * The run is bounded by an *effort* budget (≈ nonzeros touched, like
//!   [`super::fj`]) instead of wall-clock time, so results are reproducible.
//! * The paper's robustness bonus (Def. 10, `w_i` per strictly satisfied
//!   row) is left out, as in the reference code, whose bonus is the
//!   breakthrough term alone.
//! * Moves are applied to the start point instead of a zero start: this is
//!   an improver, so it begins in Improve mode with the start as incumbent.
//! * The returned point is re-checked with [`MipProblem::violation`] and the
//!   exact objective; on incremental drift past the tolerances the run
//!   reports nothing rather than an unverified point.

use crate::problem::MipProblem;
use rusolve_lp::rng::XorShift64;

const NONE: usize = usize::MAX;
/// Deltas at or below this are no move (reference `zero_tolerance`).
const ZERO_TOL: f64 = 1e-9;

/// Per-polish effort budget for the tree search, as a multiple of FJ's full
/// cap ([`super::fj::effort_cap`], the same effort currency): an
/// "about 3× FJ's" ratio (Spoorendonk's shipped 3.287 vs 0.3317).
pub const POLISH_EFFORT_MULT: usize = 3;
/// Effort the improver may spend over a whole search, in multiples of FJ's
/// cap. TODO(perf): tune with the `--localmip` A/B gate table.
pub const TOTAL_EFFORT_MULT: usize = 24;

#[derive(Debug, Clone, Copy)]
pub struct LocalMipConfig {
    /// A row counts as satisfied when `aᵀx − b ≤ feastol`.
    pub feastol: f64,
    /// A new point must beat the incumbent objective by this much.
    pub opt_tol: f64,
    /// A move forbids its reverse for `tabu_base + rand(tabu_variation)`
    /// steps.
    pub tabu_base: usize,
    pub tabu_variation: usize,
    /// Restart from the best point after this many steps without
    /// improvement (0 = never).
    pub restart_steps: usize,
    /// Probability of smoothing instead of bumping weights at a local
    /// optimum (0 = the reference's monotone default; the paper uses 3e-4).
    pub smooth_probability: f64,
    /// BMS sample sizes: violated rows and their ops, satisfied rows and
    /// their ops, binary flips, "easy" bound moves, ops of the random
    /// violated row.
    pub bms_unsat_cons: usize,
    pub bms_unsat_ops: usize,
    pub bms_sat_cons: usize,
    pub bms_sat_ops: usize,
    pub bms_flip_ops: usize,
    pub bms_easy_ops: usize,
    pub bms_random_ops: usize,
    /// Recompute all activities from scratch after this many moves.
    pub activity_period: usize,
}

impl Default for LocalMipConfig {
    fn default() -> Self {
        Self {
            feastol: 1e-6,
            opt_tol: 1e-4,
            tabu_base: 4,
            tabu_variation: 7,
            restart_steps: 1_000_000,
            smooth_probability: 0.0,
            bms_unsat_cons: 10,
            bms_unsat_ops: 2250,
            bms_sat_cons: 1,
            bms_sat_ops: 80,
            bms_flip_ops: 0,
            bms_easy_ops: 5,
            bms_random_ops: 250,
            activity_period: 100_000,
        }
    }
}

/// An improved feasible point and its (internal, minimization) objective.
#[derive(Debug, Clone)]
pub struct LocalMipResult {
    pub x: Vec<f64>,
    pub objective: f64,
    /// Effort units spent (≈ nonzeros touched).
    pub effort: usize,
    /// Search steps (moves tried) and improving points recorded.
    pub steps: usize,
    pub improvements: usize,
}

/// Improves the feasible point `start` of `p` by Local-MIP within
/// `effort_budget` effort units, with the default configuration. Returns a
/// feasible point with a strictly better objective, or `None` (no
/// improvement, or `start` itself is not feasible).
pub fn local_mip(
    p: &MipProblem,
    start: &[f64],
    effort_budget: usize,
    seed: u64,
) -> Option<LocalMipResult> {
    local_mip_with(p, start, effort_budget, seed, &LocalMipConfig::default())
}

/// [`local_mip`] with an explicit configuration (ablations).
pub fn local_mip_with(
    p: &MipProblem,
    start: &[f64],
    effort_budget: usize,
    seed: u64,
    cfg: &LocalMipConfig,
) -> Option<LocalMipResult> {
    let n = p.num_col;
    if start.len() != n || n == 0 {
        return None;
    }
    let v = p.violation(start);
    if !(v.bound <= cfg.feastol && v.row <= cfg.feastol && v.integrality <= cfg.feastol) {
        return None;
    }
    if p.cost.iter().all(|&c| c == 0.0) {
        return None;
    }
    let mut s = Search::new(p, start, seed, *cfg)?;
    s.run(effort_budget);
    if s.improvements == 0 {
        return None;
    }
    let x = s.best_x;
    let v = p.violation(&x);
    let objective = p.objective(&x);
    let ok = v.bound <= cfg.feastol && v.row <= cfg.feastol && v.integrality <= cfg.feastol;
    (ok && objective < p.objective(start)).then_some(LocalMipResult {
        x,
        objective,
        effort: s.effort,
        steps: s.step,
        improvements: s.improvements,
    })
}

/// Swap-remove set of indices with O(1) membership via a position array.
struct IndexSet {
    items: Vec<usize>,
    pos: Vec<usize>,
}

impl IndexSet {
    fn new(n: usize) -> Self {
        Self {
            items: Vec::new(),
            pos: vec![NONE; n],
        }
    }

    fn contains(&self, i: usize) -> bool {
        self.pos[i] != NONE
    }

    fn insert(&mut self, i: usize) {
        if self.pos[i] == NONE {
            self.pos[i] = self.items.len();
            self.items.push(i);
        }
    }

    fn remove(&mut self, i: usize) {
        let p = self.pos[i];
        if p == NONE {
            return;
        }
        let last = *self.items.last().expect("non-empty");
        self.items.swap_remove(p);
        if last != i {
            self.pos[last] = p;
        }
        self.pos[i] = NONE;
    }

    fn clear(&mut self) {
        for &i in &self.items {
            self.pos[i] = NONE;
        }
        self.items.clear();
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Neighbourhood {
    /// Mixed tight moves in sampled violated rows plus breakthrough moves.
    UnsatTightBreak,
    /// Mixed tight moves in sampled satisfied rows.
    SatTight,
    /// Binary flips.
    Flip,
    /// Variables to 0 / their bounds / their midpoint.
    Easy,
    /// All tight moves of one random violated row plus breakthrough moves,
    /// with only the latest move tabu; scored after a weight update and
    /// accepted even with a non-positive score.
    RandomUnsat,
}

const NEIGHBOURHOODS: [Neighbourhood; 5] = [
    Neighbourhood::UnsatTightBreak,
    Neighbourhood::SatTight,
    Neighbourhood::Flip,
    Neighbourhood::Easy,
    Neighbourhood::RandomUnsat,
];

struct Search {
    cfg: LocalMipConfig,
    // Variables.
    lb: Vec<f64>,
    ub: Vec<f64>,
    integer: Vec<bool>,
    cost: Vec<f64>,
    /// Columns with a nonzero cost, and each column's position there.
    obj_vars: Vec<usize>,
    obj_pos: Vec<usize>,
    binaries: Vec<usize>,
    movable: Vec<usize>,
    /// Column-wise (row, coefficient) cells.
    col_ptr: Vec<usize>,
    col_cells: Vec<(usize, f64)>,
    // Rows `aᵀx ≤ rhs`.
    rhs: Vec<f64>,
    row_ptr: Vec<usize>,
    row_cells: Vec<(usize, f64)>,
    // State.
    x: Vec<f64>,
    act: Vec<f64>,
    obj: f64,
    weight: Vec<u64>,
    obj_weight: u64,
    unsat: IndexSet,
    sat: IndexSet,
    best_x: Vec<f64>,
    best_obj: f64,
    /// `best_obj − opt_tol`: the objective row's right-hand side.
    cutoff: f64,
    step: usize,
    last_improve_step: usize,
    improvements: usize,
    allow_inc: Vec<usize>,
    allow_dec: Vec<usize>,
    last_inc: Vec<usize>,
    last_dec: Vec<usize>,
    /// Lift delta per objective column, valid while `keep_feas`.
    lift_delta: Vec<f64>,
    keep_feas: bool,
    moves_since_refresh: usize,
    // Best operation of the current selection.
    best_var: usize,
    best_delta: f64,
    best_score: i64,
    best_bonus: i64,
    best_age: usize,
    // Scratch.
    ops: Vec<(usize, f64)>,
    stamp: Vec<u32>,
    stamp_token: u32,
    picks: Vec<usize>,
    remap: Vec<(usize, usize)>,
    touched: Vec<usize>,
    touched_mark: Vec<bool>,
    rng: XorShift64,
    effort: usize,
}

impl Search {
    fn new(p: &MipProblem, start: &[f64], seed: u64, cfg: LocalMipConfig) -> Option<Self> {
        let n = p.num_col;
        let mut lb = p.col_lower.clone();
        let mut ub = p.col_upper.clone();
        let mut x = start.to_vec();
        for j in 0..n {
            if p.integer[j] {
                lb[j] = (lb[j] - cfg.feastol).ceil();
                ub[j] = (ub[j] + cfg.feastol).floor();
                x[j] = x[j].round();
            }
            if lb[j] > ub[j] {
                return None;
            }
            x[j] = x[j].clamp(lb[j], ub[j]);
        }

        // Rows: one `≤` row per finite side, `≥` sides negated.
        let csr = p.matrix.to_csr();
        let (rptr, ridx, rval) = (csr.row_ptr(), csr.col_idx(), csr.values());
        let mut rhs = Vec::new();
        let mut row_ptr = vec![0];
        let mut row_cells = Vec::new();
        for i in 0..p.num_row {
            let nz = || {
                (rptr[i]..rptr[i + 1])
                    .filter(|&k| rval[k] != 0.0)
                    .map(|k| (ridx[k] as usize, rval[k]))
            };
            if nz().next().is_none() {
                continue;
            }
            if p.row_lower[i].is_finite() {
                rhs.push(-p.row_lower[i]);
                row_cells.extend(nz().map(|(j, a)| (j, -a)));
                row_ptr.push(row_cells.len());
            }
            if p.row_upper[i].is_finite() {
                rhs.push(p.row_upper[i]);
                row_cells.extend(nz());
                row_ptr.push(row_cells.len());
            }
        }
        let m = rhs.len();
        let mut col_ptr = vec![0usize; n + 1];
        for &(j, _) in &row_cells {
            col_ptr[j + 1] += 1;
        }
        for j in 0..n {
            col_ptr[j + 1] += col_ptr[j];
        }
        let mut col_cells = vec![(0usize, 0.0); row_cells.len()];
        let mut fill = col_ptr.clone();
        for c in 0..m {
            for &(j, a) in &row_cells[row_ptr[c]..row_ptr[c + 1]] {
                col_cells[fill[j]] = (c, a);
                fill[j] += 1;
            }
        }

        let obj_vars: Vec<usize> = (0..n).filter(|&j| p.cost[j] != 0.0).collect();
        let mut obj_pos = vec![NONE; n];
        for (t, &j) in obj_vars.iter().enumerate() {
            obj_pos[j] = t;
        }
        let binaries = (0..n)
            .filter(|&j| p.integer[j] && lb[j] == 0.0 && ub[j] == 1.0)
            .collect();
        let movable = (0..n).filter(|&j| lb[j] < ub[j]).collect();
        let num_obj = obj_vars.len();
        let mut s = Self {
            cfg,
            lb,
            ub,
            integer: p.integer.clone(),
            cost: p.cost.clone(),
            obj_vars,
            obj_pos,
            binaries,
            movable,
            col_ptr,
            col_cells,
            rhs,
            row_ptr,
            row_cells,
            x,
            act: vec![0.0; m],
            obj: 0.0,
            weight: vec![1; m],
            obj_weight: 1,
            unsat: IndexSet::new(m),
            sat: IndexSet::new(m),
            best_x: Vec::new(),
            best_obj: f64::INFINITY,
            cutoff: f64::INFINITY,
            // Step 0 is reserved for "never moved" in the age arrays.
            step: 1,
            last_improve_step: 1,
            improvements: 0,
            allow_inc: vec![0; n],
            allow_dec: vec![0; n],
            last_inc: vec![0; n],
            last_dec: vec![0; n],
            lift_delta: vec![0.0; num_obj],
            keep_feas: false,
            moves_since_refresh: 0,
            best_var: NONE,
            best_delta: 0.0,
            best_score: 0,
            best_bonus: 0,
            best_age: 0,
            ops: Vec::new(),
            stamp: vec![0; n],
            stamp_token: 0,
            picks: Vec::new(),
            remap: Vec::new(),
            touched: Vec::new(),
            touched_mark: vec![false; n],
            rng: XorShift64::new(seed),
            effort: 0,
        };
        s.refresh();
        if !s.unsat.items.is_empty() {
            return None;
        }
        s.best_x = s.x.clone();
        s.best_obj = s.obj;
        s.cutoff = s.obj - cfg.opt_tol;
        Some(s)
    }

    /// Recomputes activities, the objective and the row status from `x`.
    fn refresh(&mut self) {
        self.act.iter_mut().for_each(|a| *a = 0.0);
        for (j, &v) in self.x.iter().enumerate() {
            for &(c, a) in &self.col_cells[self.col_ptr[j]..self.col_ptr[j + 1]] {
                self.act[c] += a * v;
            }
        }
        self.obj = self.cost.iter().zip(&self.x).map(|(c, v)| c * v).sum();
        self.unsat.clear();
        self.sat.clear();
        for c in 0..self.rhs.len() {
            if self.act[c] - self.rhs[c] <= self.cfg.feastol {
                self.sat.insert(c);
            } else {
                self.unsat.insert(c);
            }
        }
        self.moves_since_refresh = 0;
        self.effort += self.col_cells.len() + self.x.len();
    }

    fn run(&mut self, budget: usize) {
        while self.effort < budget {
            self.effort += 1;
            if self.cfg.restart_steps > 0
                && self.step > self.last_improve_step + self.cfg.restart_steps
            {
                self.restart();
            }
            if self.unsat.items.is_empty() {
                if self.obj <= self.cutoff {
                    self.record_best();
                }
                let lifted = self.lift_move();
                self.step += 1;
                if lifted {
                    continue;
                }
            }
            self.select_neighbour_move();
            if self.best_var != NONE {
                self.apply(self.best_var, self.best_delta);
            }
            self.keep_feas = false;
            self.step += 1;
        }
    }

    fn record_best(&mut self) {
        self.best_x.copy_from_slice(&self.x);
        self.best_obj = self.obj;
        self.cutoff = self.obj - self.cfg.opt_tol;
        self.last_improve_step = self.step;
        self.improvements += 1;
        self.effort += self.x.len();
    }

    /// Back to the best point with fresh weights and tabu state
    /// (reference `best_restart`).
    fn restart(&mut self) {
        self.x.copy_from_slice(&self.best_x);
        self.weight.iter_mut().for_each(|w| *w = 1);
        self.obj_weight = 1;
        for v in [
            &mut self.allow_inc,
            &mut self.allow_dec,
            &mut self.last_inc,
            &mut self.last_dec,
        ] {
            v.iter_mut().for_each(|s| *s = 0);
        }
        self.last_improve_step = self.step;
        self.keep_feas = false;
        self.refresh();
    }

    // ---- Lift moves (feasible point) ----

    /// The objective-improving end of `j`'s local feasible domain, as a
    /// delta: the interval where every row keeps its status (Local-ILP §5.1),
    /// intersected with the bounds; tight rows pin the variable.
    fn lift_delta_of(&mut self, j: usize) -> f64 {
        let xj = self.x[j];
        let (mut lo, mut hi) = (self.lb[j] - xj, self.ub[j] - xj);
        let range = self.col_ptr[j]..self.col_ptr[j + 1];
        self.effort += range.len();
        for &(c, a) in &self.col_cells[range] {
            let gap = self.act[c] - self.rhs[c];
            if gap >= 0.0 {
                if a > 0.0 {
                    hi = 0.0;
                } else {
                    lo = 0.0;
                }
            } else {
                let d = -gap / a;
                if a > 0.0 {
                    hi = hi.min(if self.integer[j] { d.floor() } else { d });
                } else {
                    lo = lo.max(if self.integer[j] { d.ceil() } else { d });
                }
            }
            if lo >= hi {
                break;
            }
        }
        if self.cost[j] > 0.0 { lo } else { hi }
    }

    /// Applies the best lift move; false when no objective column can move.
    fn lift_move(&mut self) -> bool {
        if !self.keep_feas {
            for t in 0..self.obj_vars.len() {
                self.lift_delta[t] = self.lift_delta_of(self.obj_vars[t]);
            }
        }
        let (mut best_var, mut best_delta) = (NONE, 0.0);
        let (mut best_score, mut best_age) = (0.0, usize::MAX);
        self.effort += self.obj_vars.len();
        for (t, &j) in self.obj_vars.iter().enumerate() {
            let d = self.lift_delta[t];
            if d.abs() <= ZERO_TOL || !d.is_finite() {
                continue;
            }
            let score = -self.cost[j] * d;
            let age = self.last_inc[j].max(self.last_dec[j]);
            if best_score + self.cfg.opt_tol < score || (best_score <= score && age < best_age) {
                (best_var, best_delta, best_score, best_age) = (j, d, score, age);
            }
        }
        if best_var == NONE {
            self.keep_feas = false;
            return false;
        }
        self.apply(best_var, best_delta);
        self.keep_feas = true;
        // Only objective columns sharing a row with the moved one change.
        let mut touched = std::mem::take(&mut self.touched);
        touched.push(best_var);
        self.touched_mark[best_var] = true;
        for k in self.col_ptr[best_var]..self.col_ptr[best_var + 1] {
            let c = self.col_cells[k].0;
            self.effort += self.row_ptr[c + 1] - self.row_ptr[c];
            for &(j, _) in &self.row_cells[self.row_ptr[c]..self.row_ptr[c + 1]] {
                if !self.touched_mark[j] {
                    self.touched_mark[j] = true;
                    touched.push(j);
                }
            }
        }
        for &j in &touched {
            self.touched_mark[j] = false;
            let t = self.obj_pos[j];
            if t != NONE {
                self.lift_delta[t] = self.lift_delta_of(j);
            }
        }
        touched.clear();
        self.touched = touched;
        true
    }

    // ---- Neighbourhood moves (infeasible point, or no lift left) ----

    fn reset_best(&mut self, require_positive: bool) {
        self.stamp_token = self.stamp_token.wrapping_add(1);
        if self.stamp_token == 0 {
            self.stamp.iter_mut().for_each(|s| *s = 0);
            self.stamp_token = 1;
        }
        self.best_var = NONE;
        self.best_delta = 0.0;
        self.best_score = if require_positive { 0 } else { i64::MIN };
        self.best_bonus = i64::MIN;
        self.best_age = usize::MAX;
    }

    /// Walks the neighbourhood list until one yields a positive score; the
    /// last one runs after a weight update and takes its best op whatever
    /// the score (Local-MIP Alg. 2).
    fn select_neighbour_move(&mut self) {
        self.reset_best(true);
        for (k, &nb) in NEIGHBOURHOODS.iter().enumerate() {
            self.ops.clear();
            if k + 1 == NEIGHBOURHOODS.len() {
                self.reset_best(false);
                self.update_weights();
            }
            self.collect(nb);
            for i in 0..self.ops.len() {
                let (j, d) = self.ops[i];
                self.score(j, d);
            }
            if self.best_score > 0 {
                break;
            }
        }
    }

    fn update_weights(&mut self) {
        let sp = self.cfg.smooth_probability;
        if sp > 0.0 && self.rng.next_f64() < sp {
            self.effort += self.sat.items.len();
            for &c in &self.sat.items {
                self.weight[c] = self.weight[c].saturating_sub(1);
            }
            if self.obj <= self.cutoff {
                self.obj_weight = self.obj_weight.saturating_sub(1);
            }
        } else {
            self.effort += self.unsat.items.len();
            for &c in &self.unsat.items {
                self.weight[c] += 1;
            }
            if self.unsat.items.is_empty() {
                self.obj_weight += 1;
            }
        }
    }

    fn tabu(&self, j: usize, d: f64) -> bool {
        (d < 0.0 && self.step < self.allow_dec[j]) || (d > 0.0 && self.step < self.allow_inc[j])
    }

    /// Only the immediate reversal of the previous step is forbidden.
    fn tabu_latest(&self, j: usize, d: f64) -> bool {
        (d < 0.0 && self.last_inc[j] > 0 && self.step == self.last_inc[j] + 1)
            || (d > 0.0 && self.last_dec[j] > 0 && self.step == self.last_dec[j] + 1)
    }

    fn push_op(&mut self, j: usize, d: f64, latest_only: bool) {
        let tabu = if latest_only {
            self.tabu_latest(j, d)
        } else {
            self.tabu(j, d)
        };
        if d.abs() > ZERO_TOL && d.is_finite() && !tabu {
            self.ops.push((j, d));
        }
    }

    fn out_of_bounds(&self, j: usize, v: f64) -> bool {
        v < self.lb[j] || v > self.ub[j]
    }

    /// Mixed tight move of `x_j` in row `c` (coefficient `a`): the delta
    /// that makes the row tight, rounded towards keeping it satisfied; past a
    /// bound, the bound in the direction that helps a violated row (or that
    /// uses up a satisfied row's slack). Local-MIP Def. 5.
    fn tight_delta(&self, c: usize, j: usize, a: f64) -> f64 {
        let gap = self.act[c] - self.rhs[c];
        let mut d = -gap / a;
        if self.integer[j] {
            d = if a > 0.0 { d.floor() } else { d.ceil() };
        }
        if self.out_of_bounds(j, self.x[j] + d) {
            let towards_lower = (a > 0.0) == self.unsat.contains(c);
            d = if towards_lower {
                self.lb[j]
            } else {
                self.ub[j]
            } - self.x[j];
        }
        d
    }

    /// Breakthrough move: the delta that takes the objective just below the
    /// incumbent's cutoff, or as far as the bound allows. Local-MIP Def. 4.
    fn breakthrough_delta(&self, j: usize) -> f64 {
        let c = self.cost[j];
        let mut d = -(self.obj - self.cutoff) / c;
        if self.integer[j] {
            d = if c > 0.0 { d.floor() } else { d.ceil() };
        }
        if self.out_of_bounds(j, self.x[j] + d) {
            d = if c > 0.0 { self.lb[j] } else { self.ub[j] } - self.x[j];
        }
        d
    }

    fn push_row_tight_ops(&mut self, c: usize, latest_only: bool) {
        self.effort += self.row_ptr[c + 1] - self.row_ptr[c];
        for k in self.row_ptr[c]..self.row_ptr[c + 1] {
            let (j, a) = self.row_cells[k];
            let d = self.tight_delta(c, j, a);
            self.push_op(j, d, latest_only);
        }
    }

    fn push_breakthrough_ops(&mut self, latest_only: bool) {
        if self.obj <= self.cutoff {
            return;
        }
        self.effort += self.obj_vars.len();
        for t in 0..self.obj_vars.len() {
            let j = self.obj_vars[t];
            let d = self.breakthrough_delta(j);
            self.push_op(j, d, latest_only);
        }
    }

    /// Fills `picks` with up to `k` distinct uniform positions in `0..len`
    /// (sparse partial Fisher–Yates; BMS, Cai 2015).
    fn sample_positions(&mut self, len: usize, k: usize) {
        self.picks.clear();
        if len <= k {
            self.picks.extend(0..len);
            return;
        }
        self.remap.clear();
        let lookup =
            |remap: &[(usize, usize)], i: usize| remap.iter().find(|e| e.0 == i).map_or(i, |e| e.1);
        for s in 0..k {
            let avail = len - s;
            let r = (self.rng.next_u64() % avail as u64) as usize;
            let pick = lookup(&self.remap, r);
            let last = lookup(&self.remap, avail - 1);
            match self.remap.iter_mut().find(|e| e.0 == r) {
                Some(e) => e.1 = last,
                None => self.remap.push((r, last)),
            }
            self.picks.push(pick);
        }
    }

    /// Keeps a uniform sample of at most `k` collected ops.
    fn sample_ops(&mut self, k: usize) {
        let len = self.ops.len();
        if len <= k {
            return;
        }
        for i in 0..k {
            let r = i + (self.rng.next_u64() % (len - i) as u64) as usize;
            self.ops.swap(i, r);
        }
        self.ops.truncate(k);
    }

    fn collect(&mut self, nb: Neighbourhood) {
        let cfg = self.cfg;
        match nb {
            Neighbourhood::UnsatTightBreak => {
                if cfg.bms_unsat_cons == 0 || cfg.bms_unsat_ops == 0 {
                    return;
                }
                self.sample_positions(self.unsat.items.len(), cfg.bms_unsat_cons);
                for i in 0..self.picks.len() {
                    let c = self.unsat.items[self.picks[i]];
                    self.push_row_tight_ops(c, false);
                }
                self.push_breakthrough_ops(false);
                self.sample_ops(cfg.bms_unsat_ops);
            }
            Neighbourhood::SatTight => {
                if cfg.bms_sat_cons == 0 || cfg.bms_sat_ops == 0 {
                    return;
                }
                self.sample_positions(self.sat.items.len(), cfg.bms_sat_cons);
                for i in 0..self.picks.len() {
                    let c = self.sat.items[self.picks[i]];
                    self.push_row_tight_ops(c, false);
                }
                self.sample_ops(cfg.bms_sat_ops);
            }
            Neighbourhood::Flip => {
                self.sample_positions(self.binaries.len(), cfg.bms_flip_ops);
                for i in 0..self.picks.len() {
                    let j = self.binaries[self.picks[i]];
                    let d = if self.x[j] > 0.5 { -1.0 } else { 1.0 };
                    self.push_op(j, d, false);
                }
            }
            Neighbourhood::Easy => {
                self.sample_positions(self.movable.len(), cfg.bms_easy_ops);
                for i in 0..self.picks.len() {
                    let j = self.movable[self.picks[i]];
                    let (l, u, xj) = (self.lb[j], self.ub[j], self.x[j]);
                    let zero = if l > 0.0 {
                        l
                    } else if u < 0.0 {
                        u
                    } else {
                        0.0
                    };
                    self.push_easy(j, zero - xj);
                    if !self.integer[j] && l.is_finite() && u.is_finite() {
                        self.push_easy(j, 0.5 * (l + u) - xj);
                    }
                    if l.is_finite() && l < 0.0 {
                        self.push_easy(j, l - xj);
                    }
                    if u.is_finite() && u > 0.0 {
                        self.push_easy(j, u - xj);
                    }
                }
            }
            Neighbourhood::RandomUnsat => {
                if !self.unsat.items.is_empty() {
                    let r = (self.rng.next_u64() % self.unsat.items.len() as u64) as usize;
                    let c = self.unsat.items[r];
                    self.push_row_tight_ops(c, true);
                }
                self.push_breakthrough_ops(true);
                self.sample_ops(cfg.bms_random_ops);
            }
        }
    }

    fn push_easy(&mut self, j: usize, d: f64) {
        if d.abs() > self.cfg.feastol && !self.tabu(j, d) {
            self.ops.push((j, d));
        }
    }

    /// Two-level score of moving `x_j` by `d` (Local-MIP Defs. 7–10 as in
    /// the reference `progress`): per row +2w when it becomes satisfied, −2w
    /// when it becomes violated, ±w when a violated row gets less / more
    /// violated; ±2w_obj for the objective direction. Ties go to a larger
    /// breakthrough bonus (w_obj if the new objective beats the incumbent),
    /// then to the variable moved longest ago.
    fn score(&mut self, j: usize, d: f64) {
        if self.integer[j] && self.lb[j] == 0.0 && self.ub[j] == 1.0 {
            // A binary has a single move; score it once per selection.
            if self.stamp[j] == self.stamp_token {
                return;
            }
            self.stamp[j] = self.stamp_token;
        }
        let range = self.col_ptr[j]..self.col_ptr[j + 1];
        if range.is_empty() && self.cost[j] == 0.0 {
            return;
        }
        self.effort += range.len() + 1;
        let tol = self.cfg.feastol;
        let (mut score, mut bonus) = (0i64, 0i64);
        if self.cost[j] != 0.0 {
            let w = self.obj_weight as i64;
            let new_obj = self.obj + self.cost[j] * d;
            score += if new_obj < self.obj { 2 * w } else { -2 * w };
            if new_obj < self.best_obj {
                bonus += w;
            }
        }
        for &(c, a) in &self.col_cells[range] {
            let w = self.weight[c] as i64;
            let pre = self.act[c] - self.rhs[c];
            let new = pre + a * d;
            score += match (pre <= tol, new <= tol) {
                (false, true) => 2 * w,
                (true, false) => -2 * w,
                (false, false) if new < pre => w,
                (false, false) => -w,
                (true, true) => 0,
            };
        }
        let age = self.last_inc[j].max(self.last_dec[j]);
        let better = score > self.best_score
            || (score == self.best_score
                && (bonus > self.best_bonus || (bonus == self.best_bonus && age < self.best_age)));
        if better {
            self.best_var = j;
            self.best_delta = d;
            self.best_score = score;
            self.best_bonus = bonus;
            self.best_age = age;
        }
    }

    /// Moves `x_j` by `d` (clamped to the bounds), updates activities, row
    /// status and tabu state.
    fn apply(&mut self, j: usize, d: f64) {
        let mut v = (self.x[j] + d).clamp(self.lb[j], self.ub[j]);
        if self.integer[j] {
            v = v.round();
        }
        let d = v - self.x[j];
        if d == 0.0 {
            return;
        }
        self.x[j] = v;
        self.obj += self.cost[j] * d;
        let tol = self.cfg.feastol;
        let range = self.col_ptr[j]..self.col_ptr[j + 1];
        self.effort += range.len();
        for k in range {
            let (c, a) = self.col_cells[k];
            self.act[c] += a * d;
            if self.act[c] - self.rhs[c] <= tol {
                if self.unsat.contains(c) {
                    self.unsat.remove(c);
                    self.sat.insert(c);
                }
            } else if self.sat.contains(c) {
                self.sat.remove(c);
                self.unsat.insert(c);
            }
        }
        let tenure = self.cfg.tabu_base
            + (self.rng.next_u64() % self.cfg.tabu_variation.max(1) as u64) as usize;
        if d > 0.0 {
            self.last_inc[j] = self.step;
            self.allow_dec[j] = self.step + tenure;
        } else {
            self.last_dec[j] = self.step;
            self.allow_inc[j] = self.step + tenure;
        }
        self.moves_since_refresh += 1;
        if self.moves_since_refresh >= self.cfg.activity_period {
            self.refresh();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bnb::{MipConfig, MipStatus, solve};
    use rusolve_core::CscMatrix;

    fn mip(
        rows: &[(&[f64], f64, f64)],
        lower: &[f64],
        upper: &[f64],
        integer: &[bool],
        cost: &[f64],
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
            cost: cost.to_vec(),
            maximize: false,
            col_lower: lower.to_vec(),
            col_upper: upper.to_vec(),
            integer: integer.to_vec(),
        }
    }

    fn assert_feasible(p: &MipProblem, x: &[f64]) {
        let v = p.violation(x);
        assert!(
            v.bound <= 1e-9 && v.row <= 1e-6 && v.integrality <= 1e-9,
            "{v:?} at {x:?}"
        );
    }

    #[test]
    fn lift_moves_alone_reach_the_knapsack_optimum_from_zero() {
        // max 5a + 4b + 3c s.t. 2a + 3b + c ≤ 5, binaries: optimum a = c = 1
        // and b = 0 is one lift away per variable from the all-zero start.
        let p = mip(
            &[(&[2.0, 3.0, 1.0], f64::NEG_INFINITY, 5.0)],
            &[0.0; 3],
            &[1.0; 3],
            &[true; 3],
            &[-5.0, -4.0, -3.0],
        );
        let r = local_mip(&p, &[0.0; 3], 10_000, 1).expect("improves");
        assert_feasible(&p, &r.x);
        assert_eq!(r.objective, -9.0, "{:?}", r.x);
    }

    #[test]
    fn breakthrough_escapes_a_lift_local_optimum() {
        // max 3x + 2y s.t. x + y ≤ 4, x ≤ 3 (row), x, y ∈ [0, 4] integer.
        // Start (0, 4): y is at its bound and the shared row is tight, so no
        // lift exists; the optimum (3, 1) needs moves through infeasibility.
        let inf = f64::INFINITY;
        let p = mip(
            &[(&[1.0, 1.0], -inf, 4.0), (&[1.0, 0.0], -inf, 3.0)],
            &[0.0, 0.0],
            &[4.0, 4.0],
            &[true, true],
            &[-3.0, -2.0],
        );
        let r = local_mip(&p, &[0.0, 4.0], 100_000, 3).expect("improves");
        assert_feasible(&p, &r.x);
        assert_eq!(r.objective, -11.0, "{:?}", r.x);
    }

    #[test]
    fn continuous_columns_and_equalities() {
        // min −x − 2z s.t. x + y = 3.5 (y continuous ∈ [0, 2]), z − x ≤ 0.5,
        // x, z ∈ [0, 5] integer: optimum x = 3, z = 3, y = 0.5 → −9.
        let inf = f64::INFINITY;
        let p = mip(
            &[(&[1.0, 1.0, 0.0], 3.5, 3.5), (&[-1.0, 0.0, 1.0], -inf, 0.5)],
            &[0.0, 0.0, 0.0],
            &[5.0, 2.0, 5.0],
            &[true, false, true],
            &[-1.0, 0.0, -2.0],
        );
        let start = [2.0, 1.5, 0.0];
        assert_feasible(&p, &start);
        let r = local_mip(&p, &start, 200_000, 5).expect("improves");
        assert_feasible(&p, &r.x);
        assert!(r.objective <= -9.0 + 1e-9, "{} at {:?}", r.objective, r.x);
    }

    #[test]
    fn infeasible_start_or_optimal_start_gives_none() {
        let p = mip(
            &[(&[1.0, 1.0], 1.0, f64::INFINITY)],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[true, true],
            &[1.0, 2.0],
        );
        assert!(local_mip(&p, &[0.0, 0.0], 10_000, 0).is_none());
        assert!(local_mip(&p, &[1.0, 0.0], 10_000, 0).is_none());
        let r = local_mip(&p, &[1.0, 1.0], 10_000, 0).expect("drops y");
        assert_eq!(r.objective, 1.0);
    }

    /// Random small MIP: a pure-integer box like the bnb tests, with a few
    /// columns turned continuous on request.
    fn random_mip(rng: &mut XorShift64, n: usize, m: usize, continuous: bool) -> MipProblem {
        let pick =
            |rng: &mut XorShift64, xs: &[f64]| xs[(rng.next_u64() % xs.len() as u64) as usize];
        let mut trip = Vec::new();
        for i in 0..m {
            for j in 0..n {
                let a = pick(rng, &[0.0, 0.0, 1.0, 2.0, -1.0, 3.0, -2.0, 0.5]);
                if a != 0.0 {
                    trip.push((i, j, a));
                }
            }
        }
        let inf = f64::INFINITY;
        let (mut row_lower, mut row_upper) = (Vec::new(), Vec::new());
        for _ in 0..m {
            let b = pick(rng, &[-2.0, 0.0, 1.0, 3.0, 4.5, 6.0]);
            let (lo, hi) = match rng.next_u64() % 4 {
                0 => (-inf, b),
                1 => (b, inf),
                2 => (b, b),
                _ => (b - 2.0, b + 1.5),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        let col_lower: Vec<f64> = (0..n).map(|_| pick(rng, &[-2.0, -1.0, 0.0, 0.0])).collect();
        let col_upper = col_lower
            .iter()
            .map(|&l| l + pick(rng, &[1.0, 2.0, 3.0, 4.0]))
            .collect();
        let integer = (0..n)
            .map(|_| !continuous || !rng.next_u64().is_multiple_of(3))
            .collect();
        MipProblem {
            num_col: n,
            num_row: m,
            matrix: CscMatrix::from_triplets(m, n, &trip),
            row_lower,
            row_upper,
            cost: (0..n)
                .map(|_| pick(rng, &[-3.0, -1.0, 0.0, 1.0, 2.0, 1.5]))
                .collect(),
            maximize: false,
            col_lower,
            col_upper,
            integer,
        }
    }

    /// All feasible points of a pure-integer box, worst objective first.
    fn feasible_points(p: &MipProblem) -> Vec<Vec<f64>> {
        let mut x: Vec<f64> = p.col_lower.clone();
        let mut out = Vec::new();
        loop {
            if p.violation(&x).row <= 1e-9 {
                out.push(x.clone());
            }
            let mut j = 0;
            loop {
                if j == p.num_col {
                    out.sort_by(|a, b| p.objective(b).total_cmp(&p.objective(a)));
                    return out;
                }
                if x[j] < p.col_upper[j] {
                    x[j] += 1.0;
                    break;
                }
                x[j] = p.col_lower[j];
                j += 1;
            }
        }
    }

    #[test]
    fn random_pure_integer_against_enumeration_and_bnb() {
        // From the worst and a random feasible point: the result is feasible,
        // strictly better than the start when returned, never better than
        // the enumerated optimum (which bnb must also match), and reaches
        // the optimum on most instances.
        let mut rng = XorShift64::new(20261022);
        let (mut runs, mut optimal, mut improvable) = (0, 0, 0);
        for case in 0..300 {
            let n = 2 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_mip(&mut rng, n, m, false);
            let points = feasible_points(&p);
            if points.is_empty() {
                continue;
            }
            let z = p.objective(points.last().unwrap());
            let sol = solve(&p, &MipConfig::default()).unwrap();
            assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
            assert!((sol.objective - z).abs() <= 1e-6, "case {case}");
            let random = (rng.next_u64() % points.len() as u64) as usize;
            for start in [&points[0], &points[random]] {
                runs += 1;
                let z0 = p.objective(start);
                improvable += usize::from(z0 > z + 1e-9);
                match local_mip(&p, start, 20_000, case) {
                    Some(r) => {
                        assert_feasible(&p, &r.x);
                        assert!((p.objective(&r.x) - r.objective).abs() <= 1e-9);
                        assert!(r.objective < z0, "case {case}: {} vs {z0}", r.objective);
                        assert!(r.objective >= z - 1e-9, "case {case}: below optimum");
                        optimal += usize::from(r.objective <= z + 1e-9);
                    }
                    None => optimal += usize::from(z0 <= z + 1e-9),
                }
            }
        }
        assert!(improvable > 200, "{improvable}/{runs}");
        assert!(optimal * 20 >= runs * 19, "optimal {optimal}/{runs}");
    }

    #[test]
    fn random_mixed_against_bnb() {
        // Continuous columns: bnb gives the optimum; the start is bnb's
        // point with the objective sign flipped (a feasible, usually poor
        // point). Without an LP the continuous part is only searched by
        // tight moves, so the optimum is reached less often than in the
        // pure-integer case.
        let mut rng = XorShift64::new(7);
        let (mut runs, mut improved, mut optimal) = (0, 0, 0);
        for case in 0..300 {
            let n = 2 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 4) as usize;
            let p = random_mip(&mut rng, n, m, true);
            let sol = solve(&p, &MipConfig::default()).unwrap();
            if sol.status != MipStatus::Optimal {
                continue;
            }
            let worst = MipProblem {
                cost: p.cost.iter().map(|c| -c).collect(),
                ..p.clone()
            };
            let w = solve(&worst, &MipConfig::default()).unwrap();
            if w.status != MipStatus::Optimal || p.violation(&w.x).row > 1e-6 {
                continue;
            }
            runs += 1;
            let z0 = p.objective(&w.x);
            if let Some(r) = local_mip(&p, &w.x, 50_000, case) {
                improved += 1;
                assert_feasible(&p, &r.x);
                assert!(r.objective < z0, "case {case}");
                assert!(
                    r.objective >= sol.objective - 1e-6,
                    "case {case}: below optimum"
                );
                optimal += usize::from(r.objective <= sol.objective + 1e-6);
            }
        }
        assert!(runs > 100 && improved * 2 > runs, "{improved}/{runs}");
        assert!(optimal * 4 >= runs * 3, "optimal {optimal}/{runs}");
    }
}
