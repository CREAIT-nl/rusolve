//! Feasibility Jump: an LP-free weighted local search for a first incumbent.
//!
//! Every variable owns one candidate move, its *jump value*: the value that
//! minimizes the weighted violation of its constraints with all other
//! variables held fixed. The search repeatedly applies the best-scoring move
//! among a small random sample of variables with a positive score; at a local
//! minimum it bumps the weights of all violated constraints (or the objective
//! weight, once feasible) — a Lagrangian-style penalty update — and restarts
//! from a variable of a random violated row.
//!
//! Provenance:
//!
//! * B. Luteberget & G. Sartor, "Feasibility Jump: an LP-free Lagrangian MIP
//!   heuristic", Math. Prog. Comp. 15:365–388 (2023) — the algorithm.
//! * HiGHS `mip/feasibilityjump.hh` (MIT; HiGHS' vendored copy of the authors'
//!   reference code) read as a reference text and reimplemented, not
//!   translated: constraint score `-max(0, violation)` (:75), incremental
//!   move-score updates on LHS changes (`modifyMove` :299), the jump-value
//!   sweep over sorted breakpoints (`JumpMove::updateValue` :327),
//!   variable selection constants — random good var with p = 0.001, sample of
//!   25 good vars, random cell of a violated row with p = 0.01 (:461–470),
//!   weight bump / 1e20 rescale (`updateWeights` :640), relative objective
//!   improvement 1e-4 to record a new solution (:475, :534).
//! * HiGHS `mip/HighsFeasibilityJump.cpp` (MIT): problem setup — integer
//!   bounds rounded inward with feastol slack (:51), start at the finite lower
//!   bound else the finite upper else 0 (:70), a `≥` and a `≤` constraint per
//!   finite row side (:86), effort limits `nnz·2¹⁰` total and `nnz·2⁸` since
//!   the last improvement (:107), equality tolerance = HiGHS `epsilon` (1e-9),
//!   violation tolerance = feastol (:33).
//!
//! Deviations (deliberate):
//!
//! * The jump-value sweep swaps the valid-range endpoints for negative
//!   coefficients. The HiGHS copy computes `[(lo−r)/a, (hi−r)/a]` without the
//!   swap (feasibilityjump.hh:362), which makes every row with `a < 0` an
//!   empty range there; the paper defines the jump value over all rows.
//! * Effort limits are checked every step instead of every 500 000 effort
//!   units (HiGHS `CALLBACK_EFFORT` :774), so tiny models stop at their
//!   `nnz`-scaled budget instead of a fixed floor.
//! * The objective term of a move score is `−w_obj·c_j·Δ` (cheaper ⇒
//!   positive, like the `−violation` row terms). The HiGHS copy adds
//!   `+w_obj·c_j·Δ` (feasibilityjump.hh:654, :755), which rewards *more*
//!   expensive moves once feasible. Both bugs also reported by Spoorendonk,
//!   arXiv 2609.22938.

use crate::problem::MipProblem;
use rusolve_lp::rng::XorShift64;

const RANDOM_VAR_PROBABILITY: f64 = 0.001;
const RANDOM_CELL_PROBABILITY: f64 = 0.01;
const MAX_MOVES_TO_EVALUATE: usize = 25;
const MIN_RELATIVE_OBJECTIVE_IMPROVEMENT: f64 = 1e-4;
const WEIGHT_RESCALE_THRESHOLD: f64 = 1e20;
const NONE: usize = usize::MAX;

#[derive(Debug, Clone, Copy)]
pub struct FjConfig {
    /// Integer-snapping and "same value" tolerance (HiGHS `epsilon`).
    pub equality_tol: f64,
    /// A row counts as violated when its violation exceeds this.
    pub violation_tol: f64,
    pub seed: u64,
    /// Weight increments grow by `1/decay` per bump (1.0 = constant).
    pub weight_update_decay: f64,
    /// Total effort limit; `None` = HiGHS's `nnz·2¹⁰` ([`effort_cap`]).
    pub max_effort: Option<usize>,
    /// Wall-clock stop (the MIP time limit). Polled only every
    /// [`DEADLINE_EFFORT`] effort units, like HiGHS's FJ termination
    /// callback (feasibilityjump.hh:774, `CALLBACK_EFFORT = 500000`), so it
    /// costs nothing in the move loop and never fires before the limit.
    pub deadline: Option<std::time::Instant>,
}

/// Effort units between two deadline polls (HiGHS `CALLBACK_EFFORT`).
const DEADLINE_EFFORT: usize = 500_000;

fn past(deadline: Option<std::time::Instant>) -> bool {
    deadline.is_some_and(|d| std::time::Instant::now() >= d)
}

impl FjConfig {
    pub fn new(feastol: f64) -> Self {
        Self {
            equality_tol: 1e-9,
            violation_tol: feastol,
            seed: 0,
            weight_update_decay: 1.0,
            max_effort: None,
            deadline: None,
        }
    }
}

/// The best feasible point the search found and its (internal) objective.
#[derive(Debug, Clone)]
pub struct FjResult {
    pub x: Vec<f64>,
    pub objective: f64,
    /// Effort units spent (≈ nonzeros touched).
    pub effort: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Sense {
    /// `aᵀx ≤ rhs`
    Le,
    /// `aᵀx ≥ rhs`
    Ge,
}

#[derive(Clone, Copy)]
struct Move {
    value: f64,
    score: f64,
}

struct Constraint {
    sense: Sense,
    rhs: f64,
    weight: f64,
    lhs: f64,
}

impl Constraint {
    /// 0 if `lhs` satisfies the row, `−violation` otherwise.
    fn score(&self, lhs: f64) -> f64 {
        match self.sense {
            Sense::Le => -(lhs - self.rhs).max(0.0),
            Sense::Ge => -(self.rhs - lhs).max(0.0),
        }
    }
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
}

struct Solver {
    cfg: FjConfig,
    // Variables.
    lb: Vec<f64>,
    ub: Vec<f64>,
    integer: Vec<bool>,
    obj: Vec<f64>,
    /// Column-wise (constraint, coefficient) cells.
    col_ptr: Vec<usize>,
    col_cells: Vec<(usize, f64)>,
    // Constraints.
    cons: Vec<Constraint>,
    row_ptr: Vec<usize>,
    row_cells: Vec<(usize, f64)>,
    // State.
    x: Vec<f64>,
    objective: f64,
    moves: Vec<Move>,
    violated: IndexSet,
    good: IndexSet,
    objective_weight: f64,
    weight_increment: f64,
    rng: XorShift64,
    effort: usize,
    /// Scratch for the jump-value sweep: (breakpoint, slope change).
    breakpoints: Vec<(f64, f64)>,
}

/// Runs Feasibility Jump on `p` restricted to the box `lower..upper`
/// (integer bounds are rounded inward here). `start` overrides the initial
/// assignment (clamped to the box). Returns the best feasible point found
/// within the effort budget, or `None`.
pub fn feasibility_jump(
    p: &MipProblem,
    lower: &[f64],
    upper: &[f64],
    start: Option<&[f64]>,
    cfg: FjConfig,
) -> Option<FjResult> {
    let n = p.num_col;
    let mut lb = lower.to_vec();
    let mut ub = upper.to_vec();
    for j in 0..n {
        if p.integer[j] {
            lb[j] = (lb[j] - cfg.violation_tol).ceil();
            ub[j] = (ub[j] + cfg.violation_tol).floor();
        }
        let legal = lb[j] <= ub[j] && lb[j] < f64::INFINITY && ub[j] > f64::NEG_INFINITY;
        if !legal {
            return None;
        }
    }
    let x0: Vec<f64> = (0..n)
        .map(|j| match start {
            Some(s) if s[j].is_finite() => s[j].clamp(lb[j], ub[j]),
            _ if lb[j].is_finite() => lb[j],
            _ if ub[j].is_finite() => ub[j],
            _ => 0.0,
        })
        .collect();

    // Rows: one `≥` and/or one `≤` constraint per finite side.
    let csr = p.matrix.to_csr();
    let (rptr, ridx, rval) = (csr.row_ptr(), csr.col_idx(), csr.values());
    let mut cons = Vec::new();
    let mut row_ptr = vec![0];
    let mut row_cells = Vec::new();
    for i in 0..p.num_row {
        let cells: Vec<(usize, f64)> = (rptr[i]..rptr[i + 1])
            .filter(|&k| rval[k] != 0.0)
            .map(|k| (ridx[k] as usize, rval[k]))
            .collect();
        if cells.is_empty() {
            // Constant row: nothing to search; the caller's check decides.
            continue;
        }
        for (sense, rhs) in [(Sense::Ge, p.row_lower[i]), (Sense::Le, p.row_upper[i])] {
            if rhs.is_finite() {
                cons.push(Constraint {
                    sense,
                    rhs,
                    weight: 1.0,
                    lhs: 0.0,
                });
                row_cells.extend_from_slice(&cells);
                row_ptr.push(row_cells.len());
            }
        }
    }
    let max_total = cfg.max_effort.unwrap_or_else(|| effort_cap(p));
    let nnz = rval.iter().filter(|v| **v != 0.0).count();
    // Column view of the constraint cells.
    let mut count = vec![0usize; n + 1];
    for &(j, _) in &row_cells {
        count[j + 1] += 1;
    }
    for j in 0..n {
        count[j + 1] += count[j];
    }
    let col_ptr = count.clone();
    let mut col_cells = vec![(0usize, 0.0); row_cells.len()];
    let mut fill = count;
    for c in 0..cons.len() {
        for &(j, a) in &row_cells[row_ptr[c]..row_ptr[c + 1]] {
            col_cells[fill[j]] = (c, a);
            fill[j] += 1;
        }
    }

    let num_cons = cons.len();
    let mut s = Solver {
        cfg,
        lb,
        ub,
        integer: p.integer.clone(),
        obj: p.cost.clone(),
        col_ptr,
        col_cells,
        cons,
        row_ptr,
        row_cells,
        x: x0,
        objective: 0.0,
        moves: vec![
            Move {
                value: 0.0,
                score: 0.0
            };
            n
        ],
        violated: IndexSet::new(num_cons),
        good: IndexSet::new(n),
        objective_weight: 0.0,
        weight_increment: 1.0,
        rng: XorShift64::new(cfg.seed),
        effort: 0,
        breakpoints: Vec::new(),
    };
    let max_since_improvement = nnz.max(1) << 8;
    s.run(max_total, max_since_improvement)
}

/// HiGHS's default total FJ effort, `nnz·2¹⁰` (HighsFeasibilityJump.cpp:107).
pub fn effort_cap(p: &MipProblem) -> usize {
    p.matrix
        .values()
        .iter()
        .filter(|v| **v != 0.0)
        .count()
        .max(1)
        << 10
}

/// Restarts Feasibility Jump from its best point `best` (fresh weights, next
/// seed after `cfg.seed`) while `restart_effort · effort_cap` effort units
/// remain; returns the best point over all runs, `effort` = restarts' total.
///
/// Measured on our benchmark set: at equal effort, restarting FJ from its best point closes 3–6% of the gap
/// to best-known on gap, 12–30% on multi_knapsack and 13–50% on set_cover,
/// where a Local-MIP improver closes 0–2.5%.
pub fn restart_from_best(
    p: &MipProblem,
    lower: &[f64],
    upper: &[f64],
    mut best: FjResult,
    cfg: FjConfig,
    restart_effort: f64,
) -> FjResult {
    let mut spent = 0;
    let mut left = (restart_effort.max(0.0) * effort_cap(p) as f64) as usize;
    let mut seed = cfg.seed;
    while left > 0 && !past(cfg.deadline) {
        seed = seed.wrapping_add(1);
        let c = FjConfig {
            seed,
            max_effort: Some(left),
            ..cfg
        };
        let Some(r) = feasibility_jump(p, lower, upper, Some(&best.x), c) else {
            break;
        };
        spent += r.effort;
        left = left.saturating_sub(r.effort.max(1));
        if r.objective < best.objective {
            best = r;
        }
    }
    best.effort = spent;
    best
}

impl Solver {
    fn run(&mut self, max_total: usize, max_since_improvement: usize) -> Option<FjResult> {
        self.init();
        let n = self.x.len();
        let mut best: Option<FjResult> = None;
        let mut best_violated = usize::MAX;
        let mut effort_at_improvement = 0;
        let mut effort_at_poll = 0;
        loop {
            if self.effort - effort_at_improvement > max_since_improvement
                || self.effort > max_total
            {
                break;
            }
            if self.effort - effort_at_poll > DEADLINE_EFFORT {
                effort_at_poll = self.effort;
                if past(self.cfg.deadline) {
                    break;
                }
            }
            if self.violated.items.len() < best_violated {
                effort_at_improvement = self.effort;
                best_violated = self.violated.items.len();
            }
            if self.violated.items.is_empty() {
                let best_obj = best.as_ref().map_or(f64::INFINITY, |b| b.objective);
                if self.objective < best_obj {
                    let rel = (best_obj - self.objective) / self.objective.abs().max(1.0);
                    if rel > MIN_RELATIVE_OBJECTIVE_IMPROVEMENT {
                        effort_at_improvement = self.effort;
                        best = Some(FjResult {
                            x: self.x.clone(),
                            objective: self.objective,
                            effort: self.effort,
                        });
                    }
                }
            }
            if n == 0 {
                break;
            }
            let var = self.select_variable();
            self.do_move(var);
        }
        if let Some(b) = best.as_mut() {
            b.effort = self.effort;
        }
        best
    }

    fn init(&mut self) {
        self.objective = self.obj.iter().zip(&self.x).map(|(c, v)| c * v).sum();
        for c in 0..self.cons.len() {
            let lhs = self.row_cells[self.row_ptr[c]..self.row_ptr[c + 1]]
                .iter()
                .map(|&(j, a)| a * self.x[j])
                .sum();
            self.cons[c].lhs = lhs;
            if self.cons[c].score(lhs) < -self.cfg.violation_tol {
                self.violated.insert(c);
            }
        }
        self.effort += self.row_cells.len();
        for j in 0..self.x.len() {
            self.reset_move(j);
        }
    }

    fn cells(&self, j: usize) -> &[(usize, f64)] {
        &self.col_cells[self.col_ptr[j]..self.col_ptr[j + 1]]
    }

    fn select_variable(&mut self) -> usize {
        if !self.good.items.is_empty() {
            let len = self.good.items.len();
            if self.rng.next_f64() < RANDOM_VAR_PROBABILITY {
                let k = self.rand_below(len);
                return self.good.items[k];
            }
            let sample = MAX_MOVES_TO_EVALUATE.min(len);
            self.effort += sample;
            let mut best = (f64::NEG_INFINITY, NONE);
            for _ in 0..sample {
                let k = self.rand_below(len);
                let j = self.good.items[k];
                if self.moves[j].score > best.0 {
                    best = (self.moves[j].score, j);
                }
            }
            return best.1;
        }
        // Local minimum.
        self.update_weights();
        if !self.violated.items.is_empty() {
            let k = self.rand_below(self.violated.items.len());
            let c = self.violated.items[k];
            let (from, to) = (self.row_ptr[c], self.row_ptr[c + 1]);
            if self.rng.next_f64() < RANDOM_CELL_PROBABILITY {
                let k = self.rand_below(to - from);
                return self.row_cells[from + k].0;
            }
            let mut best = (f64::NEG_INFINITY, self.row_cells[from].0);
            for &(j, _) in &self.row_cells[from..to] {
                if self.moves[j].score > best.0 {
                    best = (self.moves[j].score, j);
                }
            }
            return best.1;
        }
        self.rand_below(self.x.len())
    }

    fn rand_below(&mut self, n: usize) -> usize {
        (self.rng.next_u64() % n as u64) as usize
    }

    /// Bumps violated-row weights (or the objective weight when feasible)
    /// and updates the affected move scores in place.
    fn update_weights(&mut self) {
        let inc = self.weight_increment;
        let mut rescale = false;
        if self.violated.items.is_empty() {
            self.objective_weight += inc;
            rescale |= self.objective_weight > WEIGHT_RESCALE_THRESHOLD;
            self.effort += self.x.len();
            for j in 0..self.x.len() {
                // Minimization: a cheaper move scores positive.
                self.moves[j].score -= inc * self.obj[j] * (self.moves[j].value - self.x[j]);
            }
        } else {
            let violated = self.violated.items.clone();
            for c in violated {
                self.cons[c].weight += inc;
                rescale |= self.cons[c].weight > WEIGHT_RESCALE_THRESHOLD;
                let (from, to) = (self.row_ptr[c], self.row_ptr[c + 1]);
                self.effort += to - from;
                for k in from..to {
                    let (j, a) = self.row_cells[k];
                    let con = &self.cons[c];
                    let cand = con.lhs + a * (self.moves[j].value - self.x[j]);
                    self.moves[j].score += inc * (con.score(cand) - con.score(con.lhs));
                    self.update_good(j);
                }
            }
        }
        self.weight_increment /= self.cfg.weight_update_decay;
        if rescale {
            self.weight_increment *= 1e-20;
            self.objective_weight *= 1e-20;
            for c in &mut self.cons {
                c.weight *= 1e-20;
            }
            self.effort += self.cons.len();
            for j in 0..self.x.len() {
                self.reset_move(j);
            }
        }
    }

    fn do_move(&mut self, j: usize) {
        let new_value = self.moves[j].value;
        let delta = new_value - self.x[j];
        self.x[j] = new_value;
        self.objective += self.obj[j] * delta;
        let tol = self.cfg.violation_tol;
        for k in self.col_ptr[j]..self.col_ptr[j + 1] {
            let (c, a) = self.col_cells[k];
            let old_lhs = self.cons[c].lhs;
            let new_lhs = old_lhs + a * delta;
            self.cons[c].lhs = new_lhs;
            if self.cons[c].score(new_lhs) < -tol {
                self.violated.insert(c);
            } else {
                self.violated.remove(c);
            }
            let (from, to) = (self.row_ptr[c], self.row_ptr[c + 1]);
            self.effort += to - from;
            let con = &self.cons[c];
            for kk in from..to {
                let (v, coeff) = self.row_cells[kk];
                if v == j {
                    continue;
                }
                // Re-score v's move against the shifted LHS (modifyMove).
                let shift = coeff * (self.moves[v].value - self.x[v]);
                let old_term = con.score(old_lhs + shift) - con.score(old_lhs);
                let new_term = con.score(new_lhs + shift) - con.score(new_lhs);
                self.moves[v].score += con.weight * (new_term - old_term);
                let good = self.moves[v].score > 0.0;
                if good {
                    self.good.insert(v);
                } else {
                    self.good.remove(v);
                }
            }
        }
        self.reset_move(j);
    }

    fn update_good(&mut self, j: usize) {
        if self.moves[j].score > 0.0 {
            self.good.insert(j);
        } else {
            self.good.remove(j);
        }
    }

    /// Recomputes `j`'s jump value and its score from scratch.
    fn reset_move(&mut self, j: usize) {
        self.effort += self.col_ptr[j + 1] - self.col_ptr[j];
        let value = self.jump_value(j);
        let xj = self.x[j];
        let mut score = -self.objective_weight * self.obj[j] * (value - xj);
        for &(c, a) in self.cells(j) {
            let con = &self.cons[c];
            score += con.weight * (con.score(con.lhs + a * (value - xj)) - con.score(con.lhs));
        }
        self.moves[j] = Move { value, score };
        self.update_good(j);
    }

    /// Minimizer of `Σ_c w_c · dist(x_j, valid range of c)` over `[lb, ub]`,
    /// preferring values different from the current one; ties keep the
    /// smallest breakpoint.
    fn jump_value(&mut self, j: usize) -> f64 {
        let eq_tol = self.cfg.equality_tol;
        let xj = self.x[j];
        let (lower, upper) = (self.lb[j], self.ub[j]);
        let mut slope = 0.0;
        let mut bps = std::mem::take(&mut self.breakpoints);
        bps.clear();
        for &(c, a) in &self.col_cells[self.col_ptr[j]..self.col_ptr[j + 1]] {
            let con = &self.cons[c];
            let (lo, hi) = match con.sense {
                Sense::Le => (f64::NEG_INFINITY, con.rhs),
                Sense::Ge => (con.rhs, f64::INFINITY),
            };
            let residual = con.lhs - a * xj;
            let (mut first, mut second) = ((lo - residual) / a, (hi - residual) / a);
            if a < 0.0 {
                std::mem::swap(&mut first, &mut second);
            }
            if self.integer[j] {
                first = (first - eq_tol).ceil();
                second = (second + eq_tol).floor();
            }
            if first > second {
                continue;
            }
            if first > lower {
                slope -= con.weight;
                if first < upper {
                    bps.push((first, con.weight));
                }
            }
            if second <= lower {
                slope += con.weight;
            } else if second < upper {
                bps.push((second, con.weight));
            }
        }
        if lower.is_finite() {
            bps.push((lower, 0.0));
        }
        if upper.is_finite() {
            bps.push((upper, 0.0));
        }
        bps.sort_by(|a, b| a.0.total_cmp(&b.0));
        let eq = |a: f64, b: f64| (a - b).abs() < eq_tol;
        let mut value = bps.first().map_or(xj, |b| b.0);
        let (mut score, mut best_score, mut best_value) = (0.0, 0.0, value);
        for &(bp, dslope) in &bps {
            score += (bp - value) * slope;
            slope += dslope;
            value = bp;
            if eq(best_value, xj) || (!eq(value, xj) && score < best_score) {
                best_score = score;
                best_value = value;
            }
            // Convex: once the slope is non-negative the minimum is behind us.
            if !eq(best_value, xj) && slope >= 0.0 {
                break;
            }
        }
        self.breakpoints = bps;
        best_value
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

    fn run(p: &MipProblem) -> Option<FjResult> {
        feasibility_jump(p, &p.col_lower, &p.col_upper, None, FjConfig::new(1e-6))
    }

    fn assert_feasible(p: &MipProblem, x: &[f64]) {
        let v = p.violation(x);
        assert!(
            v.bound <= 1e-9 && v.row <= 1e-6 && v.integrality <= 1e-9,
            "{v:?} at {x:?}"
        );
    }

    #[test]
    fn set_cover_rows_with_positive_coefficients() {
        // Cover 4 elements with 4 sets; start at all-zero (every row violated).
        let inf = f64::INFINITY;
        let p = mip(
            &[
                (&[1.0, 1.0, 0.0, 0.0], 1.0, inf),
                (&[0.0, 1.0, 1.0, 0.0], 1.0, inf),
                (&[0.0, 0.0, 1.0, 1.0], 1.0, inf),
                (&[1.0, 0.0, 0.0, 1.0], 1.0, inf),
            ],
            &[0.0; 4],
            &[1.0; 4],
            &[true; 4],
            &[1.0; 4],
        );
        let r = run(&p).expect("cover is easy");
        assert_feasible(&p, &r.x);
        assert!(r.objective <= 2.0 + 1e-9, "{}", r.objective);
    }

    #[test]
    fn negative_coefficients_drive_the_jump() {
        // −x ≤ −3 (x ≥ 3) and −x − y ≤ −7 over integers in [0, 10]: only
        // negative coefficients, so a jump value ignoring them never moves.
        let inf = f64::INFINITY;
        let p = mip(
            &[(&[-1.0, 0.0], -inf, -3.0), (&[-1.0, -1.0], -inf, -7.0)],
            &[0.0, 0.0],
            &[10.0, 10.0],
            &[true, true],
            &[1.0, 2.0],
        );
        let r = run(&p).expect("feasible");
        assert_feasible(&p, &r.x);
    }

    #[test]
    fn equality_with_continuous_and_free_integer() {
        // 2x + y = 6.5, x free integer, y ∈ [0, 1] continuous; x − z ≤ 1.
        let inf = f64::INFINITY;
        let p = mip(
            &[(&[2.0, 1.0, 0.0], 6.5, 6.5), (&[1.0, 0.0, -1.0], -inf, 1.0)],
            &[-inf, 0.0, 0.0],
            &[inf, 1.0, 5.0],
            &[true, false, true],
            &[0.0, 0.0, 1.0],
        );
        let r = run(&p).expect("feasible");
        assert_feasible(&p, &r.x);
    }

    #[test]
    fn objective_weight_pushes_towards_cheaper_points() {
        // x + y ≥ 1 over integers in [0, 5], cost (3, 1): the first feasible
        // point is x = 1 or y = 1 depending on the seed; the objective weight
        // must then pull the search to y = 1 (cost 1), never reward the
        // costlier side.
        let p = mip(
            &[(&[1.0, 1.0], 1.0, f64::INFINITY)],
            &[0.0, 0.0],
            &[5.0, 5.0],
            &[true, true],
            &[3.0, 1.0],
        );
        let mut cheapest = 0;
        for seed in 0..50 {
            let cfg = FjConfig {
                seed,
                ..FjConfig::new(1e-6)
            };
            let r = feasibility_jump(&p, &p.col_lower, &p.col_upper, None, cfg).unwrap();
            assert_feasible(&p, &r.x);
            cheapest += usize::from(r.objective == 1.0);
        }
        assert!(cheapest >= 45, "{cheapest}/50");
    }

    /// Past the deadline no restart starts; the first point stands.
    #[test]
    fn restarts_stop_at_the_deadline() {
        let w = [3.0, 5.0, 2.0, 7.0, 4.0, 6.0];
        let p = mip(
            &[(&w, f64::NEG_INFINITY, 11.0)],
            &[0.0; 6],
            &[1.0; 6],
            &[true; 6],
            &[-4.0, -6.0, -3.0, -8.0, -5.0, -7.0],
        );
        let first = run(&p).unwrap();
        let cfg = FjConfig {
            deadline: Some(std::time::Instant::now()),
            ..FjConfig::new(1e-6)
        };
        let b = restart_from_best(&p, &p.col_lower, &p.col_upper, first.clone(), cfg, 3.0);
        assert_eq!(b.effort, 0);
        assert_eq!(b.x, first.x);
        let free = restart_from_best(
            &p,
            &p.col_lower,
            &p.col_upper,
            first,
            FjConfig::new(1e-6),
            3.0,
        );
        assert!(free.effort > 0);
    }

    #[test]
    fn infeasible_parity_finds_nothing() {
        // 2x − 2y = 1 over integers.
        let p = mip(
            &[(&[2.0, -2.0], 1.0, 1.0)],
            &[0.0, 0.0],
            &[10.0, 10.0],
            &[true, true],
            &[1.0, 1.0],
        );
        assert!(run(&p).is_none());
    }

    #[test]
    fn empty_integer_domain_is_rejected() {
        let p = mip(&[], &[0.2], &[0.8], &[true], &[1.0]);
        assert!(run(&p).is_none());
    }

    #[test]
    fn random_knapsacks_feasible_points_verify() {
        // Every returned point must pass the independent violation check.
        let mut rng = XorShift64::new(11);
        let mut found = 0;
        for _ in 0..50 {
            let n = 8;
            let w: Vec<f64> = (0..n).map(|_| 1.0 + (rng.next_u64() % 9) as f64).collect();
            let v: Vec<f64> = (0..n).map(|_| -1.0 - (rng.next_u64() % 9) as f64).collect();
            let cap = w.iter().sum::<f64>() * 0.4;
            let p = mip(
                &[(&w, f64::NEG_INFINITY, cap), (&v, f64::NEG_INFINITY, -6.0)],
                &[0.0; 8],
                &[1.0; 8],
                &[true; 8],
                &v,
            );
            if let Some(r) = run(&p) {
                found += 1;
                assert_feasible(&p, &r.x);
                assert!((p.objective(&r.x) - r.objective).abs() < 1e-9);
            }
        }
        assert!(found >= 40, "{found}");
    }

    #[test]
    fn restarts_never_lose_the_first_point_and_stay_feasible() {
        let mut rng = XorShift64::new(12);
        let mut improved = 0;
        for _ in 0..50 {
            let n = 12;
            let w: Vec<f64> = (0..n).map(|_| 1.0 + (rng.next_u64() % 9) as f64).collect();
            let v: Vec<f64> = (0..n).map(|_| -1.0 - (rng.next_u64() % 9) as f64).collect();
            let cap = w.iter().sum::<f64>() * 0.4;
            let p = mip(
                &[(&w, f64::NEG_INFINITY, cap)],
                &[0.0; 12],
                &[1.0; 12],
                &[true; 12],
                &v,
            );
            let cfg = FjConfig::new(1e-6);
            let first = run(&p).unwrap();
            let b = restart_from_best(&p, &p.col_lower, &p.col_upper, first.clone(), cfg, 3.0);
            assert_feasible(&p, &b.x);
            assert!((p.objective(&b.x) - b.objective).abs() < 1e-9);
            assert!(b.objective <= first.objective);
            assert!(b.effort > 0 && b.effort <= 3 * effort_cap(&p) + effort_cap(&p));
            improved += usize::from(b.objective < first.objective);
        }
        // Restarts do find better points (6/50 with this seed).
        assert!(improved > 0);
    }
}
