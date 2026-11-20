//! Native-engine Feasibility Jump on the port's presolved model (lever `mip_fj_native`, default off).
//!
//! On decomp2 the port's HiGHS 1.16 FJ stalls at -135..-144 (best -160) and the root sub-MIP needs
//! 7-19 s, while this search, run on the same presolved model, reaches -160 in 0.1-0.4 s with 0 LP
//! iterations (measured on our benchmark set). It is the FJ of `rusolve-mip`
//! (`crates/rusolve-mip/src/heuristics/fj.rs`), copied here
//! because `rusolve-mip` depends on this crate.
//!
//! Provenance: B. Luteberget & G. Sartor, "Feasibility Jump: an LP-free Lagrangian MIP heuristic",
//! Math. Prog. Comp. 15:365-388 (2023); HiGHS `mip/feasibilityjump.hh` (MIT) as a reference text,
//! reimplemented. Deviations from the HiGHS copy are the ones listed in the native module: jump-value
//! range endpoints swapped for negative coefficients, objective term `-w*c*delta`, effort checked
//! every step.

use crate::lp_data::hconst::{HighsVarType, kHighsInf};
use crate::lp_data::highs_lp::HighsLp;
use crate::util::sparse_matrix::HighsSparseMatrix;
use std::time::Instant;

const RANDOM_VAR_PROBABILITY: f64 = 0.001;
const RANDOM_CELL_PROBABILITY: f64 = 0.01;
const MAX_MOVES_TO_EVALUATE: usize = 25;
const MIN_RELATIVE_OBJECTIVE_IMPROVEMENT: f64 = 1e-4;
const WEIGHT_RESCALE_THRESHOLD: f64 = 1e20;
const DEADLINE_EFFORT: usize = 500_000;
const NONE: usize = usize::MAX;

struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                seed
            },
        }
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

struct FjConfig {
    equality_tol: f64,
    violation_tol: f64,
    weight_update_decay: f64,
    deadline: Option<Instant>,
}

fn past(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|d| Instant::now() >= d)
}

/// The best feasible point found, as column values.
pub struct FjNativeResult {
    pub x: Vec<f64>,
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

impl Solver {
    fn run(&mut self, max_total: usize, max_since_improvement: usize) -> Option<Vec<f64>> {
        self.init();
        let n = self.x.len();
        let mut best: Option<Vec<f64>> = None;
        let mut best_objective = f64::INFINITY;
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
                let best_obj = best_objective;
                if self.objective < best_obj {
                    let rel = (best_obj - self.objective) / self.objective.abs().max(1.0);
                    if rel > MIN_RELATIVE_OBJECTIVE_IMPROVEMENT {
                        effort_at_improvement = self.effort;
                        best_objective = self.objective;
                        best = Some(self.x.clone());
                    }
                }
            }
            if n == 0 {
                break;
            }
            let var = self.select_variable();
            self.do_move(var);
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

/// Runs the search on `model` (minimising `sense * cost`), starting from `start` clamped to the
/// bounds, and returns the best feasible point within the `nnz * 2^10` effort budget.
pub fn feasibility_jump_native(
    model: &HighsLp,
    start: Option<&[f64]>,
    feastol: f64,
    epsilon: f64,
    seed: u64,
    deadline: Option<Instant>,
) -> Option<FjNativeResult> {
    let n = model.num_col_ as usize;
    let sense = model.sense_ as i32 as f64;
    let cfg = FjConfig {
        equality_tol: epsilon,
        violation_tol: feastol,
        weight_update_decay: 1.0,
        deadline,
    };
    let mut lb = vec![0.0; n];
    let mut ub = vec![0.0; n];
    let mut integer = vec![false; n];
    for j in 0..n {
        let (mut l, mut u) = (model.col_lower_[j], model.col_upper_[j]);
        if model.integrality_[j] != HighsVarType::kContinuous {
            integer[j] = true;
            l = (l - feastol).ceil();
            u = (u + feastol).floor();
        }
        if !(l <= u && l < kHighsInf && u > -kHighsInf) || l.is_nan() || u.is_nan() {
            return None;
        }
        lb[j] = l;
        ub[j] = u;
    }
    let x0: Vec<f64> = (0..n)
        .map(|j| match start {
            Some(s) if s[j].is_finite() => s[j].clamp(lb[j], ub[j]),
            _ if lb[j].is_finite() => lb[j],
            _ if ub[j].is_finite() => ub[j],
            _ => 0.0,
        })
        .collect();
    let mut a = HighsSparseMatrix::default();
    a.create_rowwise(&model.a_matrix_);
    let mut cons = Vec::new();
    let mut row_ptr = vec![0usize];
    let mut row_cells: Vec<(usize, f64)> = Vec::new();
    for i in 0..model.num_row_ as usize {
        let (s, e) = (a.start_[i] as usize, a.start_[i + 1] as usize);
        let cells: Vec<(usize, f64)> = (s..e)
            .filter(|&k| a.value_[k] != 0.0)
            .map(|k| (a.index_[k] as usize, a.value_[k]))
            .collect();
        if cells.is_empty() {
            continue;
        }
        for (sn, rhs) in [
            (Sense::Ge, model.row_lower_[i]),
            (Sense::Le, model.row_upper_[i]),
        ] {
            if rhs.is_finite() {
                cons.push(Constraint {
                    sense: sn,
                    rhs,
                    weight: 1.0,
                    lhs: 0.0,
                });
                row_cells.extend_from_slice(&cells);
                row_ptr.push(row_cells.len());
            }
        }
    }
    let nnz = row_cells.len().max(1);
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
        for &(j, v) in &row_cells[row_ptr[c]..row_ptr[c + 1]] {
            col_cells[fill[j]] = (c, v);
            fill[j] += 1;
        }
    }
    let num_cons = cons.len();
    let mut s = Solver {
        cfg,
        lb,
        ub,
        integer,
        obj: model.col_cost_.iter().map(|c| sense * c).collect(),
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
        rng: XorShift64::new(seed),
        effort: 0,
        breakpoints: Vec::new(),
    };
    s.run(nnz << 10, nnz << 8).map(|x| FjNativeResult { x })
}
