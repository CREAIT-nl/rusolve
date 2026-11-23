//! Primal heuristics: cheap ways to find an
//! incumbent before or alongside the tree search, so the cutoff bites early.
//!
//! Every heuristic here only *proposes* points; the tree search accepts one
//! only after the same bound / row / integrality check every LP incumbent
//! passes (`bnb::Search::try_incumbent`).
//!
//! * [`trivial`] — all-zero / all-lower / all-upper / lock points for pure
//!   integer models (HiGHS `HighsMipSolverData::trivialHeuristics`,
//!   reference/HiGHS/highs/mip/HighsMipSolverData.cpp:265–410, MIT).
//! * [`fj`] — Feasibility Jump before the root LP (Luteberget & Sartor 2023).
//! * [`localmip`] — Local-MIP improver of a feasible point (Lin, Zou & Cai
//!   2024).
//! * [`rounding`] — LP-solution rounding: simple, randomized, ZI round,
//!   shifting.
//! * [`diving`] — fractional and coefficient LP diving with propagation.
//!
//! Motivation: on our generated families HiGHS has an incumbent before the
//! root LP from Feasibility Jump (8/16 families) or the trivial heuristics
//! (6/16).

pub mod diving;
pub mod fj;
pub mod localmip;
pub mod rounding;

use crate::problem::MipProblem;

/// Where an incumbent came from (statistics / ablation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    /// An integral node-LP solution.
    Branching,
    Trivial,
    FeasibilityJump,
    SimpleRounding,
    RandomizedRounding,
    ZiRound,
    Shifting,
    Diving,
    /// Local-MIP polish of an existing incumbent ([`localmip`]).
    LocalMip,
    /// RINS / RENS sub-MIP ([`crate::submip`]).
    SubMip,
    /// Found by another racer ([`crate::race`]).
    Racer,
}

/// Which heuristics run. Defaults follow measurement (our
/// medium benchmark set, same build with heuristics off vs on, 2026-10-22):
///
/// * ZI round and shifting are off, as in HiGHS (`mip_heuristic_run_zi_round`
///   and `mip_heuristic_run_shifting` default false, HighsOptions.h). They
///   never changed the incumbent, node count or time on any medium family.
/// * LP diving is off for now. On the cold-start engine every dive step is
///   a full LP solve: about 2x the LP iterations on set_cover, +15% time on
///   cap_facility_location, and node/time gains only on lot_sizing.
///   TODO(perf): revisit once node LPs warm-start.
/// * FJ restarts from its best point with one extra FJ budget
///   (`fj_restart_effort = 1`). Overall time ratio 0.89 against a single
///   FJ run. bin_packing s1 drops from 139 nodes to 1. No family moves
///   outside the 0.89–1.06 noise band in the slow direction.
#[derive(Debug, Clone, Copy)]
pub struct HeuristicConfig {
    pub trivial: bool,
    pub feasibility_jump: bool,
    pub simple_rounding: bool,
    pub randomized_rounding: bool,
    pub zi_round: bool,
    pub shifting: bool,
    pub diving: bool,
    /// Extra Feasibility Jump effort, in multiples of its own cap
    /// (`fj::effort_cap`), spent restarting it from its best point before
    /// the root LP (`fj::restart_from_best`); 0 = single run as in HiGHS.
    pub fj_restart_effort: f64,
    /// Share of LP iterations heuristics may use once the initial offset of
    /// 10 000 is spent (HiGHS `mip_heuristic_effort` default 0.05,
    /// HighsOptions.h:1224; budget rule HighsMipSolverData.cpp:683–687).
    pub effort: f64,
}

impl Default for HeuristicConfig {
    fn default() -> Self {
        Self {
            trivial: true,
            feasibility_jump: true,
            simple_rounding: true,
            randomized_rounding: true,
            zi_round: false,
            shifting: false,
            diving: false,
            fj_restart_effort: 1.0,
            effort: 0.05,
        }
    }
}

impl HeuristicConfig {
    /// Every heuristic on.
    pub fn all() -> Self {
        Self {
            trivial: true,
            feasibility_jump: true,
            simple_rounding: true,
            randomized_rounding: true,
            zi_round: true,
            shifting: true,
            diving: true,
            fj_restart_effort: 1.0,
            effort: 0.05,
        }
    }

    pub fn none() -> Self {
        Self {
            trivial: false,
            feasibility_jump: false,
            simple_rounding: false,
            randomized_rounding: false,
            zi_round: false,
            shifting: false,
            diving: false,
            fj_restart_effort: 0.0,
            effort: 0.0,
        }
    }

    /// Switches off the heuristics named in a comma-separated list
    /// (`trivial,fj,simple,randomized,zi,shifting,diving`, or `all`).
    pub fn disable(&mut self, names: &str) -> Result<(), String> {
        self.set(names, false)
    }

    /// Switches on the heuristics named like in [`Self::disable`].
    pub fn enable(&mut self, names: &str) -> Result<(), String> {
        self.set(names, true)
    }

    fn set(&mut self, names: &str, on: bool) -> Result<(), String> {
        for name in names.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            match name {
                "all" if on => *self = Self::all(),
                "all" => *self = Self::none(),
                "trivial" => self.trivial = on,
                "fj" => self.feasibility_jump = on,
                "simple" => self.simple_rounding = on,
                "randomized" => self.randomized_rounding = on,
                "zi" => self.zi_round = on,
                "shifting" => self.shifting = on,
                "diving" => self.diving = on,
                other => return Err(format!("unknown heuristic '{other}'")),
            }
        }
        if on && self.effort == 0.0 {
            self.effort = Self::all().effort;
        }
        Ok(())
    }
}

/// Rounding locks: `up[j]` counts the rows that may become violated when
/// `x_j` increases (`a > 0` with a finite upper side, `a < 0` with a finite
/// lower side), `down[j]` the converse. A column without down-locks can
/// always be rounded down without hurting any row. HiGHS
/// HighsMipSolverData.cpp:922–943 (MIT); Achterberg 2007 §3.3.3.
#[derive(Debug, Clone)]
pub struct Locks {
    pub up: Vec<u32>,
    pub down: Vec<u32>,
}

impl Locks {
    pub fn new(p: &MipProblem) -> Self {
        let mut up = vec![0; p.num_col];
        let mut down = vec![0; p.num_col];
        let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
        for j in 0..p.num_col {
            for k in ptr[j]..ptr[j + 1] {
                let (i, a) = (idx[k] as usize, val[k]);
                if a == 0.0 {
                    continue;
                }
                if p.row_lower[i].is_finite() {
                    if a < 0.0 {
                        up[j] += 1;
                    } else {
                        down[j] += 1;
                    }
                }
                if p.row_upper[i].is_finite() {
                    if a < 0.0 {
                        down[j] += 1;
                    } else {
                        up[j] += 1;
                    }
                }
            }
        }
        Self { up, down }
    }
}

/// Candidate points of the trivial heuristics, for pure integer models only
/// (HiGHS skips them with continuous columns, whose values would need an LP):
/// all-zero, all-lower, all-upper, and the *lock point* — each column at the
/// bound its coefficient signs favour (lower if it has more positive than
/// negative coefficients). Integer bounds are rounded inward first; points
/// with infinite entries are skipped. The caller checks feasibility.
pub fn trivial(p: &MipProblem, lower: &[f64], upper: &[f64]) -> Vec<Vec<f64>> {
    if p.integer.iter().any(|&int| !int) {
        return vec![];
    }
    let lo: Vec<f64> = lower.iter().map(|l| l.ceil()).collect();
    let hi: Vec<f64> = upper.iter().map(|u| u.floor()).collect();
    if lo.iter().zip(&hi).any(|(l, h)| l > h) {
        return vec![];
    }
    let mut out = Vec::new();
    let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
    if lo.iter().all(|&l| l <= 0.0) && hi.iter().all(|&h| h >= 0.0) {
        out.push(vec![0.0; p.num_col]);
    }
    if finite(&lo) && lo.iter().any(|&l| l != 0.0) {
        out.push(lo.clone());
    }
    if finite(&hi) {
        out.push(hi.clone());
    }
    if finite(&lo) && finite(&hi) {
        let (ptr, val) = (p.matrix.col_ptr(), p.matrix.values());
        let lock_point = (0..p.num_col)
            .map(|j| {
                let pos = val[ptr[j]..ptr[j + 1]].iter().filter(|&&a| a > 0.0).count();
                let neg = val[ptr[j]..ptr[j + 1]].iter().filter(|&&a| a < 0.0).count();
                if pos > neg { lo[j] } else { hi[j] }
            })
            .collect();
        out.push(lock_point);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_core::CscMatrix;

    #[test]
    fn locks_follow_row_sides_and_signs() {
        // r0: x0 − x1 ≤ 4 ; r1: 1 ≤ x0 + x1 ≤ 3 ; r2: x1 ≥ 0 (free above).
        let inf = f64::INFINITY;
        let trip = [
            (0, 0, 1.0),
            (0, 1, -1.0),
            (1, 0, 1.0),
            (1, 1, 1.0),
            (2, 1, 1.0),
        ];
        let p = MipProblem {
            num_col: 2,
            num_row: 3,
            matrix: CscMatrix::from_triplets(3, 2, &trip),
            row_lower: vec![-inf, 1.0, 0.0],
            row_upper: vec![4.0, 3.0, inf],
            cost: vec![0.0; 2],
            maximize: false,
            col_lower: vec![0.0; 2],
            col_upper: vec![5.0; 2],
            integer: vec![true; 2],
        };
        let l = Locks::new(&p);
        assert_eq!(l.up, vec![2, 1]);
        assert_eq!(l.down, vec![1, 3]);
    }

    #[test]
    fn trivial_points_skip_infinite_and_mixed() {
        let inf = f64::INFINITY;
        let mut p = MipProblem {
            num_col: 2,
            num_row: 1,
            matrix: CscMatrix::from_triplets(1, 2, &[(0, 0, 1.0), (0, 1, -1.0)]),
            row_lower: vec![-inf],
            row_upper: vec![1.0],
            cost: vec![0.0; 2],
            maximize: false,
            col_lower: vec![-1.0, 0.0],
            col_upper: vec![2.0, inf],
            integer: vec![true; 2],
        };
        let pts = trivial(&p, &p.col_lower, &p.col_upper);
        // zero, lower; upper/lock skipped (x1 unbounded above).
        assert_eq!(pts, vec![vec![0.0, 0.0], vec![-1.0, 0.0]]);
        p.integer[1] = false;
        assert!(trivial(&p, &p.col_lower, &p.col_upper).is_empty());
    }
}
