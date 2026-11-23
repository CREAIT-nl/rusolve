//! LP-solution rounding heuristics: simple rounding, randomized rounding,
//! ZI round, shifting, and the shared fix-propagate-resolve step
//! (`try_rounded_point`).
//!
//! Provenance:
//!
//! * Simple rounding — T. Achterberg, *Constraint Integer Programming*, PhD
//!   thesis TU Berlin (2007), §9.1 (rounding heuristics): a fractional
//!   column without locks in one direction can be rounded that way without
//!   violating any row. (HiGHS has no separate simple rounding.)
//! * Randomized rounding — HiGHS `HighsPrimalHeuristics::randomizedRounding`
//!   (reference/HiGHS/highs/mip/HighsPrimalHeuristics.cpp:1117, MIT): round
//!   lock-free columns towards their free side, others at a random threshold
//!   in [0.1, 0.9]; fix + propagate column by column; LP over the continuous
//!   columns when there are any.
//! * ZI round — C. Wallace, "ZI round, a MIP rounding heuristic", J.
//!   Heuristics 16:715–722 (2010); HiGHS `ziRound` (:1465): move each
//!   fractional column within its row slacks to the side that reduces its
//!   fractionality most; ≤ 6 passes.
//! * Shifting — Achterberg 2007 §9.1 (SCIP "shifting"); HiGHS `shifting`
//!   (:1204): repair violated rows by shifting columns (fractional first,
//!   lock-scored), with the 1.1^(shift−t) anti-cycling score.
//! * `try_rounded_point` — HiGHS `tryRoundedPoint` (:971).
//!
//! All of it is reimplemented from those descriptions, not translated.
//! Candidates are only proposed (`Ctx::found`); the tree search verifies
//! each before accepting it.

use super::{Locks, Source};
use crate::engine::{LpEngine, LpStatus};
use crate::problem::MipProblem;
use crate::propagate::Propagator;
use rusolve_lp::rng::XorShift64;

/// ZI round's pass limit (HiGHS `max_loop_count = 5`, loop runs while
/// `count <= 5`, i.e. at most 6 passes).
const ZI_MAX_PASSES: usize = 6;
/// Shifting gives up after this many passes without progress (HiGHS
/// `maxIterationsWithoutReductions = 5`, `<=` test).
const SHIFT_MAX_STALL: usize = 5;
/// Base of the anti-cycling penalty for recently shifted columns.
const SHIFT_RECENCY_BASE: f64 = 1.1;

/// What a rounding heuristic needs from the search.
pub struct Ctx<'a, E: LpEngine> {
    pub p: &'a MipProblem,
    pub locks: &'a Locks,
    pub prop: &'a mut Propagator,
    pub engine: &'a mut E,
    /// The global domain fixings start from.
    pub lower: &'a [f64],
    pub upper: &'a [f64],
    /// Objective cutoff for propagation (the search's upper limit).
    pub cutoff: f64,
    pub feastol: f64,
    pub rng: &'a mut XorShift64,
    /// LP iterations spent by heuristic LPs.
    pub lp_iterations: u64,
    /// Proposed points, to be verified by the caller.
    pub found: Vec<(Vec<f64>, Source)>,
    /// Wall-clock end of the search: a heuristic that sees it passed gives
    /// up without proposing a point:
    /// one randomized-rounding pass over nw04's 87k columns (a propagation
    /// per fixing over 17k-entry rows) took 80 s, past a 60 s limit. Polled
    /// every [`DEADLINE_STRIDE`] steps; HiGHS polls its time limit inside
    /// heuristics the same coarse way (e.g. FJ, feasibilityjump.hh:774).
    pub deadline: Option<std::time::Instant>,
}

/// Heuristic steps (fixings, shift moves) between two deadline polls.
const DEADLINE_STRIDE: usize = 64;

fn fractionality(v: f64) -> f64 {
    (v - v.round()).abs()
}

fn has_continuous(p: &MipProblem) -> bool {
    p.integer.iter().any(|&int| !int)
}

impl<E: LpEngine> Ctx<'_, E> {
    /// `true` once the search's deadline has passed (reads the clock).
    pub fn past_deadline(&self) -> bool {
        self.deadline
            .is_some_and(|d| std::time::Instant::now() >= d)
    }

    /// Fractional column with a lock-free side → round there; any column
    /// locked both ways ⇒ fail. Continuous columns keep their LP values, so
    /// from an LP-feasible `x` the result is feasible by construction.
    pub fn simple_rounding(&mut self, x: &[f64]) {
        let mut y = x.to_vec();
        for j in 0..self.p.num_col {
            if !self.p.integer[j] || fractionality(x[j]) <= self.feastol {
                continue;
            }
            y[j] = if self.locks.down[j] == 0 {
                x[j].floor()
            } else if self.locks.up[j] == 0 {
                x[j].ceil()
            } else {
                return;
            };
        }
        self.found.push((y, Source::SimpleRounding));
    }

    /// Solves the LP over the continuous columns with integers fixed in
    /// `lower`/`upper`; `None` unless Optimal.
    fn resolve(&mut self, lower: &[f64], upper: &[f64]) -> Option<Vec<f64>> {
        let out = self.engine.solve(lower, upper).ok()?;
        self.lp_iterations += out.iterations;
        (out.status == LpStatus::Optimal).then_some(out.x)
    }

    /// Fixes integer `j` to `v` (clamped into the local box) and propagates.
    /// `false` iff the local domain became infeasible.
    fn fix_and_propagate(&mut self, j: usize, v: f64, lo: &mut [f64], hi: &mut [f64]) -> bool {
        let v = v.clamp(lo[j], hi[j]);
        lo[j] = v;
        hi[j] = v;
        self.prop.propagate(lo, hi, Some(j), self.cutoff)
    }

    pub fn randomized_rounding(&mut self, x: &[f64]) {
        let (mut lo, mut hi) = (self.lower.to_vec(), self.upper.to_vec());
        let mut fixed = 0usize;
        for (j, &xj) in x.iter().enumerate() {
            if !self.p.integer[j] {
                continue;
            }
            fixed += 1;
            if fixed.is_multiple_of(DEADLINE_STRIDE) && self.past_deadline() {
                return;
            }
            let v = if self.locks.up[j] == 0 {
                (xj - self.feastol).ceil()
            } else if self.locks.down[j] == 0 {
                (xj + self.feastol).floor()
            } else {
                (xj + self.rng.uniform(0.1, 0.9)).floor()
            };
            if !self.fix_and_propagate(j, v, &mut lo, &mut hi) {
                return;
            }
        }
        if has_continuous(self.p) {
            if let Some(y) = self.resolve(&lo, &hi) {
                self.found.push((y, Source::RandomizedRounding));
            }
        } else {
            self.found.push((lo, Source::RandomizedRounding));
        }
    }

    /// Fixes the integral integer entries of `point`, propagates, and
    /// completes the continuous part by LP; fractional entries stay free in
    /// that LP and its solution goes through ZI round.
    pub fn try_rounded_point(&mut self, point: &[f64], source: Source) {
        let (mut lo, mut hi) = (self.lower.to_vec(), self.upper.to_vec());
        let mut integral = true;
        for (j, &v) in point.iter().enumerate() {
            if !self.p.integer[j] {
                continue;
            }
            if fractionality(v) > self.feastol {
                integral = false;
                continue;
            }
            if !self.fix_and_propagate(j, v.round(), &mut lo, &mut hi) {
                return;
            }
        }
        if has_continuous(self.p) {
            let Some(y) = self.resolve(&lo, &hi) else {
                return;
            };
            if !integral {
                self.zi_round(&y);
            }
            self.found.push((y, source));
        } else {
            self.found.push((lo, source));
        }
    }

    /// ZI round: per fractional integer column, the largest moves up/down
    /// that keep every row within its bounds (and the column within its
    /// bounds); take the one that lowers the fractionality more (cheaper
    /// objective on ties). Proposes the result.
    pub fn zi_round(&mut self, x: &[f64]) {
        let p = self.p;
        let ft = self.feastol;
        let zi = |v: f64| ((v - ft).ceil() - v).min(v - (v + ft).floor());
        let mut y = x.to_vec();
        let mut total: f64 = (0..p.num_col)
            .filter(|&j| p.integer[j])
            .map(|j| zi(y[j]))
            .sum();
        if total <= ft {
            return;
        }
        let mut act = p.activities(&y);
        let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
        let mut passes = 0;
        let mut improvement = f64::INFINITY;
        while total > ft && improvement > ft && passes < ZI_MAX_PASSES {
            if self.past_deadline() {
                return;
            }
            let previous = total;
            passes += 1;
            for j in 0..p.num_col {
                let v = y[j];
                if !p.integer[j] || fractionality(v) <= ft {
                    continue;
                }
                let (mut up_room, mut down_room) = (f64::INFINITY, f64::INFINITY);
                for k in ptr[j]..ptr[j + 1] {
                    let (i, a) = (idx[k] as usize, val[k]);
                    if a == 0.0 {
                        continue;
                    }
                    let slack_up = p.row_upper[i] - act[i];
                    let slack_down = act[i] - p.row_lower[i];
                    if a > 0.0 {
                        up_room = up_room.min(slack_up / a);
                        down_room = down_room.min(slack_down / a);
                    } else {
                        up_room = up_room.min(-slack_down / a);
                        down_room = down_room.min(-slack_up / a);
                    }
                }
                let up = (p.col_upper[j] - v).min(up_room);
                let down = (v - p.col_lower[j]).min(down_room);
                let (z_up, z_down, z_now) = (zi(v + up), zi(v - down), zi(v));
                let step = if (z_up - z_down).abs() <= ft && z_up < z_now {
                    // Equal gain both ways: the cheaper objective side.
                    if p.cost[j] * (v + up) <= p.cost[j] * (v - down) {
                        up
                    } else {
                        -down
                    }
                } else if z_up < z_down && z_up < z_now {
                    up
                } else if z_up > z_down && z_down < z_now {
                    -down
                } else {
                    continue;
                };
                if !step.is_finite() || step == 0.0 {
                    continue;
                }
                y[j] = v + step;
                total += zi(y[j]) - z_now;
                for k in ptr[j]..ptr[j + 1] {
                    act[idx[k] as usize] += val[k] * step;
                }
            }
            improvement = previous - total;
        }
        self.found.push((y, Source::ZiRound));
    }

    /// Shifting: alternately repair a violated row (shift one column in the
    /// repairing direction, preferring fractional integers, then columns not
    /// shifted recently) or, with all rows satisfied, round the fractional
    /// integer with the most locks on its rounding side. Ends in
    /// `try_rounded_point` (rows still violated), ZI round (fractional
    /// integers left) or a direct proposal.
    pub fn shifting(&mut self, x: &[f64]) {
        let p = self.p;
        let ft = self.feastol;
        let num_int = p.integer.iter().filter(|&&b| b).count();
        let mut y = x.to_vec();
        let mut fractional: Vec<usize> = (0..p.num_col)
            .filter(|&j| p.integer[j] && fractionality(y[j]) > ft)
            .collect();
        let mut is_frac = vec![false; p.num_col];
        for &j in &fractional {
            is_frac[j] = true;
        }
        let csr = p.matrix.to_csr();
        let (rptr, ridx, rval) = (csr.row_ptr(), csr.col_idx(), csr.values());
        let infeasible_rows = |y: &[f64]| -> Vec<(usize, i8, f64)> {
            let act = p.activities(y);
            let mut out = Vec::new();
            for (i, &a) in act.iter().enumerate() {
                if a > p.row_upper[i] + ft {
                    out.push((i, 1, a - p.row_upper[i]));
                }
                if a < p.row_lower[i] - ft {
                    out.push((i, -1, p.row_lower[i] - a));
                }
            }
            out
        };
        // TODO(perf): activities are recomputed per pass (O(nnz) each, ≤
        // #integers passes); fine at mip_diff sizes, incremental later.
        let mut rows = infeasible_rows(&y);
        let mut prev_rows = rows.len();
        // Signed iteration numbers of past shifts per column (+ up, − down).
        let mut shifts: Vec<Vec<i64>> = vec![Vec::new(); p.num_col];
        let mut stall = 0;
        let mut t: i64 = 0;
        while (!fractional.is_empty() || !rows.is_empty())
            && stall <= SHIFT_MAX_STALL
            && t <= num_int as i64
        {
            t += 1;
            if (t as usize).is_multiple_of(DEADLINE_STRIDE) && self.past_deadline() {
                return;
            }
            stall += 1;
            let mut reduced = false;
            if !rows.is_empty() {
                // First violated row touching a fractional integer, else a
                // random violated row.
                let pick = rows
                    .iter()
                    .position(|&(r, _, _)| {
                        (rptr[r]..rptr[r + 1]).any(|k| is_frac[ridx[k] as usize])
                    })
                    .unwrap_or_else(|| (self.rng.next_u64() % rows.len() as u64) as usize);
                let (r, sense, infeas) = rows[pick];
                // (score, column, coefficient, direction)
                let mut best: Option<(f64, usize, f64, i8)> = None;
                for k in rptr[r]..rptr[r + 1] {
                    let (j, a) = (ridx[k] as usize, rval[k]);
                    if p.col_lower[j] == p.col_upper[j] || a == 0.0 {
                        continue;
                    }
                    for dir in [1i8, -1] {
                        // Must move the activity towards the violated side.
                        if f64::from(dir) * a * f64::from(sense) > 0.0 {
                            continue;
                        }
                        let at_bound = if dir > 0 {
                            (p.col_upper[j] - y[j]).abs() <= ft
                        } else {
                            (y[j] - p.col_lower[j]).abs() <= ft
                        };
                        if at_bound {
                            continue;
                        }
                        let score = if is_frac[j] {
                            // More locks on the opposite side = the move is
                            // safe for more rows (HiGHS passes downlocks for
                            // up-moves, :1349).
                            let locks = if dir > 0 {
                                self.locks.down[j]
                            } else {
                                self.locks.up[j]
                            };
                            -1.0 + 1.0 / f64::from(locks + 1)
                        } else {
                            let d = f64::from(dir);
                            let mut s = if shifts[j].is_empty() {
                                d * p.cost[j]
                            } else {
                                shifts[j]
                                    .iter()
                                    .filter(|&&s| d * s as f64 > 0.0)
                                    .map(|&s| SHIFT_RECENCY_BASE.powf(d * s as f64 - t as f64))
                                    .sum()
                            };
                            if p.integer[j] {
                                s += 1.0;
                            }
                            s
                        };
                        if best.is_none_or(|b| score < b.0) {
                            best = Some((score, j, a, dir));
                        }
                    }
                }
                if let Some((_, j, a, dir)) = best {
                    let xj = y[j];
                    if is_frac[j] {
                        is_frac[j] = false;
                        fractional.retain(|&c| c != j);
                        reduced = true;
                        y[j] = if dir > 0 {
                            (xj - ft).ceil()
                        } else {
                            (xj + ft).floor()
                        };
                    } else {
                        y[j] = if p.integer[j] {
                            xj + f64::from(dir)
                        } else if dir > 0 {
                            (xj + infeas / a.abs()).min(p.col_upper[j])
                        } else {
                            (xj - infeas / a.abs()).max(p.col_lower[j])
                        };
                        shifts[j].push(i64::from(dir) * t);
                    }
                }
            } else {
                // All rows satisfied: round the fractional integer whose
                // rounding direction carries the most locks (cheapest
                // objective change on ties).
                let mut best: Option<(u32, f64, usize, f64)> = None;
                for (pos, &j) in fractional.iter().enumerate() {
                    let v = y[j];
                    for (locks, r) in [
                        (self.locks.up[j], (v + ft).floor()),
                        (self.locks.down[j], (v - ft).ceil()),
                    ] {
                        let dc = p.cost[j] * (r - v);
                        if best.is_none_or(|b| locks > b.0 || (locks == b.0 && dc < b.1)) {
                            best = Some((locks, dc, pos, r));
                        }
                    }
                }
                if let Some((_, _, pos, r)) = best {
                    y[fractional[pos]] = r;
                    is_frac[fractional[pos]] = false;
                    fractional.remove(pos);
                    reduced = true;
                }
            }
            rows = infeasible_rows(&y);
            if rows.len() < prev_rows || reduced {
                stall = 0;
            }
            prev_rows = rows.len();
        }
        if !rows.is_empty() {
            self.try_rounded_point(&y, Source::Shifting);
        } else if !fractional.is_empty() {
            self.zi_round(&y);
        } else {
            self.found.push((y, Source::Shifting));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::RusolveEngine;
    use rusolve_core::CscMatrix;
    use rusolve_lp::SimplexConfig;

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

    fn with_ctx<R>(p: &MipProblem, f: impl FnOnce(&mut Ctx<'_, RusolveEngine>) -> R) -> R {
        let locks = Locks::new(p);
        let mut prop = Propagator::new(p, 1e-6);
        let cfg = SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        };
        let mut engine = RusolveEngine::new(p, cfg);
        let mut rng = XorShift64::new(1);
        let mut ctx = Ctx {
            p,
            locks: &locks,
            prop: &mut prop,
            engine: &mut engine,
            lower: &p.col_lower,
            upper: &p.col_upper,
            cutoff: f64::INFINITY,
            feastol: 1e-6,
            rng: &mut rng,
            lp_iterations: 0,
            found: Vec::new(),
            deadline: None,
        };
        f(&mut ctx)
    }

    fn feasible(p: &MipProblem, x: &[f64]) -> bool {
        let v = p.violation(x);
        v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6
    }

    #[test]
    fn simple_rounding_uses_lock_free_side() {
        // Packing row x0 + x1 ≤ 1.5: no down-locks ⇒ round down.
        let p = mip(
            &[(&[1.0, 1.0], f64::NEG_INFINITY, 1.5)],
            &[0.0; 2],
            &[1.0; 2],
            &[true; 2],
            &[-1.0; 2],
        );
        let found = with_ctx(&p, |c| {
            c.simple_rounding(&[0.75, 0.75]);
            c.found.clone()
        });
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].0, vec![0.0, 0.0]);
        // Equality rows lock both ways ⇒ nothing.
        let q = mip(
            &[(&[1.0, 1.0], 1.0, 1.0)],
            &[0.0; 2],
            &[1.0; 2],
            &[true; 2],
            &[1.0; 2],
        );
        assert!(with_ctx(&q, |c| {
            c.simple_rounding(&[0.5, 0.5]);
            c.found.is_empty()
        }));
    }

    #[test]
    fn randomized_rounding_completes_continuous_part_by_lp() {
        // x int ∈ [0,3], y cont ∈ [0,10]: x + y ≥ 2.5, x − y ≤ 0.
        let p = mip(
            &[
                (&[1.0, 1.0], 2.5, f64::INFINITY),
                (&[1.0, -1.0], f64::NEG_INFINITY, 0.0),
            ],
            &[0.0, 0.0],
            &[3.0, 10.0],
            &[true, false],
            &[1.0, 1.0],
        );
        let found = with_ctx(&p, |c| {
            c.randomized_rounding(&[1.25, 1.25]);
            c.found.clone()
        });
        assert_eq!(found.len(), 1);
        assert!(feasible(&p, &found[0].0), "{:?}", found[0].0);
    }

    #[test]
    fn zi_round_moves_within_slack() {
        // x0 + x1 ≤ 3 with LP point (0.4, 1.6): both can reach integers.
        let p = mip(
            &[(&[1.0, 1.0], f64::NEG_INFINITY, 3.0)],
            &[0.0; 2],
            &[2.0; 2],
            &[true; 2],
            &[1.0, 1.0],
        );
        let found = with_ctx(&p, |c| {
            c.zi_round(&[0.4, 1.6]);
            c.found.clone()
        });
        assert_eq!(found.len(), 1);
        assert!(feasible(&p, &found[0].0), "{:?}", found[0].0);
    }

    #[test]
    fn shifting_repairs_a_violated_cover_row() {
        // Cover rows: x0 + x1 ≥ 1, x1 + x2 ≥ 1, from LP point (0.5, 0.5, 0.5).
        let inf = f64::INFINITY;
        let p = mip(
            &[(&[1.0, 1.0, 0.0], 1.0, inf), (&[0.0, 1.0, 1.0], 1.0, inf)],
            &[0.0; 3],
            &[1.0; 3],
            &[true; 3],
            &[1.0; 3],
        );
        let found = with_ctx(&p, |c| {
            c.shifting(&[0.5, 0.5, 0.5]);
            c.found.clone()
        });
        assert!(!found.is_empty());
        assert!(found.iter().any(|(x, _)| feasible(&p, x)), "{found:?}");
    }
}
