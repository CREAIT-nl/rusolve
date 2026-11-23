//! LP diving: from a node's LP solution, repeatedly bound one fractional
//! integer column towards a rounding, propagate, and re-solve the LP, until
//! the LP point is integral (a candidate) or the dive fails. One-level
//! backtracking: when a bound change makes the domain or LP infeasible, the
//! opposite rounding of the same column is tried before giving up.
//!
//! Provenance: T. Achterberg, *Constraint Integer Programming*, PhD thesis
//! TU Berlin (2007), §9.2 (diving heuristics, single backtrack); T. Berthold,
//! "Primal heuristics for mixed integer programs", Diploma thesis TU Berlin
//! (2006), §3.1 — selection rules:
//!
//! * *fractional diving*: the column with the smallest fractionality,
//!   rounded to the nearest integer;
//! * *coefficient diving*: the column whose rounding direction with fewer
//!   locks has the fewest locks (≥ 1; lock-free columns are left to simple
//!   rounding), rounded that way; ties by smaller fractionality.
//!
//! Every intermediate LP point also goes through simple rounding, as in the
//! SCIP diving loop (Achterberg 2007 §9.2). HiGHS has no LP diving (it uses
//! RENS/RINS sub-MIPs instead), so no HiGHS anchor here.

use super::Source;
use super::rounding::Ctx;
use crate::engine::{LpEngine, LpStatus};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiveRule {
    Fractional,
    Coefficient,
}

impl<E: LpEngine> Ctx<'_, E> {
    /// Dives from `x`, the LP solution under the local box `lower`/`upper`
    /// with objective `obj`; stops after `max_lp_iterations` LP iterations.
    /// Returns the dive's final LP point when it reached one with no
    /// fractional integer column (the integral point the caller can use as
    /// a sub-MIP fixing reference), `None` when the dive failed or hit
    /// its budget first.
    pub fn dive(
        &mut self,
        rule: DiveRule,
        lower: &[f64],
        upper: &[f64],
        x: &[f64],
        max_lp_iterations: u64,
    ) -> Option<Vec<f64>> {
        let (mut lo, mut hi) = (lower.to_vec(), upper.to_vec());
        let mut x = x.to_vec();
        let start_iters = self.lp_iterations;
        let num_int = self.p.integer.iter().filter(|&&b| b).count();
        for _ in 0..=num_int {
            let Some((j, up)) = self.select(rule, &x) else {
                self.found.push((x.clone(), Source::Diving));
                return Some(x);
            };
            if self.lp_iterations - start_iters > max_lp_iterations {
                return None;
            }
            let mut next = None;
            for dir_up in [up, !up] {
                let (mut l, mut h) = (lo.clone(), hi.clone());
                if dir_up {
                    l[j] = x[j].ceil();
                } else {
                    h[j] = x[j].floor();
                }
                if l[j] > h[j] || !self.prop.propagate(&mut l, &mut h, Some(j), self.cutoff) {
                    continue;
                }
                let Ok(out) = self.engine.solve(&l, &h) else {
                    return None;
                };
                self.lp_iterations += out.iterations;
                if out.status == LpStatus::Optimal && out.objective <= self.cutoff {
                    next = Some((l, h, out.x));
                    break;
                }
            }
            let (l, h, y) = next?;
            (lo, hi, x) = (l, h, y);
            self.simple_rounding(&x);
        }
        None
    }

    /// Branching column and direction (`true` = up) for the dive, `None`
    /// when `x` is integral.
    fn select(&self, rule: DiveRule, x: &[f64]) -> Option<(usize, bool)> {
        let ft = self.feastol;
        let mut best: Option<(u32, f64, usize, bool)> = None;
        for (j, &v) in x.iter().enumerate() {
            if !self.p.integer[j] {
                continue;
            }
            let f = v - v.floor();
            if f <= ft || f >= 1.0 - ft {
                continue;
            }
            let (locks, frac, up) = match rule {
                DiveRule::Fractional => (0, f.min(1.0 - f), f >= 0.5),
                DiveRule::Coefficient => {
                    let (lu, ld) = (self.locks.up[j], self.locks.down[j]);
                    if lu == 0 || ld == 0 {
                        // Trivially roundable: simple rounding's job.
                        continue;
                    }
                    if lu < ld || (lu == ld && f >= 0.5) {
                        (lu, 1.0 - f, true)
                    } else {
                        (ld, f, false)
                    }
                }
            };
            if best.is_none_or(|b| (locks, frac) < (b.0, b.1)) {
                best = Some((locks, frac, j, up));
            }
        }
        match best {
            Some((_, _, j, up)) => Some((j, up)),
            // Coefficient diving skips lock-free columns: fall back to the
            // nearest rounding so the dive still terminates on them.
            None if rule == DiveRule::Coefficient => self.select(DiveRule::Fractional, x),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::RusolveEngine;
    use crate::heuristics::Locks;
    use crate::problem::MipProblem;
    use crate::propagate::Propagator;
    use rusolve_core::CscMatrix;
    use rusolve_lp::SimplexConfig;
    use rusolve_lp::rng::XorShift64;

    #[test]
    fn dives_reach_integral_points_on_a_knapsack_pair() {
        // max 5a + 4b + 3c s.t. 2a + 3b + c ≤ 4.5, 4a + b + 2c ≥ 3, binaries.
        let inf = f64::INFINITY;
        let trip = [
            (0, 0, 2.0),
            (0, 1, 3.0),
            (0, 2, 1.0),
            (1, 0, 4.0),
            (1, 1, 1.0),
            (1, 2, 2.0),
        ];
        let p = MipProblem {
            num_col: 3,
            num_row: 2,
            matrix: CscMatrix::from_triplets(2, 3, &trip),
            row_lower: vec![-inf, 3.0],
            row_upper: vec![4.5, inf],
            cost: vec![-5.0, -4.0, -3.0],
            maximize: false,
            col_lower: vec![0.0; 3],
            col_upper: vec![1.0; 3],
            integer: vec![true; 3],
        };
        let locks = Locks::new(&p);
        let cfg = SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        };
        for rule in [DiveRule::Fractional, DiveRule::Coefficient] {
            let mut prop = Propagator::new(&p, 1e-6);
            let mut engine = RusolveEngine::new(&p, cfg);
            let root = engine.solve(&p.col_lower, &p.col_upper).unwrap();
            let mut rng = XorShift64::new(1);
            let mut ctx = Ctx {
                p: &p,
                locks: &locks,
                prop: &mut prop,
                engine: &mut engine,
                lower: &p.col_lower,
                upper: &p.col_upper,
                cutoff: inf,
                feastol: 1e-6,
                rng: &mut rng,
                lp_iterations: 0,
                found: Vec::new(),
                deadline: None,
            };
            ctx.dive(rule, &p.col_lower, &p.col_upper, &root.x, 1000);
            let dived: Vec<_> = ctx
                .found
                .iter()
                .filter(|(_, s)| *s == Source::Diving)
                .collect();
            assert_eq!(dived.len(), 1, "{rule:?}");
            let v = p.violation(&dived[0].0);
            assert!(v.row <= 1e-6 && v.integrality <= 1e-6, "{rule:?} {v:?}");
        }
    }

    #[test]
    fn dive_returns_the_integral_point_it_reached() {
        // Same knapsack pair: the returned point is the found Diving point,
        // integral on every integer column; a dive from an already-integral
        // point returns it unchanged.
        let inf = f64::INFINITY;
        let trip = [
            (0, 0, 2.0),
            (0, 1, 3.0),
            (0, 2, 1.0),
            (1, 0, 4.0),
            (1, 1, 1.0),
            (1, 2, 2.0),
        ];
        let p = MipProblem {
            num_col: 3,
            num_row: 2,
            matrix: CscMatrix::from_triplets(2, 3, &trip),
            row_lower: vec![-inf, 3.0],
            row_upper: vec![4.5, inf],
            cost: vec![-5.0, -4.0, -3.0],
            maximize: false,
            col_lower: vec![0.0; 3],
            col_upper: vec![1.0; 3],
            integer: vec![true; 3],
        };
        let locks = Locks::new(&p);
        let cfg = SimplexConfig {
            presolve: false,
            ..SimplexConfig::default()
        };
        let mut prop = Propagator::new(&p, 1e-6);
        let mut engine = RusolveEngine::new(&p, cfg);
        let root = engine.solve(&p.col_lower, &p.col_upper).unwrap();
        let mut rng = XorShift64::new(1);
        let mut ctx = Ctx {
            p: &p,
            locks: &locks,
            prop: &mut prop,
            engine: &mut engine,
            lower: &p.col_lower,
            upper: &p.col_upper,
            cutoff: inf,
            feastol: 1e-6,
            rng: &mut rng,
            lp_iterations: 0,
            found: Vec::new(),
            deadline: None,
        };
        let point = ctx
            .dive(
                DiveRule::Coefficient,
                &p.col_lower,
                &p.col_upper,
                &root.x,
                1000,
            )
            .expect("dive reaches an integral point");
        assert_eq!(point.len(), p.num_col);
        assert!(point.iter().all(|v| (v - v.round()).abs() <= 1e-6));
        // An integral start needs no step at all.
        let again = ctx.dive(
            DiveRule::Fractional,
            &p.col_lower,
            &p.col_upper,
            &point,
            1000,
        );
        assert_eq!(again.as_deref(), Some(point.as_slice()));
    }
}
