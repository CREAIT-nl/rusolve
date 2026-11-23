//! HiGHS source: ipm/ipx/ipm.cc (lines 450-895)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::basis::Basis;
use super::super::control::{fixed, format_int, format_str, scientific, time_fmt};
use super::super::iterate::State;
use super::super::model::Model;
use super::super::sparse_matrix::dot_column;
use super::super::timer::Timer;
use super::super::utils::all_finite;
use super::part1::step_to_boundary;
use super::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

impl IPM {
    /// HiGHS: ipx/ipm.cc::IPM::Centring (lines 450-487)
    pub(super) fn centring(&mut self, r: &mut IpmRefs<'_>, step: &mut Step, mu_to_use: f64) {
        let model_rc = r.iterate.model_rc();
        let model: &Model = &model_rc;
        let m = model.rows();
        let n = model.cols();
        let nm = n + m;
        let mut sl = Vector::new(nm);
        let mut su = Vector::new(nm);

        // Set sigma to 1 for pure centring
        let sigma = 1.0;

        {
            let it: &Iterate = &*r.iterate;
            let (xl, xu, zl, zu) = (it.xl(), it.xu(), it.zl(), it.zu());
            // sl = -xl.*zl + sigma*mu
            for j in 0..nm {
                let ju = j as usize;
                sl[ju] = if it.has_barrier_lb(j) {
                    -xl[ju] * zl[ju] + sigma * mu_to_use
                } else {
                    0.0
                };
            }
            debug_assert!(all_finite(&sl));

            // su = -xu.*zu + sigma*mu
            for j in 0..nm {
                let ju = j as usize;
                su[ju] = if it.has_barrier_ub(j) {
                    -xu[ju] * zu[ju] + sigma * mu_to_use
                } else {
                    0.0
                };
            }
            debug_assert!(all_finite(&su));
        }
        let rb = r.iterate.rb().clone();
        let rc = r.iterate.rc().clone();
        let rl = r.iterate.rl().clone();
        let ru = r.iterate.ru().clone();
        self.solve_newton_system(
            &mut *r.kkt,
            &*r.iterate,
            &mut *r.info,
            &rb,
            Some(&rc[..]),
            &rl,
            &ru,
            &sl,
            &su,
            step,
        );
    }

    /// HiGHS: ipx/ipm.cc::IPM::AssessCentrality (lines 489-544)
    ///
    /// C++ default arguments: `print = true`.
    pub(super) fn assess_centrality(
        &mut self,
        r: &mut IpmRefs<'_>,
        xl: &Vector,
        xu: &Vector,
        zl: &Vector,
        zu: &Vector,
        mu: f64,
        print: bool,
    ) {
        // The function computes the ratio
        //  min(x_j * z_j) / max(x_j * z_j)
        // and the number of products x_j * z_j that are not in the interval
        //  [0.1 * mu, 10 * mu]
        // and prints information to screen if print is on.
        let m = r.iterate.model().rows();
        let n = r.iterate.model().cols();

        let mut minxz = f64::INFINITY; // kHighsInf
        let mut maxxz = 0.0;

        let gamma = 0.1;
        self.bad_products = 0;

        for j in 0..n + m {
            let ju = j as usize;
            if r.iterate.has_barrier_lb(j) {
                let product = xl[ju] * zl[ju];
                if product < gamma * mu || product > mu / gamma {
                    self.bad_products += 1;
                }
                minxz = cpp_min(minxz, product);
                maxxz = cpp_max(maxxz, product);
            }
        }

        for j in 0..n + m {
            let ju = j as usize;
            if r.iterate.has_barrier_ub(j) {
                let product = xu[ju] * zu[ju];
                if product < gamma * mu || product > mu / gamma {
                    self.bad_products += 1;
                }
                minxz = cpp_min(minxz, product);
                maxxz = cpp_max(maxxz, product);
            }
        }

        maxxz = cpp_max(maxxz, mu);
        minxz = cpp_min(minxz, mu);

        self.centring_ratio = maxxz / minxz;

        if print {
            let mut s = format!(
                "\txj*zj in [ {}, {}]; Ratio = {}; (xj*zj / mu) not_in [0.1, 10]: {}\n",
                scientific(minxz / mu, 8, 2),
                scientific(maxxz / mu, 8, 2),
                scientific(self.centring_ratio, 8, 2),
                self.bad_products
            );
            self.control_.h_log(&mut s);
        }
    }

    /// HiGHS: ipx/ipm.cc::IPM::EvaluateCentringStep (lines 546-600)
    pub(super) fn evaluate_centring_step(
        &mut self,
        r: &mut IpmRefs<'_>,
        step: &Step,
        prev_ratio: f64,
        prev_bad: Int,
    ) -> bool {
        // The function returns true is the step is to be accepted.
        // The step is accepted if the ratio of the new point is not worse
        // than the previous one times centringRatioReduction or if the
        // number of outliers products is reduced.
        self.step_sizes(r, step, true);

        let n = r.iterate.model().cols();
        let m = r.iterate.model().rows();

        let mut xl_temp = r.iterate.xl().clone();
        let mut xu_temp = r.iterate.xu().clone();
        let mut zl_temp = r.iterate.zl().clone();
        let mut zu_temp = r.iterate.zu().clone();

        // perform temporary step
        for j in 0..n + m {
            let ju = j as usize;
            if r.iterate.has_barrier_lb(j) {
                xl_temp[ju] += self.step_primal_ * step.xl[ju];
            }
            if r.iterate.has_barrier_ub(j) {
                xu_temp[ju] += self.step_primal_ * step.xu[ju];
            }
            if r.iterate.has_barrier_lb(j) {
                zl_temp[ju] += self.step_dual_ * step.zl[ju];
            }
            if r.iterate.has_barrier_ub(j) {
                zu_temp[ju] += self.step_dual_ * step.zu[ju];
            }
        }

        // compute temporary mu
        let mut mu_temp = 0.0;
        let mut num_finite: Int = 0;
        for j in 0..n + m {
            let ju = j as usize;
            if r.iterate.has_barrier_lb(j) {
                mu_temp += xl_temp[ju] * zl_temp[ju];
                num_finite += 1;
            }
            if r.iterate.has_barrier_ub(j) {
                mu_temp += xu_temp[ju] * zu_temp[ju];
                num_finite += 1;
            }
        }
        mu_temp /= num_finite as f64;

        // assess quality of temporary point
        self.assess_centrality(r, &xl_temp, &xu_temp, &zl_temp, &zu_temp, mu_temp, false);

        // accept the step if the new ratio is not more than centringRatioReduction
        // times the previous one, or if the new point has fewer outliers
        self.centring_ratio < self.control_.centring_ratio_reduction() * prev_ratio
            || self.bad_products < prev_bad
    }

    /// HiGHS: ipx/ipm.cc::IPM::StepSizes (lines 602-690)
    ///
    /// C++ default arguments: `isCentring = false`.
    pub(super) fn step_sizes(&mut self, r: &mut IpmRefs<'_>, step: &Step, isCentring: bool) {
        let model_rc = r.iterate.model_rc();
        let m = model_rc.rows();
        let n = model_rc.cols();
        let nm = n + m;
        let it: &Iterate = &*r.iterate;
        let (xl, xu, zl, zu) = (it.xl(), it.xu(), it.zl(), it.zu());
        let (dxl, dxu, dzl, dzu) = (&step.xl, &step.xu, &step.zl, &step.zu);
        let gammaf = 0.9;
        let gammaa = 1.0 / (1.0 - gammaf);

        let mut block_xl: Int = -1;
        let mut block_xu: Int = -1;
        let mut block_zl: Int = -1;
        let mut block_zu: Int = -1;
        let step_xl = step_to_boundary(xl, dxl, Some(&mut block_xl), 1.0);
        let step_xu = step_to_boundary(xu, dxu, Some(&mut block_xu), 1.0);
        let step_zl = step_to_boundary(zl, dzl, Some(&mut block_zl), 1.0);
        let step_zu = step_to_boundary(zu, dzu, Some(&mut block_zu), 1.0);
        let maxp = f64::min(step_xl, step_xu);
        let maxd = f64::min(step_zl, step_zu);
        let mut mufull = 0.0;
        let mut num_finite: Int = 0;
        for j in 0..nm {
            let ju = j as usize;
            if it.has_barrier_lb(j) {
                debug_assert!(xl[ju].is_finite());
                debug_assert!(xl[ju] != 0.0);
                mufull += (xl[ju] + maxp * dxl[ju]) * (zl[ju] + maxd * dzl[ju]);
                num_finite += 1;
            }
            if it.has_barrier_ub(j) {
                debug_assert!(xu[ju].is_finite());
                debug_assert!(xu[ju] != 0.0);
                mufull += (xu[ju] + maxp * dxu[ju]) * (zu[ju] + maxd * dzu[ju]);
                num_finite += 1;
            }
        }
        debug_assert!(mufull.is_finite());
        mufull /= num_finite as f64;
        mufull /= gammaa;

        let mut alphap = 1.0;
        let mut alphad = 1.0;
        let mut blockp: Int = -1;
        let mut blockd: Int = -1;
        if maxp < 1.0 {
            if step_xl <= step_xu {
                blockp = block_xl;
                let bp = blockp as usize;
                let buffer = mufull / (zl[bp] + maxd * dzl[bp]);
                alphap = (xl[bp] - buffer) / (-dxl[bp]);
            } else {
                blockp = block_xu;
                let bp = blockp as usize;
                let buffer = mufull / (zu[bp] + maxd * dzu[bp]);
                alphap = (xu[bp] - buffer) / (-dxu[bp]);
            }
            alphap = cpp_max(alphap, gammaf * maxp);
            alphap = cpp_min(alphap, 1.0);
            debug_assert!(blockp >= 0);
        }
        if maxd < 1.0 {
            if step_zl <= step_zu {
                blockd = block_zl;
                let bd = blockd as usize;
                let buffer = mufull / (xl[bd] + maxp * dxl[bd]);
                alphad = (zl[bd] - buffer) / (-dzl[bd]);
            } else {
                blockd = block_zu;
                let bd = blockd as usize;
                let buffer = mufull / (xu[bd] + maxp * dxu[bd]);
                alphad = (zu[bd] - buffer) / (-dzu[bd]);
            }
            alphad = cpp_max(alphad, gammaf * maxd);
            alphad = cpp_min(alphad, 1.0);
            debug_assert!(blockd >= 0);
        }
        let _ = (blockp, blockd);
        self.step_primal_ = cpp_min(alphap, 1.0 - 1e-6);
        self.step_dual_ = cpp_min(alphad, 1.0 - 1e-6);

        if isCentring {
            // When computing stepsizes for a centring step, reduce them
            // by centringAlphaScaling. This ensures that the point is
            // well centred and does not get too close to the boundary.
            self.step_primal_ = alphap * self.control_.centring_alpha_scaling();
            self.step_dual_ = alphad * self.control_.centring_alpha_scaling();
        }
    }

    /// HiGHS: ipx/ipm.cc::IPM::MakeStep (lines 692-704)
    /// HiGHS 1.16 (#3401): ipm/ipx/ipm.cc::IPM::MakeStep (lines 694-738)
    ///
    /// C++ default arguments: `isCentring = false`.
    pub(super) fn make_step(&mut self, r: &mut IpmRefs<'_>, step: &Step, isCentring: bool) {
        self.step_sizes(r, step, isCentring);
        // HiGHS 1.16 #3401 (h16_ipx_stall): ipm/ipx/ipm.cc lines 697-699 (1.16.0). The
        // residuals are read only on the flag path: presidual() and friends
        // evaluate the iterate lazily, so an extra call is a side effect.
        let h16_ipx_stall = self.control_.h16_ipx_stall();
        let (h16_presidual, h16_dresidual, h16_complementarity) = if h16_ipx_stall {
            (
                r.iterate.presidual(),
                r.iterate.dresidual(),
                r.iterate.complementarity(),
            )
        } else {
            (0.0, 0.0, 0.0)
        };
        r.iterate.update(
            self.step_primal_,
            Some(&step.x[..]),
            Some(&step.xl[..]),
            Some(&step.xu[..]),
            self.step_dual_,
            Some(&step.y[..]),
            Some(&step.zl[..]),
            Some(&step.zu[..]),
        );
        if !isCentring {
            // HiGHS 1.16 #3401 (h16_ipx_stall): ipm/ipx/ipm.cc lines 705-731 (1.16.0).
            // A long step that barely reduces a residual, or a long step that
            // leaves the complementarity gap unchanged, counts as a bad iteration.
            let bad_iter = if h16_ipx_stall {
                let large_step_thresh = 0.9;
                let reduction_thresh = 0.9;

                let model = r.iterate.model();
                let tol = r.iterate.feasibility_tol();
                let presidual_stalled = self.step_primal_ >= large_step_thresh
                    && r.iterate.presidual() > reduction_thresh * h16_presidual
                    && r.iterate.presidual() > tol * (1.0 + model.norm_bounds());
                let dresidual_stalled = self.step_dual_ >= large_step_thresh
                    && r.iterate.dresidual() > reduction_thresh * h16_dresidual
                    && r.iterate.dresidual() > tol * (1.0 + model.norm_c());
                let complementarity_stalled = cpp_min(self.step_primal_, self.step_dual_)
                    >= large_step_thresh
                    && r.iterate.complementarity() > reduction_thresh * h16_complementarity;

                cpp_min(self.step_primal_, self.step_dual_) < 0.05
                    || presidual_stalled
                    || dresidual_stalled
                    || complementarity_stalled
            } else {
                cpp_min(self.step_primal_, self.step_dual_) < 0.05
            };
            if bad_iter {
                self.num_bad_iter_ += 1;
            } else {
                self.num_bad_iter_ = 0;
            }
            self.best_complementarity_ =
                cpp_min(self.best_complementarity_, r.iterate.complementarity());
        }
    }

    /// HiGHS: ipx/ipm.cc::IPM::SolveNewtonSystem (lines 706-819)
    ///
    /// PORT NOTE: takes the three `IpmRefs` members separately, with
    /// `iterate` shared, because the callers pass `&iterate_->rb()[0]` etc.
    /// (a `Ref` borrowed from the iterate) as right-hand sides:
    /// `let IpmRefs { kkt, iterate, info } = r; let rb = iterate.rb(); ...`.
    /// Only `rc` is null-checked in C++ (line 727).
    pub(super) fn solve_newton_system(
        &mut self,
        kkt: &mut dyn KKTSolver,
        iterate: &Iterate,
        info: &mut Info,
        rb: &[f64],
        rc: Option<&[f64]>,
        rl: &[f64],
        ru: &[f64],
        sl: &[f64],
        su: &[f64],
        lhs: &mut Step,
    ) {
        let model_rc = iterate.model_rc();
        let m = model_rc.rows();
        let n = model_rc.cols();
        let nm = n + m;
        let ai = model_rc.AI();
        let xl = iterate.xl();
        let xu = iterate.xu();
        let zl = iterate.zl();
        let zu = iterate.zu();
        let Step {
            x: dx,
            xl: dxl,
            xu: dxu,
            y: dy,
            zl: dzl,
            zu: dzu,
        } = lhs;

        // Build RHS for KKT system.
        let mut rhs1 = Vector::new(nm);
        let mut rhs2 = Vector::new(m);
        if let Some(rc) = rc {
            for j in 0..nm as usize {
                rhs1[j] = -rc[j];
            }
        }
        for j in 0..nm {
            let ju = j as usize;
            let rlj = rl[ju];
            let ruj = ru[ju];
            if iterate.has_barrier_lb(j) {
                rhs1[ju] += (sl[ju] + zl[ju] * rlj) / xl[ju];
            }
            if iterate.has_barrier_ub(j) {
                rhs1[ju] -= (su[ju] - zu[ju] * ruj) / xu[ju];
            }
            if iterate.state_of(j) == State::fixed {
                rhs1[ju] = 0.0;
            }
        }
        debug_assert!(all_finite(&rhs1));
        for i in 0..m as usize {
            rhs2[i] = rb[i];
        }

        // Solve KKT system.
        let tol = self.control_.kkt_tol() * iterate.mu().sqrt();
        kkt.solve(&rhs1, &rhs2, tol, dx, dy, info);
        if info.errflag != 0 {
            return;
        }

        // Recover solution to Newton system.
        dy.iter_mut().for_each(|v| *v *= -1.0);
        for j in 0..nm {
            let ju = j as usize;
            match iterate.state_of(j) {
                State::fixed => {
                    debug_assert!(dx[ju] == 0.0);
                    dxl[ju] = 0.0;
                    dzl[ju] = 0.0;
                }
                State::free => {
                    debug_assert!(rl[ju] == 0.0);
                    dxl[ju] = 0.0;
                    dzl[ju] = 0.0;
                }
                _ => {
                    let rlj = rl[ju];
                    dxl[ju] = dx[ju] - rlj;
                    dzl[ju] = (sl[ju] - zl[ju] * dxl[ju]) / xl[ju];
                }
            }
        }
        for j in 0..nm {
            let ju = j as usize;
            match iterate.state_of(j) {
                State::fixed => {
                    debug_assert!(dx[ju] == 0.0);
                    dxu[ju] = 0.0;
                    dzu[ju] = 0.0;
                }
                State::free => {
                    debug_assert!(ru[ju] == 0.0);
                    dxu[ju] = 0.0;
                    dzu[ju] = 0.0;
                }
                _ => {
                    let ruj = ru[ju];
                    dxu[ju] = ruj - dx[ju];
                    dzu[ju] = (su[ju] - zu[ju] * dxu[ju]) / xu[ju];
                }
            }
        }
        debug_assert!(all_finite(dxl));
        debug_assert!(all_finite(dxu));
        debug_assert!(all_finite(dzl));
        debug_assert!(all_finite(dzu));

        // Shift residual to the last two block equations.
        for j in 0..nm {
            let ju = j as usize;
            if iterate.state_of(j) == State::barrier {
                debug_assert!(xl[ju].is_finite() || xu[ju].is_finite());
                let atdy = dot_column(ai, j, dy);
                let rcj = match rc {
                    Some(rc) => rc[ju],
                    None => 0.0,
                };
                if xl[ju].is_finite() && xu[ju].is_finite() {
                    if zl[ju] * xu[ju] >= zu[ju] * xl[ju] {
                        dzl[ju] = rcj + dzu[ju] - atdy;
                    } else {
                        dzu[ju] = -rcj + dzl[ju] + atdy;
                    }
                } else if xl[ju].is_finite() {
                    dzl[ju] = rcj + dzu[ju] - atdy;
                } else {
                    dzu[ju] = -rcj + dzl[ju] + atdy;
                }
            }
        }

        // Check solution for free and fixed variables (debug builds only in C++).
        for j in 0..nm {
            let ju = j as usize;
            if iterate.state_of(j) == State::fixed {
                debug_assert!(dx[ju] == 0.0);
                debug_assert!(dxl[ju] == 0.0 && dxu[ju] == 0.0);
                debug_assert!(dzl[ju] == 0.0 && dzu[ju] == 0.0);
            }
            if iterate.state_of(j) == State::free {
                debug_assert!(dzl[ju] == 0.0 && dzu[ju] == 0.0);
            }
        }
    }

    /// HiGHS: ipx/ipm.cc::IPM::PrintHeader (lines 821-841)
    pub fn print_header(&mut self) {
        let mut s = String::new();
        s += &format!(
            " {}  {}  {}  {}  {}  {}",
            format_str("Iter", 4),
            format_str("primal obj", 15),
            format_str("dual obj", 15),
            format_str("pinf", 9),
            format_str("dinf", 9),
            format_str("gap", 8)
        );
        if !self.control_.timeless_log() {
            s += &format!("   {}", format_str("time", 7));
        }
        self.control_.h_log(&mut s);
        self.control_.debug(
            1,
            &format!(
                "  {}  {} {}  {} {}",
                format_str("stepsizes", 9),
                format_str("pivots", 7),
                format_str("kktiter", 7),
                format_str("P.fixed", 7),
                format_str("D.fixed", 7)
            ),
        );
        self.control_
            .debug(4, &format!("  {}", format_str("svdmin(B)", 9)));
        self.control_
            .debug(4, &format!("  {}", format_str("density", 8)));
        self.control_.h_log_str("\n");
    }

    /// HiGHS: ipx/ipm.cc::IPM::PrintOutput (lines 843-895)
    pub(super) fn print_output(&mut self, r: &mut IpmRefs<'_>) {
        let ipm_optimal = r.iterate.feasible() && r.iterate.optimal();

        let logging_pobj = r.iterate.pobjective_after_postproc();
        let logging_dobj = r.iterate.dobjective_after_postproc();
        let mut logging_presidual = r.iterate.presidual();
        let mut logging_dresidual = r.iterate.dresidual();

        // Now logging relative primal and dual infeasibility, and also
        // the relative primal dual objective gap
        logging_presidual /= r.iterate.bounds_measure_;
        logging_dresidual /= r.iterate.costs_measure_;
        let logging_gap =
            (logging_pobj - logging_dobj).abs() / (1.0 + 0.5 * (logging_pobj + logging_dobj).abs());

        let mut s = format!(
            " {}{}  {}  {}  {}  {}  {}",
            format_int(r.info.iter, 3),
            if ipm_optimal { "*" } else { " " },
            scientific(logging_pobj, 15, 8),
            scientific(logging_dobj, 15, 8),
            scientific(logging_presidual, 9, 2),
            scientific(logging_dresidual, 9, 2),
            scientific(logging_gap, 8, 2)
        );
        if !self.control_.timeless_log() {
            s += &format!("  {}", time_fmt(self.control_.elapsed()));
        }
        self.control_.h_log(&mut s);
        self.control_.debug(
            1,
            &format!(
                "  {} {}  {} {}",
                fixed(self.step_primal_, 4, 2),
                fixed(self.step_dual_, 4, 2),
                format_int(r.kkt.basis_changes(), 7),
                format_int(r.kkt.iter(), 7)
            ),
        );
        self.control_.debug(
            1,
            &format!(
                "  {} {}",
                format_int(r.info.dual_dropped, 7),
                format_int(r.info.primal_dropped, 7)
            ),
        );

        let basis: Option<&Basis> = r.kkt.basis();
        match basis {
            Some(basis) => {
                if self.control_.debug_enabled(4) {
                    self.control_.debug(
                        4,
                        &format!("  {}", scientific(basis.min_singular_value(), 9, 2)),
                    );
                    let timer = Timer::default();
                    let density = basis.density_inverse();
                    r.info.time_symb_invert += timer.elapsed();
                    self.control_
                        .debug(4, &format!("  {}", scientific(density, 8, 2)));
                }
            }
            None => {
                self.control_.debug(4, &format!("  {}", format_str("-", 9)));
                self.control_.debug(4, &format!("  {}", format_str("-", 8)));
            }
        }
        self.control_.h_log_str("\n");
    }
}

/// `h16_ipx_stall` tests (HiGHS 1.16 #3401, ported from check/TestIpm.cpp test-3397).
#[cfg(test)]
mod h16_ipx_stall_tests {
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::*;

    /// TestIpm.cpp test-3397 (diff lines 13-25): min x + y s.t. y - x >= 3;
    /// x, y free is unbounded.
    fn unbounded_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![1.0, 1.0],
            col_lower_: vec![-kHighsInf, -kHighsInf],
            col_upper_: vec![kHighsInf, kHighsInf],
            row_lower_: vec![3.0],
            row_upper_: vec![kHighsInf],
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![-1.0, 1.0];
        lp
    }

    /// TestIpm.cpp test-3397 (diff lines 32-46): LP with optimal objective -185.875.
    fn optimal_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 6,
            num_row_: 8,
            col_cost_: vec![-5.0, 0.0, 8.0, 0.0, -4.0, -7.0],
            col_lower_: vec![3.0, -kHighsInf, -5.0, -kHighsInf, -5.0, -5.0],
            col_upper_: vec![kHighsInf, kHighsInf, -1.0, kHighsInf, -4.0, kHighsInf],
            row_lower_: vec![
                4.0, 11.0, -16.875, -kHighsInf, -kHighsInf, -34.0, 1.5, -kHighsInf,
            ],
            row_upper_: vec![kHighsInf, 11.0, -16.875, 16.0, 0.0, -33.5, 1.5, -2.0],
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 6;
        lp.a_matrix_.num_row_ = 8;
        lp.a_matrix_.start_ = vec![0, 3, 6, 6, 10, 12, 14];
        lp.a_matrix_.index_ = vec![0, 3, 5, 1, 3, 6, 2, 4, 5, 7, 5, 6, 2, 4];
        lp.a_matrix_.value_ = vec![
            5.0, -1.0, 1.0, 2.0, 5.0, 1.0, -1.0, 5.0, 7.0, 1.0, 4.0, 1.0, 5.0, -1.0,
        ];
        lp
    }

    /// The TestIpm.cpp option settings (diff lines 19-22), plus the HiGHS 1.16 flag.
    fn ipx_highs(lp: HighsLp, h16_ipx_stall: bool) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", kIpxString),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", "off"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_double("time_limit", 10.0),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_ipx_stall", h16_ipx_stall),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    /// Runs the unbounded LP; returns the model status and the IPM iteration count.
    fn run_unbounded(h16_ipx_stall: bool) -> (HighsModelStatus, HighsInt) {
        let mut highs = ipx_highs(unbounded_lp(), h16_ipx_stall);
        highs.run();
        (
            highs.get_model_status(),
            highs.get_info().ipm_iteration_count,
        )
    }

    /// TestIpm.cpp diff lines 50-57: optimal with objective -185.875 for both
    /// crossover settings.
    fn check_optimal_lp(h16_ipx_stall: bool) {
        let mut highs = ipx_highs(optimal_lp(), h16_ipx_stall);
        for crossover in ["on", "off"] {
            assert_eq!(
                highs.set_option_value_str("run_crossover", crossover),
                HighsStatus::kOk
            );
            highs.clear_solver();
            highs.run();
            assert_eq!(
                highs.get_model_status(),
                HighsModelStatus::kOptimal,
                "{crossover}"
            );
            let objective = highs.get_info().objective_function_value;
            assert!(
                (objective + 185.875).abs() < 1e-6,
                "{crossover}: {objective}"
            );
        }
    }

    /// Flag on: the 1.16 test-3397 expectations (kUnbounded; optimal -185.875).
    /// The stall tests fire on this LP: the IPM stops in fewer iterations than
    /// with the 1.15.1 rule, which shows that the new branch is reached.
    #[test]
    fn h16_ipx_stall_test_3397_on() {
        let (status_on, iter_on) = run_unbounded(true);
        assert_eq!(status_on, HighsModelStatus::kUnbounded);
        let (_, iter_off) = run_unbounded(false);
        assert!(iter_on < iter_off, "{iter_on} vs {iter_off}");
        check_optimal_lp(true);
    }

    /// Flag off: the 1.15.1 rule. On this port it already reports kUnbounded
    /// (at IPM iteration 9) and solves the second LP to optimality.
    #[test]
    fn h16_ipx_stall_test_3397_off() {
        assert_eq!(run_unbounded(false), (HighsModelStatus::kUnbounded, 9));
        check_optimal_lp(false);
    }
}
