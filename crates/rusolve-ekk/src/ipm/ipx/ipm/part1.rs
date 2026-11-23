//! HiGHS source: ipm/ipx/ipm.cc (lines 11-448)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::model::Model;
use super::super::oracle_trace;
use super::super::sparse_matrix::{multiply_add, scatter_column};
use super::super::utils::{all_finite, infnorm, twonorm};
use super::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// `valarray += valarray` (elementwise, same size).
fn va_add_assign(a: &mut Vector, b: &Vector) {
    for (ai, bi) in a.iter_mut().zip(b.iter()) {
        *ai += *bi;
    }
}

impl Step {
    /// HiGHS: ipx/ipm.cc::IPM::Step::Step (line 12)
    pub fn new(m: Int, n: Int) -> Self {
        Step {
            x: Vector::new(n + m),
            xl: Vector::new(n + m),
            xu: Vector::new(n + m),
            y: Vector::new(m),
            zl: Vector::new(n + m),
            zu: Vector::new(n + m),
        }
    }
    /// HiGHS: ipx/ipm.cc::IPM::Step::operator+= (lines 14-18)
    pub fn add_assign(&mut self, rhs: &Step) {
        va_add_assign(&mut self.x, &rhs.x);
        va_add_assign(&mut self.xl, &rhs.xl);
        va_add_assign(&mut self.xu, &rhs.xu);
        va_add_assign(&mut self.y, &rhs.y);
        va_add_assign(&mut self.zl, &rhs.zl);
        va_add_assign(&mut self.zu, &rhs.zu);
    }
}

impl IPM {
    /// HiGHS: ipx/ipm.cc::IPM::IPM (line 21)
    pub fn new(control: Rc<Control>) -> Self {
        IPM {
            control_: control,
            step_primal_: 0.0,
            step_dual_: 0.0,
            num_bad_iter_: 0,
            best_complementarity_: 0.0,
            maxiter_: -1,
            centring_ratio: 0.0,
            bad_products: 0,
        }
    }

    /// HiGHS: ipx/ipm.cc::IPM::StartingPoint (lines 23-43)
    pub fn starting_point(
        &mut self,
        kkt: &mut dyn KKTSolver,
        iterate: &mut Iterate,
        info: &mut Info,
    ) {
        let mut r = IpmRefs { kkt, iterate, info };
        self.print_header();
        self.compute_starting_point(&mut r);
        if r.info.errflag == 0 {
            self.print_output(&mut r);
            oracle_trace::ipx_iter(
                "start",
                r.iterate,
                r.info,
                &*r.kkt,
                self.step_primal_,
                self.step_dual_,
            );
        }
        // Set status_ipm.
        if r.info.errflag == IPX_ERROR_user_interrupt {
            r.info.errflag = 0;
            r.info.status_ipm = IPX_STATUS_user_interrupt;
        } else if r.info.errflag == IPX_ERROR_time_interrupt {
            r.info.errflag = 0;
            r.info.status_ipm = IPX_STATUS_time_limit;
        } else if r.info.errflag != 0 {
            r.info.status_ipm = IPX_STATUS_failed;
        } else {
            r.info.status_ipm = IPX_STATUS_not_run;
        }
    }

    /// HiGHS: ipx/ipm.cc::IPM::Driver (lines 45-192)
    pub fn driver(&mut self, kkt: &mut dyn KKTSolver, iterate: &mut Iterate, info: &mut Info) {
        let (m, n) = {
            let model = iterate.model();
            (model.rows(), model.cols())
        };
        let mut step = Step::new(m, n);
        self.num_bad_iter_ = 0;
        let mut r = IpmRefs { kkt, iterate, info };

        loop {
            if r.iterate.term_crit_reached() {
                r.info.status_ipm = IPX_STATUS_optimal;
                break;
            }
            if self.num_bad_iter_ >= 5
                || (self.best_complementarity_ > 0.0
                    && r.iterate.complementarity() > Self::kDivergeTol * self.best_complementarity_)
            {
                // No progress in reducing the complementarity gap.
                // Check if model seems to be primal or dual infeasible.
                let dualized = r.iterate.model().dualized();
                let pobjective = r.iterate.pobjective_after_postproc();
                let dobjective = r.iterate.dobjective_after_postproc();
                // HiGHS 1.16 #3513 (h16_ipx_scaled): ipm/ipx/ipm.cc lines 68-69, 74-75 (1.16.0)
                let h16_ipx_scaled = self.control_.h16_ipx_scaled();
                if if h16_ipx_scaled {
                    r.iterate.presidual()
                        > r.iterate.feasibility_tol() * (1.0 + r.iterate.model().norm_bounds())
                        && dobjective > cpp_max(10.0 * pobjective.abs(), 1.0)
                } else {
                    dobjective > cpp_max(10.0 * pobjective.abs(), 1.0)
                } {
                    // Dual objective tends to positive infinity. Looks like the
                    // model is dual unbounded, i.e. primal infeasible.
                    r.info.status_ipm = if dualized {
                        IPX_STATUS_dual_infeas
                    } else {
                        IPX_STATUS_primal_infeas
                    };
                } else if if h16_ipx_scaled {
                    r.iterate.dresidual()
                        > r.iterate.feasibility_tol() * (1.0 + r.iterate.model().norm_c())
                        && pobjective < -cpp_max(10.0 * dobjective.abs(), 1.0)
                } else {
                    pobjective < -cpp_max(10.0 * dobjective.abs(), 1.0)
                } {
                    // Primal objective tends to negative infinity. Looks like the
                    // model is primal unbounded, i.e. dual infeasible.
                    r.info.status_ipm = if dualized {
                        IPX_STATUS_primal_infeas
                    } else {
                        IPX_STATUS_dual_infeas
                    };
                } else {
                    r.info.status_ipm = IPX_STATUS_no_progress;
                }
                break;
            }
            if r.info.iter >= self.maxiter_ {
                r.info.status_ipm = IPX_STATUS_iter_limit;
                break;
            }
            r.info.errflag = self.control_.interrupt_check(r.info.iter);
            if r.info.errflag != 0 {
                break;
            }
            r.kkt.factorize(Some(&mut *r.iterate), &mut *r.info);
            if r.info.errflag != 0 {
                break;
            }
            self.predictor(&mut r, &mut step);
            if r.info.errflag != 0 {
                break;
            }
            self.add_corrector(&mut r, &mut step);
            if r.info.errflag != 0 {
                break;
            }
            self.make_step(&mut r, &step, false);
            r.info.iter += 1;
            self.print_output(&mut r);
            oracle_trace::ipx_iter(
                "iter",
                r.iterate,
                r.info,
                &*r.kkt,
                self.step_primal_,
                self.step_dual_,
            );
        }

        // Set status_ipm if errflag terminated IPM.
        if r.info.errflag != 0 {
            if r.info.errflag == IPX_ERROR_user_interrupt {
                r.info.errflag = 0;
                r.info.status_ipm = IPX_STATUS_user_interrupt;
            } else if r.info.errflag == IPX_ERROR_time_interrupt {
                r.info.errflag = 0;
                r.info.status_ipm = IPX_STATUS_time_limit;
            } else {
                r.info.status_ipm = IPX_STATUS_failed;
            }
        }

        // Centrality of a point is evaluated by the quantities
        //  min (xj * zj) / mu
        //  max (xj * zj) / mu
        // Ideally, they are in the interval [0.1,10.0]. As soon as the ratio
        //  max (xj * zj) / min (xj * zj)
        // is below centringRatioTolerance, the point is considered centred.
        // If the new point after centring has a ratio that is lower than the
        // previous ratio times centringRatioReduction, then the step is
        // accepted. Otherwise, the step is rejected and no more centring steps
        // are performed.
        if self.control_.run_centring() != 0
            && r.info.status_ipm == IPX_STATUS_optimal
            && r.info.centring_tried == 0
        {
            // If IPM is optimal and centring has not yet run, run centring
            // (to avoid running it twice during initial IPM and main IPM).
            self.control_.h_log_str("Performing centring steps...\n");

            // freeze mu to its current value
            let mu_frozen = r.iterate.mu();

            // assess and print centrality of current point
            let (xl, xu, zl, zu) = (
                r.iterate.xl().clone(),
                r.iterate.xu().clone(),
                r.iterate.zl().clone(),
                r.iterate.zu().clone(),
            );
            let mu_now = r.iterate.mu();
            self.assess_centrality(&mut r, &xl, &xu, &zl, &zu, mu_now, true);
            let mut prev_ratio = self.centring_ratio;
            let mut prev_bad_products = self.bad_products;

            r.info.centring_success = 0;
            // if ratio is below tolerance, point is centred
            if prev_ratio < self.control_.centring_ratio_tolerance() {
                self.control_.h_log_str("\tPoint is now centred\n");
                r.info.centring_success = 1;
            } else {
                // perform centring steps
                let mut centring_complete = false;
                let max_steps = self.control_.max_centring_steps();
                let mut ii: Int = 0;
                while ii < max_steps {
                    // compute centring step
                    self.centring(&mut r, &mut step, mu_frozen);

                    // assess whether to take the step
                    let accept =
                        self.evaluate_centring_step(&mut r, &step, prev_ratio, prev_bad_products);
                    if !accept {
                        self.control_
                            .h_log_str("\tPoint cannot be centred further\n");
                        centring_complete = true;
                        break;
                    }

                    // take the step and print output
                    self.make_step(&mut r, &step, true);
                    r.info.iter += 1;
                    self.print_output(&mut r);
                    oracle_trace::ipx_iter(
                        "centring",
                        r.iterate,
                        r.info,
                        &*r.kkt,
                        self.step_primal_,
                        self.step_dual_,
                    );
                    let (xl, xu, zl, zu) = (
                        r.iterate.xl().clone(),
                        r.iterate.xu().clone(),
                        r.iterate.zl().clone(),
                        r.iterate.zu().clone(),
                    );
                    let mu_now = r.iterate.mu();
                    self.assess_centrality(&mut r, &xl, &xu, &zl, &zu, mu_now, true);
                    prev_ratio = self.centring_ratio;
                    prev_bad_products = self.bad_products;

                    // if ratio is below tolerance, point is centred
                    if prev_ratio < self.control_.centring_ratio_tolerance() {
                        self.control_.h_log_str("\tPoint is now centred\n");
                        r.info.centring_success = 1;
                        centring_complete = true;
                        break;
                    }
                    ii += 1;
                }
                if !centring_complete {
                    let mut s = String::new();
                    s += &format!(
                        "\tPoint could not be centred within {} iterations\n",
                        max_steps
                    );
                    self.control_.h_log(&mut s);
                }
            }
            r.info.centring_tried = 1;
        }
        oracle_trace::ipx_info("driver_end", r.info);
    }

    /// HiGHS: ipx/ipm.cc::IPM::ComputeStartingPoint (lines 194-328)
    pub(super) fn compute_starting_point(&mut self, r: &mut IpmRefs<'_>) {
        let model_rc = r.iterate.model_rc();
        let model: &Model = &model_rc;
        let m = model.rows();
        let n = model.cols();
        let nm = n + m;
        let ai = model.AI();
        let b = model.b();
        let c = model.c();
        let lb = model.lb();
        let ub = model.ub();
        let mut x = Vector::new(nm);
        let mut xl = Vector::new(nm);
        let mut xu = Vector::new(nm);
        let mut y = Vector::new(m);
        let mut zl = Vector::new(nm);
        let mut zu = Vector::new(nm);
        // Factorize the KKT matrix with the identity matrix in the (1,1) block.
        r.kkt.factorize(None, &mut *r.info);
        if r.info.errflag != 0 {
            return;
        }

        // Set x within its bounds and compute the minimum norm solution dx to
        // AI*dx = (b-AI*x). Then update x := x + dx to obtain a feasible point.
        let mut rb = b.clone(); // workspace
        for j in 0..nm {
            let ju = j as usize;
            let mut xj = 0.0;
            if xj < lb[ju] {
                xj = lb[ju];
            }
            if xj > ub[ju] {
                xj = ub[ju];
            }
            x[ju] = xj;
            if xj != 0.0 {
                scatter_column(ai, j, -xj, &mut rb);
            }
        }
        let tol = 0.1 * infnorm(&rb);
        zl.fill(0.0);
        r.kkt.solve(&zl, &rb, tol, &mut xl, &mut y, &mut *r.info);
        if r.info.errflag != 0 {
            return;
        }
        va_add_assign(&mut x, &xl);

        // Compute xl, xu and shift to become positive.
        let mut xinfeas = 0.0;
        for j in 0..nm as usize {
            xl[j] = x[j] - lb[j];
            xinfeas = cpp_max(xinfeas, -xl[j]);
            xu[j] = ub[j] - x[j];
            xinfeas = cpp_max(xinfeas, -xu[j]);
        }
        let xshift1 = 1.0 + 1.5 * xinfeas;
        xl.iter_mut().for_each(|v| *v += xshift1);
        xu.iter_mut().for_each(|v| *v += xshift1);

        let cnorm = twonorm(c);
        if cnorm == 0.0 {
            // Special treatment for zero objective.
            for j in 0..nm as usize {
                zl[j] = if lb[j].is_finite() { 1.0 } else { 0.0 };
                zu[j] = if ub[j].is_finite() { 1.0 } else { 0.0 };
            }
        } else {
            // Compute y as the least-squares solution to AI'*y=c.
            // Recompute zl = c-AI'*y because the KKT system is solved
            // approximately with a residual in the first block equation.
            rb.fill(0.0);
            let tol_c = 0.1 * infnorm(c);
            r.kkt.solve(c, &rb, tol_c, &mut zl, &mut y, &mut *r.info);
            if r.info.errflag != 0 {
                return;
            }
            zl = c.clone();
            multiply_add(ai, &y, -1.0, &mut zl, b'T');

            // When c lies in range(AI'), then the dual slack variables are (close
            // to) zero, and the initial point would be almost complementary but
            // usually not primal feasible. To prevent this from happening, add
            // a fraction of the objective to zl and adjust y. In exact computation
            // this does not affect dual feasibility.
            let znorm = twonorm(&zl);
            let rho = 0.05;
            if znorm < rho * cnorm {
                for j in 0..nm as usize {
                    zl[j] += rho * c[j];
                }
                y.iter_mut().for_each(|v| *v *= 1.0 - rho);
            }

            // Split dual slack solution into zl, zu and shift to become positive.
            let mut zinfeas = 0.0;
            for j in 0..nm as usize {
                let zval = zl[j];
                zl[j] = 0.0;
                zu[j] = 0.0;
                if lb[j].is_finite() && ub[j].is_finite() {
                    zl[j] = 0.5 * zval;
                    zu[j] = -0.5 * zval;
                } else if lb[j].is_finite() {
                    zl[j] = zval;
                } else if ub[j].is_finite() {
                    zu[j] = -zval;
                }
                zinfeas = cpp_max(zinfeas, -zl[j]);
                zinfeas = cpp_max(zinfeas, -zu[j]);
            }
            let zshift1 = 1.0 + 1.5 * zinfeas;
            for j in 0..nm as usize {
                if lb[j].is_finite() {
                    zl[j] += zshift1;
                }
                if ub[j].is_finite() {
                    zu[j] += zshift1;
                }
            }
        }

        // Level pairwise complementarity products.
        let mut xsum = 1.0;
        let mut zsum = 1.0;
        let mut mu = 1.0;
        for j in 0..nm as usize {
            if lb[j].is_finite() {
                xsum += xl[j];
                zsum += zl[j];
                mu += xl[j] * zl[j];
            }
            if ub[j].is_finite() {
                xsum += xu[j];
                zsum += zu[j];
                mu += xu[j] * zu[j];
            }
        }
        let xshift2 = 0.5 * mu / zsum;
        let zshift2 = 0.5 * mu / xsum;
        xl.iter_mut().for_each(|v| *v += xshift2);
        xu.iter_mut().for_each(|v| *v += xshift2);
        for j in 0..nm as usize {
            if lb[j].is_finite() {
                zl[j] += zshift2;
            }
            if ub[j].is_finite() {
                zu[j] += zshift2;
            }
        }
        r.iterate.initialize(&x, &xl, &xu, &y, &zl, &zu);
        self.best_complementarity_ = r.iterate.complementarity();
    }

    /// HiGHS: ipx/ipm.cc::IPM::Predictor (lines 353-383)
    pub(super) fn predictor(&mut self, r: &mut IpmRefs<'_>, step: &mut Step) {
        let model_rc = r.iterate.model_rc();
        let n = model_rc.cols();
        let m = model_rc.rows();
        let nm = n + m;
        let mut sl = Vector::new(nm);
        let mut su = Vector::new(nm);
        {
            let it: &Iterate = &*r.iterate;
            let (xl, zl) = (it.xl(), it.zl());
            // sl = -xl.*zl
            for j in 0..nm {
                let ju = j as usize;
                sl[ju] = if it.has_barrier_lb(j) {
                    -xl[ju] * zl[ju]
                } else {
                    0.0
                };
            }
            debug_assert!(all_finite(&sl));
            // su = -xu.*zu
            let (xu, zu) = (it.xu(), it.zu());
            for j in 0..nm {
                let ju = j as usize;
                su[ju] = if it.has_barrier_ub(j) {
                    -xu[ju] * zu[ju]
                } else {
                    0.0
                };
            }
            debug_assert!(all_finite(&su));
        }
        // Clone the residuals: the Ref guards must not outlive the solve.
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

    /// HiGHS: ipx/ipm.cc::IPM::AddCorrector (lines 385-448)
    pub(super) fn add_corrector(&mut self, r: &mut IpmRefs<'_>, step: &mut Step) {
        let model_rc = r.iterate.model_rc();
        let n = model_rc.cols();
        let m = model_rc.rows();
        let nm = n + m;
        let mu = r.iterate.mu();
        let mut sl = Vector::new(nm);
        let mut su = Vector::new(nm);
        {
            let it: &Iterate = &*r.iterate;
            let (xl, xu, zl, zu) = (it.xl(), it.xu(), it.zl(), it.zu());
            let (dxl, dxu, dzl, dzu) = (&step.xl, &step.xu, &step.zl, &step.zu);

            // Choose centering parameter.
            let step_xl = step_to_boundary(xl, dxl, None, 1.0);
            let step_xu = step_to_boundary(xu, dxu, None, 1.0);
            let step_zl = step_to_boundary(zl, dzl, None, 1.0);
            let step_zu = step_to_boundary(zu, dzu, None, 1.0);
            let maxp = cpp_min(step_xl, step_xu);
            let maxd = cpp_min(step_zl, step_zu);
            let mut muaff = 0.0;
            let mut num_finite: Int = 0;
            for j in 0..nm {
                let ju = j as usize;
                if it.has_barrier_lb(j) {
                    debug_assert!(xl[ju].is_finite());
                    debug_assert!(xl[ju] != 0.0);
                    muaff += (xl[ju] + maxp * dxl[ju]) * (zl[ju] + maxd * dzl[ju]);
                    num_finite += 1;
                }
                if it.has_barrier_ub(j) {
                    debug_assert!(xu[ju].is_finite());
                    debug_assert!(xu[ju] != 0.0);
                    muaff += (xu[ju] + maxp * dxu[ju]) * (zu[ju] + maxd * dzu[ju]);
                    num_finite += 1;
                }
            }
            debug_assert!(muaff.is_finite());
            muaff /= num_finite as f64;
            let ratio = muaff / mu;
            let sigma = ratio * ratio * ratio;

            // sl = -xl.*zl + sigma*mu - dxl.*dzl
            for j in 0..nm {
                let ju = j as usize;
                sl[ju] = if it.has_barrier_lb(j) {
                    -xl[ju] * zl[ju] + sigma * mu - dxl[ju] * dzl[ju]
                } else {
                    0.0
                };
            }
            debug_assert!(all_finite(&sl));

            // su = -xu.*zu + sigma*mu - dxu.*dzu
            for j in 0..nm {
                let ju = j as usize;
                su[ju] = if it.has_barrier_ub(j) {
                    -xu[ju] * zu[ju] + sigma * mu - dxu[ju] * dzu[ju]
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
}

/// Computes maximum alpha such that x + alpha*dx >= 0.
/// The blocking index is returned in blocking_index if not NULL.
///
/// HiGHS: ipx/ipm.cc::StepToBoundary (static, lines 332-351)
///
/// C++ default arguments: `alpha = 1.0`. `blocking_index` may be `nullptr`.
pub(super) fn step_to_boundary(
    x: &Vector,
    dx: &Vector,
    blocking_index: Option<&mut Int>,
    alpha: f64,
) -> f64 {
    let n = x.size();
    let damp = 1.0 - f64::EPSILON;
    debug_assert!(damp < 1.0);
    let mut alpha = alpha;

    let mut iblock: Int = -1;
    for i in 0..n as usize {
        debug_assert!(x[i] >= 0.0);
        if x[i] + alpha * dx[i] < 0.0 {
            alpha = -(x[i] * damp) / dx[i];
            debug_assert!(x[i] + alpha * dx[i] >= 0.0);
            iblock = i as Int;
        }
    }
    debug_assert!(alpha >= 0.0);
    if let Some(bi) = blocking_index {
        *bi = iblock;
    }
    alpha
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// Hand-computed: x0 + 1*dx0 = -1 < 0 blocks at index 0 with
    /// alpha = -(x0*damp)/dx0 = damp/2, where damp = 1 - eps. Index 1 does
    /// not block (2 + 1 >= 0).
    #[test]
    fn step_to_boundary_damped_block() {
        let x = Vector::from_slice(&[1.0, 2.0]);
        let dx = Vector::from_slice(&[-2.0, 1.0]);
        let mut block: Int = -7;
        let alpha = step_to_boundary(&x, &dx, Some(&mut block), 1.0);
        assert_eq!(block, 0);
        assert_eq!(alpha, (1.0 - f64::EPSILON) / 2.0);
    }

    /// No blocking: the default alpha is returned and the index is -1.
    #[test]
    fn step_to_boundary_no_block() {
        let x = Vector::from_slice(&[1.0, 2.0]);
        let dx = Vector::from_slice(&[0.5, 1.0]);
        let mut block: Int = 9;
        let alpha = step_to_boundary(&x, &dx, Some(&mut block), 1.0);
        assert_eq!(block, -1);
        assert_eq!(alpha, 1.0);
    }

    /// Builds a minimisation LP from row-wise triplets. Each row is
    /// `(lower, upper, [(col, value)])`.
    fn h16_ipx_scaled_lp(
        cost: &[f64],
        col_bounds: &[(f64, f64)],
        rows: &[(f64, f64, &[(usize, f64)])],
    ) -> HighsLp {
        let num_col = cost.len();
        let mut start: Vec<HighsInt> = vec![0];
        let mut index: Vec<HighsInt> = Vec::new();
        let mut value: Vec<f64> = Vec::new();
        for j in 0..num_col {
            for (i, row) in rows.iter().enumerate() {
                for &(col, v) in row.2 {
                    if col == j {
                        index.push(i as HighsInt);
                        value.push(v);
                    }
                }
            }
            start.push(index.len() as HighsInt);
        }
        let mut lp = HighsLp {
            num_col_: num_col as HighsInt,
            num_row_: rows.len() as HighsInt,
            col_cost_: cost.to_vec(),
            col_lower_: col_bounds.iter().map(|b| b.0).collect(),
            col_upper_: col_bounds.iter().map(|b| b.1).collect(),
            row_lower_: rows.iter().map(|r| r.0).collect(),
            row_upper_: rows.iter().map(|r| r.1).collect(),
            sense_: ObjSense::kMinimize,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = num_col as HighsInt;
        lp.a_matrix_.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
        lp
    }

    /// check/instances/3512-feasible.lp (TestIpm.cpp test-3512, objective -24).
    fn h16_ipx_scaled_feasible() -> HighsLp {
        let cost = [
            322898756.0,
            -341262841.0,
            485421114.0,
            131623047.0,
            -598680075.0,
            -545725019.0,
            -52955059.0,
            726068482.0,
            -403169727.0,
            -341262841.0,
            526021831.0,
            116115101.0,
            -25092767.0,
            -259011453.0,
            -414101970.0,
            454702685.0,
            -481061453.0,
            -44960378.0,
            744432569.0,
        ];
        let mut bounds = Vec::new();
        bounds.extend([(-1e9, -999999989.0); 5]);
        bounds.extend([(999999999.0, 1e9); 8]);
        bounds.extend([(-1e9, -999999990.0); 6]);
        let rows: [(f64, f64, &[(usize, f64)]); 10] = [
            (0.0, 0.0, &[(7, -1.0), (8, 1.0), (13, 1.0), (18, -1.0)]),
            (0.0, 0.0, &[(16, -1.0), (17, 1.0)]),
            (0.0, 0.0, &[(0, -1.0), (1, 1.0), (8, -1.0), (9, 1.0)]),
            (
                0.0,
                0.0,
                &[
                    (1, -1.0),
                    (2, 1.0),
                    (9, -1.0),
                    (10, 1.0),
                    (17, -1.0),
                    (18, 1.0),
                ],
            ),
            (0.0, 0.0, &[(10, -1.0), (11, 1.0), (15, -1.0), (16, 1.0)]),
            (0.0, 0.0, &[(3, -1.0), (4, 1.0), (5, 1.0), (12, -1.0)]),
            (0.0, 0.0, &[(2, -1.0), (3, 1.0), (13, -1.0), (14, 1.0)]),
            (0.0, 0.0, &[(5, -1.0), (6, 1.0), (14, -1.0), (15, 1.0)]),
            (0.0, 0.0, &[(0, 1.0), (4, -1.0), (6, -1.0), (7, 1.0)]),
            (0.0, 0.0, &[(11, -1.0), (12, 1.0)]),
        ];
        h16_ipx_scaled_lp(&cost, &bounds, &rows)
    }

    /// check/instances/3512-solve-error.lp (TestIpm.cpp test-3512, objective -499999991).
    fn h16_ipx_scaled_solve_error() -> HighsLp {
        let cost = [
            -555720215.0,
            0.0,
            -10.0,
            -16723359.0,
            -1.0,
            10.0,
            100.0,
            538996856.0,
            -10.0,
            16723358.0,
        ];
        let bounds = [
            (0.0, 1e9),
            (0.0, 5e8),
            (0.0, 1000.0),
            (0.0, 1e9),
            (0.0, 1000.0),
            (0.0, 500000.0),
            (0.0, 1e6),
            (0.0, 1e9),
            (0.0, 1e6),
            (0.0, 5e8),
        ];
        let rows: [(f64, f64, &[(usize, f64)]); 7] = [
            (0.0, 0.0, &[(0, -1.0), (3, -1.0), (9, 1.0)]),
            (-1.0, -1.0, &[(4, 1.0), (5, -1.0)]),
            (-1.0, -1.0, &[(1, 1.0), (4, -1.0), (6, -1.0)]),
            (0.0, 0.0, &[(3, 1.0), (7, 1.0), (9, -1.0)]),
            (1.0, 1.0, &[(1, -1.0), (5, 1.0), (8, -1.0)]),
            (1.0, 1.0, &[(2, 1.0), (6, 1.0), (8, 1.0)]),
            (0.0, 0.0, &[(0, 1.0), (2, -1.0), (7, -1.0)]),
        ];
        h16_ipx_scaled_lp(&cost, &bounds, &rows)
    }

    /// check/instances/3512-dual-cancellation.lp (TestIpm.cpp test-3512, objective -33).
    fn h16_ipx_scaled_dual_cancellation() -> HighsLp {
        let cost = [
            448490616.0,
            403849445.0,
            -852340062.0,
            -348413366.0,
            219745604.0,
            128667762.0,
            -469649735.0,
            517246712.0,
            -47596979.0,
            2.0,
            -51178168.0,
            176755522.0,
            787563265.0,
            -711518122.0,
            -201622496.0,
            -294845557.0,
            294845558.0,
            109003890.0,
            -337748902.0,
            128330511.0,
            -292894215.0,
            710415490.0,
            -317106778.0,
            710415492.0,
            -417521278.0,
            -148967297.0,
            358385685.0,
            -502312603.0,
        ];
        let mut bounds = Vec::new();
        bounds.extend([(999999998.0, 1e9); 6]);
        bounds.extend([(-1e9, -999999997.0); 3]);
        bounds.push((999999990.0, 1e9));
        bounds.extend([(999999998.0, 1e9); 5]);
        bounds.extend([(-1e9, -899999999.0); 2]);
        bounds.extend([(499999999.0, 500000001.0); 6]);
        bounds.extend([(999999999.0, 999999999.0); 5]);
        let rows: [(f64, f64, &[(usize, f64)]); 15] = [
            (0.0, 0.0, &[(25, -1.0), (26, 1.0)]),
            (
                0.0,
                0.0,
                &[
                    (6, 1.0),
                    (8, -1.0),
                    (19, -1.0),
                    (20, 1.0),
                    (24, -1.0),
                    (25, 1.0),
                ],
            ),
            (0.0, 0.0, &[(21, -1.0), (22, 1.0), (23, -1.0), (24, 1.0)]),
            (
                0.0,
                0.0,
                &[
                    (11, -1.0),
                    (12, 1.0),
                    (20, -1.0),
                    (21, 1.0),
                    (23, 1.0),
                    (27, -1.0),
                ],
            ),
            (0.0, 0.0, &[(1, -1.0), (2, 1.0), (15, 1.0), (16, -1.0)]),
            (0.0, 0.0, &[(0, -1.0), (1, 1.0), (17, 1.0), (22, -1.0)]),
            (0.0, 0.0, &[(6, -1.0), (7, 1.0), (10, -1.0), (11, 1.0)]),
            (0.0, 0.0, &[(0, 1.0), (2, -1.0), (3, -1.0), (4, 1.0)]),
            (0.0, 0.0, &[(13, -1.0), (14, 1.0)]),
            (0.0, 0.0, &[(4, -1.0), (5, 1.0), (18, -1.0), (19, 1.0)]),
            (0.0, 0.0, &[(12, -1.0), (13, 1.0)]),
            (0.0, 0.0, &[(7, -1.0), (8, 1.0)]),
            (0.0, 0.0, &[(3, 1.0), (5, -1.0)]),
            (0.0, 0.0, &[(10, 1.0), (14, -1.0)]),
            (
                0.0,
                0.0,
                &[
                    (15, -1.0),
                    (16, 1.0),
                    (17, -1.0),
                    (18, 1.0),
                    (26, -1.0),
                    (27, 1.0),
                ],
            ),
        ];
        h16_ipx_scaled_lp(&cost, &bounds, &rows)
    }

    /// Runs the LP through HiGHS with the IPX solver (as TestIpm.cpp does) and the
    /// HiGHS 1.16 #3513 option set to `scaled`. ipm_iteration_limit = 100 bounds the runs:
    /// 3512-feasible and 3512-dual-cancellation never terminate without a limit in
    /// C++ 1.15.1, nor in 1.15.1 + the #3513 IPX hunks (their 1.16.0 optimal result
    /// needs other 1.16 changes).
    fn h16_ipx_scaled_run(lp: HighsLp, scaled: bool) -> (HighsModelStatus, f64, HighsInt) {
        let mut highs = crate::lp_data::highs_run::Highs::default();
        assert_eq!(
            highs.set_option_value_bool("log_to_console", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", "ipx"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_int("ipm_iteration_limit", 100),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_ipx_scaled", scaled),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        let _ = highs.run();
        (
            highs.get_model_status(),
            highs.get_info().objective_function_value,
            highs.get_info().ipm_iteration_count,
        )
    }

    /// HiGHS 1.16 #3513, flag on: check/instances/3512-solve-error.lp no longer ends in a
    /// solve error. The CR solve fails in double precision at IPM iteration 13 and is
    /// restarted with extended precision (kkt_solver_basis.cc lines 149-158).
    /// Expected values come from a C++ oracle (1.15.1 + the #3513 IPX hunks, built
    /// with -ffp-contract=off): IPM optimal after 17 iterations, crossover optimal,
    /// model status Optimal, objective -499999991 (the TestIpm value).
    #[test]
    #[ignore = "h16_ipx_scaled: slow"]
    fn h16_ipx_scaled_solve_error_on_optimal() {
        let (status, objective, ipm_iter) = h16_ipx_scaled_run(h16_ipx_scaled_solve_error(), true);
        assert_eq!(status, HighsModelStatus::kOptimal);
        assert_eq!(ipm_iter, 17);
        let error = (objective - -499999991.0).abs() / (1.0 + objective.abs());
        assert!(error < 1e-8, "objective {objective}");
    }

    /// HiGHS 1.16 #3513, flag on: check/instances/3512-feasible.lp and
    /// 3512-dual-cancellation.lp are not reported infeasible; as in the C++ oracle
    /// (1.15.1 + the #3513 IPX hunks), both stop at the IPM iteration limit.
    #[test]
    #[ignore = "h16_ipx_scaled: slow"]
    fn h16_ipx_scaled_feasible_on_iteration_limit() {
        let (feasible, _, _) = h16_ipx_scaled_run(h16_ipx_scaled_feasible(), true);
        assert_eq!(feasible, HighsModelStatus::kIterationLimit);
        let (dual, _, _) = h16_ipx_scaled_run(h16_ipx_scaled_dual_cancellation(), true);
        assert_eq!(dual, HighsModelStatus::kIterationLimit);
    }

    /// HiGHS 1.16 #3513, flag off: the 1.15.1 statuses are kept (observed with C++ 1.15.1
    /// at ipm_iteration_limit = 100): 3512-solve-error ends in a solve error, and
    /// 3512-feasible and 3512-dual-cancellation stop at the IPM iteration limit.
    /// (C++ 1.15.1 built with -ffp-contract=off.)
    #[test]
    #[ignore = "h16_ipx_scaled: slow"]
    fn h16_ipx_scaled_off_keeps_1_15_1_statuses() {
        let (solve_error, _, _) = h16_ipx_scaled_run(h16_ipx_scaled_solve_error(), false);
        assert_eq!(solve_error, HighsModelStatus::kSolveError);
        let (feasible, _, _) = h16_ipx_scaled_run(h16_ipx_scaled_feasible(), false);
        assert_eq!(feasible, HighsModelStatus::kIterationLimit);
        let (dual, _, _) = h16_ipx_scaled_run(h16_ipx_scaled_dual_cancellation(), false);
        assert_eq!(dual, HighsModelStatus::kIterationLimit);
    }
}
