//! HiGHS source: ipm/ipx/model.cc (lines 982-1526)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::sparse_matrix::{dot_column, scatter_column};
use super::*;
use crate::util::highs_utils::cpp_max;

// C-shaped index loops and long parameter lists are kept as in model.cc.
#[allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::explicit_counter_loop,
    clippy::manual_memcpy
)]
impl Model {
    /// HiGHS: ipx/model.cc::Model::ScalePoint(x, slack, y, z) (lines 982-995)
    pub(super) fn scale_point(
        &self,
        x: &mut Vector,
        slack: &mut Vector,
        y: &mut Vector,
        z: &mut Vector,
    ) {
        if self.colscale_.size() > 0 {
            div_assign(x, &self.colscale_);
            mul_assign(z, &self.colscale_);
        }
        if self.rowscale_.size() > 0 {
            div_assign(y, &self.rowscale_);
            mul_assign(slack, &self.rowscale_);
        }
        for &j in self.flipped_vars_.iter() {
            let j = j as usize;
            x[j] *= -1.0;
            z[j] *= -1.0;
        }
    }
    /// HiGHS: ipx/model.cc::Model::ScalePoint(x, xl, xu, slack, y, zl, zu) (lines 997-1017)
    pub(super) fn scale_point_bounds(
        &self,
        x: &mut Vector,
        xl: &mut Vector,
        xu: &mut Vector,
        slack: &mut Vector,
        y: &mut Vector,
        zl: &mut Vector,
        zu: &mut Vector,
    ) {
        if self.colscale_.size() > 0 {
            div_assign(x, &self.colscale_);
            div_assign(xl, &self.colscale_);
            div_assign(xu, &self.colscale_);
            mul_assign(zl, &self.colscale_);
            mul_assign(zu, &self.colscale_);
        }
        if self.rowscale_.size() > 0 {
            div_assign(y, &self.rowscale_);
            mul_assign(slack, &self.rowscale_);
        }
        for &j in self.flipped_vars_.iter() {
            let j = j as usize;
            x[j] *= -1.0;
            xl[j] = xu[j];
            xu[j] = f64::INFINITY;
            zl[j] = zu[j];
            zu[j] = 0.0;
        }
    }
    /// HiGHS: ipx/model.cc::Model::ScaleBackInteriorSolution (lines 1019-1043)
    pub(super) fn scale_back_interior_solution(
        &self,
        x: &mut Vector,
        xl: &mut Vector,
        xu: &mut Vector,
        slack: &mut Vector,
        y: &mut Vector,
        zl: &mut Vector,
        zu: &mut Vector,
    ) {
        if self.colscale_.size() > 0 {
            mul_assign(x, &self.colscale_);
            mul_assign(xl, &self.colscale_);
            mul_assign(xu, &self.colscale_);
            div_assign(zl, &self.colscale_);
            div_assign(zu, &self.colscale_);
        }
        if self.rowscale_.size() > 0 {
            mul_assign(y, &self.rowscale_);
            div_assign(slack, &self.rowscale_);
        }
        for &j in self.flipped_vars_.iter() {
            let j = j as usize;
            debug_assert!(xl[j].is_finite());
            debug_assert!(xu[j].is_infinite());
            debug_assert!(zu[j] == 0.0);
            x[j] *= -1.0;
            xu[j] = xl[j];
            xl[j] = f64::INFINITY;
            zu[j] = zl[j];
            zl[j] = 0.0;
        }
    }
    /// HiGHS: ipx/model.cc::Model::ScaleBackResiduals (lines 1045-1060)
    pub(super) fn scale_back_residuals(
        &self,
        rb: &mut Vector,
        rc: &mut Vector,
        rl: &mut Vector,
        ru: &mut Vector,
    ) {
        if self.colscale_.size() > 0 {
            div_assign(rc, &self.colscale_);
            mul_assign(rl, &self.colscale_);
            mul_assign(ru, &self.colscale_);
        }
        if self.rowscale_.size() > 0 {
            div_assign(rb, &self.rowscale_);
        }
        for &j in self.flipped_vars_.iter() {
            let j = j as usize;
            rc[j] *= -1.0;
            debug_assert!(ru[j] == 0.0);
            ru[j] = -rl[j];
            rl[j] = 0.0;
        }
    }
    /// HiGHS: ipx/model.cc::Model::ScaleBackBasicSolution (lines 1062-1076)
    pub(super) fn scale_back_basic_solution(
        &self,
        x: &mut Vector,
        slack: &mut Vector,
        y: &mut Vector,
        z: &mut Vector,
    ) {
        if self.colscale_.size() > 0 {
            mul_assign(x, &self.colscale_);
            div_assign(z, &self.colscale_);
        }
        if self.rowscale_.size() > 0 {
            mul_assign(y, &self.rowscale_);
            div_assign(slack, &self.rowscale_);
        }
        for &j in self.flipped_vars_.iter() {
            let j = j as usize;
            x[j] *= -1.0;
            z[j] *= -1.0;
        }
    }
    /// HiGHS: ipx/model.cc::Model::ScaleBackBasis (lines 1078-1085)
    pub(super) fn scale_back_basis(&self, cbasis: &mut [Int], vbasis: &mut [Int]) {
        // C++ parameter `cbasis` is unused here.
        let _ = cbasis;
        for &j in self.flipped_vars_.iter() {
            let j = j as usize;
            debug_assert!(vbasis[j] != IPX_nonbasic_ub);
            if vbasis[j] == IPX_nonbasic_lb {
                vbasis[j] = IPX_nonbasic_ub;
            }
        }
    }
    /// HiGHS: ipx/model.cc::Model::DualizeBasicSolution (lines 1087-1136)
    pub(super) fn dualize_basic_solution(
        &self,
        x_user: &Vector,
        slack_user: &Vector,
        y_user: &Vector,
        z_user: &Vector,
        x_solver: &mut Vector,
        y_solver: &mut Vector,
        z_solver: &mut Vector,
    ) {
        let m = self.rows();
        let n = self.cols();
        let nu = n as usize;
        let cu = self.num_constr_ as usize;
        let mu = m as usize;

        if self.dualized_ {
            debug_assert_eq!(self.num_var_, m);
            debug_assert_eq!(cu + self.boxed_vars_.len(), nu);

            // Build dual solver variables from primal user variables.
            y_solver.assign(&negated(&x_user[..]));
            for i in 0..cu {
                z_solver[i] = -slack_user[i];
            }
            for k in 0..self.boxed_vars_.len() {
                let j = self.boxed_vars_[k] as usize;
                z_solver[cu + k] = self.c_at(self.num_constr_ + k as Int) + y_solver[j];
            }
            for i in 0..mu {
                z_solver[nu + i] = self.c_at(n + i as Int) - y_solver[i];
            }

            // Build primal solver variables from dual user variables.
            copy_n(&y_user[..], self.num_constr_, &mut x_solver[..]);
            copy_n(&z_user[..], self.num_var_, &mut x_solver[nu..]);
            for k in 0..self.boxed_vars_.len() {
                let j = self.boxed_vars_[k] as usize;
                if x_solver[nu + j] < 0.0 {
                    // j is a boxed variable and z_user[j] < 0
                    x_solver[cu + k] = -x_solver[nu + j];
                    x_solver[nu + j] = 0.0;
                } else {
                    x_solver[cu + k] = 0.0;
                }
            }
        } else {
            debug_assert_eq!(self.num_constr_, m);
            debug_assert_eq!(self.num_var_, n);
            copy_n(&x_user[..], n, &mut x_solver[..]);
            copy_n(&slack_user[..], m, &mut x_solver[nu..]);
            copy_n(&y_user[..], m, &mut y_solver[..]);
            copy_n(&z_user[..], n, &mut z_solver[..]);
            for i in 0..mu {
                z_solver[nu + i] = self.c_at(n + i as Int) - y_solver[i];
            }
        }
    }
    /// HiGHS: ipx/model.cc::Model::DualizeIPMStartingPoint (lines 1138-1199)
    #[allow(clippy::assertions_on_constants)]
    pub(super) fn dualize_ipm_starting_point(
        &self,
        x_user: &Vector,
        xl_user: &Vector,
        xu_user: &Vector,
        slack_user: &Vector,
        y_user: &Vector,
        zl_user: &Vector,
        zu_user: &Vector,
        x_solver: &mut Vector,
        xl_solver: &mut Vector,
        xu_solver: &mut Vector,
        y_solver: &mut Vector,
        zl_solver: &mut Vector,
        zu_solver: &mut Vector,
    ) {
        let m = self.rows();
        let n = self.cols();
        let nu = n as usize;

        if self.dualized_ {
            // Not implemented.
            debug_assert!(
                false,
                "DualizeIPMStartingPoint: not implemented for dualized models"
            );
        } else {
            debug_assert_eq!(self.num_constr_, m);
            debug_assert_eq!(self.num_var_, n);

            copy_n(&x_user[..], self.num_var_, &mut x_solver[..]);
            copy_n(&slack_user[..], self.num_constr_, &mut x_solver[nu..]);
            copy_n(&xl_user[..], self.num_var_, &mut xl_solver[..]);
            copy_n(&xu_user[..], self.num_var_, &mut xu_solver[..]);
            copy_n(&y_user[..], self.num_constr_, &mut y_solver[..]);
            copy_n(&zl_user[..], self.num_var_, &mut zl_solver[..]);
            copy_n(&zu_user[..], self.num_var_, &mut zu_solver[..]);

            for i in 0..m as usize {
                match self.constr_type_[i] {
                    b'=' => {
                        debug_assert!(self.lb_[nu + i] == 0.0 && self.ub_[nu + i] == 0.0);
                        // For a fixed slack variable xl, xu, zl and zu won't be used
                        // by the IPM. Just put them to zero.
                        xl_solver[nu + i] = 0.0;
                        xu_solver[nu + i] = 0.0;
                        zl_solver[nu + i] = 0.0;
                        zu_solver[nu + i] = 0.0;
                    }
                    b'<' => {
                        debug_assert!(self.lb_[nu + i] == 0.0 && self.ub_[nu + i] == f64::INFINITY);
                        xl_solver[nu + i] = slack_user[i];
                        xu_solver[nu + i] = f64::INFINITY;
                        zl_solver[nu + i] = -y_user[i];
                        zu_solver[nu + i] = 0.0;
                    }
                    b'>' => {
                        debug_assert!(
                            self.lb_[nu + i] == f64::NEG_INFINITY && self.ub_[nu + i] == 0.0
                        );
                        xl_solver[nu + i] = f64::INFINITY;
                        xu_solver[nu + i] = -slack_user[i];
                        zl_solver[nu + i] = 0.0;
                        zu_solver[nu + i] = y_user[i];
                    }
                    _ => {}
                }
            }
        }
    }
    /// HiGHS: ipx/model.cc::Model::DualizeBackInteriorSolution (lines 1201-1344)
    pub(super) fn dualize_back_interior_solution(
        &self,
        x_solver: &Vector,
        xl_solver: &Vector,
        xu_solver: &Vector,
        y_solver: &Vector,
        zl_solver: &Vector,
        zu_solver: &Vector,
        x_user: &mut Vector,
        xl_user: &mut Vector,
        xu_user: &mut Vector,
        slack_user: &mut Vector,
        y_user: &mut Vector,
        zl_user: &mut Vector,
        zu_user: &mut Vector,
    ) {
        let m = self.rows();
        let n = self.cols();
        let nu = n as usize;
        let cu = self.num_constr_ as usize;
        let mu = m as usize;

        if self.dualized_ {
            debug_assert_eq!(self.num_var_, m);
            debug_assert_eq!(cu + self.boxed_vars_.len(), nu);
            x_user.assign(&negated(&y_solver[..]));

            // If the solution from the solver would be exact, we could copy the
            // first num_constr_ entries from x_solver into y_user. However, to
            // satisfy the sign condition on y_user even if the solution is not
            // exact, we have to use the xl_solver and xu_solver entries for
            // inequality constraints.
            for i in 0..cu {
                match self.constr_type_[i] {
                    b'=' => y_user[i] = x_solver[i],
                    b'<' => y_user[i] = -xu_solver[i],
                    b'>' => y_user[i] = xl_solver[i],
                    _ => {}
                }
                debug_assert!(y_user[i].is_finite());
            }

            // Dual variables associated with lbuser <= x in the scaled user model
            // are the slack variables from the solver. For an exact solution we
            // would have x_solver[n+1:n+m] == xl_solver[n+1:n+m]. Using xl_solver
            // guarantees that zl_user >= 0 in any case. If variable j has no lower
            // bound in the scaled user model (i.e. is free), then the j-th slack
            // variable was fixed at zero in the solver model, but the IPM solution
            // may not satisfy this. Hence we must set zl_user[j] = 0 explicitly.
            copy_n(&xl_solver[nu..], self.num_var_, &mut zl_user[..]);
            for j in 0..self.num_var_ as usize {
                if !self.scaled_lbuser_[j].is_finite() {
                    zl_user[j] = 0.0;
                }
            }

            // Dual variables associated with x <= ubuser in the scaled user model
            // are the primal variables that were added for boxed variables in the
            // solver model.
            zu_user.fill(0.0);
            let mut k = self.num_constr_;
            for &j in self.boxed_vars_.iter() {
                zu_user[j as usize] = xl_solver[k as usize];
                k += 1;
            }
            debug_assert_eq!(k, n);

            // xl in the scaled user model is zl[n+1:n+m] in the solver model or
            // infinity.
            for i in 0..mu {
                if self.scaled_lbuser_[i].is_finite() {
                    xl_user[i] = zl_solver[nu + i];
                } else {
                    xl_user[i] = f64::INFINITY;
                }
            }

            // xu in the scaled user model are the entries in zl for columns of the
            // negative identity matrix (that were added for boxed variables).
            xu_user.fill(f64::INFINITY);
            k = self.num_constr_;
            for &j in self.boxed_vars_.iter() {
                xu_user[j as usize] = zl_solver[k as usize];
                k += 1;
            }
            debug_assert_eq!(k, n);

            for i in 0..cu {
                match self.constr_type_[i] {
                    b'=' => slack_user[i] = 0.0,
                    b'<' => slack_user[i] = zu_solver[i],
                    b'>' => slack_user[i] = -zl_solver[i],
                    _ => {}
                }
            }
        } else {
            debug_assert_eq!(self.num_constr_, m);
            debug_assert_eq!(self.num_var_, n);
            copy_n(&x_solver[..], self.num_var_, &mut x_user[..]);

            // Instead of copying y_solver into y_user, we use the entries from
            // zl_solver and zu_solver for inequality constraints, so that the sign
            // condition on y_user is satisfied.
            for i in 0..mu {
                debug_assert!(self.lb_[nu + i] == 0.0 || self.lb_[nu + i] == f64::NEG_INFINITY);
                debug_assert!(self.ub_[nu + i] == 0.0 || self.ub_[nu + i] == f64::INFINITY);
                debug_assert!(self.lb_[nu + i] == 0.0 || self.ub_[nu + i] == 0.0);
                match self.constr_type_[i] {
                    b'=' => y_user[i] = y_solver[i],
                    b'<' => y_user[i] = -zl_solver[nu + i],
                    b'>' => y_user[i] = zu_solver[nu + i],
                    _ => {}
                }
                debug_assert!(y_user[i].is_finite());
            }
            copy_n(&zl_solver[..], self.num_var_, &mut zl_user[..]);
            copy_n(&zu_solver[..], self.num_var_, &mut zu_user[..]);
            copy_n(&xl_solver[..], self.num_var_, &mut xl_user[..]);
            copy_n(&xu_solver[..], self.num_var_, &mut xu_user[..]);

            // If the solution would be exact, slack_user were given by the entries
            // of x_solver corresponding to slack columns. To satisfy the sign
            // condition in any case, we build the slack for inequality constraints
            // from xl_solver and xu_solver and set the slack for equality
            // constraints to zero.
            for i in 0..mu {
                match self.constr_type_[i] {
                    b'=' => slack_user[i] = 0.0,
                    b'<' => slack_user[i] = xl_solver[nu + i],
                    b'>' => slack_user[i] = -xu_solver[nu + i],
                    _ => {}
                }
                debug_assert!(slack_user[i].is_finite());
            }
        }
    }
    /// HiGHS: ipx/model.cc::Model::DualizeBackBasicSolution (lines 1346-1378)
    pub(super) fn dualize_back_basic_solution(
        &self,
        x_solver: &Vector,
        y_solver: &Vector,
        z_solver: &Vector,
        x_user: &mut Vector,
        slack_user: &mut Vector,
        y_user: &mut Vector,
        z_user: &mut Vector,
    ) {
        let m = self.rows();
        let n = self.cols();
        let nu = n as usize;
        let cu = self.num_constr_ as usize;

        if self.dualized_ {
            debug_assert_eq!(self.num_var_, m);
            debug_assert_eq!(cu + self.boxed_vars_.len(), nu);
            x_user.assign(&negated(&y_solver[..]));
            for i in 0..cu {
                slack_user[i] = -z_solver[i];
            }
            copy_n(&x_solver[..], self.num_constr_, &mut y_user[..]);
            copy_n(&x_solver[nu..], self.num_var_, &mut z_user[..]);
            let mut k = self.num_constr_;
            for &j in self.boxed_vars_.iter() {
                z_user[j as usize] -= x_solver[k as usize];
                k += 1;
            }
            debug_assert_eq!(k, n);
        } else {
            debug_assert_eq!(self.num_constr_, m);
            debug_assert_eq!(self.num_var_, n);
            copy_n(&x_solver[..], self.num_var_, &mut x_user[..]);
            copy_n(&x_solver[nu..], self.num_constr_, &mut slack_user[..]);
            copy_n(&y_solver[..], self.num_constr_, &mut y_user[..]);
            copy_n(&z_solver[..], self.num_var_, &mut z_user[..]);
        }
    }
    /// HiGHS: ipx/model.cc::Model::DualizeBackBasis (lines 1380-1425)
    pub(super) fn dualize_back_basis(
        &self,
        basic_status_solver: &[Int],
        cbasis_user: &mut [Int],
        vbasis_user: &mut [Int],
    ) {
        let m = self.rows();
        let n = self.cols();
        let nu = n as usize;
        let cu = self.num_constr_ as usize;

        if self.dualized_ {
            debug_assert_eq!(self.num_var_, m);
            debug_assert_eq!(cu + self.boxed_vars_.len(), nu);
            for i in 0..cu {
                if basic_status_solver[i] == IPX_basic {
                    cbasis_user[i] = IPX_nonbasic;
                } else {
                    cbasis_user[i] = IPX_basic;
                }
            }
            for j in 0..self.num_var_ as usize {
                // slack cannot be superbasic
                debug_assert!(basic_status_solver[nu + j] != IPX_superbasic);
                // IPX quirk: compares against the literal 0 (equal to IPX_basic).
                if basic_status_solver[nu + j] == 0 {
                    vbasis_user[j] = if self.scaled_lbuser_[j].is_finite() {
                        IPX_nonbasic_lb
                    } else {
                        IPX_superbasic
                    };
                } else {
                    vbasis_user[j] = IPX_basic;
                }
            }
            let mut k = self.num_constr_;
            for &j in self.boxed_vars_.iter() {
                let status = basic_status_solver[k as usize];
                k += 1;
                if status == IPX_basic {
                    debug_assert!(vbasis_user[j as usize] == IPX_basic);
                    vbasis_user[j as usize] = IPX_nonbasic_ub;
                }
            }
        } else {
            debug_assert_eq!(self.num_constr_, m);
            debug_assert_eq!(self.num_var_, n);
            for i in 0..cu {
                // slack cannot be superbasic
                debug_assert!(basic_status_solver[nu + i] != IPX_superbasic);
                if basic_status_solver[nu + i] == IPX_basic {
                    cbasis_user[i] = IPX_basic;
                } else {
                    cbasis_user[i] = IPX_nonbasic;
                }
            }
            for j in 0..self.num_var_ as usize {
                vbasis_user[j] = basic_status_solver[j];
            }
        }
    }
    /// HiGHS: ipx/model.cc::Model::CorrectScaledBasicSolution (lines 1427-1445)
    ///
    /// C++ takes `cbasis`/`vbasis` by value (copies); they are only read.
    pub(super) fn correct_scaled_basic_solution(
        &self,
        x: &mut Vector,
        slack: &mut Vector,
        y: &mut Vector,
        z: &mut Vector,
        cbasis: &[Int],
        vbasis: &[Int],
    ) {
        for j in 0..self.num_var_ as usize {
            if vbasis[j] == IPX_nonbasic_lb {
                x[j] = self.scaled_lbuser_[j];
            }
            if vbasis[j] == IPX_nonbasic_ub {
                x[j] = self.scaled_ubuser_[j];
            }
            if vbasis[j] == IPX_basic {
                z[j] = 0.0;
            }
        }
        for i in 0..self.num_constr_ as usize {
            if cbasis[i] == IPX_nonbasic {
                slack[i] = 0.0;
            }
            if cbasis[i] == IPX_basic {
                y[i] = 0.0;
            }
        }
    }
    /// HiGHS: ipx/model.cc::Model::MultiplyWithScaledMatrix (lines 1447-1469)
    ///
    /// `trans`: `b't'` or `b'T'` for multiplication with A'.
    pub(super) fn multiply_with_scaled_matrix(
        &self,
        rhs: &Vector,
        alpha: f64,
        lhs: &mut Vector,
        trans: u8,
    ) {
        if trans == b't' || trans == b'T' {
            debug_assert_eq!(rhs.size(), self.num_constr_);
            debug_assert_eq!(lhs.size(), self.num_var_);
            if self.dualized() {
                for i in 0..self.num_constr_ {
                    scatter_column(&self.AI_, i, alpha * rhs[i as usize], &mut lhs[..]);
                }
            } else {
                for j in 0..self.num_var_ as usize {
                    lhs[j] += alpha * dot_column(&self.AI_, j as Int, &rhs[..]);
                }
            }
        } else {
            debug_assert_eq!(rhs.size(), self.num_var_);
            debug_assert_eq!(lhs.size(), self.num_constr_);
            if self.dualized() {
                for i in 0..self.num_constr_ as usize {
                    lhs[i] += alpha * dot_column(&self.AI_, i as Int, &rhs[..]);
                }
            } else {
                for j in 0..self.num_var_ {
                    scatter_column(&self.AI_, j, alpha * rhs[j as usize], &mut lhs[..]);
                }
            }
        }
    }
}

/// HiGHS: ipx/model.cc::PrimalInfeasibility (lines 1471-1482)
#[allow(clippy::needless_range_loop)]
pub fn primal_infeasibility(model: &Model, x: &Vector) -> f64 {
    let lb = model.lb();
    let ub = model.ub();
    debug_assert_eq!(x.size(), lb.size());

    let mut infeas = 0.0;
    for j in 0..x.len() {
        infeas = cpp_max(infeas, lb[j] - x[j]);
        infeas = cpp_max(infeas, x[j] - ub[j]);
    }
    infeas
}

/// HiGHS: ipx/model.cc::DualInfeasibility (lines 1484-1499)
#[allow(clippy::needless_range_loop)]
pub fn dual_infeasibility(model: &Model, x: &Vector, z: &Vector) -> f64 {
    let lb = model.lb();
    let ub = model.ub();
    debug_assert_eq!(x.size(), lb.size());
    debug_assert_eq!(z.size(), lb.size());

    let mut infeas = 0.0;
    for j in 0..x.len() {
        if x[j] > lb[j] {
            infeas = cpp_max(infeas, z[j]);
        }
        if x[j] < ub[j] {
            infeas = cpp_max(infeas, -z[j]);
        }
    }
    infeas
}

/// HiGHS: ipx/model.cc::PrimalResidual (lines 1501-1512)
#[allow(clippy::needless_range_loop)]
pub fn primal_residual(model: &Model, x: &Vector) -> f64 {
    let ait = model.AIt();
    let bvec = model.b();
    debug_assert_eq!(x.size(), ait.rows());

    let mut res = 0.0;
    for i in 0..bvec.len() {
        let r = bvec[i] - dot_column(ait, i as Int, &x[..]);
        res = cpp_max(res, r.abs());
    }
    res
}

/// HiGHS: ipx/model.cc::DualResidual (lines 1514-1526)
#[allow(clippy::needless_range_loop)]
pub fn dual_residual(model: &Model, y: &Vector, z: &Vector) -> f64 {
    let ai = model.AI();
    let cvec = model.c();
    debug_assert_eq!(y.size(), ai.rows());
    debug_assert_eq!(z.size(), ai.cols());

    let mut res = 0.0;
    for j in 0..cvec.len() {
        let r = cvec[j] - z[j] - dot_column(ai, j as Int, &y[..]);
        res = cpp_max(res, r.abs());
    }
    res
}
