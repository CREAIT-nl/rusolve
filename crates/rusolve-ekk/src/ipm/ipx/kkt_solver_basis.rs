//! HiGHS source: ipm/ipx/kkt_solver_basis.h, kkt_solver_basis.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! C++ `Basis& basis_` is `&'a mut Basis` (`LpSolver::RunMainIPM` builds the
//! solver on the stack from `*basis_`).
//!
//! Aliased `Basis::SolveDense(v, v, trans)` calls (C++ `basis_.SolveDense(work,
//! work, 'T')`) pass a copy of the right-hand side: `lu_solve_dense` copies `rhs`
//! into its workspace before it writes `lhs` (basiclu lu_solve_dense.c), so the
//! copy gives the same result as the C++ alias.

use super::basis::{BASIC, BASIC_FREE, Basis, NONBASIC};
use super::conjugate_residuals::ConjugateResiduals;
use super::control::{Control, sci2};
use super::indexed_vector::{IndexedVector, for_each_nonzero};
use super::iterate::Iterate;
use super::kkt_solver::KKTSolver;
use super::maxvolume::Maxvolume;
use super::model::Model;
use super::sparse_matrix::{dot_column, scatter_column};
use super::splitted_normal_matrix::SplittedNormalMatrix;
use super::*;
use std::rc::Rc;

/// HiGHS: ipx/kkt_solver_basis.h::KKTSolverBasis
pub struct KKTSolverBasis<'a> {
    control_: Rc<Control>,
    model_: Rc<Model>,
    basis_: &'a mut Basis,
    splitted_normal_matrix_: SplittedNormalMatrix,
    colscale_: Vector,
    factorized_: bool,
    maxiter_: Int,
    iter_: Int,
    basis_changes_: Int,
}

impl<'a> KKTSolverBasis<'a> {
    pub const kPivotZeroTol: f64 = 1e-7;

    /// HiGHS: ipx/kkt_solver_basis.cc::KKTSolverBasis::KKTSolverBasis (lines 10-16)
    pub fn new(control: Rc<Control>, basis: &'a mut Basis) -> Self {
        let model = basis.model_rc();
        let m = model.rows();
        let n = model.cols();
        let splitted_normal_matrix_ = SplittedNormalMatrix::new(Rc::clone(&model));
        KKTSolverBasis {
            control_: control,
            model_: model,
            basis_: basis,
            splitted_normal_matrix_,
            colscale_: Vector::new(n + m),
            factorized_: false,
            maxiter_: -1,
            iter_: 0,
            basis_changes_: 0,
        }
    }
    pub fn maxiter(&self) -> Int {
        self.maxiter_
    }
    /// HiGHS: ipx/kkt_solver_basis.h::KKTSolverBasis::maxiter(Int)
    pub fn set_maxiter(&mut self, new_maxiter: Int) {
        self.maxiter_ = new_maxiter;
    }
    /// HiGHS: ipx/kkt_solver_basis.cc::KKTSolverBasis::DropPrimal (lines 194-291)
    #[allow(clippy::needless_range_loop)]
    fn drop_primal(&mut self, iterate: &mut Iterate, info: &mut Info) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let mut btran = IndexedVector::new(m);
        let mut row = IndexedVector::new(n + m);
        let drop_primal = self.control_.ipm_drop_primal();
        let volume_tol = 2.0;
        info.errflag = 0;

        let mut candidates: Vec<Int> = Vec::new();
        for p in 0..m {
            let jb = self.basis_.at(p);
            if self.basis_.status_of(jb) != BASIC {
                continue; // ignore free variables
            }
            let (xl, xu) = (iterate.xl_at(jb), iterate.xu_at(jb));
            debug_assert!(xl.is_finite() || xu.is_finite());
            debug_assert!(xl > 0.0);
            debug_assert!(xu > 0.0);
            // choose which bound is nearer
            let (xj, zj) = if xl <= xu {
                (xl, iterate.zl_at(jb))
            } else {
                (xu, iterate.zu_at(jb))
            };
            if xj < 0.01 * zj && xj <= drop_primal {
                candidates.push(jb);
            }
        }
        if candidates.is_empty() {
            return;
        }

        // Maintain a copy of the inverse scaling factors of basic variables for
        // faster access.
        let mut invscale_basic = Vector::new(m);
        for p in 0..m {
            let j = self.basis_.at(p);
            invscale_basic[p as usize] = 1.0 / self.colscale_[j as usize];
            debug_assert!(invscale_basic[p as usize].is_finite());
            debug_assert!(invscale_basic[p as usize] >= 0.0);
        }

        while let Some(&jb) = candidates.last() {
            let p = self.basis_.position_of(jb);
            debug_assert!(p >= 0);
            let pu = p as usize;
            // Pivot jb out of the basis if the volume increase sufficiently.
            let s = invscale_basic[pu];
            self.basis_.tableau_row(jb, &mut btran, &mut row, true);
            let mut jmax: Int = -1;
            let mut vmax = volume_tol;
            {
                let colscale = &self.colscale_;
                let search_pivot = |j: Int, pivot: f64| {
                    let pivot = pivot.abs();
                    if pivot > Self::kPivotZeroTol {
                        let v = pivot * colscale[j as usize] * s;
                        if v > vmax {
                            vmax = v;
                            jmax = j;
                        }
                    }
                };
                for_each_nonzero(&row, search_pivot);
            }
            if jmax >= 0 {
                // Pivot jb out of the basis.
                let pivot = row.get(jmax);
                if pivot.abs() < 1e-3 && self.control_.debug_enabled(3) {
                    self.control_.debug(
                        3,
                        &format!(
                            " |pivot| = {} (primal basic variable close to bound)\n",
                            sci2(pivot.abs())
                        ),
                    );
                }
                debug_assert!(self.basis_.status_of(jmax) == NONBASIC);
                let mut exchanged = false;
                info.errflag = self
                    .basis_
                    .exchange_if_stable(jb, jmax, pivot, 1, &mut exchanged);
                if info.errflag != 0 {
                    return;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                invscale_basic[pu] = 1.0 / self.colscale_[jmax as usize];
                debug_assert!(invscale_basic[pu].is_finite());
                debug_assert!(invscale_basic[pu] >= 0.0);
                info.updates_ipm += 1;
                self.basis_changes_ += 1;
            } else {
                // Make variable jb "implied" at a bound.
                if iterate.zl_at(jb) / iterate.xl_at(jb) > iterate.zu_at(jb) / iterate.xu_at(jb) {
                    iterate.make_implied_lb(jb);
                } else {
                    iterate.make_implied_ub(jb);
                }
                self.basis_.free_basic_variable(jb);
                invscale_basic[pu] = 0.0;
                self.colscale_[jb as usize] = f64::INFINITY;
                info.primal_dropped += 1;
            }
            candidates.pop();
        }
    }
    /// HiGHS: ipx/kkt_solver_basis.cc::KKTSolverBasis::DropDual (lines 293-385)
    #[allow(clippy::needless_range_loop)]
    fn drop_dual(&mut self, iterate: &mut Iterate, info: &mut Info) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let mut ftran = IndexedVector::new(m);
        let drop_dual = self.control_.ipm_drop_dual();
        let soft = self.control_.ipm_drop_dual_soft();
        let volume_tol = 2.0;
        info.errflag = 0;

        let mut candidates: Vec<Int> = Vec::new();
        for jn in 0..n + m {
            if self.basis_.status_of(jn) != NONBASIC {
                continue;
            }
            let (xl, xu) = (iterate.xl_at(jn), iterate.xu_at(jn));
            let (zl, zu) = (iterate.zl_at(jn), iterate.zu_at(jn));
            debug_assert!(xl.is_finite() || xu.is_finite());
            debug_assert!(xl > 0.0);
            debug_assert!(xu > 0.0);
            debug_assert!(zl > 0.0 || zu > 0.0);
            // choose larger dual variable
            let (xj, zj) = if zl >= zu { (xl, zl) } else { (xu, zu) };
            if zj < 0.01 * xj && zj <= drop_dual {
                // rusolve lp_ipx_soft_drop: dropping zeroes the dual on both bounds, which
                // costs z*x of dual objective per finite bound
                if soft {
                    let mut comp = 0.0;
                    if xl.is_finite() {
                        comp += zl * xl;
                    }
                    if xu.is_finite() {
                        comp += zu * xu;
                    }
                    if comp > drop_dual {
                        continue;
                    }
                }
                candidates.push(jn);
            }
        }
        if candidates.is_empty() {
            return;
        }

        // Maintain a copy of the inverse scaling factors of basic variables for
        // faster access.
        let mut invscale_basic = Vector::new(m);
        for p in 0..m {
            let j = self.basis_.at(p);
            invscale_basic[p as usize] = 1.0 / self.colscale_[j as usize];
            debug_assert!(invscale_basic[p as usize].is_finite());
            debug_assert!(invscale_basic[p as usize] >= 0.0);
        }

        while let Some(&jn) = candidates.last() {
            // Pivot jn into the basis if volume increases sufficiently.
            let s = self.colscale_[jn as usize];
            self.basis_.solve_for_update_lhs(jn, &mut ftran);
            let mut pmax: Int = -1;
            let mut vmax = volume_tol;
            {
                let invscale = &invscale_basic;
                let search_pivot = |p: Int, pivot: f64| {
                    let pivot = pivot.abs();
                    if pivot > Self::kPivotZeroTol {
                        let v = pivot * invscale[p as usize] * s;
                        if v > vmax {
                            vmax = v;
                            pmax = p;
                        }
                    }
                };
                for_each_nonzero(&ftran, search_pivot);
            }
            if pmax >= 0 {
                let pivot = ftran.get(pmax);
                if pivot.abs() < 1e-3 && self.control_.debug_enabled(3) {
                    self.control_.debug(
                        3,
                        &format!(
                            " |pivot| = {} (dual nonbasic variable close to zero)\n",
                            sci2(pivot.abs())
                        ),
                    );
                }
                let jb = self.basis_.at(pmax);
                // Pivot jn into the basis.
                debug_assert!(self.basis_.status_of(jb) == BASIC);
                let mut exchanged = false;
                info.errflag = self
                    .basis_
                    .exchange_if_stable(jb, jn, pivot, -1, &mut exchanged);
                if info.errflag != 0 {
                    return;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                let pmu = pmax as usize;
                invscale_basic[pmu] = 1.0 / self.colscale_[jn as usize];
                debug_assert!(invscale_basic[pmu].is_finite());
                debug_assert!(invscale_basic[pmu] >= 0.0);
                info.updates_ipm += 1;
                self.basis_changes_ += 1;
            } else {
                // Make variable jn "fixed" at its current value.
                iterate.make_fixed(jn);
                self.basis_.fix_nonbasic_variable(jn);
                self.colscale_[jn as usize] = 0.0;
                info.dual_dropped += 1;
            }
            candidates.pop();
        }
    }
}

impl KKTSolver for KKTSolverBasis<'_> {
    /// HiGHS: ipx/kkt_solver_basis.cc::KKTSolverBasis::_Factorize (lines 18-65)
    fn _factorize(&mut self, iterate: Option<&mut Iterate>, info: &mut Info) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        info.errflag = 0;
        self.factorized_ = false;
        self.iter_ = 0;
        self.basis_changes_ = 0;

        // In the call to Factorize() @iterate must not be NULL.
        let iterate = iterate.expect("KKTSolverBasis::_Factorize: iterate must not be NULL");
        for j in 0..n + m {
            self.colscale_[j as usize] = iterate.scaling_factor(j);
        }

        // Remove degenerate variables unless the primal objective is smaller than
        // the dual objective. In the latter case the model might be infeasible or
        // unbounded. In order not to affect an infeasibility test in the IPM (not
        // yet implemented) we do not remove variables if the model looks
        // infeasible.
        if iterate.pobjective() >= iterate.dobjective() {
            self.drop_primal(iterate, info);
            if info.errflag != 0 {
                return;
            }
            self.drop_dual(iterate, info);
            if info.errflag != 0 {
                return;
            }
        }

        // Run maxvolume ("Russian algorithm").
        let mut maxvol = Maxvolume::new(Rc::clone(&self.control_));
        if self.control_.update_heuristic() == 0 {
            info.errflag = maxvol.run_sequential(&self.colscale_, &mut *self.basis_);
        } else {
            info.errflag = maxvol.run_heuristic(&self.colscale_, &mut *self.basis_);
        }
        info.updates_ipm += maxvol.updates();
        info.time_maxvol += maxvol.time();
        self.basis_changes_ += maxvol.updates();
        if info.errflag != 0 {
            return;
        }

        // Refactorize and build preconditioned normal matrix.
        if !self.basis_.factorization_is_fresh() {
            info.errflag = self.basis_.factorize();
            if info.errflag != 0 {
                return;
            }
        }
        self.splitted_normal_matrix_
            .prepare(&*self.basis_, &self.colscale_);

        self.factorized_ = true;
    }

    /// HiGHS: ipx/kkt_solver_basis.cc::KKTSolverBasis::_Solve (lines 73-192)
    ///
    /// Reduces the KKT system to preconditioned normal equations while taking
    /// free variables (which must be basic) into account. See [1, Section 6.4].
    ///
    /// [1] L. Schork, "Basis Preconditioning in Interior Point Methods", PhD
    ///     thesis (2018)
    #[allow(clippy::needless_range_loop)]
    fn _solve(
        &mut self,
        a: &Vector,
        b: &Vector,
        tol: f64,
        x: &mut Vector,
        y: &mut Vector,
        info: &mut Info,
    ) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let ai = model.AI();
        let mut rhs = Vector::new(m); // unpermuted right-hand side
        let mut work = Vector::new(m);
        info.errflag = 0;
        debug_assert!(self.factorized_);

        // Compute work = inverse(B')*v, where v[p] = a[basis[p]] if variable
        // basis[p] is free, and v[p] = 0 otherwise.
        let mut num_free: Int = 0;
        work.fill(0.0);
        for p in 0..m {
            let j = self.basis_.at(p);
            if self.basis_.status_of(j) == BASIC_FREE {
                work[p as usize] = a[j as usize];
                num_free += 1;
            }
        }
        if num_free > 0 {
            let work_copy = work.clone();
            self.basis_.solve_dense(&work_copy, &mut work, b'T');
        }

        // Compute rhs = inverse(B)*(N*D2[nonbasic]*(a[nonbasic]-N'*work)).
        rhs.fill(0.0);
        if num_free > 0 {
            for j in 0..n + m {
                if self.basis_.status_of(j) == NONBASIC {
                    let ju = j as usize;
                    let d2 = self.colscale_[ju] * self.colscale_[ju];
                    let mut alpha = a[ju] - dot_column(ai, j, &work);
                    alpha *= d2;
                    debug_assert!(alpha.is_finite());
                    scatter_column(ai, j, alpha, &mut rhs);
                }
            }
        } else {
            for j in 0..n + m {
                if self.basis_.status_of(j) == NONBASIC {
                    let ju = j as usize;
                    let d2 = self.colscale_[ju] * self.colscale_[ju];
                    let alpha = d2 * a[ju];
                    debug_assert!(alpha.is_finite());
                    scatter_column(ai, j, alpha, &mut rhs);
                }
            }
        }
        let rhs_copy = rhs.clone();
        self.basis_.solve_dense(&rhs_copy, &mut rhs, b'N');

        // Compute work = inverse(B)*b.
        self.basis_.solve_dense(b, &mut work, b'N');

        // Build rhs[p] = (rhs[p]-work[p])/D[j] + D[j]*a[j], where j = basis[p]
        // is not a free variable, and rhs[p] = 0 otherwise.
        for p in 0..m {
            let j = self.basis_.at(p);
            let pu = p as usize;
            if self.basis_.status_of(j) == BASIC {
                let ju = j as usize;
                let d = self.colscale_[ju];
                rhs[pu] = (rhs[pu] - work[pu]) / d + a[ju] * d;
                debug_assert!(rhs[pu].is_finite());
            } else {
                debug_assert!(self.basis_.status_of(j) == BASIC_FREE);
                rhs[pu] = 0.0;
            }
        }

        // Build permuted rhs in work.
        {
            let colperm = self.splitted_normal_matrix_.colperm();
            for k in 0..m {
                work[k as usize] = rhs[colperm[k as usize] as usize];
            }
        }

        // Solve normal equations.
        self.splitted_normal_matrix_.reset_time();
        let mut lhs = rhs; // don't need rhs any more
        lhs.fill(0.0);
        let mut cr = ConjugateResiduals::new(Rc::clone(&self.control_));
        cr.solve(
            &mut self.splitted_normal_matrix_,
            &work,
            tol,
            None,
            self.maxiter_,
            &mut lhs,
        );
        // HiGHS 1.16 #3513 (h16_ipx_scaled): ipm/ipx/kkt_solver_basis.cc lines 149-158 (1.16.0)
        if self.control_.h16_ipx_scaled() {
            let restart = (cr.errflag() == IPX_ERROR_cr_matrix_not_posdef
                || cr.errflag() == IPX_ERROR_cr_iter_limit)
                && self.splitted_normal_matrix_.use_extended_precision();
            if restart {
                // Restart the recurrence; its earlier products used double precision.
                info.kktiter2 += cr.iter();
                info.time_cr2 += cr.time();
                self.iter_ += cr.iter();
                lhs.fill(0.0);
                cr.solve(
                    &mut self.splitted_normal_matrix_,
                    &work,
                    tol,
                    None,
                    self.maxiter_,
                    &mut lhs,
                );
            }
        }
        info.errflag = cr.errflag();
        info.kktiter2 += cr.iter();
        info.time_cr2 += cr.time();
        info.time_cr2_NNt += self.splitted_normal_matrix_.time_NNt();
        info.time_cr2_B += self.splitted_normal_matrix_.time_B();
        info.time_cr2_Bt += self.splitted_normal_matrix_.time_Bt();
        self.iter_ += cr.iter();

        // Permute back solution to normal equations.
        {
            let colperm = self.splitted_normal_matrix_.colperm();
            for k in 0..m {
                y[colperm[k as usize] as usize] = lhs[k as usize];
            }
        }

        // Recover dual solution to KKT system.
        for p in 0..m {
            let j = self.basis_.at(p);
            let pu = p as usize;
            if self.basis_.status_of(j) == BASIC {
                y[pu] /= self.colscale_[j as usize];
                debug_assert!(y[pu].is_finite());
            } else {
                debug_assert!(self.basis_.status_of(j) == BASIC_FREE);
                debug_assert!(y[pu] == 0.0);
                y[pu] = a[j as usize]; // slot in solution to free basic variable
            }
        }
        let y_copy = y.clone();
        self.basis_.solve_dense(&y_copy, y, b'T');

        // Compute x[nonbasic] and work = b - N*x[nonbasic].
        work.assign(b);
        for j in 0..n + m {
            let ju = j as usize;
            let mut xj = 0.0;
            if self.basis_.status_of(j) == NONBASIC {
                xj = a[ju] - dot_column(ai, j, y);
                xj *= self.colscale_[ju] * self.colscale_[ju];
                debug_assert!(xj.is_finite());
                scatter_column(ai, j, -xj, &mut work);
            }
            x[ju] = xj;
        }

        // Compute x[basic].
        let work_copy = work.clone();
        self.basis_.solve_dense(&work_copy, &mut work, b'N');
        for p in 0..m {
            x[self.basis_.at(p) as usize] = work[p as usize];
        }
    }
    fn _iter(&self) -> Int {
        self.iter_
    }
    fn _basis_changes(&self) -> Int {
        self.basis_changes_
    }
    fn _basis(&self) -> Option<&Basis> {
        Some(self.basis_)
    }
}
