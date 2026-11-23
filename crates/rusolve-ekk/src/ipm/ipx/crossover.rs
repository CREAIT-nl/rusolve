//! HiGHS source: ipm/ipx/crossover.h, crossover.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::basis::{Basis, copy_basic};
use super::control::{Control, format_int, sci2, textline};
use super::indexed_vector::{IndexedVector, for_each_nonzero};
use super::model::{
    Model, dual_infeasibility, dual_residual, primal_infeasibility, primal_residual,
};
use super::timer::Timer;
use super::utils::sortperm;
use super::*;
use crate::util::highs_utils::{cpp_max, cpp_min};
use std::rc::Rc;

/// HiGHS: ipx/crossover.h::Crossover
pub struct Crossover {
    control_: Rc<Control>,
    primal_pushes_: Int,
    dual_pushes_: Int,
    primal_pivots_: Int,
    dual_pivots_: Int,
    time_primal_: f64,
    time_dual_: f64,
}

impl Crossover {
    pub const kPivotZeroTol: f64 = 1e-5;

    /// HiGHS: ipx/crossover.cc::Crossover::Crossover (line 11)
    pub fn new(control: Rc<Control>) -> Self {
        Crossover {
            control_: control,
            primal_pushes_: 0,
            dual_pushes_: 0,
            primal_pivots_: 0,
            dual_pivots_: 0,
            time_primal_: 0.0,
            time_dual_: 0.0,
        }
    }

    /// HiGHS: ipx/crossover.cc::Crossover::PushAll (lines 13-71)
    ///
    /// `weights` may be `NULL` (lp_solver.cc:590).
    pub fn push_all(
        &mut self,
        basis: &mut Basis,
        x: &mut Vector,
        y: &mut Vector,
        z: &mut Vector,
        weights: Option<&[f64]>,
        info: &mut Info,
    ) {
        let model_rc = basis.model_rc();
        let model: &Model = &model_rc;
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        let perm = sortperm(n + m, weights, false);

        let mut log = String::new();
        log += &textline("Primal residual before push phase:");
        log += &sci2(primal_residual(model, x));
        log += "\n";
        log += &textline("Dual residual before push phase:");
        log += &sci2(dual_residual(model, y, z));
        log += "\n";
        self.control_.h_log(&mut log);

        // Run dual push phase.
        let mut dual_superbasics: Vec<Int> = Vec::new();
        for &j in perm.iter() {
            if basis.is_basic(j) && z[j as usize] != 0.0 {
                dual_superbasics.push(j);
            }
        }
        log += &textline("Number of dual pushes required:");
        log += &format!("{}\n", dual_superbasics.len());
        self.control_.h_log(&mut log);
        self.push_dual_x(basis, y, z, &dual_superbasics, x, info);
        debug_assert!(dual_infeasibility(model, x, z) == 0.0);
        if info.status_crossover != IPX_STATUS_optimal {
            return;
        }

        // Run primal push phase. Because z[j]==0 for all basic variables, none
        // of the primal variables is fixed at its bound.
        let mut primal_superbasics: Vec<Int> = Vec::new();
        for p in (0..perm.len()).rev() {
            let j = perm[p];
            let ju = j as usize;
            if basis.is_nonbasic(j)
                && x[ju] != lb[ju]
                && x[ju] != ub[ju]
                && !(lb[ju].is_infinite() && ub[ju].is_infinite() && x[ju] == 0.0)
            {
                primal_superbasics.push(j);
            }
        }
        log += &textline("Number of primal pushes required:");
        log += &format!("{}\n", primal_superbasics.len());
        self.control_.h_log(&mut log);
        self.push_primal(basis, x, &primal_superbasics, None, info);
        debug_assert!(primal_infeasibility(model, x) == 0.0);
        if info.status_crossover != IPX_STATUS_optimal {
            return;
        }

        let mut dbg = String::new();
        dbg += &textline("Primal residual after push phase:");
        dbg += &sci2(primal_residual(model, x));
        dbg += "\n";
        dbg += &textline("Dual residual after push phase:");
        dbg += &sci2(dual_residual(model, y, z));
        dbg += "\n";
        self.control_.debug(1, &dbg);
        info.status_crossover = IPX_STATUS_optimal;
    }

    /// HiGHS: ipx/crossover.cc::Crossover::PushPrimal(.., const bool* fixed_at_bound, ..) (lines 73-226)
    ///
    /// `fixed_at_bound` may be `nullptr` (crossover.cc:60).
    ///
    /// PORT NOTE: the C++ `throw std::logic_error` checks become `panic!`
    /// (the signatures carry no error channel).
    pub fn push_primal(
        &mut self,
        basis: &mut Basis,
        x: &mut Vector,
        variables: &[Int],
        fixed_at_bound: Option<&[bool]>,
        info: &mut Info,
    ) {
        let timer = Timer::default();
        let model_rc = basis.model_rc();
        let model: &Model = &model_rc;
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();
        let mut ftran = IndexedVector::new(m);
        let feastol = if model.dualized() {
            self.control_.dfeasibility_tol()
        } else {
            self.control_.pfeasibility_tol()
        };
        self.primal_pushes_ = 0;
        self.primal_pivots_ = 0;

        // Check that variables are nonbasic and that x satisfies bound condition.
        for &j in variables {
            if !basis.is_nonbasic(j) {
                panic!("invalid variable in Crossover::PushPrimal");
            }
        }
        for j in 0..n + m {
            let ju = j as usize;
            if x[ju] < lb[ju] || x[ju] > ub[ju] {
                panic!("bound condition violated in Crossover::PushPrimal");
            }
            if let Some(fixed) = fixed_at_bound {
                if fixed[ju] && x[ju] != lb[ju] && x[ju] != ub[ju] {
                    panic!("bound condition violated in Crossover::PushPrimal");
                }
            }
        }

        // Maintain a copy of primal basic variables and their bounds for faster
        // ratio test. Fixed-at-bound variables are handled by setting their bounds
        // equal.
        let mut xbasic = copy_basic(x, basis);
        let mut lbbasic = copy_basic(lb, basis);
        let mut ubbasic = copy_basic(ub, basis);
        if let Some(fixed) = fixed_at_bound {
            for p in 0..m {
                let j = basis.at(p) as usize;
                if fixed[j] {
                    lbbasic[p as usize] = x[j];
                    ubbasic[p as usize] = x[j];
                }
            }
        }

        self.control_.reset_print_interval();
        let batch = self.control_.crossover_batch_push() && fixed_at_bound.is_none();
        let mut batch_k: usize = 4;
        let mut rhs_idx: Vec<Int> = Vec::new();
        let mut rhs_val: Vec<f64> = Vec::new();
        let mut rhs_pos: Vec<Int> = vec![-1; m as usize];
        let mut next: usize = 0;
        while next < variables.len() {
            info.errflag = self.control_.interrupt_check(-1);
            if info.errflag != 0 {
                break;
            }

            // rusolve `lp_crossover_batch_push`: try to push the
            // next `batch_k` superbasics together; accept when the summed move keeps every basic
            // variable within its bounds, so no pivot would have happened for any of them.
            if batch && batch_k >= 2 {
                let ai = model.AI();
                let mut members: Vec<(usize, f64)> = Vec::with_capacity(batch_k);
                let mut q = next;
                while q < variables.len() && members.len() < batch_k {
                    let ju = variables[q] as usize;
                    q += 1;
                    if x[ju] == lb[ju]
                        || x[ju] == ub[ju]
                        || (x[ju] == 0.0 && lb[ju].is_infinite() && ub[ju].is_infinite())
                    {
                        continue;
                    }
                    let mt = if lb[ju].is_finite() && ub[ju].is_finite() {
                        if x[ju] - lb[ju] <= ub[ju] - x[ju] {
                            lb[ju]
                        } else {
                            ub[ju]
                        }
                    } else if lb[ju].is_finite() {
                        lb[ju]
                    } else if ub[ju].is_finite() {
                        ub[ju]
                    } else {
                        0.0
                    };
                    members.push((ju, mt));
                }
                if members.len() >= 2 {
                    rhs_idx.clear();
                    rhs_val.clear();
                    for &(ju, mt) in &members {
                        let step = x[ju] - mt;
                        let (b, e) = (ai.begin(ju as Int), ai.end(ju as Int));
                        for t in b as usize..e as usize {
                            let r = ai.rowidx()[t] as usize;
                            if rhs_pos[r] < 0 {
                                rhs_pos[r] = rhs_idx.len() as Int;
                                rhs_idx.push(r as Int);
                                rhs_val.push(0.0);
                            }
                            rhs_val[rhs_pos[r] as usize] += step * ai.values()[t];
                        }
                    }
                    for &r in &rhs_idx {
                        rhs_pos[r as usize] = -1;
                    }
                    basis.solve_rhs_lhs(&rhs_idx, &rhs_val, &mut ftran);
                    let mut ok = true;
                    for_each_nonzero(&ftran, |p, d| {
                        let pu = p as usize;
                        let v = xbasic[pu] + d;
                        if d.abs() > Self::kPivotZeroTol
                            && (v < lbbasic[pu] - feastol || v > ubbasic[pu] + feastol)
                        {
                            ok = false;
                        }
                    });
                    if ok {
                        for_each_nonzero(&ftran, |p, d| {
                            let pu = p as usize;
                            xbasic[pu] += d;
                            xbasic[pu] = cpp_max(xbasic[pu], lbbasic[pu]);
                            xbasic[pu] = cpp_min(xbasic[pu], ubbasic[pu]);
                        });
                        for &(ju, mt) in &members {
                            x[ju] = mt;
                        }
                        self.primal_pushes_ += members.len() as Int;
                        next = q;
                        batch_k = (batch_k * 2).min(512);
                        continue;
                    }
                    batch_k = (batch_k / 2).max(1);
                    continue;
                }
            }
            let jn = variables[next];
            let jnu = jn as usize;
            if x[jnu] == lb[jnu]
                || x[jnu] == ub[jnu]
                || (x[jnu] == 0.0 && lb[jnu].is_infinite() && ub[jnu].is_infinite())
            {
                // nothing to do
                next += 1;
                continue;
            }
            // Choose bound to push to. If the variable has two finite bounds, move
            // to the nearer. If it has none, move to zero.
            let mut move_to = 0.0;
            if lb[jnu].is_finite() && ub[jnu].is_finite() {
                move_to = if x[jnu] - lb[jnu] <= ub[jnu] - x[jnu] {
                    lb[jnu]
                } else {
                    ub[jnu]
                };
            } else if lb[jnu].is_finite() {
                move_to = lb[jnu];
            } else if ub[jnu].is_finite() {
                move_to = ub[jnu];
            }

            // A full step is such that x[jn]-step is at its bound.
            let mut step = x[jnu] - move_to;

            basis.solve_for_update_lhs(jn, &mut ftran);
            let mut block_at_lb = false;
            let pblock = self.primal_ratio_test(
                &xbasic,
                &ftran,
                &lbbasic,
                &ubbasic,
                step,
                feastol,
                &mut block_at_lb,
            );
            let jb = if pblock >= 0 { basis.at(pblock) } else { -1 };

            // If step was blocked, update basis and compute step size.
            if pblock >= 0 {
                let pb = pblock as usize;
                let pivot = ftran.get(pblock);
                debug_assert!(pivot != 0.0);
                if pivot.abs() < 1e-4 {
                    self.control_
                        .debug(3, &format!(" |pivot| = {}\n", sci2(pivot.abs())));
                }
                let mut exchanged = false;
                info.errflag = basis.exchange_if_stable(jb, jn, pivot, -1, &mut exchanged);
                if info.errflag == IPX_ERROR_user_interrupt
                    || info.errflag == IPX_ERROR_time_interrupt
                {
                    break;
                }
                if info.errflag != 0 {
                    self.control_.debug(
                        1,
                        &format!(
                            "{}{}\n",
                            textline("Minimum singular value of basis matrix:"),
                            sci2(basis.min_singular_value())
                        ),
                    );
                    break;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                self.primal_pivots_ += 1;
                // We must use lbbasic[pblock] and ubbasic[pblock] (and not lb[jb]
                // and ub[jb]) so that step is 0.0 if a fixed-at-bound variable
                // blocked.
                if block_at_lb {
                    step = (lbbasic[pb] - xbasic[pb]) / ftran.get(pblock);
                } else {
                    step = (ubbasic[pb] - xbasic[pb]) / ftran.get(pblock);
                }
            }
            // Update solution.
            if step != 0.0 {
                for_each_nonzero(&ftran, |p, pivot| {
                    let pu = p as usize;
                    xbasic[pu] += step * pivot;
                    xbasic[pu] = cpp_max(xbasic[pu], lbbasic[pu]);
                    xbasic[pu] = cpp_min(xbasic[pu], ubbasic[pu]);
                });
                x[jnu] -= step;
            }
            if pblock >= 0 {
                let pb = pblock as usize;
                // make clean
                x[jb as usize] = if block_at_lb {
                    lbbasic[pb]
                } else {
                    ubbasic[pb]
                };
                debug_assert!(x[jb as usize].is_finite());
                // Update copy of basic variables and bounds. Note: jn cannot be
                // a fixed-at-bound variable since it was pushed to a bound.
                xbasic[pb] = x[jnu];
                lbbasic[pb] = lb[jnu];
                ubbasic[pb] = ub[jnu];
            } else {
                x[jnu] = move_to; // make clean
                debug_assert!(x[jnu].is_finite());
            }

            self.primal_pushes_ += 1;
            next += 1;
            if batch && batch_k < 2 && pblock < 0 {
                batch_k = 2;
            }
            let mut log = String::new();
            log += &format!(
                " {} primal pushes remaining",
                format_int((variables.len() - next) as Int, 8)
            );
            log += &format!(" ({} pivots)\n", format_int(self.primal_pivots_, 7));
            self.control_.h_interval_log(&mut log);
        }
        for p in 0..m {
            x[basis.at(p) as usize] = xbasic[p as usize];
        }

        // Set status flag.
        if info.errflag == IPX_ERROR_user_interrupt {
            info.errflag = 0;
            // IPX quirk: PushPrimal reports a user interrupt in status_ipm.
            info.status_ipm = IPX_STATUS_user_interrupt;
        } else if info.errflag == IPX_ERROR_time_interrupt {
            info.errflag = 0;
            info.status_crossover = IPX_STATUS_time_limit;
        } else if info.errflag != 0 {
            info.status_crossover = IPX_STATUS_failed;
        } else {
            info.status_crossover = IPX_STATUS_optimal;
        }
        self.time_primal_ = timer.elapsed();
    }

    /// HiGHS: ipx/crossover.cc::Crossover::PushPrimal(.., const Vector& z, ..) (lines 228-233)
    pub fn push_primal_z(
        &mut self,
        basis: &mut Basis,
        x: &mut Vector,
        variables: &[Int],
        z: &Vector,
        info: &mut Info,
    ) {
        // std::valarray<bool> bound_restrict = z != 0.0;
        let bound_restrict: Vec<bool> = z.iter().map(|&zj| zj != 0.0).collect();
        // An empty valarray has no data; pass NULL (no variable is fixed).
        let fixed = if bound_restrict.is_empty() {
            None
        } else {
            Some(&bound_restrict[..])
        };
        self.push_primal(basis, x, variables, fixed, info);
    }

    /// HiGHS: ipx/crossover.cc::Crossover::PushDual(.., const int sign_restrict[], ..) (lines 235-352)
    ///
    /// PORT NOTE: the C++ `throw std::logic_error` checks become `panic!`.
    pub fn push_dual(
        &mut self,
        basis: &mut Basis,
        y: &mut Vector,
        z: &mut Vector,
        variables: &[Int],
        sign_restrict: &[i32],
        info: &mut Info,
    ) {
        let timer = Timer::default();
        let model_rc = basis.model_rc();
        let model: &Model = &model_rc;
        let m = model.rows();
        let n = model.cols();
        let mut btran = IndexedVector::new(m);
        let mut row = IndexedVector::new(n + m);
        let feastol = if model.dualized() {
            self.control_.pfeasibility_tol()
        } else {
            self.control_.dfeasibility_tol()
        };
        self.dual_pushes_ = 0;
        self.dual_pivots_ = 0;

        // Check that variables are basic and that z satisfies sign condition.
        for &j in variables {
            if !basis.is_basic(j) {
                panic!("invalid variable in Crossover::PushDual");
            }
        }
        for j in 0..n + m {
            let ju = j as usize;
            if ((sign_restrict[ju] & 1) != 0 && z[ju] < 0.0)
                || ((sign_restrict[ju] & 2) != 0 && z[ju] > 0.0)
            {
                panic!("sign condition violated in Crossover::PushDual");
            }
        }

        self.control_.reset_print_interval();
        let mut next: usize = 0;
        while next < variables.len() {
            info.errflag = self.control_.interrupt_check(-1);
            if info.errflag != 0 {
                break;
            }

            let jb = variables[next];
            let jbu = jb as usize;
            if z[jbu] == 0.0 {
                // nothing to do
                next += 1;
                continue;
            }
            // The update operation applied below is
            // y := y + step*btran, z := z - step*row, z[jb] := z[jb] - step,
            // where row is the tableau row for variable jb. In exact arithmetic
            // this leaves A'y+z unchanged.
            basis.tableau_row(jb, &mut btran, &mut row, false);
            let mut step = z[jbu];
            let jn = self.dual_ratio_test(z, &row, sign_restrict, step, feastol);

            // If step was blocked, update basis and compute step size.
            if jn >= 0 {
                let jnu = jn as usize;
                debug_assert!(basis.is_nonbasic(jn));
                let pivot = row.get(jn);
                debug_assert!(pivot != 0.0);
                if pivot.abs() < 1e-4 {
                    self.control_
                        .debug(3, &format!(" |pivot| = {}\n", sci2(pivot.abs())));
                }
                let mut exchanged = false;
                info.errflag = basis.exchange_if_stable(jb, jn, pivot, 1, &mut exchanged);
                if info.errflag == IPX_ERROR_user_interrupt
                    || info.errflag == IPX_ERROR_time_interrupt
                {
                    break;
                }
                if info.errflag != 0 {
                    self.control_.debug(
                        1,
                        &format!(
                            "{}{}\n",
                            textline("Minimum singular value of basis matrix:"),
                            sci2(basis.min_singular_value())
                        ),
                    );
                    break;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                self.dual_pivots_ += 1;
                step = z[jnu] / row.get(jn);
                // Update must move z[jb] toward zero.
                if (sign_restrict[jbu] & 1) != 0 {
                    debug_assert!(step >= 0.0);
                }
                if (sign_restrict[jbu] & 2) != 0 {
                    debug_assert!(step <= 0.0);
                }
            }
            // Update solution.
            if step != 0.0 {
                for_each_nonzero(&btran, |i, xv| {
                    y[i as usize] += step * xv;
                });
                for_each_nonzero(&row, |j, pivot| {
                    let ju = j as usize;
                    z[ju] -= step * pivot;
                    if (sign_restrict[ju] & 1) != 0 {
                        z[ju] = cpp_max(z[ju], 0.0);
                    }
                    if (sign_restrict[ju] & 2) != 0 {
                        z[ju] = cpp_min(z[ju], 0.0);
                    }
                });
                z[jbu] -= step;
            }
            if jn >= 0 {
                z[jn as usize] = 0.0; // make clean
            } else {
                debug_assert!(z[jbu] == 0.0);
            }

            self.dual_pushes_ += 1;
            next += 1;
            let mut log = String::new();
            log += &format!(
                " {} dual pushes remaining",
                format_int((variables.len() - next) as Int, 8)
            );
            log += &format!(" ({} pivots)\n", format_int(self.dual_pivots_, 7));
            self.control_.h_interval_log(&mut log);
        }

        // Set status flag.
        if info.errflag == IPX_ERROR_user_interrupt {
            info.errflag = 0;
            // IPX quirk: PushDual reports a user interrupt in status_ipm.
            info.status_ipm = IPX_STATUS_user_interrupt;
        } else if info.errflag == IPX_ERROR_time_interrupt {
            info.errflag = 0;
            info.status_crossover = IPX_STATUS_time_limit;
        } else if info.errflag != 0 {
            info.status_crossover = IPX_STATUS_failed;
        } else {
            info.status_crossover = IPX_STATUS_optimal;
        }
        self.time_dual_ = timer.elapsed();
    }

    /// HiGHS: ipx/crossover.cc::Crossover::PushDual(.., const Vector& x, ..) (lines 354-369)
    pub fn push_dual_x(
        &mut self,
        basis: &mut Basis,
        y: &mut Vector,
        z: &mut Vector,
        variables: &[Int],
        x: &Vector,
        info: &mut Info,
    ) {
        let model_rc = basis.model_rc();
        let model: &Model = &model_rc;
        let m = model.rows();
        let n = model.cols();
        let lb = model.lb();
        let ub = model.ub();

        let mut sign_restrict = vec![0i32; (n + m) as usize];
        for j in 0..sign_restrict.len() {
            if x[j] != ub[j] {
                sign_restrict[j] |= 1;
            }
            if x[j] != lb[j] {
                sign_restrict[j] |= 2;
            }
        }
        self.push_dual(basis, y, z, variables, &sign_restrict, info);
    }

    pub fn primal_pushes(&self) -> Int {
        self.primal_pushes_
    }
    pub fn dual_pushes(&self) -> Int {
        self.dual_pushes_
    }
    pub fn primal_pivots(&self) -> Int {
        self.primal_pivots_
    }
    pub fn dual_pivots(&self) -> Int {
        self.dual_pivots_
    }
    pub fn time_primal(&self) -> f64 {
        self.time_primal_
    }
    pub fn time_dual(&self) -> f64 {
        self.time_dual_
    }
    /// HiGHS: ipx/crossover.cc::Crossover::PrimalRatioTest (lines 371-428)
    fn primal_ratio_test(
        &mut self,
        xbasic: &Vector,
        ftran: &IndexedVector,
        lbbasic: &Vector,
        ubbasic: &Vector,
        step: f64,
        feastol: f64,
        block_at_lb: &mut bool,
    ) -> Int {
        let mut step = step;
        let mut pblock: Int = -1; // return value
        *block_at_lb = true;

        // First pass: determine maximum step size exploiting feasibility tol.
        for_each_nonzero(ftran, |p, pivot| {
            let pu = p as usize;
            if pivot.abs() > Self::kPivotZeroTol {
                // test block at lower bound
                if xbasic[pu] + step * pivot < lbbasic[pu] - feastol {
                    step = (lbbasic[pu] - xbasic[pu] - feastol) / pivot;
                    pblock = p;
                    *block_at_lb = true;
                }
                // test block at upper bound
                if xbasic[pu] + step * pivot > ubbasic[pu] + feastol {
                    step = (ubbasic[pu] - xbasic[pu] + feastol) / pivot;
                    pblock = p;
                    *block_at_lb = false;
                }
            }
        });

        // If the step was not blocked, we are done.
        if pblock < 0 {
            return pblock;
        }

        // Second pass: choose maximum pivot among all that block within step.
        pblock = -1;
        let mut max_pivot = Self::kPivotZeroTol;
        for_each_nonzero(ftran, |p, pivot| {
            let pu = p as usize;
            if pivot.abs() > max_pivot {
                // test block at lower bound
                if step * pivot < 0.0 {
                    let step_p = (lbbasic[pu] - xbasic[pu]) / pivot;
                    if step_p.abs() <= step.abs() {
                        pblock = p;
                        *block_at_lb = true;
                        max_pivot = pivot.abs();
                    }
                }
                // test block at upper bound
                if step * pivot > 0.0 {
                    let step_p = (ubbasic[pu] - xbasic[pu]) / pivot;
                    if step_p.abs() <= step.abs() {
                        pblock = p;
                        *block_at_lb = false;
                        max_pivot = pivot.abs();
                    }
                }
            }
        });
        debug_assert!(pblock >= 0);
        pblock
    }

    /// HiGHS: ipx/crossover.cc::Crossover::DualRatioTest (lines 430-477)
    fn dual_ratio_test(
        &mut self,
        z: &Vector,
        row: &IndexedVector,
        sign_restrict: &[i32],
        step: f64,
        feastol: f64,
    ) -> Int {
        let mut step = step;
        let mut jblock: Int = -1; // return value

        // First pass: determine maximum step size exploiting feasibility tol.
        for_each_nonzero(row, |j, pivot| {
            let ju = j as usize;
            if pivot.abs() > Self::kPivotZeroTol {
                if (sign_restrict[ju] & 1) != 0 && z[ju] - step * pivot < -feastol {
                    step = (z[ju] + feastol) / pivot;
                    jblock = j;
                    debug_assert!(z[ju] >= 0.0);
                    debug_assert!(step * pivot > 0.0);
                }
                if (sign_restrict[ju] & 2) != 0 && z[ju] - step * pivot > feastol {
                    step = (z[ju] - feastol) / pivot;
                    jblock = j;
                    debug_assert!(z[ju] <= 0.0);
                    debug_assert!(step * pivot < 0.0);
                }
            }
        });

        // If step was not block, we are done.
        if jblock < 0 {
            return jblock;
        }

        // Second pass: choose maximum pivot among all that block within step.
        jblock = -1;
        let mut max_pivot = Self::kPivotZeroTol;
        for_each_nonzero(row, |j, pivot| {
            let ju = j as usize;
            if pivot.abs() > max_pivot && (z[ju] / pivot).abs() <= step.abs() {
                if (sign_restrict[ju] & 1) != 0 && step * pivot > 0.0 {
                    jblock = j;
                    max_pivot = pivot.abs();
                }
                if (sign_restrict[ju] & 2) != 0 && step * pivot < 0.0 {
                    jblock = j;
                    max_pivot = pivot.abs();
                }
            }
        });
        debug_assert!(jblock >= 0);
        jblock
    }
}
