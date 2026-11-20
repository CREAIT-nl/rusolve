//! HiGHS source: ipm/ipx/basis.cc (lines 484-965)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::control::{format_int, format_str, ostream_double, sci2, textline, time_fmt};
use super::super::guess_basis::guess_basis;
use super::super::indexed_vector::{IndexedVector, dot_indexed, for_each_nonzero};
use super::super::sparse_matrix::scatter_column;
use super::super::timer::Timer;
use super::super::utils::{all_finite, find_max_abs};
use super::*;
use crate::ipm::basiclu::lu_internal::lu_time;
use crate::util::c_printf::{c_fmt_f_width, c_fmt_g_width, c_pad_right};
use crate::util::highs_utils::cpp_max;

impl Basis {
    /// HiGHS: ipx/basis.cc::Basis::AdaptToSingularFactorization (lines 484-506)
    pub(super) fn adapt_to_singular_factorization(&mut self) -> Int {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let mut rowperm: Vec<Int> = vec![0; m as usize];
        let mut colperm: Vec<Int> = vec![0; m as usize];
        let mut dependent_cols: Vec<Int> = Vec::new();

        self.lu_.borrow_mut().get_factors(
            None,
            None,
            Some(rowperm.as_mut_slice()),
            Some(colperm.as_mut_slice()),
            Some(&mut dependent_cols),
        );
        for &k in &dependent_cols {
            // Column p of the basis matrix was replaced by the i-th unit
            // column. Insert the corresponding slack variable jn into
            // position p of the basis.
            let p = colperm[k as usize];
            let i = rowperm[k as usize];
            let jb = self.basis_[p as usize];
            let jn = n + i;
            debug_assert!(self.map2basis_[jn as usize] < 0);
            self.basis_[p as usize] = jn;
            self.map2basis_[jn as usize] = p; // now BASIC at position p
            if jb >= 0 {
                self.map2basis_[jb as usize] = -1; // now NONBASIC
            }
        }
        dependent_cols.len() as Int
    }

    /// HiGHS: ipx/basis.cc::Basis::TightenLuPivotTol (lines 508-523)
    pub(super) fn tighten_lu_pivot_tol(&mut self) -> bool {
        let tol = self.lu_.borrow().pivottol();
        if tol <= 0.05 {
            self.lu_.borrow_mut().set_pivottol(0.1);
        } else if tol <= 0.25 {
            self.lu_.borrow_mut().set_pivottol(0.3);
        } else if tol <= 0.5 {
            self.lu_.borrow_mut().set_pivottol(0.9);
        } else {
            return false;
        }
        let pivottol = self.lu_.borrow().pivottol();
        let mut h = format!(
            " LU pivot tolerance tightened to {}\n",
            ostream_double(pivottol)
        );
        self.control_.h_log(&mut h);
        true
    }

    /// HiGHS: ipx/basis.cc::Basis::CrashBasis (lines 525-551)
    pub(super) fn crash_basis(&mut self, colweights: &[f64], interrupt: &mut bool) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();

        // Make a guess for a basis. Then use LU factorization with a strict
        // absolute pivot tolerance to remove dependent columns. This is not a
        // rank revealing factorization, but it detects many dependencies in
        // practice.
        let cols_guessed: Vec<Int> = guess_basis(&self.control_, &model, colweights, interrupt);
        if *interrupt {
            return;
        }
        debug_assert!((cols_guessed.len() as Int) <= m);
        debug_assert!((cols_guessed.len() as Int) == m); // at the moment

        // Initialize the Basis object and factorize the (partial) basis. If
        // basis_[p] is negative, the p-th column of the basis matrix is zero,
        // and a slack column will be inserted by CrashFacorize().
        for x in self.basis_.iter_mut() {
            *x = -1;
        }
        for x in self.map2basis_.iter_mut() {
            *x = -1;
        }
        for (k, &col) in cols_guessed.iter().enumerate() {
            self.basis_[k] = col;
            debug_assert!(self.map2basis_[col as usize] == -1); // must not have duplicates
            self.map2basis_[col as usize] = k as Int;
        }
        let mut num_dropped: Int = 0;
        self.crash_factorize(Some(&mut num_dropped), interrupt);
        if *interrupt {
            return;
        }
        self.control_.debug(
            1,
            &format!(
                "{}{}\n",
                textline("Number of columns dropped from guessed basis:"),
                num_dropped
            ),
        );
    }

    /// HiGHS: ipx/basis.cc::Basis::Repair (lines 590-624)
    pub(super) fn repair(&mut self, info: &mut Info, interrupt: &mut bool) {
        let model = Rc::clone(&self.model_);
        let n = model.cols();
        let m = model.rows();
        let mut work = Vector::new(m);
        info.basis_repairs = 0;

        loop {
            let (pmax, imax, pivot) = inverse_search(self, &mut work);
            if pmax < 0 || imax < 0 || !pivot.is_finite() {
                info.basis_repairs = -1;
                break;
            }
            if pivot.abs() < Self::kBasisRepairThreshold {
                break;
            }
            let jb = self.basis_[pmax as usize];
            let jn = n + imax;
            if !self.is_nonbasic(jn) {
                info.basis_repairs = -2;
                break;
            }
            if info.basis_repairs >= Self::kMaxBasisRepair {
                info.basis_repairs = -3;
                break;
            }
            self.solve_for_update(jb);
            self.solve_for_update(jn);
            self.crash_exchange(jb, jn, pivot, 0, None, interrupt);
            if *interrupt {
                return;
            }
            info.basis_repairs += 1;
            if self.control_.debug_enabled(3) {
                self.control_.debug(
                    3,
                    &format!(" basis repair: |pivot| = {}\n", sci2(pivot.abs())),
                );
            }
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::CrashFactorize (lines 626-664)
    ///
    /// `num_dropped` may be `nullptr` (null-checked at line 653).
    pub(super) fn crash_factorize(&mut self, num_dropped: Option<&mut Int>, interrupt: &mut bool) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let ai = model.AI();
        let timer = Timer::new(0.0);

        // Build column pointers for passing to LU factorization. A negative index
        // in basis_ means an empty slot. This option is kept for use by the crash
        // procedure, if an incomplete basis was constructed. The zero column in the
        // basis matrix causes a singularity in the LU factorization, so that the
        // empty slot will be replaced by a slack variable below.
        let mut begin = vec![0; m as usize];
        let mut end = vec![0; m as usize];
        for i in 0..m as usize {
            if self.basis_[i] >= 0 {
                begin[i] = ai.begin(self.basis_[i]);
                end[i] = ai.end(self.basis_[i]);
            } else {
                begin[i] = 0;
                end[i] = 0;
            }
        }
        // rusolve BASICLU time limit (not in the HiGHS reference): hand the remaining time to BASICLU
        let highs_time_limit = self.control_.time_limit();
        debug_assert!(highs_time_limit >= 0.0);
        self.lu_.borrow_mut().time_start(lu_time());
        let mut basiclu_time_limit = f64::INFINITY;
        if highs_time_limit < f64::INFINITY {
            basiclu_time_limit = highs_time_limit - self.control_.elapsed();
            if basiclu_time_limit <= 0.0 {
                *interrupt = true;
                return;
            }
        }
        self.lu_.borrow_mut().time_limit(basiclu_time_limit);
        let flag = self
            .lu_
            .borrow_mut()
            .factorize(&begin, &end, ai.rowidx(), ai.values(), true);
        if flag == IPX_ERROR_time_interrupt {
            *interrupt = true;
            return;
        }
        self.num_factorizations_ += 1;
        let fill = self.lu_.borrow().fill_factor();
        self.fill_factors_.push(fill);
        let mut ndropped: Int = 0;
        if (flag & 2) != 0 {
            ndropped = self.adapt_to_singular_factorization();
        }
        if let Some(nd) = num_dropped {
            *nd = ndropped;
        }

        self.time_factorize_ += timer.elapsed();
        self.factorization_is_fresh_ = true;

        if cfg!(debug_assertions) {
            // All empty slots must have been replaced by slack variables.
            for i in 0..m as usize {
                debug_assert!(self.basis_[i] >= 0);
            }
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::CrashExchange (lines 666-694)
    ///
    /// `num_dropped` is `nullptr` from `Repair` (line 619).
    pub(super) fn crash_exchange(
        &mut self,
        jb: Int,
        jn: Int,
        tableau_entry: f64,
        sys: i32,
        num_dropped: Option<&mut Int>,
        interrupt: &mut bool,
    ) {
        debug_assert!(self.is_basic(jb));
        debug_assert!(self.is_nonbasic(jn));
        if sys > 0 {
            // forward system needs to be solved
            self.solve_for_update(jn);
        }
        if sys < 0 {
            // transposed system needs to be solved
            self.solve_for_update(jb);
        }

        // Update basis.
        let ib = self.position_of(jb);
        debug_assert!(self.basis_[ib as usize] == jb);
        self.basis_[ib as usize] = jn;
        self.map2basis_[jn as usize] = ib; // status now BASIC
        self.map2basis_[jb as usize] = -1; // status now NONBASIC
        self.num_updates_ += 1;
        self.factorization_is_fresh_ = false;

        // Update factorization.
        let mut num_dropped = num_dropped;
        if let Some(nd) = num_dropped.as_deref_mut() {
            *nd = 0;
        }
        let timer = Timer::new(0.0);
        let err = self.lu_.borrow_mut().update(tableau_entry);
        self.time_update_ += timer.elapsed();
        let need_refactor = err != 0 || self.lu_.borrow_mut().need_fresh_factorization();
        if need_refactor {
            if self.control_.debug_enabled(3) {
                self.control_
                    .debug(3, " refactorization required in CrashExchange()\n");
            }
            self.crash_factorize(num_dropped, interrupt);
        }
    }

    /// HiGHS: ipx/basis.cc::Basis::PivotFreeVariablesIntoBasis (lines 696-803)
    pub(super) fn pivot_free_variables_into_basis(&mut self, colweights: &[f64], info: &mut Info) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let cost = model.c();
        let mut ftran = IndexedVector::new(m);
        let dependency_tol = cpp_max(0.0, self.control_.dependency_tol());
        info.errflag = 0;
        info.dependent_cols = 0;
        let mut stability_pivots: Int = 0;

        // Maintain stack of free nonbasic variables.
        let mut remaining: Vec<Int> = Vec::new();
        for j in 0..n + m {
            if colweights[j as usize].is_infinite() && self.map2basis_[j as usize] < 0 {
                remaining.push(j);
            }
        }
        self.control_.debug(
            1,
            &format!(
                "{}{}\n",
                textline("Number of free variables nonbasic:"),
                remaining.len()
            ),
        );

        self.control_.reset_print_interval();
        while let Some(&jn) = remaining.last() {
            debug_assert!(colweights[jn as usize].is_infinite());
            debug_assert!(self.map2basis_[jn as usize] < 0);
            info.errflag = self.control_.interrupt_check(-1);
            if info.errflag != 0 {
                return;
            }

            self.solve_for_update_lhs(jn, &mut ftran);
            let mut pmax: Int = -1;
            let mut pmax_nonfree: Int = -1;
            let mut fmax: f64 = 0.0;
            let mut fmax_nonfree: f64 = 0.0;
            {
                let basis_ = &self.basis_;
                for_each_nonzero(&ftran, |p, f| {
                    let f = f.abs();
                    if f > fmax {
                        fmax = f;
                        pmax = p;
                    }
                    if !colweights[basis_[p as usize] as usize].is_infinite() && f > fmax_nonfree {
                        fmax_nonfree = f;
                        pmax_nonfree = p;
                    }
                });
            }

            if fmax > 4.0 && fmax_nonfree < 1.0 {
                let jb = self.basis_[pmax as usize];
                debug_assert!(colweights[jb as usize].is_infinite());
                let mut exchanged = false;
                info.errflag = self.exchange_if_stable(jb, jn, ftran.get(pmax), -1, &mut exchanged);
                if info.errflag != 0 {
                    return;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                remaining.pop();
                remaining.push(jb);
                info.updates_start += 1;
                stability_pivots += 1;
            } else if fmax_nonfree <= dependency_tol {
                // jn cannot be pivoted into the basis. If we do not have an
                // unbounded primal ray yet, then test if column jn yields one.
                // Compute the change in the primal objective that is caused by
                // a unit increase of x[jn] with corresponding adjustment of
                // free basic variables.
                if info.cols_inconsistent == 0 {
                    let mut delta_obj = cost[jn as usize];
                    {
                        let basis_ = &self.basis_;
                        for_each_nonzero(&ftran, |p, f| {
                            let j = basis_[p as usize];
                            if colweights[j as usize].is_infinite() {
                                delta_obj -= cost[j as usize] * f;
                            }
                        });
                    }
                    if delta_obj.abs() > dependency_tol {
                        self.control_.debug(
                            1,
                            &format!(
                                "{}{}\n",
                                textline("Unbounded primal ray with objective change:"),
                                sci2(delta_obj)
                            ),
                        );
                        info.cols_inconsistent = 1;
                    }
                }
                info.dependent_cols += 1;
                remaining.pop();
            } else {
                let jb = self.basis_[pmax_nonfree as usize];
                let mut exchanged = false;
                info.errflag =
                    self.exchange_if_stable(jb, jn, ftran.get(pmax_nonfree), -1, &mut exchanged);
                if info.errflag != 0 {
                    return;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                remaining.pop();
                info.updates_start += 1;
            }
            let mut h = format!(" {} free variables remaining\n", remaining.len());
            self.control_.h_interval_log(&mut h);
        }
        self.control_.debug(
            1,
            &format!(
                "{}{}\n",
                textline("Number of free variables swapped for stability:"),
                stability_pivots
            ),
        );
    }

    /// HiGHS: ipx/basis.cc::Basis::PivotFixedVariablesOutOfBasis (lines 805-942)
    pub(super) fn pivot_fixed_variables_out_of_basis(
        &mut self,
        colweights: &[f64],
        info: &mut Info,
    ) {
        let model = Rc::clone(&self.model_);
        let m = model.rows();
        let n = model.cols();
        let ai = model.AI();
        let lb = model.lb();
        let ub = model.ub();
        let mut btran = IndexedVector::new(m);
        let mut row = IndexedVector::new(n + m);
        let dependency_tol = cpp_max(0.0, self.control_.dependency_tol());
        info.errflag = 0;
        info.dependent_rows = 0;
        let mut stability_pivots: Int = 0;

        // Build model right-hand side after subtracting fixed columns.
        // Needed for dual unboundedness test.
        let mut b_minus_fixed_columns: Vector = model.b().clone();
        for j in 0..n + m {
            let ju = j as usize;
            if lb[ju] == ub[ju] && lb[ju] != 0.0 {
                scatter_column(ai, j, -lb[ju], &mut b_minus_fixed_columns[..]);
            }
        }

        // Maintain stack of fixed basic variables.
        let mut remaining: Vec<Int> = Vec::new();
        for j in 0..n {
            // Structural variables with zero weight must not be in the crash basis.
            if colweights[j as usize] == 0.0 {
                debug_assert!(self.map2basis_[j as usize] < 0);
            }
        }
        for j in n..n + m {
            if colweights[j as usize] == 0.0 && self.map2basis_[j as usize] >= 0 {
                remaining.push(j);
            }
        }
        self.control_.debug(
            1,
            &format!(
                "{}{}\n",
                textline("Number of fixed variables basic:"),
                remaining.len()
            ),
        );

        self.control_.reset_print_interval();
        while let Some(&jb) = remaining.last() {
            debug_assert!(colweights[jb as usize] == 0.0);
            debug_assert!(self.map2basis_[jb as usize] >= 0);
            info.errflag = self.control_.interrupt_check(-1);
            if info.errflag != 0 {
                return;
            }

            self.tableau_row(jb, &mut btran, &mut row, false);
            let mut jmax: Int = -1;
            // C++ `jmax_nonfixed` is written in the scan but never read; it is
            // not carried over.
            let mut rmax: f64 = 0.0;
            let mut rmax_nonfixed: f64 = 0.0;
            for_each_nonzero(&row, |j, r| {
                // Ignore structural variables with zero weight.
                if j >= n || colweights[j as usize] != 0.0 {
                    let r = r.abs();
                    if r > rmax {
                        rmax = r;
                        jmax = j;
                    }
                    if colweights[j as usize] != 0.0 && r > rmax_nonfixed {
                        rmax_nonfixed = r;
                    }
                }
            });

            if rmax > 4.0 && rmax_nonfixed < 1.0 {
                debug_assert!(colweights[jmax as usize] == 0.0);
                let mut exchanged = false;
                info.errflag = self.exchange_if_stable(jb, jmax, row.get(jmax), 1, &mut exchanged);
                if info.errflag != 0 {
                    return;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                remaining.pop();
                remaining.push(jmax);
                info.updates_start += 1;
                stability_pivots += 1;
            } else if rmax_nonfixed <= dependency_tol {
                // jb cannot be pivoted out of the basis. If we do not have an
                // unbounded dual ray yet, then test if row jb-n yields one.
                // Compute the change in the dual objective that is caused by a
                // unit increase of y[jb-n] with corresponding adjustment of the
                // remaining y[i].
                if info.rows_inconsistent == 0 {
                    // Fix for #280 replaces this line by next
                    // double delta_obj = Dot(btran, model_.b());
                    let delta_obj = dot_indexed(&btran, &b_minus_fixed_columns[..]);
                    if delta_obj.abs() > dependency_tol {
                        self.control_.debug(
                            1,
                            &format!(
                                "{}{}\n",
                                textline("Unbounded dual ray with objective change:"),
                                sci2(delta_obj)
                            ),
                        );
                        info.rows_inconsistent = 1;
                    }
                }
                info.dependent_rows += 1;
                remaining.pop();
            } else {
                // jb can be exchanged by a non-fixed variable.
                // Among all numerically stable pivots, choose the one that
                // maximizes volume of the basis matrix.
                let mut jmax_scaled: Int = -1;
                let mut rmax_scaled: f64 = 0.0;
                for_each_nonzero(&row, |j, r| {
                    let r = r.abs();
                    let rscaled = r * colweights[j as usize];
                    if r >= 0.1 * rmax_nonfixed && rscaled > rmax_scaled {
                        rmax_scaled = rscaled;
                        jmax_scaled = j;
                    }
                });
                debug_assert!(jmax_scaled >= 0);
                let pivot = row.get(jmax_scaled);
                let mut exchanged = false;
                info.errflag = self.exchange_if_stable(jb, jmax_scaled, pivot, 1, &mut exchanged);
                if info.errflag != 0 {
                    return;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                remaining.pop();
                info.updates_start += 1;
            }
            let mut h = format!(
                "{} fixed variables remaining{}{}\n",
                format_int(remaining.len() as Int, 9),
                format_str("", 38),
                time_fmt(self.control_.elapsed())
            );
            self.control_.h_interval_log(&mut h);
        }
        self.control_.debug(
            1,
            &format!(
                "{}{}\n",
                textline("Number of fixed variables swapped for stability:"),
                stability_pivots
            ),
        );
    }

    /// HiGHS: ipx/basis.cc::Basis::reportBasisData (lines 952-965)
    pub fn report_basis_data(&self) {
        println!("\nBasis data");
        println!("    Num factorizations = {}", self.factorizations());
        println!("    Num updates = {}", self.updates_total());
        if self.num_ftran_ != 0 {
            println!(
                "    Average density of {} FTRANs is {}; sparse proportion = {}",
                c_pad_right(&self.num_ftran_.to_string(), 7),
                c_fmt_f_width(self.sum_ftran_density_ / self.num_ftran_ as f64, 4, 6),
                c_fmt_f_width(self.frac_ftran_sparse(), 4, 6)
            );
        }
        if self.num_btran_ != 0 {
            println!(
                "    Average density of {} BTRANs is {}; sparse proportion = {}",
                c_pad_right(&self.num_btran_.to_string(), 7),
                c_fmt_f_width(self.sum_btran_density_ / self.num_btran_ as f64, 4, 6),
                c_fmt_f_width(self.frac_btran_sparse(), 4, 6)
            );
        }
        println!(
            "    Mean fill-in {}",
            c_fmt_g_width(self.mean_fill(), 4, 11)
        );
        println!("    Max  fill-in {}", c_fmt_g_width(self.max_fill(), 4, 11));
    }
}

/// HiGHS: ipx/basis.cc::InverseSearch (static, lines 561-588)
///
/// Returns C++ `std::tuple<Int,Int,double>` as a Rust tuple.
pub(super) fn inverse_search(basis: &Basis, work: &mut Vector) -> (Int, Int, f64) {
    let mut inverse_max: f64 = 0.0;

    for (i, w) in work.iter_mut().enumerate() {
        *w = 1.0 / (i + 1) as f64;
    }

    loop {
        let rhs = work.clone();
        basis.solve_dense(&rhs, work, b'N');
        if !all_finite(&work[..]) {
            break;
        }
        let pmax = find_max_abs(&work[..]);
        work.fill(0.0);
        work[pmax as usize] = 1.0;
        let rhs = work.clone();
        basis.solve_dense(&rhs, work, b'T');
        if !all_finite(&work[..]) {
            break;
        }
        let imax = find_max_abs(&work[..]);
        let inverse_entry = work[imax as usize];
        if inverse_entry.abs() <= 2.0 * inverse_max {
            return (pmax, imax, inverse_entry);
        }
        inverse_max = inverse_entry.abs();
        work.fill(0.0);
        work[imax as usize] = 1.0;
    }
    (-1, -1, f64::INFINITY) // failure
}

/// HiGHS: ipx/basis.cc::CopyBasic (lines 944-950)
pub fn copy_basic(x: &Vector, basis: &Basis) -> Vector {
    let m = basis.model().rows();
    let mut xbasic = Vector::new(m);
    for p in 0..m {
        xbasic[p as usize] = x[basis.at(p) as usize];
    }
    xbasic
}
