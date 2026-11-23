//! HiGHS source: presolve/HPresolve.cpp (lines 6047-6474)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Limits, column bounds, and the run entry point. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::detach;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::removeSlacks (lines 6047-6095, 49 lines)
    pub fn remove_slacks(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        // SingletonColumns data structure appears not to be retained
        // throughout presolve
        for i_col in 0..self.model.num_col_ {
            if self.col_deleted[i_col as usize] != 0 {
                continue;
            }
            if self.colsize[i_col as usize] != 1 {
                continue;
            }
            if self.model.integrality_[i_col as usize] == HighsVarType::kInteger {
                continue;
            }
            let coliter = self.colhead[i_col as usize];
            let i_row = self.a_row[coliter as usize];
            debug_assert!(self.a_col[coliter as usize] == i_col);
            debug_assert!(self.row_deleted[i_row as usize] == 0);
            if !self.is_equation(i_row) {
                continue;
            }
            let lower = self.model.col_lower_[i_col as usize];
            let upper = self.model.col_upper_[i_col as usize];
            let cost = self.model.col_cost_[i_col as usize];
            let rhs = self.model.row_lower_[i_row as usize];
            let coeff = self.a_value[coliter as usize];
            debug_assert!(coeff != 0.0);
            // Slack is s = (rhs - a^Tx)/coeff
            //
            // Constraint bounds become:
            //
            // For coeff > 0 [rhs - coeff * upper, rhs - coeff * lower]
            //
            // For coeff < 0 [rhs - coeff * lower, rhs - coeff * upper]
            self.model.row_lower_[i_row as usize] = if coeff > 0.0 {
                rhs - coeff * upper
            } else {
                rhs - coeff * lower
            };
            self.model.row_upper_[i_row as usize] = if coeff > 0.0 {
                rhs - coeff * lower
            } else {
                rhs - coeff * upper
            };
            if cost != 0.0 {
                // Cost is (cost * rhs / coeff) + (col_cost - (cost/coeff) row_values)^Tx
                let multiplier = cost / coeff;
                let mut it = self.get_row_vector(i_row);
                while let Some(nonzero) = it.next_nz(self) {
                    let local_i_col = nonzero.index();
                    let local_value = nonzero.value();
                    self.model.col_cost_[local_i_col as usize] -= multiplier * local_value;
                }
                self.model.offset_ += multiplier * rhs;
            }
            //
            postsolve_stack.slack_col_substitution(
                i_row,
                i_col,
                rhs,
                self.get_row_vector(i_row).iter(self),
            );

            self.mark_col_deleted(i_col);

            self.unlink(coliter);
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::checkTimeLimit (lines 6097-6102, 6 lines)
    ///
    /// `timer->read()` -> `self.timer.as_ref().unwrap().read()`.
    pub fn check_time_limit(&mut self) -> HPresolveResult {
        debug_assert!(self.timer.is_some());
        if self.options.time_limit < kHighsInf
            && self.timer.as_ref().unwrap().read() >= self.options.time_limit
        {
            return HPresolveResult::kStopped;
        }
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::checkLimits (lines 6104-6151, 48 lines)
    ///
    /// `debugGetCheckCol()` / `debugGetCheckRow()` return -1 (see `sparsify.rs`).
    pub fn check_limits(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        let numreductions: usize = postsolve_stack.num_reductions();

        let mut debug_report = false;
        let check_col = self.debug_get_check_col();
        let check_row = self.debug_get_check_row();
        let mut col_bound_change = false;
        let mut row_bound_change = false;
        if check_col >= 0 || check_row >= 0 {
            if check_col >= 0 {
                col_bound_change = numreductions == 1
                    || postsolve_stack.debug_prev_col_lower
                        != self.model.col_lower_[check_col as usize]
                    || postsolve_stack.debug_prev_col_upper
                        != self.model.col_upper_[check_col as usize];
                postsolve_stack.debug_prev_col_lower = self.model.col_lower_[check_col as usize];
                postsolve_stack.debug_prev_col_upper = self.model.col_upper_[check_col as usize];
            }
            if check_row >= 0 {
                row_bound_change = numreductions == 1
                    || postsolve_stack.debug_prev_row_lower
                        != self.model.row_lower_[check_row as usize]
                    || postsolve_stack.debug_prev_row_upper
                        != self.model.row_upper_[check_row as usize];
                postsolve_stack.debug_prev_row_lower = self.model.row_lower_[check_row as usize];
                postsolve_stack.debug_prev_row_upper = self.model.row_upper_[check_row as usize];
            }
            debug_report = numreductions > postsolve_stack.debug_prev_numreductions;
        }
        // Debug output only: C++ `%11.4g` is rendered with Rust's default
        // float Display. Dead at the defaults (both check indices are -1).
        if check_col >= 0 && col_bound_change && debug_report {
            println!(
                "After reduction {:4}: col = {:4}[{}] has bounds [{}, {}]",
                numreductions - 1,
                check_col,
                self.model.col_names_[check_col as usize],
                self.model.col_lower_[check_col as usize],
                self.model.col_upper_[check_col as usize]
            );
            postsolve_stack.debug_prev_numreductions = numreductions;
        }
        if check_row >= 0 && row_bound_change && debug_report {
            println!(
                "After reduction {:4}: row = {:4}[{}] has bounds [{}, {}]",
                numreductions - 1,
                check_row,
                self.model.row_names_[check_row as usize],
                self.model.row_lower_[check_row as usize],
                self.model.row_upper_[check_row as usize]
            );
            postsolve_stack.debug_prev_numreductions = numreductions;
        }

        if (numreductions & 1023) == 0 {
            hpresolve_checked_call!(self.check_time_limit());
        }

        if numreductions >= self.reduction_limit {
            HPresolveResult::kStopped
        } else {
            HPresolveResult::kOk
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::storeCurrentProblemSize (lines 6153-6156, 4 lines)
    pub fn store_current_problem_size(&mut self) {
        self.old_num_col = self.model.num_col_ - self.num_deleted_cols;
        self.old_num_row = self.model.num_row_ - self.num_deleted_rows;
    }

    /// HiGHS: HPresolve.cpp::HPresolve::problemSizeReduction (lines 6158-6169, 12 lines)
    pub fn problem_size_reduction(&self) -> f64 {
        let col_reduction = 100.0
            * ((self.old_num_col - (self.model.num_col_ - self.num_deleted_cols)) as f64)
            / (self.old_num_col as f64);
        let row_reduction = 100.0
            * ((self.old_num_row - (self.model.num_row_ - self.num_deleted_rows)) as f64)
            / (self.old_num_row as f64);

        cpp_max(row_reduction, col_reduction)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeColBounds (lines 6171-6310, 140 lines)
    ///
    /// HiGHS 1.16 (#3411): HPresolve.cpp::computeColBounds (lines 4334-4490), rounding in computeBound (lines 4440-4455), gate h16_int_bound_rounding
    ///
    /// C++ default arguments: `boundCol = -1`, `boundColValue = -kHighsInf`, `boundColCoeffPattern = 0`,
    /// all four pointers `nullptr` (`None`).
    pub fn compute_col_bounds(
        &mut self,
        col: HighsInt,
        bound_col: HighsInt,
        bound_col_value: f64,
        bound_col_coeff_pattern: HighsInt,
        lower_bound: Option<&mut f64>,
        upper_bound: Option<&mut f64>,
        worst_case_lower_bound: Option<&mut f64>,
        worst_case_upper_bound: Option<&mut f64>,
    ) {
        debug_assert!(self.col_deleted[col as usize] == 0);
        debug_assert!(bound_col == -1 || self.col_deleted[bound_col as usize] == 0);
        debug_assert!(col != bound_col);
        // HiGHS 1.16 #3411 (h16_int_bound_rounding): HPresolve.cpp lines 4440-4455 (1.16.0)
        let h16_int_bound_rounding = self.options.h16_int_bound_rounding;

        // return if nothing to do
        if lower_bound.is_none()
            && upper_bound.is_none()
            && worst_case_lower_bound.is_none()
            && worst_case_upper_bound.is_none()
        {
            return;
        }

        // lambda for skipping non-zeros. The pointer tests are fixed for the
        // whole call, so they are evaluated once here.
        let lower_requested = lower_bound.is_some() || worst_case_lower_bound.is_some();
        let upper_requested = upper_bound.is_some() || worst_case_upper_bound.is_some();
        let skip_non_zero = |this: &HPresolve, row: HighsInt, val: f64| -> bool {
            let has_lower = this.yields_implied_lower_bound(row, val);
            let has_upper = this.yields_implied_upper_bound(row, val);
            (!upper_requested && lower_requested && !has_lower)
                || (!lower_requested && upper_requested && !has_upper)
                || (upper_requested && lower_requested && !has_lower && !has_upper)
        };

        // struct for storing non-zeros
        struct NonZeros {
            row: HighsInt,
            jval: f64,
            kval: f64,
        }

        let mut nzs: Vec<NonZeros> = Vec::with_capacity(self.colsize[col as usize] as usize);

        let store_triplet = |nzs: &mut Vec<NonZeros>, row: HighsInt, jval: f64, kval: f64| {
            // are coefficients required to have the same or opposite signs?
            // otherColCoeffPattern =  0: accept any coefficient pair
            // otherColCoeffPattern =  1: accept coefficients with same signs
            // otherColCoeffPattern = -1: accept coefficients with opposite signs
            if bound_col_coeff_pattern != 0
                && jval.is_sign_negative()
                    != ((bound_col_coeff_pattern as f64) * kval).is_sign_negative()
            {
                return;
            }
            nzs.push(NonZeros { row, jval, kval });
        };

        // store triplets (row, nonzero, nonzero) in a vector to speed up bound
        // computation. The findNonzero splay makes the row tree mutable, so
        // the column walks use the cursor form.
        if bound_col != -1 {
            if self.colsize[col as usize] < self.colsize[bound_col as usize] {
                let mut it = self.get_column_vector(col);
                while let Some(col_nz) = it.next_nz(self) {
                    // skip non-zero if it does not yield requested bounds
                    if skip_non_zero(self, col_nz.index(), col_nz.value()) {
                        continue;
                    }
                    let nz_pos = self.find_nonzero(col_nz.index(), bound_col);
                    if nz_pos == -1 {
                        continue;
                    }
                    store_triplet(
                        &mut nzs,
                        col_nz.index(),
                        col_nz.value(),
                        self.a_value[nz_pos as usize],
                    );
                }
            } else {
                let mut it = self.get_column_vector(bound_col);
                while let Some(col_nz) = it.next_nz(self) {
                    let nz_pos = self.find_nonzero(col_nz.index(), col);
                    if nz_pos == -1 {
                        continue;
                    }
                    // skip non-zero if it does not yield requested bounds
                    if skip_non_zero(self, col_nz.index(), self.a_value[nz_pos as usize]) {
                        continue;
                    }
                    store_triplet(
                        &mut nzs,
                        col_nz.index(),
                        self.a_value[nz_pos as usize],
                        col_nz.value(),
                    );
                }
            }
        } else {
            for col_nz in self.get_column_vector(col).iter(self) {
                // skip non-zero if it does not yield requested bounds
                if skip_non_zero(self, col_nz.index(), col_nz.value()) {
                    continue;
                }
                nzs.push(NonZeros {
                    row: col_nz.index(),
                    jval: col_nz.value(),
                    kval: -kHighsInf,
                });
            }
        }

        // initialise bounds
        let mut lower_bound = lower_bound;
        let mut upper_bound = upper_bound;
        let mut worst_case_lower_bound = worst_case_lower_bound;
        let mut worst_case_upper_bound = worst_case_upper_bound;
        if let Some(v) = lower_bound.as_deref_mut() {
            *v = -kHighsInf;
        }
        if let Some(v) = upper_bound.as_deref_mut() {
            *v = kHighsInf;
        }
        if let Some(v) = worst_case_lower_bound.as_deref_mut() {
            *v = -kHighsInf;
        }
        if let Some(v) = worst_case_upper_bound.as_deref_mut() {
            *v = kHighsInf;
        }

        // From here on HPresolve is only read.
        let this: &HPresolve = self;

        // lambda for actual bound computation
        let compute_bound = |triplet: &NonZeros,
                             rhs: f64,
                             direction: HighsInt,
                             is_lower_bound: bool,
                             is_worst_case_bound: bool|
         -> f64 {
            let bound = if (direction > 0 && !is_worst_case_bound)
                || (direction < 0 && is_worst_case_bound)
            {
                let residual = HighsCDouble::from(if is_worst_case_bound {
                    this.implied_row_bounds.get_residual_sum_lower(
                        &col_bound_arrays!(this),
                        triplet.row,
                        col,
                        triplet.jval,
                        bound_col,
                        triplet.kval,
                        bound_col_value,
                    )
                } else {
                    this.implied_row_bounds.get_residual_sum_lower_orig(
                        &col_bound_arrays!(this),
                        triplet.row,
                        col,
                        triplet.jval,
                        bound_col,
                        triplet.kval,
                        bound_col_value,
                    )
                });
                if residual == -kHighsInf {
                    return kHighsInf.copysign(triplet.jval);
                }
                f64::from((HighsCDouble::from(rhs) - residual) / triplet.jval)
            } else {
                let residual = HighsCDouble::from(if is_worst_case_bound {
                    this.implied_row_bounds.get_residual_sum_upper(
                        &col_bound_arrays!(this),
                        triplet.row,
                        col,
                        triplet.jval,
                        bound_col,
                        triplet.kval,
                        bound_col_value,
                    )
                } else {
                    this.implied_row_bounds.get_residual_sum_upper_orig(
                        &col_bound_arrays!(this),
                        triplet.row,
                        col,
                        triplet.jval,
                        bound_col,
                        triplet.kval,
                        bound_col_value,
                    )
                });
                if residual == kHighsInf {
                    return -(kHighsInf.copysign(triplet.jval));
                }
                f64::from((HighsCDouble::from(rhs) - residual) / triplet.jval)
            };
            // HiGHS 1.16 #3411: round bounds on integer columns (1.16.0 lines 4442-4454). Implied
            // bounds use the primal feasibility tolerance; worst-case bounds use it
            // divided by max(1, |jval|).
            if h16_int_bound_rounding
                && this.model.integrality_[col as usize] != HighsVarType::kContinuous
            {
                let tol = if is_worst_case_bound {
                    this.primal_feastol / cpp_max(1.0, triplet.jval.abs())
                } else {
                    this.primal_feastol
                };
                return if is_lower_bound {
                    (bound - tol).ceil()
                } else {
                    (bound + tol).floor()
                };
            }
            bound
        };

        // lambda for updating tightest bounds
        let mut update_bounds = |triplet: &NonZeros, rhs: f64, direction: HighsInt| {
            if direction as f64 * rhs == kHighsInf {
                return;
            }
            if (direction as f64) * triplet.jval < 0.0 {
                // lower bounds
                if let Some(lb) = lower_bound.as_deref_mut() {
                    *lb = cpp_max(*lb, compute_bound(triplet, rhs, direction, true, false));
                }
                if let Some(wlb) = worst_case_lower_bound.as_deref_mut() {
                    if *wlb != kHighsInf {
                        *wlb = cpp_max(*wlb, compute_bound(triplet, rhs, direction, true, true));
                    }
                }
            } else {
                // upper bounds
                if let Some(ub) = upper_bound.as_deref_mut() {
                    *ub = cpp_min(*ub, compute_bound(triplet, rhs, direction, false, false));
                }
                if let Some(wub) = worst_case_upper_bound.as_deref_mut() {
                    if *wub != -kHighsInf {
                        *wub = cpp_min(*wub, compute_bound(triplet, rhs, direction, false, true));
                    }
                }
            }
        };

        // compute bounds
        for triplet in nzs.iter() {
            // compute bounds using right-hand side (direction = 1) and left-hand side
            // (direction = -1)
            update_bounds(triplet, this.model.row_upper_[triplet.row as usize], 1);
            update_bounds(triplet, this.model.row_lower_[triplet.row as usize], -1);
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeDynamism (lines 6312-6324, 13 lines)
    ///
    /// C++ template over `HighsMatrixSlice<storageFormat>`; iterate `vector.iter(self)`.
    pub fn compute_dynamism<S: HighsTripletSlice>(&self, vector: S) -> HighsCDouble {
        let mut min_abs_coef = kHighsInf;
        let mut max_abs_coef = -kHighsInf;
        for nonzero in vector.iter(self) {
            let abs_coef = nonzero.value().abs();
            min_abs_coef = cpp_min(min_abs_coef, abs_coef);
            max_abs_coef = cpp_max(max_abs_coef, abs_coef);
        }
        HighsCDouble::from(max_abs_coef) / HighsCDouble::from(min_abs_coef)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::silentLog (lines 6326-6328, 3 lines)
    ///
    /// `mipsolver && mipsolver->mipdata_->numRestarts > 0`. Read through
    /// `as_ref()`: callers (`aggregator.rs`, `presolve_loop.rs`) run outside
    /// `with_mip_model!`, so `mipsolver` is present whenever it is a MIP.
    pub fn silent_log(&self) -> bool {
        self.mipsolver
            .as_ref()
            .is_some_and(|ms| ms.mipdata_.num_restarts > 0)
    }

    /// HiGHS 1.16 (#3250): HPresolve.cpp::HPresolve::moveCutsToPool (1.16.0 lines
    /// 4496-4558), under `mip_presolve_row_append`: the cut rows are those typed
    /// `kCut` in the postsolve stack, so rows appended by FM presolve stay in the model.
    fn move_cuts_to_pool(&mut self, postsolve_stack: &mut HighsPostsolveStack) {
        if self.mipsolver.is_none() {
            return;
        }
        let mut cutinds: Vec<HighsInt> = Vec::with_capacity(self.model.num_col_ as usize);
        let mut cutvals: Vec<f64> = Vec::with_capacity(self.model.num_col_ as usize);
        let mut numcuts: HighsInt = 0;
        for i in postsolve_stack.get_cut_rows() {
            numcuts += 1;
            self.store_row(i);
            cutinds.clear();
            cutvals.clear();
            for k in 0..self.rowpositions.len() {
                let j = self.rowpositions[k];
                cutinds.push(self.a_col[j as usize]);
                cutvals.push(self.a_value[j as usize]);
            }

            let rhs = self.model.row_upper_[i as usize];
            let integral = self.rowsize_integer[i as usize] + self.rowsize_impl_int[i as usize]
                == self.rowsize[i as usize]
                && self.row_coefficients_integral(i, 1.0);
            let len = cutinds.len() as HighsInt;
            with_mip_model!(self, ms => {
                detach!(ms.mipdata_.cutpools[0], |p| p.add_cut(
                    &mut ms,
                    &mut cutinds,
                    &mut cutvals,
                    len,
                    rhs,
                    integral,
                    true,
                    false,
                    false
                ));
            });

            self.mark_row_deleted(i);
            let rowpositions = std::mem::take(&mut self.rowpositions);
            for &j in rowpositions.iter() {
                self.unlink(j);
            }
            self.rowpositions = rowpositions;
        }

        if numcuts == 0 {
            return;
        }

        // Compact deleted cut rows. shrinkProblem must not be used here
        // because it replaces the cutpool with a new empty one, destroying
        // the cuts that were just added above.
        let old_num_row = self.model.num_row_;
        let mut new_row_index: Vec<HighsInt> = vec![0; old_num_row as usize];
        let mut new_num_row: HighsInt = 0;
        for i in 0..old_num_row as usize {
            if self.row_deleted[i] != 0 {
                new_row_index[i] = -1;
            } else {
                new_row_index[i] = new_num_row;
                new_num_row += 1;
            }
        }
        self.model.num_row_ = new_num_row;

        for i in 0..old_num_row as usize {
            let ni = new_row_index[i];
            if ni == -1 || ni == i as HighsInt {
                continue;
            }
            self.model.row_lower_[ni as usize] = self.model.row_lower_[i];
            self.model.row_upper_[ni as usize] = self.model.row_upper_[i];
        }
        let num_row = self.model.num_row_ as usize;
        self.model.row_lower_.truncate(num_row);
        self.model.row_upper_.truncate(num_row);
        self.model.row_names_.truncate(num_row);

        for i in 0..self.a_value.len() {
            if self.a_value[i] == 0.0 {
                continue;
            }
            debug_assert!(new_row_index[self.a_row[i] as usize] != -1);
            self.a_row[i] = new_row_index[self.a_row[i] as usize];
        }

        let nr = self.model.num_row_;
        self.mip_mut()
            .mipdata_
            .implications
            .reindex_rows(nr, &new_row_index);
    }

    /// HiGHS: HPresolve.cpp::HPresolve::run (lines 6330-6455, 126 lines)
    ///
    /// The MIP branches (cut pool, `mipsolver->mipdata_`) run in MIP mode (module doc).
    /// `toCSC(model->a_matrix_...)`: take the three vectors (module doc).
    pub fn run(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HighsModelStatus {
        self.presolve_status_ = HighsPresolveStatus::kNotSet;
        self.shrink_problem_enabled = true;
        postsolve_stack.debug_prev_numreductions = 0;
        postsolve_stack.debug_prev_col_lower = 0.0;
        postsolve_stack.debug_prev_col_upper = 0.0;
        postsolve_stack.debug_prev_row_lower = 0.0;
        postsolve_stack.debug_prev_row_upper = 0.0;
        // Presolve should only be called with a model that has a non-empty
        // constraint matrix unless it has no rows
        debug_assert!(self.model.a_matrix_.num_nz() != 0 || self.model.num_row_ == 0);
        let report_reductions = |this: &HPresolve, postsolve_stack: &HighsPostsolveStack| {
            if this.options.presolve != kHighsOffString && this.reduction_limit < kHighsSize_tInf {
                crate::highs_log_user!(
                    &this.options.log_options,
                    HighsLogType::kInfo,
                    "Presolve performed {} of {} permitted reductions\n",
                    postsolve_stack.num_reductions(),
                    this.reduction_limit
                );
            }
        };
        match self.presolve(postsolve_stack) {
            HPresolveResult::kStopped | HPresolveResult::kOk => {}
            HPresolveResult::kPrimalInfeasible => {
                self.presolve_status_ = HighsPresolveStatus::kInfeasible;
                report_reductions(self, postsolve_stack);
                return HighsModelStatus::kInfeasible;
            }
            HPresolveResult::kDualInfeasible => {
                self.presolve_status_ = HighsPresolveStatus::kUnboundedOrInfeasible;
                report_reductions(self, postsolve_stack);
                return HighsModelStatus::kUnboundedOrInfeasible;
            }
        }
        report_reductions(self, postsolve_stack);

        // HiGHS 1.16 #3167 (h16_presolve_rules): HPresolve.cpp::shrinkProblem line 1058 (1.16.0):
        // nothing to shrink before the presolve data structures exist.
        if !self.options.h16_presolve_rules || self.h16_has_presolve_data_structures() {
            self.shrink_problem(postsolve_stack);
        }

        if self.mipsolver.is_some() {
            let nnz = self.num_nonzeros();
            {
                let x1b = self.options.mip_presolve_clique_table;
                let ms = self.mip_mut();
                // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 7219-7220 at 7b890c255
                // (removed: the flag is scoped to runProbing / enumerateSolutions)
                if !x1b {
                    ms.mipdata_.cliquetable.set_presolve_flag(false);
                }
                ms.mipdata_.cliquetable.set_max_entries(nnz);
            }
            // C++ getCutPool() / getConflictPool() are cutpools[0] / conflictpools[0].
            // The propagators size their arrays by `numCol()` (the presolved model), so
            // they are built inside `with_mip_model!`.
            with_mip_model!(self, ms => {
                detach!(ms.mipdata_.domains[0], |d| d.add_cutpool(&mut ms, 0));
                detach!(ms.mipdata_.domains[0], |d| d.add_conflict_pool(&mut ms, 0));
            });

            if self.mip().mipdata_.num_restarts != 0 && self.options.mip_presolve_row_append {
                // HiGHS #3250: HPresolve.cpp line 4636 (1.16.0)
                self.move_cuts_to_pool(postsolve_stack);
            } else if self.mip().mipdata_.num_restarts != 0 {
                let mut cutinds: Vec<HighsInt> = Vec::with_capacity(self.model.num_col_ as usize);
                let mut cutvals: Vec<f64> = Vec::with_capacity(self.model.num_col_ as usize);
                let mut numcuts: HighsInt = 0;
                let mut i: HighsInt = self.model.num_row_ - 1;
                while i >= 0 {
                    // check if we already reached the original rows
                    if postsolve_stack.get_orig_row_index(i) < self.mip().orig_model_.num_row_ {
                        break;
                    }

                    // row is a cut, remove it from matrix but add to cutpool
                    numcuts += 1;
                    self.store_row(i);
                    cutinds.clear();
                    cutvals.clear();
                    for k in 0..self.rowpositions.len() {
                        let j = self.rowpositions[k];
                        cutinds.push(self.a_col[j as usize]);
                        cutvals.push(self.a_value[j as usize]);
                    }
                    let rhs = self.model.row_upper_[i as usize];
                    let integral = self.rowsize_integer[i as usize]
                        + self.rowsize_impl_int[i as usize]
                        == self.rowsize[i as usize]
                        && self.row_coefficients_integral(i, 1.0);
                    let len = cutinds.len() as HighsInt;
                    with_mip_model!(self, ms => {
                        detach!(ms.mipdata_.cutpools[0], |p| p.add_cut(
                            &mut ms,
                            &mut cutinds,
                            &mut cutvals,
                            len,
                            rhs,
                            integral,
                            true,
                            false,
                            false
                        ));
                    });

                    self.mark_row_deleted(i);
                    let rowpositions = std::mem::take(&mut self.rowpositions);
                    for &j in rowpositions.iter() {
                        self.unlink(j);
                    }
                    self.rowpositions = rowpositions;
                    i -= 1;
                }

                self.model.num_row_ -= numcuts;
                // C++ resize() to a smaller size: truncation.
                let num_row = self.model.num_row_ as usize;
                self.model.row_lower_.truncate(num_row);
                self.model.row_upper_.truncate(num_row);
                self.model.row_names_.truncate(num_row);
            }
        }

        // HiGHS 1.16 #3167 (h16_presolve_rules): HPresolve.cpp::run lines 4639-4643 (1.16.0): only
        // populate the model matrix if the presolve data structures were built.
        let h16_presolve_rules = self.options.h16_presolve_rules;
        if !h16_presolve_rules || self.h16_has_presolve_data_structures() {
            let mut a_value = std::mem::take(&mut self.model.a_matrix_.value_);
            let mut a_index = std::mem::take(&mut self.model.a_matrix_.index_);
            let mut a_start = std::mem::take(&mut self.model.a_matrix_.start_);
            self.to_csc(&mut a_value, &mut a_index, &mut a_start);
            self.model.a_matrix_.value_ = a_value;
            self.model.a_matrix_.index_ = a_index;
            self.model.a_matrix_.start_ = a_start;
        }

        if self.model.num_col_ == 0 {
            // Reduced to empty
            if let Some(ms) = self.mipsolver.as_mut() {
                if self.model.offset_ > ms.mipdata_.upper_limit {
                    self.presolve_status_ = HighsPresolveStatus::kInfeasible;
                    return HighsModelStatus::kInfeasible;
                }
                ms.mipdata_.lower_bound = 0.0;
            } else {
                // An LP with no columns must have no rows, unless the reduction
                // limit has been reached. As exposed when studying 2326
                let num_row_ok = self.model.num_row_ == 0
                    || postsolve_stack.num_reductions() >= self.reduction_limit;
                debug_assert!(num_row_ok);
                if !num_row_ok {
                    self.presolve_status_ = HighsPresolveStatus::kNotPresolved;
                    return HighsModelStatus::kNotset;
                }
            }
            self.presolve_status_ = HighsPresolveStatus::kReducedToEmpty;
            // Make sure that zero row activity from the columnless model is
            // consistent with the bounds
            return if self.zero_row_activity_feasible() {
                HighsModelStatus::kOptimal
            } else {
                HighsModelStatus::kInfeasible
            };
        } else if postsolve_stack.num_reductions() > 0 {
            // Reductions performed
            self.presolve_status_ = HighsPresolveStatus::kReduced;
        } else {
            // No reductions performed
            self.presolve_status_ = HighsPresolveStatus::kNotReduced;
        }

        if self.mipsolver.is_none() && self.options.use_implied_bounds_from_presolve {
            self.set_relaxed_implied_bounds();
        }

        debug_assert!(self.presolve_status_ != HighsPresolveStatus::kNotSet);
        HighsModelStatus::kNotset
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeIntermediateMatrix (lines 6457-6474, 18 lines)
    ///
    /// Only called by the unported `HPresolve::debug`; translate anyway.
    pub fn compute_intermediate_matrix(
        &mut self,
        flag_row: &mut Vec<HighsInt>,
        flag_col: &mut Vec<HighsInt>,
        numreductions: &mut usize,
    ) {
        self.shrink_problem_enabled = false;
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(flag_row.len() as HighsInt, flag_col.len() as HighsInt);
        self.set_reduction_limit(*numreductions);
        self.presolve(&mut stack);
        *numreductions = stack.num_reductions();

        let mut a_value = std::mem::take(&mut self.model.a_matrix_.value_);
        let mut a_index = std::mem::take(&mut self.model.a_matrix_.index_);
        let mut a_start = std::mem::take(&mut self.model.a_matrix_.start_);
        self.to_csc(&mut a_value, &mut a_index, &mut a_start);
        self.model.a_matrix_.value_ = a_value;
        self.model.a_matrix_.index_ = a_index;
        self.model.a_matrix_.start_ = a_start;

        for i in 0..self.model.num_row_ {
            flag_row[i as usize] = 1 - self.row_deleted[i as usize] as HighsInt;
        }
        for i in 0..self.model.num_col_ {
            flag_col[i as usize] = 1 - self.col_deleted[i as usize] as HighsInt;
        }
    }
}

#[cfg(test)]
mod x1b_psr_move_cuts_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;

    /// HiGHS #3250 (mip_presolve_row_append): HPresolve.cpp::moveCutsToPool (1.16.0 lines
    /// 4496-4558). Rows 0-1 are original, row 2 is a cut and row 3 is a row appended by
    /// presolve. Only the cut goes to the cut pool; the appended row stays in the model and
    /// moves to index 2, and the implications are reindexed to the new rows.
    #[test]
    fn x1b_psr_move_cuts_to_pool() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 4;
        lp.col_cost_ = vec![-1.0, -1.0, -1.0];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.integrality_ = vec![HighsVarType::kInteger; 3];
        // row 0: x0 + x1 <= 1, row 1: x1 + x2 <= 1, row 2 (cut): x0 + x2 <= 1,
        // row 3 (appended): 2 <= 2 x0 + x1 + x2
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf, -kHighsInf, 2.0];
        lp.row_upper_ = vec![1.0, 1.0, 1.0, kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 6, 9];
        lp.a_matrix_.index_ = vec![0, 2, 3, 0, 1, 3, 1, 2, 3];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_row_append = true;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.timer_.start();
        ms.new_data();
        ms.init();
        ms.setup_domain_propagation();
        // VLBs on x0 with origins in the cut row and in the appended row
        detach!(ms.mipdata_.implications, |imp| {
            imp.row_to_var_bounds.resize_with(4, HighsHashTree::default);
            imp.vlbs.resize_with(3, HighsHashTree::default);
            imp.vubs.resize_with(3, HighsHashTree::default);
            let vb = |origin| HighsImplications_VarBound {
                coef: 1.0,
                constant: 0.0,
                origin,
            };
            imp.vlbs[0].insert(HighsHashTableEntry::new(1, vb(2)));
            imp.vlbs[0].insert(HighsHashTableEntry::new(2, vb(3)));
        });

        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 3);
        stack.append_cuts_to_model_x1b(1);
        stack.append_rows_to_model(1);
        assert_eq!(stack.get_cut_rows(), vec![2]);

        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), -1));
        let ncuts0 = presolve.mip().mipdata_.cutpools[0].get_num_cuts();
        presolve.move_cuts_to_pool(&mut stack);

        assert_eq!(
            presolve.mip().mipdata_.cutpools[0].get_num_cuts(),
            ncuts0 + 1
        );
        assert_eq!(presolve.model.num_row_, 3);
        assert_eq!(presolve.row_deleted[2], 1);
        assert_eq!(presolve.model.row_lower_, vec![-kHighsInf, -kHighsInf, 2.0]);
        assert_eq!(presolve.model.row_upper_, vec![1.0, 1.0, kHighsInf]);
        // the nonzeros of the cut are unlinked; those of row 3 now live in row 2
        let mut rows: Vec<(HighsInt, HighsInt, f64)> = (0..presolve.a_value.len())
            .filter(|&k| presolve.a_value[k] != 0.0)
            .map(|k| (presolve.a_row[k], presolve.a_col[k], presolve.a_value[k]))
            .collect();
        rows.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(
            rows,
            vec![
                (0, 0, 1.0),
                (0, 1, 1.0),
                (1, 1, 1.0),
                (1, 2, 1.0),
                (2, 0, 2.0),
                (2, 1, 1.0),
                (2, 2, 1.0)
            ]
        );
        let imp = &presolve.mip().mipdata_.implications;
        assert_eq!(imp.vlbs[0].find(&1).unwrap().origin, -1);
        assert_eq!(imp.vlbs[0].find(&2).unwrap().origin, 2);
        assert_eq!(imp.row_to_var_bounds.len(), 3);
        assert_eq!(imp.row_to_var_bounds[2].find(&0), Some(&2));
    }
}
