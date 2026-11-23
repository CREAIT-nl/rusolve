//! HiGHS source: lp_data/Highs.cpp (Highs::getIis), lp_data/HighsInterface.cpp
//! (Highs::getIisInterfaceReturn, Highs::getIisInterface, Highs::extractIis)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `writeIisModel` (Highs.cpp 880-884) is file-writer code and is
//! deferred. The `computeIis` C API is out of scope; `get_iis` is its
//! solver-side entry point.

use crate::prelude::*;
use crate::util::c_printf::{c_fmt_f, c_fmt_f_width, c_fmt_g_width, c_pad_right};

impl Highs {
    /// HiGHS: Highs.cpp::Highs::getIis (lines 2238-2242)
    pub fn get_iis(&mut self, iis: &mut HighsIis) -> HighsStatus {
        let return_status = self.get_iis_interface();
        *iis = self.iis_.clone();
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::getIisInterfaceReturn (lines 1833-1927)
    pub(crate) fn get_iis_interface_return(
        &mut self,
        return_status: HighsStatus,
        original_options: &HighsOptions,
        original_callback_active: &[bool],
    ) -> HighsStatus {
        // Restore options and callbacks
        self.options_ = original_options.clone();
        for i in kCallbackMin..=kCallbackMax {
            if original_callback_active[i as usize] {
                let _ = self.start_callback(i);
            }
        }

        // Exit early if there was an error
        if return_status == HighsStatus::kError {
            // Report final results
            self.iis_.report_final(original_options);
            return return_status;
        }
        debug_assert!(self.iis_.valid_);

        // A valid HighsIis instance is one for which the information is
        // known to be correct
        //
        // In the case of a system that's feasible this will give empty
        // HighsIis::col_index_ and HighsIis::row_index_
        //
        // If the system is infeasible, then it's possible that only an IS
        // (not an IIS) has been formed. This will be identified by the
        // HighsIis::col_bound_ and HighsIis::row_bound_ of the columns and
        // rows in HighsIis::col_index_ and HighsIis::row_index_ being
        // kIisStatusMaybeInConflict, rather than kIisStatusInConflict
        //
        // The columns and rows not in HighsIis::col_index_ and
        // HighsIis::row_index_ have HighsIis::col_bound_ and
        // HighsIis::row_bound_ values set to kIisStatusNotInConflict

        // If the IIS process has identified infeasibility, then set
        if self.iis_.status_ >= kIisModelStatusTimeLimit {
            self.model_status_ = HighsModelStatus::kInfeasible;
        }

        let lp = &self.model_lp_;
        let has_is = !self.iis_.col_index_.is_empty() || !self.iis_.row_index_.is_empty();
        let has_iis = self.iis_.status_ == kIisModelStatusIrreducible;

        let mut opts = self.options_.clone();
        opts.time_limit = kHighsInf;
        opts.simplex_iteration_limit = kHighsIInf;
        opts.objective_bound = kHighsInf;

        if has_iis {
            debug_assert!(has_is);
        }
        if has_is {
            // Construct the HighsIis LP
            self.iis_.set_lp(lp);
            // Check that the LP data are OK (correspond to original model
            // reduced to HighsIis col/row and bound data).
            let lp_data_ok = self.iis_.lp_data_ok(lp, &opts);
            debug_assert!(lp_data_ok);
            if !lp_data_ok {
                // Report final results
                self.iis_.report_final(original_options);
                return HighsStatus::kError;
            }
            // Check that the HighsIis LP has the property of being infeasible
            // and, if a true IIS is claimed, optimal/unbounded if any bound
            // is relaxed
            let lp_ok = self.iis_.lp_ok(&opts);
            if !lp_ok {
                // If we fail to prove infeasibility mark candidate set as invalid and
                // return with a warning
                self.iis_.valid_ = false;
                if self.iis_.status_ != kIisModelStatusTimeLimit {
                    self.iis_.status_ = kIisModelStatusUnknown;
                }
                // Report final results
                self.iis_.report_final(original_options);
                return HighsStatus::kWarning;
            }
        } else {
            debug_assert!(self.iis_.status_ <= kIisModelStatusReducible);
            debug_assert!(self.iis_.col_bound_.is_empty());
            debug_assert!(self.iis_.row_bound_.is_empty());
        }
        // Construct the ISS status vectors for cols and rows of original
        // model
        self.iis_.set_status(lp);
        // Check consistency of the col/row_index_ and col/row_status_
        let index_status_ok = self.iis_.index_status_ok(lp);
        debug_assert!(index_status_ok);
        if !index_status_ok {
            self.iis_.report_final(original_options);
            return HighsStatus::kError;
        }

        // Report final results
        self.iis_.report_final(original_options);
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::getIisInterface (lines 1929-2153)
    ///
    /// PORT NOTE: `original_options` is a clone of `options_` and
    /// `original_callback_active` a clone of `callback_.active`, as the C++
    /// copies. Calls that pass the live `options_` as `original_options`
    /// (lines 1940-1970) pass `&self.options_.clone()`.
    pub(crate) fn get_iis_interface(&mut self) -> HighsStatus {
        if self.model_status_ == HighsModelStatus::kOptimal
            || self.model_status_ == HighsModelStatus::kUnbounded
        {
            // Strange to call getIis for a model that's known to be feasible
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Calling Highs::getIis for a model that is known to be feasible\n",
            );
            self.iis_.clear();
            // No IIS exists, so validate the empty HighsIis instance
            self.iis_.valid_ = true;
            self.iis_.status_ = kIisModelStatusFeasible;
            let opts_clone = self.options_.clone();
            let cb_active = self.callback_.active.clone();
            return self.get_iis_interface_return(HighsStatus::kOk, &opts_clone, &cb_active);
        }
        // Early exit for existing valid IIS
        if self.iis_.valid_ {
            let opts_clone = self.options_.clone();
            let cb_active = self.callback_.active.clone();
            return self.get_iis_interface_return(HighsStatus::kOk, &opts_clone, &cb_active);
        }
        // Clear IIS
        self.iis_.clear();
        // Check for trivial IIS: empty infeasible row or inconsistent bounds
        let lp = &self.model_lp_;
        if self.iis_.trivial(lp, &self.options_) {
            let opts_clone = self.options_.clone();
            let cb_active = self.callback_.active.clone();
            return self.get_iis_interface_return(HighsStatus::kOk, &opts_clone, &cb_active);
        }
        let num_row = lp.num_row_;
        if num_row == 0 {
            // For an LP with no rows, the only scope for infeasibility is
            // inconsistent columns bounds - which has already been assessed,
            // so validate the empty HighsIis instance
            self.iis_.valid_ = true;
            let opts_clone = self.options_.clone();
            let cb_active = self.callback_.active.clone();
            return self.get_iis_interface_return(HighsStatus::kOk, &opts_clone, &cb_active);
        }
        // Look for infeasible rows based on row value bounds
        if self.iis_.row_value_bounds(lp, &self.options_) {
            let opts_clone = self.options_.clone();
            let cb_active = self.callback_.active.clone();
            return self.get_iis_interface_return(HighsStatus::kOk, &opts_clone, &cb_active);
        }
        // Don't continue with more expensive techniques if using the IIS
        // light strategy
        if self.options_.iis_strategy == kIisStrategyLight {
            let opts_clone = self.options_.clone();
            let cb_active = self.callback_.active.clone();
            return self.get_iis_interface_return(HighsStatus::kOk, &opts_clone, &cb_active);
        }
        // Clear IIS
        self.iis_.clear();
        // Save original options
        let original_options = self.options_.clone();
        // Save original active callbacks and disable all except for
        // kCallbackLogging and kCallbackSimplexInterrupt
        let original_callback_active = self.callback_.active.clone();
        for i in kCallbackMin..=kCallbackMax {
            if i != kCallbackLogging
                && i != kCallbackSimplexInterrupt
                && self.callback_.active[i as usize]
            {
                let _ = self.stop_callback(i);
            }
        }
        // Zero out all clocks and set time_limit to iis_time_limit
        self.zero_all_clocks();
        let _ = self.set_option_value_double("time_limit", self.options_.iis_time_limit);
        // Disallow kIterationLimit and kObjectiveBound model statuses
        let _ = self.set_option_value_int("simplex_iteration_limit", kHighsIInf);
        let _ = self.set_option_value_double("objective_bound", kHighsInf);
        // Disallow kUnboundedOrInfeasible model status
        let _ = self.set_option_value_bool("allow_unbounded_or_infeasible", false);
        // Run model with new options if needed
        if self.model_status_ == HighsModelStatus::kNotset
            || self.model_status_ == HighsModelStatus::kIterationLimit
            || self.model_status_ == HighsModelStatus::kObjectiveBound
            || self.model_status_ == HighsModelStatus::kInterrupt
            || self.model_status_ == HighsModelStatus::kTimeLimit
            || self.model_status_ == HighsModelStatus::kUnboundedOrInfeasible
        {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Model status is {}. Resolving to establish infeasibility before computing IIS\n\n",
                self.model_status_to_string(self.model_status_),
            );
            let mut simplex_time = -self.get_run_time();
            let mut simplex_iterations = -self.info_.simplex_iteration_count;
            let _run_status = self.optimize_model();
            simplex_time += self.get_run_time();
            simplex_iterations += self.info_.simplex_iteration_count;
            self.iis_.info_.update(simplex_time, simplex_iterations);
            crate::highs_log_user!(&self.options_.log_options, HighsLogType::kInfo, "\n");
        }
        if self.model_status_ == HighsModelStatus::kOptimal
            || self.model_status_ == HighsModelStatus::kUnbounded
        {
            // Model became feasible
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Model became feasible\n",
            );
            // No IIS exists, so validate the empty HighsIis instance
            self.iis_.valid_ = true;
            self.iis_.status_ = kIisModelStatusFeasible;
            self.iis_.strategy_ = self.options_.iis_strategy;
            return self.get_iis_interface_return(
                HighsStatus::kOk,
                &original_options,
                &original_callback_active,
            );
        } else if self.model_status_ == HighsModelStatus::kTimeLimit {
            // Time limit reached
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Time limit reached prior to establishing infeasibility\n",
            );
            self.iis_.status_ = kIisModelStatusTimeLimit;
            self.iis_.strategy_ = self.options_.iis_strategy;
            return self.get_iis_interface_return(
                HighsStatus::kError,
                &original_options,
                &original_callback_active,
            );
        } else if self.model_status_ != HighsModelStatus::kInfeasible {
            // kModelEmpty, kModelError, kSolveError, kMemoryLimit, kUnknown model
            // status
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Can not compute IIS for a model with status {}\n",
                self.model_status_to_string(self.model_status_),
            );
            self.iis_.strategy_ = self.options_.iis_strategy;
            return self.get_iis_interface_return(
                HighsStatus::kError,
                &original_options,
                &original_callback_active,
            );
        }
        debug_assert!(self.model_status_ == HighsModelStatus::kInfeasible);

        // Construct an IS from ray or lp
        // Line 2044: commented-out `kIisStrategyFromRay & options_.iis_strategy` dropped.
        let ray_option = false;
        let lp_option = (kIisStrategyFromLp & self.options_.iis_strategy) != 0;
        if ray_option && !self.ekk_instance_.status_.has_invert {
            // Model is known to be infeasible, and a dual ray option is
            // chosen, but it has no INVERT, presumably because infeasibility
            // detected in presolve, so solve without presolve
            let presolve = self.options_.presolve.clone();
            self.options_.presolve = kHighsOffString.to_string();

            let mut simplex_time = -self.get_run_time();
            let mut simplex_iterations = -self.info_.simplex_iteration_count;
            let _run_status = self.optimize_model();
            self.options_.presolve = presolve;
            simplex_time += self.get_run_time();
            simplex_iterations += self.info_.simplex_iteration_count;
            self.iis_.info_.update(simplex_time, simplex_iterations);

            // Model should remain infeasible!
            if self.model_status_ == HighsModelStatus::kTimeLimit {
                // Time limit reached
                self.iis_.status_ = kIisModelStatusTimeLimit;
                return self.get_iis_interface_return(
                    HighsStatus::kError,
                    &original_options,
                    &original_callback_active,
                );
            } else if self.model_status_ != HighsModelStatus::kInfeasible {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kError,
                    "Model status has switched from {} to {} when solving without presolve\n",
                    self.model_status_to_string(HighsModelStatus::kInfeasible),
                    self.model_status_to_string(self.model_status_),
                );
                return self.get_iis_interface_return(
                    HighsStatus::kError,
                    &original_options,
                    &original_callback_active,
                );
            }
        }
        let has_dual_ray = self.ekk_instance_.dual_ray_record_.index != kNoRayIndex;
        if ray_option && !has_dual_ray {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kWarning,
                "No known dual ray from which to compute IIS\n",
            );
        }
        let ray_option = ray_option && has_dual_ray;
        let mut return_status = HighsStatus::kOk;
        if ray_option {
            debug_assert!(has_dual_ray);
            // Compute the dual ray to identify an infeasible subset of rows
            debug_assert!(self.ekk_instance_.status_.has_invert);
            debug_assert!(!self.model_lp_.is_moved_);
            let mut rhs = vec![0.0f64; num_row as usize];
            let i_row = self.ekk_instance_.dual_ray_record_.index;
            rhs[i_row as usize] = 1.0;
            let mut dual_ray_value = vec![0.0f64; num_row as usize];
            self.basis_solve_interface(&rhs, &mut dual_ray_value, None, None, true);
            for i_row in 0..self.model_lp_.num_row_ {
                if dual_ray_value[i_row as usize] != 0.0 {
                    self.iis_.row_index_.push(i_row);
                }
            }
        } else if lp_option {
            // Full LP option chosen or no dual ray to use
            // Apply the elasticity filter to the whole model in order to
            // determine an infeasible subset of rows
            return_status = self.elasticity_filter(-1.0, -1.0, 1.0, None, None, None, true);
        }
        // Do not continue if not using kIisStrategyIrreducible
        if (kIisStrategyIrreducible & self.options_.iis_strategy) == 0 {
            return self.get_iis_interface_return(
                return_status,
                &original_options,
                &original_callback_active,
            );
        }
        // If neither ray and lp options were requested or if they fail to produce a
        // valid IS, make one consisting of all constraints
        if !self.iis_.valid_ || self.iis_.status_ != kIisModelStatusReducible {
            self.iis_.valid_ = true;
            self.iis_.status_ = kIisModelStatusReducible;
            for i_row in 0..self.model_lp_.num_row_ {
                self.iis_.row_index_.push(i_row);
            }
        }

        debug_assert!(!self.iis_.row_index_.is_empty());
        // A subset of infeasible rows, so at least have a reducible
        // infeasibility set

        //  Attempt to compute a true IIS
        //
        //  To do this the matrix must be column-wise
        self.model_lp_.a_matrix_.ensure_colwise();
        return_status = self.iis_.deduce(
            &self.model_lp_,
            &self.options_,
            &self.callback_,
            &self.basis_,
        );

        // Analyse the LP solution data
        let num_lp_solved = self.iis_.info_.num_lp_solved;
        let min_time = self.iis_.info_.min_simplex_time;
        let sum_time = self.iis_.info_.sum_simplex_times;
        let max_time = self.iis_.info_.max_simplex_time;
        let min_iterations = self.iis_.info_.min_simplex_iteration_count;
        let sum_iterations = self.iis_.info_.sum_simplex_iteration_counts;
        let max_iterations = self.iis_.info_.max_simplex_iteration_count;

        if kIisDevReport {
            // C `% 6d` / `% 6.2f`: the space flag prefixes a space when no
            // sign is written, then the field is right-justified to width 6.
            let space_flag = |f: String| -> String {
                if f.starts_with('-') {
                    f
                } else {
                    format!(" {}", f)
                }
            };
            let average_iterations = if num_lp_solved > 0 {
                (1.0 * sum_iterations as f64) / num_lp_solved as f64
            } else {
                0.0
            };
            let average_time = if num_lp_solved > 0 {
                sum_time / num_lp_solved as f64
            } else {
                0.0
            };
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                " {} cols, {} rows, {} LPs solved (min / average / max) iteration count ({:>6} / {} / {}) and time ({} / {} / {}) \n",
                self.iis_.col_index_.len() as i32,
                self.iis_.row_index_.len() as i32,
                num_lp_solved,
                min_iterations,
                c_fmt_g_width(average_iterations, 2, 6),
                c_pad_right(&space_flag(format!("{}", max_iterations)), 6),
                c_fmt_f_width(min_time, 2, 6),
                c_fmt_f_width(average_time, 2, 6),
                c_pad_right(&space_flag(c_fmt_f(max_time, 2)), 6)
            );
        }
        self.get_iis_interface_return(return_status, &original_options, &original_callback_active)
    }

    /// HiGHS: HighsInterface.cpp::Highs::extractIis (lines 2841-2861)
    ///
    /// PORT NOTE: nullable output arrays are `Option<&mut [HighsInt]>`.
    pub fn extract_iis(
        &self,
        num_iis_col: &mut HighsInt,
        num_iis_row: &mut HighsInt,
        iis_col_index: Option<&mut [HighsInt]>,
        iis_row_index: Option<&mut [HighsInt]>,
        iis_col_bound: Option<&mut [HighsInt]>,
        iis_row_bound: Option<&mut [HighsInt]>,
    ) -> HighsStatus {
        debug_assert!(self.iis_.valid_);
        *num_iis_col = self.iis_.col_index_.len() as HighsInt;
        *num_iis_row = self.iis_.row_index_.len() as HighsInt;
        let mut iis_col_index = iis_col_index;
        let mut iis_col_bound = iis_col_bound;
        if iis_col_index.is_some() || iis_col_bound.is_some() {
            for i_col in 0..*num_iis_col {
                if let Some(arr) = iis_col_index.as_deref_mut() {
                    arr[i_col as usize] = self.iis_.col_index_[i_col as usize];
                }
                if let Some(arr) = iis_col_bound.as_deref_mut() {
                    arr[i_col as usize] = self.iis_.col_bound_[i_col as usize];
                }
            }
        }
        let mut iis_row_index = iis_row_index;
        let mut iis_row_bound = iis_row_bound;
        if iis_row_index.is_some() || iis_row_bound.is_some() {
            for i_row in 0..*num_iis_row {
                if let Some(arr) = iis_row_index.as_deref_mut() {
                    arr[i_row as usize] = self.iis_.row_index_[i_row as usize];
                }
                if let Some(arr) = iis_row_bound.as_deref_mut() {
                    arr[i_row as usize] = self.iis_.row_bound_[i_row as usize];
                }
            }
        }
        HighsStatus::kOk
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// Column-wise LP builder: `cols[j] = (cost, lower, upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    fn quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    fn iis_lp() -> HighsLp {
        lp(
            &[
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, 1.0)]),
                (0.0, 0.0, kHighsInf, &[(0, 1.0), (1, -1.0)]),
            ],
            &[(-kHighsInf, 1.0), (2.0, kHighsInf)],
        )
    }

    #[test]
    fn get_iis_from_lp_strategy() {
        let mut highs = quiet(iis_lp());
        assert_eq!(
            highs
                .set_option_value_int("iis_strategy", kIisStrategyFromLp | kIisStrategyIrreducible),
            HighsStatus::kOk
        );
        let mut iis = HighsIis::default();
        assert_eq!(highs.get_iis(&mut iis), HighsStatus::kOk);
        assert!(iis.valid_);
        assert_eq!(iis.status_, kIisModelStatusIrreducible);
        // Expected values are the C++ HiGHS output on this LP (Highs::getIis
        // with iis_strategy = kIisStrategyFromLp | kIisStrategyIrreducible):
        // the IIS LP keeps x0 with its bound dropped (kIisBoundStatusFree),
        // and HighsIis::setStatus marks every column in col_index_ as
        // kIisStatusInConflict.
        let mut rows = iis.row_index_.clone();
        rows.sort();
        assert_eq!(rows, vec![0, 1]);
        assert_eq!(iis.col_index_, vec![0, 1]);
        assert_eq!(
            iis.col_bound_,
            vec![kIisBoundStatusFree, kIisBoundStatusLower]
        );
        assert_eq!(iis.row_index_, vec![0, 1]);
        assert_eq!(
            iis.row_bound_,
            vec![kIisBoundStatusUpper, kIisBoundStatusLower]
        );
        assert_eq!(
            iis.col_status_,
            vec![kIisStatusInConflict, kIisStatusInConflict]
        );
        assert_eq!(
            iis.row_status_,
            vec![kIisStatusInConflict, kIisStatusInConflict]
        );
    }

    #[test]
    fn feasible_lp_has_empty_iis() {
        let mut l = iis_lp();
        l.row_lower_[1] = -kHighsInf;
        let mut highs = quiet(l);
        let mut iis = HighsIis::default();
        assert_eq!(highs.get_iis(&mut iis), HighsStatus::kOk);
        assert!(iis.col_index_.is_empty() && iis.row_index_.is_empty());
        // C++ HiGHS on this LP (default iis_strategy, the light strategy):
        // valid_ true, status_ kIisModelStatusUnknown (feasibility is not
        // established), so HighsIis::nonIsStatus gives kIisStatusMaybeInConflict.
        assert!(iis.valid_);
        assert_eq!(iis.status_, kIisModelStatusUnknown);
        assert_eq!(
            iis.col_status_,
            vec![kIisStatusMaybeInConflict, kIisStatusMaybeInConflict]
        );
        assert_eq!(
            iis.row_status_,
            vec![kIisStatusMaybeInConflict, kIisStatusMaybeInConflict]
        );
    }
}
