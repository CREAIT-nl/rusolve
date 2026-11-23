//! HiGHS source: lp_data/HighsInterface.cpp (feasibleWrtBounds,
//! changeColBoundsInterface, setNonbasicStatusInterface)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice. Column-bound changes on the
//! incumbent LP. This is the TRUE hot start of the MIP node loop:
//! `ekk_instance_.update_status(LpAction::kNewBounds)` (clear.rs:286) only
//! drops `has_fresh_rebuild` and the objective-value flags, so invert, basis
//! and DSE weights survive into the next `run`.
//!
//! Dropped: `doubleUserDataNotNull` (slices are never NULL).
//! `changeRowBoundsInterface` (lines 1009-1058) is not ported: nothing in
//! the EkkEngine slice changes row bounds.

use crate::lp_data::highs_run::Highs;
use crate::lp_data::lp_assess::assess_bounds;
use crate::prelude::*;

impl Highs {
    /// HiGHS: HighsInterface.cpp::Highs::feasibleWrtBounds (lines 940-956, 17 lines)
    /// C++ default arguments: columns = true
    ///
    /// C++ copies the vectors; borrowing them is equivalent (read only).
    pub(crate) fn feasible_wrt_bounds(&self, columns: bool) -> bool {
        if self.info_.primal_solution_status != kSolutionStatusFeasible {
            return false;
        }
        let lp: &HighsLp = &self.model_lp_;
        let primal_feasibility_tolerance: f64 = self.options_.primal_feasibility_tolerance;
        let value: &[f64] = if columns {
            &self.solution_.col_value
        } else {
            &self.solution_.row_value
        };
        let lower: &[f64] = if columns {
            &lp.col_lower_
        } else {
            &lp.row_lower_
        };
        let upper: &[f64] = if columns {
            &lp.col_upper_
        } else {
            &lp.row_upper_
        };
        let dim: HighsInt = if columns { lp.num_col_ } else { lp.num_row_ };
        for i_x in 0..dim {
            let ix: usize = i_x as usize;
            if value[ix] < lower[ix] - primal_feasibility_tolerance {
                return false;
            }
            if value[ix] > upper[ix] + primal_feasibility_tolerance {
                return false;
            }
        }
        true
    }

    /// HiGHS: HighsInterface.cpp::Highs::changeColBoundsInterface (lines 958-1007, 50 lines)
    /// HiGHS 1.16 (#3386): HighsInterface.cpp::changeColBoundsInterface (lines 961-1010)
    ///
    /// * `num_col_bounds = dataSize(index_collection)` (`data_size`,
    ///   util/highs_utils.rs); `<= 0` returns kOk.
    /// * Local copies `local_col_lower/upper` of the first `num_col_bounds`
    ///   entries; for a set, `sort_set_data(set_num_entries_, &mut
    ///   index_collection.set_, Some(col_lower), Some(col_upper), None,
    ///   Some(&mut local_col_lower), Some(&mut local_col_upper), None)`
    ///   (util/highs_sort_set.rs). Note that C++ sorts the collection's own
    ///   `set_` in place, hence `&mut`.
    /// * `assess_bounds(&options, "col", 0, ic, &mut local_lower, &mut
    ///   local_upper, options.infinite_bound, None)`.
    /// * `change_lp_col_bounds` (highs_lp_modify.rs), then
    ///   `set_nonbasic_status_interface(ic, true)`, then the
    ///   `!basis_.useful && feasible_wrt_bounds(true)` choice between
    ///   `invalidate_model_status_and_info` and
    ///   `invalidate_model_status_solution_and_info` (highs_api_support.rs),
    ///   then `ekk_instance_.update_status(LpAction::kNewBounds)`.
    pub(crate) fn change_col_bounds_interface(
        &mut self,
        index_collection: &mut HighsIndexCollection,
        col_lower: &[f64],
        col_upper: &[f64],
    ) -> HighsStatus {
        let num_col_bounds: HighsInt = data_size(index_collection);
        // If a non-positive number of costs (may) need changing nothing needs to be
        // done
        if num_col_bounds <= 0 {
            return HighsStatus::kOk;
        }
        let n: usize = num_col_bounds as usize;
        // Take a copy of the cost that can be normalised
        let mut local_col_lower: Vec<f64> = col_lower[..n].to_vec();
        let mut local_col_upper: Vec<f64> = col_upper[..n].to_vec();
        // If changing the bounds for a set of columns, ensure that the
        // set and data are in ascending order
        if index_collection.is_set_ {
            sort_set_data(
                index_collection.set_num_entries_,
                &mut index_collection.set_,
                Some(col_lower),
                Some(col_upper),
                None,
                Some(local_col_lower.as_mut_slice()),
                Some(local_col_upper.as_mut_slice()),
                None,
            );
        }
        let return_status: HighsStatus = interpret_call_status(
            &self.options_.log_options,
            assess_bounds(
                &self.options_,
                "col",
                0,
                index_collection,
                &mut local_col_lower,
                &mut local_col_upper,
                self.options_.infinite_bound,
                None,
            ),
            HighsStatus::kOk,
            "assessBounds",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }

        change_lp_col_bounds(
            &mut self.model_lp_,
            index_collection,
            &local_col_lower,
            &local_col_upper,
        );
        // Update HiGHS basis status and (any) simplex move status of
        // nonbasic variables whose bounds have changed
        self.set_nonbasic_status_interface(index_collection, true);
        // Deduce the consequences of new col bounds
        if !self.basis_.useful && self.feasible_wrt_bounds(true) {
            // Retain the solution if there's no basis, and the solution is
            // feasible
            self.invalidate_model_status_and_info();
        } else {
            // Invalidate the solution
            self.invalidate_model_status_solution_and_info();
        }
        // Determine any implications for simplex data
        self.ekk_instance_.update_status(LpAction::kNewBounds);
        // HiGHS 1.16 #3386 (h16_add_no_modify_on_error): HighsInterface.cpp line 1009 (1.16.0)
        if self.options_.h16_add_no_modify_on_error {
            return_status
        } else {
            HighsStatus::kOk
        }
    }

    /// HiGHS: HighsInterface.cpp::Highs::setNonbasicStatusInterface (lines 1187-1327, 141 lines)
    ///
    /// Returns at once unless `basis_.valid`. For every nonbasic entry in
    /// the collection, recomputes the HiGHS status from the (new) bounds and,
    /// iff `ekk_instance_.status_.has_basis`, writes `nonbasic_flag_ =
    /// kNonbasicFlagTrue` and `nonbasic_move_` into `ekk_instance_.basis_`.
    /// Columns and rows use OPPOSITE move signs (columns: lower -> Up,
    /// rows: lower -> Dn); copy each branch verbatim. `fabs` is `f64::abs`,
    /// `highs_isInfinity` is `highs_is_infinity`.
    pub(crate) fn set_nonbasic_status_interface(
        &mut self,
        index_collection: &HighsIndexCollection,
        columns: bool,
    ) {
        if !self.basis_.valid {
            return;
        }
        let has_simplex_basis: bool = self.ekk_instance_.status_.has_basis;
        let highs_basis: &mut HighsBasis = &mut self.basis_;
        let simplex_basis: &mut SimplexBasis = &mut self.ekk_instance_.basis_;
        let lp: &HighsLp = &self.model_lp_;

        debug_assert!(ok(index_collection));
        let mut from_k: HighsInt = 0;
        let mut to_k: HighsInt = 0;
        limits(index_collection, &mut from_k, &mut to_k);
        let ix_dim: HighsInt = if columns { lp.num_col_ } else { lp.num_row_ };
        // Surely this is checked elsewhere
        debug_assert!(0 <= from_k && to_k < ix_dim);
        debug_assert!(from_k <= to_k);
        let mut set_from_ix: HighsInt = 0;
        let mut set_to_ix: HighsInt = 0;
        let mut ignore_from_ix: HighsInt = 0;
        let mut ignore_to_ix: HighsInt = -1;
        let mut current_set_entry: HighsInt = 0;
        // Given a basic-nonbasic partition, all status settings are defined
        // by the bounds unless boxed, in which case any definitive (ie not
        // just kNonbasic) existing status is retained. Otherwise, set to
        // bound nearer to zero
        for _k in from_k..=to_k {
            update_out_in_index(
                index_collection,
                &mut set_from_ix,
                &mut set_to_ix,
                &mut ignore_from_ix,
                &mut ignore_to_ix,
                &mut current_set_entry,
            );
            debug_assert!(set_to_ix < ix_dim);
            debug_assert!(ignore_to_ix < ix_dim);
            if columns {
                for i_col in set_from_ix..=set_to_ix {
                    let ic: usize = i_col as usize;
                    if highs_basis.col_status[ic] == HighsBasisStatus::kBasic {
                        continue;
                    }
                    // Nonbasic column
                    let lower: f64 = lp.col_lower_[ic];
                    let upper: f64 = lp.col_upper_[ic];
                    let mut status: HighsBasisStatus = highs_basis.col_status[ic];
                    let move_: i8;
                    if lower == upper {
                        if status == HighsBasisStatus::kNonbasic {
                            status = HighsBasisStatus::kLower;
                        }
                        move_ = kNonbasicMoveZe;
                    } else if !highs_is_infinity(-lower) {
                        // Finite lower bound so boxed or lower
                        if !highs_is_infinity(upper) {
                            // Finite upper bound so boxed
                            if status == HighsBasisStatus::kNonbasic {
                                // No definitive status, so set to bound nearer to zero
                                if lower.abs() < upper.abs() {
                                    status = HighsBasisStatus::kLower;
                                    move_ = kNonbasicMoveUp;
                                } else {
                                    status = HighsBasisStatus::kUpper;
                                    move_ = kNonbasicMoveDn;
                                }
                            } else if status == HighsBasisStatus::kLower {
                                move_ = kNonbasicMoveUp;
                            } else {
                                move_ = kNonbasicMoveDn;
                            }
                        } else {
                            // Lower (since upper bound is infinite)
                            status = HighsBasisStatus::kLower;
                            move_ = kNonbasicMoveUp;
                        }
                    } else if !highs_is_infinity(upper) {
                        // Upper
                        status = HighsBasisStatus::kUpper;
                        move_ = kNonbasicMoveDn;
                    } else {
                        // FREE
                        status = HighsBasisStatus::kZero;
                        move_ = kNonbasicMoveZe;
                    }
                    highs_basis.col_status[ic] = status;
                    if has_simplex_basis {
                        debug_assert!(move_ != kIllegalMoveValue);
                        simplex_basis.nonbasic_flag_[ic] = kNonbasicFlagTrue;
                        simplex_basis.nonbasic_move_[ic] = move_;
                    }
                }
            } else {
                for i_row in set_from_ix..=set_to_ix {
                    let ir: usize = i_row as usize;
                    if highs_basis.row_status[ir] == HighsBasisStatus::kBasic {
                        continue;
                    }
                    // Nonbasic column
                    let lower: f64 = lp.row_lower_[ir];
                    let upper: f64 = lp.row_upper_[ir];
                    let mut status: HighsBasisStatus = highs_basis.row_status[ir];
                    let move_: i8;
                    if lower == upper {
                        if status == HighsBasisStatus::kNonbasic {
                            status = HighsBasisStatus::kLower;
                        }
                        move_ = kNonbasicMoveZe;
                    } else if !highs_is_infinity(-lower) {
                        // Finite lower bound so boxed or lower
                        if !highs_is_infinity(upper) {
                            // Finite upper bound so boxed
                            if status == HighsBasisStatus::kNonbasic {
                                // No definitive status, so set to bound nearer to zero
                                if lower.abs() < upper.abs() {
                                    status = HighsBasisStatus::kLower;
                                    move_ = kNonbasicMoveDn;
                                } else {
                                    status = HighsBasisStatus::kUpper;
                                    move_ = kNonbasicMoveUp;
                                }
                            } else if status == HighsBasisStatus::kLower {
                                move_ = kNonbasicMoveDn;
                            } else {
                                move_ = kNonbasicMoveUp;
                            }
                        } else {
                            // Lower (since upper bound is infinite)
                            status = HighsBasisStatus::kLower;
                            move_ = kNonbasicMoveDn;
                        }
                    } else if !highs_is_infinity(upper) {
                        // Upper
                        status = HighsBasisStatus::kUpper;
                        move_ = kNonbasicMoveUp;
                    } else {
                        // FREE
                        status = HighsBasisStatus::kZero;
                        move_ = kNonbasicMoveZe;
                    }
                    highs_basis.row_status[ir] = status;
                    if has_simplex_basis {
                        debug_assert!(move_ != kIllegalMoveValue);
                        let flag_ix: usize = (lp.num_col_ + i_row) as usize;
                        simplex_basis.nonbasic_flag_[flag_ix] = kNonbasicFlagTrue;
                        simplex_basis.nonbasic_move_[flag_ix] = move_;
                    }
                }
            }
            if ignore_to_ix >= ix_dim - 1 {
                break;
            }
        }
    }
}
