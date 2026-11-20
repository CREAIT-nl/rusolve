//! HiGHS source: presolve/HPresolveInitialSweep.cpp (1.16.0 lines 1-446),
//! presolve/HPresolveUtils.cpp::computeSingletonRowBounds (1.16.0 lines 23-97),
//! presolve/HPresolve.cpp (1.16.0 lines 3994-4019)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS 1.16 LP initial sweep (#3167 #3210 #3442), behind `h16_presolve_rules`.
//! Before the presolve data structures are built, one pass over
//! the columns of the `HighsLp` removes empty and fixed columns, and one pass over the rows
//! removes empty, singleton and redundant rows, compressing the LP in place. It runs only
//! when `mipsolver` is unset, so MIP presolve never sees it.
//!
//! `HighsLp::fme_obj_col_` (Fourier–Motzkin, #3218 #3255) is updated as in 1.16; it is -1
//! unless MIP FM presolve reformulated the objective, so the updates are inert on LPs.
//! `HighsPostsolveStack::origRowType` is compressed inside `compress_index_maps`.

use super::*;
use crate::prelude::*;

/// HiGHS 1.16: HPresolveUtils.h::SingletonRowResult (lines 51-55)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // the C++ enumerator names
pub(super) enum H16SingletonRowResult {
    kRedundant,
    kPrimalInfeasible,
    kBoundsTightened,
}

/// HiGHS 1.16: HPresolveUtils.cpp::computeSingletonRowBounds (lines 23-97)
///
/// Returns `(result, lb, ub, lower_tightened, upper_tightened)`; the C++ out-parameters are
/// left unset on `kRedundant`, here they are then the column bounds and `false`.
#[allow(clippy::too_many_arguments)]
pub(super) fn h16_compute_singleton_row_bounds(
    val: f64,
    row_lower: f64,
    row_upper: f64,
    col_lower: f64,
    col_upper: f64,
    primal_feastol: f64,
    max_abs_col_val: f64,
    is_integral: bool,
) -> (H16SingletonRowResult, f64, f64, bool, bool) {
    if val > 0.0 {
        if col_upper * val <= row_upper + primal_feastol
            && col_lower * val >= row_lower - primal_feastol
        {
            return (
                H16SingletonRowResult::kRedundant,
                col_lower,
                col_upper,
                false,
                false,
            );
        }
    } else if col_lower * val <= row_upper + primal_feastol
        && col_upper * val >= row_lower - primal_feastol
    {
        return (
            H16SingletonRowResult::kRedundant,
            col_lower,
            col_upper,
            false,
            false,
        );
    }

    let mut new_col_upper = kHighsInf;
    let mut new_col_lower = -kHighsInf;
    if val > 0.0 {
        if row_upper != kHighsInf {
            new_col_upper = row_upper / val;
        }
        if row_lower != -kHighsInf {
            new_col_lower = row_lower / val;
        }
    } else {
        if row_upper != kHighsInf {
            new_col_lower = row_upper / val;
        }
        if row_lower != -kHighsInf {
            new_col_upper = row_lower / val;
        }
    }

    // use either the primal feasibility tolerance for the bound constraint or
    // for the singleton row including scaling, whichever is tighter.
    let bound_tol = (primal_feastol / 1.0f64.max(val.abs())).max(f64::EPSILON);

    let mut lower_tightened = new_col_lower > col_lower + bound_tol;
    let mut upper_tightened = new_col_upper < col_upper - bound_tol;

    let mut lb;
    let mut ub;
    if lower_tightened {
        if is_integral {
            new_col_lower = (new_col_lower - bound_tol).ceil();
        }
        lb = new_col_lower;
    } else {
        lb = col_lower;
    }

    if upper_tightened {
        if is_integral {
            new_col_upper = (new_col_upper + bound_tol).floor();
        }
        ub = new_col_upper;
    } else {
        ub = col_upper;
    }

    // check whether the bounds are equal in tolerances
    if ub <= lb + primal_feastol {
        // bounds could be infeasible or equal in tolerances, first check infeasible
        if ub < lb - primal_feastol {
            return (
                H16SingletonRowResult::kPrimalInfeasible,
                lb,
                ub,
                lower_tightened,
                upper_tightened,
            );
        }

        // bounds are equal in tolerances, if they have a slight infeasibility below
        // those tolerances or they have a slight numerical distance which changes
        // the largest contribution below feasibility tolerance then we can safely
        // set the bound to one of the values. To heuristically get rid of numerical
        // errors we choose the bound that was not tightened, or the midpoint if
        // both where tightened.
        if ub < lb || (ub > lb && (ub - lb) * val.abs().max(max_abs_col_val) <= primal_feastol) {
            if lower_tightened && upper_tightened {
                ub = 0.5 * (ub + lb);
                lb = ub;
                lower_tightened = lb > col_lower;
                upper_tightened = ub < col_upper;
            } else if lower_tightened {
                lb = ub;
                lower_tightened = lb > col_lower;
            } else {
                ub = lb;
                upper_tightened = ub < col_upper;
            }
        }
    }

    (
        H16SingletonRowResult::kBoundsTightened,
        lb,
        ub,
        lower_tightened,
        upper_tightened,
    )
}

/// HiGHS 1.16: HPresolveInitialSweep.h::HPresolveInitialSweep::Result
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // the C++ enumerator names
pub(super) enum H16SweepResult {
    kOk,
    kPrimalInfeasible,
    kDualInfeasible,
}

macro_rules! sweep_checked_call {
    ($call:expr) => {{
        let result = $call;
        if result != H16SweepResult::kOk {
            return result;
        }
    }};
}

/// HiGHS 1.16: HPresolveInitialSweep.h::HPresolveInitialSweep
pub(super) struct H16InitialSweep<'a> {
    lp: &'a mut HighsLp,
    options: &'a HighsOptions,
    allow_rule: &'a [bool],
    primal_feastol: f64,
    pub num_deleted_rows: HighsInt,
    pub num_deleted_cols: HighsInt,
}

impl<'a> H16InitialSweep<'a> {
    /// HiGHS 1.16: HPresolveInitialSweep.cpp::HPresolveInitialSweep (lines 24-32)
    pub(super) fn new(
        lp: &'a mut HighsLp,
        options: &'a HighsOptions,
        allow_rule: &'a [bool],
        primal_feastol: f64,
    ) -> Self {
        H16InitialSweep {
            lp,
            options,
            allow_rule,
            primal_feastol,
            num_deleted_rows: 0,
            num_deleted_cols: 0,
        }
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::checkColBounds (lines 34-57)
    fn check_col_bounds(&self, col: HighsInt, is_fixed: &mut bool) -> H16SweepResult {
        let ci = col as usize;
        let bound_diff = self.lp.col_upper_[ci] - self.lp.col_lower_[ci];
        let mut max_abs_col_value: f64 = 0.0;
        let a = &self.lp.a_matrix_;
        for i_el in a.start_[ci]..a.start_[ci + 1] {
            max_abs_col_value = a.value_[i_el as usize].abs().max(max_abs_col_value);
        }
        *is_fixed = false;
        // Check for simple infeasibility in the original LP should
        // already have been carried out in
        // HPresolve::checkOriginalLpBounds()
        debug_assert!(bound_diff >= 0.0);
        if bound_diff <= self.primal_feastol
            && (bound_diff <= self.options.small_matrix_value
                || max_abs_col_value * bound_diff <= self.primal_feastol)
        {
            // check for unboundedness
            if self.lp.col_lower_[ci].abs() == kHighsInf {
                return H16SweepResult::kDualInfeasible;
            }
            // column is fixed
            *is_fixed = true;
        }
        H16SweepResult::kOk
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::emptyCol (lines 59-93)
    fn empty_col(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> H16SweepResult {
        let ci = col as usize;
        let col_nnz = self.lp.a_matrix_.start_[ci + 1] - self.lp.a_matrix_.start_[ci];
        debug_assert!(col_nnz == 0);
        let mut cost = self.lp.col_cost_[ci];
        let lower = self.lp.col_lower_[ci];
        let upper = self.lp.col_upper_[ci];

        if (cost > 0.0 && lower == -kHighsInf) || (cost < 0.0 && upper == kHighsInf) {
            if cost.abs() <= self.options.dual_feasibility_tolerance {
                cost = 0.0;
            } else {
                return H16SweepResult::kDualInfeasible;
            }
        }
        let fixval;
        if cost > 0.0 {
            fixval = lower;
            if fixval == -kHighsInf {
                return H16SweepResult::kDualInfeasible;
            }
        } else if cost < 0.0 || upper.abs() < lower.abs() {
            fixval = upper;
            if fixval == kHighsInf {
                return H16SweepResult::kDualInfeasible;
            }
        } else if lower != -kHighsInf {
            fixval = lower;
        } else {
            fixval = 0.0;
        }
        debug_assert!(fixval != kHighsInf);
        // removedModelFixedCol (HighsPostsolveStack.h lines 847-858) is removedFixedCol
        // with the column read from the LP instead of a slice.
        postsolve_stack.removed_fixed_col(
            col,
            fixval,
            cost,
            std::iter::empty::<HighsSliceNonzero>(),
        );
        self.lp.offset_ += cost * fixval;
        if col == self.lp.fme_obj_col_ {
            self.lp.fme_obj_col_ = -1;
        }
        H16SweepResult::kOk
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::removeFixedCol (lines 95-110)
    fn remove_fixed_col(&mut self, col: HighsInt) {
        let ci = col as usize;
        let fixval = self.lp.col_lower_[ci];
        if col == self.lp.fme_obj_col_ {
            self.lp.fme_obj_col_ = -1;
        }
        let (start, end) = (
            self.lp.a_matrix_.start_[ci],
            self.lp.a_matrix_.start_[ci + 1],
        );
        for i_el in start..end {
            let colrow = self.lp.a_matrix_.index_[i_el as usize] as usize;
            let colval = self.lp.a_matrix_.value_[i_el as usize];
            if self.lp.row_lower_[colrow] != -kHighsInf {
                self.lp.row_lower_[colrow] -= colval * fixval;
            }
            if self.lp.row_upper_[colrow] != kHighsInf {
                self.lp.row_upper_[colrow] -= colval * fixval;
            }
        }
        self.lp.offset_ += self.lp.col_cost_[ci] * fixval;
        debug_assert!(self.lp.offset_.is_finite());
        self.lp.col_cost_[ci] = 0.0;
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::emptyRow (lines 112-119)
    fn empty_row(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
    ) -> H16SweepResult {
        let ri = row as usize;
        if self.lp.row_upper_[ri] < -self.primal_feastol
            || self.lp.row_lower_[ri] > self.primal_feastol
        {
            return H16SweepResult::kPrimalInfeasible;
        }
        postsolve_stack.redundant_row(row);
        H16SweepResult::kOk
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::getMaxAbsColVal (lines 121-127)
    fn get_max_abs_col_val(&self, col: HighsInt) -> f64 {
        let ci = col as usize;
        let a = &self.lp.a_matrix_;
        let mut max_val: f64 = 0.0;
        for i_el in a.start_[ci]..a.start_[ci + 1] {
            max_val = a.value_[i_el as usize].abs().max(max_val);
        }
        max_val
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::singletonRow (lines 129-153)
    fn singleton_row(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
        col: HighsInt,
        val: f64,
    ) -> H16SweepResult {
        debug_assert!(val.abs() > self.options.small_matrix_value);
        let (ri, ci) = (row as usize, col as usize);

        let (sr, lb, ub, lower_tightened, upper_tightened) = h16_compute_singleton_row_bounds(
            val,
            self.lp.row_lower_[ri],
            self.lp.row_upper_[ri],
            self.lp.col_lower_[ci],
            self.lp.col_upper_[ci],
            self.primal_feastol,
            self.get_max_abs_col_val(col),
            self.lp.integrality_[ci] != HighsVarType::kContinuous,
        );
        if sr == H16SingletonRowResult::kRedundant {
            postsolve_stack.redundant_row(row);
            return H16SweepResult::kOk;
        }
        if sr == H16SingletonRowResult::kPrimalInfeasible {
            return H16SweepResult::kPrimalInfeasible;
        }

        postsolve_stack.singleton_row(row, col, val, lower_tightened, upper_tightened);

        self.lp.col_lower_[ci] = lb;
        self.lp.col_upper_[ci] = ub;
        H16SweepResult::kOk
    }

    /// HiGHS 1.16: HPresolveInitialSweep.cpp::run (lines 155-444)
    pub(super) fn run(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> H16SweepResult {
        let have_col_names = !self.lp.col_names_.is_empty();
        let have_row_names = !self.lp.row_names_.is_empty();
        let original_num_col = self.lp.num_col_;
        let original_num_row = self.lp.num_row_;

        let mut num_fixed_col: HighsInt = 0;
        let mut num_empty_col: HighsInt = 0;
        let mut num_empty_row: HighsInt = 0;
        let mut num_singleton_row: HighsInt = 0;
        let mut num_redundant_row: HighsInt = 0;
        let mut num_col: HighsInt = 0;
        let mut nnz: HighsInt = 0;
        let mut is_fixed = false;

        let lp_num_row = self.lp.num_row_ as usize;
        let mut new_col_index: Vec<HighsInt> = vec![0; self.lp.num_col_ as usize];
        let mut row_count: Vec<HighsInt> = vec![0; lp_num_row];
        // Col of row is used to identify the column containing each
        // singleton row, and val_of_row the matrix entry of the singleton
        let mut col_of_row: Vec<HighsInt> = vec![-1; lp_num_row];
        let mut val_of_row: Vec<f64> = vec![0.0; lp_num_row];
        let compute_implied_bounds = self.allow_rule[kPresolveRuleRedundantRow as usize];
        // Compute the implied bounds on rows
        let mut implied_row_lower: Vec<HighsCDouble> = Vec::new();
        let mut implied_row_upper: Vec<HighsCDouble> = Vec::new();
        if compute_implied_bounds {
            implied_row_lower = vec![HighsCDouble::from(0.0); lp_num_row];
            implied_row_upper = vec![HighsCDouble::from(0.0); lp_num_row];
        }

        // Pass through the columns, identifying any that are empty or
        // fixed, so can be removed, updating the LP in place.
        for i_col in 0..self.lp.num_col_ {
            let ic = i_col as usize;
            let col_nnz = self.lp.a_matrix_.start_[ic + 1] - self.lp.a_matrix_.start_[ic];
            sweep_checked_call!(self.check_col_bounds(i_col, &mut is_fixed));
            if self.allow_rule[kPresolveRuleEmptyCol as usize] && col_nnz == 0 {
                // Remove empty column
                sweep_checked_call!(self.empty_col(postsolve_stack, i_col));
                new_col_index[ic] = -1;
                num_empty_col += 1;
                self.num_deleted_cols += 1;
            } else if self.allow_rule[kPresolveRuleFixedCol as usize] && is_fixed {
                // Remove fixed column
                let i_el = self.lp.a_matrix_.start_[ic];
                let a = &self.lp.a_matrix_;
                let col_vec = (i_el..i_el + col_nnz)
                    .map(|k| HighsSliceNonzero::new(a.index_[k as usize], a.value_[k as usize], k));
                postsolve_stack.removed_fixed_col(
                    i_col,
                    self.lp.col_lower_[ic],
                    self.lp.col_cost_[ic],
                    col_vec,
                );
                self.remove_fixed_col(i_col);
                new_col_index[ic] = -1;
                num_fixed_col += 1;
                self.num_deleted_cols += 1;
            } else {
                // Column retained: update the LP in place by shifting the cost,
                // bounds, any names, and the matrix data. Also compute this
                // column's contribution to the implied row bounds.
                let nc = num_col as usize;
                new_col_index[ic] = num_col;
                if num_col < i_col {
                    self.lp.col_cost_[nc] = self.lp.col_cost_[ic];
                    self.lp.col_lower_[nc] = self.lp.col_lower_[ic];
                    self.lp.col_upper_[nc] = self.lp.col_upper_[ic];
                    self.lp.integrality_[nc] = self.lp.integrality_[ic];
                    if have_col_names {
                        self.lp.col_names_[nc] = std::mem::take(&mut self.lp.col_names_[ic]);
                    }
                }
                let from_os = self.lp.a_matrix_.start_[ic];
                let new_col_start = nnz;
                for i_el in 0..col_nnz {
                    let src = (from_os + i_el) as usize;
                    let i_row = self.lp.a_matrix_.index_[src];
                    let value = self.lp.a_matrix_.value_[src];
                    let ir = i_row as usize;
                    row_count[ir] += 1;
                    col_of_row[ir] = num_col;
                    val_of_row[ir] = value;
                    self.lp.a_matrix_.index_[nnz as usize] = i_row;
                    self.lp.a_matrix_.value_[nnz as usize] = value;
                    nnz += 1;
                    if compute_implied_bounds {
                        let row_lower_bnd = if value > 0.0 {
                            self.lp.col_lower_[nc]
                        } else {
                            self.lp.col_upper_[nc]
                        };
                        let row_upper_bnd = if value > 0.0 {
                            self.lp.col_upper_[nc]
                        } else {
                            self.lp.col_lower_[nc]
                        };
                        if row_lower_bnd.abs() == kHighsInf {
                            implied_row_lower[ir] = HighsCDouble::from(-kHighsInf);
                        } else if f64::from(implied_row_lower[ir]) > -kHighsInf {
                            implied_row_lower[ir] += HighsCDouble::from(value) * row_lower_bnd;
                        }
                        if row_upper_bnd.abs() == kHighsInf {
                            implied_row_upper[ir] = HighsCDouble::from(kHighsInf);
                        } else if f64::from(implied_row_upper[ir]) < kHighsInf {
                            implied_row_upper[ir] += HighsCDouble::from(value) * row_upper_bnd;
                        }
                    }
                }
                self.lp.a_matrix_.start_[nc] = new_col_start;
                num_col += 1;
            }
        }
        self.lp.a_matrix_.start_[num_col as usize] = nnz;
        let num_removed_cols = num_empty_col + num_fixed_col;
        debug_assert!(self.num_deleted_cols == num_removed_cols);
        debug_assert!(num_col + num_removed_cols == self.lp.num_col_);
        let ncu = num_col as usize;
        self.lp.col_cost_.resize(ncu, 0.0);
        self.lp.col_lower_.resize(ncu, 0.0);
        self.lp.col_upper_.resize(ncu, 0.0);
        self.lp.integrality_.resize(ncu, HighsVarType::kContinuous);
        if have_col_names {
            self.lp.col_names_.resize(ncu, String::new());
        }
        self.lp.num_col_ = num_col;
        self.lp.a_matrix_.num_col_ = num_col;
        self.lp.a_matrix_.start_.resize(ncu + 1, 0);
        self.lp.a_matrix_.index_.resize(nnz as usize, 0);
        self.lp.a_matrix_.value_.resize(nnz as usize, 0.0);
        if self.lp.fme_obj_col_ >= 0 {
            self.lp.fme_obj_col_ = new_col_index[self.lp.fme_obj_col_ as usize];
        }
        // compressColIndexMap (HighsPostsolveStack.cpp lines 52-55): origColIndex only.
        postsolve_stack.compress_index_maps(&[], &new_col_index);

        let primal_feastol = self.primal_feastol;
        let is_redundant = |lp: &HighsLp, i_row: usize| -> bool {
            f64::from(implied_row_lower[i_row]) >= lp.row_lower_[i_row] - primal_feastol
                && f64::from(implied_row_upper[i_row]) <= lp.row_upper_[i_row] + primal_feastol
        };

        let allow_empty_row = self.allow_rule[kPresolveRuleEmptyRow as usize];
        let allow_singleton_row = self.allow_rule[kPresolveRuleSingletonRow as usize];
        let allow_redundant_row = self.allow_rule[kPresolveRuleRedundantRow as usize];

        // Row pass: count empty, singleton, and redundant rows
        for i_row in 0..lp_num_row {
            if allow_empty_row && row_count[i_row] == 0 {
                num_empty_row += 1;
            } else if allow_singleton_row && row_count[i_row] == 1 {
                num_singleton_row += 1;
            } else if allow_redundant_row && is_redundant(&*self.lp, i_row) {
                num_redundant_row += 1;
            }
        }

        let num_removed_rows = num_empty_row + num_singleton_row + num_redundant_row;
        if num_empty_row != 0 || num_singleton_row != 0 || num_redundant_row != 0 {
            let mut num_row: HighsInt = 0;
            let mut has_singleton_row: Vec<bool> = vec![false; self.lp.num_col_ as usize];
            let mut new_row_index: Vec<HighsInt> = vec![0; lp_num_row];
            for i_row in 0..lp_num_row {
                let row = i_row as HighsInt;
                if allow_empty_row && row_count[i_row] == 0 {
                    // Remove empty row
                    sweep_checked_call!(self.empty_row(postsolve_stack, row));
                    new_row_index[i_row] = -1;
                    self.num_deleted_rows += 1;
                } else if allow_singleton_row && row_count[i_row] == 1 {
                    // Remove singleton row
                    has_singleton_row[col_of_row[i_row] as usize] = true;
                    debug_assert!(val_of_row[i_row] != 0.0);
                    sweep_checked_call!(self.singleton_row(
                        postsolve_stack,
                        row,
                        col_of_row[i_row],
                        val_of_row[i_row]
                    ));
                    new_row_index[i_row] = -1;
                    self.num_deleted_rows += 1;
                } else if allow_redundant_row && is_redundant(&*self.lp, i_row) {
                    // Remove redundant row
                    postsolve_stack.redundant_row(row);
                    new_row_index[i_row] = -1;
                    self.num_deleted_rows += 1;
                } else {
                    new_row_index[i_row] = num_row;
                    let nr = num_row as usize;
                    if nr < i_row {
                        self.lp.row_lower_[nr] = self.lp.row_lower_[i_row];
                        self.lp.row_upper_[nr] = self.lp.row_upper_[i_row];
                        if have_row_names {
                            self.lp.row_names_[nr] = std::mem::take(&mut self.lp.row_names_[i_row]);
                        }
                    }
                    num_row += 1;
                }
            }
            debug_assert!(self.num_deleted_rows == num_removed_rows);
            debug_assert!(num_row + num_removed_rows == self.lp.num_row_);

            let a = &mut self.lp.a_matrix_;
            if num_redundant_row == 0 {
                // Only removing entries corresponding to singleton rows so
                // there are few to remove and it can be done efficiently
                nnz = 0;
                let mut from_col: HighsInt = 0;
                // Lambda for shifting column data and updating row indices
                let shift_cols = |a: &mut HighsSparseMatrix,
                                  nnz: &mut HighsInt,
                                  from_col: HighsInt,
                                  to_col: HighsInt| {
                    for i_col in from_col..to_col {
                        let ic = i_col as usize;
                        let from_os = a.start_[ic];
                        let col_nnz = a.start_[ic + 1] - from_os;
                        let new_col_start = *nnz;
                        for i_el in 0..col_nnz {
                            let src = (from_os + i_el) as usize;
                            let i_row = a.index_[src];
                            let new_row = new_row_index[i_row as usize];
                            debug_assert!(new_row >= 0);
                            a.index_[*nnz as usize] = new_row;
                            a.value_[*nnz as usize] = a.value_[src];
                            *nnz += 1;
                        }
                        a.start_[ic] = new_col_start;
                    }
                };
                for i_col0 in 0..num_col {
                    let ic0 = i_col0 as usize;
                    if !has_singleton_row[ic0] {
                        continue;
                    }
                    // Column iCol0 contains a row singleton, so update the matrix
                    // entries for the columns since the last with a row singleton
                    shift_cols(a, &mut nnz, from_col, i_col0);
                    let from_os = a.start_[ic0];
                    let col_nnz = a.start_[ic0 + 1] - from_os;
                    let new_col_start = nnz;
                    let mut found_row_singleton = false;
                    for i_el in 0..col_nnz {
                        let src = (from_os + i_el) as usize;
                        let i_row = a.index_[src];
                        let new_row = new_row_index[i_row as usize];
                        if new_row >= 0 {
                            a.index_[nnz as usize] = new_row;
                            a.value_[nnz as usize] = a.value_[src];
                            nnz += 1;
                        } else {
                            debug_assert!(row_count[i_row as usize] == 1);
                            debug_assert!(col_of_row[i_row as usize] == i_col0);
                            debug_assert!(val_of_row[i_row as usize] == a.value_[src]);
                            found_row_singleton = true;
                        }
                    }
                    debug_assert!(found_row_singleton);
                    a.start_[ic0] = new_col_start;
                    from_col = i_col0 + 1;
                }
                // Update the matrix entries for the columns since the last with a
                // row singleton
                shift_cols(a, &mut nnz, from_col, num_col);
                a.start_[ncu] = nnz;
            } else {
                // Also removing redundant rows, so make the matrix rowwise and
                // remove rows simply above
                nnz = 0;
                let mut from_row: usize = 0;
                num_row = 0;
                // Lambda for shifting row data and updating row indices
                let shift_rows = |a: &mut HighsSparseMatrix,
                                  nnz: &mut HighsInt,
                                  num_row: &mut HighsInt,
                                  from_row: usize,
                                  to_row: usize| {
                    for i_row in from_row..to_row {
                        let new_row_start = *nnz;
                        for i_el in a.start_[i_row]..a.start_[i_row + 1] {
                            a.index_[*nnz as usize] = a.index_[i_el as usize];
                            a.value_[*nnz as usize] = a.value_[i_el as usize];
                            *nnz += 1;
                        }
                        a.start_[*num_row as usize] = new_row_start;
                        *num_row += 1;
                    }
                };
                a.ensure_rowwise();
                for i_row0 in 0..lp_num_row {
                    if new_row_index[i_row0] >= 0 {
                        continue;
                    }
                    // Row iRow0 is removed, so update the matrix entries for the
                    // rows since the last removed
                    shift_rows(a, &mut nnz, &mut num_row, from_row, i_row0);
                    from_row = i_row0 + 1;
                }
                // Update the matrix entries for the rows since the last removed
                shift_rows(a, &mut nnz, &mut num_row, from_row, lp_num_row);
                debug_assert!(num_row + num_removed_rows == lp_num_row as HighsInt);
                a.start_[num_row as usize] = nnz;
                a.num_row_ = num_row;
                a.ensure_colwise();
            }
            let nru = num_row as usize;
            self.lp.row_lower_.resize(nru, 0.0);
            self.lp.row_upper_.resize(nru, 0.0);
            if have_row_names {
                self.lp.row_names_.resize(nru, String::new());
            }
            self.lp.num_row_ = num_row;
            self.lp.a_matrix_.num_row_ = num_row;
            self.lp.a_matrix_.index_.resize(nnz as usize, 0);
            self.lp.a_matrix_.value_.resize(nnz as usize, 0.0);
            // compressRowIndexMap (HighsPostsolveStack.cpp lines 46-50), with origRowType.
            postsolve_stack.compress_index_maps(&new_row_index, &[]);
        }
        // Add doubleton equations, column singletons, variable locks

        if num_fixed_col != 0 || num_empty_col != 0 {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Initial sweep removes {} + {} = {} / {} empty + fixed columns\n",
                num_empty_col,
                num_fixed_col,
                num_removed_cols,
                original_num_col
            );
        }
        if num_empty_row != 0 || num_singleton_row != 0 || num_redundant_row != 0 {
            crate::highs_log_user!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Initial sweep identifies {} + {} + {} = {} / {} empty + singleton + redundant rows\n",
                num_empty_row,
                num_singleton_row,
                num_redundant_row,
                num_removed_rows,
                original_num_row
            );
        }

        H16SweepResult::kOk
    }
}

impl HPresolve {
    /// HiGHS 1.16: HPresolve.cpp::presolve, the initial-sweep block (lines 3994-4019).
    /// The presolve timer clocks are not ported.
    pub(super) fn h16_initial_sweep(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        if self.options.presolve != kHighsOffString
            && self.mipsolver.is_none()
            && self.options.presolve_rule_test == 0
            && self.analysis_.allow_rule_[kPresolveRuleInitialSweep as usize]
        {
            self.num_deleted_cols = 0;
            self.num_deleted_rows = 0;
            let logging_on = self.analysis_.logging_on_;
            if logging_on {
                self.analysis_.start_presolve_rule_log(
                    kPresolveRuleInitialSweep,
                    self.num_deleted_rows,
                    self.num_deleted_cols,
                    &self.model,
                );
            }
            let (sweep_result, del_cols, del_rows) = {
                let mut sweep = H16InitialSweep::new(
                    &mut self.model,
                    &self.options,
                    &self.analysis_.allow_rule_,
                    self.primal_feastol,
                );
                let r = sweep.run(postsolve_stack);
                (r, sweep.num_deleted_cols, sweep.num_deleted_rows)
            };
            self.num_deleted_cols = del_cols;
            self.num_deleted_rows = del_rows;
            self.analysis_.logging_on_ = logging_on;
            if logging_on {
                self.analysis_.stop_presolve_rule_log(
                    kPresolveRuleInitialSweep,
                    self.num_deleted_rows,
                    self.num_deleted_cols,
                    &self.model,
                );
            }
            if sweep_result == H16SweepResult::kPrimalInfeasible {
                return HPresolveResult::kPrimalInfeasible;
            }
            if sweep_result == H16SweepResult::kDualInfeasible {
                return HPresolveResult::kDualInfeasible;
            }
            // Check that time and reduction limits have not been reached
            //
            // NB Setting presolve_reduction_limit = 0 ensures that presolve
            // returns after initial sweep
            hpresolve_checked_call!(self.check_limits(postsolve_stack));
        }
        HPresolveResult::kOk
    }
}

#[cfg(test)]
mod h16_initial_sweep_tests {
    use super::*;

    /// 4 rows x 4 columns:
    /// * x0 empty (cost 1, lower 2): fixed at 2, offset += 2;
    /// * x1 fixed at 3, in rows 0 and 1: row sides shift by 3 and 6;
    /// * x2 in rows 0 and 1, x3 in rows 0, 1 and 3;
    /// * row 2 is empty, row 3 holds only x3 with 2 x3 <= 4 (singleton: x3 <= 2).
    fn lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 4;
        lp.col_cost_ = vec![1.0, 1.0, 1.0, 1.0];
        lp.col_lower_ = vec![2.0, 3.0, 0.0, 0.0];
        lp.col_upper_ = vec![5.0, 3.0, 10.0, 10.0];
        lp.integrality_ = vec![HighsVarType::kContinuous; 4];
        lp.row_lower_ = vec![-kHighsInf, 1.0, -1.0, -kHighsInf];
        lp.row_upper_ = vec![10.0, 20.0, 1.0, 4.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 0, 2, 4, 7];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1, 3];
        lp.a_matrix_.value_ = vec![1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0];
        lp
    }

    #[test]
    fn h16_initial_sweep_removes_cols_and_rows() {
        let mut lp = lp();
        let options = HighsOptions::default();
        let allow = vec![true; kH16PresolveRuleCount as usize];
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(4, 4);
        let mut sweep = H16InitialSweep::new(&mut lp, &options, &allow, 1e-7);
        assert_eq!(sweep.run(&mut stack), H16SweepResult::kOk);
        assert_eq!((sweep.num_deleted_cols, sweep.num_deleted_rows), (2, 2));
        // empty x0 at its lower bound 2 (cost 1) and fixed x1 at 3 (cost 1)
        assert_eq!(lp.offset_, 5.0);
        assert_eq!(lp.num_col_, 2);
        assert_eq!(lp.num_row_, 2);
        // rows 0 and 1 remain with sides shifted by x1 = 3 (coefficients 1 and 2)
        assert_eq!(lp.row_upper_, vec![7.0, 14.0]);
        assert_eq!(lp.row_lower_, vec![-kHighsInf, -5.0]);
        // the singleton row 3 tightened x3 <= 2
        assert_eq!(lp.col_upper_, vec![10.0, 2.0]);
        assert_eq!(lp.a_matrix_.start_, vec![0, 2, 4]);
        assert_eq!(lp.a_matrix_.index_, vec![0, 1, 0, 1]);
        assert_eq!(lp.a_matrix_.value_, vec![1.0; 4]);
        assert_eq!(stack.get_orig_cols_index(), &[2, 3]);
        assert_eq!(stack.get_orig_rows_index(), &[0, 1]);
        // empty col, fixed col, empty row, singleton row
        assert_eq!(stack.num_reductions(), 4);
    }

    /// A redundant row (x2 + x3 <= 100) sends the row removal through the rowwise copy.
    #[test]
    fn h16_initial_sweep_redundant_row_rowwise() {
        let mut lp = lp();
        lp.row_upper_[0] = 100.0;
        let options = HighsOptions::default();
        let allow = vec![true; kH16PresolveRuleCount as usize];
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(4, 4);
        let mut sweep = H16InitialSweep::new(&mut lp, &options, &allow, 1e-7);
        assert_eq!(sweep.run(&mut stack), H16SweepResult::kOk);
        assert_eq!((sweep.num_deleted_cols, sweep.num_deleted_rows), (2, 3));
        assert_eq!((lp.num_col_, lp.num_row_), (2, 1));
        assert_eq!(
            (lp.row_lower_.clone(), lp.row_upper_.clone()),
            (vec![-5.0], vec![14.0])
        );
        assert!(lp.a_matrix_.is_colwise());
        assert_eq!(lp.a_matrix_.start_, vec![0, 1, 2]);
        assert_eq!(lp.a_matrix_.index_, vec![0, 0]);
        assert_eq!(stack.get_orig_rows_index(), &[1]);
        assert_eq!(stack.num_reductions(), 5);
    }

    /// computeSingletonRowBounds: integral rounding and the redundant case.
    #[test]
    fn h16_singleton_row_bounds() {
        let (r, lb, ub, lt, ut) =
            h16_compute_singleton_row_bounds(2.0, 1.0, 9.0, 0.0, 10.0, 1e-6, 2.0, true);
        assert_eq!(r, H16SingletonRowResult::kBoundsTightened);
        assert_eq!((lb, ub, lt, ut), (1.0, 4.0, true, true));
        let (r, ..) =
            h16_compute_singleton_row_bounds(1.0, -1.0, 20.0, 0.0, 10.0, 1e-6, 1.0, false);
        assert_eq!(r, H16SingletonRowResult::kRedundant);
        let (r, ..) =
            h16_compute_singleton_row_bounds(1.0, 11.0, 20.0, 0.0, 10.0, 1e-6, 1.0, false);
        assert_eq!(r, H16SingletonRowResult::kPrimalInfeasible);
    }
}
