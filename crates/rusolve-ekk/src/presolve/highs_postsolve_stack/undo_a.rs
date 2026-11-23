//! HiGHS source: presolve/HighsPostsolveStack.cpp (lines 84-127, 149-576; 1.16.0 lines 429-438)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `undo` of FreeColSubstitution, DoubletonEquation, EqualityRowAddition(s),
//! ForcingColumn, ForcingColumnRemovedRow, SingletonRow, FixedCol,
//! RedundantRow, ImpliedEquation (1.16, #3218), ForcingRow and DuplicateRow. The file-static helpers
//! `computeRowStatus` / `computeStatus` are `super::compute_row_status`,
//! `super::compute_status_inout` (the `HighsBasisStatus&` overload) and
//! `super::compute_status` (the 2-argument overload). Fields are
//! snake_case (`colCost` -> `col_cost`, `rowType` -> `row_type`).

use super::*;
use crate::prelude::*;

impl FreeColSubstitution {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::FreeColSubstitution::undo (lines 84-127, 44 lines)
    pub fn undo(
        &mut self,
        options: &HighsOptions,
        row_values: &[Nonzero],
        col_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // compute primal values
        let mut col_coef: f64 = 0.0;
        let mut row_value = HighsCDouble::from(0.0);
        for row_val in row_values {
            if row_val.index == self.col {
                col_coef = row_val.value;
            } else {
                row_value += row_val.value * solution.col_value[row_val.index as usize];
            }
        }

        debug_assert!(col_coef != 0.0);
        // Row values aren't fully postsolved, so why do this?
        if (self.row as usize) < solution.row_value.len() {
            solution.row_value[self.row as usize] =
                f64::from(row_value + col_coef * solution.col_value[self.col as usize]);
        }
        solution.col_value[self.col as usize] = f64::from((-row_value + self.rhs) / col_coef);

        // if no dual values requested, return here
        if !solution.dual_valid {
            return;
        }

        // compute the row dual value such that reduced cost of basic column is 0
        if (self.row as usize) < solution.row_value.len() {
            solution.row_dual[self.row as usize] = 0.0;
            let mut dualval = HighsCDouble::from(self.col_cost);
            for col_val in col_values {
                if (col_val.index as usize) < solution.row_value.len() {
                    dualval -= col_val.value * solution.row_dual[col_val.index as usize];
                }
            }
            solution.row_dual[self.row as usize] = f64::from(dualval / col_coef);
        }

        solution.col_dual[self.col as usize] = 0.0;

        // set basis status if necessary
        if !basis.valid {
            return;
        }

        basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
        if (self.row as usize) < solution.row_value.len() {
            basis.row_status[self.row as usize] =
                compute_row_status(solution.row_dual[self.row as usize], self.row_type);
        }
    }
}

impl DoubletonEquation {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::DoubletonEquation::undo (lines 149-229, 81 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        col_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // retrieve the row and column index, the row side and the two
        // coefficients then compute the primal values
        solution.col_value[self.col_subst as usize] = f64::from(
            (-(HighsCDouble::from(self.coef) * solution.col_value[self.col as usize]) + self.rhs)
                / self.coef_subst,
        );

        // can only do primal postsolve, stop here
        if self.row == -1 || !solution.dual_valid {
            return;
        }

        let col_status: HighsBasisStatus = if !basis.valid {
            compute_status(
                solution.col_dual[self.col as usize],
                options.dual_feasibility_tolerance,
            )
        } else {
            compute_status_inout(
                solution.col_dual[self.col as usize],
                &mut basis.col_status[self.col as usize],
                options.dual_feasibility_tolerance,
            )
        };

        // compute the current dual values of the row and the substituted column
        // before deciding on which column becomes basic
        // for each entry in a row i of the substituted column we added the
        // doubleton equation row with scale -a_i/substCoef. Therefore the dual
        // multiplier of this row i implicitly increases the dual multiplier of this
        // doubleton equation row with that scale.
        let mut row_dual = HighsCDouble::from(0.0);
        if (self.row as usize) < solution.row_value.len() {
            solution.row_dual[self.row as usize] = 0.0;
            for col_val in col_values {
                if (col_val.index as usize) < solution.row_value.len() {
                    row_dual -= col_val.value * solution.row_dual[col_val.index as usize];
                }
            }
            row_dual /= self.coef_subst;
            solution.row_dual[self.row as usize] = f64::from(row_dual);
        }
        // the equation was also added to the objective, so the current values need
        // to be changed
        solution.col_dual[self.col_subst as usize] = self.subst_cost;
        solution.col_dual[self.col as usize] += self.subst_cost * self.coef / self.coef_subst;

        if (self.upper_tightened && col_status == HighsBasisStatus::kUpper)
            || (self.lower_tightened && col_status == HighsBasisStatus::kLower)
        {
            // column must get zero reduced cost as the current bound cannot be used
            // so alter the dual multiplier of the row to make the dual multiplier of
            // column zero
            let row_dual_delta: f64 = solution.col_dual[self.col as usize] / self.coef;
            if (self.row as usize) < solution.row_value.len() {
                solution.row_dual[self.row as usize] = f64::from(row_dual + row_dual_delta);
            }
            solution.col_dual[self.col as usize] = 0.0;
            solution.col_dual[self.col_subst as usize] = f64::from(
                HighsCDouble::from(solution.col_dual[self.col_subst as usize])
                    - row_dual_delta * self.coef_subst,
            );

            if basis.valid {
                if (self.coef.is_sign_negative() == self.coef_subst.is_sign_negative()
                    && basis.col_status[self.col as usize] == HighsBasisStatus::kUpper)
                    || (self.coef.is_sign_negative() != self.coef_subst.is_sign_negative()
                        && basis.col_status[self.col as usize] == HighsBasisStatus::kLower)
                {
                    basis.col_status[self.col_subst as usize] = HighsBasisStatus::kLower;
                } else {
                    basis.col_status[self.col_subst as usize] = HighsBasisStatus::kUpper;
                }
                basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
            }
        } else {
            // otherwise make the reduced cost of the substituted column zero and make
            // that column basic
            let row_dual_delta: f64 = solution.col_dual[self.col_subst as usize] / self.coef_subst;
            if (self.row as usize) < solution.row_value.len() {
                solution.row_dual[self.row as usize] = f64::from(row_dual + row_dual_delta);
            }
            solution.col_dual[self.col_subst as usize] = 0.0;
            solution.col_dual[self.col as usize] = f64::from(
                HighsCDouble::from(solution.col_dual[self.col as usize])
                    - row_dual_delta * self.coef,
            );
            if basis.valid {
                basis.col_status[self.col_subst as usize] = HighsBasisStatus::kBasic;
            }
        }

        if !basis.valid {
            return;
        }

        if (self.row as usize) < solution.row_value.len() {
            basis.row_status[self.row as usize] =
                compute_row_status(solution.row_dual[self.row as usize], self.row_type);
        }
    }
}

impl EqualityRowAddition {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::EqualityRowAddition::undo (lines 231-248, 18 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        eq_row_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // (removed) cuts may have been used in this reduction.
        if (self.row as usize) >= solution.row_value.len()
            || (self.added_eq_row as usize) >= solution.row_value.len()
        {
            return;
        }

        // nothing more to do if the row is zero in the dual solution or there is
        // no dual solution
        if !solution.dual_valid || solution.row_dual[self.row as usize] == 0.0 {
            return;
        }

        // the dual multiplier of the row implicitly increases the dual multiplier
        // of the equation with the scale the equation was added with
        solution.row_dual[self.added_eq_row as usize] = f64::from(
            HighsCDouble::from(self.eq_row_scale) * solution.row_dual[self.row as usize]
                + solution.row_dual[self.added_eq_row as usize],
        );

        debug_assert!(!basis.valid);
    }
}

impl EqualityRowAdditions {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::EqualityRowAdditions::undo (lines 250-273, 24 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        eq_row_values: &[Nonzero],
        target_rows: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // a (removed) cut may have been used in this reduction.
        if (self.added_eq_row as usize) >= solution.row_value.len() {
            return;
        }

        // nothing more to do if the row is zero in the dual solution or there is
        // no dual solution
        if !solution.dual_valid {
            return;
        }

        // the dual multiplier of the rows where the eq row was added implicitly
        // increases the dual multiplier of the equation with the scale that was
        // used for adding the equation
        let mut eq_row_dual = HighsCDouble::from(solution.row_dual[self.added_eq_row as usize]);
        for target_row in target_rows {
            if (target_row.index as usize) < solution.row_value.len() {
                eq_row_dual += HighsCDouble::from(target_row.value)
                    * solution.row_dual[target_row.index as usize];
            }
        }
        solution.row_dual[self.added_eq_row as usize] = f64::from(eq_row_dual);

        debug_assert!(!basis.valid);
    }
}

impl ForcingColumn {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::ForcingColumn::undo (lines 275-337, 63 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        col_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        let mut nonbasic_row: HighsInt = -1;
        let mut nonbasic_row_status = HighsBasisStatus::kNonbasic;
        let mut col_val_from_nonbasic_row: f64 = self.col_bound;

        let mut debug_num_use_row_value: HighsInt = 0;
        let debug_report = false;

        let mut compute_col_val = |direction: HighsInt| {
            // choose solution value that makes all rows feasible
            for col_val in col_values {
                // Row values aren't fully postsolved, so how can this work?
                debug_num_use_row_value += 1;
                if (col_val.index as usize) < solution.row_value.len() {
                    let col_val_from_row =
                        solution.row_value[col_val.index as usize] / col_val.value;
                    if (direction as f64) * col_val_from_row
                        > (direction as f64) * col_val_from_nonbasic_row
                    {
                        nonbasic_row = col_val.index;
                        col_val_from_nonbasic_row = col_val_from_row;
                        nonbasic_row_status = if (direction as f64) * col_val.value > 0.0 {
                            HighsBasisStatus::kLower
                        } else {
                            HighsBasisStatus::kUpper
                        };
                    }
                }
            }
            // round solution value if column is integer-constrained
            if nonbasic_row != -1 && self.col_integral {
                col_val_from_nonbasic_row = (direction as f64)
                    * ((direction as f64) * col_val_from_nonbasic_row
                        - options.mip_feasibility_tolerance)
                        .ceil();
            }
        };

        if self.at_infinite_upper {
            // choose largest value as then all rows are feasible
            compute_col_val(1);
        } else {
            // choose smallest value, as then all rows are feasible
            compute_col_val(-1);
        }
        if debug_num_use_row_value != 0 && debug_report {
            println!(
                "HighsPostsolveStack::ForcingColumn::undo Using {} unknown row activit{}",
                debug_num_use_row_value,
                if debug_num_use_row_value > 1 {
                    "ies"
                } else {
                    "y"
                }
            );
        }

        solution.col_value[self.col as usize] = col_val_from_nonbasic_row;

        if !solution.dual_valid {
            return;
        }

        solution.col_dual[self.col as usize] = 0.0;

        if !basis.valid {
            return;
        }
        if nonbasic_row == -1 {
            basis.col_status[self.col as usize] = if self.at_infinite_upper {
                HighsBasisStatus::kLower
            } else {
                HighsBasisStatus::kUpper
            };
        } else {
            basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
            basis.row_status[nonbasic_row as usize] = nonbasic_row_status;
        }
    }
}

impl ForcingColumnRemovedRow {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::ForcingColumnRemovedRow::undo (lines 339-356, 18 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        row_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // a (removed) cut may have been used in this reduction.
        if (self.row as usize) >= solution.row_value.len() {
            return;
        }

        // we use the row value as storage for the scaled value implied on the
        // column dual
        let mut val = HighsCDouble::from(self.rhs);
        for row_val in row_values {
            val -= row_val.value * solution.col_value[row_val.index as usize];
        }

        // Row values aren't fully postsolved, so why do this?
        solution.row_value[self.row as usize] = f64::from(val);

        if solution.dual_valid {
            solution.row_dual[self.row as usize] = 0.0;
        }
        if basis.valid {
            basis.row_status[self.row as usize] = HighsBasisStatus::kBasic;
        }
    }
}

impl SingletonRow {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::SingletonRow::undo (lines 358-418, 61 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // nothing to do if the rows dual value is zero in the dual solution or
        // there is no dual solution
        if !solution.dual_valid {
            return;
        }

        let col_status: HighsBasisStatus = if !basis.valid {
            compute_status(
                solution.col_dual[self.col as usize],
                options.dual_feasibility_tolerance,
            )
        } else {
            compute_status_inout(
                solution.col_dual[self.col as usize],
                &mut basis.col_status[self.col as usize],
                options.dual_feasibility_tolerance,
            )
        };

        if (!self.col_lower_tightened || col_status != HighsBasisStatus::kLower)
            && (!self.col_upper_tightened || col_status != HighsBasisStatus::kUpper)
        {
            // the tightened bound is not used in the basic solution
            // hence we simply make the row basic and give it a dual multiplier of 0
            if (self.row as usize) < solution.row_value.len() {
                if basis.valid {
                    basis.row_status[self.row as usize] = HighsBasisStatus::kBasic;
                }
                solution.row_dual[self.row as usize] = 0.0;
            }
            return;
        }

        // choose the row dual value such that the columns reduced cost becomes
        // zero
        if (self.row as usize) < solution.row_value.len() {
            solution.row_dual[self.row as usize] = solution.col_dual[self.col as usize] / self.coef;
        }
        solution.col_dual[self.col as usize] = 0.0;

        if !basis.valid {
            return;
        }

        if (self.row as usize) < solution.row_value.len() {
            match col_status {
                HighsBasisStatus::kLower => {
                    debug_assert!(self.col_lower_tightened);
                    if self.coef > 0.0 {
                        // tightened lower bound comes from row lower bound
                        basis.row_status[self.row as usize] = HighsBasisStatus::kLower;
                    } else {
                        // tightened lower bound comes from row upper bound
                        basis.row_status[self.row as usize] = HighsBasisStatus::kUpper;
                    }
                }
                HighsBasisStatus::kUpper => {
                    if self.coef > 0.0 {
                        // tightened upper bound comes from row upper bound
                        basis.row_status[self.row as usize] = HighsBasisStatus::kUpper;
                    } else {
                        // tightened upper bound comes from row lower bound
                        basis.row_status[self.row as usize] = HighsBasisStatus::kLower;
                    }
                }
                _ => debug_assert!(false),
            }
        }

        // column becomes basic
        basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
    }
}

impl FixedCol {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::FixedCol::undo (lines 421-448, 28 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        col_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // set solution value
        solution.col_value[self.col as usize] = self.fix_value;

        if !solution.dual_valid {
            return;
        }

        // compute reduced cost

        let mut reduced_cost = HighsCDouble::from(self.col_cost);
        for col_val in col_values {
            if (col_val.index as usize) < solution.row_value.len() {
                reduced_cost -= col_val.value * solution.row_dual[col_val.index as usize];
            }
        }

        solution.col_dual[self.col as usize] = f64::from(reduced_cost);

        // set basis status
        if basis.valid {
            basis.col_status[self.col as usize] = self.fix_type;
            if basis.col_status[self.col as usize] == HighsBasisStatus::kNonbasic {
                basis.col_status[self.col as usize] = if solution.col_dual[self.col as usize] >= 0.0
                {
                    HighsBasisStatus::kLower
                } else {
                    HighsBasisStatus::kUpper
                };
            }
        }
    }
}

impl RedundantRow {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::RedundantRow::undo (lines 450-462, 13 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // a (removed) cut may have been used in this reduction.
        if (self.row as usize) >= solution.row_value.len() {
            return;
        }

        // set row dual to zero if dual solution requested
        if !solution.dual_valid {
            return;
        }

        solution.row_dual[self.row as usize] = 0.0;

        if basis.valid {
            basis.row_status[self.row as usize] = HighsBasisStatus::kBasic;
        }
    }
}

impl ImpliedEquation {
    /// HiGHS 1.16 (#3218): HighsPostsolveStack.cpp::ImpliedEquation::undo (lines 429-438)
    ///
    /// rowPresolve turned an inequality into an equation because complementary
    /// slackness implied it; a postsolved dual with the wrong sign for the original
    /// side is moved back onto the columns.
    pub fn undo(&self, row_values: &[Nonzero], solution: &mut HighsSolution) {
        if !solution.dual_valid {
            return;
        }
        let old_dual = solution.row_dual[self.row as usize];
        if if self.at_lower {
            old_dual < 0.0
        } else {
            old_dual > 0.0
        } {
            solution.row_dual[self.row as usize] = 0.0;
            for nz in row_values {
                solution.col_dual[nz.index as usize] += nz.value * old_dual;
            }
        }
    }
}

impl ForcingRow {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::ForcingRow::undo (lines 464-503, 40 lines)
    pub fn undo(
        &self,
        options: &HighsOptions,
        row_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        if !solution.dual_valid {
            return;
        }

        // compute the row dual multiplier and determine the new basic column
        let mut basic_col: HighsInt = -1;
        let mut dual_delta: f64 = 0.0;
        let direction: HighsInt = if self.row_type == RowType::kLeq {
            1
        } else {
            -1
        };
        for row_val in row_values {
            let col_dual = solution.col_dual[row_val.index as usize] - row_val.value * dual_delta;
            if (direction as f64) * col_dual * row_val.value < 0.0 {
                // column is dual infeasible, decrease the row dual such that its
                // reduced cost become zero and remember this column as the new basic
                // column for this row
                dual_delta = solution.col_dual[row_val.index as usize] / row_val.value;
                basic_col = row_val.index;
            }
        }

        if basic_col != -1 {
            if (self.row as usize) < solution.row_value.len() {
                solution.row_dual[self.row as usize] =
                    solution.row_dual[self.row as usize] + dual_delta;
            }
            for row_val in row_values {
                solution.col_dual[row_val.index as usize] = f64::from(
                    -(HighsCDouble::from(dual_delta) * row_val.value)
                        + solution.col_dual[row_val.index as usize],
                );
            }
            solution.col_dual[basic_col as usize] = 0.0;

            if basis.valid {
                if (self.row as usize) < solution.row_value.len() {
                    basis.row_status[self.row as usize] = if self.row_type == RowType::kGeq {
                        HighsBasisStatus::kLower
                    } else {
                        HighsBasisStatus::kUpper
                    };
                }

                basis.col_status[basic_col as usize] = HighsBasisStatus::kBasic;
            }
        }
    }
}

impl DuplicateRow {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::DuplicateRow::undo (lines 505-576, 72 lines)
    /// HiGHS 1.16 (#3313): HighsPostsolveStack.cpp::DuplicateRow::undo (lines 479-540)
    pub fn undo(
        &self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        // (removed) cuts may have been used in this reduction.
        if (self.row as usize) >= solution.row_value.len() {
            return;
        }

        if !solution.dual_valid {
            return;
        }
        if !self.row_upper_tightened && !self.row_lower_tightened {
            // simple case of row2 being redundant, in which case it just gets a
            // dual multiplier of 0 and is made basic
            if (self.duplicate_row as usize) < solution.row_value.len() {
                solution.row_dual[self.duplicate_row as usize] = 0.0;
                if basis.valid {
                    basis.row_status[self.duplicate_row as usize] = HighsBasisStatus::kBasic;
                }
            }
            return;
        }

        // HiGHS 1.16 #3313 (h16_dup_row_undo): HighsPostsolveStack.cpp lines 498-536 (1.16.0)
        let h16_dup_row_undo = options.h16_dup_row_undo;

        let row_status: HighsBasisStatus = if !basis.valid {
            compute_status(
                solution.row_dual[self.row as usize],
                options.dual_feasibility_tolerance,
            )
        } else {
            compute_status_inout(
                solution.row_dual[self.row as usize],
                &mut basis.row_status[self.row as usize],
                options.dual_feasibility_tolerance,
            )
        };

        let compute_row_dual_and_status =
            |solution: &mut HighsSolution,
             basis: &mut HighsBasis,
             tighened: bool,
             dir: HighsInt| {
                if tighened {
                    if (self.duplicate_row as usize) < solution.row_value.len() {
                        solution.row_dual[self.duplicate_row as usize] = if h16_dup_row_undo {
                            solution.row_dual[self.row as usize] * self.duplicate_row_scale
                        } else {
                            solution.row_dual[self.row as usize] / self.duplicate_row_scale
                        };
                        if basis.valid {
                            let upper = if h16_dup_row_undo {
                                (dir as f64) * self.duplicate_row_scale > 0.0
                            } else {
                                self.duplicate_row_scale > 0.0
                            };
                            if upper {
                                basis.row_status[self.duplicate_row as usize] =
                                    HighsBasisStatus::kUpper;
                            } else {
                                basis.row_status[self.duplicate_row as usize] =
                                    HighsBasisStatus::kLower;
                            }
                        }
                    }
                    solution.row_dual[self.row as usize] = 0.0;
                    if basis.valid {
                        basis.row_status[self.row as usize] = HighsBasisStatus::kBasic;
                    }
                } else if (self.duplicate_row as usize) < solution.row_value.len() {
                    solution.row_dual[self.duplicate_row as usize] = 0.0;
                    if basis.valid {
                        basis.row_status[self.duplicate_row as usize] = HighsBasisStatus::kBasic;
                    }
                }
            };

        // at least one bound of the row was tightened by using the bound of the
        // scaled parallel row, hence we might need to make the parallel row
        // nonbasic and the row basic

        match row_status {
            HighsBasisStatus::kBasic => {
                // if row is basic the parallel row is also basic
                if (self.duplicate_row as usize) < solution.row_value.len() {
                    solution.row_dual[self.duplicate_row as usize] = 0.0;
                    if basis.valid {
                        basis.row_status[self.duplicate_row as usize] = HighsBasisStatus::kBasic;
                    }
                }
            }
            HighsBasisStatus::kUpper => {
                // if row sits on its upper bound, and the row upper bound was
                // tightened using the parallel row we make the row basic and
                // transfer its dual value to the parallel row with the proper scale
                compute_row_dual_and_status(solution, basis, self.row_upper_tightened, 1);
            }
            HighsBasisStatus::kLower => {
                compute_row_dual_and_status(solution, basis, self.row_lower_tightened, -1);
            }
            _ => debug_assert!(false),
        }
    }
}

#[cfg(test)]
mod h16_dup_row_undo_tests {
    use super::*;
    use crate::prelude::*;

    /// The LP of 1.16.0 check/TestPresolve.cpp `test-duplicate-row` (#3313):
    /// rows r0 and r1 are parallel, with r1 = -4 * r0 in the first matrix.
    /// `parallel` selects the row bounds of the first pass of that test;
    /// `values` is the column-wise matrix (x, y, z per row, see the C++).
    fn dup_row_lp(parallel: bool, values: [f64; 6]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        lp.col_cost_ = vec![2.0, -1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, 2.0, kHighsInf];
        lp.a_matrix_.start_ = vec![0, 2, 4, 6];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1];
        lp.a_matrix_.value_ = values.to_vec();
        if parallel {
            lp.row_lower_ = vec![-kHighsInf, -4.0];
            lp.row_upper_ = vec![0.0, 0.0];
        } else {
            lp.row_lower_ = vec![0.0, 0.0];
            lp.row_upper_ = vec![kHighsInf, 4.0];
        }
        lp
    }

    /// The two matrix value sets of the C++ test (the second is the sign flip).
    const VALUES_A: [f64; 6] = [1.0, -4.0, -1.0, 4.0, 3.0, -12.0];
    const VALUES_B: [f64; 6] = [-1.0, 4.0, 1.0, -4.0, -3.0, 12.0];

    fn solve(lp: HighsLp, presolve: bool, flag: bool) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_dup_row_undo", flag),
            HighsStatus::kOk
        );
        let presolve_value = if presolve { "on" } else { "off" };
        assert_eq!(
            highs.set_option_value_str("presolve", presolve_value),
            HighsStatus::kOk
        );
        // As in the C++ test: IPX with crossover, which gives the basis
        // that exposes the wrong nonbasic bound.
        assert_eq!(
            highs.set_option_value_str("solver", "ipx"),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("run_crossover", "on"),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        highs
    }

    /// Coefficient a(i, j) of the column-wise matrix.
    fn coef(lp: &HighsLp, i: HighsInt, j: usize) -> f64 {
        let m = &lp.a_matrix_;
        for k in m.start_[j] as usize..m.start_[j + 1] as usize {
            if m.index_[k] == i {
                return m.value_[k];
            }
        }
        0.0
    }

    /// Checks the returned solution against the original LP: optimality,
    /// row activities within their bounds, the reduced-cost identity
    /// `col_dual = c - A' row_dual`, and that the basis statuses agree with
    /// the primal values and duals. Returns the first violation found.
    fn kkt_violation(lp: &HighsLp, highs: &Highs) -> Option<String> {
        let tol = 1e-7;
        if highs.get_model_status() != HighsModelStatus::kOptimal {
            return Some("model status is not optimal".to_string());
        }
        let sol = highs.get_solution();
        let basis = highs.get_basis();
        if !sol.dual_valid || !basis.valid {
            return Some("no valid duals or basis".to_string());
        }
        let x = &sol.col_value;
        let y = &sol.row_dual;
        for j in 0..3 {
            let d = lp.col_cost_[j] - (0..2).map(|i| coef(lp, i, j) * y[i as usize]).sum::<f64>();
            if (d - sol.col_dual[j]).abs() > tol {
                return Some(format!(
                    "reduced cost of col {j}: {d} vs {}",
                    sol.col_dual[j]
                ));
            }
            match basis.col_status[j] {
                HighsBasisStatus::kBasic => {
                    if sol.col_dual[j].abs() > tol {
                        return Some(format!("basic col {j} has dual {}", sol.col_dual[j]));
                    }
                }
                HighsBasisStatus::kLower => {
                    if (x[j] - lp.col_lower_[j]).abs() > tol {
                        return Some(format!("col {j} at lower status but value {}", x[j]));
                    }
                }
                HighsBasisStatus::kUpper if (x[j] - lp.col_upper_[j]).abs() > tol => {
                    return Some(format!("col {j} at upper status but value {}", x[j]));
                }
                _ => {}
            }
        }
        for i in 0..2 {
            let act: f64 = (0..3).map(|j| coef(lp, i as HighsInt, j) * x[j]).sum();
            let (lo, up) = (lp.row_lower_[i], lp.row_upper_[i]);
            if act < lo - tol || act > up + tol {
                return Some(format!("row {i} activity {act} outside [{lo}, {up}]"));
            }
            match basis.row_status[i] {
                HighsBasisStatus::kBasic => {
                    if y[i].abs() > tol {
                        return Some(format!("basic row {i} has dual {}", y[i]));
                    }
                }
                HighsBasisStatus::kLower => {
                    if !lo.is_finite() || (act - lo).abs() > tol {
                        return Some(format!("row {i} at lower status but activity {act}"));
                    }
                }
                HighsBasisStatus::kUpper if (!up.is_finite() || (act - up).abs() > tol) => {
                    return Some(format!("row {i} at upper status but activity {act}"));
                }
                _ => {}
            }
        }
        None
    }

    fn objective(highs: &Highs) -> f64 {
        highs.get_info().objective_function_value
    }

    /// Presolve on with the flag on: the postsolved duals and basis pass the
    /// KKT check, and the objective matches the presolve-off solve.
    #[test]
    fn h16_dup_row_undo_on_kkt() {
        for parallel in [true, false] {
            for values in [VALUES_A, VALUES_B] {
                let lp = dup_row_lp(parallel, values);
                let reference = solve(lp.clone(), false, false);
                assert_eq!(kkt_violation(&lp, &reference), None);
                let on = solve(lp.clone(), true, true);
                assert_eq!(kkt_violation(&lp, &on), None, "parallel={parallel}");
                assert!((objective(&on) - objective(&reference)).abs() <= 1e-7);
            }
        }
    }

    /// Flag off: the same LPs solve to the same optimal objective as the
    /// presolve-off reference (1.15.1 path, unchanged by the new gate).
    #[test]
    fn h16_dup_row_undo_off_objective() {
        for parallel in [true, false] {
            for values in [VALUES_A, VALUES_B] {
                let lp = dup_row_lp(parallel, values);
                let reference = solve(lp.clone(), false, false);
                let off = solve(lp.clone(), true, false);
                assert_eq!(off.get_model_status(), HighsModelStatus::kOptimal);
                assert!((objective(&off) - objective(&reference)).abs() <= 1e-7);
            }
        }
    }

    /// Direct call of `DuplicateRow::undo` in the state the C++ comment
    /// describes: row 0 is nonbasic at its lower bound, its lower bound was
    /// tightened by the parallel row 1 with scale -0.25. The flag selects
    /// the 1.16 (#3313) or the 1.15.1 wrong-bound result.
    #[test]
    fn h16_dup_row_undo_direct_gate() {
        let reduction = DuplicateRow {
            duplicate_row_scale: -0.25,
            duplicate_row: 1,
            row: 0,
            row_lower_tightened: true,
            row_upper_tightened: false,
        };
        for flag in [false, true] {
            let mut options = HighsOptions::default();
            options.h16_dup_row_undo = flag;
            let mut solution = HighsSolution::default();
            solution.dual_valid = true;
            solution.row_value = vec![0.0, 0.0];
            solution.row_dual = vec![0.5, 0.0];
            let mut basis = HighsBasis::default();
            basis.valid = true;
            basis.row_status = vec![HighsBasisStatus::kBasic, HighsBasisStatus::kBasic];
            reduction.undo(&options, &mut solution, &mut basis);
            assert_eq!(solution.row_dual[0], 0.0);
            assert!(matches!(basis.row_status[0], HighsBasisStatus::kBasic));
            if flag {
                assert_eq!(solution.row_dual[1], 0.5 * -0.25);
                assert!(matches!(basis.row_status[1], HighsBasisStatus::kUpper));
            } else {
                assert_eq!(solution.row_dual[1], 0.5 / -0.25);
                assert!(matches!(basis.row_status[1], HighsBasisStatus::kLower));
            }
        }
    }
}
