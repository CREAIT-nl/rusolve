//! HiGHS source: presolve/HighsPostsolveStack.cpp (lines 578-1012)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `DuplicateColumn::undo` and `DuplicateColumn::okMerge`. Fields are
//! snake_case (`colScale` -> `col_scale`, `duplicateColLower` ->
//! `duplicate_col_lower`, ...).

use super::*;
use crate::prelude::*;

impl DuplicateColumn {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::DuplicateColumn::undo (lines 578-867, 290 lines)
    /// HiGHS 1.16 (#3285): HighsPostsolveStack.cpp::DuplicateColumn::undo (lines 675-720)
    /// HiGHS 1.16 (#3580): HighsPostsolveStack.cpp::DuplicateColumn::undo (lines 749-756)
    pub fn undo(
        &self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        let debug_report = false;
        // HiGHS 1.16 #3285 (h16_dominated_cols) and #3580 (h16_dupcol_warn): flags read once.
        let h16_dominated_cols = options.h16_dominated_cols;
        let h16_dupcol_warn = options.h16_dupcol_warn;
        let merge_val = solution.col_value[self.col as usize];

        let ok_residual = |x: f64, y: f64| -> bool {
            let check_merge_val = x + self.col_scale * y;
            let residual = (check_merge_val - merge_val).abs();
            let ok = residual <= options.primal_feasibility_tolerance;
            if !ok && debug_report {
                println!(
                    "HighsPostsolveStack::DuplicateColumn::undo {} + {}.{} = {} != {}: residual = {}",
                    x, self.col_scale, y, check_merge_val, merge_val, residual
                );
            }
            ok
        };

        let is_at_bound = |value: f64, bound: f64| -> bool {
            if value < bound - options.primal_feasibility_tolerance {
                return false;
            }
            if value <= bound + options.primal_feasibility_tolerance {
                return true;
            }
            false
        };

        // const bool ok_merge = okMerge(options.mip_feasibility_tolerance);
        // assert(ok_merge);
        //
        // the column dual of the duplicate column is easily computed by scaling
        // since col * colScale yields the coefficient values and cost of the
        // duplicate column.
        if solution.dual_valid {
            solution.col_dual[self.duplicate_col as usize] =
                solution.col_dual[self.col as usize] * self.col_scale;
        }

        if basis.valid {
            // do postsolve using basis status if a basis is available: if the
            // merged column is nonbasic, we can just set both columns to
            // appropriate nonbasic status and value
            //
            // Undoing z = x + a.y
            //
            // Since x became z, its basis status is unchanged
            //
            // For a > 0, z in [x_l + a.y_l, x_u + a.y_u]
            //
            // If z is nonbasic at its lower (upper) bound, set y to be
            // nonbasic at its lower (upper) bound
            //
            // For a < 0, z in [x_l + a.y_u, x_u + a.y_l]
            //
            // If z is nonbasic at lower (upper) bound, set y to be nonbasic
            // at its upper (lower) bounds
            //
            // Check for perturbations
            match basis.col_status[self.col as usize] {
                HighsBasisStatus::kLower => {
                    solution.col_value[self.col as usize] = self.col_lower;
                    if self.col_scale > 0.0 {
                        basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kLower;
                        solution.col_value[self.duplicate_col as usize] = self.duplicate_col_lower;
                    } else {
                        basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kUpper;
                        solution.col_value[self.duplicate_col as usize] = self.duplicate_col_upper;
                    }
                    // nothing else to do
                    debug_assert!(ok_residual(
                        solution.col_value[self.col as usize],
                        solution.col_value[self.duplicate_col as usize]
                    ));
                    return;
                }
                HighsBasisStatus::kUpper => {
                    solution.col_value[self.col as usize] = self.col_upper;
                    if self.col_scale > 0.0 {
                        basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kUpper;
                        solution.col_value[self.duplicate_col as usize] = self.duplicate_col_upper;
                    } else {
                        basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kLower;
                        solution.col_value[self.duplicate_col as usize] = self.duplicate_col_lower;
                    }
                    // nothing else to do
                    debug_assert!(ok_residual(
                        solution.col_value[self.col as usize],
                        solution.col_value[self.duplicate_col as usize]
                    ));
                    return;
                }
                HighsBasisStatus::kZero => {
                    solution.col_value[self.col as usize] = 0.0;
                    basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kZero;
                    solution.col_value[self.duplicate_col as usize] = 0.0;
                    // nothing else to do
                    debug_assert!(ok_residual(
                        solution.col_value[self.col as usize],
                        solution.col_value[self.duplicate_col as usize]
                    ));
                    return;
                }
                HighsBasisStatus::kBasic | HighsBasisStatus::kNonbasic => {}
            }
            // Nonbasic cases should have been considered; basic case
            // considered later
            debug_assert!(basis.col_status[self.col as usize] == HighsBasisStatus::kBasic);
        }

        // There is either no basis for postsolve, or column status is
        // basic. One of the two columns must become nonbasic. In case of
        // integrality it is simpler to choose col, since it has a
        // coefficient of +1 in the equation z = col + colScale *
        // duplicateCol where the merged column is z and is currently using
        // the index of col. The duplicateCol can have a positive or
        // negative coefficient. So for postsolve, we first start out with
        // col sitting at the lower bound and compute the corresponding
        // value for the duplicate column as (z - colLower)/colScale. Then
        // the following things might happen:
        //
        // - case 1: the value computed for duplicateCol is within the bounds
        //   - case 1.1: duplicateCol is continuous -> accept value, make col nonbasic
        //     at lower and duplicateCol basic
        //   - case 1.2: duplicateCol is integer -> accept value if integer feasible,
        //     otherwise round down and compute value of col as
        //     col = z - colScale * duplicateCol
        // - case 2: the value for duplicateCol violates the column bounds: make it
        //   sit at the bound that is violated and compute the value of col as col = z - colScale *
        //   duplicateCol for basis postsolve col is basic and duplicateCol
        //   nonbasic at lower/upper depending on which bound is violated.

        if self.col_lower != -kHighsInf {
            solution.col_value[self.col as usize] = self.col_lower;
        } else {
            solution.col_value[self.col as usize] = cpp_min(0.0, self.col_upper);
        }
        // Determine a prospective value for col_value[duplicateCol]
        solution.col_value[self.duplicate_col as usize] = f64::from(
            (HighsCDouble::from(merge_val) - solution.col_value[self.col as usize])
                / self.col_scale,
        );

        let mut recompute_col = false;
        // HiGHS 1.16 #3285 (h16_dominated_cols): HighsPostsolveStack.cpp lines 675, 686 (1.16.0)
        let mut duplicate_col_at_upper = false;

        // Set any basis status for duplicateCol to kNonbasic to check that
        // it is set
        if basis.valid {
            basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kNonbasic;
        }

        if solution.col_value[self.duplicate_col as usize] > self.duplicate_col_upper {
            // Prospective value exceeds the upper bound, so trim it to the
            // upper bound and force recalculation of col_value[col]
            solution.col_value[self.duplicate_col as usize] = self.duplicate_col_upper;
            recompute_col = true;
            if h16_dominated_cols {
                // HiGHS 1.16 #3285 (h16_dominated_cols): HighsPostsolveStack.cpp line 686 (1.16.0)
                duplicate_col_at_upper = true;
            }
            if basis.valid {
                basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kUpper;
            }
        } else if solution.col_value[self.duplicate_col as usize] < self.duplicate_col_lower {
            // Prospective value exceeds the lower bound, so trim it to the
            // lower bound and force recalculation of col_value[col]
            solution.col_value[self.duplicate_col as usize] = self.duplicate_col_lower;
            recompute_col = true;
            if basis.valid {
                basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kLower;
            }
        } else if self.duplicate_col_integral {
            debug_assert!(!basis.valid);
            // fractionality(x) = abs(x - round(x)), HighsUtils.h
            let v = solution.col_value[self.duplicate_col as usize];
            if (v - v.round()).abs() > options.mip_feasibility_tolerance {
                // Prospective value is not integer, so trim it to the integer
                // value and force recalculation of col_value[col]
                solution.col_value[self.duplicate_col as usize] = v.floor();
                recompute_col = true;
            }
        }

        if recompute_col {
            // Determine the value of col_value[col] corresponding to
            // col_value[duplicateCol]
            solution.col_value[self.col as usize] =
                merge_val - self.col_scale * solution.col_value[self.duplicate_col as usize];
            if !self.duplicate_col_integral && self.col_integral {
                // If column is integral and duplicateCol is not we need to make sure
                // we split the values into an integral one for col
                debug_assert!(!basis.valid);
                if h16_dominated_cols {
                    // HiGHS 1.16 #3285 (h16_dominated_cols): HighsPostsolveStack.cpp lines 715-720 (1.16.0)
                    solution.col_value[self.col as usize] = if duplicate_col_at_upper
                        == (self.col_scale > 0.0)
                    {
                        (solution.col_value[self.col as usize] - options.mip_feasibility_tolerance)
                            .ceil()
                    } else {
                        (solution.col_value[self.col as usize] + options.mip_feasibility_tolerance)
                            .floor()
                    };
                } else {
                    solution.col_value[self.col as usize] = (solution.col_value[self.col as usize]
                        - options.mip_feasibility_tolerance)
                        .ceil();
                }
                solution.col_value[self.duplicate_col as usize] = f64::from(
                    (HighsCDouble::from(merge_val) - solution.col_value[self.col as usize])
                        / self.col_scale,
                );
            }
        } else {
            // Have set col_value[col] to its lower bound (if finite),
            // otherwise min(0, upper), and if there is a basis, then
            // col_status[col] is kBasic.
            if basis.valid {
                // Now set the basis status of col and duplicateCol.
                //
                // Set col_status[col] to be nonbasic corresponding to its
                // value. This is fine if it's at its lower or upper bound, or
                // at zero and free, but what if it's at zero with a finite
                // upper bound?
                //
                // The basis is maintained by setting duplicateCol to be basic
                debug_assert!(basis.col_status[self.col as usize] == HighsBasisStatus::kBasic);
                basis.col_status[self.duplicate_col as usize] = basis.col_status[self.col as usize];
                if self.col_lower != -kHighsInf {
                    debug_assert!(solution.col_value[self.col as usize] == self.col_lower);
                    basis.col_status[self.col as usize] = HighsBasisStatus::kLower;
                } else if self.col_upper <= 0.0 {
                    debug_assert!(solution.col_value[self.col as usize] == self.col_upper);
                    basis.col_status[self.col as usize] = HighsBasisStatus::kUpper;
                } else {
                    debug_assert!(solution.col_value[self.col as usize] == 0.0);
                    basis.col_status[self.col as usize] = HighsBasisStatus::kZero;
                    if self.col_upper < kHighsInf {
                        // Nonbasic at zero with bounds (-inf, colUpper)
                        if h16_dupcol_warn {
                            // HiGHS 1.16 #3580 (h16_dupcol_warn): HighsPostsolveStack.cpp lines 749-756 (1.16.0)
                            crate::highs_log_user!(
                                &options.log_options,
                                HighsLogType::kWarning,
                                "HighsPostsolveStack::DuplicateColumn::undo Col is nonbasic at zero with upper bound of {}\n",
                                crate::util::c_printf::c_fmt_g(self.col_upper, 6)
                            );
                        } else {
                            println!(
                                "HighsPostsolveStack::DuplicateColumn::undo Col is nonbasic at zero with upper bound of {}",
                                self.col_upper
                            );
                            if cfg!(debug_assertions) {
                                panic!("assert(111 == 679)");
                            }
                        }
                    }
                }
                debug_assert!(
                    basis.col_status[self.duplicate_col as usize] == HighsBasisStatus::kBasic
                );
            }
        }
        // Check that any basis status for duplicateCol has been set
        if basis.valid {
            debug_assert!(
                basis.col_status[self.duplicate_col as usize] != HighsBasisStatus::kNonbasic
            );
        }

        let has_error = |sol: &HighsSolution| -> bool {
            let illegal_duplicate_col_lower = sol.col_value[self.duplicate_col as usize]
                < self.duplicate_col_lower - options.mip_feasibility_tolerance;
            let illegal_duplicate_col_upper = sol.col_value[self.duplicate_col as usize]
                > self.duplicate_col_upper + options.mip_feasibility_tolerance;
            let illegal_col_lower = sol.col_value[self.col as usize]
                < self.col_lower - options.mip_feasibility_tolerance;
            let illegal_col_upper = sol.col_value[self.col as usize]
                > self.col_upper + options.mip_feasibility_tolerance;
            let illegal_residual = !ok_residual(
                sol.col_value[self.col as usize],
                sol.col_value[self.duplicate_col as usize],
            );
            illegal_duplicate_col_lower
                || illegal_duplicate_col_upper
                || illegal_col_lower
                || illegal_col_upper
                || illegal_residual
        };

        if has_error(solution) {
            if debug_report {
                println!(
                    "DuplicateColumn::undo error: col = {}({}), duplicateCol = {}({})\n{}\n{}\n{} {} {}\n{} {} {}",
                    self.col,
                    solution.col_value[self.col as usize],
                    self.duplicate_col,
                    solution.col_value[self.duplicate_col as usize],
                    merge_val,
                    self.col_scale,
                    self.col_lower,
                    self.col_upper,
                    self.col_integral as i32,
                    self.duplicate_col_lower,
                    self.duplicate_col_upper,
                    self.duplicate_col_integral as i32
                );
            }
            // Fix error due to undo
            self.undo_fix(options, solution, merge_val);
        } else {
            return;
        }
        let allow_assert = false;
        if allow_assert {
            debug_assert!(!has_error(solution));
        }
        // Following undoFix, set any basis status, ideally keeping col basic
        if basis.valid {
            let mut duplicate_col_basic = false;
            if self.duplicate_col_lower <= -kHighsInf && self.duplicate_col_upper >= kHighsInf {
                // duplicateCol is free, so may be zero
                if solution.col_value[self.duplicate_col as usize] == 0.0 {
                    basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
                    basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kZero;
                } else {
                    duplicate_col_basic = true;
                }
            } else if is_at_bound(
                solution.col_value[self.duplicate_col as usize],
                self.duplicate_col_lower,
            ) {
                basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
                basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kLower;
            } else if is_at_bound(
                solution.col_value[self.duplicate_col as usize],
                self.duplicate_col_upper,
            ) {
                basis.col_status[self.col as usize] = HighsBasisStatus::kBasic;
                basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kUpper;
            } else {
                // duplicateCol is not free or at a bound, so must be basic
                duplicate_col_basic = true;
            }
            if duplicate_col_basic {
                // duplicateCol must be basic
                basis.col_status[self.duplicate_col as usize] = HighsBasisStatus::kBasic;
                // Hopefully col can be nonbasic
                if is_at_bound(solution.col_value[self.col as usize], self.col_lower) {
                    basis.col_status[self.col as usize] = HighsBasisStatus::kLower;
                } else if is_at_bound(solution.col_value[self.col as usize], self.col_upper) {
                    basis.col_status[self.col as usize] = HighsBasisStatus::kUpper;
                } else {
                    basis.col_status[self.col as usize] = HighsBasisStatus::kNonbasic;
                    if debug_report {
                        println!("When demerging, neither col nor duplicateCol can be nonbasic");
                    }
                    if kAllowDeveloperAssert && cfg!(debug_assertions) {
                        panic!("assert(666 == 999)");
                    }
                }
            }
        }
    }

    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::DuplicateColumn::okMerge (lines 869-1012, 144 lines)
    pub fn ok_merge(&self, tolerance: f64) -> bool {
        // When merging x and y to x+a.y, not all values of a are permitted,
        // since it must be possible to map back onto feasible values of x
        // and y.
        //
        // Assume WLOG that a > 0, x in [x_l, x_u], y in [y_l, y_u]
        //
        // Let z = x + a.y
        //
        // Range for z is [x_l+a.y_l, x_u+a.y_u]
        //
        // * If x and y are both integer:
        //
        // z will be integer and x+a.y must generate all integer values in
        // [x_l+a.y_l, x_u+a.y_u]. Hence a must be an integer. If a >=
        // (x_u-x_l)+2 then, since [a.y_l, a.y_u] contains integer multiples
        // of a, some of the intervening integers don't correspond to a
        // value of x. Hence a must be an integer and a <= (x_u-x_l)+1
        //
        // For example, if x and y are binary, then x+a.y is [0, 1, a,
        // 1+a]. For this to be a continuous sequence of integers, we must
        // have a <= 2.
        //
        // * If x is integer and y is continuous:
        //
        // z will be continuous and x+a.y must generate all values in
        // [x_l+a.y_l, x_u+a.y_u]. Since [x_l, x_u] are integers, [a.y_l,
        // a.y_u] = a[y_l, y_u] must be of length at least 1. Hence a must
        // be at least 1/(y_u-y_l) in magnitude.
        //
        // * If x is continuous and y is integer:
        //
        // z will be continuous and x+a.y must generate all values in
        // [x_l+a.y_l, x_u+a.y_u]. Since [a.y_l, a.y_u] contains integer
        // multiples of a, the gaps between them must not exceed the length
        // of [x_l, x_u]. Hence a must be at most x_u-x_l in
        // magnitude.
        //
        // Observe that this is equivalent to requiring 1/a to be at least
        // 1/(x_u-x_l) in magnitude, the symmetric result corresponding to
        // the merge (1/a)x+y.
        //
        // * If x and y are both continuous
        //
        // z will be continuous and x+a.y naturally generates all values in
        // [x_l+a.y_l, x_u+a.y_u].

        let scale = self.col_scale;
        let x_int = self.col_integral;
        let y_int = self.duplicate_col_integral;
        let x_lo = if x_int {
            (self.col_lower - tolerance).ceil()
        } else {
            self.col_lower
        };
        let x_up = if x_int {
            (self.col_upper + tolerance).floor()
        } else {
            self.col_upper
        };
        let y_lo = if y_int {
            (self.duplicate_col_lower - tolerance).ceil()
        } else {
            self.duplicate_col_lower
        };
        let y_up = if y_int {
            (self.duplicate_col_upper + tolerance).floor()
        } else {
            self.duplicate_col_upper
        };
        let x_len = x_up - x_lo;
        let y_len = y_up - y_lo;
        let mut newline: &str = "\n";
        let mut ok_merge = true;
        let debug_report = false;
        if scale == 0.0 {
            if debug_report {
                println!(
                    "{}DuplicateColumn::checkMerge: Scale cannot be zero",
                    newline
                );
            }
            newline = "";
            ok_merge = false;
        }
        let abs_scale = scale.abs();
        if x_int {
            if y_int {
                // Scale must be integer and not exceed (x_u-x_l)+1 in magnitude
                // fractionality(scale) = abs(scale - round(scale)), HighsUtils.h
                if (scale - scale.round()).abs() > tolerance {
                    if debug_report {
                        println!(
                            "{}DuplicateColumn::checkMerge: scale must be integer, but is {}",
                            newline, scale
                        );
                    }
                    newline = "";
                    ok_merge = false;
                }
                let scale_limit = x_len + 1.0 + tolerance;
                if abs_scale > scale_limit {
                    if debug_report {
                        println!(
                            "{}DuplicateColumn::checkMerge: scale = {}, but |scale| cannot exceed {} since x is [{}, {}]",
                            newline, scale, scale_limit, x_lo, x_up
                        );
                    }
                    newline = "";
                    ok_merge = false;
                }
            } else {
                // y is continuous
                if debug_report {
                    println!("DuplicateColumn::checkMerge: x-integer; y-continuous");
                }
                // Scale must be at least 1/(y_u-y_l) in magnitude
                if y_len == 0.0 {
                    if debug_report {
                        println!(
                            "{}DuplicateColumn::checkMerge: scale = {} is too small in magnitude, as y is [{}, {}]",
                            newline, scale, y_lo, y_up
                        );
                    }
                    newline = "";
                    ok_merge = false;
                } else {
                    let scale_limit = 1.0 / y_len;
                    if abs_scale < scale_limit {
                        if debug_report {
                            println!(
                                "{}DuplicateColumn::checkMerge: scale = {}, but |scale| must be at least {} since y is [{}, {}]",
                                newline, scale, scale_limit, y_lo, y_up
                            );
                        }
                        newline = "";
                        ok_merge = false;
                    }
                }
            }
        } else if y_int {
            if debug_report {
                println!("DuplicateColumn::checkMerge: x-continuous; y-integer");
            }
            // Scale must be at most (x_u-x_l) in magnitude
            let scale_limit = x_len;
            if abs_scale > scale_limit {
                if debug_report {
                    println!(
                        "{}DuplicateColumn::checkMerge: scale = {}, but |scale| must be at most {} since x is [{}, {}]",
                        newline, scale, scale_limit, x_lo, x_up
                    );
                }
                newline = "";
                ok_merge = false;
            }
        } else {
            // x and y are continuous: no restriction on the scale
        }
        // `newline` is only read by the debug_report prints (const false);
        // this read keeps those C++ assignments live for the compiler.
        let _ = newline;
        ok_merge
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// HiGHS 1.16 #3580 test model (3580.diff 5-67, `presolve-free-integer-column`); optimum -468.
    fn h16_dupcol_undo_mip_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 6,
            num_row_: 2,
            col_cost_: vec![-50.0, 0.0, 0.0, -3.0, 5.0, 0.0],
            col_lower_: vec![-kHighsInf, -kHighsInf, -kHighsInf, 0.0, 0.0, -3.0],
            col_upper_: vec![9.0, kHighsInf, kHighsInf, 6.0, 1.0, 6.0],
            row_lower_: vec![-kHighsInf, -2.0],
            row_upper_: vec![-6.0, -2.0],
            integrality_: vec![
                HighsVarType::kContinuous,
                HighsVarType::kContinuous,
                HighsVarType::kInteger,
                HighsVarType::kInteger,
                HighsVarType::kInteger,
                HighsVarType::kInteger,
            ],
            ..Default::default()
        };
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 6;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 1, 3, 4, 5, 6, 7];
        lp.a_matrix_.index_ = vec![1, 0, 1, 1, 1, 1, 0];
        lp.a_matrix_.value_ = vec![-1.0, -1.0, -1.0, -3.0, 1.0, 4.0, -1.0];
        lp
    }

    /// Runs the #3580 MIP with both HiGHS 1.16 undo flags (#3285 and #3580) set to `on`
    /// and returns the objective. Both flags are needed to reach the 1.16 optimum.
    fn h16_dupcol_undo_mip_objective(on: bool) -> f64 {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("h16_dupcol_warn", on),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_dominated_cols", on),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(h16_dupcol_undo_mip_lp()), HighsStatus::kOk);
        highs.run();
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        highs.get_info().objective_function_value
    }

    /// HiGHS 1.16 #3285/#3580 on path: the optimum -468 (checked by hand: x = (9, 2, -1, 6, 0, 4)).
    #[test]
    fn h16_dupcol_undo_mip_optimum_on() {
        let objective = h16_dupcol_undo_mip_objective(true);
        assert!((objective + 468.0).abs() < 1e-8, "objective = {objective}");
    }

    /// Off path: the 1.15.1 result on this model, pinned. It is suboptimal (-418
    /// against the optimum -468 above); #3285 fixes it only behind the flags.
    #[test]
    fn h16_dupcol_undo_mip_optimum_off() {
        let objective = h16_dupcol_undo_mip_objective(false);
        assert!((objective + 418.0).abs() < 1e-8, "objective = {objective}");
    }

    /// Basic col at zero with a finite upper bound: the branch HiGHS 1.16 #3580 changes.
    /// Only the log call differs, so the postsolved values are the same in both flags.
    fn h16_dupcol_undo_zero_basis_case(dupcol_warn: bool) -> (HighsSolution, HighsBasis) {
        let mut options = HighsOptions::default();
        options.h16_dupcol_warn = dupcol_warn;
        let dupcol = DuplicateColumn {
            col_scale: 1.0,
            col_lower: -kHighsInf,
            col_upper: 5.0,
            duplicate_col_lower: 0.0,
            duplicate_col_upper: 10.0,
            col: 0,
            duplicate_col: 1,
            col_integral: false,
            duplicate_col_integral: false,
        };
        let mut solution = HighsSolution {
            col_value: vec![3.0, 0.0],
            ..Default::default()
        };
        let mut basis = HighsBasis {
            valid: true,
            col_status: vec![HighsBasisStatus::kBasic, HighsBasisStatus::kBasic],
            ..Default::default()
        };
        dupcol.undo(&options, &mut solution, &mut basis);
        (solution, basis)
    }

    #[test]
    fn h16_dupcol_undo_zero_basis_warn_on() {
        let (solution, basis) = h16_dupcol_undo_zero_basis_case(true);
        assert_eq!(solution.col_value, vec![0.0, 3.0]);
        assert_eq!(
            basis.col_status,
            vec![HighsBasisStatus::kZero, HighsBasisStatus::kBasic]
        );
    }

    #[test]
    #[cfg_attr(debug_assertions, should_panic(expected = "assert(111 == 679)"))]
    fn h16_dupcol_undo_zero_basis_warn_off() {
        let (solution, basis) = h16_dupcol_undo_zero_basis_case(false);
        assert_eq!(solution.col_value, vec![0.0, 3.0]);
        assert_eq!(
            basis.col_status,
            vec![HighsBasisStatus::kZero, HighsBasisStatus::kBasic]
        );
    }

    /// Duplicate column merged at its upper bound (3285.diff 326-356). Integer
    /// col, continuous duplicate col, col_scale = `scale`, merged value `z`.
    fn h16_dupcol_undo_at_upper_case(
        dominated_cols: bool,
        col_scale: f64,
        merge_val: f64,
    ) -> HighsSolution {
        let mut options = HighsOptions::default();
        options.h16_dominated_cols = dominated_cols;
        let dupcol = DuplicateColumn {
            col_scale,
            col_lower: -kHighsInf,
            col_upper: 10.0,
            duplicate_col_lower: 0.0,
            duplicate_col_upper: 1.0,
            col: 0,
            duplicate_col: 1,
            col_integral: true,
            duplicate_col_integral: false,
        };
        let mut solution = HighsSolution {
            col_value: vec![merge_val, 0.0],
            ..Default::default()
        };
        let mut basis = HighsBasis::default();
        dupcol.undo(&options, &mut solution, &mut basis);
        solution
    }

    /// With a positive scale the flag does not change the rounding direction:
    /// both flags give ceil(1.5) = 2 and duplicate col 0.75.
    #[test]
    fn h16_dupcol_undo_at_upper_positive_scale_same() {
        for dominated_cols in [false, true] {
            let solution = h16_dupcol_undo_at_upper_case(dominated_cols, 2.0, 3.5);
            assert_eq!(solution.col_value, vec![2.0, 0.75], "{dominated_cols}");
        }
    }

    /// With a negative scale and the duplicate col at its upper bound, the on path
    /// rounds col up/down with floor(col + tol) = floor(-1.5) = -2 and gives
    /// duplicate col (z - col) / scale = 0.75, which is feasible.
    #[test]
    fn h16_dupcol_undo_at_upper_negative_scale_on() {
        let solution = h16_dupcol_undo_at_upper_case(true, -2.0, -3.5);
        assert_eq!(solution.col_value, vec![-2.0, 0.75]);
    }

    /// Off path: 1.15.1 rounds with ceil(-1.5 - tol) = -1, which leaves the duplicate
    /// col above its upper bound, so undo_fix repairs it. The result must be feasible.
    #[test]
    fn h16_dupcol_undo_at_upper_negative_scale_off() {
        let solution = h16_dupcol_undo_at_upper_case(false, -2.0, -3.5);
        let (col, dup) = (solution.col_value[0], solution.col_value[1]);
        assert!(dup >= -1e-6 && dup <= 1.0 + 1e-6, "dup = {dup}");
        assert!(col <= 10.0 + 1e-6, "col = {col}");
        assert!((col - 2.0 * dup + 3.5).abs() < 1e-6, "residual");
    }
}
