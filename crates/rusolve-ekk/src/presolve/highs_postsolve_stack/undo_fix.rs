//! HiGHS source: presolve/HighsPostsolveStack.cpp (lines 1014-1393)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `DuplicateColumn::undoFix`, `DuplicateColumn::transformToPresolvedSpace`
//! and `SlackColSubstitution::undo`. The file-static helpers are in
//! `super` (see `undo_a.rs`).

use super::*;
use crate::prelude::*;

impl DuplicateColumn {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::DuplicateColumn::undoFix (lines 1014-1320, 307 lines)
    pub fn undo_fix(&self, options: &HighsOptions, solution: &mut HighsSolution, merge_value: f64) {
        let mip_feasibility_tolerance = options.mip_feasibility_tolerance;
        let primal_feasibility_tolerance = options.primal_feasibility_tolerance;
        let col_value = &mut solution.col_value;
        let allow_assert = false;
        let debug_report = false;
        //=============================================================================================

        let is_integer = |v: f64| -> bool { (v - v.round()).abs() <= mip_feasibility_tolerance };

        let is_feasible = |l: f64, v: f64, u: f64| -> bool {
            v >= l - primal_feasibility_tolerance && v <= u + primal_feasibility_tolerance
        };
        let scale = self.col_scale;
        let x_int = self.col_integral;
        let y_int = self.duplicate_col_integral;
        let x_ix = self.col;
        let y_ix = self.duplicate_col;
        let x_lo = if x_int {
            (self.col_lower - mip_feasibility_tolerance).ceil()
        } else {
            self.col_lower
        };
        let x_up = if x_int {
            (self.col_upper + mip_feasibility_tolerance).floor()
        } else {
            self.col_upper
        };
        let y_lo = if y_int {
            (self.duplicate_col_lower - mip_feasibility_tolerance).ceil()
        } else {
            self.duplicate_col_lower
        };
        let y_up = if y_int {
            (self.duplicate_col_upper + mip_feasibility_tolerance).floor()
        } else {
            self.duplicate_col_upper
        };
        if kAllowDeveloperAssert {
            debug_assert!(scale != 0.0);
        }
        let mut x_v = merge_value;
        let mut y_v = -kHighsInf;

        let check_int_var = |z_lo: f64,
                             z_up: f64,
                             z_0: &mut f64,
                             z_d: &mut f64,
                             z_1: &mut f64,
                             z_free: &mut bool| {
            let value_max: f64 = 1000.0;
            *z_0 = 0.0;
            *z_d = 0.0;
            *z_1 = 0.0;
            *z_free = false;
            if z_lo <= -kHighsInf {
                if z_up >= kHighsInf {
                    // z is free
                    *z_free = true;
                    *z_0 = 0.0;
                    *z_d = 1.0;
                    *z_1 = value_max;
                } else {
                    // z is (-int, u]
                    *z_0 = z_up;
                    *z_d = -1.0;
                    *z_1 = -value_max;
                }
            } else {
                // z is [l, u] or [l, inf) respectively
                *z_0 = z_lo;
                *z_d = 1.0;
                *z_1 = if z_up >= kHighsInf { value_max } else { z_up };
            }
        };

        let compute_value =
            |value: f64| -> f64 { f64::from((HighsCDouble::from(merge_value) - value) / scale) };

        let compute_inv_value = |value: f64| -> f64 {
            f64::from(HighsCDouble::from(merge_value) - HighsCDouble::from(value) * scale)
        };

        let find_value = |z_value: &mut f64,
                          z_0: f64,
                          z_1: f64,
                          z_delta: f64,
                          other_value: &mut f64,
                          other_lower: f64,
                          other_upper: f64,
                          other_int: bool|
         -> bool {
            let eps = 1e-8;
            // C++ `for (z_value = z_0;; z_value += z_delta)`: the loop never
            // exits through its condition, only through the returns below.
            *z_value = z_0;
            loop {
                *other_value = compute_value(*z_value);
                if is_feasible(other_lower, *other_value, other_upper)
                    && (!other_int || is_integer(*other_value))
                {
                    return true;
                }
                if z_delta > 0.0 && *z_value + z_delta >= z_1 + eps {
                    return false;
                }
                if z_delta < 0.0 && *z_value + z_delta <= z_1 - eps {
                    return false;
                }
                *z_value += z_delta;
            }
        };

        // C++ lambdas setXValue / setYValue capture x_v and y_v by reference;
        // here they take them as parameters (see crate docs, rule (f)).
        let set_x_value = |x_v: &mut f64, y_v: &mut f64, value: f64| {
            *x_v = value;
            *y_v = compute_value(value);
        };

        let set_y_value = |x_v: &mut f64, y_v: &mut f64, value: f64| {
            *y_v = value;
            *x_v = compute_inv_value(value);
        };

        if x_int {
            let mut x_0 = 0.0;
            let mut x_d = 0.0;
            let mut x_1 = 0.0;
            let mut x_free = false;
            check_int_var(x_lo, x_up, &mut x_0, &mut x_d, &mut x_1, &mut x_free);
            // x is integer, so look through its possible values to find a
            // suitable y
            if x_free && debug_report {
                println!("DuplicateColumn::undo x is free");
            }
            if debug_report {
                println!("DuplicateColumn::undo Using x ({}; {}; {})", x_0, x_d, x_1);
            }
            let found_y = find_value(&mut x_v, x_0, x_1, x_d, &mut y_v, y_lo, y_up, y_int);
            if allow_assert {
                debug_assert!(found_y);
            }
        } else if y_int {
            let mut y_0 = 0.0;
            let mut y_d = 0.0;
            let mut y_1 = 0.0;
            let mut y_free = false;
            check_int_var(y_lo, y_up, &mut y_0, &mut y_d, &mut y_1, &mut y_free);
            // y is integer, so look through its possible values to find a
            // suitable x
            if y_free && debug_report {
                println!("DuplicateColumn::undo y is free");
            }
            if debug_report {
                println!("DuplicateColumn::undo Using y ({}; {}; {})", y_0, y_d, y_1);
            }
            let found_x = find_value(&mut y_v, y_0, y_1, y_d, &mut x_v, x_lo, x_up, x_int);
            if allow_assert {
                debug_assert!(found_x);
            }
        } else {
            // x and y are both continuous
            let v_m_a_ylo = if y_lo <= -kHighsInf {
                if scale > 0.0 { kHighsInf } else { -kHighsInf }
            } else {
                compute_inv_value(y_lo)
            };
            let v_m_a_yup = if y_up >= kHighsInf {
                if scale > 0.0 { -kHighsInf } else { kHighsInf }
            } else {
                compute_inv_value(y_up)
            };
            // Need to ensure that y puts x in [x_l, x_u]
            if scale > 0.0 {
                if debug_report {
                    println!(
                        "DuplicateColumn::undo [V-a(y_u), V-a(y_l)] == [{}, {}]",
                        v_m_a_yup, v_m_a_ylo
                    );
                }
                // V-ay is in [V-a(y_u), V-a(y_l)] == [v_m_a_yup, v_m_a_ylo]
                if y_up < kHighsInf {
                    // If v_m_a_yup is right of x_up+eps then [v_m_a_yup, v_m_a_ylo] is
                    // right of [x_lo-eps, x_up+eps] so there's no solution. [Could
                    // try v_m_a_ylo computed from y_lo-eps.]
                    if kAllowDeveloperAssert {
                        debug_assert!(x_up + primal_feasibility_tolerance >= v_m_a_yup);
                    }
                    // This assignment is OK unless x_v < x_lo-eps
                    set_y_value(&mut x_v, &mut y_v, y_up);
                    if x_v < x_lo - primal_feasibility_tolerance {
                        // Try y_v corresponding to x_lo
                        set_x_value(&mut x_v, &mut y_v, x_lo);
                        if y_v < y_lo - primal_feasibility_tolerance {
                            // Very tight: use x_v on its margin and hope!
                            set_x_value(&mut x_v, &mut y_v, x_lo - primal_feasibility_tolerance);
                        }
                    }
                } else if y_lo > -kHighsInf {
                    // If v_m_a_ylo is left of x_lo-eps then [v_m_a_yup, v_m_a_ylo] is
                    // left of [x_lo-eps, x_up+eps] so there's no solution. [Could
                    // try v_m_a_yup computed from y_up+eps.]
                    if kAllowDeveloperAssert {
                        debug_assert!(x_lo - primal_feasibility_tolerance <= v_m_a_ylo);
                    }
                    // This assignment is OK unless x_v > x_up-eps
                    set_y_value(&mut x_v, &mut y_v, y_lo);
                    if x_v > x_up + primal_feasibility_tolerance {
                        // Try y_v corresponding to x_up
                        set_x_value(&mut x_v, &mut y_v, x_up);
                        if y_v > y_up + primal_feasibility_tolerance {
                            // Very tight: use x_v on its margin and hope!
                            set_x_value(&mut x_v, &mut y_v, x_up + primal_feasibility_tolerance);
                        }
                    }
                } else {
                    // y is free, so use x_v = max(0, x_lo)
                    set_x_value(&mut x_v, &mut y_v, cpp_max(0.0, x_lo));
                }
            } else {
                // scale < 0
                if debug_report {
                    println!(
                        "DuplicateColumn::undo [V-a(y_l), V-a(y_u)] == [{}, {}]",
                        v_m_a_ylo, v_m_a_yup
                    );
                }
                // V-ay is in [V-a(y_l), V-a(y_u)] == [v_m_a_ylo, v_m_a_yup]
                if y_lo > -kHighsInf {
                    // If v_m_a_ylo is right of x_up+eps then [v_m_a_ylo, v_m_a_yup] is
                    // right of [x_lo-eps, x_up+eps] so there's no solution. [Could
                    // try v_m_a_ylo computed from y_up+eps.]
                    if kAllowDeveloperAssert {
                        debug_assert!(x_up + primal_feasibility_tolerance >= v_m_a_ylo);
                    }
                    // This assignment is OK unless x_v < x_lo-eps
                    set_y_value(&mut x_v, &mut y_v, y_lo);
                    if x_v < x_lo - primal_feasibility_tolerance {
                        // Try y_v corresponding to x_lo
                        set_x_value(&mut x_v, &mut y_v, x_lo);
                        if y_v > y_up + primal_feasibility_tolerance {
                            // Very tight: use x_v on its margin and hope!
                            set_x_value(&mut x_v, &mut y_v, x_lo - primal_feasibility_tolerance);
                        }
                    }
                } else if y_up < kHighsInf {
                    // If v_m_a_yup is left of x_lo-eps then [v_m_a_ylo, v_m_a_yup] is
                    // left of [x_lo-eps, x_up+eps] so there's no solution. [Could
                    // try v_m_a_yup computed from y_lo-eps.]
                    if kAllowDeveloperAssert {
                        debug_assert!(x_lo - primal_feasibility_tolerance <= v_m_a_yup);
                    }
                    // This assignment is OK unless x_v < x_lo-eps
                    set_y_value(&mut x_v, &mut y_v, y_up);
                    if x_v > x_up + primal_feasibility_tolerance {
                        // Try y_v corresponding to x_up
                        set_x_value(&mut x_v, &mut y_v, x_up);
                        if y_v < y_lo - primal_feasibility_tolerance {
                            // Very tight: use x_v on its margin and hope!
                            set_x_value(&mut x_v, &mut y_v, x_up + primal_feasibility_tolerance);
                        }
                    }
                } else {
                    // y is free, so use x_v = max(0, x_lo)
                    set_x_value(&mut x_v, &mut y_v, cpp_max(0.0, x_lo));
                }
            }
        }
        let residual_tolerance = 1e-12;
        let residual = f64::from(HighsCDouble::from(x_v) - compute_inv_value(y_v)).abs();
        let x_y_ok = is_feasible(x_lo, x_v, x_up)
            && is_feasible(y_lo, y_v, y_up)
            && (!x_int || is_integer(x_v))
            && (!y_int || is_integer(y_v))
            && (x_v.abs() < kHighsInf)
            && (y_v.abs() < kHighsInf)
            && (residual <= residual_tolerance);

        let mut check = is_feasible(x_lo, x_v, x_up);
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: isFeasible(x_lo, x_v, x_up) is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = is_feasible(y_lo, y_v, y_up);
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: isFeasible(y_lo, y_v, y_up) is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = !x_int || is_integer(x_v);
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: !x_int || isInteger(x_v) is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = !y_int || is_integer(y_v);
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: !y_int || isInteger(y_v) is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = x_v.abs() < kHighsInf;
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: std::fabs(x_v) < kHighsInf is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = y_v.abs() < kHighsInf;
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: std::fabs(y_v) < kHighsInf is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = residual <= residual_tolerance;
        if !check {
            if debug_report {
                println!("DuplicateColumn::undo error: residual <= residual_tolerance is false");
            }
            if allow_assert {
                debug_assert!(check);
            }
        }
        check = residual <= residual_tolerance;
        if debug_report {
            print!(
                "DuplicateColumn::undo{} x = {}; y = {} to give x + ({})y = {}",
                if x_y_ok { "" } else { " ERROR" },
                x_v,
                y_v,
                scale,
                merge_value
            );
        }
        if x_y_ok {
            if debug_report {
                println!(": FIXED");
            }
        } else if check {
            if debug_report {
                println!();
            }
        } else if debug_report {
            println!(": residual = {}", residual);
        }
        //=============================================================================================
        if x_y_ok {
            col_value[x_ix as usize] = x_v;
            col_value[y_ix as usize] = y_v;
        }
    }

    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::DuplicateColumn::transformToPresolvedSpace (lines 1322-1325, 4 lines)
    pub fn transform_to_presolved_space(&self, primal_sol: &mut [f64]) {
        primal_sol[self.col as usize] = primal_sol[self.col as usize]
            + self.col_scale * primal_sol[self.duplicate_col as usize];
    }
}

impl SlackColSubstitution {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::SlackColSubstitution::undo (lines 1327-1393, 67 lines)
    pub fn undo(
        &mut self,
        options: &HighsOptions,
        row_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        let debug_print = false;
        let row = self.row;
        let col = self.col;
        // May have to determine row dual and basis status unless doing
        // primal-only transformation in MIP solver, in which case row may
        // no longer exist if it corresponds to a removed cut, so have to
        // avoid exceeding array bounds of solution.row_value

        // compute primal values
        let mut col_coef = 0.0;
        let mut row_value = HighsCDouble::from(0.0);
        for row_val in row_values {
            if row_val.index == col {
                col_coef = row_val.value;
            } else {
                row_value += row_val.value * solution.col_value[row_val.index as usize];
            }
        }

        debug_assert!(col_coef != 0.0);
        // Row values aren't fully postsolved, so why do this?
        if (row as usize) < solution.row_value.len() {
            solution.row_value[row as usize] =
                f64::from(row_value + col_coef * solution.col_value[col as usize]);
        }

        solution.col_value[col as usize] = f64::from((self.rhs - row_value) / col_coef);

        // If no dual values requested, return here
        if !solution.dual_valid {
            return;
        }

        // Row retains its dual value, and column has this dual value scaled by coeff
        if (row as usize) < solution.row_value.len() {
            solution.col_dual[col as usize] = -solution.row_dual[row as usize] / col_coef;
        }

        // Set basis status if necessary
        if !basis.valid {
            return;
        }

        // If row is basic, then slack is basic, otherwise row retains its status
        if (row as usize) < solution.row_value.len() {
            let save_row_basis_status = basis.row_status[row as usize];
            if basis.row_status[row as usize] == HighsBasisStatus::kBasic {
                basis.col_status[col as usize] = HighsBasisStatus::kBasic;
                basis.row_status[row as usize] =
                    compute_row_status(solution.row_dual[row as usize], RowType::kEq);
            } else if basis.row_status[row as usize] == HighsBasisStatus::kLower {
                basis.col_status[col as usize] = if col_coef > 0.0 {
                    HighsBasisStatus::kUpper
                } else {
                    HighsBasisStatus::kLower
                };
            } else {
                basis.col_status[col as usize] = if col_coef > 0.0 {
                    HighsBasisStatus::kLower
                } else {
                    HighsBasisStatus::kUpper
                };
            }
            if debug_print {
                println!(
                    "HighsPostsolveStack::SlackColSubstitution::undo OgRowStatus = {}; RowStatus = {}; ColStatus = {}",
                    util_basis_status_to_string(save_row_basis_status),
                    util_basis_status_to_string(basis.row_status[row as usize]),
                    util_basis_status_to_string(basis.col_status[col as usize])
                );
            }
            if basis.col_status[col as usize] == HighsBasisStatus::kLower {
                debug_assert!(
                    solution.col_dual[col as usize] > -options.dual_feasibility_tolerance
                );
            } else if basis.col_status[col as usize] == HighsBasisStatus::kUpper {
                debug_assert!(solution.col_dual[col as usize] < options.dual_feasibility_tolerance);
            }
        } else {
            basis.col_status[col as usize] = HighsBasisStatus::kNonbasic;
        }
    }
}
