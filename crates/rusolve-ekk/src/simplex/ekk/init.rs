//! HiGHS source: simplex/HEkk.cpp (lines 2327-2772)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Nonbasic moves, work arrays, cost/bound initialisation including the
//! numerically-sensitive cost and bound perturbation.

use crate::highs_log_dev;
use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::setNonbasicMove (lines 2327-2396, 70 lines)
    pub fn set_nonbasic_move(&mut self) {
        let have_solution = false;
        // Don't have a simplex basis since nonbasicMove is not set up.

        // Assign nonbasicMove using as much information as is available
        let mut lower: f64;
        let mut upper: f64;
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        self.basis_.nonbasic_move_.resize(num_tot as usize, 0);

        for i_var in 0..num_tot {
            if self.basis_.nonbasic_flag_[i_var as usize] == 0 {
                // Basic variable
                self.basis_.nonbasic_move_[i_var as usize] = kNonbasicMoveZe;
                continue;
            }
            // Nonbasic variable
            if i_var < self.lp_.num_col_ {
                lower = self.lp_.col_lower_[i_var as usize];
                upper = self.lp_.col_upper_[i_var as usize];
            } else {
                let i_row = i_var - self.lp_.num_col_;
                lower = -self.lp_.row_upper_[i_row as usize];
                upper = -self.lp_.row_lower_[i_row as usize];
            }
            let mut mv: i8 = kIllegalMoveValue;
            if lower == upper {
                // Fixed
                mv = kNonbasicMoveZe;
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    //
                    // Determine the bound to set the value to according to, in order of
                    // priority
                    //
                    // 1. Any solution value
                    if have_solution {
                        let midpoint = 0.5 * (lower + upper);
                        let value = self.info_.work_value_[i_var as usize];
                        if value < midpoint {
                            mv = kNonbasicMoveUp;
                        } else {
                            mv = kNonbasicMoveDn;
                        }
                    }
                    // 2. Bound of original LP that is closer to zero
                    if mv == kIllegalMoveValue {
                        if lower.abs() < upper.abs() {
                            mv = kNonbasicMoveUp;
                        } else {
                            mv = kNonbasicMoveDn;
                        }
                    }
                } else {
                    // Lower (since upper bound is infinite)
                    mv = kNonbasicMoveUp;
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                mv = kNonbasicMoveDn;
            } else {
                // FREE
                mv = kNonbasicMoveZe;
            }
            debug_assert!(mv != kIllegalMoveValue);
            self.basis_.nonbasic_move_[i_var as usize] = mv;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::allocateWorkAndBaseArrays (lines 2398-2419, 22 lines)
    pub fn allocate_work_and_base_arrays(&mut self) {
        let num_tot = (self.lp_.num_col_ + self.lp_.num_row_) as usize;
        self.info_.work_cost_.resize(num_tot, 0.0);
        self.info_.work_dual_.resize(num_tot, 0.0);
        self.info_.work_shift_.resize(num_tot, 0.0);

        self.info_.work_lower_.resize(num_tot, 0.0);
        self.info_.work_upper_.resize(num_tot, 0.0);
        self.info_.work_range_.resize(num_tot, 0.0);
        self.info_.work_value_.resize(num_tot, 0.0);
        self.info_.work_lower_shift_.resize(num_tot, 0.0);
        self.info_.work_upper_shift_.resize(num_tot, 0.0);

        // Feel that it should be possible to resize this with in dual
        // solver, and only if Devex is being used, but a pointer to it
        // needs to be set up when constructing HDual
        self.info_.devex_index_.resize(num_tot, 0);

        let num_row = self.lp_.num_row_ as usize;
        self.info_.base_lower_.resize(num_row, 0.0);
        self.info_.base_upper_.resize(num_row, 0.0);
        self.info_.base_value_.resize(num_row, 0.0);
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseLpColBound (lines 2421-2429, 9 lines)
    pub fn initialise_lp_col_bound(&mut self) {
        for i_col in 0..self.lp_.num_col_ as usize {
            self.info_.work_lower_[i_col] = self.lp_.col_lower_[i_col];
            self.info_.work_upper_[i_col] = self.lp_.col_upper_[i_col];
            self.info_.work_range_[i_col] =
                self.info_.work_upper_[i_col] - self.info_.work_lower_[i_col];
            self.info_.work_lower_shift_[i_col] = 0.0;
            self.info_.work_upper_shift_[i_col] = 0.0;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseLpRowBound (lines 2431-2440, 10 lines)
    pub fn initialise_lp_row_bound(&mut self) {
        for i_row in 0..self.lp_.num_row_ as usize {
            let i_col = self.lp_.num_col_ as usize + i_row;
            self.info_.work_lower_[i_col] = -self.lp_.row_upper_[i_row];
            self.info_.work_upper_[i_col] = -self.lp_.row_lower_[i_row];
            self.info_.work_range_[i_col] =
                self.info_.work_upper_[i_col] - self.info_.work_lower_[i_col];
            self.info_.work_lower_shift_[i_col] = 0.0;
            self.info_.work_upper_shift_[i_col] = 0.0;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseCost (lines 2442-2569, 128 lines)
    /// HiGHS 1.16 (#3226): HEkk.cpp::HEkk::initialiseCost (lines 2472-2577)
    /// MUST-PORT (determinism, random stream): cost perturbation reads `num_tot_random_value_`; keep the expression order exactly.
    /// C++ default arguments: perturb = false
    pub fn initialise_cost(
        &mut self,
        algorithm: SimplexAlgorithm,
        solve_phase: HighsInt,
        perturb: bool,
    ) {
        let _ = solve_phase;
        // Copy the cost
        self.initialise_lp_col_cost();
        self.initialise_lp_row_cost();
        self.info_.costs_shifted = false;
        self.info_.costs_perturbed = false;
        self.analysis_.net_num_single_cost_shift = 0;
        // Primal simplex costs are either from the LP or set specially in phase 1
        if algorithm == SimplexAlgorithm::kPrimal {
            return;
        }
        // Dual simplex costs are either from the LP or perturbed
        if !perturb || self.info_.dual_simplex_cost_perturbation_multiplier == 0.0 {
            return;
        }
        // Perturb the original costs, scale down if is too big
        // HiGHS 1.16 #3226 (h16_cost_perturb): HEkk.cpp::initialiseCost lines 2483-2516 (1.16.0)
        let h16_cost_perturb: bool = self.options_.h16_cost_perturb;
        let report_cost_perturbation = self.options_.output_flag; // && analysis_.analyse_simplex_runtime_data;
        let mut num_original_nonzero_cost: HighsInt = 0;
        let mut min_abs_cost = kHighsInf;
        let mut max_abs_cost: f64 = 0.0;
        let mut sum_abs_cost: f64 = 0.0;
        if h16_cost_perturb {
            highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Cost perturbation for {}\n",
                self.lp_.model_name_
            );
            for i in 0..self.lp_.num_col_ as usize {
                let abs_cost = self.info_.work_cost_[i].abs();
                if abs_cost != 0.0 {
                    num_original_nonzero_cost += 1;
                    min_abs_cost = cpp_min(min_abs_cost, abs_cost);
                }
                sum_abs_cost += abs_cost;
                max_abs_cost = cpp_max(max_abs_cost, abs_cost);
            }
            let pct0 = (100 * num_original_nonzero_cost) / self.lp_.num_col_;
            highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "   Initially have {} nonzero costs ({:3}%)",
                num_original_nonzero_cost,
                pct0
            );
            if num_original_nonzero_cost != 0 {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    " with min / average / max = {} / {} / {}\n",
                    min_abs_cost,
                    sum_abs_cost / num_original_nonzero_cost as f64,
                    max_abs_cost
                );
            } else {
                max_abs_cost = 1.0;
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    " but perturb as if max cost was 1\n"
                );
            }
        } else {
            if report_cost_perturbation {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Cost perturbation for {}\n",
                    self.lp_.model_name_
                );
            }
            for i in 0..self.lp_.num_col_ as usize {
                let abs_cost = self.info_.work_cost_[i].abs();
                if report_cost_perturbation {
                    if abs_cost != 0.0 {
                        num_original_nonzero_cost += 1;
                        min_abs_cost = cpp_min(min_abs_cost, abs_cost);
                    }
                    sum_abs_cost += abs_cost;
                }
                max_abs_cost = cpp_max(max_abs_cost, abs_cost);
            }
            let pct0 = (100 * num_original_nonzero_cost) / self.lp_.num_col_;
            let mut average_abs_cost: f64 = 0.0;
            if report_cost_perturbation {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "   Initially have {} nonzero costs ({:3}%)",
                    num_original_nonzero_cost,
                    pct0
                );
                if num_original_nonzero_cost != 0 {
                    average_abs_cost = sum_abs_cost / num_original_nonzero_cost as f64;
                    highs_log_dev!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        " with min / average / max = {} / {} / {}\n",
                        min_abs_cost,
                        average_abs_cost,
                        max_abs_cost
                    );
                } else {
                    min_abs_cost = 1.0;
                    max_abs_cost = 1.0;
                    average_abs_cost = 1.0;
                    highs_log_dev!(
                        &self.options_.log_options,
                        HighsLogType::kInfo,
                        " but perturb as if max cost was 1\n"
                    );
                }
            }
            let _ = (min_abs_cost, average_abs_cost);
        }
        if max_abs_cost > 100.0 {
            max_abs_cost = max_abs_cost.sqrt().sqrt();
            if h16_cost_perturb {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "   Large so set max_abs_cost = sqrt(sqrt(max_abs_cost)) = {}\n",
                    max_abs_cost
                );
            } else if report_cost_perturbation {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "   Large so set max_abs_cost = sqrt(sqrt(max_abs_cost)) = {}\n",
                    max_abs_cost
                );
            }
        }

        // If there are few boxed variables, we will just use simple perturbation
        let mut boxed_rate: f64 = 0.0;
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        for i in 0..num_tot as usize {
            boxed_rate += (self.info_.work_range_[i] < 1e30) as i32 as f64;
        }
        boxed_rate /= num_tot as f64;
        if boxed_rate < 0.01 {
            max_abs_cost = cpp_min(max_abs_cost, 1.0);
            if h16_cost_perturb {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "   Small boxedRate ({}) so set max_abs_cost = min(max_abs_cost, 1.0) = {}\n",
                    boxed_rate,
                    max_abs_cost
                );
            } else if report_cost_perturbation {
                highs_log_dev!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "   Small boxedRate ({}) so set max_abs_cost = min(max_abs_cost, 1.0) = {}\n",
                    boxed_rate,
                    max_abs_cost
                );
            }
        }
        // Determine the perturbation base
        self.cost_perturbation_max_abs_cost_ = max_abs_cost;
        self.cost_perturbation_base_ =
            self.info_.dual_simplex_cost_perturbation_multiplier * 5e-7 * max_abs_cost;
        if h16_cost_perturb {
            highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "   Perturbation column base = {}\n",
                self.cost_perturbation_base_
            );
        } else if report_cost_perturbation {
            highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "   Perturbation column base = {}\n",
                self.cost_perturbation_base_
            );
        }

        // Now do the perturbation
        for i in 0..self.lp_.num_col_ as usize {
            let lower = self.lp_.col_lower_[i];
            let upper = self.lp_.col_upper_[i];
            let xpert = (1.0 + self.info_.num_tot_random_value_[i])
                * (self.info_.work_cost_[i].abs() + 1.0)
                * self.cost_perturbation_base_;
            if lower <= -kHighsInf && upper >= kHighsInf {
                // Free - no perturb
            } else if upper >= kHighsInf {
                // Lower
                self.info_.work_cost_[i] += xpert;
            } else if lower <= -kHighsInf {
                // Upper
                self.info_.work_cost_[i] += -xpert;
            } else if lower != upper {
                // Boxed
                self.info_.work_cost_[i] += if self.info_.work_cost_[i] >= 0.0 {
                    xpert
                } else {
                    -xpert
                };
            } else {
                // Fixed - no perturb
            }
        }
        let row_cost_perturbation_base_ =
            self.info_.dual_simplex_cost_perturbation_multiplier * 1e-12;
        if h16_cost_perturb {
            highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "   Perturbation row    base = {}\n",
                row_cost_perturbation_base_
            );
        } else if report_cost_perturbation {
            highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "   Perturbation row    base = {}\n",
                row_cost_perturbation_base_
            );
        }
        for i in self.lp_.num_col_ as usize..num_tot as usize {
            let perturbation2 =
                (0.5 - self.info_.num_tot_random_value_[i]) * row_cost_perturbation_base_;
            self.info_.work_cost_[i] += perturbation2;
        }
        self.info_.costs_perturbed = true;
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseBound (lines 2571-2697, 127 lines)
    /// MUST-PORT (determinism, random stream): bound perturbation reads `num_tot_random_value_`; keep the expression order exactly.
    /// C++ default arguments: perturb = false
    pub fn initialise_bound(
        &mut self,
        algorithm: SimplexAlgorithm,
        solve_phase: HighsInt,
        perturb: bool,
    ) {
        self.initialise_lp_col_bound();
        self.initialise_lp_row_bound();
        self.info_.bounds_shifted = false;
        self.info_.bounds_perturbed = false;
        // Primal simplex bounds are either from the LP or perturbed
        if algorithm == SimplexAlgorithm::kPrimal {
            if !perturb || self.info_.primal_simplex_bound_perturbation_multiplier == 0.0 {
                return;
            }
            // Perturb the bounds
            // Determine the smallest and largest finite lower/upper bounds
            let num_col = self.lp_.num_col_;
            let num_row = self.lp_.num_row_;
            let num_tot = num_col + num_row;
            let mut min_abs_lower = kHighsInf;
            let mut max_abs_lower: f64 = -1.0;
            let mut min_abs_upper = kHighsInf;
            let mut max_abs_upper: f64 = -1.0;
            for i_var in 0..num_tot as usize {
                let abs_lower = self.info_.work_lower_[i_var].abs();
                let abs_upper = self.info_.work_upper_[i_var].abs();
                if abs_lower != 0.0 && abs_lower < kHighsInf {
                    min_abs_lower = cpp_min(abs_lower, min_abs_lower);
                    max_abs_lower = cpp_max(abs_lower, max_abs_lower);
                }
                if abs_upper != 0.0 && abs_upper < kHighsInf {
                    min_abs_upper = cpp_min(abs_upper, min_abs_upper);
                    max_abs_upper = cpp_max(abs_upper, max_abs_upper);
                }
            }
            let _ = (min_abs_lower, max_abs_lower, min_abs_upper, max_abs_upper);

            let base = self.info_.primal_simplex_bound_perturbation_multiplier * 5e-7;
            for i_var in 0..num_tot as usize {
                let mut lower = self.info_.work_lower_[i_var];
                let mut upper = self.info_.work_upper_[i_var];
                let fixed = lower == upper;
                // Don't perturb bounds of nonbasic fixed variables as they stay nonbasic
                if self.basis_.nonbasic_flag_[i_var] == kNonbasicFlagTrue && fixed {
                    continue;
                }
                let random_value = self.info_.num_tot_random_value_[i_var];
                if lower > -kHighsInf {
                    if lower < -1.0 {
                        lower -= random_value * base * (-lower);
                    } else if lower < 1.0 {
                        lower -= random_value * base;
                    } else {
                        lower -= random_value * base * lower;
                    }
                    self.info_.work_lower_[i_var] = lower;
                }
                if upper < kHighsInf {
                    if upper < -1.0 {
                        upper += random_value * base * (-upper);
                    } else if upper < 1.0 {
                        upper += random_value * base;
                    } else {
                        upper += random_value * base * upper;
                    }
                    self.info_.work_upper_[i_var] = upper;
                }
                self.info_.work_range_[i_var] =
                    self.info_.work_upper_[i_var] - self.info_.work_lower_[i_var];
                if self.basis_.nonbasic_flag_[i_var] == kNonbasicFlagFalse {
                    continue;
                }
                // Set values of nonbasic variables
                if self.basis_.nonbasic_move_[i_var] > 0 {
                    self.info_.work_value_[i_var] = lower;
                } else if self.basis_.nonbasic_move_[i_var] < 0 {
                    self.info_.work_value_[i_var] = upper;
                }
            }
            for i_row in 0..num_row as usize {
                let i_var = self.basis_.basic_index_[i_row] as usize;
                self.info_.base_lower_[i_row] = self.info_.work_lower_[i_var];
                self.info_.base_upper_[i_row] = self.info_.work_upper_[i_var];
            }
            self.info_.bounds_perturbed = true;
            return;
        }
        // Dual simplex bounds are either from the LP or set to special values in
        // phase 1
        debug_assert!(algorithm == SimplexAlgorithm::kDual);
        if solve_phase == kSolvePhase2 {
            return;
        }

        // The dual objective is the sum of products of primal and dual
        // values for nonbasic variables. For dual simplex phase 1, the
        // primal bounds are set so that when the dual value is feasible, the
        // primal value is set to zero. Otherwise the value is +1/-1
        // according to the required sign of the dual, except for free
        // variables, where the bounds are [-1000, 1000]. Hence the dual
        // objective is the negation of the sum of infeasibilities, unless there are
        // free In Phase 1: change to dual phase 1 bound.
        let inf = kHighsInf;
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        for i_col in 0..num_tot as usize {
            if self.info_.work_lower_[i_col] == -inf && self.info_.work_upper_[i_col] == inf {
                // Phase 1 bounds were not modified for free rows in hsol.
                self.info_.work_lower_[i_col] = -1000.0;
                self.info_.work_upper_[i_col] = 1000.0; // FREE
            } else if self.info_.work_lower_[i_col] == -inf {
                self.info_.work_lower_[i_col] = -1.0;
                self.info_.work_upper_[i_col] = 0.0; // UPPER
            } else if self.info_.work_upper_[i_col] == inf {
                self.info_.work_lower_[i_col] = 0.0;
                self.info_.work_upper_[i_col] = 1.0; // LOWER
            } else {
                self.info_.work_lower_[i_col] = 0.0;
                self.info_.work_upper_[i_col] = 0.0; // BOXED or FIXED
            }
            self.info_.work_range_[i_col] =
                self.info_.work_upper_[i_col] - self.info_.work_lower_[i_col];
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseLpColCost (lines 2699-2706, 8 lines)
    pub fn initialise_lp_col_cost(&mut self) {
        let cost_scale_factor = 2.0f64.powf(self.options_.cost_scale_factor as f64);
        for i_col in 0..self.lp_.num_col_ as usize {
            self.info_.work_cost_[i_col] = (self.lp_.sense_ as HighsInt) as f64
                * cost_scale_factor
                * self.lp_.col_cost_[i_col];
            self.info_.work_shift_[i_col] = 0.0;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseLpRowCost (lines 2708-2714, 7 lines)
    pub fn initialise_lp_row_cost(&mut self) {
        for i_col in self.lp_.num_col_..self.lp_.num_col_ + self.lp_.num_row_ {
            self.info_.work_cost_[i_col as usize] = 0.0;
            self.info_.work_shift_[i_col as usize] = 0.0;
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::initialiseNonbasicValueAndMove (lines 2716-2772, 57 lines)
    pub fn initialise_nonbasic_value_and_move(&mut self) {
        // Initialise workValue and nonbasicMove from nonbasicFlag and
        // bounds, except for boxed variables when nonbasicMove is used to
        // set workValue=workLower/workUpper
        let num_tot = self.lp_.num_col_ + self.lp_.num_row_;
        for i_var in 0..num_tot as usize {
            if self.basis_.nonbasic_flag_[i_var] == 0 {
                // Basic variable
                self.basis_.nonbasic_move_[i_var] = kNonbasicMoveZe;
                continue;
            }
            // Nonbasic variable
            let lower = self.info_.work_lower_[i_var];
            let upper = self.info_.work_upper_[i_var];
            let original_move = self.basis_.nonbasic_move_[i_var];
            let value: f64;
            let mv: i8;
            if lower == upper {
                // Fixed
                value = lower;
                mv = kNonbasicMoveZe;
            } else if !highs_is_infinity(-lower) {
                // Finite lower bound so boxed or lower
                if !highs_is_infinity(upper) {
                    // Finite upper bound so boxed
                    if original_move == kNonbasicMoveUp {
                        // Set at lower
                        value = lower;
                        mv = kNonbasicMoveUp;
                    } else if original_move == kNonbasicMoveDn {
                        // Set at upper
                        value = upper;
                        mv = kNonbasicMoveDn;
                    } else {
                        // Invalid nonbasicMove: correct and set value at lower
                        value = lower;
                        mv = kNonbasicMoveUp;
                    }
                } else {
                    // Lower
                    value = lower;
                    mv = kNonbasicMoveUp;
                }
            } else if !highs_is_infinity(upper) {
                // Upper
                value = upper;
                mv = kNonbasicMoveDn;
            } else {
                // FREE
                value = 0.0;
                mv = kNonbasicMoveZe;
            }
            debug_assert!(mv != kIllegalMoveValue);
            self.basis_.nonbasic_move_[i_var] = mv;
            self.info_.work_value_[i_var] = value;
        }
    }
}

#[cfg(test)]
mod h16_cost_perturb_tests {
    use crate::prelude::*;

    /// Three columns (boxed, lower-only, upper-only) and one row with a finite upper bound.
    /// `col_cost` is the original objective; the random stream is fixed at 0.25.
    /// HiGHS 1.16 (#3226): HEkk.cpp::HEkk::initialiseCost (lines 2472-2577)
    fn h16_cost_perturb_ekk(col_cost: Vec<f64>, output_flag: bool, flag: bool) -> HEkk {
        let mut ekk = HEkk::default();
        ekk.options_.output_flag = output_flag;
        ekk.options_.h16_cost_perturb = flag;
        ekk.lp_.num_col_ = 3;
        ekk.lp_.num_row_ = 1;
        ekk.lp_.col_cost_ = col_cost;
        ekk.lp_.col_lower_ = vec![0.0, 0.0, -kHighsInf];
        ekk.lp_.col_upper_ = vec![10.0, kHighsInf, 5.0];
        ekk.lp_.row_lower_ = vec![-kHighsInf];
        ekk.lp_.row_upper_ = vec![4.0];
        ekk.allocate_work_and_base_arrays();
        ekk.initialise_lp_col_bound();
        ekk.initialise_lp_row_bound();
        ekk.info_.dual_simplex_cost_perturbation_multiplier = 1.0;
        ekk.info_.num_tot_random_value_ = vec![0.25; 4];
        ekk
    }

    /// HiGHS 1.16 (#3226): HEkk.cpp::HEkk::initialiseCost (lines 2472-2577)
    #[test]
    fn h16_cost_perturb_zero_cost_on() {
        let mut ekk = h16_cost_perturb_ekk(vec![0.0, 0.0, 0.0], false, true);
        ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, true);
        // No nonzero cost: the 1.16 path perturbs as if max cost was 1.
        assert_eq!(ekk.cost_perturbation_max_abs_cost_, 1.0);
        let expected = (1.0 + 0.25) * (0.0 + 1.0) * 5e-7;
        assert!((ekk.info_.work_cost_[0] - expected).abs() < 1e-18);
        assert!((ekk.info_.work_cost_[1] - expected).abs() < 1e-18);
        assert!((ekk.info_.work_cost_[2] + expected).abs() < 1e-18);
        assert!(ekk.info_.costs_perturbed);
    }

    /// HiGHS 1.16 (#3226): HEkk.cpp::HEkk::initialiseCost (lines 2472-2577)
    #[test]
    fn h16_cost_perturb_zero_cost_off() {
        let mut ekk = h16_cost_perturb_ekk(vec![0.0, 0.0, 0.0], false, false);
        ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, true);
        // 1.15.1 with output_flag off: max_abs_cost stays 0, so no column perturbation.
        assert_eq!(ekk.cost_perturbation_max_abs_cost_, 0.0);
        assert_eq!(ekk.cost_perturbation_base_, 0.0);
        assert_eq!(ekk.info_.work_cost_[0], 0.0);
        assert_eq!(ekk.info_.work_cost_[1], 0.0);
        assert_eq!(ekk.info_.work_cost_[2], 0.0);
        assert!(ekk.info_.costs_perturbed);
    }

    /// HiGHS 1.16 (#3226): HEkk.cpp::HEkk::initialiseCost (lines 2472-2577)
    /// The 1.16 path and the 1.15.1 path give bit-identical costs whenever the
    /// output flag is on, and whenever the LP has a nonzero cost. (Zero costs with the
    /// output flag off are the case the change alters, see the tests above.)
    #[test]
    fn h16_cost_perturb_output_agrees() {
        for (col_cost, output_flag) in [
            (vec![0.0, 0.0, 0.0], true),
            (vec![3.0, -4.0, 0.0], false),
            (vec![3.0, -4.0, 0.0], true),
        ] {
            let mut on = h16_cost_perturb_ekk(col_cost.clone(), output_flag, true);
            let mut off = h16_cost_perturb_ekk(col_cost.clone(), output_flag, false);
            on.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, true);
            off.initialise_cost(SimplexAlgorithm::kDual, kSolvePhase2, true);
            assert_eq!(
                on.cost_perturbation_max_abs_cost_, off.cost_perturbation_max_abs_cost_,
                "cost {col_cost:?} output_flag {output_flag}"
            );
            assert_eq!(
                on.info_.work_cost_, off.info_.work_cost_,
                "cost {col_cost:?} output_flag {output_flag}"
            );
        }
    }
}
