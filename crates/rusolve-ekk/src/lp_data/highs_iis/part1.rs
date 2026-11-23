//! HiGHS source: lp_data/HighsIis.cpp (lines 14-366)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;
use crate::prelude::*;
use crate::util::c_printf::{c_fmt_f, c_fmt_g, c_fmt_g_width, c_pad_left, c_pad_right};

impl HighsIis {
    /// HiGHS: HighsIis.cpp::HighsIis::clearData (lines 14-25)
    pub fn clear_data(&mut self) {
        self.valid_ = false;
        self.status_ = kIisModelStatusUnknown;
        self.strategy_ = kIisStrategyMin;
        self.col_index_.clear();
        self.row_index_.clear();
        self.col_bound_.clear();
        self.row_bound_.clear();
        self.col_status_.clear();
        self.row_status_.clear();
        self.model_.clear();
    }

    /// HiGHS: HighsIis.cpp::HighsIis::clear (lines 27-30)
    pub fn clear(&mut self) {
        self.clear_data();
        self.info_.clear();
    }

    /// HiGHS: HighsIis.cpp::HighsIis::clearLogInfo (lines 32-35)
    pub fn clear_log_info(&mut self) {
        self.info_.iis_last_disptime = -kHighsInf;
        self.info_.iis_num_disp_lines = 0;
    }

    /// HiGHS: HighsIis.cpp::HighsIis::invalid (lines 37-41)
    pub fn invalid(&mut self, lp: &HighsLp) {
        self.clear();
        self.col_status_ = vec![kIisStatusMaybeInConflict; lp.num_col_ as usize];
        self.row_status_ = vec![kIisStatusMaybeInConflict; lp.num_row_ as usize];
    }

    /// HiGHS: HighsIis.cpp::HighsIis::iisBoundStatusToString (lines 43-51)
    pub fn iis_bound_status_to_string(&self, bound_status: HighsInt) -> String {
        if bound_status == kIisBoundStatusDropped {
            return "Dropped".to_string();
        }
        if bound_status == kIisBoundStatusNull {
            return "   Null".to_string();
        }
        if bound_status == kIisBoundStatusFree {
            return "   Free".to_string();
        }
        if bound_status == kIisBoundStatusLower {
            return "  Lower".to_string();
        }
        if bound_status == kIisBoundStatusUpper {
            return "  Upper".to_string();
        }
        if bound_status == kIisBoundStatusBoxed {
            return "  Boxed".to_string();
        }
        "*****".to_string()
    }

    /// HiGHS: HighsIis.cpp::HighsIis::iisModelStatusToString (lines 53-60)
    pub fn iis_model_status_to_string(&self, model_status: HighsInt) -> String {
        if model_status == kIisModelStatusFeasible {
            return "Feasible".to_string();
        }
        if model_status == kIisModelStatusUnknown {
            return "Unknown".to_string();
        }
        if model_status == kIisModelStatusTimeLimit {
            return "Time limit reached".to_string();
        }
        if model_status == kIisModelStatusReducible {
            return "Reducible".to_string();
        }
        if model_status == kIisModelStatusIrreducible {
            return "Irreducible".to_string();
        }
        "*****".to_string()
    }

    /// HiGHS: HighsIis.cpp::HighsIis::report (lines 62-85)
    pub fn report(&self, message: &str, lp: &HighsLp) {
        let num_iis_col = self.col_index_.len() as HighsInt;
        let num_iis_row = self.row_index_.len() as HighsInt;
        if num_iis_col > 10 || num_iis_row > 10 {
            return;
        }
        print!("\nIIS {}\n===\n", message);
        print!("Column: ");
        for i_col in 0..num_iis_col {
            print!("{:>9} ", i_col);
        }
        print!("\nStatus: ");
        for i_col in 0..num_iis_col {
            print!(
                "{:>9} ",
                self.iis_bound_status_to_string(self.col_bound_[i_col as usize])
            );
        }
        print!("\nLower:  ");
        for i_col in 0..num_iis_col {
            print!("{} ", c_fmt_g_width(lp.col_lower_[i_col as usize], 2, 9));
        }
        print!("\nUpper:  ");
        for i_col in 0..num_iis_col {
            print!("{} ", c_fmt_g_width(lp.col_upper_[i_col as usize], 2, 9));
        }
        println!();
        println!("Row:    Status     Lower     Upper");
        for i_row in 0..num_iis_row {
            println!(
                "{:>2}   {:>9} {} {}",
                i_row,
                self.iis_bound_status_to_string(self.row_bound_[i_row as usize]),
                c_fmt_g_width(lp.row_lower_[i_row as usize], 2, 9),
                c_fmt_g_width(lp.row_upper_[i_row as usize], 2, 9)
            );
        }
        println!();
    }

    /// HiGHS: HighsIis.cpp::HighsIis::reportIteration (lines 87-114)
    pub fn report_iteration(
        &mut self,
        options: &HighsOptions,
        iter: HighsInt,
        num_rows_remaining: HighsInt,
        force: bool,
    ) {
        let output_flag = options.log_options.output_flag;
        if !output_flag {
            return;
        }
        let min_interval = 5.0;
        let runtime = self.info_.sum_simplex_times;
        if !force
            && self.info_.iis_last_disptime > -0.5 * kHighsInf
            && runtime - self.info_.iis_last_disptime < min_interval
        {
            return;
        }
        // Update last time only when we actually print a line
        self.info_.iis_last_disptime = runtime;
        let gap = 17;
        let w_iter = 9 + 2; // strlen("Iteration") + 2
        let w_rows = 4 + gap; // strlen("Rows") + gap
        let w_time = 7 + gap; // strlen("Runtime") + gap

        // Print header every 20 lines (and on first line)
        if self.info_.iis_num_disp_lines % 20 == 0 {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "{}{}{}\n",
                c_pad_right("Iteration", w_iter as usize),
                c_pad_right("Rows", w_rows as usize),
                c_pad_right("Runtime", w_time as usize)
            );
        }
        self.info_.iis_num_disp_lines += 1;
        let time_string = format!("{}s", c_fmt_f(runtime, 2));
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "{}{}{}\n",
            c_pad_right(&iter.to_string(), w_iter as usize),
            c_pad_right(&num_rows_remaining.to_string(), w_rows as usize),
            c_pad_right(&time_string, w_time as usize)
        );
    }

    /// HiGHS: HighsIis.cpp::HighsIis::reportFinal (lines 116-131)
    pub fn report_final(&self, options: &HighsOptions) {
        let output_flag = options.log_options.output_flag;
        if !output_flag {
            return;
        }
        let log_options = &options.log_options;
        // Align colons by padding labels to the same width
        let k_label_width = 19;
        crate::highs_log_user!(log_options, HighsLogType::kInfo, "\n");
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "{} : {}\n",
            c_pad_left("IIS status", k_label_width as usize),
            self.iis_model_status_to_string(self.status_)
        );
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "{} : {}\n",
            c_pad_left("Rows", k_label_width as usize),
            self.row_index_.len()
        );
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "{} : {}\n",
            c_pad_left("Columns", k_label_width as usize),
            self.col_index_.len()
        );
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "{} : {}s\n",
            c_pad_left("HiGHS run time", k_label_width as usize),
            c_fmt_f(self.info_.sum_simplex_times, 2)
        );
    }

    /// HiGHS: HighsIis.cpp::HighsIis::addCol (lines 133-136)
    /// C++ default arguments: status = kIisBoundStatusNull
    pub fn add_col(&mut self, col: HighsInt, status: HighsInt) {
        self.col_index_.push(col);
        self.col_bound_.push(status);
    }

    /// HiGHS: HighsIis.cpp::HighsIis::addRow (lines 138-141)
    /// C++ default arguments: status = kIisBoundStatusNull
    pub fn add_row(&mut self, row: HighsInt, status: HighsInt) {
        self.row_index_.push(row);
        self.row_bound_.push(status);
    }

    /// HiGHS: HighsIis.cpp::HighsIis::removeCol (lines 143-148)
    pub fn remove_col(&mut self, col: HighsInt) {
        let num_col = self.col_index_.len() as HighsInt;
        debug_assert!(col < num_col);
        self.col_index_[col as usize] = self.col_index_[(num_col - 1) as usize];
        self.col_index_.resize((num_col - 1) as usize, 0);
    }

    /// HiGHS: HighsIis.cpp::HighsIis::removeRow (lines 150-155)
    pub fn remove_row(&mut self, row: HighsInt) {
        let num_row = self.row_index_.len() as HighsInt;
        debug_assert!(row < num_row);
        self.row_index_[row as usize] = self.row_index_[(num_row - 1) as usize];
        self.row_index_.resize((num_row - 1) as usize, 0);
    }

    /// HiGHS: HighsIis.cpp::HighsIis::trivial (lines 157-226)
    pub fn trivial(&mut self, lp: &HighsLp, options: &HighsOptions) -> bool {
        self.clear();
        let col_priority = (kIisStrategyColPriority & options.iis_strategy) != 0;
        for k in 0..2 {
            if (col_priority && k == 0) || (!col_priority && k == 1) {
                // Loop over columns first
                for i_col in 0..lp.num_col_ {
                    if lp.col_lower_[i_col as usize] - lp.col_upper_[i_col as usize]
                        > 2.0 * options.primal_feasibility_tolerance
                    {
                        self.add_col(i_col, kIisBoundStatusBoxed);
                        break;
                    }
                }
                if self.col_index_.len() > 0 {
                    break;
                }
            } else {
                // Loop over rows first
                for i_row in 0..lp.num_row_ {
                    if lp.row_lower_[i_row as usize] - lp.row_upper_[i_row as usize]
                        > 2.0 * options.primal_feasibility_tolerance
                    {
                        self.add_row(i_row, kIisBoundStatusBoxed);
                        break;
                    }
                }
                if self.row_index_.len() > 0 {
                    break;
                }
            }
        }
        let num_iis_col = self.col_index_.len() as HighsInt;
        let num_iis_row = self.row_index_.len() as HighsInt;
        // If one is found then we've found an IIS
        if num_iis_col + num_iis_row > 0 {
            // Should have found exactly 1
            debug_assert!(num_iis_col + num_iis_row == 1);
            self.valid_ = true;
            self.status_ = kIisModelStatusIrreducible;
            self.strategy_ = options.iis_strategy;
            return true;
        }
        // Now look for empty rows that cannot have zero activity
        let mut count: Vec<HighsInt> = Vec::new();
        // Get the row counts
        if lp.a_matrix_.is_colwise() {
            count = vec![0; lp.num_row_ as usize];
            for i_col in 0..lp.num_col_ {
                for i_el in
                    lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
                {
                    count[lp.a_matrix_.index_[i_el as usize] as usize] += 1;
                }
            }
        } else {
            for i_row in 0..lp.num_row_ {
                count.push(
                    lp.a_matrix_.start_[(i_row + 1) as usize] - lp.a_matrix_.start_[i_row as usize],
                );
            }
        }
        debug_assert!(self.row_index_.len() == 0);
        for i_row in 0..lp.num_row_ {
            if count[i_row as usize] > 0 {
                continue;
            }
            if lp.row_lower_[i_row as usize] > options.primal_feasibility_tolerance {
                self.add_row(i_row, kIisBoundStatusLower);
            } else if lp.row_upper_[i_row as usize] < -options.primal_feasibility_tolerance {
                self.add_row(i_row, kIisBoundStatusUpper);
            }
            if self.row_index_.len() > 0 {
                // If one is found then we've found an IIS
                self.valid_ = true;
                self.status_ = kIisModelStatusIrreducible;
                self.strategy_ = options.iis_strategy;
                return true;
            }
        }
        false
    }

    /// HiGHS: HighsIis.cpp::HighsIis::rowValueBounds (lines 228-366)
    pub fn row_value_bounds(&mut self, lp: &HighsLp, options: &HighsOptions) -> bool {
        // Look for infeasible rows based on row value bounds
        self.clear();
        let mut lower_value: Vec<f64> = Vec::new();
        let mut upper_value: Vec<f64> = Vec::new();
        if lp.a_matrix_.is_colwise() {
            lower_value = vec![0.0; lp.num_row_ as usize];
            upper_value = vec![0.0; lp.num_row_ as usize];
            for i_col in 0..lp.num_col_ {
                let lower = lp.col_lower_[i_col as usize];
                let upper = lp.col_upper_[i_col as usize];
                for i_el in
                    lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
                {
                    let i_row = lp.a_matrix_.index_[i_el as usize];
                    let value = lp.a_matrix_.value_[i_el as usize];
                    if value > 0.0 {
                        lower_value[i_row as usize] += value * lower;
                        upper_value[i_row as usize] += value * upper;
                    } else {
                        lower_value[i_row as usize] += value * upper;
                        upper_value[i_row as usize] += value * lower;
                    }
                }
            }
        } else {
            for i_row in 0..lp.num_row_ {
                let mut lower_row_value = 0.0;
                let mut upper_row_value = 0.0;
                for i_el in
                    lp.a_matrix_.start_[i_row as usize]..lp.a_matrix_.start_[(i_row + 1) as usize]
                {
                    let i_col = lp.a_matrix_.index_[i_el as usize];
                    let lower = lp.col_lower_[i_col as usize];
                    let upper = lp.col_upper_[i_col as usize];
                    let value = lp.a_matrix_.value_[i_el as usize];
                    if value > 0.0 {
                        lower_row_value += value * lower;
                        upper_row_value += value * upper;
                    } else {
                        lower_row_value += value * upper;
                        upper_row_value += value * lower;
                    }
                }
                lower_value.push(lower_row_value);
                upper_value.push(upper_row_value);
            }
        }
        let mut below_lower = false;
        let mut above_upper = false; // C++ leaves it uninitialised; read only after a row set it
        for i_row in 0..lp.num_row_ {
            below_lower = upper_value[i_row as usize]
                < lp.row_lower_[i_row as usize] - options.primal_feasibility_tolerance;
            above_upper = lower_value[i_row as usize]
                > lp.row_upper_[i_row as usize] + options.primal_feasibility_tolerance;
            if below_lower || above_upper {
                self.row_index_.push(i_row);
                if below_lower {
                    self.row_bound_.push(kIisBoundStatusLower);
                } else {
                    self.row_bound_.push(kIisBoundStatusUpper);
                }
                break;
            }
        }
        if self.row_index_.len() == 0 {
            // Nothing found, but IIS data still valid
            self.clear();
            self.valid_ = true;
            self.status_ = kIisModelStatusUnknown;
            self.strategy_ = options.iis_strategy;
            return false;
        }
        debug_assert!(below_lower || above_upper);
        debug_assert!(!(below_lower && above_upper));
        // Found an infeasible row
        let i_row = self.row_index_[0];
        let row_name_string = if lp.row_names_.len() > 0 {
            format!("({})", lp.row_names_[i_row as usize])
        } else {
            String::new()
        };
        if below_lower {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "LP row {} {}has maximum row value of {}, below lower bound of {}\n",
                i_row,
                row_name_string,
                c_fmt_g(upper_value[i_row as usize], 6),
                c_fmt_g(lp.row_lower_[i_row as usize], 6)
            );
        } else {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "LP row {} {}has minimum row value of {}, above upper bound of {}\n",
                i_row,
                row_name_string,
                c_fmt_g(lower_value[i_row as usize], 6),
                c_fmt_g(lp.row_upper_[i_row as usize], 6)
            );
        }
        let mut value;
        // C++ lambda setColBound (lines 319-333)
        let set_col_bound = |below_lower: bool, value: f64, col_bound: &mut Vec<HighsInt>| {
            if below_lower {
                if value > 0.0 {
                    col_bound.push(kIisBoundStatusUpper);
                } else {
                    col_bound.push(kIisBoundStatusLower);
                }
            } else {
                if value > 0.0 {
                    col_bound.push(kIisBoundStatusLower);
                } else {
                    col_bound.push(kIisBoundStatusUpper);
                }
            }
        };

        if lp.a_matrix_.is_colwise() {
            for i_col in 0..lp.num_col_ {
                for i_el in
                    lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
                {
                    value = lp.a_matrix_.value_[i_el as usize];
                    if lp.a_matrix_.index_[i_el as usize] == i_row && value != 0.0 {
                        self.col_index_.push(i_col);
                        set_col_bound(below_lower, value, &mut self.col_bound_);
                    }
                }
            }
        } else {
            for i_el in
                lp.a_matrix_.start_[i_row as usize]..lp.a_matrix_.start_[(i_row + 1) as usize]
            {
                let i_col = lp.a_matrix_.index_[i_el as usize];
                value = lp.a_matrix_.value_[i_el as usize];
                if value != 0.0 {
                    self.col_index_.push(i_col);
                    set_col_bound(below_lower, value, &mut self.col_bound_);
                }
            }
        }

        // There must be at least one column in the IIS
        debug_assert!(self.col_index_.len() > 0);
        debug_assert!(self.col_index_.len() == self.col_bound_.len());
        debug_assert!(self.row_index_.len() == self.row_bound_.len());
        self.valid_ = true;
        self.status_ = kIisModelStatusIrreducible;
        self.strategy_ = options.iis_strategy;
        true
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

    #[test]
    fn strings_and_index_lists() {
        let mut iis = HighsIis::default();
        assert_eq!(
            iis.iis_bound_status_to_string(kIisBoundStatusLower),
            "  Lower"
        );
        assert_eq!(iis.iis_bound_status_to_string(99), "*****");
        assert_eq!(
            iis.iis_model_status_to_string(kIisModelStatusIrreducible),
            "Irreducible"
        );
        iis.add_col(3, kIisBoundStatusUpper);
        iis.add_row(1, kIisBoundStatusBoxed);
        iis.add_row(4, kIisBoundStatusLower);
        assert_eq!(
            (iis.col_index_.clone(), iis.col_bound_.clone()),
            (vec![3], vec![kIisBoundStatusUpper])
        );
        iis.remove_row(0);
        assert_eq!(iis.row_index_, vec![4]);
        iis.clear();
        assert!(iis.col_index_.is_empty() && iis.row_index_.is_empty() && !iis.valid_);
        assert_eq!(iis.status_, kIisModelStatusUnknown);
    }

    #[test]
    fn trivial_detects_inconsistent_bounds() {
        let options = HighsOptions::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![2.0];
        lp.col_upper_ = vec![1.0];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.start_ = vec![0, 0];
        let mut iis = HighsIis::default();
        assert!(iis.trivial(&lp, &options));
        assert_eq!(iis.col_index_, vec![0]);
        assert_eq!(iis.col_bound_, vec![kIisBoundStatusBoxed]);
    }

    /// HighsIis.cpp lines 194-224: an empty row whose lower bound is positive.
    #[test]
    fn trivial_detects_empty_row() {
        let options = HighsOptions::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![1.0];
        lp.row_lower_ = vec![0.0, 1.0];
        lp.row_upper_ = vec![1.0, 2.0];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let mut iis = HighsIis::default();
        assert!(iis.trivial(&lp, &options));
        assert!(iis.col_index_.is_empty());
        assert_eq!(iis.row_index_, vec![1]);
        assert_eq!(iis.row_bound_, vec![kIisBoundStatusLower]);
        assert!(iis.valid_);
        assert_eq!(iis.status_, kIisModelStatusIrreducible);
    }

    /// HighsIis.cpp lines 233-251, 305-310, 335-345: column-wise matrix,
    /// row 1 (x0 - x1 >= 2, x in [0,1]) has maximum activity 1.
    #[test]
    fn row_value_bounds_colwise_below_lower() {
        let options = HighsOptions::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0];
        lp.row_lower_ = vec![-kHighsInf, 2.0];
        lp.row_upper_ = vec![kHighsInf, kHighsInf];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 3];
        lp.a_matrix_.index_ = vec![0, 1, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, -1.0];
        let mut iis = HighsIis::default();
        assert!(iis.row_value_bounds(&lp, &options));
        assert_eq!(iis.row_index_, vec![1]);
        assert_eq!(iis.row_bound_, vec![kIisBoundStatusLower]);
        assert_eq!(iis.col_index_, vec![0, 1]);
        assert_eq!(
            iis.col_bound_,
            vec![kIisBoundStatusUpper, kIisBoundStatusLower]
        );
        assert!(iis.valid_);
        assert_eq!(iis.status_, kIisModelStatusIrreducible);
    }

    /// HighsIis.cpp lines 252-273, 311-317, 346-356: row-wise matrix,
    /// x0 - 2 x1 + 0 x2 <= -7 with x0 in [1,2], x1 in [0,3] has minimum
    /// activity -5; the zero entry is not added to the IIS.
    #[test]
    fn row_value_bounds_rowwise_above_upper() {
        let options = HighsOptions::default();
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0, 0.0, 0.0];
        lp.col_lower_ = vec![1.0, 0.0, 0.0];
        lp.col_upper_ = vec![2.0, 3.0, 1.0];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![-7.0];
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 3];
        lp.a_matrix_.index_ = vec![0, 1, 2];
        lp.a_matrix_.value_ = vec![1.0, -2.0, 0.0];
        let mut iis = HighsIis::default();
        assert!(iis.row_value_bounds(&lp, &options));
        assert_eq!(iis.row_index_, vec![0]);
        assert_eq!(iis.row_bound_, vec![kIisBoundStatusUpper]);
        assert_eq!(iis.col_index_, vec![0, 1]);
        assert_eq!(
            iis.col_bound_,
            vec![kIisBoundStatusLower, kIisBoundStatusUpper]
        );
    }

    /// HighsIis.cpp lines 291-298: nothing found, IIS data still valid.
    #[test]
    fn row_value_bounds_nothing_found() {
        let mut options = HighsOptions::default();
        options.iis_strategy = kIisStrategyColPriority;
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![1.0];
        lp.row_lower_ = vec![0.5];
        lp.row_upper_ = vec![0.5];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        let mut iis = HighsIis::default();
        assert!(!iis.row_value_bounds(&lp, &options));
        assert!(iis.row_index_.is_empty() && iis.col_index_.is_empty());
        assert!(iis.valid_);
        assert_eq!(iis.status_, kIisModelStatusUnknown);
        assert_eq!(iis.strategy_, kIisStrategyColPriority);
    }

    /// HighsIis.h lines 60-74.
    #[test]
    fn info_update_tracks_extremes() {
        let mut info = HighsIisInfo::default();
        info.update(2.0, 5);
        info.update(1.0, 7);
        assert_eq!(info.num_lp_solved, 2);
        assert_eq!(info.sum_simplex_times, 3.0);
        assert_eq!(info.min_simplex_time, 1.0);
        assert_eq!(info.max_simplex_time, 2.0);
        assert_eq!(info.sum_simplex_iteration_counts, 12);
        assert_eq!(info.min_simplex_iteration_count, 5);
        assert_eq!(info.max_simplex_iteration_count, 7);
    }
}
