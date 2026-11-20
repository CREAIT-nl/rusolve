//! HiGHS source: lp_data/HighsIis.cpp (lines 981-1290)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g;

impl HighsIis {
    /// HiGHS: HighsIis.cpp::HighsIis::indexStatusOk (lines 981-1040)
    pub fn index_status_ok(&self, lp: &HighsLp) -> bool {
        let num_col = lp.num_col_;
        let num_row = lp.num_row_;
        let col_status_size_ok = self.col_status_.len() == num_col as usize;
        let row_status_size_ok = self.row_status_.len() == num_row as usize;
        debug_assert!(col_status_size_ok);
        debug_assert!(row_status_size_ok);
        if !col_status_size_ok {
            return self.index_status_ok_return(false);
        }
        if !row_status_size_ok {
            return self.index_status_ok_return(false);
        }
        let num_iis_col = self.col_index_.len() as HighsInt;
        let num_iis_row = self.row_index_.len() as HighsInt;
        // Determine whether this is an IIS or just an IS
        let mut true_iis = false;
        for i_col in 0..num_col {
            if self.col_status_[i_col as usize] == kIisStatusInConflict {
                true_iis = true;
                break;
            }
        }
        if !true_iis {
            for i_row in 0..num_row {
                if self.row_status_[i_row as usize] == kIisStatusInConflict {
                    true_iis = true;
                    break;
                }
            }
        }
        // Now check that cols and rows in the IIS are kIisStatusInConflict
        // or kIisStatusMaybeInConflict, according to true_iis, and that all
        // other cols and rows are kIisStatusNotConflict
        let mut col_status = self.col_status_.clone();
        let mut row_status = self.row_status_.clone();
        let illegal_status: HighsInt = -99;
        for i_x in 0..num_iis_col {
            let i_col = self.col_index_[i_x as usize];
            if self.col_status_[i_col as usize]
                != if true_iis {
                    kIisStatusInConflict
                } else {
                    kIisStatusMaybeInConflict
                }
            {
                return self.index_status_ok_return(false);
            }
            col_status[i_col as usize] = illegal_status;
        }
        for i_x in 0..num_iis_row {
            let i_row = self.row_index_[i_x as usize];
            if self.row_status_[i_row as usize]
                != if true_iis {
                    kIisStatusInConflict
                } else {
                    kIisStatusMaybeInConflict
                }
            {
                return self.index_status_ok_return(false);
            }
            row_status[i_row as usize] = illegal_status;
        }
        let non_is_status = self.non_is_status();
        for i_col in 0..num_col {
            if col_status[i_col as usize] > illegal_status
                && col_status[i_col as usize] != non_is_status
            {
                return self.index_status_ok_return(false);
            }
        }
        for i_row in 0..num_row {
            if row_status[i_row as usize] > illegal_status
                && row_status[i_row as usize] != non_is_status
            {
                return self.index_status_ok_return(false);
            }
        }
        self.index_status_ok_return(true)
    }

    /// HiGHS: HighsIis.cpp::HighsIis::lpDataOk (lines 1042-1197)
    pub fn lp_data_ok(&self, lp: &HighsLp, options: &HighsOptions) -> bool {
        let iis_lp = &self.model_.lp_;
        let iis_num_col = self.col_index_.len() as HighsInt;
        let iis_num_row = self.row_index_.len() as HighsInt;
        if !(iis_lp.num_col_ == iis_num_col) {
            return self.lp_data_ok_return(false);
        }
        if !(iis_lp.num_row_ == iis_num_row) {
            return self.lp_data_ok_return(false);
        }

        let colwise = lp.a_matrix_.is_colwise();

        let illegal_index: HighsInt = -1;
        let illegal_value: f64 = kHighsInf;
        // iis_row/col give the row/col in the IIS for each row/col in the
        // LP, or an illegal index if the LP row/col isn't in the IIS
        let mut iis_row: Vec<HighsInt>;
        iis_row = vec![illegal_index; lp.num_row_ as usize];
        let mut iis_col: Vec<HighsInt>;
        iis_col = vec![illegal_index; lp.num_col_ as usize];
        let mut bound: f64;
        for iis_row_idx in 0..iis_num_row {
            let i_row = self.row_index_[iis_row_idx as usize];
            iis_row[i_row as usize] = iis_row_idx;
            let row_bound = self.row_bound_[iis_row_idx as usize];
            bound = if row_bound == kIisBoundStatusLower || row_bound == kIisBoundStatusBoxed {
                lp.row_lower_[i_row as usize]
            } else {
                -kHighsInf
            };
            if iis_lp.row_lower_[iis_row_idx as usize] != bound {
                return self.lp_data_ok_return(false);
            }
            bound = if row_bound == kIisBoundStatusUpper || row_bound == kIisBoundStatusBoxed {
                lp.row_upper_[i_row as usize]
            } else {
                kHighsInf
            };
            if iis_lp.row_upper_[iis_row_idx as usize] != bound {
                return self.lp_data_ok_return(false);
            }
        }

        // Work through the LP columns checking the zero costs and bounds
        for iis_col_idx in 0..iis_num_col {
            let i_col = self.col_index_[iis_col_idx as usize];
            iis_col[i_col as usize] = iis_col_idx;
            if iis_lp.col_cost_[iis_col_idx as usize] != 0.0 {
                return self.lp_data_ok_return(false);
            }
            let col_bound = self.col_bound_[iis_col_idx as usize];
            bound = if col_bound == kIisBoundStatusLower || col_bound == kIisBoundStatusBoxed {
                lp.col_lower_[i_col as usize]
            } else {
                -kHighsInf
            };
            if iis_lp.col_lower_[iis_col_idx as usize] != bound {
                return self.lp_data_ok_return(false);
            }
            bound = if col_bound == kIisBoundStatusUpper || col_bound == kIisBoundStatusBoxed {
                lp.col_upper_[i_col as usize]
            } else {
                kHighsInf
            };
            if iis_lp.col_upper_[iis_col_idx as usize] != bound {
                return self.lp_data_ok_return(false);
            }
        }
        let mut index: Vec<HighsInt>;
        let mut value: Vec<f64>;
        // Work through the LP matrix, checking the matrix index/value
        if colwise {
            for iis_col_idx in 0..iis_num_col {
                let i_col = self.col_index_[iis_col_idx as usize];
                // Use index/value to scatter the IIS matrix column
                index = vec![illegal_index; iis_num_row as usize];
                value = vec![illegal_value; iis_num_row as usize];
                for i_el in iis_lp.a_matrix_.start_[iis_col_idx as usize]
                    ..iis_lp.a_matrix_.start_[(iis_col_idx + 1) as usize]
                {
                    let iis_row_idx = iis_lp.a_matrix_.index_[i_el as usize];
                    let i_row = self.row_index_[iis_row_idx as usize];
                    index[iis_row_idx as usize] = i_row;
                    value[iis_row_idx as usize] = iis_lp.a_matrix_.value_[i_el as usize];
                }
                for i_el in
                    lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
                {
                    let i_row = lp.a_matrix_.index_[i_el as usize];
                    let iis_row_idx = iis_row[i_row as usize];
                    if iis_row_idx >= 0 {
                        if index[iis_row_idx as usize] != i_row {
                            return self.lp_data_ok_return(false);
                        }
                        if value[iis_row_idx as usize] != lp.a_matrix_.value_[i_el as usize] {
                            return self.lp_data_ok_return(false);
                        }
                        index[iis_row_idx as usize] = illegal_index;
                        value[iis_row_idx as usize] = illegal_value;
                    }
                }
            }
        } else {
            for iis_row_idx in 0..iis_num_row {
                let i_row = self.row_index_[iis_row_idx as usize];
                // Use index/value to scatter the IIS matrix row
                index = vec![illegal_index; iis_num_col as usize];
                value = vec![illegal_value; iis_num_col as usize];
                for i_el in iis_lp.a_matrix_.start_[iis_row_idx as usize]
                    ..iis_lp.a_matrix_.start_[(iis_row_idx + 1) as usize]
                {
                    let iis_col_idx = iis_lp.a_matrix_.index_[i_el as usize];
                    let i_col = self.col_index_[iis_col_idx as usize];
                    index[iis_col_idx as usize] = i_col;
                    value[iis_col_idx as usize] = iis_lp.a_matrix_.value_[i_el as usize];
                }
                for i_el in
                    lp.a_matrix_.start_[i_row as usize]..lp.a_matrix_.start_[(i_row + 1) as usize]
                {
                    let i_col = lp.a_matrix_.index_[i_el as usize];
                    let iis_col_idx = iis_col[i_col as usize];
                    if iis_col_idx >= 0 {
                        if index[iis_col_idx as usize] != i_col {
                            return self.lp_data_ok_return(false);
                        }
                        if value[iis_col_idx as usize] != lp.a_matrix_.value_[i_el as usize] {
                            return self.lp_data_ok_return(false);
                        }
                        index[iis_col_idx as usize] = illegal_index;
                        value[iis_col_idx as usize] = illegal_value;
                    }
                }
            }
        }
        // Work through the IIS LP matrix, making sure that the index/value
        // are correct
        if colwise {
            for iis_col_idx in 0..iis_num_col {
                let i_col = self.col_index_[iis_col_idx as usize];
                // Use index/value to scatter the LP matrix column
                index = vec![illegal_index; lp.num_row_ as usize];
                value = vec![illegal_value; lp.num_row_ as usize];
                for i_el in
                    lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
                {
                    let i_row = lp.a_matrix_.index_[i_el as usize];
                    let iis_row_idx = iis_row[i_row as usize];
                    index[i_row as usize] = iis_row_idx;
                    value[i_row as usize] = lp.a_matrix_.value_[i_el as usize];
                }
                for i_el in iis_lp.a_matrix_.start_[iis_col_idx as usize]
                    ..iis_lp.a_matrix_.start_[(iis_col_idx + 1) as usize]
                {
                    let iis_row_idx = iis_lp.a_matrix_.index_[i_el as usize];
                    let i_row = self.row_index_[iis_row_idx as usize];
                    if index[i_row as usize] != iis_row_idx {
                        return self.lp_data_ok_return(false);
                    }
                    if value[i_row as usize] != iis_lp.a_matrix_.value_[i_el as usize] {
                        return self.lp_data_ok_return(false);
                    }
                }
            }
        } else {
            for iis_row_idx in 0..iis_num_row {
                let i_row = self.row_index_[iis_row_idx as usize];
                // Use index/value to scatter the LP matrix row
                index = vec![illegal_index; lp.num_col_ as usize];
                value = vec![illegal_value; lp.num_col_ as usize];
                for i_el in
                    lp.a_matrix_.start_[i_row as usize]..lp.a_matrix_.start_[(i_row + 1) as usize]
                {
                    let i_col = lp.a_matrix_.index_[i_el as usize];
                    let iis_col_idx = iis_col[i_col as usize];
                    index[i_col as usize] = iis_col_idx;
                    value[i_col as usize] = lp.a_matrix_.value_[i_el as usize];
                }
                for i_el in iis_lp.a_matrix_.start_[iis_row_idx as usize]
                    ..iis_lp.a_matrix_.start_[(iis_row_idx + 1) as usize]
                {
                    let iis_col_idx = iis_lp.a_matrix_.index_[i_el as usize];
                    let i_col = self.col_index_[iis_col_idx as usize];
                    if index[i_col as usize] != iis_col_idx {
                        return self.lp_data_ok_return(false);
                    }
                    if value[i_col as usize] != iis_lp.a_matrix_.value_[i_el as usize] {
                        return self.lp_data_ok_return(false);
                    }
                }
            }
        }
        self.lp_data_ok_return(true)
    }

    /// HiGHS: HighsIis.cpp::HighsIis::lpOk (lines 1199-1290)
    pub fn lp_ok(&self, options: &HighsOptions) -> bool {
        // Check that the IIS LP is OK (infeasible and optimal if
        // any bound is relaxed)
        if !self.valid_ {
            return self.lp_ok_return(true);
        }
        let num_iis_col = self.col_index_.len() as HighsInt;
        let num_iis_row = self.row_index_.len() as HighsInt;
        // If an LP has a row with inconsistent bounds, or an empty row with
        // a positive lower bound or negative upper bound, then it is
        // infeasible, but the IIS contains no columns
        if num_iis_col == 0 {
            return true;
        }
        let log_options = &options.log_options;
        let iis_lp = &self.model_.lp_;
        debug_assert!(iis_lp.num_col_ == num_iis_col);
        debug_assert!(iis_lp.num_row_ == num_iis_row);
        let mut h = Highs::default();
        h.pass_options(options);
        h.set_option_value_bool("output_flag", false);
        h.pass_model(iis_lp.clone());
        // HighsIis.cpp:1217 h.writeModel("") dropped (file writers not ported).
        let status = h.optimize_model();
        if status != HighsStatus::kOk || h.get_model_status() != HighsModelStatus::kInfeasible {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kWarning,
                "HighsIis: Failed to prove infeasibility for IIS LP\n"
            );
            return self.lp_ok_return(false);
        }
        if !(self.status_ == kIisModelStatusIrreducible) {
            return self.lp_ok_return(true);
        }
        let optimal = |h: &mut Highs| -> bool {
            if options.log_dev_level > 0 {
                // HighsIis.cpp:1227 h.writeModel("") dropped (file writers not ported).
            }
            h.optimize_model();
            h.get_model_status() == HighsModelStatus::kOptimal
        };
        for iis_col in 0..num_iis_col {
            let i_col = self.col_index_[iis_col as usize];
            if self.col_bound_[iis_col as usize] == kIisBoundStatusLower {
                h.change_col_bounds(iis_col, -kHighsInf, iis_lp.col_upper_[iis_col as usize]);
                if !optimal(&mut h) {
                    crate::highs_log_user!(
                        log_options,
                        HighsLogType::kWarning,
                        "HighsIis: IIS column {} (LP column {}): relaxing lower bound of {} yield IIS LP with status {}\n",
                        iis_col,
                        i_col,
                        c_fmt_g(iis_lp.col_lower_[iis_col as usize], 6),
                        h.model_status_to_string(h.get_model_status())
                    );
                    return self.lp_ok_return(false);
                }
                h.change_col_bounds(
                    iis_col,
                    iis_lp.col_lower_[iis_col as usize],
                    iis_lp.col_upper_[iis_col as usize],
                );
            }
            if self.col_bound_[iis_col as usize] == kIisBoundStatusUpper {
                h.change_col_bounds(iis_col, iis_lp.col_lower_[iis_col as usize], kHighsInf);
                if !optimal(&mut h) {
                    crate::highs_log_user!(
                        log_options,
                        HighsLogType::kWarning,
                        "HighsIis: IIS column {} (LP column {}): relaxing upper bound of {} yield IIS LP with status {}\n",
                        iis_col,
                        i_col,
                        c_fmt_g(iis_lp.col_upper_[iis_col as usize], 6),
                        h.model_status_to_string(h.get_model_status())
                    );
                    return self.lp_ok_return(false);
                }
                h.change_col_bounds(
                    iis_col,
                    iis_lp.col_lower_[iis_col as usize],
                    iis_lp.col_upper_[iis_col as usize],
                );
            }
        }
        for iis_row in 0..num_iis_row {
            let i_row = self.row_index_[iis_row as usize];
            if self.row_bound_[iis_row as usize] == kIisBoundStatusLower {
                h.change_row_bounds(iis_row, -kHighsInf, iis_lp.row_upper_[iis_row as usize]);
                if !optimal(&mut h) {
                    crate::highs_log_user!(
                        log_options,
                        HighsLogType::kError,
                        "HighsIis: IIS row {} (LP row {}): relaxing lower bound of {} yield IIS LP with status {}\n",
                        iis_row,
                        i_row,
                        c_fmt_g(iis_lp.row_lower_[iis_row as usize], 6),
                        h.model_status_to_string(h.get_model_status())
                    );
                    return self.lp_ok_return(false);
                }
                h.change_row_bounds(
                    iis_row,
                    iis_lp.row_lower_[iis_row as usize],
                    iis_lp.row_upper_[iis_row as usize],
                );
            }
            if self.row_bound_[iis_row as usize] == kIisBoundStatusUpper {
                h.change_row_bounds(iis_row, iis_lp.row_lower_[iis_row as usize], kHighsInf);
                if !optimal(&mut h) {
                    crate::highs_log_user!(
                        log_options,
                        HighsLogType::kError,
                        "HighsIis: IIS row {} (LP row {}): relaxing upper bound of {} yield IIS LP with status {}\n",
                        iis_row,
                        i_row,
                        c_fmt_g(iis_lp.row_upper_[iis_row as usize], 6),
                        h.model_status_to_string(h.get_model_status())
                    );
                    return self.lp_ok_return(false);
                }
                h.change_row_bounds(
                    iis_row,
                    iis_lp.row_lower_[iis_row as usize],
                    iis_lp.row_upper_[iis_row as usize],
                );
            }
        }
        self.lp_ok_return(true)
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

    fn lp21() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp
    }

    #[test]
    fn index_status_consistency() {
        let lp = lp21();
        let mut iis = HighsIis::default();
        iis.status_ = kIisModelStatusIrreducible;
        iis.add_col(0, kIisBoundStatusLower);
        iis.add_row(0, kIisBoundStatusUpper);
        iis.col_status_ = vec![kIisStatusInConflict, kIisStatusNotInConflict];
        iis.row_status_ = vec![kIisStatusInConflict];
        assert!(iis.index_status_ok(&lp));
        iis.col_status_[1] = kIisStatusMaybeInConflict;
        assert!(!iis.index_status_ok(&lp));
    }

    fn index_status_wrong_size() -> bool {
        let lp = lp21();
        let mut iis = HighsIis::default();
        iis.status_ = kIisModelStatusIrreducible;
        iis.col_status_ = vec![kIisStatusInConflict];
        iis.row_status_ = vec![kIisStatusInConflict];
        iis.index_status_ok(&lp)
    }

    /// The size checks are C++ `assert`s before the early return (HighsIis.cpp:988-991).
    #[test]
    #[cfg(not(debug_assertions))]
    fn index_status_wrong_size_is_not_ok() {
        assert!(!index_status_wrong_size());
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "col_status_size_ok")]
    fn index_status_wrong_size_trips_cpp_assert() {
        index_status_wrong_size();
    }

    // LP: 2 cols x 2 rows, A = [1 0; 2 3]; IIS = {col 0 (lower), row 1 (upper)}
    fn lp_data_case(colwise: bool) -> (HighsLp, HighsIis) {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![1.0, 0.0];
        lp.col_upper_ = vec![4.0, 5.0];
        lp.row_lower_ = vec![-1.0, -2.0];
        lp.row_upper_ = vec![6.0, 1.0];
        let mut iis = HighsIis::default();
        iis.add_col(0, kIisBoundStatusLower);
        iis.add_row(1, kIisBoundStatusUpper);
        let mut iis_lp = HighsLp::default();
        iis_lp.num_col_ = 1;
        iis_lp.num_row_ = 1;
        iis_lp.col_cost_ = vec![0.0];
        iis_lp.col_lower_ = vec![1.0];
        iis_lp.col_upper_ = vec![kHighsInf];
        iis_lp.row_lower_ = vec![-kHighsInf];
        iis_lp.row_upper_ = vec![1.0];
        if colwise {
            lp.a_matrix_.start_ = vec![0, 2, 3];
            lp.a_matrix_.index_ = vec![0, 1, 1];
            lp.a_matrix_.value_ = vec![1.0, 2.0, 3.0];
        } else {
            lp.a_matrix_.format_ = MatrixFormat::kRowwise;
            lp.a_matrix_.start_ = vec![0, 1, 3];
            lp.a_matrix_.index_ = vec![0, 0, 1];
            lp.a_matrix_.value_ = vec![1.0, 2.0, 3.0];
            iis_lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        }
        iis_lp.a_matrix_.start_ = vec![0, 1];
        iis_lp.a_matrix_.index_ = vec![0];
        iis_lp.a_matrix_.value_ = vec![2.0];
        iis.model_.lp_ = iis_lp;
        (lp, iis)
    }

    #[test]
    fn lp_data_consistency() {
        let options = HighsOptions::default();
        for colwise in [true, false] {
            let (lp, mut iis) = lp_data_case(colwise);
            assert!(iis.lp_data_ok(&lp, &options));
            // Wrong matrix value
            iis.model_.lp_.a_matrix_.value_[0] = 3.0;
            assert!(!iis.lp_data_ok(&lp, &options));
            iis.model_.lp_.a_matrix_.value_[0] = 2.0;
            // Nonzero cost (C++ `if (iis_lp.col_cost_[iisCol])`)
            iis.model_.lp_.col_cost_[0] = -0.5;
            assert!(!iis.lp_data_ok(&lp, &options));
            iis.model_.lp_.col_cost_[0] = 0.0;
            // Upper bound kept although only the lower bound is in the IIS
            iis.model_.lp_.col_upper_[0] = 4.0;
            assert!(!iis.lp_data_ok(&lp, &options));
            iis.model_.lp_.col_upper_[0] = kHighsInf;
            // Row lower bound kept although only the upper bound is in the IIS
            iis.model_.lp_.row_lower_[0] = -2.0;
            assert!(!iis.lp_data_ok(&lp, &options));
            iis.model_.lp_.row_lower_[0] = -kHighsInf;
            // IIS dimensions inconsistent with col_index_
            iis.model_.lp_.num_col_ = 2;
            assert!(!iis.lp_data_ok(&lp, &options));
        }
    }

    #[test]
    fn lp_ok_trivial_returns() {
        let options = HighsOptions::default();
        let mut iis = HighsIis::default();
        // Not valid: lpOkReturn(true) without solving
        iis.add_col(0, kIisBoundStatusLower);
        assert!(iis.lp_ok(&options));
        // Valid but no IIS columns: true without solving
        let mut iis = HighsIis::default();
        iis.valid_ = true;
        iis.add_row(0, kIisBoundStatusLower);
        assert!(iis.lp_ok(&options));
    }
}
