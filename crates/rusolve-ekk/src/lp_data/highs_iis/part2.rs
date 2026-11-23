//! HiGHS source: lp_data/HighsIis.cpp (lines 368-591)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;
use crate::prelude::*;

impl HighsIis {
    /// HiGHS: HighsIis.cpp::HighsIis::deduce (lines 368-435)
    pub fn deduce(
        &mut self,
        lp: &HighsLp,
        options: &HighsOptions,
        callback: &HighsCallback,
        basis: &HighsBasis,
    ) -> HighsStatus {
        // The number of infeasible rows must be positive
        debug_assert!(!self.row_index_.is_empty());
        // Identify the LP corresponding to the set of infeasible rows
        let mut from_row = self.row_index_.clone();
        let mut from_col: Vec<HighsInt> = Vec::new();
        let mut to_row: Vec<HighsInt>;
        to_row = vec![-1; lp.num_row_ as usize];
        // Only uses this->row_index_ to initialise from_row, so can clear
        self.clear_data();
        // ToDo Exploit the known col_index_ and row_bound_ HighsIis
        // information
        //
        // To get the IIS data needs the matrix to be column-wise
        debug_assert!(lp.a_matrix_.is_colwise());
        // Determine how to detect whether a row is in from_row and (then)
        // gather information about it
        for i_x in 0..from_row.len() as HighsInt {
            to_row[from_row[i_x as usize] as usize] = i_x;
        }
        // Identify the columns (from_col) with nonzeros in the infeasible
        // rows
        for i_col in 0..lp.num_col_ {
            let mut use_col = false;
            for i_el in
                lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
            {
                use_col = use_col || to_row[lp.a_matrix_.index_[i_el as usize] as usize] >= 0;
            }
            if use_col {
                from_col.push(i_col);
            }
        }
        let to_num_col = from_col.len() as HighsInt;
        let to_num_row = from_row.len() as HighsInt;
        let mut to_lp = HighsLp::default();
        to_lp.num_col_ = to_num_col;
        to_lp.num_row_ = to_num_row;
        to_lp.a_matrix_.num_col_ = to_lp.num_col_;
        to_lp.a_matrix_.num_row_ = to_lp.num_row_;
        let has_col_names = !lp.col_names_.is_empty();
        for i_col in 0..to_num_col {
            to_lp.col_cost_.push(0.0);
            to_lp
                .col_lower_
                .push(lp.col_lower_[from_col[i_col as usize] as usize]);
            to_lp
                .col_upper_
                .push(lp.col_upper_[from_col[i_col as usize] as usize]);
            if has_col_names {
                to_lp
                    .col_names_
                    .push(lp.col_names_[from_col[i_col as usize] as usize].clone());
            }
            for i_el in lp.a_matrix_.start_[from_col[i_col as usize] as usize]
                ..lp.a_matrix_.start_[(from_col[i_col as usize] + 1) as usize]
            {
                let i_row = lp.a_matrix_.index_[i_el as usize];
                if to_row[i_row as usize] >= 0 {
                    to_lp.a_matrix_.index_.push(to_row[i_row as usize]);
                    to_lp
                        .a_matrix_
                        .value_
                        .push(lp.a_matrix_.value_[i_el as usize]);
                }
            }
            to_lp
                .a_matrix_
                .start_
                .push(to_lp.a_matrix_.index_.len() as HighsInt);
        }
        let has_row_names = !lp.row_names_.is_empty();
        for i_row in 0..to_num_row {
            to_lp
                .row_lower_
                .push(lp.row_lower_[from_row[i_row as usize] as usize]);
            to_lp
                .row_upper_
                .push(lp.row_upper_[from_row[i_row as usize] as usize]);
            if has_row_names {
                to_lp
                    .row_names_
                    .push(lp.row_names_[from_row[i_row as usize] as usize].clone());
            }
        }
        let return_status = self.compute(&to_lp, options, callback, None);
        // Indirect the values into the original LP
        for colindex in self.col_index_.iter_mut() {
            *colindex = from_col[*colindex as usize];
        }
        for rowindex in self.row_index_.iter_mut() {
            *rowindex = from_row[*rowindex as usize];
        }
        if kIisDevReport {
            self.report("On exit", lp);
        }
        return_status
    }

    /// HiGHS: HighsIis.cpp::HighsIis::setLp (lines 437-533)
    pub fn set_lp(&mut self, lp: &HighsLp) {
        let iis_lp = &mut self.model_.lp_;
        iis_lp.clear();
        let iis_num_col = self.col_index_.len() as HighsInt;
        let iis_num_row = self.row_index_.len() as HighsInt;
        let colwise = lp.a_matrix_.is_colwise();
        // Scatter the IIS rows (cols) into a full-length vector to identify
        // IIS rows (cols) with LP rows (cols) according to whether the
        // incumbent matrix is col-wise or row-wise
        let mut iis_row: Vec<HighsInt> = Vec::new();
        let mut iis_col: Vec<HighsInt> = Vec::new();
        if colwise {
            iis_row = vec![-1; lp.num_row_ as usize];
            for iis_row_idx in 0..iis_num_row {
                iis_row[self.row_index_[iis_row_idx as usize] as usize] = iis_row_idx;
            }
        } else {
            iis_col = vec![-1; lp.num_col_ as usize];
            for iis_col_idx in 0..iis_num_col {
                iis_col[self.col_index_[iis_col_idx as usize] as usize] = iis_col_idx;
            }
        }
        let mut bound: f64;

        let has_row_name = !lp.row_names_.is_empty();
        for iis_row_idx in 0..iis_num_row {
            let i_row = self.row_index_[iis_row_idx as usize];
            if has_row_name {
                iis_lp
                    .row_names_
                    .push(lp.row_names_[i_row as usize].clone());
            }
            let row_bound = self.row_bound_[iis_row_idx as usize];
            debug_assert!(
                row_bound == kIisBoundStatusLower
                    || row_bound == kIisBoundStatusUpper
                    || row_bound == kIisBoundStatusBoxed
            );
            bound = if row_bound == kIisBoundStatusLower || row_bound == kIisBoundStatusBoxed {
                lp.row_lower_[i_row as usize]
            } else {
                -kHighsInf
            };
            iis_lp.row_lower_.push(bound);
            bound = if row_bound == kIisBoundStatusUpper || row_bound == kIisBoundStatusBoxed {
                lp.row_upper_[i_row as usize]
            } else {
                kHighsInf
            };
            iis_lp.row_upper_.push(bound);
            if !colwise {
                for i_el in
                    lp.a_matrix_.start_[i_row as usize]..lp.a_matrix_.start_[(i_row + 1) as usize]
                {
                    let i_col = lp.a_matrix_.index_[i_el as usize];
                    let iis_col_idx = iis_col[i_col as usize];
                    if iis_col_idx >= 0 {
                        iis_lp.a_matrix_.index_.push(iis_col_idx);
                        iis_lp
                            .a_matrix_
                            .value_
                            .push(lp.a_matrix_.value_[i_el as usize]);
                    }
                }
            }
        }

        let has_col_name = !lp.col_names_.is_empty();
        for iis_col_idx in 0..iis_num_col {
            let i_col = self.col_index_[iis_col_idx as usize];
            // Costs in the IIS LP are zero since they play no role in IIS,
            // and when dropping bounds, optimality is the only valid model
            // status
            iis_lp.col_cost_.push(0.0);
            if has_col_name {
                iis_lp
                    .col_names_
                    .push(lp.col_names_[i_col as usize].clone());
            }
            let col_bound = self.col_bound_[iis_col_idx as usize];
            debug_assert!(
                col_bound == kIisBoundStatusLower
                    || col_bound == kIisBoundStatusUpper
                    || col_bound == kIisBoundStatusBoxed
                    || col_bound == kIisBoundStatusFree
            );
            bound = if col_bound == kIisBoundStatusLower || col_bound == kIisBoundStatusBoxed {
                lp.col_lower_[i_col as usize]
            } else {
                -kHighsInf
            };
            iis_lp.col_lower_.push(bound);
            bound = if col_bound == kIisBoundStatusUpper || col_bound == kIisBoundStatusBoxed {
                lp.col_upper_[i_col as usize]
            } else {
                kHighsInf
            };
            iis_lp.col_upper_.push(bound);
            if colwise {
                for i_el in
                    lp.a_matrix_.start_[i_col as usize]..lp.a_matrix_.start_[(i_col + 1) as usize]
                {
                    let i_row = lp.a_matrix_.index_[i_el as usize];
                    let iis_row_idx = iis_row[i_row as usize];
                    if iis_row_idx >= 0 {
                        iis_lp.a_matrix_.index_.push(iis_row_idx);
                        iis_lp
                            .a_matrix_
                            .value_
                            .push(lp.a_matrix_.value_[i_el as usize]);
                    }
                }
            }
            iis_lp
                .a_matrix_
                .start_
                .push(iis_lp.a_matrix_.index_.len() as HighsInt);
        }
        iis_lp.num_col_ = iis_lp.col_cost_.len() as HighsInt;
        iis_lp.num_row_ = iis_lp.row_lower_.len() as HighsInt;
        // The IIS LP matrix will have the same format as the incumbent LP
        iis_lp.a_matrix_.format_ = lp.a_matrix_.format_;
        iis_lp.a_matrix_.num_col_ = iis_lp.num_col_;
        iis_lp.a_matrix_.num_row_ = iis_lp.num_row_;
        iis_lp.model_name_ = format!("{}_IIS", lp.model_name_);
    }

    /// HiGHS: HighsIis.cpp::HighsIis::nonIsStatus (lines 535-548)
    pub fn non_is_status(&self) -> HighsInt {
        let is_feasible = self.status_ == kIisModelStatusFeasible;
        let has_is = !self.col_index_.is_empty() || !self.row_index_.is_empty();
        // If the model is known to be feasible, then there should be no IS,
        // and all columns and rows are kIisStatusNotInConflict
        if is_feasible {
            debug_assert!(!has_is);
        }
        if has_is {
            debug_assert!(self.status_ >= kIisModelStatusTimeLimit);
        }
        // If there is an IS, then all columns and rows not in the IS are
        // kIisStatusNotInConflict
        if is_feasible || has_is {
            kIisStatusNotInConflict
        } else {
            kIisStatusMaybeInConflict
        }
    }

    /// HiGHS: HighsIis.cpp::HighsIis::setStatus (lines 550-564)
    pub fn set_status(&mut self, lp: &HighsLp) {
        if !self.valid_ {
            return;
        }
        let non_is_status = self.non_is_status();
        let in_is_status = if self.status_ == kIisModelStatusIrreducible {
            kIisStatusInConflict
        } else {
            kIisStatusMaybeInConflict
        };
        self.col_status_ = vec![non_is_status; lp.num_col_ as usize];
        self.row_status_ = vec![non_is_status; lp.num_row_ as usize];
        let iis_num_col = self.col_index_.len() as HighsInt;
        let iis_num_row = self.row_index_.len() as HighsInt;
        for iis_col in 0..iis_num_col {
            self.col_status_[self.col_index_[iis_col as usize] as usize] = in_is_status;
        }
        for iis_row in 0..iis_num_row {
            self.row_status_[self.row_index_[iis_row as usize] as usize] = in_is_status;
        }
    }

    /// HiGHS: HighsIis.cpp::HighsIis::determineBoundStatus (lines 566-591)
    pub fn determine_bound_status(&self, lower: f64, upper: f64, is_row: bool) -> HighsInt {
        let iss_bound_status: HighsInt;
        if lower <= -kHighsInf {
            if upper >= kHighsInf {
                if is_row {
                    // Free rows can be dropped
                    iss_bound_status = kIisBoundStatusDropped;
                } else {
                    // Free columns can only be dropped if they are empty
                    iss_bound_status = kIisBoundStatusFree;
                }
            } else {
                iss_bound_status = kIisBoundStatusUpper;
            }
        } else if upper >= kHighsInf {
            iss_bound_status = kIisBoundStatusLower;
        } else {
            // FX or BX
            iss_bound_status = kIisBoundStatusBoxed;
        }
        debug_assert!(iss_bound_status != kIisBoundStatusNull);
        iss_bound_status
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
    fn bound_status_classes() {
        let iis = HighsIis::default();
        assert_eq!(
            iis.determine_bound_status(-kHighsInf, kHighsInf, true),
            kIisBoundStatusDropped
        );
        assert_eq!(
            iis.determine_bound_status(-kHighsInf, kHighsInf, false),
            kIisBoundStatusFree
        );
        assert_eq!(
            iis.determine_bound_status(-kHighsInf, 1.0, false),
            kIisBoundStatusUpper
        );
        assert_eq!(
            iis.determine_bound_status(0.0, kHighsInf, true),
            kIisBoundStatusLower
        );
        assert_eq!(
            iis.determine_bound_status(1.0, 1.0, false),
            kIisBoundStatusBoxed
        );
    }

    #[test]
    fn non_is_status_with_an_is() {
        let mut iis = HighsIis::default();
        iis.status_ = kIisModelStatusIrreducible;
        iis.add_row(0, kIisBoundStatusUpper);
        assert_eq!(iis.non_is_status(), kIisStatusNotInConflict);
    }
}
