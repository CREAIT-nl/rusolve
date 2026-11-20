//! HiGHS source: lp_data/HighsRanging.h, lp_data/HighsRanging.cpp (lines 14-68),
//! lp_data/Highs.cpp (Highs::getRanging), lp_data/HighsInterface.cpp
//! (Highs::getRangingInterface)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `writeRangingFile` (HighsRanging.cpp 603-end) is file-writer code and is
//! deferred. `getRangingData` lives in `get_ranging_data.rs`.

pub mod get_ranging_data;
pub use get_ranging_data::*;

use crate::prelude::*;

/// HiGHS: HighsRanging.h::HighsRangingRecord (lines 17-22)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HighsRangingRecord {
    pub value_: Vec<f64>,
    pub objective_: Vec<f64>,
    pub in_var_: Vec<HighsInt>,
    pub ou_var_: Vec<HighsInt>,
}

/// HiGHS: HighsRanging.h::HighsRanging (lines 24-34)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HighsRanging {
    pub valid: bool,
    pub col_cost_up: HighsRangingRecord,
    pub col_cost_dn: HighsRangingRecord,
    pub col_bound_up: HighsRangingRecord,
    pub col_bound_dn: HighsRangingRecord,
    pub row_bound_up: HighsRangingRecord,
    pub row_bound_dn: HighsRangingRecord,
}

impl HighsRanging {
    /// HiGHS: HighsRanging.cpp::HighsRanging::invalidate (line 22)
    pub fn invalidate(&mut self) {
        self.valid = false;
    }

    /// HiGHS: HighsRanging.cpp::HighsRanging::clear (lines 24-50)
    pub fn clear(&mut self) {
        self.invalidate();
        for record in [
            &mut self.col_cost_up,
            &mut self.col_cost_dn,
            &mut self.col_bound_up,
            &mut self.col_bound_dn,
            &mut self.row_bound_up,
            &mut self.row_bound_dn,
        ] {
            record.value_.clear();
            record.objective_.clear();
            record.in_var_.clear();
            record.ou_var_.clear();
        }
    }
}

/// HiGHS: HighsRanging.cpp::infProduct (lines 52-59)
pub(crate) fn inf_product(value: f64) -> f64 {
    // Multiplying value and kHighsInf
    if value == 0.0 { 0.0 } else { value * kHighsInf }
}

/// HiGHS: HighsRanging.cpp::possInfProduct (lines 61-68)
pub(crate) fn poss_inf_product(poss_inf: f64, value: f64) -> f64 {
    // Multiplying something that could be infinite and value
    if value == 0.0 { 0.0 } else { poss_inf * value }
}

impl Highs {
    /// HiGHS: Highs.cpp::Highs::getRanging (lines 2191-2199)
    ///
    /// PORT NOTE: the single-threaded profiling scaffolding
    /// (`initializeSingleThreadedProfiling` / `clearProfiling`) is a no-op in
    /// this port, as elsewhere.
    pub fn get_ranging(&mut self, ranging: &mut HighsRanging) -> HighsStatus {
        // HighsProfiling profiling; (line 2192) is a no-op in this port.
        // const bool already_profiling = this->profiling_; (line 2193) is a
        // no-op in this port.
        // if (!already_profiling) this->initializeSingleThreadedProfiling(&profiling);
        // (line 2194) is a no-op in this port.
        let return_status = self.get_ranging_interface();
        // if (!already_profiling) this->clearProfiling(); (line 2196) is a
        // no-op in this port.
        *ranging = self.ranging_.clone();
        return_status
    }

    /// HiGHS: HighsInterface.cpp::Highs::getRangingInterface (lines 1825-1831)
    ///
    /// Build a `HighsLpSolverObject` over `model_lp_` exactly as
    /// `call_solve_lp` does, copy `model_status_` into it, and call
    /// `get_ranging_data(&mut self.ranging_, &mut solver_object, &mut
    /// self.ekk_instance_)` (the `HEkk` is passed separately, lib.rs rule e).
    pub(crate) fn get_ranging_interface(&mut self) -> HighsStatus {
        let mut solver_object = HighsLpSolverObject::new(
            &mut self.model_lp_,
            &mut self.basis_,
            &mut self.solution_,
            &mut self.info_,
            &mut self.callback_,
            &mut self.options_,
            &mut self.timer_,
        );
        // solver_object.setProfiling(this->profiling_); (line 1828) is a no-op
        // in this port.
        solver_object.model_status_ = self.model_status_;
        get_ranging_data(
            &mut self.ranging_,
            &mut solver_object,
            &mut self.ekk_instance_,
        )
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

    fn ranging_lp() -> HighsLp {
        lp(
            &[
                (1.0, 0.0, kHighsInf, &[(0, 1.0)]),
                (2.0, 0.0, kHighsInf, &[(0, 1.0)]),
            ],
            &[(1.0, kHighsInf)],
        )
    }

    #[test]
    fn ranging_needs_a_solved_lp() {
        let mut highs = quiet(ranging_lp());
        let mut r = HighsRanging::default();
        assert_eq!(highs.get_ranging(&mut r), HighsStatus::kError);
    }

    #[test]
    fn cost_ranging_two_columns() {
        let mut highs = quiet(ranging_lp());
        assert_eq!(highs.run(), HighsStatus::kOk);
        let mut r = HighsRanging::default();
        assert_eq!(highs.get_ranging(&mut r), HighsStatus::kOk);
        assert!(r.valid);
        // C++ HighsRanging.cpp:109 numTotal = numCol + numRow (2 cols + 1 row);
        // the cost records are full-length vectors over numTotal, not numCol.
        assert_eq!(r.col_cost_up.value_.len(), 3);
        assert_eq!(r.row_bound_up.value_.len(), 1);
        // x0 basic: optimal until its cost reaches that of x1
        assert!((r.col_cost_up.value_[0] - 2.0).abs() <= 1e-9);
        // x1 nonbasic at lower: cost may rise without bound, fall to 1
        assert_eq!(r.col_cost_up.value_[1], kHighsInf);
        assert!((r.col_cost_dn.value_[1] - 1.0).abs() <= 1e-9);
    }
}
