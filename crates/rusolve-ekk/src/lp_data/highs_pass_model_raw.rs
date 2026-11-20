//! HiGHS source: lp_data/Highs.cpp (the two raw-array passModel overloads),
//! lp_data/HighsInterface.cpp (aFormatOk)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Building a model from raw CSC / CSR arrays, as the C API and
//! the oracle driver do. The `HighsLp` / `HighsModel` overloads are ported
//! (`pass_model`, highs_run.rs; `pass_highs_model`, highs_qp.rs) and
//! `qFormatOk` is highs_qp.rs. `logHeader` is dropped. Each nullable
//! `const T*` -> `Option<&[T]>` (an absent `costs`, bounds or `integrality`
//! is legal in C++ when the dimension is zero; `integrality` may be absent
//! for an LP).

use crate::prelude::*;

impl Highs {
    /// HiGHS: Highs.cpp::Highs::passModel (lines 533-640), the overload with
    /// a Hessian (`q_*`). `a_format` / `q_format` are `MatrixFormat` /
    /// `HessianFormat` as `HighsInt`; `sense` is `ObjSense` as `HighsInt`.
    #[allow(clippy::too_many_arguments)]
    pub fn pass_model_qp_raw(
        &mut self,
        num_col: HighsInt,
        num_row: HighsInt,
        a_num_nz: HighsInt,
        q_num_nz: HighsInt,
        a_format: HighsInt,
        q_format: HighsInt,
        sense: HighsInt,
        offset: f64,
        costs: Option<&[f64]>,
        col_lower: Option<&[f64]>,
        col_upper: Option<&[f64]>,
        row_lower: Option<&[f64]>,
        row_upper: Option<&[f64]>,
        a_start: Option<&[HighsInt]>,
        a_index: Option<&[HighsInt]>,
        a_value: Option<&[f64]>,
        q_start: Option<&[HighsInt]>,
        q_index: Option<&[HighsInt]>,
        q_value: Option<&[f64]>,
        integrality: Option<&[HighsInt]>,
    ) -> HighsStatus {
        // logHeader dropped (print-only)
        let mut model = HighsModel::default();
        let lp = &mut model.lp_;

        // Check that the formats of the constraint matrix and Hessian are valid
        if !self.a_format_ok(a_num_nz, a_format) {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Model has illegal constraint matrix format\n"
            );
            return HighsStatus::kError;
        }
        if !self.q_format_ok(q_num_nz, q_format) {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Model has illegal Hessian matrix format\n"
            );
            return HighsStatus::kError;
        }
        let a_rowwise = if a_num_nz > 0 {
            a_format == MatrixFormat::kRowwise as HighsInt
        } else {
            false
        };
        //  if (num_nz) a_rowwise = a_format == (HighsInt)MatrixFormat::kRowwise;

        lp.num_col_ = num_col;
        lp.num_row_ = num_row;
        if num_col > 0 {
            debug_assert!(costs.is_some());
            debug_assert!(col_lower.is_some());
            debug_assert!(col_upper.is_some());
            lp.col_cost_ = costs.unwrap()[..num_col as usize].to_vec();
            lp.col_lower_ = col_lower.unwrap()[..num_col as usize].to_vec();
            lp.col_upper_ = col_upper.unwrap()[..num_col as usize].to_vec();
        }
        if num_row > 0 {
            debug_assert!(row_lower.is_some());
            debug_assert!(row_upper.is_some());
            lp.row_lower_ = row_lower.unwrap()[..num_row as usize].to_vec();
            lp.row_upper_ = row_upper.unwrap()[..num_row as usize].to_vec();
        }
        if a_num_nz > 0 {
            debug_assert!(num_col > 0);
            debug_assert!(num_row > 0);
            debug_assert!(a_start.is_some());
            debug_assert!(a_index.is_some());
            debug_assert!(a_value.is_some());
            if a_rowwise {
                lp.a_matrix_.start_ = a_start.unwrap()[..num_row as usize].to_vec();
            } else {
                lp.a_matrix_.start_ = a_start.unwrap()[..num_col as usize].to_vec();
            }
            lp.a_matrix_.index_ = a_index.unwrap()[..a_num_nz as usize].to_vec();
            lp.a_matrix_.value_ = a_value.unwrap()[..a_num_nz as usize].to_vec();
        }
        if a_rowwise {
            lp.a_matrix_.start_.resize(num_row as usize + 1, 0);
            lp.a_matrix_.start_[num_row as usize] = a_num_nz;
            lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        } else {
            lp.a_matrix_.start_.resize(num_col as usize + 1, 0);
            lp.a_matrix_.start_[num_col as usize] = a_num_nz;
            lp.a_matrix_.format_ = MatrixFormat::kColwise;
        }
        if sense == ObjSense::kMaximize as HighsInt {
            lp.sense_ = ObjSense::kMaximize;
        } else {
            lp.sense_ = ObjSense::kMinimize;
        }
        lp.offset_ = offset;
        if let Some(integrality) = integrality.filter(|_| num_col > 0) {
            lp.integrality_ = vec![HighsVarType::kContinuous; num_col as usize];
            for i_col in 0..num_col {
                let integrality_status = integrality[i_col as usize];
                let legal_integrality_status = integrality_status
                    == HighsVarType::kContinuous as HighsInt
                    || integrality_status == HighsVarType::kInteger as HighsInt
                    || integrality_status == HighsVarType::kSemiContinuous as HighsInt
                    || integrality_status == HighsVarType::kSemiInteger as HighsInt
                    || integrality_status == HighsVarType::kImplicitInteger as HighsInt;
                if !legal_integrality_status {
                    crate::highs_log_user!(
                        &self.options_.log_options,
                        HighsLogType::kError,
                        "Model has illegal integer value of {} (type {}) for integrality[{}]\n",
                        integrality_status,
                        match integrality_status {
                            0 => "continuous",
                            1 => "integer",
                            2 => "semi continuous",
                            3 => "semi integer",
                            4 => "implicit integer",
                            _ => "unknown",
                        },
                        i_col
                    );
                    return HighsStatus::kError;
                }
                lp.integrality_[i_col as usize] = match integrality_status {
                    1 => HighsVarType::kInteger,
                    2 => HighsVarType::kSemiContinuous,
                    3 => HighsVarType::kSemiInteger,
                    4 => HighsVarType::kImplicitInteger,
                    _ => HighsVarType::kContinuous,
                };
            }
        }
        if q_num_nz > 0 {
            debug_assert!(num_col > 0);
            debug_assert!(q_start.is_some());
            debug_assert!(q_index.is_some());
            debug_assert!(q_value.is_some());
            let hessian = &mut model.hessian_;
            hessian.dim_ = num_col;
            hessian.format_ = if q_format == HessianFormat::kSquare as HighsInt {
                HessianFormat::kSquare
            } else {
                HessianFormat::kTriangular
            };
            hessian.start_ = q_start.unwrap()[..num_col as usize].to_vec();
            hessian.start_.resize(num_col as usize + 1, 0);
            hessian.start_[num_col as usize] = q_num_nz;
            hessian.index_ = q_index.unwrap()[..q_num_nz as usize].to_vec();
            hessian.value_ = q_value.unwrap()[..q_num_nz as usize].to_vec();
        }
        self.pass_highs_model(model)
    }

    /// HiGHS: Highs.cpp::Highs::passModel (lines 642-653), the LP overload:
    /// `pass_model_qp_raw` with `q_num_nz = 0`, `q_format = 0` and no
    /// Hessian arrays.
    #[allow(clippy::too_many_arguments)]
    pub fn pass_model_raw(
        &mut self,
        num_col: HighsInt,
        num_row: HighsInt,
        num_nz: HighsInt,
        a_format: HighsInt,
        sense: HighsInt,
        offset: f64,
        costs: Option<&[f64]>,
        col_lower: Option<&[f64]>,
        col_upper: Option<&[f64]>,
        row_lower: Option<&[f64]>,
        row_upper: Option<&[f64]>,
        a_start: Option<&[HighsInt]>,
        a_index: Option<&[HighsInt]>,
        a_value: Option<&[f64]>,
        integrality: Option<&[HighsInt]>,
    ) -> HighsStatus {
        self.pass_model_qp_raw(
            num_col,
            num_row,
            num_nz,
            0,
            a_format,
            0,
            sense,
            offset,
            costs,
            col_lower,
            col_upper,
            row_lower,
            row_upper,
            a_start,
            a_index,
            a_value,
            None,
            None,
            None,
            integrality,
        )
    }

    /// HiGHS: HighsInterface.cpp::Highs::aFormatOk (lines 2863-2873).
    /// `assert(ok_format)` -> `debug_assert!`, as in `q_format_ok`.
    pub(crate) fn a_format_ok(&self, num_nz: HighsInt, format: HighsInt) -> bool {
        if num_nz == 0 {
            return true;
        }
        let ok_format = format == MatrixFormat::kColwise as HighsInt
            || format == MatrixFormat::kRowwise as HighsInt;
        if !ok_format {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Non-empty Constraint matrix has illegal format = {}\n",
                format
            );
        }
        debug_assert!(ok_format);
        ok_format
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn rowwise_raw_lp_is_stored_colwise_equivalent() {
        // min x0 + x1, rows: x0 + 2 x1 >= 1, 3 x0 >= 0 (row-wise input)
        let mut highs = Highs::default();
        let status = highs.pass_model_raw(
            2,
            2,
            3,
            MatrixFormat::kRowwise as HighsInt,
            ObjSense::kMinimize as HighsInt,
            0.0,
            Some(&[1.0, 1.0]),
            Some(&[0.0, 0.0]),
            Some(&[kHighsInf, kHighsInf]),
            Some(&[1.0, 0.0]),
            Some(&[kHighsInf, kHighsInf]),
            Some(&[0, 2]),
            Some(&[0, 1, 0]),
            Some(&[1.0, 2.0, 3.0]),
            None,
        );
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(highs.model_lp_.num_col_, 2);
        assert_eq!(highs.model_lp_.a_matrix_.num_nz(), 3);
    }

    #[test]
    fn only_the_first_num_col_and_num_row_entries_are_read() {
        // C++ assigns [p, p + n): longer input arrays are truncated.
        let mut highs = Highs::default();
        let status = highs.pass_model_raw(
            1,
            1,
            0,
            MatrixFormat::kColwise as HighsInt,
            ObjSense::kMinimize as HighsInt,
            0.0,
            Some(&[1.0, 7.0]),
            Some(&[0.0, 7.0]),
            Some(&[1.0, 7.0]),
            Some(&[0.0, 7.0]),
            Some(&[1.0, 7.0]),
            None,
            None,
            None,
            None,
        );
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(highs.model_lp_.col_cost_, vec![1.0]);
        assert_eq!(highs.model_lp_.col_upper_, vec![1.0]);
        assert_eq!(highs.model_lp_.row_lower_, vec![0.0]);
    }

    #[test]
    fn empty_matrix_format_is_not_checked() {
        let highs = Highs::default();
        assert!(highs.a_format_ok(0, 99));
    }

    fn pass_illegal_matrix_format() -> HighsStatus {
        let mut highs = Highs::default();
        highs.pass_model_raw(
            1,
            1,
            1,
            99,
            ObjSense::kMinimize as HighsInt,
            0.0,
            Some(&[1.0]),
            Some(&[0.0]),
            Some(&[1.0]),
            Some(&[0.0]),
            Some(&[1.0]),
            Some(&[0, 1]),
            Some(&[0]),
            Some(&[1.0]),
            None,
        )
    }

    /// `aFormatOk` logs and then `assert(ok_format)`s (HighsInterface.cpp:2909), so the
    /// error return is only reachable with debug assertions off.
    #[test]
    #[cfg(not(debug_assertions))]
    fn illegal_matrix_format_is_an_error() {
        assert_eq!(pass_illegal_matrix_format(), HighsStatus::kError);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "ok_format")]
    fn illegal_matrix_format_trips_cpp_assert() {
        pass_illegal_matrix_format();
    }
}
