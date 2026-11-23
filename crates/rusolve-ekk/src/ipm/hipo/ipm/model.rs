//! HiGHS source: ipm/hipo/ipm/Model.h, ipm/hipo/ipm/Model.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The reformulated problem (`min c'x + x'Qx/2, Ax = b, x - xl = lower,
//! x + xu = upper`). Fields are `pub` because the `Preprocess*` actions are
//! C++ friends that rewrite them. `lp_orig_` is an owned clone of the
//! original `HighsLp` (C++: `const HighsLp*`). `offset_`, `AS_nz_`,
//! `NE_nz_lb_`, `NE_nz_ub_` and the norms are uninitialised in C++ and set
//! before use; they start at zero. `constraints_*` are `char` vectors (`u8`).

use super::iterate::Iterate;
use super::pre_process::Preprocessor;
use super::pre_process::PreprocessorPoint;
use crate::ipm::hipo::auxiliary::logger::{Logger, sci, textline};
use crate::ipm::hipo::auxiliary::vector_operations::infNorm;
use crate::ipm::hipo::ipm::parameters::kFreeVarsCloseRatio;
use crate::ipm::hipo::ipm::status::kStatusBadModel;
use crate::ipm::hipo::ipm::status::kStatusError;
use crate::ipm::hipo::{Int, Int64};
use crate::ipm::ipx::lp_solver::LpSolver;
use crate::ipm::ipx_wrapper::{fill_in_ipx_data, vec_data};
use crate::lp_data::hconst::ObjSense;
use crate::lp_data::hconst::kHighsInf;
use crate::lp_data::highs_lp::HighsLp;
use crate::model::highs_hessian::HighsHessian;
use crate::model::highs_hessian_utils::normalise::complete_hessian;
use crate::util::highs_utils::{cpp_max, cpp_min};
use crate::util::sparse_matrix::HighsSparseMatrix;

/// HiGHS: Model.h `class Model`
#[derive(Default)]
pub struct Model {
    // data of original problem
    pub n_orig_: Int,
    pub m_orig_: Int,
    pub lp_orig_: Option<HighsLp>,
    pub offset_: f64,
    pub sense_: ObjSense,
    pub rhs_orig_: Vec<f64>,
    pub constraints_orig_: Vec<u8>,

    // data of reformulated problem
    pub n_: Int,
    pub m_: Int,
    pub c_: Vec<f64>,
    pub b_: Vec<f64>,
    pub lower_: Vec<f64>,
    pub upper_: Vec<f64>,
    pub A_: HighsSparseMatrix,
    pub Q_: HighsHessian,
    pub constraints_: Vec<u8>,

    pub AS_nz_: Int64,
    pub NE_nz_lb_: Int64,
    pub NE_nz_ub_: Int64,

    pub colscale_: Vec<f64>,
    pub rowscale_: Vec<f64>,

    pub norm_unscaled_rhs_: f64,
    pub norm_scaled_rhs_: f64,
    pub norm_unscaled_obj_: f64,
    pub norm_scaled_obj_: f64,

    pub ready_: bool,

    pub preprocessor_: Preprocessor,

    // norms of rows and cols of A
    pub one_norm_cols_: Vec<f64>,
    pub one_norm_rows_: Vec<f64>,
    pub inf_norm_cols_: Vec<f64>,
    pub inf_norm_rows_: Vec<f64>,

    pub free_variables_: Vec<Int>,

    // rusolve (HiGHS 1.16 backport, default off): inert flags, default false; set by
    // Solver::load before init (preprocessing runs inside init)
    pub h16_ipx_free_rows_: bool,
    pub h16_ipm_gap_offset_: bool,
    pub h16_hipo_empty_rows_: bool,
}

impl Model {
    /// HiGHS: Model.cpp::Model::init (lines 11-42)
    pub fn init(&mut self, lp: &HighsLp, Q: &HighsHessian) -> Int {
        fill_in_ipx_data(
            lp,
            &mut self.n_,
            &mut self.m_,
            &mut self.offset_,
            &mut self.c_,
            &mut self.lower_,
            &mut self.upper_,
            &mut self.A_.start_,
            &mut self.A_.index_,
            &mut self.A_.value_,
            &mut self.b_,
            &mut self.constraints_,
            self.h16_ipx_free_rows_,
        );
        self.rhs_orig_ = self.b_.clone();
        self.constraints_orig_ = self.constraints_.clone();
        self.Q_ = Q.clone();
        if self.qp() {
            complete_hessian(self.n_, &mut self.Q_);
        }
        self.sense_ = lp.sense_;

        if self.checkData() != 0 {
            return kStatusBadModel;
        }

        self.lp_orig_ = Some(lp.clone());
        self.n_orig_ = self.n_;
        self.m_orig_ = self.m_;
        self.A_.num_col_ = self.n_;
        self.A_.num_row_ = self.m_;

        self.preprocess();

        self.nzBounds();
        self.computeNorms();

        // double transpose to sort indices of each column
        self.A_.ensure_rowwise();
        self.A_.ensure_colwise();

        if self.checkData() != 0 {
            return kStatusBadModel;
        }

        self.ready_ = true;

        0
    }

    /// HiGHS: Model.cpp::Model::nzBounds (lines 44-69)
    pub fn nzBounds(&mut self) {
        // compute lower and upper bounds for the number of nonzeros in normal
        // equations.
        let mut mark = vec![false; self.m_ as usize];
        self.NE_nz_lb_ = self.A_.num_row_ as Int64;
        self.NE_nz_ub_ = self.A_.num_row_ as Int64;
        let mut col: Int = 0;
        while col < self.A_.num_col_ {
            let mut used: Int64 = 0;
            let mut unused: Int64 = 0;
            let mut el: Int = self.A_.start_[col as usize];
            while el < self.A_.start_[(col + 1) as usize] {
                let row: Int = self.A_.index_[el as usize];
                if mark[row as usize] {
                    used += 1;
                } else {
                    mark[row as usize] = true;
                    unused += 1;
                }
                el += 1;
            }

            self.NE_nz_ub_ += (used + unused) * (used + unused - 1) / 2;
            self.NE_nz_lb_ += unused * (unused - 1) / 2 + used * unused;
            col += 1;
        }
        self.NE_nz_ub_ = cpp_min(
            self.NE_nz_ub_,
            (self.A_.num_row_ as Int64) * (self.A_.num_row_ + 1) as Int64 / 2,
        );

        self.AS_nz_ = self.A_.num_nz() as Int64
            + self.A_.num_row_ as Int64
            + if self.qp() {
                self.Q_.num_nz() as Int64
            } else {
                self.A_.num_col_ as Int64
            };
    }

    /// HiGHS: Model.cpp::Model::checkData (lines 71-114)
    pub fn checkData(&self) -> Int {
        // Check if model provided by the user is ok.
        // Return kStatusBadModel if something is wrong.

        // Dimensions are valid
        if self.n_ <= 0 || self.m_ < 0 {
            return kStatusBadModel;
        }

        // Vectors are of correct size
        if self.c_.len() as Int != self.n_
            || self.b_.len() as Int != self.m_
            || self.lower_.len() as Int != self.n_
            || self.upper_.len() as Int != self.n_
            || self.constraints_.len() as Int != self.m_
            || self.A_.start_.len() as Int != self.n_ + 1
            || self.A_.index_.len() as Int != self.A_.start_[self.A_.start_.len() - 1]
            || self.A_.value_.len() as Int != self.A_.start_[self.A_.start_.len() - 1]
        {
            return kStatusBadModel;
        }

        // Hessian is ok, for QPs only
        if self.qp()
            && (self.Q_.dim_ != self.n_
                || self.Q_.format_ != crate::lp_data::hconst::HessianFormat::kTriangular)
        {
            return kStatusBadModel;
        }

        // Vectors are valid
        let mut i: Int = 0;
        while i < self.n_ {
            if !self.c_[i as usize].is_finite() {
                return kStatusBadModel;
            }
            i += 1;
        }
        let mut i: Int = 0;
        while i < self.m_ {
            if !self.b_[i as usize].is_finite() {
                return kStatusBadModel;
            }
            i += 1;
        }
        let mut i: Int = 0;
        while i < self.n_ {
            if !self.lower_[i as usize].is_finite() && self.lower_[i as usize] != -f64::INFINITY {
                return kStatusBadModel;
            }
            if !self.upper_[i as usize].is_finite() && self.upper_[i as usize] != f64::INFINITY {
                return kStatusBadModel;
            }
            if self.lower_[i as usize] > self.upper_[i as usize] {
                return kStatusBadModel;
            }
            i += 1;
        }
        let mut i: Int = 0;
        while i < self.m_ {
            if self.constraints_[i as usize] != b'<'
                && self.constraints_[i as usize] != b'='
                && self.constraints_[i as usize] != b'>'
            {
                return kStatusBadModel;
            }
            i += 1;
        }

        // Matrix is valid
        let mut i: Int = 0;
        while i < self.A_.start_[self.n_ as usize] {
            if !self.A_.value_[i as usize].is_finite() {
                return kStatusBadModel;
            }
            i += 1;
        }

        0
    }

    /// HiGHS: Model.cpp::Model::preprocess (lines 116-116)
    /// `preprocessor_.apply(*this)`: take `preprocessor_` out with
    /// `std::mem::take`, apply it to `self`, put it back.
    pub fn preprocess(&mut self) {
        let mut preprocessor = std::mem::take(&mut self.preprocessor_);
        preprocessor.apply(self);
        self.preprocessor_ = preprocessor;
    }

    /// HiGHS: Model.cpp::Model::postprocess (lines 118-124)
    pub fn postprocess(
        &self,
        x: &mut Vec<f64>,
        xl: &mut Vec<f64>,
        xu: &mut Vec<f64>,
        slack: &mut Vec<f64>,
        y: &mut Vec<f64>,
        zl: &mut Vec<f64>,
        zu: &mut Vec<f64>,
        it: &Iterate,
    ) {
        let mut point = PreprocessorPoint {
            x,
            xl,
            xu,
            slack,
            y,
            zl,
            zu,
        };
        self.preprocessor_.undo(&mut point, self, it);
    }

    /// HiGHS: Model.cpp::Model::computeNorms (lines 126-178)
    pub fn computeNorms(&mut self) {
        self.norm_scaled_obj_ = infNorm(&self.c_);

        self.norm_unscaled_obj_ = 0.0;
        let mut i: Int = 0;
        while i < self.n_ {
            let mut val = self.c_[i as usize].abs();
            if self.scaled() {
                val /= self.colScale(i);
            }
            self.norm_unscaled_obj_ = cpp_max(self.norm_unscaled_obj_, val);
            i += 1;
        }

        self.norm_scaled_rhs_ = infNorm(&self.b_);
        for d in self.lower_.iter() {
            if d.is_finite() {
                self.norm_scaled_rhs_ = cpp_max(self.norm_scaled_rhs_, d.abs());
            }
        }
        for d in self.upper_.iter() {
            if d.is_finite() {
                self.norm_scaled_rhs_ = cpp_max(self.norm_scaled_rhs_, d.abs());
            }
        }

        self.norm_unscaled_rhs_ = 0.0;
        let mut i: Int = 0;
        while i < self.m_ {
            let mut val = self.b_[i as usize].abs();
            if self.scaled() {
                val /= self.rowScale(i);
            }
            self.norm_unscaled_rhs_ = cpp_max(self.norm_unscaled_rhs_, val);
            i += 1;
        }
        let mut i: Int = 0;
        while i < self.n_ {
            if self.lower_[i as usize].is_finite() {
                let mut val = self.lower_[i as usize].abs();
                if self.scaled() {
                    val *= self.colScale(i);
                }
                self.norm_unscaled_rhs_ = cpp_max(self.norm_unscaled_rhs_, val);
            }
            if self.upper_[i as usize].is_finite() {
                let mut val = self.upper_[i as usize].abs();
                if self.scaled() {
                    val *= self.colScale(i);
                }
                self.norm_unscaled_rhs_ = cpp_max(self.norm_unscaled_rhs_, val);
            }
            i += 1;
        }

        // norms of rows and cols of A
        self.one_norm_cols_.resize(self.n_ as usize, 0.0);
        self.one_norm_rows_.resize(self.m_ as usize, 0.0);
        self.inf_norm_cols_.resize(self.n_ as usize, 0.0);
        self.inf_norm_rows_.resize(self.m_ as usize, 0.0);
        let mut col: Int = 0;
        while col < self.n_ {
            let mut el: Int = self.A_.start_[col as usize];
            while el < self.A_.start_[(col + 1) as usize] {
                let row: Int = self.A_.index_[el as usize];
                let val: f64 = self.A_.value_[el as usize];
                self.one_norm_cols_[col as usize] += val.abs();
                self.one_norm_rows_[row as usize] += val.abs();
                self.inf_norm_rows_[row as usize] =
                    cpp_max(self.inf_norm_rows_[row as usize], val.abs());
                self.inf_norm_cols_[col as usize] =
                    cpp_max(self.inf_norm_cols_[col as usize], val.abs());
                el += 1;
            }
            col += 1;
        }
    }

    /// HiGHS: Model.cpp::Model::print (lines 180-353)
    pub fn print(&self, logger: &Logger) {
        let mut log_stream = String::new();

        log_stream.push_str(&textline("Rows:"));
        log_stream.push_str(&sci(self.m_ as f64, 0, 1));
        log_stream.push('\n');
        log_stream.push_str(&textline("Cols:"));
        log_stream.push_str(&sci(self.n_ as f64, 0, 1));
        log_stream.push('\n');
        log_stream.push_str(&textline("Nnz A:"));
        log_stream.push_str(&sci(self.A_.num_nz() as f64, 0, 1));
        log_stream.push('\n');

        if self.qp() {
            log_stream.push_str(&textline("Nnz Q:"));
            log_stream.push_str(&sci(self.Q_.num_nz() as f64, 0, 1));
            if self.nonSeparableQp() {
                log_stream.push_str(", non-separable\n");
            } else {
                log_stream.push_str(", separable\n");
            }
        }

        // compute max and min entry of A in absolute value
        let mut Amin: f64 = kHighsInf;
        let mut Amax: f64 = 0.0;
        let mut norm_inf_cols = vec![0.0f64; self.n_ as usize];
        let mut norm_inf_rows = vec![0.0f64; self.m_ as usize];
        let mut norm_one_cols = vec![0.0f64; self.n_ as usize];
        let mut norm_one_rows = vec![0.0f64; self.m_ as usize];
        let mut col: Int = 0;
        while col < self.n_ {
            let mut el: Int = self.A_.start_[col as usize];
            while el < self.A_.start_[(col + 1) as usize] {
                let row: Int = self.A_.index_[el as usize];
                let value: f64 = self.A_.value_[el as usize];
                if value != 0.0 {
                    Amin = cpp_min(Amin, value.abs());
                    Amax = cpp_max(Amax, value.abs());
                    norm_inf_cols[col as usize] = cpp_max(norm_inf_cols[col as usize], value.abs());
                    norm_inf_rows[row as usize] = cpp_max(norm_inf_rows[row as usize], value.abs());
                    norm_one_cols[col as usize] += value.abs();
                    norm_one_rows[row as usize] += value.abs();
                }
                el += 1;
            }
            col += 1;
        }
        if Amin.is_infinite() {
            Amin = 0.0;
        }

        let mut norm_inf_col_min: f64 = kHighsInf;
        let mut norm_inf_col_max: f64 = 0.0;
        let mut norm_inf_row_min: f64 = kHighsInf;
        let mut norm_inf_row_max: f64 = 0.0;
        for d in norm_inf_cols.iter() {
            norm_inf_col_min = cpp_min(norm_inf_col_min, *d);
            norm_inf_col_max = cpp_max(norm_inf_col_max, *d);
        }
        for d in norm_inf_rows.iter() {
            norm_inf_row_min = cpp_min(norm_inf_row_min, *d);
            norm_inf_row_max = cpp_max(norm_inf_row_max, *d);
        }

        let mut norm_one_col_min: f64 = kHighsInf;
        let mut norm_one_col_max: f64 = 0.0;
        let mut norm_one_row_min: f64 = kHighsInf;
        let mut norm_one_row_max: f64 = 0.0;
        for d in norm_one_cols.iter() {
            norm_one_col_min = cpp_min(norm_one_col_min, *d);
            norm_one_col_max = cpp_max(norm_one_col_max, *d);
        }
        for d in norm_one_rows.iter() {
            norm_one_row_min = cpp_min(norm_one_row_min, *d);
            norm_one_row_max = cpp_max(norm_one_row_max, *d);
        }

        // compute max and min entry of c
        let mut cmin: f64 = kHighsInf;
        let mut cmax: f64 = 0.0;
        let mut i: Int = 0;
        while i < self.n_ {
            if self.c_[i as usize] != 0.0 {
                cmin = cpp_min(cmin, self.c_[i as usize].abs());
                cmax = cpp_max(cmax, self.c_[i as usize].abs());
            }
            i += 1;
        }
        if cmin.is_infinite() {
            cmin = 0.0;
        }

        // compute max and min entry of b
        let mut bmin: f64 = kHighsInf;
        let mut bmax: f64 = 0.0;
        let mut i: Int = 0;
        while i < self.m_ {
            if self.b_[i as usize] != 0.0 {
                bmin = cpp_min(bmin, self.b_[i as usize].abs());
                bmax = cpp_max(bmax, self.b_[i as usize].abs());
            }
            i += 1;
        }
        if bmin.is_infinite() {
            bmin = 0.0;
        }

        // compute max and min entry of Q in absolute value
        let mut Qmin: f64 = kHighsInf;
        let mut Qmax: f64 = 0.0;
        for val in self.Q_.value_.iter() {
            if *val != 0.0 {
                Qmin = cpp_min(Qmin, val.abs());
                Qmax = cpp_max(Qmax, val.abs());
            }
        }
        if Qmin.is_infinite() {
            Qmin = 0.0;
        }

        // compute max and min for bounds intervals
        let mut boundmin: f64 = kHighsInf;
        let mut boundmax: f64 = 0.0;
        let mut i: Int = 0;
        while i < self.n_ {
            if self.lower_[i as usize].is_finite() && self.upper_[i as usize].is_finite() {
                let diff: f64 = (self.upper_[i as usize] - self.lower_[i as usize]).abs();
                boundmin = cpp_min(boundmin, diff);
                boundmax = cpp_max(boundmax, diff);
            }
            i += 1;
        }
        if boundmin.is_infinite() {
            boundmin = 0.0;
        }

        // compute max and min scaling
        let mut scalemin: f64 = kHighsInf;
        let mut scalemax: f64 = 0.0;
        if self.scaled() {
            let mut i: Int = 0;
            while i < self.n_ {
                scalemin = cpp_min(scalemin, self.colScale(i));
                scalemax = cpp_max(scalemax, self.colScale(i));
                i += 1;
            }
            let mut i: Int = 0;
            while i < self.m_ {
                scalemin = cpp_min(scalemin, self.rowScale(i));
                scalemax = cpp_max(scalemax, self.rowScale(i));
                i += 1;
            }
        }
        if scalemin.is_infinite() {
            scalemin = 0.0;
        }

        // print ranges
        log_stream.push_str(&textline("Range of A:"));
        log_stream.push('[');
        log_stream.push_str(&sci(Amin, 5, 1));
        log_stream.push_str(", ");
        log_stream.push_str(&sci(Amax, 5, 1));
        log_stream.push_str("]\n");

        if logger.debug(1) {
            log_stream.push_str(&textline("Inf-norm rows:"));
            log_stream.push('[');
            log_stream.push_str(&sci(norm_inf_row_min, 5, 1));
            log_stream.push_str(", ");
            log_stream.push_str(&sci(norm_inf_row_max, 5, 1));
            log_stream.push_str("]\n");

            log_stream.push_str(&textline("Inf-norm cols:"));
            log_stream.push('[');
            log_stream.push_str(&sci(norm_inf_col_min, 5, 1));
            log_stream.push_str(", ");
            log_stream.push_str(&sci(norm_inf_col_max, 5, 1));
            log_stream.push_str("]\n");

            log_stream.push_str(&textline("One-norm rows:"));
            log_stream.push('[');
            log_stream.push_str(&sci(norm_one_row_min, 5, 1));
            log_stream.push_str(", ");
            log_stream.push_str(&sci(norm_one_row_max, 5, 1));
            log_stream.push_str("]\n");

            log_stream.push_str(&textline("One-norm cols:"));
            log_stream.push('[');
            log_stream.push_str(&sci(norm_one_col_min, 5, 1));
            log_stream.push_str(", ");
            log_stream.push_str(&sci(norm_one_col_max, 5, 1));
            log_stream.push_str("]\n");
        }

        log_stream.push_str(&textline("Range of b:"));
        log_stream.push('[');
        log_stream.push_str(&sci(bmin, 5, 1));
        log_stream.push_str(", ");
        log_stream.push_str(&sci(bmax, 5, 1));
        log_stream.push_str("]\n");

        log_stream.push_str(&textline("Range of c:"));
        log_stream.push('[');
        log_stream.push_str(&sci(cmin, 5, 1));
        log_stream.push_str(", ");
        log_stream.push_str(&sci(cmax, 5, 1));
        log_stream.push_str("]\n");

        if self.qp() {
            log_stream.push_str(&textline("Range of Q:"));
            log_stream.push('[');
            log_stream.push_str(&sci(Qmin, 5, 1));
            log_stream.push_str(", ");
            log_stream.push_str(&sci(Qmax, 5, 1));
            log_stream.push_str("]\n");
        }

        log_stream.push_str(&textline("Range of bound intervals:"));
        log_stream.push('[');
        log_stream.push_str(&sci(boundmin, 5, 1));
        log_stream.push_str(", ");
        log_stream.push_str(&sci(boundmax, 5, 1));
        log_stream.push_str("]\n");

        log_stream.push_str(&textline("Scaling coefficients:"));
        log_stream.push('[');
        log_stream.push_str(&sci(scalemin, 5, 1));
        log_stream.push_str(", ");
        log_stream.push_str(&sci(scalemax, 5, 1));
        log_stream.push_str("]\n");

        if logger.debug(1) {
            self.preprocessor_.print(&mut log_stream);

            log_stream.push_str("Expected nnz: ");
            log_stream.push_str("AS ");
            log_stream.push_str(&sci(self.AS_nz_ as f64, 0, 1));
            log_stream.push_str("; ");
            log_stream.push_str("NE [");
            log_stream.push_str(&sci(self.NE_nz_lb_ as f64, 0, 1));
            log_stream.push_str(", ");
            log_stream.push_str(&sci(self.NE_nz_ub_ as f64, 0, 1));
            log_stream.push_str("]\n");
        }

        logger.print(&log_stream);
    }

    /// HiGHS: Model.cpp::Model::loadIntoIpx (lines 355-374)
    pub fn loadIntoIpx(&self, lps: &mut LpSolver) -> Int {
        let mut ipx_m: Int = 0;
        let mut ipx_n: Int = 0;
        let mut ipx_b: Vec<f64> = Vec::new();
        let mut ipx_c: Vec<f64> = Vec::new();
        let mut ipx_lower: Vec<f64> = Vec::new();
        let mut ipx_upper: Vec<f64> = Vec::new();
        let mut ipx_A_vals: Vec<f64> = Vec::new();
        let mut ipx_A_ptr: Vec<Int> = Vec::new();
        let mut ipx_A_rows: Vec<Int> = Vec::new();
        let mut ipx_constraints: Vec<u8> = Vec::new();
        let mut ipx_offset: f64 = 0.0;

        let lp_orig = match &self.lp_orig_ {
            Some(lp) => lp,
            None => return kStatusError,
        };

        fill_in_ipx_data(
            lp_orig,
            &mut ipx_n,
            &mut ipx_m,
            &mut ipx_offset,
            &mut ipx_c,
            &mut ipx_lower,
            &mut ipx_upper,
            &mut ipx_A_ptr,
            &mut ipx_A_rows,
            &mut ipx_A_vals,
            &mut ipx_b,
            &mut ipx_constraints,
            self.h16_ipx_free_rows_,
        );

        // PORT NOTE: `std::vector::data()` of an empty vector is null (ipx_wrapper.rs
        // "Call surface"), hence `vec_data`.
        lps.load_model(
            ipx_n,
            ipx_offset,
            vec_data(&ipx_c),
            vec_data(&ipx_lower),
            vec_data(&ipx_upper),
            ipx_m,
            vec_data(&ipx_A_ptr),
            vec_data(&ipx_A_rows),
            vec_data(&ipx_A_vals),
            vec_data(&ipx_b),
            vec_data(&ipx_constraints),
        )
    }

    /// HiGHS: Model.cpp::Model::printDense (lines 376-418)
    pub fn printDense(&self) {
        let mut Adense = vec![vec![0.0f64; self.n_ as usize]; self.m_ as usize];
        let mut Qdense = vec![vec![0.0f64; self.n_ as usize]; self.n_ as usize];

        let mut col: Int = 0;
        while col < self.n_ {
            let mut el: Int = self.A_.start_[col as usize];
            while el < self.A_.start_[(col + 1) as usize] {
                let row: Int = self.A_.index_[el as usize];
                let val: f64 = self.A_.value_[el as usize];
                Adense[row as usize][col as usize] = val;
                el += 1;
            }
            col += 1;
        }
        let mut col: Int = 0;
        while col < self.n_ {
            let mut el: Int = self.Q_.start_[col as usize];
            while el < self.Q_.start_[(col + 1) as usize] {
                let row: Int = self.Q_.index_[el as usize];
                let val: f64 = self.Q_.value_[el as usize];
                Qdense[row as usize][col as usize] = val;
                el += 1;
            }
            col += 1;
        }

        // C++ lines 395-417: printf-based dense printing.
        // Skipped: IO is out of scope.
        let _ = &Adense;
        let _ = &Qdense;
    }

    /// HiGHS: Model.cpp::Model::adjustFreeVars (lines 420-442)
    pub fn adjustFreeVars(
        &mut self,
        x: &mut Vec<f64>,
        xl: &mut Vec<f64>,
        xu: &mut Vec<f64>,
        logger: &Logger,
    ) {
        for k in 0..self.free_variables_.len() {
            let i: Int = self.free_variables_[k];
            if x[i as usize] < self.lower_[i as usize] * kFreeVarsCloseRatio {
                // getting close to lower bound
                let new_lower: f64 = x[i as usize] / kFreeVarsCloseRatio;
                logger.printDetailed(&format!(
                    "Free var {:5} is at {} with lb {}, lb changed to {}\n",
                    i,
                    sci(x[i as usize], 8, 1),
                    sci(self.lower_[i as usize], 8, 1),
                    sci(new_lower, 8, 1)
                ));
                self.lower_[i as usize] = new_lower;
                xl[i as usize] = x[i as usize] - self.lower_[i as usize];
            }
            if x[i as usize] > self.upper_[i as usize] * kFreeVarsCloseRatio {
                // getting close to upper bound
                let new_upper: f64 = x[i as usize] / kFreeVarsCloseRatio;
                logger.printDetailed(&format!(
                    "Free var {:5} is at {} with ub {}, ub changed to {}\n",
                    i,
                    sci(x[i as usize], 8, 1),
                    sci(self.upper_[i as usize], 8, 1),
                    sci(new_upper, 8, 1)
                ));
                self.upper_[i as usize] = new_upper;
                xu[i as usize] = self.upper_[i as usize] - x[i as usize];
            }
        }
    }

    // Inline getters of Model.h (complete).
    pub fn normScaledRhs(&self) -> f64 {
        self.norm_scaled_rhs_
    }
    pub fn normScaledObj(&self) -> f64 {
        self.norm_scaled_obj_
    }
    pub fn normUnscaledObj(&self) -> f64 {
        self.norm_unscaled_obj_
    }
    pub fn normUnscaledRhs(&self) -> f64 {
        self.norm_unscaled_rhs_
    }
    pub fn hasLb(&self, j: Int) -> bool {
        self.lower_[j as usize].is_finite()
    }
    pub fn hasUb(&self, j: Int) -> bool {
        self.upper_[j as usize].is_finite()
    }
    pub fn m(&self) -> Int {
        self.m_
    }
    pub fn n(&self) -> Int {
        self.n_
    }
    pub fn n_orig(&self) -> Int {
        self.n_orig_
    }
    pub fn m_orig(&self) -> Int {
        self.m_orig_
    }
    pub fn lpOrig(&self) -> Option<&HighsLp> {
        self.lp_orig_.as_ref()
    }
    pub fn rhsOrig(&self) -> &[f64] {
        &self.rhs_orig_
    }
    pub fn constraintsOrig(&self) -> &[u8] {
        &self.constraints_orig_
    }
    pub fn qp(&self) -> bool {
        !self.Q_.empty()
    }
    pub fn nonSeparableQp(&self) -> bool {
        self.qp() && !self.Q_.is_diagonal()
    }
    pub fn sense(&self) -> f64 {
        (self.sense_ as Int) as f64
    }
    pub fn A(&self) -> &HighsSparseMatrix {
        &self.A_
    }
    pub fn Q(&self) -> &HighsHessian {
        &self.Q_
    }
    pub fn b(&self) -> &[f64] {
        &self.b_
    }
    pub fn c(&self) -> &[f64] {
        &self.c_
    }
    pub fn lb(&self, i: Int) -> f64 {
        self.lower_[i as usize]
    }
    pub fn ub(&self, i: Int) -> f64 {
        self.upper_[i as usize]
    }
    pub fn constraint(&self, i: Int) -> u8 {
        self.constraints_[i as usize]
    }
    pub fn colScale(&self, i: Int) -> f64 {
        if self.scaled() {
            self.colscale_[i as usize]
        } else {
            1.0
        }
    }
    pub fn rowScale(&self, i: Int) -> f64 {
        if self.scaled() {
            self.rowscale_[i as usize]
        } else {
            1.0
        }
    }
    pub fn ready(&self) -> bool {
        self.ready_
    }
    pub fn scaled(&self) -> bool {
        !self.colscale_.is_empty()
    }
    pub fn offset(&self) -> f64 {
        self.offset_
    }
    pub fn oneNormRows(&self, i: Int) -> f64 {
        self.one_norm_rows_[i as usize]
    }
    pub fn oneNormCols(&self, i: Int) -> f64 {
        self.one_norm_cols_[i as usize]
    }
    pub fn infNormRows(&self, i: Int) -> f64 {
        self.inf_norm_rows_[i as usize]
    }
    pub fn infNormCols(&self, i: Int) -> f64 {
        self.inf_norm_cols_[i as usize]
    }
    pub fn nzAS(&self) -> Int64 {
        self.AS_nz_
    }
    pub fn nzNElb(&self) -> Int64 {
        self.NE_nz_lb_
    }
    pub fn nzNEub(&self) -> Int64 {
        self.NE_nz_ub_
    }
    pub fn freeVars(&self) -> &[Int] {
        &self.free_variables_
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn default_model_is_unscaled_lp() {
        let m = Model::default();
        assert!(!m.ready() && !m.scaled() && !m.qp());
        assert_eq!(m.sense(), 1.0);
    }

    use crate::io::highs_io::{HighsLogOptions, kHighsLogDevLevelDetailed, kHighsLogDevLevelInfo};
    use crate::lp_data::highs_callback::{HighsCallbackInput, HighsCallbackOutput};
    use std::sync::{Arc, Mutex};

    /// Row 0 general-bounded (slack column), row 1 free (dropped), row 2
    /// equality: the IPX form has n = 3, m = 2 (ipx_wrapper.rs `tiny_lp`).
    fn tiny_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 3,
            col_cost_: vec![3.0, -4.0],
            col_lower_: vec![0.0, -kHighsInf],
            col_upper_: vec![10.0, kHighsInf],
            row_lower_: vec![1.0, -kHighsInf, 2.0],
            row_upper_: vec![4.0, kHighsInf, 2.0],
            sense_: ObjSense::kMaximize,
            offset_: 0.5,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 1, 2];
        lp.a_matrix_.value_ = vec![2.0, 3.0, 5.0, 7.0];
        lp
    }

    /// Unscaled LP model filled by hand (no `init`): m = 1, n = 2.
    fn hand_model() -> Model {
        let mut m = Model {
            n_: 2,
            m_: 1,
            c_: vec![-3.0, 0.0],
            b_: vec![1.5],
            lower_: vec![-1.0, -kHighsInf],
            upper_: vec![2.0, kHighsInf],
            constraints_: vec![b'='],
            ..Default::default()
        };
        m.A_.num_col_ = 2;
        m.A_.num_row_ = 1;
        m.A_.start_ = vec![0, 1, 2];
        m.A_.index_ = vec![0, 0];
        m.A_.value_ = vec![-1.0, 3.0];
        m
    }

    fn capture_logger(dev_level: Int) -> (Logger, Arc<Mutex<String>>) {
        let out = Arc::new(Mutex::new(String::new()));
        let sink = out.clone();
        let cb = crate::lp_data::highs_callback::HighsCallbackFunctionType::new(
            move |_t: Int, msg: &str, _o: &HighsCallbackOutput, _i: &mut HighsCallbackInput| {
                sink.lock().unwrap().push_str(msg);
            },
        );
        let opts = HighsLogOptions {
            log_to_console: false,
            log_dev_level: dev_level,
            user_callback: Some(cb),
            user_callback_active: true,
            ..Default::default()
        };
        (
            Logger {
                log_options_: Some(opts),
            },
            out,
        )
    }

    #[test]
    fn init_tiny_lp() {
        let lp = tiny_lp();
        let mut m = Model::default();
        assert_eq!(m.init(&lp, &HighsHessian::default()), 0);
        assert!(m.ready());
        assert!(!m.qp());
        assert_eq!(m.n_orig(), 3);
        assert_eq!(m.m_orig(), 2);
        assert_eq!(m.rhsOrig(), &[0.0, 2.0]);
        assert_eq!(m.constraintsOrig(), b"==");
        assert_eq!(m.sense(), -1.0);
        assert_eq!(m.offset(), -0.5);
        assert!(m.lpOrig().is_some());
        assert_eq!(m.checkData(), 0);
        assert_eq!(m.A().num_col_, m.n());
        assert_eq!(m.A().num_row_, m.m());
        // Model.cpp 34-35: the double transpose leaves the indices of each column sorted.
        for col in 0..m.n() as usize {
            let (s, e) = (m.A_.start_[col] as usize, m.A_.start_[col + 1] as usize);
            assert!(m.A_.index_[s..e].windows(2).all(|w| w[0] < w[1]));
        }
        assert_eq!(m.one_norm_cols_.len(), m.n() as usize);
        assert_eq!(m.inf_norm_rows_.len(), m.m() as usize);
        assert!(m.nzNElb() <= m.nzNEub());
    }

    #[test]
    fn init_rejects_bad_models() {
        // lower > upper (Model.cpp 102): rejected by the first checkData,
        // before lp_orig_ is set.
        let mut lp = tiny_lp();
        lp.col_lower_[0] = 11.0;
        let mut m = Model::default();
        assert_eq!(m.init(&lp, &HighsHessian::default()), kStatusBadModel);
        assert!(!m.ready() && m.lpOrig().is_none() && m.n_orig() == 0);
        // Non-finite cost (Model.cpp 94).
        let mut lp = tiny_lp();
        lp.col_cost_[1] = f64::NAN;
        let mut m = Model::default();
        assert_eq!(m.init(&lp, &HighsHessian::default()), kStatusBadModel);
        // Non-finite matrix value (Model.cpp 111).
        let mut lp = tiny_lp();
        lp.a_matrix_.value_[3] = f64::INFINITY;
        let mut m = Model::default();
        assert_eq!(m.init(&lp, &HighsHessian::default()), kStatusBadModel);
    }

    #[test]
    fn check_data_branches() {
        assert_eq!(hand_model().checkData(), 0);
        // n_ <= 0 (Model.cpp 76).
        assert_eq!(Model::default().checkData(), kStatusBadModel);
        // Wrong size (Model.cpp 79-86).
        let mut m = hand_model();
        m.b_.push(0.0);
        assert_eq!(m.checkData(), kStatusBadModel);
        let mut m = hand_model();
        m.A_.value_.pop();
        assert_eq!(m.checkData(), kStatusBadModel);
        // NaN bounds are rejected, infinite ones of the right sign are not (98-101).
        let mut m = hand_model();
        m.lower_[1] = f64::NAN;
        assert_eq!(m.checkData(), kStatusBadModel);
        let mut m = hand_model();
        m.upper_[0] = -f64::INFINITY;
        assert_eq!(m.checkData(), kStatusBadModel);
        let mut m = hand_model();
        m.lower_[0] = f64::INFINITY;
        assert_eq!(m.checkData(), kStatusBadModel);
        // NaN rhs (96) and bad constraint type (105-107).
        let mut m = hand_model();
        m.b_[0] = f64::NAN;
        assert_eq!(m.checkData(), kStatusBadModel);
        let mut m = hand_model();
        m.constraints_[0] = b'x';
        assert_eq!(m.checkData(), kStatusBadModel);
        // QP: Hessian dimension must equal n_ (89-90).
        let mut m = hand_model();
        m.Q_.dim_ = 1;
        m.Q_.start_ = vec![0, 1];
        m.Q_.index_ = vec![0];
        m.Q_.value_ = vec![1.0];
        assert!(m.qp());
        assert_eq!(m.checkData(), kStatusBadModel);
        m.Q_.dim_ = 2;
        m.Q_.start_ = vec![0, 1, 1];
        assert_eq!(m.checkData(), 0);
    }

    #[test]
    fn nz_bounds() {
        // m = 3: col 0 = rows {0, 1}, col 1 = rows {1, 2}.
        let mut m = Model {
            m_: 3,
            ..Default::default()
        };
        m.A_.num_col_ = 2;
        m.A_.num_row_ = 3;
        m.A_.start_ = vec![0, 2, 4];
        m.A_.index_ = vec![0, 1, 1, 2];
        m.A_.value_ = vec![1.0; 4];
        m.nzBounds();
        // col 0: unused 2 -> ub += 1, lb += 1; col 1: used 1, unused 1 -> ub += 1, lb += 1.
        assert_eq!((m.nzNElb(), m.nzNEub()), (5, 5));
        assert_eq!(m.nzAS(), 4 + 3 + 2);

        // m = 2, both columns on rows {0, 1}: ub = 4 is capped at 2 * 3 / 2 = 3.
        let mut m = Model {
            m_: 2,
            ..Default::default()
        };
        m.A_.num_col_ = 2;
        m.A_.num_row_ = 2;
        m.A_.start_ = vec![0, 2, 4];
        m.A_.index_ = vec![0, 1, 0, 1];
        m.A_.value_ = vec![1.0; 4];
        m.nzBounds();
        assert_eq!((m.nzNElb(), m.nzNEub()), (3, 3));
        // QP: AS counts Q's nonzeros instead of n.
        m.Q_.dim_ = 2;
        m.Q_.start_ = vec![0, 1, 3];
        m.Q_.index_ = vec![0, 0, 1];
        m.Q_.value_ = vec![1.0, 2.0, 3.0];
        m.nzBounds();
        assert_eq!(m.nzAS(), 4 + 2 + 3);
    }

    #[test]
    fn compute_norms_unscaled_and_scaled() {
        let mut m = hand_model();
        m.computeNorms();
        assert_eq!(m.normScaledObj(), 3.0);
        assert_eq!(m.normUnscaledObj(), 3.0);
        assert_eq!(m.normScaledRhs(), 2.0);
        assert_eq!(m.normUnscaledRhs(), 2.0);
        assert_eq!((m.oneNormCols(0), m.oneNormCols(1)), (1.0, 3.0));
        assert_eq!((m.infNormCols(0), m.infNormCols(1)), (1.0, 3.0));
        assert_eq!((m.oneNormRows(0), m.infNormRows(0)), (4.0, 3.0));

        let mut m = hand_model();
        m.c_ = vec![-3.0, 8.0];
        m.colscale_ = vec![2.0, 4.0];
        m.rowscale_ = vec![0.5];
        assert!(m.scaled());
        m.computeNorms();
        assert_eq!(m.normScaledObj(), 8.0);
        // max(3 / 2, 8 / 4)
        assert_eq!(m.normUnscaledObj(), 2.0);
        // max(|b| = 1.5, |lower_0| = 1, |upper_0| = 2)
        assert_eq!(m.normScaledRhs(), 2.0);
        // max(1.5 / 0.5, 1 * 2, 2 * 2)
        assert_eq!(m.normUnscaledRhs(), 4.0);
    }

    #[test]
    fn adjust_free_vars() {
        let mut m = Model {
            n_: 3,
            lower_: vec![-10.0, -10.0, -10.0],
            upper_: vec![10.0, 10.0, 10.0],
            free_variables_: vec![0, 2],
            ..Default::default()
        };
        // x_0 near lb, x_1 near lb but not free, x_2 near ub.
        let mut x = vec![-6.0, -6.0, 7.0];
        let mut xl = vec![4.0, 4.0, 17.0];
        let mut xu = vec![16.0, 16.0, 3.0];
        let (logger, out) = capture_logger(kHighsLogDevLevelDetailed);
        m.adjustFreeVars(&mut x, &mut xl, &mut xu, &logger);
        assert_eq!(m.lower_, vec![-6.0 / kFreeVarsCloseRatio, -10.0, -10.0]);
        assert_eq!(m.upper_, vec![10.0, 10.0, 7.0 / kFreeVarsCloseRatio]);
        assert_eq!(xl, vec![-6.0 - (-6.0 / kFreeVarsCloseRatio), 4.0, 17.0]);
        assert_eq!(xu, vec![16.0, 16.0, 7.0 / kFreeVarsCloseRatio - 7.0]);
        assert_eq!(x, vec![-6.0, -6.0, 7.0]);
        assert_eq!(
            out.lock().unwrap().as_str(),
            "Free var     0 is at -6.0e+00 with lb -1.0e+01, lb changed to -1.2e+01\n\
             Free var     2 is at  7.0e+00 with ub  1.0e+01, ub changed to  1.4e+01\n"
        );
    }

    #[test]
    fn load_into_ipx_without_lp_is_error() {
        let m = Model::default();
        let mut lps = LpSolver::default();
        assert_eq!(m.loadIntoIpx(&mut lps), kStatusError);
    }

    #[test]
    fn print_lp_lines() {
        let mut m = hand_model();
        m.nzBounds();
        let (logger, out) = capture_logger(0);
        m.print(&logger);
        let expected = format!(
            "{:<32}1.0e+00\n{:<32}2.0e+00\n{:<32}2.0e+00\n\
             {:<32}[1.0e+00, 3.0e+00]\n{:<32}[1.5e+00, 1.5e+00]\n\
             {:<32}[3.0e+00, 3.0e+00]\n{:<32}[3.0e+00, 3.0e+00]\n\
             {:<32}[0.0e+00, 0.0e+00]\n",
            "Rows:",
            "Cols:",
            "Nnz A:",
            "Range of A:",
            "Range of b:",
            "Range of c:",
            "Range of bound intervals:",
            "Scaling coefficients:"
        );
        assert_eq!(out.lock().unwrap().as_str(), expected);

        // debug(1): norms, Preprocessor::print (empty stack) and the nnz estimates.
        let (logger, out) = capture_logger(kHighsLogDevLevelInfo);
        m.print(&logger);
        let s = out.lock().unwrap();
        assert!(s.contains(&format!("{:<32}[3.0e+00, 3.0e+00]\n", "Inf-norm rows:")));
        assert!(s.contains(&format!("{:<32}[1.0e+00, 3.0e+00]\n", "One-norm cols:")));
        assert!(s.ends_with("Expected nnz: AS 5.0e+00; NE [1.0e+00, 1.0e+00]\n"));
    }
}
