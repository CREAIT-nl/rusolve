//! HiGHS source: ipm/ipx/model.cc (lines 501-980)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::super::control::{Control, ostream_double, scientific, textline};
use super::super::sparse_matrix::{scale_column, transpose};
use super::super::utils::infnorm;
use super::*;
use crate::util::highs_cdouble::ldexp_f64;
use crate::util::highs_hash::frexp;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: ipx/model.cc::CheckVectors (static, lines 501-522)
///
/// C++ returns `int`.
pub(super) fn check_vectors(
    m: Int,
    n: Int,
    rhs: &[f64],
    constr_type: &[u8],
    obj: &[f64],
    lb: &[f64],
    ub: &[f64],
) -> i32 {
    for i in 0..m as usize {
        if !rhs[i].is_finite() {
            return -1;
        }
    }
    for j in 0..n as usize {
        if !obj[j].is_finite() {
            return -2;
        }
    }
    for j in 0..n as usize {
        if !lb[j].is_finite() && lb[j] != f64::NEG_INFINITY {
            return -3;
        }
        if !ub[j].is_finite() && ub[j] != f64::INFINITY {
            return -3;
        }
        if lb[j] > ub[j] {
            return -3;
        }
    }
    for i in 0..m as usize {
        if constr_type[i] != b'=' && constr_type[i] != b'<' && constr_type[i] != b'>' {
            return -4;
        }
    }
    0
}

/// HiGHS: ipx/model.cc::CheckMatrix (static, lines 526-548)
pub(super) fn check_matrix(m: Int, n: Int, Ap: &[Int], Ai: &[Int], Ax: &[f64]) -> Int {
    if Ap[0] != 0 {
        return -5;
    }
    for j in 0..n as usize {
        if Ap[j] > Ap[j + 1] {
            return -5;
        }
    }
    for p in 0..Ap[n as usize] as usize {
        if !Ax[p].is_finite() {
            return -6;
        }
    }
    // Test for out of bound indices and duplicates.
    let mut marked: Vec<Int> = vec![-1; m as usize];
    for j in 0..n as usize {
        for p in Ap[j] as usize..Ap[j + 1] as usize {
            let i = Ai[p];
            if i < 0 || i >= m {
                return -7;
            }
            if marked[i as usize] == j as Int {
                return -8;
            }
            marked[i as usize] = j as Int;
        }
    }
    0
}

/// HiGHS: ipx/model.cc::EquilibrationFactor (static, lines 737-743)
///
/// C++ `int` arguments.
pub(super) fn equilibration_factor(expmin: i32, expmax: i32, exp: i32) -> f64 {
    if exp < expmin {
        return ldexp_f64(1.0, (expmin - exp + 1) / 2);
    }
    if exp > expmax {
        return ldexp_f64(1.0, -((exp - expmax + 1) / 2));
    }
    1.0
}

/// `std::minmax_element` values for a nonempty range: the first minimum and
/// the last maximum. Only the values are used by IPX, so ties do not matter.
fn minmax_values(v: &[f64]) -> (f64, f64) {
    let mut mn = v[0];
    let mut mx = v[0];
    for &x in &v[1..] {
        if x < mn {
            mn = x;
        }
        if !(x < mx) {
            mx = x;
        }
    }
    (mn, mx)
}

impl Model {
    /// HiGHS: ipx/model.cc::Model::CopyInput (lines 550-595)
    ///
    /// Returns `IPX_ERROR_argument_null` if any argument is `None`.
    pub(super) fn copy_input(
        &mut self,
        num_constr: Int,
        num_var: Int,
        Ap: Option<&[Int]>,
        Ai: Option<&[Int]>,
        Ax: Option<&[f64]>,
        rhs: Option<&[f64]>,
        constr_type: Option<&[u8]>,
        offset: f64,
        obj: Option<&[f64]>,
        lbuser: Option<&[f64]>,
        ubuser: Option<&[f64]>,
    ) -> Int {
        let (
            Some(ap),
            Some(ai),
            Some(ax),
            Some(rhs),
            Some(constr_type),
            Some(obj),
            Some(lbuser),
            Some(ubuser),
        ) = (Ap, Ai, Ax, rhs, constr_type, obj, lbuser, ubuser)
        else {
            return IPX_ERROR_argument_null;
        };
        if num_constr < 0 || num_var <= 0 {
            return IPX_ERROR_invalid_dimension;
        }
        if check_vectors(num_constr, num_var, rhs, constr_type, obj, lbuser, ubuser) != 0 {
            return IPX_ERROR_invalid_vector;
        }
        if check_matrix(num_constr, num_var, ap, ai, ax) != 0 {
            return IPX_ERROR_invalid_matrix;
        }
        let nc = num_constr as usize;
        let nv = num_var as usize;
        self.num_constr_ = num_constr;
        self.num_eqconstr_ = constr_type[..nc].iter().filter(|&&ch| ch == b'=').count() as Int;
        self.num_var_ = num_var;
        self.num_entries_ = ap[nv];
        self.num_free_var_ = 0;
        self.boxed_vars_.clear();
        for j in 0..nv {
            if lbuser[j].is_infinite() && ubuser[j].is_infinite() {
                self.num_free_var_ += 1;
            }
            if lbuser[j].is_finite() && ubuser[j].is_finite() {
                self.boxed_vars_.push(j as Int);
            }
        }
        self.constr_type_ = constr_type[..nc].to_vec();
        self.offset_ = offset;
        self.scaled_obj_ = Vector::from_slice(&obj[..nv]);
        self.scaled_rhs_ = Vector::from_slice(&rhs[..nc]);
        self.scaled_lbuser_ = Vector::from_slice(&lbuser[..nv]);
        self.scaled_ubuser_ = Vector::from_slice(&ubuser[..nv]);
        self.A_
            .load_from_arrays(num_constr, num_var, &ap[..nv], &ap[1..nv + 1], ai, ax);
        self.norm_obj_ = infnorm(&self.scaled_obj_);
        self.norm_rhs_ = infnorm(&self.scaled_rhs_);
        for &x in self.scaled_lbuser_.iter() {
            if x.is_finite() {
                self.norm_rhs_ = cpp_max(self.norm_rhs_, x.abs());
            }
        }
        for &x in self.scaled_ubuser_.iter() {
            if x.is_finite() {
                self.norm_rhs_ = cpp_max(self.norm_rhs_, x.abs());
            }
        }
        0
    }

    /// HiGHS: ipx/model.cc::Model::ScaleModel (lines 597-626)
    pub(super) fn scale_model(&mut self, control: &Control) {
        self.flipped_vars_.clear();
        for j in 0..self.num_var_ as usize {
            if self.scaled_ubuser_[j].is_finite() && self.scaled_lbuser_[j].is_infinite() {
                self.scaled_lbuser_[j] = -self.scaled_ubuser_[j];
                self.scaled_ubuser_[j] = f64::INFINITY;
                scale_column(&mut self.A_, j as Int, -1.0);
                self.scaled_obj_[j] *= -1.0;
                self.flipped_vars_.push(j as Int);
            }
        }
        self.colscale_.resize(0);
        self.rowscale_.resize(0);

        // Choose scaling method.
        if control.scale() >= 1 {
            self.equilibrate_matrix();
        }

        // Apply scaling to vectors.
        if self.colscale_.size() > 0 {
            debug_assert_eq!(self.colscale_.size(), self.num_var_);
            mul_assign(&mut self.scaled_obj_, &self.colscale_);
            div_assign(&mut self.scaled_lbuser_, &self.colscale_);
            div_assign(&mut self.scaled_ubuser_, &self.colscale_);
        }
        if self.rowscale_.size() > 0 {
            debug_assert_eq!(self.rowscale_.size(), self.num_constr_);
            mul_assign(&mut self.scaled_rhs_, &self.rowscale_);
        }
    }

    /// HiGHS: ipx/model.cc::Model::LoadPrimal (lines 628-666)
    pub(super) fn load_primal(&mut self) {
        self.num_rows_ = self.num_constr_;
        self.num_cols_ = self.num_var_;
        self.dualized_ = false;

        // Copy A and append identity matrix.
        self.AI_ = self.A_.clone();
        for i in 0..self.num_constr_ {
            self.AI_.push_back(i, 1.0);
            self.AI_.add_column();
        }
        debug_assert_eq!(self.AI_.cols(), self.num_var_ + self.num_constr_);

        // Copy vectors and set bounds on slack variables.
        self.b_.assign(&self.scaled_rhs_);
        self.c_.resize(self.num_var_ + self.num_constr_);
        self.c_.fill(0.0);
        copy_n(&self.scaled_obj_[..], self.num_var_, &mut self.c_[..]);
        self.lb_.resize(self.num_rows_ + self.num_cols_);
        copy_n(&self.scaled_lbuser_[..], self.num_var_, &mut self.lb_[..]);
        self.ub_.resize(self.num_rows_ + self.num_cols_);
        copy_n(&self.scaled_ubuser_[..], self.num_var_, &mut self.ub_[..]);
        let nv = self.num_var_ as usize;
        for i in 0..self.num_constr_ as usize {
            match self.constr_type_[i] {
                b'=' => {
                    self.lb_[nv + i] = 0.0;
                    self.ub_[nv + i] = 0.0;
                }
                b'<' => {
                    self.lb_[nv + i] = 0.0;
                    self.ub_[nv + i] = f64::INFINITY;
                }
                b'>' => {
                    self.lb_[nv + i] = f64::NEG_INFINITY;
                    self.ub_[nv + i] = 0.0;
                }
                _ => {}
            }
        }
    }

    /// HiGHS: ipx/model.cc::Model::LoadDual (lines 668-733)
    pub(super) fn load_dual(&mut self) {
        self.num_rows_ = self.num_var_;
        self.num_cols_ = self.num_constr_ + self.boxed_vars_.len() as Int;
        self.dualized_ = true;

        // Check that every variable with finite scaled_ubuser_ has finite
        // scaled_lbuser_ (must be the case after scaling).
        for j in 0..self.num_var_ as usize {
            if self.scaled_ubuser_[j].is_finite() {
                debug_assert!(self.scaled_lbuser_[j].is_finite());
            }
        }

        // Build AI.
        self.AI_ = transpose(&self.A_);
        for j in 0..self.num_var_ {
            if self.scaled_ubuser_[j as usize].is_finite() {
                self.AI_.push_back(j, -1.0);
                self.AI_.add_column();
            }
        }
        debug_assert_eq!(self.AI_.cols(), self.num_cols_);
        for i in 0..self.num_rows_ {
            self.AI_.push_back(i, 1.0);
            self.AI_.add_column();
        }

        // Build vectors.
        self.b_.assign(&self.scaled_obj_);
        self.c_.resize(self.num_cols_ + self.num_rows_);
        let mut put: usize = 0;
        for &x in self.scaled_rhs_.iter() {
            self.c_[put] = -x;
            put += 1;
        }
        for &x in self.scaled_ubuser_.iter() {
            if x.is_finite() {
                self.c_[put] = x;
                put += 1;
            }
        }
        debug_assert_eq!(put as Int, self.num_cols_);
        for &x in self.scaled_lbuser_.iter() {
            // If x is negative infinity, then the variable will be fixed and we can
            // give it any (finite) cost.
            self.c_[put] = if x.is_finite() { -x } else { 0.0 };
            put += 1;
        }
        self.lb_.resize(self.num_cols_ + self.num_rows_);
        self.ub_.resize(self.num_cols_ + self.num_rows_);
        for i in 0..self.num_constr_ as usize {
            match self.constr_type_[i] {
                b'=' => {
                    self.lb_[i] = f64::NEG_INFINITY;
                    self.ub_[i] = f64::INFINITY;
                }
                b'<' => {
                    self.lb_[i] = f64::NEG_INFINITY;
                    self.ub_[i] = 0.0;
                }
                b'>' => {
                    self.lb_[i] = 0.0;
                    self.ub_[i] = f64::INFINITY;
                }
                _ => {}
            }
        }
        for j in self.num_constr_ as usize..self.num_cols_ as usize {
            self.lb_[j] = 0.0;
            self.ub_[j] = f64::INFINITY;
        }
        for j in 0..self.num_var_ as usize {
            self.lb_[(self.num_cols_ as usize) + j] = 0.0;
            self.ub_[(self.num_cols_ as usize) + j] = if self.scaled_lbuser_[j].is_finite() {
                f64::INFINITY
            } else {
                0.0
            };
        }
    }

    /// HiGHS: ipx/model.cc::Model::EquilibrateMatrix (lines 745-836)
    pub(super) fn equilibrate_matrix(&mut self) {
        let m = self.A_.rows();
        let n = self.A_.cols();
        let (ap, ai, ax) = self.A_.arrays_mut();

        self.colscale_.resize(0);
        self.rowscale_.resize(0);

        // The absolute value of each nonzero entry of AI can be written as x*2^exp
        // with x in the range [0.5,1). We consider AI well scaled if each entry is
        // such that expmin <= exp <= expmax for parameters expmin and expmax.
        // For example,
        //
        //  expmin  expmax  allowed range
        //  -----------------------------
        //       0       2    [0.5,    4)
        //      -1       3    [0.25,   8)
        //       1       7    [1.0,  128)
        //
        // If AI is not well scaled, a recursive row and column equilibration is
        // applied. In each iteration, the scaling factors are roughly 1/sqrt(max),
        // where max is the maximum row or column entry. However, the factors are
        // truncated to powers of 2, so that no round-off errors occur.

        let expmin: i32 = 0;
        let expmax: i32 = 3;
        let maxround: Int = 10;

        // Quick return if entries are within the target range.
        let mut out_of_range = false;
        for p in 0..ap[n as usize] as usize {
            let (_, exp) = frexp(ax[p].abs());
            if exp < expmin || exp > expmax {
                out_of_range = true;
                break;
            }
        }
        if !out_of_range {
            return;
        }

        self.colscale_.resize(n);
        self.rowscale_.resize(m);
        self.colscale_.fill(1.0);
        self.rowscale_.fill(1.0);
        let mut colmax = Vector::new(n);
        let mut rowmax = Vector::new(m);

        for _round in 0..maxround {
            // Compute infinity norm of each row and column.
            rowmax.fill(0.0);
            for j in 0..n as usize {
                colmax[j] = 0.0;
                for p in ap[j] as usize..ap[j + 1] as usize {
                    let i = ai[p] as usize;
                    let xa = ax[p].abs();
                    colmax[j] = cpp_max(colmax[j], xa);
                    rowmax[i] = cpp_max(rowmax[i], xa);
                }
            }
            // Replace rowmax and colmax entries by scaling factors from this round.
            let mut out_of_range = false;
            for i in 0..m as usize {
                let (_, exp) = frexp(rowmax[i]);
                rowmax[i] = equilibration_factor(expmin, expmax, exp);
                if rowmax[i] != 1.0 {
                    out_of_range = true;
                    self.rowscale_[i] *= rowmax[i];
                }
            }
            for j in 0..n as usize {
                let (_, exp) = frexp(colmax[j]);
                colmax[j] = equilibration_factor(expmin, expmax, exp);
                if colmax[j] != 1.0 {
                    out_of_range = true;
                    self.colscale_[j] *= colmax[j];
                }
            }
            if !out_of_range {
                break;
            }
            // Rescale A.
            for j in 0..n as usize {
                for p in ap[j] as usize..ap[j + 1] as usize {
                    ax[p] *= colmax[j]; // column scaling
                    ax[p] *= rowmax[ai[p] as usize]; // row scaling
                }
            }
        }
    }

    /// HiGHS: ipx/model.cc::Model::FindDenseColumns (lines 838-861)
    pub(super) fn find_dense_columns(&mut self) {
        self.num_dense_cols_ = 0;
        self.nz_dense_ = self.rows() + 1;

        let mut colcount: Vec<Int> = vec![0; self.num_cols_ as usize];
        for j in 0..self.num_cols_ {
            colcount[j as usize] = self.AI_.end(j) - self.AI_.begin(j);
        }
        // std::less<Int>: pdqsort (non-branchless variant) on plain integers.
        pdqsort_by(&mut colcount[..], |a, b| a < b);

        for j in 1..self.num_cols_ as usize {
            if colcount[j] > cpp_max(40, 10 * colcount[j - 1]) {
                // j is the first dense column
                self.num_dense_cols_ = self.num_cols_ - j as Int;
                self.nz_dense_ = colcount[j];
                break;
            }
        }

        if self.num_dense_cols_ > 1000 {
            self.num_dense_cols_ = 0;
            self.nz_dense_ = self.rows() + 1;
        }
    }

    /// HiGHS: ipx/model.cc::Model::PrintCoefficientRange (lines 863-938)
    pub(super) fn print_coefficient_range(&self, control: &Control) {
        let mut amin = f64::INFINITY;
        let mut amax = 0.0;
        for j in 0..self.A_.cols() {
            for p in self.A_.begin(j)..self.A_.end(j) {
                let x = self.A_.value(p);
                if x != 0.0 {
                    amin = cpp_min(amin, x.abs());
                    amax = cpp_max(amax, x.abs());
                }
            }
        }
        if amin == f64::INFINITY {
            // no nonzero entries in A_
            amin = 0.0;
        }
        let mut logging = String::new();
        logging.push_str(&format!(
            "{}[{}, {}]\n",
            textline("Matrix range:"),
            scientific(amin, 5, 0),
            scientific(amax, 5, 0)
        ));
        control.h_log(&mut logging);

        let mut rhsmin = f64::INFINITY;
        let mut rhsmax = 0.0;
        for &x in self.scaled_rhs_.iter() {
            if x != 0.0 {
                rhsmin = cpp_min(rhsmin, x.abs());
                rhsmax = cpp_max(rhsmax, x.abs());
            }
        }
        if rhsmin == f64::INFINITY {
            // no nonzero entries in rhs
            rhsmin = 0.0;
        }
        logging.push_str(&format!(
            "{}[{}, {}]\n",
            textline("RHS range:"),
            scientific(rhsmin, 5, 0),
            scientific(rhsmax, 5, 0)
        ));
        control.h_log(&mut logging);

        let mut objmin = f64::INFINITY;
        let mut objmax = 0.0;
        for &x in self.scaled_obj_.iter() {
            if x != 0.0 {
                objmin = cpp_min(objmin, x.abs());
                objmax = cpp_max(objmax, x.abs());
            }
        }
        if objmin == f64::INFINITY {
            // no nonzero entries in obj
            objmin = 0.0;
        }
        logging.push_str(&format!(
            "{}[{}, {}]\n",
            textline("Objective range:"),
            scientific(objmin, 5, 0),
            scientific(objmax, 5, 0)
        ));
        control.h_log(&mut logging);

        let mut boundmin = f64::INFINITY;
        let mut boundmax = 0.0;
        for &x in self.scaled_lbuser_.iter() {
            if x != 0.0 && x.is_finite() {
                boundmin = cpp_min(boundmin, x.abs());
                boundmax = cpp_max(boundmax, x.abs());
            }
        }
        for &x in self.scaled_ubuser_.iter() {
            if x != 0.0 && x.is_finite() {
                boundmin = cpp_min(boundmin, x.abs());
                boundmax = cpp_max(boundmax, x.abs());
            }
        }
        if boundmin == f64::INFINITY {
            // no finite nonzeros entries in bounds
            boundmin = 0.0;
        }
        logging.push_str(&format!(
            "{}[{}, {}]\n",
            textline("Bounds range:"),
            scientific(boundmin, 5, 0),
            scientific(boundmax, 5, 0)
        ));
        control.h_log(&mut logging);
    }

    /// HiGHS: ipx/model.cc::Model::PrintPreprocessingLog (lines 940-980)
    pub(super) fn print_preprocessing_log(&self, control: &Control) {
        // Find the minimum and maximum scaling factor.
        let mut minscale = f64::INFINITY;
        let mut maxscale = 0.0;
        if self.colscale_.size() > 0 {
            let (mn, mx) = minmax_values(&self.colscale_);
            minscale = cpp_min(minscale, mn);
            maxscale = cpp_max(maxscale, mx);
        }
        if self.rowscale_.size() > 0 {
            let (mn, mx) = minmax_values(&self.rowscale_);
            minscale = cpp_min(minscale, mn);
            maxscale = cpp_max(maxscale, mx);
        }
        if minscale == f64::INFINITY {
            minscale = 1.0;
        }
        if maxscale == 0.0 {
            maxscale = 1.0;
        }

        let mut logging = String::new();
        logging.push_str(&format!(
            "Preprocessing\n{}{}\n{}{}\n",
            textline("Dualized model:"),
            if self.dualized() { "yes" } else { "no" },
            textline("Number of dense columns:"),
            self.num_dense_cols()
        ));
        control.h_log(&mut logging);
        if control.scale() > 0 {
            logging.push_str(&format!(
                "{}[{}, {}]\n",
                textline("Range of scaling factors:"),
                scientific(minscale, 8, 2),
                scientific(maxscale, 8, 2)
            ));
            control.h_log(&mut logging);
        }
        logging.push_str(&format!(
            "{}{}\n{}{}\n",
            textline("Scaled cost norm:   "),
            ostream_double(self.norm_c_),
            textline("Scaled bounds norm: "),
            ostream_double(self.norm_bounds_)
        ));
        control.h_log(&mut logging);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pins the C integer division `(expmin-exp+1)/2` and `-((exp-expmax+1)/2)`
    /// (truncation toward zero) for exponents just outside [expmin, expmax].
    #[test]
    fn equilibration_factor_truncates_like_c() {
        assert_eq!(equilibration_factor(0, 3, -1), 2.0); // (0+1+1)/2 = 1
        assert_eq!(equilibration_factor(0, 3, 0), 1.0);
        assert_eq!(equilibration_factor(0, 3, 3), 1.0);
        assert_eq!(equilibration_factor(0, 3, 4), 0.5); // -((4-3+1)/2) = -1
        assert_eq!(equilibration_factor(0, 3, 5), 0.5); // -((5-3+1)/2) = -1
        assert_eq!(equilibration_factor(0, 3, 6), 0.25); // -((6-3+1)/2) = -2
    }

    /// A duplicated row index inside one column is rejected with -8.
    #[test]
    fn check_matrix_rejects_duplicate_row() {
        let ap = [0, 2];
        let ai = [0, 0];
        let ax = [1.0, 2.0];
        assert_eq!(check_matrix(2, 1, &ap, &ai, &ax), -8);
        let ai_ok = [0, 1];
        assert_eq!(check_matrix(2, 1, &ap, &ai_ok, &ax), 0);
    }
}
