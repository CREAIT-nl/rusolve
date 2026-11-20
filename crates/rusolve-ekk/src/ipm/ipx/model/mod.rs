//! HiGHS source: ipm/ipx/model.h, model.cc (class Model)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `model.cc` (1528 lines) is split: [`part1`] lines 10-497 (load, info,
//! pre/postsolve API), [`part2`] lines 501-980 (input checks, scaling,
//! primal/dual computational form, logging), [`part3`] lines 982-1526
//! (scale/dualize helpers and the free residual functions). Header inline
//! members are implemented here.
//!
//! # Null pointers
//!
//! Every C pointer argument that the C++ checks against `nullptr` is an
//! `Option`. `CopyInput` returns `IPX_ERROR_argument_null` if any input
//! pointer is null; the IpxWrapper passes `std::vector::data()`, which is
//! `nullptr` for an empty libc++ vector (for example an LP with no rows).
//! The wrapper must pass `None` exactly where libc++ `data()` would be null.

use super::sparse_matrix::SparseMatrix;
use super::*;

mod part1;
mod part2;
mod part3;

pub use part3::{dual_infeasibility, dual_residual, primal_infeasibility, primal_residual};

/// HiGHS: ipx/model.h::Model (lines 31-393)
///
/// PORT NOTE: `double offset_;` has no initializer in C++; `Model()` is
/// `= default` and `clear()` does not reset it, so it is only read after
/// `CopyInput` set it. Rust starts it at `0.0`.
#[derive(Debug, Clone, Default)]
pub struct Model {
    // Computational form model.
    /// model was dualized in preprocessing?
    dualized_: bool,
    /// # rows of AI
    num_rows_: Int,
    /// # structural columns of AI
    num_cols_: Int,
    /// # columns classified as dense
    num_dense_cols_: Int,
    /// minimum # nonzeros in a dense column
    nz_dense_: Int,
    /// matrix AI columnwise
    AI_: SparseMatrix,
    /// matrix AI rowwise
    AIt_: SparseMatrix,
    b_: Vector,
    c_: Vector,
    lb_: Vector,
    ub_: Vector,
    /// infinity norm of [b;lb;ub]
    norm_bounds_: f64,
    /// infinity norm of c
    norm_c_: f64,

    // User model after scaling.
    /// # constraints
    num_constr_: Int,
    /// # equality constraints
    num_eqconstr_: Int,
    /// # variables
    num_var_: Int,
    /// # free variables
    num_free_var_: Int,
    /// # entries in input matrix
    num_entries_: Int,
    /// indices of boxed variables
    boxed_vars_: Vec<Int>,
    /// C++ `std::vector<char>`: `'='`, `'<'` or `'>'`
    constr_type_: Vec<u8>,
    /// Infnorm(obj) as given by user
    norm_obj_: f64,
    /// Infnorm(rhs,lb,ub) as given by user
    norm_rhs_: f64,
    offset_: f64,
    scaled_obj_: Vector,
    scaled_rhs_: Vector,
    scaled_lbuser_: Vector,
    scaled_ubuser_: Vector,
    /// is cleared after preprocessing
    A_: SparseMatrix,

    // Data from ScaleModel() that is required by ScaleBack*().
    flipped_vars_: Vec<Int>,
    colscale_: Vector,
    rowscale_: Vector,
}

impl Model {
    /// HiGHS: ipx/model.h::Model::Model (line 40, `= default`)
    pub fn new() -> Self {
        Self::default()
    }
    /// HiGHS: ipx/model.h::Model::empty (line 69)
    pub fn empty(&self) -> bool {
        self.cols() == 0
    }
    /// HiGHS: ipx/model.h::Model::rows (line 75)
    pub fn rows(&self) -> Int {
        self.num_rows_
    }
    /// HiGHS: ipx/model.h::Model::cols (line 79)
    pub fn cols(&self) -> Int {
        self.num_cols_
    }
    /// HiGHS: ipx/model.h::Model::num_dense_cols (line 82)
    pub fn num_dense_cols(&self) -> Int {
        self.num_dense_cols_
    }
    /// HiGHS: ipx/model.h::Model::IsDenseColumn (lines 85-87)
    pub fn is_dense_column(&self, j: Int) -> bool {
        self.AI_.entries_col(j) >= self.nz_dense_
    }
    /// HiGHS: ipx/model.h::Model::dualized (line 90)
    pub fn dualized(&self) -> bool {
        self.dualized_
    }
    /// HiGHS: ipx/model.h::Model::AI (line 93)
    pub fn AI(&self) -> &SparseMatrix {
        &self.AI_
    }
    /// HiGHS: ipx/model.h::Model::AIt (line 94)
    pub fn AIt(&self) -> &SparseMatrix {
        &self.AIt_
    }
    /// HiGHS: ipx/model.h::Model::offset (line 97)
    pub fn offset(&self) -> f64 {
        self.offset_
    }
    /// HiGHS: ipx/model.h::Model::b() (line 100)
    pub fn b(&self) -> &Vector {
        &self.b_
    }
    /// HiGHS: ipx/model.h::Model::c() (line 101)
    pub fn c(&self) -> &Vector {
        &self.c_
    }
    /// HiGHS: ipx/model.h::Model::lb() (line 102)
    pub fn lb(&self) -> &Vector {
        &self.lb_
    }
    /// HiGHS: ipx/model.h::Model::ub() (line 103)
    pub fn ub(&self) -> &Vector {
        &self.ub_
    }
    /// HiGHS: ipx/model.h::Model::b(Int) (line 106)
    pub fn b_at(&self, i: Int) -> f64 {
        self.b_[i as usize]
    }
    /// HiGHS: ipx/model.h::Model::c(Int) (line 107)
    pub fn c_at(&self, j: Int) -> f64 {
        self.c_[j as usize]
    }
    /// HiGHS: ipx/model.h::Model::lb(Int) (line 108)
    pub fn lb_at(&self, j: Int) -> f64 {
        self.lb_[j as usize]
    }
    /// HiGHS: ipx/model.h::Model::ub(Int) (line 109)
    pub fn ub_at(&self, j: Int) -> f64 {
        self.ub_[j as usize]
    }
    /// HiGHS: ipx/model.h::Model::norm_bounds (line 112)
    pub fn norm_bounds(&self) -> f64 {
        self.norm_bounds_
    }
    /// HiGHS: ipx/model.h::Model::norm_c (line 115)
    pub fn norm_c(&self) -> f64 {
        self.norm_c_
    }

    /// Test hook (not in HiGHS): every member, for the golden replay; the
    /// oracle reads them with `#define private public`
    /// (tools/highs-oracle/src/ipx.cpp modelObj).
    #[doc(hidden)]
    pub fn oracle_view(&self) -> ModelOracleView<'_> {
        ModelOracleView {
            dualized: self.dualized_,
            num_rows: self.num_rows_,
            num_cols: self.num_cols_,
            num_dense_cols: self.num_dense_cols_,
            nz_dense: self.nz_dense_,
            AI: &self.AI_,
            AIt: &self.AIt_,
            b: &self.b_,
            c: &self.c_,
            lb: &self.lb_,
            ub: &self.ub_,
            norm_bounds: self.norm_bounds_,
            norm_c: self.norm_c_,
            num_constr: self.num_constr_,
            num_eqconstr: self.num_eqconstr_,
            num_var: self.num_var_,
            num_free_var: self.num_free_var_,
            num_entries: self.num_entries_,
            boxed_vars: &self.boxed_vars_,
            constr_type: &self.constr_type_,
            norm_obj: self.norm_obj_,
            norm_rhs: self.norm_rhs_,
            offset: self.offset_,
            scaled_obj: &self.scaled_obj_,
            scaled_rhs: &self.scaled_rhs_,
            scaled_lbuser: &self.scaled_lbuser_,
            scaled_ubuser: &self.scaled_ubuser_,
            A: &self.A_,
            flipped_vars: &self.flipped_vars_,
            colscale: &self.colscale_,
            rowscale: &self.rowscale_,
        }
    }
}

/// C `std::copy_n(src, count, dst)`: copies the first `count` elements.
pub(super) fn copy_n<T: Copy>(src: &[T], count: Int, dst: &mut [T]) {
    let count = count as usize;
    dst[..count].copy_from_slice(&src[..count]);
}

/// C `std::copy(std::begin(src), std::end(src), dst)`: all of `src`.
pub(super) fn copy_all<T: Copy>(src: &[T], dst: &mut [T]) {
    dst[..src.len()].copy_from_slice(src);
}

/// valarray `-v`: elementwise negation (`-x`, so `-0.0` for `+0.0`).
pub(super) fn negated(v: &[f64]) -> Vec<f64> {
    v.iter().map(|&x| -x).collect()
}

/// valarray `x *= s`: `x[i] = x[i] * s[i]` for every `i`.
pub(super) fn mul_assign(x: &mut [f64], s: &[f64]) {
    debug_assert_eq!(x.len(), s.len());
    for (xi, si) in x.iter_mut().zip(s) {
        *xi *= *si;
    }
}

/// valarray `x /= s`: `x[i] = x[i] / s[i]` for every `i`.
pub(super) fn div_assign(x: &mut [f64], s: &[f64]) {
    debug_assert_eq!(x.len(), s.len());
    for (xi, si) in x.iter_mut().zip(s) {
        *xi /= *si;
    }
}

/// See [`Model::oracle_view`].
#[doc(hidden)]
#[allow(non_snake_case)]
pub struct ModelOracleView<'a> {
    pub dualized: bool,
    pub num_rows: Int,
    pub num_cols: Int,
    pub num_dense_cols: Int,
    pub nz_dense: Int,
    pub AI: &'a SparseMatrix,
    pub AIt: &'a SparseMatrix,
    pub b: &'a Vector,
    pub c: &'a Vector,
    pub lb: &'a Vector,
    pub ub: &'a Vector,
    pub norm_bounds: f64,
    pub norm_c: f64,
    pub num_constr: Int,
    pub num_eqconstr: Int,
    pub num_var: Int,
    pub num_free_var: Int,
    pub num_entries: Int,
    pub boxed_vars: &'a [Int],
    pub constr_type: &'a [u8],
    pub norm_obj: f64,
    pub norm_rhs: f64,
    pub offset: f64,
    pub scaled_obj: &'a Vector,
    pub scaled_rhs: &'a Vector,
    pub scaled_lbuser: &'a Vector,
    pub scaled_ubuser: &'a Vector,
    pub A: &'a SparseMatrix,
    pub flipped_vars: &'a [Int],
    pub colscale: &'a Vector,
    pub rowscale: &'a Vector,
}
