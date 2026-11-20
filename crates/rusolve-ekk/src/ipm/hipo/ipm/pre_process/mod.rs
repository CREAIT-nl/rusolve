//! HiGHS source: ipm/hipo/ipm/PreProcess.h, ipm/hipo/ipm/PreProcess.cpp (lines 1-287)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Model preprocessing actions, applied in order by `Preprocessor::apply`
//! and undone in reverse by `Preprocessor::undo`. `PreprocessAction` is the
//! base struct [`PreprocessActionBase`] (`n_pre`, `m_pre`, `n_post`,
//! `m_post`, uninitialised in C++, zero here) plus the trait
//! [`PreprocessAction`]. `std::map<Int, FixedVarsData>` is a `BTreeMap`
//! (same ascending key order). `std::stringstream&` is `&mut String`.
//! Empty rows and fixed variables are here; scaling, formulation, free
//! variables and the `Preprocessor` driver are in `scaling.rs`.

pub mod scaling;

pub use scaling::{PreprocessFormulation, PreprocessFreeVars, PreprocessScaling, Preprocessor};

use super::iterate::Iterate;
use super::model::Model;
use crate::ipm::hipo::Int;
use crate::lp_data::hconst::kHighsInf;
use std::collections::BTreeMap;

/// HiGHS: PreProcess.h `struct PreprocessorPoint` (seven references)
pub struct PreprocessorPoint<'a> {
    pub x: &'a mut Vec<f64>,
    pub xl: &'a mut Vec<f64>,
    pub xu: &'a mut Vec<f64>,
    pub slack: &'a mut Vec<f64>,
    pub y: &'a mut Vec<f64>,
    pub zl: &'a mut Vec<f64>,
    pub zu: &'a mut Vec<f64>,
}

impl PreprocessorPoint<'_> {
    /// HiGHS: PreProcess.cpp::PreprocessorPoint::assertConsistency (lines 9-17)
    pub fn assertConsistency(&self, n: Int, m: Int) {
        debug_assert!(self.x.len() as Int == n);
        debug_assert!(self.xl.len() as Int == n);
        debug_assert!(self.xu.len() as Int == n);
        debug_assert!(self.slack.len() as Int == m);
        debug_assert!(self.y.len() as Int == m);
        debug_assert!(self.zl.len() as Int == n);
        debug_assert!(self.zu.len() as Int == n);
    }
}

/// HiGHS: PreProcess.h `struct PreprocessAction` (data members)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PreprocessActionBase {
    pub n_pre: Int,
    pub m_pre: Int,
    pub n_post: Int,
    pub m_post: Int,
}

/// HiGHS: PreProcess.h `struct PreprocessAction` (virtual methods)
pub trait PreprocessAction {
    fn base(&self) -> &PreprocessActionBase;
    fn base_mut(&mut self) -> &mut PreprocessActionBase;
    fn apply(&mut self, model: &mut Model);
    fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate);
    fn print(&self, stream: &mut String);
}

/// HiGHS: PreProcess.h `struct PreprocessEmptyRows`
/// HiGHS 1.16 (#3376): ipm/hipo/ipm/PreProcess.h::PreprocessEmptyRows (line 41, `b_pre`)
#[derive(Clone, Debug, Default)]
pub struct PreprocessEmptyRows {
    pub base: PreprocessActionBase,
    pub rows_shift: Vec<Int>,
    /// HiGHS 1.16 #3376 (h16_hipo_empty_rows): inert when the flag is off (never read there).
    pub b_pre: Vec<f64>,
    pub empty_rows: Int,
}

impl PreprocessAction for PreprocessEmptyRows {
    fn base(&self) -> &PreprocessActionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut PreprocessActionBase {
        &mut self.base
    }

    /// HiGHS: PreProcess.cpp::PreprocessEmptyRows::apply (lines 19-78)
    /// HiGHS 1.16 (#3376): ipm/hipo/ipm/PreProcess.cpp::PreprocessEmptyRows::apply (lines 19-79; line 28)
    fn apply(&mut self, model: &mut Model) {
        let n = model.n_;
        let m = model.m_;
        let A = &mut model.A_;
        let b = &mut model.b_;
        let constraints = &mut model.constraints_;

        self.base.n_pre = n;
        self.base.m_pre = m;
        // HiGHS 1.16 #3376 (h16_hipo_empty_rows): ipm/hipo/ipm/PreProcess.cpp line 28 (1.16.0)
        if model.h16_hipo_empty_rows_ {
            self.b_pre = b.clone();
        }

        // find empty rows
        let mut entries_per_row = vec![0 as Int; m as usize];
        for col in 0..n {
            let mut el = A.start_[col as usize];
            while el < A.start_[(col + 1) as usize] {
                let row = A.index_[el as usize];
                entries_per_row[row as usize] += 1;
                el += 1;
            }
        }

        self.empty_rows = 0;
        for &i in &entries_per_row {
            if i == 0 {
                self.empty_rows += 1;
            }
        }

        if self.empty_rows > 0 {
            self.rows_shift.clear();
            self.rows_shift.resize(m as usize, 0 as Int);
            for i in 0..m {
                if entries_per_row[i as usize] == 0 {
                    // count how many empty rows there are before a given row
                    for j in (i + 1)..m {
                        self.rows_shift[j as usize] += 1;
                    }
                    self.rows_shift[i as usize] = -1;
                }
            }

            // shift each row index by the number of empty rows before it
            for col in 0..n {
                let mut el = A.start_[col as usize];
                while el < A.start_[(col + 1) as usize] {
                    let row = A.index_[el as usize];
                    A.index_[el as usize] -= self.rows_shift[row as usize];
                    el += 1;
                }
            }
            A.num_row_ -= self.empty_rows;

            // shift entries in b and constraints
            for i in 0..m {
                // ignore entries to be removed
                if self.rows_shift[i as usize] == -1 {
                    continue;
                }

                let shifted_pos = i - self.rows_shift[i as usize];
                b[shifted_pos as usize] = b[i as usize];
                constraints[shifted_pos as usize] = constraints[i as usize];
            }
            b.resize(A.num_row_ as usize, 0.0);
            constraints.resize(A.num_row_ as usize, 0);

            model.m_ = A.num_row_;
        }

        self.base.n_post = model.n_;
        self.base.m_post = model.m_;
    }

    /// HiGHS: PreProcess.cpp::PreprocessEmptyRows::undo (lines 80-107)
    /// HiGHS 1.16 (#3376): ipm/hipo/ipm/PreProcess.cpp::PreprocessEmptyRows::undo (lines 81-111; lines 95-98)
    fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate) {
        point.assertConsistency(self.base.n_post, self.base.m_post);
        if self.empty_rows > 0 {
            // Add Lagrange multiplier for empty rows that were removed
            // Add slack for constraints that were removed

            let mut new_y = vec![0.0; self.base.m_pre as usize];
            let mut new_slack = vec![0.0; self.base.m_pre as usize];

            // position to read from y and slack
            let mut pos: Int = 0;

            for i in 0..self.base.m_pre {
                // ignore shift of empty rows, they will receive a value of 0
                if self.rows_shift[i as usize] == -1 {
                    // HiGHS 1.16 #3376 (h16_hipo_empty_rows): ipm/hipo/ipm/PreProcess.cpp lines 95-98 (1.16.0)
                    if model.h16_hipo_empty_rows_ {
                        new_slack[i as usize] = self.b_pre[i as usize];
                    }
                    continue;
                }

                // re-align value of y and slack, considering empty rows
                new_y[(pos + self.rows_shift[i as usize]) as usize] = point.y[pos as usize];
                new_slack[(pos + self.rows_shift[i as usize]) as usize] = point.slack[pos as usize];
                pos += 1;
            }

            *point.y = new_y;
            *point.slack = new_slack;
        }
        point.assertConsistency(self.base.n_pre, self.base.m_pre);
    }

    /// HiGHS: PreProcess.cpp::PreprocessEmptyRows::print (lines 109-111)
    fn print(&self, stream: &mut String) {
        if self.empty_rows > 0 {
            stream.push_str(&format!("Removed {} empty rows\n", self.empty_rows));
        }
    }
}

/// HiGHS: PreProcess.h `struct PreprocessFixedVars::FixedVarsData`
/// (`c` is uninitialised in C++ until `apply` writes it).
#[derive(Clone, Debug, Default)]
pub struct FixedVarsData {
    pub c: f64,
    pub indA: Vec<Int>,
    pub indQ: Vec<Int>,
    pub valA: Vec<f64>,
    pub valQ: Vec<f64>,
}

/// HiGHS: PreProcess.h `struct PreprocessFixedVars`
#[derive(Clone, Debug, Default)]
pub struct PreprocessFixedVars {
    pub base: PreprocessActionBase,
    pub fixed_vars: Int,
    pub fixed_at: Vec<f64>,
    pub data: BTreeMap<Int, FixedVarsData>,
}

impl PreprocessAction for PreprocessFixedVars {
    fn base(&self) -> &PreprocessActionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut PreprocessActionBase {
        &mut self.base
    }

    /// HiGHS: PreProcess.cpp::PreprocessFixedVars::apply (lines 113-213)
    fn apply(&mut self, model: &mut Model) {
        let n = model.n_;
        let m = model.m_;
        // `model.qp()` is `!Q_.empty()`; Q_ is not modified before the
        // deleteCols below, and afterwards keeps n - fixed_vars > 0 columns.
        let qp = model.qp();
        let A = &mut model.A_;
        let b = &mut model.b_;
        let c = &mut model.c_;
        let lower = &mut model.lower_;
        let upper = &mut model.upper_;
        let Q = &mut model.Q_;
        let offset = &mut model.offset_;

        self.base.n_pre = n;
        self.base.m_pre = m;

        // See "Preprocessing for quadratic programming", Gould, Toint, Math Program
        self.fixed_vars = 0;
        for j in 0..n {
            if lower[j as usize] == upper[j as usize] {
                self.fixed_vars += 1;
            }
            self.data.entry(j).or_default().c = c[j as usize];
        }

        // cannot remove all variables
        if self.fixed_vars == n {
            self.fixed_vars = 0;
        }

        if self.fixed_vars > 0 {
            self.fixed_at.clear();
            self.fixed_at.resize(n as usize, kHighsInf);
            let mut index_to_remove: Vec<Int> = Vec::new();
            for j in 0..n {
                if lower[j as usize] == upper[j as usize] {
                    self.fixed_at[j as usize] = lower[j as usize];
                    index_to_remove.push(j);
                    let xcol = self.fixed_at[j as usize];

                    let dataj = self.data.get_mut(&j).unwrap();

                    *offset += dataj.c * xcol;
                    if qp {
                        *offset += 0.5 * Q.diag(j) * xcol * xcol;
                    }

                    let mut el = A.start_[j as usize];
                    while el < A.start_[(j + 1) as usize] {
                        let row = A.index_[el as usize];
                        let val = A.value_[el as usize];
                        b[row as usize] -= val * xcol;

                        dataj.indA.push(row);
                        dataj.valA.push(val);
                        el += 1;
                    }

                    if qp {
                        for colQ in 0..j {
                            let mut el = Q.start_[colQ as usize];
                            while el < Q.start_[(colQ + 1) as usize] {
                                let rowQ = Q.index_[el as usize];
                                if rowQ == j {
                                    c[colQ as usize] += Q.value_[el as usize] * xcol;

                                    dataj.indQ.push(colQ);
                                    dataj.valQ.push(Q.value_[el as usize]);
                                }
                                el += 1;
                            }
                        }
                        let mut el = Q.start_[j as usize];
                        while el < Q.start_[(j + 1) as usize] {
                            let rowQ = Q.index_[el as usize];
                            c[rowQ as usize] += Q.value_[el as usize] * xcol;

                            dataj.indQ.push(rowQ);
                            dataj.valQ.push(Q.value_[el as usize]);
                            el += 1;
                        }
                    }
                }
            }

            let mut index_collection = crate::util::highs_utils::HighsIndexCollection::default();
            crate::util::highs_utils::create_set(
                &mut index_collection,
                index_to_remove.len() as Int,
                &index_to_remove,
                n,
            );
            A.delete_cols(&index_collection);
            if qp {
                Q.delete_cols(&index_collection);
            }

            let mut next: Int = 0;
            let mut copy_to: Int = 0;
            for i in 0..n {
                if next < index_to_remove.len() as Int && i == index_to_remove[next as usize] {
                    next += 1;
                    continue;
                } else {
                    c[copy_to as usize] = c[i as usize];
                    lower[copy_to as usize] = lower[i as usize];
                    upper[copy_to as usize] = upper[i as usize];
                    copy_to += 1;
                }
            }

            model.n_ -= self.fixed_vars;
            debug_assert!(A.num_col_ == model.n_);
            debug_assert!(!qp || Q.dim_ == model.n_);
            c.resize(model.n_ as usize, 0.0);
            lower.resize(model.n_ as usize, 0.0);
            upper.resize(model.n_ as usize, 0.0);
        }

        self.base.n_post = model.n_;
        self.base.m_post = model.m_;
    }

    /// HiGHS: PreProcess.cpp::PreprocessFixedVars::undo (lines 215-282)
    fn undo(&self, point: &mut PreprocessorPoint<'_>, model: &Model, it: &Iterate) {
        point.assertConsistency(self.base.n_post, self.base.m_post);

        if self.fixed_vars > 0 {
            // Add primal and dual variables for fixed variables

            let mut new_x = vec![0.0; self.base.n_pre as usize];
            let mut new_xl = vec![0.0; self.base.n_pre as usize];
            let mut new_xu = vec![0.0; self.base.n_pre as usize];
            let mut new_zl = vec![0.0; self.base.n_pre as usize];
            let mut new_zu = vec![0.0; self.base.n_pre as usize];

            let mut pos: Int = 0;
            for j in 0..self.base.n_pre {
                if self.fixed_at[j as usize].is_finite() {
                    new_x[j as usize] = self.fixed_at[j as usize];
                    new_xl[j as usize] = 0.0;
                    new_xu[j as usize] = 0.0;
                } else {
                    new_x[j as usize] = point.x[pos as usize];
                    new_xl[j as usize] = point.xl[pos as usize];
                    new_xu[j as usize] = point.xu[pos as usize];
                    new_zl[j as usize] = point.zl[pos as usize];
                    new_zu[j as usize] = point.zu[pos as usize];
                    pos += 1;
                }
            }

            for j in 0..self.base.n_pre {
                if self.fixed_at[j as usize].is_finite() {
                    // compute dual variables so that they are dual feasible
                    // z = c - A^T * y + Q * x
                    // need to do this after all x have been computed, due to Q*x term
                    let dataj = self.data.get(&j).unwrap();
                    let mut z = dataj.c;
                    for i in 0..dataj.indA.len() as Int {
                        let row = dataj.indA[i as usize];
                        let val = dataj.valA[i as usize];
                        z -= val * point.y[row as usize];
                    }
                    if model.qp() {
                        for i in 0..dataj.indQ.len() as Int {
                            let row = dataj.indQ[i as usize];
                            let val = dataj.valQ[i as usize];
                            z += val * new_x[row as usize];
                        }
                    }

                    if z >= 0.0 {
                        new_zl[j as usize] = z;
                        new_zu[j as usize] = 0.0;
                    } else {
                        new_zl[j as usize] = 0.0;
                        new_zu[j as usize] = -z;
                    }
                }
            }

            *point.x = new_x;
            *point.xl = new_xl;
            *point.xu = new_xu;
            *point.zl = new_zl;
            *point.zu = new_zu;
        }

        point.assertConsistency(self.base.n_pre, self.base.m_pre);
    }

    /// HiGHS: PreProcess.cpp::PreprocessFixedVars::print (lines 284-287)
    fn print(&self, stream: &mut String) {
        if self.fixed_vars > 0 {
            stream.push_str(&format!("Removed {} fixed variables\n", self.fixed_vars));
        }
    }
}

#[cfg(test)]
mod h16_tests {
    //! HiGHS 1.16 #3376 tests. Port of check/TestHipo.cpp "issue-3369" at the model level:
    //! `min x` s.t. `x >= 1` and an empty row `0 <= 5`. The empty row's slack after undo
    //! is its rhs with the flag on, and 0 with it off (1.15.1).
    use super::*;
    use crate::util::sparse_matrix::HighsSparseMatrix;

    fn h16_empty_row_model(flag: bool) -> Model {
        Model {
            n_: 1,
            m_: 2,
            b_: vec![1.0, 5.0],
            constraints_: vec![0, 0],
            A_: HighsSparseMatrix {
                num_col_: 1,
                num_row_: 2,
                start_: vec![0, 1],
                index_: vec![0],
                value_: vec![1.0],
                ..Default::default()
            },
            h16_hipo_empty_rows_: flag,
            ..Default::default()
        }
    }

    /// Applies the empty-row preprocessing, then undoes it on a point with
    /// slack = [0.25] and y = [0.5]. Returns the slack after undo (length m_pre = 2).
    fn h16_apply_undo_slack(flag: bool) -> Vec<f64> {
        let mut model = h16_empty_row_model(flag);
        let mut pre = PreprocessEmptyRows::default();
        pre.apply(&mut model);
        assert_eq!(pre.empty_rows, 1);
        assert_eq!(model.m_, 1);
        let it = Iterate::new(&model);
        let mut x = vec![0.0];
        let mut xl = vec![0.0];
        let mut xu = vec![0.0];
        let mut zl = vec![0.0];
        let mut zu = vec![0.0];
        let mut slack = vec![0.25];
        let mut y = vec![0.5];
        let mut point = PreprocessorPoint {
            x: &mut x,
            xl: &mut xl,
            xu: &mut xu,
            slack: &mut slack,
            y: &mut y,
            zl: &mut zl,
            zu: &mut zu,
        };
        pre.undo(&mut point, &model, &it);
        slack
    }

    #[test]
    fn h16_hipo_empty_rows_on_postsolves_empty_row_slack() {
        // Flag on: the empty row (index 1, rhs 5) gets slack = b; row 0 keeps its slack.
        assert_eq!(h16_apply_undo_slack(true), vec![0.25, 5.0]);
    }

    #[test]
    fn h16_hipo_empty_rows_off_keeps_zero_slack() {
        // Flag off: 1.15.1 behaviour, the empty row's slack stays 0.
        assert_eq!(h16_apply_undo_slack(false), vec![0.25, 0.0]);
    }
}
