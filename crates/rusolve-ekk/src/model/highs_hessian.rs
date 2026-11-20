//! HiGHS source: model/HighsHessian.h, model/HighsHessian.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `HighsHessian()` calls `clear()`, so `Default` is the cleared
//! state: `dim_ = 0`, `format_ = kTriangular`, `start_ = [0]`.

use crate::prelude::*;

/// HiGHS: HighsHessian.h::HighsHessian (lines 21-52)
#[derive(Clone, Debug)]
pub struct HighsHessian {
    pub dim_: HighsInt,
    pub format_: HessianFormat,
    pub start_: Vec<HighsInt>,
    pub index_: Vec<HighsInt>,
    pub value_: Vec<f64>,
}

impl Default for HighsHessian {
    /// HiGHS: HighsHessian.h::HighsHessian::HighsHessian (line 23: `clear()`)
    fn default() -> Self {
        let mut hessian = Self {
            dim_: 0,
            format_: HessianFormat::kTriangular,
            start_: Vec::new(),
            index_: Vec::new(),
            value_: Vec::new(),
        };
        hessian.clear();
        hessian
    }
}

impl PartialEq for HighsHessian {
    /// HiGHS: HighsHessian.cpp::HighsHessian::operator== (lines 173-180)
    fn eq(&self, hessian: &Self) -> bool {
        let mut equal = true;
        equal = self.dim_ == hessian.dim_ && equal;
        equal = self.start_ == hessian.start_ && equal;
        equal = self.index_ == hessian.index_ && equal;
        equal = self.value_ == hessian.value_ && equal;
        equal
    }
}

impl HighsHessian {
    /// HiGHS: HighsHessian.cpp::HighsHessian::clear (lines 16-23)
    ///
    /// Architect-filled: `Highs::clearModel` and `HighsIis::clearData` reach it
    /// on the LP path.
    pub fn clear(&mut self) {
        self.dim_ = 0;
        self.start_.clear();
        self.index_.clear();
        self.value_.clear();
        self.format_ = HessianFormat::kTriangular;
        // C++: this->start_.assign(1, 0) — replaces with [0]
        self.start_.clear();
        self.start_.push(0);
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::exactResize (lines 25-34)
    pub fn exact_resize(&mut self) {
        if self.dim_ != 0 {
            self.start_.resize((self.dim_ + 1) as usize, 0);
            let num_nz = self.start_[self.dim_ as usize];
            self.index_.resize(num_nz as usize, 0);
            self.value_.resize(num_nz as usize, 0.0);
        } else {
            self.clear();
        }
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::deleteCols (lines 36-125)
    pub fn delete_cols(&mut self, index_collection: &HighsIndexCollection) {
        if self.dim_ == 0 {
            return;
        }
        // Can't handle non-triangular matrices yet
        debug_assert!(self.format_ == HessianFormat::kTriangular);
        debug_assert!(ok(index_collection));
        let mut from_k: HighsInt = 0;
        let mut to_k: HighsInt = 0;
        limits(index_collection, &mut from_k, &mut to_k);
        if from_k > to_k {
            return;
        }
        let mut delete_from_col: HighsInt = 0;
        let mut delete_to_col: HighsInt = 0;
        let mut keep_from_col: HighsInt = 0;
        let mut keep_to_col: HighsInt = -1;
        let mut current_set_entry: HighsInt = 0;

        // Initial pass creates a look-up to for the new index of columns
        // being retained, and -1 for columns being deleted
        let mut new_index: Vec<HighsInt>;
        new_index = vec![-1; self.dim_ as usize];
        let mut new_dim: HighsInt = 0;
        let mut k: HighsInt = from_k;
        while k <= to_k {
            update_out_in_index(
                index_collection,
                &mut delete_from_col,
                &mut delete_to_col,
                &mut keep_from_col,
                &mut keep_to_col,
                &mut current_set_entry,
            );
            if k == from_k {
                // Account for the initial columns being kept
                let mut i_col: HighsInt = 0;
                while i_col < delete_from_col {
                    new_index[i_col as usize] = new_dim;
                    new_dim += 1;
                    i_col += 1;
                }
            }
            let mut i_col: HighsInt = keep_from_col;
            while i_col <= keep_to_col {
                new_index[i_col as usize] = new_dim;
                new_dim += 1;
                i_col += 1;
            }
            // When using a mask, to_k = this->dim_, but consecutive
            // keep/delete entries are accumulated, so may not need all passes
            if keep_to_col >= self.dim_ - 1 {
                break;
            }
            k += 1;
        }
        debug_assert!(new_dim < self.dim_);
        // Now perform the pass that deletes rows/columns
        keep_to_col = -1;
        current_set_entry = 0;
        // Have to accumulate new number of entries from the outset, as
        // entries may be lost from columns being kept before any are
        // deleted. Also keep a count of the number of nonzeros, in case the new
        let check_new_dim: HighsInt = new_dim;
        new_dim = 0;
        let mut new_num_nz: HighsInt = 0;
        let mut new_num_entries: HighsInt = 0;
        let save_start: Vec<HighsInt> = self.start_.clone();
        let mut k: HighsInt = from_k;
        while k <= to_k {
            update_out_in_index(
                index_collection,
                &mut delete_from_col,
                &mut delete_to_col,
                &mut keep_from_col,
                &mut keep_to_col,
                &mut current_set_entry,
            );
            if k == from_k {
                // Account for the initial columns being kept
                let mut i_col: HighsInt = 0;
                while i_col < delete_from_col {
                    debug_assert!(new_index[i_col as usize] >= 0);
                    let mut i_el: HighsInt = save_start[i_col as usize];
                    while i_el < save_start[(i_col + 1) as usize] {
                        let i_row: HighsInt = new_index[self.index_[i_el as usize] as usize];
                        if i_row < 0 {
                            i_el += 1;
                            continue;
                        }
                        self.index_[new_num_entries as usize] = i_row;
                        self.value_[new_num_entries as usize] = self.value_[i_el as usize];
                        if self.value_[new_num_entries as usize] != 0.0 {
                            new_num_nz += 1;
                        }
                        new_num_entries += 1;
                        i_el += 1;
                    }
                    new_dim += 1;
                    self.start_[new_dim as usize] = new_num_entries;
                    i_col += 1;
                }
                debug_assert!(new_dim == delete_from_col);
            }
            let mut i_col: HighsInt = keep_from_col;
            while i_col <= keep_to_col {
                debug_assert!(new_index[i_col as usize] >= 0);
                let mut i_el: HighsInt = save_start[i_col as usize];
                while i_el < save_start[(i_col + 1) as usize] {
                    let i_row: HighsInt = new_index[self.index_[i_el as usize] as usize];
                    if i_row < 0 {
                        i_el += 1;
                        continue;
                    }
                    self.index_[new_num_entries as usize] = i_row;
                    self.value_[new_num_entries as usize] = self.value_[i_el as usize];
                    if self.value_[new_num_entries as usize] != 0.0 {
                        new_num_nz += 1;
                    }
                    new_num_entries += 1;
                    i_el += 1;
                }
                new_dim += 1;
                self.start_[new_dim as usize] = new_num_entries;
                i_col += 1;
            }
            if keep_to_col >= self.dim_ - 1 {
                break;
            }
            k += 1;
        }
        debug_assert!(new_dim == check_new_dim);
        self.dim_ = new_dim;
        if new_num_nz == 0 {
            self.clear();
        } else {
            self.exact_resize();
        }
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::scaleOk (lines 127-138)
    pub fn scale_ok(
        &self,
        hessian_scale: HighsInt,
        small_matrix_value: f64,
        large_matrix_value: f64,
    ) -> bool {
        if self.dim_ == 0 {
            return true;
        }
        // C++ std::pow(2, hessian_scale): double pow, never powi
        let hessian_scale_value: f64 = 2.0_f64.powf(hessian_scale as f64);
        let mut i_el: HighsInt = 0;
        while i_el < self.start_[self.dim_ as usize] {
            let abs_new_value: f64 = (self.value_[i_el as usize] * hessian_scale_value).abs();
            if abs_new_value >= large_matrix_value {
                return false;
            }
            if abs_new_value <= small_matrix_value {
                return false;
            }
            i_el += 1;
        }
        true
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::numNz (lines 140-144)
    pub fn num_nz(&self) -> HighsInt {
        debug_assert!(self.format_ok());
        debug_assert!((self.start_.len() as HighsInt) > self.dim_);
        self.start_[self.dim_ as usize]
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::print (lines 146-172)
    ///
    /// C++ default arguments: message = ""
    pub fn print(&self, message: &str) {
        let num_nz: HighsInt = self.num_nz();
        let format_str = if self.format_ == HessianFormat::kTriangular {
            "Triangular"
        } else {
            "Square"
        };
        println!(
            "{} Hessian of dimension {} and {} entries: {}",
            format_str, self.dim_, num_nz, message
        );
        println!(
            "Start; Index; Value of sizes {}; {}; {}",
            self.start_.len(),
            self.index_.len(),
            self.value_.len()
        );
        if self.dim_ <= 0 {
            return;
        }
        print!(" Row|");
        let mut i_col: i32 = 0;
        while i_col < self.dim_ {
            print!(" {:>4}", i_col);
            i_col += 1;
        }
        println!();
        print!("-----");
        let mut i_col: i32 = 0;
        while i_col < self.dim_ {
            print!("-----");
            i_col += 1;
        }
        println!();
        let mut col: Vec<String>;
        col = vec![String::new(); self.dim_ as usize];
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let mut i_el: HighsInt = self.start_[i_col as usize];
            while i_el < self.start_[(i_col + 1) as usize] {
                col[self.index_[i_el as usize] as usize] = format!(
                    "{:>4}",
                    crate::util::c_printf::c_fmt_g(self.value_[i_el as usize], 6)
                );
                i_el += 1;
            }
            print!("{:>4}|", i_col);
            let mut i_row: i32 = 0;
            while i_row < self.dim_ {
                print!(" {:>4}", col[i_row as usize]);
                i_row += 1;
            }
            println!();
            let mut i_el: HighsInt = self.start_[i_col as usize];
            while i_el < self.start_[(i_col + 1) as usize] {
                col[self.index_[i_el as usize] as usize] = String::new();
                i_el += 1;
            }
            i_col += 1;
        }
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::product (lines 182-196)
    pub fn product(&self, solution: &[f64], product: &mut Vec<f64>) {
        if self.dim_ <= 0 {
            return;
        }
        product.clear();
        product.resize(self.dim_ as usize, 0.0);
        let triangular: bool = self.format_ == HessianFormat::kTriangular;
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let mut i_el: HighsInt = self.start_[i_col as usize];
            while i_el < self.start_[(i_col + 1) as usize] {
                let i_row: HighsInt = self.index_[i_el as usize];
                product[i_row as usize] += self.value_[i_el as usize] * solution[i_col as usize];
                if triangular && i_row != i_col {
                    product[i_col as usize] +=
                        self.value_[i_el as usize] * solution[i_row as usize];
                }
                i_el += 1;
            }
            i_col += 1;
        }
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::alphaProductPlusY (lines 198-213)
    pub fn alpha_product_plus_y(&self, alpha: f64, x: &[f64], y: &mut Vec<f64>) {
        if self.dim_ <= 0 {
            return;
        }
        let triangular: bool = self.format_ == HessianFormat::kTriangular;
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let mut i_el: HighsInt = self.start_[i_col as usize];
            while i_el < self.start_[(i_col + 1) as usize] {
                let i_row: HighsInt = self.index_[i_el as usize];
                y[i_row as usize] += alpha * self.value_[i_el as usize] * x[i_col as usize];
                if triangular && i_row != i_col {
                    y[i_col as usize] += alpha * self.value_[i_el as usize] * x[i_row as usize];
                }
                i_el += 1;
            }
            i_col += 1;
        }
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::objectiveValue (lines 215-229)
    pub fn objective_value(&self, solution: &[f64]) -> f64 {
        let mut objective_function_value: f64 = 0.0;
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let i_el: HighsInt = self.start_[i_col as usize];
            // Assumes Hessian format is triangular
            debug_assert!(self.index_[i_el as usize] == i_col);
            objective_function_value += 0.5
                * solution[i_col as usize]
                * self.value_[i_el as usize]
                * solution[i_col as usize];
            let mut i_el: HighsInt = self.start_[i_col as usize] + 1;
            while i_el < self.start_[(i_col + 1) as usize] {
                objective_function_value += solution[i_col as usize]
                    * self.value_[i_el as usize]
                    * solution[self.index_[i_el as usize] as usize];
                i_el += 1;
            }
            i_col += 1;
        }
        objective_function_value
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::objectiveCDoubleValue (lines 231-246)
    pub fn objective_cdouble_value(&self, solution: &[f64]) -> HighsCDouble {
        let mut objective_function_value: HighsCDouble = HighsCDouble::from(0.0_f64);
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let i_el: HighsInt = self.start_[i_col as usize];
            // Assumes Hessian format is triangular
            debug_assert!(self.index_[i_el as usize] == i_col);
            objective_function_value += 0.5
                * solution[i_col as usize]
                * self.value_[i_el as usize]
                * solution[i_col as usize];
            let mut i_el: HighsInt = self.start_[i_col as usize] + 1;
            while i_el < self.start_[(i_col + 1) as usize] {
                objective_function_value += solution[i_col as usize]
                    * self.value_[i_el as usize]
                    * solution[self.index_[i_el as usize] as usize];
                i_el += 1;
            }
            i_col += 1;
        }
        objective_function_value
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::empty (line 248)
    pub fn empty(&self) -> bool {
        self.dim_ <= 0
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::isDiagonal (lines 249-258)
    pub fn is_diagonal(&self) -> bool {
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let mut i_el: HighsInt = self.start_[i_col as usize];
            while i_el < self.start_[(i_col + 1) as usize] {
                let i_row: HighsInt = self.index_[i_el as usize];
                if i_row != i_col {
                    return false;
                }
                i_el += 1;
            }
            i_col += 1;
        }
        true
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::diag (lines 260-264)
    pub fn diag(&self, i: HighsInt) -> f64 {
        debug_assert!(i < self.dim_);
        debug_assert!(self.index_[self.start_[i as usize] as usize] == i);
        self.value_[self.start_[i as usize] as usize]
    }

    /// HiGHS: HighsHessian.cpp::HighsHessian::toSquare (lines 266-319)
    pub fn to_square(&self) -> HighsHessian {
        if self.format_ == HessianFormat::kSquare {
            return self.clone();
        }
        debug_assert!(self.format_ == HessianFormat::kTriangular);
        let mut iwork: Vec<HighsInt> = vec![0; self.dim_ as usize];
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let mut i_el: HighsInt = self.start_[i_col as usize];
            let i_row: HighsInt = self.index_[i_el as usize];
            debug_assert!(i_row == i_col);
            i_el += 1;
            while i_el < self.start_[(i_col + 1) as usize] {
                let i_row: HighsInt = self.index_[i_el as usize];
                debug_assert!(i_row > i_col);
                iwork[i_row as usize] += 1;
                i_el += 1;
            }
            i_col += 1;
        }
        let triangular_hessian_off_diagonal: HighsInt = self.num_nz() - self.dim_;
        let mut square_hessian: HighsHessian = HighsHessian::default();
        square_hessian.format_ = HessianFormat::kSquare;
        square_hessian.dim_ = self.dim_;
        square_hessian.start_[0] = 0;
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let square_hessian_col_nnz: HighsInt = iwork[i_col as usize]
                + self.start_[(i_col + 1) as usize]
                - self.start_[i_col as usize];
            let next_start = square_hessian.start_[i_col as usize] + square_hessian_col_nnz;
            square_hessian.start_.push(next_start);
            iwork[i_col as usize] = square_hessian.start_[i_col as usize] + 1;
            i_col += 1;
        }
        let square_hessian_nnz: HighsInt = square_hessian.num_nz();
        debug_assert!(square_hessian_nnz == self.dim_ + 2 * triangular_hessian_off_diagonal);
        square_hessian.index_.resize(square_hessian_nnz as usize, 0);
        square_hessian
            .value_
            .resize(square_hessian_nnz as usize, 0.0);
        let mut i_col: HighsInt = 0;
        while i_col < self.dim_ {
            let mut i_el: HighsInt = self.start_[i_col as usize];
            let i_row: HighsInt = self.index_[i_el as usize];
            let square_i_el: HighsInt = square_hessian.start_[i_col as usize];
            square_hessian.index_[square_i_el as usize] = i_col;
            square_hessian.value_[square_i_el as usize] = self.value_[i_el as usize];
            i_el += 1;
            while i_el < self.start_[(i_col + 1) as usize] {
                let i_row: HighsInt = self.index_[i_el as usize];
                let value: f64 = self.value_[i_el as usize];
                let mut square_i_el: HighsInt = iwork[i_col as usize];
                square_hessian.index_[square_i_el as usize] = i_row;
                square_hessian.value_[square_i_el as usize] = value;
                iwork[i_col as usize] += 1;
                square_i_el = iwork[i_row as usize];
                square_hessian.index_[square_i_el as usize] = i_col;
                square_hessian.value_[square_i_el as usize] = value;
                iwork[i_row as usize] += 1;
                i_el += 1;
            }
            i_col += 1;
        }
        square_hessian
    }

    /// HiGHS: HighsHessian.h::HighsHessian::formatOk (lines 40-43)
    pub fn format_ok(&self) -> bool {
        self.format_ == HessianFormat::kTriangular || self.format_ == HessianFormat::kSquare
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

    /// [[2, 1], [1, 4]] in triangular (lower, column-wise) format.
    fn h2() -> HighsHessian {
        HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 2, 3],
            index_: vec![0, 1, 1],
            value_: vec![2.0, 1.0, 4.0],
        }
    }

    #[test]
    fn product_objective_2x2() {
        let h = h2();
        assert_eq!(h.num_nz(), 3);
        assert!(h.format_ok());
        assert!(!h.empty());
        assert!(!h.is_diagonal());
        assert_eq!(h.diag(1), 4.0);
        let mut p = Vec::new();
        h.product(&[1.0, 1.0], &mut p);
        assert_eq!(p, vec![3.0, 5.0]);
        // y += 2 Q e0 = (4, 2)
        let mut y = vec![1.0, 1.0];
        h.alpha_product_plus_y(2.0, &[1.0, 0.0], &mut y);
        assert_eq!(y, vec![5.0, 3.0]);
        // 0.5*2 + 1*1 + 0.5*4
        assert_eq!(h.objective_value(&[1.0, 1.0]), 4.0);
        assert!(h == h.clone());
        let sq = h.to_square();
        assert_eq!(sq.format_, HessianFormat::kSquare);
        assert_eq!(sq.num_nz(), 4);
    }

    #[test]
    fn exact_resize_trims() {
        let mut h = h2();
        h.index_.push(7);
        h.value_.push(7.0);
        h.exact_resize();
        assert_eq!(h.start_.len(), 3);
        assert_eq!(h.index_.len(), 3);
        assert_eq!(h.value_.len(), 3);
    }

    #[test]
    fn empty_hessian() {
        let h = HighsHessian::default();
        assert!(h.empty());
        assert_eq!(h.num_nz(), 0);
        assert!(h.is_diagonal());
        assert_eq!(h.objective_value(&[]), 0.0);
        let mut p = Vec::new();
        h.product(&[], &mut p);
        assert!(p.is_empty());
    }

    #[test]
    fn scale_ok_basic() {
        let h = h2();
        assert!(h.scale_ok(0, 1e-10, 1e10));
        assert!(!h.scale_ok(100, 1e-10, 1e10));
    }

    #[test]
    fn diagonal_hessian() {
        let h = HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 1, 2],
            index_: vec![0, 1],
            value_: vec![3.0, 5.0],
        };
        assert!(h.is_diagonal());
        assert_eq!(h.diag(0), 3.0);
        assert_eq!(h.diag(1), 5.0);
        assert_eq!(h.objective_value(&[1.0, 1.0]), 0.5 * 3.0 + 0.5 * 5.0);
    }

    /// [[2, 1, 0], [1, 4, 3], [0, 3, 6]] in triangular format.
    fn h3() -> HighsHessian {
        HighsHessian {
            dim_: 3,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 2, 4, 5],
            index_: vec![0, 1, 1, 2, 2],
            value_: vec![2.0, 1.0, 4.0, 3.0, 6.0],
        }
    }

    #[test]
    fn delete_cols_interval() {
        let mut h = h3();
        let ic = HighsIndexCollection {
            dimension_: 3,
            is_interval_: true,
            from_: 1,
            to_: 1,
            ..Default::default()
        };
        h.delete_cols(&ic);
        assert_eq!(h.dim_, 2);
        assert_eq!(h.start_, vec![0, 1, 2]);
        assert_eq!(h.index_, vec![0, 1]);
        assert_eq!(h.value_, vec![2.0, 6.0]);
    }

    #[test]
    fn delete_cols_set() {
        let mut h = h3();
        let ic = HighsIndexCollection {
            dimension_: 3,
            is_set_: true,
            set_num_entries_: 2,
            set_: vec![0, 2],
            ..Default::default()
        };
        h.delete_cols(&ic);
        assert_eq!(h.dim_, 1);
        assert_eq!(h.start_, vec![0, 1]);
        assert_eq!(h.index_, vec![0]);
        assert_eq!(h.value_, vec![4.0]);
    }

    #[test]
    fn delete_cols_mask() {
        let mut h = h3();
        let ic = HighsIndexCollection {
            dimension_: 3,
            is_mask_: true,
            mask_: vec![1, 0, 0],
            ..Default::default()
        };
        h.delete_cols(&ic);
        assert_eq!(h.dim_, 2);
        assert_eq!(h.start_, vec![0, 2, 3]);
        assert_eq!(h.index_, vec![0, 1, 1]);
        assert_eq!(h.value_, vec![4.0, 3.0, 6.0]);
    }

    #[test]
    fn delete_cols_all_zero_clears() {
        // Remaining column holds only an explicit zero: new_num_nz == 0 -> clear().
        let mut h = HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 1, 2],
            index_: vec![0, 1],
            value_: vec![5.0, 0.0],
        };
        let ic = HighsIndexCollection {
            dimension_: 2,
            is_interval_: true,
            from_: 0,
            to_: 0,
            ..Default::default()
        };
        h.delete_cols(&ic);
        assert_eq!(h.dim_, 0);
        assert_eq!(h.start_, vec![0]);
        assert!(h.index_.is_empty());
        assert!(h.value_.is_empty());
    }

    #[test]
    fn to_square_3x3() {
        let sq = h3().to_square();
        assert_eq!(sq.dim_, 3);
        assert_eq!(sq.start_, vec![0, 2, 5, 7]);
        assert_eq!(sq.index_, vec![0, 1, 1, 0, 2, 2, 1]);
        assert_eq!(sq.value_, vec![2.0, 1.0, 4.0, 1.0, 3.0, 6.0, 3.0]);
    }

    #[test]
    fn objective_cdouble_matches_double() {
        let h = h3();
        let x = [1.0, 2.0, 3.0];
        // 0.5*2 + 1*1*2 + 0.5*2*4*2 + 2*3*3 + 0.5*3*6*3 = 1 + 2 + 8 + 18 + 27
        assert_eq!(h.objective_value(&x), 56.0);
        assert_eq!(f64::from(h.objective_cdouble_value(&x)), 56.0);
    }

    #[test]
    fn to_square_already_square() {
        let mut h = h2();
        h.format_ = HessianFormat::kSquare;
        let sq = h.to_square();
        assert_eq!(sq.format_, HessianFormat::kSquare);
        assert_eq!(sq.dim_, h.dim_);
    }
}
