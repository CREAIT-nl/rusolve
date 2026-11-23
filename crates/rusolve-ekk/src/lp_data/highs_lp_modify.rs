//! HiGHS source: lp_data/HighsLp.cpp (deleteRowsFromVectors, deleteRows),
//! lp_data/HighsLpUtils.cpp (deleteScale, changeLpCosts, changeLpColBounds, changeBounds)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice: the in-place LP edits behind
//! `Highs::deleteRows` and `Highs::changeColsBounds`. All index walks use
//! `limits` + `update_out_in_index` from util/highs_utils.rs exactly as C++
//! does (the `keep_to_* = -1` / `current_set_entry = 0` initialisation and
//! both early `break`s matter for masks and sets).

use crate::prelude::*;

impl HighsLp {
    /// HiGHS: HighsLp.cpp::HighsLp::deleteRowsFromVectors (lines 373-412, 40 lines)
    ///
    /// `row_names_` is compacted iff it is non-empty (`have_names`), as in C++.
    pub fn delete_rows_from_vectors(
        &mut self,
        new_num_row: &mut HighsInt,
        index_collection: &HighsIndexCollection,
    ) {
        debug_assert!(ok(index_collection));
        let mut from_k: HighsInt = 0;
        let mut to_k: HighsInt = 0;
        limits(index_collection, &mut from_k, &mut to_k);
        // Initialise new_num_row in case none is removed due to from_k > to_k
        *new_num_row = self.num_row_;
        if from_k > to_k {
            return;
        }

        let mut delete_from_row: HighsInt = 0;
        let mut delete_to_row: HighsInt = 0;
        let mut keep_from_row: HighsInt = 0;
        let mut keep_to_row: HighsInt = -1;
        let mut current_set_entry: HighsInt = 0;

        let row_dim: HighsInt = self.num_row_;
        *new_num_row = 0;
        let have_names: bool = self.row_names_.len() > 0;
        for k in from_k..=to_k {
            update_out_in_index(
                index_collection,
                &mut delete_from_row,
                &mut delete_to_row,
                &mut keep_from_row,
                &mut keep_to_row,
                &mut current_set_entry,
            );
            if k == from_k {
                // Account for the initial rows being kept
                *new_num_row = delete_from_row;
            }
            if delete_to_row >= row_dim - 1 {
                break;
            }
            debug_assert!(delete_to_row < row_dim);
            for row in keep_from_row..=keep_to_row {
                self.row_lower_[*new_num_row as usize] = self.row_lower_[row as usize];
                self.row_upper_[*new_num_row as usize] = self.row_upper_[row as usize];
                if have_names {
                    self.row_names_[*new_num_row as usize] = self.row_names_[row as usize].clone();
                }
                *new_num_row += 1;
            }
            if keep_to_row >= row_dim - 1 {
                break;
            }
        }
        self.row_lower_.resize(*new_num_row as usize, 0.0);
        self.row_upper_.resize(*new_num_row as usize, 0.0);
        if have_names {
            self.row_names_.resize(*new_num_row as usize, String::new());
        }
    }

    /// HiGHS: HighsLp.cpp::HighsLp::deleteRows (lines 421-426, 6 lines)
    ///
    /// `this->a_matrix_.deleteRows(index_collection)` is the existing
    /// `HighsSparseMatrix::delete_rows` (util/sparse_matrix/modify.rs:402).
    pub fn delete_rows(&mut self, index_collection: &HighsIndexCollection) {
        let mut new_num_row: HighsInt = 0;
        self.delete_rows_from_vectors(&mut new_num_row, index_collection);
        self.a_matrix_.delete_rows(index_collection);
        self.num_row_ = new_num_row;
    }
}

/// HiGHS: HighsLpUtils.cpp::deleteScale (lines 1636-1666, 31 lines)
///
/// Compacts `scale` in place; it does NOT resize (the caller,
/// `deleteRowsInterface`, resizes afterwards).
pub fn delete_scale(scale: &mut [f64], index_collection: &HighsIndexCollection) {
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

    let col_dim: HighsInt = index_collection.dimension_;
    let mut new_num_col: HighsInt = 0;
    for k in from_k..=to_k {
        update_out_in_index(
            index_collection,
            &mut delete_from_col,
            &mut delete_to_col,
            &mut keep_from_col,
            &mut keep_to_col,
            &mut current_set_entry,
        );
        // Account for the initial columns being kept
        if k == from_k {
            new_num_col = delete_from_col;
        }
        if delete_to_col >= col_dim - 1 {
            break;
        }
        debug_assert!(delete_to_col < col_dim);
        for col in keep_from_col..=keep_to_col {
            debug_assert!((scale.len() as HighsInt) > new_num_col);
            scale[new_num_col as usize] = scale[col as usize];
            new_num_col += 1;
        }
        if keep_to_col >= col_dim - 1 {
            break;
        }
    }
}

/// HiGHS: HighsLpUtils.cpp::changeLpCosts (lines 1764-1799, 36 lines)
///
/// `usr_col` starts at -1 and is incremented for an interval, set to `k`
/// otherwise; a mask entry of 0 is skipped (`continue`), as in `change_bounds`.
pub fn change_lp_costs(
    lp: &mut HighsLp,
    index_collection: &HighsIndexCollection,
    new_col_cost: &[f64],
    infinite_cost: f64,
) {
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    if from_k > to_k {
        return;
    }

    let interval: bool = index_collection.is_interval_;
    let mask: bool = index_collection.is_mask_;
    let col_set: &[HighsInt] = &index_collection.set_;
    let col_mask: &[HighsInt] = &index_collection.mask_;

    // Change the costs to the user-supplied costs, according to the technique
    let mut usr_col: HighsInt = -1;
    for k in from_k..to_k + 1 {
        let lp_col: HighsInt = if interval || mask {
            k
        } else {
            col_set[k as usize]
        };
        let col: HighsInt = lp_col;
        if interval {
            usr_col += 1;
        } else {
            usr_col = k;
        }
        if mask && col_mask[col as usize] == 0 {
            continue;
        }
        lp.col_cost_[col as usize] = new_col_cost[usr_col as usize];
    }
    // Check whether the LP still has an infinite cost
    if lp.has_infinite_cost_ {
        lp.has_infinite_cost_ = lp.has_infinite_cost(infinite_cost);
    }
}

/// HiGHS: HighsLpUtils.cpp::changeLpColBounds (lines 1801-1807, 7 lines)
pub fn change_lp_col_bounds(
    lp: &mut HighsLp,
    index_collection: &HighsIndexCollection,
    new_col_lower: &[f64],
    new_col_upper: &[f64],
) {
    change_bounds(
        &mut lp.col_lower_,
        &mut lp.col_upper_,
        index_collection,
        new_col_lower,
        new_col_upper,
    );
}

/// HiGHS: HighsLpUtils.cpp::changeBounds (lines 1817-1851, 35 lines)
///
/// `usr_ix` starts at -1 and is incremented for an interval, set to `k`
/// otherwise; a mask entry of 0 is skipped (`continue`).
pub fn change_bounds(
    lower: &mut [f64],
    upper: &mut [f64],
    index_collection: &HighsIndexCollection,
    new_lower: &[f64],
    new_upper: &[f64],
) {
    debug_assert!(ok(index_collection));
    let mut from_k: HighsInt = 0;
    let mut to_k: HighsInt = 0;
    limits(index_collection, &mut from_k, &mut to_k);
    if from_k > to_k {
        return;
    }

    let interval: bool = index_collection.is_interval_;
    let mask: bool = index_collection.is_mask_;
    let ix_set: &[HighsInt] = &index_collection.set_;
    let ix_mask: &[HighsInt] = &index_collection.mask_;

    // Change the bounds to the user-supplied bounds, according to the technique
    let mut usr_ix: HighsInt = -1;
    for k in from_k..to_k + 1 {
        let lp_ix: HighsInt = if interval || mask {
            k
        } else {
            ix_set[k as usize]
        };
        let ix: HighsInt = lp_ix;
        if interval {
            usr_ix += 1;
        } else {
            usr_ix = k;
        }
        if mask && ix_mask[ix as usize] == 0 {
            continue;
        }
        lower[ix as usize] = new_lower[usr_ix as usize];
        upper[ix as usize] = new_upper[usr_ix as usize];
    }
}
