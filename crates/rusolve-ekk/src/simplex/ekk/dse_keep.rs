//! rusolve `mip_lp_keep_dse_weights`: keep dual steepest-edge weights across a basis or row
//! change instead of recomputing all of them. Not in HiGHS; default off.
//!
//! The dual steepest-edge weight of basic row `r` is `w_r = ||e_r^T B^{-1}||^2` (Forrest and
//! Goldfarb, "Steepest-edge simplex algorithms for linear programming", Math. Programming 57
//! (1992) 341-374; Koberstein, PhD thesis, Paderborn 2005, section 3.3). HiGHS keeps the weights
//! while only bounds change, but `setBasis`, `addRows` and `deleteRows` drop them, and the next
//! dual solve from a non-logical basis computes every weight again with one BTRAN per row
//! (`HEkk::computeDualSteepestEdgeWeights`). In a MIP search these events happen at almost every
//! node LP that is not a plain child: a sibling or a jump restores a stored basis, separation and
//! the LP heuristics add and remove cut rows.
//!
//! This lever saves the weights by variable just before they are dropped (the same scatter by
//! variable that `HEkk::getNonsingularInverse` uses to carry weights across an INVERT that
//! permutes `basic_index_`). The next dual solve that would compute all weights takes the saved
//! weight of every variable that is still basic and runs a BTRAN only for the rows whose basic
//! variable has no saved weight.
//!
//! With the value 2, `HighsLpRelaxation::storeBasis` also stores the weights by variable with
//! the basis it stores for a search node, and `recoverBasis` hands them back, so a sibling or a
//! backtrack starts from the weights of the very basis it restores instead of the weights of the
//! last LP solved. HiGHS has the same idea inside one LP: `Highs::putIterate` and `getIterate`
//! save and restore the basis, the INVERT and the weights together (used by the MIP
//! `HighsLpRelaxation::Playground`).
//!
//! Exactness: if cut rows are appended with their slacks basic, `B' = [B 0; R_B I]` and the old
//! rows of `B'^{-1}` are `[e_r^T B^{-1}, 0]`, so their weights do not change; deleting rows whose
//! slacks are basic is the reverse. The weights of the new rows are computed exactly. After a
//! basis restore that changes `k` basic variables, the kept weights are an approximation of the
//! exact weights of the restored basis. That is valid for pricing: at each iteration the dual
//! simplex computes the exact weight of the chosen row from `row_ep` and rechooses when the kept
//! weight is less than a quarter of it (`HEkkDual::acceptDualSteepestEdgeWeight`). Weights never
//! affect feasibility or the optimality test, only the choice of the leaving row.

use crate::prelude::*;

/// Weights saved by variable, for the variable indexing of the LP they were saved from with the
/// row deletions applied since.
#[derive(Clone, Debug)]
pub struct DseKeep {
    /// Weight per variable (`num_col` columns then the row slacks), 0 when the variable has none.
    pub weight: Vec<f64>,
    /// Number of columns of the saved indexing; -1 when nothing is saved.
    pub num_col: HighsInt,
}

impl Default for DseKeep {
    fn default() -> Self {
        Self {
            weight: Vec::new(),
            num_col: -1,
        }
    }
}

impl HEkk {
    /// Scatters the current weights by variable into `out` (cleared first) and returns the
    /// number of columns, or None when there are no valid weights.
    fn dse_scatter(&self, out: &mut Vec<f64>) -> Option<HighsInt> {
        if !self.status_.has_dual_steepest_edge_weights || !self.status_.has_basis {
            return None;
        }
        let num_row = self.basis_.basic_index_.len();
        let num_tot = self.basis_.nonbasic_flag_.len();
        if num_row == 0 || num_tot <= num_row {
            return None;
        }
        // `Highs::addRows` appends the new basic slacks to the simplex basis before it reaches
        // `HEkk::addRows`, so the basis can have more rows than there are weights: those rows
        // get no saved weight.
        let num_weight = num_row.min(self.dual_edge_weight_.len());
        out.clear();
        out.resize(num_tot, 0.0);
        for r in 0..num_weight {
            let var = self.basis_.basic_index_[r] as usize;
            out[var] = self.dual_edge_weight_[r];
        }
        Some((num_tot - num_row) as HighsInt)
    }

    /// Saves the current weights by variable, when the lever is on and the weights are valid.
    /// Keeps any earlier save otherwise (it is still indexed by variable).
    pub fn dse_keep_save(&mut self) {
        if self.options_.mip_lp_keep_dse_weights <= 0 {
            return;
        }
        let mut weight = std::mem::take(&mut self.dse_keep_.weight);
        // Without valid weights the earlier save goes back untouched.
        if let Some(num_col) = self.dse_scatter(&mut weight) {
            self.dse_keep_.num_col = num_col;
        }
        self.dse_keep_.weight = weight;
    }

    /// The current weights by variable (lever value 2, stored with a search basis).
    pub fn dse_weights_by_variable(&self) -> Option<Vec<f64>> {
        let mut weight = Vec::new();
        self.dse_scatter(&mut weight).map(|_| weight)
    }

    /// Replaces the save by weights stored with a basis of `num_col` columns.
    pub fn dse_keep_set(&mut self, weight: &[f64], num_col: usize) {
        self.dse_keep_.weight.clear();
        self.dse_keep_.weight.extend_from_slice(weight);
        self.dse_keep_.num_col = num_col as HighsInt;
    }

    /// Drops the saved weights.
    pub fn dse_keep_drop(&mut self) {
        self.dse_keep_.num_col = -1;
        self.dse_keep_.weight.clear();
    }

    /// Renumbers the saved slack weights for a row deletion. Only a mask collection (the one
    /// the MIP uses) is renumbered; any other collection drops the save.
    pub fn dse_keep_delete_rows(&mut self, index_collection: &HighsIndexCollection) {
        let keep = &mut self.dse_keep_;
        if keep.num_col < 0 {
            return;
        }
        let num_col = keep.num_col as usize;
        let old_num_row = keep.weight.len() - num_col;
        if !index_collection.is_mask_ || index_collection.mask_.len() < old_num_row {
            self.dse_keep_drop();
            return;
        }
        let mut new_row = num_col;
        for row in 0..old_num_row {
            if index_collection.mask_[row] == 0 {
                keep.weight[new_row] = keep.weight[num_col + row];
                new_row += 1;
            }
        }
        keep.weight.truncate(new_row);
    }

    /// Sets the weights of a non-logical basis from the save, computing only the missing ones.
    /// Returns false (nothing done) when there is no usable save; the caller then computes all
    /// weights as HiGHS does. Consumes the save.
    pub fn dse_keep_apply(&mut self) -> bool {
        if self.options_.mip_lp_keep_dse_weights <= 0 || self.dse_keep_.num_col < 0 {
            return false;
        }
        let num_col = self.lp_.num_col_;
        let num_row = self.lp_.num_row_;
        let num_tot = (num_col + num_row) as usize;
        // Rows can only have been appended since the save (deletions are renumbered).
        if self.dse_keep_.num_col != num_col || self.dse_keep_.weight.len() > num_tot {
            self.dse_keep_drop();
            return false;
        }
        let mut weight = std::mem::take(&mut self.dse_keep_.weight);
        weight.resize(num_tot, 0.0);
        debug_assert!(self.dual_edge_weight_.len() as HighsInt >= num_row);
        let mut row_ep = self.hvector_pool_.take();
        row_ep.setup(num_row);
        for i_row in 0..num_row {
            let var = self.basis_.basic_index_[i_row as usize] as usize;
            let w = weight[var];
            self.dual_edge_weight_[i_row as usize] = if w > 0.0 {
                w
            } else {
                self.compute_dual_steepest_edge_weight(i_row, &mut row_ep)
            };
        }
        self.hvector_pool_.give(row_ep);
        // Consumed: the solve updates the weights from here, and the next save takes them.
        self.dse_keep_.weight = weight;
        self.dse_keep_drop();
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn saved(num_col: usize, weight: &[f64]) -> HEkk {
        let mut ekk = HEkk::default();
        ekk.options_.mip_lp_keep_dse_weights = 1;
        ekk.dse_keep_set(weight, num_col);
        ekk
    }

    #[test]
    fn delete_rows_renumbers_the_kept_slack_weights() {
        let mut ekk = saved(2, &[1.0, 2.0, 10.0, 20.0, 30.0]);
        let mut ic = HighsIndexCollection::default();
        ic.is_mask_ = true;
        ic.mask_ = vec![0, 1, 0];
        ekk.dse_keep_delete_rows(&ic);
        assert_eq!(ekk.dse_keep_.num_col, 2);
        assert_eq!(ekk.dse_keep_.weight, vec![1.0, 2.0, 10.0, 30.0]);
    }

    #[test]
    fn delete_rows_by_set_drops_the_save() {
        let mut ekk = saved(2, &[1.0, 2.0, 10.0]);
        let mut ic = HighsIndexCollection::default();
        ic.is_set_ = true;
        ic.set_ = vec![0];
        ic.set_num_entries_ = 1;
        ekk.dse_keep_delete_rows(&ic);
        assert_eq!(ekk.dse_keep_.num_col, -1);
        assert!(ekk.dse_keep_.weight.is_empty());
    }

    #[test]
    fn lever_off_applies_nothing() {
        let mut ekk = saved(2, &[1.0, 2.0, 10.0]);
        ekk.options_.mip_lp_keep_dse_weights = 0;
        assert!(!ekk.dse_keep_apply());
    }

    #[test]
    fn a_save_for_other_columns_is_dropped() {
        let mut ekk = saved(3, &[1.0, 2.0, 3.0, 10.0]);
        ekk.lp_.num_col_ = 2;
        ekk.lp_.num_row_ = 1;
        assert!(!ekk.dse_keep_apply());
        assert_eq!(ekk.dse_keep_.num_col, -1);
    }
}
