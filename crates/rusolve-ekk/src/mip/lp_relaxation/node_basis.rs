//! rusolve `mip_node_install_basis`: a node taken from the node queue starts its LP from the
//! basis its parent stored, instead of the basis the LP relaxation happens to hold.
//!
//! HiGHS installs a queued node with no basis (`HighsSearch::installNode`), so the node's first
//! LP starts from the basis of whatever the worker solved last, often a node far away in the
//! tree. SCIP keeps the LP state of each fork node and loads it when a node below the fork is
//! focused (`SCIPtreeLoadLPState` in src/scip/tree.c, SCIP 11.0.0, Apache-2.0).
//!
//! Here the open node keeps a shared reference to the basis its parent stored (the same
//! `Arc` that siblings already restore in `backtrack`). The cut rows of the LP change while the
//! node waits in the queue, so `storeBasis` also records the identity of each cut row, and the
//! stored statuses are matched to the current rows on install:
//! - model rows keep their position and status;
//! - a cut row still in the LP takes its stored status;
//! - a cut row added since the basis was stored is basic (its slack).
//!
//! A stored cut row that has left the LP drops out with its status. If its slack was basic the
//! count of basic variables still equals the row count. If it was nonbasic (the slack became
//! basic later and the cut was aged out) there are more basic variables than rows, and the
//! stored basis is not used: the LP keeps the basis it holds, as in HiGHS. (Passing such a
//! basis as an alien basis, which `Highs::setBasis` completes to a nonsingular one, was
//! measured and made the first LP of a node longer on two of six test instances.)

use super::*;
use std::collections::HashMap;

impl HighsLpRelaxation {
    /// Number of model rows: the leading rows of origin `kModel` (cut rows always follow them).
    fn num_leading_model_rows(&self) -> usize {
        self.lprows
            .partition_point(|r| r.origin == HighsLpRelaxation_LpRow::kModel)
    }

    /// `[origin, cutpoolindex, index]` of each LP row after the model rows.
    pub(crate) fn cut_row_keys(&self) -> Vec<[HighsInt; 3]> {
        let nmodel = self.num_leading_model_rows();
        self.lprows[nmodel..]
            .iter()
            .map(|r| [r.origin, r.cutpoolindex, r.index])
            .collect()
    }

    /// Install `basis` (stored by `store_basis` with row keys) on the current LP, its cut rows
    /// matched by identity. Returns false, leaving the LP basis unchanged, when the basis has
    /// no row keys, does not fit the model, or the matched basis has more basic variables than
    /// rows.
    pub(crate) fn restore_node_basis(&mut self, basis: &HighsBasis) -> bool {
        let Some(keys) = basis.row_keys.as_deref() else {
            return false;
        };
        let nrow = self.lprows.len();
        let nmodel = self.num_leading_model_rows();
        if !basis.valid
            || basis.col_status.len() != self.lpsolver.get_num_col() as usize
            || basis.row_status.len() != nmodel + keys.len()
        {
            return false;
        }
        let current = self.cut_row_keys();
        let row_status = if current.as_slice() == keys {
            basis.row_status.clone()
        } else {
            let stored: HashMap<[HighsInt; 3], HighsBasisStatus> = keys
                .iter()
                .copied()
                .zip(basis.row_status[nmodel..].iter().copied())
                .collect();
            let mut row_status = Vec::with_capacity(nrow);
            row_status.extend_from_slice(&basis.row_status[..nmodel]);
            row_status.extend(
                current
                    .iter()
                    .map(|k| stored.get(k).copied().unwrap_or(HighsBasisStatus::kBasic)),
            );
            row_status
        };
        let nbasic = basis
            .col_status
            .iter()
            .chain(row_status.iter())
            .filter(|&&s| s == HighsBasisStatus::kBasic)
            .count();
        if nbasic != nrow {
            return false;
        }
        let restored = HighsBasis {
            valid: true,
            alien: basis.alien,
            useful: basis.useful,
            was_alien: basis.was_alien,
            debug_origin_name: "HighsLpRelaxation::restoreNodeBasis".to_string(),
            col_status: basis.col_status.clone(),
            row_status,
            ..Default::default()
        };
        self.set_stored_basis(Some(Arc::new(restored)));
        self.recover_basis();
        true
    }
}
