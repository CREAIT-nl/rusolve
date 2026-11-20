//! HiGHS source: presolve/HPresolve.cpp (1.16.0 lines 964-1030, 2647-2662)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiGHS 1.16 row-bound updates, behind `h16_row_duals` (HiGHS #3270 with the
//! #3318 loosening-only reset). In 1.16 every presolve
//! write of a row side goes through `changeRowLower` / `changeRowUpper`, which
//! * return early when the side does not change (no `markChangedRow`),
//! * keep the row-dual bounds consistent with the side's finiteness and reset the row-dual
//!   implied bounds derived through the row's columns when a side becomes finite,
//! * reset the column implied bounds derived from the row only when the side loosens.
//!
//! Sums of a side and a shift are taken in `HighsCDouble` (`addToRowLower/Upper`).
//!
//! The call sites read `h16_row_duals` once and keep the 1.15.1 direct writes when it is off.
//! These are the only copies: `mip_probing.rs` (#3318) and `sparsify.rs` (#3264)
//! each carried a private copy of the `skipRowDualUpdate = false` form and now call
//! `h16_change_row_lower/Upper(row, v, false)` (deduplicated, byte-identical bodies,
//! no behaviour change). `mip_probing.rs` likewise wrapped `addToRowLower/Upper` as
//! `add_to_row_lower/upper`; those wrappers are deleted too (again byte-identical),
//! so `finalise_probing` calls `h16_add_to_row_lower/Upper` like every other site.

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS 1.16 (#3270, #3318): HPresolve.cpp::changeRowLower (lines 964-991)
    pub fn h16_change_row_lower(
        &mut self,
        row: HighsInt,
        new_lower: f64,
        skip_row_dual_update: bool,
    ) {
        let old_lower = self.model.row_lower_[row as usize];
        if old_lower == new_lower {
            return;
        }
        self.model.row_lower_[row as usize] = new_lower;

        if !skip_row_dual_update {
            if old_lower == -kHighsInf && new_lower != -kHighsInf {
                // row gained a finite lower bound: dual upper bound loosens from 0
                // to kHighsInf, so implied dual bounds derived through columns in
                // this row may be stale (too tight) and must be reset
                self.change_row_dual_upper(row, kHighsInf);
                let mut cols = Vec::new();
                let mut it = self.get_row_vector(row);
                while let Some(nz) = it.next_nz(self) {
                    cols.push(nz.index());
                }
                for col in cols {
                    self.reset_row_dual_implied_bounds_derived_from_col(col);
                }
            } else if old_lower != -kHighsInf && new_lower == -kHighsInf {
                // row lost its finite lower bound: dual upper bound tightens from
                // kHighsInf to 0; existing implied dual bounds remain valid and new
                // tightening opportunities are picked up by updateRowDualImpliedBounds
                self.change_row_dual_upper(row, 0.0);
            }
        }

        // tightening preserves validity of column implied bounds; only reset on
        // loosening
        if new_lower < old_lower {
            self.reset_col_implied_bounds_derived_from_row(row);
        }
        if self.mipsolver.is_some() {
            self.x1_row_modified(row);
        }
        self.mark_changed_row(row);
    }

    /// HiGHS 1.16 (#3270, #3318): HPresolve.cpp::changeRowUpper (lines 993-1020)
    pub fn h16_change_row_upper(
        &mut self,
        row: HighsInt,
        new_upper: f64,
        skip_row_dual_update: bool,
    ) {
        let old_upper = self.model.row_upper_[row as usize];
        if old_upper == new_upper {
            return;
        }
        self.model.row_upper_[row as usize] = new_upper;

        if !skip_row_dual_update {
            if old_upper == kHighsInf && new_upper != kHighsInf {
                // row gained a finite upper bound: dual lower bound loosens from 0
                // to -kHighsInf, so implied dual bounds derived through columns in
                // this row may be stale (too tight) and must be reset
                self.change_row_dual_lower(row, -kHighsInf);
                let mut cols = Vec::new();
                let mut it = self.get_row_vector(row);
                while let Some(nz) = it.next_nz(self) {
                    cols.push(nz.index());
                }
                for col in cols {
                    self.reset_row_dual_implied_bounds_derived_from_col(col);
                }
            } else if old_upper != kHighsInf && new_upper == kHighsInf {
                // row lost its finite upper bound: dual lower bound tightens from
                // -kHighsInf to 0; existing implied dual bounds remain valid and new
                // tightening opportunities are picked up by updateRowDualImpliedBounds
                self.change_row_dual_lower(row, 0.0);
            }
        }

        // tightening preserves validity of column implied bounds; only reset on
        // loosening
        if new_upper > old_upper {
            self.reset_col_implied_bounds_derived_from_row(row);
        }
        if self.mipsolver.is_some() {
            self.x1_row_modified(row);
        }
        self.mark_changed_row(row);
    }

    /// HiGHS 1.16 (#3270): HPresolve.cpp::addToRowLower (lines 1022-1025)
    ///
    /// Call sites: `aggregator.rs`, `csx.rs`, `sparsify.rs`, and `finalise_probing` in
    /// `mip_probing.rs` (#3318, under `mip_coef_change_impl_bounds`; that
    /// site's former `add_to_row_lower` wrapper was deduplicated into this one, so the C++
    /// `skipRowDualUpdate` argument is passed explicitly as false).
    pub fn h16_add_to_row_lower(&mut self, row: HighsInt, delta: HighsCDouble) {
        if self.model.row_lower_[row as usize] == -kHighsInf {
            return;
        }
        let new_lower = f64::from(self.model.row_lower_[row as usize] + delta);
        self.h16_change_row_lower(row, new_lower, false);
    }

    /// HiGHS 1.16 (#3270): HPresolve.cpp::addToRowUpper (lines 1027-1030)
    ///
    /// Same call sites as `h16_add_to_row_lower`.
    pub fn h16_add_to_row_upper(&mut self, row: HighsInt, delta: HighsCDouble) {
        if self.model.row_upper_[row as usize] == kHighsInf {
            return;
        }
        let new_upper = f64::from(self.model.row_upper_[row as usize] + delta);
        self.h16_change_row_upper(row, new_upper, false);
    }

    /// HiGHS 1.16 (#3218): HPresolve.cpp::rowPresolve (lines 2653, 2660), the
    /// `postsolve_stack.impliedEquation` record of an implied-equation conversion.
    /// It only acts on the dual postsolve, so it rides with the FM postsolve changes
    /// (`h16_fourier_motzkin`); with that flag off the conversion stays unrecorded as in 1.15.1.
    pub fn h16_record_implied_equation(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        row: HighsInt,
        at_lower: bool,
    ) {
        if !self.options.h16_fourier_motzkin {
            return;
        }
        postsolve_stack.implied_equation(row, at_lower, self.get_row_vector(row).iter(self));
    }
}

#[cfg(test)]
mod h16_row_duals_tests {
    use super::*;

    /// 2 rows x 3 columns: row 0 = x0 + x1 <= 4, row 1 = x1 + x2 in [-inf, 6];
    /// x0 and x2 are singleton columns, x1 sits in both rows.
    fn setup() -> HPresolve {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0, 1.0, 1.0];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![10.0; 3];
        lp.integrality_ = vec![HighsVarType::kContinuous; 3];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![4.0, 6.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 1, 3, 4];
        lp.a_matrix_.index_ = vec![0, 0, 1, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0];
        let options = HighsOptions::default();
        let mut presolve = HPresolve::default();
        presolve.analysis_.setup(&lp, &options, true);
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
        presolve
    }

    /// changeRowLower (1.16 lines 964-991): an unchanged side is a no-op (no
    /// markChangedRow); a side that becomes finite opens the row dual's upper bound.
    #[test]
    fn h16_change_row_lower_finiteness_and_noop() {
        let mut p = setup();
        assert_eq!(p.row_dual_upper[0], 0.0);
        p.changed_row_indices.clear();
        p.changed_row_flag[0] = 0;
        p.h16_change_row_lower(0, -kHighsInf, false);
        assert!(p.changed_row_indices.is_empty());
        p.h16_change_row_lower(0, 1.0, false);
        assert_eq!(p.model.row_lower_[0], 1.0);
        assert_eq!(p.row_dual_upper[0], kHighsInf);
        assert_eq!(p.changed_row_indices, vec![0]);
        p.h16_change_row_lower(0, -kHighsInf, false);
        assert_eq!(p.row_dual_upper[0], 0.0);
        // skipRowDualUpdate leaves the row dual bound alone
        p.h16_change_row_lower(0, 2.0, true);
        assert_eq!(p.row_dual_upper[0], 0.0);
    }

    /// addToRowUpper (1.16 lines 1027-1030) skips an infinite side and adds in HighsCDouble.
    #[test]
    fn h16_add_to_row_sides() {
        let mut p = setup();
        p.h16_add_to_row_lower(0, HighsCDouble::from(3.0));
        assert_eq!(p.model.row_lower_[0], -kHighsInf);
        p.h16_add_to_row_upper(0, HighsCDouble::from(0.1) * 3.0);
        assert_eq!(
            p.model.row_upper_[0],
            f64::from(HighsCDouble::from(4.0) + HighsCDouble::from(0.1) * 3.0)
        );
    }

    /// #3367: rowsizeSingleton counts singleton columns per row and follows unlink.
    #[test]
    fn h16_rowsize_singleton_tracks_unlink() {
        let mut p = setup();
        assert_eq!(p.rowsize_singleton, vec![1, 1]);
        // remove x1 from row 0: x1 becomes a singleton column of row 1
        let pos = p.find_nonzero(0, 1);
        p.unlink(pos);
        assert_eq!(p.rowsize_singleton, vec![1, 2]);
        // remove x2 from row 1: x2 is now empty
        let pos = p.find_nonzero(1, 2);
        p.unlink(pos);
        assert_eq!(p.rowsize_singleton, vec![1, 1]);
    }

    /// #3218: ImpliedEquation::undo moves a wrong-sign row dual onto the columns.
    #[test]
    fn h16_implied_equation_undo() {
        let mut sol = HighsSolution::default();
        sol.dual_valid = true;
        sol.row_dual = vec![-1.5];
        sol.col_dual = vec![0.0, 1.0];
        let red = ImpliedEquation {
            row: 0,
            at_lower: true,
        };
        red.undo(&[Nonzero::new(1, 2.0)], &mut sol);
        assert_eq!(sol.row_dual[0], 0.0);
        assert_eq!(sol.col_dual[1], 1.0 + 2.0 * -1.5);
        // right sign for the lower side: unchanged
        sol.row_dual[0] = 0.5;
        red.undo(&[Nonzero::new(1, 2.0)], &mut sol);
        assert_eq!(sol.row_dual[0], 0.5);
    }
}
