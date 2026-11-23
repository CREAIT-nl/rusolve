//! HiGHS source: presolve/HPresolveZeroCostSingleton.cpp (1.16.0 lines 12-72)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! zeroCostSingleton: relax a zero-cost singleton continuous column out of a
//! ranged row, recording the reduction for postsolve. Gated by
//! `h16_zero_cost_singleton` at the call site in `singletons.rs`.

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS 1.16 (#2962, #3455): HPresolveZeroCostSingleton.cpp::zeroCostSingleton (lines 12-72)
    ///
    /// Row bounds are assigned directly (the #2962 form, not the later
    /// `changeRowLower`/`changeRowUpper` of #3270). No rule logging: the
    /// zero-cost-singleton rule has no `kPresolveRule*` entry in this port.
    pub fn zero_cost_singleton(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        col: HighsInt,
    ) -> HPresolveResult {
        let ci = col as usize;
        // For a double-sided row b_0 <= a^Tx + cs <= b_1,
        // where s is a singleton continuous column with 0 cost,
        // relax s out as its value can be determined in postsolve.
        // The row may now admit additional reductions afterwards.
        // Dual fixing already handles single-sided row case with fixings.
        if self.model.integrality_[ci] == HighsVarType::kInteger
            || self.model.col_cost_[ci] != 0.0
            || self.colsize[ci] != 1
        {
            return HPresolveResult::kOk;
        }
        debug_assert!(self.col_deleted[ci] == 0);

        let nz_pos = self.colhead[ci];
        let row = self.a_row[nz_pos as usize];
        let ri = row as usize;
        debug_assert!(self.row_deleted[ri] == 0);
        // Row must be ranged, but why can't it be an equation?
        let was_equation = self.is_equation(row);
        if !self.is_ranged(row) {
            return HPresolveResult::kOk;
        }

        let coef = self.a_value[nz_pos as usize];

        if coef.abs() == kHighsInf {
            return HPresolveResult::kOk;
        }

        self.store_row(row);

        let lb = self.model.col_lower_[ci];
        let ub = self.model.col_upper_[ci];
        let change_from_col_lb = coef * lb;
        let change_from_col_ub = coef * ub;

        let new_row_lower =
            self.model.row_lower_[ri] - cpp_max(change_from_col_lb, change_from_col_ub);
        let new_row_upper =
            self.model.row_upper_[ri] - cpp_min(change_from_col_lb, change_from_col_ub);

        postsolve_stack.zero_cost_singleton(
            row,
            col,
            self.model.row_lower_[ri],
            self.model.row_upper_[ri],
            new_row_lower,
            new_row_upper,
            lb,
            ub,
            coef,
            self.get_stored_row().iter(self),
        );

        if self.options.h16_row_duals {
            // HiGHS 1.16 #3270 (h16_row_duals): HPresolveZeroCostSingleton.cpp lines 56-57 (1.16.0)
            self.h16_change_row_lower(row, new_row_lower, false);
            self.h16_change_row_upper(row, new_row_upper, false);
        } else {
            self.model.row_lower_[ri] = new_row_lower;
            self.model.row_upper_[ri] = new_row_upper;
        }
        if was_equation && new_row_lower != new_row_upper {
            if let Some(key) = self.eqiters[ri].take() {
                self.equations.remove(&key);
            }
        }

        // Delete the singleton column
        self.mark_col_deleted(col);
        self.unlink(nz_pos);

        // HiGHS 1.16 #3300 (h16_presolve_rules): HPresolveZeroCostSingleton.cpp line 71 (1.16.0)
        // returns checkLimits(postsolve_stack); #2962 as merged returned kOk.
        if self.options.h16_presolve_rules {
            self.check_limits(postsolve_stack)
        } else {
            HPresolveResult::kOk
        }
    }
}

#[cfg(test)]
mod h16_zero_cost_tests {
    use super::*;
    use crate::prelude::*;

    /// One row `row_lower <= x0 + 2 s <= row_upper`, x0 in [0, 10] with cost 1, s in
    /// [0, 3] with cost 0 and the given integrality. Calls zero_cost_singleton on s and
    /// returns (s deleted, new row bounds, number of postsolve reductions).
    fn run_zcs(s_type: HighsVarType, row_upper: f64) -> (bool, f64, f64, usize) {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![1.0, 0.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![10.0, 3.0],
            row_lower_: vec![1.0],
            row_upper_: vec![row_upper],
            ..Default::default()
        };
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 2.0];
        let options = HighsOptions::default();
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, options.presolve_reduction_limit, None));
        // ok_set_input resets an LP's integrality to continuous; restore s's type so the
        // integrality test of zero_cost_singleton is exercised directly
        presolve.model.integrality_[1] = s_type;
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(1, 2);
        assert_eq!(
            presolve.zero_cost_singleton(&mut stack, 1),
            HPresolveResult::kOk
        );
        (
            presolve.col_deleted[1] != 0,
            presolve.model.row_lower_[0],
            presolve.model.row_upper_[0],
            stack.num_reductions(),
        )
    }

    /// HiGHS 1.16 #2962: s is relaxed out of 1 <= x0 + 2 s <= 4; with coef * [lb, ub] = [0, 6]
    /// the C++ gives row bounds [1 - 6, 4 - 0] = [-5, 4] and one postsolve reduction.
    #[test]
    fn h16_zero_cost_continuous_removed() {
        assert_eq!(
            run_zcs(HighsVarType::kContinuous, 4.0),
            (true, -5.0, 4.0, 1)
        );
    }

    /// HiGHS 1.16 #3455: the 1.16 test is `== kInteger`, so an implicit integer is relaxed too
    /// (the #2962 `!= kContinuous` test kept it).
    #[test]
    fn h16_zero_cost_implied_integer_removed() {
        assert_eq!(
            run_zcs(HighsVarType::kImplicitInteger, 4.0),
            (true, -5.0, 4.0, 1)
        );
    }

    /// HiGHS 1.16 #3455: an integer singleton is kept.
    #[test]
    fn h16_zero_cost_integer_kept() {
        assert_eq!(run_zcs(HighsVarType::kInteger, 4.0), (false, 1.0, 4.0, 0));
    }

    /// HiGHS 1.16 #2962: a one-sided row is not ranged, so the column is kept.
    #[test]
    fn h16_zero_cost_one_sided_row_kept() {
        assert_eq!(
            run_zcs(HighsVarType::kContinuous, kHighsInf),
            (false, 1.0, kHighsInf, 0)
        );
    }

    /// HiGHS 1.16 #3300 (h16_presolve_rules): with the 1.16 rule slots the final return is
    /// `check_limits(postsolve_stack)` (HPresolveZeroCostSingleton.cpp line 71, 1.16.0).
    /// With the default `presolve_reduction_limit` of -1 (`reduction_limit` is
    /// `kHighsSize_tInf`) one reduction cannot reach the limit, so `check_limits` returns
    /// kOk and the reduction still happens; with the flag off the #2962 kOk return is
    /// unchanged (run_zcs above).
    #[test]
    fn h16_presolve_rules_zero_cost_returns_check_limits() {
        use crate::util::highs_timer::HighsTimer;
        let mut options = HighsOptions::default();
        options.h16_presolve_rules = true;
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![1.0, 0.0],
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![10.0, 3.0],
            row_lower_: vec![1.0],
            row_upper_: vec![4.0],
            ..Default::default()
        };
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 2.0];
        let mut presolve = HPresolve::default();
        // A timer is passed so the time-limit check of check_limits stays safe even though
        // it is not reached (one reduction: 1 & 1023 != 0).
        assert!(presolve.ok_set_input(
            lp,
            &options,
            options.presolve_reduction_limit,
            Some(HighsTimer::new()),
        ));
        // With h16_presolve_rules on, ok_set_input returns via h16_set_input without
        // building the link structures (presolve_h16.rs: "presolve_h16.rs builds them
        // inside presolve"); the 1.16 driver calls this before any rule runs.
        assert!(presolve.h16_ok_setup_presolve_data_structures());
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(1, 2);
        // check_limits under default limits returns kOk for one reduction:
        // numreductions = 1 < reduction_limit = kHighsSize_tInf.
        assert_eq!(
            presolve.zero_cost_singleton(&mut stack, 1),
            HPresolveResult::kOk
        );
        // The reduction happened exactly as in the flag-off case: s relaxed out of
        // 1 <= x0 + 2 s <= 4 with coef * [lb, ub] = [0, 6] gives row bounds [-5, 4].
        assert_ne!(presolve.col_deleted[1], 0);
        assert_eq!(presolve.model.row_lower_[0], -5.0);
        assert_eq!(presolve.model.row_upper_[0], 4.0);
        assert_eq!(stack.num_reductions(), 1);
    }
}
