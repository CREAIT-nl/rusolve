//! HiGHS 1.16 #2962 (h16_zero_cost_singleton): postsolve of a zero-cost singleton
//! continuous column removed from a ranged row.
//!
//! HiGHS source: presolve/HighsPostsolveStack.cpp::ZeroObjSingletonContinuousCol::undo
//! (1.16.0 lines 1675-1762). The `dual_tol` local of the C++ is unused and omitted.

use super::*;
use crate::prelude::*;

impl ZeroObjSingletonContinuousCol {
    /// HiGHS 1.16 (#2962): HighsPostsolveStack.cpp::ZeroObjSingletonContinuousCol::undo (lines 1675-1762)
    pub fn undo(
        &self,
        options: &HighsOptions,
        row_values: &[Nonzero],
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
    ) {
        debug_assert!(self.orig_row_lower != -kHighsInf && self.orig_row_upper != kHighsInf);

        let primal_tol = options.primal_feasibility_tolerance;
        let row = self.row as usize;
        let col = self.col as usize;
        let coef = self.coef;
        let lb = self.lb;
        let ub = self.ub;

        // Get activity of row without removed singleton
        let mut act = HighsCDouble::from(0.0);
        solution.col_value[col] = 0.0;
        for row_val in row_values {
            act += row_val.value * solution.col_value[row_val.index as usize];
        }

        // Find a suitable bound within the domain
        let mut col_value = kHighsInf;
        // Detemine whether the row was at its lower or upper bound by
        // virtue of basis status or value
        let mut row_at_lower = (basis.valid && basis.row_status[row] == HighsBasisStatus::kLower)
            || f64::from(act - primal_tol) <= self.new_row_lb;
        let row_at_upper = (basis.valid && basis.row_status[row] == HighsBasisStatus::kUpper)
            || f64::from(act + primal_tol) >= self.new_row_ub;
        let mut col_at_lower = false;
        let mut col_at_upper = false;
        if row_at_lower {
            if coef > 0.0 {
                col_at_upper = true;
            } else {
                col_at_lower = true;
            }
        } else if row_at_upper {
            if coef > 0.0 {
                col_at_lower = true;
            } else {
                col_at_upper = true;
            }
        } else {
            // Determine whether the column can be at its lower or upper bound
            // with the row between its bound basic
            let act_with_col_at_lower = f64::from(act + coef * lb);
            let act_with_col_at_upper = f64::from(act + coef * ub);
            if act_with_col_at_lower >= self.orig_row_lower - primal_tol
                && act_with_col_at_lower <= self.orig_row_upper + primal_tol
            {
                col_at_lower = true;
            } else if act_with_col_at_upper >= self.orig_row_lower - primal_tol
                && act_with_col_at_upper <= self.orig_row_upper + primal_tol
            {
                col_at_upper = true;
            } else {
                // Otherwise, there will be a column value for which the row is at a bound
                let col_value_for_lower = f64::from((self.orig_row_lower - act) / coef);
                let col_value_for_upper = f64::from((self.orig_row_upper - act) / coef);
                if col_value_for_lower >= lb - primal_tol && col_value_for_lower <= ub + primal_tol
                {
                    // Can set column basic at this value
                    col_value = col_value_for_lower;
                    row_at_lower = true;
                } else if col_value_for_upper >= lb - primal_tol
                    && col_value_for_upper <= ub + primal_tol
                {
                    // C++ also sets row_at_upper = true here, but it is never read again
                    col_value = col_value_for_upper;
                }
            }
        }
        if col_at_lower || col_at_upper {
            col_value = if col_at_lower { lb } else { ub };
            if solution.dual_valid {
                solution.col_dual[col] = -coef * solution.row_dual[row];
            }
            if basis.valid {
                basis.col_status[col] = if col_at_lower {
                    HighsBasisStatus::kLower
                } else {
                    HighsBasisStatus::kUpper
                };
            }
        } else {
            // Column takes value between its bounds and inherits any basic status from
            // the row
            if solution.dual_valid {
                solution.col_dual[col] = 0.0;
                solution.row_dual[row] = 0.0;
            }
            if basis.valid {
                basis.row_status[row] = if row_at_lower {
                    HighsBasisStatus::kLower
                } else {
                    HighsBasisStatus::kUpper
                };
                basis.col_status[col] = HighsBasisStatus::kBasic;
            }
        }
        solution.col_value[col] = col_value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Row 0: 1 <= x0 + 2 s <= 4 with s in [0, 3], after relaxing s the row is -5 <= x0 <= 4.
    fn reduction() -> ZeroObjSingletonContinuousCol {
        ZeroObjSingletonContinuousCol {
            orig_row_lower: 1.0,
            orig_row_upper: 4.0,
            new_row_lb: -5.0,
            new_row_ub: 4.0,
            lb: 0.0,
            ub: 3.0,
            coef: 2.0,
            col: 1,
            row: 0,
        }
    }

    fn undo_with_x0(x0: f64) -> f64 {
        let row_values = vec![Nonzero::new(0, 1.0), Nonzero::new(1, 2.0)];
        let mut solution = HighsSolution {
            col_value: vec![x0, 0.0],
            row_value: vec![0.0],
            ..Default::default()
        };
        let mut basis = HighsBasis::default();
        let options = HighsOptions::default();
        reduction().undo(&options, &row_values, &mut solution, &mut basis);
        solution.col_value[1]
    }

    /// HiGHS 1.16 #2962 (h16_zero_cost_singleton): the row is not at a bound, so s takes the
    /// value that puts the row activity on its original lower bound.
    #[test]
    fn h16_zero_cost_undo_interior_value() {
        assert!((undo_with_x0(0.0) - 0.5).abs() < 1e-12);
    }

    /// HiGHS 1.16 #2962: the row activity with s at its lower bound already lies in the
    /// original row range, so s takes its lower bound.
    #[test]
    fn h16_zero_cost_undo_column_at_lower_bound() {
        assert_eq!(undo_with_x0(2.0), 0.0);
    }

    /// HiGHS 1.16 #2962: with x0 = -3 only s = ub puts the activity inside [1, 4].
    #[test]
    fn h16_zero_cost_undo_column_at_upper_bound() {
        assert_eq!(undo_with_x0(-3.0), 3.0);
    }

    /// test-dual-fixing-primal-dual-postsolve (3342.diff): the C++ comment's step 1 is a
    /// zero objective singleton continuous column removing x from w - x - y - z = -5.
    /// Unique optimum w = 7, x = 9, y = 3, z = 0, objective -56.
    fn dual_fixing_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 5,
            num_row_: 3,
            col_cost_: vec![0.0, -8.0, 0.0, 0.0, 1.0],
            col_lower_: vec![-kHighsInf, -kHighsInf, 0.0, -kHighsInf, 0.0],
            col_upper_: vec![kHighsInf, 7.0, 9.0, 3.0, 1.0],
            row_lower_: vec![-kHighsInf, -5.0, 0.0],
            row_upper_: vec![8.0, -5.0, kHighsInf],
            sense_: ObjSense::kMinimize,
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 5;
        lp.a_matrix_.num_row_ = 3;
        lp.a_matrix_.start_ = vec![0, 2, 4, 5, 7, 8];
        lp.a_matrix_.index_ = vec![0, 2, 1, 2, 1, 0, 1, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, -1.0, -1.0, -2.0, -1.0, -1.0];
        lp
    }

    fn new_highs(solver: &str, crossover: &str, zcs: bool, presolve: &str) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("solver", solver),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("run_crossover", crossover),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_str("presolve", presolve),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_zero_cost_singleton", zcs),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(dual_fixing_lp()), HighsStatus::kOk);
        highs
    }

    fn num_zcs_reductions(solver: &str, crossover: &str, zcs: bool) -> usize {
        let mut highs = new_highs(solver, crossover, zcs, "on");
        let _ = highs.presolve();
        highs
            .presolve_
            .data_
            .post_solve_stack
            .reductions
            .iter()
            .filter(|r| r.0 == ReductionType::kZeroObjSingletonContinuousCol)
            .count()
    }

    /// HiGHS 1.16 #2962 (h16_zero_cost_singleton): presolve pushes the new reduction only with
    /// the option on, and not under IPM without crossover (mayRequireDualPostsolve).
    #[test]
    fn h16_zero_cost_presolve_pushes_reduction() {
        assert!(num_zcs_reductions("simplex", "on", true) > 0);
        assert_eq!(num_zcs_reductions("simplex", "on", false), 0);
        assert_eq!(num_zcs_reductions("ipm", "off", true), 0);
    }

    /// HiGHS 1.16 #2962/#3306: the postsolved primal and dual solution with the reduction
    /// matches the solve without presolve.
    #[test]
    fn h16_zero_cost_postsolve_matches_presolve_off() {
        let mut on = new_highs("simplex", "on", true, "on");
        let _ = on.run();
        let mut reference = new_highs("simplex", "on", false, "off");
        let _ = reference.run();
        for h in [&on, &reference] {
            assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
            assert!((h.get_info().objective_function_value + 56.0).abs() < 1e-6);
            assert_eq!(h.get_info().num_primal_infeasibilities, 0);
            assert_eq!(h.get_info().num_dual_infeasibilities, 0);
        }
        let (a, b) = (on.get_solution(), reference.get_solution());
        for j in 1..5 {
            assert!(
                (a.col_value[j] - b.col_value[j]).abs() < 1e-6,
                "{a:?} {b:?}"
            );
            assert!((a.col_dual[j] - b.col_dual[j]).abs() < 1e-6, "{a:?} {b:?}");
        }
    }
}
