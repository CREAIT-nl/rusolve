//! HiGHS source: presolve/HPresolve.cpp (HiGHS latest, #3312, merge 7b890c255, lines 9842-9898)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiGHS #3312: the presolve clique-table updates for a
//! fixed or a substituted column. Called only from `mip_presolve_clique_table` gates.
//! Read `presolve/hpresolve/mod.rs` first (ownership, MIP mode).

use super::*;
use crate::mip::clique_table::HighsCliqueTable_CliqueVar;
use crate::prelude::*;

impl HPresolve {
    /// `mipsolver != nullptr && mipsolver->mipdata_->cliquesExtracted`
    fn x1b_cliques_extracted(&self) -> bool {
        match &self.mipsolver {
            Some(ms) => ms.mipdata_.cliques_extracted,
            None => false,
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::updateCliqueTableFixedCol (lines 9842-9864)
    /// C++: `Result updateCliqueTableFixedCol(HighsInt col, double val)`
    pub fn update_clique_table_fixed_col(&mut self, col: HighsInt, val: f64) -> HPresolveResult {
        if !self.x1b_cliques_extracted()
            || self.model.integrality_[col as usize] == HighsVarType::kContinuous
            || (val != 0.0 && val != 1.0)
        {
            return HPresolveResult::kOk;
        }
        let mut implied_fixings: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
        let ms = self.mipsolver.as_mut().unwrap();
        if !self.presolve_clique_table.fix_col(
            &mut ms.mipdata_.cliquetable,
            col,
            val != 0.0,
            &mut implied_fixings,
        ) {
            return HPresolveResult::kPrimalInfeasible;
        }

        for fixing in &implied_fixings {
            if self.col_deleted[fixing.col() as usize] != 0 {
                continue;
            }
            let fixval = fixing.val() as f64;
            if fixval < self.model.col_lower_[fixing.col() as usize] - self.primal_feastol
                || fixval > self.model.col_upper_[fixing.col() as usize] + self.primal_feastol
            {
                return HPresolveResult::kPrimalInfeasible;
            }
            hpresolve_checked_call!(self.change_col_bounds(
                fixing.col() as HighsInt,
                fixval,
                fixval
            ));
        }
        HPresolveResult::kOk
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::updateCliqueTableSubstituteCol (lines 9866-9898)
    /// C++: `Result updateCliqueTableSubstituteCol(HighsInt substCol, HighsInt stayCol, double offset, double scale)`
    pub fn update_clique_table_substitute_col(
        &mut self,
        subst_col: HighsInt,
        stay_col: HighsInt,
        offset: f64,
        scale: f64,
    ) -> HPresolveResult {
        if !self.x1b_cliques_extracted() {
            return HPresolveResult::kOk;
        }
        let s = subst_col as usize;
        let t = stay_col as usize;
        let is_binary = self.model.integrality_[s] != HighsVarType::kContinuous
            && self.model.integrality_[t] != HighsVarType::kContinuous
            && self.model.col_lower_[s] == 0.0
            && self.model.col_lower_[t] == 0.0
            && self.model.col_upper_[s] == 1.0
            && self.model.col_upper_[t] == 1.0
            && ((offset == 0.0 && scale == 1.0) || (offset == 1.0 && scale == -1.0));

        if !is_binary {
            return HPresolveResult::kOk;
        }

        let replacement =
            HighsCliqueTable_CliqueVar::new(stay_col, if scale == 1.0 { 1 } else { 0 });
        let mut implied_fixings: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
        let ms = self.mipsolver.as_mut().unwrap();
        if !self.presolve_clique_table.substitute_col(
            &mut ms.mipdata_.cliquetable,
            subst_col,
            replacement,
            &mut implied_fixings,
        ) {
            return HPresolveResult::kPrimalInfeasible;
        }

        for v in &implied_fixings {
            let col = v.col() as HighsInt;
            let val = v.val() as f64;
            if self.col_deleted[col as usize] != 0 {
                continue;
            }
            if val < self.model.col_lower_[col as usize] - self.primal_feastol
                || val > self.model.col_upper_[col as usize] + self.primal_feastol
            {
                return HPresolveResult::kPrimalInfeasible;
            }
            hpresolve_checked_call!(self.change_col_bounds(col, val, val));
        }
        HPresolveResult::kOk
    }
}

#[cfg(test)]
mod x1b_clqp_d_tests {
    use super::*;
    use crate::mip::clique_table::HighsCliqueTable;

    /// LP with a doubleton equation (r0), a singleton row (r2), a fixed column (x4) and a
    /// column (x3) that dual fixing can move to a bound; solved with the option `flag`.
    fn solve_lp(flag: bool) -> (Vec<f64>, Vec<f64>, HighsInt, f64) {
        let mut lp = HighsLp::default();
        lp.num_col_ = 5;
        lp.num_row_ = 4;
        lp.col_cost_ = vec![1.0, 2.0, 3.0, -1.0, 0.5];
        lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0, 1.0];
        lp.col_upper_ = vec![10.0, 10.0, 10.0, 10.0, 1.0];
        lp.row_lower_ = vec![2.0, 1.0, -kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![2.0, kHighsInf, 3.0, 5.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        // x0: r0 r3; x1: r0 r1; x2: r1 r2 r3; x3: r1 r3; x4: r1
        lp.a_matrix_.start_ = vec![0, 2, 4, 7, 9, 10];
        lp.a_matrix_.index_ = vec![0, 3, 0, 1, 1, 2, 3, 1, 3, 1];
        lp.a_matrix_.value_ = vec![1.0; 10];
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_bool("mip_presolve_clique_table", flag);
        h.pass_model(lp);
        h.run();
        assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
        let sol = h.get_solution().clone();
        let iters = h.get_info().simplex_iteration_count;
        let obj = h.get_info().objective_function_value;
        (sol.col_value, sol.row_dual, iters, obj)
    }

    /// LP presolve has no MIP solver, so every gated site runs the 1.15.1 statements:
    /// the option on and off give the same solution, duals and iteration count.
    #[test]
    fn x1b_clqp_d_lp_unchanged() {
        let off = solve_lp(false);
        let on = solve_lp(true);
        assert_eq!(off, on);
    }

    /// Clique `x0 + x1 + x2 <= 1` (all binary) in a MIP-mode HPresolve with the cliques
    /// extracted; returns the upper bounds after `changeColLower(x0, 1)`.
    fn fix_x0(flag: bool) -> Vec<f64> {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 1;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3];
        lp.a_matrix_.index_ = vec![0, 0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0];
        lp.col_cost_ = vec![-1.0, -1.0, -1.0];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.integrality_ = vec![HighsVarType::kInteger; 3];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![1.0];

        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_clique_table = flag;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        let mut table = HighsCliqueTable::new(3);
        let clq = [
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(1, 1),
            HighsCliqueTable_CliqueVar::new(2, 1),
        ];
        table.do_add_clique(&clq, false, kHighsIInf);
        *ms.mipdata_.cliquetable = table;
        ms.mipdata_.cliques_extracted = true;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), options.presolve_reduction_limit));
        let ms = presolve.mipsolver.as_ref().unwrap();
        presolve
            .presolve_clique_table
            .rebuild(&ms.mipdata_.cliquetable);
        assert_eq!(presolve.change_col_lower(0, 1.0), HPresolveResult::kOk);
        presolve.model.col_upper_.clone()
    }

    /// On (HPresolve.cpp 2733-2737 at 7b890c255): fixing x0 = 1 fixes its clique partners
    /// to 0. Off (1.15.1): only x0's bound changes.
    #[test]
    fn x1b_clqp_d_fixed_col_updates_cliques() {
        assert_eq!(fix_x0(true), vec![1.0, 0.0, 0.0]);
        assert_eq!(fix_x0(false), vec![1.0, 1.0, 1.0]);
    }
}
