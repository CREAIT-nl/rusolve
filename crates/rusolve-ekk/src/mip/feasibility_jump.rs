//! HiGHS source: mip/HighsFeasibilityJump.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::mip::feasibilityjump as fj;
use crate::prelude::*;

impl HighsMipSolver {
    /// rusolve mip_first_inc: called when the pre-root FJ left no incumbent.
    /// Native FJ (Luteberget and Sartor 2023, own budget nnz * 2^10), then the FPR and Local-MIP
    /// presolve chain that is otherwise off by default.
    pub fn first_incumbent_fallback(&mut self) {
        if !self.options_mip_.mip_fj_native {
            let model: &HighsLp = mip_model!(self);
            let r = crate::mip::fj_native::feasibility_jump_native(
                model,
                None,
                self.mipdata_.feastol,
                self.mipdata_.epsilon,
                self.options_mip_.random_seed as u64,
                None,
            );
            if let Some(r) = r {
                let _ = self.try_solution(&r.x, kSolutionSourceFeasibilityJump);
            }
        }
        if self.mipdata_.upper_limit == kHighsInf
            && !(self.options_mip_.mip_heuristic_run_fpr
                || self.options_mip_.mip_heuristic_run_local_mip)
        {
            let (fpr, lm) = (
                self.options_mip_.mip_heuristic_run_fpr,
                self.options_mip_.mip_heuristic_run_local_mip,
            );
            self.options_mip_.mip_heuristic_run_fpr = true;
            self.options_mip_.mip_heuristic_run_local_mip = true;
            let _ = crate::mip::fpr_local_mip::driver::s5_run_presolve(self);
            self.options_mip_.mip_heuristic_run_fpr = fpr;
            self.options_mip_.mip_heuristic_run_local_mip = lm;
        }
    }

    /// HiGHS: HighsFeasibilityJump.cpp::HighsMipSolverData::feasibilityJump (lines 12-139)
    /// C++: `HighsModelStatus feasibilityJump()`
    pub fn feasibility_jump(&mut self) -> HighsModelStatus {
        // This is the (presolved) model being solved
        let feastol = self.mipdata_.feastol;
        let epsilon = self.mipdata_.epsilon;
        let model: &HighsLp = mip_model!(self);
        let sense_multiplier = model.sense_ as i32 as f64;

        // The C++ `#ifdef HIGHSINT64` branch (skip FJ) is not taken: HighsInt is 32-bit here.

        if self.options_mip_.mip_fj_native {
            // rusolve mip_fj_native: the native engine's FJ (`fj_native.rs`).
            let start: Option<Vec<f64>> =
                (!self.mipdata_.incumbent.is_empty()).then(|| self.mipdata_.incumbent.clone());
            let r = crate::mip::fj_native::feasibility_jump_native(
                model,
                start.as_deref(),
                feastol,
                epsilon,
                self.options_mip_.random_seed as u64,
                None,
            );
            if let Some(r) = r {
                let _ = self.try_solution(&r.x, kSolutionSourceFeasibilityJump);
            }
            return HighsModelStatus::kNotset;
        }
        let mut found_integer_feasible_solution = false;
        let mut col_value: Vec<f64> = vec![0.0; model.num_col_ as usize];

        let use_incumbent = !self.mipdata_.incumbent.is_empty();

        // Configure Feasibility Jump and pass it the problem
        let mut solver = fj::FeasibilityJumpSolver::new(
            &self.options_mip_.log_options,
            self.options_mip_.random_seed,
            epsilon,
            feastol,
            1.0,
        );
        // HiGHS 1.16 #3310 / #3269 (h16_fj): option plumbing into the solver and its JumpMove (no
        // HighsOptions in feasibilityjump.hh); the default `false` runs the 1.15.1 statements.
        solver.h16_fj = self.options_mip_.h16_fj;
        solver.jump_move.h16_fj = self.options_mip_.h16_fj;

        for col in 0..model.num_col_ {
            let cu = col as usize;
            let mut lower = model.col_lower_[cu];
            let mut upper = model.col_upper_[cu];

            debug_assert!(
                model.integrality_[cu] == HighsVarType::kContinuous
                    || model.integrality_[cu] == HighsVarType::kInteger
                    || model.integrality_[cu] == HighsVarType::kImplicitInteger
            );
            let fj_var_type: fj::VarType;
            if model.integrality_[cu] == HighsVarType::kContinuous {
                fj_var_type = fj::Continuous;
            } else {
                fj_var_type = fj::Integer;
                lower = (lower - feastol).ceil();
                upper = (upper + feastol).floor();
            }

            let legal_bounds = lower <= upper
                && lower < kHighsInf
                && upper > -kHighsInf
                && !lower.is_nan()
                && !upper.is_nan();
            if !legal_bounds {
                crate::highs_log_user!(
                    &self.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "HighsMipSolverData::feasibilityJump() has detected infeasible/illegal bounds [{}, {}] for column {}: MIP is infeasible\n",
                    lower,
                    upper,
                    col
                );
                debug_assert!(legal_bounds);
                return HighsModelStatus::kInfeasible;
            }
            solver.add_var(
                fj_var_type,
                lower,
                upper,
                sense_multiplier * model.col_cost_[cu],
            );

            let mut initial_assignment = 0.0f64;
            if use_incumbent && self.mipdata_.incumbent[cu].is_finite() {
                initial_assignment = cpp_max(lower, cpp_min(upper, self.mipdata_.incumbent[cu]));
            } else if lower.is_finite() {
                initial_assignment = lower;
            } else if upper.is_finite() {
                initial_assignment = upper;
            }
            col_value[cu] = initial_assignment;
        }

        let mut a_matrix = HighsSparseMatrix::default();
        a_matrix.create_rowwise(&model.a_matrix_);

        for row in 0..model.num_row_ {
            let ru = row as usize;
            let has_finite_lower = model.row_lower_[ru].is_finite();
            let has_finite_upper = model.row_upper_[ru].is_finite();
            if has_finite_lower || has_finite_upper {
                let row_num_nz = a_matrix.start_[ru + 1] - a_matrix.start_[ru];
                let s = a_matrix.start_[ru] as usize;
                let e = a_matrix.start_[ru + 1] as usize;
                if has_finite_lower {
                    solver.add_constraint(
                        fj::Gte,
                        model.row_lower_[ru],
                        row_num_nz,
                        &mut a_matrix.index_[s..e],
                        &mut a_matrix.value_[s..e],
                        0,
                    );
                }
                if has_finite_upper {
                    solver.add_constraint(
                        fj::Lte,
                        model.row_upper_[ru],
                        row_num_nz,
                        &mut a_matrix.index_[s..e],
                        &mut a_matrix.value_[s..e],
                        0,
                    );
                }
            }
        }

        let nnz = a_matrix.num_nz();
        let max_total_effort: usize = (nnz as usize) << 10;
        let max_effort_since_last_improvement: usize = (nnz as usize) << 8;

        // C++ `col_value` is reassigned by the callback; `solve` gets a copy of the initial
        // values (solve only reads them in `init`, before the first callback).
        let mut initial_values: Vec<f64> = col_value.clone();
        let ls_t0 = std::time::Instant::now();
        // PORT NOTE: the C++ callback also computes `objective_function_value` (never read).
        let _ = solver.solve(&mut initial_values, &mut |status| {
            if let Some(solution) = status.solution {
                found_integer_feasible_solution = true;
                col_value = solution[..status.num_vars as usize].to_vec();
            }
            if status.effort_since_last_improvement > max_effort_since_last_improvement
                || status.total_effort > max_total_effort
            {
                fj::Terminate
            } else {
                fj::Continue
            }
        });
        let fj_stats = &crate::mip::fpr_local_mip::stats::FJ;
        fj_stats.add_moves(solver.ls_evals, solver.ls_steps);
        fj_stats.finish_call(ls_t0, solver.total_effort as u64);
        crate::work::add_fj_units(solver.total_effort as u64); // reporting-only

        if found_integer_feasible_solution {
            // Initial assignments that violate integrality or column bounds can lead to
            // infeasible results. Even if those initial assignments should not occur,
            // use trySolution rather than addIncumbent for an explicit check.
            let _ = self.try_solution(&col_value, kSolutionSourceFeasibilityJump);
        }
        // rusolve (mip_ls_int_lp_repair, after Rigal 2020 / Gardi et al. 2014). Polish the FJ
        // solution's continuous part, or repair FJ's infeasible end state, by an LP with
        // the integers fixed.
        if self.options_mip_.mip_ls_int_lp_repair {
            let point: &[f64] = if found_integer_feasible_solution {
                &col_value
            } else {
                &solver.problem.incumbent_assignment
            };
            let _ = crate::mip::ls_int_lp_repair::repair(
                self,
                point,
                kSolutionSourceFeasibilityJump,
                found_integer_feasible_solution,
            );
        }
        HighsModelStatus::kNotset
    }
}
