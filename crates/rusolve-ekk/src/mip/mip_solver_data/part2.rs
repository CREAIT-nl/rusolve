//! HiGHS source: mip/HighsMipSolverData.cpp (lines 831-1841)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::mip::scip_heuristics::{S6Origin, S6State};
use crate::prelude::*;
use crate::util::c_printf::{c_fmt_f, c_fmt_f_width, c_fmt_g};

/// SCIP 11.0.0 (32dd425): primal.c::primalAddSol (lines 859-1100) pool insert, D4
fn s6_pool_insert_incumbent(
    s6: &mut S6State,
    obj: f64,
    x: Vec<f64>,
    solution_source: i32,
    nodenum: i64,
) {
    let origin = if s6.cur_origin == S6Origin::Highs(kSolutionSourceNone) {
        S6Origin::Highs(solution_source)
    } else {
        s6.cur_origin
    };
    if let Some((pos, _)) = s6
        .pool
        .insert(obj, x, origin, nodenum, &mut s6.next_sol_index)
    {
        if origin == S6Origin::Crossover && pos == 0 {
            s6.xover_nbestsolsfound += 1;
        }
    }
}

impl HighsMipSolver {
    /// HiGHS 1.16 (#3523): HighsMipSolverData.cpp::checkAddSolution (lines 837-931)
    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::runSetup (lines 851-901)
    /// C++: `void checkAddSolution()`
    ///
    /// Extracted from `runSetup` (1.15.1 lines 851-901), verbatim on the off path. Called
    /// unconditionally from `run_setup` (and, on the on path, from `cleanup_solve`).
    pub fn check_add_solution(&mut self) {
        // HiGHS 1.16 #3523 (h16_mip_start): read the option once; off = the 1.15.1 block below
        let h16_mip_start = self.options_mip_.h16_mip_start;
        if self.solution_objective_ != kHighsInf {
            // Assigning new incumbent
            self.mipdata_.incumbent = self
                .mipdata_
                .post_solve_stack
                .get_reduced_primal_solution(&self.solution_);
            // return the objective value in the transformed space
            let mut solobj = self.solution_objective_ * (self.orig_model_.sense_ as i32) as f64
                - mip_model!(self).offset_;
            let feasible = self.bound_violation_ <= self.options_mip_.mip_feasibility_tolerance
                && self.integrality_violation_ <= self.options_mip_.mip_feasibility_tolerance
                && self.row_violation_ <= self.options_mip_.mip_feasibility_tolerance;
            if self.mipdata_.num_restarts == 0 {
                if h16_mip_start {
                    // HiGHS 1.16 #3523 (h16_mip_start): HighsMipSolverData.cpp checkAddSolution lines 854-897 (1.16.0)
                    if feasible {
                        crate::highs_log_user!(
                            &self.options_mip_.log_options,
                            HighsLogType::kInfo,
                            "\nMIP start solution is feasible: objective value is {}\n",
                            c_fmt_g(self.solution_objective_, 12)
                        );
                    } else {
                        // The initial solution shouldn't be infeasible in the
                        // original space, since it's tested externally
                        crate::highs_log_user!(
                            &self.options_mip_.log_options,
                            HighsLogType::kWarning,
                            concat!(
                                "\nMIP start solution is infeasible: objective value is ",
                                "{} with violations:",
                                " variable bound = {}; constraint bound = {}; ",
                                "integrality = {}\n"
                            ),
                            c_fmt_g(self.solution_objective_, 12),
                            c_fmt_g(self.bound_violation_, 4),
                            c_fmt_g(self.row_violation_, 4),
                            c_fmt_g(self.integrality_violation_, 4)
                        );
                    }
                    // Since the supplied solution (mipsolver.solution_) may not
                    // correspond to postsolve applied to incumbent (#3510), compute
                    // the objective locally using incumbent, saving the value of
                    // mipsolver.solution_objective_ to identify changes
                    let previous_mipsolver_solution_objective = self.solution_objective_;
                    let mipsolver_solution_objective = f64::from(
                        mip_model!(self).objective_c_double_value(&self.mipdata_.incumbent),
                    );
                    self.solution_objective_ =
                        (self.orig_model_.sense_ as i32) as f64 * mipsolver_solution_objective;
                    let has_start_solution_objective_difference =
                        (self.solution_objective_ - previous_mipsolver_solution_objective).abs()
                            / (1.0 + previous_mipsolver_solution_objective.abs())
                            > self.options_mip_.mip_feasibility_tolerance;
                    // This affects solobj - which is use used to determine bounds
                    // etc (without the objective offset)
                    solobj = mipsolver_solution_objective - mip_model!(self).offset_;
                    if has_start_solution_objective_difference {
                        let reduced = self.mipdata_.post_solve_stack.num_reductions() > 0;
                        crate::highs_log_user!(
                            &self.options_mip_.log_options,
                            HighsLogType::kInfo,
                            "MIP start solution is {}{}: objective value is {}\n",
                            if feasible { "feasible" } else { "infeasible" },
                            if reduced { " in the reduced space" } else { "" },
                            c_fmt_g(self.solution_objective_, 12)
                        );
                    }
                } else {
                    crate::highs_log_user!(
                        &self.options_mip_.log_options,
                        HighsLogType::kInfo,
                        "\nMIP start solution is {}, objective value is {}\n",
                        if feasible { "feasible" } else { "infeasible" },
                        c_fmt_g(self.solution_objective_, 12)
                    );
                }
            }
            if feasible && solobj < self.mipdata_.upper_bound {
                let prev_upper_bound = self.mipdata_.upper_bound;

                self.mipdata_.upper_bound = solobj;

                let bound_change = self.mipdata_.upper_bound != prev_upper_bound;
                if !self.submip && bound_change {
                    let lower_bound = self.mipdata_.lower_bound;
                    let upper_bound = self.mipdata_.upper_bound;
                    self.update_primal_dual_integral(
                        lower_bound,
                        lower_bound,
                        prev_upper_bound,
                        upper_bound,
                        true,
                        true,
                    );
                }

                let new_upper_limit = self.compute_new_upper_limit(solobj, 0.0, 0.0);

                self.save_report_mip_solution(new_upper_limit);
                if new_upper_limit < self.mipdata_.upper_limit {
                    self.mipdata_.upper_limit = new_upper_limit;
                    self.mipdata_.optimality_limit = self.compute_new_upper_limit(
                        solobj,
                        self.options_mip_.mip_abs_gap,
                        self.options_mip_.mip_rel_gap,
                    );
                    let optimality_limit = self.mipdata_.optimality_limit;
                    self.mipdata_
                        .nodequeue
                        .set_optimality_limit(optimality_limit);
                }
            }
            if !self.submip
                && feasible
                && self.callback_.user_callback.is_some()
                && self.callback_.active[kCallbackMipSolution as usize]
            {
                debug_assert!(!self.submip);
                self.callback_.clear_highs_callback_output();
                self.callback_.data_out.mip_solution = self.solution_.clone();
                let interrupt = self.interrupt_from_callback_with_data(
                    kCallbackMipSolution,
                    self.solution_objective_,
                    "Feasible solution",
                );
                debug_assert!(!interrupt);
            }
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::runSetup (lines 831-1157)
    /// C++: `void runSetup()`
    ///
    /// PORT NOTE: the C++ `const HighsLp& model = *mipsolver.model_` alias is not held across
    /// mutation of `presolvedModel`; `mip_model!(self)` is re-borrowed at each use, which reads
    /// the same object the alias would at that point. The HIGHS_DEBUGSOL block is compiled out
    /// (the oracle build does not define it).
    #[allow(clippy::needless_range_loop)]
    pub fn run_setup(&mut self) {
        // Indicate that the first LP has not been solved
        self.get_lp().set_solved_first_lp(false);

        self.mipdata_.last_disptime = -kHighsInf;
        self.mipdata_.disptime = 0.0;

        // Transform the reference of the objective limit and lower/upper
        // bounds from the original model to the current model, undoing the
        // transformation done before restart so that the offset change due
        // to presolve is incorporated. Bound changes are transitory, so no
        // real gap change, and no update to P-D integral is necessary
        self.mipdata_.upper_limit -= mip_model!(self).offset_;
        self.mipdata_.optimality_limit -= mip_model!(self).offset_;

        self.mipdata_.lower_bound -= mip_model!(self).offset_;
        self.mipdata_.upper_bound -= mip_model!(self).offset_;

        self.check_add_solution();

        // S6 (mip_heuristic_run_oneopt|crossover): SCIP heur_oneopt.c lines 454-468, heur_crossover.c lines 835-864 at 32dd425
        let s6_on = self.options_mip_.mip_alns
            || self.options_mip_.mip_proximity
            || self.options_mip_.mip_improve_early
            || self.options_mip_.mip_primal_stall
            || self.options_mip_.mip_local_branching
            || self.options_mip_.mip_heuristic_run_oneopt
            || self.options_mip_.mip_heuristic_run_crossover;
        if s6_on {
            self.mipdata_.s6 = Some(Box::new(S6State::new(self)));
            if self.options_mip_.mip_alns && !self.submip {
                self.mipdata_.alns =
                    Some(Box::new(crate::mip::scip_heuristics::alns::AlnsState::new(
                        self.options_mip_.random_seed,
                    )));
            }
            if (self.options_mip_.mip_proximity
                || self.options_mip_.mip_improve_early
                || self.options_mip_.mip_primal_stall)
                && !self.submip
            {
                self.mipdata_.proximity = Some(Box::default());
            }
            if self.options_mip_.mip_local_branching && !self.submip {
                self.mipdata_.local_branching = Some(Box::default());
            }
            // insert the reduced incumbent assigned above if it was feasible. solobj and
            // feasible are recomputed with the same expressions as the HiGHS block; nothing in
            // that block writes the fields they read.
            if self.solution_objective_ != kHighsInf {
                let solobj = self.solution_objective_ * (self.orig_model_.sense_ as i32) as f64
                    - mip_model!(self).offset_;
                let feasible = self.bound_violation_ <= self.options_mip_.mip_feasibility_tolerance
                    && self.integrality_violation_ <= self.options_mip_.mip_feasibility_tolerance
                    && self.row_violation_ <= self.options_mip_.mip_feasibility_tolerance;
                if feasible {
                    let nodenum = self.mipdata_.num_nodes;
                    let x = self.mipdata_.incumbent.clone();
                    if let Some(s6) = self.mipdata_.s6.as_mut() {
                        s6_pool_insert_incumbent(
                            s6,
                            solobj,
                            x,
                            kSolutionSourceHighsSolution,
                            nodenum,
                        );
                    }
                }
            }
        }

        if self.num_col() == 0 {
            self.add_incumbent(&[], 0.0, kSolutionSourceEmptyMip, true, false);
        }

        *self.mipdata_.redcostfixing = HighsRedcostFixing::default();
        let pseudocost = HighsPseudocost::new(self);
        *self.get_pseudo_cost() = pseudocost;
        let num_col = self.num_col();
        self.mipdata_.nodequeue.set_num_col(num_col);
        let optimality_limit = self.mipdata_.optimality_limit;
        self.mipdata_
            .nodequeue
            .set_optimality_limit(optimality_limit);

        self.mipdata_.continuous_cols.clear();
        self.mipdata_.integer_cols.clear();
        self.mipdata_.implint_cols.clear();
        self.mipdata_.integral_cols.clear();

        self.mipdata_.row_matrix_set = false;
        if !self.mipdata_.row_matrix_set {
            self.mipdata_.row_matrix_set = true;
            highs_sparse_transpose(
                mip_model!(self).num_row_,
                mip_model!(self).num_col_,
                &mip_model!(self).a_matrix_.start_,
                &mip_model!(self).a_matrix_.index_,
                &mip_model!(self).a_matrix_.value_,
                &mut self.mipdata_.a_rstart_,
                &mut self.mipdata_.a_rindex_,
                &mut self.mipdata_.a_rvalue_,
            );
            // (re-)initialize number of uplocks and downlocks
            let model_num_col = mip_model!(self).num_col_ as usize;
            self.mipdata_.uplocks.clear();
            self.mipdata_.uplocks.resize(model_num_col, 0);
            self.mipdata_.downlocks.clear();
            self.mipdata_.downlocks.resize(model_num_col, 0);
            let model = mip_model!(self);
            for i in 0..model_num_col {
                let start = model.a_matrix_.start_[i] as usize;
                let end = model.a_matrix_.start_[i + 1] as usize;
                for j in start..end {
                    let row = model.a_matrix_.index_[j] as usize;

                    if model.row_lower_[row] != -kHighsInf {
                        if model.a_matrix_.value_[j] < 0.0 {
                            self.mipdata_.uplocks[i] += 1;
                        } else {
                            self.mipdata_.downlocks[i] += 1;
                        }
                    }
                    if model.row_upper_[row] != kHighsInf {
                        if model.a_matrix_.value_[j] < 0.0 {
                            self.mipdata_.downlocks[i] += 1;
                        } else {
                            self.mipdata_.uplocks[i] += 1;
                        }
                    }
                }
            }
        }

        let num_row = self.num_row() as usize;
        self.mipdata_.rowintegral.resize(num_row, 0);

        // compute the maximal absolute coefficients to filter propagation
        self.mipdata_.max_abs_row_coef.resize(num_row, 0.0);
        for i in 0..num_row {
            let mut maxabsval = 0.0_f64;

            let start = self.mipdata_.a_rstart_[i] as usize;
            let end = self.mipdata_.a_rstart_[i + 1] as usize;
            let mut integral = true;
            for j in start..end {
                integral = integral
                    && self.is_col_integral(self.mipdata_.a_rindex_[j])
                    && fractionality(self.mipdata_.a_rvalue_[j], None) <= self.mipdata_.epsilon;

                maxabsval = cpp_max(maxabsval, self.mipdata_.a_rvalue_[j].abs());
            }

            if integral {
                if self.mipdata_.presolved_model.row_lower_[i] != -kHighsInf {
                    self.mipdata_.presolved_model.row_lower_[i] =
                        (self.mipdata_.presolved_model.row_lower_[i] - self.mipdata_.feastol)
                            .ceil();
                }

                if self.mipdata_.presolved_model.row_upper_[i] != kHighsInf {
                    self.mipdata_.presolved_model.row_upper_[i] =
                        (self.mipdata_.presolved_model.row_upper_[i] + self.mipdata_.feastol)
                            .floor();
                }
            }

            self.mipdata_.rowintegral[i] = integral as u8;
            self.mipdata_.max_abs_row_coef[i] = maxabsval;
        }

        // compute row activities and propagate all rows once
        // The clique partition reads the model HighsObjectiveFunction captured at
        // construction (`mipsolver.model_` then): the original model for the object built
        // by the HighsMipSolverData ctor, the presolved one after HPresolve::shrinkProblem
        // rebuilt it (HPresolve.cpp:1016).
        let objective_model = match self.mipdata_.objective_function.model {
            HighsMipModel::Orig => &self.orig_model_,
            HighsMipModel::Presolved => &self.mipdata_.presolved_model,
        };
        self.mipdata_.objective_function.setup_clique_partition(
            objective_model,
            &self.mipdata_.domains[0],
            &mut self.mipdata_.cliquetable,
        );
        detach!(self.mipdata_.domains[0], |dom| dom
            .setup_objective_propagation(self));
        detach!(self.mipdata_.domains[0], |dom| dom
            .compute_row_activities(self));
        detach!(self.mipdata_.domains[0], |dom| dom.propagate(self));
        if self.mipdata_.domains[0].infeasible() {
            self.modelstatus_ = HighsModelStatus::kInfeasible;

            self.update_lower_bound(kHighsInf, true, true);

            self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
            return;
        }

        if mip_model!(self).num_col_ == 0 {
            self.modelstatus_ = HighsModelStatus::kOptimal;
            return;
        }

        if self.check_limits(0) {
            return;
        }
        // extract cliques if they have not been extracted before

        // HiGHS quirk: C++ iterates the live changed-cols vector; cleanupVarbounds never
        // changes it (that would invalidate the range-for), so a snapshot is equivalent.
        let changed_cols: Vec<HighsInt> = self.mipdata_.domains[0].get_changed_cols().to_vec();
        for col in changed_cols {
            detach!(self.mipdata_.implications, |imp| imp
                .cleanup_varbounds(self, col));
        }
        self.mipdata_.domains[0].clear_changed_cols();

        self.get_lp()
            .get_lp_solver()
            .set_option_value_str("presolve", kHighsOffString);

        self.check_obj_integrality();
        self.mipdata_.rootlpsol.clear();
        self.mipdata_.firstlpsol.clear();
        let mut num_binary: HighsInt = 0;
        let mut num_domain_fixed: HighsInt = 0;
        self.mipdata_.max_tree_size_log2 = 0;
        let model_col_count = self.num_col();
        for i in 0..model_col_count {
            match self.variable_type_2(i) {
                HighsVarType::kContinuous => {
                    if self.mipdata_.domains[0].is_fixed(i) {
                        num_domain_fixed += 1;
                        continue;
                    }
                    self.mipdata_.continuous_cols.push(i);
                }
                HighsVarType::kImplicitInteger => {
                    if self.mipdata_.domains[0].is_fixed(i) {
                        num_domain_fixed += 1;
                        continue;
                    }
                    self.mipdata_.implint_cols.push(i);
                    self.mipdata_.integral_cols.push(i);
                }
                HighsVarType::kInteger => {
                    if self.mipdata_.domains[0].is_fixed(i) {
                        num_domain_fixed += 1;
                        if fractionality(self.mipdata_.domains[0].col_lower_[i as usize], None)
                            > self.mipdata_.feastol
                        {
                            // integer variable is fixed to a fractional value -> infeasible
                            self.modelstatus_ = HighsModelStatus::kInfeasible;

                            self.update_lower_bound(kHighsInf, true, true);

                            self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                            return;
                        }
                        continue;
                    }
                    self.mipdata_.integer_cols.push(i);
                    self.mipdata_.integral_cols.push(i);
                    let col_lower = mip_model!(self).col_lower_[i as usize];
                    let col_upper = mip_model!(self).col_upper_[i as usize];
                    self.mipdata_.max_tree_size_log2 +=
                        cpp_min(1024.0, 1.0 + col_upper - col_lower).log2().ceil() as HighsInt;
                    // NB Since this is for counting the number of times the condition is true,
                    // a bitwise operator avoids a conditional branch (C++ `&` on the two
                    // comparisons), so it is kept as bitwise here.
                    num_binary += ((col_lower == 0.0) & (col_upper == 1.0)) as HighsInt;
                }
                HighsVarType::kSemiContinuous | HighsVarType::kSemiInteger => {
                    crate::highs_log_user!(
                        &self.options_mip_.log_options,
                        HighsLogType::kError,
                        "Semicontinuous or semiinteger variables should have been reformulated away before HighsMipSolverData::runSetup() is called."
                    );
                    panic!("Unexpected variable type");
                }
            }
        }

        self.basis_transfer();

        self.mipdata_.numintegercols = self.mipdata_.integer_cols.len() as HighsInt;
        self.mipdata_.detect_symmetries = self.mipdata_.detect_symmetries && num_binary > 0;
        self.mipdata_.num_clique_entries_after_presolve =
            self.mipdata_.cliquetable.get_num_entries();
        let num_col = self.num_col();
        let num_general_integer = self.mipdata_.numintegercols - num_binary;
        let num_implied_integer = self.mipdata_.implint_cols.len() as HighsInt;
        let num_continuous = self.mipdata_.continuous_cols.len() as HighsInt;
        debug_assert!(
            num_col
                == num_continuous
                    + num_binary
                    + num_general_integer
                    + num_implied_integer
                    + num_domain_fixed
        );
        let num_row = self.num_row();
        let num_nonzero = self.num_nonzero();
        if self.mipdata_.num_restarts == 0 {
            self.mipdata_.num_clique_entries_after_first_presolve =
                self.mipdata_.cliquetable.get_num_entries();
            // Hardware concurrency is machine-dependent and reads 0 when unknown, as
            // `std::thread::hardware_concurrency()` does.
            let num_threads: HighsInt = crate::parallel::highs_parallel::num_threads();
            let hardware_threads = std::thread::available_parallelism()
                .map(|n| n.get() as HighsInt)
                .unwrap_or(0);
            let max_workers = self.get_max_num_workers();
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "\nSolving MIP model with:\n   {} row{}\n   {} col{} ({} binary, {} integer, {} implied int., {} continuous, {} domain fixed)\n   {} nonzero{}\n   Thread count {} (of {} threads). Using {} max workers. Parallel search {}\n",
                num_row,
                if num_row == 1 { "" } else { "s" },
                num_col,
                if num_col == 1 { "" } else { "s" },
                num_binary,
                num_general_integer,
                num_implied_integer,
                num_continuous,
                num_domain_fixed,
                num_nonzero,
                if num_nonzero == 1 { "" } else { "s" },
                num_threads,
                hardware_threads,
                max_workers,
                if max_workers > 1 { "on" } else { "off" }
            );
        } else {
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "Model after restart has {} row{}, {} col{} ({} bin., {} int., {} impl., {} cont., {} dom.fix.), and {} nonzero{}\n",
                num_row,
                if num_row == 1 { "" } else { "s" },
                num_col,
                if num_col == 1 { "" } else { "s" },
                num_binary,
                num_general_integer,
                num_implied_integer,
                num_continuous,
                num_domain_fixed,
                num_nonzero,
                if num_nonzero == 1 { "" } else { "s" }
            );
        }

        detach!(self.mipdata_.heuristics, |h| h.setup_int_cols(self));

        // HIGHS_DEBUGSOL block (mip_debug_solution_file): HighsMipSolverData.cpp lines 1132-1145
        HighsDebugSol::run_setup_block(self);

        if self.mipdata_.upper_limit == kHighsInf {
            self.mipdata_.analytic_center_computed = false;
        }
        self.mipdata_.analytic_center_status = HighsModelStatus::kNotset;
        self.mipdata_.analytic_center.clear();

        self.mipdata_.symmetries.clear();
        // rusolve mip_symmetry_dynamic: clear dynamic symmetry data
        if self.options_mip_.mip_symmetry_dynamic {
            self.mipdata_.symmetry_dynamic = None;
        }

        if self.mipdata_.num_restarts != 0 {
            crate::highs_log_user!(&self.options_mip_.log_options, HighsLogType::kInfo, "\n");
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::transformNewIntegerFeasibleSolution (lines 1159-1308)
    /// HiGHS 1.16 (#3181): HighsMipSolverData.cpp::transformNewIntegerFeasibleSolution (lines 1288-1299)
    /// C++: `double transformNewIntegerFeasibleSolution(const std::vector<double> & sol, const bool possibly_store_as_new_incumbent = true)`
    /// C++ default arguments: `possibly_store_as_new_incumbent = true`
    pub fn transform_new_integer_feasible_solution(
        &mut self,
        sol: &[f64],
        possibly_store_as_new_incumbent: bool,
    ) -> f64 {
        let mut solution = HighsSolution::default();
        solution.col_value = sol.to_vec();
        solution.value_valid = true;
        // Perform primal postsolve to get the original column values
        self.mipdata_
            .post_solve_stack
            .undo_primal(&self.options_mip_, &mut solution, -1, false);
        // Determine the row values, as they aren't computed in primal
        // postsolve
        let return_status = crate::lp_data::presolve_support::calculate_row_values_quad(
            &self.orig_model_,
            &mut solution,
            -1,
        );
        if kAllowDeveloperAssert {
            debug_assert!(return_status == HighsStatus::kOk);
        }
        let mut allow_try_again = true;
        // `try_again:` label: the `continue` below re-runs the body from the top
        loop {
            // compute the objective value in the original space
            let mut bound_violation = 0.0_f64;
            let mut row_violation = 0.0_f64;
            let mut integrality_violation = 0.0_f64;
            let mut mipsolver_quad_objective_value = HighsCDouble::from(0.0);
            let feasible = self.solution_feasible(
                Some(&self.orig_model_),
                &solution.col_value,
                Some(&solution.row_value),
                &mut bound_violation,
                &mut row_violation,
                &mut integrality_violation,
                &mut mipsolver_quad_objective_value,
            );
            let mipsolver_objective_value = f64::from(mipsolver_quad_objective_value);
            if !feasible && allow_try_again {
                let mut fixed_model = self.orig_model_.clone();
                fixed_model.integrality_.clear();
                for i in 0..self.orig_model_.num_col_ as usize {
                    if self.orig_model_.integrality_[i] == HighsVarType::kInteger {
                        let solval = solution.col_value[i].round();
                        fixed_model.col_lower_[i] = cpp_max(fixed_model.col_lower_[i], solval);
                        fixed_model.col_upper_[i] = cpp_min(fixed_model.col_upper_[i], solval);
                    }
                }
                self.mipdata_.total_repair_lp += 1;
                let time_available =
                    cpp_max(self.options_mip_.time_limit - self.timer_.read(), 0.1);
                // Highs instantiation
                let mut tmp_solver = Highs::default();
                tmp_solver.set_profiling();
                // HiGHS quirk: the `debug_report` branch (log_dev_level / highs_analysis_level)
                // is dead (`const bool debug_report = false`) and is not ported.
                tmp_solver.set_option_value_bool("output_flag", false);
                // h16_sub_solver_flags: the C++ `Highs tmpSolver` runs 1.16.0 code
                if self.options_mip_.h16_sub_solver_flags {
                    tmp_solver
                        .options_
                        .copy_version_flags_from(&self.options_mip_);
                }
                tmp_solver.set_option_value_double("time_limit", time_available);
                // Set primal feasibility tolerance for LP solves according to
                // mip_feasibility_tolerance. Interestingly, dual feasibility
                // tolerance not set to smaller tolerance as in
                // HighsLpRelaxationconstructor.
                let mip_primal_feasibility_tolerance = self.options_mip_.mip_feasibility_tolerance;
                tmp_solver.set_option_value_double(
                    "primal_feasibility_tolerance",
                    mip_primal_feasibility_tolerance,
                );
                // check if only root presolve is allowed
                let use_presolve = !self.options_mip_.mip_root_presolve_only;
                let presolve = if use_presolve {
                    kHighsChooseString
                } else {
                    kHighsOffString
                };
                tmp_solver.set_option_value_str("presolve", presolve);
                tmp_solver.pass_model(fixed_model);
                // Until a good decision can be made on whether to use simplex,
                // HiPO or IPX to solve an LP without a basis, use simplex
                tmp_solver.set_option_value_str("solver", kSimplexString);
                tmp_solver.optimize_lp();
                // HiGHS 1.16 #3181 (h16_repair_lp): HighsMipSolverData.cpp lines 1288-1292 (1.16.0)
                let h16_repair_lp = self.options_mip_.h16_repair_lp;
                if h16_repair_lp
                    && use_presolve
                    && tmp_solver.get_model_status() == HighsModelStatus::kSolveError
                {
                    tmp_solver.set_option_value_str("presolve", kHighsOffString);
                    tmp_solver.optimize_lp();
                }
                self.mipdata_.total_repair_lp_iterations +=
                    i64::from(tmp_solver.get_info().simplex_iteration_count);
                if tmp_solver.get_info().primal_solution_status == kSolutionStatusFeasible {
                    self.mipdata_.total_repair_lp_feasible += 1;
                    solution = tmp_solver.get_solution().clone();
                    // HiGHS 1.16 #3181 (h16_repair_lp): HighsMipSolverData.cpp lines 1297-1299 (1.16.0)
                    if h16_repair_lp {
                        let _ = crate::lp_data::presolve_support::calculate_row_values_quad_vec(
                            &self.orig_model_,
                            &solution.col_value,
                            &mut solution.row_value,
                            -1,
                        );
                    }
                    allow_try_again = false;
                    continue; // goto try_again
                }
            }

            let transformed_solobj = f64::from(
                HighsCDouble::from(self.orig_model_.sense_ as i32 as f64)
                    * mipsolver_quad_objective_value
                    - mip_model!(self).offset_,
            );

            // Possible MIP solution callback
            if !self.submip
                && feasible
                && self.callback_.user_callback.is_some()
                && self.callback_.active[kCallbackMipSolution as usize]
            {
                self.callback_.clear_highs_callback_output();
                self.callback_.data_out.mip_solution = solution.col_value.clone();
                let interrupt = self.interrupt_from_callback_with_data(
                    kCallbackMipSolution,
                    mipsolver_objective_value,
                    "Feasible solution",
                );
                debug_assert!(!interrupt);
            }

            // Catch the case where the repaired solution now has worse objective
            // than the current stored solution
            if transformed_solobj >= self.mipdata_.upper_bound && !sol.is_empty() {
                return transformed_solobj;
            }

            if possibly_store_as_new_incumbent {
                // Store the solution as incumbent in the original space if there
                // is no solution or if it is feasible
                if feasible {
                    // store
                    self.row_violation_ = row_violation;
                    self.bound_violation_ = bound_violation;
                    self.integrality_violation_ = integrality_violation;
                    self.solution_ = std::mem::take(&mut solution.col_value);
                    self.solution_objective_ = mipsolver_objective_value;
                } else {
                    let current_feasible = self.solution_objective_ != kHighsInf
                        && self.bound_violation_ <= self.options_mip_.mip_feasibility_tolerance
                        && self.integrality_violation_
                            <= self.options_mip_.mip_feasibility_tolerance
                        && self.row_violation_ <= self.options_mip_.mip_feasibility_tolerance;
                    crate::highs_log_user!(
                        &self.options_mip_.log_options,
                        HighsLogType::kWarning,
                        "Solution with objective {} has untransformed violations: bound = {}; integrality = {}; row = {}\n",
                        c_fmt_g(mipsolver_objective_value, 6),
                        c_fmt_g(bound_violation, 4),
                        c_fmt_g(integrality_violation, 4),
                        c_fmt_g(row_violation, 4)
                    );
                    if !current_feasible {
                        // if the current incumbent is non existent or also not feasible we
                        // still store the new one
                        self.row_violation_ = row_violation;
                        self.bound_violation_ = bound_violation;
                        self.integrality_violation_ = integrality_violation;
                        self.solution_ = std::mem::take(&mut solution.col_value);
                        self.solution_objective_ = mipsolver_objective_value;
                    }

                    // return infinity so that it is not used for bounding
                    return kHighsInf;
                }
            }

            // return the objective value in the transformed space
            return transformed_solobj;
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::percentageInactiveIntegers (lines 1310-1315)
    /// C++: `double percentageInactiveIntegers() const`
    pub fn percentage_inactive_integers(&self) -> f64 {
        // HiGHS quirk: the size_t difference wraps when there are more substitutions than
        // integer columns (libc++ size_t arithmetic), so `wrapping_sub` is kept.
        let inactive = self
            .mipdata_
            .integer_cols
            .len()
            .wrapping_sub(self.mipdata_.cliquetable.get_substitutions_2().len());
        100.0 * (1.0 - (inactive as f64) / (self.mipdata_.numintegercols as f64))
    }

    /// rusolve `mip_cutoff_bigm`: cutoff-driven big-M tightening.
    ///
    /// Savelsbergh (1994) coefficient tightening, "Preprocessing and probing techniques for
    /// mixed integer programming problems", on the *global* bounds the objective cutoff has
    /// propagated. For a row `a x <= b` with a binary `z` and `a_z` positive, let
    /// `d = b - (maxact - a_z)`; when `d` is positive the row is slack at `z = 0` and the
    /// coefficient can drop to `a_z - d` with right-hand side `b - d`. For `a_z` negative the
    /// test is `d = b - (maxact + a_z)` and the coefficient rises to `a_z + d`. HiGHS presolve
    /// does this once, but the cutoff arrives after presolve and a big-M row only tightens when
    /// the start variables have a finite upper bound, which only the cutoff gives them. Returns
    /// true (restart, so presolve re-applies the tightening) when at least 1 % of the rows would
    /// lose at least 5 % of a big coefficient and the cutoff fell since the last such restart.
    pub fn cutoff_bigm_restart_due(&mut self) -> bool {
        if !self.options_mip_.mip_cutoff_bigm
            || self.submip
            || !self.options_mip_.mip_allow_restart
            || self.options_mip_.presolve == kHighsOffString
            || self.mipdata_.cutoff_bigm_restarts >= 8
        {
            return false;
        }
        let limit = self.mipdata_.upper_limit;
        if !limit.is_finite() {
            return false;
        }
        if let Some(last) = self.mipdata_.cutoff_bigm_last {
            if limit > last - 0.01 * last.abs().max(1.0) {
                return false;
            }
        }
        let num_row = self.num_row() as usize;
        let num_col = self.num_col() as usize;
        if num_row == 0 || self.mipdata_.a_rstart_.len() < num_row + 1 {
            return false;
        }
        let dom = &self.mipdata_.domains[0];
        if dom.col_lower_.len() < num_col {
            return false;
        }
        let model = mip_model!(self);
        let mut hits = 0usize;
        for i in 0..num_row {
            let (start, end) = (
                self.mipdata_.a_rstart_[i] as usize,
                self.mipdata_.a_rstart_[i + 1] as usize,
            );
            for (sign, bound) in [(1.0, model.row_upper_[i]), (-1.0, -model.row_lower_[i])] {
                if !bound.is_finite() {
                    continue;
                }
                let mut maxact = 0.0f64;
                let mut finite = true;
                for k in start..end {
                    let col = self.mipdata_.a_rindex_[k] as usize;
                    let a = sign * self.mipdata_.a_rvalue_[k];
                    let v = if a > 0.0 {
                        dom.col_upper_[col]
                    } else {
                        dom.col_lower_[col]
                    };
                    if !v.is_finite() {
                        finite = false;
                        break;
                    }
                    maxact += a * v;
                }
                if !finite {
                    continue;
                }
                let mut best = 0.0f64;
                for k in start..end {
                    let col = self.mipdata_.a_rindex_[k] as usize;
                    if !(self.is_col_integral(col as HighsInt)
                        && dom.col_lower_[col] == 0.0
                        && dom.col_upper_[col] == 1.0)
                    {
                        continue;
                    }
                    let a = sign * self.mipdata_.a_rvalue_[k];
                    let d = if a > 0.0 {
                        bound - (maxact - a)
                    } else {
                        bound - (maxact + a)
                    };
                    if d > 1e-6 && d < a.abs() {
                        best = best.max(d / a.abs());
                    }
                }
                if best >= 0.05 {
                    hits += 1;
                    break;
                }
            }
        }
        if (hits as f64) < 0.01 * num_row as f64 || hits == 0 {
            return false;
        }
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "\n{} of {} rows have big-M coefficients tightenable against the cutoff {}, restarting\n",
            hits,
            num_row,
            limit
        );
        self.mipdata_.cutoff_bigm_last = Some(limit);
        self.mipdata_.cutoff_bigm_restarts += 1;
        true
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::performRestart (lines 1317-1467)
    /// C++: `void performRestart()`
    pub fn perform_restart(&mut self) {
        // HiGHS #3218/#3250 (mip_presolve_row_append): cut rows are typed in the
        // postsolve stack and the root basis is expanded to its next row/col index
        // (HighsMipSolverData.cpp lines 1393 and 1409-1410, 1.16.0).
        let x1b_row_append = self.options_mip_.mip_presolve_row_append;
        let pscostinit = HighsPseudocostInitialization::new_2(
            &self.mipdata_.pseudocosts[0],
            self.options_mip_.mip_pscost_minreliable,
            &self.mipdata_.post_solve_stack,
            x1b_row_append,
        );

        self.pscostinit = Some(Box::new(pscostinit));
        self.mipdata_.num_restarts += 1;
        self.mipdata_.num_leaves_before_run = self.mipdata_.num_leaves;
        self.mipdata_.num_nodes_before_run = self.mipdata_.num_nodes;
        self.mipdata_.total_lp_iterations_before_run = self.mipdata_.total_lp_iterations;
        self.mipdata_.heuristic_lp_iterations_before_run = self.mipdata_.heuristic_lp_iterations;
        self.mipdata_.sepa_lp_iterations_before_run = self.mipdata_.sepa_lp_iterations;
        self.mipdata_.sb_lp_iterations_before_run = self.mipdata_.sb_lp_iterations;
        // rusolve mip_node_bigm: the restart takes the LP relaxation as the new model, so the
        // node-local big-M rows of the last node must go first, or they become global rows.
        for lp in self.mipdata_.lps.iter_mut() {
            lp.restore_local_mods();
        }
        let num_lp_rows = self.mipdata_.lps[0].get_lp().num_row_;
        let num_model_rows = self.num_row();
        let num_cuts = num_lp_rows - num_model_rows;
        if x1b_row_append {
            self.mipdata_
                .post_solve_stack
                .append_cuts_to_model_x1b(num_cuts);
        } else if num_cuts > 0 {
            self.mipdata_
                .post_solve_stack
                .append_cuts_to_model(num_cuts);
        }
        let integrality = std::mem::take(&mut self.mipdata_.presolved_model.integrality_);
        let offset = self.mipdata_.presolved_model.offset_;
        self.mipdata_.presolved_model = self.mipdata_.lps[0].get_lp().clone();
        self.mipdata_.presolved_model.offset_ = offset;
        self.mipdata_.presolved_model.integrality_ = integrality;
        // HIGHS_DEBUGSOL (mip_debug_solution_file): HighsMipSolverData.cpp line 1342
        self.mipdata_.debug_solution.suspend();

        // `root_basis` is the C++ stack local; `rootbasis_is_local` records the pointer
        // identity `rootbasis == &root_basis` checked at the end.
        let mut rootbasis_is_local = false;
        let basis = &self.mipdata_.firstrootbasis;
        if basis.valid {
            // if we have a basis after solving the root LP, we expand it to the
            // original space so that it can be used for constructing a starting basis
            // for the presolved model after the restart
            let mut root_basis = HighsBasis::default();
            let (basis_num_col, basis_num_row) = if x1b_row_append {
                (
                    self.mipdata_.post_solve_stack.get_next_col_index(),
                    self.mipdata_.post_solve_stack.get_next_row_index(),
                )
            } else {
                (
                    self.mipdata_.post_solve_stack.get_orig_num_col(),
                    self.mipdata_.post_solve_stack.get_orig_num_row(),
                )
            };
            root_basis
                .col_status
                .resize(basis_num_col as usize, HighsBasisStatus::kLower);
            root_basis
                .row_status
                .resize(basis_num_row as usize, HighsBasisStatus::kBasic);
            root_basis.valid = true;
            root_basis.useful = true;

            for i in 0..self.num_col() {
                let orig = self.mipdata_.post_solve_stack.get_orig_col_index(i) as usize;
                root_basis.col_status[orig] = basis.col_status[i as usize];
            }

            let num_row = basis.row_status.len() as HighsInt;
            for i in 0..num_row {
                let orig = self.mipdata_.post_solve_stack.get_orig_row_index(i) as usize;
                root_basis.row_status[orig] = basis.row_status[i as usize];
            }

            self.rootbasis = Some(Box::new(root_basis));
            rootbasis_is_local = true;
        }

        // Transform the reference of the objective limit and lower/upper
        // bounds to the original model, since offset will generally change
        // in presolve. Bound changes are transitory, so no real gap change,
        // and no update to P-D integral is necessary
        self.mipdata_.upper_limit += mip_model!(self).offset_;
        self.mipdata_.optimality_limit += mip_model!(self).offset_;

        self.mipdata_.upper_bound += mip_model!(self).offset_;
        self.mipdata_.lower_bound += mip_model!(self).offset_;

        // remove the current incumbent. Any incumbent is already transformed into the
        // original space and kept there
        self.mipdata_.incumbent.clear();
        self.mipdata_.pruned_treeweight = HighsCDouble::from(0.0);
        self.mipdata_.nodequeue.clear();
        self.mipdata_.global_orbits = None;
        // rusolve mip_symmetry_dynamic: clear dynamic symmetry data
        if self.options_mip_.mip_symmetry_dynamic {
            self.mipdata_.symmetry_dynamic = None;
        }

        // Need to be able to set presolve reduction limit separately when
        // restarting - so that bugs in presolve restart can be investigated
        // independently (see #1553)
        //
        // However, when restarting, presolve is (naturally) applied to the
        // presolved problem, so have to control the number of _further_
        // presolve reductions
        //
        // The number of further presolve reductions must be positive,
        // otherwise the MIP solver cycles, hence
        // restart_presolve_reduction_limit cannot be zero
        //
        // Although postSolveStack.numReductions() is size_t, it makes no
        // sense to use presolve_reduction_limit when the number of
        // reductions is vast
        let num_reductions = self.mipdata_.post_solve_stack.num_reductions() as HighsInt;
        let restart_presolve_reduction_limit = self.options_mip_.restart_presolve_reduction_limit;
        debug_assert!(restart_presolve_reduction_limit != 0);
        let further_presolve_reduction_limit = if restart_presolve_reduction_limit >= 0 {
            num_reductions + restart_presolve_reduction_limit
        } else {
            -1
        };
        // C++ `runMipPresolve` here is HighsMipSolverData's (HighsMipSolverData.cpp:1408),
        // not the HighsMipSolver wrapper that rebuilds `mipdata_` (HighsMipSolver.cpp:1213).
        self.run_mip_presolve_data(further_presolve_reduction_limit);

        if self.modelstatus_ != HighsModelStatus::kNotset {
            // transform the objective limit to the current model
            self.mipdata_.upper_limit -= mip_model!(self).offset_;
            self.mipdata_.optimality_limit -= mip_model!(self).offset_;

            if self.modelstatus_ == HighsModelStatus::kOptimal {
                self.mipdata_.upper_bound = 0.0;
                self.transform_new_integer_feasible_solution(&[], true);
            } else {
                self.mipdata_.upper_bound -= mip_model!(self).offset_;
            }

            // lower_bound still relates to the original model, and the offset
            // is never applied, since MIP solving is complete, and
            // lower_bound is set to upper_bound, so apply the offset now, so
            // that housekeeping in updatePrimalDualIntegral is correct
            self.mipdata_.lower_bound -= mip_model!(self).offset_;

            // There must be a gap change, since it's now zero, so always call
            // updatePrimalDualIntegral (unless solving a sub-MIP)
            //
            // Surely there must be a lower bound change
            let upper_bound = self.mipdata_.upper_bound;
            let check_prev_data = self.modelstatus_ != HighsModelStatus::kOptimal;
            self.update_lower_bound(upper_bound, true, check_prev_data);
            if self.solution_objective_ != kHighsInf
                && self.modelstatus_ == HighsModelStatus::kInfeasible
            {
                self.modelstatus_ = HighsModelStatus::kOptimal;
            }
            // HiGHS quirk: this early return does not clear `rootbasis` or `pscostinit`, which
            // point into this frame in C++ (dangling). The port keeps the owned copies.
            return;
        }
        // Bounds are currently in the original space since presolve will have
        // changed offset_
        // HIGHS_DEBUGSOL (mip_debug_solution_file): HighsMipSolverData.cpp line 1443
        self.mipdata_.debug_solution.resume();
        self.run_setup();
        if self.terminate() {
            return;
        }

        // HiGHS #3250 (mip_presolve_row_append): HighsMipSolverData.cpp line 1508 (1.16.0)
        if x1b_row_append {
            self.mipdata_
                .post_solve_stack
                .remove_cuts_from_model_x1b(num_cuts);
        } else {
            self.mipdata_
                .post_solve_stack
                .remove_cuts_from_model(num_cuts);
        }

        // HighsNodeQueue oldNodeQueue; std::swap(nodequeue, oldNodeQueue); is commented out in C++.

        // Ensure master worker is pointing to the correct cut and conflict pools
        if !self.mipdata_.workers.is_empty() {
            self.mipdata_.workers[0].set_cut_pool(0);
            self.mipdata_.workers[0].set_conflict_pool(0);
            self.mipdata_.workers[0].set_global_domain(0);
            self.mipdata_.workers[0].set_pseudocost(0);
            self.mipdata_.workers[0].upper_bound = self.mipdata_.upper_bound;
            self.mipdata_.workers[0].upper_limit = self.mipdata_.upper_limit;
            self.mipdata_.workers[0].optimality_limit = self.mipdata_.optimality_limit;
        }

        // remove the pointer into the stack-space of this function
        if rootbasis_is_local {
            self.rootbasis = None;
        }
        self.pscostinit = None;
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::basisTransfer (lines 1469-1493)
    /// C++: `void basisTransfer()`
    pub fn basis_transfer(&mut self) {
        // if a root basis is given, construct a basis for the root LP from
        // in the reduced problem space after presolving
        let num_row = self.num_row();
        let num_col = self.num_col();
        if let Some(rootbasis) = self.rootbasis.as_deref() {
            self.mipdata_.firstrootbasis.col_status =
                vec![HighsBasisStatus::kNonbasic; num_col as usize];
            self.mipdata_.firstrootbasis.row_status =
                vec![HighsBasisStatus::kNonbasic; num_row as usize];
            self.mipdata_.firstrootbasis.valid = true;
            self.mipdata_.firstrootbasis.alien = true;
            self.mipdata_.firstrootbasis.useful = true;

            // HiGHS #3250 (mip_presolve_row_append): HighsMipSolverData.cpp lines
            // 1541-1547 (1.16.0) loop over getOrigRows() / getOrigCols().
            let (rows, cols): (Vec<HighsInt>, Vec<HighsInt>) =
                if self.options_mip_.mip_presolve_row_append {
                    (
                        self.mipdata_.post_solve_stack.get_orig_rows(),
                        self.mipdata_.post_solve_stack.get_orig_cols(),
                    )
                } else {
                    ((0..num_row).collect(), (0..num_col).collect())
                };
            for i in rows {
                let orig = self.mipdata_.post_solve_stack.get_orig_row_index(i) as usize;
                let status = rootbasis.row_status[orig];
                self.mipdata_.firstrootbasis.row_status[i as usize] = status;
            }

            for i in cols {
                let orig = self.mipdata_.post_solve_stack.get_orig_col_index(i) as usize;
                let status = rootbasis.col_status[orig];
                self.mipdata_.firstrootbasis.col_status[i as usize] = status;
            }
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::getSolution (lines 1495-1497)
    /// C++: `const std::vector<double> & getSolution() const`
    pub fn get_solution(&self) -> &[f64] {
        &self.mipdata_.incumbent
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::addIncumbent (lines 1499-1607)
    /// C++: `bool addIncumbent(const std::vector<double> & sol, double solobj, const int solution_source, const bool print_display_line = true, const bool is_user_solution = false)`
    /// C++ default arguments: `print_display_line = true`, `is_user_solution = false`
    pub fn add_incumbent(
        &mut self,
        sol: &[f64],
        solobj: f64,
        solution_source: i32,
        print_display_line: bool,
        is_user_solution: bool,
    ) -> bool {
        debug_assert!(!self.parallel_lock_active());
        if self.callback_.user_callback.is_some() {
            // let an improving-solution callback name the source
            self.callback_.last_solution_source = solution_source;
        }
        let execute_mip_solution_callback = !is_user_solution
            && !self.submip
            && self.callback_.user_callback.is_some()
            && self.callback_.active[kCallbackMipSolution as usize];
        // Happens if solobj improves on the upper bound or the MIP solution
        // callback is active
        let mut solobj = solobj;
        let possibly_store_as_new_incumbent = solobj < self.mipdata_.upper_bound;
        let get_transformed_solution =
            possibly_store_as_new_incumbent || execute_mip_solution_callback;
        let transformed_solobj = if get_transformed_solution {
            self.transform_new_integer_feasible_solution(sol, possibly_store_as_new_incumbent)
        } else {
            0.0
        };
        // HiGHS quirk: `highs_solution_report` is the constant `false` in
        // HighsMipSolverData.cpp, so its debug report block is dead code and is not ported.
        if possibly_store_as_new_incumbent {
            solobj = transformed_solobj;
            if solobj >= self.mipdata_.upper_bound {
                return false;
            }

            let prev_upper_bound = self.mipdata_.upper_bound;

            self.mipdata_.upper_bound = solobj;
            let upper_bound = self.mipdata_.upper_bound;
            for worker in self.mipdata_.workers.iter_mut() {
                worker.upper_bound = upper_bound;
            }

            let bound_change = self.mipdata_.upper_bound != prev_upper_bound;
            if !self.submip && bound_change {
                let lower_bound = self.mipdata_.lower_bound;
                self.update_primal_dual_integral(
                    lower_bound,
                    lower_bound,
                    prev_upper_bound,
                    self.mipdata_.upper_bound,
                    true,
                    true,
                );
            }

            // Assigning new incumbent
            self.mipdata_.incumbent = sol.to_vec();
            self.race_publish(solobj);
            let c5_time = self.timer_.read();
            self.stats_.record_incumbent(solution_source, c5_time); // reporting-only
            // S6 (mip_heuristic_run_oneopt|crossover): SCIP primal.c lines 859-1100 at 32dd425
            let nodenum = self.mipdata_.num_nodes;
            if let Some(s6) = self.mipdata_.s6.as_mut() {
                s6_pool_insert_incumbent(s6, solobj, sol.to_vec(), solution_source, nodenum);
            }
            let new_upper_limit = self.compute_new_upper_limit(solobj, 0.0, 0.0);

            if !is_user_solution && !self.submip {
                self.save_report_mip_solution(new_upper_limit);
            }
            if new_upper_limit < self.mipdata_.upper_limit {
                self.mipdata_.num_improving_sols += 1;
                self.mipdata_.upper_limit = new_upper_limit;
                self.mipdata_.optimality_limit = self.compute_new_upper_limit(
                    solobj,
                    self.options_mip_.mip_abs_gap,
                    self.options_mip_.mip_rel_gap,
                );
                let optimality_limit = self.mipdata_.optimality_limit;
                self.mipdata_
                    .nodequeue
                    .set_optimality_limit(optimality_limit);
                let upper_limit = self.mipdata_.upper_limit;
                for worker in self.mipdata_.workers.iter_mut() {
                    worker.upper_limit = upper_limit;
                    worker.optimality_limit = optimality_limit;
                }
                let ul = self.mipdata_.upper_limit;
                self.mipdata_.debug_solution.new_incumbent_found(ul);
                detach!(self.mipdata_.domains[0], |dom| dom.propagate(self));
                if !self.mipdata_.domains[0].infeasible() {
                    detach!(self.mipdata_.redcostfixing, |rc| rc
                        .propagate_root_redcost(self));
                }

                // Two calls to printDisplayLine added for completeness,
                // ensuring that when the root node has an integer solution, a
                // logging line is issued

                if self.mipdata_.domains[0].infeasible() {
                    self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                    self.mipdata_.nodequeue.clear();
                    if print_display_line {
                        self.print_display_line(solution_source); // Added for completeness
                    }
                    return true;
                }
                detach!(self.mipdata_.cliquetable, |ct| ct.extract_obj_cliques(self));
                if self.mipdata_.domains[0].infeasible() {
                    self.mipdata_.pruned_treeweight = HighsCDouble::from(1.0);
                    self.mipdata_.nodequeue.clear();
                    if print_display_line {
                        self.print_display_line(solution_source); // Added for completeness
                    }
                    return true;
                }
                let upper_limit = self.mipdata_.upper_limit;
                let bounded = self.mipdata_.nodequeue.perform_bounding(upper_limit);
                self.mipdata_.pruned_treeweight += bounded;
                self.print_display_line(solution_source);
            }
        } else if self.mipdata_.incumbent.is_empty() {
            // Assigning new incumbent
            self.mipdata_.incumbent = sol.to_vec();
        }

        // S6 (mip_heuristic_run_oneopt|crossover): SCIP primal.c lines 859-1100 at 32dd425
        if !possibly_store_as_new_incumbent {
            let nodenum = self.mipdata_.num_nodes;
            if let Some(s6) = self.mipdata_.s6.as_mut() {
                s6_pool_insert_incumbent(s6, solobj, sol.to_vec(), solution_source, nodenum);
            }
        }

        true
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::printSolutionSourceKey (lines 1674-1713)
    /// C++: `void printSolutionSourceKey() const`
    pub fn print_solution_source_key(&self) {
        // Last MipSolutionSource enum is kSolutionSourceCleanup - which is
        // not a solution source, but used to force the last logging line to
        // be printed
        let last_enum = kSolutionSourceCount - 1;
        // Set the index of the last solution source to be printed in each
        // line of the key. Four or five can be printed, depending on the
        // lengths of the solution source strings in that line
        let limits: [i32; 4] = [4, 9, 14, last_enum];
        debug_assert!(last_enum > limits[limits.len() - 2]);

        let mut ss = String::new();
        for k in 0..limits[0] {
            if k == 0 {
                ss.push_str("\nSrc: ");
            } else {
                ss.push_str("; ");
            }
            ss.push_str(&format!(
                "{} => {}",
                self.solution_source_to_string(k, true),
                self.solution_source_to_string(k, false)
            ));
        }
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "{};\n",
            ss
        );
        let to_line = limits.len() as i32 - 1;
        for line in 0..to_line {
            let mut ss = String::new();
            for k in limits[line as usize]..limits[line as usize + 1] {
                if k == limits[line as usize] {
                    ss.push_str("     ");
                } else {
                    ss.push_str("; ");
                }
                ss.push_str(&format!(
                    "{} => {}",
                    self.solution_source_to_string(k, true),
                    self.solution_source_to_string(k, false)
                ));
            }
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "{}{}\n",
                ss,
                if line < to_line - 1 { ";" } else { "" }
            );
        }
    }

    /// HiGHS: HighsMipSolverData.cpp::HighsMipSolverData::printDisplayLine (lines 1715-1840)
    /// C++: `void printDisplayLine(const int solution_source = kSolutionSourceNone)`
    /// C++ default arguments: `solution_source = kSolutionSourceNone`
    pub fn print_display_line(&mut self, solution_source: i32) {
        // MIP logging method
        //
        // Note that if the original problem is a maximization, the cost
        // coefficients are negated so that the MIP solver only solves a
        // minimization. Hence, in preparing to print the display line, the
        // dual bound (lb) is always less than the primal bound (ub). When
        // printed, the sense of the optimization is applied so that the
        // values printed correspond to the original objective.

        // No point in computing all the logging values if logging is off
        let output_flag = self.options_mip_.log_options.output_flag;
        if !output_flag {
            return;
        }

        let timeless_log = self.options_mip_.timeless_log;
        self.mipdata_.disptime = if timeless_log {
            self.mipdata_.disptime + 1.0
        } else {
            self.timer_.read()
        };
        if solution_source == kSolutionSourceNone
            && self.mipdata_.disptime - self.mipdata_.last_disptime
                < self.options_mip_.mip_min_logging_interval
        {
            return;
        }
        self.mipdata_.last_disptime = self.mipdata_.disptime;
        let time_string = if timeless_log {
            String::new()
        } else {
            format!(" {}s", c_fmt_f_width(self.mipdata_.disptime, 1, 7))
        };

        if self.mipdata_.num_disp_lines % 20 == 0 {
            if self.mipdata_.num_disp_lines == 0 {
                self.print_solution_source_key();
            }
            let work_string0 = if timeless_log {
                "   Work"
            } else {
                "      Work      "
            };
            let work_string1 = if timeless_log {
                "LpIters"
            } else {
                "LpIters     Time"
            };
            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                "\n        Nodes      |    B&B Tree     |            Objective Bounds              |  Dynamic Constraints | {}\nSrc  Proc. InQueue |  Leaves   Expl. | BestBound       BestSol              Gap |   Cuts   InLp Confl. | {}\n\n",
                work_string0,
                work_string1
            );
        }

        self.mipdata_.num_disp_lines += 1;

        let print_nodes = convert_to_print_string(self.mipdata_.num_nodes);
        let queue_nodes = convert_to_print_string(self.mipdata_.nodequeue.num_active_nodes());
        let print_leaves =
            convert_to_print_string(self.mipdata_.num_leaves - self.mipdata_.num_leaves_before_run);

        let explored = 100.0 * f64::from(self.mipdata_.pruned_treeweight);

        let mut lb = 0.0_f64;
        let mut ub = 0.0_f64;
        let mut gap = self.limits_to_gap(
            self.mipdata_.lower_bound,
            self.mipdata_.upper_bound,
            &mut lb,
            &mut ub,
        );
        gap *= 1e2;
        // HiGHS quirk: after this clamp `ub <= objective_bound`, so the `objective_bound < ub`
        // tests on ub below are never true and the "*" suffix is never printed.
        if self.options_mip_.objective_bound < ub {
            ub = self.options_mip_.objective_bound;
        }

        let print_lp_iters = convert_to_print_string(self.mipdata_.total_lp_iterations);
        let dynamic_constraints_in_lp: HighsInt = {
            let lp = &self.mipdata_.lps[0];
            if lp.num_rows() > 0 {
                lp.num_rows() - lp.get_num_model_rows(self)
            } else {
                0
            }
        };
        let sense = f64::from(self.orig_model_.sense_ as i32);
        if self.mipdata_.upper_bound != kHighsInf {
            let gap_string = if gap >= 9999.0 {
                "Large".to_string()
            } else {
                format!("{}%", c_fmt_f(gap, 2))
            };

            let ub_string = if self.options_mip_.objective_bound < ub {
                convert_to_print_string_2(sense * ub, b"*")
            } else {
                convert_to_print_string_2(sense * ub, b"")
            };

            let lb_string = convert_to_print_string_2(sense * lb, b"");

            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                " {} {:>7} {:>7}   {:>7} {:>6}%   {:<15} {:<15} {:>8}   {:>6} {:>6} {:>6}   {:>7}{}\n",
                self.solution_source_to_string(solution_source, true),
                print_buf_str(&print_nodes),
                print_buf_str(&queue_nodes),
                print_buf_str(&print_leaves),
                c_fmt_f(explored, 2),
                print_buf_str(&lb_string),
                print_buf_str(&ub_string),
                gap_string,
                self.mipdata_.cutpools[0].get_num_cuts(),
                dynamic_constraints_in_lp,
                self.mipdata_.conflictpools[0].get_num_conflicts(),
                print_buf_str(&print_lp_iters),
                time_string
            );
        } else {
            let ub_string = if self.options_mip_.objective_bound < ub {
                convert_to_print_string_2(sense * ub, b"*")
            } else {
                convert_to_print_string_2(sense * ub, b"")
            };

            let lb_string = convert_to_print_string_2(sense * lb, b"");

            crate::highs_log_user!(
                &self.options_mip_.log_options,
                HighsLogType::kInfo,
                " {} {:>7} {:>7}   {:>7} {:>6}%   {:<15} {:<15} {:>8}   {:>6} {:>6} {:>6}   {:>7}{}\n",
                self.solution_source_to_string(solution_source, true),
                print_buf_str(&print_nodes),
                print_buf_str(&queue_nodes),
                print_buf_str(&print_leaves),
                c_fmt_f(explored, 2),
                print_buf_str(&lb_string),
                print_buf_str(&ub_string),
                c_fmt_f_width(gap, 2, 8),
                self.mipdata_.cutpools[0].get_num_cuts(),
                dynamic_constraints_in_lp,
                self.mipdata_.conflictpools[0].get_num_conflicts(),
                print_buf_str(&print_lp_iters),
                time_string
            );
        }
        // Check that limitsToBounds yields the same values for the
        // dual_bound, primal_bound (modulo optimization sense) and
        // mip_rel_gap
        let mut dual_bound = 0.0_f64;
        let mut primal_bound = 0.0_f64;
        let mut mip_rel_gap = 0.0_f64;
        self.limits_to_bounds(&mut dual_bound, &mut primal_bound, &mut mip_rel_gap);
        mip_rel_gap *= 1e2;
        debug_assert!(dual_bound == sense * lb);
        debug_assert!(primal_bound == sense * ub);
        debug_assert!(gap == mip_rel_gap);

        // Possibly interrupt from MIP logging callback
        self.callback_.clear_highs_callback_output();
        let solution_objective = self.solution_objective_;
        let interrupt = self.interrupt_from_callback_with_data(
            kCallbackMipLogging,
            solution_objective,
            "MIP logging",
        );
        debug_assert!(!interrupt);
    }
}

/// HiGHS: HighsMipSolverData.cpp::convertToPrintString (lines 1609-1633)
/// C++: `static std::array<char, 22> convertToPrintString(int64_t val)`
pub fn convert_to_print_string(val: i64) -> [u8; 22] {
    let l = cpp_max(1.0, val as f64).log10();
    let print = match l as i32 {
        0..=5 => format!("{val}"),
        6..=8 => format!("{}k", val / 1000),
        _ => format!("{}m", val / 1000000),
    };
    to_print_buf(&print)
}

/// HiGHS: HighsMipSolverData.cpp::convertToPrintString (lines 1635-1672)
/// C++: `static std::array<char, 22> convertToPrintString(double val, const char * trailingStr = "")`
/// C++ default arguments: `trailing_str = ""`
pub fn convert_to_print_string_2(val: f64, trailing_str: &[u8]) -> [u8; 22] {
    let trailing = String::from_utf8_lossy(trailing_str);
    let l = if val.abs() == kHighsInf {
        0.0
    } else {
        cpp_max(1e-6, val.abs()).log10()
    };
    let body = match l as i32 {
        0..=3 => c_fmt_g(val, 10),
        4 => c_fmt_g(val, 11),
        5 => c_fmt_g(val, 12),
        6..=10 => c_fmt_g(val, 13),
        _ => c_fmt_g(val, 9),
    };
    to_print_buf(&format!("{body}{trailing}"))
}

/// HiGHS: HighsUtils.cpp::highsSparseTranspose (lines 60-88)
/// C++: `void highsSparseTranspose(HighsInt numRow, HighsInt numCol, const std::vector<HighsInt>& Astart, const std::vector<HighsInt>& Aindex, const std::vector<double>& Avalue, std::vector<HighsInt>& ARstart, std::vector<HighsInt>& ARindex, std::vector<double>& ARvalue)`
///
/// Like C++ `resize`, the output vectors keep their existing entries where the new length allows.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
fn highs_sparse_transpose(
    num_row: HighsInt,
    num_col: HighsInt,
    astart: &[HighsInt],
    aindex: &[HighsInt],
    avalue: &[f64],
    ar_start: &mut Vec<HighsInt>,
    ar_index: &mut Vec<HighsInt>,
    ar_value: &mut Vec<f64>,
) {
    let num_row = num_row as usize;
    // Make a AR copy
    let mut iwork: Vec<HighsInt> = vec![0; num_row];
    ar_start.resize(num_row + 1, 0);
    let a_countx = aindex.len();
    ar_index.resize(a_countx, 0);
    ar_value.resize(a_countx, 0.0);
    for &row in aindex {
        debug_assert!((row as usize) < num_row);
        iwork[row as usize] += 1;
    }
    for i in 1..=num_row {
        ar_start[i] = ar_start[i - 1] + iwork[i - 1];
    }
    iwork[..num_row].copy_from_slice(&ar_start[..num_row]);
    for i_col in 0..num_col as usize {
        for k in astart[i_col] as usize..astart[i_col + 1] as usize {
            let i_row = aindex[k] as usize;
            let i_put = iwork[i_row] as usize;
            iwork[i_row] += 1;
            ar_index[i_put] = i_col as HighsInt;
            ar_value[i_put] = avalue[k];
        }
    }
}

/// C `snprintf` into a `char[22]`: at most 21 bytes are kept, then a NUL.
fn to_print_buf(s: &str) -> [u8; 22] {
    let mut buf = [0u8; 22];
    let bytes = s.as_bytes();
    let n = bytes.len().min(21);
    buf[..n].copy_from_slice(&bytes[..n]);
    buf
}

/// C `%s` of a `char[22]`: the bytes up to the first NUL.
fn print_buf_str(buf: &[u8; 22]) -> &str {
    let end = buf.iter().position(|&b| b == 0).unwrap_or(buf.len());
    std::str::from_utf8(&buf[..end]).unwrap_or("")
}
#[cfg(test)]
mod s6_hook_tests {
    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    fn s6_hook_knapsack_lp() -> HighsLp {
        // 2-row binary knapsack: maximize 3x1+2x2+4x3+1x4 subject to
        //   2x1+1x2+3x3+1x4 <= 4
        //   1x1+2x2+1x3+2x4 <= 5
        //   x_i in {0,1}
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 2;
        lp.sense_ = ObjSense::kMaximize;
        lp.offset_ = 0.0;
        lp.col_cost_ = vec![3.0, 2.0, 4.0, 1.0];
        lp.col_lower_ = vec![0.0; 4];
        lp.col_upper_ = vec![1.0; 4];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![4.0, 5.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4, 6, 8];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![2.0, 1.0, 1.0, 2.0, 3.0, 1.0, 1.0, 2.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 4];
        lp
    }

    fn s6_hook_solve<R: Send + 'static>(
        on: bool,
        f: impl FnOnce(&mut HighsMipSolver) -> R + Send + 'static,
    ) -> R {
        std::thread::spawn(move || {
            initialize_scheduler(1);
            let mut options = HighsOptions::default();
            options.output_flag = false;
            options.mip_heuristic_run_oneopt = on;
            options.mip_heuristic_run_crossover = false;
            let lp = s6_hook_knapsack_lp();
            let solution = HighsSolution::default();
            let mut callback = HighsCallback::default();
            let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
            ms.run();
            let result = f(&mut ms);
            drop(ms);
            HighsTaskExecutor::shutdown(true);
            result
        })
        .join()
        .unwrap()
    }

    #[test]
    fn s6_hook_off_pool_is_none() {
        s6_hook_solve(false, |ms| {
            assert!(ms.mipdata_.s6.is_none());
        });
    }

    #[test]
    fn s6_hook_on_pool_holds_incumbent() {
        s6_hook_solve(true, |ms| {
            let s6 = ms.mipdata_.s6.as_ref().expect("s6 should be Some");
            assert!(s6.pool.sols.len() >= 1);
            assert!(
                (s6.pool.sols[0].obj - ms.mipdata_.upper_bound).abs() < 1e-6,
                "pool[0].obj {} vs upper_bound {}",
                s6.pool.sols[0].obj,
                ms.mipdata_.upper_bound
            );
            for w in s6.pool.sols.windows(2) {
                assert!(w[0].obj <= w[1].obj + 1e-12);
            }
            let mut seen = std::collections::HashSet::new();
            for sol in &s6.pool.sols {
                assert!(seen.insert(sol.index), "duplicate index {}", sol.index);
                assert!(sol.index < s6.next_sol_index);
            }
        });
    }

    #[test]
    fn s6_hook_restart_resets_pool() {
        s6_hook_solve(true, |ms| {
            ms.run_setup();
            let s6 = ms
                .mipdata_
                .s6
                .as_ref()
                .expect("s6 should be Some after restart");
            assert_eq!(
                s6.pool.sols.len(),
                1,
                "pool should hold only the reduced incumbent after restart"
            );
        });
    }
}

#[cfg(test)]
mod h16_repair_lp_tests {
    use super::*;
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::presolve_support::calculate_row_values_quad_vec;
    use crate::mip::mip_solver::HighsMipSolver;

    /// Two columns: x0 binary (cost 0, no entries), y in [0, 10] (cost 0),
    /// and one row `y >= 1`. The incumbent `(0, 0)` violates the row, and
    /// the repair LP with x0 fixed to 0 is feasible at y = 1.
    fn repair_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 10.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![kHighsInf];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kContinuous];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        lp
    }

    /// Runs the repair branch of `transformNewIntegerFeasibleSolution` with the
    /// option set to `flag`, and returns the solver for inspection.
    fn run_repair(flag: bool) -> HighsMipSolver {
        let mut callback = HighsCallback::default();
        let mut options = HighsOptions::default();
        options.h16_repair_lp = flag;
        let lp = repair_lp();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.upper_bound = kHighsInf;
        // No presolve reductions: identity index maps for the postsolve undo.
        ms.mipdata_.post_solve_stack.initialize_index_maps(1, 2);
        ms.transform_new_integer_feasible_solution(&[0.0, 0.0], true);
        ms
    }

    fn check_repaired(ms: &HighsMipSolver) {
        assert_eq!(ms.mipdata_.total_repair_lp_feasible, 1);
        assert_eq!(ms.solution_.len(), 2);
        assert_eq!(ms.solution_[0], 0.0);
        assert!((ms.solution_[1] - 1.0).abs() < 1e-9);
        assert_eq!(ms.solution_objective_, 0.0);
    }

    #[test]
    fn h16_repair_lp_on_repairs_incumbent() {
        check_repaired(&run_repair(true));
    }

    #[test]
    fn h16_repair_lp_off_repairs_incumbent_as_1_15_1() {
        // The repair branch exists on both paths; with the flag off the result
        // is the 1.15.1 one, with no row-value recomputation after the repair LP.
        check_repaired(&run_repair(false));
    }

    #[test]
    fn h16_repair_lp_on_row_values_use_quad_sum() {
        // The on path recomputes row values of the repaired solution with
        // calculate_row_values_quad. Row 0 is 1e16*x0 + 1*x1 - 1e16*x2 at
        // x = (1, 1, 1): the quad sum is exactly 1, a double sum gives 0.
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 1;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3];
        lp.a_matrix_.index_ = vec![0, 0, 0];
        lp.a_matrix_.value_ = vec![1e16, 1.0, -1e16];
        let col_value = vec![1.0, 1.0, 1.0];
        let mut row_value = Vec::new();
        let status = calculate_row_values_quad_vec(&lp, &col_value, &mut row_value, -1);
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(row_value, vec![1.0]);
    }
}
