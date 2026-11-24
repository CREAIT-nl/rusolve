//! HiGHS source: simplex/HEkkDual.cpp (lines 28-571)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Dual simplex driver and per-solve initialisation.
//! See `HEkkDual` docs for the C++ alias -> Rust expression table.

use super::*;
use crate::prelude::*;

impl HEkkDual {
    /// HiGHS: HEkkDual.cpp::HEkkDual::solve (lines 28-393, 366 lines)
    /// C++ default arguments: pass_force_phase2 = false
    pub fn solve(&mut self, ekk: &mut HEkk, pass_force_phase2: bool) -> HighsStatus {
        // Initialise control data for a particular solve
        self.initialise_solve(ekk);

        if self.debug_dual_simplex(ekk, "Initialise", true) == HighsDebugStatus::kLogicalError {
            return ekk.return_from_solve(HighsStatus::kError);
        }
        // Assumes that the LP has a positive number of rows
        if ekk.is_unconstrained_lp() {
            return ekk.return_from_solve(HighsStatus::kError);
        }

        if !self.dual_info_ok(&*ekk, &ekk.lp_) {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kError,
                "HPrimalDual::solve has error in dual information\n"
            );
            return ekk.return_from_solve(HighsStatus::kError);
        }

        // Possibly use Li dual steepest edge weights by not storing squared
        // primal infeasibilities
        self.possibly_use_li_dual_steepest_edge(ekk);

        debug_assert!(ekk.status_.has_invert);
        if !ekk.status_.has_invert {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kError,
                "HDual:: Should enter solve with INVERT\n"
            );
            return ekk.return_from_solve(HighsStatus::kError);
        }

        // Determine the duals without cost perturbation
        ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhaseUnknown, false);
        ekk.compute_dual();
        ekk.compute_simplex_dual_infeasible();
        // Record whether the solution with unperturbed costs is dual feasible
        let dual_feasible_with_unperturbed_costs = ekk.info_.num_dual_infeasibilities == 0;
        // Force phase 2 if dual infeasibilities without cost perturbation
        // involved fixed variables or were (at most) small
        self.force_phase2 = pass_force_phase2
            || ekk.info_.max_dual_infeasibility * ekk.info_.max_dual_infeasibility
                < ekk.options_.dual_feasibility_tolerance;
        // Within the MIP solver, unless the basis supplied was alien, the
        // simplex solver should be able to start from dual feasibility, so
        // possibly debug this property. Note that debug_dual_feasible is
        // set in HEkk::setBasis, and is false if kDebugMipNodeDualFeasible
        // is false
        if ekk.debug_dual_feasible && !dual_feasible_with_unperturbed_costs {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kWarning,
                "Basis should be dual feasible, but duals without cost perturbation \
                 have num / max / sum = {:4} / {} / {} infeasibilities",
                ekk.info_.num_dual_infeasibilities,
                ekk.info_.max_dual_infeasibility,
                ekk.info_.sum_dual_infeasibilities
            );
            if !self.force_phase2 {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    " !!Not forcing phase 2!! basis Id = {}; update count = {}; \
                     name = {}\n",
                    ekk.basis_.debug_id,
                    ekk.basis_.debug_update_count,
                    ekk.basis_.debug_origin_name
                );
            } else {
                crate::highs_log_dev!(&ekk.options_.log_options, HighsLogType::kWarning, "\n");
            }
        }
        // Determine whether the solution is near-optimal. Values 1000 and
        // 1e-3 (ensuring sum<1) are unimportant, as the sum of primal
        // infeasibilities for near-optimal solutions is typically many
        // orders of magnitude smaller than 1, and the sum of primal
        // infeasibilities will be very much larger for non-trivial LPs
        // that are dual feasible for a logical or crash basis.
        //
        // Consider there to be no dual infeasibilities if there are none,
        // or if phase 2 is forced, in which case any dual infeasibilities
        // will be shifted
        let no_simplex_dual_infeasibilities =
            dual_feasible_with_unperturbed_costs || self.force_phase2;
        let near_optimal = no_simplex_dual_infeasibilities
            && ekk.info_.num_primal_infeasibilities < 1000
            && ekk.info_.max_primal_infeasibility < 1e-3;
        if near_optimal {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Dual feasible with unperturbed costs and num / max / sum primal infeasibilities of \
                 {} / {} / {}, so near-optimal\n",
                ekk.info_.num_primal_infeasibilities,
                ekk.info_.max_primal_infeasibility,
                ekk.info_.sum_primal_infeasibilities
            );
        }

        // Perturb costs according to whether the solution is near-optimal
        let perturb_costs = !near_optimal;
        if !perturb_costs {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kDetailed,
                "Near-optimal, so don't use cost perturbation\n"
            );
        }
        ekk.initialise_cost(SimplexAlgorithm::kDual, kSolvePhaseUnknown, perturb_costs);
        // Check whether the time/iteration limit has been reached. First
        // point at which a non-error return can occur
        if ekk.bailout() {
            return ekk.return_from_solve(HighsStatus::kWarning);
        }

        // Consider initialising edge weights
        if ekk.status_.has_dual_steepest_edge_weights {
            // Dual steepest edge weights are known, so possibly check
            debug_assert!(ekk.dual_edge_weight_.len() as HighsInt >= self.solver_num_row);
            debug_assert!(ekk.scattered_dual_edge_weight_.len() as HighsInt >= self.solver_num_tot);
            ekk.dev_debug_dual_steepest_edge_weights("before solve");
        } else {
            // Set up edge weights
            //
            // Assign unit weights - necessary for Dantzig and Devex, and
            // correct for steepest edge when B=I so, for clarity, do it for
            // all. Also ensure that the scattering vector is the right size
            ekk.dual_edge_weight_.clear();
            ekk.dual_edge_weight_
                .resize(self.solver_num_row as usize, 1.0);
            ekk.scattered_dual_edge_weight_
                .resize(self.solver_num_tot as usize, 0.0);
            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                // Intending to using dual steepest edge weights
                //
                // Exact DSE weights need to be computed if the basis contains
                // structurals
                if ekk.logical_basis() {
                    // Unit weights already set up for B=I
                    ekk.status_.has_dual_steepest_edge_weights = true;
                } else {
                    // Non-logical basis
                    if near_optimal || ekk.options_.lp_dual_devex_start {
                        // Use Devex rather than compute steepest edge weights
                        crate::highs_log_dev!(
                            &ekk.options_.log_options,
                            HighsLogType::kDetailed,
                            "Basis is not logical, but near-optimal, so use Devex rather \
                             than compute steepest edge weights\n"
                        );
                        self.edge_weight_mode = EdgeWeightMode::kDevex;
                        debug_assert!(!ekk.status_.has_dual_steepest_edge_weights);
                    } else {
                        // Compute steepest edge weights
                        crate::highs_log_dev!(
                            &ekk.options_.log_options,
                            HighsLogType::kDetailed,
                            "Basis is not logical, so compute steepest edge weights\n"
                        );
                        // rusolve lever `mip_lp_keep_dse_weights` (not in HiGHS): take the
                        // weights saved by variable before the last basis or row change and
                        // compute only the missing ones (simplex/ekk/dse_keep.rs).
                        if !ekk.dse_keep_apply() {
                            ekk.compute_dual_steepest_edge_weights(true);
                        }
                        ekk.status_.has_dual_steepest_edge_weights = true;
                    }
                }
            }
            if self.edge_weight_mode == EdgeWeightMode::kDevex {
                self.initialise_devex_framework(ekk);
            }
            // Check on consistency between edge_weight_mode and
            // status.has_dual_steepest_edge_weights
            if self.edge_weight_mode == EdgeWeightMode::kSteepestEdge {
                debug_assert!(ekk.status_.has_dual_steepest_edge_weights);
            } else {
                debug_assert!(!ekk.status_.has_dual_steepest_edge_weights);
            }
        }
        // Resize the copy of scattered edge weights for backtracking
        ekk.info_
            .backtracking_basis_edge_weight_
            .resize(self.solver_num_tot as usize, 0.0);

        if perturb_costs {
            // Compute the dual values with perturbed costs
            ekk.compute_dual();
            // Determine the number of dual infeasibilities after fixed
            // variable flips
            self.compute_dual_infeasibilities_with_fixed_variable_flips(ekk);
            self.dual_infeas_count = ekk.info_.num_dual_infeasibilities;
        }

        // Determine the solve phase
        if self.force_phase2 {
            // Dual infeasibilities without cost perturbation involved
            // fixed variables or were (at most) small, so can easily be
            // removed by flips for fixed variables and shifts for the rest
            self.solve_phase = kSolvePhase2;
        } else {
            // Phase depends on the number of dual infeasibilities after fixed
            // variable flips
            self.solve_phase = if self.dual_infeas_count > 0 {
                kSolvePhase1
            } else {
                kSolvePhase2
            };
        }
        if ekk.debug_ok_for_solve(SimplexAlgorithm::kDual, self.solve_phase)
            == HighsDebugStatus::kLogicalError
        {
            return ekk.return_from_solve(HighsStatus::kError);
        }
        //
        // The major solving loop
        //
        while self.solve_phase != 0 {
            let it0 = ekk.iteration_count_;
            // When starting a new phase the (updated) dual objective function
            // value isn't known. Indicate this so that when the value
            // computed from scratch in rebuild() isn't checked against the
            // the updated value
            ekk.status_.has_dual_objective_value = false;
            if self.solve_phase == kSolvePhaseUnknown {
                // Reset the phase 2 bounds so that true number of dual
                // infeasibilities can be determined
                ekk.initialise_bound(SimplexAlgorithm::kDual, kSolvePhaseUnknown, false);
                ekk.initialise_nonbasic_value_and_move();
                // Determine the number of unavoidable dual infeasibilities, and
                // hence the solve phase
                self.compute_dual_infeasibilities_with_fixed_variable_flips(ekk);
                self.dual_infeas_count = ekk.info_.num_dual_infeasibilities;
                self.solve_phase = if self.dual_infeas_count > 0 {
                    kSolvePhase1
                } else {
                    kSolvePhase2
                };
                if ekk.info_.backtracking_ {
                    // Backtracking, so set the bounds and primal values
                    ekk.initialise_bound(SimplexAlgorithm::kDual, self.solve_phase, false);
                    ekk.initialise_nonbasic_value_and_move();
                    // Can now forget that we might have been backtracking
                    ekk.info_.backtracking_ = false;
                }
            }
            debug_assert!(self.solve_phase == kSolvePhase1 || self.solve_phase == kSolvePhase2);
            if self.solve_phase == kSolvePhase1 {
                // Phase 1
                ekk.analysis_.simplex_timer_start(SimplexDualPhase1Clock, 0);
                self.solve_phase1(ekk);
                ekk.analysis_.simplex_timer_stop(SimplexDualPhase1Clock, 0);
                ekk.info_.dual_phase1_iteration_count += ekk.iteration_count_ - it0;
            } else if self.solve_phase == kSolvePhase2 {
                // Phase 2
                ekk.analysis_.simplex_timer_start(SimplexDualPhase2Clock, 0);
                self.solve_phase2(ekk);
                ekk.analysis_.simplex_timer_stop(SimplexDualPhase2Clock, 0);
                ekk.info_.dual_phase2_iteration_count += ekk.iteration_count_ - it0;
            } else {
                // Should only be kSolvePhase1 or kSolvePhase2
                ekk.model_status_ = HighsModelStatus::kSolveError;
                return ekk.return_from_solve(HighsStatus::kError);
            }
            // Return if bailing out from solve
            if ekk.solve_bailout_ {
                return ekk.return_from_solve(HighsStatus::kWarning);
            }
            // Can have all possible cases of solve_phase
            debug_assert!(self.solve_phase >= kSolvePhaseMin && self.solve_phase <= kSolvePhaseMax);
            // Look for scenarios when the major solving loop ends
            if self.solve_phase == kSolvePhaseTabooBasis {
                // Only basis change is taboo
                ekk.model_status_ = HighsModelStatus::kUnknown;
                return ekk.return_from_solve(HighsStatus::kWarning);
            }
            if self.solve_phase == kSolvePhaseError {
                // Solver error so return HighsStatus::kError
                debug_assert!(matches!(ekk.model_status_, HighsModelStatus::kSolveError));
                return ekk.return_from_solve(HighsStatus::kError);
            }
            if self.solve_phase == kSolvePhaseExit {
                // LP identified as not having an optimal solution
                debug_assert!(matches!(
                    ekk.model_status_,
                    HighsModelStatus::kUnboundedOrInfeasible | HighsModelStatus::kInfeasible
                ));
                break;
            }
            if self.solve_phase == kSolvePhaseOptimalCleanup
                || self.solve_phase == kSolvePhasePrimalInfeasibleCleanup
            {
                // Dual infeasibilities after phase 2 which ends either
                //
                // primal feasible with dual infeasibilities, so use primal
                // simplex to clean up expecting to identify optimality
                //
                // primal infeasible with dual infeasibilities so use primal
                // simplex to clean up expecting to identify (primal)
                // infeasibility
                break;
            }
            // If solve_phase == kSolvePhaseOptimal == 0 then major solving
            // loop ends naturally since solve_phase is false
        }
        // If bailing out, should have returned already
        debug_assert!(!ekk.solve_bailout_);
        // Should only have these cases
        debug_assert!(
            self.solve_phase == kSolvePhaseExit
                || self.solve_phase == kSolvePhaseUnknown
                || self.solve_phase == kSolvePhaseOptimal
                || self.solve_phase == kSolvePhaseOptimalCleanup
                || self.solve_phase == kSolvePhasePrimalInfeasibleCleanup
        );
        // Can't be solve_phase == kSolvePhase1 since this requires simplex
        // solver to have continued after identifying dual infeasibility.
        if self.solve_phase == kSolvePhaseOptimalCleanup
            || self.solve_phase == kSolvePhasePrimalInfeasibleCleanup
        {
            ekk.dual_simplex_cleanup_level_ += 1;
            if self.solve_phase == kSolvePhasePrimalInfeasibleCleanup {
                // Primal and dual infeasibilities aren't known when cleaning up
                // after suspected unboundedness in dual phase 2. All that's
                // known is that cost shifting was required to get dual
                // feasibility after removing cost perturbations, and dual
                // simplex iterations may also have been done. This is unlike
                // clean-up of dual infeasibilities after suspected optimality,
                // when no shifting and dual simplex iterations are done after
                // removing cost perturbations.
                //
                // Determine the primal and dual infeasibilities
                ekk.compute_simplex_infeasible();
            }
            if ekk.dual_simplex_cleanup_level_ > ekk.options_.max_dual_simplex_cleanup_level {
                // No clean up. Dual simplex was optimal or unbounded with
                // unperturbed costs, so say that the scaled LP has been solved
                // optimally. Optimality or infeasibility for the unscaled LP
                // are unlikely but will still be assessed honestly, so leave it
                // to the user to decide whether the solution can be accepted.
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "HEkkDual:: Cannot use level {} primal simplex cleanup for {} dual infeasibilities\n",
                    ekk.dual_simplex_cleanup_level_,
                    ekk.info_.num_dual_infeasibilities
                );
                if self.solve_phase == kSolvePhaseOptimalCleanup {
                    ekk.model_status_ = HighsModelStatus::kOptimal;
                } else {
                    ekk.model_status_ = HighsModelStatus::kInfeasible;
                }
            } else {
                // Use primal simplex to clean up. This usually yields
                // optimality or infeasibility (according to whether solve_phase
                // is kSolvePhaseOptimalCleanup or
                // kSolvePhasePrimalInfeasibleCleanup) but can yield
                // unboundedness. Time/iteration limit return is, of course,
                // possible, as are solver error
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kInfo,
                    "HEkkDual:: Using primal simplex to try to clean up num / max / sum = {} / {} / {} dual infeasibilities\n",
                    ekk.info_.num_dual_infeasibilities,
                    ekk.info_.max_dual_infeasibility,
                    ekk.info_.sum_dual_infeasibilities
                );
                let mut return_status = HighsStatus::kOk;
                ekk.analysis_
                    .simplex_timer_start(SimplexPrimalPhase2Clock, 0);
                // Switch off any bound perturbation
                let save_primal_simplex_bound_perturbation_multiplier =
                    ekk.info_.primal_simplex_bound_perturbation_multiplier;
                ekk.info_.primal_simplex_bound_perturbation_multiplier = 0.0;
                let mut primal_solver = HEkkPrimal::new(ekk);
                let call_status = primal_solver.solve(ekk, true);
                // Restore any bound perturbation
                ekk.info_.primal_simplex_bound_perturbation_multiplier =
                    save_primal_simplex_bound_perturbation_multiplier;
                ekk.analysis_
                    .simplex_timer_stop(SimplexPrimalPhase2Clock, 0);
                debug_assert!(ekk.called_return_from_solve_);
                return_status = interpret_call_status(
                    &ekk.options_.log_options,
                    call_status,
                    return_status,
                    "HEkkPrimal::solve",
                );
                // Reset called_return_from_solve_ to be false, since it's
                // called for this solve
                ekk.called_return_from_solve_ = false;
                if return_status != HighsStatus::kOk {
                    return ekk.return_from_solve(return_status);
                }
                if ekk.model_status_ == HighsModelStatus::kOptimal
                    && (ekk.info_.num_primal_infeasibilities + ekk.info_.num_dual_infeasibilities)
                        != 0
                {
                    crate::highs_log_dev!(
                        &ekk.options_.log_options,
                        HighsLogType::kWarning,
                        "HEkkDual:: Primal simplex clean up yields optimality, but with {} (max {}) primal infeasibilities and {} (max {}) dual infeasibilities\n",
                        ekk.info_.num_primal_infeasibilities,
                        ekk.info_.max_primal_infeasibility,
                        ekk.info_.num_dual_infeasibilities,
                        ekk.info_.max_dual_infeasibility
                    );
                }
            }
        }
        debug_assert!(matches!(
            ekk.model_status_,
            HighsModelStatus::kOptimal
                | HighsModelStatus::kInfeasible
                | HighsModelStatus::kUnbounded
                | HighsModelStatus::kUnboundedOrInfeasible
        ));
        if ekk.debug_ok_for_solve(SimplexAlgorithm::kDual, self.solve_phase)
            == HighsDebugStatus::kLogicalError
        {
            return ekk.return_from_solve(HighsStatus::kError);
        }
        ekk.return_from_solve(HighsStatus::kOk)
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::initialiseInstance (lines 395-430, 36 lines)
    pub fn initialise_instance(&mut self, ekk: &mut HEkk) {
        // Called in constructor for HEkkDual class
        // Copy size, matrix and simplex NLA

        self.solver_num_col = ekk.lp_.num_col_;
        self.solver_num_row = ekk.lp_.num_row_;
        self.solver_num_tot = self.solver_num_col + self.solver_num_row;
        self.inv_solver_num_row = 1.0 / (self.solver_num_row as f64);

        // Pointer copies (a_matrix, simplex_nla, analysis, jMove, workDual,
        // workValue, workRange, baseLower, baseUpper, baseValue) are removed
        // by the ownership pattern; see the HEkkDual struct docs.

        // Setup local vectors
        self.col_dse.setup(self.solver_num_row);
        self.col_bfrt.setup(self.solver_num_row);
        self.col_aq.setup(self.solver_num_row);
        self.row_ep.setup(self.solver_num_row);
        self.row_ap.setup(self.solver_num_col);

        self.dev_row_ep.setup(self.solver_num_row);
        self.dev_col_dse.setup(self.solver_num_row);

        // Setup other buffers
        self.dual_row.setup(ekk);
        self.dual_rhs.setup(ekk);
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::initialiseInstanceParallel (lines 432-480, 49 lines)
    /// SIP / PAMI setup (slices, `multi_choice` / `multi_finish` vectors); never reached with
    /// kSimplexStrategyDualPlain.
    // C++ lines 453-456 keep their two if-assignments (rule 1), not `clamp`.
    #[allow(clippy::manual_clamp)]
    pub fn initialise_instance_parallel(&mut self, ekk: &mut HEkk) {
        // No need to call this with kSimplexStrategyDualPlain
        if ekk.info_.simplex_strategy == kSimplexStrategyDualPlain {
            return;
        }

        // Identify the (current) number of HiGHS tasks to be used
        let num_concurrency = ekk.info_.num_concurrency;

        let pass_num_slice: HighsInt;
        if ekk.info_.simplex_strategy == kSimplexStrategyDualTasks {
            // Initialize for tasks
            pass_num_slice = num_concurrency - 2;
            debug_assert!(pass_num_slice > 0);
            if pass_num_slice <= 0 {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "SIP trying to use using {} slices due to concurrency ({}) being too small: results unpredictable\n",
                    pass_num_slice,
                    num_concurrency
                );
            }
        } else {
            // Initialize for multi
            self.multi_num = num_concurrency;
            if self.multi_num < 1 {
                self.multi_num = 1;
            }
            if self.multi_num > kSimplexConcurrencyLimit {
                self.multi_num = kSimplexConcurrencyLimit;
            }
            for i in 0..self.multi_num {
                self.multi_choice[i as usize].row_out = -1;
                self.multi_choice[i as usize]
                    .row_ep
                    .setup(self.solver_num_row);
                self.multi_choice[i as usize]
                    .col_aq
                    .setup(self.solver_num_row);
                self.multi_choice[i as usize]
                    .col_bfrt
                    .setup(self.solver_num_row);
            }
            pass_num_slice = cpp_max(self.multi_num - 1, 1 as HighsInt);
            debug_assert!(pass_num_slice > 0);
            if pass_num_slice <= 0 {
                crate::highs_log_dev!(
                    &ekk.options_.log_options,
                    HighsLogType::kWarning,
                    "PAMI trying to use using {} slices due to concurrency ({}) being too small: results unpredictable\n",
                    pass_num_slice,
                    num_concurrency
                );
            }
        }
        // Create the multiple HEkkDualRow instances: one for each column
        // slice
        for _i in 0..pass_num_slice {
            self.slice_dual_row.push(HEkkDualRow::default());
        }
        // Initialise the column slices
        self.init_slice(ekk, pass_num_slice);
        self.multi_iteration = 0;
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::initSlice (lines 482-538, 57 lines)
    /// Column partition of the matrix for the parallel PRICE / CHUZC slices.
    pub fn init_slice(&mut self, ekk: &mut HEkk, initial_num_slice: HighsInt) {
        // Number of slices
        self.slice_num = initial_num_slice;
        if self.slice_num < 1 {
            self.slice_num = 1;
        }
        debug_assert!(self.slice_num <= kHighsSlicedLimit);
        if self.slice_num > kHighsSlicedLimit {
            crate::highs_log_dev!(
                &ekk.options_.log_options,
                HighsLogType::kWarning,
                "WARNING: {} = slice_num > kHighsSlicedLimit = {} so truncating slice_num\n",
                self.slice_num,
                kHighsSlicedLimit
            );
            self.slice_num = kHighsSlicedLimit;
        }

        // Size the Vec fields that correspond to C++ fixed-size arrays
        self.slice_start.resize(kHighsSlicedLimit as usize + 1, 0);
        self.slice_a_matrix
            .resize(kHighsSlicedLimit as usize, HighsSparseMatrix::default());
        self.slice_ar_matrix
            .resize(kHighsSlicedLimit as usize, HighsSparseMatrix::default());
        self.slice_row_ap
            .resize(kHighsSlicedLimit as usize, HVector::default());

        // Alias to the matrix
        let acount_x = ekk.lp_.a_matrix_.start_[self.solver_num_col as usize];

        // Figure out partition weight
        let sliced_count_x = (acount_x as f64) / (self.slice_num as f64);
        self.slice_start[0] = 0;
        for i in 0..(self.slice_num - 1) {
            let mut end_column = self.slice_start[i as usize] + 1; // At least one column
            let mut end_x = ekk.lp_.a_matrix_.start_[end_column as usize];
            let stop_x = ((i + 1) as f64 * sliced_count_x) as HighsInt;
            while end_x < stop_x {
                end_column += 1;
                end_x = ekk.lp_.a_matrix_.start_[end_column as usize];
            }
            self.slice_start[(i + 1) as usize] = end_column;
            if end_column >= self.solver_num_col {
                self.slice_num = i; // SHRINK
                break;
            }
        }
        self.slice_start[self.slice_num as usize] = self.solver_num_col;

        // Partition the matrix, row_ap and related packet
        // C++ lines 520-529: sliced_Astart is written but never read; omitted.
        for i in 0..self.slice_num {
            // The matrix
            let from_col = self.slice_start[i as usize];
            let to_col = self.slice_start[(i + 1) as usize] - 1;
            let slice_num_col = self.slice_start[(i + 1) as usize] - from_col;
            // C++ lines 527-529: sliced_Astart.resize and fill are dead code; omitted.
            self.slice_a_matrix[i as usize].create_slice(&ekk.lp_.a_matrix_, from_col, to_col);
            self.slice_ar_matrix[i as usize].create_rowwise(&self.slice_a_matrix[i as usize]);

            // The row_ap and its packages
            self.slice_row_ap[i as usize].setup(slice_num_col);
            self.slice_dual_row[i as usize].setup_slice(ekk, slice_num_col);
        }
    }

    /// HiGHS: HEkkDual.cpp::HEkkDual::initialiseSolve (lines 540-571, 32 lines)
    pub fn initialise_solve(&mut self, ekk: &mut HEkk) {
        // Copy values of simplex solver options to dual simplex options
        self.primal_feasibility_tolerance = ekk.options_.primal_feasibility_tolerance;
        self.dual_feasibility_tolerance = ekk.options_.dual_feasibility_tolerance;
        self.objective_bound = ekk.options_.objective_bound;

        // Copy tolerances
        // ToDo: Eliminate these horribly-named unnecessary copies!
        self.tp = self.primal_feasibility_tolerance;
        self.td = self.dual_feasibility_tolerance;

        self.initial_basis_is_logical_ = true;
        for i_row in 0..self.solver_num_row {
            if ekk.basis_.basic_index_[i_row as usize] < self.solver_num_col {
                self.initial_basis_is_logical_ = false;
                break;
            }
        }

        let dual_edge_weight_strategy = ekk.info_.dual_edge_weight_strategy;
        self.interpret_dual_edge_weight_strategy(ekk, dual_edge_weight_strategy);

        // Initialise model and run status values
        ekk.model_status_ = HighsModelStatus::kNotset;
        ekk.solve_bailout_ = false;
        ekk.called_return_from_solve_ = false;
        ekk.exit_algorithm_ = SimplexAlgorithm::kDual;

        self.rebuild_reason = kRebuildReasonNo;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;
    use crate::prelude::*;

    #[test]
    fn w5_sip_init_slice_partition() {
        // `initialiseInstanceParallel` / `initSlice` (HEkkDual.cpp 432-538) give the C++ `slice_num`, `slice_start` and slice matrices for a small LP at `simplex_strategy` 2 and 3.
        let mut ekk = HEkk::default();
        ekk.lp_.num_col_ = 4;
        ekk.lp_.num_row_ = 3;
        // Colwise matrix with 4 columns of nonzero counts 2,1,3,2 (total 8)
        // start_ = [0, 2, 3, 6, 8]
        // index_ = [0,1, 0, 0,1,2, 1,2]
        // value_ = [1,1, 1, 1,1,1, 1,1]
        ekk.lp_.a_matrix_.num_col_ = 4;
        ekk.lp_.a_matrix_.num_row_ = 3;
        ekk.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
        ekk.lp_.a_matrix_.start_ = vec![0, 2, 3, 6, 8];
        ekk.lp_.a_matrix_.index_ = vec![0, 1, 0, 0, 1, 2, 1, 2];
        ekk.lp_.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

        ekk.info_.simplex_strategy = kSimplexStrategyDualTasks;
        ekk.info_.num_concurrency = 4;

        let mut dual = HEkkDual::default();
        dual.solver_num_col = 4;
        dual.solver_num_row = 3;
        dual.solver_num_tot = 7;

        dual.initialise_instance_parallel(&mut ekk);

        // pass_num_slice = num_concurrency - 2 = 4 - 2 = 2
        // slice_num = 2
        // AcountX = Astart[4] = 8
        // sliced_countX = 8 / 2 = 4.0
        // i=0: endColumn = slice_start[0]+1 = 1; endX = Astart[1] = 2; stopX = (0+1)*4.0 = 4
        //   2 < 4: endColumn=2, endX=Astart[2]=3; 3 < 4: endColumn=3, endX=Astart[3]=6; 6 >= 4 stop
        //   slice_start[1] = 3; endColumn(3) < solver_num_col(4), no shrink
        // slice_start[2] = solver_num_col = 4
        // So slice_start = [0, 3, 4]
        assert_eq!(dual.slice_num, 2);
        assert_eq!(dual.slice_start[0], 0);
        assert_eq!(dual.slice_start[1], 3);
        assert_eq!(dual.slice_start[2], 4);

        // slice 0: from_col=0, to_col=2, slice_num_col=3
        // slice 1: from_col=3, to_col=3, slice_num_col=1
        assert_eq!(dual.slice_a_matrix[0].num_col_, 3);
        assert_eq!(dual.slice_a_matrix[1].num_col_, 1);
        assert_eq!(dual.slice_dual_row.len(), 2);
        assert_eq!(dual.slice_dual_row[0].work_size, 3);
        assert_eq!(dual.slice_dual_row[1].work_size, 1);
        assert_eq!(dual.slice_row_ap[0].size, 3);
        assert_eq!(dual.slice_row_ap[1].size, 1);
        assert_eq!(dual.multi_iteration, 0);
    }

    fn w5_sip_init_lp() -> (HEkk, HEkkDual) {
        let mut ekk = HEkk::default();
        ekk.lp_.num_col_ = 4;
        ekk.lp_.num_row_ = 3;
        ekk.lp_.a_matrix_.num_col_ = 4;
        ekk.lp_.a_matrix_.num_row_ = 3;
        ekk.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
        ekk.lp_.a_matrix_.start_ = vec![0, 2, 3, 6, 8];
        ekk.lp_.a_matrix_.index_ = vec![0, 1, 0, 0, 1, 2, 1, 2];
        ekk.lp_.a_matrix_.value_ = vec![1.0; 8];
        let mut dual = HEkkDual::default();
        dual.solver_num_col = 4;
        dual.solver_num_row = 3;
        dual.solver_num_tot = 7;
        (ekk, dual)
    }

    #[test]
    fn w5_sip_init_multi_partition() {
        // kSimplexStrategyDualMulti, num_concurrency 3 (HEkkDual.cpp 453-463):
        // multi_num = 3, pass_num_slice = max(3 - 1, 1) = 2, same partition as SIP above.
        let (mut ekk, mut dual) = w5_sip_init_lp();
        ekk.info_.simplex_strategy = kSimplexStrategyDualMulti;
        ekk.info_.num_concurrency = 3;
        dual.initialise_instance_parallel(&mut ekk);
        assert_eq!(dual.multi_num, 3);
        for i in 0..3 {
            assert_eq!(dual.multi_choice[i].row_out, -1);
            assert_eq!(dual.multi_choice[i].row_ep.size, 3);
            assert_eq!(dual.multi_choice[i].col_aq.size, 3);
            assert_eq!(dual.multi_choice[i].col_bfrt.size, 3);
        }
        assert_eq!(dual.slice_dual_row.len(), 2);
        assert_eq!(dual.slice_num, 2);
        assert_eq!(&dual.slice_start[..3], &[0, 3, 4]);
    }

    #[test]
    fn w5_sip_init_multi_clamp_and_shrink() {
        // num_concurrency 10 > kSimplexConcurrencyLimit: multi_num = 8, pass_num_slice = 7.
        // initSlice: sliced_countX = 8 / 7.0; stopX = 1, 2, 3, 4 for i = 0..3, giving
        // slice_start[1..=4] = 1, 2, 3, 4; at i = 3 endColumn 4 >= solver_num_col, so
        // slice_num = 3 (SHRINK, line 513) and slice_start[3] = solver_num_col = 4.
        let (mut ekk, mut dual) = w5_sip_init_lp();
        ekk.info_.simplex_strategy = kSimplexStrategyDualMulti;
        ekk.info_.num_concurrency = 10;
        dual.initialise_instance_parallel(&mut ekk);
        assert_eq!(dual.multi_num, kSimplexConcurrencyLimit);
        assert_eq!(dual.slice_dual_row.len(), 7);
        assert_eq!(dual.slice_num, 3);
        assert_eq!(&dual.slice_start[..5], &[0, 1, 2, 4, 4]);
        assert_eq!(dual.slice_a_matrix[0].num_col_, 1);
        assert_eq!(dual.slice_a_matrix[1].num_col_, 1);
        assert_eq!(dual.slice_a_matrix[2].num_col_, 2);
        assert_eq!(dual.slice_dual_row[2].work_size, 2);
        // Untouched slices keep their default work_size (only slice_num of them are set up).
        assert_eq!(dual.slice_dual_row[3].work_size, -1);
    }

    #[test]
    fn w5_sip_init_plain_returns_early() {
        let (mut ekk, mut dual) = w5_sip_init_lp();
        ekk.info_.simplex_strategy = kSimplexStrategyDualPlain;
        ekk.info_.num_concurrency = 4;
        dual.multi_iteration = 5;
        dual.initialise_instance_parallel(&mut ekk);
        assert!(dual.slice_dual_row.is_empty());
        assert_eq!(dual.multi_iteration, 5);
    }
}
