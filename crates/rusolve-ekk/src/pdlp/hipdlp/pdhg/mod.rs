//! HiGHS source: pdlp/hipdlp/pdhg.hpp, pdlp/hipdlp/pdhg.cc (lines 1-492)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `class PDLPSolver`, CPU build. `pdhg.cc` is 2910 lines, split by
//! contiguous line range (one `impl` block per file):
//!
//! | Rust module    | lines                    | contents |
//! |----------------|--------------------------|----------|
//! | this file      | 1-492                    | struct, getters, `printConstraintInfo`, `preprocessLp`, `postprocess` |
//! | [`solve`]      | 494-739                  | `solve`, `computeFixedPointError` |
//! | [`iterate`]    | 784-1136                 | `runConvergenceCheck`, `checkRestartCriteria`, `performPdhgStep`, `performHalpernPdhgStep`, `accumulateAverages`, `prepareNextIteration`, `solveReturn` |
//! | [`converge`]   | 1138-1527                | `initialize`, averages, feasibility / slack / gap / objective, `checkConvergence` |
//! | [`power`]      | 1529-1670, 1783-1901     | `powerMethod`, `setup`, `scaleProblem`, `unscaleSolution`, `logSummary` |
//! | [`update`]     | 1944-2359                | step sizes, primal weight, `updateX` / `updateY`, `updateIterates*`, stability helpers |
//!
//! `pdhg.cc` 1903-1938 (the `PrimalDualParams` / `MalitskyPockParams` /
//! `AdaptiveLinesearchParams::initialise` methods) is in
//! [`super::defs`]. Not ported: GPU code 742-781, 1087-1125, 1673-1780,
//! 2395-2909; the profiling / debug helpers 2362-2388 (`PDLP_PROFILE == 0`);
//! the `#if PDLP_DEBUG_LOG` / `#if PDLP_PROFILE` members; `hipdlpTimer*`
//! calls. `performHalpernStep` and `applyHalpernAveraging` are declared in
//! `pdhg.hpp` but never defined; they are not ported.
//!
//! # Borrowing
//!
//! * `original_lp_` is `Option<&'a HighsLp>` (`passLp`), `highs_timer_p_`
//!   is `Option<&'a HighsTimer>` (`setup`). Both are only read.
//! * A method whose C++ parameters alias members (`checkConvergence(...,
//!   x_current_, ..., dSlackPos_, dSlackNeg_)`, `computeAverageIterate(
//!   Ax_avg_, ATy_avg_)`, `updateAverageIterates(x_current_, ...)`) takes
//!   slices; the caller `std::mem::take`s the member vectors, calls, and
//!   puts them back, so the callee sees the same values.
//! * `computeDualSlacks` tests pointer identity `&dSlackPos == &dSlackPos_
//!   && &dSlackNeg == &dSlackNeg_`; Rust passes that fact as
//!   `slacks_are_members: bool`, threaded through `computeDualFeasibility`
//!   and `checkConvergence` (true exactly for the `"[L]"` calls in
//!   `runConvergenceCheck`).
//! * `scaling_.passLp(&processed_lp)` / `passParams` are not stored; the
//!   scaling methods take `&mut self.lp_` and `&self.params_` (disjoint
//!   field borrows).
//! * `restart_scheme_.updateBeta(b)` is `self.restart_scheme_.update_beta(b,
//!   &mut self.params_)`.

pub mod converge;
pub mod iterate;
pub mod power;
pub mod solve;
pub mod update;

use super::defs::*;
use super::linalg;
use super::logger::Logger;
use super::restart::RestartScheme;
use super::scaling::Scaling;
use super::solver_results::*;
use crate::prelude::*;

/// HiGHS: pdhg.hpp `class PDLPSolver`
#[derive(Clone, Debug)]
pub struct PDLPSolver<'a> {
    // --- Data Members ---
    pub(crate) lp_: HighsLp,
    pub(crate) original_lp_: Option<&'a HighsLp>,
    pub(crate) unscaled_processed_lp_: HighsLp,

    pub(crate) params_: PrimalDualParams,
    pub(crate) stepsize_: StepSizeConfig,
    pub(crate) logger_: Logger,
    pub(crate) log_options_: HighsLogOptions,
    pub(crate) results_: SolverResults,

    pub(crate) scaling_: Scaling,
    pub(crate) restart_scheme_: RestartScheme,

    // Problem dimensions & metadata
    pub(crate) original_num_col_: HighsInt,
    pub(crate) num_eq_rows_: HighsInt,
    pub(crate) sense_origin_: HighsInt,
    pub(crate) unscaled_rhs_norm_: f64,
    pub(crate) unscaled_c_norm_: f64,
    pub(crate) is_equality_row_: Vec<bool>,
    pub(crate) constraint_new_idx_: Vec<HighsInt>,
    /// `std::vector<ConstraintType>`: HConst.h plain enum `EQ`..`FREE`.
    pub(crate) constraint_types_: Vec<HighsInt>,

    // Iteration State Vectors
    pub(crate) x_current_: Vec<f64>,
    pub(crate) y_current_: Vec<f64>,
    pub(crate) x_next_: Vec<f64>,
    pub(crate) y_next_: Vec<f64>,
    pub(crate) x_avg_: Vec<f64>,
    pub(crate) y_avg_: Vec<f64>,
    pub(crate) x_sum_: Vec<f64>,
    pub(crate) y_sum_: Vec<f64>,
    pub(crate) x_at_last_restart_: Vec<f64>,
    pub(crate) y_at_last_restart_: Vec<f64>,
    /// For over-relaxed Halpern
    pub(crate) reflected_x_: Vec<f64>,
    pub(crate) reflected_y_: Vec<f64>,
    /// For Halpern
    pub(crate) x_anchor_: Vec<f64>,
    pub(crate) y_anchor_: Vec<f64>,

    // Caches
    pub(crate) Ax_cache_: Vec<f64>,
    pub(crate) ATy_cache_: Vec<f64>,
    pub(crate) Ax_next_: Vec<f64>,
    pub(crate) ATy_next_: Vec<f64>,
    // Average Caches
    pub(crate) Ax_avg_: Vec<f64>,
    pub(crate) ATy_avg_: Vec<f64>,
    // Residual Caches
    pub(crate) K_times_x_diff_: Vec<f64>,
    pub(crate) dSlackPos_: Vec<f64>,
    pub(crate) dSlackNeg_: Vec<f64>,
    pub(crate) dSlackPosAvg_: Vec<f64>,
    pub(crate) dSlackNegAvg_: Vec<f64>,
    pub(crate) halpern_dual_slack_next_: Vec<f64>,
    pub(crate) halpern_dual_slack_next_valid_: bool,
    /// For Halpern restart
    pub(crate) initial_fpe_: f64,
    pub(crate) fpe_: f64,

    // Scalars
    pub(crate) final_iter_count_: HighsInt,
    pub(crate) num_rejected_steps_: HighsInt,
    pub(crate) halpern_iteration_: HighsInt,
    pub(crate) sum_weights_: f64,
    pub(crate) current_eta_: f64,
    pub(crate) ratio_last_two_step_sizes_: f64,

    // PID
    pub(crate) primal_weight_: f64,
    pub(crate) best_primal_weight_: f64,
    pub(crate) primal_weight_error_sum_: f64,
    pub(crate) primal_weight_last_error_: f64,
    pub(crate) best_primal_dual_residual_gap_: f64,

    pub(crate) highs_timer_p_: Option<&'a HighsTimer>,
}

impl Default for PDLPSolver<'_> {
    /// The C++ default member initialisers; `HighsLp` members are
    /// `HighsLp::default()` (a cleared LP).
    fn default() -> Self {
        Self {
            lp_: HighsLp::default(),
            original_lp_: None,
            unscaled_processed_lp_: HighsLp::default(),
            params_: PrimalDualParams::default(),
            stepsize_: StepSizeConfig::default(),
            logger_: Logger::default(),
            log_options_: HighsLogOptions::default(),
            results_: SolverResults::default(),
            scaling_: Scaling::default(),
            restart_scheme_: RestartScheme::default(),
            original_num_col_: 0,
            num_eq_rows_: 0,
            sense_origin_: 1,
            unscaled_rhs_norm_: 0.0,
            unscaled_c_norm_: 0.0,
            is_equality_row_: Vec::new(),
            constraint_new_idx_: Vec::new(),
            constraint_types_: Vec::new(),
            x_current_: Vec::new(),
            y_current_: Vec::new(),
            x_next_: Vec::new(),
            y_next_: Vec::new(),
            x_avg_: Vec::new(),
            y_avg_: Vec::new(),
            x_sum_: Vec::new(),
            y_sum_: Vec::new(),
            x_at_last_restart_: Vec::new(),
            y_at_last_restart_: Vec::new(),
            reflected_x_: Vec::new(),
            reflected_y_: Vec::new(),
            x_anchor_: Vec::new(),
            y_anchor_: Vec::new(),
            Ax_cache_: Vec::new(),
            ATy_cache_: Vec::new(),
            Ax_next_: Vec::new(),
            ATy_next_: Vec::new(),
            Ax_avg_: Vec::new(),
            ATy_avg_: Vec::new(),
            K_times_x_diff_: Vec::new(),
            dSlackPos_: Vec::new(),
            dSlackNeg_: Vec::new(),
            dSlackPosAvg_: Vec::new(),
            dSlackNegAvg_: Vec::new(),
            halpern_dual_slack_next_: Vec::new(),
            halpern_dual_slack_next_valid_: false,
            initial_fpe_: 0.0,
            fpe_: 0.0,
            final_iter_count_: 0,
            num_rejected_steps_: 0,
            halpern_iteration_: 0,
            sum_weights_: 0.0,
            current_eta_: 0.0,
            ratio_last_two_step_sizes_: 1.0,
            primal_weight_: 1.0,
            best_primal_weight_: 1.0,
            primal_weight_error_sum_: 0.0,
            primal_weight_last_error_: 0.0,
            best_primal_dual_residual_gap_: f64::INFINITY,
            highs_timer_p_: None,
        }
    }
}

impl<'a> PDLPSolver<'a> {
    /// HiGHS: pdhg.hpp::PDLPSolver::passLp
    pub fn passLp(&mut self, lp: &'a HighsLp) {
        self.original_lp_ = Some(lp);
    }

    /// HiGHS: pdhg.hpp::PDLPSolver::getTerminationCode
    pub fn getTerminationCode(&self) -> TerminationStatus {
        self.results_.term_code
    }

    /// HiGHS: pdhg.hpp::PDLPSolver::getIterationCount
    pub fn getIterationCount(&self) -> HighsInt {
        self.final_iter_count_
    }

    /// HiGHS: pdhg.hpp::PDLPSolver::getNumCol
    pub fn getNumCol(&self) -> HighsInt {
        self.lp_.num_col_
    }

    /// HiGHS: pdhg.hpp::PDLPSolver::getNumRow
    pub fn getNumRow(&self) -> HighsInt {
        self.lp_.num_row_
    }

    /// HiGHS: pdhg.cc::PDLPSolver::printConstraintInfo (lines 37-150)
    ///
    /// Counts constraint and bound types and prints them through
    /// `logger_.detailed`; nothing it computes is stored. Print only.
    pub(crate) fn printConstraintInfo(&self) {}

    /// HiGHS: pdhg.cc::PDLPSolver::preprocessLp (lines 152-357)
    /// HiGHS 1.16 (#3469): pdhg.cc::PDLPSolver::preprocessLp (lines 150-365)
    /// With `h16_hipdlp_max`: costs are multiplied by `sense_origin_`, the
    /// offset by `sense_origin_` too, and `FREE` rows are removed (no row,
    /// no slack column, `constraint_new_idx_` = -1).
    ///
    /// The cuPDLP formulation of `original_lp_` into `lp_`: row
    /// classification (`EQ` / `BOUND` / `GEQ` / `LEQ` / `FREE` against
    /// `kHighsInf`), row permutation `EQ` / `BOUND` / `FREE` first, a slack
    /// column per `BOUND` / `FREE` row (bounds = the row bounds, row becomes
    /// `= 0`), `LEQ` rows negated, `is_equality_row_`,
    /// `constraint_new_idx_`, `num_eq_rows_`, `sense_origin_` (costs are
    /// copied unchanged; only `postprocess` applies the sense),
    /// `original_num_col_`. Each column's `(new_row, value)` pairs are
    /// sorted with `std::sort` on the pair (row first; rows are distinct).
    /// Then
    /// `unscaled_processed_lp_ = lp_` and the unscaled norms
    /// (`linalg::vectorNorm(col_cost_)`, `linalg::vectorNorm(row_lower_)`,
    /// both `p = 2`). `scaling_.passLp` is not stored.
    pub fn preprocessLp(&mut self) {
        // pdhg.cc 153-155: #if PDLP_PROFILE hipdlpTimerStart(kHipdlpClockPreprocess); dropped (PDLP_PROFILE == 0).
        self.logger_
            .detailed("Preprocessing LP using cupdlp formulation (slack variables for bounds)...");

        let original_lp = self.original_lp_.unwrap();
        let nRows_orig = original_lp.num_row_;
        let nCols_orig = original_lp.num_col_;
        // HiGHS 1.16 #3469 (h16_hipdlp_max): read once; HiGHS 1.16.0 pdhg.cc lines 150-365
        let h16_hipdlp_max = self.params_.h16_hipdlp_max;

        if original_lp.a_matrix_.is_rowwise() {
            self.logger_
                .info("Original LP matrix must be in column-wise format,");
        }

        // HiGHS 1.16 #3469 (h16_hipdlp_max): HiGHS 1.16.0 pdhg.cc lines 168-169 set the offset after sense_origin_
        if !h16_hipdlp_max {
            self.lp_.offset_ = original_lp.offset_;
        }

        let mut num_new_cols: HighsInt = 0;
        let mut nEqs: HighsInt = 0;
        // HiGHS 1.16 #3469 (h16_hipdlp_max): HiGHS 1.16.0 pdhg.cc line 167
        let mut num_free_rows: HighsInt = 0;
        self.sense_origin_ = if original_lp.sense_ == ObjSense::kMinimize {
            1
        } else {
            -1
        };
        // HiGHS 1.16 #3469 (h16_hipdlp_max): HiGHS 1.16.0 pdhg.cc line 169
        if h16_hipdlp_max {
            self.lp_.offset_ = original_lp.offset_ * self.sense_origin_ as f64;
        }
        self.constraint_types_.resize(nRows_orig as usize, 0);

        // 1. First pass: Classify constraints and count slack variables needed
        for i in 0..nRows_orig {
            let has_lower = original_lp.row_lower_[i as usize] > -kHighsInf;
            let has_upper = original_lp.row_upper_[i as usize] < kHighsInf;

            if has_lower && has_upper {
                if original_lp.row_lower_[i as usize] == original_lp.row_upper_[i as usize] {
                    self.constraint_types_[i as usize] = EQ;
                    nEqs += 1;
                } else {
                    self.constraint_types_[i as usize] = BOUND;
                    num_new_cols += 1;
                    nEqs += 1;
                }
            } else if has_lower {
                self.constraint_types_[i as usize] = GEQ;
            } else if has_upper {
                self.constraint_types_[i as usize] = LEQ;
            } else {
                self.constraint_types_[i as usize] = FREE;
                // HiGHS 1.16 #3469 (h16_hipdlp_max): HiGHS 1.16.0 pdhg.cc lines 189-192
                if h16_hipdlp_max {
                    num_free_rows += 1;
                } else {
                    num_new_cols += 1;
                    nEqs += 1;
                }
            }
        }

        // 2. Set new dimensions
        self.lp_.num_col_ = nCols_orig + num_new_cols;
        // HiGHS 1.16 #3469 (h16_hipdlp_max): HiGHS 1.16.0 pdhg.cc line 199 (free rows are removed)
        self.lp_.num_row_ = if h16_hipdlp_max {
            nRows_orig - num_free_rows
        } else {
            nRows_orig
        };
        self.original_num_col_ = nCols_orig;
        self.constraint_new_idx_.resize(nRows_orig as usize, 0);

        self.lp_.col_cost_.resize(self.lp_.num_col_ as usize, 0.0);
        self.lp_.col_lower_.resize(self.lp_.num_col_ as usize, 0.0);
        self.lp_.col_upper_.resize(self.lp_.num_col_ as usize, 0.0);
        self.lp_.row_lower_.resize(self.lp_.num_row_ as usize, 0.0);
        self.lp_.row_upper_.resize(self.lp_.num_row_ as usize, 0.0);

        // 3. Determine row permutation: EQ/BOUND/FREE first, then LEQ/GEQ
        // (HiGHS 1.16 #3469: EQ/BOUND first, FREE rows get -1, HiGHS 1.16.0 pdhg.cc lines 211-221)
        let mut eq_idx: HighsInt = 0;
        let mut ineq_idx: HighsInt = nEqs;
        for i in 0..nRows_orig {
            let is_eq_block = if h16_hipdlp_max {
                self.constraint_types_[i as usize] == EQ
                    || self.constraint_types_[i as usize] == BOUND
            } else {
                self.constraint_types_[i as usize] == EQ
                    || self.constraint_types_[i as usize] == BOUND
                    || self.constraint_types_[i as usize] == FREE
            };
            if is_eq_block {
                self.constraint_new_idx_[i as usize] = eq_idx;
                eq_idx += 1;
            } else if h16_hipdlp_max && self.constraint_types_[i as usize] == FREE {
                self.constraint_new_idx_[i as usize] = -1;
            } else {
                self.constraint_new_idx_[i as usize] = ineq_idx;
                ineq_idx += 1;
            }
        }

        // 4. Create is_equality_row_ vector based on NEW ordering
        // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc lines 223-230 (free rows skipped, size num_row_)
        let is_equality_len = if h16_hipdlp_max {
            self.lp_.num_row_
        } else {
            nRows_orig
        };
        self.is_equality_row_
            .resize(is_equality_len as usize, false);
        for i in 0..nRows_orig {
            if h16_hipdlp_max && self.constraint_types_[i as usize] == FREE {
                continue;
            }
            let new_row_idx = self.constraint_new_idx_[i as usize];
            self.is_equality_row_[new_row_idx as usize] = if h16_hipdlp_max {
                self.constraint_types_[i as usize] == EQ
                    || self.constraint_types_[i as usize] == BOUND
            } else {
                self.constraint_types_[i as usize] == EQ
                    || self.constraint_types_[i as usize] == BOUND
                    || self.constraint_types_[i as usize] == FREE
            };
        }

        // 5. Populate costs and bounds for original and new slack variables
        for i in 0..nCols_orig {
            // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc line 234 (cost times sense_origin_)
            self.lp_.col_cost_[i as usize] = if h16_hipdlp_max {
                original_lp.col_cost_[i as usize] * self.sense_origin_ as f64
            } else {
                original_lp.col_cost_[i as usize]
            };
            self.lp_.col_lower_[i as usize] = original_lp.col_lower_[i as usize];
            self.lp_.col_upper_[i as usize] = original_lp.col_upper_[i as usize];
        }

        let mut current_slack_col: HighsInt = nCols_orig;
        for i in 0..nRows_orig {
            // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc lines 240-241 (free rows get no slack column)
            let slack_row = if h16_hipdlp_max {
                self.constraint_types_[i as usize] == BOUND
            } else {
                self.constraint_types_[i as usize] == BOUND
                    || self.constraint_types_[i as usize] == FREE
            };
            if slack_row {
                self.lp_.col_cost_[current_slack_col as usize] = 0.0;
                self.lp_.col_lower_[current_slack_col as usize] =
                    original_lp.row_lower_[i as usize];
                self.lp_.col_upper_[current_slack_col as usize] =
                    original_lp.row_upper_[i as usize];
                current_slack_col += 1;
            }
        }

        // 6. Set the new RHS (row bounds) in PERMUTED order
        for i in 0..nRows_orig {
            // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc line 251 (free rows skipped)
            if h16_hipdlp_max && self.constraint_types_[i as usize] == FREE {
                continue;
            }
            let new_row_idx = self.constraint_new_idx_[i as usize];
            // HiGHS 1.16 #3469: the BOUND arm takes FREE rows only on the 1.15.1 path
            // (1.16.0 pdhg.cc lines 267-273)
            let bound_arm = if h16_hipdlp_max {
                self.constraint_types_[i as usize] == BOUND
            } else {
                self.constraint_types_[i as usize] == BOUND
                    || self.constraint_types_[i as usize] == FREE
            };
            match self.constraint_types_[i as usize] {
                EQ => {
                    self.lp_.row_lower_[new_row_idx as usize] = original_lp.row_lower_[i as usize];
                    self.lp_.row_upper_[new_row_idx as usize] = original_lp.row_upper_[i as usize];
                }
                GEQ => {
                    self.lp_.row_lower_[new_row_idx as usize] = original_lp.row_lower_[i as usize];
                    self.lp_.row_upper_[new_row_idx as usize] = kHighsInf;
                }
                LEQ => {
                    // Becomes -ax >= -b
                    self.lp_.row_lower_[new_row_idx as usize] = -original_lp.row_upper_[i as usize];
                    self.lp_.row_upper_[new_row_idx as usize] = kHighsInf;
                }
                _ if bound_arm => {
                    // Becomes ax - z = 0
                    self.lp_.row_lower_[new_row_idx as usize] = 0.0;
                    self.lp_.row_upper_[new_row_idx as usize] = 0.0;
                }
                _ => {}
            }
        }

        // 7. Build the constraint matrix with proper row ordering
        // First, get the original matrix in column-wise format
        let mut original_matrix = original_lp.a_matrix_.clone();
        original_matrix.ensure_colwise();

        // Create new matrix
        self.lp_.a_matrix_.clear();
        self.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
        self.lp_.a_matrix_.num_col_ = self.lp_.num_col_;
        self.lp_.a_matrix_.num_row_ = self.lp_.num_row_;

        // Reserve space
        self.lp_
            .a_matrix_
            .start_
            .resize((self.lp_.num_col_ + 1) as usize, 0);
        self.lp_
            .a_matrix_
            .index_
            .reserve((original_matrix.num_nz() + num_new_cols) as usize);
        self.lp_
            .a_matrix_
            .value_
            .reserve((original_matrix.num_nz() + num_new_cols) as usize);

        // Build column by column
        self.lp_.a_matrix_.start_[0] = 0;
        for col in 0..nCols_orig {
            // For each column, add entries in the new row order
            let mut col_entries: Vec<(HighsInt, f64)> = Vec::new();

            let mut el = original_matrix.start_[col as usize];
            while el < original_matrix.start_[(col + 1) as usize] {
                let old_row = original_matrix.index_[el as usize];
                // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc lines 305-308 (skip free rows;
                // the C++ `continue` also advances the loop, so `el += 1` first)
                if h16_hipdlp_max && self.constraint_types_[old_row as usize] == FREE {
                    el += 1;
                    continue;
                }
                let new_row = self.constraint_new_idx_[old_row as usize];
                let mut val = original_matrix.value_[el as usize];

                // Apply scaling for LEQ constraints
                if self.constraint_types_[old_row as usize] == LEQ {
                    val = -val;
                }

                col_entries.push((new_row, val));
                el += 1;
            }

            // Sort by row index to maintain proper sparse matrix format
            // PORT NOTE: std::sort on std::pair<HighsInt, double> compares
            // (row, value) lexicographically; tuple `partial_cmp` is the same
            // order. Rows in a column are distinct, so stability is moot.
            col_entries.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            // Add to matrix
            for entry in &col_entries {
                self.lp_.a_matrix_.index_.push(entry.0);
                self.lp_.a_matrix_.value_.push(entry.1);
            }

            self.lp_.a_matrix_.start_[(col + 1) as usize] =
                self.lp_.a_matrix_.index_.len() as HighsInt;
        }

        // Add slack variable columns
        current_slack_col = nCols_orig;
        for i in 0..nRows_orig {
            // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc lines 336 (BOUND only)
            let slack_row = if h16_hipdlp_max {
                self.constraint_types_[i as usize] == BOUND
            } else {
                self.constraint_types_[i as usize] == BOUND
                    || self.constraint_types_[i as usize] == FREE
            };
            if slack_row {
                let row_idx = self.constraint_new_idx_[i as usize];
                self.lp_.a_matrix_.index_.push(row_idx);
                self.lp_.a_matrix_.value_.push(-1.0);
                self.lp_.a_matrix_.start_[(current_slack_col + 1) as usize] =
                    self.lp_.a_matrix_.index_.len() as HighsInt;
                current_slack_col += 1;
            }
        }

        self.num_eq_rows_ = nEqs;
        // pdhg.cc:339 scaling_.passLp(&processed_lp); dropped (Scaling::passLp not stored).
        self.unscaled_processed_lp_ = self.lp_.clone();

        // 8. Compute and store norms of unscaled cost and rhs
        self.unscaled_c_norm_ = linalg::vectorNorm(&self.lp_.col_cost_, 2.0);
        self.unscaled_rhs_norm_ = linalg::vectorNorm(&self.lp_.row_lower_, 2.0);

        self.logger_.detailed(&format!(
            "Preprocessing complete. New dimensions: {} rows, {} cols.",
            self.lp_.num_row_, self.lp_.num_col_
        ));
        self.logger_.detailed(&format!(
            "Unscaled norms: ||c|| = {}, ||b|| = {}",
            self.unscaled_c_norm_, self.unscaled_rhs_norm_
        ));

        self.printConstraintInfo();
        // pdhg.cc 354-356: #if PDLP_PROFILE hipdlpTimerStop(kHipdlpClockPreprocess); dropped (PDLP_PROFILE == 0).
    }

    /// HiGHS: pdhg.cc::PDLPSolver::postprocess (lines 359-492)
    /// HiGHS 1.16 (#3469): pdhg.cc::PDLPSolver::postprocess (lines 367-506)
    /// With `h16_hipdlp_max`: an empty `y_current_` is allowed when the LP
    /// has no rows, and the row duals are multiplied by `sense_origin_`
    /// (`FREE` rows get 0).
    ///
    /// Returns before writing when `original_lp_` is `None`
    /// (`DIMENSION_MISMATCH`), `x_current_` / `y_current_` empty
    /// (`INVALID_SOLUTION`) or of the wrong size; `NUMERICAL_ERROR` on a
    /// non-finite primal value. Writes `results_.primal_obj` and divides
    /// `unscaled_processed_lp_` cost / matrix by the scale factors when
    /// scaled (state kept even though nothing reads it later).
    /// `orig_matrix.ensureColwise()` acts on a copy.
    pub fn postprocess(&mut self, solution: &mut HighsSolution) -> PostSolveRetcode {
        self.logger_.detailed("Post-solving the solution...");

        // pdhg.cc 362-364: #if PDLP_PROFILE hipdlpTimerStart(kHipdlpClockPostprocess); dropped (PDLP_PROFILE == 0).

        let original_lp = match self.original_lp_ {
            Some(lp) => lp,
            None => return PostSolveRetcode::DIMENSION_MISMATCH,
        };

        // HiGHS 1.16 #3469 (h16_hipdlp_max): HiGHS 1.16.0 pdhg.cc line 377
        let h16_hipdlp_max = self.params_.h16_hipdlp_max;
        let y_bad = if h16_hipdlp_max {
            self.y_current_.is_empty() && self.lp_.num_row_ > 0
        } else {
            self.y_current_.is_empty()
        };
        if self.x_current_.is_empty() || y_bad {
            return PostSolveRetcode::INVALID_SOLUTION;
        }

        if self.x_current_.len() != self.lp_.num_col_ as usize
            || self.y_current_.len() != self.lp_.num_row_ as usize
        {
            self.logger_.info(&format!(
                "Solution dimension mismatch: x_current size={} vs expected={}, y_current size={} vs expected={}",
                self.x_current_.len(),
                self.lp_.num_col_,
                self.y_current_.len(),
                self.lp_.num_row_
            ));
            return PostSolveRetcode::DIMENSION_MISMATCH;
        }

        // 2. Resize solution object to original dimensions
        solution
            .col_value
            .resize(self.original_num_col_ as usize, 0.0);
        solution
            .row_value
            .resize(original_lp.num_row_ as usize, 0.0);
        solution
            .col_dual
            .resize(self.original_num_col_ as usize, 0.0);
        solution.row_dual.resize(original_lp.num_row_ as usize, 0.0);

        // 3. Recover Primal Column Values (x)
        for i in 0..original_lp.num_col_ {
            if i >= self.x_current_.len() as HighsInt {
                self.logger_.info(&format!(
                    "Index {} out of bounds for x_current_ of size {}",
                    i,
                    self.x_current_.len()
                ));
                return PostSolveRetcode::DIMENSION_MISMATCH;
            }

            solution.col_value[i as usize] = self.x_current_[i as usize];

            if !solution.col_value[i as usize].is_finite() {
                self.logger_.info(&format!(
                    "Non-finite primal variable value at index {}: {}",
                    i, solution.col_value[i as usize]
                ));
                return PostSolveRetcode::NUMERICAL_ERROR;
            }
        }

        let mut final_primal_objective = original_lp.offset_;
        for i in 0..self.original_num_col_ {
            final_primal_objective +=
                original_lp.col_cost_[i as usize] * solution.col_value[i as usize];
        }
        self.results_.primal_obj = final_primal_objective;

        // 4. Recover Dual Row Values (y)
        let y_reordered = self.y_current_.clone();
        for orig_row in 0..original_lp.num_row_ {
            let reordered_row = self.constraint_new_idx_[orig_row as usize];

            // HiGHS 1.16 #3469: HiGHS 1.16.0 pdhg.cc lines 428-431 (free rows have zero duals)
            if h16_hipdlp_max && self.constraint_types_[orig_row as usize] == FREE {
                solution.row_dual[orig_row as usize] = 0.0;
                continue;
            }

            // Get the dual value from the reordered position
            let dual_value = y_reordered[reordered_row as usize];

            // Apply sign correction for LEQ constraints
            // HiGHS 1.16 #3469: the sense factor is applied on the 1.16.0 path (lines 437-441)
            if self.constraint_types_[orig_row as usize] == LEQ {
                solution.row_dual[orig_row as usize] = if h16_hipdlp_max {
                    -dual_value * self.sense_origin_ as f64
                } else {
                    -dual_value
                };
            } else {
                solution.row_dual[orig_row as usize] = if h16_hipdlp_max {
                    dual_value * self.sense_origin_ as f64
                } else {
                    dual_value
                };
            }
        }

        // 5. Recover Primal Row Values (ax)
        let row_scale = self.scaling_.getRowScaling().to_vec();
        let col_scale = self.scaling_.getColScaling().to_vec();
        if self.scaling_.isScaled() {
            // Unscale matrix, costs, and rhs
            for iCol in 0..self.unscaled_processed_lp_.num_col_ {
                self.unscaled_processed_lp_.col_cost_[iCol as usize] /= col_scale[iCol as usize];
                let mut iEl = self.unscaled_processed_lp_.a_matrix_.start_[iCol as usize];
                while iEl < self.unscaled_processed_lp_.a_matrix_.start_[(iCol + 1) as usize] {
                    let iRow = self.unscaled_processed_lp_.a_matrix_.index_[iEl as usize];
                    self.unscaled_processed_lp_.a_matrix_.value_[iEl as usize] /=
                        row_scale[iRow as usize] * col_scale[iCol as usize];
                    iEl += 1;
                }
            }
        }

        let mut ax_original = vec![0.0; original_lp.num_row_ as usize];

        // Get original matrix in column-wise format
        let mut orig_matrix = original_lp.a_matrix_.clone();
        orig_matrix.ensure_colwise();

        // Compute ax using only the original columns (not slack variables)
        for col in 0..self.original_num_col_ {
            let x_val = self.x_current_[col as usize]; // Use unscaled x values

            let mut el = orig_matrix.start_[col as usize];
            while el < orig_matrix.start_[(col + 1) as usize] {
                let row = orig_matrix.index_[el as usize];
                let a_val = orig_matrix.value_[el as usize];
                ax_original[row as usize] += a_val * x_val;
                el += 1;
            }
        }

        // Now ax_original contains the correct row activity values
        for orig_row in 0..original_lp.num_row_ {
            solution.row_value[orig_row as usize] = ax_original[orig_row as usize];
        }

        // 6. Recover Dual Column Values (Reduced Costs)
        // The duals on the variable bounds l <= x <= u are the reduced costs.
        // In the PDLP framework, these are given by dSlackPos - dSlackNeg.
        for i in 0..self.original_num_col_ {
            if i >= self.dSlackPos_.len() as HighsInt || i >= self.dSlackNeg_.len() as HighsInt {
                self.logger_.info(&format!(
                    "Index {} out of bounds for dSlackPos/Neg of size {}",
                    i,
                    self.dSlackPos_.len()
                ));
                return PostSolveRetcode::DIMENSION_MISMATCH;
            }
            solution.col_dual[i as usize] =
                self.dSlackPos_[i as usize] - self.dSlackNeg_[i as usize];
            solution.col_dual[i as usize] *= self.sense_origin_ as f64;
        }

        solution.value_valid = true; // to do
        solution.dual_valid = true;
        self.logger_.detailed("Post-solve complete.");

        // pdhg.cc 488-490: #if PDLP_PROFILE hipdlpTimerStop(kHipdlpClockPostprocess); dropped (PDLP_PROFILE == 0).
        PostSolveRetcode::OK
    }
}

#[cfg(test)]
pub(crate) mod tests {
    //! Unit tests.
    use super::*;

    /// min x0 + x1 s.t. r0: x0 + x1 >= 1 (GEQ), r1: x0 - x1 <= 2 (LEQ),
    /// r2: 0 <= x0 <= 3 (BOUND), r3: x1 = 0.5 (EQ); 0 <= x <= 10.
    pub(crate) fn four_row_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 4;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.row_lower_ = vec![1.0, -kHighsInf, 0.0, 0.5];
        lp.row_upper_ = vec![kHighsInf, 2.0, 3.0, 0.5];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 6];
        lp.a_matrix_.index_ = vec![0, 1, 2, 0, 1, 3];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0, -1.0, 1.0];
        lp
    }

    #[test]
    fn preprocess_orders_equalities_first_and_adds_slack() {
        let lp = four_row_lp();
        let mut s = PDLPSolver::default();
        s.passLp(&lp);
        s.preprocessLp();
        // EQ/BOUND rows (r2, r3) first, one slack column for the BOUND row.
        assert_eq!(s.num_eq_rows_, 2);
        assert_eq!(s.lp_.num_row_, 4);
        assert_eq!(s.lp_.num_col_, 3);
        assert_eq!(s.original_num_col_, 2);
        assert_eq!(s.constraint_types_, vec![GEQ, LEQ, BOUND, EQ]);
        assert_eq!(s.sense_origin_, 1);
    }

    #[test]
    fn postprocess_without_lp_is_dimension_mismatch() {
        let mut s = PDLPSolver::default();
        let mut sol = HighsSolution::default();
        assert_eq!(
            s.postprocess(&mut sol),
            PostSolveRetcode::DIMENSION_MISMATCH
        );
    }

    /// Runs `solve_lp_hipdlp` on `lp` with `h16_hipdlp_max` set to `h16`.
    fn h16_hipdlp_max_solve(lp: &HighsLp, h16: bool) -> (HighsModelStatus, HighsSolution) {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.h16_hipdlp_max = h16;
        let mut timer = HighsTimer::default();
        let mut basis = HighsBasis::default();
        let mut solution = HighsSolution::default();
        let mut status = HighsModelStatus::kNotset;
        let mut info = HighsInfo::default();
        let mut callback = HighsCallback::default();
        let rc = crate::pdlp::hipdlp_wrapper::solve_lp_hipdlp(
            &options,
            &mut timer,
            lp,
            &mut basis,
            &mut solution,
            &mut status,
            &mut info,
            &mut callback,
        );
        assert_eq!(rc, HighsStatus::kOk);
        (status, solution)
    }

    /// `offset + c . x` over the original columns of `lp`.
    fn h16_hipdlp_max_objective(lp: &HighsLp, solution: &HighsSolution) -> f64 {
        let mut objective = lp.offset_;
        for i in 0..lp.num_col_ as usize {
            objective += lp.col_cost_[i] * solution.col_value[i];
        }
        objective
    }

    /// max 3x + 2y s.t. x + y <= 4, x + 3y <= 6, x, y >= 0. Optimum 12 at (4, 0).
    fn h16_hipdlp_max_lp(offset: f64) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.sense_ = ObjSense::kMaximize;
        lp.offset_ = offset;
        lp.col_cost_ = vec![3.0, 2.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![4.0, 6.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 3.0];
        lp
    }

    /// min x0 + 2 x1 + 3 x2 with row 1 free (x0 - x2), as in TestHiPdlp
    /// "hi-pdlp-free-rows". Minimum 6 (at (2, 2, 0)), maximum 11 (at (0, 1, 3)).
    /// With x0 = 4 - x1 - x2 the objective is 4 + x1 + 2 x2, and row 0 caps x2
    /// at 3 and x1 + x2 at 4, so the maximum is 4 + 1 + 6 = 11.
    fn h16_hipdlp_max_free_row_lp(sense: ObjSense) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 4;
        lp.sense_ = sense;
        lp.col_cost_ = vec![1.0, 2.0, 3.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf, kHighsInf];
        lp.row_lower_ = vec![1.0, -kHighsInf, 2.0, 4.0];
        lp.row_upper_ = vec![kHighsInf, kHighsInf, kHighsInf, 4.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 6, 9];
        lp.a_matrix_.index_ = vec![0, 1, 3, 0, 2, 3, 1, 2, 3];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0];
        lp
    }

    /// Maximise on: the maximise LP reaches its optimum 12, not the minimum 0.
    #[test]
    fn h16_hipdlp_max_maximise_on() {
        let lp = h16_hipdlp_max_lp(0.0);
        let (status, solution) = h16_hipdlp_max_solve(&lp, true);
        assert_eq!(status, HighsModelStatus::kOptimal);
        let objective = h16_hipdlp_max_objective(&lp, &solution);
        assert!(
            (objective - 12.0).abs() < 1e-3 * (1.0 + 12.0),
            "{objective}"
        );
    }

    /// Maximise off: the 1.15.1 path minimises the un-negated costs of the same
    /// LP, so the reported "maximum" is the minimum 0 (wrong; recorded).
    #[test]
    fn h16_hipdlp_max_maximise_off_records_1_15_1() {
        let lp = h16_hipdlp_max_lp(0.0);
        let (status, solution) = h16_hipdlp_max_solve(&lp, false);
        assert_eq!(status, HighsModelStatus::kOptimal);
        let objective = h16_hipdlp_max_objective(&lp, &solution);
        assert!(objective.abs() < 1e-3, "{objective}");
    }

    /// The offset is reported on the original scale: 12 + 5 = 17.
    #[test]
    fn h16_hipdlp_max_offset_on() {
        let lp = h16_hipdlp_max_lp(5.0);
        let (status, solution) = h16_hipdlp_max_solve(&lp, true);
        assert_eq!(status, HighsModelStatus::kOptimal);
        let objective = h16_hipdlp_max_objective(&lp, &solution);
        assert!(
            (objective - 17.0).abs() < 1e-3 * (1.0 + 17.0),
            "{objective}"
        );
    }

    /// Free rows are removed: both senses reach their optimum, and the free
    /// row has a zero dual.
    #[test]
    fn h16_hipdlp_max_free_rows_on() {
        for (sense, expected) in [(ObjSense::kMinimize, 6.0), (ObjSense::kMaximize, 11.0)] {
            let lp = h16_hipdlp_max_free_row_lp(sense);
            let (status, solution) = h16_hipdlp_max_solve(&lp, true);
            assert_eq!(status, HighsModelStatus::kOptimal);
            let objective = h16_hipdlp_max_objective(&lp, &solution);
            assert!(
                (objective - expected).abs() < 1e-3 * (1.0 + expected),
                "{sense:?}: {objective}"
            );
            assert_eq!(solution.row_dual.len(), 4);
            assert_eq!(solution.row_dual[1], 0.0);
        }
    }

    /// Free rows off: the 1.15.1 path keeps the free row (as a slack column). The
    /// minimum is still 6. The maximum comes out as the minimum 6 (wrong; recorded).
    #[test]
    fn h16_hipdlp_max_free_rows_off_records_1_15_1() {
        let lp = h16_hipdlp_max_free_row_lp(ObjSense::kMinimize);
        let (status, solution) = h16_hipdlp_max_solve(&lp, false);
        assert_eq!(status, HighsModelStatus::kOptimal);
        let objective = h16_hipdlp_max_objective(&lp, &solution);
        assert!((objective - 6.0).abs() < 1e-3 * (1.0 + 6.0), "{objective}");

        let lp = h16_hipdlp_max_free_row_lp(ObjSense::kMaximize);
        let (status, solution) = h16_hipdlp_max_solve(&lp, false);
        assert_eq!(status, HighsModelStatus::kOptimal);
        let objective = h16_hipdlp_max_objective(&lp, &solution);
        assert!((objective - 6.0).abs() < 1e-3 * (1.0 + 6.0), "{objective}");
    }

    /// Every row is free: no rows are left after removal, so the LP is solved
    /// with zero rows. min x0 + 2 x1 on boxes gives x = (0, 1), objective 2.
    #[test]
    fn h16_hipdlp_max_all_free_rows_on() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![1.0, 2.0];
        lp.col_lower_ = vec![0.0, 1.0];
        lp.col_upper_ = vec![3.0, 2.0];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, kHighsInf];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 4];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, -1.0];
        let (status, solution) = h16_hipdlp_max_solve(&lp, true);
        assert_eq!(status, HighsModelStatus::kOptimal);
        let objective = h16_hipdlp_max_objective(&lp, &solution);
        assert!((objective - 2.0).abs() < 1e-3, "{objective}");
        assert_eq!(solution.row_dual, vec![0.0, 0.0]);
    }
}
