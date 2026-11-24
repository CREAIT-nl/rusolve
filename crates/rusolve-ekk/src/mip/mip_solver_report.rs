//! HiGHS source: mip/HighsMipSolver.cpp (cleanupSolve, solvingReport, callbackGetCutPool,
//! getGapString, solutionFeasible)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//! Split out of `mip_solver.rs` to keep the driver file readable.

use crate::io::highs_io::highs_double_to_string;
use crate::lp_data::presolve_support::{calculate_row_values_quad_vec, is_basis_consistent};
use crate::prelude::*;
use crate::util::c_printf::{c_fmt_f, c_fmt_g};

/// Copy a C string into the fixed `[u8; 128]` buffer `snprintf` writes (NUL-terminated,
/// truncated to 127 bytes).
fn c_buffer(s: &str) -> [u8; 128] {
    let mut buf = [0u8; 128];
    let bytes = s.as_bytes();
    let n = bytes.len().min(127);
    buf[..n].copy_from_slice(&bytes[..n]);
    buf
}

/// HiGHS: HighsMipSolver.cpp::getGapString (lines 1252-1285)
/// C++: `static std::array<char, 128> getGapString(const double gap_, const double primal_bound_, const HighsOptions * options_mip_)`
///
/// `options_mip_` is always non-null in C++ (it is dereferenced), so `None` panics.
pub fn get_gap_string(
    gap_: f64,
    primal_bound_: f64,
    options_mip_: Option<&HighsOptions>,
) -> [u8; 128] {
    let options_mip_ =
        options_mip_.expect("get_gap_string: options_mip_ is null (C++ dereferences it)");
    if gap_ == kHighsInf {
        return c_buffer("inf");
    }
    let mut print_tol = cpp_max(cpp_min(1e-2, 1e-1 * gap_), 1e-6);
    let gap_val_string = highs_double_to_string(100.0 * gap_, print_tol);
    let mut gap_tol = options_mip_.mip_rel_gap;

    if options_mip_.mip_abs_gap > options_mip_.mip_feasibility_tolerance {
        gap_tol = if primal_bound_ == 0.0 {
            kHighsInf
        } else {
            cpp_max(gap_tol, options_mip_.mip_abs_gap / primal_bound_.abs())
        };
    }

    if gap_tol == 0.0 {
        c_buffer(&format!("{}%", gap_val_string))
    } else if gap_tol != kHighsInf {
        print_tol = cpp_max(cpp_min(1e-2, 1e-1 * gap_tol), 1e-6);
        let gap_tol_string = highs_double_to_string(100.0 * gap_tol, print_tol);
        c_buffer(&format!(
            "{}% (tolerance: {}%)",
            gap_val_string, gap_tol_string
        ))
    } else {
        c_buffer(&format!("{}% (tolerance: inf)", gap_val_string))
    }
}

/// Clock names of `profiling_->name[clock]` for the three clocks `solvingReport` reports
/// (HighsTimer.h / HighsInterface.cpp:4299-4301).
const CLOCK_NAMES: [&str; 3] = ["Presolve", "Solve", "Postsolve"];

impl HighsMipSolver {
    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::cleanupSolve (lines 1025-1121)
    /// C++: `void cleanupSolve()`
    ///
    pub fn cleanup_solve(&mut self) {
        for w in 0..self.mipdata_.workers.len() {
            debug_assert!(self.mipdata_.workers[w].solutions_.is_empty());
        }
        if self.terminator_active() {
            if self.terminator_terminated() {
                // Indicate that this instance has been interrupted
                self.modelstatus_ = HighsModelStatus::kHighsInterrupt;
            } else if !self.submip {
                // When sub-MIPs call cleanupSolve(), they generally don't have a
                // termination criterion for the whole MIP solver.
                self.terminator_terminate();
            }
        }

        // Force a final logging line
        self.print_display_line(kSolutionSourceCleanup);
        // Stop the solve clock - which won't be running if presolve determines the model status
        if self.profiling_.running(kSolveTime) {
            self.profiling_.stop(kSolveTime);
        }

        // Need to complete the calculation of P-D integral, checking for NO gap change
        let lb = self.mipdata_.lower_bound;
        let ub = self.mipdata_.upper_bound;
        self.update_primal_dual_integral(lb, lb, ub, ub, false, true);

        self.profiling_.start(kPostsolveTime, false);
        // HiGHS 1.16 #3135 (h16_mip_start): read the option once; off = the 1.15.1 statements
        let h16_mip_start = self.options_mip_.h16_mip_start;
        let havesolution = self.solution_objective_ != kHighsInf;
        let feasible = if havesolution {
            let feasible = self.bound_violation_ <= self.options_mip_.mip_feasibility_tolerance
                && self.integrality_violation_ <= self.options_mip_.mip_feasibility_tolerance
                && self.row_violation_ <= self.options_mip_.mip_feasibility_tolerance;
            if h16_mip_start && feasible && self.mipdata_.upper_bound == kHighsInf {
                // HiGHS 1.16 #3135: HighsMipSolver.cpp cleanupSolve lines 1112-1125 (1.16.0)
                self.check_add_solution();
                if self.modelstatus_ == HighsModelStatus::kInfeasible
                    && self.mipdata_.upper_bound != kHighsInf
                {
                    self.mipdata_.lower_bound = self.mipdata_.upper_bound;
                }
            }
            feasible
        } else {
            false
        };

        self.dual_bound_ = self.mipdata_.lower_bound;
        if self.mipdata_.objective_function.is_integral() {
            let integral_scale = self.mipdata_.objective_function.integral_scale();
            let rounded_lower_bound =
                (self.mipdata_.lower_bound * integral_scale - self.mipdata_.feastol).ceil()
                    / integral_scale;
            self.dual_bound_ = cpp_max(self.dual_bound_, rounded_lower_bound);
        }
        let offset = mip_model!(self).offset_;
        self.dual_bound_ += offset;
        self.primal_bound_ = self.mipdata_.upper_bound + offset;
        self.node_count_ = self.mipdata_.num_nodes;
        self.total_lp_iterations_ = self.mipdata_.total_lp_iterations;
        self.dual_bound_ = cpp_min(self.dual_bound_, self.primal_bound_);
        self.primal_dual_integral_ = self.mipdata_.primal_dual_integral.value;

        // adjust objective sense in case of maximization problem
        if self.orig_model_.sense_ == ObjSense::kMaximize {
            self.dual_bound_ = -self.dual_bound_;
            self.primal_bound_ = -self.primal_bound_;
        }

        if self.modelstatus_ == HighsModelStatus::kNotset
            || self.modelstatus_ == HighsModelStatus::kInfeasible
        {
            if feasible && havesolution {
                self.modelstatus_ = HighsModelStatus::kOptimal;
            } else {
                self.modelstatus_ = HighsModelStatus::kInfeasible;
            }
        }

        self.profiling_.stop(kPostsolveTime);
        self.timer_.stop();

        let p = self.primal_bound_;
        let d = self.dual_bound_;
        self.gap_ = (p - d).abs();
        if p == 0.0 {
            self.gap_ = if d == 0.0 { 0.0 } else { kHighsInf };
        } else if p != kHighsInf {
            self.gap_ = (p - d).abs() / p.abs();
        } else {
            self.gap_ = kHighsInf;
        }

        if self.options_mip_.output_flag {
            // solutionstatus is only used in solvingReport, but depends on the values of
            // havesolution and feasible that are local to cleanupSolve
            let solutionstatus = if havesolution {
                if feasible { "feasible" } else { "infeasible" }
            } else {
                "-"
            };
            self.solving_report(solutionstatus);
        }

        // fclose(improving_solution_file_): dropping the last handle closes the file
        self.improving_solution_file_ = None;
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::solvingReport (lines 1123-1210)
    /// C++: `void solvingReport(const std::string & solutionstatus) const`
    ///
    /// PORT NOTE: `profiling_->read(clock, record)` / `numCall` are not available on the
    /// profiling stub (mip/mip_timer.rs), so both record reads are 0 and the call counts
    /// are 0: only the `"%.2f (%s)"` line per clock is printed, with 0.00.
    pub fn solving_report(&self, solutionstatus: &str) {
        let gap_string = get_gap_string(self.gap_, self.primal_bound_, Some(&self.options_mip_));
        let gap_end = gap_string.iter().position(|&b| b == 0).unwrap_or(128);
        let gap_str = String::from_utf8_lossy(&gap_string[..gap_end]).into_owned();

        let timeless_log = self.options_mip_.timeless_log;
        // HiGHS 1.16 #3272 (h16_mip_start): read the option once; off = the 1.15.1 statements
        let h16_mip_start = self.options_mip_.h16_mip_start;
        // HiGHS 1.16 #3521 #3578 (h16_mip_workers): read the option once; off = the 1.15.1 statements
        let h16_mip_workers = self.options_mip_.h16_mip_workers;
        let pd_integral = if h16_mip_workers {
            // HiGHS 1.16 #3521 (h16_mip_workers): HighsMipSolver.cpp solvingReport line 1335 (1.16.0)
            c_fmt_f(self.mipdata_.primal_dual_integral.value, 2)
        } else {
            c_fmt_g(self.mipdata_.primal_dual_integral.value, 12)
        };
        let opts = &self.options_mip_.log_options;
        crate::highs_log_user!(opts, HighsLogType::kInfo, "\nSolving report\n");
        if !self.orig_model_.model_name_.is_empty() {
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Model             {}\n",
                self.orig_model_.model_name_
            );
        }
        if h16_mip_start {
            // HiGHS 1.16 #3272 (h16_mip_start): HighsMipSolver.cpp solvingReport lines 1282-1300 (1.16.0)
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Status            {}\n",
                util_model_status_to_string(self.modelstatus_)
            );
            // HiGHS 1.16 #3578 (h16_mip_workers): HighsMipSolver.cpp solvingReport lines 1324-1326 (1.16.0)
            let report_bounds = if h16_mip_workers {
                self.modelstatus_ != HighsModelStatus::kSolveError
                    && self.modelstatus_ != HighsModelStatus::kMemoryLimit
                    && self.modelstatus_ != HighsModelStatus::kPresolveError
            } else {
                self.modelstatus_ != HighsModelStatus::kSolveError
            };
            if report_bounds {
                crate::highs_log_user!(
                    opts,
                    HighsLogType::kInfo,
                    "  Primal bound      {}\n  Dual bound        {}\n  Gap               {}\n",
                    c_fmt_g(self.primal_bound_, 12),
                    c_fmt_g(self.dual_bound_, 12),
                    gap_str
                );
                if !timeless_log {
                    crate::highs_log_user!(
                        opts,
                        HighsLogType::kInfo,
                        "  P-D integral      {}\n",
                        pd_integral
                    );
                }
            }
        } else {
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Status            {}\n  Primal bound      {}\n  Dual bound        {}\n  Gap               {}\n",
                util_model_status_to_string(self.modelstatus_),
                c_fmt_g(self.primal_bound_, 12),
                c_fmt_g(self.dual_bound_, 12),
                gap_str
            );
            if !timeless_log {
                crate::highs_log_user!(
                    opts,
                    HighsLogType::kInfo,
                    "  P-D integral      {}\n",
                    pd_integral
                );
            }
        }
        if self.options_mip_.mip_shift_and_propagate && !self.submip {
            let s = &self.stats_;
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Shift-and-prop    {} calls, {} found, {} probings, {} cutoffs, {} rows violated, {} s\n",
                s.shiftprop_calls,
                s.shiftprop_found,
                s.shiftprop_probings,
                s.shiftprop_cutoffs,
                s.shiftprop_violations_left,
                c_fmt_g(s.shiftprop_s, 4)
            );
        }
        if self.options_mip_.mip_node_rounding > 0 && !self.submip {
            let s = &self.stats_;
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Node rounding     {} calls, {} / {} / {} found (simple / rounding / shifting), {} polished, {} improving, {} s\n",
                s.noderound_calls,
                s.noderound_found[0],
                s.noderound_found[1],
                s.noderound_found[2],
                s.noderound_polished,
                s.noderound_improving,
                c_fmt_g(s.noderound_s, 4)
            );
        }
        if self.options_mip_.mip_alns && !self.submip {
            let s = &self.stats_;
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  ALNS              {} calls, {} improving, {} nodes, {} work\n",
                s.alns_calls,
                s.alns_found,
                s.alns_nodes,
                s.alns_work
            );
        }
        if (self.options_mip_.mip_proximity || self.options_mip_.mip_primal_stall) && !self.submip {
            let s = &self.stats_;
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Proximity         {} calls, {} improving, {} nodes, {} work\n",
                s.prox_calls,
                s.prox_found,
                s.prox_nodes,
                s.prox_work
            );
        }
        if self.options_mip_.mip_local_branching && !self.submip {
            let s = &self.stats_;
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Local branching   {} calls, {} improving, {} nodes, {} work\n",
                s.lb_calls,
                s.lb_found,
                s.lb_nodes,
                s.lb_work
            );
        }
        crate::highs_log_user!(
            opts,
            HighsLogType::kInfo,
            "  Solution status   {}\n",
            solutionstatus
        );
        if solutionstatus != "-" {
            if h16_mip_start {
                // HiGHS 1.16 #3272 (h16_mip_start): HighsMipSolver.cpp solvingReport lines 1305-1322 (1.16.0)
                crate::highs_log_user!(
                    opts,
                    HighsLogType::kInfo,
                    "                    {} (objective)\n",
                    c_fmt_g(self.solution_objective_, 12)
                );
                if self.bound_violation_ != 0.0 {
                    crate::highs_log_user!(
                        opts,
                        HighsLogType::kInfo,
                        "                    {} (variable bound violation)\n",
                        c_fmt_g(self.bound_violation_, 12)
                    );
                }
                if self.row_violation_ != 0.0 {
                    crate::highs_log_user!(
                        opts,
                        HighsLogType::kInfo,
                        "                    {} (constraint bound violation)\n",
                        c_fmt_g(self.row_violation_, 12)
                    );
                }
                if self.integrality_violation_ != 0.0 {
                    crate::highs_log_user!(
                        opts,
                        HighsLogType::kInfo,
                        "                    {} (integrality violation)\n",
                        c_fmt_g(self.integrality_violation_, 12)
                    );
                }
            } else {
                crate::highs_log_user!(
                    opts,
                    HighsLogType::kInfo,
                    "                    {} (objective)\n                    {} (bound viol.)\n                    {} (int. viol.)\n                    {} (row viol.)\n",
                    c_fmt_g(self.solution_objective_, 12),
                    c_fmt_g(self.bound_violation_, 12),
                    c_fmt_g(self.integrality_violation_, 12),
                    c_fmt_g(self.row_violation_, 12)
                );
            }
        }
        if !timeless_log {
            let total = self.timer_.read();
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Timing            {}\n",
                c_fmt_f(total, 2)
            );
            for clock in [kPresolveTime, kSolveTime, kPostsolveTime] {
                // C++ callRecord(clock): mip_time / submip_time / call counts from the profiler
                let mip_time = self.profiling_.read(clock);
                let submip_time = 0.0;
                let mip_calls: HighsInt = self.profiling_.num_call(clock);
                let submip_calls: HighsInt = 0;
                let total_time = mip_time + submip_time;
                // Only log postsolve if it's written as nonzero
                crate::highs_log_user!(
                    opts,
                    HighsLogType::kInfo,
                    "                    {} ({})\n",
                    c_fmt_f(total_time, 2),
                    CLOCK_NAMES[clock as usize]
                );
                if mip_calls > 1 || submip_calls > 0 {
                    crate::highs_log_user!(
                        opts,
                        HighsLogType::kInfo,
                        "                        MIP    time [calls] = {} [{}]\n",
                        c_fmt_f(mip_time, 2),
                        mip_calls
                    );
                    if submip_calls > 0 {
                        crate::highs_log_user!(
                            opts,
                            HighsLogType::kInfo,
                            "                        subMIP time [calls] = {} [{}]\n",
                            c_fmt_f(submip_time, 2),
                            submip_calls
                        );
                    }
                }
            }
        }
        crate::highs_log_user!(
            opts,
            HighsLogType::kInfo,
            "  Max sub-MIP depth {}\n  Nodes             {}\n",
            self.max_submip_level,
            self.mipdata_.num_nodes
        );
        if self.mipdata_.total_repair_lp != 0 {
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "  Repair LPs        {} ({} feasible; {} iterations)\n",
                self.mipdata_.total_repair_lp,
                self.mipdata_.total_repair_lp_feasible,
                self.mipdata_.total_repair_lp_iterations
            );
        } else {
            crate::highs_log_user!(opts, HighsLogType::kInfo, "  Repair LPs        0\n");
        }
        crate::highs_log_user!(
            opts,
            HighsLogType::kInfo,
            "  LP iterations     {}\n",
            self.mipdata_.total_lp_iterations
        );
        if self.mipdata_.total_lp_iterations != 0 {
            crate::highs_log_user!(
                opts,
                HighsLogType::kInfo,
                "                    {} (strong br.)\n                    {} (separation)\n                    {} (heuristics)\n",
                self.mipdata_.sb_lp_iterations,
                self.mipdata_.sepa_lp_iterations,
                self.mipdata_.heuristic_lp_iterations
            );
        }
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::callbackGetCutPool (lines 1231-1250)
    /// C++: `void callbackGetCutPool() const`
    ///
    /// `&mut self`: writes `callback_->data_out` (C++ `const`, through the pointer).
    pub fn callback_get_cut_pool(&mut self) {
        debug_assert!(self.callback_.callback_active(kCallbackMipGetCutPool));
        self.callback_.clear_highs_callback_output();

        let mut cutpool_num_col: HighsInt = 0;
        let mut cutpool_num_cut: HighsInt = 0;
        let mut cutpool_lower: Vec<f64> = Vec::new();
        let mut cutpool_upper: Vec<f64> = Vec::new();
        let mut cut_matrix = HighsSparseMatrix::default();
        self.mipdata_.lps[0].get_cut_pool(
            self,
            &mut cutpool_num_col,
            &mut cutpool_num_cut,
            &mut cutpool_lower,
            &mut cutpool_upper,
            &mut cut_matrix,
        );

        let data_out = &mut self.callback_.data_out;
        data_out.cutpool_num_col = cutpool_num_col;
        data_out.cutpool_num_cut = cutpool_num_cut;
        data_out.cutpool_lower = cutpool_lower;
        data_out.cutpool_upper = cutpool_upper;
        // take ownership
        data_out.cutpool_start = std::mem::take(&mut cut_matrix.start_);
        data_out.cutpool_index = std::mem::take(&mut cut_matrix.index_);
        data_out.cutpool_value = std::mem::take(&mut cut_matrix.value_);

        let interrupt = self.interrupt_from_callback_with_data(
            kCallbackMipGetCutPool,
            self.solution_objective_,
            "MIP cut pool",
        );
        debug_assert!(!interrupt);
    }

    /// HiGHS: HighsMipSolver.cpp::HighsMipSolver::solutionFeasible (lines 1287-1363)
    /// C++: `bool solutionFeasible(const HighsLp * lp, const std::vector<double> & col_value, const std::vector<double> * pass_row_value, double & bound_violation, double & row_violation, double & integrality_violation, HighsCDouble & obj) const`
    ///
    /// PORT NOTE: `lp == None` stands for the C++ `orig_model_` argument (the constructor
    /// `HighsMipSolver::HighsMipSolver` is the only caller that passes the original model;
    /// it passes `None` here because `orig_model_` is `self`).
    pub fn solution_feasible(
        &self,
        lp: Option<&HighsLp>,
        col_value: &[f64],
        pass_row_value: Option<&Vec<f64>>,
        bound_violation: &mut f64,
        row_violation: &mut f64,
        integrality_violation: &mut f64,
        obj: &mut HighsCDouble,
    ) -> bool {
        *bound_violation = 0.0;
        *row_violation = 0.0;
        *integrality_violation = 0.0;
        let mip_feasibility_tolerance = self.options_mip_.mip_feasibility_tolerance;

        let lp = lp.unwrap_or(&self.orig_model_);

        *obj = HighsCDouble::from(lp.offset_);

        // HiGHS 1.16 #3149 (h16_cdouble): HighsMipSolver.cpp::solutionFeasible objective loop (1.16.0)
        let h16_cdouble_obj = self.options_mip_.h16_cdouble;
        for i in 0..lp.num_col_ as usize {
            let value = col_value[i];
            if h16_cdouble_obj {
                *obj += HighsCDouble::from(lp.col_cost_[i]) * value;
            } else {
                *obj += lp.col_cost_[i] * value;
            }

            if lp.integrality_[i] == HighsVarType::kInteger {
                *integrality_violation =
                    cpp_max(fractionality(value, None), *integrality_violation);
            }

            let lower = lp.col_lower_[i];
            let upper = lp.col_upper_[i];
            let primal_infeasibility = if value < lower - mip_feasibility_tolerance {
                lower - value
            } else if value > upper + mip_feasibility_tolerance {
                value - upper
            } else {
                continue;
            };

            *bound_violation = cpp_max(*bound_violation, primal_infeasibility);
        }

        // Check row feasibility if there are a positive number of rows.
        if lp.num_row_ > 0 {
            let mut row_value: Vec<f64> = Vec::new();
            if pass_row_value.is_none() {
                calculate_row_values_quad_vec(lp, col_value, &mut row_value, -1);
            }
            let row_value_p: &[f64] = match pass_row_value {
                Some(v) => v.as_slice(),
                None => &row_value,
            };

            for i in 0..lp.num_row_ as usize {
                let value = row_value_p[i];
                let lower = lp.row_lower_[i];
                let upper = lp.row_upper_[i];

                let primal_infeasibility = if value < lower - mip_feasibility_tolerance {
                    lower - value
                } else if value > upper + mip_feasibility_tolerance {
                    value - upper
                } else {
                    continue;
                };

                *row_violation = cpp_max(*row_violation, primal_infeasibility);
            }
        }

        *bound_violation <= mip_feasibility_tolerance
            && *integrality_violation <= mip_feasibility_tolerance
            && *row_violation <= mip_feasibility_tolerance
    }
}

#[cfg(test)]
mod h16_cdouble_obj_tests {
    //! HiGHS 1.16 #3149 (h16_cdouble) objective accumulation in `solutionFeasible`.

    use crate::lp_data::hconst::HighsVarType;
    use crate::prelude::*;

    /// Two columns: `x0 = 1 + 2^-30` with cost `1 + 2^-30`, and `x1 = 1` with cost
    /// `-(1 + 2^-29)`. The exact objective is `2^-60`. Rounding the `x0` product to a
    /// double first gives `1 + 2^-29`, which cancels exactly and loses the `2^-60` term.
    fn cancel_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 0;
        lp.offset_ = 0.0;
        lp.col_cost_ = vec![1.0 + 2f64.powi(-30), -(1.0 + 2f64.powi(-29))];
        lp.col_lower_ = vec![0.0; 2];
        lp.col_upper_ = vec![10.0; 2];
        lp.integrality_ = vec![HighsVarType::kContinuous; 2];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 0, 0];
        lp
    }

    /// `solution_feasible` objective (as `f64`) for `col_value` with `h16_cdouble` set.
    fn objective_of(h16_cdouble: bool, col_value: &[f64]) -> f64 {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.h16_cdouble = h16_cdouble;
        let lp = cancel_lp();
        let solution = HighsSolution::default();
        let mut callback = HighsCallback::default();
        let ms = HighsMipSolver::new(&mut callback, &options, &lp, &solution, false, 0);
        let mut bound_violation = 0.0;
        let mut row_violation = 0.0;
        let mut integrality_violation = 0.0;
        let mut obj = HighsCDouble::from(0.0);
        let feasible = ms.solution_feasible(
            None,
            col_value,
            None,
            &mut bound_violation,
            &mut row_violation,
            &mut integrality_violation,
            &mut obj,
        );
        assert!(feasible);
        f64::from(obj)
    }

    #[test]
    fn h16_cdouble_obj_on_keeps_low_product_bits() {
        let got = objective_of(true, &[1.0 + 2f64.powi(-30), 1.0]);
        assert_eq!(got.to_bits(), 2f64.powi(-60).to_bits());
    }

    #[test]
    fn h16_cdouble_obj_off_rounds_double_product() {
        let got = objective_of(false, &[1.0 + 2f64.powi(-30), 1.0]);
        assert_eq!(got.to_bits(), 0.0f64.to_bits());
    }
}
