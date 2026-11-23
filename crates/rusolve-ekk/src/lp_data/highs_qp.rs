//! HiGHS source: lp_data/Highs.cpp (Highs::callSolveQp, Highs::passModel(HighsModel),
//! Highs::passHessian), lp_data/HighsInterface.cpp (Highs::qFormatOk)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The `Highs` glue for QP. `Highs` keeps the HiGHS
//! `model_` as `model_lp_` + `model_hessian_`.
//!
//! * `call_solve_qp` builds the QUASS `Instance` from `model_lp_` /
//!   `model_hessian_`, fills `Settings`, subscribes the four log
//!   closures (each captures a clone of `options_.log_options` and of
//!   `options_.timeless_log`, as the C++ lambdas capture `this`), calls
//!   `qpsolver::a_quass::solveqp`, then `model_objective_value` /
//!   `get_kkt_failures_model` (model/highs_model.rs) and
//!   `check_optimality("QP")`. `use_hipo` reads the non-upstream switch
//!   `HighsOptions::hipo_available` for `isAvailable<HighsExtras::hipo>()`
//!   (default `false`, so `"ipm"`/`"hipo"` QPs go to QUASS as in the IPX and QP
//!   oracle builds; the HiPO tests set it and calls `crate::ipm::hipo_wrapper::solve_hipo`
//!   with the Hessian); the profiling calls are no-ops.
//! * `pass_highs_model` is the master `passModel(HighsModel)`. The port
//!   also rewires `Highs::pass_model(HighsLp)` (highs_run.rs) to
//!   `passModel(HighsLp)` (Highs.cpp 525-529), i.e. to call
//!   `pass_highs_model` with an empty Hessian, and adds the Hessian lines to
//!   `report_model_stats` (HighsInterface.cpp 25-125: `q_num_nz`, "QP" /
//!   "MIQP" problem type). Until then `pass_model` keeps its LP-only body.
//! * `analyseLp`, `writeLpMatrixPicToFile`, `writeHessianPicToFile`
//!   (analysis / image-file options, off by default) and `logHeader` are
//!   not ported.

use crate::lp_data::lp_assess::{assess_lp, lp_dimensions_ok};
use crate::prelude::*;
use crate::qpsolver::{
    Instance, PricingStrategy, QpAsmStatus, QpModelStatus, Settings, Statistics,
    qp_model_status_to_string, solveqp,
};
use crate::util::c_printf::{c_fmt_f, c_fmt_f_width, c_fmt_g_width};

/// rusolve (convex QP): which engine `call_solve_qp` ran last, `true` for HiPO. The
/// rusolve `lp_solver_auto` lever switches `solver` to `ipm` only around `optimize_highs` and
/// restores it, so the options after `run` do not say what ran. Read by the CLI's `qp.engine`.
pub static LAST_QP_ENGINE_HIPO: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

impl Highs {
    /// HiGHS: Highs.cpp::Highs::callSolveQp (lines 4139-4325)
    pub(crate) fn call_solve_qp(&mut self) -> HighsStatus {
        // Check that the model is column-wise
        debug_assert!(self.model_lp_.a_matrix_.is_colwise());
        debug_assert!(self.model_hessian_.format_ == HessianFormat::kTriangular);
        if self.model_hessian_.dim_ > self.model_lp_.num_col_ {
            crate::highs_log_dev!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Hessian dimension = {} is incompatible with matrix dimension = {}\n",
                self.model_hessian_.dim_,
                self.model_lp_.num_col_
            );
            self.model_status_ = HighsModelStatus::kModelError;
            self.solution_.value_valid = false;
            self.solution_.dual_valid = false;
            return HighsStatus::kError;
        }

        let return_status;

        // Choose solver
        let use_hipo = (self.options_.solver == kHipoString || self.options_.solver == kIpmString)
            && self.options_.hipo_available;

        LAST_QP_ENGINE_HIPO.store(use_hipo, std::sync::atomic::Ordering::Relaxed);
        if use_hipo {
            // profiling no-op
            return_status = crate::ipm::hipo_wrapper::solve_hipo(
                &self.options_,
                &mut self.timer_,
                &self.model_lp_,
                &self.model_hessian_,
                &mut self.basis_,
                &mut self.solution_,
                &mut self.model_status_,
                &mut self.info_,
                &mut self.callback_,
            );
            // profiling no-op
            if return_status == HighsStatus::kError {
                return return_status;
            }
        } else {
            //
            // Run the QP solver
            // profiling no-op

            let mut instance = Instance::new(self.model_lp_.num_col_, self.model_lp_.num_row_);

            instance.sense = self.model_lp_.sense_ as HighsInt;
            instance.num_con = self.model_lp_.num_row_;
            instance.num_var = self.model_lp_.num_col_;

            instance.A.mat.num_col = self.model_lp_.num_col_;
            instance.A.mat.num_row = self.model_lp_.num_row_;
            instance.A.mat.start = self.model_lp_.a_matrix_.start_.clone();
            instance.A.mat.index = self.model_lp_.a_matrix_.index_.clone();
            instance.A.mat.value = self.model_lp_.a_matrix_.value_.clone();
            instance.c.value = self.model_lp_.col_cost_.clone();
            instance.offset = self.model_lp_.offset_;
            instance.con_lo = self.model_lp_.row_lower_.clone();
            instance.con_up = self.model_lp_.row_upper_.clone();
            instance.var_lo = self.model_lp_.col_lower_.clone();
            instance.var_up = self.model_lp_.col_upper_.clone();
            instance.Q.mat.num_col = self.model_lp_.num_col_;
            instance.Q.mat.num_row = self.model_lp_.num_col_;
            triangular_to_square_hessian(
                &self.model_hessian_,
                &mut instance.Q.mat.start,
                &mut instance.Q.mat.index,
                &mut instance.Q.mat.value,
            );

            for i in 0..instance.c.value.len() as HighsInt {
                if instance.c.value[i as usize] != 0.0 {
                    instance.c.index[instance.c.num_nz as usize] = i;
                    instance.c.num_nz += 1;
                }
            }

            if self.model_lp_.sense_ == ObjSense::kMaximize {
                // Negate the vector and Hessian
                for i in instance.c.value.iter_mut() {
                    *i *= -1.0;
                }
                for i in instance.Q.mat.value.iter_mut() {
                    *i *= -1.0;
                }
            }

            let mut settings = Settings::default();
            let mut stats = Statistics::default();

            settings.reportingfequency = 100;
            if self.options_.qp_iteration_limit <= 10 {
                settings.reportingfequency = 1;
            } else if self.options_.qp_iteration_limit <= 100 {
                settings.reportingfequency = 10;
            }
            // Setting qp_update_limit = 10 leads to error with lpHighs3
            let qp_update_limit: HighsInt = 1000; // 1000; // default
            if qp_update_limit != settings.reinvertfrequency {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Changing QP reinversion frequency from {} to {}\n",
                    settings.reinvertfrequency,
                    qp_update_limit
                );
                settings.reinvertfrequency = qp_update_limit;
            }

            settings.allow_hot_start = self.options_.qp_allow_hot_start;
            settings.iteration_limit = self.options_.qp_iteration_limit;
            settings.nullspace_limit = self.options_.qp_nullspace_limit;
            // HiGHS 1.16 #3153 (h16_qp_bounded_start): flag carriers. The tolerance is copied only
            // when the flag is on (HiGHS 1.16 QpAsmWrapper.cpp:184 copies it unconditionally;
            // nothing reads it off-path).
            settings.h16_qp_bounded_start = self.options_.h16_qp_bounded_start;
            settings.h16_qp_convexity = self.options_.h16_qp_convexity;
            settings.h16_hfactor_pivot = self.options_.h16_hfactor_pivot;
            debug_assert!(settings.hessian_regularization_value == kHessianRegularizationValue);
            settings.hessian_regularization_value = self.options_.qp_regularization_value;
            if self.options_.h16_qp_bounded_start {
                settings.primal_feasibility_tolerance = self.options_.primal_feasibility_tolerance;
            }

            // Define the QP model status logging function
            let log_opts = self.options_.log_options.clone();
            settings
                .qp_model_status_log
                .subscribe(Rc::new(move |qp_model_status| {
                    if *qp_model_status == QpModelStatus::kUndetermined
                        || *qp_model_status == QpModelStatus::kLargeNullspace
                        || *qp_model_status == QpModelStatus::kNonConvex
                        || *qp_model_status == QpModelStatus::kError
                        || *qp_model_status == QpModelStatus::kNotset
                    {
                        crate::highs_log_user!(
                            &log_opts,
                            HighsLogType::kInfo,
                            "QP solver model status: {}\n",
                            qp_model_status_to_string(*qp_model_status)
                        );
                    }
                }));

            // Define the QP solver iteration logging function
            let log_opts = self.options_.log_options.clone();
            let timeless_log = self.options_.timeless_log;
            settings.iteration_log.subscribe(Rc::new(move |stats| {
                let rep = stats.iteration.len() - 1;
                let time_string = if timeless_log {
                    String::new()
                } else {
                    format!(" {}s", c_fmt_f_width(stats.time[rep], 2, 9))
                };
                crate::highs_log_user!(
                    &log_opts,
                    HighsLogType::kInfo,
                    "{:11}  {}           {:6}{}\n",
                    stats.iteration[rep],
                    c_fmt_g_width(stats.objval[rep], 8, 15),
                    stats.nullspacedimension[rep],
                    time_string
                );
            }));

            // Define the QP nullspace limit logging function
            let log_opts = self.options_.log_options.clone();
            settings
                .nullspace_limit_log
                .subscribe(Rc::new(move |nullspace_limit| {
                    crate::highs_log_user!(
                        &log_opts,
                        HighsLogType::kError,
                        "QP solver has exceeded nullspace limit of {}\n",
                        *nullspace_limit
                    );
                }));

            // Define the degeneracy failure logging function
            let log_opts = self.options_.log_options.clone();
            settings
                .degeneracy_fail_log
                .subscribe(Rc::new(move |degeneracy_fail_data| {
                    crate::highs_log_user!(
                        &log_opts,
                        HighsLogType::kError,
                        "QP solver has failed due to degeneracy: \
                         cannot find non-active constraint to leave basis. \
                         max: log(d[{}]) = {}\n",
                        degeneracy_fail_data.0,
                        c_fmt_f(degeneracy_fail_data.1, 6)
                    );
                }));

            settings.time_limit = self.options_.time_limit;
            settings.lambda_zero_threshold = self.options_.dual_feasibility_tolerance;

            match self.options_.simplex_primal_edge_weight_strategy {
                0 => {
                    settings.pricing = PricingStrategy::DantzigWolfe;
                }
                1 => {
                    settings.pricing = PricingStrategy::Devex;
                }
                2 => {
                    settings.pricing = PricingStrategy::SteepestEdge;
                }
                _ => {
                    settings.pricing = PricingStrategy::Devex;
                }
            }

            // print header for QP solver output
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "  Iteration        Objective     NullspaceDim\n"
            );

            let status = solveqp(
                &mut instance,
                &mut settings,
                &mut stats,
                &mut self.model_status_,
                &mut self.basis_,
                &mut self.solution_,
                &mut self.timer_,
            );
            // profiling no-op

            // QP solver can fail, so should return something other than
            // QpAsmStatus::kOk
            if status == QpAsmStatus::kError {
                return HighsStatus::kError;
            }

            debug_assert!(status == QpAsmStatus::kOk || status == QpAsmStatus::kWarning);
            return_status = if status == QpAsmStatus::kWarning {
                HighsStatus::kWarning
            } else {
                HighsStatus::kOk
            };

            // Set the QP-specific values of info_
            self.info_.simplex_iteration_count += stats.phase1_iterations;
            self.info_.qp_iteration_count += stats.num_iterations;
        }

        // Get the objective and any KKT failures
        self.info_.objective_function_value = model_objective_value(
            &self.model_lp_,
            &self.model_hessian_,
            &self.solution_.col_value,
            self.options_.h16_obj_cdouble,
        );
        get_kkt_failures_model(
            &self.options_,
            &self.model_lp_,
            &self.model_hessian_,
            &self.solution_,
            &self.basis_,
            &mut self.info_,
        );
        self.info_.valid = true;
        if self.model_status_ == HighsModelStatus::kOptimal {
            return self.check_optimality("QP");
        }
        return_status
    }

    /// HiGHS: Highs.cpp::Highs::passModel(HighsModel) (lines 443-523)
    pub fn pass_highs_model(&mut self, model: HighsModel) -> HighsStatus {
        // This is the "master" Highs::passModel, in that all the others
        // (and readModel) eventually call it
        // logHeader not ported.
        // Possibly analyse the LP data
        // analyseLp not ported.
        let mut return_status = HighsStatus::kOk;
        // Clear the incumbent model and any associated data
        self.clear_model();
        // Move the model's LP and Hessian to the internal LP and Hessian
        self.model_lp_ = model.lp_;
        self.model_hessian_ = model.hessian_;
        self.model_lp_.origin_name_ = "Original".to_string();
        debug_assert!(self.model_lp_.a_matrix_.format_ok());
        if self.model_lp_.num_col_ == 0 || self.model_lp_.num_row_ == 0 {
            // Model constraint matrix has either no columns or no
            // rows. Clearly the matrix is empty, so may have no orientation
            // or starts assigned. HiGHS assumes that such a model will have
            // null starts, so make it column-wise
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Model has either no columns or no rows, so ignoring user \
                 constraint matrix data and initialising empty matrix\n"
            );
            self.model_lp_.a_matrix_.format_ = MatrixFormat::kColwise;
            self.model_lp_.a_matrix_.start_.clear();
            self.model_lp_
                .a_matrix_
                .start_
                .resize((self.model_lp_.num_col_ + 1) as usize, 0);
            self.model_lp_.a_matrix_.index_.clear();
            self.model_lp_.a_matrix_.value_.clear();
        } else {
            // Matrix has rows and columns, so a_matrix format must be valid
            if !self.model_lp_.a_matrix_.format_ok() {
                return HighsStatus::kError;
            }
        }
        // Dimensions in a_matrix_ may not be set, so take them from lp.
        self.model_lp_.set_matrix_dimensions();
        // Residual scale factors may be present. ToDo Allow user-defined
        // scale factors to be passed
        debug_assert!(!self.model_lp_.is_scaled_);
        debug_assert!(!self.model_lp_.is_moved_);
        self.model_lp_.reset_scale();
        // Check that the LP array dimensions are valid
        if !lp_dimensions_ok("passModel", &self.model_lp_, &self.options_.log_options) {
            return HighsStatus::kError;
        }
        // Check that the Hessian format is valid
        if !self.model_hessian_.format_ok() {
            return HighsStatus::kError;
        }
        // Check validity of the LP, normalising its values
        return_status = interpret_call_status(
            &self.options_.log_options,
            assess_lp(&mut self.model_lp_, &self.options_),
            return_status,
            "assessLp",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        // Now legality of matrix is established, ensure that it is
        // column-wise
        self.model_lp_.ensure_colwise();
        // Check validity of any Hessian, normalising its entries
        return_status = interpret_call_status(
            &self.options_.log_options,
            assess_hessian(&mut self.model_hessian_, &self.options_),
            return_status,
            "assessHessian",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        if self.model_hessian_.dim_ != 0 {
            // Clear any zero Hessian
            if self.model_hessian_.num_nz() == 0 {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Hessian has dimension {} but no nonzeros, so is ignored\n",
                    self.model_hessian_.dim_
                );
                self.model_hessian_.clear();
            }
        }
        // Ensure that any non-zero Hessian of dimension less than the
        // number of columns in the model is completed
        if self.model_hessian_.dim_ != 0 {
            complete_hessian(self.model_lp_.num_col_, &mut self.model_hessian_);
        }
        // Possibly create image files of the constraint matrix and Hessian
        // writeLpMatrixPicToFile / writeHessianPicToFile not ported.

        // Clear solver status, solution, basis and info associated with any
        // previous model; clear any HiGHS model object; create a HiGHS
        // model object for this LP
        let clear_solver_status = self.clear_solver();
        return_status = interpret_call_status(
            &self.options_.log_options,
            clear_solver_status,
            return_status,
            "clearSolver",
        );
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::passHessian(HighsHessian) (lines 655-682)
    pub fn pass_hessian(&mut self, hessian: HighsHessian) -> HighsStatus {
        // logHeader not ported.
        let mut return_status = HighsStatus::kOk;
        self.model_hessian_ = hessian;
        // Check validity of any Hessian, normalising its entries
        return_status = interpret_call_status(
            &self.options_.log_options,
            assess_hessian(&mut self.model_hessian_, &self.options_),
            return_status,
            "assessHessian",
        );
        if return_status == HighsStatus::kError {
            return return_status;
        }
        if self.model_hessian_.dim_ != 0 {
            // Clear any zero Hessian
            if self.model_hessian_.num_nz() == 0 {
                crate::highs_log_user!(
                    &self.options_.log_options,
                    HighsLogType::kInfo,
                    "Hessian has dimension {} but no nonzeros, so is ignored\n",
                    self.model_hessian_.dim_
                );
                self.model_hessian_.clear();
            }
        }
        // Ensure that any non-zero Hessian of dimension less than the
        // number of columns in the model is completed
        if self.model_hessian_.dim_ != 0 {
            complete_hessian(self.model_lp_.num_col_, &mut self.model_hessian_);
        }

        let clear_solver_status = self.clear_solver();
        return_status = interpret_call_status(
            &self.options_.log_options,
            clear_solver_status,
            return_status,
            "clearSolver",
        );
        self.return_from_highs(return_status)
    }

    /// HiGHS: Highs.cpp::Highs::passHessian(dim, num_nz, format, start, index,
    /// value) (lines 684-711)
    ///
    /// `start` has `dim` entries (no end sentinel), as in C++.
    pub fn pass_hessian_raw(
        &mut self,
        dim: HighsInt,
        num_nz: HighsInt,
        format: HighsInt,
        start: &[HighsInt],
        index: &[HighsInt],
        value: &[f64],
    ) -> HighsStatus {
        // logHeader not ported.
        let mut hessian = HighsHessian::default();
        if !self.q_format_ok(num_nz, format) {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Model has illegal Hessian matrix format\n"
            );
            return HighsStatus::kError;
        }
        let num_col = self.model_lp_.num_col_;
        if dim != num_col {
            return HighsStatus::kError;
        }
        hessian.dim_ = num_col;
        // C++ static_cast<HessianFormat>(format); only the two valid values
        // are representable here (qFormatOk rejects others when num_nz > 0).
        hessian.format_ = if format == HessianFormat::kSquare as HighsInt {
            HessianFormat::kSquare
        } else {
            HessianFormat::kTriangular
        };
        if dim > 0 {
            debug_assert!(!start.is_empty());
            hessian.start_ = start[..num_col as usize].to_vec();
            hessian.start_.resize((num_col + 1) as usize, 0);
            hessian.start_[num_col as usize] = num_nz;
        }
        if num_nz > 0 {
            debug_assert!(!index.is_empty());
            debug_assert!(!value.is_empty());
            hessian.index_ = index[..num_nz as usize].to_vec();
            hessian.value_ = value[..num_nz as usize].to_vec();
        }
        self.pass_hessian(hessian)
    }

    /// HiGHS: HighsInterface.cpp::Highs::qFormatOk (lines 2875-2886)
    pub(crate) fn q_format_ok(&self, num_nz: HighsInt, format: HighsInt) -> bool {
        if num_nz == 0 {
            return true;
        }
        let ok_format = format == HessianFormat::kTriangular as HighsInt
            || format == HessianFormat::kSquare as HighsInt;
        if !ok_format {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kError,
                "Non-empty Hessian matrix has illegal format = {}\n",
                format
            );
        }
        debug_assert!(ok_format);
        ok_format
    }
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    /// Column-wise LP builder: `cols[j] = (cost, lower, upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    fn quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    #[test]
    fn highs_runs_box_qp() {
        let lp = lp(
            &[
                (-2.0, 0.0, 10.0, &[(0, 1.0)]),
                (-4.0, 0.0, 10.0, &[(0, 1.0)]),
            ],
            &[(-kHighsInf, 10.0)],
        );
        let mut highs = quiet(lp);
        let hessian = HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 1, 2],
            index_: vec![0, 1],
            value_: vec![2.0, 4.0],
        };
        assert_eq!(highs.pass_hessian(hessian), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        // The default qp_regularization_value (1e-7) is added to the Hessian
        // diagonal, so the optimum is x_j = c_j / (h_j + reg) exactly, i.e.
        // 2/(2+1e-7) and 4/(4+1e-7), not 1.0 (the C++ result, not a bug).
        let reg = 1e-7;
        let x = &highs.get_solution().col_value;
        assert!((x[0] - 2.0 / (2.0 + reg)).abs() <= 1e-12);
        assert!((x[1] - 4.0 / (4.0 + reg)).abs() <= 1e-12);
        assert!((highs.get_info().objective_function_value + 3.0).abs() <= 1e-6);
    }

    #[test]
    fn non_convex_qp_is_rejected() {
        let lp = lp(&[(0.0, 0.0, 1.0, &[]), (0.0, 0.0, 1.0, &[])], &[]);
        let mut highs = quiet(lp);
        let hessian = HighsHessian {
            dim_: 2,
            format_: HessianFormat::kTriangular,
            start_: vec![0, 1, 2],
            index_: vec![0, 1],
            value_: vec![-1.0, 1.0],
        };
        assert_eq!(highs.pass_hessian(hessian), HighsStatus::kOk);
        assert_eq!(highs.run(), HighsStatus::kError);
    }
}
