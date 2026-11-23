//! HiGHS source: ipm/ipx/lp_solver.h, lp_solver.cc (class LpSolver)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `control_` and `model_` are `Rc` (see "Ownership" in `ipx/mod.rs`): `Iterate`, `Basis`, the KKT
//! solvers and `Crossover` keep clones. `LoadModel` and `ClearModel` mutate
//! `model_` through `Rc::get_mut`; both run `ClearModel`, which drops
//! `iterate_`/`basis_` (the only clones that outlive `Solve`) first, so the
//! `expect` cannot fire. `SetParameters` and `SetCallback` go through the
//! `RefCell`s inside `Control` and work at any time, as in C++.
//!
//! `SetCallback(HighsCallback*)` moves the callback in; the wrapper gets it
//! back with [`LpSolver::take_callback`] after `Solve`.

use super::basis::Basis;
use super::control::{Control, fix2, sci2, sci8, textline};
use super::crossover::Crossover;
use super::info::{status_string, write_info};
use super::ipm::IPM;
use super::iterate::{Iterate, State};
use super::kkt_solver_basis::KKTSolverBasis;
use super::kkt_solver_diag::KKTSolverDiag;
use super::model::{Model, dual_infeasibility, primal_infeasibility};
use super::timer::Timer;
use super::*;
use crate::lp_data::highs_callback::HighsCallback;
use std::rc::Rc;

/// HiGHS: ipx/lp_solver.h::LpSolver (lines 17-207)
pub struct LpSolver {
    control_: Rc<Control>,
    info_: Info,
    model_: Rc<Model>,
    iterate_: Option<Box<Iterate>>,
    basis_: Option<Box<Basis>>,

    // Basic solution computed by crossover and basic status of each variable
    // (one of IPX_nonbasic_lb, IPX_nonbasic_ub, IPX_basic, IPX_superbasic).
    x_crossover_: Vector,
    y_crossover_: Vector,
    z_crossover_: Vector,
    crossover_weights_: Vector,
    basic_statuses_: Vec<Int>,

    // IPM starting point provided by user (presolved).
    x_start_: Vector,
    xl_start_: Vector,
    xu_start_: Vector,
    y_start_: Vector,
    zl_start_: Vector,
    zu_start_: Vector,
}

impl Default for LpSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl LpSolver {
    /// HiGHS: ipx/lp_solver.h::LpSolver (implicit default constructor)
    pub fn new() -> Self {
        Self {
            control_: Rc::new(Control::new()),
            info_: Info::default(),
            model_: Rc::new(Model::new()),
            iterate_: None,
            basis_: None,
            x_crossover_: Vector::default(),
            y_crossover_: Vector::default(),
            z_crossover_: Vector::default(),
            crossover_weights_: Vector::default(),
            basic_statuses_: Vec::new(),
            x_start_: Vector::default(),
            xl_start_: Vector::default(),
            xu_start_: Vector::default(),
            y_start_: Vector::default(),
            zl_start_: Vector::default(),
            zu_start_: Vector::default(),
        }
    }

    /// HiGHS: ipx/lp_solver.cc::LpSolver::LoadModel (lines 15-25)
    ///
    /// Nullable C pointers are `Option` (`Model::load`).
    #[allow(clippy::too_many_arguments)]
    pub fn load_model(
        &mut self,
        num_var: Int,
        offset: f64,
        obj: Option<&[f64]>,
        lb: Option<&[f64]>,
        ub: Option<&[f64]>,
        num_constr: Int,
        Ap: Option<&[Int]>,
        Ai: Option<&[Int]>,
        Ax: Option<&[f64]>,
        rhs: Option<&[f64]>,
        constr_type: Option<&[u8]>,
    ) -> Int {
        self.clear_model();
        let errflag = Rc::get_mut(&mut self.model_)
            .expect("LpSolver::load_model: model_ is shared")
            .load(
                &self.control_,
                num_constr,
                num_var,
                Ap,
                Ai,
                Ax,
                rhs,
                constr_type,
                offset,
                obj,
                lb,
                ub,
            );
        self.model_.get_info(&mut self.info_);
        errflag
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::LoadIPMStartingPoint (lines 27-49)
    #[allow(clippy::too_many_arguments)]
    pub fn load_ipm_starting_point(
        &mut self,
        x: Option<&[f64]>,
        xl: Option<&[f64]>,
        xu: Option<&[f64]>,
        slack: Option<&[f64]>,
        y: Option<&[f64]>,
        zl: Option<&[f64]>,
        zu: Option<&[f64]>,
    ) -> Int {
        let m = self.model_.rows();
        let n = self.model_.cols();
        self.x_start_.resize(n + m);
        self.xl_start_.resize(n + m);
        self.xu_start_.resize(n + m);
        self.y_start_.resize(m);
        self.zl_start_.resize(n + m);
        self.zu_start_.resize(n + m);
        let errflag = self.model_.presolve_ipm_starting_point(
            x,
            xl,
            xu,
            slack,
            y,
            zl,
            zu,
            &mut self.x_start_,
            &mut self.xl_start_,
            &mut self.xu_start_,
            &mut self.y_start_,
            &mut self.zl_start_,
            &mut self.zu_start_,
        );
        if errflag != 0 {
            self.clear_ipm_starting_point();
            return errflag;
        }
        self.make_ipm_starting_point_valid();
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::Solve (lines 51-126)
    ///
    /// PORT NOTE: the C++ `try { ... } catch (const std::bad_alloc&)` and
    /// `catch (const std::exception&)` blocks are not reproduced: Rust has no
    /// `bad_alloc` to catch (an allocation failure aborts) and any internal
    /// panic propagates to the caller. Consequently `IPX_STATUS_out_of_memory`
    /// and `IPX_STATUS_internal_error` are never set here.
    pub fn solve(&mut self) -> Int {
        if self.model_.empty() {
            self.info_.status = IPX_STATUS_no_model;
            return self.info_.status;
        }
        self.clear_solution();
        self.control_.reset_timer();
        self.control_.open_logfile();
        self.control_.h_log_str("IPX version 1.0\n");

        self.interior_point_solve();
        let run_crossover_on = self.control_.run_crossover() == 1;
        let run_crossover_choose = self.control_.run_crossover() == -1;
        let run_crossover_not_off = run_crossover_choose || run_crossover_on;
        let run_crossover = (self.info_.status_ipm == IPX_STATUS_optimal && run_crossover_on)
            || (self.info_.status_ipm == IPX_STATUS_imprecise && run_crossover_not_off);
        //        if ((info_.status_ipm == IPX_STATUS_optimal ||
        //             info_.status_ipm == IPX_STATUS_imprecise) && run_crossover_on) {
        if run_crossover {
            if run_crossover_on {
                self.control_.h_log_str("Running crossover as requested\n");
            } else if run_crossover_choose {
                debug_assert!(self.info_.status_ipm == IPX_STATUS_imprecise);
                self.control_
                    .h_log_str("Running crossover since IPX is imprecise\n");
            } else {
                debug_assert!(run_crossover_on || run_crossover_choose);
            }
            self.build_crossover_starting_point();
            self.run_crossover();
        }
        if let Some(basis) = self.basis_.as_deref() {
            self.info_.ftran_sparse = basis.frac_ftran_sparse();
            self.info_.btran_sparse = basis.frac_btran_sparse();
            self.info_.time_lu_invert = basis.time_factorize();
            self.info_.time_lu_update = basis.time_update();
            self.info_.time_ftran = basis.time_ftran();
            self.info_.time_btran = basis.time_btran();
            self.info_.mean_fill = basis.mean_fill();
            self.info_.max_fill = basis.max_fill();
        }
        if self.info_.status_ipm == IPX_STATUS_primal_infeas
            || self.info_.status_ipm == IPX_STATUS_dual_infeas
            || self.info_.status_crossover == IPX_STATUS_primal_infeas
            || self.info_.status_crossover == IPX_STATUS_dual_infeas
        {
            // When IPM or crossover detect the model to be infeasible
            // (currently only the former is implemented), then the problem is
            // solved.
            self.info_.status = IPX_STATUS_solved;
        } else {
            let method_status = if run_crossover {
                self.info_.status_crossover
            } else {
                self.info_.status_ipm
            };
            if method_status == IPX_STATUS_optimal || method_status == IPX_STATUS_imprecise {
                self.info_.status = IPX_STATUS_solved;
            } else {
                self.info_.status = IPX_STATUS_stopped;
            }
        }
        self.print_summary();

        self.info_.time_total = self.control_.elapsed();
        super::oracle_trace::ipx_info("solve_end", &self.info_);
        let mut info_text = String::new();
        write_info(&mut info_text, &self.info_);
        self.control_.debug(2, &info_text);
        self.control_.close_logfile();
        if self.control_.report_basis_data() {
            // PORT NOTE: C++ dereferences `basis_` unconditionally here; it is
            // null when BuildStartingBasis returned before creating it, so the
            // Rust port skips the report in that case.
            if let Some(basis) = self.basis_.as_deref() {
                basis.report_basis_data();
            }
        }
        self.info_.status
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetInfo (lines 128-130)
    pub fn get_info(&self) -> Info {
        self.info_
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetInteriorSolution (lines 132-142)
    ///
    /// Each `None` output is not returned.
    #[allow(clippy::too_many_arguments)]
    pub fn get_interior_solution(
        &self,
        x: Option<&mut [f64]>,
        xl: Option<&mut [f64]>,
        xu: Option<&mut [f64]>,
        slack: Option<&mut [f64]>,
        y: Option<&mut [f64]>,
        zl: Option<&mut [f64]>,
        zu: Option<&mut [f64]>,
    ) -> Int {
        let Some(iterate) = self.iterate_.as_deref() else {
            return -1;
        };
        self.model_.postsolve_interior_solution(
            iterate.x(),
            iterate.xl(),
            iterate.xu(),
            iterate.y(),
            iterate.zl(),
            iterate.zu(),
            x,
            xl,
            xu,
            slack,
            y,
            zl,
            zu,
        );
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetBasicSolution (lines 144-152)
    ///
    /// Each `None` output is not returned.
    pub fn get_basic_solution(
        &self,
        x: Option<&mut [f64]>,
        slack: Option<&mut [f64]>,
        y: Option<&mut [f64]>,
        z: Option<&mut [f64]>,
        cbasis: Option<&mut [Int]>,
        vbasis: Option<&mut [Int]>,
    ) -> Int {
        if self.basic_statuses_.is_empty() {
            return -1;
        }
        self.model_.postsolve_basic_solution(
            &self.x_crossover_,
            &self.y_crossover_,
            &self.z_crossover_,
            &self.basic_statuses_,
            x,
            slack,
            y,
            z,
        );
        self.model_
            .postsolve_basis(&self.basic_statuses_, cbasis, vbasis);
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetParameters (lines 154-156)
    pub fn get_parameters(&self) -> Parameters {
        self.control_.parameters().clone()
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::SetParameters (lines 158-160)
    ///
    /// `Control::parameters_` is a `RefCell`, so this works while `basis_`/
    /// `iterate_` share `control_` (after `Solve`), as in C++.
    pub fn set_parameters(&mut self, new_parameters: Parameters) {
        self.control_.set_parameters(&new_parameters);
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::SetCallback (lines 162-164)
    ///
    /// PORT NOTE: moves the callback into `control_`; see
    /// [`LpSolver::take_callback`].
    pub fn set_callback(&mut self, callback: Option<HighsCallback>) {
        self.control_.callback(callback);
    }
    /// PORT NOTE: returns the callback moved in by
    /// [`LpSolver::set_callback`] (no C++ counterpart; the C++ wrapper keeps
    /// ownership of `HighsCallback` and IPX only holds a pointer).
    pub fn take_callback(&mut self) -> Option<HighsCallback> {
        self.control_.take_callback()
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::ClearModel (lines 166-170)
    ///
    /// PORT NOTE: C++ calls `model_.clear()` before `ClearSolution()`. Here
    /// `iterate_`/`basis_` hold `Rc` clones of `model_`, so they are dropped
    /// first; `ClearSolution()` drops them anyway, and nothing in their
    /// destructors touches the model, so the observable order is unchanged.
    pub fn clear_model(&mut self) {
        self.iterate_ = None;
        self.basis_ = None;
        Rc::get_mut(&mut self.model_)
            .expect("LpSolver::clear_model: model_ is shared")
            .clear();
        self.clear_solution();
        self.clear_ipm_starting_point();
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::ClearIPMStartingPoint (lines 172-179)
    pub fn clear_ipm_starting_point(&mut self) {
        self.x_start_.resize(0);
        self.xl_start_.resize(0);
        self.xu_start_.resize(0);
        self.y_start_.resize(0);
        self.zl_start_.resize(0);
        self.zu_start_.resize(0);
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::CrossoverFromStartingPoint (lines 181-245)
    ///
    /// Each `None` vector is taken as all zeros.
    #[allow(clippy::needless_range_loop)]
    pub fn crossover_from_starting_point(
        &mut self,
        x_start: Option<&[f64]>,
        slack_start: Option<&[f64]>,
        y_start: Option<&[f64]>,
        z_start: Option<&[f64]>,
    ) -> Int {
        let m = self.model_.rows();
        let n = self.model_.cols();

        self.clear_solution();
        self.control_.h_log_str("Crossover from starting point\n");

        self.x_crossover_.resize(n + m);
        self.y_crossover_.resize(m);
        self.z_crossover_.resize(n + m);
        self.crossover_weights_.resize(0);
        self.model_.presolve_starting_point(
            x_start,
            slack_start,
            y_start,
            z_start,
            &mut self.x_crossover_,
            &mut self.y_crossover_,
            &mut self.z_crossover_,
        );

        // Check that starting point is complementary and satisfies bound and
        // sign conditions.
        {
            let lb = self.model_.lb();
            let ub = self.model_.ub();
            let x = &self.x_crossover_;
            let z = &self.z_crossover_;
            for j in 0..(n + m) as usize {
                if x[j] < lb[j] || x[j] > ub[j] {
                    return IPX_ERROR_invalid_vector;
                }
                if x[j] != lb[j] && z[j] > 0.0 {
                    return IPX_ERROR_invalid_vector;
                }
                if x[j] != ub[j] && z[j] < 0.0 {
                    return IPX_ERROR_invalid_vector;
                }
            }
        }

        // Construct starting basis.
        self.basis_ = Some(Box::new(Basis::new(
            self.control_.clone(),
            self.model_.clone(),
        )));
        if self.control_.crash_basis() != 0 {
            // Take columns in the following order of priority:
            // - free columns
            // - columns between their bounds, in increasing number of nonzeros
            // - columns with zero dual, in increasing number of nonzeros
            // - Fixed columns and those with nonzero dual
            let timer = Timer::new(0.0);
            let ipm_crash = self.control_.crossover_ipm_crash();
            let colweight = {
                let lb = self.model_.lb();
                let ub = self.model_.ub();
                let ai = self.model_.AI();
                let x = &self.x_crossover_;
                let z = &self.z_crossover_;
                let mut colweight = Vector::new(n + m);
                for j in 0..(n + m) as usize {
                    let nz = ai.entries_col(j as Int);
                    colweight[j] = if lb[j] == ub[j] {
                        0.0
                    } else if lb[j].is_infinite() && ub[j].is_infinite() {
                        f64::INFINITY
                    } else if z[j] != 0.0 {
                        0.0
                    } else if x[j] != lb[j] && x[j] != ub[j] {
                        if ipm_crash {
                            // rusolve lp_crossover_ipm_crash (Bixby and
                            // Saltzman 1994: basic candidates are the variables farthest from
                            // their bounds): rank by distance to the nearest finite bound,
                            // then by sparsity
                            let d = (x[j] - lb[j]).min(ub[j] - x[j]);
                            let t = d / (1.0 + d);
                            (m + 1) as f64
                                + t * m as f64
                                + 1e-9 * (m - nz + 1) as f64 / (m + 1) as f64
                        } else {
                            (m + (m - nz + 1)) as f64
                        }
                    } else {
                        (m - nz + 1) as f64
                    };
                }
                colweight
            };
            self.basis_
                .as_deref_mut()
                .expect("basis_ was just constructed")
                .construct_basis_from_weights(&colweight, &mut self.info_);
            self.info_.time_starting_basis += timer.elapsed();
            if self.info_.errflag == IPX_ERROR_time_interrupt {
                // rusolve BASICLU time limit
                self.clear_solution();
                self.info_.status_crossover = IPX_STATUS_time_limit;
                return 0;
            }
            if self.info_.errflag != 0 {
                self.clear_solution();
                return 0;
            }
        }

        self.run_crossover();
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetIterate (lines 247-264)
    #[allow(clippy::too_many_arguments)]
    pub fn get_iterate(
        &mut self,
        x: Option<&mut [f64]>,
        y: Option<&mut [f64]>,
        zl: Option<&mut [f64]>,
        zu: Option<&mut [f64]>,
        xl: Option<&mut [f64]>,
        xu: Option<&mut [f64]>,
    ) -> Int {
        let Some(iterate) = self.iterate_.as_deref() else {
            return -1;
        };
        if let Some(out) = x {
            copy_out_f64(iterate.x(), out);
        }
        if let Some(out) = y {
            copy_out_f64(iterate.y(), out);
        }
        if let Some(out) = zl {
            copy_out_f64(iterate.zl(), out);
        }
        if let Some(out) = zu {
            copy_out_f64(iterate.zu(), out);
        }
        if let Some(out) = xl {
            copy_out_f64(iterate.xl(), out);
        }
        if let Some(out) = xu {
            copy_out_f64(iterate.xu(), out);
        }
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetBasis (lines 289-299)
    pub fn get_basis(&mut self, cbasis: Option<&mut [Int]>, vbasis: Option<&mut [Int]>) -> Int {
        let Some(basis) = self.basis_.as_deref() else {
            return -1;
        };
        if !self.basic_statuses_.is_empty() {
            // crossover provides basic statuses
            self.model_
                .postsolve_basis(&self.basic_statuses_, cbasis, vbasis);
        } else {
            let statuses = build_basic_statuses(basis);
            self.model_.postsolve_basis(&statuses, cbasis, vbasis);
        }
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::GetKKTMatrix (lines 301-333)
    ///
    /// The matrix is returned only if all of `AIp`, `AIi`, `AIx` are `Some`.
    pub fn get_kkt_matrix(
        &mut self,
        AIp: Option<&mut [Int]>,
        AIi: Option<&mut [Int]>,
        AIx: Option<&mut [f64]>,
        g: Option<&mut [f64]>,
    ) -> Int {
        let Some(iterate) = self.iterate_.as_deref() else {
            return -1;
        };
        if let (Some(aip), Some(aii), Some(aix)) = (AIp, AIi, AIx) {
            let ai = self.model_.AI();
            copy_out_int(&ai.colptr()[..(ai.cols() + 1) as usize], aip);
            copy_out_int(&ai.rowidx()[..ai.entries() as usize], aii);
            copy_out_f64(&ai.values()[..ai.entries() as usize], aix);
        }
        if let Some(g) = g {
            let m = self.model_.rows();
            let n = self.model_.cols();
            for j in 0..(n + m) as usize {
                g[j] = match iterate.state_of(j as Int) {
                    State::fixed => f64::INFINITY,
                    State::free => 0.0,
                    State::barrier => {
                        let value = iterate.zl_at(j as Int) / iterate.xl_at(j as Int)
                            + iterate.zu_at(j as Int) / iterate.xu_at(j as Int);
                        debug_assert!(value.is_finite());
                        debug_assert!(value > 0.0);
                        value
                    }
                };
            }
        }
        0
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::SymbolicInvert (lines 335-340)
    pub fn symbolic_invert(
        &mut self,
        rowcounts: Option<&mut [Int]>,
        colcounts: Option<&mut [Int]>,
    ) -> Int {
        let Some(basis) = self.basis_.as_deref() else {
            return -1;
        };
        basis.symbolic_invert(rowcounts, colcounts);
        0
    }
    /// HiGHS: ipx/lp_solver.h::LpSolver::setTimerOffset (line 167)
    pub fn set_timer_offset(&self, offset: f64) {
        self.control_.set_timer_offset(offset);
    }

    // private

    /// HiGHS: ipx/lp_solver.cc::LpSolver::ClearSolution (lines 342-355)
    fn clear_solution(&mut self) {
        self.iterate_ = None;
        self.basis_ = None;
        self.x_crossover_.resize(0);
        self.y_crossover_.resize(0);
        self.z_crossover_.resize(0);
        self.crossover_weights_.resize(0);
        self.basic_statuses_.clear();
        self.basic_statuses_.shrink_to_fit();
        self.info_ = Info::default();
        // Restore info entries that belong to model.
        self.model_.get_info(&mut self.info_);
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::InteriorPointSolve (lines 357-388)
    fn interior_point_solve(&mut self) {
        if self.control_.run_centring() != 0 {
            self.control_
                .h_log_str("Interior point solve for analytic centre\n");
        } else {
            self.control_.h_log_str("Interior point solve\n");
        }

        // Allocate new iterate and set tolerances for IPM termination test.
        let mut iterate = Iterate::new(self.model_.clone());
        iterate.set_feasibility_tol(self.control_.ipm_feasibility_tol());
        iterate.set_optimality_tol(self.control_.ipm_optimality_tol());
        // rusolve (HiGHS 1.16 backport, default off): passed to the IPX iterate (no IPX code reads it yet)
        iterate.set_h16_gap_offset(self.control_.h16_ipm_gap_offset());
        if self.control_.run_crossover() != 0 {
            iterate.set_start_crossover_tol(self.control_.start_crossover_tol());
        }
        self.iterate_ = Some(Box::new(iterate));

        self.run_ipm();

        let iterate = self
            .iterate_
            .as_deref_mut()
            .expect("iterate_ is set by run_ipm's caller");
        iterate.postprocess();
        iterate.evaluate_postsolved(&mut self.info_);

        // Declare status_ipm "imprecise" if the IPM terminated optimal but the
        // solution after postprocessing/postsolve does not satisfy tolerances.
        if self.info_.status_ipm == IPX_STATUS_optimal
            && (self.info_.rel_objgap.abs() > self.control_.ipm_optimality_tol()
                || self.info_.rel_presidual > self.control_.ipm_feasibility_tol()
                || self.info_.rel_dresidual > self.control_.ipm_feasibility_tol())
        {
            self.info_.status_ipm = IPX_STATUS_imprecise;
        }

        // rusolve: an IPM that stalls within a loose distance of
        // the optimum is "imprecise", so crossover runs on it (Andersen & Andersen 2000 start
        // crossover from an approximate optimum) instead of a cold simplex solve.
        if self.control_.stall_to_crossover()
            && self.info_.status_ipm == IPX_STATUS_no_progress
            && self.info_.rel_objgap.abs() <= 1e-5
            && self.info_.rel_presidual <= 10.0 * self.control_.ipm_feasibility_tol()
            && self.info_.rel_dresidual <= 10.0 * self.control_.ipm_feasibility_tol()
        {
            self.info_.status_ipm = IPX_STATUS_imprecise;
        }

        // Assess the success of analytic centre calculation
        if self.info_.centring_tried != 0 {
            self.info_.status_ipm = if self.info_.centring_success != 0 {
                IPX_STATUS_optimal
            } else {
                IPX_STATUS_imprecise
            };
        }
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::RunIPM (lines 390-413)
    fn run_ipm(&mut self) {
        let mut ipm = IPM::new(self.control_.clone());
        self.info_.centring_tried = 0;
        self.info_.centring_success = 0;

        if self.x_start_.size() != 0 {
            self.control_.h_log_str(
                " Using starting point provided by user. Skipping initial iterations.\n",
            );
            self.iterate_
                .as_deref_mut()
                .expect("iterate_ is set by interior_point_solve")
                .initialize(
                    &self.x_start_,
                    &self.xl_start_,
                    &self.xu_start_,
                    &self.y_start_,
                    &self.zl_start_,
                    &self.zu_start_,
                );
        } else {
            self.compute_starting_point(&mut ipm);
            if self.info_.status_ipm != IPX_STATUS_not_run {
                return;
            }
            self.run_initial_ipm(&mut ipm);
            if self.info_.status_ipm != IPX_STATUS_not_run {
                return;
            }
        }
        self.build_starting_basis();
        if self.info_.status_ipm != IPX_STATUS_not_run || self.info_.centring_tried != 0 {
            return;
        }
        self.run_main_ipm(&mut ipm);
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::MakeIPMStartingPointValid (lines 415-468)
    #[allow(clippy::needless_range_loop)]
    fn make_ipm_starting_point_valid(&mut self) {
        let m = self.model_.rows();
        let n = self.model_.cols();
        let lb = self.model_.lb();
        let ub = self.model_.ub();
        let xl = &mut self.xl_start_;
        let xu = &mut self.xu_start_;
        let zl = &mut self.zl_start_;
        let zu = &mut self.zu_start_;

        let mut num_complementarity_products: Int = 0;
        let mut sum_complementarity_products: f64 = 0.0;
        for j in 0..(n + m) as usize {
            if xl[j] > 0.0 && zl[j] > 0.0 {
                sum_complementarity_products += xl[j] * zl[j];
                num_complementarity_products += 1;
            }
            if xu[j] > 0.0 && zu[j] > 0.0 {
                sum_complementarity_products += xu[j] * zu[j];
                num_complementarity_products += 1;
            }
        }
        let mu = if num_complementarity_products != 0 {
            sum_complementarity_products / num_complementarity_products as f64
        } else {
            1.0
        };

        for j in 0..(n + m) as usize {
            if lb[j].is_finite() {
                debug_assert!(xl[j].is_finite() && xl[j] >= 0.0);
                debug_assert!(zl[j].is_finite() && zl[j] >= 0.0);
                if xl[j] == 0.0 && zl[j] == 0.0 {
                    let s = mu.sqrt();
                    xl[j] = s;
                    zl[j] = s;
                } else if xl[j] == 0.0 {
                    xl[j] = mu / zl[j];
                } else if zl[j] == 0.0 {
                    zl[j] = mu / xl[j];
                }
            } else {
                debug_assert!(xl[j] == f64::INFINITY);
                debug_assert!(zl[j] == 0.0);
            }
            if ub[j].is_finite() {
                debug_assert!(xu[j].is_finite() && xu[j] >= 0.0);
                debug_assert!(zu[j].is_finite() && zu[j] >= 0.0);
                if xu[j] == 0.0 && zu[j] == 0.0 {
                    let s = mu.sqrt();
                    xu[j] = s;
                    zu[j] = s;
                } else if xu[j] == 0.0 {
                    xu[j] = mu / zu[j];
                } else if zu[j] == 0.0 {
                    zu[j] = mu / xu[j];
                }
            } else {
                debug_assert!(xu[j] == f64::INFINITY);
                debug_assert!(zu[j] == 0.0);
            }
        }
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::ComputeStartingPoint (lines 470-479)
    fn compute_starting_point(&mut self, ipm: &mut IPM) {
        let timer = Timer::new(0.0);
        let mut kkt = KKTSolverDiag::new(self.control_.clone(), self.model_.clone());

        // If the starting point procedure fails, then iterate_ remains as
        // initialized by the constructor, which is a valid state for
        // postprocessing/postsolving.
        ipm.starting_point(
            &mut kkt,
            self.iterate_
                .as_deref_mut()
                .expect("iterate_ is set by interior_point_solve"),
            &mut self.info_,
        );
        self.info_.time_ipm1 += timer.elapsed();
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::RunInitialIPM (lines 481-517)
    fn run_initial_ipm(&mut self, ipm: &mut IPM) {
        let timer = Timer::new(0.0);
        let mut kkt = KKTSolverDiag::new(self.control_.clone(), self.model_.clone());

        let switchiter = self.control_.switchiter();
        if switchiter < 0 {
            // Switch iteration not specified by user. Run as long as KKT solver
            // converges within min(500,10+m/20) iterations.
            let m = self.model_.rows();
            kkt.set_maxiter((10 + m / 20).min(500));
            ipm.set_maxiter(self.control_.ipm_maxiter());
        } else {
            ipm.set_maxiter(switchiter.min(self.control_.ipm_maxiter()));
        }
        ipm.driver(
            &mut kkt,
            self.iterate_
                .as_deref_mut()
                .expect("iterate_ is set by interior_point_solve"),
            &mut self.info_,
        );
        match self.info_.status_ipm {
            IPX_STATUS_optimal => {
                // If the IPM reached its termination criterion in the initial
                // iterations (happens rarely), we still call the IPM again with
                // basis preconditioning. In exact arithmetic it would terminate
                // without an additional iteration. A starting basis is then
                // available for crossover.
                self.info_.status_ipm = IPX_STATUS_not_run;
            }
            IPX_STATUS_no_progress => {
                self.info_.status_ipm = IPX_STATUS_not_run;
            }
            IPX_STATUS_failed => {
                self.info_.status_ipm = IPX_STATUS_not_run;
                self.info_.errflag = 0;
            }
            // stopped at switchiter (other iter_limit stops fall through, as in C++)
            IPX_STATUS_iter_limit if self.info_.iter < self.control_.ipm_maxiter() => {
                self.info_.status_ipm = IPX_STATUS_not_run;
            }
            _ => {}
        }
        self.info_.time_ipm1 += timer.elapsed();
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::BuildStartingBasis (lines 519-555)
    fn build_starting_basis(&mut self) {
        if self.control_.stop_at_switch() < 0 {
            self.info_.status_ipm = IPX_STATUS_debug;
            return;
        }
        self.control_.h_log_str(" Constructing starting basis...\n");
        self.basis_ = Some(Box::new(Basis::new(
            self.control_.clone(),
            self.model_.clone(),
        )));
        super::starting_basis::starting_basis(
            self.iterate_
                .as_deref_mut()
                .expect("iterate_ is set by interior_point_solve"),
            self.basis_
                .as_deref_mut()
                .expect("basis_ was just constructed"),
            &mut self.info_,
        );
        if self.info_.errflag == IPX_ERROR_user_interrupt {
            self.info_.errflag = 0;
            self.info_.status_ipm = IPX_STATUS_user_interrupt;
            return;
        } else if self.info_.errflag == IPX_ERROR_time_interrupt {
            self.info_.errflag = 0;
            self.info_.status_ipm = IPX_STATUS_time_limit;
            return;
        } else if self.info_.errflag != 0 {
            self.info_.status_ipm = IPX_STATUS_failed;
            return;
        }
        if self.model_.dualized() {
            std::mem::swap(
                &mut self.info_.dependent_rows,
                &mut self.info_.dependent_cols,
            );
            std::mem::swap(
                &mut self.info_.rows_inconsistent,
                &mut self.info_.cols_inconsistent,
            );
        }
        if self.control_.stop_at_switch() > 0 {
            self.info_.status_ipm = IPX_STATUS_debug;
            return;
        }
        if self.info_.rows_inconsistent != 0 {
            self.info_.status_ipm = IPX_STATUS_primal_infeas;
            return;
        }
        if self.info_.cols_inconsistent != 0 {
            self.info_.status_ipm = IPX_STATUS_dual_infeas;
        }
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::RunMainIPM (lines 557-564)
    fn run_main_ipm(&mut self, ipm: &mut IPM) {
        let control = self.control_.clone();
        let basis = self
            .basis_
            .as_deref_mut()
            .expect("basis_ is set by build_starting_basis");
        let mut kkt = KKTSolverBasis::new(control, basis);
        let timer = Timer::new(0.0);
        ipm.print_header();
        ipm.set_maxiter(self.control_.ipm_maxiter());
        ipm.driver(
            &mut kkt,
            self.iterate_
                .as_deref_mut()
                .expect("iterate_ is set by interior_point_solve"),
            &mut self.info_,
        );
        self.info_.time_ipm2 = timer.elapsed();
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::BuildCrossoverStartingPoint (lines 566-582)
    #[allow(clippy::needless_range_loop)]
    fn build_crossover_starting_point(&mut self) {
        let m = self.model_.rows();
        let n = self.model_.cols();

        // Construct a complementary primal-dual point from the final IPM iterate.
        // This usually increases the residuals to Ax=b and A'y+z=c.
        self.x_crossover_.resize(n + m);
        self.y_crossover_.resize(m);
        self.z_crossover_.resize(n + m);
        self.iterate_
            .as_deref()
            .expect("iterate_ is set by interior_point_solve")
            .drop_to_complementarity(
                &mut self.x_crossover_,
                &mut self.y_crossover_,
                &mut self.z_crossover_,
            );

        // Perform dual pushes in increasing order and primal pushes in decreasing
        // order of the scaling factors from the final IPM iterate.
        self.crossover_weights_.resize(n + m);
        let iterate = self
            .iterate_
            .as_deref()
            .expect("iterate_ is set by interior_point_solve");
        for j in 0..(n + m) as usize {
            self.crossover_weights_[j] = iterate.scaling_factor(j as Int);
        }
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::RunCrossover (lines 584-647)
    #[allow(clippy::needless_range_loop)]
    fn run_crossover(&mut self) {
        debug_assert!(self.basis_.is_some());
        let m = self.model_.rows();
        let n = self.model_.cols();
        self.basic_statuses_.clear();

        let weights: Option<&[f64]> = if self.crossover_weights_.size() != 0 {
            debug_assert!(self.crossover_weights_.size() == n + m);
            Some(&self.crossover_weights_[..])
        } else {
            None
        };

        let mut crossover = Crossover::new(self.control_.clone());
        crossover.push_all(
            self.basis_
                .as_deref_mut()
                .expect("basis_ is set before crossover"),
            &mut self.x_crossover_,
            &mut self.y_crossover_,
            &mut self.z_crossover_,
            weights,
            &mut self.info_,
        );
        self.info_.time_crossover = crossover.time_primal() + crossover.time_dual();
        self.info_.updates_crossover = crossover.primal_pivots() + crossover.dual_pivots();
        if self.info_.status_crossover != IPX_STATUS_optimal {
            // Crossover failed. Discard solution.
            self.x_crossover_.resize(0);
            self.y_crossover_.resize(0);
            self.z_crossover_.resize(0);
            return;
        }

        let basis = self
            .basis_
            .as_deref()
            .expect("basis_ is set before crossover");
        // Recompute vertex solution and set basic statuses.
        basis.compute_basic_solution(
            &mut self.x_crossover_,
            &mut self.y_crossover_,
            &mut self.z_crossover_,
        );
        self.basic_statuses_.resize((n + m) as usize, 0);
        let lb = self.model_.lb();
        let ub = self.model_.ub();
        let x = &self.x_crossover_;
        let z = &self.z_crossover_;
        for j in 0..self.basic_statuses_.len() {
            self.basic_statuses_[j] = if basis.is_basic(j as Int) {
                IPX_basic
            } else if lb[j] == ub[j] {
                if z[j] >= 0.0 {
                    IPX_nonbasic_lb
                } else {
                    IPX_nonbasic_ub
                }
            } else if x[j] == lb[j] {
                IPX_nonbasic_lb
            } else if x[j] == ub[j] {
                IPX_nonbasic_ub
            } else {
                IPX_superbasic
            };
        }
        let mut bound_text = String::new();
        bound_text.push_str(&textline("Bound violation of basic solution:"));
        bound_text.push_str(&sci2(primal_infeasibility(
            &self.model_,
            &self.x_crossover_,
        )));
        bound_text.push('\n');
        bound_text.push_str(&textline("Dual sign violation of basic solution:"));
        bound_text.push_str(&sci2(dual_infeasibility(
            &self.model_,
            &self.x_crossover_,
            &self.z_crossover_,
        )));
        bound_text.push('\n');
        self.control_.debug(1, &bound_text);
        let mut sv_text = String::new();
        sv_text.push_str(&textline("Minimum singular value of basis matrix:"));
        sv_text.push_str(&sci2(basis.min_singular_value()));
        sv_text.push('\n');
        self.control_.debug(1, &sv_text);

        // Declare crossover status "imprecise" if the vertex solution defined by
        // the final basis does not satisfy tolerances.
        self.model_.evaluate_basic_solution(
            &self.x_crossover_,
            &self.y_crossover_,
            &self.z_crossover_,
            &self.basic_statuses_,
            &mut self.info_,
        );
        if self.info_.primal_infeas > self.control_.pfeasibility_tol()
            || self.info_.dual_infeas > self.control_.dfeasibility_tol()
        {
            self.info_.status_crossover = IPX_STATUS_imprecise;
        }
    }
    /// HiGHS: ipx/lp_solver.cc::LpSolver::PrintSummary (lines 649-684)
    fn print_summary(&mut self) {
        let mut h = String::from("Summary\n");
        if !self.control_.timeless_log() {
            h.push_str(&textline("Runtime:"));
            h.push_str(&fix2(self.control_.elapsed()));
            h.push_str("s\n");
        }
        h.push_str(&textline("Status interior point solve:"));
        h.push_str(&status_string(self.info_.status_ipm));
        h.push('\n');
        h.push_str(&textline("Status crossover:"));
        h.push_str(&status_string(self.info_.status_crossover));
        h.push('\n');
        self.control_.h_log(&mut h);

        if self.info_.status_ipm == IPX_STATUS_optimal
            || self.info_.status_ipm == IPX_STATUS_imprecise
        {
            h.push_str(&textline("objective value:"));
            h.push_str(&sci8(self.info_.pobjval));
            h.push('\n');
            h.push_str(&textline("interior solution primal residual (abs/rel):"));
            h.push_str(&sci2(self.info_.abs_presidual));
            h.push_str(" / ");
            h.push_str(&sci2(self.info_.rel_presidual));
            h.push('\n');
            h.push_str(&textline("interior solution dual residual (abs/rel):"));
            h.push_str(&sci2(self.info_.abs_dresidual));
            h.push_str(" / ");
            h.push_str(&sci2(self.info_.rel_dresidual));
            h.push('\n');
            h.push_str(&textline("interior solution objective gap (abs/rel):"));
            h.push_str(&sci2(self.info_.pobjval - self.info_.dobjval));
            h.push_str(" / ");
            h.push_str(&sci2(self.info_.rel_objgap));
            h.push('\n');
            self.control_.h_log(&mut h);
        }
        if self.info_.status_crossover == IPX_STATUS_optimal
            || self.info_.status_crossover == IPX_STATUS_imprecise
        {
            h.push_str(&textline("basic solution primal infeasibility:"));
            h.push_str(&sci2(self.info_.primal_infeas));
            h.push('\n');
            h.push_str(&textline("basic solution dual infeasibility:"));
            h.push_str(&sci2(self.info_.dual_infeas));
            h.push('\n');
            self.control_.h_log(&mut h);
        }
    }
}

/// `std::copy(src.begin(), src.end(), out)`: writes `src.len()` elements at
/// the front of `out`.
fn copy_out_f64(src: &[f64], out: &mut [f64]) {
    out[..src.len()].copy_from_slice(src);
}

/// `std::copy_n` for `Int` arrays: the caller passes the element count via
/// the length of `src`.
fn copy_out_int(src: &[Int], out: &mut [Int]) {
    out[..src.len()].copy_from_slice(src);
}

/// Returns a vector of basic statuses that is consistent with the basis and
/// the bounds from the model.
/// HiGHS: ipx/lp_solver.cc::BuildBasicStatuses (static, lines 268-287)
fn build_basic_statuses(basis: &Basis) -> Vec<Int> {
    let model = basis.model();
    let m = model.rows();
    let n = model.cols();
    let lb = model.lb();
    let ub = model.ub();
    let mut basic_statuses = vec![0; (n + m) as usize];
    for j in 0..(n + m) {
        let ju = j as usize;
        basic_statuses[ju] = if basis.is_basic(j) {
            IPX_basic
        } else if lb[ju].is_finite() {
            IPX_nonbasic_lb
        } else if ub[ju].is_finite() {
            IPX_nonbasic_ub
        } else {
            IPX_superbasic
        };
    }
    basic_statuses
}

#[cfg(test)]
mod tests {
    use super::*;

    /// min x0 + x1 s.t. x0 + x1 >= 1 (`>`), x0 - x1 = 0, 0 <= x <= 10.
    fn load(lps: &mut LpSolver) -> Int {
        lps.load_model(
            2,
            0.0,
            Some(&[1.0, 1.0]),
            Some(&[0.0, 0.0]),
            Some(&[10.0, 10.0]),
            2,
            Some(&[0, 2, 4]),
            Some(&[0, 1, 0, 1]),
            Some(&[1.0, 1.0, 1.0, -1.0]),
            Some(&[1.0, 0.0]),
            Some(b">="),
        )
    }

    /// Review item L10: `SetParameters` / `SetCallback` / `LoadModel` after
    /// a `Solve` (while `basis_`/`iterate_` hold `Rc` clones of `control_`
    /// and `model_`) must not panic, and a reused solver gives the same
    /// result as a fresh one.
    #[test]
    fn reuse_after_solve() {
        let quiet = Parameters {
            display: 0,
            ..Parameters::default()
        };
        let mut fresh = LpSolver::new();
        fresh.set_parameters(quiet.clone());
        fresh.set_callback(Some(HighsCallback::default()));
        assert_eq!(load(&mut fresh), 0);
        let s0 = fresh.solve();
        let info0 = fresh.get_info();

        let mut lps = LpSolver::new();
        lps.set_parameters(quiet.clone());
        lps.set_callback(Some(HighsCallback::default()));
        assert_eq!(load(&mut lps), 0);
        lps.solve();
        assert!(lps.basis_.is_some() || lps.iterate_.is_some());
        lps.set_parameters(quiet.clone());
        let cb = lps.take_callback();
        assert!(cb.is_some());
        lps.set_callback(cb);
        assert_eq!(lps.get_parameters().display, 0);
        assert_eq!(load(&mut lps), 0);
        let s1 = lps.solve();
        let info1 = lps.get_info();
        assert_eq!(s0, s1);
        assert_eq!(info0.status_ipm, info1.status_ipm);
        assert_eq!(info0.status_crossover, info1.status_crossover);
        assert_eq!(info0.iter, info1.iter);
        assert_eq!(info0.objval.to_bits(), info1.objval.to_bits());
    }

    /// Option plumbing: the HiGHS 1.16 backport Parameters fields reach Control and the IPX Iterate.
    #[test]
    fn h16_lp_plumb_ipx_fields() {
        for on in [false, true] {
            let params = Parameters {
                display: 0,
                run_crossover: 0,
                h16_ipx_stall: on,
                h16_ipx_scaled: on,
                h16_ipm_gap_offset: on,
                ..Parameters::default()
            };
            let mut lps = LpSolver::new();
            // Control::interruptCheck asserts a callback is set (control.cc), as the wrapper does.
            lps.set_callback(Some(HighsCallback::default()));
            lps.set_parameters(params);
            assert_eq!(lps.control_.h16_ipx_stall(), on);
            assert_eq!(lps.control_.h16_ipx_scaled(), on);
            assert_eq!(lps.control_.h16_ipm_gap_offset(), on);
            assert_eq!(load(&mut lps), 0);
            lps.solve();
            let iterate = lps.iterate_.as_deref().expect("iterate after IPM solve");
            assert_eq!(iterate.h16_gap_offset(), on);
        }
    }
}
