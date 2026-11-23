//! HiGHS source: ipm/ipx/ipx_internal.h, ipx_config.h, ipx_info.h,
//! ipx_parameters.h, ipx_status.h (shared types and constants), plus one
//! Rust module per IPX `.cc` file
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # Ownership
//!
//! IPX objects keep `const Model&` and `const Control&` members for their
//! whole life (`Iterate`, `Basis`, `KKTSolver*`, `NormalMatrix`, `Crossover`,
//! ...). The crate rule "the owner keeps the data, borrowers get it as a
//! parameter" would thread `&Model`/`&Control` through every IPX call, so
//! IPX is the documented exception: these members are `Rc<Model>` /
//! `Rc<Control>`, cloned from `LpSolver::model_` / `LpSolver::control_`.
//! `Model` is immutable while shared: `Model::Load` / `Model::clear` run
//! from `LpSolver::LoadModel` / `ClearModel` after `ClearModel` has dropped
//! every clone (`iterate_`, `basis_`), and go through
//! `Rc::get_mut(..).expect(..)`. All mutable state of `Control`
//! (`parameters_`, `timer_`, `interval_`, `callback_`, the log streams) uses
//! `Cell`/`RefCell`, because `LpSolver::SetParameters` / `SetCallback` may
//! run while `basis_` still holds a clone (C++ mutates the shared object).
//!
//! Mutable references (`Basis&` in `KKTSolverBasis`, `KKTSolver*`/
//! `Iterate*`/`Info*` stored by `IPM::Driver`) follow the crate rule: a
//! lifetime-bound `&'a mut` field, or a parameter bundle.

use std::ops::{Deref, DerefMut};

pub mod basiclu_kernel;
pub mod basiclu_wrapper;
pub mod basis;
pub mod conjugate_residuals;
pub mod control;
pub mod crossover;
pub mod diagonal_precond;
pub mod forrest_tomlin;
pub mod guess_basis;
pub mod indexed_vector;
pub mod info;
pub mod ipm;
pub mod iterate;
pub mod kkt_solver;
pub mod kkt_solver_basis;
pub mod kkt_solver_diag;
pub mod linear_operator;
pub mod lp_solver;
pub mod lu_factorization;
pub mod lu_update;
pub mod maxvolume;
pub mod model;
pub mod normal_matrix;
pub mod oracle_trace;
pub mod power_method;
pub mod sparse_matrix;
pub mod sparse_utils;
pub mod splitted_normal_matrix;
pub mod starting_basis;
pub mod symbolic_invert;
pub mod timer;
pub mod utils;

/// HiGHS: ipx/ipx_internal.h `using Int = ipxint` (ipx_config.h: `typedef HighsInt ipxint`)
pub type Int = crate::HighsInt;

// ipx_internal.h
pub const kHypersparseThreshold: f64 = 0.1;
pub const kLuDependencyTol: f64 = 1e-3;
pub const kLuStabilityThreshold: f64 = 1e-12;
pub const kFtDiagErrorTol: f64 = 1e-8;

// ipx_status.h: for status
pub const IPX_STATUS_not_run: Int = 0;
pub const IPX_STATUS_solved: Int = 1000;
pub const IPX_STATUS_stopped: Int = 1005;
pub const IPX_STATUS_no_model: Int = 1006;
pub const IPX_STATUS_out_of_memory: Int = 1003;
pub const IPX_STATUS_internal_error: Int = 1004;
// for status_ipm and status_crossover
pub const IPX_STATUS_optimal: Int = 1;
pub const IPX_STATUS_imprecise: Int = 2;
pub const IPX_STATUS_primal_infeas: Int = 3;
pub const IPX_STATUS_dual_infeas: Int = 4;
pub const IPX_STATUS_user_interrupt: Int = 5;
pub const IPX_STATUS_time_limit: Int = 6;
pub const IPX_STATUS_iter_limit: Int = 7;
pub const IPX_STATUS_no_progress: Int = 8;
pub const IPX_STATUS_failed: Int = 9;
pub const IPX_STATUS_debug: Int = 10;
// error flags for invalid input
pub const IPX_ERROR_argument_null: Int = 102;
pub const IPX_ERROR_invalid_dimension: Int = 103;
pub const IPX_ERROR_invalid_matrix: Int = 104;
pub const IPX_ERROR_invalid_vector: Int = 105;
pub const IPX_ERROR_invalid_basis: Int = 107;
// error flags CR method
pub const IPX_ERROR_cr_iter_limit: Int = 201;
pub const IPX_ERROR_cr_matrix_not_posdef: Int = 202;
pub const IPX_ERROR_cr_precond_not_posdef: Int = 203;
pub const IPX_ERROR_cr_no_progress: Int = 204;
pub const IPX_ERROR_cr_inf_or_nan: Int = 205;
// error flags basis factorization
pub const IPX_ERROR_basis_singular: Int = 301;
pub const IPX_ERROR_basis_almost_singular: Int = 302;
pub const IPX_ERROR_basis_update_singular: Int = 303;
pub const IPX_ERROR_basis_repair_overflow: Int = 304;
pub const IPX_ERROR_basis_repair_search: Int = 305;
pub const IPX_ERROR_basis_too_ill_conditioned: Int = 306;
pub const IPX_ERROR_lapack_chol: Int = 401;
pub const IPX_ERROR_not_implemented: Int = 901;
pub const IPX_ERROR_user_interrupt: Int = 998;
pub const IPX_ERROR_time_interrupt: Int = 999;
// basic status
pub const IPX_basic: Int = 0;
pub const IPX_nonbasic: Int = -1;
pub const IPX_nonbasic_lb: Int = -1;
pub const IPX_nonbasic_ub: Int = -2;
pub const IPX_superbasic: Int = -3;

/// HiGHS: ipx/ipx_info.h::ipx_info + ipx_internal.h::Info (`memset` to 0)
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Info {
    pub status: Int,
    pub status_ipm: Int,
    pub status_crossover: Int,
    pub errflag: Int,
    pub num_var: Int,
    pub num_constr: Int,
    pub num_entries: Int,
    pub num_rows_solver: Int,
    pub num_cols_solver: Int,
    pub num_entries_solver: Int,
    pub dualized: Int,
    pub dense_cols: Int,
    pub centring_tried: Int,
    pub centring_success: Int,
    pub dependent_rows: Int,
    pub dependent_cols: Int,
    pub rows_inconsistent: Int,
    pub cols_inconsistent: Int,
    pub primal_dropped: Int,
    pub dual_dropped: Int,
    pub abs_presidual: f64,
    pub abs_dresidual: f64,
    pub rel_presidual: f64,
    pub rel_dresidual: f64,
    pub pobjval: f64,
    pub dobjval: f64,
    pub rel_objgap: f64,
    pub complementarity: f64,
    pub normx: f64,
    pub normy: f64,
    pub normz: f64,
    pub objval: f64,
    pub primal_infeas: f64,
    pub dual_infeas: f64,
    pub iter: Int,
    pub kktiter1: Int,
    pub kktiter2: Int,
    pub basis_repairs: Int,
    pub updates_start: Int,
    pub updates_ipm: Int,
    pub updates_crossover: Int,
    pub time_total: f64,
    pub time_ipm1: f64,
    pub time_ipm2: f64,
    pub time_starting_basis: f64,
    pub time_crossover: f64,
    pub time_kkt_factorize: f64,
    pub time_kkt_solve: f64,
    pub time_maxvol: f64,
    pub time_cr1: f64,
    pub time_cr1_AAt: f64,
    pub time_cr1_pre: f64,
    pub time_cr2: f64,
    pub time_cr2_NNt: f64,
    pub time_cr2_B: f64,
    pub time_cr2_Bt: f64,
    pub ftran_sparse: f64,
    pub btran_sparse: f64,
    pub time_ftran: f64,
    pub time_btran: f64,
    pub time_lu_invert: f64,
    pub time_lu_update: f64,
    pub mean_fill: f64,
    pub max_fill: f64,
    pub time_symb_invert: f64,
    pub maxvol_updates: Int,
    pub maxvol_skipped: Int,
    pub maxvol_passes: Int,
    pub tbl_nnz: Int,
    pub tbl_max: f64,
    pub frobnorm_squared: f64,
    pub lambdamax: f64,
    pub volume_increase: f64,
}

/// HiGHS: ipx/ipx_parameters.h::ipx_parameters + ipx_internal.h::Parameters
///
/// `logfile` (`const char*`) is `Option<String>`; `log_options`
/// (`const HighsLogOptions*`) is a clone of the HiGHS log options (IPX
/// never writes through the pointer).
#[derive(Debug, Clone)]
pub struct Parameters {
    // Solver control
    pub display: Int,
    pub logfile: Option<String>,
    pub print_interval: f64,
    pub time_limit: f64,
    pub analyse_basis_data: bool,
    // Preprocessing
    pub dualize: Int,
    pub scale: Int,
    // Interior point method
    pub ipm_maxiter: Int,
    pub ipm_feasibility_tol: f64,
    pub ipm_optimality_tol: f64,
    pub ipm_drop_primal: f64,
    pub ipm_drop_dual: f64,
    pub ipm_drop_dual_soft: bool,
    // Linear solver
    pub kkt_tol: f64,
    // Basis construction in IPM
    pub crash_basis: Int,
    pub dependency_tol: f64,
    pub volume_tol: f64,
    pub rows_per_slice: Int,
    pub maxskip_updates: Int,
    // LU factorization
    pub lu_kernel: Int,
    pub lu_pivottol: f64,
    // Crossover
    pub run_crossover: Int,
    pub start_crossover_tol: f64,
    pub pfeasibility_tol: f64,
    pub dfeasibility_tol: f64,
    // Debugging
    pub debug: Int,
    pub switchiter: Int,
    pub stop_at_switch: Int,
    pub update_heuristic: Int,
    pub maxpasses: Int,
    // Centring
    pub run_centring: Int,
    pub max_centring_steps: Int,
    pub centring_ratio_tolerance: f64,
    pub centring_ratio_reduction: f64,
    pub centring_alpha_scaling: f64,
    pub bad_products_tolerance: Int,
    // HiGHS logging parameters
    pub highs_logging: bool,
    pub timeless_log: bool,
    pub log_options: Option<crate::io::highs_io::HighsLogOptions>,
    // rusolve (HiGHS 1.16 backport, default off): inert flags, default false.
    // Set from HighsOptions in ipm/ipx_wrapper.rs solve_lp_ipx; read via Control accessors.
    pub h16_ipx_stall: bool,
    pub h16_ipx_scaled: bool,
    pub h16_ipm_gap_offset: bool,
    // rusolve: `lp_ipm_failure_policy`; an IPM that
    // stops with no progress near the optimum is declared imprecise, so crossover runs on it
    pub stall_to_crossover: bool,
    // rusolve `lp_crossover_ipm_crash`: crash basis ranks candidates by
    // the IPM point's distance to the bounds
    pub crossover_ipm_crash: bool,
    // rusolve `lp_crossover_batch_push`: group the unblocked pushes
    pub crossover_batch_push: bool,
}

impl Default for Parameters {
    /// HiGHS: ipx/ipx_internal.h::Parameters::Parameters (defaults)
    fn default() -> Self {
        Self {
            display: 1,
            logfile: None,
            print_interval: 60.0,
            analyse_basis_data: false,
            time_limit: -1.0,
            dualize: -1,
            scale: 1,
            ipm_maxiter: 300,
            ipm_feasibility_tol: 1e-6,
            ipm_optimality_tol: 1e-8,
            ipm_drop_primal: 1e-9,
            ipm_drop_dual: 1e-9,
            ipm_drop_dual_soft: false,
            kkt_tol: 0.3,
            crash_basis: 1,
            dependency_tol: 1e-6,
            volume_tol: 2.0,
            rows_per_slice: 10000,
            maxskip_updates: 10,
            lu_kernel: 0,
            lu_pivottol: 0.0625,
            run_crossover: 1,
            start_crossover_tol: 1e-8,
            pfeasibility_tol: 1e-7,
            dfeasibility_tol: 1e-7,
            debug: 0,
            switchiter: -1,
            stop_at_switch: 0,
            update_heuristic: 1,
            maxpasses: -1,
            run_centring: 0,
            max_centring_steps: 5,
            centring_ratio_tolerance: 100.0,
            centring_ratio_reduction: 1.5,
            centring_alpha_scaling: 0.5,
            bad_products_tolerance: 3,
            highs_logging: false,
            timeless_log: false,
            log_options: None,
            h16_ipx_stall: false,
            h16_ipx_scaled: false,
            h16_ipm_gap_offset: false,
            stall_to_crossover: false,
            crossover_ipm_crash: false,
            crossover_batch_push: false,
        }
    }
}

/// HiGHS: ipx/ipx_internal.h `using Vector = std::valarray<double>`
///
/// Derefs to `[f64]` (not `Vec<f64>`) so that `Vec::resize` cannot be
/// called by mistake: `valarray::resize(n)` sets *every* element to zero.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Vector(pub Vec<f64>);

impl Vector {
    /// `Vector v(n)`: `n` zeros.
    pub fn new(n: Int) -> Self {
        Vector(vec![0.0; n as usize])
    }
    /// `Vector v(value, n)` (note the valarray argument order).
    pub fn from_value(value: f64, n: Int) -> Self {
        Vector(vec![value; n as usize])
    }
    /// `Vector v(ptr, n)`: copy of `n` elements.
    pub fn from_slice(x: &[f64]) -> Self {
        Vector(x.to_vec())
    }
    /// `v.size()`
    pub fn size(&self) -> Int {
        self.0.len() as Int
    }
    /// `v.resize(n)` (libc++): destroys the contents, `n` zeros.
    pub fn resize(&mut self, n: Int) {
        self.0.clear();
        self.0.resize(n as usize, 0.0);
    }
    /// `v.resize(n, value)` (libc++): `n` copies of `value`.
    pub fn resize_value(&mut self, n: Int, value: f64) {
        self.0.clear();
        self.0.resize(n as usize, value);
    }
    /// `v = scalar`: fill, size unchanged.
    pub fn fill(&mut self, value: f64) {
        for e in self.0.iter_mut() {
            *e = value;
        }
    }
    /// `v = w` (libc++ `valarray::operator=`): resize to `w.size()` if the
    /// sizes differ, then copy.
    pub fn assign(&mut self, w: &[f64]) {
        if self.0.len() != w.len() {
            self.0.clear();
            self.0.resize(w.len(), 0.0);
        }
        self.0.copy_from_slice(w);
    }
    /// `v.sum()` (libc++): `0.0` if empty, else `v[0] + v[1] + ...` in order.
    pub fn sum(&self) -> f64 {
        let mut it = self.0.iter();
        match it.next() {
            None => 0.0,
            Some(&first) => {
                let mut r = first;
                for &e in it {
                    r += e;
                }
                r
            }
        }
    }
    /// `v.max()` (libc++): first maximal element by `<` (empty: UB in C++, panics here).
    pub fn max(&self) -> f64 {
        let mut r = self.0[0];
        for &e in &self.0[1..] {
            if r < e {
                r = e;
            }
        }
        r
    }
    /// `v.min()` (libc++): first minimal element by `<` (empty: UB in C++, panics here).
    pub fn min(&self) -> f64 {
        let mut r = self.0[0];
        for &e in &self.0[1..] {
            if e < r {
                r = e;
            }
        }
        r
    }
}

impl Deref for Vector {
    type Target = [f64];
    fn deref(&self) -> &[f64] {
        &self.0
    }
}

impl DerefMut for Vector {
    fn deref_mut(&mut self) -> &mut [f64] {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_valarray_semantics() {
        let mut v = Vector::from_slice(&[1.0, 2.0, 3.0]);
        v.resize(2);
        assert_eq!(v.0, vec![0.0, 0.0]);
        v.assign(&[4.0, 5.0, 6.0]);
        assert_eq!(v.size(), 3);
        assert_eq!(v.sum(), 15.0);
        v.fill(1.5);
        assert_eq!(v.0, vec![1.5; 3]);
        assert_eq!(Vector::default().sum(), 0.0);
        assert_eq!(Vector::from_value(2.0, 2).0, vec![2.0, 2.0]);
    }
}
