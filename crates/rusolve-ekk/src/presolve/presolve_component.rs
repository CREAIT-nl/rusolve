//! HiGHS source: presolve/PresolveComponent.h, presolve/PresolveComponent.cpp,
//! util/HighsComponent.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The thin wrapper `Highs` keeps around one LP presolve: it owns the
//! reduced LP, the postsolve stack and the recovered solution/basis.
//!
//! # Design: owned options/timer
//!
//! C++ stores `const HighsOptions* options_` (set by `Highs::runPresolve`
//! to `&options_`) and `HighsTimer* timer` (set by `init`). Rust stores
//! clones: `options_: HighsOptions` (assigned by `run_presolve` as
//! `self.presolve_.options_ = self.options_.clone()`, which is the
//! C++ pointer assignment) and `timer: HighsTimer` (cloned in `init`; a
//! clone shares `start_time`, so `read()` is the same clock).
//!
//! The `HighsComponent` base classes only contribute `is_valid = false`
//! (data and info); they are flattened into the two structs.

use crate::prelude::*;

/// HiGHS: PresolveComponent.h::PresolveComponentData
#[derive(Clone, Debug, Default)]
pub struct PresolveComponentData {
    /// HighsComponentData::is_valid
    pub is_valid: bool,
    pub reduced_lp_: HighsLp,
    pub post_solve_stack: HighsPostsolveStack,
    pub recovered_solution_: HighsSolution,
    pub recovered_basis_: HighsBasis,
    pub presolve_log_: HighsPresolveLog,
}

impl PresolveComponentData {
    /// HiGHS: PresolveComponent.h::PresolveComponentData::clear (lines 40-48)
    ///
    /// C++ does not clear `presolve_log_`.
    pub fn clear(&mut self) {
        self.is_valid = false;

        self.post_solve_stack = HighsPostsolveStack::default();
        self.reduced_lp_.clear();
        self.recovered_solution_.clear();
        self.recovered_basis_.clear();
    }
}

/// HiGHS: PresolveComponent.h::PresolveComponentInfo
///
/// Not referenced by the LP path apart from construction.
#[derive(Clone, Debug, Default)]
pub struct PresolveComponentInfo {
    /// HighsComponentInfo::is_valid
    pub is_valid: bool,
    pub n_rows_removed: HighsInt,
    pub n_cols_removed: HighsInt,
    pub n_nnz_removed: HighsInt,

    pub presolve_time: f64,
    pub postsolve_time: f64,
}

/// HiGHS: PresolveComponent.h::PresolveComponent
#[derive(Clone, Debug)]
pub struct PresolveComponent {
    pub info_: PresolveComponentInfo,
    pub data_: PresolveComponentData,
    /// C++ `const HighsOptions* options_` (see module doc).
    pub options_: HighsOptions,
    /// C++ `HighsTimer* timer` (see module doc).
    pub timer: HighsTimer,

    pub presolve_status_: HighsPresolveStatus,
    pub postsolve_status_: HighsPostsolveStatus,
}

impl Default for PresolveComponent {
    fn default() -> Self {
        Self {
            info_: PresolveComponentInfo::default(),
            data_: PresolveComponentData::default(),
            options_: HighsOptions::default(),
            timer: HighsTimer::default(),
            presolve_status_: HighsPresolveStatus::kNotPresolved,
            postsolve_status_: HighsPostsolveStatus::kNotPresolved,
        }
    }
}

impl PresolveComponent {
    /// HiGHS: PresolveComponent.cpp::PresolveComponent::init (lines 16-22, 7 lines)
    ///
    /// C++ default arguments: `mip = false`.
    pub fn init(&mut self, lp: &HighsLp, timer: &HighsTimer, mip: bool) -> HighsStatus {
        let _ = mip;
        self.data_
            .post_solve_stack
            .initialize_index_maps(lp.num_row_, lp.num_col_);
        self.data_.reduced_lp_ = lp.clone();
        self.timer = timer.clone();
        HighsStatus::kOk
    }

    /// HiGHS: PresolveComponent.h::PresolveComponent::getReducedProblem
    pub fn get_reduced_problem(&self) -> &HighsLp {
        &self.data_.reduced_lp_
    }

    /// HiGHS: PresolveComponent.h::PresolveComponent::getPresolveLog
    pub fn get_presolve_log(&self) -> &HighsPresolveLog {
        &self.data_.presolve_log_
    }

    /// HiGHS: PresolveComponent.cpp::PresolveComponent::negateReducedLpColDuals (lines 24-29, 6 lines)
    pub fn negate_reduced_lp_col_duals(&mut self) {
        for col in 0..self.data_.reduced_lp_.num_col_ {
            self.data_.recovered_solution_.col_dual[col as usize] =
                -self.data_.recovered_solution_.col_dual[col as usize];
        }
    }

    /// HiGHS: PresolveComponent.cpp::PresolveComponent::run (lines 31-43, 13 lines)
    ///
    /// C++ passes `data_.reduced_lp_` by reference into `okSetInput`, which
    /// keeps `HighsLp* model` pointing at it. Rust moves the LP into the
    /// HPresolve (which owns `model: HighsLp`) and moves it back afterwards,
    /// on both paths.
    pub fn run(&mut self) -> HighsPresolveStatus {
        let mut presolve = HPresolve::default();
        let lp = std::mem::take(&mut self.data_.reduced_lp_);
        if !presolve.ok_set_input(
            lp,
            &self.options_,
            self.options_.presolve_reduction_limit,
            Some(self.timer.clone()),
        ) {
            self.data_.reduced_lp_ = presolve.take_model();
            self.presolve_status_ = HighsPresolveStatus::kOutOfMemory;
            return self.presolve_status_;
        }

        presolve.run(&mut self.data_.post_solve_stack);
        self.data_.reduced_lp_ = presolve.take_model();
        self.data_.presolve_log_ = presolve.get_presolve_log().clone();
        self.presolve_status_ = presolve.get_presolve_status();
        self.presolve_status_
    }

    /// HiGHS: PresolveComponent.cpp::PresolveComponent::clear (line 45, 1 line)
    pub fn clear(&mut self) {
        self.data_.clear();
    }
}
