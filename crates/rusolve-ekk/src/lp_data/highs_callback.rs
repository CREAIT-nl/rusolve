//! HiGHS source: lp_data/HighsCallback.h, lp_data/HighsCallback.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The LP path only needs the surface that `HEkk::bailout` uses
//! (`user_callback`, `active[]`, `data_out`, `clearHighsCallbackOutput`,
//! `callbackAction`); with no user callback set the C++ branch is never
//! taken. The MIP port fills in the MIP fields
//! of `HighsCallbackOutput`/`HighsCallbackInput` and the remaining methods.
//!
//! Ownership: C++ `std::function<void(int, const std::string&, const
//! HighsCallbackOutput*, HighsCallbackInput*, void*)> user_callback` is an
//! `Option<HighsCallbackFunctionType>` (an `Arc<Mutex>` closure, so `HighsCallback`
//! stays `Clone` and `Send`); `void* user_callback_data` is captured by the closure.
//! The C API (`c_callback`, `HighsCallbackDataOut/In`) is not ported.

use crate::HighsInt;
use crate::lp_data::hconst::*;
use crate::lp_data::highs_status::HighsStatus;
use std::sync::{Arc, Mutex, MutexGuard};

/// HiGHS: HighsCallback.h::ExternalMipSolutionQueryOrigin (plain enum)
pub type ExternalMipSolutionQueryOrigin = HighsInt;
pub const kExternalMipSolutionQueryOriginAfterSetup: ExternalMipSolutionQueryOrigin = 0;
pub const kExternalMipSolutionQueryOriginBeforeDive: ExternalMipSolutionQueryOrigin = 1;
pub const kExternalMipSolutionQueryOriginEvaluateRootNode0: ExternalMipSolutionQueryOrigin = 2;
pub const kExternalMipSolutionQueryOriginEvaluateRootNode1: ExternalMipSolutionQueryOrigin = 3;
pub const kExternalMipSolutionQueryOriginEvaluateRootNode2: ExternalMipSolutionQueryOrigin = 4;
pub const kExternalMipSolutionQueryOriginEvaluateRootNode3: ExternalMipSolutionQueryOrigin = 5;
pub const kExternalMipSolutionQueryOriginEvaluateRootNode4: ExternalMipSolutionQueryOrigin = 6;

/// HiGHS: HighsCallback.h::HighsCallbackOutput (lines 36-59).
/// `Highs* highs` is dropped (the callback closure captures what it needs);
/// `HighsLogType log_type` is kept as its `HighsInt` value.
#[derive(Clone, Debug, Default)]
pub struct HighsCallbackOutput {
    pub log_type: HighsInt,
    pub running_time: f64,
    pub simplex_iteration_count: HighsInt,
    pub ipm_iteration_count: HighsInt,
    pub pdlp_iteration_count: HighsInt,
    pub objective_function_value: f64,
    pub mip_node_count: i64,
    pub mip_total_lp_iterations: i64,
    pub mip_primal_bound: f64,
    pub mip_dual_bound: f64,
    pub mip_gap: f64,
    pub mip_solution: Vec<f64>,
    pub cutpool_num_col: HighsInt,
    pub cutpool_num_cut: HighsInt,
    pub cutpool_start: Vec<HighsInt>,
    pub cutpool_index: Vec<HighsInt>,
    pub cutpool_value: Vec<f64>,
    pub cutpool_lower: Vec<f64>,
    pub cutpool_upper: Vec<f64>,
    pub external_solution_query_origin: ExternalMipSolutionQueryOrigin,
    /// rusolve A4 (not in HiGHS): the `kSolutionSource*` code of the incumbent an
    /// improving-solution callback reports (`-1` = unknown). Set by
    /// `HighsMipSolverData::set_callback_data_out` from `HighsCallback::last_solution_source`;
    /// read-only information for the callback, no search logic reads it.
    pub mip_solution_source: HighsInt,
}

/// HiGHS: HighsCallback.h::HighsCallbackInput (lines 61-76). `Highs* highs`
/// is replaced by the column count where `setSolution`/`repairSolution`
/// need it.
#[derive(Clone, Debug, Default)]
pub struct HighsCallbackInput {
    pub user_interrupt: bool,
    pub user_has_solution: bool,
    pub user_solution: Vec<f64>,
}

impl HighsCallbackInput {
    /// HiGHS: HighsCallback.cpp::HighsCallbackInput::setSolution (dense overload)
    ///
    /// `num_col` is `user_solution.len()`: `clear_highs_callback_input` sizes the
    /// buffer to the model's column count. The C++ error log (through
    /// `highs->getOptions().log_options`) is not written, because this struct
    /// carries no log options; the status is the same.
    pub fn set_solution(&mut self, num_entries: HighsInt, value: &[f64]) -> HighsStatus {
        let num_col = self.user_solution.len() as HighsInt;
        if num_entries <= num_col {
            let n = num_entries.max(0) as usize;
            self.user_solution[..n].copy_from_slice(&value[..n]);
            self.user_has_solution = true;
            HighsStatus::kOk
        } else {
            HighsStatus::kError
        }
    }

    /// HiGHS: HighsCallback.cpp::HighsCallbackInput::setSolution (sparse overload)
    ///
    /// PORT-BLOCKED: the body reads `lp.num_col_`, `lp.col_lower_`,
    /// `lp.col_upper_`, `options.primal_feasibility_tolerance` and the log
    /// options of the owning `Highs`. None of these is reachable from
    /// `HighsCallbackInput`, and the signature cannot carry them.
    ///
    /// Left `todo!()` on purpose: this is user-callback API only. No caller
    /// exists under reference/HiGHS/highs/mip (the MIP solver reads the dense
    /// `data_in.user_solution` directly), so no MIP path reaches it.
    pub fn set_solution_sparse(
        &mut self,
        num_entries: HighsInt,
        index: &[HighsInt],
        value: &[f64],
    ) -> HighsStatus {
        todo!()
    }

    /// HiGHS: HighsCallback.cpp::HighsCallbackInput::repairSolution
    ///
    /// PORT-BLOCKED: the body clones the owning `Highs`, passes its model,
    /// fixes columns with `changeColBounds`, and runs it under a callback. The
    /// owning `Highs` is not reachable from `HighsCallbackInput`, and the
    /// signature cannot carry it.
    ///
    /// Left `todo!()` on purpose: this is user-callback API only, with no
    /// caller under reference/HiGHS/highs/mip, so no MIP path reaches it.
    pub fn repair_solution(&mut self) -> HighsStatus {
        todo!()
    }
}

/// C++ `HighsCallbackFunctionType`; the `void* user_callback_data` argument
/// is whatever the closure captures.
///
/// Rust-only `Send` change (CB-SEND, manifest open item 5): the holder is an
/// `Arc<Mutex<dyn FnMut + Send>>` instead of `Rc<RefCell<dyn FnMut>>`, so
/// `HighsOptions`, `HEkk` and `HighsMipSolver` can move into a task. `borrow_mut()`
/// returns the lock guard, so the call sites `cb.borrow_mut()(...)` are unchanged.
#[derive(Clone)]
pub struct HighsCallbackFunctionType(
    Arc<Mutex<dyn FnMut(HighsInt, &str, &HighsCallbackOutput, &mut HighsCallbackInput) + Send>>,
);

impl HighsCallbackFunctionType {
    /// Wraps a user closure. The closure must be `Send` (it runs under a mutex).
    pub fn new<F>(callback: F) -> Self
    where
        F: FnMut(HighsInt, &str, &HighsCallbackOutput, &mut HighsCallbackInput) + Send + 'static,
    {
        Self(Arc::new(Mutex::new(callback)))
    }

    /// Locks the callback for one call. A poisoned lock means an earlier user
    /// callback panicked, which is a bug in the user closure; it is not recovered.
    pub fn borrow_mut(
        &self,
    ) -> MutexGuard<
        '_,
        dyn FnMut(HighsInt, &str, &HighsCallbackOutput, &mut HighsCallbackInput) + Send + 'static,
    > {
        self.0.lock().expect("user callback mutex poisoned")
    }
}

/// HiGHS: HighsCallback.h::HighsCallback (lines 82-100)
#[derive(Clone)]
pub struct HighsCallback {
    /// C++: `HighsCallbackFunctionType user_callback = nullptr`
    pub user_callback: Option<HighsCallbackFunctionType>,
    pub active: Vec<bool>,
    pub data_out: HighsCallbackOutput,
    pub data_in: HighsCallbackInput,
    /// optional MIP hooks (`mip/hooks.rs`); `None` leaves the search untouched.
    pub hooks: Option<crate::mip::hooks::MipHooksHandle>,
    /// counters of hook activity and gate rejections.
    pub hook_stats: crate::mip::hooks::HookStatsHandle,
    /// source code of the incumbent being added (`kSolutionSource*`, `-1` none),
    /// recorded by `add_incumbent` only while a user callback is installed.
    pub last_solution_source: HighsInt,
}

impl std::fmt::Debug for HighsCallback {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsCallback")
            .field("user_callback", &self.user_callback.is_some())
            .field("active", &self.active)
            .field("data_out", &self.data_out)
            .field("data_in", &self.data_in)
            .finish()
    }
}

impl Default for HighsCallback {
    fn default() -> Self {
        Self {
            user_callback: None,
            active: vec![false; kNumCallbackType as usize],
            data_out: HighsCallbackOutput::default(),
            data_in: HighsCallbackInput::default(),
            hooks: None,
            hook_stats: Default::default(),
            last_solution_source: -1,
        }
    }
}

impl HighsCallback {
    /// HiGHS: HighsCallback.cpp::HighsCallback::clearHighsCallbackOutput (lines 21-40)
    pub fn clear_highs_callback_output(&mut self) {
        // Field by field, as C++: `mip_solution` and the cut-pool vectors are
        // cleared, not reallocated, and the fields C++ does not touch keep
        // their values.
        self.data_out.log_type = HighsLogType::kInfo as HighsInt;
        self.data_out.running_time = -1.0;
        self.data_out.simplex_iteration_count = -1;
        self.data_out.ipm_iteration_count = -1;
        self.data_out.pdlp_iteration_count = -1;
        self.data_out.objective_function_value = -kHighsInf;
        self.data_out.mip_node_count = -1;
        self.data_out.mip_primal_bound = kHighsInf;
        self.data_out.mip_dual_bound = -kHighsInf;
        self.data_out.mip_gap = -1.0;
        self.data_out.mip_solution.clear();
        self.data_out.cutpool_start.clear();
        self.data_out.cutpool_index.clear();
        self.data_out.cutpool_value.clear();
        self.data_out.cutpool_lower.clear();
        self.data_out.cutpool_upper.clear();
        self.data_out.external_solution_query_origin = kExternalMipSolutionQueryOriginAfterSetup;
    }

    /// HiGHS: HighsCallback.cpp::HighsCallback::clearHighsCallbackInput (lines 42-53)
    pub fn clear_highs_callback_input(&mut self, num_col: usize) {
        // make sure buffer size is correct and reset the contents if previously used
        if self.data_in.user_has_solution || num_col != self.data_in.user_solution.len() {
            self.data_in.user_solution.clear();
            self.data_in.user_solution.resize(num_col, kHighsUndefined);
        }
        self.data_in.user_interrupt = false;
        self.data_in.user_has_solution = false;
    }

    /// HiGHS: HighsCallback.cpp::HighsCallback::clear (lines 55-61)
    pub fn clear(&mut self, num_col: usize) {
        self.user_callback = None;
        self.active.clear();
        self.active.resize(kNumCallbackType as usize, false);
        self.clear_highs_callback_output();
        self.clear_highs_callback_input(num_col);
    }

    /// HiGHS: HighsCallback.cpp::HighsCallback::callbackActive (lines 63-75)
    pub fn callback_active(&self, callback_type: HighsInt) -> bool {
        // Check that callback function has been defined
        if self.user_callback.is_none() {
            return false;
        }
        // Check that callback_type is within range
        let callback_type_ok = callback_type >= kCallbackMin && callback_type <= kCallbackMax;
        debug_assert!(callback_type_ok);
        if !callback_type_ok {
            return false;
        }
        // Don't call callback if it is not active
        debug_assert!(!self.active.is_empty());
        if !self.active[callback_type as usize] {
            return false;
        }
        true
    }

    /// HiGHS: HighsCallback.cpp::HighsCallback::callbackAction (lines 77-93)
    pub fn callback_action(&mut self, callback_type: HighsInt, message: &str) -> bool {
        if !self.callback_active(callback_type) {
            return false;
        }
        if let Some(user_callback) = self.user_callback.clone() {
            let mut callback = user_callback.borrow_mut();
            (*callback)(callback_type, message, &self.data_out, &mut self.data_in);
        }
        // Assess any action
        let action = self.data_in.user_interrupt;
        // Check for no action if case not handled internally
        if callback_type == kCallbackMipImprovingSolution
            || callback_type == kCallbackMipSolution
            || callback_type == kCallbackMipLogging
            || callback_type == kCallbackMipGetCutPool
            || callback_type == kCallbackMipDefineLazyConstraints
            || callback_type == kCallbackMipUserSolution
        {
            debug_assert!(!action);
        }
        action
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;

    #[test]
    fn w5_cb_send_user_callback_is_send() {
        // Every user-callback holder (`HighsCallback`, the log / MIP-solution sinks) is `Send`, so `HEkk` and `HighsMipSolver` can move into a task; the serial callback order and payloads are unchanged.
        fn assert_send<T: Send>() {}
        assert_send::<HighsCallback>();
        assert_send::<HighsCallbackFunctionType>();
        assert_send::<crate::io::highs_io::HighsLogOptions>();
        assert_send::<crate::lp_data::highs_options::HighsOptions>();
        assert_send::<crate::simplex::ekk::HEkk>();
    }
}
