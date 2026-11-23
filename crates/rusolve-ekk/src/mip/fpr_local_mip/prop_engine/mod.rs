//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! prop_engine.h: `IndexedMinHeap`, `PropResult` and the `PropEngine` struct
//! with its inline accessors. The out-of-line methods are in
//! `activity.rs` and `propagate.rs`.
//!
//! The C++ engine holds raw pointers into problem-owned arrays; here they are
//! slices borrowed for `'a`. The pointer accessors return those slices, so
//! pointer identity (fpr_core.cpp acquire_engine) is `std::ptr::eq` on
//! `as_ptr()`.

pub(crate) mod activity;
pub(crate) mod propagate;

use super::{Deadline, VarState};
use crate::prelude::*;

/// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap::Entry (lines 56-59)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct HeapEntry {
    pub(crate) key: f64,
    pub(crate) var: HighsInt,
}

/// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap (lines 27-74)
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct IndexedMinHeap {
    pub(crate) heap: Vec<HeapEntry>,
    /// pos[var] = heap index, or kNotPresent.
    pub(crate) pos: Vec<HighsInt>,
}

impl IndexedMinHeap {
    /// mip-heuristics prop_engine.h:61
    pub(crate) const kNotPresent: HighsInt = -1;

    /// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap::empty (line 37)
    pub(crate) fn empty(&self) -> bool {
        self.heap.is_empty()
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap::size (line 38)
    pub(crate) fn size(&self) -> usize {
        self.heap.len()
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap::contains (lines 39-41)
    pub(crate) fn contains(&self, var: HighsInt) -> bool {
        var >= 0 && (var as usize) < self.pos.len() && self.pos[var as usize] != Self::kNotPresent
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap::top_var (line 44)
    pub(crate) fn top_var(&self) -> HighsInt {
        self.heap[0].var
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::IndexedMinHeap::entry_less (lines 65-67)
    pub(crate) fn entry_less(a: &HeapEntry, b: &HeapEntry) -> bool {
        a.key < b.key || (a.key == b.key && a.var < b.var)
    }
}

/// mip-heuristics (36f7d9d): prop_engine.h::PropResult (lines 93-102)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)] // C++ k-prefixed enumerator names
pub(crate) enum PropResult {
    kFixpoint,
    kInfeasible,
    kBudgetExhausted,
    kDeadlineExpired,
}

/// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::ActUndo (lines 358-364)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct ActUndo {
    pub(crate) row: HighsInt,
    pub(crate) old_min_ninf: HighsInt,
    pub(crate) old_max_ninf: HighsInt,
    pub(crate) old_min_finite: f64,
    pub(crate) old_max_finite: f64,
}

/// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::PqUndo (lines 373-379)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct PqUndo {
    pub(crate) var: HighsInt,
    /// domain size before change
    pub(crate) old_dom: f64,
    /// domain size after change
    pub(crate) new_dom: f64,
    /// was the var in the PQ before this change?
    pub(crate) was_present: bool,
    /// is the var in the PQ after this change?
    pub(crate) is_present: bool,
}

/// mip-heuristics (36f7d9d): prop_engine.h::PropEngine (lines 107-382)
#[derive(Clone, Debug)]
pub(crate) struct PropEngine<'a> {
    pub(crate) ncol_: HighsInt,
    pub(crate) nrow_: HighsInt,
    pub(crate) ar_start_: &'a [HighsInt],
    pub(crate) ar_index_: &'a [HighsInt],
    pub(crate) ar_value_: &'a [f64],
    pub(crate) col_start_: &'a [HighsInt],
    pub(crate) col_row_: &'a [HighsInt],
    pub(crate) col_val_: &'a [f64],
    pub(crate) col_lb_: &'a [f64],
    pub(crate) col_ub_: &'a [f64],
    pub(crate) row_lo_: &'a [f64],
    pub(crate) row_hi_: &'a [f64],
    pub(crate) integrality_: &'a [HighsVarType],
    pub(crate) feastol_: f64,
    pub(crate) nnz_: usize,
    pub(crate) vs_: Vec<VarState>,
    pub(crate) solution_: Vec<f64>,
    pub(crate) vs_undo_: Vec<(HighsInt, VarState)>,
    pub(crate) sol_undo_: Vec<(HighsInt, f64)>,
    pub(crate) prop_worklist_: Vec<HighsInt>,
    pub(crate) prop_in_wl_: Vec<u8>,
    /// C++ initialiser: 0 (prop_engine.h:334)
    pub(crate) prop_work_: usize,
    pub(crate) deadline_: Deadline,
    pub(crate) min_activity_: Vec<f64>,
    pub(crate) max_activity_: Vec<f64>,
    pub(crate) min_finite_: Vec<f64>,
    pub(crate) max_finite_: Vec<f64>,
    pub(crate) min_ninf_: Vec<HighsInt>,
    pub(crate) max_ninf_: Vec<HighsInt>,
    pub(crate) act_undo_: Vec<ActUndo>,
    /// C++ initialiser: false (prop_engine.h:371)
    pub(crate) pq_active_: bool,
    pub(crate) domain_pq_: IndexedMinHeap,
    pub(crate) pq_undo_: Vec<PqUndo>,
}

impl<'a> PropEngine<'a> {
    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::var (line 230)
    pub(crate) fn var(&self, j: HighsInt) -> &VarState {
        &self.vs_[j as usize]
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::var (line 230), non-const overload
    pub(crate) fn var_mut(&mut self, j: HighsInt) -> &mut VarState {
        &mut self.vs_[j as usize]
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::sol (line 233)
    pub(crate) fn sol(&self, j: HighsInt) -> f64 {
        self.solution_[j as usize]
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::sol (line 232), non-const overload
    pub(crate) fn sol_mut(&mut self, j: HighsInt) -> &mut f64 {
        &mut self.solution_[j as usize]
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::sol_data (line 237)
    pub(crate) fn sol_data(&self) -> &[f64] {
        &self.solution_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::sol_data (line 236), non-const overload
    pub(crate) fn sol_data_mut(&mut self) -> &mut [f64] {
        &mut self.solution_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::row_min_activity (line 246)
    pub(crate) fn row_min_activity(&self, i: HighsInt) -> f64 {
        self.min_activity_[i as usize]
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::row_max_activity (line 247)
    pub(crate) fn row_max_activity(&self, i: HighsInt) -> f64 {
        self.max_activity_[i as usize]
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::min_activity_data (line 248)
    pub(crate) fn min_activity_data(&self) -> &[f64] {
        &self.min_activity_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::max_activity_data (line 249)
    pub(crate) fn max_activity_data(&self) -> &[f64] {
        &self.max_activity_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::activities_initialized (line 250)
    pub(crate) fn activities_initialized(&self) -> bool {
        !self.min_activity_.is_empty()
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::pq_initialized (line 258)
    pub(crate) fn pq_initialized(&self) -> bool {
        self.pq_active_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::set_deadline (line 276)
    pub(crate) fn set_deadline(&mut self, d: &Deadline) {
        self.deadline_ = *d;
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::deadline (line 277)
    pub(crate) fn deadline(&self) -> &Deadline {
        &self.deadline_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::effort (line 280)
    pub(crate) fn effort(&self) -> usize {
        self.prop_work_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::add_effort (line 281)
    pub(crate) fn add_effort(&mut self, e: usize) {
        self.prop_work_ += e;
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::ncol (line 283)
    pub(crate) fn ncol(&self) -> HighsInt {
        self.ncol_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::nrow (line 284)
    pub(crate) fn nrow(&self) -> HighsInt {
        self.nrow_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::feastol (line 285)
    pub(crate) fn feastol(&self) -> f64 {
        self.feastol_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::is_int (lines 286-288)
    pub(crate) fn is_int(&self, j: HighsInt) -> bool {
        self.integrality_[j as usize] != HighsVarType::kContinuous
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::col_lb (line 291)
    pub(crate) fn col_lb(&self) -> &'a [f64] {
        self.col_lb_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::col_ub (line 292)
    pub(crate) fn col_ub(&self) -> &'a [f64] {
        self.col_ub_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::row_lo (line 293)
    pub(crate) fn row_lo(&self) -> &'a [f64] {
        self.row_lo_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::row_hi (line 294)
    pub(crate) fn row_hi(&self) -> &'a [f64] {
        self.row_hi_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::integrality (line 295)
    pub(crate) fn integrality(&self) -> &'a [HighsVarType] {
        self.integrality_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::ar_start (line 296)
    pub(crate) fn ar_start(&self) -> &'a [HighsInt] {
        self.ar_start_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::ar_index (line 297)
    pub(crate) fn ar_index(&self) -> &'a [HighsInt] {
        self.ar_index_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::ar_value (line 298)
    pub(crate) fn ar_value(&self) -> &'a [f64] {
        self.ar_value_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::csc_start (line 299)
    pub(crate) fn csc_start(&self) -> &'a [HighsInt] {
        self.col_start_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::csc_row (line 300)
    pub(crate) fn csc_row(&self) -> &'a [HighsInt] {
        self.col_row_
    }

    /// mip-heuristics (36f7d9d): prop_engine.h::PropEngine::csc_val (line 301)
    pub(crate) fn csc_val(&self) -> &'a [f64] {
        self.col_val_
    }
}
