//! HiGHS source: mip/ (the MIP solver; serial path) plus presolve/HighsSymmetry
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # The HiGHS 1.15.1 MIP solver, faithful port
//!
//! Every struct carries the full C++ field list; every method has a
//! `/// HiGHS: File.cpp::Class::method (lines a-b)` doc line naming
//! the C++ it translates. The goal is node-for-node
//! equality with HiGHS run with `threads=1` (oracle:
//! `tools/highs-oracle/mip_trace.cpp`).
//!
//! # Ownership design (the one rule set for C++ back-pointers)
//!
//! The crate-wide rules (a)-(f) of `lib.rs` still hold. The MIP solver adds:
//!
//! **M1. `HighsMipSolver` is the context.** No MIP component stores a
//! pointer or reference to the solver, its `mipdata_`, a worker, an LP, a
//! domain or a pool. A method whose C++ body reads `mipsolver` /
//! `mipsolver_` / `mipdata_` takes `ms: &mut HighsMipSolver` (the
//! translator may downgrade to `&HighsMipSolver` when the body only reads).
//! `HighsMipSolverData` methods are `impl HighsMipSolver` (C++ `this->x` /
//! `x` -> `self.mipdata_.x`, `mipsolver.y` -> `self.y`); a name that
//! exists on both classes gets the suffix `_data` on the data version.
//! The model is read with [`mip_model!`]`(ms)` (C++ `mipsolver.model_`).
//!
//! **M2. Detach-on-call.** Objects that live inside `ms` and whose methods
//! need `ms` are held in a [`Slot<T>`]: the caller detaches the object for
//! the duration of the call and puts it back afterwards ([`detach!`]):
//!
//! ```ignore
//! // C++: mipdata_->domain.propagate();
//! detach!(ms.mipdata_.domains[0], |dom| dom.propagate(ms));
//! ```
//!
//! While detached, an access through the slot panics with a message naming
//! the slot: the detached object must be reached through the parameter the
//! callee received (e.g. `ctx.lp`), exactly as C++ reaches it through `this`.
//! When the callee re-enters the detached object through `ms` (a bound
//! change inside a clique-table method reaching `addImplications`), the
//! method lends `self` back for that call with [`lend!`]; the slot's stable
//! address ([`Slot::holds`]) gives C++ pointer identity
//! (`this == &mipdata_->getDomain()`).
//! The `std::deque` members of `HighsMipSolverData` (`lps`, `cutpools`,
//! `conflictpools`, `domains`, `pseudocosts`, `workers`) are
//! `Vec<Slot<T>>`; the singletons that call back into the solver
//! (`heuristics`, `cliquetable`, `implications`, `redcostfixing`,
//! `objective_function`, `nodequeue`) are `Slot<T>`. The global domain is
//! `domains[0]` (C++ `getDomain()`), the master LP `lps[0]`.
//! `HighsMipWorker`'s pointers (`lp_`, `globaldom_`, `cutpool_`,
//! `conflictpool_`, `pseudocost_`) are `usize` indices into those
//! vectors (not always 0: `createNewWorkers` /
//! `constructAdditionalWorkerData` append entries).
//! `HighsSearch` and `HighsSeparation` methods take one bundled context
//! ([`SearchCtx`] / [`SepaCtx`]) holding the solver, the worker and the
//! objects the C++ class pointed at (`lp`, `pseudocost`), because the
//! heuristics run searches on their own local LP and pseudocost.
//! Parallel workers: workers stay attached in `ms.mipdata_.workers` while a search runs
//! (`SearchCtx::worker` is an index; `addIncumbent` loops over all workers);
//! `SepaCtx::worker` is detached by its caller.
//!
//! **M3. Pool observers via an exact lazy event log.** C++ pools
//! (`HighsCutPool`, `HighsConflictPool`) keep raw pointers to every
//! domain's `CutpoolPropagation` / `ConflictPoolPropagation` and call
//! `cutAdded` / `cutDeleted` / `conflictAdded` / `conflictDeleted` on them
//! synchronously. In Rust a pool appends a [`PoolEvent`] (with a snapshot of
//! the row so later index reuse cannot change it) to its [`PoolEventLog`];
//! each propagation object keeps a `cursor` into the log and replays the
//! pending events, in log order, at the **entry of every `HighsDomain`
//! method that reads or writes domain or propagation state**, before any
//! bound change (activities are maintained incrementally in
//! `HighsCDouble`, so replay order is part of the arithmetic). Because no
//! domain state changes between the C++ notification and the next domain
//! method call, the replayed calls see the same domain state the C++
//! calls saw: the result is exact. Registration (`addPropagationDomain`)
//! records a cursor keyed by [`DomainKey`]; the domain copy constructor /
//! `operator=` inherit the source's cursor. `HighsDomain::release(ms)`
//! replaces the destructor's `removePropagationDomain` (a missing release
//! only delays log GC; replay through an unregistered key panics).
//! *Rejected alternative:* eager notification with explicit observer lists
//! threaded through every pool signature (kept as the fallback should the
//! lazy replay ever prove inexact).

use crate::prelude::*;
use std::ops::{Deref, DerefMut};
use std::ptr::NonNull;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

/// M2: detach a [`Slot`] for the duration of one call (rule f, with a
/// panic instead of a silent default while detached).
///
/// `detach!(ms.mipdata_.lps[0], |lp| lp.resolve_lp(ms, None))`
#[macro_export]
macro_rules! detach {
    ($slot:expr, |$x:ident| $body:expr) => {{
        let mut __detached = $slot.take();
        let __result = {
            let $x = &mut *__detached;
            $body
        };
        $slot.put(__detached);
        __result
    }};
}

/// M2 re-entrancy: run `$body` with the detached value `$this` (the
/// running method's `self`) lent back to `$slot`, so that callees which
/// re-enter through `ms` (C++ `mipdata_->cliquetable` / `implications` /
/// `cutpools[i]` reached from `HighsDomain::changeBound` or `propagate`)
/// see it. `$body` must not touch `$this` (it is `T::default()` while lent);
/// reach the object as `ms.mipdata_.x` inside the body instead. A no-op
/// when `$this` is not the value detached from `$slot`.
///
/// `lend!(ms.mipdata_.cliquetable, self, globaldom.fix_col(ms, col, val, reason))`
#[macro_export]
macro_rules! lend {
    ($slot:expr, $this:expr, $body:expr) => {{
        let __lent = $slot.lend($this);
        let __result = $body;
        $slot.reclaim($this, __lent);
        __result
    }};
}

/// M2 re-entrancy over a `std::deque` member (`lps`, `domains`, `cutpools`,
/// `workers`, ...): like [`lend!`], but finds the slot of `$slots` that holds
/// `$this` ([`slot_position`]). A no-op (the body still runs) when no slot
/// holds `$this`, e.g. a search's local domain or a heuristic's local LP copy.
/// Use it when the running method does not know which slot its `self` (or a
/// `&mut` parameter) came from.
///
/// `lend_in!(ms.mipdata_.domains, propdomain, detach!(ms.mipdata_.redcostfixing, |rc| rc.propagate_root_redcost(ms)))`
#[macro_export]
macro_rules! lend_in {
    ($slots:expr, $this:expr, $body:expr) => {{
        let __pos = $crate::mip::slot_position(&$slots, &*$this);
        let __lent = match __pos {
            Some(__i) => $slots[__i].lend($this),
            None => false,
        };
        let __result = $body;
        if let Some(__i) = __pos {
            $slots[__i].reclaim($this, __lent);
        }
        __result
    }};
}

/// M1: C++ `mipsolver.model_` / `mipsolver->model_` (a `const HighsLp*`
/// that points at either the original or the presolved model).
#[macro_export]
macro_rules! mip_model {
    ($ms:expr) => {
        match $ms.model_ {
            $crate::mip::HighsMipModel::Orig => &$ms.orig_model_,
            $crate::mip::HighsMipModel::Presolved => &$ms.mipdata_.presolved_model,
        }
    };
}

pub mod activation_dff;
pub mod cdcl_dive;
pub mod cdcl_uip;
pub mod clique_table;
pub mod conflict_pool;
pub mod conflict_resolution;
pub mod cut_generation;
pub mod cut_pool;
pub mod cut_selection;
pub mod debug_sol;
pub mod dff_separator;
pub mod domain;
pub mod domain_change;
pub mod dynamic_row_matrix;
pub mod feasibility_jump;
pub mod feasibilityjump;
pub mod fj_native;
pub mod fpr_local_mip;
pub mod gfk_solve;
pub mod hooks;
pub mod implications;
pub mod ipm_replay;
pub mod knapsack_cover_separator;
pub mod lazy;
pub mod lp_aggregator;
pub mod lp_relaxation;
pub mod ls_int_lp_repair;
#[cfg(test)]
mod ls_int_lp_repair_tests;
pub mod machine_sched_separator;
pub mod mip_probe;
pub mod mip_solver;
pub mod mip_solver_data;
pub mod mip_solver_report;
pub mod mip_stats;
pub mod mip_timer;
pub mod mip_worker;
pub mod modk_separator;
pub mod netdesign_separator;
pub mod node_queue;
pub mod objective_function;
pub mod path_separator;
pub mod primal_heuristics;
pub mod pseudocost;
pub mod race;
pub mod redcost_fixing;
pub mod scip_heuristics;
pub mod search;
pub mod separation;
pub mod separator;
pub mod symmetry_dynamic;
pub mod tableau_separator;
pub mod ti_precedence_separator;
pub mod transformed_lp;

pub use clique_table::*;
pub use conflict_pool::*;
pub use cut_generation::*;
pub use cut_pool::*;
pub use debug_sol::*;
pub use domain::*;
pub use domain_change::*;
pub use dynamic_row_matrix::*;
pub use feasibility_jump::*;
pub use gfk_solve::*;
pub use implications::*;
pub use ipm_replay::*;
pub use lp_aggregator::*;
pub use lp_relaxation::*;
pub use mip_solver::*;
pub use mip_solver_data::*;
pub use mip_timer::*;
pub use mip_worker::*;
pub use modk_separator::*;
pub use node_queue::*;
pub use objective_function::*;
pub use path_separator::*;
pub use primal_heuristics::*;
pub use pseudocost::*;
pub use redcost_fixing::*;
pub use search::*;
pub use separation::*;
pub use separator::*;
pub use tableau_separator::*;
pub use transformed_lp::*;
// `feasibilityjump` (external_feasibilityjump namespace: `Var`, `Problem`,
// `RowType`, ...) is reached by path: `crate::mip::feasibilityjump::*`.

// ---------------------------------------------------------------- M2

/// M2: a detachable owned value (`Option<Box<T>>` that is always `Some`
/// except while a method of the value runs).
///
/// `addr` is the address of the value this slot owns, kept while it is
/// detached (the `Box` moves out, its heap allocation does not), so that
/// C++ identity tests `&x == &mipdata_->x` stay answerable from inside the
/// detached method ([`Slot::holds`]). It changes only when the slot owns a
/// different allocation (`new`, `put` of another box, `lend`, `clone`).
///
/// `spare` (perf) caches the placeholder that
/// [`Slot::lend`] leaves behind in the lender's `*this`: the first `lend` builds one
/// `T::default()`, later ones swap it in and out, so `lend!` costs two memcpys instead
/// of a `Box` allocation plus a `T::default()` (hash tables, a whole `Highs` for an LP
/// relaxation) and its drop. The placeholder is never read (the `lend!` contract).
///
/// `pinned` (parallel workers): while `HighsMipSolver::run_task` runs worker tasks
/// concurrently, the slots of the objects the C++ workers share (`heuristics`,
/// `cliquetable`, `implications`, `redcostfixing`, `objective_function`, `nodequeue`,
/// `domains[0]`, `cutpools[0]`, `conflictpools[0]`, `pseudocosts[0]`) are pinned: the value
/// stays in place and [`Slot::take`] hands out a pointer to it ([`Detached::Pinned`])
/// instead of moving it out, so every worker reaches the one object in place, as the C++
/// workers reach `mipdata_->x` (see [`Slot::set_pinned`] for the safety argument).
#[derive(Debug)]
pub struct Slot<T> {
    val: Option<Box<T>>,
    /// Atomic (relaxed) because other workers read it ([`slot_position`] scans a whole
    /// deque) while the owning worker's `lend` / `put` writes it.
    addr: AtomicUsize,
    spare: Option<Box<T>>,
    pinned: bool,
}

fn box_addr<T>(b: &T) -> usize {
    b as *const T as usize
}

/// M2: the value [`Slot::take`] hands out: the moved-out box, or (pinned slot) a pointer
/// to the value left in place. Dereferences to the value either way; give it back with
/// [`Slot::put`].
#[derive(Debug)]
pub enum Detached<T> {
    Owned(Box<T>),
    Pinned(NonNull<T>),
}

impl<T> Deref for Detached<T> {
    type Target = T;
    fn deref(&self) -> &T {
        match self {
            Detached::Owned(b) => b,
            // SAFETY: see `Slot::set_pinned`: the pointee lives in its pinned slot, which is
            // neither moved nor dropped until it is unpinned after the worker tasks joined.
            Detached::Pinned(p) => unsafe { p.as_ref() },
        }
    }
}

impl<T> DerefMut for Detached<T> {
    fn deref_mut(&mut self) -> &mut T {
        match self {
            Detached::Owned(b) => b,
            // SAFETY: as in `deref`; concurrent workers only read the pinned object (or use
            // its atomics) under the C++ parallel-lock protocol (`Slot::set_pinned`).
            Detached::Pinned(p) => unsafe { p.as_mut() },
        }
    }
}

impl<T: Default> Default for Slot<T> {
    fn default() -> Self {
        Slot::new(T::default())
    }
}

impl<T: Clone> Clone for Slot<T> {
    fn clone(&self) -> Self {
        let val = self.val.clone();
        let addr = val.as_deref().map_or(0, box_addr);
        Slot {
            val,
            addr: AtomicUsize::new(addr),
            spare: None,
            pinned: false,
        }
    }
}

impl<T> Slot<T> {
    pub fn new(value: T) -> Self {
        let val = Box::new(value);
        let addr = box_addr(&*val);
        Slot {
            val: Some(val),
            addr: AtomicUsize::new(addr),
            spare: None,
            pinned: false,
        }
    }

    fn set_addr(&self, addr: usize) {
        // Write only on change: a pinned or re-put value keeps its address, and other
        // workers may be reading it.
        if self.addr.load(Ordering::Relaxed) != addr {
            self.addr.store(addr, Ordering::Relaxed);
        }
    }

    /// Detach the value (panics if it is already detached). A pinned slot
    /// ([`Slot::set_pinned`]) keeps its value and hands out a pointer to it.
    pub fn take(&mut self) -> Detached<T> {
        if self.pinned {
            let val = self
                .val
                .as_deref_mut()
                .expect("Slot::take: pinned slot is empty");
            return Detached::Pinned(NonNull::from(val));
        }
        Detached::Owned(
            self.val.take().expect(
                "Slot::take: value already detached (re-entrant call: pass it as a parameter)",
            ),
        )
    }

    /// Put a detached value back.
    pub fn put(&mut self, value: Detached<T>) {
        match value {
            Detached::Owned(value) => {
                debug_assert!(self.val.is_none(), "Slot::put: slot is not detached");
                self.set_addr(box_addr(&*value));
                self.val = Some(value);
            }
            Detached::Pinned(p) => {
                debug_assert!(
                    self.pinned && self.holds(p.as_ptr()),
                    "Slot::put: pointer is not this pinned slot's value"
                );
            }
        }
    }

    /// Parallel workers: pin (or unpin) the slot; see the `pinned` field. Only
    /// `HighsMipSolver::run_task` calls it, serially, before spawning the worker tasks
    /// and after joining them; the slot must be attached.
    ///
    /// SAFETY argument (the C++ protocol, HighsMipSolver.cpp `runTask` /
    /// `parallelLockActive()`): while the tasks run, the shared objects stay in
    /// `mipdata_` and every worker reaches them through `mipdata_->x` at once. The C++
    /// worker code only reads them (`getDomain()`, `getCutPool().isDuplicate`, the clique
    /// table and implication queries, the heuristics' statistics) or writes through
    /// atomics (`HighsCutPool::numLps_`, `ageResetWhileLocked_`, the pool event cursors);
    /// every other write is guarded by `!parallelLockActive()` or `thread_safe` and moved
    /// to worker-local state (worker cut / conflict pools, `worker.randgen`, local
    /// `intcols_`), and the port keeps those branches. A pinned value is never moved or
    /// dropped (`take` does not move it, `lend` is a no-op on an attached slot), so the
    /// pointers handed out stay valid until the tasks are joined.
    pub fn set_pinned(&mut self, pinned: bool) {
        debug_assert!(self.val.is_some(), "Slot::set_pinned: slot is detached");
        self.pinned = pinned;
    }

    pub fn is_pinned(&self) -> bool {
        self.pinned
    }

    pub fn is_detached(&self) -> bool {
        self.val.is_none()
    }

    /// C++ `p == &slotted_object`: true iff `p` is the value this slot owns,
    /// attached, detached or lent (identity only; the value is not read).
    pub fn holds(&self, p: *const T) -> bool {
        self.addr.load(Ordering::Relaxed) == p as usize
    }

    /// Re-entrancy: while a method
    /// of the detached value runs as `this`, hand the value back to the slot
    /// so that a callee re-entering through `ms` (C++ `mipdata_->x`) finds it.
    /// Moves `*this` into the slot (leaving a `T::default()` placeholder, cached in `spare`, behind) and
    /// returns true; returns false and does nothing when `this` is not the
    /// value detached from this slot (e.g. a local copy outside `ms`).
    /// Undo with [`Slot::reclaim`]; use through [`lend!`].
    pub fn lend(&mut self, this: &mut T) -> bool
    where
        T: Default,
    {
        if self.val.is_some() || !self.holds(this) {
            return false;
        }
        let mut lent = self.spare.take().unwrap_or_default();
        std::mem::swap(&mut *lent, this);
        self.set_addr(box_addr(&*lent));
        self.val = Some(lent);
        true
    }

    /// Undo [`Slot::lend`]: move the (possibly modified) value back into
    /// `*this` and leave the slot detached again, owning `this`.
    pub fn reclaim(&mut self, this: &mut T, lent: bool) {
        if !lent {
            return;
        }
        let mut lent = self
            .val
            .take()
            .expect("Slot::reclaim: value already detached");
        std::mem::swap(&mut *lent, this);
        self.spare = Some(lent);
        self.set_addr(box_addr(&*this));
    }
}

impl<T> Deref for Slot<T> {
    type Target = T;
    fn deref(&self) -> &T {
        self.val.as_deref().expect(
            "Slot: value is detached by a running method; reach it through that method's parameter",
        )
    }
}

impl<T> DerefMut for Slot<T> {
    fn deref_mut(&mut self) -> &mut T {
        self.val.as_deref_mut().expect(
            "Slot: value is detached by a running method; reach it through that method's parameter",
        )
    }
}

/// C++ `const T*` to an object of a *parent* solver (a sub-MIP's `clqtableinit` /
/// `implicinit`): the sub-MIP reads the parent's object in place, as the C++ does,
/// instead of moving it out of its (possibly pinned, see [`Slot::set_pinned`]) slot.
pub struct ParentPtr<T>(NonNull<T>);

impl<T> ParentPtr<T> {
    pub fn new(value: &T) -> Self {
        ParentPtr(NonNull::from(value))
    }
}

impl<T> Clone for ParentPtr<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for ParentPtr<T> {}

impl<T> std::fmt::Debug for ParentPtr<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("ParentPtr").field(&self.0).finish()
    }
}

impl<T> Deref for ParentPtr<T> {
    type Target = T;
    fn deref(&self) -> &T {
        // SAFETY: `solveSubMip` sets the pointer to the parent's object, runs the sub-MIP
        // and clears it before returning; the parent object stays attached (or pinned) in
        // its slot meanwhile, and the sub-MIP only reads it (C++ `const` pointer).
        unsafe { self.0.as_ref() }
    }
}

// SAFETY: the pointee is only read through it (see `Deref`), possibly by sub-MIPs of
// several parallel workers at once, as the C++ `const` pointer is.
unsafe impl<T> Send for ParentPtr<T> {}
unsafe impl<T> Sync for ParentPtr<T> {}

/// C++ `&x == &mipdata_->deque[i]` searched over the whole deque: the index of
/// the slot of `slots` that holds `p` ([`Slot::holds`]), or `None`.
pub fn slot_position<T>(slots: &[Slot<T>], p: &T) -> Option<usize> {
    let p = p as *const T;
    slots.iter().position(|s| s.holds(p))
}

/// M1: which model C++ `HighsMipSolver::model_` points at.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HighsMipModel {
    /// `model_ = orig_model_`
    #[default]
    Orig,
    /// `model_ = &mipdata_->presolvedModel`
    Presolved,
}

/// M2: context of every `HighsSearch` method that touches the solver.
/// C++ members `mipsolver`, `mipworker`, `lp` and `pseudocost` of
/// `HighsSearch` are these fields (the heuristics pass their local LP and
/// pseudocost, as C++ constructs `HighsSearch` with them).
///
/// Parallel workers:
/// * `worker` is an index into `ms.mipdata_.workers`; the worker stays
///   **attached** (C++ `mipworker.x` is `ctx.worker().x` /
///   `ctx.worker_mut().x`). The worker's own search runs taken out of
///   `search_ptr_` ([`HighsSearch::with_worker_search`]).
/// * `lp` is either detached from `ms.mipdata_.lps[worker.lp_]` (the worker's
///   own search) or a heuristic's local `heurlp` / `branch()`'s `lpCopy`.
/// * `pseudocost` is either detached from
///   `ms.mipdata_.pseudocosts[worker.pseudocost_]` (the worker's own search)
///   or a heuristic's local `pscost` copy. C++ `mipworker.getPseudocost()`
///   (distinct from the search's `pseudocost`) is
///   [`SearchCtx::with_worker_pseudocost`].
pub struct SearchCtx<'a> {
    pub ms: &'a mut HighsMipSolver,
    pub worker: usize,
    pub lp: &'a mut HighsLpRelaxation,
    pub pseudocost: &'a mut HighsPseudocost,
}

impl<'a> SearchCtx<'a> {
    pub fn new(
        ms: &'a mut HighsMipSolver,
        worker: usize,
        lp: &'a mut HighsLpRelaxation,
        pseudocost: &'a mut HighsPseudocost,
    ) -> Self {
        SearchCtx {
            ms,
            worker,
            lp,
            pseudocost,
        }
    }

    /// C++ `mipworker` (read).
    pub fn worker(&self) -> &HighsMipWorker {
        &self.ms.mipdata_.workers[self.worker]
    }

    /// C++ `mipworker` (write).
    pub fn worker_mut(&mut self) -> &mut HighsMipWorker {
        &mut self.ms.mipdata_.workers[self.worker]
    }

    /// C++ `mipworker.getPseudocost()` for a callee that also takes `ms` (and
    /// the search's LP): the search's own `pseudocost` when it is the worker's
    /// (`pseudocosts[k].holds`), else `pseudocosts[k]` detached for the call
    /// (a heuristic search runs on a local copy while the worker's stays
    /// attached).
    pub fn with_worker_pseudocost<R>(
        &mut self,
        f: impl FnOnce(&mut HighsMipSolver, &mut HighsLpRelaxation, &mut HighsPseudocost) -> R,
    ) -> R {
        let k = self.ms.mipdata_.workers[self.worker].pseudocost_;
        if self.ms.mipdata_.pseudocosts[k].holds(&*self.pseudocost) {
            f(self.ms, self.lp, self.pseudocost)
        } else {
            let ms = &mut *self.ms;
            let lp = &mut *self.lp;
            detach!(ms.mipdata_.pseudocosts[k], |pc| f(ms, lp, pc))
        }
    }

    /// Re-entry through `ms` (`addIncumbent` / `trySolution` read every LP;
    /// `getLp()` / pseudocost readers): lend `lp` and `pseudocost` back to the
    /// slots they were detached from for the call ([`lend_in!`]; a no-op for
    /// a heuristic's local copies). Inside `f` both hold `T::default()`:
    /// clone what the call reads from them (the LP solution) **before**.
    pub fn lend_back<R>(&mut self, f: impl FnOnce(&mut HighsMipSolver) -> R) -> R {
        let ms = &mut *self.ms;
        let lp = &mut *self.lp;
        let pc = &mut *self.pseudocost;
        lend_in!(
            ms.mipdata_.lps,
            lp,
            lend_in!(ms.mipdata_.pseudocosts, pc, f(ms))
        )
    }
}

/// M2: context of every `HighsSeparation` method that touches the solver
/// (C++ members `mipworker_` and `lp`).
///
/// `worker` is detached from `ms.mipdata_.workers[..]` and `lp` from
/// `ms.mipdata_.lps[HighsSeparation::lp]` by the caller (C++ `lp` is a pointer
/// into `mipdata_->lps`; the worker's `lp_` is the same index). Around a call
/// that re-enters through `ms` (`resolveLp` -> `addIncumbent` loops over all
/// workers and prints `getLp().numRows()`; `addRootRedcost` runs
/// `getLp().computeBasicDegenerateDuals`), the method lends them back with
/// [`lend_in!`].
pub struct SepaCtx<'a> {
    pub ms: &'a mut HighsMipSolver,
    pub worker: &'a mut HighsMipWorker,
    pub lp: &'a mut HighsLpRelaxation,
}

// ---------------------------------------------------------------- M3

/// M3: identity of one propagation observer (a `CutpoolPropagation` /
/// `ConflictPoolPropagation` inside one `HighsDomain`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DomainKey(pub u64);

/// Process-wide source of [`DomainKey`]s: keys stay unique across the worker threads of one
/// solve (a thread-local counter would hand the same key to two threads).
static NEXT_DOMAIN_KEY: AtomicU64 = AtomicU64::new(1);

impl DomainKey {
    /// A key never handed out before (identity only; values do not affect the solve).
    pub fn fresh() -> Self {
        DomainKey(NEXT_DOMAIN_KEY.fetch_add(1, Ordering::Relaxed))
    }
}

/// M3: one C++ pool -> observer notification, with the data the observer
/// reads from the pool at notification time.
#[derive(Clone, Debug)]
pub enum PoolEvent {
    /// `CutpoolPropagation::cutAdded(cut, propagate)` (HighsCutPool.cpp:627).
    CutAdded {
        cut: HighsInt,
        propagate: bool,
        /// row of the cut at notification time (`getMatrix().getRow*`)
        inds: Vec<HighsInt>,
        vals: Vec<f64>,
        rhs: f64,
    },
    /// `CutpoolPropagation::cutDeleted(cut, deletedOnlyForPropagation)`.
    CutDeleted { cut: HighsInt, only_for_prop: bool },
    /// `ConflictPoolPropagation::conflictAdded(conflict)`.
    ConflictAdded {
        conflict: HighsInt,
        entries: Vec<HighsDomainChange>,
    },
    /// `ConflictPoolPropagation::conflictDeleted(conflict)`.
    ConflictDeleted { conflict: HighsInt },
}

/// A pool's append-only notification log plus registered observer
/// cursors (C++ `propagationDomains`). Event `i` of the pool's history is
/// `events[i - base]`; events every registered cursor has passed can be
/// dropped (`gc`).
#[derive(Clone, Debug, Default)]
pub struct PoolEventLog {
    pub events: Vec<PoolEvent>,
    /// absolute index of `events[0]`
    pub base: u64,
    /// registered observers and the absolute index of their next unread event
    pub observers: Vec<(DomainKey, CursorCell)>,
    /// `push` runs `gc` once `events.len()` reaches this (0 = the initial
    /// threshold). After each gc it is set to twice the retained length, so
    /// the gc cost is amortized O(1) per event. Memory only: gc never drops an
    /// event that a registered cursor has not passed.
    pub gc_at: usize,
}

/// An observer's cursor (C++ `propagationDomains` entry). Atomic so that a `PoolEventLog`
/// can be shared across threads; `advance` takes `&self`.
#[derive(Debug)]
pub struct CursorCell(AtomicU64);

impl CursorCell {
    pub fn new(cursor: u64) -> Self {
        Self(AtomicU64::new(cursor))
    }
    pub fn get(&self) -> u64 {
        self.0.load(Ordering::Acquire)
    }
    pub fn set(&self, cursor: u64) {
        self.0.store(cursor, Ordering::Release);
    }
}

impl Clone for CursorCell {
    fn clone(&self) -> Self {
        Self::new(self.get())
    }
}

/// First `PoolEventLog::push` gc threshold (events).
const POOL_EVENT_LOG_GC_MIN: usize = 1024;

impl PoolEventLog {
    /// Absolute index one past the last event (the cursor of a fresh observer).
    pub fn end(&self) -> u64 {
        self.base + self.events.len() as u64
    }
    /// C++ `for (p : propagationDomains) p->cutAdded(..)` etc.: append one event.
    ///
    /// Also runs the amortized [`PoolEventLog::gc`]. Nothing else calls gc, so
    /// without this the log (with its row snapshots) would grow for the whole
    /// solve. Invariant this relies on: only registered observers hold cursors
    /// (a domain copy takes its cursor from a registered source).
    pub fn push(&mut self, event: PoolEvent) {
        self.events.push(event);
        if self.events.len() >= cpp_max(self.gc_at, POOL_EVENT_LOG_GC_MIN) {
            self.gc();
            self.gc_at = 2 * self.events.len();
        }
    }
    /// C++ `addPropagationDomain(domain)`: register `key` with `cursor`.
    pub fn register(&mut self, key: DomainKey, cursor: u64) {
        debug_assert!(
            self.observers.iter().all(|(k, _)| *k != key),
            "PoolEventLog::register: domain {key:?} already registered"
        );
        debug_assert!(cursor >= self.base && cursor <= self.end());
        self.observers.push((key, CursorCell::new(cursor)));
    }
    /// C++ `removePropagationDomain(domain)`.
    ///
    /// Returns whether `key` was registered. An unknown key writes nothing (C++ erases only a
    /// matching entry): a copy made under `parallelLockActive()` never registered with pool 0,
    /// and its destruction in a worker task must not touch the shared log.
    pub fn unregister(&mut self, key: DomainKey) -> bool {
        match self.observers.iter().position(|(k, _)| *k == key) {
            Some(pos) => {
                self.observers.remove(pos);
                true
            }
            None => false,
        }
    }
    /// Events from absolute index `cursor` to the end (panics if GC dropped some).
    pub fn pending(&self, cursor: u64) -> &[PoolEvent] {
        assert!(
            cursor >= self.base,
            "PoolEventLog::pending: cursor {cursor} precedes the GC base {}",
            self.base
        );
        let start = (cursor - self.base) as usize;
        assert!(
            start <= self.events.len(),
            "PoolEventLog::pending: cursor {cursor} is past the end {}",
            self.end()
        );
        &self.events[start..]
    }
    /// Record that observer `key` has consumed up to `cursor`.
    pub fn advance(&self, key: DomainKey, cursor: u64) {
        let (_, cell) = self
            .observers
            .iter()
            .find(|(k, _)| *k == key)
            .unwrap_or_else(|| panic!("PoolEventLog::advance: domain {key:?} is not registered"));
        debug_assert!(cursor >= cell.get() && cursor <= self.end());
        cell.set(cursor);
    }
    /// Drop events every registered observer has consumed.
    pub fn gc(&mut self) {
        let min = self
            .observers
            .iter()
            .map(|(_, cursor)| cursor.get())
            .min()
            .unwrap_or_else(|| self.end());
        let dropped = (min - self.base) as usize;
        self.events.drain(..dropped);
        self.base = min;
    }
}

// ---------------------------------------------------------------- shared type aliases

/// C++ `std::multimap<double, HighsInt>` (HighsCliqueTable / HighsImplications
/// use `emplace`, ordered iteration and `erase(iterator)`): per key, values
/// in insertion order (libstdc++/libc++ `multimap::emplace` appends to the
/// equal range).
pub type CppMultimapF64 = BTreeMap<OrderedF64, Vec<HighsInt>>;

/// C++ `HighsNodeQueue::NodeSet::iterator` (= `std::set<std::pair<double,
/// int64_t>>::iterator`): stored as the key itself.
pub type NodeSetKey = (OrderedF64, i64);

/// C++ `HighsNodeQueue::NodeSet` (the pool allocator is dropped).
pub type NodeSet = BTreeSet<NodeSetKey>;

/// C++ `std::unique_ptr<HighsSeparator>` in `HighsSeparation::separators`:
/// the virtual `separateLpSolution` dispatches on this enum.
#[derive(Clone, Debug)]
pub enum HighsSeparatorBox {
    Tableau(HighsTableauSeparator),
    Path(HighsPathSeparator),
    Modk(HighsModkSeparator),
    MachineSched(crate::mip::machine_sched_separator::HighsMachineSchedSeparator),
    KnapsackCover(crate::mip::knapsack_cover_separator::HighsKnapsackCoverSeparator),
    TiPrecedence(crate::mip::ti_precedence_separator::HighsTiPrecedenceSeparator),
    Netdesign(crate::mip::netdesign_separator::HighsNetdesignSeparator),
    Dff(crate::mip::dff_separator::HighsDffSeparator),
}

/// C++ `std::mt19937` (used by the feasibility-jump heuristic only). Must
/// be bit-exact with the 32-bit Mersenne Twister, including
/// `std::uniform_real_distribution<double>(0, 1)(rng)` as libc++ computes
/// it (`generate_canonical<double, 53>`: two draws). Job M03.
#[derive(Clone, Debug)]
pub struct Mt19937 {
    /// state words (624)
    pub mt: Vec<u32>,
    /// position in the state
    pub mti: usize,
}

/// C++ `std::mt19937` default constructor: seed 5489.
impl Default for Mt19937 {
    fn default() -> Self {
        Mt19937::new(5489)
    }
}

/// MT19937 parameters (C++ `mersenne_twister_engine<uint_fast32_t, 32, 624, 397,
/// 31, 0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253>`).
const MT_N: usize = 624;
const MT_M: usize = 397;
const MT_A: u32 = 0x9908_b0df;
const MT_UPPER_MASK: u32 = 0x8000_0000;
const MT_LOWER_MASK: u32 = 0x7fff_ffff;
const MT_F: u32 = 1_812_433_253;

impl Mt19937 {
    /// C++ `std::mt19937(seed)` (default-constructed = seed 5489).
    pub fn new(seed: u32) -> Self {
        // libc++ `mersenne_twister_engine::seed`: x[0] = s, x[i] = f * (x[i-1] ^ (x[i-1] >> 30)) + i.
        let mut mt = vec![0u32; MT_N];
        mt[0] = seed;
        for i in 1..MT_N {
            let prev = mt[i - 1];
            mt[i] = MT_F
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        Mt19937 { mt, mti: 0 }
    }
    /// C++ `rng()`: libc++ `mersenne_twister_engine::operator()` (lazy one-word twist).
    pub fn next_u32(&mut self) -> u32 {
        let i = self.mti;
        let j = (i + 1) % MT_N;
        let yp = (self.mt[i] & MT_UPPER_MASK) | (self.mt[j] & MT_LOWER_MASK);
        let k = (i + MT_M) % MT_N;
        self.mt[i] = self.mt[k] ^ (yp >> 1) ^ (MT_A * (yp & 1));
        // libc++: z = x[i] ^ ((x[i] >> u) & d), with u = 11 and d = 0xffffffff.
        let mut z = self.mt[i] ^ (self.mt[i] >> 11);
        self.mti = j;
        z ^= (z << 7) & 0x9d2c_5680;
        z ^= (z << 15) & 0xefc6_0000;
        z ^ (z >> 18)
    }
    /// C++ `std::uniform_real_distribution<double>(a, b)(rng)` (libc++):
    /// `(b - a) * generate_canonical<double, 53>(rng) + a`.
    pub fn uniform_real(&mut self, a: f64, b: f64) -> f64 {
        (b - a) * self.generate_canonical_f64() + a
    }

    /// libc++ `generate_canonical<double, 53>` over a 32-bit engine:
    /// k = 2 draws, `(g1 + g2 * 2^32) / 2^64`, with the same rounding steps.
    fn generate_canonical_f64(&mut self) -> f64 {
        // rp = static_cast<double>(max - min) + 1 = 2^32; log2(rp) = 32; k = ceil(53 / 32) = 2.
        let rp: f64 = u32::MAX as f64 + 1.0;
        let mut base = rp;
        let mut sp = self.next_u32() as f64;
        for _ in 1..2 {
            sp += (self.next_u32() as f64) * base;
            base *= rp;
        }
        sp / base
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mt19937_matches_standard_sequence() {
        let mut rng = Mt19937::new(5489);
        assert_eq!(rng.next_u32(), 3_499_211_612);
        let mut rng = Mt19937::new(5489);
        let mut last = 0;
        for _ in 0..10_000 {
            last = rng.next_u32();
        }
        // C++ standard [rand.predef]: the 10000th output of a default-constructed mt19937.
        assert_eq!(last, 4_123_659_995);
        assert_eq!(Mt19937::default().next_u32(), 3_499_211_612);
    }

    #[test]
    fn uniform_real_is_two_draw_generate_canonical() {
        let mut rng = Mt19937::new(12345);
        let mut reference = rng.clone();
        for _ in 0..1000 {
            let g1 = reference.next_u32() as f64;
            let g2 = reference.next_u32() as f64;
            let expected = (g1 + g2 * 4_294_967_296.0) / 18_446_744_073_709_551_616.0;
            let got = rng.uniform_real(0.0, 1.0);
            assert_eq!(got.to_bits(), expected.to_bits());
            assert!((0.0..1.0).contains(&got));
        }
    }

    /// Differential test against libc++ `std::mt19937` and
    /// `uniform_real_distribution<double>(0, 1)` (generator and output in
    /// `tests/fixtures/mip/unit/mt19937.{cpp,out}`): seeds as the `int` FJ receives.
    #[test]
    fn mt19937_matches_libcxx_dump() {
        let expected = include_str!("../../tests/fixtures/mip/unit/mt19937.out");
        let mut got = String::new();
        for seed in [0i32, 1, -1, 12345, 5489] {
            let mut rng = Mt19937::new(seed as u32);
            let mut h: u64 = 0xcbf2_9ce4_8422_2325;
            for i in 0..10_000 {
                let v: u64 = if i % 3 == 0 {
                    rng.uniform_real(0., 1.).to_bits()
                } else {
                    u64::from(rng.next_u32())
                };
                if i < 4 {
                    got += &format!("{seed} {i} {v:016x}\n");
                }
                for b in 0..8 {
                    h ^= (v >> (8 * b)) & 0xff;
                    h = h.wrapping_mul(0x0000_0100_0000_01b3);
                }
            }
            got += &format!("{seed} fnv {h}\n");
        }
        assert_eq!(got, expected);
    }

    fn conflict_deleted(conflict: HighsInt) -> PoolEvent {
        PoolEvent::ConflictDeleted { conflict }
    }

    fn conflict_of(event: &PoolEvent) -> HighsInt {
        match event {
            PoolEvent::ConflictDeleted { conflict } => *conflict,
            _ => panic!("unexpected event"),
        }
    }

    #[test]
    fn pool_event_log_replays_in_order_from_cursor() {
        let mut log = PoolEventLog::default();
        let key = DomainKey::fresh();
        log.register(key, log.end());
        assert!(log.pending(0).is_empty());
        log.push(conflict_deleted(1));
        log.push(conflict_deleted(2));
        assert_eq!(log.end(), 2);
        let pending = log.pending(0);
        assert_eq!(pending.len(), 2);
        assert_eq!(conflict_of(&pending[0]), 1);
        assert_eq!(conflict_of(&pending[1]), 2);
        log.advance(key, 1);
        assert_eq!(log.pending(1).len(), 1);
        assert_eq!(conflict_of(&log.pending(1)[0]), 2);
    }

    #[test]
    fn pool_event_log_gc_drops_only_events_all_observers_passed() {
        let mut log = PoolEventLog::default();
        let a = DomainKey::fresh();
        let b = DomainKey::fresh();
        log.register(a, 0);
        log.register(b, 0);
        for conflict in 0..4 {
            log.push(conflict_deleted(conflict));
        }
        log.advance(a, 3);
        log.advance(b, 2);
        log.gc();
        assert_eq!(log.base, 2);
        assert_eq!(log.events.len(), 2);
        assert_eq!(log.end(), 4);
        assert_eq!(conflict_of(&log.pending(2)[0]), 2);
        log.advance(b, 4);
        log.gc();
        assert_eq!(log.base, 3);
        assert_eq!(conflict_of(&log.pending(3)[0]), 3);
    }

    #[test]
    fn pool_event_log_unregister_releases_gc_constraint() {
        let mut log = PoolEventLog::default();
        let a = DomainKey::fresh();
        log.register(a, 0);
        log.push(conflict_deleted(7));
        log.unregister(a);
        log.gc();
        assert_eq!(log.base, 1);
        assert!(log.events.is_empty());
        assert!(log.pending(1).is_empty());
    }

    #[test]
    fn pool_event_log_push_gc_is_bounded_and_keeps_lagging_cursor() {
        let mut log = PoolEventLog::default();
        let fast = DomainKey::fresh();
        let slow = DomainKey::fresh();
        log.register(fast, 0);
        log.register(slow, 0);
        let slow_cursor: u64 = 1500;
        let mut max_len = 0usize;
        for i in 0..20000 {
            log.push(conflict_deleted(i));
            log.advance(fast, log.end());
            if log.end() == slow_cursor {
                log.advance(slow, slow_cursor);
            }
            max_len = max_len.max(log.events.len());
        }
        // The slow observer still sees every event from its cursor, in order.
        let pending = log.pending(slow_cursor);
        assert_eq!(pending.len() as u64, 20000 - slow_cursor);
        for (j, e) in pending.iter().enumerate() {
            assert_eq!(conflict_of(e), slow_cursor as HighsInt + j as HighsInt);
        }
        // Once the slow observer catches up, push-driven gc bounds the log.
        log.advance(slow, log.end());
        for i in 20000..60000 {
            log.push(conflict_deleted(i));
            log.advance(fast, log.end());
            log.advance(slow, log.end());
        }
        assert!(log.events.len() <= 2 * 20000, "len {}", log.events.len());
        assert!(log.base > slow_cursor);
        assert!(max_len >= 1024);
    }

    #[test]
    #[should_panic(expected = "precedes the GC base")]
    fn pool_event_log_pending_before_gc_base_panics() {
        let mut log = PoolEventLog::default();
        let a = DomainKey::fresh();
        log.register(a, 0);
        log.push(conflict_deleted(1));
        log.advance(a, 1);
        log.gc();
        let _ = log.pending(0);
    }

    #[test]
    #[should_panic(expected = "is not registered")]
    fn pool_event_log_advance_unregistered_panics() {
        let log = PoolEventLog::default();
        log.advance(DomainKey::fresh(), 0);
    }

    #[test]
    fn slot_detach_and_put_back() {
        let mut slot = Slot::new(41u32);
        assert!(!slot.is_detached());
        let result = detach!(slot, |x| {
            *x += 1;
            *x
        });
        assert_eq!(result, 42);
        assert_eq!(*slot, 42);
        let taken = slot.take();
        assert!(slot.is_detached());
        slot.put(taken);
        assert!(!slot.is_detached());
        assert_eq!(*slot, 42);
    }

    #[test]
    #[should_panic(expected = "detached")]
    fn slot_access_while_detached_panics() {
        let mut slot = Slot::new(1u32);
        let _taken = slot.take();
        let _ = *slot;
    }

    #[test]
    fn slot_lend_hands_detached_value_back() {
        let mut slot = Slot::new(vec![1u32]);
        let mut detached = slot.take();
        let this: &mut Vec<u32> = &mut detached;
        assert!(slot.holds(this));
        let len = lend!(slot, this, {
            // a re-entrant callee reaches the value through the slot
            slot.push(2);
            slot.len()
        });
        assert_eq!(len, 2);
        assert_eq!(*this, vec![1, 2]);
        assert!(slot.is_detached());
        assert!(slot.holds(this));
        // not the detached value: lend is a no-op
        let mut other = vec![7u32];
        assert!(lend!(slot, &mut other, slot.is_detached()));
        slot.put(detached);
        assert_eq!(*slot, vec![1, 2]);
    }

    #[test]
    fn slot_default_holds_default_value() {
        let slot: Slot<Vec<u8>> = Slot::default();
        assert!(slot.is_empty());
    }

    // Unit tests: `Send` / `Sync` of the MIP types

    fn assert_send<T: Send>() {}

    #[test]
    fn w5_mip_send_mip_state_is_send() {
        // The MIP state moves into worker tasks: `Rc` -> `Arc`, the process-wide
        // `NEXT_DOMAIN_KEY`, and atomic observer cursors (`CursorCell`) make it `Send`.
        assert_send::<HighsMipSolver>();
        assert_send::<HighsMipWorker>();
        assert_send::<HighsLpRelaxation>();
        assert_send::<HighsSearch>();
        assert_send::<HighsDomain>();
    }
}
