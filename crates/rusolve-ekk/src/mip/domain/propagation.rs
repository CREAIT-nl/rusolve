//! HiGHS source: mip/HighsDomain.cpp (lines 1-1258)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! Borrow story (re-entrant objects and HighsDomain propagators):
//! * The C++ back-pointer `domain` of the three propagator structs is gone. A propagator
//!   method that only *reads and writes* domain fields takes `domain: &mut HighsDomain` (the
//!   owning domain, with its `cutpoolpropagation` / `conflict_pool_propagation` vector
//!   `mem::take`n out by the caller, so the method must not reach those vectors through
//!   `domain`). A propagator method whose C++ body calls `domain->changeBound` /
//!   `checkChangeBound` (which re-enters every propagator of the domain) is a `HighsDomain`
//!   method instead, reading `self.conflict_pool_propagation[k]` / `self.obj_prop_` afresh at
//!   each access: `HighsDomain::propagate_conflict` and `HighsDomain::obj_prop_propagate`.
//! * The C++ `cutpool` / `conflictpool_` pointers are indices into `ms.mipdata_.cutpools` /
//!   `conflictpools` (`CutpoolPropagation::cutpool`, `ConflictPoolPropagation::conflictpool`).
//!   The pool is read through `ms: &HighsMipSolver`; it must be attached (a pool method that
//!   reaches a domain lends itself back first, see `lend!`).
//! * Pool notifications (`cutAdded`, `cutDeleted`, `conflictAdded`, `conflictDeleted`) are
//!   replayed from the pool's M3 event log by `HighsDomain::replay_pool_events` (mod.rs). The
//!   row / conflict data the C++ reads from the pool at notification time comes from the
//!   event snapshot (`inds`, `vals`, `rhs`, `entries` parameters).

use crate::prelude::*;

/// HiGHS: HighsDomain.cpp::activityContributionMin (lines 20-31)
/// C++: `static double activityContributionMin(double coef, const double & lb, const double & ub)`
pub fn activity_contribution_min(coef: f64, lb: f64, ub: f64) -> f64 {
    if coef < 0.0 {
        if ub == kHighsInf {
            -kHighsInf
        } else {
            coef * ub
        }
    } else if lb == -kHighsInf {
        -kHighsInf
    } else {
        coef * lb
    }
}

/// HiGHS: HighsDomain.cpp::activityContributionMax (lines 33-44)
/// C++: `static double activityContributionMax(double coef, const double & lb, const double & ub)`
pub fn activity_contribution_max(coef: f64, lb: f64, ub: f64) -> f64 {
    if coef < 0.0 {
        if lb == -kHighsInf {
            kHighsInf
        } else {
            coef * lb
        }
    } else if ub == kHighsInf {
        kHighsInf
    } else {
        coef * ub
    }
}

/// HiGHS latest (#3184, c5ebac196): HighsDomain.cpp::activityContributionMin (lines 20-31)
/// C++: `static HighsCDouble activityContributionMin(double coef, const double & lb, const double & ub)`
pub fn activity_contribution_min_cd(coef: f64, lb: f64, ub: f64) -> HighsCDouble {
    if coef < 0.0 {
        if ub == kHighsInf {
            HighsCDouble::from(-kHighsInf)
        } else {
            coef * HighsCDouble::from(ub)
        }
    } else if lb == -kHighsInf {
        HighsCDouble::from(-kHighsInf)
    } else {
        coef * HighsCDouble::from(lb)
    }
}

/// HiGHS latest (#3184, c5ebac196): HighsDomain.cpp::activityContributionMax (lines 33-44)
/// C++: `static HighsCDouble activityContributionMax(double coef, const double & lb, const double & ub)`
pub fn activity_contribution_max_cd(coef: f64, lb: f64, ub: f64) -> HighsCDouble {
    if coef < 0.0 {
        if lb == -kHighsInf {
            HighsCDouble::from(kHighsInf)
        } else {
            coef * HighsCDouble::from(lb)
        }
    } else if ub == kHighsInf {
        HighsCDouble::from(kHighsInf)
    } else {
        coef * HighsCDouble::from(ub)
    }
}

/// HiGHS: HighsDomain.cpp::computeDelta (lines 46-59)
/// C++: `static HighsCDouble computeDelta(double val, double oldbound, double newbound, double inf, HighsInt & numinfs)`
#[inline]
pub fn compute_delta(
    val: f64,
    oldbound: f64,
    newbound: f64,
    inf: f64,
    numinfs: &mut HighsInt,
) -> HighsCDouble {
    if oldbound == inf {
        *numinfs -= 1;
        HighsCDouble::from(newbound) * val
    } else if newbound == inf {
        *numinfs += 1;
        HighsCDouble::from(-oldbound) * val
    } else {
        (HighsCDouble::from(newbound) - HighsCDouble::from(oldbound)) * val
    }
}

/// HiGHS: HighsDomain.cpp::boundRange (lines 61-67)
/// HiGHS 1.16 (#3525): HighsDomain.cpp::boundRange (lines 61-69)
/// C++: `static double boundRange(double upper_bound, double lower_bound, double tolerance, HighsVarType var_type)`
// The 6th argument is `mip_cont_prop_improve` (HiGHS: 0.3); the 5th is `h16_bound_range_inf` (`ms.options_mip_.h16_bound_range_inf` at every caller).
pub fn bound_range(
    upper_bound: f64,
    lower_bound: f64,
    tolerance: f64,
    var_type: HighsVarType,
    h16_bound_range_inf: bool,
    cont_frac: f64,
) -> f64 {
    if h16_bound_range_inf {
        // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp lines 63-64 (1.16.0)
        if lower_bound <= -kHighsInf || upper_bound >= kHighsInf {
            return kHighsInf;
        }
    }
    let range = upper_bound - lower_bound;
    range
        - (if var_type == HighsVarType::kContinuous {
            cpp_max(cont_frac * range, 1000.0 * tolerance)
        } else {
            tolerance
        })
}

impl HighsDomain {
    /// HiGHS: HighsDomain.cpp::HighsDomain::HighsDomain (lines 69-78)
    /// C++: `HighsDomain(HighsMipSolver & mipsolver)`
    // PORT NOTE: signature changed: the solver is only read (`model_`, `numCol()`).
    pub fn new(ms: &HighsMipSolver) -> Self {
        let num_col = ms.num_col() as usize;
        let model = mip_model!(ms);
        HighsDomain {
            col_lower_: model.col_lower_.clone(),
            col_upper_: model.col_upper_.clone(),
            col_lower_pos_: vec![-1; num_col],
            col_upper_pos_: vec![-1; num_col],
            changedcolsflags_: vec![0; num_col],
            changedcols_: Vec::with_capacity(num_col),
            infeasible_reason: HighsDomain_Reason::unspecified(),
            infeasible_: false,
            x1_cdouble_act_contrib: ms.options_mip_.mip_cdouble_act_contrib,
            x1_reduced_stack_fast: ms.options_mip_.mip_reduced_stack_fast,
            ..Default::default()
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::addCutpool (lines 80-83)
    /// C++: `void addCutpool(HighsCutPool & cutpool)`
    // PORT NOTE: signature changed: the pool is `ms.mipdata_.cutpools[cutpool]` (all C++ pools
    // live in `mipdata_->cutpools`); registration needs it mutably. Caller idiom for the
    // global domain: `detach!(ms.mipdata_.domains[0], |d| d.add_cutpool(ms, 0))`.
    pub fn add_cutpool(&mut self, ms: &mut HighsMipSolver, cutpool: usize) {
        let cutpoolindex = self.cutpoolpropagation.len() as HighsInt;
        let prop = HighsDomain_CutpoolPropagation::new(ms, cutpoolindex, cutpool);
        self.cutpoolpropagation.push(prop);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::addConflictPool (lines 85-88)
    /// C++: `void addConflictPool(HighsConflictPool & conflictPool)`
    // PORT NOTE: signature changed: the pool is `ms.mipdata_.conflictpools[conflict_pool]`.
    pub fn add_conflict_pool(&mut self, ms: &mut HighsMipSolver, conflict_pool: usize) {
        let conflict_pool_index = self.conflict_pool_propagation.len() as HighsInt;
        let prop = HighsDomain_ConflictPoolPropagation::new(ms, conflict_pool_index, conflict_pool);
        self.conflict_pool_propagation.push(prop);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::propagateConflict (lines 305-376)
    /// C++: `void ConflictPoolPropagation::propagateConflict(HighsInt conflict)`
    ///
    /// PORT NOTE: a `HighsDomain` method on `self.conflict_pool_propagation[k]` (C++ `this` is
    /// that element, `domain` is `self`) because the body calls `domain->changeBound`, which
    /// updates every conflict propagator of the domain, this one included. Re-index
    /// `self.conflict_pool_propagation[k]` after each `change_bound`; the C++ `watched`
    /// pointer is only used in `case 2`, where no bound changes. The pool is
    /// `ms.mipdata_.conflictpools[self.conflict_pool_propagation[k].conflictpool]`
    /// (`resetAge` needs it mutably: read entries into locals first).
    pub fn propagate_conflict(&mut self, ms: &mut HighsMipSolver, k: usize, conflict: HighsInt) {
        // HiGHS 1.16 #3571 (h16_watch_link): read once, passed to link_watched_literal (HighsDomain.cpp lines 95-107, 1.16.0)
        let h16 = ms.options_mip_.h16_watch_link;
        let ci = conflict as usize;
        // remove propagate flag, but keep watched and deleted information
        self.conflict_pool_propagation[k].conflict_flag_[ci] &= 3 | 8;
        // if two inactive literals are watched or conflict has been deleted skip
        if self.conflict_pool_propagation[k].conflict_flag_[ci] >= 2 {
            return;
        }

        if self.infeasible_ {
            return;
        }

        let pool = self.conflict_pool_propagation[k].conflictpool;
        let start = ms.mipdata_.conflictpools[pool].get_conflict_ranges()[ci].0;
        if start == -1 {
            self.conflict_pool_propagation[k].unlink_watched_literal(2 * conflict);
            self.conflict_pool_propagation[k].unlink_watched_literal(2 * conflict + 1);
            return;
        }
        let end = ms.mipdata_.conflictpools[pool].get_conflict_ranges()[ci].1;

        let mut inactive: [HighsInt; 2] = [0; 2];
        let mut num_inactive: u8 = 0;
        // PORT NOTE (perf): the C++ `entries` reference, sliced to the conflict once, so the scan
        // does not re-reach the pool slot and bounds-check every entry.
        let entries = &ms.mipdata_.conflictpools[pool].get_conflict_entry_vector()
            [start as usize..end as usize];
        for (i, entry) in (start..end).zip(entries) {
            if self.is_active(entry) {
                continue;
            }

            inactive[num_inactive as usize] = i;
            num_inactive += 1;
            if num_inactive == 2 {
                break;
            }
        }

        self.conflict_pool_propagation[k].conflict_flag_[ci] = num_inactive;

        match num_inactive {
            0 => {
                debug_assert!(!self.infeasible_);
                ms.mipdata_.debug_solution.node_pruned(self);
                self.infeasible_ = true;
                self.infeasible_reason = HighsDomain_Reason::cut(
                    self.cutpoolpropagation.len() as HighsInt
                        + self.conflict_pool_propagation[k].conflictpoolindex,
                    conflict,
                );
                self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                crate::work::add_event(crate::work::Event::ConfPropCutoffs, 1);
                // printf("conflict propagation found infeasibility\n");
                ms.mipdata_.conflictpools[pool].reset_age(conflict);
            }
            1 => {
                let entry = ms.mipdata_.conflictpools[pool].get_conflict_entry_vector()
                    [inactive[0] as usize]
                    .clone();
                let domchg = self.flip(ms, &entry);
                if !self.is_active(&domchg) {
                    let flipped = self.flip(ms, &entry);
                    let reason = HighsDomain_Reason::cut(
                        self.cutpoolpropagation.len() as HighsInt
                            + self.conflict_pool_propagation[k].conflictpoolindex,
                        conflict,
                    );
                    self.change_bound(ms, flipped, reason);
                    crate::work::add_event(crate::work::Event::ConfPropBoundChanges, 1);
                    ms.mipdata_.conflictpools[pool].reset_age(conflict);
                }
                // printf("conflict propagation found bound change\n");
            }
            _ => {
                let e0 = ms.mipdata_.conflictpools[pool].get_conflict_entry_vector()
                    [inactive[0] as usize]
                    .clone();
                if self.conflict_pool_propagation[k].watched_literals_[2 * ci].domchg != e0 {
                    self.conflict_pool_propagation[k].unlink_watched_literal(2 * conflict);
                    self.conflict_pool_propagation[k].watched_literals_[2 * ci].domchg = e0;
                    self.conflict_pool_propagation[k].link_watched_literal(2 * conflict, h16);
                }

                let e1 = ms.mipdata_.conflictpools[pool].get_conflict_entry_vector()
                    [inactive[1] as usize]
                    .clone();
                if self.conflict_pool_propagation[k].watched_literals_[2 * ci + 1].domchg != e1 {
                    self.conflict_pool_propagation[k].unlink_watched_literal(2 * conflict + 1);
                    self.conflict_pool_propagation[k].watched_literals_[2 * ci + 1].domchg = e1;
                    self.conflict_pool_propagation[k].link_watched_literal(2 * conflict + 1, h16);
                }
            }
        }
    }
}

impl HighsDomain_ConflictPoolPropagation {
    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::linkWatchedLiteral (lines 90-104)
    /// HiGHS 1.16 (#3571): HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::linkWatchedLiteral (lines 95-107)
    /// C++: `void linkWatchedLiteral(HighsInt linkPos)`
    /// `h16` (`h16_watch_link`) selects the 1.16 placement of `head = linkPos` (outside the `if`);
    /// off keeps the 1.15.1 placement (inside the `if`, so an empty list stays unreachable).
    pub fn link_watched_literal(&mut self, link_pos: HighsInt, h16: bool) {
        let lp = link_pos as usize;
        debug_assert!(self.watched_literals_[lp].domchg.column != -1);
        let col = self.watched_literals_[lp].domchg.column as usize;
        let head: &mut HighsInt =
            if self.watched_literals_[lp].domchg.boundtype == HighsBoundType::kLower {
                &mut self.col_lower_watched_[col]
            } else {
                &mut self.col_upper_watched_[col]
            };

        self.watched_literals_[lp].prev = -1;
        self.watched_literals_[lp].next = *head;
        if h16 {
            // HiGHS 1.16 #3571 (h16_watch_link): HighsDomain.cpp lines 105-106 (1.16.0): head = linkPos always
            if *head != -1 {
                self.watched_literals_[*head as usize].prev = link_pos;
            }
            *head = link_pos;
        } else if *head != -1 {
            self.watched_literals_[*head as usize].prev = link_pos;
            // HiGHS quirk: `head = linkPos` sits inside the `if`, so linking into an empty list leaves
            // the column head at -1 and the literal is unreachable from the head.
            *head = link_pos;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::unlinkWatchedLiteral (lines 106-123)
    /// C++: `void unlinkWatchedLiteral(HighsInt linkPos)`
    pub fn unlink_watched_literal(&mut self, link_pos: HighsInt) {
        let lp = link_pos as usize;
        if self.watched_literals_[lp].domchg.column == -1 {
            return;
        }

        let col = self.watched_literals_[lp].domchg.column as usize;
        let head: &mut HighsInt =
            if self.watched_literals_[lp].domchg.boundtype == HighsBoundType::kLower {
                &mut self.col_lower_watched_[col]
            } else {
                &mut self.col_upper_watched_[col]
            };
        self.watched_literals_[lp].domchg.column = -1;
        let prev = self.watched_literals_[lp].prev;
        let next = self.watched_literals_[lp].next;
        if prev != -1 {
            self.watched_literals_[prev as usize].next = next;
        } else {
            *head = next;
        }

        if next != -1 {
            self.watched_literals_[next as usize].prev = prev;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::ConflictPoolPropagation (lines 125-134)
    /// C++: `ConflictPoolPropagation(HighsInt conflictpoolindex, HighsDomain * domain, HighsConflictPool & conflictpool_)`
    // PORT NOTE: signature changed: no back-pointer; the pool is `ms.mipdata_.conflictpools[conflictpool]`.
    // `addPropagationDomain(this)` registers a fresh `DomainKey::fresh()` at the log end and
    // sets `cursor` to it.
    pub fn new(ms: &mut HighsMipSolver, conflictpoolindex: HighsInt, conflictpool: usize) -> Self {
        let num_col = ms.num_col() as usize;
        let key = DomainKey::fresh();
        ms.mipdata_.conflictpools[conflictpool].add_propagation_domain(key);
        let cursor = ms.mipdata_.conflictpools[conflictpool].events.end();
        Self {
            conflictpoolindex,
            conflictpool,
            cursor,
            key,
            col_lower_watched_: vec![-1; num_col],
            col_upper_watched_: vec![-1; num_col],
            ..Default::default()
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::ConflictPoolPropagation (lines 136-149)
    /// C++: `ConflictPoolPropagation(const ConflictPoolPropagation & other)`
    // PORT NOTE: signature changed: the copy gets a fresh key and inherits `other.cursor`
    // (registered at that cursor, so it replays the same pending events `other` would).
    // C++ skips the registration when `parallelLockActive()` and the pool is
    // `getConflictPool()` (index 0); the copy then must not replay: set `cursor` to `u64::MAX`
    // (see `HighsDomain::replay_pool_events`).
    pub fn copy_from(ms: &mut HighsMipSolver, other: &HighsDomain_ConflictPoolPropagation) -> Self {
        let key = DomainKey::fresh();
        let mut cursor = u64::MAX;
        if !ms.parallel_lock_active() || other.conflictpool != 0 {
            let pool_end = ms.mipdata_.conflictpools[other.conflictpool].events.end();
            cursor = if other.cursor == u64::MAX {
                pool_end
            } else {
                other.cursor
            };
            ms.mipdata_.conflictpools[other.conflictpool]
                .events
                .register(key, cursor);
        }
        Self {
            conflictpoolindex: other.conflictpoolindex,
            conflictpool: other.conflictpool,
            cursor,
            key,
            col_lower_watched_: other.col_lower_watched_.clone(),
            col_upper_watched_: other.col_upper_watched_.clone(),
            conflict_flag_: other.conflict_flag_.clone(),
            propagate_conflict_inds_: other.propagate_conflict_inds_.clone(),
            watched_literals_: other.watched_literals_.clone(),
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::operator= (lines 151-167)
    /// C++: `ConflictPoolPropagation & operator=(const ConflictPoolPropagation & other)`
    // PORT NOTE: signature changed: unregister `self.key` from the old pool, copy, keep
    // `self.key`, register it with the new pool at `other.cursor`.
    pub fn op_assign(
        &mut self,
        ms: &mut HighsMipSolver,
        other: &HighsDomain_ConflictPoolPropagation,
    ) {
        if std::ptr::eq(&*self, other) {
            return;
        }
        let key = self.key;
        ms.mipdata_.conflictpools[self.conflictpool].remove_propagation_domain(key);
        self.conflictpoolindex = other.conflictpoolindex;
        self.conflictpool = other.conflictpool;
        self.col_lower_watched_
            .clone_from(&other.col_lower_watched_);
        self.col_upper_watched_
            .clone_from(&other.col_upper_watched_);
        self.conflict_flag_.clone_from(&other.conflict_flag_);
        self.propagate_conflict_inds_
            .clone_from(&other.propagate_conflict_inds_);
        self.watched_literals_.clone_from(&other.watched_literals_);
        debug_assert!(!ms.parallel_lock_active());
        let pool_end = ms.mipdata_.conflictpools[self.conflictpool].events.end();
        self.cursor = if other.cursor == u64::MAX {
            pool_end
        } else {
            other.cursor
        };
        ms.mipdata_.conflictpools[self.conflictpool]
            .events
            .register(key, self.cursor);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::~ConflictPoolPropagation (lines 169-171)
    /// C++: `~ConflictPoolPropagation()`
    // PORT NOTE: explicit destructor (M3): `remove_propagation_domain(self.key)` on the pool.
    pub fn release(&mut self, ms: &mut HighsMipSolver) {
        ms.mipdata_.conflictpools[self.conflictpool].remove_propagation_domain(self.key);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::conflictDeleted (lines 173-177)
    /// HiGHS 1.16 (#3571): HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::conflictDeleted (lines 176-181)
    /// C++: `void conflictDeleted(HighsInt conflict)`
    pub fn conflict_deleted(&mut self, conflict: HighsInt, h16: bool) {
        let ci = conflict as usize;
        // HiGHS 1.16 #3571 (h16_watch_link): HighsDomain.cpp line 177 (1.16.0), debug-only assert
        if h16 {
            debug_assert!(ci < self.conflict_flag_.len());
        }
        self.conflict_flag_[ci] |= 8;
        self.unlink_watched_literal(2 * conflict);
        self.unlink_watched_literal(2 * conflict + 1);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::conflictAdded (lines 179-259)
    /// HiGHS 1.16 (#3571): HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::conflictAdded (lines 183-266)
    /// C++: `void conflictAdded(HighsInt conflict)`
    // PORT NOTE: signature changed: `entries` is the M3 snapshot of
    // `conflictEntries[start..end]`; C++ index `i` in `start..end` is `entries[i - start]`.
    // `domain` is the owning domain (read only).
    #[allow(clippy::collapsible_match)] // keeps the C++ `switch (numActive)` shape
    pub fn conflict_added(
        &mut self,
        domain: &HighsDomain,
        conflict: HighsInt,
        entries: &[HighsDomainChange],
        h16: bool,
    ) {
        let ci = conflict as usize;
        // `entries` is conflictEntries[start..end]: C++ index i is entries[i - start].
        if self.conflict_flag_.len() as HighsInt <= conflict {
            self.watched_literals_
                .resize_with(2 * ci + 2, Default::default);
            self.conflict_flag_.resize(ci + 1, 0);
        }

        let mut num_watched: HighsInt = 0;
        for i in 0..entries.len() {
            if domain.is_active(&entries[i]) {
                continue;
            }

            let watch_pos = 2 * conflict + num_watched;
            self.watched_literals_[watch_pos as usize].domchg = entries[i].clone();
            self.link_watched_literal(watch_pos, h16);
            num_watched += 1;
            if num_watched == 2 {
                break;
            }
        }

        match num_watched {
            0 => {
                let mut latest_active: [(HighsInt, usize); 2] = [(0, 0); 2];
                let mut num_active: HighsInt = 0;
                for i in 0..entries.len() {
                    let entry = &entries[i];
                    let pos = if entry.boundtype == HighsBoundType::kLower {
                        domain.col_lower_pos_[entry.column as usize]
                    } else {
                        domain.col_upper_pos_[entry.column as usize]
                    };
                    match num_active {
                        0 => {
                            latest_active[0] = (pos, i);
                            num_active = 1;
                        }
                        1 => {
                            latest_active[1] = (pos, i);
                            num_active = 2;
                            if latest_active[0].0 < latest_active[1].0 {
                                latest_active.swap(0, 1);
                            }
                        }
                        2 => {
                            if pos > latest_active[1].0 {
                                latest_active[1] = (pos, i);
                                if latest_active[0].0 < latest_active[1].0 {
                                    latest_active.swap(0, 1);
                                }
                            }
                        }
                        _ => {}
                    }
                }
                for i in 0..num_active {
                    let watch_pos = 2 * conflict + i;
                    self.watched_literals_[watch_pos as usize].domchg =
                        entries[latest_active[i as usize].1].clone();
                    self.link_watched_literal(watch_pos, h16);
                }
            }
            1 => {
                let mut latest_active: HighsInt = -1;
                let mut latest_pos: HighsInt = -1;
                for i in 0..entries.len() {
                    // HiGHS 1.16 #3571 (h16_watch_link): HighsDomain.cpp lines 242-244 (1.16.0): the inactive
                    // literal is already watched, so only active ones are considered
                    if h16 && !domain.is_active(&entries[i]) {
                        continue;
                    }
                    let entry = &entries[i];
                    let pos = if entry.boundtype == HighsBoundType::kLower {
                        domain.col_lower_pos_[entry.column as usize]
                    } else {
                        domain.col_upper_pos_[entry.column as usize]
                    };
                    if pos > latest_pos {
                        latest_active = i as HighsInt;
                        latest_pos = pos;
                    }
                }
                if latest_active != -1 {
                    let watch_pos = 2 * conflict + 1;
                    self.watched_literals_[watch_pos as usize].domchg =
                        entries[latest_active as usize].clone();
                    self.link_watched_literal(watch_pos, h16);
                }
            }
            _ => {}
        }

        self.conflict_flag_[ci] = (num_watched as u8) | (self.conflict_flag_[ci] & 4);
        self.mark_propagate_conflict(conflict);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::markPropagateConflict (lines 261-267)
    /// C++: `void markPropagateConflict(HighsInt conflict)`
    pub fn mark_propagate_conflict(&mut self, conflict: HighsInt) {
        let ci = conflict as usize;
        if self.conflict_flag_[ci] < 2 {
            self.propagate_conflict_inds_.push(conflict);
            self.conflict_flag_[ci] |= 4;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::updateActivityLbChange (lines 269-285)
    /// C++: `void updateActivityLbChange(HighsInt col, double oldbound, double newbound)`
    // PORT NOTE: `domain` only feeds the `assert(!domain->infeasible_)`.
    pub fn update_activity_lb_change(
        &mut self,
        domain: &HighsDomain,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
    ) {
        debug_assert!(!domain.infeasible_);

        let mut i = self.col_lower_watched_[col as usize];
        while i != -1 {
            let conflict = i >> 1;

            // PORT NOTE (perf): the literal is read once (boundval and next) and the flag is
            // indexed once, with `markPropagateConflict` applied to the updated flag in place.
            let literal = &self.watched_literals_[i as usize];
            let boundval = literal.domchg.boundval;
            let next = literal.next;
            // C++ uint8_t arithmetic: the int difference -1 wraps to 255, and `+=` wraps.
            let num_inactive_delta: u8 =
                (((boundval > newbound) as i32) - ((boundval > oldbound) as i32)) as u8;
            if num_inactive_delta != 0 {
                let flag = &mut self.conflict_flag_[conflict as usize];
                let updated = flag.wrapping_add(num_inactive_delta);
                if updated < 2 {
                    *flag = updated | 4;
                    self.propagate_conflict_inds_.push(conflict);
                } else {
                    *flag = updated;
                }
            }
            i = next;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictPoolPropagation::updateActivityUbChange (lines 287-303)
    /// C++: `void updateActivityUbChange(HighsInt col, double oldbound, double newbound)`
    pub fn update_activity_ub_change(
        &mut self,
        domain: &HighsDomain,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
    ) {
        debug_assert!(!domain.infeasible_);

        let mut i = self.col_upper_watched_[col as usize];
        while i != -1 {
            let conflict = i >> 1;

            // PORT NOTE (perf): the literal is read once (boundval and next) and the flag is
            // indexed once, with `markPropagateConflict` applied to the updated flag in place.
            let literal = &self.watched_literals_[i as usize];
            let boundval = literal.domchg.boundval;
            let next = literal.next;
            // C++ uint8_t arithmetic: the int difference -1 wraps to 255, and `+=` wraps.
            let num_inactive_delta: u8 =
                (((boundval < newbound) as i32) - ((boundval < oldbound) as i32)) as u8;
            if num_inactive_delta != 0 {
                let flag = &mut self.conflict_flag_[conflict as usize];
                let updated = flag.wrapping_add(num_inactive_delta);
                if updated < 2 {
                    *flag = updated | 4;
                    self.propagate_conflict_inds_.push(conflict);
                } else {
                    *flag = updated;
                }
            }
            i = next;
        }
    }
}

impl HighsDomain_CutpoolPropagation {
    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::CutpoolPropagation (lines 378-383)
    /// C++: `CutpoolPropagation(HighsInt cutpoolindex, HighsDomain * domain, HighsCutPool & cutpool_)`
    // PORT NOTE: signature changed: no back-pointer; the pool is `ms.mipdata_.cutpools[cutpool]`.
    // Registers a fresh key at the log end (`add_propagation_domain`).
    pub fn new(ms: &mut HighsMipSolver, cutpoolindex: HighsInt, cutpool: usize) -> Self {
        let key = DomainKey::fresh();
        ms.mipdata_.cutpools[cutpool].add_propagation_domain(key);
        let cursor = ms.mipdata_.cutpools[cutpool].events.end();
        Self {
            cutpoolindex,
            cutpool,
            cursor,
            key,
            ..Default::default()
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::CutpoolPropagation (lines 385-398)
    /// C++: `CutpoolPropagation(const CutpoolPropagation & other)`
    // PORT NOTE: as `ConflictPoolPropagation::copy_from` (fresh key, inherited cursor, the
    // `parallelLockActive() && cutpool == 0` exception sets `cursor = u64::MAX`).
    pub fn copy_from(ms: &mut HighsMipSolver, other: &HighsDomain_CutpoolPropagation) -> Self {
        let key = DomainKey::fresh();
        let mut cursor = u64::MAX;
        if !ms.parallel_lock_active() || other.cutpool != 0 {
            let pool_end = ms.mipdata_.cutpools[other.cutpool].events.end();
            cursor = if other.cursor == u64::MAX {
                pool_end
            } else {
                other.cursor
            };
            ms.mipdata_.cutpools[other.cutpool]
                .events
                .register(key, cursor);
        }
        Self {
            cutpoolindex: other.cutpoolindex,
            cutpool: other.cutpool,
            cursor,
            key,
            activitycuts_: other.activitycuts_.clone(),
            activitycutsinf_: other.activitycutsinf_.clone(),
            propagatecutflags_: other.propagatecutflags_.clone(),
            propagatecutinds_: other.propagatecutinds_.clone(),
            capacity_threshold_: other.capacity_threshold_.clone(),
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::operator= (lines 406-421)
    /// C++: `CutpoolPropagation & operator=(const CutpoolPropagation & other)`
    pub fn op_assign(&mut self, ms: &mut HighsMipSolver, other: &HighsDomain_CutpoolPropagation) {
        if std::ptr::eq(&*self, other) {
            return;
        }
        let key = self.key;
        ms.mipdata_.cutpools[self.cutpool].remove_propagation_domain(key);
        self.cutpoolindex = other.cutpoolindex;
        self.cutpool = other.cutpool;
        self.activitycuts_.clone_from(&other.activitycuts_);
        self.activitycutsinf_.clone_from(&other.activitycutsinf_);
        self.propagatecutflags_
            .clone_from(&other.propagatecutflags_);
        self.propagatecutinds_.clone_from(&other.propagatecutinds_);
        self.capacity_threshold_
            .clone_from(&other.capacity_threshold_);
        debug_assert!(!ms.parallel_lock_active());
        let pool_end = ms.mipdata_.cutpools[self.cutpool].events.end();
        self.cursor = if other.cursor == u64::MAX {
            pool_end
        } else {
            other.cursor
        };
        ms.mipdata_.cutpools[self.cutpool]
            .events
            .register(key, self.cursor);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::~CutpoolPropagation (lines 423-425)
    /// C++: `~CutpoolPropagation()`
    pub fn release(&mut self, ms: &mut HighsMipSolver) {
        ms.mipdata_.cutpools[self.cutpool].remove_propagation_domain(self.key);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::recomputeCapacityThreshold (lines 427-446)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::HighsDomain::CutpoolPropagation::recomputeCapacityThreshold (lines 434-453)
    /// C++: `void recomputeCapacityThreshold(HighsInt cut)`
    // PORT NOTE: signature changed: the cut row is passed as slices (`getRowStart..getRowEnd` of
    // the pool matrix in `propagate`, the event snapshot in `cut_added`); `ms` for
    // `domain->feastol()` / `variableType`.
    pub fn recompute_capacity_threshold(
        &mut self,
        ms: &HighsMipSolver,
        domain: &HighsDomain,
        cut: HighsInt,
        arindex: &[HighsInt],
        arvalue: &[f64],
    ) {
        let ci = cut as usize;
        self.capacity_threshold_[ci] = -domain.feastol(ms);
        for i in 0..arindex.len() {
            let cl = arindex[i] as usize;
            if domain.col_upper_[cl] == domain.col_lower_[cl] {
                continue;
            }

            let threshold = arvalue[i].abs()
                * bound_range(
                    domain.col_upper_[cl],
                    domain.col_lower_[cl],
                    domain.feastol(ms),
                    domain.variable_type(ms, arindex[i]),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );

            // C++ std::max({a, b, c}) is max_element with `<`: cpp_max(cpp_max(a, b), c).
            self.capacity_threshold_[ci] = cpp_max(
                cpp_max(self.capacity_threshold_[ci], threshold),
                domain.feastol(ms),
            );
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::cutAdded (lines 448-486)
    /// C++: `void cutAdded(HighsInt cut, bool propagate)`
    // PORT NOTE: signature changed: replayed from `PoolEvent::CutAdded` (row snapshot `inds`,
    // `vals`, `rhs`). `domain != &getDomain()` is `!domain.is_mip_global_domain(ms)`;
    // `computeMinActivity(start, end, ...)` becomes `compute_min_activity(0, len, inds, vals, ..)`;
    // `markPropagateCut(cut)` becomes `mark_propagate_cut(rhs, cut)`.
    pub fn cut_added(
        &mut self,
        ms: &HighsMipSolver,
        domain: &HighsDomain,
        cut: HighsInt,
        propagate: bool,
        inds: &[HighsInt],
        vals: &[f64],
        rhs: f64,
    ) {
        let ci = cut as usize;
        if !propagate {
            if !domain.is_mip_global_domain(ms) {
                return;
            }
        }

        if self.activitycuts_.len() as HighsInt <= cut {
            self.activitycuts_.resize(ci + 1, HighsCDouble::default());
            self.activitycutsinf_.resize(ci + 1, 0);
            self.propagatecutflags_.resize(ci + 1, 2);
            self.capacity_threshold_.resize(ci + 1, 0.0);
        }

        self.propagatecutflags_[ci] &= !2u8;
        domain.compute_min_activity(
            0,
            inds.len() as HighsInt,
            inds,
            vals,
            &mut self.activitycutsinf_[ci],
            &mut self.activitycuts_[ci],
        );

        if propagate {
            self.recompute_capacity_threshold(ms, domain, cut, inds, vals);
            self.mark_propagate_cut(rhs, cut, ms.options_mip_.h16_cdouble);
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::cutDeleted (lines 488-497)
    /// C++: `void cutDeleted(HighsInt cut, bool deletedOnlyForPropagation = false)`
    /// C++ default arguments: `deleted_only_for_propagation = false`
    // PORT NOTE: `domain == &getDomain()` is `domain.is_mip_global_domain(ms)`.
    pub fn cut_deleted(
        &mut self,
        ms: &HighsMipSolver,
        domain: &HighsDomain,
        cut: HighsInt,
        deleted_only_for_propagation: bool,
    ) {
        if deleted_only_for_propagation && domain.is_mip_global_domain(ms) {
            debug_assert!(domain.branch_pos_.is_empty());
            return;
        }

        if cut < self.propagatecutflags_.len() as HighsInt {
            self.propagatecutflags_[cut as usize] |= 2;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::markPropagateCut (lines 499-507)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::CutpoolPropagation::markPropagateCut (lines 507-515)
    /// C++: `void markPropagateCut(HighsInt cut)`
    // PORT NOTE: signature changed: `cutpool->getRhs()[cut]` is passed as `rhs` (the live pool
    // value from `update_activity_*` / `HighsDomain::mark_propagate_cut`, the event snapshot
    // from `cut_added`). `h16_cdouble` is `ms.options_mip_.h16_cdouble` (the callers have `ms`).
    pub fn mark_propagate_cut(&mut self, rhs: f64, cut: HighsInt, h16_cdouble: bool) {
        let ci = cut as usize;
        if self.propagatecutflags_[ci] == 0
            && (self.activitycutsinf_[ci] == 1 || {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 510 (1.16.0)
                if h16_cdouble {
                    rhs - self.activitycuts_[ci] <= self.capacity_threshold_[ci]
                } else {
                    rhs - f64::from(self.activitycuts_[ci]) <= self.capacity_threshold_[ci]
                }
            })
        {
            self.propagatecutinds_.push(cut);
            self.propagatecutflags_[ci] |= 1;
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::updateActivityLbChange (lines 509-572)
    /// C++: `void updateActivityLbChange(HighsInt col, double oldbound, double newbound, bool threshold, bool activity, bool infeasdomain)`
    // PORT NOTE: signature changed: `domain` is the owner (its `cutpoolpropagation` is taken out
    // by the caller); `cutpool->getMatrix().forEach*ColumnEntry` reads
    // `ms.mipdata_.cutpools[self.cutpool]`; the closures write `self.activitycuts_[row]` etc.,
    // `domain.infeasible_*` and read `ms.mipdata_.feastol`.
    pub fn update_activity_lb_change(
        &mut self,
        ms: &HighsMipSolver,
        domain: &mut HighsDomain,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
        threshold: bool,
        activity: bool,
        infeasdomain: bool,
    ) {
        if !infeasdomain {
            debug_assert!(!domain.infeasible_);
        }

        let cutpool = self.cutpool;
        if newbound < oldbound && threshold {
            ms.mipdata_.cutpools[cutpool]
                .get_matrix()
                .for_each_negative_column_entry(col, &mut |row: HighsInt, val: f64| {
                    domain.update_threshold_lb_change(
                        ms,
                        col,
                        newbound,
                        val,
                        &mut self.capacity_threshold_[row as usize],
                    );
                    true
                });
        }

        if !activity {
            return;
        }

        if !infeasdomain {
            ms.mipdata_.cutpools[cutpool]
                .get_matrix()
                .for_each_positive_column_entry(col, &mut |row: HighsInt, val: f64| {
                    debug_assert!(val > 0.0);
                    let ri = row as usize;
                    let deltamin: HighsCDouble = compute_delta(
                        val,
                        oldbound,
                        newbound,
                        -kHighsInf,
                        &mut self.activitycutsinf_[ri],
                    );
                    self.activitycuts_[ri] += deltamin;

                    if deltamin <= 0.0 {
                        domain.update_threshold_lb_change(
                            ms,
                            col,
                            newbound,
                            val,
                            &mut self.capacity_threshold_[ri],
                        );
                        return true;
                    }

                    if self.activitycutsinf_[ri] == 0
                        && self.activitycuts_[ri] - ms.mipdata_.cutpools[cutpool].get_rhs()[ri]
                            > ms.mipdata_.feastol
                    {
                        ms.mipdata_.debug_solution.node_pruned(&*domain);
                        domain.infeasible_ = true;
                        domain.infeasible_pos = domain.domchgstack_.len() as HighsInt;
                        domain.infeasible_reason = HighsDomain_Reason::cut(self.cutpoolindex, row);
                        return false;
                    }

                    let rhs = ms.mipdata_.cutpools[cutpool].get_rhs()[ri];
                    self.mark_propagate_cut(rhs, row, ms.options_mip_.h16_cdouble);

                    true
                });
        }

        if domain.infeasible_ {
            // assert(domain->infeasible_reason.type == cutpoolindex);
            debug_assert!(domain.infeasible_reason.index >= 0);
            let (oldbound, newbound) = (newbound, oldbound);
            ms.mipdata_.cutpools[cutpool]
                .get_matrix()
                .for_each_positive_column_entry(col, &mut |row: HighsInt, val: f64| {
                    debug_assert!(val > 0.0);
                    let ri = row as usize;
                    self.activitycuts_[ri] += compute_delta(
                        val,
                        oldbound,
                        newbound,
                        -kHighsInf,
                        &mut self.activitycutsinf_[ri],
                    );

                    if domain.infeasible_reason.index == row
                        && domain.infeasible_reason.r#type == self.cutpoolindex
                    {
                        return false;
                    }

                    true
                });
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::CutpoolPropagation::updateActivityUbChange (lines 574-638)
    /// C++: `void updateActivityUbChange(HighsInt col, double oldbound, double newbound, bool threshold, bool activity, bool infeasdomain)`
    pub fn update_activity_ub_change(
        &mut self,
        ms: &HighsMipSolver,
        domain: &mut HighsDomain,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
        threshold: bool,
        activity: bool,
        infeasdomain: bool,
    ) {
        if !infeasdomain {
            debug_assert!(!domain.infeasible_);
        }

        let cutpool = self.cutpool;
        if newbound > oldbound && threshold {
            ms.mipdata_.cutpools[cutpool]
                .get_matrix()
                .for_each_positive_column_entry(col, &mut |row: HighsInt, val: f64| {
                    domain.update_threshold_ub_change(
                        ms,
                        col,
                        newbound,
                        val,
                        &mut self.capacity_threshold_[row as usize],
                    );
                    true
                });
        }

        if !activity {
            return;
        }

        if !infeasdomain {
            ms.mipdata_.cutpools[cutpool]
                .get_matrix()
                .for_each_negative_column_entry(col, &mut |row: HighsInt, val: f64| {
                    debug_assert!(val < 0.0);
                    let ri = row as usize;
                    let deltamin: HighsCDouble = compute_delta(
                        val,
                        oldbound,
                        newbound,
                        kHighsInf,
                        &mut self.activitycutsinf_[ri],
                    );

                    self.activitycuts_[ri] += deltamin;

                    if deltamin <= 0.0 {
                        domain.update_threshold_ub_change(
                            ms,
                            col,
                            newbound,
                            val,
                            &mut self.capacity_threshold_[ri],
                        );
                        return true;
                    }

                    if self.activitycutsinf_[ri] == 0
                        && self.activitycuts_[ri] - ms.mipdata_.cutpools[cutpool].get_rhs()[ri]
                            > ms.mipdata_.feastol
                    {
                        ms.mipdata_.debug_solution.node_pruned(&*domain);
                        domain.infeasible_ = true;
                        domain.infeasible_pos = domain.domchgstack_.len() as HighsInt;
                        domain.infeasible_reason = HighsDomain_Reason::cut(self.cutpoolindex, row);
                        return false;
                    }

                    let rhs = ms.mipdata_.cutpools[cutpool].get_rhs()[ri];
                    self.mark_propagate_cut(rhs, row, ms.options_mip_.h16_cdouble);

                    true
                });
        }

        if domain.infeasible_ {
            // assert(domain->infeasible_reason.type == cutpoolindex);
            debug_assert!(domain.infeasible_reason.index >= 0);
            let (oldbound, newbound) = (newbound, oldbound);
            ms.mipdata_.cutpools[cutpool]
                .get_matrix()
                .for_each_negative_column_entry(col, &mut |row: HighsInt, val: f64| {
                    debug_assert!(val < 0.0);
                    let ri = row as usize;
                    self.activitycuts_[ri] += compute_delta(
                        val,
                        oldbound,
                        newbound,
                        kHighsInf,
                        &mut self.activitycutsinf_[ri],
                    );

                    if domain.infeasible_reason.index == row
                        && domain.infeasible_reason.r#type == self.cutpoolindex
                    {
                        return false;
                    }

                    true
                });
        }
    }
}

/// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation (lines 673-1258)
///
/// PORT NOTE: every `ObjectivePropagation` member function that reads the domain is a
/// `HighsDomain` method named `obj_prop_<name>` operating on `self.obj_prop_` (C++ `this`)
/// and `self` (C++ `domain`). `objFunc` is `ms.mipdata_.objective_function`, `cost` is
/// `mip_model!(ms).col_cost_`. Borrow the tree through the field, not through `self`:
/// `ObjectiveContributionTree::new(&mut self.obj_prop_, i)` next to reads of
/// `self.col_lower_` is a disjoint-field borrow and compiles; a `&self` helper call while
/// the tree is alive does not (read the values into locals first).
impl HighsDomain {
    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::ObjectivePropagation (lines 673-744)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::ObjectivePropagation::ObjectivePropagation (lines 921-994)
    /// C++: `ObjectivePropagation(HighsDomain * domain)`
    // PORT NOTE: builds `self.obj_prop_` in place (C++ `objProp_ = ObjectivePropagation(this)`;
    // the constructor reads only the domain bounds, not the old `objProp_`). Sets
    // `obj_prop_.active = true` (C++ `domain != nullptr`).
    pub fn obj_prop_init(&mut self, ms: &HighsMipSolver) {
        let objfunc = &ms.mipdata_.objective_function;
        let cost = &mip_model!(ms).col_cost_;
        let obj_nonzeros = objfunc.get_objective_nonzeros();
        let partition_starts = objfunc.get_clique_partition_starts();
        let num_partitions = objfunc.get_num_clique_partitions();

        // Built as a local so that the contribution tree (which borrows it) does not overlap `self`.
        let mut op = HighsDomain_ObjectivePropagation {
            active: true,
            is_propagated: false,
            capacity_threshold: kHighsInf,
            ..Default::default()
        };
        if num_partitions != 0 {
            op.propagation_cons_buffer = objfunc.get_objective_values_packed().to_vec();
            op.partition_clique_data = vec![Default::default(); num_partitions as usize];
        }
        op.objective_lower_contributions =
            vec![Default::default(); partition_starts[num_partitions as usize] as usize];
        op.contribution_partition_sets = vec![(-1, -1); num_partitions as usize];

        // for each clique partition first set up all columns as contributing with their largest
        // possible value and then remove the largest contribution of each clique partition. For these
        // columns we do not need to check for finite bounds as they must be binary.
        let mut objective_lower = HighsCDouble::from(0.0);
        for i in 0..num_partitions {
            let mut rhs: HighsInt = 1;
            {
                let mut tree =
                    HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(&mut op, i);
                for j in partition_starts[i as usize]..partition_starts[i as usize + 1] {
                    let cl = obj_nonzeros[j as usize];
                    tree.nodes[j as usize].col = cl;
                    tree.nodes[j as usize].partition = i;
                    if cost[cl as usize] > 0.0 {
                        objective_lower += cost[cl as usize];
                        tree.nodes[j as usize].contribution = cost[cl as usize];
                        rhs -= 1;
                        if self.col_lower_[cl as usize] == 0.0 {
                            tree.link(j);
                        }
                    } else {
                        tree.nodes[j as usize].contribution = -cost[cl as usize];
                        if self.col_upper_[cl as usize] == 1.0 {
                            tree.link(j);
                        }
                    }
                }

                let worst_pos = tree.first();
                if worst_pos != -1 {
                    objective_lower -= tree.nodes[worst_pos as usize].contribution;
                }
            }
            op.partition_clique_data[i as usize].rhs = rhs;
        }

        // add contribution of remaining objective nonzeros
        let num_obj_nz = obj_nonzeros.len();
        for i in partition_starts[num_partitions as usize] as usize..num_obj_nz {
            let cl = obj_nonzeros[i] as usize;
            if cost[cl] > 0.0 {
                if self.col_lower_[cl] == -kHighsInf {
                    op.num_inf_obj_lower += 1;
                } else if ms.options_mip_.h16_cdouble {
                    // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 982 (1.16.0)
                    objective_lower += HighsCDouble::from(self.col_lower_[cl]) * cost[cl];
                } else {
                    objective_lower += self.col_lower_[cl] * cost[cl];
                }
            } else if self.col_upper_[cl] == kHighsInf {
                op.num_inf_obj_lower += 1;
            } else if ms.options_mip_.h16_cdouble {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 988 (1.16.0)
                objective_lower += HighsCDouble::from(self.col_upper_[cl]) * cost[cl];
            } else {
                objective_lower += self.col_upper_[cl] * cost[cl];
            }
        }

        op.objective_lower = objective_lower;
        self.obj_prop_ = op;
        self.obj_prop_recompute_capacity_threshold(ms);
        self.obj_prop_debug_check_objective_lower(ms);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::getPropagationConstraint (lines 746-791)
    /// C++: `void getPropagationConstraint(HighsInt domchgStackSize, const double *& vals, const HighsInt *& inds, HighsInt & len, double & rhs, HighsInt domchgCol = - 1)`
    /// C++ default arguments: `domchg_col = - 1`
    // PORT NOTE: the out-slices borrow either `ms` (`objFunc` data) or `self.obj_prop_`
    // (`propagationConsBuffer`), so both are borrowed for `'s`. Callers that keep using the
    // domain copy them out (`to_vec()`), as `extract_obj_cliques` does.
    pub fn obj_prop_get_propagation_constraint<'s>(
        &'s mut self,
        ms: &'s HighsMipSolver,
        domchg_stack_size: HighsInt,
        vals: &mut &'s [f64],
        inds: &mut &'s [HighsInt],
        len: &mut HighsInt,
        rhs: &mut f64,
        domchg_col: HighsInt,
    ) {
        // HiGHS #3223 (mip_cdouble_proofs): HighsDomain.cpp line 1034 at 7db1c8180
        let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;

        let num_partitions = ms.mipdata_.objective_function.get_num_clique_partitions();
        *inds = ms.mipdata_.objective_function.get_objective_nonzeros();
        *len = inds.len() as HighsInt;
        if num_partitions == 0 {
            *vals = ms.mipdata_.objective_function.get_objective_values_packed();
            *rhs = ms.mipdata_.upper_limit;
            return;
        }
        let cost = &mip_model!(ms).col_cost_;
        let partition_starts = ms.mipdata_.objective_function.get_clique_partition_starts();
        let packed_obj_vals = ms.mipdata_.objective_function.get_objective_values_packed();

        let mut tmp_rhs = HighsCDouble::from(ms.mipdata_.upper_limit);
        for i in 0..num_partitions as usize {
            let start = partition_starts[i];
            let end = partition_starts[i + 1];
            let mut largest = 0.0f64;
            for j in start..end {
                let cl = inds[j as usize];
                // skip the column we might want to explain a bound change for and take the second
                // largest column instead.
                if cl == domchg_col {
                    continue;
                }
                let mut pos: HighsInt = 0;
                if cost[cl as usize] > 0.0 {
                    let lb = self.get_col_lower_pos(cl, domchg_stack_size, &mut pos);
                    if lb < 1.0 {
                        largest = cpp_max(largest, cost[cl as usize]);
                    }
                } else {
                    let ub = self.get_col_upper_pos(cl, domchg_stack_size, &mut pos);
                    if ub > 0.0 {
                        largest = cpp_max(largest, -cost[cl as usize]);
                    }
                }
            }

            if x1b_cdouble {
                tmp_rhs += HighsCDouble::from(largest)
                    * (self.obj_prop_.partition_clique_data[i].rhs as f64);
            } else {
                tmp_rhs += largest * (self.obj_prop_.partition_clique_data[i].rhs as f64);
            }
            if self.obj_prop_.partition_clique_data[i].multiplier != largest {
                self.obj_prop_.partition_clique_data[i].multiplier = largest;
                for j in start..end {
                    let pv = packed_obj_vals[j as usize];
                    self.obj_prop_.propagation_cons_buffer[j as usize] = pv - largest.copysign(pv);
                }
            }
        }

        *vals = self.obj_prop_.propagation_cons_buffer.as_slice();
        *rhs = f64::from(tmp_rhs);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::recomputeCapacityThreshold (lines 793-824)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::ObjectivePropagation::recomputeCapacityThreshold (lines 1043-1074)
    /// C++: `void recomputeCapacityThreshold()`
    pub fn obj_prop_recompute_capacity_threshold(&mut self, ms: &HighsMipSolver) {
        let objfunc = &ms.mipdata_.objective_function;
        let cost = &mip_model!(ms).col_cost_;
        let obj_nonzeros = objfunc.get_objective_nonzeros();
        let partition_starts = objfunc.get_clique_partition_starts();
        let num_partitions = objfunc.get_num_clique_partitions();
        let feastol = self.feastol(ms);

        self.obj_prop_.capacity_threshold = -feastol;
        for i in 0..num_partitions {
            let worst_pos = obj_tree_first(&mut self.obj_prop_, i);
            if worst_pos == -1 {
                continue;
            }
            let worst_col = self.obj_prop_.objective_lower_contributions[worst_pos as usize].col;
            if self.is_fixed(worst_col) {
                continue;
            }

            let mut contribution =
                self.obj_prop_.objective_lower_contributions[worst_pos as usize].contribution;
            let best_pos = obj_tree_last(&mut self.obj_prop_, i);
            if best_pos != worst_pos {
                contribution -=
                    self.obj_prop_.objective_lower_contributions[best_pos as usize].contribution;
            }

            self.obj_prop_.capacity_threshold = cpp_max(
                self.obj_prop_.capacity_threshold,
                contribution * (1.0 - feastol),
            );
        }

        for i in partition_starts[num_partitions as usize] as usize..obj_nonzeros.len() {
            let cl = obj_nonzeros[i] as usize;
            // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1071 (1.16.0)
            let range = bound_range(
                self.col_upper_[cl],
                self.col_lower_[cl],
                feastol,
                self.variable_type(ms, cl as HighsInt),
                ms.options_mip_.h16_bound_range_inf,
                ms.options_mip_.mip_cont_prop_improve,
            );
            self.obj_prop_.capacity_threshold =
                cpp_max(self.obj_prop_.capacity_threshold, cost[cl].abs() * range);
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::updateActivityLbChange (lines 826-945)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::ObjectivePropagation::updateActivityLbChange (lines 1076-1195)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::ObjectivePropagation::updateActivityLbChange (lines 1076-1195)
    /// C++: `void updateActivityLbChange(HighsInt col, double oldbound, double newbound)`
    pub fn obj_prop_update_activity_lb_change(
        &mut self,
        ms: &HighsMipSolver,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
    ) {
        let cl = col as usize;
        let cost = mip_model!(ms).col_cost_[cl];
        if cost <= 0.0 {
            if cost != 0.0 && newbound < oldbound {
                let feastol = self.feastol(ms);
                // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1082 (1.16.0)
                let range = bound_range(
                    self.col_upper_[cl],
                    newbound,
                    feastol,
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );
                self.obj_prop_.capacity_threshold =
                    cpp_max(self.obj_prop_.capacity_threshold, -cost * range);
                self.obj_prop_.is_propagated = false;
            }
            self.obj_prop_debug_check_objective_lower(ms);
            return;
        }

        self.obj_prop_.is_propagated = false;

        let upper_limit = ms.mipdata_.upper_limit;
        let partition_pos = ms.mipdata_.objective_function.get_col_clique_partition(col);
        if partition_pos == -1 {
            if oldbound == -kHighsInf {
                self.obj_prop_.num_inf_obj_lower -= 1;
            } else if ms.options_mip_.h16_cdouble {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1098 (1.16.0)
                self.obj_prop_.objective_lower -= HighsCDouble::from(oldbound) * cost;
            } else {
                self.obj_prop_.objective_lower -= oldbound * cost;
            }

            if newbound == -kHighsInf {
                self.obj_prop_.num_inf_obj_lower += 1;
            } else if ms.options_mip_.h16_cdouble {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1103 (1.16.0)
                self.obj_prop_.objective_lower += HighsCDouble::from(newbound) * cost;
            } else {
                self.obj_prop_.objective_lower += newbound * cost;
            }

            self.obj_prop_debug_check_objective_lower(ms);

            if newbound < oldbound {
                let feastol = self.feastol(ms);
                // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1110 (1.16.0)
                let range = bound_range(
                    self.col_upper_[cl],
                    self.col_lower_[cl],
                    feastol,
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );
                self.obj_prop_.capacity_threshold =
                    cpp_max(self.obj_prop_.capacity_threshold, cost * range);
            } else if self.obj_prop_.num_inf_obj_lower == 0
                && self.obj_prop_.objective_lower > upper_limit
            {
                self.infeasible_ = true;
                self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                self.infeasible_reason = HighsDomain_Reason::objective();
                self.obj_prop_update_activity_lb_change(ms, col, newbound, oldbound);
            }
        } else if newbound == 0.0 {
            // binary lower bound of variable in clique partition is relaxed to 0
            debug_assert!(oldbound == 1.0);
            self.obj_prop_relax_clique_column(ms, partition_pos);
        } else {
            // binary lower bound of variable in clique partition is tightened to 1
            debug_assert!(oldbound == 0.0);
            debug_assert!(newbound == 1.0);
            self.obj_prop_tighten_clique_column(ms, partition_pos);

            if self.obj_prop_.num_inf_obj_lower == 0 && self.obj_prop_.objective_lower > upper_limit
            {
                self.infeasible_ = true;
                self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                self.infeasible_reason = HighsDomain_Reason::objective();
                self.obj_prop_update_activity_lb_change(ms, col, newbound, oldbound);
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::updateActivityUbChange (lines 947-1064)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::ObjectivePropagation::updateActivityUbChange (lines 1197-1314)
    /// HiGHS 1.16 (#3525): HighsDomain.cpp::ObjectivePropagation::updateActivityUbChange (lines 1197-1314)
    /// C++: `void updateActivityUbChange(HighsInt col, double oldbound, double newbound)`
    pub fn obj_prop_update_activity_ub_change(
        &mut self,
        ms: &HighsMipSolver,
        col: HighsInt,
        oldbound: f64,
        newbound: f64,
    ) {
        let cl = col as usize;
        let cost = mip_model!(ms).col_cost_[cl];
        if cost >= 0.0 {
            if cost != 0.0 && newbound > oldbound {
                let feastol = self.feastol(ms);
                // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1203 (1.16.0)
                let range = bound_range(
                    newbound,
                    self.col_lower_[cl],
                    feastol,
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );
                self.obj_prop_.capacity_threshold =
                    cpp_max(self.obj_prop_.capacity_threshold, cost * range);
                self.obj_prop_.is_propagated = false;
            }
            self.obj_prop_debug_check_objective_lower(ms);
            return;
        }

        self.obj_prop_.is_propagated = false;

        let upper_limit = ms.mipdata_.upper_limit;
        let partition_pos = ms.mipdata_.objective_function.get_col_clique_partition(col);
        if partition_pos == -1 {
            if oldbound == kHighsInf {
                self.obj_prop_.num_inf_obj_lower -= 1;
            } else if ms.options_mip_.h16_cdouble {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1218 (1.16.0)
                self.obj_prop_.objective_lower -= HighsCDouble::from(oldbound) * cost;
            } else {
                self.obj_prop_.objective_lower -= oldbound * cost;
            }

            if newbound == kHighsInf {
                self.obj_prop_.num_inf_obj_lower += 1;
            } else if ms.options_mip_.h16_cdouble {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1223 (1.16.0)
                self.obj_prop_.objective_lower += HighsCDouble::from(newbound) * cost;
            } else {
                self.obj_prop_.objective_lower += newbound * cost;
            }

            self.obj_prop_debug_check_objective_lower(ms);

            if newbound > oldbound {
                let feastol = self.feastol(ms);
                // HiGHS 1.16 #3525 (h16_bound_range_inf): HighsDomain.cpp line 1230 (1.16.0)
                let range = bound_range(
                    self.col_upper_[cl],
                    self.col_lower_[cl],
                    feastol,
                    self.variable_type(ms, col),
                    ms.options_mip_.h16_bound_range_inf,
                    ms.options_mip_.mip_cont_prop_improve,
                );
                self.obj_prop_.capacity_threshold =
                    cpp_max(self.obj_prop_.capacity_threshold, -cost * range);
            } else if self.obj_prop_.num_inf_obj_lower == 0
                && self.obj_prop_.objective_lower > upper_limit
            {
                self.infeasible_ = true;
                self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                self.infeasible_reason = HighsDomain_Reason::objective();
                self.obj_prop_update_activity_ub_change(ms, col, newbound, oldbound);
            }
        } else if newbound == 1.0 {
            // binary upper bound of variable in clique partition is relaxed to 1
            debug_assert!(oldbound == 0.0);
            self.obj_prop_relax_clique_column(ms, partition_pos);
        } else {
            // binary upper bound of variable in clique partition is tightened to 0
            debug_assert!(oldbound == 1.0);
            debug_assert!(newbound == 0.0);
            self.obj_prop_tighten_clique_column(ms, partition_pos);

            if self.obj_prop_.num_inf_obj_lower == 0 && self.obj_prop_.objective_lower > upper_limit
            {
                self.infeasible_ = true;
                self.infeasible_pos = self.domchgstack_.len() as HighsInt;
                self.infeasible_reason = HighsDomain_Reason::objective();
                self.obj_prop_update_activity_ub_change(ms, col, newbound, oldbound);
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::shouldBePropagated (lines 1066-1075)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::ObjectivePropagation::shouldBePropagated (lines 1316-1325)
    /// C++: `bool shouldBePropagated() const`
    pub fn obj_prop_should_be_propagated(&self, ms: &HighsMipSolver) -> bool {
        if self.obj_prop_.is_propagated {
            return false;
        }
        if self.obj_prop_.num_inf_obj_lower > 1 {
            return false;
        }
        if self.infeasible_ {
            return false;
        }
        let upper_limit = ms.mipdata_.upper_limit;
        if upper_limit == kHighsInf {
            return false;
        }
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1322 (1.16.0)
        let over_capacity = if ms.options_mip_.h16_cdouble {
            upper_limit - self.obj_prop_.objective_lower > self.obj_prop_.capacity_threshold
        } else {
            upper_limit - f64::from(self.obj_prop_.objective_lower)
                > self.obj_prop_.capacity_threshold
        };
        if over_capacity {
            return false;
        }

        true
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::debugCheckObjectiveLower (lines 1077-1122)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::ObjectivePropagation::debugCheckObjectiveLower (lines 1327-1375)
    /// C++: `void debugCheckObjectiveLower() const`
    pub fn obj_prop_debug_check_objective_lower(&self, ms: &HighsMipSolver) {
        if !cfg!(debug_assertions) || self.infeasible_ {
            return;
        }
        let objfunc = &ms.mipdata_.objective_function;
        let cost = &mip_model!(ms).col_cost_;
        let num_partitions = objfunc.get_num_clique_partitions();
        let partition_starts = objfunc.get_clique_partition_starts();
        let obj_nonzeros = objfunc.get_objective_nonzeros();

        let mut lower_from_scratch = HighsCDouble::from(0.0);
        let mut num_inf: HighsInt = 0;
        for i in 0..num_partitions as usize {
            let start = partition_starts[i];
            let end = partition_starts[i + 1];
            let mut largest = 0.0f64;
            for j in start..end {
                let cl = obj_nonzeros[j as usize] as usize;
                if cost[cl] > 0.0 {
                    lower_from_scratch += cost[cl];

                    if self.col_lower_[cl] < 1.0 {
                        largest = cpp_max(largest, cost[cl]);
                    }
                } else if self.col_upper_[cl] > 0.0 {
                    largest = cpp_max(largest, -cost[cl]);
                }
            }
            lower_from_scratch -= largest;
        }

        for i in partition_starts[num_partitions as usize] as usize..obj_nonzeros.len() {
            let cl = obj_nonzeros[i] as usize;
            if cost[cl] > 0.0 {
                if self.col_lower_[cl] > -kHighsInf {
                    // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1360 (1.16.0)
                    if ms.options_mip_.h16_cdouble {
                        lower_from_scratch += HighsCDouble::from(self.col_lower_[cl]) * cost[cl];
                    } else {
                        lower_from_scratch += self.col_lower_[cl] * cost[cl];
                    }
                } else {
                    num_inf += 1;
                }
            } else if self.col_upper_[cl] < kHighsInf {
                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1367 (1.16.0)
                if ms.options_mip_.h16_cdouble {
                    lower_from_scratch += HighsCDouble::from(self.col_upper_[cl]) * cost[cl];
                } else {
                    lower_from_scratch += self.col_upper_[cl] * cost[cl];
                }
            } else {
                num_inf += 1;
            }
        }
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 1372 (1.16.0)
        if ms.options_mip_.h16_cdouble {
            assert!(
                (lower_from_scratch - self.obj_prop_.objective_lower).abs() <= self.feastol(ms)
            );
        } else {
            assert!(
                f64::from(lower_from_scratch - self.obj_prop_.objective_lower).abs()
                    <= self.feastol(ms)
            );
        }
        assert_eq!(num_inf, self.obj_prop_.num_inf_obj_lower);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::propagate (lines 1124-1258)
    /// C++: `void propagate()`
    // PORT NOTE: `domain->changeBound` / `checkChangeBound` update `self.obj_prop_` (through
    // `obj_prop_update_activity_*`) while this loop runs: re-read every `self.obj_prop_`
    // field after each bound change, never hold a reference (or an `ObjectiveContributionTree`)
    // across one. Copy `objFunc` arrays you iterate (`getObjectiveNonzeros`,
    // `getCliquePartitionStarts`) into locals or re-index `ms.mipdata_.objective_function`.
    pub fn obj_prop_propagate(&mut self, ms: &mut HighsMipSolver) {
        // Reporting only: tally runs, bound changes and cutoffs.
        if !self.obj_prop_should_be_propagated(ms) {
            return;
        }
        let n0 = self.domchgstack_.len();
        let was_infeasible = self.infeasible_;
        self.obj_prop_propagate_body(ms);
        crate::work::add_event(crate::work::Event::ObjPropRuns, 1);
        crate::work::add_event(
            crate::work::Event::ObjPropBoundChanges,
            self.domchgstack_.len().saturating_sub(n0) as u64,
        );
        if self.infeasible_ && !was_infeasible {
            crate::work::add_event(crate::work::Event::ObjPropCutoffs, 1);
        }
    }

    fn obj_prop_propagate_body(&mut self, ms: &mut HighsMipSolver) {
        if !self.obj_prop_should_be_propagated(ms) {
            return;
        }

        self.obj_prop_debug_check_objective_lower(ms);

        let upper_limit = ms.mipdata_.upper_limit;
        if self.obj_prop_.num_inf_obj_lower == 0 && self.obj_prop_.objective_lower > upper_limit {
            self.infeasible_ = true;
            self.infeasible_pos = self.domchgstack_.len() as HighsInt;
            self.infeasible_reason = HighsDomain_Reason::objective();
            return;
        }

        // HiGHS #3223 (mip_cdouble_proofs): HighsDomain.cpp lines 1495-1507 at 7db1c8180
        let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;

        let mut capacity: HighsCDouble = upper_limit - self.obj_prop_.objective_lower;
        if self.obj_prop_.num_inf_obj_lower == 1 {
            // Scan non-binary columns for infinite bound contribution until the one column that
            // contributes with an infinite bound is found, which is the only column that can be propagated.
            let num_col = ms
                .mipdata_
                .objective_function
                .get_objective_nonzeros()
                .len() as HighsInt;
            let start = ms
                .mipdata_
                .objective_function
                .get_num_binaries_in_objective();
            for i in start..num_col {
                let cl = ms.mipdata_.objective_function.get_objective_nonzeros()[i as usize];
                let cost = mip_model!(ms).col_cost_[cl as usize];
                if (cost > 0.0 && self.col_lower_[cl as usize] != -kHighsInf)
                    || (cost < 0.0 && self.col_upper_[cl as usize] != kHighsInf)
                {
                    continue;
                }
                let boundtype = if cost > 0.0 {
                    HighsBoundType::kUpper
                } else {
                    HighsBoundType::kLower
                };
                self.check_change_bound(
                    ms,
                    boundtype,
                    cl,
                    capacity / cost,
                    HighsDomain_Reason::objective(),
                );
                break;
            }
        } else {
            let num_partitions = ms.mipdata_.objective_function.get_num_clique_partitions();
            loop {
                let mut num_bound_changes: HighsInt = 0;
                for i in 0..num_partitions {
                    let worst = obj_tree_first(&mut self.obj_prop_, i);
                    if worst == -1 {
                        continue;
                    }

                    let mut contribution =
                        self.obj_prop_.objective_lower_contributions[worst as usize].contribution;
                    let second_worst = obj_tree_successor(&mut self.obj_prop_, i, worst);
                    if second_worst != -1 {
                        contribution -= self.obj_prop_.objective_lower_contributions
                            [second_worst as usize]
                            .contribution;
                    }

                    // the upper limit already uses a tolerance, so we can do a hard cutoff
                    if contribution > capacity {
                        let cl = self.obj_prop_.objective_lower_contributions[worst as usize].col;
                        if mip_model!(ms).col_cost_[cl as usize] > 0.0 {
                            if self.col_upper_[cl as usize] > 0.0 {
                                num_bound_changes += 1;
                                self.change_bound_2(
                                    ms,
                                    HighsBoundType::kUpper,
                                    cl,
                                    0.0,
                                    HighsDomain_Reason::objective(),
                                );
                                if self.infeasible_ {
                                    break;
                                }
                            }
                        } else if self.col_lower_[cl as usize] < 1.0 {
                            num_bound_changes += 1;
                            self.change_bound_2(
                                ms,
                                HighsBoundType::kLower,
                                cl,
                                1.0,
                                HighsDomain_Reason::objective(),
                            );
                            if self.infeasible_ {
                                break;
                            }
                        }
                    } else if second_worst != -1 {
                        // it might be that we can fix the column with the lowest possible contribution
                        // to the bound value that yields the highest objective contribution.
                        let mut best = obj_tree_last(&mut self.obj_prop_, i);
                        while best != obj_tree_first(&mut self.obj_prop_, i) {
                            let first = obj_tree_first(&mut self.obj_prop_, i);
                            if self.obj_prop_.objective_lower_contributions[first as usize]
                                .contribution
                                - self.obj_prop_.objective_lower_contributions[best as usize]
                                    .contribution
                                > capacity
                            {
                                let cl =
                                    self.obj_prop_.objective_lower_contributions[best as usize].col;
                                if mip_model!(ms).col_cost_[cl as usize] > 0.0 {
                                    debug_assert!(self.col_lower_[cl as usize] < 1.0);
                                    num_bound_changes += 1;
                                    self.change_bound_2(
                                        ms,
                                        HighsBoundType::kLower,
                                        cl,
                                        1.0,
                                        HighsDomain_Reason::objective(),
                                    );
                                    if self.infeasible_ {
                                        break;
                                    }
                                } else {
                                    debug_assert!(self.col_upper_[cl as usize] > 0.0);
                                    num_bound_changes += 1;
                                    self.change_bound_2(
                                        ms,
                                        HighsBoundType::kUpper,
                                        cl,
                                        0.0,
                                        HighsDomain_Reason::objective(),
                                    );
                                    if self.infeasible_ {
                                        break;
                                    }
                                }
                            } else {
                                break;
                            }

                            best = obj_tree_last(&mut self.obj_prop_, i);
                        }

                        if self.infeasible_ {
                            break;
                        }
                    }
                }

                if self.infeasible_ {
                    break;
                }

                let num_obj_nzs = ms
                    .mipdata_
                    .objective_function
                    .get_objective_nonzeros()
                    .len();
                let start = ms.mipdata_.objective_function.get_clique_partition_starts()
                    [num_partitions as usize] as usize;
                for i in start..num_obj_nzs {
                    let cl = ms.mipdata_.objective_function.get_objective_nonzeros()[i] as usize;
                    let cost = mip_model!(ms).col_cost_[cl];
                    if cost > 0.0 {
                        if x1b_cdouble {
                            if self.check_change_bound(
                                ms,
                                HighsBoundType::kUpper,
                                cl as HighsInt,
                                (capacity + HighsCDouble::from(self.col_lower_[cl]) * cost) / cost,
                                HighsDomain_Reason::objective(),
                            ) {
                                num_bound_changes += 1;
                            }
                        } else {
                            if self.check_change_bound(
                                ms,
                                HighsBoundType::kUpper,
                                cl as HighsInt,
                                (capacity + self.col_lower_[cl] * cost) / cost,
                                HighsDomain_Reason::objective(),
                            ) {
                                num_bound_changes += 1;
                            }
                        }
                    } else {
                        if x1b_cdouble {
                            if self.check_change_bound(
                                ms,
                                HighsBoundType::kLower,
                                cl as HighsInt,
                                (capacity + HighsCDouble::from(self.col_upper_[cl]) * cost) / cost,
                                HighsDomain_Reason::objective(),
                            ) {
                                num_bound_changes += 1;
                            }
                        } else {
                            if self.check_change_bound(
                                ms,
                                HighsBoundType::kLower,
                                cl as HighsInt,
                                (capacity + self.col_upper_[cl] * cost) / cost,
                                HighsDomain_Reason::objective(),
                            ) {
                                num_bound_changes += 1;
                            }
                        }
                    }
                    if self.infeasible_ {
                        break;
                    }
                }
                if self.infeasible_ {
                    break;
                }

                if num_bound_changes == 0 {
                    break;
                }
                capacity = upper_limit - self.obj_prop_.objective_lower;
            }
        }

        self.obj_prop_recompute_capacity_threshold(ms);
        self.obj_prop_.is_propagated = true;
    }
}

/// C++ `ObjectiveContributionTree contributionTree(this, partition)` queries. The C++ view is a
/// pair of references into `contributionPartitionSets[partition]` plus `nodes`, so a fresh view
/// observes the same state as the C++ object held across a `changeBound` (see `propagate`).
fn obj_tree_first(
    obj_prop: &mut HighsDomain_ObjectivePropagation,
    partition: HighsInt,
) -> HighsInt {
    HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(obj_prop, partition).first()
}

fn obj_tree_last(obj_prop: &mut HighsDomain_ObjectivePropagation, partition: HighsInt) -> HighsInt {
    HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(obj_prop, partition).last()
}

fn obj_tree_successor(
    obj_prop: &mut HighsDomain_ObjectivePropagation,
    partition: HighsInt,
    node: HighsInt,
) -> HighsInt {
    HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(obj_prop, partition)
        .successor(node)
}

impl HighsDomain {
    /// The "binary bound relaxed to its non-restrictive value" branch of
    /// `ObjectivePropagation::updateActivityLbChange` (lines 862-890) and
    /// `updateActivityUbChange` (lines 1000-1024), which are textually identical in C++.
    fn obj_prop_relax_clique_column(&mut self, ms: &HighsMipSolver, partition_pos: HighsInt) {
        let feastol = self.feastol(ms);
        let partition =
            self.obj_prop_.objective_lower_contributions[partition_pos as usize].partition;
        let (old_contribution, is_first, contribution, best_pos) = {
            let mut tree = HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(
                &mut self.obj_prop_,
                partition,
            );
            let curr_first = tree.first();

            tree.link(partition_pos);

            let old_contribution = if curr_first != -1 {
                tree.nodes[curr_first as usize].contribution
            } else {
                0.0
            };
            let is_first = partition_pos == tree.first();
            let contribution = tree.nodes[partition_pos as usize].contribution;
            let best_pos = if is_first && contribution != old_contribution {
                tree.last()
            } else {
                -1
            };
            (old_contribution, is_first, contribution, best_pos)
        };

        if is_first && contribution != old_contribution {
            self.obj_prop_.objective_lower += old_contribution;
            self.obj_prop_.objective_lower -= contribution;

            // update the capacity threshold with the difference of the new highest
            // contribution position to the lowest contribution as the column with
            // the lowest contribution can be fixed to its bound that yields the
            // highest objective value.
            let mut delta = contribution;
            if best_pos != partition_pos {
                delta -=
                    self.obj_prop_.objective_lower_contributions[best_pos as usize].contribution;
            }
            self.obj_prop_.capacity_threshold =
                cpp_max(delta * (1.0 - feastol), self.obj_prop_.capacity_threshold);
        } else {
            // the new linked column could be the one with the new lowest contribution so update
            // the capacity threshold to ensure propagation runs when it can be fixed to the bound
            // that yields the highest objective value
            self.obj_prop_.capacity_threshold = cpp_max(
                (old_contribution - contribution) * (1.0 - feastol),
                self.obj_prop_.capacity_threshold,
            );
        }

        self.obj_prop_debug_check_objective_lower(ms);
    }

    /// The "binary bound tightened to its restrictive value" branch of
    /// `ObjectivePropagation::updateActivityLbChange` (lines 893-919) and
    /// `updateActivityUbChange` (lines 1027-1052), which are textually identical in C++.
    fn obj_prop_tighten_clique_column(&mut self, ms: &HighsMipSolver, partition_pos: HighsInt) {
        let partition =
            self.obj_prop_.objective_lower_contributions[partition_pos as usize].partition;
        let was_first = {
            let tree = HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(
                &mut self.obj_prop_,
                partition,
            );
            tree.first() == partition_pos
        };
        if was_first {
            let contrib =
                self.obj_prop_.objective_lower_contributions[partition_pos as usize].contribution;
            self.obj_prop_.objective_lower += contrib;
        }

        {
            let mut tree = HighsDomain_ObjectivePropagation_ObjectiveContributionTree::new(
                &mut self.obj_prop_,
                partition,
            );
            tree.unlink(partition_pos);
        }

        if was_first {
            let new_worst = obj_tree_first(&mut self.obj_prop_, partition);
            if new_worst != -1 {
                let contrib =
                    self.obj_prop_.objective_lower_contributions[new_worst as usize].contribution;
                self.obj_prop_.objective_lower -= contrib;
            }
        }

        self.obj_prop_debug_check_objective_lower(ms);
    }
}

#[cfg(test)]
mod h16_watch_link_tests {
    use crate::prelude::*;

    /// Two lower-bound watches on column 0 (literals 0 and 1), no list linked yet.
    fn pool_with_two_literals_on_col0() -> HighsDomain_ConflictPoolPropagation {
        let mut prop = HighsDomain_ConflictPoolPropagation {
            col_lower_watched_: vec![-1; 2],
            col_upper_watched_: vec![-1; 2],
            watched_literals_: vec![Default::default(); 4],
            ..Default::default()
        };
        for lp in 0..2 {
            prop.watched_literals_[lp].domchg = HighsDomainChange {
                boundval: 1.0,
                column: 0,
                boundtype: HighsBoundType::kLower,
            };
        }
        prop
    }

    /// Links literal 0 into the empty list of column 0, then literal 1 on top of it.
    fn link_two(h16: bool) -> HighsDomain_ConflictPoolPropagation {
        let mut prop = pool_with_two_literals_on_col0();
        prop.link_watched_literal(0, h16);
        prop.link_watched_literal(1, h16);
        prop
    }

    #[test]
    fn h16_watch_link_on_empty_list_head_points_to_literal() {
        // 1.16.0 HighsDomain.cpp lines 105-106: head = linkPos is unconditional.
        let mut prop = pool_with_two_literals_on_col0();
        prop.link_watched_literal(0, true);
        assert_eq!(prop.col_lower_watched_[0], 0);
        assert_eq!(prop.watched_literals_[0].prev, -1);
        assert_eq!(prop.watched_literals_[0].next, -1);
    }

    #[test]
    fn h16_watch_link_on_second_literal_becomes_head() {
        let prop = link_two(true);
        assert_eq!(prop.col_lower_watched_[0], 1);
        assert_eq!(prop.watched_literals_[1].next, 0);
        assert_eq!(prop.watched_literals_[0].prev, 1);
        assert_eq!(prop.watched_literals_[0].next, -1);
    }

    #[test]
    fn h16_watch_link_off_empty_list_stays_unreachable() {
        // 1.15.1 placement: head = linkPos only inside `if (head != -1)`, so the literal is
        // not reachable from the column head.
        let mut prop = pool_with_two_literals_on_col0();
        prop.link_watched_literal(0, false);
        assert_eq!(prop.col_lower_watched_[0], -1);
        assert_eq!(prop.watched_literals_[0].prev, -1);
        assert_eq!(prop.watched_literals_[0].next, -1);
    }

    #[test]
    fn h16_watch_link_off_second_literal_also_unreachable() {
        // Off: literal 1 is linked into the list that literal 0 left empty, so the same quirk
        // applies again and the column head stays -1 (1.15.1 HighsDomain.cpp lines 100-103).
        let prop = link_two(false);
        assert_eq!(prop.col_lower_watched_[0], -1);
        assert_eq!(prop.watched_literals_[1].prev, -1);
        assert_eq!(prop.watched_literals_[1].next, -1);
        assert_eq!(prop.watched_literals_[0].prev, -1);
        assert_eq!(prop.watched_literals_[0].next, -1);
    }

    /// Conflict 0 = {x0 >= 5 (inactive: lb 0), x1 >= 0 (active: lb 0)}; x0's bound change is the
    /// later one on the stack (colLowerPos_ 3 vs 1).
    fn conflict_added_case1(h16: bool) -> HighsDomain_ConflictPoolPropagation {
        let domain = HighsDomain {
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![10.0, 10.0],
            col_lower_pos_: vec![3, 1],
            col_upper_pos_: vec![-1, -1],
            ..Default::default()
        };
        let entries = vec![
            HighsDomainChange {
                boundval: 5.0,
                column: 0,
                boundtype: HighsBoundType::kLower,
            },
            HighsDomainChange {
                boundval: 0.0,
                column: 1,
                boundtype: HighsBoundType::kLower,
            },
        ];
        let mut prop = HighsDomain_ConflictPoolPropagation {
            col_lower_watched_: vec![-1; 2],
            col_upper_watched_: vec![-1; 2],
            ..Default::default()
        };
        prop.conflict_added(&domain, 0, &entries, h16);
        prop
    }

    #[test]
    fn h16_watch_link_on_conflict_added_watches_active_literal() {
        // 1.16.0 HighsDomain.cpp conflictAdded: first loop watches the inactive x0 >= 5 at slot 0
        // (numWatched = 1); case 1 skips inactive literals (lines 242-244), so slot 1 watches the
        // active x1 >= 0. linkWatchedLiteral sets both column heads (lines 105-106).
        let prop = conflict_added_case1(true);
        assert_eq!(prop.watched_literals_[0].domchg.column, 0);
        assert_eq!(prop.watched_literals_[1].domchg.column, 1);
        assert_eq!(prop.col_lower_watched_[0], 0);
        assert_eq!(prop.col_lower_watched_[1], 1);
        assert_eq!(prop.conflict_flag_[0] & 3, 1);
    }

    #[test]
    fn h16_watch_link_off_conflict_added_rewatches_inactive_literal() {
        // 1.15.1 conflictAdded case 1 scans every literal and takes the largest position, which is
        // the already-watched x0 >= 5 (pos 3 > 1); both heads stay -1 (quirk placement).
        let prop = conflict_added_case1(false);
        assert_eq!(prop.watched_literals_[0].domchg.column, 0);
        assert_eq!(prop.watched_literals_[1].domchg.column, 0);
        assert_eq!(prop.col_lower_watched_[0], -1);
        assert_eq!(prop.col_lower_watched_[1], -1);
        assert_eq!(prop.conflict_flag_[0] & 3, 1);
    }
}
