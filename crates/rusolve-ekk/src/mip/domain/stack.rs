//! HiGHS source: mip/HighsDomain.cpp (lines 2119-2753)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port: `set_domain_change_stack*`, `backtrack*`, `propagate` (.cpp
//! 2119-2562) and `get_col_lower_pos` .. `feastol` (.cpp 2563-2753).
//!
//! Every `&mut self` method here that changes bounds or reads propagator state starts with
//! `self.replay_pool_events(ms)`: `set_domain_change_stack*`, `backtrack_to_global`,
//! `backtrack`, `propagate`, `conflict_analysis*`, `conflict_analyze_reconvergence`.

use crate::prelude::*;
use std::cell::Cell;

thread_local! {
    /// PORT NOTE (perf, not a C++ member): reusable storage for the
    /// `std::unique_ptr<HighsDomainChange[]> changedbounds(new HighsDomainChange[n])` that
    /// `HighsDomain::propagate()` allocates per call. The C++ array is default-initialised,
    /// i.e. left uninitialised for this POD, and `propagate` only reads slots it wrote in
    /// the same pass, so stale contents are unobservable; one grow-only buffer shared by
    /// every domain (the heuristics' fresh local copies included) skips the per-call
    /// `vec![default; n]` zero-fill. Per thread, because
    /// parallel MIP workers propagate their own domains concurrently (a buffer in
    /// `mipdata_` would be moved out by two workers at once).
    static CHANGEDBOUNDS_BUF: Cell<Vec<HighsDomainChange>> = const { Cell::new(Vec::new()) };
}

/// C++ `havePropagationRows` lambda (HighsDomain.cpp:2368-2381).
fn have_propagation_rows(domain: &HighsDomain, ms: &HighsMipSolver) -> bool {
    if !domain.propagateinds_.is_empty() {
        return true;
    }

    if domain.obj_prop_.is_active() && domain.obj_prop_should_be_propagated(ms) {
        return true;
    }

    for cutpoolprop in domain.cutpoolpropagation.iter() {
        if !cutpoolprop.propagatecutinds_.is_empty() {
            return true;
        }
    }

    for conflictprop in domain.conflict_pool_propagation.iter() {
        if !conflictprop.propagate_conflict_inds_.is_empty() {
            return true;
        }
    }

    // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp line 2644 at 8a5f9c1d8
    if !domain.infeasible_ && domain.dual_fix_probing_active && domain.dual_fix_probing_.is_active()
    {
        return true;
    }

    false
}

/// HiGHS 1.16 (#3575): HighsDomain.cpp propagateIndex lambda, upper and lower tests (lines 2683, 2695)
/// C++ `propagateIndex` lambda over the model rows (HighsDomain.cpp:2412-2445): `k` indexes
/// `propagateinds` and `prop_row_num_changed_bounds_`; the row's changed bounds go to
/// `changedbounds[2 * ARstart_[i] ..]`.
fn propagate_model_row_index(
    domain: &mut HighsDomain,
    ms: &HighsMipSolver,
    propagateinds: &[HighsInt],
    changedbounds: &mut [HighsDomainChange],
    k: usize,
) {
    let i = propagateinds[k];
    let start = ms.mipdata_.a_rstart_[i as usize];
    let end = ms.mipdata_.a_rstart_[i as usize + 1];
    let rlen = end - start;
    let rindex = &ms.mipdata_.a_rindex_[start as usize..end as usize];
    let rvalue = &ms.mipdata_.a_rvalue_[start as usize..end as usize];
    let feastol = ms.mipdata_.feastol;
    let mut recompute_cap_threshold = false;

    // HiGHS 1.16 #3575 (h16_row_redundant): HighsDomain.cpp lines 2683 and 2695 (1.16.0)
    // The lower test is evaluated after the upper body, as in C++ (renormalize() runs first).
    let h16_row_redundant = ms.options_mip_.h16_row_redundant;
    let upper_active = if h16_row_redundant {
        !domain.is_row_upper_redundant(ms, i)
    } else {
        ms.row_upper_2(i) != kHighsInf
            && (domain.activitymaxinf_[i as usize] != 0
                || domain.activitymax_[i as usize] > ms.row_upper_2(i) + feastol)
    };

    if upper_active {
        // computeMinActivity(...) is the incrementally maintained activitymin_[i].
        domain.activitymin_[i as usize].renormalize();
        let rupper = ms.row_upper_2(i);
        let ninfmin = domain.activitymininf_[i as usize];
        let changed = domain.propagate_row_upper(
            ms,
            rindex,
            rvalue,
            rlen,
            rupper,
            &domain.activitymin_[i as usize],
            ninfmin,
            &mut changedbounds[(2 * start) as usize..],
        );
        domain.prop_row_num_changed_bounds_[k].0 = changed;

        recompute_cap_threshold = true;
    }

    let lower_active = if h16_row_redundant {
        !domain.is_row_lower_redundant(ms, i)
    } else {
        ms.row_lower_2(i) != -kHighsInf
            && (domain.activitymininf_[i as usize] != 0
                || domain.activitymin_[i as usize] < ms.row_lower_2(i) - feastol)
    };

    if lower_active {
        // computeMaxActivity(...) is the incrementally maintained activitymax_[i].
        domain.activitymax_[i as usize].renormalize();
        let rlower = ms.row_lower_2(i);
        let ninfmax = domain.activitymaxinf_[i as usize];
        let offset = 2 * start + domain.prop_row_num_changed_bounds_[k].0;
        let changed = domain.propagate_row_lower(
            ms,
            rindex,
            rvalue,
            rlen,
            rlower,
            &domain.activitymax_[i as usize],
            ninfmax,
            &mut changedbounds[offset as usize..],
        );
        domain.prop_row_num_changed_bounds_[k].1 = changed;

        recompute_cap_threshold = true;
    }

    if recompute_cap_threshold {
        domain.recompute_capacity_threshold(ms, i);
    }
}

/// C++ `propagateIndex` lambda over the cut rows of `cutpoolpropagation[cutpool]`
/// (HighsDomain.cpp:2486-2507). `cp` is the propagator index (C++ `cutpool`, the `Reason::cut`
/// type), `k` indexes `propagateinds` and `prop_row_num_changed_bounds_`.
fn propagate_cut_row_index(
    domain: &mut HighsDomain,
    ms: &HighsMipSolver,
    cp: usize,
    propagateinds: &[HighsInt],
    changedbounds: &mut [HighsDomainChange],
    k: usize,
) {
    let i = propagateinds[k];
    let pool = domain.cutpoolpropagation[cp].cutpool;

    // first check if cut is marked as deleted
    if (domain.cutpoolpropagation[cp].propagatecutflags_[i as usize] & 2) != 0 {
        return;
    }

    let mut rlen: HighsInt = 0;
    let mut rindex: &[HighsInt] = &[];
    let mut rvalue: &[f64] = &[];
    ms.mipdata_.cutpools[pool].get_cut(i, &mut rlen, &mut rindex, &mut rvalue);
    domain.cutpoolpropagation[cp].activitycuts_[i as usize].renormalize();
    let activity = domain.cutpoolpropagation[cp].activitycuts_[i as usize];
    let ninfmin = domain.cutpoolpropagation[cp].activitycutsinf_[i as usize];
    let rhs = ms.mipdata_.cutpools[pool].get_rhs()[i as usize];
    let row_start = ms.mipdata_.cutpools[pool].get_matrix().get_row_start(i);

    let changed = domain.propagate_row_upper(
        ms,
        rindex,
        rvalue,
        rlen,
        rhs,
        &activity,
        ninfmin,
        &mut changedbounds[row_start as usize..],
    );
    domain.prop_row_num_changed_bounds_[k].0 = changed;

    // cutpoolprop.recomputeCapacityThreshold(i): the propagator is taken out of the domain
    // while it borrows the domain (see the propagate PORT NOTE).
    let mut prop = std::mem::take(&mut domain.cutpoolpropagation[cp]);
    prop.recompute_capacity_threshold(ms, domain, i, rindex, rvalue);
    domain.cutpoolpropagation[cp] = prop;
}

impl HighsDomain {
    /// HiGHS: HighsDomain.cpp::HighsDomain::setDomainChangeStack (lines 2120-2151)
    /// C++: `void setDomainChangeStack(const std::vector<HighsDomainChange> & domchgstack)`
    pub fn set_domain_change_stack(
        &mut self,
        ms: &mut HighsMipSolver,
        domchgstack: &[HighsDomainChange],
    ) {
        self.replay_pool_events(ms);
        self.infeasible_ = false;
        // HiGHS: debugSolution.resetDomain(*this) dropped (debug-only, not ported).
        if !self.domchgstack_.is_empty() {
            for domchg in self.domchgstack_.iter() {
                if domchg.boundtype == HighsBoundType::kLower {
                    self.col_lower_pos_[domchg.column as usize] = -1;
                } else {
                    self.col_upper_pos_[domchg.column as usize] = -1;
                }
            }
        }

        self.prevboundval_.clear();
        self.domchgstack_.clear();
        self.domchgreason_.clear();
        self.branch_pos_.clear();
        for domchg in domchgstack.iter() {
            if domchg.boundtype == HighsBoundType::kLower
                && domchg.boundval <= self.col_lower_[domchg.column as usize]
            {
                continue;
            }
            if domchg.boundtype == HighsBoundType::kUpper
                && domchg.boundval >= self.col_upper_[domchg.column as usize]
            {
                continue;
            }

            self.change_bound(ms, domchg.clone(), HighsDomain_Reason::unspecified());

            if self.infeasible_ {
                break;
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::setDomainChangeStack (lines 2153-2232)
    /// C++: `void setDomainChangeStack(const std::vector<HighsDomainChange> & domchgstack, const std::vector<HighsInt> & branchingPositions)`
    pub fn set_domain_change_stack_2(
        &mut self,
        ms: &mut HighsMipSolver,
        domchgstack: &[HighsDomainChange],
        branching_positions: &[HighsInt],
    ) {
        self.replay_pool_events(ms);
        self.infeasible_ = false;
        // HiGHS: debugSolution.resetDomain(*this) dropped (debug-only, not ported).
        // HiGHS #3357 (mip_keep_branch_reason): HighsDomain.cpp lines 2418-2468 at 696d3fbb1
        let x1_keep_branch_reason = ms.options_mip_.mip_keep_branch_reason;
        let no_symmetry = !ms.mipdata_.detect_symmetries;
        if !self.domchgstack_.is_empty() {
            for domchg in self.domchgstack_.iter() {
                if domchg.boundtype == HighsBoundType::kLower {
                    self.col_lower_pos_[domchg.column as usize] = -1;
                } else {
                    self.col_upper_pos_[domchg.column as usize] = -1;
                }
            }
        }

        self.prevboundval_.clear();
        self.domchgstack_.clear();
        self.domchgreason_.clear();
        self.branch_pos_.clear();
        let stacksize = domchgstack.len() as HighsInt;
        // `k` persists across the outer loop, as in C++. The increment happens before the body
        // so that `continue` advances `k` exactly as the C++ `for (; k < branchPos; ++k)` does.
        let mut k: HighsInt = 0;
        for &branch_pos in branching_positions.iter() {
            while k < branch_pos {
                let kk = k as usize;
                k += 1;
                let domchg = &domchgstack[kk];
                if domchg.boundtype == HighsBoundType::kLower
                    && domchg.boundval <= self.col_lower_[domchg.column as usize]
                {
                    continue;
                }
                if domchg.boundtype == HighsBoundType::kUpper
                    && domchg.boundval >= self.col_upper_[domchg.column as usize]
                {
                    continue;
                }

                self.change_bound(ms, domchg.clone(), HighsDomain_Reason::unspecified());
                if !self.infeasible_ {
                    self.propagate(ms);
                }
                if self.infeasible_ {
                    return;
                }
            }

            if k == stacksize {
                return;
            }

            // For redundant branching bound changes we need to be more careful due to
            // symmetry handling. If these bound changes are redundant simply because
            // the corresponding subtree was enumerated and hence the global bound
            // updated, then we still need to keep their status as branching variables
            // for computing correct stabilizers. They can, however, be safely dropped
            // if they are either strictly redundant in the global domain, or if there
            // is already a local bound change that makes the branching change
            // redundant.
            // HiGHS #3357 (mip_keep_branch_reason): HighsDomain.cpp lines 2448-2468 at 696d3fbb1
            if !x1_keep_branch_reason || no_symmetry {
                let domchg = &domchgstack[k as usize];
                if domchg.boundtype == HighsBoundType::kLower {
                    if domchg.boundval <= self.col_lower_[domchg.column as usize] {
                        if domchg.boundval < self.col_lower_[domchg.column as usize] {
                            continue;
                        }
                        if self.col_lower_pos_[domchg.column as usize] != -1 {
                            continue;
                        }
                    }
                } else if domchg.boundval >= self.col_upper_[domchg.column as usize] {
                    if domchg.boundval > self.col_upper_[domchg.column as usize] {
                        continue;
                    }
                    if self.col_upper_pos_[domchg.column as usize] != -1 {
                        continue;
                    }
                }
            }

            self.change_bound(
                ms,
                domchgstack[k as usize].clone(),
                HighsDomain_Reason::branching(),
            );
            if !self.infeasible_ {
                self.propagate(ms);
            }
            if self.infeasible_ {
                return;
            }
        }

        for kk in k as usize..stacksize as usize {
            let domchg = &domchgstack[kk];
            if domchg.boundtype == HighsBoundType::kLower
                && domchg.boundval <= self.col_lower_[domchg.column as usize]
            {
                continue;
            }
            if domchg.boundtype == HighsBoundType::kUpper
                && domchg.boundval >= self.col_upper_[domchg.column as usize]
            {
                continue;
            }

            // HiGHS: HighsDomain.cpp:2467-2468 (1.16.0) / 2225-2226 (1.15.1) —
            //   debugSolution.boundChangeAdded(*this, domchgstack[k], true);
            // The `branching == true` argument makes this a suppression rather than a
            // report: HighsDebugSol.cpp:151-154 (1.16.0) / 148-151 (1.15.1) inserts the
            // change into conflictingBounds and returns, so the boundChangeAdded that
            // changeBound makes for this same change with `branching == false` (it is
            // replayed with Reason::unspecified()) takes the same silent path instead of the
            // assert. Without it the port reported replayed branching changes of a node.
            ms.mipdata_
                .debug_solution
                .bound_change_added(&*self, domchg, true, None);
            ms.mipdata_.debug_solution.begin_branching_replay();
            self.change_bound(ms, domchg.clone(), HighsDomain_Reason::unspecified());
            ms.mipdata_.debug_solution.end_branching_replay();
            if !self.infeasible_ {
                self.propagate(ms);
            }
            if self.infeasible_ {
                break;
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::backtrackToGlobal (lines 2234-2292)
    /// C++: `void backtrackToGlobal()`
    pub fn backtrack_to_global(&mut self, ms: &mut HighsMipSolver) {
        self.replay_pool_events(ms);
        let mut k: HighsInt = self.domchgstack_.len() as HighsInt - 1;
        let old_infeasible = self.infeasible_;
        let old_reason = self.infeasible_reason;

        if self.infeasible_ && self.infeasible_pos == self.domchgstack_.len() as HighsInt {
            debug_assert!(old_infeasible);
            debug_assert!(k == self.domchgstack_.len() as HighsInt - 1);
            self.infeasible_ = false;
            self.infeasible_reason = HighsDomain_Reason::unspecified();
        }

        while k >= 0 {
            let prevbound = self.prevboundval_[k as usize].0;
            let prevpos = self.prevboundval_[k as usize].1;
            debug_assert!(prevpos < k);

            // HiGHS: debugSolution.boundChangeRemoved(*this, domchgstack_[k]) dropped
            // (debug-only, not ported).
            let chg = self.domchgstack_[k as usize].clone();
            if chg.boundtype == HighsBoundType::kLower {
                debug_assert!(self.col_lower_pos_[chg.column as usize] == k);
                self.col_lower_pos_[chg.column as usize] = prevpos;
            } else {
                debug_assert!(self.col_upper_pos_[chg.column as usize] == k);
                self.col_upper_pos_[chg.column as usize] = prevpos;
            }

            if prevbound != chg.boundval {
                // change back to global bound
                self.do_change_bound(
                    ms,
                    &HighsDomainChange {
                        boundval: prevbound,
                        column: chg.column,
                        boundtype: chg.boundtype,
                    },
                );
            }

            if self.infeasible_ && self.infeasible_pos == k {
                debug_assert!(old_infeasible);
                debug_assert!(k == self.domchgstack_.len() as HighsInt - 1);
                self.infeasible_ = false;
                self.infeasible_reason = HighsDomain_Reason::unspecified();
            }

            k -= 1;
        }

        if old_infeasible {
            self.mark_propagate_cut(ms, old_reason);
            self.infeasible_reason = HighsDomain_Reason::unspecified();
            self.infeasible_ = false;
        }

        let numreason = self.domchgreason_.len() as HighsInt;
        for i in (k + 1)..numreason {
            self.mark_propagate_cut(ms, self.domchgreason_[i as usize]);
        }

        self.domchgstack_.clear();
        self.prevboundval_.clear();
        self.domchgreason_.clear();
        self.branch_pos_.clear();
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::backtrack (lines 2294-2364)
    /// C++: `HighsDomainChange backtrack()`
    pub fn backtrack(&mut self, ms: &mut HighsMipSolver) -> HighsDomainChange {
        self.replay_pool_events(ms);
        let mut k: HighsInt = self.domchgstack_.len() as HighsInt - 1;
        let old_infeasible = self.infeasible_;
        let old_reason = self.infeasible_reason;

        if self.infeasible_ && self.infeasible_pos == self.domchgstack_.len() as HighsInt {
            debug_assert!(old_infeasible);
            debug_assert!(k == self.domchgstack_.len() as HighsInt - 1);
            self.infeasible_ = false;
            self.infeasible_reason = HighsDomain_Reason::unspecified();
        }
        while k >= 0 {
            let prevbound = self.prevboundval_[k as usize].0;
            let prevpos = self.prevboundval_[k as usize].1;
            debug_assert!(prevpos < k);

            // HiGHS: debugSolution.boundChangeRemoved(*this, domchgstack_[k]) dropped
            // (debug-only, not ported).
            let chg = self.domchgstack_[k as usize].clone();
            if chg.boundtype == HighsBoundType::kLower {
                debug_assert!(self.col_lower_pos_[chg.column as usize] == k);
                self.col_lower_pos_[chg.column as usize] = prevpos;
            } else {
                debug_assert!(self.col_upper_pos_[chg.column as usize] == k);
                self.col_upper_pos_[chg.column as usize] = prevpos;
            }

            // change back to global bound
            self.do_change_bound(
                ms,
                &HighsDomainChange {
                    boundval: prevbound,
                    column: chg.column,
                    boundtype: chg.boundtype,
                },
            );

            if self.infeasible_ && self.infeasible_pos == k {
                debug_assert!(old_infeasible);
                debug_assert!(k == self.domchgstack_.len() as HighsInt - 1);
                self.infeasible_ = false;
                self.infeasible_reason = HighsDomain_Reason::unspecified();
            }

            if self.domchgreason_[k as usize].r#type == HighsDomain_Reason::kBranching {
                // HiGHS quirk: `branchPos_.pop_back()` on an empty vector is undefined behaviour
                // in C++; `pop` is a no-op there.
                self.branch_pos_.pop();
                break;
            }

            k -= 1;
        }

        if old_infeasible {
            self.mark_propagate_cut(ms, old_reason);
            self.infeasible_reason = HighsDomain_Reason::unspecified();
            self.infeasible_ = false;
        }

        let numreason = self.domchgreason_.len() as HighsInt;
        for i in (k + 1)..numreason {
            self.mark_propagate_cut(ms, self.domchgreason_[i as usize]);
        }

        if k < 0 {
            self.domchgstack_.clear();
            self.prevboundval_.clear();
            self.domchgreason_.clear();
            self.branch_pos_.clear();
            return HighsDomainChange {
                boundval: 0.0,
                column: -1,
                boundtype: HighsBoundType::kLower,
            };
        }

        let backtrackboundchg = self.domchgstack_[k as usize].clone();
        self.domchgstack_.truncate(k as usize);
        self.domchgreason_.truncate(k as usize);
        self.prevboundval_.truncate(k as usize);

        backtrackboundchg
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::propagate (lines 2366-2561)
    /// C++: `bool propagate()`
    ///
    /// PORT NOTE: `auto& cutpoolprop = cutpoolpropagation[cutpool]` is held across
    /// `changeBound` in C++ (deque references are stable); here re-index
    /// `self.cutpoolpropagation[cutpool]` / `self.conflict_pool_propagation[k]` at every use
    /// and never hold a borrow across `change_bound`. `propagateinds.swap(...)` is
    /// `std::mem::swap(&mut propagateinds, &mut self.cutpoolpropagation[cutpool].propagatecutinds_)`.
    /// `cutpoolprop.cutpool->getCut(i, ..)` copies the row out of
    /// `ms.mipdata_.cutpools[pool]` before `propagate_row_upper` (which takes `ms` by `&`) and
    /// `cutpoolprop.recomputeCapacityThreshold(i)` is
    /// `recompute_capacity_threshold(ms, self, i, inds, vals)` on the taken-out element.
    /// `conflictprop.propagateConflict(c)` is `self.propagate_conflict(ms, k, c)`;
    /// `objProp_.propagate()` is `self.obj_prop_propagate(ms)`; `resetAge` goes to
    /// `ms.mipdata_.cutpools[pool].reset_age(i, ms.mipdata_.parallel_lock_active())`
    /// (read the flag first).
    #[allow(clippy::needless_range_loop)]
    pub fn propagate(&mut self, ms: &mut HighsMipSolver) -> bool {
        self.replay_pool_events(ms);
        let mut propagateinds: Vec<HighsInt> = Vec::new();

        if !have_propagation_rows(self, ms) {
            return false;
        }

        let mut changedboundsize = 2 * ms.mipdata_.a_rvalue_.len();
        for cutpoolprop in self.cutpoolpropagation.iter() {
            changedboundsize = cpp_max(
                changedboundsize,
                ms.mipdata_.cutpools[cutpoolprop.cutpool]
                    .get_matrix()
                    .nonzero_capacity(),
            );
        }

        // `new HighsDomainChange[changedboundsize]`: taken out of the thread's
        // [`CHANGEDBOUNDS_BUF`] (so a nested propagate gets its own array) and grown, never
        // cleared.
        let mut changedbounds = CHANGEDBOUNDS_BUF.take();
        if changedbounds.len() < changedboundsize {
            changedbounds.resize(changedboundsize, HighsDomainChange::default());
        }

        while have_propagation_rows(self, ms) {
            if self.obj_prop_.is_active() {
                self.obj_prop_propagate(ms);
            }

            let num_conflict_pools = self.conflict_pool_propagation.len();
            for conflict_pool in 0..num_conflict_pools {
                while !self.conflict_pool_propagation[conflict_pool]
                    .propagate_conflict_inds_
                    .is_empty()
                {
                    std::mem::swap(
                        &mut propagateinds,
                        &mut self.conflict_pool_propagation[conflict_pool].propagate_conflict_inds_,
                    );

                    for &conflict in propagateinds.iter() {
                        self.propagate_conflict(ms, conflict_pool, conflict);
                    }

                    propagateinds.clear();
                }
            }

            if !self.propagateinds_.is_empty() {
                std::mem::swap(&mut propagateinds, &mut self.propagateinds_);

                let numproprows = propagateinds.len();
                for i in 0..numproprows {
                    let row = propagateinds[i];
                    self.propagateflags_[row as usize] = 0;
                }

                if !self.infeasible_ {
                    self.prop_row_num_changed_bounds_.clear();
                    self.prop_row_num_changed_bounds_
                        .resize(numproprows, (0, 0));

                    for k in 0..numproprows {
                        propagate_model_row_index(self, ms, &propagateinds, &mut changedbounds, k);
                    }

                    for k in 0..numproprows {
                        let i = propagateinds[k];

                        if self.prop_row_num_changed_bounds_[k].0 != 0 {
                            let start = 2 * ms.mipdata_.a_rstart_[i as usize];
                            let end = start + self.prop_row_num_changed_bounds_[k].0;
                            let mut j = start;
                            while j != end && !self.infeasible_ {
                                let chg = changedbounds[j as usize].clone();
                                self.change_bound(ms, chg, HighsDomain_Reason::model_row_upper(i));
                                j += 1;
                            }

                            if self.infeasible_ {
                                break;
                            }
                        }
                        if self.prop_row_num_changed_bounds_[k].1 != 0 {
                            let start = 2 * ms.mipdata_.a_rstart_[i as usize]
                                + self.prop_row_num_changed_bounds_[k].0;
                            let end = start + self.prop_row_num_changed_bounds_[k].1;
                            let mut j = start;
                            while j != end && !self.infeasible_ {
                                let chg = changedbounds[j as usize].clone();
                                self.change_bound(ms, chg, HighsDomain_Reason::model_row_lower(i));
                                j += 1;
                            }

                            if self.infeasible_ {
                                break;
                            }
                        }
                    }
                }

                propagateinds.clear();
            }

            let numpools = self.cutpoolpropagation.len();
            for cutpool in 0..numpools {
                if !self.cutpoolpropagation[cutpool]
                    .propagatecutinds_
                    .is_empty()
                {
                    std::mem::swap(
                        &mut propagateinds,
                        &mut self.cutpoolpropagation[cutpool].propagatecutinds_,
                    );

                    let numproprows = propagateinds.len();
                    for i in 0..numproprows {
                        let cut = propagateinds[i];
                        self.cutpoolpropagation[cutpool].propagatecutflags_[cut as usize] &= 2;
                    }

                    if !self.infeasible_ {
                        self.prop_row_num_changed_bounds_.clear();
                        self.prop_row_num_changed_bounds_
                            .resize(numproprows, (0, 0));

                        for k in 0..numproprows {
                            propagate_cut_row_index(
                                self,
                                ms,
                                cutpool,
                                &propagateinds,
                                &mut changedbounds,
                                k,
                            );
                        }

                        for k in 0..numproprows {
                            let i = propagateinds[k];
                            if self.prop_row_num_changed_bounds_[k].0 != 0 {
                                let pool = self.cutpoolpropagation[cutpool].cutpool;
                                // `parallelLockActive()` is read before the pool is borrowed
                                // mutably for `resetAge`.
                                let thread_safe = ms.parallel_lock_active();
                                ms.mipdata_.cutpools[pool].reset_age(i, thread_safe);
                                let start =
                                    ms.mipdata_.cutpools[pool].get_matrix().get_row_start(i);
                                let end = start + self.prop_row_num_changed_bounds_[k].0;
                                let mut j = start;
                                while j != end && !self.infeasible_ {
                                    let chg = changedbounds[j as usize].clone();
                                    self.change_bound(
                                        ms,
                                        chg,
                                        HighsDomain_Reason::cut(cutpool as HighsInt, i),
                                    );
                                    j += 1;
                                }
                            }

                            if self.infeasible_ {
                                break;
                            }
                        }
                    }

                    propagateinds.clear();
                }
            }

            // HiGHS #3244 (dual_fix_probing): HighsDomain.cpp lines 2822-2828 at 8a5f9c1d8
            // The C++ block is unguarded; `enabled_` is only set by `beginProbing`, which runs
            // only while `dualFixProbingActive_` is true, so the rule-F gate is equivalent.
            if self.dual_fix_probing_active {
                if !self.infeasible_ && self.dual_fix_probing_.is_active() {
                    self.dfp_propagate(ms);
                }
                if !self.infeasible_
                    && self.dual_fix_probing_.is_enabled()
                    && !have_propagation_rows(self, ms)
                {
                    self.dfp_propagate_zero_costs(ms);
                }
            }
        }

        CHANGEDBOUNDS_BUF.set(changedbounds);
        true
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::getColLowerPos (lines 2563-2572)
    /// C++: `double getColLowerPos(HighsInt col, HighsInt stackpos, HighsInt & pos) const`
    pub fn get_col_lower_pos(&self, col: HighsInt, stackpos: HighsInt, pos: &mut HighsInt) -> f64 {
        let mut lb = self.col_lower_[col as usize];
        *pos = self.col_lower_pos_[col as usize];
        while *pos > stackpos || (*pos != -1 && self.prevboundval_[*pos as usize].0 == lb) {
            lb = self.prevboundval_[*pos as usize].0;
            *pos = self.prevboundval_[*pos as usize].1;
        }
        lb
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::getColUpperPos (lines 2574-2583)
    /// C++: `double getColUpperPos(HighsInt col, HighsInt stackpos, HighsInt & pos) const`
    pub fn get_col_upper_pos(&self, col: HighsInt, stackpos: HighsInt, pos: &mut HighsInt) -> f64 {
        let mut ub = self.col_upper_[col as usize];
        *pos = self.col_upper_pos_[col as usize];
        while *pos > stackpos || (*pos != -1 && self.prevboundval_[*pos as usize].0 == ub) {
            ub = self.prevboundval_[*pos as usize].0;
            *pos = self.prevboundval_[*pos as usize].1;
        }
        ub
    }

    /// HiGHS 1.16 (#3196): HighsDomain.h::HighsDomain::walkBoundChain (lines 670-681)
    /// C++: `bool walkBoundChain(HighsInt& pos, double target, HighsInt direction) const`
    /// Only reached on the `h16_bound_chain` path. `direction` is `1` (walk while the
    /// previous bound value >= target) or `-1` (walk while it is <= target). Returns false
    /// when the chain is exhausted (`pos` becomes -1); 1.15.1 had no such check.
    pub fn walk_bound_chain(&self, pos: &mut HighsInt, target: f64, direction: HighsInt) -> bool {
        while *pos != -1
            && f64::from(direction) * self.prevboundval_[*pos as usize].0
                >= f64::from(direction) * target
        {
            *pos = self.prevboundval_[*pos as usize].1;
        }
        *pos != -1
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::conflictAnalysis (lines 2585-2598)
    /// C++: `void conflictAnalysis(HighsConflictPool & conflictPool, HighsDomain & globaldom, HighsPseudocost & pseudocost)`
    ///
    /// PORT NOTE: signature changed (domain by index): the pool is
    /// `ms.mipdata_.conflictpools[conflict_pool]` and the global domain
    /// `ms.mipdata_.domains[globaldom]` (all C++ callers pass `mipdata_` members). `&globaldom ==
    /// this` is `ms.mipdata_.domains[globaldom].holds(self)` and is tested first, before the
    /// slot is touched. `globaldom.propagate()` is
    /// `detach!(ms.mipdata_.domains[globaldom], |g| g.propagate(ms))`. The `ConflictSet` is
    /// built on `self` and driven through its `ms: &mut` entry points (conflict_set.rs).
    /// `pseudocost` is detached by the caller when it lives in `ms`.
    pub fn conflict_analysis(
        &mut self,
        ms: &mut HighsMipSolver,
        conflict_pool: usize,
        globaldom: usize,
        pseudocost: &mut HighsPseudocost,
    ) {
        // rusolve ablation lever mip_conflict_off_mask
        if crate::work::conflict_src_masked(ms.options_mip_.mip_conflict_off_mask) {
            return;
        }
        // C++ pools notify this domain's propagators eagerly; the ConflictSet reads
        // `conflictFlag_` (bit 8 = deleted), so pending events must be applied first.
        self.replay_pool_events(ms);
        if ms.mipdata_.domains[globaldom].holds(self as *const HighsDomain) {
            return;
        }
        if ms.mipdata_.domains[globaldom].infeasible() || !self.infeasible_ {
            return;
        }

        // Not sure how this should be modified for the workers.
        detach!(ms.mipdata_.domains[globaldom], |g| g.propagate(ms));
        if ms.mipdata_.domains[globaldom].infeasible() {
            return;
        }

        crate::work::conflict_event(crate::work::ConflictKind::Calls, 1);
        let mut conflict_set = HighsDomain_ConflictSet::new();
        conflict_set.conflict_analysis(ms, self, globaldom, conflict_pool, pseudocost);
        // mip_conflict_resolution: SCIP scip_conflict.c lines 735-764 at 32dd425
        let s2 = ms.options_mip_.mip_conflict_resolution
            && !(ms.options_mip_.mip_conflict_resolution_tree_only && ms.mipdata_.num_nodes == 0);
        if s2 {
            crate::mip::conflict_resolution::output::s2_conflict_analysis(
                ms,
                self,
                globaldom,
                conflict_pool,
                None,
            );
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::conflictAnalysis (lines 2600-2616)
    /// C++: `void conflictAnalysis(const HighsInt * proofinds, const double * proofvals, HighsInt prooflen, double proofrhs, HighsConflictPool & conflictPool, HighsDomain & globaldom, HighsPseudocost & pseudocost)`
    ///
    /// PORT NOTE: signature changed (domain by index): the pool is
    /// `ms.mipdata_.conflictpools[conflict_pool]` and the global domain
    /// `ms.mipdata_.domains[globaldom]` (all C++ callers pass `mipdata_` members). `&globaldom ==
    /// this` is `ms.mipdata_.domains[globaldom].holds(self)` and is tested first, before the
    /// slot is touched. `globaldom.propagate()` is
    /// `detach!(ms.mipdata_.domains[globaldom], |g| g.propagate(ms))`. The `ConflictSet` is
    /// built on `self` and driven through its `ms: &mut` entry points (conflict_set.rs).
    /// `pseudocost` is detached by the caller when it lives in `ms`.
    pub fn conflict_analysis_2(
        &mut self,
        ms: &mut HighsMipSolver,
        proofinds: &[HighsInt],
        proofvals: &[f64],
        prooflen: HighsInt,
        proofrhs: f64,
        conflict_pool: usize,
        globaldom: usize,
        pseudocost: &mut HighsPseudocost,
    ) {
        // rusolve ablation lever mip_conflict_off_mask
        if crate::work::conflict_src_masked(ms.options_mip_.mip_conflict_off_mask) {
            return;
        }
        // C++ pools notify this domain's propagators eagerly; the ConflictSet reads
        // `conflictFlag_` (bit 8 = deleted), so pending events must be applied first.
        self.replay_pool_events(ms);
        if ms.mipdata_.domains[globaldom].holds(self as *const HighsDomain) {
            return;
        }
        if ms.mipdata_.domains[globaldom].infeasible() {
            return;
        }

        detach!(ms.mipdata_.domains[globaldom], |g| g.propagate(ms));
        if ms.mipdata_.domains[globaldom].infeasible() {
            return;
        }

        crate::work::conflict_event(crate::work::ConflictKind::Calls, 1);
        let mut conflict_set = HighsDomain_ConflictSet::new();
        conflict_set.conflict_analysis_2(
            ms,
            self,
            globaldom,
            proofinds,
            proofvals,
            prooflen,
            proofrhs,
            conflict_pool,
            pseudocost,
        );
        // mip_conflict_resolution: SCIP scip_conflict.c lines 735-764 at 32dd425
        let s2 = ms.options_mip_.mip_conflict_resolution
            && !(ms.options_mip_.mip_conflict_resolution_tree_only && ms.mipdata_.num_nodes == 0);
        if s2 {
            crate::mip::conflict_resolution::output::s2_conflict_analysis(
                ms,
                self,
                globaldom,
                conflict_pool,
                Some((proofinds, proofvals, prooflen, proofrhs)),
            );
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::conflictAnalyzeReconvergence (lines 2618-2666)
    /// C++: `void conflictAnalyzeReconvergence(const HighsDomainChange & domchg, const HighsInt * proofinds, const double * proofvals, HighsInt prooflen, double proofrhs, HighsConflictPool & conflictPool, HighsDomain & globaldom, HighsPseudocost & pseudocost)`
    ///
    /// PORT NOTE: signature changed (domain by index): the pool is
    /// `ms.mipdata_.conflictpools[conflict_pool]` and the global domain
    /// `ms.mipdata_.domains[globaldom]` (all C++ callers pass `mipdata_` members). `&globaldom ==
    /// this` is `ms.mipdata_.domains[globaldom].holds(self)` and is tested first, before the
    /// slot is touched. `globaldom.propagate()` is
    /// `detach!(ms.mipdata_.domains[globaldom], |g| g.propagate(ms))`. The `ConflictSet` is
    /// built on `self` and driven through its `ms: &mut` entry points (conflict_set.rs).
    /// `pseudocost` is detached by the caller when it lives in `ms`.
    pub fn conflict_analyze_reconvergence(
        &mut self,
        ms: &mut HighsMipSolver,
        domchg: &HighsDomainChange,
        proofinds: &[HighsInt],
        proofvals: &[f64],
        prooflen: HighsInt,
        proofrhs: f64,
        conflict_pool: usize,
        globaldom: usize,
        pseudocost: &mut HighsPseudocost,
    ) {
        let _src = crate::work::conflict_src(crate::work::ConflictSrc::Degen);
        // rusolve ablation lever mip_conflict_off_mask
        if crate::work::conflict_src_masked(ms.options_mip_.mip_conflict_off_mask) {
            return;
        }
        // C++ pools notify this domain's propagators eagerly; the ConflictSet reads
        // `conflictFlag_` (bit 8 = deleted), so pending events must be applied first.
        self.replay_pool_events(ms);
        if ms.mipdata_.domains[globaldom].holds(self as *const HighsDomain) {
            return;
        }
        if ms.mipdata_.domains[globaldom].infeasible() {
            return;
        }

        detach!(ms.mipdata_.domains[globaldom], |g| g.propagate(ms));
        if ms.mipdata_.domains[globaldom].infeasible() {
            return;
        }

        crate::work::conflict_event(crate::work::ConflictKind::Calls, 1);
        let mut conflict_set = HighsDomain_ConflictSet::new();

        let mut ninfmin: HighsInt = 0;
        let mut activitymin = HighsCDouble::from(0.0);
        ms.mipdata_.domains[globaldom].compute_min_activity(
            0,
            prooflen,
            proofinds,
            proofvals,
            &mut ninfmin,
            &mut activitymin,
        );
        if ninfmin != 0 {
            return;
        }

        // C++ passes `conflictSet.reconvergenceFrontier` by reference to methods that also
        // take `conflictSet` mutably; the frontier is taken out for the calls and used as a
        // local (it is only ever read by `explainBoundChangeLeq` and written by `insert`).
        let mut frontier = std::mem::take(&mut conflict_set.reconvergence_frontier);
        let chg = HighsDomain_ConflictSet_LocalDomChg {
            pos: self.domchgstack_.len() as HighsInt,
            domchg: domchg.clone(),
        };
        if !conflict_set.explain_bound_change_leq(
            ms,
            self,
            &ms.mipdata_.domains[globaldom],
            &frontier,
            &chg,
            proofinds,
            proofvals,
            prooflen,
            proofrhs,
            f64::from(activitymin),
        ) {
            return;
        }

        if 10 * conflict_set.resolved_domain_changes.len()
            > 1000 + 3 * ms.mipdata_.integral_cols.len()
        {
            return;
        }

        crate::work::conflict_event(crate::work::ConflictKind::Explained, 1);
        frontier.extend(conflict_set.resolved_domain_changes.iter().cloned());

        let mut depth = self.branch_pos_.len() as HighsInt;
        while depth > 0 {
            let branch_pos = self.branch_pos_[(depth - 1) as usize];
            if self.domchgstack_[branch_pos as usize].boundval
                != self.prevboundval_[branch_pos as usize].0
            {
                break;
            }

            depth -= 1;
        }

        conflict_set.resolve_depth(
            ms,
            self,
            &ms.mipdata_.domains[globaldom],
            &mut frontier,
            depth,
            0,
            pseudocost,
            0,
            false,
        );
        detach!(ms.mipdata_.conflictpools[conflict_pool], |p| p
            .add_reconvergence_cut(ms, self, &frontier, domchg));
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::tightenCoefficients (lines 2668-2710)
    /// C++: `void tightenCoefficients(HighsInt * inds, double * vals, HighsInt len, double & rhs) const`
    pub fn tighten_coefficients(
        &self,
        ms: &HighsMipSolver,
        inds: &mut [HighsInt],
        vals: &mut [f64],
        len: HighsInt,
        rhs: &mut f64,
    ) {
        // HiGHS #3223 (mip_cdouble_proofs): HighsDomain.cpp lines 2958, 2962 at 7db1c8180
        let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;

        let mut maxactivity = HighsCDouble::from(0.0);

        for i in 0..len as usize {
            if vals[i] > 0.0 {
                if self.col_upper_[inds[i] as usize] == kHighsInf {
                    return;
                }

                if x1b_cdouble {
                    maxactivity += HighsCDouble::from(self.col_upper_[inds[i] as usize]) * vals[i];
                } else {
                    maxactivity += self.col_upper_[inds[i] as usize] * vals[i];
                }
            } else {
                if self.col_lower_[inds[i] as usize] == -kHighsInf {
                    return;
                }

                if x1b_cdouble {
                    maxactivity += HighsCDouble::from(self.col_lower_[inds[i] as usize]) * vals[i];
                } else {
                    maxactivity += self.col_lower_[inds[i] as usize] * vals[i];
                }
            }
        }

        let maxabscoef = maxactivity - *rhs;
        if maxabscoef > ms.mipdata_.feastol {
            let mut upper = HighsCDouble::from(*rhs);
            let mut tightened: HighsInt = 0;
            for i in 0..len as usize {
                if ms.is_col_continuous(inds[i]) {
                    continue;
                }
                if vals[i] > maxabscoef {
                    let delta = vals[i] - maxabscoef;
                    upper -= delta * self.col_upper_[inds[i] as usize];
                    vals[i] = f64::from(maxabscoef);
                    tightened += 1;
                } else if vals[i] < -maxabscoef {
                    let delta = -vals[i] - maxabscoef;
                    upper += delta * self.col_lower_[inds[i] as usize];
                    vals[i] = -f64::from(maxabscoef);
                    tightened += 1;
                }
            }

            if tightened != 0 {
                *rhs = f64::from(upper);
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::getMinCutActivity (lines 2712-2727)
    /// C++: `double getMinCutActivity(const HighsCutPool & cutpool, HighsInt cut) const`
    // PORT NOTE: signature changed: `cutpoolprop.cutpool == &cutpool` compares pool indices
    // (`cutpool` is the index into `ms.mipdata_.cutpools`). The C++ is `const` and sees
    // eagerly-notified propagators; here pending pool events are replayed first, which needs
    // `&mut self`. [`HighsDomain::get_min_cut_activity_replayed`] is the C++ body.
    pub fn get_min_cut_activity(
        &mut self,
        ms: &HighsMipSolver,
        cutpool: usize,
        cut: HighsInt,
    ) -> f64 {
        self.replay_pool_events(ms);
        self.get_min_cut_activity_replayed(cutpool, cut)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::getMinCutActivity (lines 2712-2727), C++ body.
    ///
    /// PORT NOTE: for callers that hold the domain by `&` (ConflictSet's `globaldom`): valid
    /// when no event of that pool is pending for this domain, which holds after
    /// `conflict_analysis*` ran `globaldom.propagate()` (adding conflicts does not touch cut
    /// pool logs).
    pub fn get_min_cut_activity_replayed(&self, cutpool: usize, cut: HighsInt) -> f64 {
        for cutpoolprop in self.cutpoolpropagation.iter() {
            if cutpoolprop.cutpool == cutpool {
                return if cut < cutpoolprop.propagatecutflags_.len() as HighsInt
                    && (cutpoolprop.propagatecutflags_[cut as usize] & 2) == 0
                    && cutpoolprop.activitycutsinf_[cut as usize] == 0
                {
                    f64::from(cutpoolprop.activitycuts_[cut as usize])
                } else {
                    -kHighsInf
                };
            }
        }

        -kHighsInf
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::isFixing (lines 2729-2735)
    /// C++: `bool isFixing(const HighsDomainChange & domchg) const`
    pub fn is_fixing(&self, ms: &HighsMipSolver, domchg: &HighsDomainChange) -> bool {
        let otherbound = if domchg.boundtype == HighsBoundType::kUpper {
            self.col_lower_[domchg.column as usize]
        } else {
            self.col_upper_[domchg.column as usize]
        };
        (domchg.boundval - otherbound).abs() <= ms.mipdata_.epsilon
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::flip (lines 2737-2751)
    /// C++: `HighsDomainChange flip(const HighsDomainChange & domchg) const`
    pub fn flip(&self, ms: &HighsMipSolver, domchg: &HighsDomainChange) -> HighsDomainChange {
        if domchg.boundtype == HighsBoundType::kLower {
            let mut flipped = HighsDomainChange {
                boundval: domchg.boundval - ms.mipdata_.feastol,
                column: domchg.column,
                boundtype: HighsBoundType::kUpper,
            };
            if ms.is_col_integral(domchg.column) {
                flipped.boundval = flipped.boundval.floor();
            }
            flipped
        } else {
            let mut flipped = HighsDomainChange {
                boundval: domchg.boundval + ms.mipdata_.feastol,
                column: domchg.column,
                boundtype: HighsBoundType::kLower,
            };
            if ms.is_col_integral(domchg.column) {
                flipped.boundval = flipped.boundval.ceil();
            }
            flipped
        }
    }

    /// HiGHS latest (#3223, 7db1c8180): HighsDomain.cpp::HighsDomain::termIsZero (lines 69-74)
    /// C++: `static bool termIsZero(double val, double lb, double ub, double epsilon)`
    pub fn term_is_zero(val: f64, lb: f64, ub: f64, epsilon: f64) -> bool {
        let abs_val = val.abs();
        if abs_val <= f64::MIN_POSITIVE {
            return true;
        }
        abs_val <= epsilon && abs_val * cpp_max(lb.abs(), ub.abs()) <= epsilon
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::feastol (lines 2753-2753)
    /// C++: `double feastol() const`
    pub fn feastol(&self, ms: &HighsMipSolver) -> f64 {
        ms.mipdata_.feastol
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    // C++ HighsDomain.cpp:69-74 at 7db1c8180
    fn x1b_train_a_term_is_zero() {
        assert!(HighsDomain::term_is_zero(
            f64::MIN_POSITIVE / 2.0,
            -1.0,
            1.0,
            1e-9
        ));
        assert!(!HighsDomain::term_is_zero(1e-10, -1e3, 1e3, 1e-9));
        assert!(HighsDomain::term_is_zero(1e-12, -1.0, 1.0, 1e-9));
    }

    #[test]
    // C++ HighsDomain.cpp:2958, 2962 at 7db1c8180
    fn x1b_train_a_tighten_cdouble_on_off() {
        // Option on: HighsCDouble accumulation gives exact product 3*3333333333333333 - 1e16 = -1,
        // maxabscoef = -0.5, no tightening.
        let mut ms_on = tiny_mip();
        ms_on.options_mip_.mip_cdouble_proofs = true;
        let mut dom_on = HighsDomain::new(&ms_on);
        dom_on.col_upper_ = vec![3.0, 10.0];
        dom_on.col_lower_ = vec![0.0, 1.0];
        let mut inds_on = [0 as HighsInt, 1 as HighsInt];
        let mut vals_on = [3333333333333333.0f64, -1e16];
        let mut rhs_on = -0.5f64;
        dom_on.tighten_coefficients(&ms_on, &mut inds_on, &mut vals_on, 2, &mut rhs_on);
        assert_eq!(vals_on, [3333333333333333.0f64, -1e16]);
        assert_eq!(rhs_on, -0.5);

        // Option off: 1.15.1 plain-double path. f64 product 3*3333333333333333 rounds to 1e16,
        // the sum is 0, maxabscoef = 0.5, both coefficients tightened to maxabscoef.
        let mut ms_off = tiny_mip();
        ms_off.options_mip_.mip_cdouble_proofs = false;
        let mut dom_off = HighsDomain::new(&ms_off);
        dom_off.col_upper_ = vec![3.0, 10.0];
        dom_off.col_lower_ = vec![0.0, 1.0];
        let mut inds_off = [0 as HighsInt, 1 as HighsInt];
        let mut vals_off = [3333333333333333.0f64, -1e16];
        let mut rhs_off = -0.5f64;
        dom_off.tighten_coefficients(&ms_off, &mut inds_off, &mut vals_off, 2, &mut rhs_off);
        assert_eq!(vals_off, [0.5f64, -0.5]);
    }

    #[test]
    fn x1_dom_cd_contrib_exact() {
        // Option on: set ms.options_mip_.mip_cdouble_act_contrib = true before HighsDomain::new
        let mut ms = tiny_mip();
        ms.options_mip_.mip_cdouble_act_contrib = true;
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);

        // Change x0 lower to 3.0, then propagate row upper (x0 + x1 <= 4)
        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        assert!(dom.propagate(&mut ms));
        assert!(!dom.infeasible());
        // x1 <= 4 - 3 = 1
        assert_eq!(dom.col_upper_[1], 1.0);

        // Option off: same input, 1.15.1 result
        let mut ms2 = tiny_mip();
        ms2.options_mip_.mip_cdouble_act_contrib = false;
        let mut dom2 = HighsDomain::new(&ms2);
        dom2.compute_row_activities(&ms2);
        dom2.change_bound(&mut ms2, lower(0, 3.0), HighsDomain_Reason::branching());
        assert!(dom2.propagate(&mut ms2));
        assert!(!dom2.infeasible());
        assert_eq!(dom2.col_upper_[1], 1.0);

        // A coefficient and bound whose f64 product rounds: 0.1 * 3.0 != 0.3 exactly, so the
        // HighsCDouble contribution carries a nonzero low word that the f64 one drops.
        let cd = activity_contribution_min_cd(0.1, 3.0, 10.0);
        assert_eq!(cd.hi, activity_contribution_min(0.1, 3.0, 10.0));
        assert_ne!(cd.lo, 0.0);
        let cdmax = activity_contribution_max_cd(0.1, 0.0, 3.0);
        assert_eq!(cdmax.hi, activity_contribution_max(0.1, 0.0, 3.0));
        assert_ne!(cdmax.lo, 0.0);

        // compute_min_activity accumulates the HighsCDouble contributions when the flag is on.
        let idx = [0, 1];
        let vals = [0.1, 0.1];
        dom.col_lower_ = vec![3.0, 3.0];
        let mut ninf: HighsInt = -1;
        let mut act = HighsCDouble::from(0.0);
        dom.compute_min_activity(0, 2, &idx, &vals, &mut ninf, &mut act);
        let mut expect = HighsCDouble::from(0.0);
        expect += cd;
        expect += cd;
        expect.renormalize();
        assert_eq!(ninf, 0);
        assert_eq!((act.hi, act.lo), (expect.hi, expect.lo));

        dom2.col_lower_ = vec![3.0, 3.0];
        let mut act2 = HighsCDouble::from(0.0);
        dom2.compute_min_activity(0, 2, &idx, &vals, &mut ninf, &mut act2);
        let mut expect2 = HighsCDouble::from(0.0);
        expect2 += activity_contribution_min(0.1, 3.0, 10.0);
        expect2 += activity_contribution_min(0.1, 3.0, 10.0);
        expect2.renormalize();
        assert_eq!((act2.hi, act2.lo), (expect2.hi, expect2.lo));
    }

    #[test]
    fn x1_dom_reduced_stack_fast_equals_slow() {
        let mut ms = tiny_mip();
        ms.options_mip_.mip_reduced_stack_fast = false;
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);

        // Build a domain with >=5 bound changes on one column
        dom.change_bound(&mut ms, lower(0, 1.0), HighsDomain_Reason::branching());
        dom.change_bound(&mut ms, lower(0, 2.0), HighsDomain_Reason::unspecified());
        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::unspecified());
        dom.change_bound(&mut ms, lower(0, 4.0), HighsDomain_Reason::unspecified());
        dom.change_bound(&mut ms, lower(0, 5.0), HighsDomain_Reason::unspecified());

        // Slow path (off)
        let mut bp_slow: Vec<HighsInt> = Vec::new();
        let stack_slow = dom.get_reduced_domain_change_stack(&mut bp_slow);

        // Fast path (on)
        dom.x1_reduced_stack_fast = true;
        let mut bp_fast: Vec<HighsInt> = Vec::new();
        let stack_fast = dom.get_reduced_domain_change_stack(&mut bp_fast);

        assert_eq!(stack_slow, stack_fast);
        assert_eq!(bp_slow, bp_fast);
    }

    #[test]
    fn x1_dom_keep_branch_reason_flag() {
        // The stack replays lower(0, 2.0) (unspecified), then a branching change lower(0, 2.0)
        // that is redundant with a local bound on col 0 (colLowerPos_ != -1 when checked).
        // 1.15.1 drops it (skip runs); with the option on and symmetry detected, it is kept.
        let stack = vec![lower(0, 2.0), lower(0, 2.0)];
        let bp = vec![1];

        // Option off: the 1.15.1 skip always runs, the redundant branching change is dropped.
        let mut ms = tiny_mip();
        ms.options_mip_.mip_keep_branch_reason = false;
        ms.mipdata_.detect_symmetries = true;
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);
        dom.set_domain_change_stack_2(&mut ms, &stack, &bp);
        assert_eq!(dom.branch_pos_.len(), 0);

        // Option on, symmetry detected: the skip is bypassed, the change is kept as branching.
        let mut ms2 = tiny_mip();
        ms2.options_mip_.mip_keep_branch_reason = true;
        ms2.mipdata_.detect_symmetries = true;
        let mut dom2 = HighsDomain::new(&ms2);
        dom2.compute_row_activities(&ms2);
        dom2.set_domain_change_stack_2(&mut ms2, &stack, &bp);
        assert_eq!(dom2.branch_pos_.len(), 1);

        // Option on, no symmetry: the skip runs again, as in 1.15.1.
        let mut ms3 = tiny_mip();
        ms3.options_mip_.mip_keep_branch_reason = true;
        ms3.mipdata_.detect_symmetries = false;
        let mut dom3 = HighsDomain::new(&ms3);
        dom3.compute_row_activities(&ms3);
        dom3.set_domain_change_stack_2(&mut ms3, &stack, &bp);
        assert_eq!(dom3.branch_pos_.len(), 0);
    }

    #[test]
    fn x1_dom_is_fixed_to_val() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        // col 0: [0, 10] — not fixed to 5
        assert!(!dom.is_fixed_to_val(0, 5));
        // Fix col 0 to 5
        dom.col_lower_[0] = 5.0;
        dom.col_upper_[0] = 5.0;
        assert!(dom.is_fixed_to_val(0, 5));
        assert!(!dom.is_fixed_to_val(0, 3));
    }

    /// Two integer columns `x0, x1 in [0, 10]` and one row `x0 + x1 <= 4`.
    fn tiny_mip() -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 2];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        ms.mipdata_.a_rstart_ = vec![0, 2];
        ms.mipdata_.a_rindex_ = vec![0, 1];
        ms.mipdata_.a_rvalue_ = vec![1.0, 1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms
    }

    fn lower(col: HighsInt, val: f64) -> HighsDomainChange {
        HighsDomainChange {
            boundval: val,
            column: col,
            boundtype: HighsBoundType::kLower,
        }
    }

    #[test]
    fn x1_dfp_b_redundancy_tracked_only_when_active() {
        // One row x0 + x1 <= 4 (lhs -inf), both columns in [0, 10]. Lowering both upper
        // bounds to 2 makes the row redundant for the DFP propagation (max activity 4).
        // The 1.15.1 path is unchanged; only the DFP redundant-row record depends on the flag.
        for active in [false, true] {
            let mut lp = HighsLp::default();
            lp.num_col_ = 2;
            lp.num_row_ = 1;
            lp.col_lower_ = vec![0.0, 0.0];
            lp.col_upper_ = vec![10.0, 10.0];
            lp.col_cost_ = vec![0.0, 0.0];
            lp.integrality_ = vec![HighsVarType::kContinuous; 2];
            lp.row_lower_ = vec![-kHighsInf];
            lp.row_upper_ = vec![4.0];
            lp.a_matrix_.format_ = MatrixFormat::kColwise;
            lp.a_matrix_.num_col_ = 2;
            lp.a_matrix_.num_row_ = 1;
            lp.a_matrix_.start_ = vec![0, 1, 2];
            lp.a_matrix_.index_ = vec![0, 0];
            lp.a_matrix_.value_ = vec![1.0, 1.0];

            let options = HighsOptions::default();
            let mut cb = HighsCallback::default();
            let solution = HighsSolution::default();
            let ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);

            let mut dom = HighsDomain::default();
            dom.col_lower_ = vec![0.0; 2];
            dom.col_upper_ = vec![10.0; 2];
            dom.activitymin_ = vec![HighsCDouble::from(0.0)];
            dom.activitymax_ = vec![HighsCDouble::from(20.0)];
            dom.activitymininf_ = vec![0];
            dom.activitymaxinf_ = vec![0];
            dom.propagateflags_ = vec![0];
            dom.dual_fix_probing_.redundant_row_flags_ = vec![false; 2];
            dom.capacity_threshold_ = vec![kHighsInf];
            dom.dual_fix_probing_active = active;

            // Max activity 20 -> 12: not redundant yet.
            dom.update_activity_ub_change(&ms, 0, 10.0, 2.0);
            assert!(dom.dual_fix_probing_.redundant_row_inds_.is_empty());
            // Max activity 12 -> 4 <= row_upper 4: redundant only when tracking is active.
            dom.update_activity_ub_change(&ms, 1, 10.0, 2.0);
            if active {
                assert_eq!(dom.dual_fix_probing_.redundant_row_inds_.len(), 1);
                assert!(
                    dom.dual_fix_probing_.redundant_row_flags_[HighsDomain_DualFixProbingRowSide {
                        row: 0,
                        is_rhs: true
                    }
                    .index()]
                );
            } else {
                assert!(dom.dual_fix_probing_.redundant_row_inds_.is_empty());
            }
        }
    }

    #[test]
    fn change_bound_propagate_backtrack_round_trip() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);

        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        assert_eq!(dom.branch_pos_, vec![0]);
        assert!(dom.propagate(&mut ms));
        assert!(!dom.infeasible());
        // x1 <= 4 - 3, reason: model row 0 upper side.
        assert_eq!(dom.col_upper_[1], 1.0);
        // stack: x0 >= 3 (branching), x0 <= 4 and x1 <= 1 (both from row 0 upper side)
        assert_eq!(dom.col_upper_, vec![4.0, 1.0]);
        assert_eq!(dom.domchgstack_.len(), 3);
        assert_eq!(
            dom.domchgreason_[1].r#type,
            HighsDomain_Reason::kModelRowUpper
        );
        assert_eq!(
            dom.domchgreason_[2].r#type,
            HighsDomain_Reason::kModelRowUpper
        );
        assert_eq!(dom.col_upper_pos_, vec![1, 2]);
        // nothing left to propagate
        assert!(!dom.propagate(&mut ms));

        let chg = dom.backtrack(&mut ms);
        assert_eq!((chg.column, chg.boundval), (0, 3.0));
        assert_eq!(chg.boundtype, HighsBoundType::kLower);
        assert!(dom.domchgstack_.is_empty() && dom.prevboundval_.is_empty());
        assert!(dom.branch_pos_.is_empty());
        assert_eq!(dom.col_lower_, vec![0.0, 0.0]);
        assert_eq!(dom.col_upper_, vec![10.0, 10.0]);
        assert_eq!(dom.col_lower_pos_, vec![-1, -1]);
        assert_eq!(dom.col_upper_pos_, vec![-1, -1]);
        assert_eq!(f64::from(dom.activitymin_[0]), 0.0);

        // Replaying the same stack with the branching position reproduces the node.
        dom.set_domain_change_stack_2(&mut ms, &[lower(0, 3.0)], &[0]);
        assert_eq!(dom.col_upper_, vec![4.0, 1.0]);
        assert_eq!(dom.branch_pos_, vec![0]);
    }

    #[test]
    fn conflict_infeasible_row_and_backtrack_to_global() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);

        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        dom.change_bound(&mut ms, lower(1, 2.0), HighsDomain_Reason::branching());
        dom.propagate(&mut ms);
        assert!(dom.infeasible());

        dom.backtrack_to_global(&mut ms);
        assert!(!dom.infeasible());
        assert!(dom.domchgstack_.is_empty() && dom.branch_pos_.is_empty());
        assert_eq!(dom.col_lower_, vec![0.0, 0.0]);
        assert_eq!(f64::from(dom.activitymin_[0]), 0.0);
    }

    #[test]
    fn conflict_pool_registration_and_replay() {
        let mut ms = tiny_mip();
        // Remove the row so only the conflict propagates.
        ms.orig_model_.row_upper_ = vec![kHighsInf];
        ms.mipdata_
            .conflictpools
            .push(Slot::new(HighsConflictPool::new(10, 100)));
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);
        dom.add_conflict_pool(&mut ms, 0);
        assert_eq!(dom.conflict_pool_propagation[0].cursor, 0);

        // Conflict: x0 >= 3 and x1 >= 1 cannot hold together.
        ms.mipdata_.conflictpools[0].add_conflict_from_other_pool(&[lower(0, 3.0), lower(1, 1.0)]);
        // Not applied until the domain next touches its propagators.
        assert!(dom.conflict_pool_propagation[0].conflict_flag_.is_empty());

        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        assert_eq!(dom.conflict_pool_propagation[0].conflict_flag_.len(), 1);
        assert_eq!(
            dom.conflict_pool_propagation[0].cursor,
            ms.mipdata_.conflictpools[0].events.end()
        );
        let prop = &dom.conflict_pool_propagation[0];
        // both literals inactive at registration: two watched, flag counts them
        assert_eq!(prop.conflict_flag_, vec![2]);
        assert_eq!(prop.watched_literals_[0].domchg, lower(0, 3.0));
        assert_eq!(prop.watched_literals_[1].domchg, lower(1, 1.0));
        // HiGHS quirk (linkWatchedLiteral sets the head only inside `if (head != -1)`): the
        // watch lists stay empty, so the conflict never propagates x1 <= 0.
        assert_eq!(prop.col_lower_watched_, vec![-1, -1]);
        assert!(!dom.propagate(&mut ms));
        assert_eq!(dom.col_upper_, vec![10.0, 10.0]);
        assert_eq!(dom.domchgstack_.len(), 1);

        // A copy inherits the cursor and replays a later deletion on its own.
        let mut copy = HighsDomain::copy_from(&mut ms, &dom);
        ms.mipdata_.conflictpools[0].remove_conflict(0);
        copy.backtrack(&mut ms);
        assert_ne!(copy.conflict_pool_propagation[0].conflict_flag_[0] & 8, 0);
        assert_eq!(dom.conflict_pool_propagation[0].conflict_flag_[0] & 8, 0);
        dom.backtrack(&mut ms);
        assert_ne!(dom.conflict_pool_propagation[0].conflict_flag_[0] & 8, 0);
    }
}

#[cfg(test)]
mod h16_row_redundant_stack_tests {
    use super::*;

    // One continuous column in [0, 10] with coefficient 1 in one row with the given sides.
    fn stack_mip(row_lower: f64, row_upper: f64, h16_row_redundant: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.row_lower_ = vec![row_lower];
        lp.row_upper_ = vec![row_upper];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1];
        lp.a_matrix_.index_ = vec![0];
        lp.a_matrix_.value_ = vec![1.0];
        ms.mipdata_.a_rstart_ = vec![0, 1];
        ms.mipdata_.a_rindex_ = vec![0];
        ms.mipdata_.a_rvalue_ = vec![1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.h16_row_redundant = h16_row_redundant;
        // these pin HiGHS' own 0.3 continuous-bound threshold
        ms.options_mip_.mip_cont_prop_improve = 0.3;
        ms
    }

    // Runs propagateIndex on row 0 and returns (capacity threshold, upper changes, lower changes).
    fn run_propagate_index(
        ms: &HighsMipSolver,
        activitymin: f64,
        activitymininf: HighsInt,
        activitymax: f64,
        activitymaxinf: HighsInt,
    ) -> (f64, HighsInt, HighsInt) {
        let mut dom = HighsDomain::new(ms);
        dom.col_lower_ = vec![0.0];
        dom.col_upper_ = vec![10.0];
        dom.activitymin_ = vec![HighsCDouble::from(activitymin)];
        dom.activitymininf_ = vec![activitymininf];
        dom.activitymax_ = vec![HighsCDouble::from(activitymax)];
        dom.activitymaxinf_ = vec![activitymaxinf];
        dom.capacity_threshold_ = vec![123.0];
        dom.prop_row_num_changed_bounds_ = vec![(0, 0)];
        let propagateinds = vec![0];
        let mut changedbounds = vec![HighsDomainChange::default(); 2];
        propagate_model_row_index(&mut dom, ms, &propagateinds, &mut changedbounds, 0);
        let (up, lo) = dom.prop_row_num_changed_bounds_[0];
        (dom.capacity_threshold_[0], up, lo)
    }

    // Row upper 2.0, max activity d1 = round(2.0 + 1e-6): 1.15.1-redundant (d1 <= d1) but not
    // 1.16-redundant (d1 - 2.0 > feastol in CDouble). On the 1.16 path propagateRowUpper tightens
    // the column upper bound from 10 to 2 (one change) and recomputeCapacityThreshold gives
    // max(-feastol, 1 * boundRange(10, 0) = 10 - max(3, 1e-3) = 7, feastol) = 7.
    #[test]
    fn h16_row_redundant_on_propagate_upper_recomputes() {
        let ms = stack_mip(-kHighsInf, 2.0, true);
        assert_eq!(run_propagate_index(&ms, 0.0, 0, 2.0 + 1e-6, 0), (7.0, 1, 0));
    }

    #[test]
    fn h16_row_redundant_off_propagate_upper_skipped() {
        let ms = stack_mip(-kHighsInf, 2.0, false);
        assert_eq!(
            run_propagate_index(&ms, 0.0, 0, 2.0 + 1e-6, 0),
            (123.0, 0, 0)
        );
    }

    // Mirror on the lower side: row lower -2.0, min activity -d1. 1.16 runs propagateRowLower (no
    // change: the implied lower bound -2 - (10 - 10) = -2 is below 0) and recomputes the threshold
    // to 7.
    #[test]
    fn h16_row_redundant_on_propagate_lower_recomputes() {
        let ms = stack_mip(-2.0, kHighsInf, true);
        assert_eq!(
            run_propagate_index(&ms, -(2.0 + 1e-6), 0, 10.0, 0),
            (7.0, 0, 0)
        );
    }

    #[test]
    fn h16_row_redundant_off_propagate_lower_skipped() {
        let ms = stack_mip(-2.0, kHighsInf, false);
        assert_eq!(
            run_propagate_index(&ms, -(2.0 + 1e-6), 0, 10.0, 0),
            (123.0, 0, 0)
        );
    }
}
