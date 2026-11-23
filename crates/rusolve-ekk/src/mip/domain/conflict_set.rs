//! HiGHS source: mip/HighsDomain.cpp (lines 2754-3967)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! PORT NOTE (borrow story): the C++ members `HighsDomain& localdom` / `const HighsDomain&
//! globaldom` are parameters. `localdom` is the domain running `conflictAnalysis*` (outside
//! `ms` or detached from it), so it is `&mut` / `&` next to `ms`. `globaldom` lives in
//! `ms.mipdata_.domains`:
//! * read-only methods take `ms: &HighsMipSolver, localdom, globaldom: &HighsDomain`;
//! * the three methods that add to the conflict pool (`compute_cuts`, `conflict_analysis`,
//!   `conflict_analysis_2`) take `ms: &mut HighsMipSolver, localdom, globaldom: usize,
//!   conflict_pool: usize`. They borrow `let g = &ms.mipdata_.domains[globaldom];` for the
//!   read-only calls and end that borrow before
//!   `detach!(ms.mipdata_.conflictpools[conflict_pool], |p| p.add_conflict_cut(ms, localdom, ..))`.
//!
//! `localdom.mipsolver->mipdata_->x` is `ms.mipdata_.x` (`nodequeue`, `feastol`, `epsilon`,
//! `getRow`). `localdom.cutpoolpropagation[k].cutpool->...` is
//! `ms.mipdata_.cutpools[localdom.cutpoolpropagation[k].cutpool]`.
//! `globaldom.getMinCutActivity(pool, cut)` is `globaldom.get_min_cut_activity_replayed(..)`.
//! `localdom.objProp_.getPropagationConstraint(..)` is
//! `localdom.obj_prop_get_propagation_constraint(ms, ..)`; copy `vals`/`inds` to `Vec`s before
//! the next `localdom` use. `std::set<LocalDomChg>::iterator` in `resolveQueue` is the element
//! value (erase by value; `LocalDomChg` order is by `pos`, unique per frontier).

use crate::prelude::*;

/// Module-private: the `LocalDomChg{pos, HighsDomainChange()}` probe key C++ builds for `find`,
/// `lower_bound` and `upper_bound` (the set is ordered by `pos` only).
fn lookup_key(pos: HighsInt) -> HighsDomain_ConflictSet_LocalDomChg {
    HighsDomain_ConflictSet_LocalDomChg {
        pos,
        domchg: HighsDomainChange::default(),
    }
}

impl HighsDomain_ConflictSet {
    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::ConflictSet (lines 2755-2762)
    /// C++: `ConflictSet(HighsDomain & localdom_, const HighsDomain & globaldom_)`
    // PORT NOTE: the referenced domains are method parameters (module note); nothing to store.
    pub fn new() -> Self {
        Self::default()
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeGeq (lines 2770-2876)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeGeq (lines 3022-3132)
    /// HiGHS 1.16 (#3267): HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeGeq (lines 3022-3132)
    /// C++: `bool explainBoundChangeGeq(const std::set<LocalDomChg> & currentFrontier, const LocalDomChg & domchg, const HighsInt * inds, const double * vals, HighsInt len, double rhs, double maxAct)`
    #[allow(clippy::needless_range_loop)]
    pub fn explain_bound_change_geq(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        globaldom: &HighsDomain,
        current_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        domchg: &HighsDomain_ConflictSet_LocalDomChg,
        inds: &[HighsInt],
        vals: &[f64],
        len: HighsInt,
        rhs: f64,
        max_act: f64,
    ) -> bool {
        // HiGHS 1.16 #3267 (h16_big_act): HighsDomain.cpp line 3026 (1.16.0)
        let h16_big_act = ms.options_mip_.h16_big_act;
        if h16_big_act {
            if max_act == kHighsInf || max_act.abs() > kHighsInvMacheps {
                return false;
            }
        } else if max_act == kHighsInf {
            return false;
        }

        // get the coefficient value of the column for which we want to explain
        // the bound change
        let mut domchg_val: f64 = 0.0;

        // reserve() is capacity only
        self.resolve_buffer.clear();
        let nodequeue = &ms.mipdata_.nodequeue;
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3022-3132 (1.16.0)
        let h16_cdouble = ms.options_mip_.h16_cdouble;
        let mut m = HighsCDouble::from(max_act);
        for i in 0..len as usize {
            let col = inds[i];

            if col == domchg.domchg.column {
                domchg_val = vals[i];
                continue;
            }

            let mut cand = HighsDomain_ConflictSet_ResolveCandidate::default();
            cand.value_pos = i as HighsInt;

            if vals[i] > 0.0 {
                let ub = localdom.get_col_upper_pos(col, domchg.pos, &mut cand.bound_pos);
                if globaldom.col_upper_[col as usize] <= ub || cand.bound_pos == -1 {
                    continue;
                }
                if let Some(it) = current_frontier.get(&lookup_key(cand.bound_pos)) {
                    cand.base_bound = it.domchg.boundval;
                    if cand.base_bound != globaldom.col_upper_[col as usize] {
                        if h16_cdouble {
                            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3055-3056 (1.16.0)
                            m += vals[i]
                                * (HighsCDouble::from(cand.base_bound)
                                    - globaldom.col_upper_[col as usize]);
                        } else {
                            m += vals[i] * (cand.base_bound - globaldom.col_upper_[col as usize]);
                        }
                    }
                    if cand.base_bound <= ub {
                        continue;
                    }
                } else {
                    cand.base_bound = globaldom.col_upper_[col as usize];
                }

                cand.delta = vals[i] * (ub - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    ub,
                    globaldom.col_upper_[col as usize],
                    nodequeue.num_nodes_down(col),
                );
            } else {
                let lb = localdom.get_col_lower_pos(col, domchg.pos, &mut cand.bound_pos);
                if globaldom.col_lower_[col as usize] >= lb || cand.bound_pos == -1 {
                    continue;
                }
                if let Some(it) = current_frontier.get(&lookup_key(cand.bound_pos)) {
                    cand.base_bound = it.domchg.boundval;
                    if cand.base_bound != globaldom.col_lower_[col as usize] {
                        if h16_cdouble {
                            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3073-3074 (1.16.0)
                            m += vals[i]
                                * (HighsCDouble::from(cand.base_bound)
                                    - globaldom.col_lower_[col as usize]);
                        } else {
                            m += vals[i] * (cand.base_bound - globaldom.col_lower_[col as usize]);
                        }
                    }
                    if cand.base_bound >= lb {
                        continue;
                    }
                } else {
                    cand.base_bound = globaldom.col_lower_[col as usize];
                }

                cand.delta = vals[i] * (lb - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    lb,
                    globaldom.col_lower_[col as usize],
                    nodequeue.num_nodes_up(col),
                );
            }

            self.resolve_buffer.push(cand);
        }

        if domchg_val == 0.0 {
            return false;
        }

        pdqsort_by(&mut self.resolve_buffer, |a, b| a.op_lt(b));

        // to explain the bound change we start from the bound constraint,
        // multiply it by the columns coefficient in the constraint. Then the
        // bound constraint is a0 * x0 >= b0 and the constraint
        // a0 * x0 + \sum ai * xi >= b. Let the max activity of \sum ai xi be M,
        // then the constraint yields the bound a0 * x0 >= b - M. If M is
        // sufficiently small the bound constraint is implied. Therefore we
        // decrease M by updating it with the stronger local bounds until
        // M <= b - b0 holds.
        let mut b0 = domchg.domchg.boundval;
        if !ms.is_col_continuous(domchg.domchg.column) {
            // in case of an integral variable the bound was rounded and can be
            // relaxed by 1-feastol. We use 1 - 10 * feastol for numerical safety.
            if domchg.domchg.boundtype == HighsBoundType::kLower {
                b0 -= 1.0 - 10.0 * ms.mipdata_.feastol;
            } else {
                b0 += 1.0 - 10.0 * ms.mipdata_.feastol;
            }
        } else {
            // for a continuous variable we relax the bound by epsilon to
            // accommodate for tiny rounding errors
            if domchg.domchg.boundtype == HighsBoundType::kLower {
                b0 -= ms.mipdata_.epsilon;
            } else {
                b0 += ms.mipdata_.epsilon;
            }
        }

        // now multiply the bound constraint with the coefficient value in the
        // constraint to obtain b0
        b0 *= domchg_val;

        // compute the lower bound of M that is necessary
        let mupper = rhs - b0;

        // M is the current residual activity initially
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3124-3129 (1.16.0)
        if h16_cdouble {
            if domchg_val < 0.0 {
                m -= HighsCDouble::from(domchg_val)
                    * globaldom.col_lower_[domchg.domchg.column as usize];
            } else {
                m -= HighsCDouble::from(domchg_val)
                    * globaldom.col_upper_[domchg.domchg.column as usize];
            }
        } else if domchg_val < 0.0 {
            m -= domchg_val * globaldom.col_lower_[domchg.domchg.column as usize];
        } else {
            m -= domchg_val * globaldom.col_upper_[domchg.domchg.column as usize];
        }

        self.resolve_linear_geq(ms, localdom, m, mupper, vals)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeLeq (lines 2878-2985)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeLeq (lines 3134-3245)
    /// HiGHS 1.16 (#3267): HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeLeq (lines 3134-3245)
    /// C++: `bool explainBoundChangeLeq(const std::set<LocalDomChg> & currentFrontier, const LocalDomChg & domchg, const HighsInt * inds, const double * vals, HighsInt len, double rhs, double minAct)`
    #[allow(clippy::needless_range_loop)]
    pub fn explain_bound_change_leq(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        globaldom: &HighsDomain,
        current_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        domchg: &HighsDomain_ConflictSet_LocalDomChg,
        inds: &[HighsInt],
        vals: &[f64],
        len: HighsInt,
        rhs: f64,
        min_act: f64,
    ) -> bool {
        // HiGHS 1.16 #3267 (h16_big_act): HighsDomain.cpp line 3138 (1.16.0)
        let h16_big_act = ms.options_mip_.h16_big_act;
        if h16_big_act {
            if min_act == -kHighsInf || min_act.abs() > kHighsInvMacheps {
                return false;
            }
        } else if min_act == -kHighsInf {
            return false;
        }
        // get the coefficient value of the column for which we want to explain
        // the bound change
        let mut domchg_val: f64 = 0.0;

        // reserve() is capacity only
        self.resolve_buffer.clear();
        let nodequeue = &ms.mipdata_.nodequeue;
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3134-3245 (1.16.0)
        let h16_cdouble = ms.options_mip_.h16_cdouble;
        let mut m = HighsCDouble::from(min_act);
        for i in 0..len as usize {
            let col = inds[i];

            if col == domchg.domchg.column {
                domchg_val = vals[i];
                continue;
            }

            let mut cand = HighsDomain_ConflictSet_ResolveCandidate::default();
            cand.value_pos = i as HighsInt;

            if vals[i] > 0.0 {
                let lb = localdom.get_col_lower_pos(col, domchg.pos, &mut cand.bound_pos);
                if globaldom.col_lower_[col as usize] >= lb || cand.bound_pos == -1 {
                    continue;
                }
                if let Some(it) = current_frontier.get(&lookup_key(cand.bound_pos)) {
                    cand.base_bound = it.domchg.boundval;
                    if cand.base_bound != globaldom.col_lower_[col as usize] {
                        if h16_cdouble {
                            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3167-3168 (1.16.0)
                            m += vals[i]
                                * (HighsCDouble::from(cand.base_bound)
                                    - globaldom.col_lower_[col as usize]);
                        } else {
                            m += vals[i] * (cand.base_bound - globaldom.col_lower_[col as usize]);
                        }
                    }
                    if cand.base_bound >= lb {
                        continue;
                    }
                } else {
                    cand.base_bound = globaldom.col_lower_[col as usize];
                }

                cand.delta = vals[i] * (lb - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    lb,
                    globaldom.col_lower_[col as usize],
                    nodequeue.num_nodes_up(col),
                );
            } else {
                let ub = localdom.get_col_upper_pos(col, domchg.pos, &mut cand.bound_pos);
                if globaldom.col_upper_[col as usize] <= ub || cand.bound_pos == -1 {
                    continue;
                }
                if let Some(it) = current_frontier.get(&lookup_key(cand.bound_pos)) {
                    cand.base_bound = it.domchg.boundval;
                    if cand.base_bound != globaldom.col_upper_[col as usize] {
                        if h16_cdouble {
                            // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3184-3185 (1.16.0)
                            m += vals[i]
                                * (HighsCDouble::from(cand.base_bound)
                                    - globaldom.col_upper_[col as usize]);
                        } else {
                            m += vals[i] * (cand.base_bound - globaldom.col_upper_[col as usize]);
                        }
                    }
                    if cand.base_bound <= ub {
                        continue;
                    }
                } else {
                    cand.base_bound = globaldom.col_upper_[col as usize];
                }

                cand.delta = vals[i] * (ub - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    ub,
                    globaldom.col_upper_[col as usize],
                    nodequeue.num_nodes_down(col),
                );
            }

            self.resolve_buffer.push(cand);
        }

        if domchg_val == 0.0 {
            return false;
        }

        pdqsort_by(&mut self.resolve_buffer, |a, b| a.op_lt(b));

        debug_assert!(domchg_val != 0.0);

        // to explain the bound change we start from the bound constraint,
        // multiply it by the columns coefficient in the constraint. Then the
        // bound constraint is a0 * x0 <= b0 and the constraint
        // a0 * x0 + \sum ai * xi <= b. Let the min activity of \sum ai xi be M,
        // then the constraint yields the bound a0 * x0 <= b - M. If M is
        // sufficiently large the bound constraint is implied. Therefore we
        // increase M by updating it with the stronger local bounds until
        // M >= b - b0 holds.
        let mut b0 = domchg.domchg.boundval;
        if !ms.is_col_continuous(domchg.domchg.column) {
            // in case of an integral variable the bound was rounded and can be
            // relaxed by 1-feastol. We use 1 - 10 * feastol for numerical safety
            if domchg.domchg.boundtype == HighsBoundType::kLower {
                b0 -= 1.0 - 10.0 * ms.mipdata_.feastol;
            } else {
                b0 += 1.0 - 10.0 * ms.mipdata_.feastol;
            }
        } else {
            // for a continuous variable we relax the bound by epsilon to
            // accommodate for tiny rounding errors
            if domchg.domchg.boundtype == HighsBoundType::kLower {
                b0 -= ms.mipdata_.epsilon;
            } else {
                b0 += ms.mipdata_.epsilon;
            }
        }

        // now multiply the bound constraint with the coefficient value in the
        // constraint to obtain b0
        b0 *= domchg_val;

        // compute the lower bound of M that is necessary
        let mlower = rhs - b0;

        // M is the global residual activity initially
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3237-3242 (1.16.0)
        if h16_cdouble {
            if domchg_val < 0.0 {
                m -= HighsCDouble::from(domchg_val)
                    * globaldom.col_upper_[domchg.domchg.column as usize];
            } else {
                m -= HighsCDouble::from(domchg_val)
                    * globaldom.col_lower_[domchg.domchg.column as usize];
            }
        } else if domchg_val < 0.0 {
            m -= domchg_val * globaldom.col_upper_[domchg.domchg.column as usize];
        } else {
            m -= domchg_val * globaldom.col_lower_[domchg.domchg.column as usize];
        }

        self.resolve_linear_leq(ms, localdom, m, mlower, vals)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::resolveLinearGeq (lines 2987-3092)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::HighsDomain::ConflictSet::resolveLinearGeq (lines 3247-3359)
    /// HiGHS 1.16 (#3196): HighsDomain.cpp::HighsDomain::ConflictSet::resolveLinearGeq (lines 3247-3359)
    /// C++: `bool resolveLinearGeq(HighsCDouble M, double Mupper, const double * vals)`
    #[allow(clippy::needless_range_loop)]
    pub fn resolve_linear_geq(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        m: HighsCDouble,
        mupper: f64,
        vals: &[f64],
    ) -> bool {
        let mut m = m;
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3247-3359 (1.16.0)
        let h16_cdouble = ms.options_mip_.h16_cdouble;
        // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3247-3359 (1.16.0)
        let h16_bound_chain = ms.options_mip_.h16_bound_chain;
        self.resolved_domain_changes.clear();
        let mut covered: f64 = f64::from(m - mupper);
        if covered > 0.0 {
            for reason_domchg in self.resolve_buffer.iter() {
                let locdomchg = HighsDomain_ConflictSet_LocalDomChg {
                    pos: reason_domchg.bound_pos,
                    domchg: localdom.domchgstack_[reason_domchg.bound_pos as usize].clone(),
                };
                m += reason_domchg.delta;
                self.resolved_domain_changes.push(locdomchg);
                debug_assert!(self.resolved_domain_changes.last().unwrap().pos >= 0);
                debug_assert!(
                    self.resolved_domain_changes.last().unwrap().pos
                        < localdom.domchgstack_.len() as HighsInt
                );
                covered = f64::from(m - mupper);
                if covered <= 0.0 {
                    break;
                }
            }

            if covered > 0.0 {
                return false;
            }

            if covered < -localdom.feastol(ms) {
                // there is room for relaxing bounds / dropping unneeded bound changes
                // from the explanation
                // C++ `for (k = size - 1; k >= 0; --k)`: the size is read once, and each
                // removal takes out the last index, which is never above k.
                let n = self.resolved_domain_changes.len();
                for k in (0..n).rev() {
                    let reason_domchg = self.resolve_buffer[k];
                    let i = reason_domchg.value_pos as usize;
                    let col = self.resolved_domain_changes[k].domchg.column;
                    if self.resolved_domain_changes[k].domchg.boundtype == HighsBoundType::kLower {
                        let lb = self.resolved_domain_changes[k].domchg.boundval;
                        let glb = reason_domchg.base_bound;
                        let mut relax_lb: f64 =
                            f64::from(((mupper - (m - reason_domchg.delta)) / vals[i]) + glb);
                        if !ms.is_col_continuous(col) {
                            relax_lb = relax_lb.ceil();
                        }

                        if relax_lb - lb >= -localdom.feastol(ms) {
                            continue;
                        }

                        self.resolved_domain_changes[k].domchg.boundval = relax_lb;

                        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3290-3292 (1.16.0)
                        let droppable = if h16_cdouble {
                            relax_lb - glb <= ms.mipdata_.epsilon
                                && m - reason_domchg.delta - mupper <= localdom.feastol(ms)
                        } else {
                            relax_lb - glb <= ms.mipdata_.epsilon
                        };
                        if droppable {
                            // domain change can be fully removed from conflict
                            let last = self.resolved_domain_changes.len() - 1;
                            self.resolved_domain_changes.swap(last, k);
                            self.resolved_domain_changes.truncate(last);

                            m -= reason_domchg.delta;
                        } else {
                            // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3302-3303 (1.16.0)
                            if h16_bound_chain {
                                if !localdom.walk_bound_chain(
                                    &mut self.resolved_domain_changes[k].pos,
                                    relax_lb,
                                    1,
                                ) {
                                    return false;
                                }
                            } else {
                                while relax_lb
                                    <= localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .0
                                {
                                    self.resolved_domain_changes[k].pos = localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .1;
                                }
                            }

                            // bound can be relaxed
                            if h16_cdouble {
                                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 3306 (1.16.0)
                                m += vals[i] * (HighsCDouble::from(relax_lb) - lb);
                            } else {
                                m += vals[i] * (relax_lb - lb);
                            }
                        }

                        covered = f64::from(m - mupper);
                        if covered >= -localdom.feastol(ms) {
                            break;
                        }
                    } else {
                        let ub = self.resolved_domain_changes[k].domchg.boundval;
                        let gub = reason_domchg.base_bound;
                        let mut relax_ub: f64 =
                            f64::from(((mupper - (m - reason_domchg.delta)) / vals[i]) + gub);
                        if !ms.is_col_continuous(col) {
                            relax_ub = relax_ub.floor();
                        }

                        if relax_ub - ub <= localdom.feastol(ms) {
                            continue;
                        }

                        self.resolved_domain_changes[k].domchg.boundval = relax_ub;

                        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3325-3326 (1.16.0)
                        let droppable = if h16_cdouble {
                            relax_ub - gub >= -ms.mipdata_.epsilon
                                && m - reason_domchg.delta - mupper <= localdom.feastol(ms)
                        } else {
                            relax_ub - gub >= -ms.mipdata_.epsilon
                        };
                        if droppable {
                            // domain change can be fully removed from conflict
                            let last = self.resolved_domain_changes.len() - 1;
                            self.resolved_domain_changes.swap(last, k);
                            self.resolved_domain_changes.truncate(last);

                            m -= reason_domchg.delta;
                        } else {
                            // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3337-3338 (1.16.0)
                            if h16_bound_chain {
                                if !localdom.walk_bound_chain(
                                    &mut self.resolved_domain_changes[k].pos,
                                    relax_ub,
                                    -1,
                                ) {
                                    return false;
                                }
                            } else {
                                // bound can be relaxed
                                while relax_ub
                                    >= localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .0
                                {
                                    self.resolved_domain_changes[k].pos = localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .1;
                                }
                            }

                            if h16_cdouble {
                                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 3340 (1.16.0)
                                m += vals[i] * (HighsCDouble::from(relax_ub) - ub);
                            } else {
                                m += vals[i] * (relax_ub - ub);
                            }
                        }

                        covered = f64::from(m - mupper);
                        if covered >= -localdom.feastol(ms) {
                            break;
                        }
                    }
                }

                debug_assert!(covered <= ms.mipdata_.feastol);
            }
        }

        true
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::resolveLinearLeq (lines 3094-3203)
    /// HiGHS 1.16 (#3149): HighsDomain.cpp::HighsDomain::ConflictSet::resolveLinearLeq (lines 3361-3477)
    /// HiGHS 1.16 (#3196): HighsDomain.cpp::HighsDomain::ConflictSet::resolveLinearLeq (lines 3361-3477)
    /// C++: `bool resolveLinearLeq(HighsCDouble M, double Mlower, const double * vals)`
    #[allow(clippy::needless_range_loop)]
    pub fn resolve_linear_leq(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        m: HighsCDouble,
        mlower: f64,
        vals: &[f64],
    ) -> bool {
        let mut m = m;
        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3361-3477 (1.16.0)
        let h16_cdouble = ms.options_mip_.h16_cdouble;
        // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3361-3477 (1.16.0)
        let h16_bound_chain = ms.options_mip_.h16_bound_chain;
        self.resolved_domain_changes.clear();
        let mut covered: f64 = f64::from(m - mlower);
        if covered < 0.0 {
            for reason_domchg in self.resolve_buffer.iter() {
                let locdomchg = HighsDomain_ConflictSet_LocalDomChg {
                    pos: reason_domchg.bound_pos,
                    domchg: localdom.domchgstack_[reason_domchg.bound_pos as usize].clone(),
                };
                m += reason_domchg.delta;
                self.resolved_domain_changes.push(locdomchg);
                debug_assert!(self.resolved_domain_changes.last().unwrap().pos >= 0);
                debug_assert!(
                    self.resolved_domain_changes.last().unwrap().pos
                        < localdom.domchgstack_.len() as HighsInt
                );
                covered = f64::from(m - mlower);
                if covered >= 0.0 {
                    break;
                }
            }

            if covered < 0.0 {
                // local bounds reach only a value below the required one
                return false;
            }

            if covered > localdom.feastol(ms) {
                // there is room for relaxing bounds / dropping unneeded bound changes
                // from the explanation
                // C++ `for (k = size - 1; k >= 0; --k)`: see resolve_linear_geq.
                let n = self.resolved_domain_changes.len();
                for k in (0..n).rev() {
                    let reason_domchg = self.resolve_buffer[k];
                    let i = reason_domchg.value_pos as usize;
                    let col = self.resolved_domain_changes[k].domchg.column;
                    if self.resolved_domain_changes[k].domchg.boundtype == HighsBoundType::kLower {
                        let lb = self.resolved_domain_changes[k].domchg.boundval;
                        let glb = reason_domchg.base_bound;
                        let mut relax_lb: f64 =
                            f64::from(((mlower - (m - reason_domchg.delta)) / vals[i]) + glb);
                        if !ms.is_col_continuous(col) {
                            relax_lb = relax_lb.ceil();
                        }

                        if relax_lb - lb >= -localdom.feastol(ms) {
                            continue;
                        }

                        self.resolved_domain_changes[k].domchg.boundval = relax_lb;

                        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3410-3411 (1.16.0)
                        let droppable = if h16_cdouble {
                            relax_lb - glb <= ms.mipdata_.epsilon
                                && m - reason_domchg.delta - mlower >= -localdom.feastol(ms)
                        } else {
                            relax_lb - glb <= ms.mipdata_.epsilon
                        };
                        if droppable {
                            // domain change can be fully removed from conflict
                            let last = self.resolved_domain_changes.len() - 1;
                            self.resolved_domain_changes.swap(last, k);
                            self.resolved_domain_changes.truncate(last);

                            m -= reason_domchg.delta;
                        } else {
                            // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3422-3423 (1.16.0)
                            if h16_bound_chain {
                                if !localdom.walk_bound_chain(
                                    &mut self.resolved_domain_changes[k].pos,
                                    relax_lb,
                                    1,
                                ) {
                                    return false;
                                }
                            } else {
                                // bound can be relaxed
                                while relax_lb
                                    <= localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .0
                                {
                                    self.resolved_domain_changes[k].pos = localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .1;
                                }
                            }

                            if h16_cdouble {
                                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 3425 (1.16.0)
                                m += vals[i] * (HighsCDouble::from(relax_lb) - lb);
                            } else {
                                m += vals[i] * (relax_lb - lb);
                            }
                        }

                        covered = f64::from(m - mlower);
                        if covered <= localdom.feastol(ms) {
                            break;
                        }
                    } else {
                        let ub = self.resolved_domain_changes[k].domchg.boundval;
                        let gub = reason_domchg.base_bound;
                        let mut relax_ub: f64 =
                            f64::from(((mlower - (m - reason_domchg.delta)) / vals[i]) + gub);
                        if !ms.is_col_continuous(col) {
                            relax_ub = relax_ub.floor();
                        }

                        if relax_ub - ub <= localdom.feastol(ms) {
                            continue;
                        }

                        self.resolved_domain_changes[k].domchg.boundval = relax_ub;

                        // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp lines 3444-3445 (1.16.0)
                        let droppable = if h16_cdouble {
                            relax_ub - gub >= -ms.mipdata_.epsilon
                                && m - reason_domchg.delta - mlower >= -localdom.feastol(ms)
                        } else {
                            relax_ub - gub >= -ms.mipdata_.epsilon
                        };
                        if droppable {
                            // domain change can be fully removed from conflict
                            let last = self.resolved_domain_changes.len() - 1;
                            self.resolved_domain_changes.swap(last, k);
                            self.resolved_domain_changes.truncate(last);

                            m -= reason_domchg.delta;
                        } else {
                            // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3456-3457 (1.16.0)
                            if h16_bound_chain {
                                if !localdom.walk_bound_chain(
                                    &mut self.resolved_domain_changes[k].pos,
                                    relax_ub,
                                    -1,
                                ) {
                                    return false;
                                }
                            } else {
                                // bound can be relaxed
                                while relax_ub
                                    >= localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .0
                                {
                                    self.resolved_domain_changes[k].pos = localdom.prevboundval_
                                        [self.resolved_domain_changes[k].pos as usize]
                                        .1;
                                }
                            }

                            if h16_cdouble {
                                // HiGHS 1.16 #3149 (h16_cdouble): HighsDomain.cpp line 3459 (1.16.0)
                                m += vals[i] * (HighsCDouble::from(relax_ub) - ub);
                            } else {
                                m += vals[i] * (relax_ub - ub);
                            }
                        }

                        covered = f64::from(m - mlower);
                        if covered <= localdom.feastol(ms) {
                            break;
                        }
                    }
                }

                debug_assert!(covered >= -ms.mipdata_.feastol);
            }
        }
        true
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibility (lines 3205-3335)
    /// C++: `bool explainInfeasibility()`
    #[allow(clippy::needless_return)]
    pub fn explain_infeasibility(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &mut HighsDomain,
        globaldom: &HighsDomain,
    ) -> bool {
        let reason = localdom.infeasible_reason;
        match reason.r#type {
            HighsDomain_Reason::kUnknown | HighsDomain_Reason::kBranching => return false,
            HighsDomain_Reason::kConflictingBounds => {
                self.resolved_domain_changes.clear();
                let conflicting_bound_pos = reason.index;
                let col = localdom.domchgstack_[conflicting_bound_pos as usize].column;

                self.resolved_domain_changes
                    .push(HighsDomain_ConflictSet_LocalDomChg {
                        pos: conflicting_bound_pos,
                        domchg: localdom.domchgstack_[conflicting_bound_pos as usize].clone(),
                    });

                // C++ `HighsInt otherBoundPos;` is an out-parameter of getColUpperPos /
                // getColLowerPos, which always assigns it; initialised here so it can be borrowed.
                let mut other_bound_pos: HighsInt = 0;
                if localdom.domchgstack_[conflicting_bound_pos as usize].boundtype
                    == HighsBoundType::kLower
                {
                    let ub = localdom.get_col_upper_pos(
                        col,
                        conflicting_bound_pos,
                        &mut other_bound_pos,
                    );
                    debug_assert!(
                        localdom.domchgstack_[conflicting_bound_pos as usize].boundval - ub
                            > ms.mipdata_.feastol
                    );
                } else {
                    let lb = localdom.get_col_lower_pos(
                        col,
                        conflicting_bound_pos,
                        &mut other_bound_pos,
                    );
                    debug_assert!(
                        localdom.domchgstack_[conflicting_bound_pos as usize].boundval - lb
                            < -ms.mipdata_.feastol
                    );
                }
                if other_bound_pos != -1 {
                    self.resolved_domain_changes
                        .push(HighsDomain_ConflictSet_LocalDomChg {
                            pos: other_bound_pos,
                            domchg: localdom.domchgstack_[other_bound_pos as usize].clone(),
                        });
                }
                return true;
            }
            HighsDomain_Reason::kCliqueTable => {
                debug_assert!(false, "clique-table infeasibility is not explained");
                return false;
            }
            HighsDomain_Reason::kModelRowLower => {
                let row_index = reason.index;

                // retrieve the matrix values of the cut
                let mut len: HighsInt = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                ms.get_row(row_index, &mut len, &mut inds, &mut vals);

                let max_act = globaldom.get_max_activity(row_index);

                return self.explain_infeasibility_geq(
                    ms,
                    localdom,
                    globaldom,
                    inds,
                    vals,
                    len,
                    ms.row_lower_2(row_index),
                    max_act,
                );
            }
            HighsDomain_Reason::kModelRowUpper => {
                let row_index = reason.index;

                // retrieve the matrix values of the cut
                let mut len: HighsInt = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                ms.get_row(row_index, &mut len, &mut inds, &mut vals);

                let min_act = globaldom.get_min_activity(row_index);

                return self.explain_infeasibility_leq(
                    ms,
                    localdom,
                    globaldom,
                    inds,
                    vals,
                    len,
                    ms.row_upper_2(row_index),
                    min_act,
                );
            }
            HighsDomain_Reason::kObjective => {
                // copy the constraint out of `localdom` before it is used again
                let infeasible_pos = localdom.infeasible_pos;
                let (ov, oi, olen, orhs) = {
                    let mut vals: &[f64] = &[];
                    let mut inds: &[HighsInt] = &[];
                    let mut len: HighsInt = 0;
                    let mut rhs: f64 = 0.0;
                    localdom.obj_prop_get_propagation_constraint(
                        ms,
                        infeasible_pos,
                        &mut vals,
                        &mut inds,
                        &mut len,
                        &mut rhs,
                        -1,
                    );
                    (vals.to_vec(), inds.to_vec(), len, rhs)
                };

                let mut ninfmin: HighsInt = 0;
                let mut min_act = HighsCDouble::from(0.0);
                globaldom.compute_min_activity(0, olen, &oi, &ov, &mut ninfmin, &mut min_act);
                // This case should only happen when a globally unbounded column is
                // bounded in the local domain, e.g., by a branching choice in some
                // heuristic.
                if ninfmin > 0 {
                    return false;
                }

                return self.explain_infeasibility_leq(
                    ms,
                    localdom,
                    globaldom,
                    &oi,
                    &ov,
                    olen,
                    orhs,
                    f64::from(min_act),
                );
            }
            _ => {
                debug_assert!(reason.r#type >= 0);
                debug_assert!(
                    reason.r#type
                        < (localdom.cutpoolpropagation.len()
                            + localdom.conflict_pool_propagation.len())
                            as HighsInt
                );

                if reason.r#type < localdom.cutpoolpropagation.len() as HighsInt {
                    let cutpool_index = reason.r#type as usize;
                    let cut_index = reason.index;
                    let pool = localdom.cutpoolpropagation[cutpool_index].cutpool;

                    // retrieve the matrix values of the cut
                    let mut len: HighsInt = 0;
                    let mut inds: &[HighsInt] = &[];
                    let mut vals: &[f64] = &[];
                    ms.mipdata_.cutpools[pool].get_cut(cut_index, &mut len, &mut inds, &mut vals);

                    let min_act = globaldom.get_min_cut_activity_replayed(pool, cut_index);
                    let rhs = ms.mipdata_.cutpools[pool].get_rhs()[cut_index as usize];

                    return self.explain_infeasibility_leq(
                        ms, localdom, globaldom, inds, vals, len, rhs, min_act,
                    );
                } else {
                    let conflict_pool_index =
                        (reason.r#type - localdom.cutpoolpropagation.len() as HighsInt) as usize;
                    let conflict_index = reason.index;

                    if (localdom.conflict_pool_propagation[conflict_pool_index].conflict_flag_
                        [conflict_index as usize]
                        & 8)
                        != 0
                    {
                        return false;
                    }

                    // retrieve the conflict entries
                    let pool = localdom.conflict_pool_propagation[conflict_pool_index].conflictpool;
                    let (first, second) = ms.mipdata_.conflictpools[pool].get_conflict_ranges()
                        [conflict_index as usize];
                    let conflict = &ms.mipdata_.conflictpools[pool].get_conflict_entry_vector()
                        [first as usize..second as usize];

                    return self.explain_infeasibility_conflict(ms, localdom, globaldom, conflict);
                }
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibilityConflict (lines 3337-3374)
    /// HiGHS 1.16 (#3196): HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibilityConflict (lines 3612-3639)
    /// C++: `bool explainInfeasibilityConflict(const HighsDomainChange * conflict, HighsInt len)`
    pub fn explain_infeasibility_conflict(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        globaldom: &HighsDomain,
        conflict: &[HighsDomainChange],
    ) -> bool {
        // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3612-3639 (1.16.0)
        let h16_bound_chain = ms.options_mip_.h16_bound_chain;
        self.resolved_domain_changes.clear();
        for entry in conflict.iter() {
            if globaldom.is_active(entry) {
                continue;
            }

            // C++ `HighsInt pos;` is an out-parameter of getColLowerPos / getColUpperPos.
            let mut pos: HighsInt = 0;
            if entry.boundtype == HighsBoundType::kLower {
                let lb =
                    localdom.get_col_lower_pos(entry.column, localdom.infeasible_pos, &mut pos);
                if pos == -1 || lb < entry.boundval {
                    return false;
                }

                // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3624-3625 (1.16.0)
                if h16_bound_chain {
                    if !localdom.walk_bound_chain(&mut pos, entry.boundval, 1) {
                        return false;
                    }
                } else {
                    while localdom.prevboundval_[pos as usize].0 >= entry.boundval {
                        pos = localdom.prevboundval_[pos as usize].1;
                        // since we checked that the bound value is not active globally and is
                        // active for the local domain at pos
                        // pos should never become -1
                        debug_assert!(pos != -1);
                    }
                }
            } else {
                let ub =
                    localdom.get_col_upper_pos(entry.column, localdom.infeasible_pos, &mut pos);
                if pos == -1 || ub > entry.boundval {
                    return false;
                }

                // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3631-3632 (1.16.0)
                if h16_bound_chain {
                    if !localdom.walk_bound_chain(&mut pos, entry.boundval, -1) {
                        return false;
                    }
                } else {
                    while localdom.prevboundval_[pos as usize].0 <= entry.boundval {
                        pos = localdom.prevboundval_[pos as usize].1;
                        // since we checked that the bound value is not active globally and is
                        // active for the local domain at pos
                        // pos should never become -1
                        debug_assert!(pos != -1);
                    }
                }
            }

            self.resolved_domain_changes
                .push(HighsDomain_ConflictSet_LocalDomChg {
                    pos,
                    domchg: entry.clone(),
                });
        }

        let _ = ms;
        true
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibilityGeq (lines 3376-3421)
    /// HiGHS 1.16 (#3267): HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibilityGeq (lines 3641-3685)
    /// C++: `bool explainInfeasibilityGeq(const HighsInt * inds, const double * vals, HighsInt len, double rhs, double maxAct)`
    #[allow(clippy::needless_range_loop)]
    pub fn explain_infeasibility_geq(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        globaldom: &HighsDomain,
        inds: &[HighsInt],
        vals: &[f64],
        len: HighsInt,
        rhs: f64,
        max_act: f64,
    ) -> bool {
        // HiGHS 1.16 #3267 (h16_big_act): HighsDomain.cpp line 3645 (1.16.0)
        let h16_big_act = ms.options_mip_.h16_big_act;
        if h16_big_act {
            if max_act == kHighsInf || max_act.abs() > kHighsInvMacheps {
                return false;
            }
        } else if max_act == kHighsInf {
            return false;
        }

        let mut infeasible_pos: HighsInt = kHighsIInf;
        if localdom.infeasible_ {
            infeasible_pos = localdom.infeasible_pos;
        }

        // reserve() is capacity only
        self.resolve_buffer.clear();
        let nodequeue = &ms.mipdata_.nodequeue;
        for i in 0..len as usize {
            let col = inds[i];

            let mut cand = HighsDomain_ConflictSet_ResolveCandidate::default();
            cand.value_pos = i as HighsInt;

            if vals[i] > 0.0 {
                let ub = localdom.get_col_upper_pos(col, infeasible_pos, &mut cand.bound_pos);
                cand.base_bound = globaldom.col_upper_[col as usize];
                if cand.base_bound <= ub || cand.bound_pos == -1 {
                    continue;
                }
                cand.delta = vals[i] * (ub - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    ub,
                    globaldom.col_upper_[col as usize],
                    nodequeue.num_nodes_down(col),
                );
            } else {
                let lb = localdom.get_col_lower_pos(col, infeasible_pos, &mut cand.bound_pos);
                cand.base_bound = globaldom.col_lower_[col as usize];
                if cand.base_bound >= lb || cand.bound_pos == -1 {
                    continue;
                }
                cand.delta = vals[i] * (lb - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    lb,
                    globaldom.col_lower_[col as usize],
                    nodequeue.num_nodes_up(col),
                );
            }

            self.resolve_buffer.push(cand);
        }

        pdqsort_by(&mut self.resolve_buffer, |a, b| a.op_lt(b));

        // compute the lower bound of M that is necessary
        let mupper = rhs - cpp_max(10.0, rhs.abs()) * ms.mipdata_.feastol;

        debug_assert!(self.reason_side_frontier.is_empty());
        self.resolve_linear_geq(ms, localdom, HighsCDouble::from(max_act), mupper, vals)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibilityLeq (lines 3423-3467)
    /// HiGHS 1.16 (#3267): HighsDomain.cpp::HighsDomain::ConflictSet::explainInfeasibilityLeq (lines 3687-3730)
    /// C++: `bool explainInfeasibilityLeq(const HighsInt * inds, const double * vals, HighsInt len, double rhs, double minAct)`
    #[allow(clippy::needless_range_loop)]
    pub fn explain_infeasibility_leq(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        globaldom: &HighsDomain,
        inds: &[HighsInt],
        vals: &[f64],
        len: HighsInt,
        rhs: f64,
        min_act: f64,
    ) -> bool {
        // HiGHS 1.16 #3267 (h16_big_act): HighsDomain.cpp line 3691 (1.16.0)
        let h16_big_act = ms.options_mip_.h16_big_act;
        if h16_big_act {
            if min_act == -kHighsInf || min_act.abs() > kHighsInvMacheps {
                return false;
            }
        } else if min_act == -kHighsInf {
            return false;
        }

        let mut infeasible_pos: HighsInt = kHighsIInf;
        if localdom.infeasible_ {
            infeasible_pos = localdom.infeasible_pos;
        }

        // reserve() is capacity only
        self.resolve_buffer.clear();
        let nodequeue = &ms.mipdata_.nodequeue;
        for i in 0..len as usize {
            let col = inds[i];

            let mut cand = HighsDomain_ConflictSet_ResolveCandidate::default();
            cand.value_pos = i as HighsInt;

            if vals[i] > 0.0 {
                let lb = localdom.get_col_lower_pos(col, infeasible_pos, &mut cand.bound_pos);
                cand.base_bound = globaldom.col_lower_[col as usize];
                if cand.base_bound >= lb || cand.bound_pos == -1 {
                    continue;
                }
                cand.delta = vals[i] * (lb - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    lb,
                    globaldom.col_lower_[col as usize],
                    nodequeue.num_nodes_up(col),
                );
            } else {
                let ub = localdom.get_col_upper_pos(col, infeasible_pos, &mut cand.bound_pos);
                cand.base_bound = globaldom.col_upper_[col as usize];
                if cand.base_bound <= ub || cand.bound_pos == -1 {
                    continue;
                }
                cand.delta = vals[i] * (ub - cand.base_bound);
                cand.prio = compute_prio(
                    vals[i],
                    ub,
                    globaldom.col_upper_[col as usize],
                    nodequeue.num_nodes_down(col),
                );
            }

            self.resolve_buffer.push(cand);
        }

        pdqsort_by(&mut self.resolve_buffer, |a, b| a.op_lt(b));

        // compute the lower bound of M that is necessary
        let mlower = rhs + cpp_max(10.0, rhs.abs()) * ms.mipdata_.feastol;

        self.resolve_linear_leq(ms, localdom, HighsCDouble::from(min_act), mlower, vals)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChange (lines 3469-3600)
    /// C++: `bool explainBoundChange(const std::set<LocalDomChg> & currentFrontier, LocalDomChg domchg)`
    /// HiGHS 1.16 (#3275): HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChange (lines 3732-3861), kObjective case (lines 3794-3811)
    #[allow(clippy::needless_return)]
    pub fn explain_bound_change(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &mut HighsDomain,
        globaldom: &HighsDomain,
        current_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        domchg: HighsDomain_ConflictSet_LocalDomChg,
    ) -> bool {
        // PORT NOTE: C++ `addConflictCut` notifies `localdom`'s ConflictPoolPropagation
        // synchronously, so a conflict slot reused by a cut added earlier in this analysis
        // (deletedConflicts_ recycling) already has a fresh `conflictFlag_` when it is read
        // below. The lazy event log is replayed here to observe the same state.
        localdom.replay_pool_events(ms);
        let dpos = domchg.pos as usize;
        let reason = localdom.domchgreason_[dpos];
        match reason.r#type {
            HighsDomain_Reason::kUnknown
            | HighsDomain_Reason::kBranching
            | HighsDomain_Reason::kConflictingBounds => return false,
            HighsDomain_Reason::kCliqueTable => {
                let col = reason.index >> 1;
                let val = reason.index & 1;
                self.resolved_domain_changes.clear();
                // C++ `HighsInt boundPos;` is an out-parameter of getColLowerPos / getColUpperPos.
                let mut bound_pos: HighsInt = 0;
                if val != 0 {
                    debug_assert!(localdom.col_lower_pos_[col as usize] >= 0);
                    debug_assert!(
                        localdom.col_lower_pos_[col as usize]
                            < localdom.domchgstack_.len() as HighsInt
                    );

                    localdom.get_col_lower_pos(col, domchg.pos, &mut bound_pos);
                } else {
                    debug_assert!(localdom.col_upper_pos_[col as usize] >= 0);
                    debug_assert!(
                        localdom.col_upper_pos_[col as usize]
                            < localdom.domchgstack_.len() as HighsInt
                    );

                    localdom.get_col_upper_pos(col, domchg.pos, &mut bound_pos);
                }

                if bound_pos != -1 {
                    self.resolved_domain_changes
                        .push(HighsDomain_ConflictSet_LocalDomChg {
                            pos: bound_pos,
                            domchg: localdom.domchgstack_[bound_pos as usize].clone(),
                        });
                }

                return true;
            }
            HighsDomain_Reason::kModelRowLower => {
                let row_index = reason.index;

                // retrieve the matrix values of the cut
                let mut len: HighsInt = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                ms.get_row(row_index, &mut len, &mut inds, &mut vals);

                let max_act = globaldom.get_max_activity(row_index);

                return self.explain_bound_change_geq(
                    ms,
                    localdom,
                    globaldom,
                    current_frontier,
                    &domchg,
                    inds,
                    vals,
                    len,
                    ms.row_lower_2(row_index),
                    max_act,
                );
            }
            HighsDomain_Reason::kModelRowUpper => {
                let row_index = reason.index;

                // retrieve the matrix values of the cut
                let mut len: HighsInt = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                ms.get_row(row_index, &mut len, &mut inds, &mut vals);

                let min_act = globaldom.get_min_activity(row_index);

                return self.explain_bound_change_leq(
                    ms,
                    localdom,
                    globaldom,
                    current_frontier,
                    &domchg,
                    inds,
                    vals,
                    len,
                    ms.row_upper_2(row_index),
                    min_act,
                );
            }
            HighsDomain_Reason::kObjective => {
                // copy the constraint out of `localdom` before it is used again
                let (ov, oi, olen, orhs) = {
                    let mut vals: &[f64] = &[];
                    let mut inds: &[HighsInt] = &[];
                    let mut len: HighsInt = 0;
                    let mut rhs: f64 = 0.0;
                    localdom.obj_prop_get_propagation_constraint(
                        ms,
                        domchg.pos,
                        &mut vals,
                        &mut inds,
                        &mut len,
                        &mut rhs,
                        domchg.domchg.column,
                    );
                    (vals.to_vec(), inds.to_vec(), len, rhs)
                };

                let h16_ninfmin = ms.options_mip_.h16_ninfmin;
                let mut ninfmin: HighsInt = 0;
                let mut min_act = HighsCDouble::from(0.0);
                globaldom.compute_min_activity(0, olen, &oi, &ov, &mut ninfmin, &mut min_act);
                if h16_ninfmin {
                    // HiGHS 1.16 #3275 (h16_ninfmin): HighsDomain.cpp lines 3806-3807 (1.16.0)
                    // todo: treat case with a single infinite contribution
                    if ninfmin > 0 {
                        return false;
                    }
                } else {
                    debug_assert!(ninfmin <= 1);
                    // todo: treat case with a single infinite contribution that propagated a
                    // bound
                    if ninfmin == 1 {
                        return false;
                    }
                }

                return self.explain_bound_change_leq(
                    ms,
                    localdom,
                    globaldom,
                    current_frontier,
                    &domchg,
                    &oi,
                    &ov,
                    olen,
                    orhs,
                    f64::from(min_act),
                );
            }
            _ => {
                debug_assert!(reason.r#type >= 0);
                debug_assert!(
                    reason.r#type
                        < (localdom.cutpoolpropagation.len()
                            + localdom.conflict_pool_propagation.len())
                            as HighsInt
                );
                if reason.r#type < localdom.cutpoolpropagation.len() as HighsInt {
                    let cutpool_index = reason.r#type as usize;
                    let cut_index = reason.index;
                    let pool = localdom.cutpoolpropagation[cutpool_index].cutpool;

                    // retrieve the matrix values of the cut
                    let mut len: HighsInt = 0;
                    let mut inds: &[HighsInt] = &[];
                    let mut vals: &[f64] = &[];
                    ms.mipdata_.cutpools[pool].get_cut(cut_index, &mut len, &mut inds, &mut vals);
                    let min_act = globaldom.get_min_cut_activity_replayed(pool, cut_index);
                    let rhs = ms.mipdata_.cutpools[pool].get_rhs()[cut_index as usize];

                    return self.explain_bound_change_leq(
                        ms,
                        localdom,
                        globaldom,
                        current_frontier,
                        &domchg,
                        inds,
                        vals,
                        len,
                        rhs,
                        min_act,
                    );
                } else {
                    let conflict_pool_index =
                        (reason.r#type - localdom.cutpoolpropagation.len() as HighsInt) as usize;
                    let conflict_index = reason.index;

                    // C++ `break` out of the switch: falls through to `return false`
                    if (localdom.conflict_pool_propagation[conflict_pool_index].conflict_flag_
                        [conflict_index as usize]
                        & 8)
                        != 0
                    {
                        return false;
                    }

                    // retrieve the conflict entries
                    let pool = localdom.conflict_pool_propagation[conflict_pool_index].conflictpool;
                    let (first, second) = ms.mipdata_.conflictpools[pool].get_conflict_ranges()
                        [conflict_index as usize];
                    let conflict = &ms.mipdata_.conflictpools[pool].get_conflict_entry_vector()
                        [first as usize..second as usize];

                    return self
                        .explain_bound_change_conflict(ms, localdom, globaldom, &domchg, conflict);
                }
            }
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeConflict (lines 3602-3659)
    /// HiGHS 1.16 (#3196): HighsDomain.cpp::HighsDomain::ConflictSet::explainBoundChangeConflict (lines 3863-3910)
    /// C++: `bool explainBoundChangeConflict(const LocalDomChg & locdomchg, const HighsDomainChange * conflict, HighsInt len)`
    pub fn explain_bound_change_conflict(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &HighsDomain,
        globaldom: &HighsDomain,
        locdomchg: &HighsDomain_ConflictSet_LocalDomChg,
        conflict: &[HighsDomainChange],
    ) -> bool {
        // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3863-3910 (1.16.0)
        let h16_bound_chain = ms.options_mip_.h16_bound_chain;
        self.resolved_domain_changes.clear();
        let domchg = localdom.flip(ms, &locdomchg.domchg);
        let mut found_domchg = false;
        for entry in conflict.iter() {
            if !found_domchg && entry.column == domchg.column && entry.boundtype == domchg.boundtype
            {
                if entry.boundtype == HighsBoundType::kLower {
                    if entry.boundval <= domchg.boundval {
                        found_domchg = true;
                        continue;
                    }
                } else if entry.boundval >= domchg.boundval {
                    found_domchg = true;
                    continue;
                }
            }

            if globaldom.is_active(entry) {
                continue;
            }

            // C++ `HighsInt pos;` is an out-parameter of getColLowerPos / getColUpperPos.
            let mut pos: HighsInt = 0;
            if entry.boundtype == HighsBoundType::kLower {
                let lb = localdom.get_col_lower_pos(entry.column, locdomchg.pos - 1, &mut pos);

                if pos == -1 || lb < entry.boundval {
                    return false;
                }

                // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3893-3894 (1.16.0)
                if h16_bound_chain {
                    if !localdom.walk_bound_chain(&mut pos, entry.boundval, 1) {
                        return false;
                    }
                } else {
                    while localdom.prevboundval_[pos as usize].0 >= entry.boundval {
                        pos = localdom.prevboundval_[pos as usize].1;
                        // since we checked that the bound value is not active globally and is
                        // active for the local domain at pos
                        // pos should never become -1
                        debug_assert!(pos != -1);
                    }
                }
            } else {
                let ub = localdom.get_col_upper_pos(entry.column, locdomchg.pos - 1, &mut pos);

                if pos == -1 || ub > entry.boundval {
                    return false;
                }

                // HiGHS 1.16 #3196 (h16_bound_chain): HighsDomain.cpp lines 3901-3902 (1.16.0)
                if h16_bound_chain {
                    if !localdom.walk_bound_chain(&mut pos, entry.boundval, -1) {
                        return false;
                    }
                } else {
                    while localdom.prevboundval_[pos as usize].0 <= entry.boundval {
                        pos = localdom.prevboundval_[pos as usize].1;
                        // since we checked that the bound value is not active globally and is
                        // active for the local domain at pos
                        // pos should never become -1
                        debug_assert!(pos != -1);
                    }
                }
            }

            self.resolved_domain_changes
                .push(HighsDomain_ConflictSet_LocalDomChg {
                    pos,
                    domchg: localdom.domchgstack_[pos as usize].clone(),
                });
        }

        found_domchg
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::pushQueue (lines 3661-3668)
    /// C++: `void pushQueue(std::set<LocalDomChg>::iterator domchg)`
    pub fn push_queue(&mut self, domchg: HighsDomain_ConflictSet_LocalDomChg) {
        self.resolve_queue.push(domchg);
        cpp_push_heap_by(&mut self.resolve_queue, |a, b| a.pos < b.pos);
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::popQueue (lines 3670-3680)
    /// C++: `std::set<LocalDomChg>::iterator popQueue()`
    pub fn pop_queue(&mut self) -> HighsDomain_ConflictSet_LocalDomChg {
        debug_assert!(!self.resolve_queue.is_empty());
        cpp_pop_heap_by(&mut self.resolve_queue, |a, b| a.pos < b.pos);
        self.resolve_queue
            .pop()
            .expect("popQueue on an empty resolve queue")
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::clearQueue (lines 3682-3682)
    /// C++: `void clearQueue()`
    pub fn clear_queue(&mut self) {
        self.resolve_queue.clear();
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::queueSize (lines 3684-3686)
    /// C++: `HighsInt queueSize() const`
    pub fn queue_size(&self) -> HighsInt {
        self.resolve_queue.len() as HighsInt
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::resolvable (lines 3688-3700)
    /// C++: `bool resolvable(HighsInt domChgPos) const`
    pub fn resolvable(&self, localdom: &HighsDomain, dom_chg_pos: HighsInt) -> bool {
        debug_assert!(dom_chg_pos >= 0);
        debug_assert!(dom_chg_pos < localdom.domchgreason_.len() as HighsInt);
        !matches!(
            localdom.domchgreason_[dom_chg_pos as usize].r#type,
            HighsDomain_Reason::kBranching | HighsDomain_Reason::kUnknown
        )
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::resolveDepth (lines 3702-3783)
    /// C++: `HighsInt resolveDepth(std::set<LocalDomChg> & frontier, HighsInt depthLevel, HighsInt stopSize, HighsPseudocost & pseudocost, HighsInt minResolve = 0, bool increaseConflictScore = false)`
    /// C++ default arguments: `min_resolve = 0`, `increase_conflict_score = false`
    pub fn resolve_depth(
        &mut self,
        ms: &HighsMipSolver,
        localdom: &mut HighsDomain,
        globaldom: &HighsDomain,
        frontier: &mut BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        depth_level: HighsInt,
        stop_size: HighsInt,
        pseudocost: &mut HighsPseudocost,
        min_resolve: HighsInt,
        increase_conflict_score: bool,
    ) -> HighsInt {
        let mut depth_level = depth_level;
        self.clear_queue();
        let start_pos = HighsDomain_ConflictSet_LocalDomChg {
            pos: if depth_level == 0 {
                0
            } else {
                localdom.branch_pos_[(depth_level - 1) as usize] + 1
            },
            domchg: HighsDomainChange::default(),
        };
        while depth_level < localdom.branch_pos_.len() as HighsInt {
            let branch_pos = localdom.branch_pos_[depth_level as usize];
            if localdom.domchgstack_[branch_pos as usize].boundval
                != localdom.prevboundval_[branch_pos as usize].0
            {
                break;
            }
            // printf("skipping redundant depth\n");
            depth_level += 1;
        }

        let mut empty = true;
        if depth_level == localdom.branch_pos_.len() as HighsInt {
            for it in frontier.range(start_pos.clone()..) {
                empty = false;
                if self.resolvable(localdom, it.pos) {
                    self.push_queue(it.clone());
                }
            }
        } else {
            let end_key = lookup_key(localdom.branch_pos_[depth_level as usize]);
            for it in frontier.range(start_pos.clone()..=end_key) {
                empty = false;
                if self.resolvable(localdom, it.pos) {
                    self.push_queue(it.clone());
                }
            }
        }

        if empty {
            return -1;
        }

        let mut num_resolved: HighsInt = 0;

        while self.queue_size() > stop_size || (self.queue_size() > 0 && num_resolved < min_resolve)
        {
            // the queue holds copies of the set elements: `pos` identifies the element
            let popped = self.pop_queue();
            if !self.explain_bound_change(ms, localdom, globaldom, frontier, popped.clone()) {
                continue;
            }

            num_resolved += 1;
            frontier.remove(&popped);
            let n = self.resolved_domain_changes.len();
            for idx in 0..n {
                let i = self.resolved_domain_changes[idx].clone();
                let inserted = frontier.insert(i.clone());
                if inserted {
                    if increase_conflict_score {
                        if localdom.domchgstack_[i.pos as usize].boundtype == HighsBoundType::kLower
                        {
                            pseudocost.increase_conflict_score_up(
                                localdom.domchgstack_[i.pos as usize].column,
                            );
                        } else {
                            pseudocost.increase_conflict_score_down(
                                localdom.domchgstack_[i.pos as usize].column,
                            );
                        }
                    }
                    if i.pos >= start_pos.pos && self.resolvable(localdom, i.pos) {
                        self.push_queue(i.clone());
                    }
                } else {
                    // C++ `insertResult.first` is the element in the set, and the queue holds
                    // iterators to set elements, so the relaxed bound is seen by a queued
                    // entry too. Here the queue entry for the same `pos` is updated alongside.
                    // if (insertResult.first->domchg.boundval != i.domchg.boundval) printf(..)
                    if let Some(new_bound) = Self::merge_into_frontier(frontier, &i) {
                        for queued in self.resolve_queue.iter_mut() {
                            if queued.pos == i.pos {
                                queued.domchg.boundval = new_bound;
                            }
                        }
                    }
                }
            }
        }

        num_resolved
    }

    /// HiGHS 1.16 (#3443): HighsDomain.cpp::insertIntoFrontier (lines 4016-4029)
    /// C++: `std::pair<std::set<LocalDomChg>::iterator, bool> insertIntoFrontier(std::set<LocalDomChg> & frontier, const LocalDomChg & locdomchg)`
    /// PORT NOTE: returns only C++ `.second` (true when newly inserted). Its one caller (the
    /// flag-on loop in `conflict_analysis`) ignores `.first`. The merge of an existing entry
    /// takes the element out of the set, updates its relaxed bound and inserts it again (same
    /// ordering key, so same set contents as the C++ in-place update).
    pub fn insert_into_frontier(
        frontier: &mut BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        locdomchg: &HighsDomain_ConflictSet_LocalDomChg,
    ) -> bool {
        // the frontier is ordered by stack position only, so if an entry for the
        // same position already exists, keep the stronger of the two relaxed bounds
        if frontier.insert(locdomchg.clone()) {
            return true;
        }
        Self::merge_into_frontier(frontier, locdomchg);
        false
    }

    /// The `!insertResult.second` branch of `insertIntoFrontier` (HighsDomain.cpp lines
    /// 4021-4027): `locdomchg.pos` is already in `frontier`, and its element's relaxed bound
    /// becomes `max` (lower) / `min` (upper) of the two.
    /// PORT NOTE: C++ updates the `mutable` member in place. A `BTreeSet` element cannot be
    /// mutated, so the merged element is written back with `replace` (same `pos`, so the
    /// same slot), and only when the merge changes the bound's bits; otherwise the set is
    /// left as is. Returns the new bound when it changed, so the caller can update the
    /// queued copy of the same `pos`.
    fn merge_into_frontier(
        frontier: &mut BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        locdomchg: &HighsDomain_ConflictSet_LocalDomChg,
    ) -> Option<f64> {
        let existing = frontier
            .get(locdomchg)
            .expect("insert reported an existing element");
        let bound = existing.domchg.boundval;
        let merged = if locdomchg.domchg.boundtype == HighsBoundType::kLower {
            cpp_max(bound, locdomchg.domchg.boundval)
        } else {
            cpp_min(bound, locdomchg.domchg.boundval)
        };
        if merged.to_bits() == bound.to_bits() {
            return None;
        }
        let mut updated = existing.clone();
        updated.domchg.boundval = merged;
        frontier.replace(updated);
        Some(merged)
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::computeCuts (lines 3785-3825)
    /// C++: `HighsInt computeCuts(HighsInt depthLevel, HighsConflictPool & conflictPool, HighsPseudocost & pseudocost)`
    pub fn compute_cuts(
        &mut self,
        ms: &mut HighsMipSolver,
        localdom: &mut HighsDomain,
        globaldom: usize,
        depth_level: HighsInt,
        conflict_pool: usize,
        pseudocost: &mut HighsPseudocost,
    ) -> HighsInt {
        // C++ `resolveDepth(reasonSideFrontier, ...)`: the frontier is taken out of `self`
        // for the call (resolveDepth never reads reasonSideFrontier), then put back.
        let mut reason_side_frontier = std::mem::take(&mut self.reason_side_frontier);
        let num_resolved = self.resolve_depth(
            &*ms,
            localdom,
            &ms.mipdata_.domains[globaldom],
            &mut reason_side_frontier,
            depth_level,
            1,
            pseudocost,
            if depth_level == localdom.branch_pos_.len() as HighsInt {
                1
            } else {
                0
            },
            true,
        );
        self.reason_side_frontier = reason_side_frontier;
        if num_resolved == -1 {
            return -1;
        }
        let mut num_conflicts: HighsInt = 0;
        if num_resolved > 0 {
            // add conflict cut
            ms.mipdata_
                .debug_solution
                .check_conflict_reason_frontier(&self.reason_side_frontier, &localdom.domchgstack_);

            detach!(ms.mipdata_.conflictpools[conflict_pool], |p| p
                .add_conflict_cut(&*ms, &*localdom, &self.reason_side_frontier));
            num_conflicts += 1;
        }

        // if the queue size is 1 then we have a resolvable UIP that is not the
        // branch vertex
        if self.queue_size() == 1 {
            let uip = self.pop_queue();
            self.clear_queue();

            // compute the UIP reconvergence cut
            self.reconvergence_frontier.clear();
            self.reconvergence_frontier.insert(uip.clone());
            let mut reconvergence_frontier = std::mem::take(&mut self.reconvergence_frontier);
            let num_resolved_2 = self.resolve_depth(
                &*ms,
                localdom,
                &ms.mipdata_.domains[globaldom],
                &mut reconvergence_frontier,
                depth_level,
                0,
                pseudocost,
                0,
                false,
            );
            self.reconvergence_frontier = reconvergence_frontier;

            if num_resolved_2 > 0 && !self.reconvergence_frontier.contains(&uip) {
                ms.mipdata_
                    .debug_solution
                    .check_conflict_reconvergence_frontier(
                        &self.reconvergence_frontier,
                        &uip,
                        &localdom.domchgstack_,
                    );
                detach!(ms.mipdata_.conflictpools[conflict_pool], |p| p
                    .add_reconvergence_cut(
                        &*ms,
                        &*localdom,
                        &self.reconvergence_frontier,
                        &uip.domchg,
                    ));
                num_conflicts += 1;
            }
        }

        num_conflicts
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::conflictAnalysis (lines 3827-3890)
    /// HiGHS 1.16 (#3443): HighsDomain.cpp::ConflictSet::conflictAnalysis (lines 4074-4138), flag `h16_frontier_merge`
    /// C++: `void conflictAnalysis(HighsConflictPool & conflictPool, HighsPseudocost & pseudocost)`
    pub fn conflict_analysis(
        &mut self,
        ms: &mut HighsMipSolver,
        localdom: &mut HighsDomain,
        globaldom: usize,
        conflict_pool: usize,
        pseudocost: &mut HighsPseudocost,
    ) {
        let h16_frontier_merge = ms.options_mip_.h16_frontier_merge;
        // reserve() is capacity only
        if !self.explain_infeasibility(&*ms, localdom, &ms.mipdata_.domains[globaldom]) {
            return;
        }

        crate::work::conflict_event(crate::work::ConflictKind::Explained, 1);
        pseudocost.increase_conflict_weight();

        for locdomchg in self.resolved_domain_changes.iter() {
            if locdomchg.domchg.boundtype == HighsBoundType::kLower {
                pseudocost.increase_conflict_score_up(locdomchg.domchg.column);
            } else {
                pseudocost.increase_conflict_score_down(locdomchg.domchg.column);
            }
        }

        if 10 * self.resolved_domain_changes.len() > 1000 + 3 * ms.mipdata_.integral_cols.len() {
            return;
        }

        // HiGHS 1.16 #3443 (h16_frontier_merge): HighsDomain.cpp lines 4094-4095 (1.16.0)
        if h16_frontier_merge {
            for entry in self.resolved_domain_changes.iter() {
                Self::insert_into_frontier(&mut self.reason_side_frontier, entry);
            }
        } else {
            for entry in self.resolved_domain_changes.iter() {
                self.reason_side_frontier.insert(entry.clone());
            }
        }

        ms.mipdata_
            .debug_solution
            .check_conflict_reason_frontier(&self.reason_side_frontier, &localdom.domchgstack_);

        let mut num_conflicts: HighsInt = 0;
        let mut last_depth: HighsInt = localdom.branch_pos_.len() as HighsInt;
        let mut curr_depth: HighsInt = last_depth;
        while curr_depth >= 0 {
            if curr_depth > 0 {
                // skip redundant branching changes which are just added for symmetry
                // handling
                let branchpos = localdom.branch_pos_[(curr_depth - 1) as usize];
                if localdom.domchgstack_[branchpos as usize].boundval
                    == localdom.prevboundval_[branchpos as usize].0
                {
                    last_depth -= 1;
                    curr_depth -= 1;
                    continue;
                }
            }
            let num_new_conflicts = self.compute_cuts(
                ms,
                localdom,
                globaldom,
                curr_depth,
                conflict_pool,
                pseudocost,
            );
            // if the depth level was empty, do not consider it
            if num_new_conflicts == -1 {
                last_depth -= 1;
                curr_depth -= 1;
                continue;
            }

            num_conflicts += num_new_conflicts;
            // if no conflict was found in the first non-empty depth level we stop here
            if num_conflicts == 0 {
                break;
            }
            // if no conflict was found in this depth level and all conflicts of the
            // first 5 non-empty depth levels are generated we stop here
            if last_depth - curr_depth >= 4 && num_new_conflicts == 0 {
                break;
            }
            curr_depth -= 1;
        }

        // if we stopped in the highest non-empty depth it means no conflicts where
        // added yet. We want to at least add the current conflict frontier as it
        // means the bound change leading to infeasibility was the last branching
        // itself and hence should have been propagated in the previous depth but was
        // not, e.g. because the threshold for an integral variable was not reached.
        if curr_depth == last_depth {
            detach!(ms.mipdata_.conflictpools[conflict_pool], |p| p
                .add_conflict_cut(&*ms, &*localdom, &self.reason_side_frontier));
        }
    }

    /// HiGHS: HighsDomain.cpp::HighsDomain::ConflictSet::conflictAnalysis (lines 3892-3966)
    /// C++: `void conflictAnalysis(const HighsInt * proofinds, const double * proofvals, HighsInt prooflen, double proofrhs, HighsConflictPool & conflictPool, HighsPseudocost & pseudocost)`
    pub fn conflict_analysis_2(
        &mut self,
        ms: &mut HighsMipSolver,
        localdom: &mut HighsDomain,
        globaldom: usize,
        proofinds: &[HighsInt],
        proofvals: &[f64],
        prooflen: HighsInt,
        proofrhs: f64,
        conflict_pool: usize,
        pseudocost: &mut HighsPseudocost,
    ) {
        // reserve() is capacity only
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

        if !self.explain_infeasibility_leq(
            &*ms,
            localdom,
            &ms.mipdata_.domains[globaldom],
            proofinds,
            proofvals,
            prooflen,
            proofrhs,
            f64::from(activitymin),
        ) {
            return;
        }

        crate::work::conflict_event(crate::work::ConflictKind::Explained, 1);
        pseudocost.increase_conflict_weight();
        for locdomchg in self.resolved_domain_changes.iter() {
            if locdomchg.domchg.boundtype == HighsBoundType::kLower {
                pseudocost.increase_conflict_score_up(locdomchg.domchg.column);
            } else {
                pseudocost.increase_conflict_score_down(locdomchg.domchg.column);
            }
        }

        if 10 * self.resolved_domain_changes.len() > 1000 + 3 * ms.mipdata_.integral_cols.len() {
            return;
        }

        for entry in self.resolved_domain_changes.iter() {
            self.reason_side_frontier.insert(entry.clone());
        }

        debug_assert!(self.resolved_domain_changes.len() == self.reason_side_frontier.len());

        ms.mipdata_
            .debug_solution
            .check_conflict_reason_frontier(&self.reason_side_frontier, &localdom.domchgstack_);

        let mut num_conflicts: HighsInt = 0;
        let mut last_depth: HighsInt = localdom.branch_pos_.len() as HighsInt;
        let mut curr_depth: HighsInt = last_depth;
        while curr_depth >= 0 {
            if curr_depth > 0 {
                // skip redundant branching changes which are just added for symmetry
                // handling
                let branchpos = localdom.branch_pos_[(curr_depth - 1) as usize];
                if localdom.domchgstack_[branchpos as usize].boundval
                    == localdom.prevboundval_[branchpos as usize].0
                {
                    last_depth -= 1;
                    curr_depth -= 1;
                    continue;
                }
            }
            let num_new_conflicts = self.compute_cuts(
                ms,
                localdom,
                globaldom,
                curr_depth,
                conflict_pool,
                pseudocost,
            );
            // if the depth level was empty, do not consider it
            if num_new_conflicts == -1 {
                last_depth -= 1;
                curr_depth -= 1;
                continue;
            }

            num_conflicts += num_new_conflicts;
            // if no conflict was found in the first non-empty depth level we stop here
            if num_conflicts == 0 {
                break;
            }
            // if no conflict was found in this depth level and all conflicts of the
            // first 5 non-empty depth levels are generated we stop here
            if last_depth - curr_depth >= 4 && num_new_conflicts == 0 {
                break;
            }
            curr_depth -= 1;
        }

        // if we stopped in the highest non-empty depth it means no conflicts where
        // added yet. We want to at least add the current conflict frontier as it
        // means the bound change leading to infeasibility was the last branching
        // itself and hence should have been propagated in the previous depth but was
        // not, e.g. because the threshold for an integral variable was not reached.
        if curr_depth == last_depth {
            detach!(ms.mipdata_.conflictpools[conflict_pool], |p| p
                .add_conflict_cut(&*ms, &*localdom, &self.reason_side_frontier));
        }
    }
}

/// HiGHS: HighsDomain.cpp::computePrio (lines 2764-2768)
/// C++: `static double computePrio(double val, double bound, double globalbound, int64_t numnodes)`
pub fn compute_prio(val: f64, bound: f64, globalbound: f64, numnodes: i64) -> f64 {
    (val * (bound - globalbound) * ((1 + numnodes) as f64)).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One continuous column in [0, 10], feastol 1e-6, epsilon 1e-3, with `h16_cdouble` as given.
    fn h16_cdouble_ms(h16_cdouble: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-3;
        ms.options_mip_.h16_cdouble = h16_cdouble;
        ms
    }

    /// `resolve_linear_geq` with one reason entry: lower bound 0.5 on column 0 (coefficient -10,
    /// base 0), M = 1e-4, Mupper = 0. The 1.15.1 drop test `relax - glb <= epsilon` holds
    /// (relax = 1e-5), but the residual without the entry is 1e-4 > feastol, so the 1.16.0 guard
    /// `M - delta - Mupper <= feastol` rejects the drop.
    fn h16_cdouble_geq_case(h16_cdouble: bool) -> (bool, HighsDomain_ConflictSet) {
        let ms = h16_cdouble_ms(h16_cdouble);
        let mut dom = HighsDomain::new(&ms);
        dom.domchgstack_ = vec![HighsDomainChange {
            boundval: 0.5,
            column: 0,
            boundtype: HighsBoundType::kLower,
        }];
        dom.prevboundval_ = vec![(-kHighsInf, -1)];
        let mut cs = HighsDomain_ConflictSet::new();
        cs.resolve_buffer = vec![HighsDomain_ConflictSet_ResolveCandidate {
            delta: -5.0,
            base_bound: 0.0,
            bound_pos: 0,
            value_pos: 0,
            ..Default::default()
        }];
        let ok = cs.resolve_linear_geq(&ms, &dom, HighsCDouble::from(1e-4), 0.0, &[-10.0]);
        (ok, cs)
    }

    #[test]
    fn h16_cdouble_conflict_on_guard_rejects_drop() {
        let (ok, cs) = h16_cdouble_geq_case(true);
        assert!(ok);
        // 1.16.0: no drop; the bound is relaxed to 1e-5 and the entry stays in the explanation.
        assert_eq!(cs.resolved_domain_changes.len(), 1);
        assert_eq!(cs.resolved_domain_changes[0].pos, 0);
        assert!((cs.resolved_domain_changes[0].domchg.boundval - 1e-5).abs() < 1e-12);
    }

    /// The 1.15.1 drop leaves the conflict uncovered, which is what the C++
    /// `assert(covered <= feastol)` (1.15.1 HighsDomain.cpp:3087) catches in debug builds; the
    /// returned state is only observable with debug assertions off.
    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_cdouble_conflict_off_guard_drops_entry() {
        let (ok, cs) = h16_cdouble_geq_case(false);
        assert!(ok);
        // 1.15.1: the entry is dropped from the explanation.
        assert!(cs.resolved_domain_changes.is_empty());
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "covered <= ms.mipdata_.feastol")]
    fn h16_cdouble_conflict_off_guard_trips_cpp_assert() {
        h16_cdouble_geq_case(false);
    }

    /// Mirror of `h16_cdouble_geq_case` for `resolve_linear_leq`: lower bound 0.5 on column 0
    /// (coefficient 10, base 0, delta 5), M = -1e-4, Mlower = 0. relax = 1e-5 passes the 1.15.1
    /// drop test, but the residual without the entry is -1e-4 < -feastol, so the 1.16.0 guard
    /// `M - delta - Mlower >= -feastol` (HighsDomain.cpp lines 3410-3411) rejects the drop.
    fn h16_cdouble_leq_case(h16_cdouble: bool) -> (bool, HighsDomain_ConflictSet) {
        let ms = h16_cdouble_ms(h16_cdouble);
        let mut dom = HighsDomain::new(&ms);
        dom.domchgstack_ = vec![HighsDomainChange {
            boundval: 0.5,
            column: 0,
            boundtype: HighsBoundType::kLower,
        }];
        dom.prevboundval_ = vec![(-kHighsInf, -1)];
        let mut cs = HighsDomain_ConflictSet::new();
        cs.resolve_buffer = vec![HighsDomain_ConflictSet_ResolveCandidate {
            delta: 5.0,
            base_bound: 0.0,
            bound_pos: 0,
            value_pos: 0,
            ..Default::default()
        }];
        let ok = cs.resolve_linear_leq(&ms, &dom, HighsCDouble::from(-1e-4), 0.0, &[10.0]);
        (ok, cs)
    }

    #[test]
    fn h16_cdouble_conflict_on_leq_guard_rejects_drop() {
        let (ok, cs) = h16_cdouble_leq_case(true);
        assert!(ok);
        // 1.16.0: no drop; the bound is relaxed to 1e-5 and the entry stays in the explanation.
        assert_eq!(cs.resolved_domain_changes.len(), 1);
        assert_eq!(cs.resolved_domain_changes[0].pos, 0);
        assert!((cs.resolved_domain_changes[0].domchg.boundval - 1e-5).abs() < 1e-12);
    }

    /// As for the geq case: the C++ `assert(covered >= -feastol)` (1.15.1
    /// HighsDomain.cpp:3199) fires on the 1.15.1 drop in debug builds.
    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_cdouble_conflict_off_leq_guard_drops_entry() {
        let (ok, cs) = h16_cdouble_leq_case(false);
        assert!(ok);
        // 1.15.1: the entry is dropped from the explanation.
        assert!(cs.resolved_domain_changes.is_empty());
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "covered >= -ms.mipdata_.feastol")]
    fn h16_cdouble_conflict_off_leq_guard_trips_cpp_assert() {
        h16_cdouble_leq_case(false);
    }
}

#[cfg(test)]
mod h16_frontier_tests {
    use super::*;

    fn entry(
        pos: HighsInt,
        boundtype: HighsBoundType,
        boundval: f64,
    ) -> HighsDomain_ConflictSet_LocalDomChg {
        HighsDomain_ConflictSet_LocalDomChg {
            pos,
            domchg: HighsDomainChange {
                boundval,
                column: 2,
                boundtype,
            },
        }
    }

    /// HiGHS 1.16 #3443 (h16_frontier_merge) on path: `insertIntoFrontier` merges an entry with the
    /// same stack position, keeping the stronger relaxed bound (max for lower, min for upper).
    #[test]
    fn h16_frontier_on_merges_same_position() {
        let mut frontier = BTreeSet::new();
        assert!(HighsDomain_ConflictSet::insert_into_frontier(
            &mut frontier,
            &entry(5, HighsBoundType::kLower, 1.0)
        ));
        assert!(!HighsDomain_ConflictSet::insert_into_frontier(
            &mut frontier,
            &entry(5, HighsBoundType::kLower, 3.0)
        ));
        assert_eq!(frontier.len(), 1);
        assert_eq!(frontier.iter().next().unwrap().domchg.boundval, 3.0);
        // a weaker relaxed lower bound does not overwrite the stronger one
        assert!(!HighsDomain_ConflictSet::insert_into_frontier(
            &mut frontier,
            &entry(5, HighsBoundType::kLower, 0.5)
        ));
        assert_eq!(frontier.iter().next().unwrap().domchg.boundval, 3.0);

        assert!(HighsDomain_ConflictSet::insert_into_frontier(
            &mut frontier,
            &entry(7, HighsBoundType::kUpper, 4.0)
        ));
        assert!(!HighsDomain_ConflictSet::insert_into_frontier(
            &mut frontier,
            &entry(7, HighsBoundType::kUpper, 2.5)
        ));
        let upper = frontier.iter().find(|e| e.pos == 7).unwrap();
        assert_eq!(upper.domchg.boundval, 2.5);
        assert_eq!(frontier.len(), 2);
    }

    /// HiGHS 1.16 #3443 (h16_frontier_merge) off path: the 1.15.1 `reasonSideFrontier.insert` keeps the
    /// first entry of a duplicate stack position and does not merge the relaxed bound.
    #[test]
    fn h16_frontier_off_keeps_first_entry() {
        let mut frontier: BTreeSet<HighsDomain_ConflictSet_LocalDomChg> = BTreeSet::new();
        frontier.insert(entry(5, HighsBoundType::kLower, 1.0));
        frontier.insert(entry(5, HighsBoundType::kLower, 3.0));
        assert_eq!(frontier.len(), 1);
        assert_eq!(frontier.iter().next().unwrap().domchg.boundval, 1.0);
    }
}

#[cfg(test)]
mod h16_ninfmin_tests {
    use super::*;
    use crate::mip::objective_function::HighsObjectiveFunction;

    /// Three columns: x0, x1 free with cost 1; x2 in [0, 10] with cost 1. No rows. The
    /// objective propagation constraint is `x0 + x1 + x2 <= upper_limit` (upper_limit 5), so
    /// the min activity over the constraint has two infinite contributions (x0, x1).
    fn objective_mip() -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 3;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![1.0, 1.0, 1.0];
        lp.col_lower_ = vec![-kHighsInf, -kHighsInf, 0.0];
        lp.col_upper_ = vec![kHighsInf, kHighsInf, 10.0];
        lp.row_lower_ = vec![];
        lp.row_upper_ = vec![];
        lp.integrality_ = vec![HighsVarType::kContinuous; 3];
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0, 0, 0, 0];
        lp.a_matrix_.index_ = vec![];
        lp.a_matrix_.value_ = vec![];
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.upper_limit = 5.0;
        let objective_function = HighsObjectiveFunction::new(&ms);
        ms.mipdata_.objective_function = Slot::new(objective_function);
        ms
    }

    /// Runs `explain_bound_change` on the objective-reason bound change `x2 <= 5`.
    fn explain_objective_change(on: bool) -> bool {
        let mut ms = objective_mip();
        ms.options_mip_.h16_ninfmin = on;
        let mut localdom = HighsDomain::new(&ms);
        let globaldom = HighsDomain::new(&ms);
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 2,
                boundtype: HighsBoundType::kUpper,
            },
            HighsDomain_Reason::objective(),
        );
        let pos = localdom.domchgstack_.len() - 1;
        let domchg = HighsDomain_ConflictSet_LocalDomChg {
            pos: pos as HighsInt,
            domchg: localdom.domchgstack_[pos].clone(),
        };
        let mut conflict_set = HighsDomain_ConflictSet::new();
        let frontier = BTreeSet::new();
        conflict_set.explain_bound_change(&ms, &mut localdom, &globaldom, &frontier, domchg)
    }

    #[test]
    // HiGHS 1.16 #3275: HighsDomain.cpp lines 3806-3807 (1.16.0)
    fn h16_ninfmin_on_two_infinite_contributions_no_explanation() {
        // 1.16: `if (ninfmin > 0) return false;` with ninfmin == 2 -> no explanation.
        assert!(!explain_objective_change(true));
    }

    #[test]
    // The 1.15.1 off path asserts `ninfmin <= 1` in debug builds (C++ `assert`); this
    // ninfmin == 2 case is only reachable with debug assertions disabled (release).
    #[cfg(not(debug_assertions))]
    // HiGHS 1.16 #3275: 1.15.1 HighsDomain.cpp lines 3546-3547 (off path, unchanged).
    fn h16_ninfmin_off_two_infinite_contributions_unchanged() {
        // Off: ninfmin == 2 is not caught by `ninfmin == 1`, so the 1.15.1 explanation runs
        // and succeeds: with upper_limit 5 the bound x2 <= 5 follows from an empty set of
        // other contributions (the infinite ones are skipped by minAct).
        assert!(explain_objective_change(false));
    }
}

#[cfg(test)]
mod h16_bound_chain_tests {
    use super::*;

    /// One continuous column in [0, 10], feastol 1e-6, epsilon 1e-3, with `h16_bound_chain` as given.
    fn h16_bound_chain_ms(h16_bound_chain: bool) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.integrality_ = vec![HighsVarType::kContinuous];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-3;
        ms.options_mip_.h16_bound_chain = h16_bound_chain;
        ms
    }

    /// Lower-bound chain on column 0: local changes at stack positions 0, 1, 2 set the bound to
    /// 0.2, 0.4, 0.6. The previous-bound links are (-inf, -1), (0.2, 0), (0.4, 1).
    fn chain_domain(ms: &HighsMipSolver) -> HighsDomain {
        let mut dom = HighsDomain::new(ms);
        dom.domchgstack_ = [0.2, 0.4, 0.6]
            .iter()
            .map(|&v| HighsDomainChange {
                boundval: v,
                column: 0,
                boundtype: HighsBoundType::kLower,
            })
            .collect();
        dom.prevboundval_ = vec![(-kHighsInf, -1), (0.2, 0), (0.4, 1)];
        dom
    }

    /// `resolve_linear_geq` on the top of the chain (stack position 2, bound 0.6, coefficient -10,
    /// delta -6, base 0), M = 3 before the entry and Mupper = 0. The relaxed bound is
    /// (0 - 3) / -10 = 0.3, which lies in (0.2, 0.4], so the chain walk must stop at stack
    /// position 1. The 1.15.1 drop test (`relax - base <= epsilon`) does not hold, so the walk runs.
    fn run_chain_geq(ms: &HighsMipSolver, dom: &HighsDomain) -> (bool, HighsDomain_ConflictSet) {
        let mut cs = HighsDomain_ConflictSet::new();
        cs.resolve_buffer = vec![HighsDomain_ConflictSet_ResolveCandidate {
            delta: -6.0,
            base_bound: 0.0,
            bound_pos: 2,
            value_pos: 0,
            ..Default::default()
        }];
        let ok = cs.resolve_linear_geq(ms, dom, HighsCDouble::from(3.0), 0.0, &[-10.0]);
        (ok, cs)
    }

    /// Chain of one local change (stack position 0, bound 0.5, previous bound 0.3). The walk
    /// for a relaxed bound of 0.2 passes the previous bound 0.3 and exhausts the chain (pos -1).
    fn exhausted_domain(ms: &HighsMipSolver) -> HighsDomain {
        let mut dom = HighsDomain::new(ms);
        dom.domchgstack_ = vec![HighsDomainChange {
            boundval: 0.5,
            column: 0,
            boundtype: HighsBoundType::kLower,
        }];
        dom.prevboundval_ = vec![(0.3, -1)];
        dom
    }

    #[test]
    fn h16_bound_chain_on_walks_past_relaxed_bound() {
        let ms = h16_bound_chain_ms(true);
        let dom = chain_domain(&ms);
        let (ok, cs) = run_chain_geq(&ms, &dom);
        assert!(ok);
        // 1.16.0 walkBoundChain: from stack position 2 (prev 0.4 >= 0.3) to 1 (prev 0.2 < 0.3).
        assert_eq!(cs.resolved_domain_changes.len(), 1);
        assert_eq!(cs.resolved_domain_changes[0].pos, 1);
        assert!((cs.resolved_domain_changes[0].domchg.boundval - 0.3).abs() < 1e-12);
    }

    #[test]
    fn h16_bound_chain_off_walks_past_relaxed_bound() {
        let ms = h16_bound_chain_ms(false);
        let dom = chain_domain(&ms);
        let (ok, cs) = run_chain_geq(&ms, &dom);
        assert!(ok);
        // 1.15.1 while loop: same chain, same stop position 1.
        assert_eq!(cs.resolved_domain_changes.len(), 1);
        assert_eq!(cs.resolved_domain_changes[0].pos, 1);
        assert!((cs.resolved_domain_changes[0].domchg.boundval - 0.3).abs() < 1e-12);
    }

    #[test]
    fn h16_bound_chain_on_exhausted_chain_rejects() {
        // 1.16.0: walkBoundChain returns false when pos becomes -1, and resolveLinearGeq returns
        // false. (1.15.1 would index prevboundval_[-1] here; that path is not reachable with the
        // flag off, so there is no off-path twin for this case.)
        let ms = h16_bound_chain_ms(true);
        let dom = exhausted_domain(&ms);
        let mut cs = HighsDomain_ConflictSet::new();
        cs.resolve_buffer = vec![HighsDomain_ConflictSet_ResolveCandidate {
            delta: -5.0,
            base_bound: 0.0,
            bound_pos: 0,
            value_pos: 0,
            ..Default::default()
        }];
        let ok = cs.resolve_linear_geq(&ms, &dom, HighsCDouble::from(2.0), 0.0, &[-10.0]);
        assert!(!ok);
    }

    #[test]
    fn h16_bound_chain_on_walk_helper_directions() {
        let ms = h16_bound_chain_ms(true);
        let mut dom = HighsDomain::new(&ms);
        // Lower direction (+1): stop at the first previous bound below the target.
        dom.prevboundval_ = vec![(-kHighsInf, -1), (0.2, 0), (0.4, 1)];
        let mut pos: HighsInt = 2;
        assert!(dom.walk_bound_chain(&mut pos, 0.3, 1));
        assert_eq!(pos, 1);
        // Upper direction (-1): walk while the previous bound is <= target.
        dom.prevboundval_ = vec![(10.0, -1), (5.0, 0), (3.0, 1)];
        let mut pos: HighsInt = 2;
        assert!(dom.walk_bound_chain(&mut pos, 4.0, -1));
        assert_eq!(pos, 1);
        // Exhaustion: every previous bound is <= target, so pos becomes -1 and the walk fails.
        let mut pos: HighsInt = 2;
        assert!(!dom.walk_bound_chain(&mut pos, 20.0, -1));
        assert_eq!(pos, -1);
    }
}

#[cfg(test)]
mod h16_big_act_tests {
    use super::*;

    /// Two continuous columns in [0, 10], no rows, feastol 1e-6, epsilon 1e-9.
    fn big_act_mip() -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 2;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.row_lower_ = vec![];
        lp.row_upper_ = vec![];
        lp.integrality_ = vec![HighsVarType::kContinuous; 2];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0, 0, 0];
        lp.a_matrix_.index_ = vec![];
        lp.a_matrix_.value_ = vec![];
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.nodequeue.set_num_col(2);
        ms
    }

    /// Upper bound 5 on column 1 (pos 0), then upper bound 5 on column 0 (pos 1, the change
    /// being explained). Both columns sit at 5 with global upper 10.
    fn explain_geq_case(on: bool, max_act: f64) -> (bool, HighsDomain_ConflictSet) {
        let mut ms = big_act_mip();
        ms.options_mip_.h16_big_act = on;
        let mut localdom = HighsDomain::new(&ms);
        let globaldom = HighsDomain::new(&ms);
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 1,
                boundtype: HighsBoundType::kUpper,
            },
            HighsDomain_Reason::unspecified(),
        );
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 0,
                boundtype: HighsBoundType::kUpper,
            },
            HighsDomain_Reason::unspecified(),
        );
        let pos = localdom.domchgstack_.len() - 1;
        let domchg = HighsDomain_ConflictSet_LocalDomChg {
            pos: pos as HighsInt,
            domchg: localdom.domchgstack_[pos].clone(),
        };
        let mut cs = HighsDomain_ConflictSet::new();
        let frontier = BTreeSet::new();
        let ok = cs.explain_bound_change_geq(
            &ms,
            &localdom,
            &globaldom,
            &frontier,
            &domchg,
            &[0, 1],
            &[1.0, 1.0],
            2,
            0.0,
            max_act,
        );
        (ok, cs)
    }

    /// Mirror of `explain_geq_case`: lower bound 5 on both columns, coefficients +1 (the min
    /// activity reads lower bounds for positive coefficients), `min_act` as given.
    fn explain_leq_case(on: bool, min_act: f64) -> (bool, HighsDomain_ConflictSet) {
        let mut ms = big_act_mip();
        ms.options_mip_.h16_big_act = on;
        let mut localdom = HighsDomain::new(&ms);
        let globaldom = HighsDomain::new(&ms);
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 1,
                boundtype: HighsBoundType::kLower,
            },
            HighsDomain_Reason::unspecified(),
        );
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 0,
                boundtype: HighsBoundType::kLower,
            },
            HighsDomain_Reason::unspecified(),
        );
        let pos = localdom.domchgstack_.len() - 1;
        let domchg = HighsDomain_ConflictSet_LocalDomChg {
            pos: pos as HighsInt,
            domchg: localdom.domchgstack_[pos].clone(),
        };
        let mut cs = HighsDomain_ConflictSet::new();
        let frontier = BTreeSet::new();
        let ok = cs.explain_bound_change_leq(
            &ms,
            &localdom,
            &globaldom,
            &frontier,
            &domchg,
            &[0, 1],
            &[1.0, 1.0],
            2,
            0.0,
            min_act,
        );
        (ok, cs)
    }

    /// Column 0 upper bound 5 (pos 0); the infeasibility explanation sees it as a candidate.
    fn explain_infeasibility_geq_case(on: bool, max_act: f64) -> (bool, HighsDomain_ConflictSet) {
        let mut ms = big_act_mip();
        ms.options_mip_.h16_big_act = on;
        let mut localdom = HighsDomain::new(&ms);
        let globaldom = HighsDomain::new(&ms);
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 0,
                boundtype: HighsBoundType::kUpper,
            },
            HighsDomain_Reason::unspecified(),
        );
        let mut cs = HighsDomain_ConflictSet::new();
        let ok =
            cs.explain_infeasibility_geq(&ms, &localdom, &globaldom, &[0], &[1.0], 1, 0.0, max_act);
        (ok, cs)
    }

    /// Column 0 lower bound 5 (pos 0); the infeasibility explanation sees it as a candidate.
    fn explain_infeasibility_leq_case(on: bool, min_act: f64) -> (bool, HighsDomain_ConflictSet) {
        let mut ms = big_act_mip();
        ms.options_mip_.h16_big_act = on;
        let mut localdom = HighsDomain::new(&ms);
        let globaldom = HighsDomain::new(&ms);
        localdom.change_bound(
            &mut ms,
            HighsDomainChange {
                boundval: 5.0,
                column: 0,
                boundtype: HighsBoundType::kLower,
            },
            HighsDomain_Reason::unspecified(),
        );
        let mut cs = HighsDomain_ConflictSet::new();
        let ok =
            cs.explain_infeasibility_leq(&ms, &localdom, &globaldom, &[0], &[1.0], 1, 0.0, min_act);
        (ok, cs)
    }

    // HiGHS 1.16 #3267: HighsDomain.cpp line 3026 (1.16.0): |maxAct| > 2^52 bails out before any work.
    #[test]
    fn h16_big_act_on_geq_bails_out_on_huge_activity() {
        let (ok, cs) = explain_geq_case(true, 9007199254740992.0);
        assert!(!ok);
        assert!(cs.resolve_buffer.is_empty());
    }

    #[test]
    fn h16_big_act_off_geq_runs_1_15_1_path_on_huge_activity() {
        // 1.15.1 only bails out on +inf, so the candidate for column 1 is collected before the
        // resolution fails on the huge residual.
        let (ok, cs) = explain_geq_case(false, 9007199254740992.0);
        assert!(!ok);
        assert_eq!(cs.resolve_buffer.len(), 1);
    }

    #[test]
    fn h16_big_act_on_leq_bails_out_on_huge_activity() {
        let (ok, cs) = explain_leq_case(true, -(9007199254740992.0));
        assert!(!ok);
        assert!(cs.resolve_buffer.is_empty());
    }

    #[test]
    fn h16_big_act_off_leq_runs_1_15_1_path_on_huge_activity() {
        let (ok, cs) = explain_leq_case(false, -(9007199254740992.0));
        assert!(!ok);
        assert_eq!(cs.resolve_buffer.len(), 1);
    }

    #[test]
    fn h16_big_act_on_infeasibility_geq_bails_out_on_huge_activity() {
        let (ok, cs) = explain_infeasibility_geq_case(true, 9007199254740992.0);
        assert!(!ok);
        assert!(cs.resolve_buffer.is_empty());
    }

    #[test]
    fn h16_big_act_off_infeasibility_geq_runs_1_15_1_path_on_huge_activity() {
        let (ok, cs) = explain_infeasibility_geq_case(false, 9007199254740992.0);
        assert!(!ok);
        assert_eq!(cs.resolve_buffer.len(), 1);
    }

    #[test]
    fn h16_big_act_on_infeasibility_leq_bails_out_on_huge_activity() {
        let (ok, cs) = explain_infeasibility_leq_case(true, -(9007199254740992.0));
        assert!(!ok);
        assert!(cs.resolve_buffer.is_empty());
    }

    #[test]
    fn h16_big_act_off_infeasibility_leq_runs_1_15_1_path_on_huge_activity() {
        let (ok, cs) = explain_infeasibility_leq_case(false, -(9007199254740992.0));
        assert!(!ok);
        assert_eq!(cs.resolve_buffer.len(), 1);
    }
}
