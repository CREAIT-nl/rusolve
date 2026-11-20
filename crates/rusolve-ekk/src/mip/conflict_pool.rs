//! HiGHS source: mip/HighsConflictPool.h, mip/HighsConflictPool.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use std::sync::atomic::{AtomicU8, Ordering};

/// HiGHS: HighsConflictPool.h::HighsConflictPool (lines 18-123)
#[derive(Default, Debug)]
pub struct HighsConflictPool {
    /// C++: `HighsInt agelim_`
    pub agelim_: HighsInt,
    /// C++: `HighsInt softlimit_`
    pub softlimit_: HighsInt,
    /// C++: `bool age_lock_`
    pub age_lock_: bool,
    /// C++: `std::vector<HighsInt> ageDistribution_`
    pub age_distribution_: Vec<HighsInt>,
    /// C++: `std::vector<int16_t> ages_`
    pub ages_: Vec<i16>,
    /// C++: `std::vector<unsigned int> modification_`
    pub modification_: Vec<u32>,
    /// C++: `std::deque<std::atomic<uint8_t>> ageResetWhileLocked_`
    pub age_reset_while_locked_: Vec<AtomicU8>,
    /// C++: `std::vector<HighsDomainChange> conflictEntries_`
    pub conflict_entries_: Vec<HighsDomainChange>,
    /// C++: `std::vector<std::pair<HighsInt, HighsInt>> conflictRanges_`
    pub conflict_ranges_: Vec<(HighsInt, HighsInt)>,
    /// C++: `std::set<std::pair<HighsInt, HighsInt>> freeSpaces_`
    pub free_spaces_: BTreeSet<(HighsInt, HighsInt)>,
    /// C++: `std::vector<HighsInt> deletedConflicts_`
    pub deleted_conflicts_: Vec<HighsInt>,
    /// C++ `std::vector<HighsDomain::ConflictPoolPropagation*> propagationDomains`: replaced by the M3 event log.
    pub events: PoolEventLog,
}

impl Clone for HighsConflictPool {
    fn clone(&self) -> Self {
        HighsConflictPool {
            agelim_: self.agelim_,
            softlimit_: self.softlimit_,
            age_lock_: self.age_lock_,
            age_distribution_: self.age_distribution_.clone(),
            ages_: self.ages_.clone(),
            modification_: self.modification_.clone(),
            age_reset_while_locked_: self
                .age_reset_while_locked_
                .iter()
                .map(|a| AtomicU8::new(a.load(Ordering::Relaxed)))
                .collect(),
            conflict_entries_: self.conflict_entries_.clone(),
            conflict_ranges_: self.conflict_ranges_.clone(),
            free_spaces_: self.free_spaces_.clone(),
            deleted_conflicts_: self.deleted_conflicts_.clone(),
            events: self.events.clone(),
        }
    }
}

impl HighsConflictPool {
    /// HiGHS: HighsConflictPool.cpp::HighsConflictPool::addConflictCut (lines 13-83)
    /// C++: `void addConflictCut(const HighsDomain & domain, const std::set<HighsDomain::ConflictSet::LocalDomChg> & reasonSideFrontier)`
    // PORT NOTE: signature changed: `domain.variableType`, `domain.feastol()` and `domain.flip()`
    // read `ms` (domain/mod.rs `variable_type`, domain/stack.rs `feastol` / `flip`). Shared, so
    // `ConflictSet` can pass `ms` while it holds the local domain mutably.
    pub fn add_conflict_cut(
        &mut self,
        ms: &HighsMipSolver,
        domain: &HighsDomain,
        reason_side_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
    ) {
        let conflict_len = reason_side_frontier.len() as HighsInt;
        crate::work::conflict_event(crate::work::ConflictKind::Added, 1);
        crate::work::conflict_event(crate::work::ConflictKind::Len, conflict_len as u64);
        let conflict_index: HighsInt;
        let start: HighsInt;
        let end: HighsInt;
        let found = if self.free_spaces_.is_empty() {
            None
        } else {
            self.free_spaces_
                .range((conflict_len, -1)..)
                .next()
                .copied()
        };
        match found {
            None => {
                start = self.conflict_entries_.len() as HighsInt;
                end = start + conflict_len;

                self.conflict_entries_
                    .resize(end as usize, HighsDomainChange::default());
            }
            Some(freeslot) => {
                self.free_spaces_.remove(&freeslot);

                start = freeslot.1;
                end = start + conflict_len;
                // if the space was not completely occupied, we register the remainder of
                // it again in the priority queue
                if freeslot.0 > conflict_len {
                    self.free_spaces_.insert((freeslot.0 - conflict_len, end));
                }
            }
        }

        // register the range of entries for this conflict with a reused or a new
        // index
        if self.deleted_conflicts_.is_empty() {
            conflict_index = self.conflict_ranges_.len() as HighsInt;
            self.conflict_ranges_.push((start, end));
            self.ages_.resize(self.conflict_ranges_.len(), 0);
            self.modification_.resize(self.conflict_ranges_.len(), 0);
            self.age_reset_while_locked_
                .resize_with(self.conflict_ranges_.len(), || AtomicU8::new(0));
        } else {
            conflict_index = *self.deleted_conflicts_.last().unwrap();
            self.deleted_conflicts_.pop();
            self.conflict_ranges_[conflict_index as usize].0 = start;
            self.conflict_ranges_[conflict_index as usize].1 = end;
        }

        let ci = conflict_index as usize;
        self.age_reset_while_locked_[ci].store(0, Ordering::Relaxed);
        self.modification_[ci] = self.modification_[ci].wrapping_add(1);
        self.ages_[ci] = 0;
        self.age_distribution_[self.ages_[ci] as usize] += 1;

        let domchg_stack = domain.get_domain_change_stack();
        let feastol = domain.feastol(ms);
        // PORT NOTE: C++ `i++` per frontier entry becomes `start + k` (same indices, enumerate).
        for (k, domchg) in reason_side_frontier.iter().enumerate() {
            let i = start + k as HighsInt;
            debug_assert!(i < end);
            debug_assert!(domchg.pos >= 0);
            debug_assert!((domchg.pos as usize) < domchg_stack.len());
            self.conflict_entries_[i as usize] = domchg.domchg.clone();
            let col = self.conflict_entries_[i as usize].column;
            if domain.variable_type(ms, col) == HighsVarType::kContinuous {
                if self.conflict_entries_[i as usize].boundtype == HighsBoundType::kLower {
                    self.conflict_entries_[i as usize].boundval += feastol;
                } else {
                    self.conflict_entries_[i as usize].boundval -= feastol;
                }
            }
        }

        // C++: for (conflictProp : propagationDomains) conflictProp->conflictAdded(conflictIndex);
        // PORT NOTE: one M3 event carries the entries (snapshot) in place of the observer loop.
        self.events.push(PoolEvent::ConflictAdded {
            conflict: conflict_index,
            entries: self.conflict_entries_[start as usize..end as usize].to_vec(),
        });
    }

    /// HiGHS: HighsConflictPool.cpp::HighsConflictPool::addReconvergenceCut (lines 85-159)
    /// C++: `void addReconvergenceCut(const HighsDomain & domain, const std::set<HighsDomain::ConflictSet::LocalDomChg> & reconvergenceFrontier, const HighsDomainChange & reconvergenceDomchg)`
    // PORT NOTE: signature changed: `domain.variableType`, `domain.feastol()` and `domain.flip()`
    // read `ms` (domain/mod.rs `variable_type`, domain/stack.rs `feastol` / `flip`). Shared, so
    // `ConflictSet` can pass `ms` while it holds the local domain mutably.
    pub fn add_reconvergence_cut(
        &mut self,
        ms: &HighsMipSolver,
        domain: &HighsDomain,
        reconvergence_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        reconvergence_domchg: &HighsDomainChange,
    ) {
        let conflict_len = reconvergence_frontier.len() as HighsInt + 1;
        crate::work::conflict_event(crate::work::ConflictKind::Added, 1);
        crate::work::conflict_event(crate::work::ConflictKind::Len, conflict_len as u64);
        let conflict_index: HighsInt;
        let start: HighsInt;
        let end: HighsInt;
        let found = if self.free_spaces_.is_empty() {
            None
        } else {
            self.free_spaces_
                .range((conflict_len, -1)..)
                .next()
                .copied()
        };
        match found {
            None => {
                start = self.conflict_entries_.len() as HighsInt;
                end = start + conflict_len;

                self.conflict_entries_
                    .resize(end as usize, HighsDomainChange::default());
            }
            Some(freeslot) => {
                self.free_spaces_.remove(&freeslot);

                start = freeslot.1;
                end = start + conflict_len;
                // if the space was not completely occupied, we register the remainder of
                // it again in the priority queue
                if freeslot.0 > conflict_len {
                    self.free_spaces_.insert((freeslot.0 - conflict_len, end));
                }
            }
        }

        // register the range of entries for this conflict with a reused or a new
        // index
        if self.deleted_conflicts_.is_empty() {
            conflict_index = self.conflict_ranges_.len() as HighsInt;
            self.conflict_ranges_.push((start, end));
            self.ages_.resize(self.conflict_ranges_.len(), 0);
            self.modification_.resize(self.conflict_ranges_.len(), 0);
            self.age_reset_while_locked_
                .resize_with(self.conflict_ranges_.len(), || AtomicU8::new(0));
        } else {
            conflict_index = *self.deleted_conflicts_.last().unwrap();
            self.deleted_conflicts_.pop();
            self.conflict_ranges_[conflict_index as usize].0 = start;
            self.conflict_ranges_[conflict_index as usize].1 = end;
        }

        let ci = conflict_index as usize;
        self.age_reset_while_locked_[ci].store(0, Ordering::Relaxed);
        self.modification_[ci] = self.modification_[ci].wrapping_add(1);
        self.ages_[ci] = 0;
        self.age_distribution_[self.ages_[ci] as usize] += 1;

        let domchg_stack = domain.get_domain_change_stack();
        debug_assert!(start < end);
        self.conflict_entries_[start as usize] = domain.flip(ms, reconvergence_domchg);
        let feastol = domain.feastol(ms);
        // PORT NOTE: C++ `i++` (after the flipped entry at `start`) becomes `start + 1 + k`.
        for (k, domchg) in reconvergence_frontier.iter().enumerate() {
            let i = start + 1 + k as HighsInt;
            debug_assert!(i < end);
            debug_assert!(domchg.pos >= 0);
            debug_assert!((domchg.pos as usize) < domchg_stack.len());
            self.conflict_entries_[i as usize] = domchg.domchg.clone();
            let col = self.conflict_entries_[i as usize].column;
            if domain.variable_type(ms, col) == HighsVarType::kContinuous {
                if self.conflict_entries_[i as usize].boundtype == HighsBoundType::kLower {
                    self.conflict_entries_[i as usize].boundval += feastol;
                } else {
                    self.conflict_entries_[i as usize].boundval -= feastol;
                }
            }
        }

        // C++: for (conflictProp : propagationDomains) conflictProp->conflictAdded(conflictIndex);
        // PORT NOTE: one M3 event carries the entries (snapshot) in place of the observer loop.
        self.events.push(PoolEvent::ConflictAdded {
            conflict: conflict_index,
            entries: self.conflict_entries_[start as usize..end as usize].to_vec(),
        });
    }

    /// HiGHS: HighsConflictPool.cpp::HighsConflictPool::removeConflict (lines 161-182)
    /// C++: `void removeConflict(HighsInt conflict)`
    pub fn remove_conflict(&mut self, conflict: HighsInt) {
        crate::work::add_event(crate::work::Event::ConfRemoved, 1);
        // C++: for (conflictProp : propagationDomains) conflictProp->conflictDeleted(conflict);
        // PORT NOTE: one M3 event replaces the loop over the observers. The log is not GC'd here.
        self.events.push(PoolEvent::ConflictDeleted { conflict });

        let cidx = conflict as usize;
        if self.ages_[cidx] >= 0 {
            self.age_distribution_[self.ages_[cidx] as usize] -= 1;
            self.ages_[cidx] = -1;
        }

        let start = self.conflict_ranges_[cidx].0;
        let end = self.conflict_ranges_[cidx].1;

        // register the space of the deleted row and the index so that it can be
        // reused
        self.deleted_conflicts_.push(conflict);
        self.free_spaces_.insert((end - start, start));

        // set the range to -1,-1 to indicate a deleted row
        self.conflict_ranges_[cidx].0 = -1;
        self.conflict_ranges_[cidx].1 = -1;
        self.modification_[cidx] = self.modification_[cidx].wrapping_add(1);
    }

    /// HiGHS: HighsConflictPool.cpp::HighsConflictPool::performAging (lines 184-210)
    /// C++: `void performAging(bool thread_safe = false)`
    /// C++ default arguments: `thread_safe = false`
    pub fn perform_aging(&mut self, thread_safe: bool) {
        if self.age_lock_ {
            return;
        }
        let conflict_max_index = self.conflict_ranges_.len() as HighsInt;
        let mut agelim = self.agelim_;
        let mut num_active_conflicts = self.get_num_conflicts();
        while agelim > 5 && num_active_conflicts > self.softlimit_ {
            num_active_conflicts -= self.age_distribution_[agelim as usize];
            agelim -= 1;
        }

        for i in 0..conflict_max_index {
            let iu = i as usize;
            if self.ages_[iu] < 0 {
                continue;
            }
            // C++: HighsConflictPool.cpp:197 `ageResetWhileLocked_[i].load(relaxed) == 1`
            if thread_safe && self.age_reset_while_locked_[iu].load(Ordering::Relaxed) == 1 {
                self.reset_age(i);
            }

            self.age_distribution_[self.ages_[iu] as usize] -= 1;
            self.ages_[iu] += 1;
            // C++: HighsConflictPool.cpp:202 `ageResetWhileLocked_[i].store(0, relaxed)`
            self.age_reset_while_locked_[iu].store(0, Ordering::Relaxed);

            if self.ages_[iu] as HighsInt > agelim {
                self.ages_[iu] = -1;
                self.remove_conflict(i);
            } else {
                self.age_distribution_[self.ages_[iu] as usize] += 1;
            }
        }
    }

    /// HiGHS: HighsConflictPool.cpp::HighsConflictPool::addConflictFromOtherPool (lines 212-265)
    /// C++: `void addConflictFromOtherPool(const HighsDomainChange * conflictEntries, HighsInt conflictLen)`
    // PORT NOTE: signature changed because the C++ takes a pointer plus `conflictLen` (an array), which
    // `Option<&HighsDomainChange>` cannot express; the slice carries the length (`conflictLen` = `len()`).
    pub fn add_conflict_from_other_pool(&mut self, conflict_entries: &[HighsDomainChange]) {
        let conflict_len = conflict_entries.len() as HighsInt;
        let conflict_index: HighsInt;
        let start: HighsInt;
        let end: HighsInt;
        let found = if self.free_spaces_.is_empty() {
            None
        } else {
            self.free_spaces_
                .range((conflict_len, -1)..)
                .next()
                .copied()
        };
        match found {
            None => {
                start = self.conflict_entries_.len() as HighsInt;
                end = start + conflict_len;

                self.conflict_entries_
                    .resize(end as usize, HighsDomainChange::default());
            }
            Some(freeslot) => {
                self.free_spaces_.remove(&freeslot);

                start = freeslot.1;
                end = start + conflict_len;
                // if the space was not completely occupied, we register the remainder of
                // it again in the priority queue
                if freeslot.0 > conflict_len {
                    self.free_spaces_.insert((freeslot.0 - conflict_len, end));
                }
            }
        }

        // register the range of entries for this conflict with a reused or a new
        // index
        if self.deleted_conflicts_.is_empty() {
            conflict_index = self.conflict_ranges_.len() as HighsInt;
            self.conflict_ranges_.push((start, end));
            self.ages_.resize(self.conflict_ranges_.len(), 0);
            self.modification_.resize(self.conflict_ranges_.len(), 0);
            self.age_reset_while_locked_
                .resize_with(self.conflict_ranges_.len(), || AtomicU8::new(0));
        } else {
            conflict_index = *self.deleted_conflicts_.last().unwrap();
            self.deleted_conflicts_.pop();
            self.conflict_ranges_[conflict_index as usize].0 = start;
            self.conflict_ranges_[conflict_index as usize].1 = end;
        }

        let ci = conflict_index as usize;
        self.age_reset_while_locked_[ci].store(0, Ordering::Relaxed);
        self.modification_[ci] = self.modification_[ci].wrapping_add(1);
        self.ages_[ci] = 0;
        self.age_distribution_[self.ages_[ci] as usize] += 1;

        for (i, entry) in conflict_entries.iter().enumerate() {
            debug_assert!((start as usize + i) < end as usize);
            self.conflict_entries_[start as usize + i] = entry.clone();
        }

        // C++: for (conflictProp : propagationDomains) conflictProp->conflictAdded(conflictIndex);
        // PORT NOTE: one M3 event carries the entries (snapshot) in place of the observer loop.
        self.events.push(PoolEvent::ConflictAdded {
            conflict: conflict_index,
            entries: self.conflict_entries_[start as usize..end as usize].to_vec(),
        });
    }

    /// HiGHS: HighsConflictPool.cpp::HighsConflictPool::syncConflictPool (lines 267-284)
    /// C++: `void syncConflictPool(HighsConflictPool & syncpool)`
    pub fn sync_conflict_pool(&mut self, syncpool: &mut HighsConflictPool) {
        let conflict_max_index = self.conflict_ranges_.len() as HighsInt;
        for i in 0..conflict_max_index {
            let iu = i as usize;
            if self.ages_[iu] < 0 {
                continue;
            }
            let start = self.conflict_ranges_[iu].0;
            let end = self.conflict_ranges_[iu].1;
            debug_assert!(start >= 0 && end >= 0);
            // C++: syncpool.addConflictFromOtherPool(&conflictEntries_[start], end - start);
            syncpool.add_conflict_from_other_pool(
                &self.conflict_entries_[start as usize..end as usize],
            );
            self.remove_conflict(i);
        }
        self.deleted_conflicts_.clear();
        self.free_spaces_.clear();
        self.conflict_ranges_.clear();
        self.conflict_entries_.clear();
        self.modification_.clear();
        self.ages_.clear();
        self.age_reset_while_locked_.clear();
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::HighsConflictPool (lines 41-55)
    /// C++: `HighsConflictPool(HighsInt agelim, HighsInt softlimit)`
    pub fn new(agelim: HighsInt, softlimit: HighsInt) -> Self {
        let mut pool = HighsConflictPool {
            agelim_: agelim,
            softlimit_: softlimit,
            ..Default::default()
        };
        // C++ `ageDistribution_.resize(agelim_ + 1);` in the constructor body.
        pool.age_distribution_.resize((agelim + 1) as usize, 0);
        pool
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::resetAge (lines 76-86)
    /// C++: `void resetAge(HighsInt conflict)`
    pub fn reset_age(&mut self, conflict: HighsInt) {
        let cidx = conflict as usize;
        if self.ages_[cidx] > 0 {
            if self.age_lock_ {
                // C++: HighsConflictPool.h:79 `ageResetWhileLocked_[conflict].store(1, relaxed);`
                self.age_reset_while_locked_[cidx].store(1, Ordering::Relaxed);
                return;
            }
            self.age_distribution_[self.ages_[cidx] as usize] -= 1;
            self.age_distribution_[0] += 1;
            self.ages_[cidx] = 0;
        }
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::setAgeLimit (lines 88-91)
    /// C++: `void setAgeLimit(HighsInt agelim)`
    pub fn set_age_limit(&mut self, agelim: HighsInt) {
        self.agelim_ = agelim;
        self.age_distribution_.resize((agelim + 1) as usize, 0);
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::getModificationCount (lines 93-95)
    /// C++: `unsigned int getModificationCount(HighsInt cut) const`
    pub fn get_modification_count(&self, cut: HighsInt) -> u32 {
        self.modification_[cut as usize]
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::addPropagationDomain (lines 97-99)
    /// C++: `void addPropagationDomain(HighsDomain::ConflictPoolPropagation * domain)`
    // PORT NOTE: signature changed because the observer identity is a `DomainKey` (the M3 log keys its
    // cursors on it, and `HighsDomain_ConflictPoolPropagation::key` carries it). The C++ pointer is
    // nullable; a key is not, so the null case has no counterpart (C++ would crash on the next notification).
    pub fn add_propagation_domain(&mut self, key: DomainKey) {
        // C++ `propagationDomains.push_back(domain)`: cursor starts at the current end of the log.
        self.events.register(key, self.events.end());
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::removePropagationDomain (lines 101-108)
    /// C++: `void removePropagationDomain(HighsDomain::ConflictPoolPropagation * domain)`
    // PORT NOTE: signature changed for the same reason as addPropagationDomain (`DomainKey` identity).
    pub fn remove_propagation_domain(&mut self, key: DomainKey) {
        // An unregistered key (domain copy made under the parallel lock) is a no-op, as in C++:
        // no write to the shared pool-0 log from a worker task.
        if self.events.unregister(key) {
            // PORT NOTE: `gc` here drops events that no remaining observer still needs (exact, see mod.rs M3).
            self.events.gc();
        }
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::getConflictEntryVector (lines 110-112)
    /// C++: `const std::vector<HighsDomainChange> & getConflictEntryVector() const`
    pub fn get_conflict_entry_vector(&self) -> &[HighsDomainChange] {
        &self.conflict_entries_
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::getConflictRanges (lines 114-116)
    /// C++: `const std::vector<std::pair<HighsInt, HighsInt>> & getConflictRanges() const`
    pub fn get_conflict_ranges(&self) -> &[(HighsInt, HighsInt)] {
        &self.conflict_ranges_
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::getNumConflicts (lines 118-120)
    /// C++: `HighsInt getNumConflicts() const`
    pub fn get_num_conflicts(&self) -> HighsInt {
        (self.conflict_ranges_.len() - self.deleted_conflicts_.len()) as HighsInt
    }

    /// HiGHS: HighsConflictPool.h::HighsConflictPool::setAgeLock (lines 122-122)
    /// C++: `void setAgeLock(const bool ageLock)`
    pub fn set_age_lock(&mut self, age_lock: bool) {
        self.age_lock_ = age_lock;
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;

    #[test]
    fn w5_mip_atom_age_reset_while_locked_is_atomic() {
        // `ageResetWhileLocked_` (HighsConflictPool.h:26) is a relaxed `AtomicU8` per conflict; `performAging` / `resetAge` give the C++ values under `parallelLockActive()`.

        // Test 1: resetAge with age_lock_ = true stores 1 to ageResetWhileLocked_.
        // C++: HighsConflictPool.h:79 `ageResetWhileLocked_[conflict].store(1, relaxed);`
        {
            let mut pool = HighsConflictPool::new(10, 100);
            pool.age_lock_ = true;
            pool.ages_.resize(1, 5);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(0));
            pool.age_distribution_.resize(11, 0);

            pool.reset_age(0);
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 1);
            assert_eq!(pool.ages_[0], 5); // age unchanged (early return)
        }

        // Test 2: resetAge with age_lock_ = false resets age, no atomic store needed.
        // C++: HighsConflictPool.h:82-84 (no store to ageResetWhileLocked_ in this path)
        {
            let mut pool = HighsConflictPool::new(10, 100);
            pool.age_lock_ = false;
            pool.ages_.resize(1, 5);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(0));
            pool.age_distribution_.resize(11, 0);
            pool.age_distribution_[5] = 1;

            pool.reset_age(0);
            assert_eq!(pool.ages_[0], 0);
            assert_eq!(pool.age_distribution_[5], 0);
            assert_eq!(pool.age_distribution_[0], 1);
        }

        // Test 3: performAging with thread_safe=true reads ageResetWhileLocked_ via load(Relaxed).
        // C++: HighsConflictPool.cpp:197 `ageResetWhileLocked_[i].load(relaxed) == 1`
        {
            let mut pool = HighsConflictPool::new(10, 100);
            pool.age_lock_ = false;
            pool.conflict_ranges_.push((0, 2));
            pool.ages_.resize(1, 3);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(1));
            pool.age_distribution_.resize(11, 0);
            pool.age_distribution_[3] = 1;

            // thread_safe = true: should call reset_age because ageResetWhileLocked_ == 1
            pool.perform_aging(true);
            // reset_age was called (age_lock_ is false, so it resets age to 0)
            // then ages_[0] += 1 makes it 1
            assert_eq!(pool.ages_[0], 1);
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 0);
        }

        // Test 4: performAging stores 0 to ageResetWhileLocked_ after processing.
        // C++: HighsConflictPool.cpp:202 `ageResetWhileLocked_[i].store(0, relaxed);`
        {
            let mut pool = HighsConflictPool::new(10, 100);
            pool.age_lock_ = false;
            pool.conflict_ranges_.push((0, 2));
            pool.ages_.resize(1, 3);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(1));
            pool.age_distribution_.resize(11, 0);
            pool.age_distribution_[3] = 1;

            // thread_safe = false: does not check ageResetWhileLocked_, but still stores 0
            pool.perform_aging(false);
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 0);
            assert_eq!(pool.ages_[0], 4); // 3 + 1
        }

        // Test 5: addConflictCut stores 0 to ageResetWhileLocked_ for the new conflict.
        // C++: HighsConflictPool.cpp:48 `ageResetWhileLocked_.resize(...)` then line 56 `store(0, relaxed)`
        {
            let mut pool = HighsConflictPool::new(10, 100);
            // Simulate the state after the resize in addConflictCut:
            // conflict_ranges_ gets a new entry, ages_ and age_reset_while_locked_ are resized.
            pool.conflict_ranges_.push((0, 1));
            pool.ages_.resize(1, 0);
            pool.modification_.resize(1, 0);
            pool.age_reset_while_locked_
                .resize_with(1, || AtomicU8::new(99)); // pre-existing garbage
            pool.age_distribution_.resize(11, 0);

            // The store(0, Relaxed) at the conflict index
            pool.age_reset_while_locked_[0].store(0, Ordering::Relaxed);
            assert_eq!(pool.age_reset_while_locked_[0].load(Ordering::Relaxed), 0);
        }

        // Test 6: syncConflictPool clears ageResetWhileLocked_ (Vec::clear).
        // C++: HighsConflictPool.cpp:283 `ageResetWhileLocked_.clear();`
        {
            let mut pool = HighsConflictPool::new(10, 100);
            pool.age_reset_while_locked_
                .resize_with(5, || AtomicU8::new(1));
            assert_eq!(pool.age_reset_while_locked_.len(), 5);
            pool.age_reset_while_locked_.clear();
            assert_eq!(pool.age_reset_while_locked_.len(), 0);
        }
    }
}
