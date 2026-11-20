//! HiGHS source: presolve/ (LP subset: HPresolve, HighsPostsolveStack,
//! HPresolveAnalysis, PresolveComponent)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! LP presolve, serial, `mipsolver == nullptr`. This file holds the small shared helpers:
//!
//! * [`OrderedF64`]: key wrapper for C++ `std::set<double>` /
//!   `std::map<double, T>` (`BTreeSet<OrderedF64>` / `BTreeMap<OrderedF64, T>`).
//! * [`HighsMultiBuckets`]: the libc++ `std::unordered_multimap<u64,
//!   HighsInt>` of `detectParallelRowsAndCols`.
//! * [`HighsMipSolver`] and [`mip_out_of_scope!`] (the latter is no longer used: the
//!   MIP branches of the public `Highs::presolve` / `Highs::postsolve` API are in
//!   `lp_data/highs_presolve_api.rs`).

/// The body of any MIP-only code path (`unimplemented!("MIP presolve: out of
/// scope")`). Usable as a statement or as an expression of any type. Defined
/// before the submodules so it is in textual scope in all of them.
#[macro_export]
macro_rules! mip_out_of_scope {
    () => {
        unimplemented!("MIP presolve: out of scope")
    };
}

pub mod highs_postsolve_stack;
pub mod highs_symmetry;
pub mod hpresolve;
pub mod hpresolve_analysis;
pub mod hpresolve_clique_table;
pub mod icrash;
pub mod icrash_util;
pub mod icrash_x;
pub mod presolve_component;

use crate::prelude::*;
use std::cmp::Ordering;

/// C++ `HighsMipSolver` (the MIP port, `crate::mip`). `HPresolve::mipsolver`
/// stays `None` on the LP path and is set in MIP mode (`ok_set_input_mip`).
pub use crate::mip::HighsMipSolver;

/// `f64` key with C++ `std::less<double>` semantics for `std::set<double>` /
/// `std::map<double, T>`: `a < b` is `Less`, `b < a` is `Greater`, anything
/// else (including `-0.0` vs `0.0`, which C++ treats as the same key) is
/// `Equal`. NaN keys are UB in C++ and never occur.
#[derive(Clone, Copy, Debug, Default)]
pub struct OrderedF64(pub f64);

impl PartialEq for OrderedF64 {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for OrderedF64 {}
impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.0 < other.0 {
            Ordering::Less
        } else if other.0 < self.0 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

/// C++ `std::unordered_multimap<std::uint64_t, HighsInt> buckets`
/// (HPresolve.cpp 7345, `detectParallelRowsAndCols`), modelled on libc++.
///
/// The C++ only ever touches one equal-key range at a time:
/// `find(key)` (first element of the range), `it++` while `it->first == key`,
/// `erase(last)` and `emplace_hint(last, key, value)`. Only the order
/// *within* one key's range is observable, so each key maps to a `Vec` in
/// libc++ iteration order. libc++ `emplace_hint(hint, ...)` with `hint`
/// holding an equal key links the new node immediately **before** `hint`;
/// with `hint == end()` (no element with this key) it starts a new range.
///
/// Translation pattern (`it`, `last` are indices into the range; `None` =
/// `end()`):
///
/// ```ignore
/// let key = self.col_hashes[i];              // auto it = buckets.find(key);
/// let mut it: usize = 0;                     // decltype(it) last = it;
/// let mut last: Option<usize> = if buckets.count(key) > 0 { Some(0) } else { None };
/// while it < buckets.count(key) {            // while (it != end && it->first == key)
///     let cand = buckets.get(key, it);       //   it->second
///     last = Some(it);                       //   last = it++;
///     it += 1;
///     ...
/// }
/// buckets.erase(key, last.unwrap());         // buckets.erase(last);
/// buckets.emplace_hint(last, key, i);        // buckets.emplace_hint(last, key, i);
/// ```
/// `if (it != buckets.end())` right after `find` is `buckets.count(key) > 0`.
#[derive(Clone, Debug, Default)]
pub struct HighsMultiBuckets {
    map: HashMap<u64, Vec<HighsInt>>,
}

impl HighsMultiBuckets {
    /// Number of elements with this key (the length of the equal range).
    pub fn count(&self, key: u64) -> usize {
        self.map.get(&key).map_or(0, |v| v.len())
    }

    /// The `idx`-th element of the equal range of `key`.
    pub fn get(&self, key: u64, idx: usize) -> HighsInt {
        self.map[&key][idx]
    }

    /// `erase(it)` for the element at index `idx` of `key`'s range.
    pub fn erase(&mut self, key: u64, idx: usize) {
        let v = self.map.get_mut(&key).expect("erase: key not present");
        v.remove(idx);
        if v.is_empty() {
            self.map.remove(&key);
        }
    }

    /// libc++ `emplace_hint(hint, key, value)`: insert before index `hint`
    /// of `key`'s range, or start the range if `hint` is `None` (`end()`).
    pub fn emplace_hint(&mut self, hint: Option<usize>, key: u64, value: HighsInt) {
        let v = self.map.entry(key).or_default();
        match hint {
            Some(h) => v.insert(h, value),
            None => v.push(value),
        }
    }

    /// `clear()`.
    pub fn clear(&mut self) {
        self.map.clear();
    }
}

/// HiGHS: util/HighsMemoryAllocation.h::okResize (lines 19-28)
///
/// `okResize(v, n)` (default `value = T{}`) -> `ok_resize(&mut v, n,
/// Default::default())`. `std::bad_alloc` is not caught in Rust (allocation
/// failure aborts), so this always returns `true`. `resize` keeps existing
/// entries, as in C++.
pub fn ok_resize<T: Clone>(use_vector: &mut Vec<T>, dimension: HighsInt, value: T) -> bool {
    use_vector.resize(dimension as usize, value);
    true
}

/// HiGHS: util/HighsMemoryAllocation.h::okReserve (lines 30-39)
///
/// Capacity is not observable; always `true`. The `unordered_map` overload
/// (lines 41-50) is `map.reserve(n as usize)` at the call site.
pub fn ok_reserve<T>(use_vector: &mut Vec<T>, dimension: HighsInt) -> bool {
    let n = dimension as usize;
    if n > use_vector.len() {
        use_vector.reserve(n - use_vector.len());
    }
    true
}

/// HiGHS: util/HighsMemoryAllocation.h::okAssign (lines 52-61)
///
/// `okAssign(v, n)` (default `value = T{}`) -> `ok_assign(&mut v, n,
/// Default::default())`. Replaces every entry (`std::vector::assign`).
pub fn ok_assign<T: Clone>(use_vector: &mut Vec<T>, dimension: HighsInt, value: T) -> bool {
    use_vector.clear();
    use_vector.resize(dimension as usize, value);
    true
}
