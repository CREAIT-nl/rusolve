//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! Re-expresses the algorithms (not the text) of libstdc++ 13.3 (GCC 13.3.0; GPL-3.0 with the GCC Runtime Library Exception) distributions and std::shuffle / std::sort. See THIRD_PARTY_NOTICES.md.
//!
//! The distributions, `std::shuffle` and `std::sort` ported here are the
//! **GCC 13.3 libstdc++** algorithms (the compiler the oracle is built with).
//! A libc++ or MSVC build of upstream gives other traces: the
//! `uniform_int_distribution` downscaling, `generate_canonical`, the
//! two-at-a-time `std::shuffle` and the introsort tie-breaking are all
//! implementation-specific.
//!
//!

use crate::HighsInt;
use crate::lp_data::hconst::HighsVarType;
use crate::mip::fpr_local_mip::{
    CscMatrix, HeuristicBudget, K_BASE_SEED_OFFSET, K_BUDGET_BASE_SHIFT,
    K_PATIENCE_CEILING_DIVISOR, K_SEED_STRIDE, RestartPoint,
};
use crate::util::highs_utils::{cpp_max, cpp_min};

// ---------------------------------------------------------------------------
// S5Rng: xoshiro256++ (rng.h:1-85)
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): rng.h::Xoshiro256PlusPlus (lines 26-81)
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct S5Rng {
    s: [u64; 4],
}

impl S5Rng {
    /// mip-heuristics (36f7d9d): rng.h::Xoshiro256PlusPlus::reseed (lines 38-51)
    pub(crate) fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut s = [0u64; 4];
        s[0] = splitmix64(&mut z);
        s[1] = splitmix64(&mut z);
        s[2] = splitmix64(&mut z);
        s[3] = splitmix64(&mut z);
        if (s[0] | s[1] | s[2] | s[3]) == 0 {
            s[0] = 0x9E3779B97F4A7C15; // golden-ratio fallback
        }
        S5Rng { s }
    }

    /// mip-heuristics (36f7d9d): rng.h::Xoshiro256PlusPlus::operator() (lines 53-63)
    pub(crate) fn next_u64(&mut self) -> u64 {
        let result = rotl(self.s[0].wrapping_add(self.s[3]), 23).wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = rotl(self.s[3], 45);
        result
    }
}

/// mip-heuristics (36f7d9d): rng.h::Xoshiro256PlusPlus::rotl (lines 66-68)
fn rotl(x: u64, k: u32) -> u64 {
    x.rotate_left(k)
}

/// mip-heuristics (36f7d9d): rng.h::Xoshiro256PlusPlus::splitmix64 (lines 72-78)
fn splitmix64(z: &mut u64) -> u64 {
    *z = z.wrapping_add(0x9E3779B97F4A7C15);
    let mut x = *z;
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

// ---------------------------------------------------------------------------
// libstdc++ uniform_int_distribution (bits/uniform_int_dist.h:257-373)
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/uniform_int_dist.h::_S_nd (lines 257-281)
/// Lemire's fast random integer generation with unsigned __int128.
fn s_nd_u128(rng: &mut S5Rng, range: u64) -> u64 {
    let mut product: u128 = (rng.next_u64() as u128) * (range as u128);
    let mut low: u64 = product as u64;
    if low < range {
        let threshold: u64 = (0u64.wrapping_sub(range)) % range;
        while low < threshold {
            product = (rng.next_u64() as u128) * (range as u128);
            low = product as u64;
        }
    }
    (product >> 64) as u64
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/uniform_int_dist.h::operator() (lines 284-373)
/// Ported for the generator with min()=0, max()=u64::MAX (urngrange = u64::MAX).
/// The upscaling branch (urngrange < urange) is unreachable for this generator
/// because urngrange = u64::MAX >= urange for any urange <= u64::MAX.
fn s5_uniform_int_u64_impl(rng: &mut S5Rng, a: u64, b: u64) -> u64 {
    // __urngrange = u64::MAX - 0 = u64::MAX
    // __urange = (b as u64) - (a as u64)
    let urange: u64 = b.wrapping_sub(a);

    // Upscaling is unreachable for this generator (urngrange = u64::MAX).
    let ret: u64 = if urange < u64::MAX {
        // downscaling: urngrange > urange; __uerange = __urange + 1
        // __urngrange == UINT64_MAX branch: use _S_nd<unsigned __int128>
        s_nd_u128(rng, urange.wrapping_add(1))
    } else {
        // __urngrange == __urange: full range
        rng.next_u64()
    };

    // return __ret + __param.a()
    ret.wrapping_add(a)
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/uniform_int_dist.h::operator() for i32
pub(crate) fn s5_uniform_int_i32(rng: &mut S5Rng, a: i32, b: i32) -> i32 {
    let au = a as u64;
    let bu = b as u64;
    let ret = s5_uniform_int_u64_impl(rng, au, bu);
    ret as i32
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/uniform_int_dist.h::operator() for i64
pub(crate) fn s5_uniform_int_i64(rng: &mut S5Rng, a: i64, b: i64) -> i64 {
    let au = a as u64;
    let bu = b as u64;
    let ret = s5_uniform_int_u64_impl(rng, au, bu);
    ret as i64
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/uniform_int_dist.h::operator() for u64
pub(crate) fn s5_uniform_int_u64(rng: &mut S5Rng, a: u64, b: u64) -> u64 {
    s5_uniform_int_u64_impl(rng, a, b)
}

// ---------------------------------------------------------------------------
// libstdc++ generate_canonical (bits/random.tcc:3349-3380)
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/random.tcc::generate_canonical (lines 3349-3380)
/// For double with digits=53, urng min=0 max=u64::MAX:
///   __b = min(53, 53) = 53
///   __r = (u64::MAX - 0 + 1) as long double = 2^64
///   __log2r = log(2^64)/log(2) = 64
///   __m = max(1, (53 + 64 - 1) / 64) = max(1, 116/64) = max(1, 1) = 1
/// So one draw: sum = (x - 0) as f64, tmp = 1, then ret = sum / 2^64.
pub(crate) fn s5_generate_canonical(rng: &mut S5Rng) -> f64 {
    let x = rng.next_u64();
    let mut ret = (x as f64) / 18446744073709551616.0;
    if ret >= 1.0 {
        ret = 1.0 - f64::EPSILON / 2.0;
    }
    ret
}

// ---------------------------------------------------------------------------
// libstdc++ uniform_real_distribution (bits/random.h:1904-1910)
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/random.h::uniform_real_distribution::operator() (lines 1904-1910)
/// return (__aurng() * (__p.b() - __p.a())) + __p.a();
pub(crate) fn s5_uniform_real(rng: &mut S5Rng, a: f64, b: f64) -> f64 {
    let canonical = s5_generate_canonical(rng);
    (canonical * (b - a)) + a
}

// ---------------------------------------------------------------------------
// libstdc++ bernoulli_distribution (bits/random.h:3741-3750)
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/random.h::bernoulli_distribution::operator() (lines 3741-3750)
/// _Adaptor has min()=0, max()=1, so the test is:
///   generate_canonical < p * (1.0 - 0.0)
/// which simplifies to generate_canonical < p since p * 1.0 == p exactly.
pub(crate) fn s5_bernoulli(rng: &mut S5Rng, p: f64) -> bool {
    s5_generate_canonical(rng) < p * (1.0 - 0.0)
}

// ---------------------------------------------------------------------------
// libstdc++ std::shuffle (bits/stl_algo.h:3739-3805)
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__gen_two_uniform_ints (lines 3717-3725)
fn gen_two_uniform_ints(rng: &mut S5Rng, b0: u64, b1: u64) -> (u64, u64) {
    let x = s5_uniform_int_u64(rng, 0, b0.wrapping_mul(b1).wrapping_sub(1));
    (x / b1, x % b1)
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::shuffle (lines 3739-3805)
pub(crate) fn s5_shuffle<T>(v: &mut [T], rng: &mut S5Rng) {
    let n = v.len();
    if n == 0 {
        return;
    }

    // __urngrange = u64::MAX - 0 = u64::MAX
    // __urange = n as u64
    let urngrange: u64 = u64::MAX;
    let urange: u64 = n as u64;

    if urngrange / urange >= urange {
        // Two-at-a-time path
        let mut i: usize = 1;

        // If urange is even, do the first swap up front
        if urange.is_multiple_of(2) {
            let j = s5_uniform_int_u64(rng, 0, 1) as usize;
            v.swap(i, j);
            i += 1;
        }

        // Now (n - i) is even, process in pairs
        while i != n {
            let swap_range: u64 = (i as u64) + 1;
            let (pos_first, pos_second) = gen_two_uniform_ints(rng, swap_range, swap_range + 1);
            v.swap(i, pos_first as usize);
            i += 1;
            v.swap(i, pos_second as usize);
            i += 1;
        }

        return;
    }

    // General one-at-a-time fallback
    for i in 1..n {
        let j = s5_uniform_int_u64(rng, 0, i as u64) as usize;
        v.swap(i, j);
    }
}

// ---------------------------------------------------------------------------
// libstdc++ std::sort (bits/stl_algo.h:85-103, 1792-1952)
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__move_median_to_first (lines 85-103)
fn move_median_to_first<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    comp: &mut F,
) {
    if comp(&v[a], &v[b]) {
        if comp(&v[b], &v[c]) {
            v.swap(result, b);
        } else if comp(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if comp(&v[a], &v[c]) {
        v.swap(result, a);
    } else if comp(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__unguarded_partition (lines 1871-1887)
fn unguarded_partition<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    pivot: usize,
    comp: &mut F,
) -> usize {
    let mut first = first;
    let mut last = last;
    loop {
        while comp(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while comp(&v[pivot], &v[last]) {
            last -= 1;
        }
        if !(first < last) {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__unguarded_partition_pivot (lines 1893-1900)
fn unguarded_partition_pivot<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(v, first, first + 1, mid, last - 1, comp);
    unguarded_partition(v, first + 1, last, first, comp)
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__unguarded_linear_insert (lines 1788-1805)
/// `__val` is a clone of `*__last`; the hole moves left as in C++, with no
/// lower-bound guard (the caller guarantees a sentinel, as in C++).
fn unguarded_linear_insert<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    last: usize,
    comp: &mut F,
) {
    let val = v[last].clone();
    let mut last = last;
    let mut next = last;
    next = next.wrapping_sub(1);
    while comp(&val, &v[next]) {
        v[last] = v[next].clone();
        last = next;
        next = next.wrapping_sub(1);
    }
    v[last] = val;
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__insertion_sort (lines 1808-1829)
fn insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if first == last {
        return;
    }
    let mut i = first + 1;
    while i != last {
        if comp(&v[i], &v[first]) {
            let val = v[i].clone();
            // _GLIBCXX_MOVE_BACKWARD3(__first, __i, __i + 1)
            let mut j = i;
            while j != first {
                v[j] = v[j - 1].clone();
                j -= 1;
            }
            v[first] = val;
        } else {
            unguarded_linear_insert(v, i, comp);
        }
        i += 1;
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__unguarded_insertion_sort (lines 1836-1842)
fn unguarded_insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    let mut i = first;
    while i != last {
        unguarded_linear_insert(v, i, comp);
        i += 1;
    }
}

const S_THRESHOLD: usize = 16;

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__final_insertion_sort (lines 1854-1865)
fn final_insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if last - first > S_THRESHOLD {
        insertion_sort(v, first, first + S_THRESHOLD, comp);
        unguarded_insertion_sort(v, first + S_THRESHOLD, last, comp);
    } else {
        insertion_sort(v, first, last, comp);
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__introsort_loop (lines 1918-1935)
fn introsort_loop<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    comp: &mut F,
) {
    while last - first > S_THRESHOLD {
        if depth_limit == 0 {
            // __partial_sort(first, last, last, comp)
            // This is heap_select + sort_heap. Port from bits/stl_heap.h.
            heap_select(v, first, last, last, comp);
            sort_heap(v, first, last, comp);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(v, first, last, comp);
        introsort_loop(v, cut, last, depth_limit, comp);
        last = cut;
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__sort (lines 1942-1952)
fn std_sort_impl<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if first != last {
        let n = last - first;
        // std::__lg(n) = floor(log2(n))
        let depth_limit = 2 * (usize::BITS - 1 - n.leading_zeros()) as usize;
        introsort_loop(v, first, last, depth_limit, comp);
        final_insertion_sort(v, first, last, comp);
    }
}

// ---------------------------------------------------------------------------
// libstdc++ heap operations (bits/stl_heap.h) for __partial_sort
// ---------------------------------------------------------------------------

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__parent (the `(__holeIndex - 1) / 2` of __push_heap, C++ signed division)
fn heap_parent(hole: usize) -> usize {
    if hole == 0 { 0 } else { (hole - 1) / 2 }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__push_heap (lines 135-147)
/// `_Tp __value` is held by the caller and moved in as `value`; `first` is the
/// heap base and `hole`/`top` are distances from it, as in C++.
fn push_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    top: usize,
    value: T,
    comp: &mut F,
) {
    let mut parent = heap_parent(hole);
    while hole > top && comp(&v[first + parent], &value) {
        v[first + hole] = v[first + parent].clone();
        hole = parent;
        parent = heap_parent(hole);
    }
    v[first + hole] = value;
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__adjust_heap (lines 224-249)
fn adjust_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    comp: &mut F,
) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if comp(&v[first + second], &v[first + second - 1]) {
            second -= 1;
        }
        v[first + hole] = v[first + second].clone();
        hole = second;
    }
    if len.is_multiple_of(2) && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[first + hole] = v[first + second - 1].clone();
        hole = second - 1;
    }
    push_heap(v, first, hole, top, value, comp);
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__pop_heap (lines 254-265)
fn pop_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    result: usize,
    comp: &mut F,
) {
    let value = v[result].clone();
    v[result] = v[first].clone();
    adjust_heap(v, first, 0, last - first, value, comp);
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__make_heap (lines 340-362)
fn make_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if last - first < 2 {
        return;
    }
    let len = last - first;
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent].clone();
        adjust_heap(v, first, parent, len, value, comp);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_heap.h::__sort_heap (lines 419-427)
fn sort_heap<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    comp: &mut F,
) {
    while last - first > 1 {
        last -= 1;
        pop_heap(v, first, last, last, comp);
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::__heap_select (lines 1631-1639)
fn heap_select<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    middle: usize,
    last: usize,
    comp: &mut F,
) {
    make_heap(v, first, middle, comp);
    for i in middle..last {
        if comp(&v[i], &v[first]) {
            pop_heap(v, first, middle, i, comp);
        }
    }
}

/// libstdc++ 13.3 (GCC 13.3.0): bits/stl_algo.h::sort (lines 1942-1952)
pub(crate) fn s5_std_sort<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    let n = v.len();
    std_sort_impl(v, 0, n, &mut comp);
}

// ---------------------------------------------------------------------------
// heuristic_common.h helpers
// ---------------------------------------------------------------------------

/// mip-heuristics (36f7d9d): heuristic_common.h::heuristic_base_seed (lines 159-161)
pub(crate) fn heuristic_base_seed(random_seed: HighsInt) -> u32 {
    (random_seed as u32).wrapping_add(K_BASE_SEED_OFFSET)
}

/// mip-heuristics (36f7d9d): heuristic_context.h::worker_seed (lines 509-511)
pub(crate) fn worker_seed(base_seed: u32, w: i32) -> u32 {
    base_seed.wrapping_add((w as u32).wrapping_mul(K_SEED_STRIDE))
}

/// mip-heuristics (36f7d9d): heuristic_common.h::heuristic_effort_budget (lines 208-217)
pub(crate) fn heuristic_effort_budget(nnz: usize, effort: f64) -> usize {
    if effort <= 0.0 {
        return 0;
    }
    let budget = ((nnz << K_BUDGET_BASE_SHIFT) as f64) * effort;
    if !(budget < (usize::MAX as f64)) {
        return usize::MAX;
    }
    budget as usize
}

/// mip-heuristics (36f7d9d): heuristic_common.h::vanilla_effort_budget (lines 202-206)
pub(crate) fn vanilla_effort_budget(nnz: usize, mip_heuristic_effort: f64) -> usize {
    // kBaseRatio = 4.0, kVanillaAnchor = 0.05
    // mip_heuristic_effort * kBaseRatio / kVanillaAnchor
    heuristic_effort_budget(nnz, mip_heuristic_effort * 4.0 / 0.05)
}

/// mip-heuristics (36f7d9d): heuristic_common.h::saturating_mul (lines 231-236)
pub(crate) fn saturating_mul(a: usize, b: usize) -> usize {
    if a == 0 || b == 0 {
        return 0;
    }
    if a > usize::MAX / b {
        usize::MAX
    } else {
        a * b
    }
}

/// mip-heuristics (36f7d9d): heuristic_common.h::patience_threshold (lines 324-334)
pub(crate) fn patience_threshold(nnz: usize, per_base: f64, budget: usize) -> usize {
    if per_base <= 0.0 {
        return usize::MAX;
    }
    let threshold = std::cmp::max(heuristic_effort_budget(nnz, per_base), 1);
    if budget == 0 {
        return threshold;
    }
    let ceiling = std::cmp::max(budget / K_PATIENCE_CEILING_DIVISOR, 1);
    std::cmp::min(threshold, ceiling)
}

/// mip-heuristics (36f7d9d): heuristic_context.h::make_budget (lines 521-543)
pub(crate) fn make_budget(total: usize, num_workers: usize, stale: usize) -> HeuristicBudget {
    if total == 0 {
        return HeuristicBudget::default();
    }
    HeuristicBudget {
        total,
        per_worker: total / num_workers,
        attempt_cap: std::cmp::max(total / (num_workers * 10), 1),
        stale,
        worker_stale: std::cmp::max(stale / num_workers, 1),
        fj_stall_ends_run: false,
        restart: RestartPoint::default(),
        randomize_restart_weights: false,
    }
}

/// mip-heuristics (36f7d9d): heuristic_context.h::make_until_stopped_budget (lines 665-672)
pub(crate) fn make_until_stopped_budget(attempt_cap: usize, patience: usize) -> HeuristicBudget {
    HeuristicBudget {
        total: usize::MAX,
        per_worker: usize::MAX,
        attempt_cap,
        stale: usize::MAX,
        worker_stale: if patience == 0 { usize::MAX } else { patience },
        fj_stall_ends_run: false,
        restart: RestartPoint::default(),
        randomize_restart_weights: false,
    }
}

/// mip-heuristics (36f7d9d): heuristic_common.h::build_csc (lines 63-89)
pub(crate) fn build_csc(
    ncol: HighsInt,
    nrow: HighsInt,
    ar_start: &[HighsInt],
    ar_index: &[HighsInt],
    ar_value: &[f64],
) -> CscMatrix {
    let nnz = ar_index.len() as HighsInt;
    let mut csc = CscMatrix {
        col_start: vec![0; (ncol + 1) as usize],
        col_row: vec![0; nnz as usize],
        col_val: vec![0.0; nnz as usize],
    };
    for k in 0..nnz {
        csc.col_start[(ar_index[k as usize] + 1) as usize] += 1;
    }
    for j in 0..ncol {
        csc.col_start[(j + 1) as usize] += csc.col_start[j as usize];
    }
    {
        let mut pos: Vec<HighsInt> = csc.col_start.clone();
        for i in 0..nrow {
            let start = ar_start[i as usize] as usize;
            let end = ar_start[(i + 1) as usize] as usize;
            for k in start..end {
                let j = ar_index[k] as usize;
                let p = pos[j] as usize;
                csc.col_row[p] = i;
                csc.col_val[p] = ar_value[k];
                pos[j] += 1;
            }
        }
    }
    csc
}

/// mip-heuristics (36f7d9d): heuristic_common.h::is_integer (lines 91-93)
pub(crate) fn is_integer(integrality: &[HighsVarType], j: HighsInt) -> bool {
    integrality[j as usize] != HighsVarType::kContinuous
}

/// mip-heuristics (36f7d9d): heuristic_common.h::row_violation (lines 108-110)
pub(crate) fn row_violation(lhs: f64, lo: f64, hi: f64) -> f64 {
    cpp_max(0.0, lhs - hi) + cpp_max(0.0, lo - lhs)
}

/// mip-heuristics (36f7d9d): heuristic_common.h::is_row_violated (lines 113-115)
pub(crate) fn is_row_violated(lhs: f64, lo: f64, hi: f64, feastol: f64) -> bool {
    lhs > hi + feastol || lhs < lo - feastol
}

/// mip-heuristics (36f7d9d): heuristic_common.h::clamp_round (lines 118-123)
pub(crate) fn clamp_round(val: f64, lb: f64, ub: f64, integer: bool) -> f64 {
    let mut val = val;
    if integer {
        val = val.round();
    }
    cpp_max(lb, cpp_min(ub, val))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const RNG_FIXTURE: &str = include_str!("../../../tests/fixtures/s5/rng.txt");

    /// Returns (header, values) for one fixture line, split at the " : " separator.
    fn fixture_parts(line: &'static str) -> (&'static str, &'static str) {
        line.split_once(" : ")
            .unwrap_or_else(|| panic!("fixture line without separator: {}", line))
    }

    fn fixture_lines(prefix: &str) -> Vec<&'static str> {
        RNG_FIXTURE
            .lines()
            .filter(|l| l.starts_with(prefix))
            .collect()
    }

    fn hex_bits(s: &str) -> u64 {
        u64::from_str_radix(s.strip_prefix("0x").unwrap_or(s), 16).unwrap()
    }

    fn header_value<'a>(parts: &[&'a str], key: &str) -> &'a str {
        parts
            .iter()
            .find_map(|p| p.strip_prefix(key))
            .unwrap_or_else(|| panic!("missing {} in header", key))
    }

    #[test]
    fn s5_rng_raw_draws() {
        for seed in [0u64, 42, 1039] {
            let prefix = format!("raw seed={} :", seed);
            let line = fixture_lines(&prefix)
                .into_iter()
                .next()
                .unwrap_or_else(|| panic!("fixture line not found: {}", prefix));
            let (_, values) = fixture_parts(line);
            let expected: Vec<u64> = values
                .split_whitespace()
                .map(|s| s.parse::<u64>().unwrap())
                .collect();
            let mut rng = S5Rng::new(seed);
            let actual: Vec<u64> = (0..expected.len()).map(|_| rng.next_u64()).collect();
            assert_eq!(actual, expected, "raw draws seed={}", seed);
        }
    }

    #[test]
    fn s5_rng_uniform_int() {
        for line in RNG_FIXTURE.lines() {
            if !line.starts_with("uniform_int_") {
                continue;
            }
            // uniform_int_<t> seed=<s> a=<a> b=<b> : <vals...>
            let (header, values) = fixture_parts(line);
            let parts: Vec<&str> = header.split_whitespace().collect();
            let kind = parts[0];
            let seed: u64 = header_value(&parts, "seed=").parse().unwrap();
            let a_str = header_value(&parts, "a=");
            let b_str = header_value(&parts, "b=");
            let mut rng = S5Rng::new(seed);
            let (expected, actual): (Vec<u64>, Vec<u64>) = if kind == "uniform_int_i32" {
                let a = a_str.parse::<i32>().unwrap();
                let b = b_str.parse::<i32>().unwrap();
                let exp = values
                    .split_whitespace()
                    .map(|s| s.parse::<i64>().unwrap() as u64)
                    .collect::<Vec<u64>>();
                let act = (0..exp.len())
                    .map(|_| s5_uniform_int_i32(&mut rng, a, b) as i64 as u64)
                    .collect();
                (exp, act)
            } else if kind == "uniform_int_i64" {
                let a = a_str.parse::<i64>().unwrap();
                let b = b_str.parse::<i64>().unwrap();
                let exp = values
                    .split_whitespace()
                    .map(|s| s.parse::<i64>().unwrap() as u64)
                    .collect::<Vec<u64>>();
                let act = (0..exp.len())
                    .map(|_| s5_uniform_int_i64(&mut rng, a, b) as u64)
                    .collect();
                (exp, act)
            } else {
                let a = a_str.parse::<u64>().unwrap();
                let b = b_str.parse::<u64>().unwrap();
                let exp = values
                    .split_whitespace()
                    .map(|s| s.parse::<u64>().unwrap())
                    .collect::<Vec<u64>>();
                let act = (0..exp.len())
                    .map(|_| s5_uniform_int_u64(&mut rng, a, b))
                    .collect();
                (exp, act)
            };
            assert_eq!(actual, expected, "uniform_int line: {}", line);
        }
    }

    #[test]
    fn s5_rng_uniform_real_canonical_bernoulli() {
        for line in RNG_FIXTURE.lines() {
            if line.starts_with("uniform_real ") || line.starts_with("uniform_real_temp ") {
                // uniform_real[_temp] seed=<s> a=<hex bits> b=<hex bits> : <hex bits...>
                let (header, values) = fixture_parts(line);
                let parts: Vec<&str> = header.split_whitespace().collect();
                let seed: u64 = header_value(&parts, "seed=").parse().unwrap();
                let a = f64::from_bits(hex_bits(header_value(&parts, "a=")));
                let b = f64::from_bits(hex_bits(header_value(&parts, "b=")));
                let expected: Vec<u64> = values.split_whitespace().map(hex_bits).collect();
                let mut rng = S5Rng::new(seed);
                let actual: Vec<u64> = (0..expected.len())
                    .map(|_| s5_uniform_real(&mut rng, a, b).to_bits())
                    .collect();
                assert_eq!(actual, expected, "uniform_real: {}", line);
            } else if line.starts_with("generate_canonical ") {
                let (header, values) = fixture_parts(line);
                let parts: Vec<&str> = header.split_whitespace().collect();
                let seed: u64 = header_value(&parts, "seed=").parse().unwrap();
                let expected: Vec<u64> = values.split_whitespace().map(hex_bits).collect();
                let mut rng = S5Rng::new(seed);
                let actual: Vec<u64> = (0..expected.len())
                    .map(|_| s5_generate_canonical(&mut rng).to_bits())
                    .collect();
                assert_eq!(actual, expected, "generate_canonical: {}", line);
            } else if line.starts_with("bernoulli ") {
                // bernoulli seed=<s> p=<hex bits> : <0/1...>
                let (header, values) = fixture_parts(line);
                let parts: Vec<&str> = header.split_whitespace().collect();
                let seed: u64 = header_value(&parts, "seed=").parse().unwrap();
                let p = f64::from_bits(hex_bits(header_value(&parts, "p=")));
                let expected: Vec<u64> = values
                    .split_whitespace()
                    .map(|s| s.parse::<u64>().unwrap())
                    .collect();
                let mut rng = S5Rng::new(seed);
                let actual: Vec<u64> = (0..expected.len())
                    .map(|_| s5_bernoulli(&mut rng, p) as u64)
                    .collect();
                assert_eq!(actual, expected, "bernoulli: {}", line);
            }
        }
    }

    #[test]
    fn s5_rng_shuffle() {
        for line in fixture_lines("shuffle ") {
            // shuffle seed=42 n=<n> : <perm...> | next=<u64>
            let (header, values) = fixture_parts(line);
            let parts: Vec<&str> = header.split_whitespace().collect();
            let seed: u64 = header_value(&parts, "seed=").parse().unwrap();
            let n: usize = header_value(&parts, "n=").parse().unwrap();
            let (perm_str, next_str) = values.split_once('|').unwrap();
            let perm: Vec<usize> = perm_str
                .split_whitespace()
                .map(|s| s.parse::<usize>().unwrap())
                .collect();
            let next_val: u64 = next_str
                .trim()
                .strip_prefix("next=")
                .unwrap()
                .parse()
                .unwrap();

            let mut v: Vec<usize> = (0..n).collect();
            let mut rng = S5Rng::new(seed);
            s5_shuffle(&mut v, &mut rng);
            assert_eq!(v, perm, "shuffle n={}", n);
            assert_eq!(rng.next_u64(), next_val, "shuffle next n={}", n);
        }
    }

    #[test]
    fn s5_rng_introsort_ties() {
        for line in fixture_lines("sort ") {
            // sort seed=42 n=<n> nkeys=<k> keys=<k0,k1,...> : <id_order...>
            let (header, values) = fixture_parts(line);
            let parts: Vec<&str> = header.split_whitespace().collect();
            let seed: u64 = header_value(&parts, "seed=").parse().unwrap();
            let n: usize = header_value(&parts, "n=").parse().unwrap();
            let nkeys: u64 = header_value(&parts, "nkeys=").parse().unwrap();
            let fixture_keys: Vec<u64> = header_value(&parts, "keys=")
                .split(',')
                .map(|s| s.parse::<u64>().unwrap())
                .collect();
            let expected_ids: Vec<usize> = values
                .split_whitespace()
                .map(|s| s.parse::<usize>().unwrap())
                .collect();

            let mut rng = S5Rng::new(seed);
            let keys: Vec<u64> = (0..n).map(|_| rng.next_u64() % nkeys).collect();
            assert_eq!(keys, fixture_keys, "sort keys n={} nkeys={}", n, nkeys);
            let mut pairs: Vec<(usize, u64)> = keys.into_iter().enumerate().collect();
            s5_std_sort(&mut pairs, |a, b| a.1 > b.1);
            let actual_ids: Vec<usize> = pairs.iter().map(|p| p.0).collect();
            assert_eq!(actual_ids, expected_ids, "sort n={} nkeys={}", n, nkeys);
        }
    }

    #[test]
    fn s5_rng_budget_edges() {
        let nnz: usize = 1000;

        // heuristic_effort_budget edges
        assert_eq!(heuristic_effort_budget(nnz, 0.0), 0);
        assert_eq!(heuristic_effort_budget(nnz, -1.0), 0);
        assert_eq!(heuristic_effort_budget(nnz, 1e300), usize::MAX);

        // saturating_mul
        assert_eq!(saturating_mul(usize::MAX, 2), usize::MAX);
        assert_eq!(saturating_mul(2, usize::MAX), usize::MAX);
        assert_eq!(saturating_mul(usize::MAX, 0), 0);
        assert_eq!(saturating_mul(3, 5), 15);

        // patience_threshold: zero/negative multiplier -> usize::MAX
        assert_eq!(patience_threshold(nnz, 0.0, 1000), usize::MAX);
        assert_eq!(patience_threshold(nnz, -1.0, 1000), usize::MAX);
        assert_eq!(patience_threshold(nnz, 0.0, 0), usize::MAX);
        assert_eq!(patience_threshold(0, 0.0, 1000), usize::MAX);

        // patience_threshold: threshold does not move with budget
        let k_per_base = 0.25;
        let small = patience_threshold(nnz, k_per_base, 100_000_000);
        let large = patience_threshold(nnz, k_per_base, 100_000_000_000);
        assert_eq!(small, nnz * 256);
        assert_eq!(large, small);
        assert_eq!(
            patience_threshold(2 * nnz, k_per_base, 100_000_000),
            2 * small
        );

        // Clamped to quarter of budget
        assert_eq!(patience_threshold(nnz, k_per_base, 1000), 250);

        // Degenerate: nnz=0 gives threshold 1
        assert_eq!(patience_threshold(0, k_per_base, 1000), 1);

        // budget=0 returns threshold unclamped
        assert_eq!(patience_threshold(nnz, k_per_base, 0), nnz * 256);

        // Clamped strictly below ceiling
        let k_budget: usize = 1_000_000;
        for &per_base in &[1.0, 10.0, 1e6] {
            let threshold = patience_threshold(nnz, per_base, k_budget);
            assert_eq!(threshold, k_budget / K_PATIENCE_CEILING_DIVISOR);
            assert!(threshold < k_budget);
        }

        // Below ceiling: option used as-is
        assert_eq!(patience_threshold(nnz, 0.1, k_budget), 102_400);

        // Small budget: floor at 1
        assert_eq!(patience_threshold(nnz, 1.0, 2), 1);

        // Saturation: patience_threshold with huge nnz
        assert_eq!(patience_threshold(usize::MAX, 2.0, 0), usize::MAX);
        assert_eq!(patience_threshold(usize::MAX / 2, 4.0, 0), usize::MAX);
        assert_eq!(patience_threshold(usize::MAX, 2.0, 1000), 250);

        // make_budget(0, ...) == default
        assert_eq!(make_budget(0, 4, 7), HeuristicBudget::default());

        // make_until_stopped_budget
        let b = make_until_stopped_budget(100, 50);
        assert_eq!(b.total, usize::MAX);
        assert_eq!(b.per_worker, usize::MAX);
        assert_eq!(b.attempt_cap, 100);
        assert_eq!(b.stale, usize::MAX);
        assert_eq!(b.worker_stale, 50);

        let b0 = make_until_stopped_budget(100, 0);
        assert_eq!(b0.worker_stale, usize::MAX);

        // Base seed and worker seed
        assert_eq!(heuristic_base_seed(0), 42);
        assert_eq!(worker_seed(42, 2), 42 + 2 * 997);
    }
}
