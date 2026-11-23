//! HiGHS source: util/HighsSort.cpp (subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Only the routines reached from the simplex path: the decreasing heap used by
//! HEkkPrimal (hyper-sparse CHUZC) and the `double`-valued `maxheapsort` chain used by
//! HEkkDualRow::chooseFinal (BFRT group sort). The `HighsInt`-valued overloads,
//! `increasingSetOk` and `sortSetData` are unused there and omitted.
//! HEAPS ARE 1-BASED: entry 0 of `heap_v`/`heap_i` is unused (C++ indexes `heap_v[1..=n]`).
//! Keep the exact sift order: ties are broken by position, and the BFRT group order
//! depends on it.

use crate::prelude::*;

/// HiGHS: HighsSort.cpp::addToDecreasingHeap (lines 17-67, 51 lines)
pub fn add_to_decreasing_heap(
    n: &mut HighsInt,
    mx_n: HighsInt,
    heap_v: &mut Vec<f64>,
    heap_ix: &mut Vec<HighsInt>,
    v: f64,
    ix: HighsInt,
) {
    let mut cd_p: HighsInt;
    let mut pa_p: HighsInt;
    if *n < mx_n {
        // The heap is not full so put the new value at the bottom of the
        // heap and let it rise up to its correct level.
        *n += 1;
        cd_p = *n;
        pa_p = cd_p / 2;
        loop {
            if pa_p > 0 && v < heap_v[pa_p as usize] {
                heap_v[cd_p as usize] = heap_v[pa_p as usize];
                heap_ix[cd_p as usize] = heap_ix[pa_p as usize];
                cd_p = pa_p;
                pa_p /= 2;
                continue;
            }
            break;
        }
        heap_v[cd_p as usize] = v;
        heap_ix[cd_p as usize] = ix;
    } else if v > heap_v[1] {
        // The heap is full so replace the least value with the new value
        // and let it sink down to its correct level.
        pa_p = 1;
        cd_p = pa_p + pa_p;
        loop {
            if cd_p <= *n {
                if cd_p < *n && heap_v[cd_p as usize] > heap_v[(cd_p + 1) as usize] {
                    cd_p += 1;
                }
                if v > heap_v[cd_p as usize] {
                    heap_v[pa_p as usize] = heap_v[cd_p as usize];
                    heap_ix[pa_p as usize] = heap_ix[cd_p as usize];
                    pa_p = cd_p;
                    cd_p += cd_p;
                    continue;
                }
            }
            break;
        }
        heap_v[pa_p as usize] = v;
        heap_ix[pa_p as usize] = ix;
    }
    // Set heap_ix[0]=1 to indicate that the values form a heap.
    heap_ix[0] = 1;
}

/// HiGHS: HighsSort.cpp::sortDecreasingHeap (lines 69-124, 56 lines)
pub fn sort_decreasing_heap(n: HighsInt, heap_v: &mut Vec<f64>, heap_ix: &mut Vec<HighsInt>) {
    let mut fo_p: HighsInt;
    let mut srt_p: HighsInt;
    let mut cd_p: HighsInt;
    let mut pa_p: HighsInt;
    if n <= 1 {
        return;
    }
    if heap_ix[0] != 1 {
        // The data are assumed to be completely unordered. A heap will be formed
        // and sorted.
        fo_p = n / 2 + 1;
        srt_p = n;
    } else {
        // The data are assumed to form a heap which is to be sorted.
        fo_p = 1;
        srt_p = n;
    }
    loop {
        let v: f64;
        let ix: HighsInt;
        if fo_p > 1 {
            fo_p -= 1;
            v = heap_v[fo_p as usize];
            ix = heap_ix[fo_p as usize];
        } else {
            v = heap_v[srt_p as usize];
            ix = heap_ix[srt_p as usize];
            heap_v[srt_p as usize] = heap_v[1];
            heap_ix[srt_p as usize] = heap_ix[1];
            srt_p -= 1;
            if srt_p == 1 {
                heap_v[1] = v;
                heap_ix[1] = ix;
                return;
            }
        }
        pa_p = fo_p;
        cd_p = fo_p + fo_p;
        loop {
            if cd_p <= srt_p {
                if cd_p < srt_p && heap_v[cd_p as usize] > heap_v[(cd_p + 1) as usize] {
                    cd_p += 1;
                }
                if v > heap_v[cd_p as usize] {
                    heap_v[pa_p as usize] = heap_v[cd_p as usize];
                    heap_ix[pa_p as usize] = heap_ix[cd_p as usize];
                    pa_p = cd_p;
                    cd_p += cd_p;
                    continue;
                }
            }
            break;
        }
        heap_v[pa_p as usize] = v;
        heap_ix[pa_p as usize] = ix;
    }
}

/// HiGHS: HighsSort.cpp::maxheapsort (lines 136-139, 4 lines)
/// Overload `maxheapsort(double* heap_v, HighsInt* heap_i, HighsInt n)` (the only one used).
pub fn maxheapsort(heap_v: &mut [f64], heap_i: &mut [HighsInt], n: HighsInt) {
    build_maxheap(heap_v, heap_i, n);
    max_heapsort(heap_v, heap_i, n);
}

/// HiGHS: HighsSort.cpp::buildMaxheap (lines 155-160, 6 lines)
/// Overload `(double*, HighsInt*, HighsInt)`.
pub fn build_maxheap(heap_v: &mut [f64], heap_i: &mut [HighsInt], n: HighsInt) {
    let mut i: HighsInt = n / 2;
    while i >= 1 {
        max_heapify(heap_v, heap_i, i, n);
        i -= 1;
    }
}

/// HiGHS: HighsSort.cpp::maxHeapsort (lines 187-199, 13 lines)
/// Overload `(double*, HighsInt*, HighsInt)`.
pub fn max_heapsort(heap_v: &mut [f64], heap_i: &mut [HighsInt], n: HighsInt) {
    let mut i: HighsInt = n;
    while i >= 2 {
        // i >= 2, so the two indices are distinct and swap is equivalent to the
        // three-statement exchange in the C++.
        heap_v.swap(i as usize, 1);
        heap_i.swap(i as usize, 1);
        max_heapify(heap_v, heap_i, 1, i - 1);
        i -= 1;
    }
}

/// HiGHS: HighsSort.cpp::maxHeapify (lines 240-259, 20 lines)
/// Overload `(double*, HighsInt*, HighsInt, HighsInt)`.
pub fn max_heapify(heap_v: &mut [f64], heap_i: &mut [HighsInt], i: HighsInt, n: HighsInt) {
    let temp_v: f64 = heap_v[i as usize];
    let temp_i: HighsInt = heap_i[i as usize];
    let mut j: HighsInt = 2 * i;
    while j <= n {
        if j < n && heap_v[(j + 1) as usize] > heap_v[j as usize] {
            j += 1;
        }
        if temp_v > heap_v[j as usize] {
            break;
        } else if temp_v <= heap_v[j as usize] {
            heap_v[(j / 2) as usize] = heap_v[j as usize];
            heap_i[(j / 2) as usize] = heap_i[j as usize];
            j *= 2;
        }
    }
    heap_v[(j / 2) as usize] = temp_v;
    heap_i[(j / 2) as usize] = temp_i;
}

/// HiGHS: extern/pdqsort/pdqsort.h::pdqsort(begin, end) with the default `std::less`.
///
/// Translate `pdqsort(v.begin(), v.begin() + n)` as `pdqsort(&mut v[..n as usize])`.
/// Every simplex-path call sorts `std::pair`s whose `second` (or `first`, for
/// `HEkkDualRow::workData`) is a distinct index, or plain distinct indices, so the
/// lexicographic order is total and any correct sort yields pdqsort's permutation:
/// pdqsort itself is not ported. Rust tuple `PartialOrd` is lexicographic like
/// `std::pair::operator<` (and, like it, treats `-0.0 == 0.0`).
pub fn pdqsort<T: PartialOrd>(v: &mut [T]) {
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
}
