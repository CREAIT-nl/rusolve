//! HiGHS source: util/HighsSort.cpp (the `HighsInt`-key heapsort overloads
//! and `sortSetData`)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice. `highs_sort.rs` holds only the
//! `double`-key overloads; `sortSetData` (called by
//! `Highs::changeColsBounds(set)` and `Highs::changeColBoundsInterface`)
//! needs the `HighsInt`-key, `HighsInt`-index overloads, so they live here.
//!
//! The heap arrays are 1-based exactly as in C++: entry 0 is unused, so
//! slices have length `n + 1`. Keep `>` / `<=` as written (ties decide the
//! permutation, and the permutation decides which data entry wins when the
//! set has repeated indices).

use crate::prelude::*;

/// HiGHS: HighsSort.cpp::maxheapsort (lines 131-134, 4 lines)
/// Overload `maxheapsort(HighsInt* heap_v, HighsInt* heap_i, HighsInt n)`.
pub fn maxheapsort_int(heap_v: &mut [HighsInt], heap_i: &mut [HighsInt], n: HighsInt) {
    build_maxheap_int(heap_v, heap_i, n);
    max_heapsort_int(heap_v, heap_i, n);
}

/// HiGHS: HighsSort.cpp::buildMaxheap (lines 148-153, 6 lines)
/// Overload `buildMaxheap(HighsInt* heap_v, HighsInt* heap_i, HighsInt n)`.
pub fn build_maxheap_int(heap_v: &mut [HighsInt], heap_i: &mut [HighsInt], n: HighsInt) {
    let mut i: HighsInt = n / 2;
    while i >= 1 {
        max_heapify_int(heap_v, heap_i, i, n);
        i -= 1;
    }
}

/// HiGHS: HighsSort.cpp::maxHeapsort (lines 173-185, 13 lines)
/// Overload `maxHeapsort(HighsInt* heap_v, HighsInt* heap_i, HighsInt n)`.
pub fn max_heapsort_int(heap_v: &mut [HighsInt], heap_i: &mut [HighsInt], n: HighsInt) {
    let mut temp_v: HighsInt;
    let mut temp_i: HighsInt;
    let mut i: HighsInt = n;
    while i >= 2 {
        temp_v = heap_v[i as usize];
        heap_v[i as usize] = heap_v[1];
        heap_v[1] = temp_v;
        temp_i = heap_i[i as usize];
        heap_i[i as usize] = heap_i[1];
        heap_i[1] = temp_i;
        max_heapify_int(heap_v, heap_i, 1, i - 1);
        i -= 1;
    }
}

/// HiGHS: HighsSort.cpp::maxHeapify (lines 219-238, 20 lines)
/// Overload `maxHeapify(HighsInt* heap_v, HighsInt* heap_i, HighsInt i, HighsInt n)`.
/// Keep the redundant `else if (temp_v <= heap_v[j])` test: it is what C++ does.
pub fn max_heapify_int(heap_v: &mut [HighsInt], heap_i: &mut [HighsInt], i: HighsInt, n: HighsInt) {
    let temp_v: HighsInt = heap_v[i as usize];
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

/// HiGHS: HighsSort.cpp::sortSetData (lines 322-344, 23 lines)
/// Overload with three `double` data arrays.
///
/// `const double* dataK` (may be NULL) -> `Option<&[f64]>`; `double*
/// sorted_dataK` -> `Option<&mut [f64]>`, written iff the matching `dataK` is
/// `Some` (C++ tests only `dataK != NULL`). The two are always distinct
/// buffers in the callers ported here (C++ callers pass the user array as
/// `dataK` and a local copy as `sorted_dataK`).
pub fn sort_set_data(
    num_entries: HighsInt,
    set: &mut [HighsInt],
    data0: Option<&[f64]>,
    data1: Option<&[f64]>,
    data2: Option<&[f64]>,
    sorted_data0: Option<&mut [f64]>,
    sorted_data1: Option<&mut [f64]>,
    sorted_data2: Option<&mut [f64]>,
) {
    if num_entries <= 0 {
        return;
    }
    let n: usize = num_entries as usize;
    let mut sort_set_vec: Vec<HighsInt> = vec![0; 1 + n];
    let mut perm_vec: Vec<HighsInt> = vec![0; 1 + n];

    let sort_set: &mut [HighsInt] = &mut sort_set_vec;
    let perm: &mut [HighsInt] = &mut perm_vec;

    for ix in 0..n {
        sort_set[1 + ix] = set[ix];
        perm[1 + ix] = ix as HighsInt;
    }
    maxheapsort_int(sort_set, perm, num_entries);

    let mut sorted0 = sorted_data0;
    let mut sorted1 = sorted_data1;
    let mut sorted2 = sorted_data2;
    for ix in 0..n {
        set[ix] = sort_set[1 + ix];
        if let Some(d) = data0 {
            sorted0
                .as_deref_mut()
                .expect("sortSetData: sorted_data0 missing")[ix] = d[perm[1 + ix] as usize];
        }
        if let Some(d) = data1 {
            sorted1
                .as_deref_mut()
                .expect("sortSetData: sorted_data1 missing")[ix] = d[perm[1 + ix] as usize];
        }
        if let Some(d) = data2 {
            sorted2
                .as_deref_mut()
                .expect("sortSetData: sorted_data2 missing")[ix] = d[perm[1 + ix] as usize];
        }
    }
}
