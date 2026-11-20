//! HiGHS source: extern/pdqsort/pdqsort.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//! Altered source version of pdqsort (Copyright (c) 2021 Orson Peters, zlib License,
//! vendored by HiGHS); the zlib notice is in THIRD_PARTY_NOTICES.md.
//!
//! Pattern-defeating quicksort, the *branchy* instantiation
//! (`pdqsort_loop<Iter, Compare, false>`, [`pdqsort_by`]) and the *branchless*
//! instantiation (`pdqsort_loop<Iter, Compare, true>`, [`pdqsort_branchless_by`]).
//! HPresolve sorts with custom
//! comparators that can compare equal (e.g. singletonColStuffing's
//! `cost/val` ratio, HPresolve.cpp 5008), so the order of tied elements is
//! decided by this exact algorithm: `slice::sort_by` would be unfaithful.
//!
//! Iterators are indices into one slice `v`: C++ `Iter begin` -> `begin:
//! usize`, `*it` -> `v[it]`, `std::iter_swap(a, b)` -> `v.swap(a, b)`,
//! `T tmp = std::move(*sift)` -> `let tmp = v[sift].clone()`. `comp` is
//! called with references in the C++ argument order.
//!
//! `pdqsort(begin, end)` without a comparator (default `std::less` on an
//! arithmetic type) selects the branchless partition in C++. Its output is
//! the unique sorted sequence because equal integers are indistinguishable,
//! so call sites translate it as the existing default-compare
//! `util::highs_sort::pdqsort(&mut v)` (or `v.sort_unstable()` on integers).
//! The explicit-comparator branchless entry point `pdqsort_branchless(begin,
//! end, comp)` is ported as [`pdqsort_branchless_by`].
//!
//! The comparator entry point is named [`pdqsort_by`] (not `pdqsort`)
//! because `util::highs_sort::pdqsort` (default compare) is already in the
//! prelude: `pdqsort(v.begin(), v.end(), comp)` -> `pdqsort_by(&mut v, comp)`.

use crate::prelude::*;

/// `pdqsort_detail::insertion_sort_threshold`.
pub const PDQSORT_INSERTION_SORT_THRESHOLD: usize = 24;
/// `pdqsort_detail::ninther_threshold`.
pub const PDQSORT_NINTHER_THRESHOLD: usize = 128;
/// `pdqsort_detail::partial_insertion_sort_limit`.
pub const PDQSORT_PARTIAL_INSERTION_SORT_LIMIT: usize = 8;
/// `pdqsort_detail::block_size` (branchless partition offset block).
pub const PDQSORT_BLOCK_SIZE: usize = 64;
/// `pdqsort_detail::cacheline_size` (only pads the offset buffers in C++;
/// the Rust buffers are index-addressed, so it sizes them only).
pub const PDQSORT_CACHELINE_SIZE: usize = 64;

/// HiGHS: pdqsort.h::pdqsort_branchless(Iter, Iter, Compare) (lines 517-521, 5 lines), named `pdqsort_branchless_by`
///
/// Entry point: `if (begin == end) return;` then
/// `pdqsort_loop<..., true>(v, 0, v.len(), comp, pdqsort_log2(v.len()))`
/// (`leftmost` defaults to true in C++).
pub fn pdqsort_branchless_by<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    if v.is_empty() {
        return;
    }
    let n = v.len();
    pdqsort_loop_impl(v, 0, n, &mut comp, pdqsort_log2(n), true, true);
}

/// HiGHS: pdqsort.h::pdqsort(Iter, Iter, Compare) (lines 495-508, 14 lines), named `pdqsort_by`
///
/// Entry point: `if (begin == end) return;` then
/// `pdqsort_loop(v, 0, v.len(), comp, pdqsort_log2(v.len()), true)`.
pub fn pdqsort_by<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    if v.is_empty() {
        return;
    }
    let n = v.len();
    pdqsort_loop(v, 0, n, &mut comp, pdqsort_log2(n), true);
}

/// HiGHS: pdqsort.h::pdqsort_detail::log2 (lines 68-73, 6 lines)
pub fn pdqsort_log2(n: usize) -> i32 {
    let mut n = n;
    let mut lg = 0;
    loop {
        n >>= 1;
        if n == 0 {
            break;
        }
        lg += 1;
    }
    lg
}

/// HiGHS: pdqsort.h::pdqsort_detail::insertion_sort (lines 76-95, 20 lines)
pub fn pdqsort_insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) {
    if begin == end {
        return;
    }

    let mut cur = begin + 1;
    while cur != end {
        let mut sift = cur;
        let mut sift_1 = cur - 1;

        // Compare first so we can avoid 2 moves for an element already positioned correctly.
        if comp(&v[sift], &v[sift_1]) {
            let tmp = v[sift].clone();

            loop {
                v[sift] = v[sift_1].clone();
                sift -= 1;
                if sift == begin {
                    break;
                }
                sift_1 -= 1;
                if !comp(&tmp, &v[sift_1]) {
                    break;
                }
            }

            v[sift] = tmp;
        }

        cur += 1;
    }
}

/// HiGHS: pdqsort.h::pdqsort_detail::unguarded_insertion_sort (lines 99-118, 20 lines)
///
/// Reads `v[begin - 1]` as the sentinel (never written).
pub fn pdqsort_unguarded_insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) {
    if begin == end {
        return;
    }

    let mut cur = begin + 1;
    while cur != end {
        let mut sift = cur;
        let mut sift_1 = cur - 1;

        // Compare first so we can avoid 2 moves for an element already positioned correctly.
        if comp(&v[sift], &v[sift_1]) {
            let tmp = v[sift].clone();

            loop {
                v[sift] = v[sift_1].clone();
                sift -= 1;
                sift_1 -= 1;
                if !comp(&tmp, &v[sift_1]) {
                    break;
                }
            }

            v[sift] = tmp;
        }

        cur += 1;
    }
}

/// HiGHS: pdqsort.h::pdqsort_detail::partial_insertion_sort (lines 123-148, 26 lines)
pub fn pdqsort_partial_insertion_sort<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) -> bool {
    if begin == end {
        return true;
    }

    let mut limit: usize = 0;
    let mut cur = begin + 1;
    while cur != end {
        let mut sift = cur;
        let mut sift_1 = cur - 1;

        // Compare first so we can avoid 2 moves for an element already positioned correctly.
        if comp(&v[sift], &v[sift_1]) {
            let tmp = v[sift].clone();

            loop {
                v[sift] = v[sift_1].clone();
                sift -= 1;
                if sift == begin {
                    break;
                }
                sift_1 -= 1;
                if !comp(&tmp, &v[sift_1]) {
                    break;
                }
            }

            v[sift] = tmp;
            limit += cur - sift;
        }

        if limit > PDQSORT_PARTIAL_INSERTION_SORT_LIMIT {
            return false;
        }

        cur += 1;
    }

    true
}

/// HiGHS: pdqsort.h::pdqsort_detail::sort2 (lines 150-153, 4 lines)
pub fn pdqsort_sort2<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    a: usize,
    b: usize,
    comp: &mut F,
) {
    if comp(&v[b], &v[a]) {
        v.swap(a, b);
    }
}

/// HiGHS: pdqsort.h::pdqsort_detail::sort3 (lines 156-161, 6 lines)
pub fn pdqsort_sort3<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    a: usize,
    b: usize,
    cc: usize,
    comp: &mut F,
) {
    pdqsort_sort2(v, a, b, comp);
    pdqsort_sort2(v, b, cc, comp);
    pdqsort_sort2(v, a, b, comp);
}

/// HiGHS: pdqsort.h::pdqsort_detail::partition_right (lines 328-366, 39 lines)
///
/// Returns `(pivot_pos, already_partitioned)`. The pivot is held in a
/// local clone (`T pivot(std::move(*begin))`).
pub fn pdqsort_partition_right<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) -> (usize, bool) {
    // Move pivot into local for speed.
    let pivot = v[begin].clone();

    let mut first = begin;
    let mut last = end;

    // Find the first element greater than or equal than the pivot (the median of 3 guarantees
    // this exists).
    loop {
        first += 1;
        if !comp(&v[first], &pivot) {
            break;
        }
    }

    // Find the first element strictly smaller than the pivot. We have to guard this search if
    // there was no element before *first.
    if first - 1 == begin {
        while first < last {
            last -= 1;
            if comp(&v[last], &pivot) {
                break;
            }
        }
    } else {
        loop {
            last -= 1;
            if comp(&v[last], &pivot) {
                break;
            }
        }
    }

    // If the first pair of elements that should be swapped to partition are the same element,
    // the passed in sequence already was correctly partitioned.
    let already_partitioned = first >= last;

    // Keep swapping pairs of elements that are on the wrong side of the pivot. Previously
    // swapped pairs guard the searches, which is why the first iteration is special-cased
    // above.
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if !comp(&v[first], &pivot) {
                break;
            }
        }
        loop {
            last -= 1;
            if comp(&v[last], &pivot) {
                break;
            }
        }
    }

    // Put the pivot in the right place.
    let pivot_pos = first - 1;
    v[begin] = v[pivot_pos].clone();
    v[pivot_pos] = pivot;

    (pivot_pos, already_partitioned)
}

/// HiGHS: pdqsort.h::pdqsort_detail::swap_offsets (lines 175-195, 21 lines)
///
/// `first`/`last` are the iterator bases (`offsets_l_base`, `offsets_r_base`);
/// `offsets_l`/`offsets_r` are the already-advanced offset blocks
/// (`offsets_l + start_l`, `offsets_r + start_r`).
pub fn pdqsort_swap_offsets<T: Clone>(
    v: &mut [T],
    first: usize,
    last: usize,
    offsets_l: &[u8],
    offsets_r: &[u8],
    num: usize,
    use_swaps: bool,
) {
    if use_swaps {
        // This case is needed for the descending distribution, where we need
        // to have proper swapping for pdqsort to remain O(n).
        for i in 0..num {
            v.swap(first + offsets_l[i] as usize, last - offsets_r[i] as usize);
        }
    } else if num > 0 {
        let mut l = first + offsets_l[0] as usize;
        let mut r = last - offsets_r[0] as usize;
        let tmp = v[l].clone();
        v[l] = v[r].clone();
        for i in 1..num {
            l = first + offsets_l[i] as usize;
            v[r] = v[l].clone();
            r = last - offsets_r[i] as usize;
            v[l] = v[r].clone();
        }
        v[r] = tmp;
    }
}

/// HiGHS: pdqsort.h::pdqsort_detail::partition_right_branchless (lines 200-295, 96 lines)
///
/// Same contract as [`pdqsort_partition_right`]: returns
/// `(pivot_pos, already_partitioned)`. The offset blocks are `u8` arrays
/// indexed from 0 (C++ `align_cacheline` only shifts the base address, which
/// has no effect on values).
pub fn pdqsort_partition_right_branchless<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) -> (usize, bool) {
    // Move pivot into local for speed.
    let pivot = v[begin].clone();
    let mut first = begin;
    let mut last = end;

    // Find the first element greater than or equal than the pivot (the median of 3 guarantees
    // this exists).
    loop {
        first += 1;
        if !comp(&v[first], &pivot) {
            break;
        }
    }

    // Find the first element strictly smaller than the pivot. We have to guard this search if
    // there was no element before *first.
    if first - 1 == begin {
        while first < last {
            last -= 1;
            if comp(&v[last], &pivot) {
                break;
            }
        }
    } else {
        loop {
            last -= 1;
            if comp(&v[last], &pivot) {
                break;
            }
        }
    }

    // If the first pair of elements that should be swapped to partition are the same element,
    // the passed in sequence already was correctly partitioned.
    let already_partitioned = first >= last;
    if !already_partitioned {
        v.swap(first, last);
        first += 1;

        // The following branchless partitioning is derived from "BlockQuicksort: How Branch
        // Mispredictions don't affect Quicksort" by Stefan Edelkamp and Armin Weiss, but
        // heavily micro-optimized.
        let mut offsets_l = [0u8; PDQSORT_BLOCK_SIZE + PDQSORT_CACHELINE_SIZE];
        let mut offsets_r = [0u8; PDQSORT_BLOCK_SIZE + PDQSORT_CACHELINE_SIZE];

        let mut offsets_l_base = first;
        let mut offsets_r_base = last;
        let mut num_l: usize = 0;
        let mut num_r: usize = 0;
        let mut start_l: usize = 0;
        let mut start_r: usize = 0;

        while first < last {
            // Fill up offset blocks with elements that are on the wrong side.
            // First we determine how much elements are considered for each offset block.
            let num_unknown = last - first;
            let left_split = if num_l == 0 {
                if num_r == 0 {
                    num_unknown / 2
                } else {
                    num_unknown
                }
            } else {
                0
            };
            let right_split = if num_r == 0 {
                num_unknown - left_split
            } else {
                0
            };

            // Fill the offset blocks. The C++ unrolled-by-8 loop over
            // `block_size` is the same as a plain loop of `min(split, block_size)`.
            let left_count = left_split.min(PDQSORT_BLOCK_SIZE);
            for i in 0..left_count {
                offsets_l[num_l] = i as u8;
                num_l += usize::from(!comp(&v[first], &pivot));
                first += 1;
            }

            let right_count = right_split.min(PDQSORT_BLOCK_SIZE);
            for i in 0..right_count {
                offsets_r[num_r] = (i + 1) as u8;
                last -= 1;
                num_r += usize::from(comp(&v[last], &pivot));
            }

            // Swap elements and update block sizes and first/last boundaries.
            let num = num_l.min(num_r);
            pdqsort_swap_offsets(
                v,
                offsets_l_base,
                offsets_r_base,
                &offsets_l[start_l..],
                &offsets_r[start_r..],
                num,
                num_l == num_r,
            );
            num_l -= num;
            num_r -= num;
            start_l += num;
            start_r += num;

            if num_l == 0 {
                start_l = 0;
                offsets_l_base = first;
            }

            if num_r == 0 {
                start_r = 0;
                offsets_r_base = last;
            }
        }

        // We have now fully identified [first, last)'s proper position. Swap the last elements.
        if num_l > 0 {
            // C++: `offsets_l += start_l; while (num_l--) iter_swap(base + offsets_l[num_l], --last);`
            while num_l > 0 {
                num_l -= 1;
                last -= 1;
                v.swap(offsets_l_base + offsets_l[start_l + num_l] as usize, last);
            }
            first = last;
        }
        if num_r > 0 {
            // C++: `offsets_r += start_r; while (num_r--) iter_swap(base - offsets_r[num_r], first), ++first;`
            while num_r > 0 {
                num_r -= 1;
                v.swap(offsets_r_base - offsets_r[start_r + num_r] as usize, first);
                first += 1;
            }
            // C++: `last = first;` is a dead store here (`last` is not read again).
        }
    }

    // Put the pivot in the right place.
    let pivot_pos = first - 1;
    v[begin] = v[pivot_pos].clone();
    v[pivot_pos] = pivot;

    (pivot_pos, already_partitioned)
}

/// HiGHS: pdqsort.h::pdqsort_detail::partition_left (lines 372-396, 25 lines)
pub fn pdqsort_partition_left<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) -> usize {
    // Move pivot into local for speed.
    let pivot = v[begin].clone();
    let mut first = begin;
    let mut last = end;

    loop {
        last -= 1;
        if !comp(&pivot, &v[last]) {
            break;
        }
    }

    if last + 1 == end {
        while first < last {
            first += 1;
            if comp(&pivot, &v[first]) {
                break;
            }
        }
    } else {
        loop {
            first += 1;
            if comp(&pivot, &v[first]) {
                break;
            }
        }
    }

    while first < last {
        v.swap(first, last);
        loop {
            last -= 1;
            if !comp(&pivot, &v[last]) {
                break;
            }
        }
        loop {
            first += 1;
            if comp(&pivot, &v[first]) {
                break;
            }
        }
    }

    let pivot_pos = last;
    v[begin] = v[pivot_pos].clone();
    v[pivot_pos] = pivot;

    pivot_pos
}

/// HiGHS: pdqsort.h::pdqsort_detail::pdqsort_loop (lines 399-491, 93 lines)
///
/// `Branchless = false`: always `partition_right`. `diff_t` arithmetic on
/// indices: sizes are `usize` and never negative here. The heapsort
/// fallback (`std::make_heap` + `std::sort_heap`) is
/// [`pdqsort_heap_fallback`].
///
/// C++ default arguments: `leftmost = true`.
pub fn pdqsort_loop<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
    bad_allowed: i32,
    leftmost: bool,
) {
    pdqsort_loop_impl(v, begin, end, comp, bad_allowed, leftmost, false);
}

/// PORT NOTE: the C++ `bool Branchless` template parameter becomes the runtime
/// flag `branchless`; both instantiations are otherwise identical code.
fn pdqsort_loop_impl<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
    bad_allowed: i32,
    leftmost: bool,
    branchless: bool,
) {
    let mut begin = begin;
    let mut end = end;
    let mut bad_allowed = bad_allowed;
    let mut leftmost = leftmost;

    // Use a while loop for tail recursion elimination.
    loop {
        let size = end - begin;

        // Insertion sort is faster for small arrays.
        if size < PDQSORT_INSERTION_SORT_THRESHOLD {
            if leftmost {
                pdqsort_insertion_sort(v, begin, end, comp);
            } else {
                pdqsort_unguarded_insertion_sort(v, begin, end, comp);
            }
            return;
        }

        // Choose pivot as median of 3 or pseudomedian of 9.
        let s2 = size / 2;
        if size > PDQSORT_NINTHER_THRESHOLD {
            pdqsort_sort3(v, begin, begin + s2, end - 1, comp);
            pdqsort_sort3(v, begin + 1, begin + (s2 - 1), end - 2, comp);
            pdqsort_sort3(v, begin + 2, begin + (s2 + 1), end - 3, comp);
            pdqsort_sort3(v, begin + (s2 - 1), begin + s2, begin + (s2 + 1), comp);
            v.swap(begin, begin + s2);
        } else {
            pdqsort_sort3(v, begin + s2, begin, end - 1, comp);
        }

        // If *(begin - 1) is the end of the right partition of a previous partition operation
        // there is no element in [begin, end) that is smaller than *(begin - 1). Then if our
        // pivot compares equal to *(begin - 1) we change strategy, putting equal elements in
        // the left partition, greater elements in the right partition. We do not have to
        // recurse on the left partition, since it's sorted (all equal).
        if !leftmost && !comp(&v[begin - 1], &v[begin]) {
            begin = pdqsort_partition_left(v, begin, end, comp) + 1;
            continue;
        }

        // Partition and get results.
        let (pivot_pos, already_partitioned) = if branchless {
            pdqsort_partition_right_branchless(v, begin, end, comp)
        } else {
            pdqsort_partition_right(v, begin, end, comp)
        };

        // Check for a highly unbalanced partition.
        let l_size = pivot_pos - begin;
        let r_size = end - (pivot_pos + 1);
        let highly_unbalanced = l_size < size / 8 || r_size < size / 8;

        // If we got a highly unbalanced partition we shuffle elements to break many patterns.
        if highly_unbalanced {
            // If we had too many bad partitions, switch to heapsort to guarantee O(n log n).
            bad_allowed -= 1;
            if bad_allowed == 0 {
                pdqsort_heap_fallback(v, begin, end, comp);
                return;
            }

            if l_size >= PDQSORT_INSERTION_SORT_THRESHOLD {
                v.swap(begin, begin + l_size / 4);
                v.swap(pivot_pos - 1, pivot_pos - l_size / 4);

                if l_size > PDQSORT_NINTHER_THRESHOLD {
                    v.swap(begin + 1, begin + (l_size / 4 + 1));
                    v.swap(begin + 2, begin + (l_size / 4 + 2));
                    v.swap(pivot_pos - 2, pivot_pos - (l_size / 4 + 1));
                    v.swap(pivot_pos - 3, pivot_pos - (l_size / 4 + 2));
                }
            }

            if r_size >= PDQSORT_INSERTION_SORT_THRESHOLD {
                v.swap(pivot_pos + 1, pivot_pos + (1 + r_size / 4));
                v.swap(end - 1, end - r_size / 4);

                if r_size > PDQSORT_NINTHER_THRESHOLD {
                    v.swap(pivot_pos + 2, pivot_pos + (2 + r_size / 4));
                    v.swap(pivot_pos + 3, pivot_pos + (3 + r_size / 4));
                    v.swap(end - 2, end - (1 + r_size / 4));
                    v.swap(end - 3, end - (2 + r_size / 4));
                }
            }
        } else {
            // If we were decently balanced and we tried to sort an already partitioned
            // sequence try to use insertion sort.
            if already_partitioned
                && pdqsort_partial_insertion_sort(v, begin, pivot_pos, comp)
                && pdqsort_partial_insertion_sort(v, pivot_pos + 1, end, comp)
            {
                return;
            }
        }

        // Sort the left partition first using recursion and do tail recursion elimination for
        // the right-hand partition.
        pdqsort_loop_impl(v, begin, pivot_pos, comp, bad_allowed, leftmost, branchless);
        begin = pivot_pos + 1;
        leftmost = false;
    }
}

/// libc++ `std::make_heap(begin, end, comp); std::sort_heap(begin, end, comp);`
/// (pdqsort.h lines 445-447).
///
/// The libc++ source is not in `reference/`, so this cannot be checked
/// line by line (unfaithful spot). Implement the
/// textbook max-heap with `comp` as "less": make_heap sifts down every
/// parent from `(n - 2) / 2` to 0; sort_heap repeatedly swaps the root to
/// the back and sifts down. Only reached after `log2(n)` highly unbalanced
/// partitions (adversarial inputs).
pub fn pdqsort_heap_fallback<T: Clone, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    begin: usize,
    end: usize,
    comp: &mut F,
) {
    // Sift the element at (relative) index `root` down within [begin, begin + len),
    // treating `comp` as "less" (max-heap).
    let sift_down = |v: &mut [T], root: usize, len: usize, comp: &mut F| {
        let mut parent = root;
        let mut child = 2 * parent + 1;
        while child < len {
            if child + 1 < len && comp(&v[begin + child], &v[begin + child + 1]) {
                child += 1;
            }
            if !comp(&v[begin + parent], &v[begin + child]) {
                break;
            }
            v.swap(begin + parent, begin + child);
            parent = child;
            child = 2 * parent + 1;
        }
    };

    let n = end - begin;

    // std::make_heap: sift down every parent from (n - 2) / 2 to 0.
    if n > 1 {
        let mut parent = (n - 2) / 2;
        loop {
            sift_down(v, parent, n, comp);
            if parent == 0 {
                break;
            }
            parent -= 1;
        }
    }

    // std::sort_heap: repeatedly swap the root to the back and sift down.
    let mut len = n;
    while len > 1 {
        len -= 1;
        v.swap(begin, begin + len);
        sift_down(v, 0, len, comp);
    }
}
