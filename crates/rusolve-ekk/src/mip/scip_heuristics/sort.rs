//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIPsortRealPtr (sorttpl.c instantiation RealPtr) for the heuristics.

use super::*;

/// SCIP 11.0.0 (32dd425): sorttpl.c (line 50) — minimum input size to use ninther for pivot selection
pub(crate) const S6_SORTTPL_MINSIZENINTHER: usize = 729;

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 157-158, SORTTPL_ISBETTER / SORTTPL_ISWORSE with SORTTPL_CMP (x)-(y))
#[inline]
fn is_better(x: f64, y: f64) -> bool {
    (x - y) < 0.0
}

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 157-158, SORTTPL_ISWORSE with SORTTPL_CMP (x)-(y))
#[inline]
fn is_worse(x: f64, y: f64) -> bool {
    (x - y) > 0.0
}

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 171-243)
fn shell_sort(keys: &mut [f64], ptrs: &mut [usize], start: isize, end: isize) {
    debug_assert!(start <= end);

    // static const int incs[3] = {1, 5, 19};
    let incs: [isize; 3] = [1, 5, 19];

    for k in (0..=2).rev() {
        let h = incs[k];
        let first = h + start;

        for i in first..=end {
            let tempkey = keys[i as usize];
            let tempptr = ptrs[i as usize];

            let mut j = i;
            while j >= first && is_better(tempkey, keys[(j - h) as usize]) {
                keys[j as usize] = keys[(j - h) as usize];
                ptrs[j as usize] = ptrs[(j - h) as usize];
                j -= h;
            }

            keys[j as usize] = tempkey;
            ptrs[j as usize] = tempptr;
        }
    }
}

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 245-296)
fn median_three(keys: &[f64], ia: isize, ib: isize, ic: isize) -> isize {
    debug_assert!(ia >= 0);
    debug_assert!(ib >= 0);
    debug_assert!(ic >= 0);
    debug_assert!(ia != ib);
    debug_assert!(ib != ic);
    debug_assert!(ic != ia);

    if is_better(keys[ia as usize], keys[ib as usize]) {
        if is_better(keys[ib as usize], keys[ic as usize]) {
            ib
        } else if is_better(keys[ia as usize], keys[ic as usize]) {
            ic
        } else {
            ia
        }
    } else if is_better(keys[ib as usize], keys[ic as usize]) {
        if is_better(keys[ia as usize], keys[ic as usize]) {
            ia
        } else {
            ic
        }
    } else {
        ib
    }
}

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 298-368)
fn select_pivot_index(keys: &[f64], start: isize, end: isize) -> isize {
    // use the middle index on small arrays
    if end - start < S6_SORTTPL_SHELLSORTMAX as isize {
        (start + end) / 2
    } else if end - start + 1 < S6_SORTTPL_MINSIZENINTHER as isize {
        // select the median of the first, last, and middle element as pivot element
        let mid = (start + end) / 2;
        median_three(keys, start, mid, end)
    } else {
        // use the median of medians of nine evenly distributed elements of the key array
        let gap = (end - start + 1) / 9;

        debug_assert!(start + 8 * gap <= end);

        // collect 3 medians evenly distributed over the array
        let median1 = median_three(keys, start, start + gap, start + 2 * gap);
        let median2 = median_three(keys, start + 3 * gap, start + 4 * gap, start + 5 * gap);
        let median3 = median_three(keys, start + 6 * gap, start + 7 * gap, start + 8 * gap);

        // compute and return the median of the medians
        median_three(keys, median1, median2, median3)
    }
}

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 371-553)
fn q_sort(keys: &mut [f64], ptrs: &mut [usize], start: isize, end: isize, type_: bool) {
    let mut start = start;
    let mut end = end;
    let mut type_ = type_;

    debug_assert!(start <= end);

    // use quick-sort for long lists
    while end - start >= S6_SORTTPL_SHELLSORTMAX as isize {
        // select pivot element
        let mid = select_pivot_index(keys, start, end);
        let pivotkey = keys[mid as usize];

        // partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end]
        let mut lo = start;
        let mut hi = end;
        loop {
            if type_ {
                while lo < end && is_better(keys[lo as usize], pivotkey) {
                    lo += 1;
                }
                while hi > start && !is_better(keys[hi as usize], pivotkey) {
                    hi -= 1;
                }
            } else {
                while lo < end && !is_worse(keys[lo as usize], pivotkey) {
                    lo += 1;
                }
                while hi > start && is_worse(keys[hi as usize], pivotkey) {
                    hi -= 1;
                }
            }

            if lo >= hi {
                break;
            }

            keys.swap(lo as usize, hi as usize);
            ptrs.swap(lo as usize, hi as usize);

            lo += 1;
            hi -= 1;
        }

        debug_assert!((hi == lo - 1) || (type_ && hi == start) || (!type_ && lo == end));

        // skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)
        if type_ {
            while lo < end && !is_better(pivotkey, keys[lo as usize]) {
                lo += 1;
            }

            // make sure that we have at least one element in the smaller partition
            if lo == start {
                // everything is greater or equal than the pivot element: move pivot to the left (degenerate case)
                debug_assert!(!is_better(keys[mid as usize], pivotkey));
                debug_assert!(!is_better(pivotkey, keys[mid as usize]));
                keys.swap(lo as usize, mid as usize);
                ptrs.swap(lo as usize, mid as usize);
                lo += 1;
            }
        } else {
            while hi > start && !is_worse(pivotkey, keys[hi as usize]) {
                hi -= 1;
            }

            // make sure that we have at least one element in the smaller partition
            if hi == end {
                // everything is greater or equal than the pivot element: move pivot to the left (degenerate case)
                debug_assert!(!is_better(keys[mid as usize], pivotkey));
                debug_assert!(!is_better(pivotkey, keys[mid as usize]));
                keys.swap(hi as usize, mid as usize);
                ptrs.swap(hi as usize, mid as usize);
                hi -= 1;
            }
        }

        // sort the smaller partition by a recursive call, sort the larger part without recursion
        if hi - start <= end - lo {
            // sort [start,hi] with a recursive call
            if start < hi {
                q_sort(keys, ptrs, start, hi, !type_);
            }

            // now focus on the larger part [lo,end]
            start = lo;
        } else {
            if lo < end {
                // sort [lo,end] with a recursive call
                q_sort(keys, ptrs, lo, end, !type_);
            }

            // now focus on the larger part [start,hi]
            end = hi;
        }
        type_ = !type_;
    }

    // use shell sort on the remaining small list
    if end > start {
        shell_sort(keys, ptrs, start, end);
    }
}

/// SCIP 11.0.0 (32dd425): sorttpl.c (lines 576-599, SCIPsortRealPtr entry)
///
/// Sorts `keys` ascending and applies the same permutation to `ptrs`
/// (SCIP's `SCIPsortRealPtr(key, ptr, len)`).
pub(crate) fn s6_sort_real_ptr(keys: &mut [f64], ptrs: &mut [usize]) {
    let len = keys.len();

    // ignore the trivial cases
    if len <= 1 {
        return;
    }

    // use shell sort on the small list
    if len <= S6_SORTTPL_SHELLSORTMAX {
        shell_sort(keys, ptrs, 0, (len - 1) as isize);
    } else {
        q_sort(keys, ptrs, 0, (len - 1) as isize, true);
    }
}

// ---------------------------------------------------------------------------
// Generic sorttpl.c (shiftandpropagate needs Real/Int keys, ascending and descending)
// ---------------------------------------------------------------------------

/// sorttpl.c with a caller-supplied `SORTTPL_ISBETTER` (`better(x, y)` is `CMP(x, y) < 0`;
/// `SORTTPL_ISWORSE(x, y)` is `better(y, x)` for the numeric comparators). Descending order
/// is `better = |x, y| y - x < 0`, the template's `SORTTPL_BACKWARDS`.
pub(crate) fn sorttpl_sort<K: Copy, P>(
    keys: &mut [K],
    ptrs: &mut [P],
    better: &dyn Fn(K, K) -> bool,
) {
    let len = keys.len();
    debug_assert_eq!(len, ptrs.len());
    if len <= 1 {
        return;
    }
    if len <= S6_SORTTPL_SHELLSORTMAX {
        g_shell_sort(keys, ptrs, 0, (len - 1) as isize, better);
    } else {
        g_q_sort(keys, ptrs, 0, (len - 1) as isize, true, better);
    }
}

fn g_shell_sort<K: Copy, P>(
    keys: &mut [K],
    ptrs: &mut [P],
    start: isize,
    end: isize,
    better: &dyn Fn(K, K) -> bool,
) {
    let incs: [isize; 3] = [1, 5, 19];
    for k in (0..=2).rev() {
        let h = incs[k];
        let first = h + start;
        for i in first..=end {
            // insertion by rotation: move element i down by steps of h
            let tempkey = keys[i as usize];
            let mut j = i;
            while j >= first && better(tempkey, keys[(j - h) as usize]) {
                keys.swap(j as usize, (j - h) as usize);
                ptrs.swap(j as usize, (j - h) as usize);
                j -= h;
            }
        }
    }
}

fn g_median_three<K: Copy>(
    keys: &[K],
    ia: isize,
    ib: isize,
    ic: isize,
    better: &dyn Fn(K, K) -> bool,
) -> isize {
    let (ka, kb, kc) = (keys[ia as usize], keys[ib as usize], keys[ic as usize]);
    if better(ka, kb) {
        if better(kb, kc) {
            ib
        } else if better(ka, kc) {
            ic
        } else {
            ia
        }
    } else if better(kb, kc) {
        if better(ka, kc) { ia } else { ic }
    } else {
        ib
    }
}

fn g_select_pivot_index<K: Copy>(
    keys: &[K],
    start: isize,
    end: isize,
    better: &dyn Fn(K, K) -> bool,
) -> isize {
    if end - start < S6_SORTTPL_SHELLSORTMAX as isize {
        (start + end) / 2
    } else if end - start + 1 < S6_SORTTPL_MINSIZENINTHER as isize {
        let mid = (start + end) / 2;
        g_median_three(keys, start, mid, end, better)
    } else {
        let gap = (end - start + 1) / 9;
        let m1 = g_median_three(keys, start, start + gap, start + 2 * gap, better);
        let m2 = g_median_three(
            keys,
            start + 3 * gap,
            start + 4 * gap,
            start + 5 * gap,
            better,
        );
        let m3 = g_median_three(
            keys,
            start + 6 * gap,
            start + 7 * gap,
            start + 8 * gap,
            better,
        );
        g_median_three(keys, m1, m2, m3, better)
    }
}

fn g_q_sort<K: Copy, P>(
    keys: &mut [K],
    ptrs: &mut [P],
    start: isize,
    end: isize,
    type_: bool,
    better: &dyn Fn(K, K) -> bool,
) {
    let mut start = start;
    let mut end = end;
    let mut type_ = type_;
    let worse = |x: K, y: K| better(y, x);

    while end - start >= S6_SORTTPL_SHELLSORTMAX as isize {
        let mid = g_select_pivot_index(keys, start, end, better);
        let pivotkey = keys[mid as usize];

        let mut lo = start;
        let mut hi = end;
        loop {
            if type_ {
                while lo < end && better(keys[lo as usize], pivotkey) {
                    lo += 1;
                }
                while hi > start && !better(keys[hi as usize], pivotkey) {
                    hi -= 1;
                }
            } else {
                while lo < end && !worse(keys[lo as usize], pivotkey) {
                    lo += 1;
                }
                while hi > start && worse(keys[hi as usize], pivotkey) {
                    hi -= 1;
                }
            }
            if lo >= hi {
                break;
            }
            keys.swap(lo as usize, hi as usize);
            ptrs.swap(lo as usize, hi as usize);
            lo += 1;
            hi -= 1;
        }

        if type_ {
            while lo < end && !better(pivotkey, keys[lo as usize]) {
                lo += 1;
            }
            if lo == start {
                keys.swap(lo as usize, mid as usize);
                ptrs.swap(lo as usize, mid as usize);
                lo += 1;
            }
        } else {
            while hi > start && !worse(pivotkey, keys[hi as usize]) {
                hi -= 1;
            }
            if hi == end {
                keys.swap(hi as usize, mid as usize);
                ptrs.swap(hi as usize, mid as usize);
                hi -= 1;
            }
        }

        if hi - start <= end - lo {
            if start < hi {
                g_q_sort(keys, ptrs, start, hi, !type_, better);
            }
            start = lo;
        } else {
            if lo < end {
                g_q_sort(keys, ptrs, lo, end, !type_, better);
            }
            end = hi;
        }
        type_ = !type_;
    }

    if end > start {
        g_shell_sort(keys, ptrs, start, end, better);
    }
}

#[cfg(test)]
mod tests {
    // Expected permutations: SCIP 11.0.0 sorttpl.c compiled with gcc (NDEBUG,
    // -ffp-contract=off) as the RealPtr instantiation of misc.c:5794-5798, run
    // on the same inputs; SCIPsortRealPtr output copied verbatim.
    use super::*;

    #[test]
    fn s6_sort_n0_1_2() {
        // n=0
        let mut keys: Vec<f64> = vec![];
        let mut ptrs: Vec<usize> = vec![];
        s6_sort_real_ptr(&mut keys, &mut ptrs);
        assert!(keys.is_empty());
        assert!(ptrs.is_empty());

        // n=1
        let mut keys: Vec<f64> = vec![1.0];
        let mut ptrs: Vec<usize> = vec![0];
        s6_sort_real_ptr(&mut keys, &mut ptrs);
        assert_eq!(keys, vec![1.0]);
        assert_eq!(ptrs, vec![0]);

        // n=2 reversed
        let mut keys: Vec<f64> = vec![2.0, 1.0];
        let mut ptrs: Vec<usize> = vec![0, 1];
        s6_sort_real_ptr(&mut keys, &mut ptrs);
        assert_eq!(keys, vec![1.0, 2.0]);
        assert_eq!(ptrs, vec![1, 0]);

        // n=2 tie
        let mut keys: Vec<f64> = vec![1.0, 1.0];
        let mut ptrs: Vec<usize> = vec![0, 1];
        s6_sort_real_ptr(&mut keys, &mut ptrs);
        assert_eq!(keys, vec![1.0, 1.0]);
        assert_eq!(ptrs, vec![0, 1]);
    }

    #[test]
    fn s6_sort_shell_25_ties() {
        let n = 25;
        let mut keys: Vec<f64> = (0..n).map(|i| ((i * 7) % 5) as f64).collect();
        let mut ptrs: Vec<usize> = (0..n).collect();

        s6_sort_real_ptr(&mut keys, &mut ptrs);

        let expected_keys: Vec<f64> = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0,
            3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0,
        ];
        let expected_ptrs: Vec<usize> = vec![
            0, 20, 5, 10, 15, 3, 23, 8, 13, 18, 21, 6, 11, 16, 1, 9, 14, 19, 4, 24, 7, 12, 17, 2,
            22,
        ];

        assert_eq!(keys, expected_keys);
        assert_eq!(ptrs, expected_ptrs);
    }

    #[test]
    fn s6_sort_quick_40_ties() {
        let n = 40;
        let mut keys: Vec<f64> = (0..n).map(|i| ((i * 13) % 6) as f64).collect();
        let mut ptrs: Vec<usize> = (0..n).collect();

        s6_sort_real_ptr(&mut keys, &mut ptrs);

        let expected_keys: Vec<f64> = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
            2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
            5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
        ];
        let expected_ptrs: Vec<usize> = vec![
            0, 36, 30, 24, 18, 12, 6, 7, 13, 19, 25, 1, 37, 31, 8, 14, 20, 38, 32, 26, 2, 9, 39,
            33, 3, 15, 21, 27, 28, 34, 22, 4, 10, 16, 23, 17, 35, 5, 11, 29,
        ];

        assert_eq!(keys, expected_keys);
        assert_eq!(ptrs, expected_ptrs);
    }

    #[test]
    fn s6_sort_sorted_and_reverse() {
        // n=30 sorted
        let n = 30;
        let mut keys: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let mut ptrs: Vec<usize> = (0..n).collect();
        s6_sort_real_ptr(&mut keys, &mut ptrs);
        assert_eq!(keys, (0..n).map(|i| i as f64).collect::<Vec<f64>>());
        assert_eq!(ptrs, (0..n).collect::<Vec<usize>>());

        // n=30 reversed
        let mut keys: Vec<f64> = (0..n).map(|i| (n - 1 - i) as f64).collect();
        let mut ptrs: Vec<usize> = (0..n).collect();
        s6_sort_real_ptr(&mut keys, &mut ptrs);
        assert_eq!(keys, (0..n).map(|i| i as f64).collect::<Vec<f64>>());
        assert_eq!(ptrs, (0..n).rev().collect::<Vec<usize>>());
    }

    #[test]
    fn s6_sort_keys_ascending_ptrs_follow() {
        // keys non-decreasing, and each ptr still refers to its original key
        let check_consistency = |keys: &[f64], ptrs: &[usize], orig_keys: &[f64]| {
            for i in 1..keys.len() {
                assert!(keys[i - 1] <= keys[i], "Keys not sorted at index {}", i);
            }
            for i in 0..keys.len() {
                assert_eq!(
                    orig_keys[ptrs[i]], keys[i],
                    "Ptr {} does not match key {}",
                    ptrs[i], keys[i]
                );
            }
        };

        // n=40 case
        let n40 = 40;
        let orig_keys_40: Vec<f64> = (0..n40).map(|i| ((i * 13) % 6) as f64).collect();
        let mut keys_40 = orig_keys_40.clone();
        let mut ptrs_40: Vec<usize> = (0..n40).collect();
        s6_sort_real_ptr(&mut keys_40, &mut ptrs_40);
        check_consistency(&keys_40, &ptrs_40, &orig_keys_40);

        // n=30 reversed case
        let n30 = 30;
        let orig_keys_30: Vec<f64> = (0..n30).map(|i| (n30 - 1 - i) as f64).collect();
        let mut keys_30 = orig_keys_30.clone();
        let mut ptrs_30: Vec<usize> = (0..n30).collect();
        s6_sort_real_ptr(&mut keys_30, &mut ptrs_30);
        check_consistency(&keys_30, &ptrs_30, &orig_keys_30);
    }

    #[test]
    fn s6_sort_ninther_1000_ties() {
        // n=1000 >= SORTTPL_MINSIZENINTHER (729): the top-level pivot uses the
        // median of nine (sorttpl.c:328-365). Oracle: same gcc build as above,
        // h = fold(h * 1000003 + ptr) over the sorted ptrs, u64 wrapping.
        let n = 1000;
        let orig: Vec<f64> = (0..n).map(|i| ((i * 37) % 11) as f64).collect();
        let mut keys = orig.clone();
        let mut ptrs: Vec<usize> = (0..n).collect();
        s6_sort_real_ptr(&mut keys, &mut ptrs);

        assert_eq!(
            ptrs[..12],
            [0, 187, 451, 176, 462, 770, 792, 198, 990, 440, 473, 11]
        );
        let h = ptrs
            .iter()
            .fold(0u64, |h, &p| h.wrapping_mul(1000003).wrapping_add(p as u64));
        assert_eq!(h, 11032466779233480904);
        for i in 1..n {
            assert!(keys[i - 1] <= keys[i]);
        }
        for i in 0..n {
            assert_eq!(orig[ptrs[i]], keys[i]);
        }
    }
}

#[cfg(test)]
mod generic_tests {
    use super::*;

    #[test]
    fn generic_matches_real_ptr_and_sorts_down() {
        // same permutation as the Real/usize port on ties and a long list
        let n = 300usize;
        let k0: Vec<f64> = (0..n).map(|i| ((i * 7919) % 13) as f64).collect();
        let mut ka = k0.clone();
        let mut pa: Vec<usize> = (0..n).collect();
        s6_sort_real_ptr(&mut ka, &mut pa);
        let mut kb = k0.clone();
        let mut pb: Vec<usize> = (0..n).collect();
        sorttpl_sort(&mut kb, &mut pb, &|x: f64, y: f64| x - y < 0.0);
        assert_eq!(pa, pb);
        assert_eq!(ka, kb);
        // descending, integer keys
        let mut ki: Vec<i32> = (0..n).map(|i| ((i * 31) % 17) as i32).collect();
        let mut pi: Vec<usize> = (0..n).collect();
        let orig = ki.clone();
        sorttpl_sort(&mut ki, &mut pi, &|x: i32, y: i32| y - x < 0);
        assert!(ki.windows(2).all(|w| w[0] >= w[1]));
        assert!(pi.iter().enumerate().all(|(r, &p)| orig[p] == ki[r]));
    }
}
