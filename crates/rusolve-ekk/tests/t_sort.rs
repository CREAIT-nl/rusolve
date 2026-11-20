//! Golden test for manifest row 6, `src/util/highs_sort.rs` (HighsSort.cpp:
//! addToDecreasingHeap, sortDecreasingHeap, maxheapsort, buildMaxheap,
//! maxHeapsort, maxHeapify; the `double* heap_v, HighsInt* heap_i`
//! overloads, which are the ones the simplex calls).
//!
//! There is no oracle fixture for HighsSort (fixture gap G-SORT). The port
//! is compared with `reference`, a line-by-line transcription of those
//! functions from HighsSort.cpp (HiGHS 1.15.1, MIT License, Copyright (c)
//! HiGHS; see THIRD_PARTY_NOTICES.md). What matters is the exact final
//! permutation, including the order of equal keys (the BFRT group order in
//! HEkkDualRow and hyper-sparse CHUZC in HEkkPrimal depend on it), so the
//! inputs are full of ties, signed zeros and already-sorted runs. Arrays are
//! 1-based as in HiGHS: entry 0 is unused (for the decreasing heap,
//! `heap_ix[0] = 1` marks "is a heap").

mod common;
mod port;

use rusolve_ekk::prelude::*;

mod reference {
    //! HighsSort.cpp, transcribed. HighsInt = i32; indices are 1-based.
    pub fn add_to_decreasing_heap(
        n: &mut i32,
        mx_n: i32,
        heap_v: &mut [f64],
        heap_ix: &mut [i32],
        v: f64,
        ix: i32,
    ) {
        let mut cd_p: i32;
        let mut pa_p: i32;
        if *n < mx_n {
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
        heap_ix[0] = 1;
    }

    pub fn sort_decreasing_heap(n: i32, heap_v: &mut [f64], heap_ix: &mut [i32]) {
        let mut fo_p: i32;
        let mut srt_p: i32;
        let mut cd_p: i32;
        let mut pa_p: i32;
        let mut ix: i32;
        let mut v: f64;
        if n <= 1 {
            return;
        }
        if heap_ix[0] != 1 {
            fo_p = n / 2 + 1;
            srt_p = n;
        } else {
            fo_p = 1;
            srt_p = n;
        }
        loop {
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

    pub fn maxheapsort(heap_v: &mut [f64], heap_i: &mut [i32], n: i32) {
        build_maxheap(heap_v, heap_i, n);
        max_heapsort(heap_v, heap_i, n);
    }

    pub fn build_maxheap(heap_v: &mut [f64], heap_i: &mut [i32], n: i32) {
        let mut i = n / 2;
        while i >= 1 {
            max_heapify(heap_v, heap_i, i, n);
            i -= 1;
        }
    }

    pub fn max_heapsort(heap_v: &mut [f64], heap_i: &mut [i32], n: i32) {
        let mut i = n;
        while i >= 2 {
            heap_v.swap(i as usize, 1);
            heap_i.swap(i as usize, 1);
            max_heapify(heap_v, heap_i, 1, i - 1);
            i -= 1;
        }
    }

    pub fn max_heapify(heap_v: &mut [f64], heap_i: &mut [i32], i: i32, n: i32) {
        let temp_v = heap_v[i as usize];
        let temp_i = heap_i[i as usize];
        let mut j = 2 * i;
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
}

/// A deterministic test-local generator (not HighsRandom, so this test
/// does not depend on row 3).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// Input arrays (1-based: slot 0 is a sentinel) of every size 0..=40 and a
/// few larger ones, with values drawn from small pools so ties are common.
fn cases() -> Vec<(String, Vec<f64>)> {
    let mut out = Vec::new();
    let pools: [&[f64]; 4] = [
        &[1.0, 2.0, 3.0],
        &[0.0, -0.0, 1e-9, -1e-9, 0.5],
        &[-3.5, 7.25, 7.25, 1.0 / 3.0, 1e300, -1e300, 4.0],
        &[1.0],
    ];
    let mut g = Lcg(20261030);
    for (p, pool) in pools.iter().enumerate() {
        for n in (0..=40).chain([63, 64, 65, 200]) {
            let mut v = vec![-12345.0];
            for _ in 0..n {
                v.push(pool[(g.next() % pool.len() as u64) as usize]);
            }
            out.push((format!("pool{p} n={n}"), v));
        }
    }
    for n in [1usize, 2, 7, 33] {
        let asc: Vec<f64> = std::iter::once(0.0)
            .chain((1..=n).map(|k| k as f64))
            .collect();
        let desc: Vec<f64> = std::iter::once(0.0)
            .chain((1..=n).rev().map(|k| k as f64))
            .collect();
        out.push((format!("ascending n={n}"), asc));
        out.push((format!("descending n={n}"), desc));
    }
    out
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

fn first_diff<T: PartialEq + std::fmt::Debug>(e: &[T], a: &[T]) -> String {
    match e.iter().zip(a).position(|(x, y)| x != y) {
        Some(k) => format!(
            "first difference at [{k}]: expected {:?}, actual {:?}",
            e[k], a[k]
        ),
        None => format!("lengths differ: expected {}, actual {}", e.len(), a.len()),
    }
}

#[track_caller]
fn same(ctx: &str, ev: &[f64], ei: &[i32], av: &[f64], ai: &[i32]) {
    assert!(
        bits(ev) == bits(av),
        "{ctx}: heap_v differs from HighsSort.cpp, {}\n  expected {ev:?}\n  actual   {av:?}",
        first_diff(&bits(ev), &bits(av))
    );
    assert!(
        ei == ai,
        "{ctx}: heap_i differs from HighsSort.cpp (tie order), {}\n  expected {ei:?}\n  actual   {ai:?}",
        first_diff(ei, ai)
    );
}

#[test]
fn maxheapsort_matches_highs_sort_cpp() {
    for (name, v) in cases() {
        let n = (v.len() - 1) as i32;
        let idx: Vec<i32> = (0..=n).map(|k| 100 + k).collect();
        let (mut ev, mut ei) = (v.clone(), idx.clone());
        reference::maxheapsort(&mut ev, &mut ei, n);
        let (mut av, mut ai) = (v.clone(), idx.clone());
        maxheapsort(&mut av, &mut ai, n);
        same(&format!("maxheapsort({name})"), &ev, &ei, &av, &ai);

        let (mut ev, mut ei) = (v.clone(), idx.clone());
        reference::build_maxheap(&mut ev, &mut ei, n);
        let (mut av, mut ai) = (v.clone(), idx.clone());
        build_maxheap(&mut av, &mut ai, n);
        same(&format!("buildMaxheap({name})"), &ev, &ei, &av, &ai);

        // maxHeapsort on the (reference) heap
        let (mut hv, mut hi) = (ev.clone(), ei.clone());
        reference::max_heapsort(&mut ev, &mut ei, n);
        max_heapsort(&mut hv, &mut hi, n);
        same(&format!("maxHeapsort(heap of {name})"), &ev, &ei, &hv, &hi);

        // maxHeapify at every internal node of the raw input
        for i in 1..=n / 2 {
            let (mut ev, mut ei) = (v.clone(), idx.clone());
            reference::max_heapify(&mut ev, &mut ei, i, n);
            let (mut av, mut ai) = (v.clone(), idx.clone());
            max_heapify(&mut av, &mut ai, i, n);
            same(&format!("maxHeapify({name}, i={i})"), &ev, &ei, &av, &ai);
        }
    }
}

#[test]
fn decreasing_heap_matches_highs_sort_cpp() {
    for (name, v) in cases() {
        let total = (v.len() - 1) as i32;
        for mx_n in [1, 3, 8, total.max(1)] {
            let size = (mx_n + 1) as usize;
            let (mut en, mut ev, mut ei) = (0i32, vec![0.0; size], vec![0i32; size]);
            let (mut an, mut av, mut ai) = (0i32, vec![0.0; size], vec![0i32; size]);
            for (k, &x) in v.iter().enumerate().skip(1) {
                reference::add_to_decreasing_heap(&mut en, mx_n, &mut ev, &mut ei, x, k as i32);
                add_to_decreasing_heap(&mut an, mx_n, &mut av, &mut ai, x, k as i32);
                let ctx = format!("addToDecreasingHeap({name}, mx_n={mx_n}) after value #{k}");
                assert_eq!(en, an, "{ctx}: n");
                same(&ctx, &ev, &ei, &av, &ai);
            }
            reference::sort_decreasing_heap(en, &mut ev, &mut ei);
            sort_decreasing_heap(an, &mut av, &mut ai);
            same(
                &format!("sortDecreasingHeap(heap of {name}, mx_n={mx_n}, n={en})"),
                &ev,
                &ei,
                &av,
                &ai,
            );
        }
        // Unordered input (heap_ix[0] != 1): sortDecreasingHeap heapifies first.
        let n = total;
        let mut ev = v.clone();
        let mut ei: Vec<i32> = (0..=n).collect();
        ei[0] = 0;
        let (mut av, mut ai) = (ev.clone(), ei.clone());
        reference::sort_decreasing_heap(n, &mut ev, &mut ei);
        sort_decreasing_heap(n, &mut av, &mut ai);
        same(
            &format!("sortDecreasingHeap(unordered {name}, n={n})"),
            &ev,
            &ei,
            &av,
            &ai,
        );
    }
}
