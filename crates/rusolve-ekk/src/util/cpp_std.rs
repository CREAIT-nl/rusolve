//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The C++ standard-library algorithms whose result order is
//! implementation-defined. Where HiGHS relies on them and ties can occur,
//! the port must reproduce the library the oracle is built against: libc++
//! (Apple clang 21, macOS). Rust's `BinaryHeap`, `sort_unstable_by` and
//! `select_nth_unstable_by` do not give the same order on ties.
//!
//! Translated from libc++ 210106 (LLVM 21.1.6),
//! the `<__algorithm/*.h>` headers of the macOS SDK
//! (`/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/c++/v1/__algorithm/`):
//! `push_heap.h` (`__sift_up`), `pop_heap.h` (`__pop_heap`), `make_heap.h`,
//! `sift_down.h` (`__sift_down`, `__floyd_sift_down`), `sort_heap.h`,
//! `partial_sort.h`, and `partition.h` (the bidirectional overload, which is
//! the one random-access iterators select). libc++ is Apache-2.0 WITH
//! LLVM-exception; see THIRD_PARTY_NOTICES.md. `comp(a, b)` is C++
//! `comp(a, b)`, i.e. "a is less than b".
//!
//! Translation technique: libc++ moves the sifted value into a register and
//! leaves a hole. Here the held value stays physically at the hole, and each
//! "move the child into the hole" is a `swap`. The comparisons see the same
//! operands in the same order, and the final slot contents are identical, so
//! the resulting order matches libc++ exactly. No `Clone`, no `unsafe`.
//!
//! Plain `std::sort` is not here: HiGHS uses `pdqsort` (util/highs_pdqsort.rs).

use crate::prelude::*;

use std::sync::atomic::{AtomicU8, Ordering};

/// 0 = unread, 1 = libc++ (the default, the macOS oracle), 2 = libstdc++.
static FLAVOR: AtomicU8 = AtomicU8::new(0);

/// True when the heap operations follow libstdc++ (the Linux highspy wheel and
/// any GCC build of HiGHS) instead of libc++. Diagnostic switch: the env var
/// `RUSOLVE_STD_LIBSTDCXX=1` (read once) or `set_libstdcxx`. Default off, so
/// the identity references (libc++ oracle) are unchanged. The two libraries
/// differ on ties: libc++ picks the left child when the children compare
/// equal, libstdc++ the right one.
pub fn libstdcxx() -> bool {
    match FLAVOR.load(Ordering::Relaxed) {
        1 => false,
        2 => true,
        _ => {
            let on = std::env::var("RUSOLVE_STD_LIBSTDCXX").is_ok_and(|v| v == "1");
            FLAVOR.store(if on { 2 } else { 1 }, Ordering::Relaxed);
            on
        }
    }
}

/// Select the library whose tie order the heap operations reproduce.
pub fn set_libstdcxx(on: bool) {
    FLAVOR.store(if on { 2 } else { 1 }, Ordering::Relaxed);
}

/// libstdc++ `bits/stl_heap.h` (GCC 13): `__adjust_heap`, `__pop_heap`,
/// `__make_heap`, `__sort_heap`, `__heap_select`. libstdc++ is GPL-3.0 with the
/// GCC Runtime Library Exception. This re-expresses the control flow of those
/// functions (written from recollection of the algorithm, validated against
/// g++ 13.3 output; no header text is copied), so it follows the same terms
/// (THIRD_PARTY_NOTICES.md, same wording as lm_search.rs).
/// The held value stays physically at the hole, as in the libc++ translation.
mod gnu {
    /// `__adjust_heap(first, hole, len, value)`: `v[hole]` holds `value`.
    fn adjust_heap<T, F: FnMut(&T, &T) -> bool>(
        v: &mut [T],
        comp: &mut F,
        hole_start: usize,
        len: usize,
    ) {
        let top = hole_start;
        let mut hole = hole_start;
        let mut second = hole_start;
        while second < len.saturating_sub(1) / 2 {
            second = 2 * (second + 1);
            if comp(&v[second], &v[second - 1]) {
                second -= 1;
            }
            v.swap(hole, second);
            hole = second;
        }
        if len & 1 == 0 && second == (len - 2) / 2 {
            second = 2 * (second + 1);
            v.swap(hole, second - 1);
            hole = second - 1;
        }
        // __push_heap(first, hole, top, value)
        while hole > top {
            let parent = (hole - 1) / 2;
            if !comp(&v[parent], &v[hole]) {
                break;
            }
            v.swap(parent, hole);
            hole = parent;
        }
    }

    /// `std::pop_heap(first, last)`: `__pop_heap(first, last-1, last-1)`.
    pub fn pop_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
        let n = v.len();
        if n > 1 {
            v.swap(0, n - 1);
            adjust_heap(&mut v[..n - 1], comp, 0, n - 1);
        }
    }

    pub fn make_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
        let len = v.len();
        if len < 2 {
            return;
        }
        let mut parent = (len - 2) / 2;
        loop {
            adjust_heap(v, comp, parent, len);
            if parent == 0 {
                return;
            }
            parent -= 1;
        }
    }

    pub fn sort_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
        let mut n = v.len();
        while n > 1 {
            pop_heap(&mut v[..n], comp);
            n -= 1;
        }
    }

    /// `std::partial_sort`: `__heap_select` then `__sort_heap`.
    pub fn partial_sort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], middle: usize, comp: &mut F) {
        make_heap(&mut v[..middle], comp);
        for i in middle..v.len() {
            if comp(&v[i], &v[0]) {
                v.swap(0, i);
                adjust_heap(&mut v[..middle], comp, 0, middle);
            }
        }
        sort_heap(&mut v[..middle], comp);
    }
}

/// libc++ `__sift_up(first, last, comp, len)` with `len = v.len()`. The
/// element to sift is at index `len - 1` (C++ `*--last`). Its hole moves up
/// while `comp(parent, held)`.
fn sift_up<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
    let len = v.len();
    if len > 1 {
        let mut l = (len - 2) / 2;
        let mut hole = len - 1;
        let mut ptr = l;
        if comp(&v[ptr], &v[hole]) {
            loop {
                v.swap(ptr, hole);
                hole = ptr;
                if l == 0 {
                    break;
                }
                l = (l - 1) / 2;
                ptr = l;
                if !comp(&v[ptr], &v[hole]) {
                    break;
                }
            }
        }
    }
}

/// libc++ `__sift_down(first, comp, len, start)`. `v[0]` is `first`.
fn sift_down<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F, len: usize, start: usize) {
    let mut child = start;
    if len < 2 || (len - 2) / 2 < child {
        return;
    }
    child = 2 * child + 1;
    let mut child_i = child;
    if child + 1 < len && comp(&v[child_i], &v[child_i + 1]) {
        child_i += 1;
        child += 1;
    }
    let mut hole = start;
    if comp(&v[child_i], &v[hole]) {
        return;
    }
    loop {
        v.swap(hole, child_i);
        hole = child_i;
        if (len - 2) / 2 < child {
            break;
        }
        child = 2 * child + 1;
        child_i = child;
        if child + 1 < len && comp(&v[child_i], &v[child_i + 1]) {
            child_i += 1;
            child += 1;
        }
        if comp(&v[child_i], &v[hole]) {
            break;
        }
    }
}

/// libc++ `__floyd_sift_down(first, comp, len)`, `len >= 2`. The hole starts
/// at `v[0]`. It moves down to a leaf along the larger children without
/// comparing against the held value, and returns the leaf index.
fn floyd_sift_down<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F, len: usize) -> usize {
    let mut hole = 0usize;
    let mut child_i = 0usize;
    let mut child = 0usize;
    loop {
        child_i += child + 1;
        child = 2 * child + 1;
        if child + 1 < len && comp(&v[child_i], &v[child_i + 1]) {
            child_i += 1;
            child += 1;
        }
        v.swap(hole, child_i);
        hole = child_i;
        if child > (len - 2) / 2 {
            return hole;
        }
    }
}

/// libc++ `__pop_heap(first, last, comp, len)` with `len = v.len()`.
fn pop_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
    let len = v.len();
    if len > 1 {
        let hole = floyd_sift_down(v, comp, len);
        let last = len - 1;
        if hole != last {
            // C++: *hole = move(*last); ++hole; *last = move(top);
            // sift_up(first, hole, ...). The held top is now at `last`.
            v.swap(hole, last);
            sift_up(&mut v[..hole + 1], comp);
        }
    }
}

/// libc++ `__make_heap(first, last, comp)` with `len = v.len()`.
fn make_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
    let n = v.len();
    if n > 1 {
        let mut start = (n - 2) / 2;
        loop {
            sift_down(v, comp, n, start);
            if start == 0 {
                break;
            }
            start -= 1;
        }
    }
}

/// libc++ `__sort_heap(first, last, comp)`.
fn sort_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], comp: &mut F) {
    let mut n = v.len();
    while n > 1 {
        pop_heap(&mut v[..n], comp);
        n -= 1;
    }
}

/// libc++ `std::push_heap(v.begin(), v.end(), comp)`: sifts up the last element.
pub fn cpp_push_heap_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    // libstdc++ `__push_heap` (value held, parent compared first) is the same
    // sequence of comparisons and moves as libc++ `__sift_up`.
    sift_up(v, &mut comp);
}

/// libc++ `std::pop_heap(v.begin(), v.end(), comp)`: moves the top to the
/// back (Floyd sift-down of the hole to a leaf, then sift-up).
pub fn cpp_pop_heap_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    if libstdcxx() {
        gnu::pop_heap(v, &mut comp);
    } else {
        pop_heap(v, &mut comp);
    }
}

/// libc++ `std::make_heap(v.begin(), v.end(), comp)`.
pub fn cpp_make_heap_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    if libstdcxx() {
        gnu::make_heap(v, &mut comp);
    } else {
        make_heap(v, &mut comp);
    }
}

/// libc++ `std::sort_heap(v.begin(), v.end(), comp)`.
pub fn cpp_sort_heap_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    if libstdcxx() {
        gnu::sort_heap(v, &mut comp);
    } else {
        sort_heap(v, &mut comp);
    }
}

/// libc++ `std::partial_sort(v.begin(), v.begin() + middle, v.end(), comp)`.
pub fn cpp_partial_sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], middle: usize, mut comp: F) {
    if middle == 0 {
        return;
    }
    if libstdcxx() {
        gnu::partial_sort(v, middle, &mut comp);
        return;
    }
    make_heap(&mut v[..middle], &mut comp);
    for i in middle..v.len() {
        if comp(&v[i], &v[0]) {
            v.swap(i, 0);
            sift_down(&mut v[..middle], &mut comp, middle, 0);
        }
    }
    sort_heap(&mut v[..middle], &mut comp);
}

/// libc++ `std::partition(v.begin(), v.end(), pred)` (bidirectional
/// overload: forward scan for a false element, backward scan for a true
/// one, swap). Returns the size of the true part.
pub fn cpp_partition<T, F: FnMut(&T) -> bool>(v: &mut [T], mut pred: F) -> usize {
    let mut first = 0usize;
    let mut last = v.len();
    loop {
        loop {
            if first == last {
                return first;
            }
            if !pred(&v[first]) {
                break;
            }
            first += 1;
        }
        loop {
            last -= 1;
            if first == last {
                return first;
            }
            if pred(&v[last]) {
                break;
            }
        }
        v.swap(first, last);
        first += 1;
    }
}

/// C++ `std::priority_queue<T, std::vector<T>, Compare>`: the container
/// plus libc++ heap operations. The comparator is passed to each call (C++
/// stores it; a stored closure would cost the derives). Use it wherever
/// two queued elements can compare equal; `BinaryHeap<Reverse<T>>` is only
/// equivalent when keys are unique (e.g. the free-slot queues).
#[derive(Clone, Debug, Default)]
pub struct CppPriorityQueue<T> {
    /// C++ `c` (the underlying `std::vector`).
    pub c: Vec<T>,
}

impl<T> CppPriorityQueue<T> {
    /// C++ `push(x)`: `c.push_back(x); push_heap(c, comp)`.
    pub fn push<F: FnMut(&T, &T) -> bool>(&mut self, x: T, comp: F) {
        self.c.push(x);
        cpp_push_heap_by(self.c.as_mut_slice(), comp);
    }

    /// C++ `top()`: `c.front()`.
    pub fn top(&self) -> &T {
        &self.c[0]
    }

    /// C++ `pop()`: `pop_heap(c, comp); c.pop_back()`; returns the popped element.
    pub fn pop<F: FnMut(&T, &T) -> bool>(&mut self, comp: F) -> T {
        cpp_pop_heap_by(self.c.as_mut_slice(), comp);
        self.c
            .pop()
            .expect("CppPriorityQueue::pop on an empty queue")
    }

    /// C++ `empty()`.
    pub fn empty(&self) -> bool {
        self.c.is_empty()
    }

    /// C++ `size()`.
    pub fn size(&self) -> usize {
        self.c.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Keys of the 24-element fixture: the LCG `x = x*1103515245 + 12345`
    /// (u32 wrap), key `(x >> 16) % 4`. Ties are dense, so the order of
    /// equal keys is what the tests check. Element = id; comparison = key.
    fn keys() -> Vec<i32> {
        let mut x: u32 = 12345;
        (0..24)
            .map(|_| {
                x = x.wrapping_mul(1103515245).wrapping_add(12345);
                ((x >> 16) % 4) as i32
            })
            .collect()
    }

    fn lt(keys: &[i32]) -> impl FnMut(&usize, &usize) -> bool + '_ {
        move |a, b| keys[*a] < keys[*b]
    }

    fn base() -> Vec<usize> {
        (0..24).collect()
    }

    #[test]
    fn cpp_std_fixture_keys() {
        // Printed by the libc++ generator (see the commit message).
        assert_eq!(
            keys(),
            vec![
                0, 0, 1, 2, 3, 1, 1, 2, 2, 1, 0, 3, 2, 3, 3, 0, 1, 0, 1, 0, 3, 1, 3, 0
            ]
        );
    }

    #[test]
    fn cpp_std_hand_checked_push_heap() {
        // Values, comp = `<`. Push 3, 1, 4 into a max-heap.
        let mut v: Vec<i32> = Vec::new();
        for x in [3, 1, 4] {
            v.push(x);
            cpp_push_heap_by(&mut v, |a: &i32, b: &i32| a < b);
        }
        assert_eq!(v, vec![4, 1, 3]);
    }

    #[test]
    fn cpp_std_hand_checked_partition() {
        // pred = even. [1,2,3,4]: swap(0,3) gives [4,2,3,1], then the split is 2.
        let mut v = vec![1, 2, 3, 4];
        let split = cpp_partition(&mut v, |x: &i32| x % 2 == 0);
        assert_eq!(split, 2);
        assert_eq!(v, vec![4, 2, 3, 1]);
    }

    #[test]
    fn cpp_std_make_heap_matches_libcxx() {
        let k = keys();
        let mut v = base();
        cpp_make_heap_by(&mut v, lt(&k));
        assert_eq!(
            v,
            vec![
                20, 4, 11, 7, 22, 12, 13, 3, 8, 9, 21, 5, 2, 6, 14, 15, 16, 17, 18, 19, 1, 0, 10,
                23
            ]
        );
    }

    #[test]
    fn cpp_std_push_heap_matches_libcxx() {
        let k = keys();
        let mut v: Vec<usize> = Vec::new();
        for id in base() {
            v.push(id);
            cpp_push_heap_by(&mut v, lt(&k));
        }
        assert_eq!(
            v,
            vec![
                4, 20, 11, 7, 22, 12, 13, 16, 8, 2, 3, 0, 5, 6, 14, 15, 1, 17, 18, 19, 9, 10, 21,
                23
            ]
        );
    }

    #[test]
    fn cpp_std_pop_heap_matches_libcxx() {
        let k = keys();
        let mut v = base();
        cpp_make_heap_by(&mut v, lt(&k));
        let mut popped = Vec::new();
        while !v.is_empty() {
            cpp_pop_heap_by(&mut v, lt(&k));
            popped.push(v.pop().unwrap());
        }
        assert_eq!(
            popped,
            vec![
                20, 4, 22, 11, 13, 14, 7, 3, 8, 12, 16, 18, 9, 21, 5, 2, 6, 15, 0, 23, 17, 10, 19,
                1
            ]
        );
    }

    #[test]
    fn cpp_std_heap_after_pops_matches_libcxx() {
        let k = keys();
        let mut v = base();
        cpp_make_heap_by(&mut v, lt(&k));
        for _ in 0..5 {
            cpp_pop_heap_by(&mut v, lt(&k));
            v.pop();
        }
        assert_eq!(
            v,
            vec![
                14, 7, 12, 3, 9, 5, 6, 16, 8, 19, 21, 23, 2, 1, 0, 15, 10, 17, 18
            ]
        );
    }

    #[test]
    fn cpp_std_sort_heap_matches_libcxx() {
        let k = keys();
        let mut v = base();
        cpp_make_heap_by(&mut v, lt(&k));
        cpp_sort_heap_by(&mut v, lt(&k));
        assert_eq!(
            v,
            vec![
                1, 19, 10, 17, 23, 0, 15, 6, 2, 5, 21, 9, 18, 16, 12, 8, 3, 7, 14, 13, 11, 22, 4,
                20
            ]
        );
    }

    #[test]
    fn cpp_std_partial_sort_matches_libcxx() {
        let k = keys();
        let mut v = base();
        cpp_partial_sort_by(&mut v, 7, lt(&k));
        assert_eq!(
            v,
            vec![
                15, 1, 19, 17, 23, 10, 0, 4, 8, 3, 7, 11, 12, 13, 14, 9, 16, 5, 18, 2, 20, 21, 22,
                6
            ]
        );
    }

    #[test]
    fn cpp_std_partition_matches_libcxx() {
        let k = keys();
        let mut v = base();
        let split = cpp_partition(&mut v, |id: &usize| k[*id] < 2);
        assert_eq!(split, 14);
        assert_eq!(
            v,
            vec![
                0, 1, 2, 23, 21, 5, 6, 19, 18, 9, 10, 17, 16, 15, 14, 13, 12, 11, 8, 7, 20, 4, 22,
                3
            ]
        );
    }

    #[test]
    fn cpp_std_priority_queue_matches_libcxx() {
        let k = keys();
        let mut pq: CppPriorityQueue<usize> = CppPriorityQueue::default();
        for id in 0..12 {
            pq.push(id, lt(&k));
        }
        let mut a = Vec::new();
        for _ in 0..5 {
            a.push(*pq.top());
            pq.pop(lt(&k));
        }
        assert_eq!(a, vec![4, 11, 3, 7, 8]);
        for id in 12..24 {
            pq.push(id, lt(&k));
        }
        let mut b = Vec::new();
        while !pq.empty() {
            b.push(*pq.top());
            pq.pop(lt(&k));
        }
        assert_eq!(
            b,
            vec![
                13, 14, 20, 22, 12, 9, 21, 2, 18, 5, 16, 6, 1, 17, 0, 23, 19, 10, 15
            ]
        );
        assert_eq!(pq.size(), 0);
    }

    /// Edge lengths 0..=5 against a libc++ dump (tests/fixtures/mip/unit/cpp_std.cpp):
    /// every key array in {0,1,2}^n, element = position, comparison on key only,
    /// so each output line is the libc++ tie order.
    #[test]
    fn short_inputs_match_libcxx_dump() {
        let expected = include_str!("../../tests/fixtures/mip/unit/cpp_std.out");
        let line = |op: &str, n: usize, code: usize, extra: usize, v: &[usize]| {
            let mut l = format!("{op} {n} {code} {extra}");
            for id in v {
                l.push_str(&format!(" {id}"));
            }
            l.push('\n');
            l
        };
        let mut got = String::new();
        for n in 0..=5usize {
            let total = 3usize.pow(n as u32);
            for code in 0..total {
                let mut keys = vec![0i32; n];
                let mut rem = code;
                for k in keys.iter_mut() {
                    *k = (rem % 3) as i32;
                    rem /= 3;
                }
                let base: Vec<usize> = (0..n).collect();
                let lt = |a: &usize, b: &usize| keys[*a] < keys[*b];

                let mut v = base.clone();
                cpp_make_heap_by(&mut v, lt);
                got.push_str(&line("make", n, code, 0, &v));

                let mut v = base.clone();
                for i in 1..=n {
                    cpp_push_heap_by(&mut v[..i], lt);
                }
                got.push_str(&line("push", n, code, 0, &v));

                let mut v = base.clone();
                cpp_make_heap_by(&mut v, lt);
                for k in (1..=n).rev() {
                    cpp_pop_heap_by(&mut v[..k], lt);
                }
                got.push_str(&line("pop", n, code, 0, &v));

                let mut v = base.clone();
                cpp_make_heap_by(&mut v, lt);
                cpp_sort_heap_by(&mut v, lt);
                got.push_str(&line("sortheap", n, code, 0, &v));

                for m in 0..=n {
                    let mut v = base.clone();
                    cpp_partial_sort_by(&mut v, m, lt);
                    got.push_str(&line("psort", n, code, m, &v));
                }

                for p in 0..3i32 {
                    let mut v = base.clone();
                    let k = cpp_partition(&mut v, |e| keys[*e] == p);
                    got.push_str(&line("part", n, code, p as usize * 100 + k, &v));
                }

                let gt = |a: &usize, b: &usize| keys[*a] > keys[*b];
                let mut pq: CppPriorityQueue<usize> = CppPriorityQueue::default();
                let mut popped = Vec::new();
                for i in 0..n {
                    pq.push(base[i], gt);
                    if i % 2 == 1 {
                        popped.push(*pq.top());
                        pq.pop(gt);
                    }
                }
                while !pq.empty() {
                    popped.push(*pq.top());
                    pq.pop(gt);
                }
                got.push_str(&line("pq", n, code, 0, &popped));
            }
        }
        for (i, (g, e)) in got.lines().zip(expected.lines()).enumerate() {
            assert_eq!(g, e, "first mismatch at line {}", i + 1);
        }
        assert_eq!(got.lines().count(), expected.lines().count());
    }

    /// Expected values printed by a g++ 13.3 / libstdc++ program running the
    /// same operations on the same 24 keys. The
    /// `gnu` functions are called directly so the test does not touch the
    /// process-wide flavor switch.
    #[test]
    fn libstdcxx_heap_ops_match_gcc13() {
        let k = keys();
        let mut v = base();
        gnu::make_heap(&mut v, &mut lt(&k));
        assert_eq!(
            v,
            vec![
                14, 22, 13, 8, 4, 11, 6, 7, 3, 20, 21, 5, 12, 2, 0, 15, 16, 17, 18, 19, 9, 1, 10,
                23
            ]
        );
        let mut q = v.clone();
        gnu::sort_heap(&mut q, &mut lt(&k));
        assert_eq!(
            q,
            vec![
                17, 1, 15, 10, 23, 0, 19, 16, 18, 9, 21, 5, 2, 6, 7, 3, 8, 12, 20, 4, 22, 11, 13,
                14
            ]
        );
        let mut w = base();
        gnu::make_heap(&mut w, &mut lt(&k));
        gnu::pop_heap(&mut w, &mut lt(&k));
        let n = w.len();
        assert_eq!(
            w,
            vec![
                13, 22, 11, 8, 4, 12, 6, 7, 3, 20, 21, 5, 23, 2, 0, 15, 16, 17, 18, 19, 9, 1, 10,
                14
            ]
        );
        gnu::pop_heap(&mut w[..n - 1], &mut lt(&k));
        assert_eq!(
            w,
            vec![
                11, 22, 12, 8, 4, 5, 6, 7, 3, 20, 21, 10, 23, 2, 0, 15, 16, 17, 18, 19, 9, 1, 13,
                14
            ]
        );
        let mut s = base();
        gnu::partial_sort(&mut s, 7, &mut lt(&k));
        assert_eq!(
            s,
            vec![
                0, 1, 19, 23, 15, 17, 10, 4, 8, 3, 7, 11, 12, 13, 14, 6, 16, 2, 18, 5, 20, 21, 22,
                9
            ]
        );
        let mut s2 = base();
        gnu::partial_sort(&mut s2, 8, &mut lt(&k));
        assert_eq!(
            s2,
            vec![
                0, 1, 23, 10, 17, 19, 15, 9, 4, 7, 3, 11, 12, 13, 14, 8, 16, 6, 18, 2, 20, 21, 22,
                5
            ]
        );
        let mut odd: Vec<usize> = (0..23).collect();
        gnu::make_heap(&mut odd, &mut lt(&k));
        gnu::pop_heap(&mut odd, &mut lt(&k));
        assert_eq!(
            odd,
            vec![
                13, 22, 11, 8, 4, 12, 6, 7, 3, 20, 21, 5, 10, 2, 0, 15, 16, 17, 18, 19, 9, 1, 14
            ]
        );
    }
}
