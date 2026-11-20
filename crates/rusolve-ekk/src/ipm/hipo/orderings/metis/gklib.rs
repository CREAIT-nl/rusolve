//! HiGHS source: extern/metis/libmetis/gklib.c, extern/metis/GKlib/{gk_mkblas.h,gk_mkmemory.h,gk_mkrandom.h,random.c,error.c}
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! The GKlib template instances libmetis uses (`gklib.c`): `i`-BLAS, the
//! `?malloc` / `?set` / `?copy` family, the `irand*` RNG on top of
//! `gk_rand_r`, and the two GKlib quicksort instances. The RNG is ported
//! operation by operation: `gk_rand_r` is an `unsigned` multiply-add that
//! wraps (`u32::wrapping_*`). The quicksort is a median-of-three quicksort
//! with an explicit stack and a final insertion sort, which is not stable;
//! neither may be replaced by `slice::sort*`. Allocation helpers are complete:
//! `imalloc` returns zeros (the C memory is uninitialised and always written
//! first), `icopy` is `memmove` (callers with overlapping ranges use
//! `copy_within`).
//!
//! License note: the quicksort (`gk_sort`, behind `isortd` / `ikvsorti`) is an
//! independent implementation of the GKlib quicksort semantics, verified
//! against the reference outputs in
//! `tests/t_gksort.rs`; no LGPL code is included.

use super::{idx_t, real_t};

/// HiGHS: gklib_defs.h `GK_MKKEYVALUE_T(ikv_t, idx_t, idx_t)`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ikv_t {
    pub key: idx_t,
    pub val: idx_t,
}

/// HiGHS: gklib_defs.h `GK_MKKEYVALUE_T(rkv_t, real_t, idx_t)`
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct rkv_t {
    pub key: real_t,
    pub val: idx_t,
}

/// HiGHS: gk_mkmemory.h `imalloc` (lines 18-21)
pub fn imalloc(n: usize) -> Vec<idx_t> {
    vec![0; n]
}

/// HiGHS: gk_mkmemory.h `ismalloc` (lines 36-45)
pub fn ismalloc(n: usize, ival: idx_t) -> Vec<idx_t> {
    vec![ival; n]
}

/// HiGHS: gk_mkmemory.h `rmalloc` (lines 18-21)
pub fn rmalloc(n: usize) -> Vec<real_t> {
    vec![0.0; n]
}

/// HiGHS: gk_mkmemory.h `rsmalloc` (lines 36-45)
pub fn rsmalloc(n: usize, ival: real_t) -> Vec<real_t> {
    vec![ival; n]
}

/// HiGHS: gk_mkmemory.h `iset` (lines 51-59)
pub fn iset(n: usize, val: idx_t, x: &mut [idx_t]) {
    x[..n].fill(val);
}

/// HiGHS: gk_mkmemory.h `rset` (lines 51-59)
pub fn rset(n: usize, val: real_t, x: &mut [real_t]) {
    x[..n].fill(val);
}

/// HiGHS: gk_mkmemory.h `icopy` (lines 65-68), `memmove(b, a)`
pub fn icopy(n: usize, a: &[idx_t], b: &mut [idx_t]) {
    b[..n].copy_from_slice(&a[..n]);
}

/// HiGHS: gk_mkmemory.h `rcopy` (lines 65-68), `memmove(b, a)`
pub fn rcopy(n: usize, a: &[real_t], b: &mut [real_t]) {
    b[..n].copy_from_slice(&a[..n]);
}

/// HiGHS: gk_mkblas.h::isum (lines 18-27)
/// `OUTTYPE` = `idx_t`: the sum wraps like the C signed add would not, so
/// keep plain `+` (overflow is a bug in either language).
pub fn isum(n: usize, x: &[idx_t], incx: usize) -> idx_t {
    let mut sum: idx_t = 0;
    let mut xi = 0usize;
    for _i in 0..n {
        sum += x[xi];
        xi += incx;
    }
    sum
}

/// HiGHS: gk_mkblas.h::iaxpy (lines 33-42)
pub fn iaxpy(n: usize, alpha: idx_t, x: &[idx_t], incx: usize, y: &mut [idx_t], incy: usize) {
    let mut xi = 0usize;
    let mut yi = 0usize;
    for _i in 0..n {
        y[yi] += alpha * x[xi];
        xi += incx;
        yi += incy;
    }
}

/// HiGHS: GKlib/random.c::gk_rand_r (lines 16-21)
/// `((*rng_state * 1103515245) + 12345) & 0x7fffffff` in `unsigned`.
pub fn gk_rand_r(rng_state: &mut u32) -> i32 {
    let result = rng_state.wrapping_mul(1103515245).wrapping_add(12345) & 0x7fff_ffff;
    *rng_state = result;
    result as i32
}

/// HiGHS: GKlib/random.c::gk_randint64 (lines 24-29)
pub fn gk_randint64(rng_state: &mut u32) -> u64 {
    let piece_1 = (gk_rand_r(rng_state) as u64) << 32;
    let piece_2 = gk_rand_r(rng_state) as u64;
    piece_1 | piece_2
}

/// HiGHS: GKlib/random.c::gk_randint32 (lines 32-35)
pub fn gk_randint32(rng_state: &mut u32) -> u32 {
    gk_rand_r(rng_state) as u32
}

/// HiGHS: gk_mkrandom.h::irand (lines 27-33)
/// `sizeof(idx_t) <= sizeof(int32_t)` picks `gk_randint32`, else
/// `gk_randint64`; the cast to `idx_t` is a C truncation (`as`).
pub fn irand(rng_state: &mut u32) -> idx_t {
    if std::mem::size_of::<idx_t>() <= std::mem::size_of::<i32>() {
        gk_randint32(rng_state) as idx_t
    } else {
        gk_randint64(rng_state) as idx_t
    }
}

/// HiGHS: gk_mkrandom.h::irandInRange (lines 39-42)
pub fn irandInRange(max: idx_t, rng_state: &mut u32) -> idx_t {
    irand(rng_state) % max
}

/// HiGHS: gk_mkrandom.h::irandArrayPermute (lines 50-81)
pub fn irandArrayPermute(
    n: idx_t,
    p: &mut [idx_t],
    nshuffles: idx_t,
    flag: i32,
    rng_state: &mut u32,
) {
    if flag == 1 {
        for i in 0..n {
            p[i as usize] = i;
        }
    }

    if n < 10 {
        for _i in 0..n {
            let v = irandInRange(n, rng_state);
            let u = irandInRange(n, rng_state);
            // gk_SWAP(p[v], p[u], tmp);
            p.swap(v as usize, u as usize);
        }
    } else {
        for _i in 0..nshuffles {
            let v = irandInRange(n - 3, rng_state);
            let u = irandInRange(n - 3, rng_state);
            // gk_SWAP(p[v+0], p[u+2], tmp);
            p.swap(v as usize, (u + 2) as usize);
            // gk_SWAP(p[v+1], p[u+3], tmp);
            p.swap((v + 1) as usize, (u + 3) as usize);
            // gk_SWAP(p[v+2], p[u+0], tmp);
            p.swap((v + 2) as usize, u as usize);
            // gk_SWAP(p[v+3], p[u+1], tmp);
            p.swap((v + 3) as usize, (u + 1) as usize);
        }
    }
}

/// Independent implementation of the GKlib quicksort semantics: a
/// median-of-three quicksort with an explicit LIFO stack of ranges that leaves
/// partitions of at most `THRESH` elements unsorted, followed by a sentinel
/// insertion sort over the whole array. Not stable; the exact comparison and
/// swap sequence determines the final order of equal keys.
fn gk_sort<T: Copy>(a: &mut [T], n: usize, less: impl Fn(&T, &T) -> bool) {
    const THRESH: isize = 8;
    if n == 0 {
        return;
    }
    let n = n as isize;
    let at = |x: isize| x as usize;

    // Phase 1: partition until every remaining range has at most THRESH elements.
    if n > THRESH {
        let mut lo: isize = 0;
        let mut hi: isize = n - 1;
        let mut stack: Vec<(isize, isize)> = Vec::new();
        loop {
            let mid = lo + (hi - lo) / 2;
            if less(&a[at(mid)], &a[at(lo)]) {
                a.swap(at(mid), at(lo));
            }
            if less(&a[at(hi)], &a[at(mid)]) {
                a.swap(at(mid), at(hi));
                if less(&a[at(mid)], &a[at(lo)]) {
                    a.swap(at(mid), at(lo));
                }
            }

            let mut p = mid;
            let mut i = lo + 1;
            let mut j = hi - 1;
            loop {
                while less(&a[at(i)], &a[at(p)]) {
                    i += 1;
                }
                while less(&a[at(p)], &a[at(j)]) {
                    j -= 1;
                }
                if i < j {
                    a.swap(at(i), at(j));
                    if p == i {
                        p = j;
                    } else if p == j {
                        p = i;
                    }
                    i += 1;
                    j -= 1;
                } else if i == j {
                    i += 1;
                    j -= 1;
                    break;
                }
                if i > j {
                    break;
                }
            }

            let left = j - lo;
            let right = hi - i;
            if left <= THRESH {
                if right <= THRESH {
                    match stack.pop() {
                        Some((l, h)) => {
                            lo = l;
                            hi = h;
                        }
                        None => break,
                    }
                } else {
                    lo = i;
                }
            } else if right <= THRESH {
                hi = j;
            } else if left > right {
                stack.push((lo, j));
                lo = i;
            } else {
                stack.push((i, hi));
                hi = j;
            }
        }
    }

    // Phase 2: put the smallest of the first THRESH+1 elements at a[0] as a
    // sentinel, then insertion-sort the rest.
    let mut c: isize = 0;
    for k in 1..=THRESH.min(n - 1) {
        if less(&a[at(k)], &a[at(c)]) {
            c = k;
        }
    }
    if c != 0 {
        a.swap(0, at(c));
    }
    for r in 2..n {
        let mut t = r - 1;
        while less(&a[at(r)], &a[at(t)]) {
            t -= 1;
        }
        t += 1;
        if t != r {
            a[at(t)..=at(r)].rotate_right(1);
        }
    }
}

/// HiGHS: gklib.c::isortd
/// Decreasing sort of `base[..n]` (`i_gt` as the "less" predicate); an
/// independent implementation of the GKlib quicksort semantics (`gk_sort`).
pub fn isortd(n: usize, base: &mut [idx_t]) {
    gk_sort(base, n, |a, b| a > b);
}

/// HiGHS: gklib.c::ikvsorti
/// Increasing sort on `key` (`ikey_lt`); an independent implementation of the
/// GKlib quicksort semantics (`gk_sort`). Equal keys keep the GKlib order,
/// which the oracle's partition order depends on.
pub fn ikvsorti(n: usize, base: &mut [ikv_t]) {
    gk_sort(base, n, |a, b| a.key < b.key);
}

/// HiGHS: GKlib/error.c::gk_errexit (lines 18-21), `GK_ERREXIT` (gk_macros.h line 36)
/// The C prints only without `NDEBUG` and aborts.
pub fn gk_errexit(msg: &str) -> ! {
    panic!("{}", msg);
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn rand_r_first_values() {
        // seed 42 (kMetisSeed): ((42*1103515245)+12345) & 0x7fffffff
        let mut s: u32 = 42;
        let r = gk_rand_r(&mut s);
        let expect = (42u32.wrapping_mul(1103515245).wrapping_add(12345) & 0x7fff_ffff) as i32;
        assert_eq!(r, expect);
        assert_eq!(s, expect as u32);
    }

    #[test]
    fn qsort_matches_sorted_order() {
        let mut v: Vec<idx_t> = vec![5, 1, 9, 3, 3, 7, 0, 2, 8, 6, 4, 1];
        isortd(v.len(), &mut v);
        assert_eq!(v, vec![9, 8, 7, 6, 5, 4, 3, 3, 2, 1, 1, 0]);
        let mut kv: Vec<ikv_t> = (0..12)
            .map(|i| ikv_t {
                key: (i * 7) % 5,
                val: i,
            })
            .collect();
        ikvsorti(kv.len(), &mut kv);
        assert!(kv.windows(2).all(|w| w[0].key <= w[1].key));
    }
}
