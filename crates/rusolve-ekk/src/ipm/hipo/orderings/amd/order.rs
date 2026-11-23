//! HiGHS source: extern/amd/amd_order.c, extern/amd/amd_defaults.c, extern/amd/amd_valid.c, extern/amd/amd_preprocess.c
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! The user entry point: validate, optionally sort/deduplicate
//! (`amd_preprocess` when `AMD_OK_BUT_JUMBLED`), form the pattern of A+A'
//! (`amd_aat`), size the workspace `slen` and call `amd_1`. `Info` statistics
//! are `double` and are reproduced exactly (the oracle records `status` only).

use super::amd_int;
use super::{
    AMD_AGGRESSIVE, AMD_CONTROL, AMD_DEFAULT_AGGRESSIVE, AMD_DEFAULT_DENSE, AMD_DENSE, AMD_INFO,
    AMD_INVALID, AMD_MEMORY, AMD_N, AMD_NZ, AMD_OK, AMD_OUT_OF_MEMORY, AMD_STATUS, EMPTY,
    amd_int_max,
};

/// HiGHS: amd_order.c::Highs_amd_order (lines 21-195)
pub fn Highs_amd_order(
    n: amd_int,
    Ap: &[amd_int],
    Ai: &[amd_int],
    P: &mut [amd_int],
    Control: Option<&[f64]>,
    Info: Option<&mut [f64]>,
) -> i32 {
    // HiGHS: amd_order.c lines 31-33
    // amd_int *Len, *S, nz, i, *Pinv, info, status, *Rp, *Ri, *Cp, *Ci, ok;
    // size_t nzaat, slen;
    // double mem = 0;

    // HiGHS: amd_order.c line 36: info = Info != (double *) NULL;
    let info = Info.is_some();
    let mut Info = Info;

    // HiGHS: amd_order.c lines 37-45
    if info {
        let info_slice = Info.as_mut().unwrap();
        for i in 0..AMD_INFO {
            info_slice[i] = EMPTY as f64;
        }
        info_slice[AMD_N] = n as f64;
        info_slice[AMD_STATUS] = super::AMD_OK as f64;
    }

    // HiGHS: amd_order.c lines 48-52
    // In Rust, Ap, Ai, P are always valid slices (non-null).
    if n < 0 {
        if info {
            Info.as_mut().unwrap()[AMD_STATUS] = AMD_INVALID as f64;
        }
        return AMD_INVALID;
    }

    // HiGHS: amd_order.c lines 54-57
    if n == 0 {
        return super::AMD_OK;
    }

    // HiGHS: amd_order.c line 59: nz = Ap[n];
    let nz = Ap[n as usize];

    // HiGHS: amd_order.c lines 60-63
    if info {
        Info.as_mut().unwrap()[AMD_NZ] = nz as f64;
    }

    // HiGHS: amd_order.c lines 64-68
    if nz < 0 {
        if info {
            Info.as_mut().unwrap()[AMD_STATUS] = AMD_INVALID as f64;
        }
        return AMD_INVALID;
    }

    // HiGHS: amd_order.c lines 71-76
    // check if n or nz will cause integer overflow
    // ((size_t) n) >= amd_int_max / sizeof(amd_int)
    // ((size_t) nz) >= amd_int_max / sizeof(amd_int)
    // amd_int_max converts to size_t before the division by sizeof (amd_int).
    let limit = (amd_int_max as usize) / std::mem::size_of::<amd_int>();
    if (n as usize) >= limit || (nz as usize) >= limit {
        if info {
            Info.as_mut().unwrap()[AMD_STATUS] = AMD_OUT_OF_MEMORY as f64;
        }
        return AMD_OUT_OF_MEMORY;
    }

    // HiGHS: amd_order.c line 79: status = Highs_amd_valid(n, n, Ap, Ai);
    let status = Highs_amd_valid(n, n, Ap, Ai);

    // HiGHS: amd_order.c lines 81-85
    if status == AMD_INVALID {
        if info {
            Info.as_mut().unwrap()[AMD_STATUS] = AMD_INVALID as f64;
        }
        return AMD_INVALID;
    }

    // HiGHS: amd_order.c lines 88-100
    // allocate two size-n integer workspaces
    let nn = n as usize;
    let mut Len = vec![0 as amd_int; nn];
    let mut Pinv = vec![0 as amd_int; nn];
    // mem += n; mem += n; (mem is double, tracking memory usage)
    let mut mem: f64 = 0.0;
    mem += n as f64;
    mem += n as f64;

    // In Rust, vec! allocation doesn't return null on OOM, it panics.
    // The C code checks for null. We assume allocation succeeds or panics.
    // If we were to emulate the OOM check, we'd need try_reserve, but
    // the rules say malloc becomes vec! and NULL test disappears.

    // HiGHS: amd_order.c lines 102-132
    let mut Rp: Option<Vec<amd_int>>;
    let mut Ri: Option<Vec<amd_int>>;
    let Cp: &[amd_int];
    let Ci: &[amd_int];

    if status == super::AMD_OK_BUT_JUMBLED {
        // sort the input matrix and remove duplicate entries
        // Rp = malloc(nn+1, sizeof(amd_int));
        // Ri = malloc(nz, sizeof(amd_int));
        Rp = Some(vec![0 as amd_int; nn + 1]);
        Ri = Some(vec![0 as amd_int; nz as usize]);
        mem += (n + 1) as f64;
        // mem += MAX(nz, 1);
        // MAX is a C ternary for integers: (a > b) ? a : b
        let max_nz_1 = if nz > 1 { nz } else { 1 };
        mem += max_nz_1 as f64;

        // use Len and Pinv as workspace to create R = A'
        Highs_amd_preprocess(
            n,
            Ap,
            Ai,
            Rp.as_mut().unwrap(),
            Ri.as_mut().unwrap(),
            &mut Len,
            &mut Pinv,
        );
        Cp = Rp.as_ref().unwrap();
        Ci = Ri.as_ref().unwrap();
    } else {
        // order the input matrix as-is. No need to compute R = A' first
        // Rp = NULL; Ri = NULL;
        Cp = Ap;
        Ci = Ai;
    }

    // HiGHS: amd_order.c line 138: nzaat = Highs_amd_aat(n, Cp, Ci, Len, P, Info);
    let nzaat = super::aat::Highs_amd_aat(n, Cp, Ci, &mut Len, P, Info.as_deref_mut());

    // HiGHS: amd_order.c lines 146-160
    // S = NULL;
    // slen = nzaat;
    // ok = ((slen + nzaat/5) >= slen);
    // slen += nzaat/5;
    // for (i = 0; ok && i < 7; i++) { ok = ((slen + nn) > slen); slen += nn; }
    // mem += slen;
    // ok = ok && (slen < SIZE_T_MAX / sizeof(amd_int));
    // if (ok) { S = malloc(slen, sizeof(amd_int)); }

    // size_t arithmetic wraps in C; the overflow checks rely on it.
    let mut slen: usize = nzaat;
    let mut ok: bool = slen.wrapping_add(nzaat / 5) >= slen;
    slen = slen.wrapping_add(nzaat / 5);

    let mut i: amd_int = 0;
    while ok && i < 7 {
        ok = slen.wrapping_add(nn) > slen;
        slen = slen.wrapping_add(nn);
        i += 1;
    }

    mem += slen as f64;

    // SIZE_T_MAX / sizeof(amd_int)
    let size_t_max_limit = usize::MAX / std::mem::size_of::<amd_int>();
    ok = ok && (slen < size_t_max_limit);

    let mut S: Option<Vec<amd_int>> = None;
    if ok {
        S = Some(vec![0 as amd_int; slen]);
    }

    // HiGHS: amd_order.c lines 162-171
    if S.is_none() {
        // :: out of memory :: (or problem too large)
        // Highs_SuiteSparse_free(Rp); ... dropped
        if info {
            Info.as_mut().unwrap()[AMD_STATUS] = AMD_OUT_OF_MEMORY as f64;
        }
        return AMD_OUT_OF_MEMORY;
    }

    // HiGHS: amd_order.c lines 172-176
    if info {
        // memory usage, in bytes.
        Info.as_mut().unwrap()[AMD_MEMORY] = mem * (std::mem::size_of::<amd_int>() as f64);
    }

    // HiGHS: amd_order.c line 182: Highs_amd_1(n, Cp, Ci, P, Pinv, Len, slen, S, Control, Info);
    // Note: slen in C is size_t, but Highs_amd_1 takes amd_int slen.
    // The C code passes slen (size_t) to amd_int slen.
    // We cast slen to amd_int.
    super::aat::Highs_amd_1(
        n,
        Cp,
        Ci,
        P,
        &mut Pinv,
        &mut Len,
        slen as amd_int,
        S.as_mut().unwrap(),
        Control,
        Info.as_deref_mut(),
    );

    // HiGHS: amd_order.c lines 188-194
    // Highs_SuiteSparse_free(Rp); ... dropped
    if info {
        Info.as_mut().unwrap()[AMD_STATUS] = status as f64;
    }
    status
}

/// HiGHS: amd_defaults.c::Highs_amd_defaults (lines 21-37)
pub fn Highs_amd_defaults(Control: &mut [f64]) {
    // HiGHS: amd_defaults.c lines 28-36
    // if (Control != NULL)
    // In Rust, Control is a slice, assumed non-null.
    for i in 0..AMD_CONTROL {
        Control[i] = 0.0;
    }
    Control[AMD_DENSE] = AMD_DEFAULT_DENSE;
    Control[AMD_AGGRESSIVE] = AMD_DEFAULT_AGGRESSIVE;
}

/// HiGHS: amd_valid.c::Highs_amd_valid (lines 38-93)
pub fn Highs_amd_valid(n_row: amd_int, n_col: amd_int, Ap: &[amd_int], Ai: &[amd_int]) -> i32 {
    // HiGHS: amd_valid.c lines 47-48
    // amd_int nz, j, p1, p2, ilast, i, p;
    // int result = AMD_OK;
    let mut result: i32 = super::AMD_OK;

    // HiGHS: amd_valid.c lines 50-53
    if n_row < 0 || n_col < 0 {
        return AMD_INVALID;
    }

    // HiGHS: amd_valid.c line 54: nz = Ap[n_col];
    let nz = Ap[n_col as usize];

    // HiGHS: amd_valid.c lines 55-60
    if Ap[0] != 0 || nz < 0 {
        return AMD_INVALID;
    }

    // HiGHS: amd_valid.c lines 61-91
    for j in 0..n_col {
        let p1 = Ap[j as usize];
        let p2 = Ap[(j + 1) as usize];

        // HiGHS: amd_valid.c lines 66-71
        if p1 > p2 {
            return AMD_INVALID;
        }

        let mut ilast = EMPTY;
        let mut p = p1;
        while p < p2 {
            let i = Ai[p as usize];

            // HiGHS: amd_valid.c lines 77-82
            if i < 0 || i >= n_row {
                return AMD_INVALID;
            }

            // HiGHS: amd_valid.c lines 83-88
            if i <= ilast {
                result = super::AMD_OK_BUT_JUMBLED;
            }
            ilast = i;
            p += 1;
        }
    }

    result
}

/// HiGHS: amd_preprocess.c::Highs_amd_preprocess (lines 25-107)
pub fn Highs_amd_preprocess(
    n: amd_int,
    Ap: &[amd_int],
    Ai: &[amd_int],
    Rp: &mut [amd_int],
    Ri: &mut [amd_int],
    W: &mut [amd_int],
    Flag: &mut [amd_int],
) {
    // HiGHS: amd_preprocess.c lines 52-56
    for i in 0..n {
        W[i as usize] = 0;
        Flag[i as usize] = EMPTY;
    }

    // HiGHS: amd_preprocess.c lines 57-70
    for j in 0..n {
        let p2 = Ap[(j + 1) as usize];
        let mut p = Ap[j as usize];
        while p < p2 {
            let i = Ai[p as usize];
            if Flag[i as usize] != j {
                W[i as usize] += 1;
                Flag[i as usize] = j;
            }
            p += 1;
        }
    }

    // HiGHS: amd_preprocess.c lines 76-80
    Rp[0] = 0;
    for i in 0..n {
        Rp[(i + 1) as usize] = Rp[i as usize] + W[i as usize];
    }

    // HiGHS: amd_preprocess.c lines 81-85
    for i in 0..n {
        W[i as usize] = Rp[i as usize];
        Flag[i as usize] = EMPTY;
    }

    // HiGHS: amd_preprocess.c lines 92-105
    for j in 0..n {
        let p2 = Ap[(j + 1) as usize];
        let mut p = Ap[j as usize];
        while p < p2 {
            let i = Ai[p as usize];
            if Flag[i as usize] != j {
                Ri[W[i as usize] as usize] = j;
                W[i as usize] += 1;
                Flag[i as usize] = j;
            }
            p += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn defaults_values() {
        let mut c = [0.0; crate::ipm::hipo::orderings::amd::AMD_CONTROL];
        Highs_amd_defaults(&mut c);
        assert_eq!(
            c[crate::ipm::hipo::orderings::amd::AMD_DENSE],
            crate::ipm::hipo::orderings::amd::AMD_DEFAULT_DENSE
        );
        assert_eq!(
            c[crate::ipm::hipo::orderings::amd::AMD_AGGRESSIVE],
            crate::ipm::hipo::orderings::amd::AMD_DEFAULT_AGGRESSIVE
        );
    }

    #[test]
    fn order_arrow_matrix() {
        // arrow matrix: node 0 adjacent to all others; AMD orders it last
        let n = 5;
        let mut ap = vec![0];
        let mut ai = Vec::new();
        for j in 0..n {
            if j == 0 {
                ai.extend(1..n);
            } else {
                ai.push(0);
            }
            ap.push(ai.len() as amd_int);
        }
        let mut p = vec![0; n as usize];
        let st = Highs_amd_order(n, &ap, &ai, &mut p, None, None);
        assert_eq!(st, crate::ipm::hipo::orderings::amd::AMD_OK);
        assert_eq!(p[n as usize - 1], 0);
    }

    #[test]
    fn jumbled_matches_sorted() {
        // amd_order.c lines 102-124: an unsorted / duplicated symmetric pattern is
        // preprocessed into R = A' (sorted, deduplicated), which for a symmetric
        // pattern equals the sorted input, so the ordering is identical and the
        // status is AMD_OK_BUT_JUMBLED.
        let ap = vec![0, 3, 5, 7, 9];
        let ai = vec![1, 2, 3, 0, 2, 0, 1, 0, 3];
        let ap_j = vec![0, 4, 7, 9, 11];
        let ai_j = vec![3, 1, 2, 1, 2, 0, 0, 0, 1, 0, 3];
        let mut info = [0.0; crate::ipm::hipo::orderings::amd::AMD_INFO];
        let mut p = vec![0; 4];
        let mut p_j = vec![0; 4];
        assert_eq!(
            Highs_amd_valid(4, 4, &ap, &ai),
            crate::ipm::hipo::orderings::amd::AMD_OK
        );
        assert_eq!(
            Highs_amd_order(4, &ap, &ai, &mut p, None, None),
            crate::ipm::hipo::orderings::amd::AMD_OK
        );
        let st = Highs_amd_order(4, &ap_j, &ai_j, &mut p_j, None, Some(&mut info));
        assert_eq!(st, crate::ipm::hipo::orderings::amd::AMD_OK_BUT_JUMBLED);
        assert_eq!(
            info[crate::ipm::hipo::orderings::amd::AMD_STATUS],
            st as f64
        );
        // amd_order.c line 62 stores Ap[n] = 11, then amd_aat.c line 172 overwrites
        // it with the nz of the deduplicated R.
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_NZ], 9.0);
        assert_eq!(p, p_j);
    }

    #[test]
    fn invalid_inputs() {
        use crate::ipm::hipo::orderings::amd::{AMD_INVALID, AMD_N, AMD_OK, AMD_STATUS};
        let mut info = [0.0; crate::ipm::hipo::orderings::amd::AMD_INFO];
        let mut p = vec![0; 2];
        // n < 0 (lines 48-52)
        assert_eq!(
            Highs_amd_order(-1, &[0], &[], &mut p, None, Some(&mut info)),
            AMD_INVALID
        );
        assert_eq!(info[AMD_STATUS], AMD_INVALID as f64);
        assert_eq!(info[AMD_N], -1.0);
        // n == 0 (lines 54-57): Info stays cleared to EMPTY except N / STATUS
        assert_eq!(
            Highs_amd_order(0, &[0], &[], &mut p, None, Some(&mut info)),
            AMD_OK
        );
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_NZ], -1.0);
        // row index out of range (amd_valid.c lines 77-82)
        assert_eq!(
            Highs_amd_order(2, &[0, 1, 2], &[0, 2], &mut p, None, Some(&mut info)),
            AMD_INVALID
        );
        assert_eq!(info[AMD_STATUS], AMD_INVALID as f64);
        // Ap[0] != 0 (amd_valid.c lines 55-60)
        assert_eq!(Highs_amd_valid(2, 2, &[1, 1, 2], &[0, 1]), AMD_INVALID);
        // descending column pointers (amd_valid.c lines 66-71)
        assert_eq!(Highs_amd_valid(2, 2, &[0, 2, 1], &[0, 1]), AMD_INVALID);
    }
}
