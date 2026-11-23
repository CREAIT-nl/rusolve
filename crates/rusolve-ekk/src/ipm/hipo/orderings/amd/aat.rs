//! HiGHS source: extern/amd/amd_aat.c, extern/amd/amd_1.c
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! Pattern of A+A' (`amd_aat`: column counts and symmetry statistics) and
//! `amd_1`, which builds that pattern in the workspace `S` and hands it to
//! `amd_2`. `S` is carved into `Pe, Nv, Head, Elen, Degree, W, Iw` by
//! pointer offsets in C; the Rust splits it with successive `split_at_mut`
//! in the same order.

use super::EMPTY;
use super::amd_int;

/// HiGHS: amd_aat.c::Highs_amd_aat (lines 20-179)
pub fn Highs_amd_aat(
    n: amd_int,
    Ap: &[amd_int],
    Ai: &[amd_int],
    Len: &mut [amd_int],
    Tp: &mut [amd_int],
    Info: Option<&mut [f64]>,
) -> usize {
    // Rebind so the Option can be reborrowed mutably at both Info use sites.
    let mut Info = Info;
    // HiGHS: amd_aat.c lines 30-32
    let mut p1: amd_int;
    let mut p2: amd_int;
    let mut p: amd_int;
    let mut i: amd_int;
    let mut j: amd_int;
    let mut pj: amd_int;
    let mut pj2: amd_int;
    let mut k: amd_int;
    let mut nzdiag: amd_int;
    let mut nzboth: amd_int;
    let mut sym: f64;
    let mut nzaat: usize;

    // HiGHS: amd_aat.c lines 35-43
    if let Some(info) = Info.as_deref_mut() {
        // clear the Info array, if it exists
        for i in 0..super::AMD_INFO as amd_int {
            info[i as usize] = EMPTY as f64;
        }
        info[super::AMD_STATUS] = super::AMD_OK as f64;
    }

    // HiGHS: amd_aat.c lines 45-48
    for k in 0..n {
        Len[k as usize] = 0;
    }

    // HiGHS: amd_aat.c lines 50-52
    nzdiag = 0;
    nzboth = 0;
    let nz: amd_int = Ap[n as usize];

    // HiGHS: amd_aat.c lines 54-121
    for k in 0..n {
        p1 = Ap[k as usize];
        p2 = Ap[(k + 1) as usize];

        // construct A+A'
        p = p1;
        while p < p2 {
            // scan the upper triangular part of A
            j = Ai[p as usize];
            if j < k {
                // entry A (j,k) is in the strictly upper triangular part,
                // add both A (j,k) and A (k,j) to the matrix A+A'
                Len[j as usize] += 1;
                Len[k as usize] += 1;

                p += 1;
            } else if j == k {
                // skip the diagonal
                p += 1;
                nzdiag += 1;
                break;
            } else {
                // j > k: first entry below the diagonal
                break;
            }
            // scan lower triangular part of A, in column j until reaching
            // row k.  Start where last scan left off.

            pj2 = Ap[(j + 1) as usize];
            pj = Tp[j as usize];
            while pj < pj2 {
                i = Ai[pj as usize];
                if i < k {
                    // A (i,j) is only in the lower part, not in upper.
                    // add both A (i,j) and A (j,i) to the matrix A+A'
                    Len[i as usize] += 1;
                    Len[j as usize] += 1;

                    pj += 1;
                } else if i == k {
                    // entry A (k,j) in lower part and A (j,k) in upper
                    pj += 1;
                    nzboth += 1;
                    break;
                } else {
                    // i > k: consider this entry later, when k advances to i
                    break;
                }
            }
            Tp[j as usize] = pj;
        }
        // Tp [k] points to the entry just below the diagonal in column k
        Tp[k as usize] = p;
    }

    // HiGHS: amd_aat.c lines 123-136
    // clean up, for remaining mismatched entries
    for j in 0..n {
        pj = Tp[j as usize];
        while pj < Ap[(j + 1) as usize] {
            i = Ai[pj as usize];
            // A (i,j) is only in the lower part, not in upper.
            // add both A (i,j) and A (j,i) to the matrix A+A'
            Len[i as usize] += 1;
            Len[j as usize] += 1;

            pj += 1;
        }
    }

    // HiGHS: amd_aat.c lines 138-155
    // compute the symmetry of the nonzero pattern of A
    if nz == nzdiag {
        sym = 1.0;
    } else {
        sym = (2.0 * (nzboth as f64)) / ((nz - nzdiag) as f64);
    }

    // HiGHS: amd_aat.c lines 157-161
    nzaat = 0;
    for k in 0..n {
        nzaat += Len[k as usize] as usize;
    }

    // HiGHS: amd_aat.c lines 168-176
    if let Some(info) = Info {
        info[super::AMD_STATUS] = super::AMD_OK as f64;
        info[super::AMD_N] = n as f64;
        info[super::AMD_NZ] = nz as f64;
        info[super::AMD_SYMMETRY] = sym;
        info[super::AMD_NZDIAG] = nzdiag as f64;
        info[super::AMD_NZ_A_PLUS_AT] = nzaat as f64;
    }

    // HiGHS: amd_aat.c line 178
    nzaat
}

/// HiGHS: amd_1.c::Highs_amd_1 (lines 29-172)
pub fn Highs_amd_1(
    n: amd_int,
    Ap: &[amd_int],
    Ai: &[amd_int],
    P: &mut [amd_int],
    Pinv: &mut [amd_int],
    Len: &mut [amd_int],
    slen: amd_int,
    S: &mut [amd_int],
    Control: Option<&[f64]>,
    Info: Option<&mut [f64]>,
) {
    // HiGHS: amd_1.c lines 44-45
    let mut i: amd_int;
    let mut j: amd_int;
    let mut k: amd_int;
    let mut p: amd_int;
    let mut pfree: amd_int;
    let mut pj: amd_int;
    let mut p1: amd_int;
    let mut p2: amd_int;
    let mut pj2: amd_int;

    // HiGHS: amd_1.c lines 53-61
    // iwlen = slen - 6*n
    // s = S; Pe = s; s += n; Nv = s; s += n; Head = s; s += n;
    // Elen = s; s += n; Degree = s; s += n; W = s; s += n; Iw = s; s += iwlen;
    let iwlen: amd_int = slen - 6 * n;

    // Carve S into Pe, Nv, Head, Elen, Degree, W, Iw in C order.
    // Sp aliases Nv, Tp aliases W (C lines 66-67).
    let n_us = n as usize;
    let iwlen_us = iwlen as usize;

    let (pe_slice, rest) = S.split_at_mut(n_us);
    let (nv_slice, rest) = rest.split_at_mut(n_us);
    let (head_slice, rest) = rest.split_at_mut(n_us);
    let (elen_slice, rest) = rest.split_at_mut(n_us);
    let (degree_slice, rest) = rest.split_at_mut(n_us);
    let (w_slice, rest) = rest.split_at_mut(n_us);
    let (iw_slice, _) = rest.split_at_mut(iwlen_us);

    // Sp = Nv (C line 66), Tp = W (C line 67)
    // We use nv_slice as Sp and w_slice as Tp directly.

    // HiGHS: amd_1.c lines 68-74
    pfree = 0;
    for j in 0..n {
        pe_slice[j as usize] = pfree;
        nv_slice[j as usize] = pfree; // Sp[j] = pfree
        pfree += Len[j as usize];
    }

    // HiGHS: amd_1.c lines 82-146
    for k in 0..n {
        p1 = Ap[k as usize];
        p2 = Ap[(k + 1) as usize];

        // construct A+A'
        p = p1;
        while p < p2 {
            // scan the upper triangular part of A
            j = Ai[p as usize];

            if j < k {
                // entry A (j,k) in the strictly upper triangular part

                iw_slice[nv_slice[j as usize] as usize] = k;
                nv_slice[j as usize] += 1;
                iw_slice[nv_slice[k as usize] as usize] = j;
                nv_slice[k as usize] += 1;
                p += 1;
            } else if j == k {
                // skip the diagonal
                p += 1;
                break;
            } else {
                // j > k: first entry below the diagonal
                break;
            }
            // scan lower triangular part of A, in column j until reaching
            // row k.  Start where last scan left off.

            pj2 = Ap[(j + 1) as usize];
            pj = w_slice[j as usize]; // Tp[j]
            while pj < pj2 {
                i = Ai[pj as usize];

                if i < k {
                    // A (i,j) is only in the lower part, not in upper

                    iw_slice[nv_slice[i as usize] as usize] = j;
                    nv_slice[i as usize] += 1;
                    iw_slice[nv_slice[j as usize] as usize] = i;
                    nv_slice[j as usize] += 1;
                    pj += 1;
                } else if i == k {
                    // entry A (k,j) in lower part and A (j,k) in upper
                    pj += 1;
                    break;
                } else {
                    // i > k: consider this entry later, when k advances to i
                    break;
                }
            }
            w_slice[j as usize] = pj; // Tp[j] = pj
        }
        w_slice[k as usize] = p; // Tp[k] = p
    }

    // HiGHS: amd_1.c lines 148-161
    // clean up, for remaining mismatched entries
    for j in 0..n {
        pj = w_slice[j as usize]; // Tp[j]
        while pj < Ap[(j + 1) as usize] {
            i = Ai[pj as usize];

            // A (i,j) is only in the lower part, not in upper

            iw_slice[nv_slice[i as usize] as usize] = j;
            nv_slice[i as usize] += 1;
            iw_slice[nv_slice[j as usize] as usize] = i;
            nv_slice[j as usize] += 1;

            pj += 1;
        }
    }

    // Tp and Sp no longer needed

    // HiGHS: amd_1.c lines 170-171
    // Highs_amd_2(n, Pe, Iw, Len, iwlen, pfree, Nv, Pinv, P, Head, Elen, Degree, W, Control, Info)
    // C passes: Nv, Pinv, P as Nv, Next, Last respectively.
    super::amd_2::Highs_amd_2(
        n,
        pe_slice,
        iw_slice,
        Len,
        iwlen,
        pfree,
        nv_slice,
        Pinv,
        P,
        head_slice,
        elen_slice,
        degree_slice,
        w_slice,
        Control,
        Info,
    );
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn aat_symmetric_counts() {
        // symmetric 3x3 pattern with off-diagonals (0,1), (1,2)
        let ap = [0, 2, 5, 7];
        let ai = [0, 1, 0, 1, 2, 1, 2];
        let mut len = [0; 3];
        let mut tp = [0; 3];
        let nz = Highs_amd_aat(3, &ap, &ai, &mut len, &mut tp, None);
        assert_eq!(nz, 4);
        assert_eq!(len, [1, 2, 1]);
    }

    #[test]
    fn aat_unsymmetric_info() {
        // A(1,0) lower only (inner-loop i < k branch), A(0,2) upper only
        let ap = [0, 2, 3, 5];
        let ai = [0, 1, 1, 0, 2];
        let mut len = [0; 3];
        let mut tp = [0; 3];
        let mut info = [0.0; crate::ipm::hipo::orderings::amd::AMD_INFO];
        let nz = Highs_amd_aat(3, &ap, &ai, &mut len, &mut tp, Some(&mut info));
        assert_eq!(nz, 4);
        assert_eq!(len, [2, 1, 1]);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_STATUS], 0.0);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_N], 3.0);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_NZ], 5.0);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_SYMMETRY], 0.0);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_NZDIAG], 3.0);
        assert_eq!(
            info[crate::ipm::hipo::orderings::amd::AMD_NZ_A_PLUS_AT],
            4.0
        );
        // entries not written after the clear stay EMPTY
        for v in &info[6..] {
            assert_eq!(*v, -1.0);
        }
    }

    #[test]
    fn aat_nzboth_symmetry() {
        // A(2,0) in lower part of column 0 and A(0,2) in upper part of column 2
        let ap = [0, 2, 2, 3];
        let ai = [0, 2, 0];
        let mut len = [0; 3];
        let mut tp = [0; 3];
        let mut info = [0.0; crate::ipm::hipo::orderings::amd::AMD_INFO];
        let nz = Highs_amd_aat(3, &ap, &ai, &mut len, &mut tp, Some(&mut info));
        assert_eq!(nz, 2);
        assert_eq!(len, [1, 0, 1]);
        // sym = (2 * 1) / (3 - 1)
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_SYMMETRY], 1.0);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_NZDIAG], 1.0);
    }

    #[test]
    fn aat_cleanup_loop_and_diag_only() {
        // A(1,0) only, no diagonal: counted by the clean-up loop
        let ap = [0, 1, 1];
        let ai = [1];
        let mut len = [0; 2];
        let mut tp = [0; 2];
        let mut info = [0.0; crate::ipm::hipo::orderings::amd::AMD_INFO];
        let nz = Highs_amd_aat(2, &ap, &ai, &mut len, &mut tp, Some(&mut info));
        assert_eq!(nz, 2);
        assert_eq!(len, [1, 1]);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_SYMMETRY], 0.0);
        // diagonal only: nz == nzdiag gives sym = 1
        let ap = [0, 1, 2];
        let ai = [0, 1];
        let nz = Highs_amd_aat(2, &ap, &ai, &mut len, &mut tp, Some(&mut info));
        assert_eq!(nz, 0);
        assert_eq!(len, [0, 0]);
        assert_eq!(info[crate::ipm::hipo::orderings::amd::AMD_SYMMETRY], 1.0);
    }

    #[test]
    fn amd1_tridiagonal_permutation() {
        // tridiagonal 4x4, full symmetric pattern
        let n = 4;
        let ap = [0, 2, 5, 8, 10];
        let ai = [0, 1, 0, 1, 2, 1, 2, 3, 2, 3];
        let mut len = [0; 4];
        let mut tp = [0; 4];
        let nzaat = Highs_amd_aat(n, &ap, &ai, &mut len, &mut tp, None) as amd_int;
        assert_eq!(nzaat, 6);
        let slen = nzaat + nzaat / 5 + 8 * n;
        let mut s = vec![0; slen as usize];
        let mut p = [0; 4];
        let mut pinv = [0; 4];
        let mut info = [0.0; crate::ipm::hipo::orderings::amd::AMD_INFO];
        Highs_amd_1(
            n,
            &ap,
            &ai,
            &mut p,
            &mut pinv,
            &mut len,
            slen,
            &mut s,
            None,
            Some(&mut info),
        );
        let mut seen = [false; 4];
        for k in 0..4 {
            let pk = p[k] as usize;
            assert!(!seen[pk]);
            seen[pk] = true;
            assert_eq!(pinv[pk] as usize, k);
        }
    }
}
