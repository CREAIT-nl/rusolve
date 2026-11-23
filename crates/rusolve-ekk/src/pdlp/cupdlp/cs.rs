//! HiGHS source: pdlp/cupdlp/cupdlp_cs.h, pdlp/cupdlp/cupdlp_cs.c
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The CSparse subset cuPDLP uses for `csr2csc` / `csc2csr`
//! (`cupdlp_utils.c` 1199-1254). The raw allocator wrappers
//! (`cupdlp_dcs_malloc` / `calloc` / `free` / `realloc`, lines 5-29) are
//! `Vec` allocation; `cupdlp_dcs_spfree` (61-68) is `Drop`;
//! `cupdlp_dcs_print` (145-187) is a printer and is not ported.
//! A `NULL` result (out of memory or rejected input) is `None`.

use super::defs::*;

/// HiGHS: cupdlp_cs.h `struct cupdlp_cs_sparse` (`cupdlp_dcs`).
///
/// `x` is `None` when the matrix was allocated with `values == 0`
/// (`A->x = NULL`); every `if (A->x)` / `if (Cx)` test reads it.
#[derive(Clone, Debug, Default)]
pub struct cupdlp_dcs {
    /// maximum number of entries
    pub nzmax: i32,
    /// number of rows
    pub m: i32,
    /// number of columns
    pub n: i32,
    /// column pointers (size n+1) or col indices (size nzmax)
    pub p: Vec<i32>,
    /// row indices, size nzmax
    pub i: Vec<i32>,
    /// numerical values, size nzmax
    pub x: Option<Vec<f64>>,
    /// # of entries in triplet matrix, -1 for compressed-col
    pub nz: i32,
}

/// HiGHS: cupdlp_cs.h `IS_CSC(A)` (`A && (A->nz == -1)`).
pub fn IS_CSC(A: Option<&cupdlp_dcs>) -> bool {
    match A {
        Some(a) => a.nz == -1,
        None => false,
    }
}

/// HiGHS: cupdlp_cs.h `IS_TRIPLET(A)` (`A && (A->nz >= 0)`).
pub fn IS_TRIPLET(A: Option<&cupdlp_dcs>) -> bool {
    match A {
        Some(a) => a.nz >= 0,
        None => false,
    }
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_spalloc (lines 31-44)
pub fn cupdlp_dcs_spalloc(
    m: i32,
    n: i32,
    nzmax: i32,
    values: i32,
    triplet: i32,
) -> Option<cupdlp_dcs> {
    // cupdlp_dcs_calloc(1, sizeof(cupdlp_dcs)) is a zeroed struct; allocation
    // cannot fail in Rust, so the NULL test disappears.
    let mut A = cupdlp_dcs::default();
    // define dimensions and nzmax
    A.m = m;
    A.n = n;
    let nzmax = MAX(nzmax, 1);
    A.nzmax = nzmax;
    // allocate triplet or comp.col
    A.nz = if triplet != 0 { 0 } else { -1 };
    // cupdlp_dcs_malloc(triplet ? nzmax : n + 1, sizeof(int))
    let p_size = if triplet != 0 { nzmax } else { n + 1 };
    A.p = vec![0i32; p_size as usize];
    // cupdlp_dcs_malloc(nzmax, sizeof(int))
    A.i = vec![0i32; nzmax as usize];
    // cupdlp_dcs_malloc(nzmax, sizeof(double)) or NULL
    A.x = if values != 0 {
        Some(vec![0.0f64; nzmax as usize])
    } else {
        None
    };
    // The NULL tests on A->p, A->i, A->x cannot fire in Rust; allocation
    // succeeded, so return Some(A).
    Some(A)
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_sprealloc (lines 47-58)
pub fn cupdlp_dcs_sprealloc(A: Option<&mut cupdlp_dcs>, nzmax: i32) -> i32 {
    // if (!A) return (0);
    let A = match A {
        Some(a) => a,
        None => return 0,
    };
    // if (nzmax <= 0) nzmax = IS_CSC(A) ? (A->p[A->n]) : A->nz;
    let mut nzmax = nzmax;
    if nzmax <= 0 {
        nzmax = if IS_CSC(Some(A)) {
            A.p[A.n as usize]
        } else {
            A.nz
        };
    }
    // nzmax = MAX(nzmax, 1);
    nzmax = MAX(nzmax, 1);
    // A->i = cupdlp_dcs_realloc(A->i, nzmax, sizeof(int), &oki);
    // In Rust, resize never fails; oki = 1.
    A.i.resize(nzmax as usize, 0);
    let oki = 1;
    // if (IS_TRIPLET(A)) A->p = cupdlp_dcs_realloc(A->p, nzmax, sizeof(int), &okj);
    let mut okj = 1;
    if IS_TRIPLET(Some(A)) {
        A.p.resize(nzmax as usize, 0);
    }
    // if (A->x) A->x = cupdlp_dcs_realloc(A->x, nzmax, sizeof(double), &okx);
    let mut okx = 1;
    if let Some(x) = A.x.as_mut() {
        x.resize(nzmax as usize, 0.0);
    }
    // ok = (oki && okj && okx);
    let ok = oki != 0 && okj != 0 && okx != 0;
    // if (ok) A->nzmax = nzmax;
    if ok {
        A.nzmax = nzmax;
    }
    // return (ok);
    if ok { 1 } else { 0 }
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_done (lines 71-76)
///
/// The workspaces `w` and `x` are freed by going out of scope at the call
/// site; only `C` and `ok` remain.
pub fn cupdlp_dcs_done(C: Option<cupdlp_dcs>, ok: i32) -> Option<cupdlp_dcs> {
    // cupdlp_dcs_free(w); /* free workspace */
    // cupdlp_dcs_free(x);
    // return (ok ? C : cupdlp_dcs_spfree(C));
    // cupdlp_dcs_spfree frees C and returns NULL; in Rust, dropping C is the
    // same as returning None.
    if ok != 0 { C } else { None }
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_cumsum (lines 78-90)
pub fn cupdlp_dcs_cumsum(p: &mut [i32], c: &mut [i32], n: i32) -> f64 {
    // if (!p || !c) return (-1);
    // In Rust, slices are always non-null; the check disappears.
    let mut nz: i32 = 0;
    let mut nz2: f64 = 0.0;
    let n_usize = n as usize;
    for i in 0..n_usize {
        p[i] = nz;
        nz += c[i];
        // also in double to avoid int overflow
        nz2 += c[i] as f64;
        // also copy p[0..n-1] back into c[0..n-1]
        c[i] = p[i];
    }
    p[n_usize] = nz;
    // return sum (c [0..n-1])
    nz2
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_entry (lines 93-102)
pub fn cupdlp_dcs_entry(T: &mut cupdlp_dcs, i: i32, j: i32, x: f64) -> i32 {
    // if (!IS_TRIPLET(T) || i < 0 || j < 0) return (0);
    if !IS_TRIPLET(Some(T)) || i < 0 || j < 0 {
        return 0;
    }
    // if (T->nz >= T->nzmax && !cupdlp_dcs_sprealloc(T, 2 * (T->nzmax))) return (0);
    if T.nz >= T.nzmax {
        let new_nzmax = 2 * T.nzmax;
        let ret = cupdlp_dcs_sprealloc(Some(T), new_nzmax);
        if ret == 0 {
            return 0;
        }
    }
    // if (T->x) T->x[T->nz] = x;
    if let Some(ref mut xv) = T.x {
        xv[T.nz as usize] = x;
    }
    // T->i[T->nz] = i;
    T.i[T.nz as usize] = i;
    // T->p[T->nz++] = j;
    T.p[T.nz as usize] = j;
    T.nz += 1;
    // T->m = MAX(T->m, i + 1);
    T.m = MAX(T.m, i + 1);
    // T->n = MAX(T->n, j + 1);
    T.n = MAX(T.n, j + 1);
    // return (1);
    1
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_compress (lines 104-128)
pub fn cupdlp_dcs_compress(T: &cupdlp_dcs) -> Option<cupdlp_dcs> {
    // if (!IS_TRIPLET(T)) return (NULL);
    if !IS_TRIPLET(Some(T)) {
        return None;
    }
    let m = T.m;
    let n = T.n;
    // Ti = T->i; Tj = T->p; Tx = T->x;
    let Ti = &T.i;
    let Tj = &T.p;
    let Tx = T.x.as_ref();
    let nz = T.nz;
    // C = cupdlp_dcs_spalloc(m, n, nz, Tx != NULL, 0);
    let mut C = cupdlp_dcs_spalloc(m, n, nz, if Tx.is_some() { 1 } else { 0 }, 0);
    // w = cupdlp_dcs_calloc(n, sizeof(int));
    let mut w = vec![0i32; n as usize];
    // if (!C || !w) return (cupdlp_dcs_done(C, w, NULL, 0));
    // Allocation cannot fail in Rust; the check disappears.
    let C = match C {
        Some(c) => c,
        None => return cupdlp_dcs_done(None, 0),
    };
    let mut C = C;
    // Cp = C->p; Ci = C->i; Cx = C->x;
    // We need mutable access to C's fields while reading T's fields.
    // Since T is borrowed immutably and C is a separate owned value, this is fine.
    // However, we need to access C.p, C.i, C.x mutably.
    // We'll use indices directly on C's vectors.

    // for (k = 0; k < nz; k++) w[Tj[k]]++; /* column counts */
    for k in 0..nz as usize {
        let col = Tj[k] as usize;
        w[col] += 1;
    }
    // cupdlp_dcs_cumsum(Cp, w, n); /* column pointers */
    cupdlp_dcs_cumsum(&mut C.p, &mut w, n);
    // for (k = 0; k < nz; k++) {
    //   Ci[p = w[Tj[k]]++] = Ti[k]; /* A(i,j) is the pth entry in C */
    //   if (Cx) Cx[p] = Tx[k];
    // }
    for k in 0..nz as usize {
        let col = Tj[k] as usize;
        let p = w[col] as usize;
        w[col] += 1;
        C.i[p] = Ti[k];
        if let Some(Cx) = C.x.as_mut() {
            if let Some(Txv) = Tx {
                Cx[p] = Txv[k];
            }
        }
    }
    // return (cupdlp_dcs_done(C, w, NULL, 1));
    cupdlp_dcs_done(Some(C), 1)
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_norm (lines 130-143)
pub fn cupdlp_dcs_norm(A: &cupdlp_dcs) -> f64 {
    // if (!IS_CSC(A) || !A->x) return (-1);
    if !IS_CSC(Some(A)) || A.x.is_none() {
        return -1.0;
    }
    let n = A.n;
    let Ap = &A.p;
    let Ax = A.x.as_ref().unwrap();
    let mut nrm: f64 = 0.0;
    let n_usize = n as usize;
    for j in 0..n_usize {
        let mut s: f64 = 0.0;
        let mut p = Ap[j] as usize;
        let p_end = Ap[j + 1] as usize;
        while p < p_end {
            s += Ax[p].abs();
            p += 1;
        }
        nrm = MAX(nrm, s);
    }
    nrm
}

/// HiGHS: cupdlp_cs.c::cupdlp_dcs_transpose (lines 189-214)
pub fn cupdlp_dcs_transpose(A: &cupdlp_dcs, values: i32) -> Option<cupdlp_dcs> {
    // if (!IS_CSC(A)) return (NULL);
    if !IS_CSC(Some(A)) {
        return None;
    }
    let m = A.m;
    let n = A.n;
    let Ap = &A.p;
    let Ai = &A.i;
    let Ax = A.x.as_ref();
    // C = cupdlp_dcs_spalloc(n, m, Ap[n], values && Ax, 0);
    let has_values = values != 0 && Ax.is_some();
    let mut C = cupdlp_dcs_spalloc(n, m, Ap[n as usize], if has_values { 1 } else { 0 }, 0);
    // w = cupdlp_calloc(m, sizeof(cupdlp_int));
    let mut w = vec![0i32; m as usize];
    // if (!C || !w) return (cupdlp_dcs_done(C, w, NULL, 0));
    // Allocation cannot fail in Rust; the check disappears.
    let C = match C {
        Some(c) => c,
        None => return cupdlp_dcs_done(None, 0),
    };
    let mut C = C;
    // Cp = C->p; Ci = C->i; Cx = C->x;
    // for (p = 0; p < Ap[n]; p++) w[Ai[p]]++; /* row counts */
    let ap_n = Ap[n as usize] as usize;
    for p in 0..ap_n {
        let row = Ai[p] as usize;
        w[row] += 1;
    }
    // cupdlp_dcs_cumsum(Cp, w, m); /* row pointers */
    cupdlp_dcs_cumsum(&mut C.p, &mut w, m);
    // for (j = 0; j < n; j++) {
    //   for (p = Ap[j]; p < Ap[j + 1]; p++) {
    //     Ci[q = w[Ai[p]]++] = j; /* place A(i,j) as entry C(j,i) */
    //     if (Cx) Cx[q] = Ax[p];
    //   }
    // }
    let n_usize = n as usize;
    for j in 0..n_usize {
        let mut p = Ap[j] as usize;
        let p_end = Ap[j + 1] as usize;
        while p < p_end {
            let row = Ai[p] as usize;
            let q = w[row] as usize;
            w[row] += 1;
            C.i[q] = j as i32;
            if let Some(Cx) = C.x.as_mut() {
                if let Some(Axv) = Ax {
                    Cx[q] = Axv[p];
                }
            }
            p += 1;
        }
    }
    // return (cupdlp_dcs_done(C, w, NULL, 1));
    cupdlp_dcs_done(Some(C), 1)
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// 2x3 matrix [[1 0 2] [0 3 0]] as triplets, compressed, transposed.
    #[test]
    fn triplet_compress_transpose_round_trip() {
        let mut t = cupdlp_dcs_spalloc(0, 0, 1, 1, 1).unwrap();
        assert_eq!(t.nz, 0);
        assert_eq!(cupdlp_dcs_entry(&mut t, 0, 0, 1.0), 1);
        assert_eq!(cupdlp_dcs_entry(&mut t, 1, 1, 3.0), 1);
        assert_eq!(cupdlp_dcs_entry(&mut t, 0, 2, 2.0), 1);
        assert_eq!((t.m, t.n, t.nz), (2, 3, 3));
        let c = cupdlp_dcs_compress(&t).unwrap();
        assert!(IS_CSC(Some(&c)));
        assert_eq!(c.p, vec![0, 1, 2, 3]);
        assert_eq!(c.i, vec![0, 1, 0]);
        assert_eq!(cupdlp_dcs_norm(&c), 3.0);
        let ct = cupdlp_dcs_transpose(&c, 1).unwrap();
        assert_eq!((ct.m, ct.n), (3, 2));
        assert_eq!(ct.p, vec![0, 2, 3]);
        assert_eq!(ct.i, vec![0, 2, 1]);
        assert_eq!(ct.x.as_deref(), Some(&[1.0, 2.0, 3.0][..]));
    }

    #[test]
    fn cumsum_returns_total_and_copies_back() {
        let mut p = vec![0; 4];
        let mut c = vec![2, 0, 3];
        assert_eq!(cupdlp_dcs_cumsum(&mut p, &mut c, 3), 5.0);
        assert_eq!(p, vec![0, 2, 2, 5]);
        assert_eq!(c, vec![0, 2, 2]);
    }

    #[test]
    fn spalloc_no_values() {
        let a = cupdlp_dcs_spalloc(2, 3, 5, 0, 0).unwrap();
        assert_eq!(a.m, 2);
        assert_eq!(a.n, 3);
        assert_eq!(a.nzmax, 5);
        assert_eq!(a.nz, -1);
        assert!(a.x.is_none());
        assert_eq!(a.p.len(), 4); // n + 1
        assert_eq!(a.i.len(), 5);
    }

    #[test]
    fn spalloc_triplet() {
        let a = cupdlp_dcs_spalloc(0, 0, 1, 1, 1).unwrap();
        assert_eq!(a.nz, 0);
        assert!(a.x.is_some());
        assert_eq!(a.p.len(), 1); // nzmax
        assert_eq!(a.i.len(), 1);
    }

    #[test]
    fn spalloc_nzmax_clamped() {
        let a = cupdlp_dcs_spalloc(1, 1, 0, 0, 0).unwrap();
        assert_eq!(a.nzmax, 1);
    }

    #[test]
    fn sprealloc_null() {
        assert_eq!(cupdlp_dcs_sprealloc(None, 10), 0);
    }

    #[test]
    fn sprealloc_csc() {
        let mut a = cupdlp_dcs_spalloc(2, 2, 4, 1, 0).unwrap();
        // Set up as CSC with nz = -1
        a.nz = -1;
        a.p = vec![0, 2, 4];
        a.i = vec![0, 1, 0, 1];
        a.x = Some(vec![1.0, 2.0, 3.0, 4.0]);
        // nzmax <= 0, so nzmax = A->p[A->n] = A->p[2] = 4
        let ret = cupdlp_dcs_sprealloc(Some(&mut a), 0);
        assert_eq!(ret, 1);
        assert_eq!(a.nzmax, 4);
    }

    #[test]
    fn sprealloc_triplet() {
        let mut a = cupdlp_dcs_spalloc(0, 0, 2, 1, 1).unwrap();
        assert_eq!(a.nz, 0);
        // nzmax <= 0, so nzmax = A->nz = 0, then MAX(0, 1) = 1
        let ret = cupdlp_dcs_sprealloc(Some(&mut a), 0);
        assert_eq!(ret, 1);
        assert_eq!(a.nzmax, 1);
    }

    #[test]
    fn done_ok() {
        let c = cupdlp_dcs_spalloc(1, 1, 1, 0, 0).unwrap();
        let r = cupdlp_dcs_done(Some(c), 1);
        assert!(r.is_some());
    }

    #[test]
    fn done_not_ok() {
        let c = cupdlp_dcs_spalloc(1, 1, 1, 0, 0).unwrap();
        let r = cupdlp_dcs_done(Some(c), 0);
        assert!(r.is_none());
    }

    #[test]
    fn entry_not_triplet() {
        let mut a = cupdlp_dcs_spalloc(1, 1, 1, 1, 0).unwrap();
        assert_eq!(cupdlp_dcs_entry(&mut a, 0, 0, 1.0), 0);
    }

    #[test]
    fn entry_negative_indices() {
        let mut a = cupdlp_dcs_spalloc(0, 0, 1, 1, 1).unwrap();
        assert_eq!(cupdlp_dcs_entry(&mut a, -1, 0, 1.0), 0);
        assert_eq!(cupdlp_dcs_entry(&mut a, 0, -1, 1.0), 0);
    }

    #[test]
    fn compress_not_triplet() {
        let a = cupdlp_dcs_spalloc(1, 1, 1, 1, 0).unwrap();
        assert!(cupdlp_dcs_compress(&a).is_none());
    }

    #[test]
    fn norm_not_csc() {
        let a = cupdlp_dcs_spalloc(0, 0, 1, 1, 1).unwrap();
        assert_eq!(cupdlp_dcs_norm(&a), -1.0);
    }

    #[test]
    fn norm_no_x() {
        let mut a = cupdlp_dcs_spalloc(1, 1, 1, 0, 0).unwrap();
        a.nz = -1;
        a.p = vec![0, 0];
        assert_eq!(cupdlp_dcs_norm(&a), -1.0);
    }

    #[test]
    fn transpose_not_csc() {
        let a = cupdlp_dcs_spalloc(0, 0, 1, 1, 1).unwrap();
        assert!(cupdlp_dcs_transpose(&a, 1).is_none());
    }

    #[test]
    fn transpose_no_values() {
        let mut a = cupdlp_dcs_spalloc(2, 2, 2, 0, 0).unwrap();
        a.nz = -1;
        a.p = vec![0, 1, 2];
        a.i = vec![0, 1];
        let ct = cupdlp_dcs_transpose(&a, 0).unwrap();
        assert!(ct.x.is_none());
    }
}
