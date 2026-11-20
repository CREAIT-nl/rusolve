//! Dense BLAS-3-style kernels for the supernodal Cholesky.
//!
//! Provenance:
//! * packed, cache-blocked GEMM with a register-blocked `MR × NR` micro-kernel
//!   (pack an `MC × KC` block of A into `MR`-row micro-panels and a `KC × N`
//!   block of B into `NR`-column micro-panels; the micro-kernel keeps the
//!   `MR × NR` block of C in registers across the whole `KC` loop):
//!   K. Goto & R. van de Geijn, "Anatomy of high-performance matrix
//!   multiplication", ACM TOMS 34(3) (2008); F. Van Zee & R. van de Geijn,
//!   "BLIS: a framework for rapidly instantiating BLAS functionality",
//!   ACM TOMS 41(3) (2015);
//! * blocked right-looking dense Cholesky of a supernode panel (factor an
//!   `NB`-column block, then one rank-`NB` trailing update by the GEMM kernel):
//!   Golub & Van Loan, *Matrix Computations* (4th ed. 2013) §4.2.9 (block
//!   Cholesky); Ng & Peyton, SIAM J. Sci. Comput. 14(5) (1993) for its use on
//!   supernode panels.
//!
//! All matrices are column-major: element (i, p) of a matrix with leading
//! dimension `ld` is at `i + p·ld`. No `unsafe`: the micro-kernel is written
//! on fixed-size arrays so LLVM vectorizes it (NEON `fmla` on aarch64; on
//! x86-64 without the `fma` target feature it falls back to mul + add).

/// Micro-kernel rows (two 4-lane / four 2-lane vectors).
pub const MR: usize = 8;
/// Micro-kernel columns.
pub const NR: usize = 4;
/// Inner-product block (packed micro-panels of A and B stay in L1/L2).
const KC: usize = 256;
/// Row block of A packed at once.
const MC: usize = 128;

#[inline(always)]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(any(target_arch = "aarch64", target_feature = "fma"))]
    {
        a.mul_add(b, c)
    }
    #[cfg(not(any(target_arch = "aarch64", target_feature = "fma")))]
    {
        a * b + c
    }
}

/// `acc += Σ_p ap[p]·bp[p]ᵀ` over `kc` packed steps.
#[inline(always)]
fn micro(kc: usize, ap: &[f64], bp: &[f64]) -> [[f64; MR]; NR] {
    let mut c = [[0.0f64; MR]; NR];
    for (a, b) in ap
        .as_chunks::<MR>()
        .0
        .iter()
        .zip(bp.as_chunks::<NR>().0)
        .take(kc)
    {
        for j in 0..NR {
            let bj = b[j];
            for i in 0..MR {
                c[j][i] = fma(a[i], bj, c[j][i]);
            }
        }
    }
    c
}

/// Reusable packing buffers (one per thread).
#[derive(Default)]
pub struct GemmWork {
    ap: Vec<f64>,
    bp: Vec<f64>,
}

/// `C[i, j] += alpha · Σ_{p<k} A[i, p] · B[j, p]` for `i < m`, `j < n`
/// (C = C + α·A·Bᵀ). With `lower`, tiles lying entirely strictly above the
/// diagonal (`i < j` for every entry) are skipped; entries with `i < j` inside
/// diagonal tiles may still be written, so callers must not rely on them.
#[allow(clippy::too_many_arguments)]
pub fn gemm_nt(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    lower: bool,
    w: &mut GemmWork,
) {
    if m == 0 || n == 0 || k == 0 {
        return;
    }
    let ntile = n.div_ceil(NR);
    for p0 in (0..k).step_by(KC) {
        let kc = KC.min(k - p0);
        // Pack B[0..n, p0..p0+kc] into NR-column micro-panels (zero-padded).
        w.bp.clear();
        w.bp.resize(ntile * kc * NR, 0.0);
        for jt in 0..ntile {
            let j0 = jt * NR;
            let nj = NR.min(n - j0);
            let dst = &mut w.bp[jt * kc * NR..(jt + 1) * kc * NR];
            for p in 0..kc {
                let col = &b[(p0 + p) * ldb + j0..(p0 + p) * ldb + j0 + nj];
                dst[p * NR..p * NR + nj].copy_from_slice(col);
            }
        }
        for i0 in (0..m).step_by(MC) {
            let mc = MC.min(m - i0);
            // Under `lower`, columns j ≥ i0 + mc are entirely above the
            // diagonal for this row block.
            let ncols = if lower { n.min(i0 + mc) } else { n };
            if ncols == 0 {
                continue;
            }
            let mtile = mc.div_ceil(MR);
            w.ap.clear();
            w.ap.resize(mtile * kc * MR, 0.0);
            for it in 0..mtile {
                let r0 = i0 + it * MR;
                let mi = MR.min(i0 + mc - r0);
                let dst = &mut w.ap[it * kc * MR..(it + 1) * kc * MR];
                for p in 0..kc {
                    let col = &a[(p0 + p) * lda + r0..(p0 + p) * lda + r0 + mi];
                    dst[p * MR..p * MR + mi].copy_from_slice(col);
                }
            }
            for jt in 0..ncols.div_ceil(NR) {
                let j0 = jt * NR;
                let nj = NR.min(n - j0);
                let bp = &w.bp[jt * kc * NR..(jt + 1) * kc * NR];
                for it in 0..mtile {
                    let r0 = i0 + it * MR;
                    if lower && r0 + MR <= j0 {
                        continue;
                    }
                    let mi = MR.min(i0 + mc - r0);
                    let ap = &w.ap[it * kc * MR..(it + 1) * kc * MR];
                    let acc = micro(kc, ap, bp);
                    for (jj, accj) in acc.iter().enumerate().take(nj) {
                        let base = (j0 + jj) * ldc + r0;
                        let dst = &mut c[base..base + mi];
                        for (d, v) in dst.iter_mut().zip(accj.iter()) {
                            *d = fma(alpha, *v, *d);
                        }
                    }
                }
            }
        }
    }
}

/// Minimum flops of a trailing panel update before it is split across threads.
const PAR_MIN_FLOPS: f64 = 2e7;

/// [`gemm_nt`] split by column blocks of C over `w.len()` scoped threads
/// (`std::thread::scope`; one packing workspace per thread). Under `lower`,
/// the column cut points balance the lower-trapezoid work Σ_j (m − j), and
/// each thread only touches rows `j0..m` of its columns `j0..j1`.
#[allow(clippy::too_many_arguments)]
pub fn par_gemm_nt(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    lower: bool,
    w: &mut [GemmWork],
) {
    let t = w.len().min(n.div_ceil(NR)).max(1);
    if t == 1 || m == 0 || k == 0 {
        gemm_nt(m, n, k, alpha, a, lda, b, ldb, c, ldc, lower, &mut w[0]);
        return;
    }
    // Column cuts at multiples of NR with equal work per chunk.
    let work = |j: usize| {
        if lower {
            m.saturating_sub(j) as f64
        } else {
            m as f64
        }
    };
    let total: f64 = (0..n).map(work).sum();
    let mut cuts = vec![0usize];
    let mut acc = 0.0;
    for j in 0..n {
        acc += work(j);
        let want = total * cuts.len() as f64 / t as f64;
        if acc >= want && cuts.len() < t && (j + 1) % NR == 0 && j + 1 < n {
            cuts.push(j + 1);
        }
    }
    cuts.push(n);
    std::thread::scope(|sc| {
        let mut rest = c;
        let mut done_cols = 0usize;
        let mut ws = w.iter_mut();
        for win in cuts.windows(2) {
            let (j0, j1) = (win[0], win[1]);
            let tmp = std::mem::take(&mut rest);
            let take = ((j1 - done_cols) * ldc).min(tmp.len());
            let (chunk, tail) = tmp.split_at_mut(take);
            rest = tail;
            done_cols = j1;
            let wk = ws.next().unwrap();
            let (r0, mm) = if lower {
                (j0.min(m), m.saturating_sub(j0))
            } else {
                (0, m)
            };
            if mm == 0 {
                continue;
            }
            let a = &a[r0..];
            let b = &b[j0..];
            let mut job = move || {
                gemm_nt(
                    mm,
                    j1 - j0,
                    k,
                    alpha,
                    a,
                    lda,
                    b,
                    ldb,
                    &mut chunk[r0..],
                    ldc,
                    lower,
                    wk,
                )
            };
            if j1 == n {
                job();
            } else {
                sc.spawn(job);
            }
        }
    });
}

/// Outcome of a pivot in the dense panel Cholesky, decided by the caller.
pub enum Pivot {
    /// Use `sqrt(d)`.
    Accept,
    /// Replace the diagonal of L by this value and zero the column below it
    /// (dependent-row skip, Wright 1999).
    Skip(f64),
    /// Replace the diagonal of L by this value and scale the column by it
    /// (pivot modification: the factored matrix gains `l² − d` on this
    /// diagonal entry; Gill & Murray modified Cholesky).
    Replace(f64),
}

/// Blocked dense Cholesky of an `nr × nc` panel (column-major, ld = `nr`,
/// `nr ≥ nc`): on exit the lower trapezoid holds L (the panel's own `nc`
/// columns), rows `nc..nr` hold `L21 = A21·L11⁻ᵀ`. `pivot(j, d)` decides each
/// pivot. `tick()` is called once per column block and returns false to abort.
/// Columns are processed in blocks of `nb`; the in-block part is a left-looking
/// column update unrolled by 4, the trailing part one [`gemm_nt`] call.
#[allow(clippy::too_many_arguments)]
pub fn panel_cholesky(
    blk: &mut [f64],
    nr: usize,
    nc: usize,
    nb: usize,
    threads: usize,
    w: &mut [GemmWork],
    mut pivot: impl FnMut(usize, f64) -> Pivot,
    mut tick: impl FnMut() -> bool,
) -> bool {
    let nb = nb.max(1);
    let mut j0 = 0;
    while j0 < nc {
        if !tick() {
            return false;
        }
        let jb = nb.min(nc - j0);
        // In-block left-looking factor of columns j0..j0+jb (all rows).
        for j in j0..j0 + jb {
            let (left, right) = blk.split_at_mut(j * nr);
            let colj = &mut right[..nr];
            let mut kk = j0;
            while kk + 4 <= j {
                let c0 = &left[kk * nr..(kk + 1) * nr];
                let c1 = &left[(kk + 1) * nr..(kk + 2) * nr];
                let c2 = &left[(kk + 2) * nr..(kk + 3) * nr];
                let c3 = &left[(kk + 3) * nr..(kk + 4) * nr];
                let (t0, t1, t2, t3) = (c0[j], c1[j], c2[j], c3[j]);
                for r in j..nr {
                    colj[r] -= c0[r] * t0 + c1[r] * t1 + c2[r] * t2 + c3[r] * t3;
                }
                kk += 4;
            }
            while kk < j {
                let ck = &left[kk * nr..(kk + 1) * nr];
                let t = ck[j];
                if t != 0.0 {
                    for r in j..nr {
                        colj[r] -= t * ck[r];
                    }
                }
                kk += 1;
            }
            let d = colj[j];
            match pivot(j, d) {
                Pivot::Accept => {
                    let sd = d.sqrt();
                    colj[j] = sd;
                    let inv = 1.0 / sd;
                    for v in colj[j + 1..].iter_mut() {
                        *v *= inv;
                    }
                }
                Pivot::Skip(h) => {
                    colj[j] = h;
                    colj[j + 1..].iter_mut().for_each(|v| *v = 0.0);
                }
                Pivot::Replace(sd) => {
                    colj[j] = sd;
                    let inv = 1.0 / sd;
                    for v in colj[j + 1..].iter_mut() {
                        *v *= inv;
                    }
                }
            }
        }
        // Trailing update: rows j1..nr, columns j1..nc  -=  L[·, j0..j1]·L[j1..nc, j0..j1]ᵀ.
        let j1 = j0 + jb;
        if j1 < nc {
            let (done, rest) = blk.split_at_mut(j1 * nr);
            let a = &done[j0 * nr + j1..];
            let (mm, nn) = (nr - j1, nc - j1);
            let fl = (2 * jb) as f64 * nn as f64 * (mm as f64 - 0.5 * nn as f64);
            if threads > 1 && fl >= PAR_MIN_FLOPS {
                par_gemm_nt(mm, nn, jb, -1.0, a, nr, a, nr, &mut rest[j1..], nr, true, {
                    let tw = threads.min(w.len());
                    &mut w[..tw]
                });
            } else {
                gemm_nt(
                    mm,
                    nn,
                    jb,
                    -1.0,
                    a,
                    nr,
                    a,
                    nr,
                    &mut rest[j1..],
                    nr,
                    true,
                    &mut w[0],
                );
            }
        }
        j0 = j1;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::XorShift64;

    #[test]
    fn gemm_matches_naive() {
        let mut rng = XorShift64::new(7);
        for &(m, n, k) in &[
            (1, 1, 1),
            (9, 5, 3),
            (37, 13, 300),
            (130, 129, 17),
            (8, 4, 256),
        ] {
            let lda = m + 3;
            let ldb = n + 1;
            let ldc = m + 2;
            let a: Vec<f64> = (0..lda * k).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let b: Vec<f64> = (0..ldb * k).map(|_| rng.uniform(-1.0, 1.0)).collect();
            for lower in [false, true] {
                let c0: Vec<f64> = (0..ldc * n).map(|_| rng.uniform(-1.0, 1.0)).collect();
                let mut c = c0.clone();
                let mut w = GemmWork::default();
                gemm_nt(m, n, k, -1.0, &a, lda, &b, ldb, &mut c, ldc, lower, &mut w);
                for j in 0..n {
                    for i in 0..m {
                        if lower && i < j {
                            continue;
                        }
                        let s: f64 = (0..k).map(|p| a[i + p * lda] * b[j + p * ldb]).sum();
                        let want = c0[i + j * ldc] - s;
                        assert!(
                            (c[i + j * ldc] - want).abs() <= 1e-12 * (1.0 + s.abs()),
                            "m{m} n{n} k{k} ({i},{j})"
                        );
                    }
                }
                // padding rows of C untouched
                for j in 0..n {
                    for i in m..ldc {
                        assert_eq!(c[i + j * ldc], c0[i + j * ldc]);
                    }
                }
            }
        }
    }

    #[test]
    fn par_gemm_matches_serial() {
        let mut rng = XorShift64::new(11);
        for &(m, n, k) in &[(300, 200, 40), (97, 97, 13), (50, 9, 5)] {
            let a: Vec<f64> = (0..m * k).map(|_| rng.uniform(-1.0, 1.0)).collect();
            for lower in [false, true] {
                let c0: Vec<f64> = (0..m * n).map(|_| rng.uniform(-1.0, 1.0)).collect();
                let (mut c1, mut c2) = (c0.clone(), c0.clone());
                gemm_nt(
                    m,
                    n,
                    k,
                    -1.0,
                    &a,
                    m,
                    &a,
                    m,
                    &mut c1,
                    m,
                    lower,
                    &mut GemmWork::default(),
                );
                let mut ws: Vec<GemmWork> = (0..4).map(|_| GemmWork::default()).collect();
                par_gemm_nt(m, n, k, -1.0, &a, m, &a, m, &mut c2, m, lower, &mut ws);
                for j in 0..n {
                    for i in 0..m {
                        if lower && i < j {
                            continue;
                        }
                        assert_eq!(c1[i + j * m], c2[i + j * m], "m{m} n{n} ({i},{j})");
                    }
                }
            }
        }
    }

    #[test]
    fn panel_cholesky_matches_unblocked() {
        let mut rng = XorShift64::new(3);
        for &(nr, nc, nb) in &[(5, 5, 2), (40, 17, 4), (100, 70, 16), (64, 64, 64)] {
            // SPD leading block: B·Bᵀ + nr·I, rectangular rows below.
            let bm: Vec<f64> = (0..nr * nc).map(|_| rng.uniform(-1.0, 1.0)).collect();
            let mut p = vec![0.0; nr * nc];
            for j in 0..nc {
                for i in j..nr {
                    let s: f64 = (0..nc).map(|q| bm[i + q * nr] * bm[j + q * nr]).sum();
                    p[i + j * nr] = s + if i == j { nr as f64 } else { 0.0 };
                }
            }
            let orig = p.clone();
            let mut w = vec![GemmWork::default()];
            assert!(panel_cholesky(
                &mut p,
                nr,
                nc,
                nb,
                1,
                &mut w,
                |_, _| Pivot::Accept,
                || true
            ));
            // check L·Lᵀ = original on the lower trapezoid
            for j in 0..nc {
                for i in j..nr {
                    let s: f64 = (0..=j).map(|q| p[i + q * nr] * p[j + q * nr]).sum();
                    assert!(
                        (s - orig[i + j * nr]).abs() <= 1e-10 * (1.0 + orig[i + j * nr].abs()),
                        "nr{nr} nc{nc} nb{nb} ({i},{j})"
                    );
                }
            }
        }
    }
}
