//! HiGHS source: extern/metis/libmetis/mcutil.c, extern/metis/libmetis/util.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Multi-constraint balance helpers and the `iargmax*_nrm` utilities. With
//! `ncon == 1` (the only HiPO case) several are reached only from the
//! `Mc*` paths, but `ivecle` / `ivecaxpylez` / `BetterBalance2Way` /
//! `ComputeLoadImbalanceDiff` are on the single-constraint path too.

use super::{graph_t, idx_t, real_t};

/// HiGHS: mcutil.c::ivecle (lines 20-28)
pub fn ivecle(n: idx_t, x: &[idx_t], z: &[idx_t]) -> i32 {
    // HiGHS: mcutil.c lines 22-25
    // C: for (n--; n>=0; n--) { if (x[n] > z[n]) return 0; }
    let mut n = n;
    n -= 1;
    while n >= 0 {
        if x[n as usize] > z[n as usize] {
            return 0;
        }
        n -= 1;
    }

    // HiGHS: mcutil.c line 27
    1
}

/// HiGHS: mcutil.c::ivecaxpylez (lines 33-41)
pub fn ivecaxpylez(n: idx_t, a: idx_t, x: &[idx_t], y: &[idx_t], z: &[idx_t]) -> i32 {
    // HiGHS: mcutil.c lines 35-38
    // C: for (n--; n>=0; n--) { if (a*x[n]+y[n] > z[n]) return 0; }
    let mut n = n;
    n -= 1;
    while n >= 0 {
        // C: a*x[n]+y[n] > z[n]
        // All are idx_t (i32). C integer arithmetic.
        if a * x[n as usize] + y[n as usize] > z[n as usize] {
            return 0;
        }
        n -= 1;
    }

    // HiGHS: mcutil.c line 40
    1
}

/// HiGHS: mcutil.c::BetterVBalance (lines 47-66)
pub fn BetterVBalance(
    ncon: idx_t,
    invtvwgt: &[real_t],
    v_vwgt: &[idx_t],
    u1_vwgt: &[idx_t],
    u2_vwgt: &[idx_t],
) -> i32 {
    // HiGHS: mcutil.c lines 50-51
    // real_t sum1=0.0, sum2=0.0, diff1=0.0, diff2=0.0;
    let mut sum1: real_t = 0.0;
    let mut sum2: real_t = 0.0;
    let mut diff1: real_t = 0.0;
    let mut diff2: real_t = 0.0;

    // HiGHS: mcutil.c lines 53-56
    // for (i=0; i<ncon; i++) {
    //   sum1 += (v_vwgt[i]+u1_vwgt[i])*invtvwgt[i];
    //   sum2 += (v_vwgt[i]+u2_vwgt[i])*invtvwgt[i];
    // }
    // C: (v_vwgt[i]+u1_vwgt[i]) is idx_t + idx_t = idx_t (i32).
    // Then * invtvwgt[i] (real_t = f32). C promotes i32 to f32, multiplies in f32.
    // sum1 += result in f32.
    for i in 0..ncon {
        let t1: f32 = (v_vwgt[i as usize] + u1_vwgt[i as usize]) as f32 * invtvwgt[i as usize];
        sum1 += t1;
        let t2: f32 = (v_vwgt[i as usize] + u2_vwgt[i as usize]) as f32 * invtvwgt[i as usize];
        sum2 += t2;
    }

    // HiGHS: mcutil.c lines 57-58
    // sum1 = sum1/ncon;
    // sum2 = sum2/ncon;
    // C: real_t / idx_t converts ncon to f32 and divides in f32.
    sum1 = sum1 / (ncon as f32);
    sum2 = sum2 / (ncon as f32);

    // HiGHS: mcutil.c lines 60-63
    // for (i=0; i<ncon; i++) {
    //   diff1 += rabs(sum1 - (v_vwgt[i]+u1_vwgt[i])*invtvwgt[i]);
    //   diff2 += rabs(sum2 - (v_vwgt[i]+u2_vwgt[i])*invtvwgt[i]);
    // }
    // rabs is fabsf (REALTYPEWIDTH == 32), so the whole update stays in f32.
    for i in 0..ncon {
        let t1: f32 = (v_vwgt[i as usize] + u1_vwgt[i as usize]) as f32 * invtvwgt[i as usize];
        diff1 += (sum1 - t1).abs();

        let t2: f32 = (v_vwgt[i as usize] + u2_vwgt[i as usize]) as f32 * invtvwgt[i as usize];
        diff2 += (sum2 - t2).abs();
    }

    // HiGHS: mcutil.c line 65
    // return (diff1 - diff2 >= 0);
    // C: diff1 - diff2 is f32. >= 0 compares f32 with int (0 promoted to f32).
    // Result is int (0 or 1).
    (diff1 - diff2 >= 0.0) as i32
}

/// HiGHS: mcutil.c::BetterBalance2Way (lines 72-81)
pub fn BetterBalance2Way(n: idx_t, x: &[real_t], y: &[real_t]) -> i32 {
    // HiGHS: mcutil.c line 74
    let mut nrm1: real_t = 0.0;
    let mut nrm2: real_t = 0.0;

    // HiGHS: mcutil.c lines 76-79
    // for (--n; n>=0; n--) {
    //   if (x[n] > 0) nrm1 += x[n]*x[n];
    //   if (y[n] > 0) nrm2 += y[n]*y[n];
    // }
    let mut n = n;
    n -= 1;
    while n >= 0 {
        if x[n as usize] > 0.0 {
            nrm1 += x[n as usize] * x[n as usize];
        }
        if y[n as usize] > 0.0 {
            nrm2 += y[n as usize] * y[n as usize];
        }
        n -= 1;
    }

    // HiGHS: mcutil.c line 80
    // return nrm2 < nrm1;
    (nrm2 < nrm1) as i32
}

/// HiGHS: mcutil.c::ComputeLoadImbalance (lines 87-105)
pub fn ComputeLoadImbalance(graph: &graph_t, nparts: idx_t, pijbm: &[real_t]) -> real_t {
    // HiGHS: mcutil.c lines 89-93
    let ncon = graph.ncon;
    let pwgts = &graph.pwgts;

    // HiGHS: mcutil.c line 95
    let mut max: real_t = 1.0;

    // HiGHS: mcutil.c lines 96-102
    for i in 0..ncon {
        for j in 0..nparts {
            // C: cur = pwgts[j*ncon+i]*pijbm[j*ncon+i];
            // pwgts[...] is idx_t (i32). pijbm[...] is real_t (f32).
            // C: i32 * f32 -> f32.
            let idx = (j * ncon + i) as usize;
            let cur: real_t = (pwgts[idx] as f32) * pijbm[idx];
            if cur > max {
                max = cur;
            }
        }
    }

    // HiGHS: mcutil.c line 104
    max
}

/// HiGHS: mcutil.c::ComputeLoadImbalanceDiff (lines 114-133)
pub fn ComputeLoadImbalanceDiff(
    graph: &graph_t,
    nparts: idx_t,
    pijbm: &[real_t],
    ubvec: &[real_t],
) -> real_t {
    // HiGHS: mcutil.c lines 117-121
    let ncon = graph.ncon;
    let pwgts = &graph.pwgts;

    // HiGHS: mcutil.c line 123
    let mut max: real_t = -1.0;

    // HiGHS: mcutil.c lines 124-130
    for i in 0..ncon {
        for j in 0..nparts {
            // C: cur = pwgts[j*ncon+i]*pijbm[j*ncon+i] - ubvec[i];
            let idx = (j * ncon + i) as usize;
            let cur: real_t = (pwgts[idx] as f32) * pijbm[idx] - ubvec[i as usize];
            if cur > max {
                max = cur;
            }
        }
    }

    // HiGHS: mcutil.c line 132
    max
}

/// HiGHS: mcutil.c::ComputeLoadImbalanceDiffVec (lines 141-162)
pub fn ComputeLoadImbalanceDiffVec(
    graph: &graph_t,
    nparts: idx_t,
    pijbm: &[real_t],
    ubfactors: &[real_t],
    diffvec: &mut [real_t],
) -> real_t {
    // HiGHS: mcutil.c lines 144-148
    let ncon = graph.ncon;
    let pwgts = &graph.pwgts;

    // HiGHS: mcutil.c line 150
    // for (max=-1.0, i=0; i<ncon; i++) {
    let mut max: real_t = -1.0;
    for i in 0..ncon {
        // HiGHS: mcutil.c line 151
        // diffvec[i] = pwgts[i]*pijbm[i] - ubfactors[i];
        diffvec[i as usize] =
            (pwgts[i as usize] as f32) * pijbm[i as usize] - ubfactors[i as usize];

        // HiGHS: mcutil.c lines 152-156
        for j in 1..nparts {
            // C: cur = pwgts[j*ncon+i]*pijbm[j*ncon+i] - ubfactors[i];
            let idx = (j * ncon + i) as usize;
            let cur: real_t = (pwgts[idx] as f32) * pijbm[idx] - ubfactors[i as usize];
            if cur > diffvec[i as usize] {
                diffvec[i as usize] = cur;
            }
        }

        // HiGHS: mcutil.c lines 157-158
        if max < diffvec[i as usize] {
            max = diffvec[i as usize];
        }
    }

    // HiGHS: mcutil.c line 161
    max
}

/// HiGHS: util.c::iargmax_nrm (lines 20-28)
pub fn iargmax_nrm(n: usize, x: &[idx_t], y: &[real_t]) -> idx_t {
    // HiGHS: util.c lines 22-23
    let mut max: idx_t = 0;

    // HiGHS: util.c lines 24-25
    // for (i=1; i<n; i++)
    //   max = (x[i]*y[i] > x[max]*y[max] ? i : max);
    // C: x[i] is idx_t (i32), y[i] is real_t (f32).
    // x[i]*y[i]: i32 * f32 -> f32.
    // Comparison in f32.
    for i in 1..n {
        let lhs: f32 = (x[i] as f32) * y[i];
        let rhs: f32 = (x[max as usize] as f32) * y[max as usize];
        if lhs > rhs {
            max = i as idx_t;
        }
    }

    // HiGHS: util.c line 27
    max
}

/// HiGHS: util.c::iargmax2_nrm (lines 34-57)
pub fn iargmax2_nrm(n: usize, x: &[idx_t], y: &[real_t]) -> idx_t {
    // HiGHS: util.c lines 36-45
    let mut max1: usize;
    let mut max2: usize;

    // C: if (x[0]*y[0] > x[1]*y[1])
    let v0: f32 = (x[0] as f32) * y[0];
    let v1: f32 = (x[1] as f32) * y[1];
    if v0 > v1 {
        max1 = 0;
        max2 = 1;
    } else {
        max1 = 1;
        max2 = 0;
    }

    // HiGHS: util.c lines 47-54
    for i in 2..n {
        let vi: f32 = (x[i] as f32) * y[i];
        let v_max1: f32 = (x[max1] as f32) * y[max1];
        if vi > v_max1 {
            max2 = max1;
            max1 = i;
        } else {
            let v_max2: f32 = (x[max2] as f32) * y[max2];
            if vi > v_max2 {
                max2 = i;
            }
        }
    }

    // HiGHS: util.c line 56
    max2 as idx_t
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn ivecle_basic() {
        assert_eq!(ivecle(2, &[1, 2], &[1, 3]), 1);
        assert_eq!(ivecle(2, &[1, 4], &[1, 3]), 0);
    }
}
