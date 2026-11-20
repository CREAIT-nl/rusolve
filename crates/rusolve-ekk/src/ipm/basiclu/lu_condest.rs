//! HiGHS source: ipm/basiclu/lu_condest.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_condest.c::lu_condest (lines 24-49)
///
/// `pivot`, `norm`, `norminv` may be `NULL` in C.
pub fn lu_condest(
    m: LuInt,
    ubegin: &[LuInt],
    ui: &[LuInt],
    ux: &[f64],
    pivot: Option<&[f64]>,
    perm: &[LuInt],
    upper: i32,
    work: &mut [f64],
    norm: Option<&mut f64>,
    norminv: Option<&mut f64>,
) -> f64 {
    let mut unorm: f64 = 0.0;

    /* compute 1-norm of U */
    for j in 0..m {
        let mut colsum: f64 = match pivot {
            Some(pv) => pv[j as usize].abs(),
            None => 1.0,
        };
        let mut p = ubegin[j as usize];
        while ui[p as usize] >= 0 {
            colsum += ux[p as usize].abs();
            p += 1;
        }
        unorm = unorm.max(colsum);
    }

    /* estimate 1-norm of U^{-1} */
    let uinvnorm = lu_normest(m, ubegin, ui, ux, pivot, perm, upper, work);

    if let Some(norm) = norm {
        *norm = unorm;
    }
    if let Some(norminv) = norminv {
        *norminv = uinvnorm;
    }

    unorm * uinvnorm
}

/// HiGHS: basiclu/lu_condest.c::lu_normest (lines 74-124)
///
/// PORT NOTE: C `perm ? perm[k] : k`. Every caller passes `this.p`
/// (basiclu_factorize.c:94,97), never `NULL`, so the identity branch is not
/// reproduced and `perm` is a plain slice.
pub fn lu_normest(
    m: LuInt,
    ubegin: &[LuInt],
    ui: &[LuInt],
    ux: &[f64],
    pivot: Option<&[f64]>,
    perm: &[LuInt],
    upper: i32,
    work: &mut [f64],
) -> f64 {
    let mut x1norm: f64 = 0.0;
    let mut xinfnorm: f64 = 0.0;

    // first sweep: U'x = b, b[i] = +1/-1 chosen dynamically
    let (kbeg, kend, kinc): (LuInt, LuInt, LuInt) = if upper != 0 {
        (0, m, 1)
    } else {
        (m - 1, -1, -1)
    };
    let mut k = kbeg;
    while k != kend {
        let j = perm[k as usize];
        let mut temp: f64 = 0.0;
        let mut p = ubegin[j as usize];
        while ui[p as usize] >= 0 {
            let i = ui[p as usize];
            temp -= work[i as usize] * ux[p as usize];
            p += 1;
        }
        temp += if temp >= 0.0 { 1.0 } else { -1.0 }; /* choose b[i] = 1 or b[i] = -1 */
        if let Some(pv) = pivot {
            temp /= pv[j as usize];
        }
        work[j as usize] = temp;
        x1norm += temp.abs();
        xinfnorm = xinfnorm.max(temp.abs());
        k += kinc;
    }

    // second sweep: Uy = x
    let mut y1norm: f64 = 0.0;
    let (kbeg, kend, kinc): (LuInt, LuInt, LuInt) = if upper != 0 {
        (m - 1, -1, -1)
    } else {
        (0, m, 1)
    };
    let mut k = kbeg;
    while k != kend {
        let j = perm[k as usize];
        if let Some(pv) = pivot {
            work[j as usize] /= pv[j as usize];
        }
        let temp = work[j as usize];
        let mut p = ubegin[j as usize];
        while ui[p as usize] >= 0 {
            let i = ui[p as usize];
            work[i as usize] -= temp * ux[p as usize];
            p += 1;
        }
        y1norm += temp.abs();
        k += kinc;
    }

    (y1norm / x1norm).max(xinfnorm)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// U = [[2, 1], [0, 4]] (upper): off-diagonal entry 1 at (0,1), pivots
    /// (2, 4), identity perm. inv(U) = [[0.5, -0.125], [0, 0.25]] has 1-norm
    /// 0.5 and ||U||_1 = 5. All values are dyadic, so the results are exact.
    #[test]
    fn condest_upper_two_by_two() {
        let ubegin = [0, 1];
        let ui = [-1, 0, -1];
        let ux = [0.0, 1.0, 0.0];
        let pivot = [2.0, 4.0];
        let perm = [0, 1];
        let mut work = [0.0; 2];
        let mut norm = 0.0;
        let mut norminv = 0.0;
        let cond = lu_condest(
            2,
            &ubegin,
            &ui,
            &ux,
            Some(&pivot),
            &perm,
            1,
            &mut work,
            Some(&mut norm),
            Some(&mut norminv),
        );
        assert_eq!(norm, 5.0);
        assert_eq!(norminv, 0.5);
        assert_eq!(cond, 2.5);
    }

    /// L = [[2, 0], [1, 4]] (lower): off-diagonal entry 1 at (1,0) stored in
    /// column 0, pivots (2, 4), identity perm. inv(L) = [[0.5, 0], [-0.125,
    /// 0.25]] has 1-norm 0.625; ||L||_1 = max(2 + 1, 4) = 4.
    #[test]
    fn condest_lower_two_by_two() {
        let ubegin = [0, 2];
        let ui = [1, -1, -1];
        let ux = [1.0, 0.0, 0.0];
        let pivot = [2.0, 4.0];
        let perm = [0, 1];
        let mut work = [0.0; 2];
        let mut norm = 0.0;
        let mut norminv = 0.0;
        let cond = lu_condest(
            2,
            &ubegin,
            &ui,
            &ux,
            Some(&pivot),
            &perm,
            0,
            &mut work,
            Some(&mut norm),
            Some(&mut norminv),
        );
        assert_eq!(norm, 4.0);
        assert_eq!(norminv, 0.625);
        assert_eq!(cond, 2.5);
    }
}
