//! Row-level convexity test and factor for the SOCP / convex QCQP route:
//! `Q = L Lᵀ` for the `½ xᵀQx` of a quadratic constraint row, so the row
//! `aᵀx + ½ xᵀQx <= b` becomes the rotated cone `‖Lᵀx‖² <= 2 (b - aᵀx)`.
//!
//! Dense complete-pivoting Cholesky on the row's support (the variables with
//! a nonzero entry). It stops at the numerical rank, so a singular PSD `Q`
//! (the usual case: a few squared terms) gives a thin `L`, and it returns the
//! variable of negative curvature when `Q` is not PSD. Source: Higham 2002,
//! *Accuracy and Stability of Numerical Algorithms*, 2nd ed., section 10.3
//! (Cholesky with complete pivoting for semidefinite matrices).
//!
//! A row whose support is larger than [`MAX_SUPPORT`] is refused with a typed
//! error: the sparse factorization the design names is a later step.
//!
//! TODO: unify with the objective convexity check (`rusolve-core`
//! `qp_convexity`). This module is kept separate for now.

use rusolve_core::Triplet;

/// Largest support (variables with a `Q` entry) the dense factor accepts.
pub const MAX_SUPPORT: usize = 3000;

/// `Q = L Lᵀ` on `vars`: `l[c][k]` is the coefficient of `x_{vars[k]}` in
/// the `c`-th component of `Lᵀx`.
#[derive(Debug, Clone, PartialEq)]
pub struct Factor {
    pub vars: Vec<usize>,
    /// `rank` columns of `L`, each as a dense vector over `vars`.
    pub l: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Refusal {
    /// `Q` has negative curvature: `x_a` (and `x_b` when the pivot is an
    /// off-diagonal block) with the Schur complement value `pivot`.
    Indefinite {
        var: usize,
        other: Option<usize>,
        pivot: f64,
    },
    TooLarge {
        support: usize,
    },
}

/// Factors `q` (lower triangle `(i, j, v)`, `i >= j`, the `Q` of `½ xᵀQx`),
/// or says why it cannot. `sign` multiplies `Q` first (`-1.0` to test a
/// concave row).
pub fn factor_psd(q: &[Triplet], sign: f64) -> Result<Factor, Refusal> {
    let mut vars: Vec<usize> = q.iter().flat_map(|&(i, j, _)| [i, j]).collect();
    vars.sort_unstable();
    vars.dedup();
    let k = vars.len();
    if k > MAX_SUPPORT {
        return Err(Refusal::TooLarge { support: k });
    }
    let pos = |v: usize| vars.binary_search(&v).unwrap();
    let mut s = vec![0.0f64; k * k];
    let mut scale = 0.0f64;
    for &(i, j, v) in q {
        let (a, b) = (pos(i), pos(j));
        let v = sign * v;
        s[a * k + b] += v;
        if a != b {
            s[b * k + a] += v;
        }
        scale = scale.max(v.abs());
    }
    // Entries below `tol` are rounding noise of the elimination.
    let tol = 1e-9 * scale.max(f64::MIN_POSITIVE);
    let mut l: Vec<Vec<f64>> = Vec::new();
    let mut active: Vec<bool> = vec![true; k];
    loop {
        // complete pivoting: the largest remaining diagonal
        let (mut p, mut dp) = (usize::MAX, f64::NEG_INFINITY);
        let (mut lo, mut lo_at) = (f64::INFINITY, 0usize);
        for i in 0..k {
            if !active[i] {
                continue;
            }
            let d = s[i * k + i];
            if d > dp {
                (p, dp) = (i, d);
            }
            if d < lo {
                (lo, lo_at) = (d, i);
            }
        }
        if p == usize::MAX {
            break;
        }
        if lo < -tol {
            return Err(Refusal::Indefinite {
                var: vars[lo_at],
                other: None,
                pivot: lo,
            });
        }
        if dp <= tol {
            // Remaining diagonal is zero: a PSD Schur complement is then
            // zero, so any off-diagonal entry is negative curvature.
            for i in 0..k {
                for j in 0..i {
                    if active[i] && active[j] && s[i * k + j].abs() > tol {
                        return Err(Refusal::Indefinite {
                            var: vars[i],
                            other: Some(vars[j]),
                            pivot: -s[i * k + j].abs(),
                        });
                    }
                }
            }
            break;
        }
        let root = dp.sqrt();
        let col: Vec<f64> = (0..k)
            .map(|i| if active[i] { s[i * k + p] / root } else { 0.0 })
            .collect();
        for i in 0..k {
            if !active[i] {
                continue;
            }
            for j in 0..k {
                if active[j] {
                    s[i * k + j] -= col[i] * col[j];
                }
            }
        }
        active[p] = false;
        l.push(col);
    }
    Ok(Factor { vars, l })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rebuild(f: &Factor) -> Vec<Vec<f64>> {
        let k = f.vars.len();
        let mut m = vec![vec![0.0; k]; k];
        for c in &f.l {
            for i in 0..k {
                for j in 0..k {
                    m[i][j] += c[i] * c[j];
                }
            }
        }
        m
    }

    #[test]
    fn psd_singular_gives_a_thin_factor_that_reproduces_q() {
        // Q = [[2, 2], [2, 2]] (rank 1), plus x7²: Q77 = 3.
        let f = factor_psd(&[(0, 0, 2.0), (1, 0, 2.0), (1, 1, 2.0), (7, 7, 3.0)], 1.0).unwrap();
        assert_eq!(f.vars, vec![0, 1, 7]);
        assert_eq!(f.l.len(), 2);
        let m = rebuild(&f);
        let want = [[2.0, 2.0, 0.0], [2.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        for i in 0..3 {
            for j in 0..3 {
                assert!((m[i][j] - want[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn negative_diagonal_names_the_variable() {
        let e = factor_psd(&[(0, 0, 1.0), (3, 3, -2.0)], 1.0).unwrap_err();
        assert_eq!(
            e,
            Refusal::Indefinite {
                var: 3,
                other: None,
                pivot: -2.0
            }
        );
    }

    #[test]
    fn zero_diagonal_with_off_diagonal_is_indefinite() {
        // x y has eigenvalues +-1/2
        let e = factor_psd(&[(1, 0, 1.0)], 1.0).unwrap_err();
        assert!(matches!(e, Refusal::Indefinite { other: Some(0), .. }));
    }

    #[test]
    fn indefinite_after_elimination_is_found() {
        // [[1, 2], [2, 1]]: first pivot 1, Schur 1 - 4 = -3
        let e = factor_psd(&[(0, 0, 1.0), (1, 0, 2.0), (1, 1, 1.0)], 1.0).unwrap_err();
        assert!(matches!(e, Refusal::Indefinite { pivot, .. } if (pivot + 3.0).abs() < 1e-12));
    }

    #[test]
    fn sign_flip_tests_concave_rows() {
        assert!(factor_psd(&[(0, 0, -4.0)], -1.0).is_ok());
        assert!(factor_psd(&[(0, 0, -4.0)], 1.0).is_err());
    }
}
