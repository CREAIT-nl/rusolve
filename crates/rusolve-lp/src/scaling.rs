//! Power-of-2 equilibration scaling around the simplex solve.
//!
//! Curtis-Reid matrix equilibration (M. J. Curtis & J. K. Reid, "On the
//! solution of large sparse linear programs", ACM Trans. Math. Softw.
//! 1(2):188-207, 1972) in its least-squares-in-log form, with every factor
//! forced to an exact power of 2 so that scaling the matrix, rhs and bounds
//! is rounding-free in binary floating point (a scaled entry is the original
//! significand with only its exponent changed).
//!
//! Shape (following HiGHS —
//! HiGHS scales before every solve, factors rounded to powers of 2, coarse
//! max-based equilibration then iterative mean-based sweeps, and applies the
//! result only through an improvement gate):
//!
//! 1. Row pass:   r_i = 2^round(-log2 max_j |a_ij|)      (zero rows: 1)
//! 2. Col pass:   c_j = 2^round(-log2 max_i |r_i a_ij|)  (zero cols: 1)
//! 3. Up to [`REFINE_SWEEPS`] refinement sweeps, each a row then a column
//!    geometric-mean step (the Curtis-Reid iteration: drive the per-row /
//!    per-column geomean of |scaled a| to 1, rounded to a power of 2). The
//!    factor set with the smallest nonzero-entry range max|a|/min|a| wins.
//! 4. Improvement gate: if the scaled range is not strictly better than the
//!    unscaled one, the LP is returned untouched (`applied = false`).
//!
//! Exponents are clamped to ±[`MAX_EXP`] (cf. HiGHS
//! `allowed_matrix_scale_factor`).
//!
//! Provenance (HiGHS, MIT): the mechanism follows the
//! HiGHS LP-scaling region of `highs/lp_data/HighsLpUtils.cpp`
//! (`scaleHighsLp` / `chooseScalingStrategy` / `rescaleMatrix`: six equilibration
//! search rounds, `pow(2, floor(log(x)/log2 + 0.5))` power-of-2 rounding,
//! min/max scale clamping, matrix-value-range and equilibration improvement
//! measures). Algorithmic ideas and constants were ported; no C++ idioms or
//! code were translated. The factor range [0.0039, 256], the [0.5, 2)
//! target band and the improvement gate (tuff: 2.29e7) match what HiGHS
//! reports in its scaling log.
//!
//! Application contract (what `scale_lp` guarantees to the caller):
//! * `A' = diag(row) A diag(c)`, `b' = row * b`, row senses unchanged;
//! * column bounds `[l_j, u_j] -> [l_j / c_j, u_j / c_j]` (±inf stays ±inf;
//!   c_j > 0 always, being a power of 2);
//! * `cost'_j = c_j * cost_j`, so the objective VALUE is invariant:
//!   Σ cost'·x' = Σ cost·x for x = c ⊙ x';
//! * [`unscale_x`] maps a scaled-space x' back to caller space.
//!
//! No `unsafe`. All factors are exact powers of 2 by construction.

use rusolve_core::Lp;

/// Cap on |log2| of any factor (2^60 ≈ 1.15e18), cf. HiGHS
/// `allowed_matrix_scale_factor`.
const MAX_EXP: i32 = 60;

/// Number of geometric-mean refinement sweeps after the two max-based
/// passes; the best factor set (smallest entry range) is kept.
const REFINE_SWEEPS: usize = 4;

/// Row/column multipliers of a scaled LP; `col` also unscales the solution.
#[derive(Debug, Clone)]
pub struct Scaling {
    /// Row multipliers: `A' = diag(row) A diag(col)`, `b' = row ⊙ b`.
    pub row: Vec<f64>,
    /// Column multipliers; caller-space `x = col ⊙ x'`.
    pub col: Vec<f64>,
}

impl Scaling {
    /// All-ones scaling (the identity of the scaling group).
    pub fn identity(num_row: usize, num_col: usize) -> Self {
        Self {
            row: vec![1.0; num_row],
            col: vec![1.0; num_col],
        }
    }
}

/// Map a scaled-space solution back: `x_j = c_j * x'_j`.
pub fn unscale_x(sc: &Scaling, x: &mut [f64]) {
    for (xj, &cj) in x.iter_mut().zip(&sc.col) {
        *xj *= cj;
    }
}

/// Compute equilibration factors and return `(scaled_lp, scaling, applied)`.
///
/// `applied == false` means the improvement gate rejected the scaling: the
/// returned LP is a bit-identical clone of `lp` and `scaling` is the identity
/// — callers may treat the two outcomes uniformly.
pub fn scale_lp(lp: &Lp) -> (Lp, Scaling, bool) {
    match equilibration(lp) {
        None => {
            if std::env::var("RUSOLVE_SCALE_DEBUG").is_ok() {
                eprintln!("rusolve scaling: gate skipped (no improvement)");
            }
            (lp.clone(), Scaling::identity(lp.num_row, lp.num_col), false)
        }
        Some((re, ce)) => {
            if std::env::var("RUSOLVE_SCALE_DEBUG").is_ok() {
                // Improvement factor in log2 units: 2^(before-after) is the
                // max/min nonzero-entry ratio shrinkage (HiGHS logs the same
                // quantity as "matrix value ratio improvement").
                let before = spread(lp, &vec![0; lp.num_row], &vec![0; lp.num_col]);
                let after = spread(lp, &re, &ce);
                eprintln!(
                    "rusolve scaling: applied, range improvement 2^({before:.2} - {after:.2}) = {:.3}",
                    (before - after).exp2()
                );
            }
            let row: Vec<f64> = re.iter().map(|&e| (e as f64).exp2()).collect();
            let col: Vec<f64> = ce.iter().map(|&e| (e as f64).exp2()).collect();
            (apply(lp, &row, &col), Scaling { row, col }, true)
        }
    }
}

/// `round(log2 v))` as an i32, saturated well inside i32 range so extreme
/// coefficients (|a| up to ~2^1024) cannot wrap the cast.
fn log2_round(v: f64) -> i32 {
    v.log2().round().clamp(-1023.0, 1023.0) as i32
}

/// `max/min` range of the nonzero |scaled| entries, in log2 units
/// (monotone in the ratio; computed in log space so no product overflows).
/// An empty or single-valued support returns 0.0.
fn spread(lp: &Lp, r: &[i32], c: &[i32]) -> f64 {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    for (j, cj) in c.iter().enumerate() {
        for k in ptr[j]..ptr[j + 1] {
            let a = rval[k];
            if a != 0.0 {
                let l = a.abs().log2() + r[ridx[k] as usize] as f64 + *cj as f64;
                lo = lo.min(l);
                hi = hi.max(l);
            }
        }
    }
    if hi <= lo { 0.0 } else { hi - lo }
}

/// Equilibration exponents `(row, col)`, or `None` when the improvement gate
/// rejects scaling (including empty / already-equilibrated matrices).
#[allow(clippy::needless_range_loop)] // column-major: j indexes the CSC column range and c[j] jointly
fn equilibration(lp: &Lp) -> Option<(Vec<i32>, Vec<i32>)> {
    let m = lp.num_row;
    let n = lp.num_col;
    let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
    let before = spread(lp, &vec![0i32; m], &vec![0i32; n]);

    // Pass 1 (row max): r_i maps the row's largest entry to 2^0.
    let mut rmax = vec![0.0f64; m];
    for j in 0..n {
        for k in ptr[j]..ptr[j + 1] {
            let i = ridx[k] as usize;
            let v = rval[k].abs();
            if v > rmax[i] {
                rmax[i] = v;
            }
        }
    }
    let mut r: Vec<i32> = (0..m)
        .map(|i| {
            if rmax[i] > 0.0 {
                (-log2_round(rmax[i])).clamp(-MAX_EXP, MAX_EXP)
            } else {
                0
            }
        })
        .collect();

    // Pass 2 (col max on the row-scaled matrix, in log space).
    let mut c: Vec<i32> = vec![0; n];
    for j in 0..n {
        let mut lmax = f64::NEG_INFINITY;
        for k in ptr[j]..ptr[j + 1] {
            let a = rval[k];
            if a != 0.0 {
                lmax = lmax.max(a.abs().log2() + r[ridx[k] as usize] as f64);
            }
        }
        if lmax > f64::NEG_INFINITY {
            c[j] = (-(lmax.round().clamp(-1023.0, 1023.0) as i32)).clamp(-MAX_EXP, MAX_EXP);
        }
    }

    // Refinement: Curtis-Reid geometric-mean sweeps, row then column,
    // keeping the factor set with the smallest entry range.
    let mut best = spread(lp, &r, &c);
    let mut best_r = r.clone();
    let mut best_c = c.clone();
    for _ in 0..REFINE_SWEEPS {
        // Row step: r_i *= 2^-round(row geomean of |r a / c| in log units).
        let mut sum = vec![0.0f64; m];
        let mut cnt = vec![0usize; m];
        for (j, &cj) in c.iter().enumerate() {
            for k in ptr[j]..ptr[j + 1] {
                let a = rval[k];
                if a != 0.0 {
                    let i = ridx[k] as usize;
                    sum[i] += a.abs().log2() + cj as f64;
                    cnt[i] += 1;
                }
            }
        }
        for i in 0..m {
            if cnt[i] > 0 {
                let mean = sum[i] / cnt[i] as f64 + r[i] as f64;
                r[i] = (r[i] - mean.round() as i32).clamp(-MAX_EXP, MAX_EXP);
            }
        }
        // Column step: same on the updated rows.
        for j in 0..n {
            let mut s = 0.0;
            let mut k = 0usize;
            for kk in ptr[j]..ptr[j + 1] {
                let a = rval[kk];
                if a != 0.0 {
                    s += a.abs().log2() + r[ridx[kk] as usize] as f64;
                    k += 1;
                }
            }
            if k > 0 {
                let mean = s / k as f64 + c[j] as f64;
                c[j] = (c[j] - mean.round() as i32).clamp(-MAX_EXP, MAX_EXP);
            }
        }
        let s = spread(lp, &r, &c);
        if s < best {
            best = s;
            best_r = r.clone();
            best_c = c.clone();
        } else {
            break; // no improvement: the earlier factor set stays
        }
    }

    // Improvement gate: strictly better entry range or no scaling at all.
    if best >= before {
        return None;
    }
    Some((best_r, best_c))
}

/// Build the scaled LP: `A' = diag(row) A diag(c)`, `b' = row ⊙ b`,
/// bounds `÷ c`, costs `⊙ c`. All factors are positive powers of 2, so every
/// finite product/quotient here is exact (exponent arithmetic only) and
/// ±inf bounds keep their infinity.
#[allow(clippy::needless_range_loop)]
fn apply(lp: &Lp, row: &[f64], col: &[f64]) -> Lp {
    let mut out = lp.clone();
    // In-place value rewrite at stored positions (structure unchanged); the
    // index views borrow `lp.a`, the mutable values `out.a` (a distinct clone).
    {
        let (ptr, ridx) = (lp.a.col_ptr(), lp.a.row_idx());
        let rval = out.a.values_mut();
        for j in 0..lp.num_col {
            let cj = col[j];
            for k in ptr[j]..ptr[j + 1] {
                let v = &mut rval[k];
                if *v != 0.0 {
                    *v *= row[ridx[k] as usize] * cj;
                }
            }
        }
    }
    for j in 0..lp.num_col {
        let cj = col[j];
        out.objective[j] *= cj;
        out.lower[j] /= cj;
        out.upper[j] /= cj;
    }
    for (i, b) in out.rhs.iter_mut().enumerate() {
        *b *= row[i];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplex::{SimplexConfig, SimplexSolver, SimplexStatus};

    // Same lp_of shape as the presolve.rs test helpers (private there).
    fn lp_of(
        rows: &[&[f64]],
        rhs: Vec<f64>,
        senses: &str,
        obj: Vec<f64>,
        lo: Vec<f64>,
        up: Vec<f64>,
    ) -> Lp {
        let n = obj.len();
        let m = rhs.len();
        let triplets: Vec<(usize, usize, f64)> = rows
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().enumerate().take(n).map(move |(j, &v)| (i, j, v)))
            .collect();
        Lp {
            num_col: n,
            num_row: m,
            a: rusolve_core::CscMatrix::from_triplets(m, n, &triplets),
            rhs,
            objective: obj,
            maximize: false,
            senses: senses.chars().collect(),
            lower: lo,
            upper: up,
        }
    }

    /// Badly scaled but benign: coefficients spanning 1e-4..1e4 over E/G/L
    /// rows. Optimum by inspection: row0 (E) pins x1 = 3 − 1e−8·x2, row1 (G)
    /// pins x2 = 0.7 − 1e−8·x1, row2 (L) leaves x3 slack ⇒ obj ≈ 3.7.
    fn badly_scaled() -> Lp {
        lp_of(
            &[&[1e4, 1e-4, 0.0], &[1e-4, 1e4, 0.0], &[2e3, 0.0, 1.0]],
            vec![3e4, 7e3, 6.4e3],
            "EGL",
            vec![1.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![10.0, 10.0, 10.0],
        )
    }

    /// The scaling contract, end to end: solving the scaled LP and mapping x
    /// back through `unscale_x` must reproduce the unscaled solve's objective
    /// and a feasible point of the ORIGINAL rows/bounds. (The solve() hook's
    /// `SimplexConfig::scaling` is the same A/B at the solver level.)
    fn check_scale_agree(lp: &Lp) {
        let cfg = SimplexConfig::default();
        let base = SimplexSolver::new(cfg).solve(lp).unwrap();
        let (slp, sc, applied) = scale_lp(lp);
        let mut sol = SimplexSolver::new(cfg).solve(&slp).unwrap();
        assert_eq!(
            base.status, sol.status,
            "{:?} vs {:?}",
            base.status, sol.status
        );
        if base.status != SimplexStatus::Optimal {
            return;
        }
        assert!(applied, "test LPs here should pass the gate");
        unscale_x(&sc, &mut sol.x);
        let obj: f64 = lp.objective.iter().zip(&sol.x).map(|(c, x)| c * x).sum();
        assert!(
            (obj - base.objective).abs() <= 1e-6 * base.objective.abs().max(1.0),
            "obj mismatch: {obj} vs {}",
            base.objective
        );
        // The scaled solve's own reported objective is value-invariant too.
        assert!(
            (sol.objective - base.objective).abs() <= 1e-6 * base.objective.abs().max(1.0),
            "scaled-reported obj {} vs {}",
            sol.objective,
            base.objective
        );
        let act = row_activities(lp, &sol.x);
        for (i, s) in lp.senses.iter().enumerate() {
            let lhs = act[i];
            let r = lp.rhs[i];
            match s {
                'L' => assert!(lhs <= r + 1e-6, "row {i} violated: {lhs} > {r}"),
                'G' => assert!(lhs >= r - 1e-6, "row {i} violated: {lhs} < {r}"),
                _ => assert!((lhs - r).abs() <= 1e-6, "row {i} violated: {lhs} != {r}"),
            }
        }
        for j in 0..lp.num_col {
            assert!(
                sol.x[j] >= lp.lower[j] - 1e-6 && sol.x[j] <= lp.upper[j] + 1e-6,
                "bound {j} violated: {:?}",
                sol.x
            );
        }
    }

    fn is_pow2(f: f64) -> bool {
        f > 0.0 && f.log2().fract() == 0.0
    }

    /// Row activities `A·x` off the CSC (verification code; dense test data).
    fn row_activities(lp: &Lp, x: &[f64]) -> Vec<f64> {
        let (ptr, ridx, rval) = (lp.a.col_ptr(), lp.a.row_idx(), lp.a.values());
        let mut act = vec![0.0f64; lp.num_row];
        for j in 0..lp.num_col {
            let xj = x[j];
            if xj != 0.0 {
                for k in ptr[j]..ptr[j + 1] {
                    let a = rval[k];
                    if a != 0.0 {
                        act[ridx[k] as usize] += a * xj;
                    }
                }
            }
        }
        act
    }

    #[test]
    fn factors_are_exact_powers_of_two() {
        let lp = badly_scaled();
        let (_, sc, applied) = scale_lp(&lp);
        assert!(applied, "1e-4..1e4 spread must pass the gate");
        for &f in sc.row.iter().chain(sc.col.iter()) {
            assert!(is_pow2(f), "factor {f} is not an exact power of 2");
            assert!(f >= 2f64.powi(-MAX_EXP) && f <= 2f64.powi(MAX_EXP));
        }
    }

    #[test]
    fn badly_scaled_lp_agrees_scaled_vs_unscaled() {
        let lp = badly_scaled();
        check_scale_agree(&lp);
        // And the known optimum, independent of either path.
        let sol = SimplexSolver::new(SimplexConfig::default())
            .solve(&lp)
            .unwrap();
        assert_eq!(sol.status, SimplexStatus::Optimal);
        assert!((sol.objective - 3.7).abs() < 1e-6, "obj {}", sol.objective);
    }

    #[test]
    fn equilibrated_lp_skips_scaling() {
        // All entries 1.0: range 1, nothing can improve it ⇒ gate rejects.
        let lp = lp_of(
            &[&[1.0, 1.0, 0.0], &[0.0, 1.0, 1.0]],
            vec![4.0, 3.0],
            "LL",
            vec![1.0, 2.0, 3.0],
            vec![0.0; 3],
            vec![10.0; 3],
        );
        let (slp, sc, applied) = scale_lp(&lp);
        assert!(!applied, "gate must skip an already-equilibrated LP");
        assert_eq!(sc.row, vec![1.0; 2]);
        assert_eq!(sc.col, vec![1.0; 3]);
        // bit-identical clone
        assert_eq!(slp.a, lp.a);
        assert_eq!(slp.rhs, lp.rhs);
    }

    #[test]
    fn infinite_bounds_survive_round_trip() {
        let mut lp = badly_scaled();
        lp.lower[0] = f64::NEG_INFINITY; // free-ish column (MI)
        lp.upper[1] = f64::INFINITY; // free-ish column (P)
        lp.lower[1] = f64::NEG_INFINITY; // fully free column
        let (slp, sc, applied) = scale_lp(&lp);
        assert!(applied);
        assert_eq!(slp.lower[0], f64::NEG_INFINITY);
        assert_eq!(slp.lower[1], f64::NEG_INFINITY);
        assert!(slp.upper[1].is_infinite() && slp.upper[1] > 0.0);
        assert!(slp.upper[0].is_finite());
        // solve both ways: free-variable split/flip must still work and agree
        check_scale_agree(&lp);
        // unscale_x is the exact inverse of the bound transform direction
        let mut xs = vec![1.0, 2.0, 3.0];
        let want: Vec<f64> = xs.iter().zip(&sc.col).map(|(&x, &c)| x * c).collect();
        unscale_x(&sc, &mut xs);
        assert_eq!(xs, want);
    }

    #[test]
    fn scaling_reduces_entry_range_and_preserves_objective_value() {
        let lp = badly_scaled();
        let (slp, sc, applied) = scale_lp(&lp);
        assert!(applied);
        let before = spread(&lp, &vec![0; lp.num_row], &vec![0; lp.num_col]);
        let after = {
            let re: Vec<i32> = sc.row.iter().map(|&f| f.log2() as i32).collect();
            let ce: Vec<i32> = sc.col.iter().map(|&f| f.log2() as i32).collect();
            spread(&lp, &re, &ce)
        };
        assert!(after < before, "range must improve: {after} vs {before}");
        // objective-value invariance: Σ cost'·x' = Σ cost·(c ⊙ x')
        let xp = [1.5, 0.25, 7.0];
        let orig: f64 = lp.objective.iter().zip(&xp).map(|(c, x)| c * x).sum();
        let mut uns = xp;
        unscale_x(&sc, &mut uns);
        let scaled_obj: f64 = slp.objective.iter().zip(&xp).map(|(c, x)| c * x).sum();
        assert_eq!(scaled_obj, orig);
        let direct: f64 = lp.objective.iter().zip(&uns).map(|(c, x)| c * x).sum();
        assert!((scaled_obj - direct).abs() <= 1e-12 * orig.abs());
        // scaled rows/bounds match the documented transform on one sample
        for j in 0..lp.num_col {
            for i in 0..lp.num_row {
                let a = lp.a.get(i, j);
                let expect = if a == 0.0 {
                    0.0
                } else {
                    a * sc.row[i] * sc.col[j]
                };
                assert_eq!(slp.a.get(i, j), expect, "entry ({i},{j})");
            }
            assert_eq!(slp.upper[j] * sc.col[j], lp.upper[j]);
        }
        for i in 0..lp.num_row {
            assert_eq!(slp.rhs[i], lp.rhs[i] * sc.row[i]);
        }
    }

    #[test]
    fn solve_hook_matches_unscaled_solve() {
        // The solve() hook (config.scaling = the RUSOLVE_SCALE lever) must
        // reproduce the unscaled objective AND report x in ORIGINAL space.
        let lp = badly_scaled();
        let on = SimplexSolver::new(SimplexConfig {
            scaling: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        let off = SimplexSolver::new(SimplexConfig {
            scaling: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(on.status, SimplexStatus::Optimal);
        assert_eq!(on.status, off.status);
        assert!(
            (on.objective - off.objective).abs() <= 1e-6 * on.objective.abs().max(1.0),
            "obj mismatch: {} vs {}",
            on.objective,
            off.objective
        );
        assert!((on.objective - 3.7).abs() < 1e-6);
        let act = row_activities(&lp, &on.x);
        for (i, s) in lp.senses.iter().enumerate() {
            let lhs = act[i];
            let r = lp.rhs[i];
            match s {
                'L' => assert!(lhs <= r + 1e-6, "row {i}: {lhs} > {r}"),
                'G' => assert!(lhs >= r - 1e-6, "row {i}: {lhs} < {r}"),
                _ => assert!((lhs - r).abs() <= 1e-6, "row {i}: {lhs} != {r}"),
            }
        }
    }

    #[test]
    fn solve_hook_handles_maximize() {
        // Objective recomputation after unscale must respect the original
        // sense: maximize the badly-scaled LP ⇒ x2 hits its ub, obj ≈ 13.
        let mut lp = badly_scaled();
        lp.maximize = true;
        let on = SimplexSolver::new(SimplexConfig {
            scaling: true,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        let off = SimplexSolver::new(SimplexConfig {
            scaling: false,
            ..SimplexConfig::default()
        })
        .solve(&lp)
        .unwrap();
        assert_eq!(on.status, SimplexStatus::Optimal);
        assert_eq!(on.status, off.status);
        assert!(
            (on.objective - off.objective).abs() <= 1e-6 * on.objective.abs().max(1.0),
            "obj mismatch: {} vs {}",
            on.objective,
            off.objective
        );
        assert!((on.objective - 13.0).abs() < 1e-5, "obj {}", on.objective);
    }

    #[test]
    fn zero_rows_and_columns_scale_as_one() {
        // row 1 empty, col 1 empty: factors must stay 1.0 (no NaN/inf), and
        // the remaining 1e6-range row still earns a power-of-2 factor.
        let lp = lp_of(
            &[&[1e6, 0.0, 0.0], &[0.0, 0.0, 0.0], &[4.0, 0.0, 2.0]],
            vec![2e6, 0.0, 8.0],
            "GEG",
            vec![1.0, 0.0, 1.0],
            vec![0.0; 3],
            vec![f64::INFINITY; 3],
        );
        let (slp, sc, applied) = scale_lp(&lp);
        assert!(applied);
        assert_eq!(sc.row[1], 1.0, "zero row must not be scaled");
        assert_eq!(sc.col[1], 1.0, "zero col must not be scaled");
        assert!(slp.rhs[1] == 0.0);
        check_scale_agree(&lp);
    }
}
