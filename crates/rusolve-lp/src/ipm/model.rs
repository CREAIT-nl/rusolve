//! IPM internal LP: formulation (fixed-column substitution, slack columns,
//! artificial bounds on free columns), Ruiz scaling, postprocessing back to
//! the caller's [`SparseLp`], and the original-LP verification gate.
//!
//! Provenance:
//! * internal form `min cᵀx, Ax = b, l ≤ x ≤ u` with one slack per
//!   inequality row: HiGHS HiPO `Model::init`/`PreProcess` (MIT) and
//!   IPX `Model::Load` (MIT), reimplemented — not translated;
//! * Ruiz ∞-norm equilibration: D. Ruiz, "A scaling algorithm to
//!   equilibrate both rows and columns norms in matrices", RAL-TR-2001-034;
//!   rounding to powers of two (exact scaling): HiPO `Model::scale` (MIT);
//! * artificial bounds ±1e4 on free columns, widened as the iterate
//!   approaches them: HiPO `adjustFreeVars` (MIT);
//! * verification on the original LP (primal residual, sign-feasible
//!   Lagrangian dual, duality gap): standard LP duality, e.g. Vanderbei,
//!   *Linear Programming* (5th ed.), ch. 5 and 9.

use rusolve_core::{CscBuilder, CscMatrix, SparseLp};

/// Initial half-width of the artificial box on a free column (scaled space).
pub const FREE_BOUND: f64 = 1e4;
/// Ruiz passes and the clamp on accumulated scale factors.
const RUIZ_PASSES: usize = 10;
const SCALE_MIN: f64 = 1e-4;
const SCALE_MAX: f64 = 1e4;

/// Why the LP never reached the IPM.
#[derive(Debug, Clone, PartialEq)]
pub enum FormError {
    /// Contradictory bounds found by the formulation (lb > ub, empty row
    /// whose bounds exclude 0).
    Infeasible(String),
    /// NaN data.
    Unsupported(String),
}

/// Scaled internal LP the IPM iterates on.
#[derive(Debug, Clone)]
pub struct IpmLp {
    pub m: usize,
    /// Internal columns: structurals first, then one slack per inequality row.
    pub n: usize,
    pub n_struct: usize,
    pub a: CscMatrix,
    pub b: Vec<f64>,
    pub c: Vec<f64>,
    pub l: Vec<f64>,
    pub u: Vec<f64>,
    /// Bounds of column j are artificial (free column in the caller's LP).
    pub art: Vec<bool>,
    /// Row / column scale factors: `Â = R·A·C` (slack column scale = 1/R_i).
    pub rscale: Vec<f64>,
    pub cscale: Vec<f64>,
    /// Internal row → caller row; internal structural → caller column.
    pub row_orig: Vec<usize>,
    pub col_orig: Vec<usize>,
    /// Caller column → fixed value (substituted out).
    pub fixed: Vec<Option<f64>>,
    /// +1 minimize, −1 maximize (internal objective = sense · c).
    pub sense: f64,
    pub n_free: usize,
    pub n_dropped_rows: usize,
}

impl IpmLp {
    /// Number of finite (lower + upper) bounds = complementarity pairs.
    pub fn num_bounds(&self) -> usize {
        self.l.iter().filter(|v| v.is_finite()).count()
            + self.u.iter().filter(|v| v.is_finite()).count()
    }
}

/// Build the internal LP from the caller's `SparseLp`.
pub fn formulate(sp: &SparseLp) -> Result<IpmLp, FormError> {
    let a0 = &sp.a;
    let (m0, n0) = (a0.nrow(), a0.ncol());
    let sense = if sp.maximize { -1.0 } else { 1.0 };
    let (ptr, ri, va) = (a0.col_ptr(), a0.row_idx(), a0.values());
    for j in 0..n0 {
        let (lb, ub) = (sp.col_lower[j], sp.col_upper[j]);
        if lb.is_nan() || ub.is_nan() || sp.objective[j].is_nan() {
            return Err(FormError::Unsupported(format!("NaN data on column {j}")));
        }
        if lb > ub || lb == f64::INFINITY || ub == f64::NEG_INFINITY {
            return Err(FormError::Infeasible(format!(
                "column {j}: bounds [{lb}, {ub}]"
            )));
        }
    }
    for i in 0..m0 {
        let (rl, ru) = (sp.row_lower[i], sp.row_upper[i]);
        if rl.is_nan() || ru.is_nan() {
            return Err(FormError::Unsupported(format!("NaN bound on row {i}")));
        }
        if rl > ru {
            return Err(FormError::Infeasible(format!("row {i}: [{rl}, {ru}]")));
        }
    }
    // Fixed columns: substitute out.
    let mut fixed = vec![None; n0];
    let mut off = vec![0.0f64; m0];
    let mut row_cnt = vec![0usize; m0];
    for j in 0..n0 {
        if sp.col_lower[j] == sp.col_upper[j] {
            let v = sp.col_lower[j];
            fixed[j] = Some(v);
            for k in ptr[j]..ptr[j + 1] {
                off[ri[k] as usize] += va[k] * v;
            }
        } else {
            for k in ptr[j]..ptr[j + 1] {
                if va[k] != 0.0 {
                    row_cnt[ri[k] as usize] += 1;
                }
            }
        }
    }
    // Rows: drop free and empty rows (an empty row must contain 0).
    let mut row_new = vec![usize::MAX; m0];
    let mut row_orig = Vec::new();
    let mut n_dropped_rows = 0;
    for i in 0..m0 {
        let (rl, ru) = (sp.row_lower[i] - off[i], sp.row_upper[i] - off[i]);
        let free = rl == f64::NEG_INFINITY && ru == f64::INFINITY;
        if free {
            n_dropped_rows += 1;
            continue;
        }
        if row_cnt[i] == 0 {
            let tol = 1e-9;
            if rl > tol * (1.0 + sp.row_lower[i].abs()) || ru < -tol * (1.0 + sp.row_upper[i].abs())
            {
                return Err(FormError::Infeasible(format!(
                    "empty row {i} needs activity in [{rl}, {ru}]"
                )));
            }
            n_dropped_rows += 1;
            continue;
        }
        row_new[i] = row_orig.len();
        row_orig.push(i);
    }
    let m = row_orig.len();
    let col_orig: Vec<usize> = (0..n0).filter(|&j| fixed[j].is_none()).collect();
    let ns = col_orig.len();

    // Structural block with internal row indices (CSC, rows ascending).
    let mut sptr = Vec::with_capacity(ns + 1);
    sptr.push(0usize);
    let mut srow: Vec<usize> = Vec::new();
    let mut sval: Vec<f64> = Vec::new();
    for &j in &col_orig {
        for k in ptr[j]..ptr[j + 1] {
            let r = row_new[ri[k] as usize];
            if r != usize::MAX && va[k] != 0.0 {
                srow.push(r);
                sval.push(va[k]);
            }
        }
        sptr.push(srow.len());
    }

    // Ruiz equilibration (geometric ∞-norm passes).
    let mut rs = vec![1.0f64; m];
    let mut cs = vec![1.0f64; ns];
    let mut rmax = vec![0.0f64; m];
    for _ in 0..RUIZ_PASSES {
        rmax.iter_mut().for_each(|v| *v = 0.0);
        let mut worst = 0.0f64;
        let mut cmaxs = vec![0.0f64; ns];
        for k in 0..ns {
            let mut cmax = 0.0f64;
            for p in sptr[k]..sptr[k + 1] {
                let v = (rs[srow[p]] * sval[p] * cs[k]).abs();
                cmax = cmax.max(v);
                rmax[srow[p]] = rmax[srow[p]].max(v);
            }
            cmaxs[k] = cmax;
        }
        for k in 0..ns {
            if cmaxs[k] > 0.0 {
                worst = worst.max((1.0 - cmaxs[k]).abs());
                cs[k] = (cs[k] / cmaxs[k].sqrt()).clamp(SCALE_MIN, SCALE_MAX);
            }
        }
        for i in 0..m {
            if rmax[i] > 0.0 {
                worst = worst.max((1.0 - rmax[i]).abs());
                rs[i] = (rs[i] / rmax[i].sqrt()).clamp(SCALE_MIN, SCALE_MAX);
            }
        }
        if worst < 1e-2 {
            break;
        }
    }
    let pow2 = |v: f64| v.log2().round().exp2();
    rs.iter_mut().for_each(|v| *v = pow2(*v));
    cs.iter_mut().for_each(|v| *v = pow2(*v));

    // Internal columns.
    let mut l = Vec::with_capacity(ns + m);
    let mut u = Vec::with_capacity(ns + m);
    let mut c = Vec::with_capacity(ns + m);
    let mut art = Vec::with_capacity(ns + m);
    let mut cscale = Vec::with_capacity(ns + m);
    let mut n_free = 0;
    let mut bld = CscBuilder::new(m);
    let mut col: Vec<(usize, f64)> = Vec::new();
    for (k, &j) in col_orig.iter().enumerate() {
        col.clear();
        for p in sptr[k]..sptr[k + 1] {
            col.push((srow[p], rs[srow[p]] * sval[p] * cs[k]));
        }
        bld.push_col(&col);
        let (lb, ub) = (sp.col_lower[j], sp.col_upper[j]);
        c.push(sense * sp.objective[j] * cs[k]);
        cscale.push(cs[k]);
        if lb == f64::NEG_INFINITY && ub == f64::INFINITY {
            l.push(-FREE_BOUND);
            u.push(FREE_BOUND);
            art.push(true);
            n_free += 1;
        } else {
            l.push(lb / cs[k]);
            u.push(ub / cs[k]);
            art.push(false);
        }
    }
    let mut b = vec![0.0; m];
    for (i, &io) in row_orig.iter().enumerate() {
        let rl = (sp.row_lower[io] - off[io]) * rs[i];
        let ru = (sp.row_upper[io] - off[io]) * rs[i];
        if rl == ru {
            b[i] = rl;
            continue;
        }
        if ru.is_finite() {
            b[i] = ru;
            bld.push_col(&[(i, 1.0)]);
            u.push(ru - rl);
        } else {
            b[i] = rl;
            bld.push_col(&[(i, -1.0)]);
            u.push(f64::INFINITY);
        }
        l.push(0.0);
        c.push(0.0);
        art.push(false);
        cscale.push(1.0 / rs[i]);
    }
    let a = bld.into_matrix();
    let n = a.ncol();
    Ok(IpmLp {
        m,
        n,
        n_struct: ns,
        a,
        b,
        c,
        l,
        u,
        art,
        rscale: rs,
        cscale,
        row_orig,
        col_orig,
        fixed,
        sense,
        n_free,
        n_dropped_rows,
    })
}

/// Caller-space primal x, row duals π and reduced costs d = c − Aᵀπ.
pub fn postprocess(
    sp: &SparseLp,
    lp: &IpmLp,
    xh: &[f64],
    yh: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n0 = sp.a.ncol();
    let m0 = sp.a.nrow();
    let mut x = vec![0.0; n0];
    for (j, f) in lp.fixed.iter().enumerate() {
        if let Some(v) = f {
            x[j] = *v;
        }
    }
    for (k, &j) in lp.col_orig.iter().enumerate() {
        let v = lp.cscale[k] * xh[k];
        // Snap onto a bound the iterate sits within rounding of.
        x[j] = v.clamp(sp.col_lower[j], sp.col_upper[j]);
    }
    let mut pi = vec![0.0; m0];
    for (i, &io) in lp.row_orig.iter().enumerate() {
        pi[io] = lp.sense * lp.rscale[i] * yh[i];
    }
    let d = reduced_costs(sp, &pi);
    (x, pi, d)
}

pub fn reduced_costs(sp: &SparseLp, pi: &[f64]) -> Vec<f64> {
    let (ptr, ri, va) = (sp.a.col_ptr(), sp.a.row_idx(), sp.a.values());
    (0..sp.a.ncol())
        .map(|j| {
            let mut s = sp.objective[j];
            for k in ptr[j]..ptr[j + 1] {
                s -= va[k] * pi[ri[k] as usize];
            }
            s
        })
        .collect()
}

/// Original-LP optimality measures (all relative, ∞-norm).
#[derive(Debug, Clone, Copy, Default)]
pub struct Verify {
    /// max bound / row violation / (1 + |bound|).
    pub pinf: f64,
    /// max sign violation of (π, d) / (1 + ‖c‖∞); NaN when no duals checked.
    pub dinf: f64,
    /// |pobj − dobj| / (1 + |pobj|); NaN when no duals checked.
    pub gap: f64,
    /// cᵀx in the caller's sense.
    pub pobj: f64,
    pub dobj: f64,
}

impl Verify {
    pub fn ok(&self, tol: f64) -> bool {
        self.pinf <= tol
            && (self.dinf.is_nan() || self.dinf <= tol)
            && (self.gap.is_nan() || self.gap <= tol)
    }
}

/// Check (x, π) on the caller's LP. With `pi = None` only primal
/// feasibility is checked (presolve-ON path: duals are not postsolved).
pub fn verify(sp: &SparseLp, x: &[f64], pi: Option<&[f64]>) -> Verify {
    let (m0, n0) = (sp.a.nrow(), sp.a.ncol());
    let mut act = vec![0.0; m0];
    sp.a.matvec(x, &mut act);
    let mut pinf = 0.0f64;
    for i in 0..m0 {
        let (rl, ru) = (sp.row_lower[i], sp.row_upper[i]);
        if rl.is_finite() {
            pinf = pinf.max((rl - act[i]) / (1.0 + rl.abs()));
        }
        if ru.is_finite() {
            pinf = pinf.max((act[i] - ru) / (1.0 + ru.abs()));
        }
    }
    for j in 0..n0 {
        let (lb, ub) = (sp.col_lower[j], sp.col_upper[j]);
        if lb.is_finite() {
            pinf = pinf.max((lb - x[j]) / (1.0 + lb.abs()));
        }
        if ub.is_finite() {
            pinf = pinf.max((x[j] - ub) / (1.0 + ub.abs()));
        }
        if !x[j].is_finite() {
            pinf = f64::INFINITY;
        }
    }
    let pobj: f64 = sp.objective.iter().zip(x).map(|(c, v)| c * v).sum();
    let mut out = Verify {
        pinf: pinf.max(0.0),
        dinf: f64::NAN,
        gap: f64::NAN,
        pobj,
        dobj: f64::NAN,
    };
    let Some(pi) = pi else {
        return out;
    };
    // Minimization form: c' = s·c, π' = s·π, d' = s·d.
    let s = if sp.maximize { -1.0 } else { 1.0 };
    let d = reduced_costs(sp, pi);
    let cnorm = sp.objective.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    let mut viol = 0.0f64;
    let mut dobj = 0.0f64;
    // Lagrangian term: for a multiplier v on [lo, hi], v·lo if v > 0, v·hi
    // if v < 0; a wrong-signed v against an infinite side is a violation.
    let mut term = |v: f64, lo: f64, hi: f64| {
        if v > 0.0 {
            if lo.is_finite() {
                dobj += v * lo;
            } else {
                viol = viol.max(v);
            }
        } else if v < 0.0 {
            if hi.is_finite() {
                dobj += v * hi;
            } else {
                viol = viol.max(-v);
            }
        }
    };
    for i in 0..m0 {
        term(s * pi[i], sp.row_lower[i], sp.row_upper[i]);
    }
    for j in 0..n0 {
        term(s * d[j], sp.col_lower[j], sp.col_upper[j]);
    }
    let pobj_min = s * pobj;
    out.dinf = viol / (1.0 + cnorm);
    out.dobj = s * dobj;
    out.gap = (pobj_min - dobj).abs() / (1.0 + pobj_min.abs());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lp1() -> SparseLp {
        // min -x0 - 2x1 + 0x2 (fixed) + x3 (free)
        // r0: x0 + x1 + x2 <= 4 ; r1: 1 <= x0 - x1 + x3 <= 3 ; r2 free ; r3: x3 = 2
        let a = CscMatrix::from_triplets(
            4,
            4,
            &[
                (0, 0, 1.0),
                (1, 0, 1.0),
                (2, 0, 5.0),
                (0, 1, 1.0),
                (1, 1, -1.0),
                (0, 2, 1.0),
                (1, 3, 1.0),
                (3, 3, 1.0),
            ],
        );
        SparseLp::new(
            a,
            vec![f64::NEG_INFINITY, 1.0, f64::NEG_INFINITY, 2.0],
            vec![4.0, 3.0, f64::INFINITY, 2.0],
            vec![0.0, 0.0, 1.0, f64::NEG_INFINITY],
            vec![10.0, f64::INFINITY, 1.0, f64::INFINITY],
            vec![-1.0, -2.0, 0.0, 1.0],
            false,
        )
    }

    #[test]
    fn formulation_shapes_and_unscale_roundtrip() {
        let sp = lp1();
        let lp = formulate(&sp).unwrap();
        assert_eq!(lp.m, 3); // free row dropped
        assert_eq!(lp.n_struct, 3); // x2 fixed
        assert_eq!(lp.n, 3 + 2); // two inequality slacks
        assert_eq!(lp.n_free, 1);
        assert_eq!(lp.fixed[2], Some(1.0));
        // A known feasible point maps to Ax = b in the internal space.
        let x = [1.5, 0.5, 1.0, 2.0];
        let v = verify(&sp, &x, None);
        assert!(v.pinf <= 0.0, "{v:?}");
        let mut xh = vec![0.0; lp.n];
        for (k, &j) in lp.col_orig.iter().enumerate() {
            xh[k] = x[j] / lp.cscale[k];
        }
        // slacks: r0 (<=4): s = 4 - 2 - 1 = 1 ; r1 ([1,3]): s = 3 - (1 + 2) = 0
        let act0 = 1.5 + 0.5 + 1.0;
        let act1 = 1.5 - 0.5 + 2.0;
        xh[3] = (4.0 - act0) * lp.rscale[0];
        xh[4] = (3.0 - act1) * lp.rscale[1];
        let mut r = vec![0.0; lp.m];
        lp.a.matvec(&xh, &mut r);
        for i in 0..lp.m {
            assert!(
                (r[i] - lp.b[i]).abs() < 1e-12,
                "row {i}: {} vs {}",
                r[i],
                lp.b[i]
            );
        }
    }

    #[test]
    fn infeasible_bounds_and_empty_rows_detected() {
        let a = CscMatrix::from_triplets(1, 1, &[(0, 0, 1.0)]);
        let sp = SparseLp::new(
            a.clone(),
            vec![3.0],
            vec![3.0],
            vec![1.0],
            vec![1.0],
            vec![1.0],
            false,
        );
        assert!(matches!(formulate(&sp), Err(FormError::Infeasible(_))));
        // empty row 1 needs activity in [1, 2]
        let a2 = CscMatrix::from_triplets(2, 1, &[(0, 0, 1.0)]);
        let sp = SparseLp::new(
            a2,
            vec![0.0, 1.0],
            vec![1.0, 2.0],
            vec![0.0],
            vec![1.0],
            vec![1.0],
            false,
        );
        assert!(matches!(formulate(&sp), Err(FormError::Infeasible(_))));
    }

    #[test]
    fn verify_detects_wrong_dual_signs() {
        // min x s.t. x >= 1 (row), x free: optimum x = 1, π = 1.
        let a = CscMatrix::from_triplets(1, 1, &[(0, 0, 1.0)]);
        let sp = SparseLp::new(
            a,
            vec![1.0],
            vec![f64::INFINITY],
            vec![f64::NEG_INFINITY],
            vec![f64::INFINITY],
            vec![1.0],
            false,
        );
        let good = verify(&sp, &[1.0], Some(&[1.0]));
        assert!(good.ok(1e-9), "{good:?}");
        let bad = verify(&sp, &[1.0], Some(&[-1.0]));
        assert!(!bad.ok(1e-6), "{bad:?}");
    }
}
