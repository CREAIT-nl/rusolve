//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

use super::*;

/// SCIP def.h:219 `MAX(x, y) ((x) >= (y) ? (x) : (y))`.
#[inline]
pub(crate) fn s4_max(x: f64, y: f64) -> f64 {
    if x >= y { x } else { y }
}

/// SCIP def.h:223 `MIN(x, y) ((x) <= (y) ? (x) : (y))`.
#[inline]
pub(crate) fn s4_min(x: f64, y: f64) -> f64 {
    if x <= y { x } else { y }
}

/// SCIP def.h:215 `ABS(x) ((x) >= 0 ? (x) : -(x))`.
#[inline]
fn s4_abs(x: f64) -> f64 {
    if x >= 0.0 { x } else { -x }
}

/// SCIP 11.0.0 (32dd425): lp.c::rowCalcNorms (lines 4974-4975, 4990-5004)
/// Computes the squared Euclidean norm of a row's coefficient vector.
pub(crate) fn s4_row_sqrnorm(val: &[f64]) -> f64 {
    let mut sqrnorm = 0.0;
    for &v in val {
        // rowAddNorms lp.c:2062: sqrnorm += SQR(absval)
        let absval = v.abs();
        sqrnorm += absval * absval;
    }
    sqrnorm
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetNorm (lines 6395-6399 via rowCalcNorms)
/// Computes the Euclidean norm of a row's coefficient vector.
pub(crate) fn s4_row_norm(val: &[f64]) -> f64 {
    s4_row_sqrnorm(val).sqrt()
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwRecalcLPActivity (lines 6402-6451)
/// Computes the activity of a row given a solution vector, clamped to [-S4_INFINITY, S4_INFINITY].
pub(crate) fn s4_row_activity(idx: &[HighsInt], val: &[f64], sol: &[f64]) -> f64 {
    // lp.c:6413: activity = row->constant (pool rows have no constant).
    // lp.c:6424-6450 (unlinked columns, #ifndef NDEBUG block): not reachable, every pool
    // column is an LP column.
    let mut activity = 0.0;
    for i in 0..idx.len() {
        let col_idx = idx[i] as usize;
        activity += val[i] * sol[col_idx];
    }
    activity = s4_max(activity, -S4_INFINITY);
    activity = s4_min(activity, S4_INFINITY);
    activity
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetLPFeasibility (lines 6484-6498)
/// Computes the feasibility of a row: min(rhs - activity, activity - lhs).
pub(crate) fn s4_row_feasibility(lhs: f64, rhs: f64, act: f64) -> f64 {
    s4_min(rhs - act, act - lhs)
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetLPEfficacy (lines 7054-7092)
/// Computes the efficacy of a row: -feasibility / max(norm, sumepsilon).
pub(crate) fn s4_row_efficacy(feas: f64, norm: f64) -> f64 {
    let norm = s4_max(norm, S4_SUMEPSILON);
    -feas / norm
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetNumIntCols (lines 6965-6977), counted as in
/// rowCalcIdxsAndVals (lp.c:5037-5061). SCIPcolIsIntegral (lp.c:17455-17463) is
/// `variable_type != kContinuous`, so implicit integers count.
pub(crate) fn s4_row_num_int_cols(idx: &[HighsInt], var_type: &[HighsVarType]) -> HighsInt {
    let mut numintcols: HighsInt = 0;
    for i in 0..idx.len() {
        if var_type[idx[i] as usize] != HighsVarType::kContinuous {
            numintcols += 1;
        }
    }
    numintcols
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIPlpRecalculateObjSqrNorm (lines 18080-18107)
/// Computes the squared Euclidean norm of the objective vector, iterating from last to first.
pub(crate) fn s4_obj_sqrnorm(cost: &[f64]) -> f64 {
    let mut objsqrnorm = 0.0;
    for i in (0..cost.len()).rev() {
        objsqrnorm += cost[i] * cost[i];
    }
    s4_max(objsqrnorm, 0.0)
}

/// SCIP 11.0.0 (32dd425): lp.c::rowCalcNorms (lines 4976, 4990-5004)
/// Computes the scalar product of the row coefficients with the objective vector.
pub(crate) fn s4_row_objprod(idx: &[HighsInt], val: &[f64], cost: &[f64]) -> f64 {
    let mut objprod = 0.0;
    for i in 0..idx.len() {
        let col_idx = idx[i] as usize;
        objprod += val[i] * cost[col_idx];
    }
    objprod
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetObjParallelism (lines 8046-8076)
/// Computes the parallelism of a row with the objective function.
pub(crate) fn s4_row_obj_parallelism(
    idx: &[HighsInt],
    val: &[f64],
    cost: &[f64],
    objsqrnorm: f64,
) -> f64 {
    let sqrnorm = s4_row_sqrnorm(val);
    let objprod = s4_row_objprod(idx, val, cost);
    let prod = sqrnorm * objsqrnorm;
    // SCIPsetIsPositive(set, prod): prod > epsilon
    let mut parallelism = if prod > S4_EPSILON {
        objprod.abs() / prod.sqrt()
    } else {
        0.0
    };
    parallelism = s4_min(parallelism, 1.0);
    parallelism = s4_max(parallelism, 0.0);
    parallelism
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetScalarProduct (lines 7254-7383, easy case)
/// Computes the scalar product of two rows, assuming sorted indices and no unlinked columns.
pub(crate) fn s4_row_scalar_product(
    idx1: &[HighsInt],
    val1: &[f64],
    idx2: &[HighsInt],
    val2: &[f64],
) -> f64 {
    // lp.c:7314-7348 (#ifndef NDEBUG partition check) dropped. Easy case 1) applies:
    // every pool column is an LP column, so MAX(nlpcols, nunlinked) - 1 = len - 1.
    let mut i1 = idx1.len() as i64 - 1;
    let mut i2 = idx2.len() as i64 - 1;
    let mut scalarprod = 0.0;

    while i1 >= 0 && i2 >= 0 {
        let id1 = idx1[i1 as usize];
        let id2 = idx2[i2 as usize];
        if id1 < id2 {
            i2 -= 1;
        } else if id1 > id2 {
            i1 -= 1;
        } else {
            scalarprod += val1[i1 as usize] * val2[i2 as usize];
            i1 -= 1;
            i2 -= 1;
        }
    }
    scalarprod
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetParallelism (lines 7970-8028, case 'e')
/// Computes the parallelism between two rows using Euclidean norm.
pub(crate) fn s4_row_parallelism(
    idx1: &[HighsInt],
    val1: &[f64],
    idx2: &[HighsInt],
    val2: &[f64],
) -> f64 {
    let scalarprod = s4_row_scalar_product(idx1, val1, idx2, val2);
    if scalarprod == 0.0 {
        return 0.0;
    }
    // lp.c:7989-8012: the norm == 0 re-evaluation adds SQR(vals) to a zero sqrnorm in
    // row order, which reproduces the freshly computed sqrnorm here, so it is a no-op.
    let norm1 = s4_row_norm(val1);
    let norm2 = s4_row_norm(val2);
    scalarprod.abs() / (norm1 * norm2)
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetLPSolCutoffDistance (lines 6997-7051, direction calc)
/// Computes the normalized direction vector from lp_sol to incumbent.
pub(crate) fn s4_sol_direction(incumbent: &[f64], lp_sol: &[f64]) -> Vec<f64> {
    let n = incumbent.len();
    let mut direction = vec![0.0; n];
    let mut scale = 0.0;

    for k in 0..n {
        direction[k] = incumbent[k] - lp_sol[k];
        scale += direction[k] * direction[k];
    }

    if scale > 0.0 {
        scale = 1.0 / scale.sqrt();
        for k in 0..n {
            direction[k] *= scale;
        }
    }
    direction
}

/// SCIP 11.0.0 (32dd425): lp.c::SCIProwGetLPSolCutoffDistance (lines 7035-7050)
/// Computes the cutoff distance of a row in the given direction.
pub(crate) fn s4_row_cutoff_distance(
    idx: &[HighsInt],
    val: &[f64],
    rhs: f64,
    dir: &[f64],
    act: f64,
) -> f64 {
    let mut solcutoffdist = 0.0;
    for i in 0..idx.len() {
        let col_idx = idx[i] as usize;
        solcutoffdist += val[i] * dir[col_idx];
    }

    // SCIPsetIsSumZero: |x| <= sumepsilon
    if solcutoffdist.abs() <= S4_SUMEPSILON {
        solcutoffdist = S4_SUMEPSILON;
    }

    let feas = s4_row_feasibility(-S4_INFINITY, rhs, act);
    -feas / s4_abs(solcutoffdist)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s4_metric_efficacy_small_norm() {
        let feas = -1.0;
        let norm = 1e-10; // < S4_SUMEPSILON (1e-6)
        let eff = s4_row_efficacy(feas, norm);
        // norm becomes S4_SUMEPSILON = 1e-6
        // eff = -(-1.0) / 1e-6 = 1e6
        assert_eq!(eff, 1e6);
    }

    #[test]
    fn s4_metric_feasibility_infinite_lhs() {
        let lhs = -1e20;
        let rhs = 10.0;
        let act = 5.0;
        let feas = s4_row_feasibility(lhs, rhs, act);
        // min(10 - 5, 5 - (-1e20)) = min(5, 1e20+5) = 5
        assert_eq!(feas, 5.0);
    }

    #[test]
    fn s4_metric_obj_parallelism_zero_cost() {
        let idx = vec![0, 1];
        let val = vec![1.0, 1.0];
        let cost = vec![0.0, 0.0];
        let objsqrnorm = s4_obj_sqrnorm(&cost); // 0.0
        let par = s4_row_obj_parallelism(&idx, &val, &cost, objsqrnorm);
        assert_eq!(par, 0.0);
    }

    #[test]
    fn s4_metric_obj_parallelism_clamp() {
        // Construct a case where prod > eps but ratio > 1 due to numerical noise or specific values?
        // Actually, Cauchy-Schwarz ensures |objprod| <= sqrt(sqrnorm * objsqrnorm).
        // So ratio <= 1. The clamp is for safety.
        // Let's test a simple case: row = [1, 0], cost = [1, 0].
        // sqrnorm = 1, objsqrnorm = 1, objprod = 1. prod = 1. ratio = 1.
        let idx = vec![0];
        let val = vec![1.0];
        let cost = vec![1.0];
        let objsqrnorm = s4_obj_sqrnorm(&cost);
        let par = s4_row_obj_parallelism(&idx, &val, &cost, objsqrnorm);
        assert_eq!(par, 1.0);
    }

    #[test]
    fn s4_metric_efficacy_hand() {
        // row 3x0 + 4x1 <= 1 at x = (1, 1): act = 7, feas = min(1-7, 7+1e20) = -6,
        // norm = 5, efficacy = 6/5.
        let idx = vec![0, 1];
        let val = vec![3.0, 4.0];
        let act = s4_row_activity(&idx, &val, &[1.0, 1.0]);
        assert_eq!(act, 7.0);
        let feas = s4_row_feasibility(-S4_INFINITY, 1.0, act);
        assert_eq!(feas, -6.0);
        assert_eq!(s4_row_norm(&val), 5.0);
        assert_eq!(s4_row_efficacy(feas, s4_row_norm(&val)), 6.0 / 5.0);
    }

    #[test]
    fn s4_metric_obj_parallelism_rounding_clamp() {
        // |x*c| / sqrt((x*x)*(c*c)) evaluates to 1.0000000000000002 here; MIN(.,1.0)
        // (lp.c:8072) clamps it to 1.
        let xv = 5.004807362210215_f64;
        let cv = 4.549961541408507_f64;
        assert!((xv * cv).abs() / ((xv * xv) * (cv * cv)).sqrt() > 1.0);
        let cost = vec![cv];
        let par = s4_row_obj_parallelism(&[0], &[xv], &cost, s4_obj_sqrnorm(&cost));
        assert_eq!(par, 1.0);
    }

    #[test]
    fn s4_metric_obj_sqrnorm_reverse_order() {
        // Forward: 1 + 1 + 1e16 = 1e16 + 2. Last column first (lp.c:18095-18098):
        // 1e16 + 1 rounds to 1e16 (ties to even), then + 1 again gives 1e16.
        let cost = vec![1.0, 1.0, 1e8];
        assert_eq!(s4_obj_sqrnorm(&cost), 1e16);
    }

    #[test]
    fn s4_metric_num_int_cols() {
        let vt = vec![
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kImplicitInteger,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ];
        // implicit integers count (SCIPcolIsIntegral)
        assert_eq!(s4_row_num_int_cols(&[0, 1, 2, 3], &vt), 2);
        assert_eq!(s4_row_num_int_cols(&[1, 3], &vt), 0);
        assert_eq!(s4_row_num_int_cols(&[0, 2, 4], &vt), 3);
    }

    #[test]
    fn s4_metric_scalar_product_interleaved() {
        // idx1: 0, 2, 4
        // idx2: 1, 2, 3
        // Common: 2
        // val1[1] * val2[1]
        let idx1 = vec![0, 2, 4];
        let val1 = vec![1.0, 2.0, 3.0];
        let idx2 = vec![1, 2, 3];
        let val2 = vec![4.0, 5.0, 6.0];
        let sp = s4_row_scalar_product(&idx1, &val1, &idx2, &val2);
        assert_eq!(sp, 10.0); // 2.0 * 5.0
    }

    #[test]
    fn s4_metric_scalar_product_descending_order() {
        // Common columns 1, 3, 6 with products 1, 1e100, -1e100. Ascending order gives
        // (1 + 1e100) - 1e100 = 0; the C loop runs from the end (lp.c:7367-7381):
        // (-1e100 + 1e100) + 1 = 1.
        let idx1 = vec![0, 1, 3, 5, 6];
        let val1 = vec![7.0, 1.0, 1e100, 9.0, -1e100];
        let idx2 = vec![1, 2, 3, 4, 6];
        let val2 = vec![1.0, 8.0, 1.0, 9.0, 1.0];
        assert_eq!(s4_row_scalar_product(&idx1, &val1, &idx2, &val2), 1.0);
        assert_eq!(s4_row_scalar_product(&idx2, &val2, &idx1, &val1), 1.0);
    }

    #[test]
    fn s4_metric_parallelism_zero_sp() {
        // Orthogonal rows
        let idx1 = vec![0];
        let val1 = vec![1.0];
        let idx2 = vec![1];
        let val2 = vec![1.0];
        let par = s4_row_parallelism(&idx1, &val1, &idx2, &val2);
        assert_eq!(par, 0.0);
    }

    #[test]
    fn s4_metric_sol_direction_zero_scale() {
        let inc = vec![1.0, 2.0];
        let lp = vec![1.0, 2.0];
        let dir = s4_sol_direction(&inc, &lp);
        assert_eq!(dir, vec![0.0, 0.0]);
    }

    #[test]
    fn s4_metric_cutoff_distance_small_dot() {
        // dot product is 0
        let idx = vec![0];
        let val = vec![1.0];
        let rhs = 10.0;
        let dir = vec![0.0];
        let act = 5.0;
        // dot = 0. |dot| <= S4_SUMEPSILON -> dot = S4_SUMEPSILON
        // feas = min(10-5, 5-(-1e20)) = 5
        // dist = -5 / 1e-6 = -5e6
        let dist = s4_row_cutoff_distance(&idx, &val, rhs, &dir, act);
        assert_eq!(dist, -5e6);
    }
}
