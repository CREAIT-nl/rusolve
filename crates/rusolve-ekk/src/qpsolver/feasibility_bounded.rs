//! HiGHS source: qpsolver/feasibility_bounded.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! NOTE: the unboundedness test compares `var_lo == +inf`
//! (sic) in HiGHS; keep it.
//!
//! HiGHS 1.16 (#3153): with `settings.h16_qp_bounded_start` the 1.16.0 body is
//! used (diagonal fast path, row-wise Cholesky, unboundedness test with the
//! sign fix, projection with `primal_feasibility_tolerance`). With the flag off
//! the 1.15.1 body below runs unchanged.

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: feasibility_bounded.hpp::computeStartingPointBounded (lines 14-111)
/// HiGHS 1.16 (#3153): feasibility_bounded.hpp::computeStartingPointBounded (lines 14-169)
pub fn compute_starting_point_bounded(
    instance: &mut Instance,
    settings: &mut Settings,
    stats: &mut Statistics,
    modelstatus: &mut QpModelStatus,
    result: &mut QpHotstartInformation,
    timer: &mut HighsTimer,
) {
    // HiGHS 1.16 #3153 (h16_qp_bounded_start): flag read once at the gate
    let h16_qp_bounded_start: bool = settings.h16_qp_bounded_start;

    // compute initial feasible point for problems with bounds only (no general
    // linear constraints)

    let mut res: QpVector;
    if h16_qp_bounded_start {
        // HiGHS 1.16 (#3153): feasibility_bounded.hpp::computeStartingPointBounded (lines 22-112)
        // Solve  Qx + c = 0 --> x = -Q^-1c
        let dim: HighsInt = instance.num_var;
        res = instance.c.op_neg();
        let hessian: &MatrixBase = &instance.Q.mat;
        debug_assert!(res.dim == dim);
        if hessian.is_diagonal() {
            // Diagonal Hessian
            for i_row in 0..dim {
                let value: f64 = hessian.diagonal(i_row);
                if value <= 0.0 {
                    *modelstatus = QpModelStatus::kNonConvex;
                    return;
                }
                res.value[i_row as usize] /= value;
            }
        } else {
            // General Hessian: compute Cholesky factorization of Q
            // (C++ L.resize(dim * dim) then L.assign(dim * dim, 0): zero-filled)
            let mut l: Vec<f64> = vec![0.0; (dim * dim) as usize];
            // C++ L[row * dim + col]
            let ix = |row: HighsInt, col: HighsInt| (row * dim + col) as usize;

            // First copy the lower triangle of Q into L
            let mut num_entries: HighsInt = 0;
            let mut value: Vec<f64> = vec![0.0; dim as usize];
            let mut index: Vec<HighsInt> = vec![0; dim as usize];
            for i_col in 0..dim {
                hessian.get_column(i_col, &mut num_entries, &mut index, &mut value);
                for i_el in 0..num_entries {
                    let i_row: HighsInt = index[i_el as usize];
                    // Take the entries above or on the diagonal in column iCol of
                    // Q as the entries before or on the diagonal in row iCol of
                    // (row-wise) L
                    if i_row <= i_col {
                        l[ix(i_col, i_row)] = value[i_el as usize];
                    }
                }
            }
            // Now compute Cholesky factorization of L
            // (C++ printL() calls dropped: debug_printing is false)
            for i_row in 0..dim {
                for i_col in 0..=i_row {
                    let mut sum: f64 = 0.0;
                    for k in 0..i_col {
                        sum += l[ix(i_row, k)] * l[ix(i_col, k)];
                    }
                    if i_col < i_row {
                        let value: f64 = (l[ix(i_row, i_col)] - sum) / l[ix(i_col, i_col)];
                        l[ix(i_row, i_col)] = value;
                    } else {
                        let value: f64 = l[ix(i_col, i_col)] - sum;
                        if value <= 0.0 {
                            *modelstatus = QpModelStatus::kNonConvex;
                            return;
                        }
                        l[ix(i_col, i_col)] = value.sqrt();
                    }
                }
            }
            // Solve for Qx = -c
            // Solve Ly = -c
            for i_row in 0..dim {
                let mut sum: f64 = 0.0;
                for i_col in 0..i_row {
                    sum += res.value[i_col as usize] * l[ix(i_row, i_col)];
                }
                res.value[i_row as usize] = (res.value[i_row as usize] - sum) / l[ix(i_row, i_row)];
            }
            // Solve L^Tx = y
            let mut i_row: HighsInt = dim - 1;
            while i_row >= 0 {
                res.value[i_row as usize] /= l[ix(i_row, i_row)];
                for i_col in 0..i_row {
                    res.value[i_col as usize] -= res.value[i_row as usize] * l[ix(i_row, i_col)];
                }
                i_row -= 1;
            }
        }
    } else {
        // compute  Qx + c = 0 --> x = Q^-1c

        let mut l: Vec<f64> = vec![0.0; (instance.num_var as usize) * (instance.num_var as usize)];

        // compute cholesky factorization of Q
        for col in 0..(instance.num_var as usize) {
            for idx in
                (instance.Q.mat.start[col] as usize)..(instance.Q.mat.start[col + 1] as usize)
            {
                let mut sum: f64 = 0.0;
                let row: usize = instance.Q.mat.index[idx] as usize;
                if row == col {
                    for k in 0..row {
                        sum += l[k * (instance.num_var as usize) + row]
                            * l[k * (instance.num_var as usize) + row];
                    }
                    l[row * (instance.num_var as usize) + row] =
                        (instance.Q.mat.value[idx] - sum).sqrt();
                } else {
                    for k in 0..row {
                        sum += l[k * (instance.num_var as usize) + col]
                            * l[k * (instance.num_var as usize) + row];
                    }
                    l[row * (instance.num_var as usize) + col] = (instance.Q.mat.value[idx] - sum)
                        / l[row * (instance.num_var as usize) + row];
                }
            }
        }

        // solve for c
        res = instance.c.op_neg();
        for r in 0..res.dim {
            for j in 0..r {
                res.value[r as usize] -= res.value[j as usize]
                    * l[(j as usize) * (instance.num_var as usize) + (r as usize)];
            }
            res.value[r as usize] /= l[(r as usize) * (instance.num_var as usize) + (r as usize)];
        }

        let mut i: HighsInt = res.dim - 1;
        while i >= 0 {
            let mut sum: f64 = 0.0;
            let mut j: HighsInt = res.dim - 1;
            while j > i {
                sum += res.value[j as usize]
                    * l[(i as usize) * (instance.num_var as usize) + (j as usize)];
                j -= 1;
            }
            res.value[i as usize] = (res.value[i as usize] - sum)
                / l[(i as usize) * (instance.num_var as usize) + (i as usize)];
            i -= 1;
        }
    }

    // project solution to bounds and collect active bounds
    let mut x0 = QpVector::new(instance.num_var);
    let ra = QpVector::new(instance.num_con);
    let mut initialactive: Vec<HighsInt> = Vec::new();
    let mut initialinactive: Vec<HighsInt> = Vec::new();
    let mut atlower: Vec<BasisStatus> = Vec::new();

    for i in 0..instance.num_var {
        if h16_qp_bounded_start {
            // HiGHS 1.16 (#3153): feasibility_bounded.hpp::computeStartingPointBounded (lines 121-149)
            // Check for unboundedness so that qp-unbounded and test-qod pass:
            // assumes diagonal Hessian and
            // settings.hessian_regularization_value > 0
            if settings.hessian_regularization_value > 0.0 {
                if res.value[i as usize] > 0.5 / settings.hessian_regularization_value
                    && instance.var_up[i as usize] == kHighsInf
                    && instance.c.value[i as usize] < 0.0
                {
                    *modelstatus = QpModelStatus::kUnbounded;
                    return;
                } else if res.value[i as usize] < -0.5 / settings.hessian_regularization_value
                    && instance.var_lo[i as usize] == -kHighsInf
                    && instance.c.value[i as usize] > 0.0
                {
                    *modelstatus = QpModelStatus::kUnbounded;
                    return;
                }
            }
            if res.value[i as usize]
                < instance.var_lo[i as usize] - settings.primal_feasibility_tolerance
            {
                res.value[i as usize] = instance.var_lo[i as usize];
                initialactive.push(i + instance.num_con);
                atlower.push(BasisStatus::kActiveAtLower);
            } else if res.value[i as usize]
                > instance.var_up[i as usize] + settings.primal_feasibility_tolerance
            {
                res.value[i as usize] = instance.var_up[i as usize];
                initialactive.push(i + instance.num_con);
                atlower.push(BasisStatus::kActiveAtUpper);
            } else {
                initialinactive.push(i + instance.num_con);
            }
        } else if res.value[i as usize] > 0.5 / settings.hessian_regularization_value
            && instance.var_up[i as usize] == f64::INFINITY
            && instance.c.value[i as usize] < 0.0
        {
            *modelstatus = QpModelStatus::kUnbounded;
            return;
        } else if res.value[i as usize] < 0.5 / settings.hessian_regularization_value
            && instance.var_lo[i as usize] == f64::INFINITY
            && instance.c.value[i as usize] > 0.0
        {
            *modelstatus = QpModelStatus::kUnbounded;
            return;
        } else if res.value[i as usize] <= instance.var_lo[i as usize] {
            res.value[i as usize] = instance.var_lo[i as usize];
            initialactive.push(i + instance.num_con);
            atlower.push(BasisStatus::kActiveAtLower);
        } else if res.value[i as usize] >= instance.var_up[i as usize] {
            res.value[i as usize] = instance.var_up[i as usize];
            initialactive.push(i + instance.num_con);
            atlower.push(BasisStatus::kActiveAtUpper);
        } else {
            initialinactive.push(i + instance.num_con);
        }
        if res.value[i as usize].abs() > 1e-4 {
            x0.value[i as usize] = res.value[i as usize];
            x0.index[x0.num_nz as usize] = i;
            x0.num_nz += 1;
        }
    }

    // if no bounds are active, solution lies in the interior -> optimal
    if initialactive.len() == 0 {
        *modelstatus = QpModelStatus::kOptimal;
    }

    debug_assert!(((initialactive.len() + initialinactive.len()) as HighsInt) == instance.num_var);

    result.status = atlower;
    result.active = initialactive;
    result.inactive = initialinactive;
    result.primal = x0;
    result.rowact = ra;
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::qpsolver::*;

    /// min 0.5 (2 x0^2 + 4 x1^2) - 2 x0 - 4 x1, 0 <= x <= 10 and, when
    /// `with_row`, x0 + x1 <= 10. Optimum x = (1, 1), objective -3.
    fn box_qp(with_row: bool) -> Instance {
        let nc = if with_row { 1 } else { 0 };
        let mut inst = Instance::new(2, nc);
        inst.c = QpVector {
            num_nz: 2,
            dim: 2,
            index: vec![0, 1],
            value: vec![-2.0, -4.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 2,
            num_col: 2,
            start: vec![0, 1, 2],
            index: vec![0, 1],
            value: vec![2.0, 4.0],
        };
        inst.A.mat = if with_row {
            MatrixBase {
                num_row: 1,
                num_col: 2,
                start: vec![0, 1, 2],
                index: vec![0, 0],
                value: vec![1.0, 1.0],
            }
        } else {
            MatrixBase {
                num_row: 0,
                num_col: 2,
                start: vec![0, 0, 0],
                index: vec![],
                value: vec![],
            }
        };
        inst.con_lo = vec![-kHighsInf; nc as usize];
        inst.con_up = vec![10.0; nc as usize];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        inst
    }

    /// Both variables active at their lower bound (constraint ids
    /// `num_con + i`); the row (id 0) is not in the basis.
    fn at_lower(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kActiveAtLower,
            BasisStatus::kActiveAtLower,
        ];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    /// Both variables in the nullspace (inactive in basis): Z = I.
    fn all_free(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kInactiveInBasis,
            BasisStatus::kInactiveInBasis,
        ];
        Basis::new(rt, vec![], status, vec![1, 2])
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
    }

    #[test]
    fn bounded_start_is_within_bounds() {
        let mut inst = box_qp(false);
        let mut settings = Settings::default();
        let mut stats = Statistics::default();
        let mut status = QpModelStatus::kNotset;
        let mut result = QpHotstartInformation::new(2, 0);
        let mut timer = HighsTimer::new();
        compute_starting_point_bounded(
            &mut inst,
            &mut settings,
            &mut stats,
            &mut status,
            &mut result,
            &mut timer,
        );
        assert_ne!(status, QpModelStatus::kInfeasible);
        assert_eq!(result.active.len() + result.inactive.len(), 2);
        for i in 0..2 {
            assert!(result.primal.value[i] >= 0.0 && result.primal.value[i] <= 10.0);
        }
    }
}
