//! HiGHS source: presolve/ICrashUtil.h, presolve/ICrashUtil.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Free functions; reached by path
//! (`crate::presolve::icrash_util::*`), not in the prelude.
//! `minimizeComponentBreakpoints` is declared (ICrashUtil.h line 52) but
//! never defined: not ported.

use crate::prelude::*;

/// HiGHS: ICrashUtil.cpp::convertToMinimization (lines 23-28)
pub fn convert_to_minimization(lp: &mut HighsLp) {
    if lp.sense_ != ObjSense::kMinimize {
        for col in 0..lp.num_col_ {
            lp.col_cost_[col as usize] = -lp.col_cost_[col as usize];
        }
    }
}

/// HiGHS: ICrashUtil.cpp::isEqualityProblem (lines 30-35)
pub fn is_equality_problem(lp: &HighsLp) -> bool {
    for row in 0..lp.num_row_ {
        if lp.row_lower_[row as usize] != lp.row_upper_[row as usize] {
            return false;
        }
    }

    true
}

/// HiGHS: ICrashUtil.cpp::vectorProduct (lines 37-43)
pub fn vector_product(v1: &[f64], v2: &[f64]) -> f64 {
    debug_assert!(v1.len() == v2.len());
    let mut sum = 0.0;
    for i in 0..(v1.len() as i32) {
        sum += v1[i as usize] * v2[i as usize];
    }
    sum
}

/// HiGHS: ICrashUtil.cpp::muptiplyByTranspose (lines 45-58)
pub fn muptiply_by_transpose(lp: &HighsLp, v: &[f64], result: &mut Vec<f64>) {
    debug_assert!(result.len() as i32 == lp.num_col_);
    debug_assert!(v.len() as i32 == lp.num_row_);

    result.clear();
    result.resize(lp.num_col_ as usize, 0.0);
    for col in 0..lp.num_col_ {
        for k in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
            let row = lp.a_matrix_.index_[k as usize];
            result[col as usize] += lp.a_matrix_.value_[k as usize] * v[row as usize];
        }
    }
}

/// HiGHS: ICrashUtil.cpp::printMinorIterationDetails (lines 60-79)
pub fn print_minor_iteration_details(
    iteration: f64,
    col: f64,
    old_value: f64,
    update: f64,
    ctx: f64,
    r: &[f64],
    quadratic_objective: f64,
    options: HighsLogOptions,
) {
    let rnorm = get_norm2(r);
    // C++ lines 66-75 build a std::stringstream; every value is a double
    // streamed with the default ostream state (%g, precision 6).
    // std::endl is "\n" on the stream.
    let g = |v: f64| crate::util::c_printf::c_fmt_g(v, 6);
    let message = format!(
        "iter {}, col {}, update {}, old_value {}, new_value {}, ctx {}, r {}, quadratic_objective {}\n",
        g(iteration),
        g(col),
        g(update),
        g(old_value),
        g(old_value + update),
        g(ctx),
        g(rnorm),
        g(quadratic_objective)
    );

    // C++ line 77: // std::cout << ss.str(); dropped (commented-out code).
    highs_log_user(&options, HighsLogType::kInfo, &message);
}

/// HiGHS: ICrashUtil.cpp::initialize (lines 81-104)
pub fn initialize(lp: &HighsLp, solution: &mut HighsSolution, lambda: &mut Vec<f64>) -> bool {
    // Clear and resize primal column solution.
    solution.clear();
    solution.col_value.resize(lp.num_col_ as usize, 0.0);

    for col in 0..lp.num_col_ {
        if lp.col_lower_[col as usize] <= 0.0 && lp.col_upper_[col as usize] >= 0.0 {
            solution.col_value[col as usize] = 0.0;
        } else if lp.col_lower_[col as usize] > 0.0 {
            solution.col_value[col as usize] = lp.col_lower_[col as usize];
        } else if lp.col_upper_[col as usize] < 0.0 {
            solution.col_value[col as usize] = lp.col_upper_[col as usize];
        } else {
            // C++ line 95: printf to stdout.
            println!("ICrash error: setting initial value for column {}", col);
            return false;
        }
    }

    lambda.resize(lp.num_row_ as usize, 0.0);
    lambda.clear();
    lambda.resize(lp.num_row_ as usize, 0.0);

    true
}

/// HiGHS: ICrashUtil.cpp::minimizeComponentQP (lines 106-164)
pub fn minimize_component_qp(
    col: i32,
    mu: f64,
    lp: &HighsLp,
    objective: &mut f64,
    residual: &mut Vec<f64>,
    sol: &mut HighsSolution,
) -> f64 {
    // Minimize quadratic for column col.

    // Formulas for a and b when minimizing for x_j
    // a = (1/(2*mu)) * sum_i a_ij^2
    // b = -(1/(2*mu)) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i))
    // b / 2 = -(1/(2*mu)) sum_i (2 * a_ij
    let mut a = 0.0;
    let mut b = 0.0;

    for k in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
        let row = lp.a_matrix_.index_[k as usize];
        a += lp.a_matrix_.value_[k as usize] * lp.a_matrix_.value_[k as usize];
        // matlab but with b = b / 2
        let bracket =
            -residual[row as usize] - lp.a_matrix_.value_[k as usize] * sol.col_value[col as usize];
        // C++ lines 126-127: clp minimizing for delta_x (commented-out code dropped).
        b += lp.a_matrix_.value_[k as usize] * bracket;
    }

    a = (0.5 / mu) * a;
    b = (0.5 / mu) * b + 0.5 * lp.col_cost_[col as usize];

    let theta = -b / a;

    // matlab
    let new_x = if theta > 0.0 {
        cpp_min(theta, lp.col_upper_[col as usize])
    } else {
        cpp_max(theta, lp.col_lower_[col as usize])
    };
    let delta_x = new_x - sol.col_value[col as usize];

    // C++ lines 144-148: clp minimizing for delta_x (commented-out code dropped).

    sol.col_value[col as usize] += delta_x;

    // C++ line 152: std::cout dropped (commented-out code).

    // Update objective, row_value, residual after each component update.
    *objective += lp.col_cost_[col as usize] * delta_x;
    for k in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
        let row = lp.a_matrix_.index_[k as usize];
        sol.row_value[row as usize] += lp.a_matrix_.value_[k as usize] * delta_x;
        residual[row as usize] = (lp.row_upper_[row as usize] - sol.row_value[row as usize]).abs();
    }

    delta_x
}

/// HiGHS: ICrashUtil.cpp::minimizeComponentIca (lines 166-228)
pub fn minimize_component_ica(
    col: i32,
    mu: f64,
    lambda: &[f64],
    lp: &HighsLp,
    objective: &mut f64,
    residual: &mut Vec<f64>,
    sol: &mut HighsSolution,
) -> f64 {
    // Minimize quadratic for column col.

    // Formulas for a and b when minimizing for x_j
    // a = (1/(2*mu)) * sum_i a_ij^2
    // b = -(1/(2*mu)) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i)) + c_j (\)
    //     + sum_i a_ij * lambda_i
    // b / 2 = -(1/(2*mu)) sum_i (2 * a_ij
    let mut a = 0.0;
    let mut b = 0.0;

    for k in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
        let row = lp.a_matrix_.index_[k as usize];
        a += lp.a_matrix_.value_[k as usize] * lp.a_matrix_.value_[k as usize];
        // matlab but with b = b / 2
        let mut bracket =
            -residual[row as usize] - lp.a_matrix_.value_[k as usize] * sol.col_value[col as usize];
        bracket += lambda[row as usize];
        // C++ lines 188-189: clp minimizing for delta_x (commented-out code dropped).
        b += lp.a_matrix_.value_[k as usize] * bracket;
    }

    a = (0.5 / mu) * a;
    b = (0.5 / mu) * b + 0.5 * lp.col_cost_[col as usize];

    let theta = -b / a;

    // matlab
    let new_x = if theta > 0.0 {
        cpp_min(theta, lp.col_upper_[col as usize])
    } else {
        cpp_max(theta, lp.col_lower_[col as usize])
    };
    let delta_x = new_x - sol.col_value[col as usize];

    // C++ lines 207-211: clp minimizing for delta_x (commented-out code dropped).

    sol.col_value[col as usize] += delta_x;

    // C++ line 215: std::cout dropped (commented-out code).

    // Update objective, row_value, residual after each component update.
    *objective += lp.col_cost_[col as usize] * delta_x;
    for k in lp.a_matrix_.start_[col as usize]..lp.a_matrix_.start_[(col + 1) as usize] {
        let row = lp.a_matrix_.index_[k as usize];
        residual[row as usize] -= lp.a_matrix_.value_[k as usize] * delta_x;
        sol.row_value[row as usize] += lp.a_matrix_.value_[k as usize] * delta_x;
        // C++ line 224: residual[row] = fabs(...) dropped (commented-out code).
    }

    delta_x
}

/// HiGHS: ICrashUtil.cpp::updateResidual (lines 230-251)
pub fn update_residual(
    piecewise: bool,
    lp: &HighsLp,
    sol: &HighsSolution,
    residual: &mut Vec<f64>,
) {
    residual.clear();
    residual.resize(lp.num_row_ as usize, 0.0);

    if !piecewise {
        debug_assert!(is_equality_problem(lp));
        for row in 0..lp.num_row_ {
            residual[row as usize] =
                (lp.row_upper_[row as usize] - sol.row_value[row as usize]).abs();
        }
    } else {
        // piecewise
        for row in 0..lp.num_row_ {
            let mut value = 0.0;
            if sol.row_value[row as usize] <= lp.row_lower_[row as usize] {
                value = lp.row_lower_[row as usize] - sol.row_value[row as usize];
            } else if sol.row_value[row as usize] >= lp.row_upper_[row as usize] {
                value = sol.row_value[row as usize] - lp.row_upper_[row as usize];
            }

            residual[row as usize] = value;
        }
    }
}

/// HiGHS: ICrashUtil.cpp::updateResidualFast (lines 253-259)
pub fn update_residual_fast(lp: &HighsLp, sol: &HighsSolution, residual: &mut Vec<f64>) {
    debug_assert!(is_equality_problem(lp));
    for row in 0..lp.num_row_ {
        residual[row as usize] = (lp.row_upper_[row as usize] - sol.row_value[row as usize]).abs();
    }
}

/// HiGHS: ICrashUtil.cpp::updateResidualIca (lines 262-267)
pub fn update_residual_ica(lp: &HighsLp, sol: &HighsSolution, residual: &mut Vec<f64>) {
    debug_assert!(is_equality_problem(lp));
    for row in 0..lp.num_row_ {
        residual[row as usize] = lp.row_upper_[row as usize] - sol.row_value[row as usize];
    }
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

    /// Column-wise LP builder: `cols[j] = (cost, lower, upper, [(row, value)])`.
    #[allow(clippy::type_complexity)]
    fn lp(cols: &[(f64, f64, f64, &[(HighsInt, f64)])], rows: &[(f64, f64)]) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = cols.len() as HighsInt;
        lp.num_row_ = rows.len() as HighsInt;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = lp.num_col_;
        lp.a_matrix_.num_row_ = lp.num_row_;
        lp.a_matrix_.start_ = vec![0];
        for &(cost, lower, upper, entries) in cols {
            lp.col_cost_.push(cost);
            lp.col_lower_.push(lower);
            lp.col_upper_.push(upper);
            for &(r, v) in entries {
                lp.a_matrix_.index_.push(r);
                lp.a_matrix_.value_.push(v);
            }
            lp.a_matrix_
                .start_
                .push(lp.a_matrix_.index_.len() as HighsInt);
        }
        for &(l, u) in rows {
            lp.row_lower_.push(l);
            lp.row_upper_.push(u);
        }
        lp
    }

    fn _quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    #[test]
    fn small_vector_helpers() {
        assert_eq!(vector_product(&[1.0, 2.0], &[3.0, 4.0]), 11.0);
        // [[1, 1], [1, -1]] column-wise; A^T (1, 2) = (3, -1)
        let mut l = lp(
            &[
                (1.0, 0.0, kHighsInf, &[(0, 1.0), (1, 1.0)]),
                (2.0, 0.0, kHighsInf, &[(0, 1.0), (1, -1.0)]),
            ],
            &[(1.0, 1.0), (0.0, 0.0)],
        );
        let mut r = vec![0.0; 2];
        muptiply_by_transpose(&l, &[1.0, 2.0], &mut r);
        assert_eq!(r, vec![3.0, -1.0]);
        assert!(is_equality_problem(&l));
        l.sense_ = ObjSense::kMaximize;
        convert_to_minimization(&mut l);
        assert_eq!(l.col_cost_, vec![-1.0, -2.0]);
        l.row_upper_[1] = 1.0;
        assert!(!is_equality_problem(&l));
    }
}
