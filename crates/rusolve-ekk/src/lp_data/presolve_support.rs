//! HiGHS source: lp_data/HighsLpUtils.cpp (cleanBounds, calculateRowValuesQuad,
//! reportPresolveReductions), lp_data/HighsSolution.cpp (refineBasis,
//! isBasisConsistent, getKktFailures, getLpKktFailures, lpKktCheck,
//! reportKktFailures), lp_data/HighsSolutionDebug.cpp
//! (debugHighsSolution), util/HighsSparseMatrix.cpp (assessSmallValues)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The free functions the presolve branch of `Highs::run`
//! (`highs_presolve.rs`) calls. Everything is filled except `refine_basis`
//! (translator task) and the `debug_highs_solution` stub.

#![allow(unused_variables)]

use crate::lp_data::highs_solution::is_basis_right_size;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g_width;

/// HiGHS: HighsLpUtils.cpp::cleanBounds (lines 724-770, 47 lines)
pub fn clean_bounds(options: &HighsOptions, lp: &mut HighsLp) -> HighsStatus {
    let mut max_residual: f64 = 0.0;
    let mut num_change: HighsInt = 0;
    for i_col in 0..lp.num_col_ as usize {
        let residual = lp.col_lower_[i_col] - lp.col_upper_[i_col];
        if residual > options.primal_feasibility_tolerance {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "Column {} has inconsistent bounds [{}, {}] (residual = {}) after presolve\n",
                i_col,
                lp.col_lower_[i_col],
                lp.col_upper_[i_col],
                residual
            );
            return HighsStatus::kError;
        } else if residual > 0.0 {
            num_change += 1;
            max_residual = cpp_max(residual, max_residual);
            let mid = 0.5 * (lp.col_lower_[i_col] + lp.col_upper_[i_col]);
            lp.col_lower_[i_col] = mid;
            lp.col_upper_[i_col] = mid;
        }
    }
    for i_row in 0..lp.num_row_ as usize {
        let residual = lp.row_lower_[i_row] - lp.row_upper_[i_row];
        if residual > options.primal_feasibility_tolerance {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kError,
                "Row {} has inconsistent bounds [{}, {}] (residual = {}) after presolve\n",
                i_row,
                lp.row_lower_[i_row],
                lp.row_upper_[i_row],
                residual
            );
            return HighsStatus::kError;
        } else if residual > 0.0 {
            num_change += 1;
            max_residual = cpp_max(residual, max_residual);
            let mid = 0.5 * (lp.row_lower_[i_row] + lp.row_upper_[i_row]);
            lp.row_lower_[i_row] = mid;
            lp.row_upper_[i_row] = mid;
        }
    }
    if num_change != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "Resolved {} inconsistent bounds (maximum residual = {:9.4e}) after presolve\n",
            num_change,
            max_residual
        );
        return HighsStatus::kWarning;
    }
    HighsStatus::kOk
}

/// HiGHS: HighsLpUtils.cpp::calculateRowValuesQuad(lp, col_value, row_value,
/// report_row) (lines 3011-3047, 37 lines)
///
/// C++ default arguments: `report_row = -1`.
pub fn calculate_row_values_quad_vec(
    lp: &HighsLp,
    col_value: &[f64],
    row_value: &mut Vec<f64>,
    report_row: HighsInt,
) -> HighsStatus {
    let correct_size = col_value.len() as HighsInt == lp.num_col_;
    let is_colwise = lp.a_matrix_.is_colwise();
    let data_error = !correct_size || !is_colwise;
    debug_assert!(!data_error);
    if data_error {
        return HighsStatus::kError;
    }

    let mut row_value_quad: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); lp.num_row_ as usize];

    for col in 0..lp.num_col_ as usize {
        for i in lp.a_matrix_.start_[col] as usize..lp.a_matrix_.start_[col + 1] as usize {
            let row = lp.a_matrix_.index_[i];
            debug_assert!(row >= 0);
            debug_assert!(row < lp.num_row_);
            row_value_quad[row as usize] +=
                HighsCDouble::from(col_value[col]) * lp.a_matrix_.value_[i];
            if row == report_row {
                println!(
                    "calculateRowValuesQuad: Row {} becomes {} due to contribution of .col_value[{}] = {}",
                    row,
                    f64::from(row_value_quad[row as usize]),
                    col,
                    col_value[col]
                );
            }
        }
    }

    // assign quad values to double vector
    row_value.resize(lp.num_row_ as usize, 0.0);
    for (dst, src) in row_value.iter_mut().zip(row_value_quad.iter()) {
        *dst = f64::from(*src);
    }

    HighsStatus::kOk
}

/// HiGHS: HighsLpUtils.cpp::calculateRowValuesQuad(lp, solution, report_row)
/// (lines 3049-3053, 5 lines)
///
/// C++ default arguments: `report_row = -1`.
pub fn calculate_row_values_quad(
    lp: &HighsLp,
    solution: &mut HighsSolution,
    report_row: HighsInt,
) -> HighsStatus {
    let col_value = std::mem::take(&mut solution.col_value);
    let status = calculate_row_values_quad_vec(lp, &col_value, &mut solution.row_value, report_row);
    solution.col_value = col_value;
    status
}

/// HiGHS: HighsLpUtils.cpp::reportPresolveReductions (lines 3100-3160, 61 lines)
pub fn report_presolve_reductions(
    log_options: &HighsLogOptions,
    presolve_status: HighsPresolveStatus,
    lp: &HighsLp,
    presolved_lp: &HighsLp,
) {
    let num_col_from = lp.num_col_;
    let num_row_from = lp.num_row_;
    let num_nz_from = lp.a_matrix_.num_nz();
    let num_col_to: HighsInt;
    let num_row_to: HighsInt;
    let num_nz_to: HighsInt;
    let message: &str;

    match presolve_status {
        HighsPresolveStatus::kNotPresolved
        | HighsPresolveStatus::kInfeasible
        | HighsPresolveStatus::kUnboundedOrInfeasible => return,
        HighsPresolveStatus::kNotReduced => {
            num_col_to = num_col_from;
            num_row_to = num_row_from;
            num_nz_to = num_nz_from;
            message = "- Not reduced";
        }
        HighsPresolveStatus::kReduced | HighsPresolveStatus::kTimeout => {
            num_col_to = presolved_lp.num_col_;
            num_row_to = presolved_lp.num_row_;
            num_nz_to = presolved_lp.a_matrix_.num_nz();
            message = if presolve_status == HighsPresolveStatus::kTimeout {
                "- Timeout"
            } else {
                ""
            };
        }
        HighsPresolveStatus::kReducedToEmpty => {
            num_col_to = 0;
            num_row_to = 0;
            num_nz_to = 0;
            message = "- Reduced to empty";
        }
        _ => {
            // case HighsPresolveStatus::kOutOfMemory
            debug_assert!(presolve_status == HighsPresolveStatus::kOutOfMemory);
            return;
        }
    }
    let mut nz_sign_char = '-';
    let mut delta_nz = num_nz_from - num_nz_to;
    if num_nz_to > num_nz_from {
        delta_nz = -delta_nz;
        nz_sign_char = '+';
    }
    crate::highs_log_user!(
        log_options,
        HighsLogType::kInfo,
        "Presolve reductions: rows {}(-{}); columns {}(-{}); nonzeros {}({}{}) {}\n",
        num_row_to,
        num_row_from - num_row_to,
        num_col_to,
        num_col_from - num_col_to,
        num_nz_to,
        nz_sign_char,
        delta_nz,
        message
    );
}

/// HiGHS: HighsSolution.cpp::refineBasis (lines 1423-1498, 76 lines)
///
/// Replaces every `HighsBasisStatus::kNonbasic` entry by kLower/kUpper/kZero
/// from the bounds (and, if `solution.value_valid`, the primal value).
/// `fabs` is `f64::abs`; `highs_isInfinity` is `highs_is_infinity`.
/// rusolve `lp_postsolve_snap_basis`: postsolve may leave a nonbasic
/// boxed variable with the status of one bound and a value at the other (Primal2_1000: 568408
/// columns with range 1e-8, status upper, value at lower). The simplex then moves them to the
/// status bound and starts with a basic solution that is primal infeasible in 576664 places.
/// Give such a variable the status of the bound its value is at.
pub fn snap_postsolve_basis(lp: &HighsLp, solution: &HighsSolution, basis: &mut HighsBasis) {
    if !solution.value_valid {
        return;
    }
    let snap = |status: &mut HighsBasisStatus, lower: f64, upper: f64, value: f64| {
        if lower == upper || highs_is_infinity(-lower) || highs_is_infinity(upper) {
            return;
        }
        let (at_lower, at_upper) = ((value - lower).abs(), (value - upper).abs());
        match *status {
            HighsBasisStatus::kUpper if at_lower < at_upper => *status = HighsBasisStatus::kLower,
            HighsBasisStatus::kLower if at_upper < at_lower => *status = HighsBasisStatus::kUpper,
            _ => {}
        }
    };
    for j in 0..lp.num_col_ as usize {
        snap(
            &mut basis.col_status[j],
            lp.col_lower_[j],
            lp.col_upper_[j],
            solution.col_value[j],
        );
    }
    for i in 0..lp.num_row_ as usize {
        snap(
            &mut basis.row_status[i],
            lp.row_lower_[i],
            lp.row_upper_[i],
            solution.row_value[i],
        );
    }
}

/// rusolve `lp_postsolve_fix_tiny`: for the clean-up solve after postsolve,
/// fix every column whose range is positive and within `tol` at its lower bound, give a nonbasic
/// one the status of the lower bound, and return the original upper bounds for
/// `restore_tiny_range_cols`. The fixed LP is an inner restriction of the original (the lower
/// bound is feasible), so the cleaned-up solution is feasible for the original LP, and the
/// objective changes by at most tol times the sum of the cost magnitudes of these columns.
/// Primal2_1000: 568408 columns of range 1e-8; recovered at their lower bound but marked upper by
/// postsolve, which gave 576664 primal infeasibilities (sum 72) to the dual simplex.
pub fn fix_tiny_range_cols(
    lp: &mut HighsLp,
    basis: &mut HighsBasis,
    tol: f64,
) -> Vec<(usize, f64)> {
    let mut saved = Vec::new();
    for j in 0..lp.num_col_ as usize {
        let (lower, upper) = (lp.col_lower_[j], lp.col_upper_[j]);
        if lower < upper && upper - lower <= tol {
            saved.push((j, upper));
            lp.col_upper_[j] = lower;
            if basis.valid && basis.col_status[j] == HighsBasisStatus::kUpper {
                basis.col_status[j] = HighsBasisStatus::kLower;
            }
        }
    }
    saved
}

/// Undo `fix_tiny_range_cols` on the model bounds.
pub fn restore_tiny_range_cols(lp: &mut HighsLp, saved: &[(usize, f64)]) {
    for &(j, upper) in saved {
        lp.col_upper_[j] = upper;
    }
}

pub fn refine_basis(lp: &HighsLp, solution: &HighsSolution, basis: &mut HighsBasis) {
    debug_assert!(basis.useful);
    debug_assert!(is_basis_right_size(lp, basis));
    let have_highs_solution = solution.value_valid;

    let num_col = lp.num_col_;
    let num_row = lp.num_row_;
    for i_col in 0..num_col as usize {
        if basis.col_status[i_col] != HighsBasisStatus::kNonbasic {
            continue;
        }
        let lower = lp.col_lower_[i_col];
        let upper = lp.col_upper_[i_col];
        let status: HighsBasisStatus;
        if lower == upper {
            status = HighsBasisStatus::kLower;
        } else if !highs_is_infinity(-lower) {
            if !highs_is_infinity(upper) {
                if have_highs_solution {
                    if solution.col_value[i_col] < 0.5 * (lower + upper) {
                        status = HighsBasisStatus::kLower;
                    } else {
                        status = HighsBasisStatus::kUpper;
                    }
                } else {
                    if lower.abs() < upper.abs() {
                        status = HighsBasisStatus::kLower;
                    } else {
                        status = HighsBasisStatus::kUpper;
                    }
                }
            } else {
                status = HighsBasisStatus::kLower;
            }
        } else if !highs_is_infinity(upper) {
            status = HighsBasisStatus::kUpper;
        } else {
            status = HighsBasisStatus::kZero;
        }
        debug_assert!(status != HighsBasisStatus::kNonbasic);
        basis.col_status[i_col] = status;
    }

    for i_row in 0..num_row as usize {
        if basis.row_status[i_row] != HighsBasisStatus::kNonbasic {
            continue;
        }
        let lower = lp.row_lower_[i_row];
        let upper = lp.row_upper_[i_row];
        let status: HighsBasisStatus;
        if lower == upper {
            status = HighsBasisStatus::kLower;
        } else if !highs_is_infinity(-lower) {
            if !highs_is_infinity(upper) {
                if have_highs_solution {
                    if solution.row_value[i_row] < 0.5 * (lower + upper) {
                        status = HighsBasisStatus::kLower;
                    } else {
                        status = HighsBasisStatus::kUpper;
                    }
                } else {
                    if lower.abs() < upper.abs() {
                        status = HighsBasisStatus::kLower;
                    } else {
                        status = HighsBasisStatus::kUpper;
                    }
                }
            } else {
                status = HighsBasisStatus::kLower;
            }
        } else if !highs_is_infinity(upper) {
            status = HighsBasisStatus::kUpper;
        } else {
            status = HighsBasisStatus::kZero;
        }
        debug_assert!(status != HighsBasisStatus::kNonbasic);
        basis.row_status[i_row] = status;
    }
}

/// HiGHS: HighsSolution.cpp::isBasisConsistent (lines 1989-2002, 14 lines)
pub fn is_basis_consistent(lp: &HighsLp, basis: &HighsBasis) -> bool {
    if !is_basis_right_size(lp, basis) {
        return false;
    }

    let mut num_basic_variables: HighsInt = 0;
    for i_col in 0..lp.num_col_ as usize {
        if basis.col_status[i_col] == HighsBasisStatus::kBasic {
            num_basic_variables += 1;
        }
    }
    for i_row in 0..lp.num_row_ as usize {
        if basis.row_status[i_row] == HighsBasisStatus::kBasic {
            num_basic_variables += 1;
        }
    }
    num_basic_variables == lp.num_row_
}

impl HighsSparseMatrix {
    /// HiGHS: HighsSparseMatrix.cpp::HighsSparseMatrix::assessSmallValues (lines 815-824, 10 lines)
    ///
    /// `analyseVectorValues` is a no-op stub in this crate (report only).
    pub fn assess_small_values(&self, log_options: &HighsLogOptions, small_matrix_value: f64) {
        let mut min_value = kHighsInf;
        let num_values = self.value_.len() as HighsInt;
        for i_x in 0..num_values as usize {
            min_value = cpp_min(self.value_[i_x].abs(), min_value);
        }
        if min_value > small_matrix_value {
            return;
        }
        analyse_vector_values(
            Some(log_options),
            "Small values in matrix",
            num_values,
            &self.value_,
            false,
            "",
        );
    }
}

/// HiGHS: HighsSolutionDebug.cpp::debugHighsSolution(message, options, model,
/// solution, basis) (lines 44-68)
///
/// Debug check: out of scope; stub returns `kNotChecked`, which is what
/// HiGHS returns for `highs_debug_level < kHighsDebugLevelCheap` (the
/// default). `model` is the LP (no Hessian).
pub fn debug_highs_solution(
    message: &str,
    options: &HighsOptions,
    lp: &HighsLp,
    solution: &HighsSolution,
    basis: &HighsBasis,
) -> HighsDebugStatus {
    HighsDebugStatus::kNotChecked
}

/// HiGHS: HighsSolution.cpp::kHighsSolutionLo / No / Up (lines 24-26):
/// `const uint8_t kHighsSolutionLo = -1;` wraps to 255.
const kHighsSolutionLo: u8 = u8::MAX;
const kHighsSolutionNo: u8 = 0;
const kHighsSolutionUp: u8 = 1;

/// HiGHS: HighsSolution.cpp::printf_kkt (line 28)
const printf_kkt: bool = false;

/// HiGHS: HighsUtils.h::infeasibility (lines 219-271)
/// C++: `inline std::pair<double, double> infeasibility(const double lower, const double value, const double upper, const double tolerance)`
pub fn infeasibility(lower: f64, value: f64, upper: f64, tolerance: f64) -> (f64, f64) {
    let mut residual: f64 = 0.0;
    let mut infeasibility: f64 = 0.0;
    if value < lower - tolerance {
        infeasibility = lower - value;
    }
    if value > upper + tolerance {
        infeasibility = value - upper;
    }
    if tolerance > 0.0 {
        if value < lower {
            residual = lower - value;
        }
        if value > upper {
            residual = value - upper;
        }
    } else {
        residual = infeasibility;
    }
    if infeasibility == 0.0 {
        residual = cpp_min(residual, tolerance);
    }
    (infeasibility, residual)
}

/// HiGHS 1.16.0 (#3299, commit cd50bf7f1): the `getEffectiveCosts` call that
/// 1.16.0 `getKktFailures` makes (`std::vector<double> effective_costs =
/// getEffectiveCosts(lp, options);`, HighsSolution.cpp line 169).
///
/// C++ passes `lp` by const reference and `getEffectiveCosts` asserts that the
/// matrix is column-wise (`assert(lp.a_matrix_.isColwise())`,
/// HighsSolution.cpp line 705, 1.16.0). Here the LP is only borrowed and
/// cannot be reformatted in place, so a row-wise matrix is column-wise
/// transposed into a local copy first: the same LP data, and the case the C++
/// assert rules out. Transposition can reorder the entries within a column, so
/// the accumulated effective cost of a column lying in several selected rows
/// can differ in its last bits from a matrix that was column-wise already.
fn kkt_effective_costs(lp: &HighsLp, options: &HighsOptions) -> Vec<f64> {
    if lp.a_matrix_.is_colwise() {
        return crate::lp_data::highs_solution::get_effective_costs(lp, options);
    }
    let mut lp_colwise = lp.clone();
    lp_colwise.ensure_colwise();
    crate::lp_data::highs_solution::get_effective_costs(&lp_colwise, options)
}

/// HiGHS: HighsSolution.cpp::getKktFailures(options, is_qp, lp, gradient,
/// solution, highs_info, get_residuals) (lines 73-561)
///
/// Called from `Highs::callSolveMip` through the model overload (lines 38-51), whose
/// `getPrimalDualBasisErrors` / `getPrimalDualGlpsolErrors` write only the discarded
/// `primal_dual_errors` and are not ported.
///
/// 1.16.0 also promotes its two local norms to new `HighsInfo` fields
/// `active_cost_norm` / `active_bound_norm` (HighsSolution.cpp lines 150-151,
/// 1.16.0, read by the reports at lines 1434 and 2211) and logs them at line
/// 522. That rename leaves every measure unchanged, so it is not ported: only
/// the effective-costs use of the cost norm is, behind
/// `h16_kkt_effective_costs`.
pub fn get_kkt_failures(
    options: &HighsOptions,
    is_qp: bool,
    lp: &HighsLp,
    gradient: &[f64],
    solution: &HighsSolution,
    highs_info: &mut HighsInfo,
    get_residuals: bool,
) {
    let mut primal_feasibility_tolerance = options.primal_feasibility_tolerance;
    let mut dual_feasibility_tolerance = options.dual_feasibility_tolerance;
    let mip_feasibility_tolerance = options.mip_feasibility_tolerance;
    let mut primal_residual_tolerance = options.primal_residual_tolerance;
    let mut dual_residual_tolerance = options.dual_residual_tolerance;
    let mut optimality_tolerance = options.optimality_tolerance;
    if options.kkt_tolerance != kDefaultKktTolerance {
        primal_feasibility_tolerance = options.kkt_tolerance;
        dual_feasibility_tolerance = options.kkt_tolerance;
        primal_residual_tolerance = options.kkt_tolerance;
        dual_residual_tolerance = options.kkt_tolerance;
        optimality_tolerance = options.kkt_tolerance;
    }
    // PORT NOTE: the C++ binds references to the highs_info fields; the port writes
    // `highs_info.<field>` directly.

    let have_primal_solution = solution.value_valid;
    let have_dual_solution = solution.dual_valid;
    let have_integrality = !lp.integrality_.is_empty();

    highs_info.primal_solution_status = kSolutionStatusNone;
    highs_info.dual_solution_status = kSolutionStatusNone;

    debug_assert!(have_primal_solution || !have_dual_solution);

    highs_info.invalidate_kkt();

    // HiGHS 1.16.0 (#3299, commit cd50bf7f1): HighsSolution.cpp line 169
    // computes the effective costs of the LP, used as `effective_cost` in the
    // cost-norm measures below (lines 274, 283-284, 397, 433). Behind
    // `h16_kkt_effective_costs`; with it off `get_effective_costs` is not run at
    // all (it costs passes over the matrix) and the raw gradient cost is used,
    // which is what 1.15.1 does.
    let h16_kkt_effective_costs = options.h16_kkt_effective_costs;
    let effective_costs_owned: Vec<f64> = if h16_kkt_effective_costs {
        kkt_effective_costs(lp, options)
    } else {
        Vec::new()
    };
    // With the option off the effective cost of a column is its cost, so
    // `effective_cost` below is bit-identical to `cost`.
    let effective_costs: &[f64] = if h16_kkt_effective_costs {
        &effective_costs_owned
    } else {
        gradient
    };

    if have_primal_solution {
        debug_assert!(solution.col_value.len() as HighsInt >= lp.num_col_);
        debug_assert!(solution.row_value.len() as HighsInt >= lp.num_row_);
        highs_info.num_primal_infeasibilities = 0;
        highs_info.max_primal_infeasibility = 0.0;
        highs_info.sum_primal_infeasibilities = 0.0;
        highs_info.num_semi_infeasibilities = 0;
        highs_info.max_semi_infeasibility = 0.0;
        highs_info.sum_semi_infeasibilities = 0.0;
        highs_info.num_relative_primal_infeasibilities = 0;
        highs_info.max_relative_primal_infeasibility = 0.0;
        if get_residuals {
            highs_info.num_primal_residual_errors = 0;
            highs_info.max_primal_residual_error = 0.0;
            highs_info.num_relative_primal_residual_errors = 0;
            highs_info.max_relative_primal_residual_error = 0.0;
        }
    }
    if have_dual_solution {
        debug_assert!(solution.col_dual.len() as HighsInt >= lp.num_col_);
        debug_assert!(solution.row_dual.len() as HighsInt >= lp.num_row_);
        highs_info.num_dual_infeasibilities = 0;
        highs_info.max_dual_infeasibility = 0.0;
        highs_info.sum_dual_infeasibilities = 0.0;
        highs_info.num_relative_dual_infeasibilities = 0;
        highs_info.max_relative_dual_infeasibility = 0.0;
        if get_residuals {
            highs_info.num_dual_residual_errors = 0;
            highs_info.max_dual_residual_error = 0.0;
            highs_info.num_relative_dual_residual_errors = 0;
            highs_info.max_relative_dual_residual_error = 0.0;
        }
    }
    // Without a primal solution, nothing can be done!
    if !have_primal_solution {
        return;
    }

    let mut primal_activity: Vec<f64> = Vec::new();
    let mut dual_activity: Vec<f64> = Vec::new();
    if get_residuals {
        lp.a_matrix_
            .product_quad(&mut primal_activity, &solution.col_value, kDebugReportOff);
    }
    if get_residuals && have_dual_solution {
        lp.a_matrix_.product_transpose_quad(
            &mut dual_activity,
            &solution.row_dual,
            kDebugReportOff,
        );
    }

    let mut max_col_primal_infeasibility: f64 = 0.0;
    let mut max_col_dual_infeasibility: f64 = 0.0;
    let mut max_relative_col_primal_infeasibility: f64 = 0.0;
    let mut max_relative_col_dual_infeasibility: f64 = 0.0;

    let mut primal_infeasibility: f64 = 0.0;
    let mut dual_infeasibility: f64 = 0.0;
    let mut semi_infeasibility: f64 = 0.0;
    let mut cost: f64 = 0.0;
    // HiGHS 1.16.0: HighsSolution.cpp line 230 (only assigned for columns, as in C++)
    let mut effective_cost: f64 = 0.0;
    let mut lower: f64;
    let mut upper: f64;
    let mut value: f64;
    let mut dual: f64 = 0.0;
    let mut at_status: u8 = kHighsSolutionNo;
    let mut mid_status: u8 = kHighsSolutionNo;
    let mut integrality = HighsVarType::kContinuous;
    let mut highs_norm_bounds: f64 = 0.0;
    let mut highs_norm_costs: f64 = 0.0;

    for pass in 0..2 {
        for i_var in 0..lp.num_col_ + lp.num_row_ {
            let is_col = i_var < lp.num_col_;
            if is_col {
                let i_col = i_var as usize;
                cost = gradient[i_col];
                // HiGHS 1.16.0: HighsSolution.cpp line 274
                effective_cost = effective_costs[i_col];
                lower = lp.col_lower_[i_col];
                upper = lp.col_upper_[i_col];
                value = solution.col_value[i_col];
                if have_dual_solution {
                    dual = solution.col_dual[i_col];
                }
                if have_integrality {
                    integrality = lp.integrality_[i_col];
                }
                if pass == 0 {
                    if dual * dual < dual_feasibility_tolerance {
                        // HiGHS 1.16.0 (HighsSolution.cpp lines 283-284) takes the
                        // effective cost here; 1.15.1 uses the raw cost, which
                        // `effective_cost` is when `h16_kkt_effective_costs` is off.
                        highs_norm_costs = cpp_max(effective_cost.abs(), highs_norm_costs);
                    }
                    if get_residuals && have_dual_solution {
                        let mut q_dual_activity = HighsCDouble::from(dual_activity[i_col]);
                        q_dual_activity -= gradient[i_col];
                        dual_activity[i_col] = f64::from(q_dual_activity);
                    }
                }
            } else {
                let i_row = (i_var - lp.num_col_) as usize;
                lower = lp.row_lower_[i_row];
                upper = lp.row_upper_[i_row];
                value = solution.row_value[i_row];
                if have_dual_solution {
                    dual = solution.row_dual[i_row];
                }
                integrality = HighsVarType::kContinuous;
            }
            // Flip dual according to lp.sense_
            dual *= (lp.sense_ as HighsInt) as f64;
            get_variable_kkt_failures(
                primal_feasibility_tolerance,
                dual_feasibility_tolerance,
                mip_feasibility_tolerance,
                lower,
                upper,
                value,
                dual,
                integrality,
                &mut primal_infeasibility,
                &mut dual_infeasibility,
                &mut semi_infeasibility,
                &mut at_status,
                &mut mid_status,
            );
            if pass == 0 {
                if at_status == kHighsSolutionLo {
                    highs_norm_bounds = cpp_max(lower.abs(), highs_norm_bounds);
                } else if at_status == kHighsSolutionUp {
                    highs_norm_bounds = cpp_max(upper.abs(), highs_norm_bounds);
                }
            } else {
                if primal_infeasibility > 0.0 {
                    if primal_infeasibility > primal_feasibility_tolerance {
                        highs_info.num_primal_infeasibilities += 1;
                    }
                    if highs_info.max_primal_infeasibility < primal_infeasibility {
                        highs_info.max_primal_infeasibility = primal_infeasibility;
                    }
                    highs_info.sum_primal_infeasibilities += primal_infeasibility;
                    let mut relative_bound_measure = highs_norm_bounds;
                    if at_status == kHighsSolutionNo {
                        if mid_status == kHighsSolutionNo || mid_status == kHighsSolutionLo {
                            relative_bound_measure = cpp_max(lower.abs(), relative_bound_measure);
                        } else {
                            debug_assert!(mid_status == kHighsSolutionUp);
                            relative_bound_measure = cpp_max(upper.abs(), relative_bound_measure);
                        }
                    }
                    let relative_primal_infeasibility =
                        primal_infeasibility / (1.0 + relative_bound_measure);

                    if relative_primal_infeasibility > primal_feasibility_tolerance {
                        highs_info.num_relative_primal_infeasibilities += 1;
                    }
                    if highs_info.max_relative_primal_infeasibility < relative_primal_infeasibility
                    {
                        highs_info.max_relative_primal_infeasibility =
                            relative_primal_infeasibility;
                    }
                }
                if semi_infeasibility > 0.0 {
                    highs_info.num_semi_infeasibilities += 1;
                    if highs_info.max_semi_infeasibility < semi_infeasibility {
                        highs_info.max_semi_infeasibility = semi_infeasibility;
                    }
                    highs_info.sum_semi_infeasibilities += semi_infeasibility;
                }
                if have_dual_solution && dual_infeasibility > 0.0 {
                    if dual_infeasibility > dual_feasibility_tolerance {
                        highs_info.num_dual_infeasibilities += 1;
                    }
                    if highs_info.max_dual_infeasibility < dual_infeasibility {
                        highs_info.max_dual_infeasibility = dual_infeasibility;
                    }
                    highs_info.sum_dual_infeasibilities += dual_infeasibility;

                    let mut relative_cost_measure = highs_norm_costs;
                    if is_col && cost != 0.0 && dual * dual >= dual_feasibility_tolerance {
                        // HiGHS 1.16.0 (HighsSolution.cpp line 397) uses the effective
                        // cost here, the test staying on the raw cost as in 1.15.1
                        relative_cost_measure =
                            cpp_max(effective_cost.abs(), relative_cost_measure);
                    }
                    let relative_dual_infeasibility =
                        dual_infeasibility / (1.0 + relative_cost_measure);

                    if relative_dual_infeasibility > dual_feasibility_tolerance {
                        highs_info.num_relative_dual_infeasibilities += 1;
                    }
                    if highs_info.max_relative_dual_infeasibility < relative_dual_infeasibility {
                        highs_info.max_relative_dual_infeasibility = relative_dual_infeasibility;
                    }
                }
                if !is_col && get_residuals {
                    let i_row = (i_var - lp.num_col_) as usize;
                    let primal_residual_error =
                        (primal_activity[i_row] - solution.row_value[i_row]).abs();
                    let relative_primal_residual_error =
                        primal_residual_error / (1.0 + highs_norm_bounds);

                    if primal_residual_error > primal_residual_tolerance {
                        highs_info.num_primal_residual_errors += 1;
                    }
                    if highs_info.max_primal_residual_error < primal_residual_error {
                        highs_info.max_primal_residual_error = primal_residual_error;
                    }

                    if relative_primal_residual_error > primal_residual_tolerance {
                        highs_info.num_relative_primal_residual_errors += 1;
                    }
                    if highs_info.max_relative_primal_residual_error
                        < relative_primal_residual_error
                    {
                        highs_info.max_relative_primal_residual_error =
                            relative_primal_residual_error;
                    }
                }
                if is_col && get_residuals && have_dual_solution {
                    let i_col = i_var as usize;
                    let dual_residual_error =
                        (dual_activity[i_col] + solution.col_dual[i_col]).abs();
                    let relative_dual_residual_error =
                        dual_residual_error / (1.0 + highs_norm_costs);

                    if dual_residual_error > dual_residual_tolerance {
                        highs_info.num_dual_residual_errors += 1;
                    }
                    if highs_info.max_dual_residual_error < dual_residual_error {
                        highs_info.max_dual_residual_error = dual_residual_error;
                    }

                    if relative_dual_residual_error > dual_residual_tolerance {
                        highs_info.num_relative_dual_residual_errors += 1;
                    }
                    if highs_info.max_relative_dual_residual_error < relative_dual_residual_error {
                        highs_info.max_relative_dual_residual_error = relative_dual_residual_error;
                    }
                }
            }
            if pass == 1 && i_var == lp.num_col_ - 1 {
                max_col_primal_infeasibility = highs_info.max_primal_infeasibility;
                max_col_dual_infeasibility = highs_info.max_dual_infeasibility;

                max_relative_col_primal_infeasibility =
                    highs_info.max_relative_primal_infeasibility;
                max_relative_col_dual_infeasibility = highs_info.max_relative_dual_infeasibility;

                highs_info.max_primal_infeasibility = 0.0;
                highs_info.max_dual_infeasibility = 0.0;

                highs_info.max_relative_primal_infeasibility = 0.0;
                highs_info.max_relative_dual_infeasibility = 0.0;
            }
        }
    }

    let max_row_primal_infeasibility = highs_info.max_primal_infeasibility;
    let max_row_dual_infeasibility = highs_info.max_dual_infeasibility;

    let max_relative_row_primal_infeasibility = highs_info.max_relative_primal_infeasibility;
    let max_relative_row_dual_infeasibility = highs_info.max_relative_dual_infeasibility;

    highs_info.max_primal_infeasibility =
        cpp_max(max_col_primal_infeasibility, max_row_primal_infeasibility);
    highs_info.max_dual_infeasibility =
        cpp_max(max_col_dual_infeasibility, max_row_dual_infeasibility);

    highs_info.max_relative_primal_infeasibility = cpp_max(
        max_relative_col_primal_infeasibility,
        max_relative_row_primal_infeasibility,
    );
    highs_info.max_relative_dual_infeasibility = cpp_max(
        max_relative_col_dual_infeasibility,
        max_relative_row_dual_infeasibility,
    );

    if have_dual_solution {
        let have_values = get_complementarity_violations(
            lp,
            solution,
            optimality_tolerance,
            &mut highs_info.num_complementarity_violations,
            &mut highs_info.max_complementarity_violation,
        );
        debug_assert!(have_values);
    }

    if have_dual_solution {
        let gradient_p: Option<&[f64]> = if is_qp { Some(gradient) } else { None };
        let mut dual_objective_value: f64 = 0.0;
        let dual_objective_status = compute_dual_objective_value(
            gradient_p,
            lp,
            solution,
            &mut dual_objective_value,
            options.h16_obj_cdouble,
        );
        debug_assert!(dual_objective_status);
        let abs_objective_difference =
            (highs_info.objective_function_value - dual_objective_value).abs();
        highs_info.primal_dual_objective_error = abs_objective_difference;
        let denominator =
            1.0 + highs_info.objective_function_value.abs() + dual_objective_value.abs();
        highs_info.primal_dual_objective_error /= denominator;
    }

    if printf_kkt || options.log_dev_level > 0 {
        let g = |x: f64| c_fmt_g_width(x, 3, 8);
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "getKktFailures:: cost norm = {}; bound norm = {}\n",
            g(highs_norm_costs),
            g(highs_norm_bounds)
        );
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "getKktFailures:                      LP  (abs / rel)         Col (abs / rel)         Row (abs / rel)\n"
        );
        crate::highs_log_dev!(
            &options.log_options,
            HighsLogType::kInfo,
            "getKktFailures: primal infeasibility {} / {}     {} / {}     {} / {}\n",
            g(highs_info.max_primal_infeasibility),
            g(highs_info.max_relative_primal_infeasibility),
            g(max_col_primal_infeasibility),
            g(max_relative_col_primal_infeasibility),
            g(max_row_primal_infeasibility),
            g(max_relative_row_primal_infeasibility)
        );
        if have_dual_solution {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "getKktFailures:   dual infeasibility {} / {}     {} / {}     {} / {}\n",
                g(highs_info.max_dual_infeasibility),
                g(highs_info.max_relative_dual_infeasibility),
                g(max_col_dual_infeasibility),
                g(max_relative_col_dual_infeasibility),
                g(max_row_dual_infeasibility),
                g(max_relative_row_dual_infeasibility)
            );
        }
        if get_residuals {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "getKktFailures: primal residual      {} / {}\n",
                g(highs_info.max_primal_residual_error),
                g(highs_info.max_relative_primal_residual_error)
            );
            if have_dual_solution {
                crate::highs_log_dev!(
                    &options.log_options,
                    HighsLogType::kInfo,
                    "getKktFailures:   dual residual      {} / {}\n",
                    g(highs_info.max_dual_residual_error),
                    g(highs_info.max_relative_dual_residual_error)
                );
            }
        }
        if !is_qp && have_dual_solution {
            crate::highs_log_dev!(
                &options.log_options,
                HighsLogType::kInfo,
                "getKktFailures: objective gap        {}\n",
                g(highs_info.primal_dual_objective_error)
            );
        }
    }

    if highs_info.num_primal_infeasibilities + highs_info.num_semi_infeasibilities != 0 {
        highs_info.primal_solution_status = kSolutionStatusInfeasible;
    } else {
        highs_info.primal_solution_status = kSolutionStatusFeasible;
    }
    if have_dual_solution {
        if highs_info.num_dual_infeasibilities != 0 {
            highs_info.dual_solution_status = kSolutionStatusInfeasible;
        } else {
            highs_info.dual_solution_status = kSolutionStatusFeasible;
        }
    }
}

/// HiGHS: HighsSolution.cpp::getVariableKktFailures (lines 568-670)
/// (`index` is an unused defaulted trailing parameter in C++.)
#[allow(clippy::too_many_arguments)]
pub fn get_variable_kkt_failures(
    primal_feasibility_tolerance: f64,
    dual_feasibility_tolerance: f64,
    mip_feasibility_tolerance: f64,
    lower: f64,
    upper: f64,
    value: f64,
    dual: f64,
    integrality: HighsVarType,
    primal_infeasibility: &mut f64,
    dual_infeasibility: &mut f64,
    semi_infeasibility: &mut f64,
    at_status: &mut u8,
    mid_status: &mut u8,
) {
    let infeasibility_residual = infeasibility(lower, value, upper, primal_feasibility_tolerance);
    *primal_infeasibility = infeasibility_residual.1;
    *semi_infeasibility = 0.0;
    // Determine whether this value is close to a bound
    *at_status = kHighsSolutionNo;
    let mut bound_residual = (lower - value).abs();
    if bound_residual * bound_residual <= primal_feasibility_tolerance {
        *at_status = kHighsSolutionLo;
    } else {
        bound_residual = (value - upper).abs();
        if bound_residual * bound_residual <= primal_feasibility_tolerance {
            *at_status = kHighsSolutionUp;
        }
    }
    *mid_status = kHighsSolutionNo;
    if lower < upper {
        let length = upper - lower;
        if lower <= -kHighsInf && upper >= kHighsInf {
            // Free variable
            *dual_infeasibility = dual.abs();
        } else if length * length > primal_feasibility_tolerance {
            let middle = (lower + upper) * 0.5;
            if value < middle {
                *mid_status = kHighsSolutionLo;
                *dual_infeasibility = cpp_max(-dual, 0.);
            } else {
                *mid_status = kHighsSolutionUp;
                *dual_infeasibility = cpp_max(dual, 0.);
            }
        } else {
            *dual_infeasibility = 0.0;
        }
    } else {
        // Fixed variable
        *dual_infeasibility = 0.0;
    }
    // Account for semi-variables
    let semi_variable =
        integrality == HighsVarType::kSemiContinuous || integrality == HighsVarType::kSemiInteger;
    let mid = cpp_max(primal_feasibility_tolerance, lower * 0.5);

    if semi_variable {
        *primal_infeasibility = 0.0;
        *semi_infeasibility = 0.0;
        if value < mid {
            *semi_infeasibility = value.abs();
            if *semi_infeasibility < mip_feasibility_tolerance {
                *semi_infeasibility = 0.0;
            }
        } else {
            if value < lower {
                *semi_infeasibility = lower - value;
            } else if value > upper {
                *semi_infeasibility = value - upper;
            }
            debug_assert!(*semi_infeasibility >= 0.0);
            if *semi_infeasibility < primal_feasibility_tolerance {
                *semi_infeasibility = 0.0;
            }
        }
    }
}

/// HiGHS: HighsSolution.cpp::getComplementarityViolations (lines 996-1033)
pub fn get_complementarity_violations(
    lp: &HighsLp,
    solution: &HighsSolution,
    optimality_tolerance: f64,
    num_complementarity_violation: &mut HighsInt,
    max_complementarity_violation: &mut f64,
) -> bool {
    *num_complementarity_violation = kHighsIllegalComplementarityCount;
    *max_complementarity_violation = kHighsIllegalComplementarityViolation;
    if !solution.dual_valid {
        return false;
    }

    *num_complementarity_violation = 0;
    *max_complementarity_violation = 0.0;
    let mut primal_residual: f64;
    for i_var in 0..lp.num_col_ + lp.num_row_ {
        let is_col = i_var < lp.num_col_;
        let i_row = i_var - lp.num_col_;
        let (primal, dual, lower, upper) = if is_col {
            let ic = i_var as usize;
            (
                solution.col_value[ic],
                solution.col_dual[ic],
                lp.col_lower_[ic],
                lp.col_upper_[ic],
            )
        } else {
            let ir = i_row as usize;
            (
                solution.row_value[ir],
                solution.row_dual[ir],
                lp.row_lower_[ir],
                lp.row_upper_[ir],
            )
        };
        if lower <= -kHighsInf && upper >= kHighsInf {
            // Free
            primal_residual = 1.0;
        } else {
            let mid = (lower + upper) * 0.5;
            primal_residual = if primal < mid {
                (lower - primal).abs()
            } else {
                (upper - primal).abs()
            };
        }
        let dual_residual = dual.abs();
        let complementarity_violation = primal_residual * dual_residual;
        if complementarity_violation > optimality_tolerance {
            *num_complementarity_violation += 1;
        }
        *max_complementarity_violation =
            cpp_max(complementarity_violation, *max_complementarity_violation);
    }
    true
}

/// HiGHS: HighsSolution.cpp::computeDualObjectiveValue(gradient, lp, solution,
/// dual_objective_value) (lines 1345-1389); `gradient` is the C++ nullable pointer.
/// HiGHS 1.16 (#3513): HighsSolution.cpp::computeDualObjectiveValue (lines 1468-1514)
pub fn compute_dual_objective_value(
    gradient: Option<&[f64]>,
    lp: &HighsLp,
    solution: &HighsSolution,
    dual_objective_value: &mut f64,
    h16_obj_cdouble: bool,
) -> bool {
    *dual_objective_value = 0.0;
    if !solution.dual_valid {
        return false;
    }
    debug_assert!(solution.col_value.len() == lp.num_col_ as usize);
    debug_assert!(solution.col_dual.len() == lp.num_col_ as usize);
    debug_assert!(solution.row_value.len() == lp.num_row_ as usize);
    debug_assert!(solution.row_dual.len() == lp.num_row_ as usize);

    // HiGHS 1.16 #3513 (h16_obj_cdouble): HighsSolution.cpp lines 1479-1512 (1.16.0)
    if h16_obj_cdouble {
        return compute_dual_objective_value_h16(gradient, lp, solution, dual_objective_value);
    }
    *dual_objective_value = lp.offset_;
    if let Some(gradient) = gradient {
        let mut quad_value: f64 = 0.0;
        for i_col in 0..lp.num_col_ as usize {
            quad_value += (lp.col_cost_[i_col] - gradient[i_col]) * solution.col_value[i_col];
        }
        *dual_objective_value += 0.5 * quad_value;
    }
    let mut bound: f64;
    for i_var in 0..lp.num_col_ + lp.num_row_ {
        let is_col = i_var < lp.num_col_;
        let i_row = i_var - lp.num_col_;
        let (primal, dual, lower, upper) = if is_col {
            let ic = i_var as usize;
            (
                solution.col_value[ic],
                solution.col_dual[ic],
                lp.col_lower_[ic],
                lp.col_upper_[ic],
            )
        } else {
            let ir = i_row as usize;
            (
                solution.row_value[ir],
                solution.row_dual[ir],
                lp.row_lower_[ir],
                lp.row_upper_[ir],
            )
        };
        if lower <= -kHighsInf && upper >= kHighsInf {
            // Free
            bound = 1.0;
        } else {
            let mid = (lower + upper) * 0.5;
            bound = if primal < mid { lower } else { upper };
        }
        *dual_objective_value += bound * dual;
    }
    true
}

/// HiGHS 1.16 (#3513): HighsSolution.cpp::computeDualObjectiveValue (lines 1479-1512),
/// the HighsCDouble accumulation. Called only from
/// [`compute_dual_objective_value`] when `h16_obj_cdouble` is on, after the
/// shared `dual_valid` check (lines 1471-1477).
fn compute_dual_objective_value_h16(
    gradient: Option<&[f64]>,
    lp: &HighsLp,
    solution: &HighsSolution,
    dual_objective_value: &mut f64,
) -> bool {
    let mut extended_dual_objective_value = HighsCDouble::from(lp.offset_);
    if let Some(gradient) = gradient {
        let mut extended_quad_value = HighsCDouble::from(0.0);
        for i_col in 0..lp.num_col_ as usize {
            let extended_coefficient =
                HighsCDouble::from(lp.col_cost_[i_col]) - HighsCDouble::from(gradient[i_col]);
            extended_quad_value += extended_coefficient * solution.col_value[i_col];
        }
        extended_dual_objective_value += 0.5 * extended_quad_value;
    }
    let mut bound: f64;
    for i_var in 0..lp.num_col_ + lp.num_row_ {
        let is_col = i_var < lp.num_col_;
        let i_row = i_var - lp.num_col_;
        let (primal, dual, lower, upper) = if is_col {
            let ic = i_var as usize;
            (
                solution.col_value[ic],
                solution.col_dual[ic],
                lp.col_lower_[ic],
                lp.col_upper_[ic],
            )
        } else {
            let ir = i_row as usize;
            (
                solution.row_value[ir],
                solution.row_dual[ir],
                lp.row_lower_[ir],
                lp.row_upper_[ir],
            )
        };
        if lower <= -kHighsInf && upper >= kHighsInf {
            // Free
            bound = 1.0;
        } else {
            let mid = (lower + upper) * 0.5;
            bound = if primal < mid { lower } else { upper };
        }
        extended_dual_objective_value += HighsCDouble::from(bound) * dual;
    }
    *dual_objective_value = f64::from(extended_dual_objective_value);
    true
}

/// HiGHS: HighsSolution.cpp::getLpKktFailures(options, lp, solution, basis,
/// highs_info, primal_dual_errors, get_residuals) (lines 53-64).
///
/// `getPrimalDualBasisErrors` / `getPrimalDualGlpsolErrors` write only the
/// caller's `primal_dual_errors`, which `lpKktCheck` discards; not ported.
pub fn get_lp_kkt_failures(
    options: &HighsOptions,
    lp: &HighsLp,
    solution: &HighsSolution,
    highs_info: &mut HighsInfo,
    get_residuals: bool,
) {
    get_kkt_failures(
        options,
        false,
        lp,
        &lp.col_cost_,
        solution,
        highs_info,
        get_residuals,
    );
}

/// C `%4.0e`.
fn fmt_tol(v: f64) -> String {
    crate::util::c_printf::c_pad_right(&crate::util::c_printf::c_fmt_e(v, 0), 4)
}

/// C `%8.3g`.
fn fmt_g83(v: f64) -> String {
    c_fmt_g_width(v, 3, 8)
}

/// HiGHS: HighsSolution.cpp::lpKktCheck (lines 1043-1327).
pub fn lp_kkt_check(
    model_status: &mut HighsModelStatus,
    info: &mut HighsInfo,
    lp: &HighsLp,
    solution: &HighsSolution,
    basis: &HighsBasis,
    options: &HighsOptions,
    message: &str,
) {
    if !solution.value_valid {
        return;
    }
    let log_options = &options.log_options;
    let mut primal_feasibility_tolerance = options.primal_feasibility_tolerance;
    let mut dual_feasibility_tolerance = options.dual_feasibility_tolerance;
    let mut primal_residual_tolerance = options.primal_residual_tolerance;
    let mut dual_residual_tolerance = options.dual_residual_tolerance;
    let mut optimality_tolerance = options.optimality_tolerance;
    if options.kkt_tolerance != kDefaultKktTolerance {
        primal_feasibility_tolerance = options.kkt_tolerance;
        dual_feasibility_tolerance = options.kkt_tolerance;
        primal_residual_tolerance = options.kkt_tolerance;
        dual_residual_tolerance = options.kkt_tolerance;
        optimality_tolerance = options.kkt_tolerance;
    }
    info.objective_function_value =
        lp.objective_value(&solution.col_value, options.h16_obj_cdouble);
    let get_residuals = !basis.valid;
    get_lp_kkt_failures(options, lp, solution, info, get_residuals);
    if *model_status == HighsModelStatus::kOptimal {
        report_kkt_failures(lp, options, info, message);
    }
    // get_residuals is false when there is a valid basis, since
    // residual errors are assumed to be small, so
    // info.num_primal_residual_errors = -1, since they aren't
    // known. Hence don't consider this in identifying unboundedness
    // from HighsModelStatus::kUnboundedOrInfeasible
    if *model_status == HighsModelStatus::kUnboundedOrInfeasible
        && info.num_primal_infeasibilities == 0
        && (!get_residuals || info.num_primal_residual_errors == 0)
    {
        *model_status = HighsModelStatus::kUnbounded;
    }
    let was_optimal = *model_status == HighsModelStatus::kOptimal;
    let mut kkt_ok = true;
    let mut written_optimality_error_header = false;

    let mut found_optimality_error = || {
        kkt_ok = false;
        if !was_optimal || written_optimality_error_header {
            return;
        }
        crate::highs_log_user!(
            log_options,
            HighsLogType::kWarning,
            "LP solver claims optimality, but with\n"
        );
        written_optimality_error_header = true;
    };

    let mut max_primal_tolerance_relative_violation: f64 = 0.0;
    let mut max_dual_tolerance_relative_violation: f64 = 0.0;
    let mut primal_dual_objective_tolerance_relative_violation: f64 = 0.0;
    let max_allowed_tolerance_relative_violation: f64 = 1e2;
    if basis.valid {
        if info.num_primal_infeasibilities > 0 {
            max_primal_tolerance_relative_violation = cpp_max(
                info.max_primal_infeasibility / primal_feasibility_tolerance,
                max_primal_tolerance_relative_violation,
            );
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max/sum {:6} / {} / {} primal infeasibilities       (tolerance = {})\n",
                    info.num_primal_infeasibilities,
                    fmt_g83(info.max_primal_infeasibility),
                    fmt_g83(info.sum_primal_infeasibilities),
                    fmt_tol(primal_feasibility_tolerance)
                );
            }
        }
        if info.num_dual_infeasibilities > 0 {
            max_dual_tolerance_relative_violation = cpp_max(
                info.max_dual_infeasibility / dual_feasibility_tolerance,
                max_dual_tolerance_relative_violation,
            );
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max/sum {:6} / {} / {}   dual infeasibilities       (tolerance = {})\n",
                    info.num_dual_infeasibilities,
                    fmt_g83(info.max_dual_infeasibility),
                    fmt_g83(info.sum_dual_infeasibilities),
                    fmt_tol(dual_feasibility_tolerance)
                );
            }
        }
        // An optimal basic solution has no complementarity violations
        // by construction, and can be assumed to have no relative
        // primal or dual residual errors or meaningful primal dual
        // objective error
        let mut unexpected_error_if_optimal = info.num_complementarity_violations != 0;
        let mut local_dual_objective: f64 = 0.0;
        if info.primal_dual_objective_error > optimality_tolerance {
            // Ignore primal-dual objective errors if both objectives are small
            let ok_dual_objective = compute_dual_objective_value(
                None,
                lp,
                solution,
                &mut local_dual_objective,
                options.h16_obj_cdouble,
            );
            debug_assert!(ok_dual_objective);
            if info.objective_function_value * info.objective_function_value > optimality_tolerance
                && local_dual_objective * local_dual_objective > optimality_tolerance
            {
                unexpected_error_if_optimal = true;
            }
        }
        let have_residual_errors = info.num_primal_residual_errors != kHighsIllegalResidualCount;
        if have_residual_errors {
            unexpected_error_if_optimal = unexpected_error_if_optimal
                || info.num_relative_primal_residual_errors != 0
                || info.num_relative_dual_residual_errors != 0;
            max_primal_tolerance_relative_violation = cpp_max(
                info.max_relative_primal_residual_error / primal_residual_tolerance,
                max_primal_tolerance_relative_violation,
            );
            max_dual_tolerance_relative_violation = cpp_max(
                info.max_relative_dual_residual_error / dual_residual_tolerance,
                max_dual_tolerance_relative_violation,
            );
        }
        primal_dual_objective_tolerance_relative_violation =
            info.primal_dual_objective_error / optimality_tolerance;

        if was_optimal && unexpected_error_if_optimal {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kWarning,
                "Optimal basic solution has {} complementarity violations and {} primal dual objective error from primal (dual) objective = {} ({})\n",
                info.num_complementarity_violations,
                crate::util::c_printf::c_fmt_g(info.primal_dual_objective_error, 6),
                crate::util::c_printf::c_fmt_g(info.objective_function_value, 6),
                crate::util::c_printf::c_fmt_g(local_dual_objective, 6)
            );
            if have_residual_errors {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max {:6} / {}  relative primal residual errors         (tolerance = {})\n",
                    info.num_relative_primal_residual_errors,
                    fmt_g83(info.max_relative_primal_residual_error),
                    fmt_tol(primal_residual_tolerance)
                );
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max {:6} / {}  relative   dual residual errors         (tolerance = {})\n",
                    info.num_relative_dual_residual_errors,
                    fmt_g83(info.max_relative_dual_residual_error),
                    fmt_tol(dual_residual_tolerance)
                );
            }
            // C++ asserts (compiled out in the NDEBUG reference build).
        }
        // Infeasibility of the primal and dual solutions based on number
        // of primal/dual infeasibilities should have been set in
        // getKktFailures, but qualify this if the residuals are
        // meaningful
        if info.num_primal_infeasibilities != 0 {
            debug_assert!(info.primal_solution_status == kSolutionStatusInfeasible);
        } else {
            info.primal_solution_status = kSolutionStatusFeasible;
        }
        if info.num_dual_infeasibilities != 0 {
            debug_assert!(info.dual_solution_status == kSolutionStatusInfeasible);
        } else {
            info.dual_solution_status = kSolutionStatusFeasible;
        }
        // Overrule feasibility if large relative tolerance failures have
        // occurred - pretty inconceivable since absolute residuals should
        // be small with a basis
        if max_primal_tolerance_relative_violation > max_allowed_tolerance_relative_violation {
            info.primal_solution_status = kSolutionStatusInfeasible;
        }
        if max_dual_tolerance_relative_violation > max_allowed_tolerance_relative_violation {
            info.dual_solution_status = kSolutionStatusInfeasible;
        }
    } else {
        // A solution without a basis may have primal or dual residual
        // errors, and complementarity errors - due to the convergence
        // being based on relative primal-dual objective error, so test
        // the latter
        let mut tolerance_relative_violation =
            info.max_relative_primal_infeasibility / primal_feasibility_tolerance;
        max_primal_tolerance_relative_violation = cpp_max(
            tolerance_relative_violation,
            max_primal_tolerance_relative_violation,
        );
        if info.num_relative_primal_infeasibilities > 0 {
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max {:6} / {} relative primal infeasibilities (tolerance = {})\n",
                    info.num_relative_primal_infeasibilities,
                    fmt_g83(info.max_relative_primal_infeasibility),
                    fmt_tol(primal_feasibility_tolerance)
                );
            }
        }
        tolerance_relative_violation =
            info.max_relative_dual_infeasibility / dual_feasibility_tolerance;
        max_dual_tolerance_relative_violation = cpp_max(
            tolerance_relative_violation,
            max_dual_tolerance_relative_violation,
        );
        if info.num_relative_dual_infeasibilities > 0 {
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max {:6} / {} relative   dual infeasibilities (tolerance = {})\n",
                    info.num_relative_dual_infeasibilities,
                    fmt_g83(info.max_relative_dual_infeasibility),
                    fmt_tol(dual_feasibility_tolerance)
                );
            }
        }
        tolerance_relative_violation =
            info.max_relative_primal_residual_error / primal_residual_tolerance;
        max_primal_tolerance_relative_violation = cpp_max(
            tolerance_relative_violation,
            max_primal_tolerance_relative_violation,
        );
        if info.num_relative_primal_residual_errors > 0 {
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max {:6} / {} relative primal residual errors (tolerance = {})\n",
                    info.num_relative_primal_residual_errors,
                    fmt_g83(info.max_relative_primal_residual_error),
                    fmt_tol(primal_residual_tolerance)
                );
            }
        }
        tolerance_relative_violation =
            info.max_relative_dual_residual_error / dual_residual_tolerance;
        max_dual_tolerance_relative_violation = cpp_max(
            tolerance_relative_violation,
            max_dual_tolerance_relative_violation,
        );
        if info.num_relative_dual_residual_errors > 0 {
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "   num/max {:6} / {} relative   dual residual errors (tolerance = {})\n",
                    info.num_relative_dual_residual_errors,
                    fmt_g83(info.max_relative_dual_residual_error),
                    fmt_tol(dual_residual_tolerance)
                );
            }
        }
        if info.primal_dual_objective_error > optimality_tolerance {
            primal_dual_objective_tolerance_relative_violation =
                info.primal_dual_objective_error / optimality_tolerance;
            found_optimality_error();
            if was_optimal {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kWarning,
                    "                    {} relative P-D objective error    (tolerance = {})\n",
                    fmt_g83(info.primal_dual_objective_error),
                    fmt_tol(optimality_tolerance)
                );
            }
        }
        // Set the primal and dual solution status according to tolerance failure
        if max_primal_tolerance_relative_violation > max_allowed_tolerance_relative_violation {
            info.primal_solution_status = kSolutionStatusInfeasible;
        } else {
            info.primal_solution_status = kSolutionStatusFeasible;
        }
        if max_dual_tolerance_relative_violation > max_allowed_tolerance_relative_violation {
            info.dual_solution_status = kSolutionStatusInfeasible;
        } else {
            info.dual_solution_status = kSolutionStatusFeasible;
        }
    }
    let _ = kkt_ok;
    let mut max_tolerance_relative_violation = primal_dual_objective_tolerance_relative_violation;
    max_tolerance_relative_violation = cpp_max(
        max_primal_tolerance_relative_violation,
        max_tolerance_relative_violation,
    );
    max_tolerance_relative_violation = cpp_max(
        max_dual_tolerance_relative_violation,
        max_tolerance_relative_violation,
    );
    //
    // Now see whether optimality is compromised or permitted given the tolerance
    // failures
    if *model_status == HighsModelStatus::kOptimal {
        if max_tolerance_relative_violation > max_allowed_tolerance_relative_violation {
            *model_status = HighsModelStatus::kUnknown;
            crate::highs_log_user!(
                log_options,
                HighsLogType::kWarning,
                "Model status changed from \"Optimal\" to \"Unknown\" since relative violation of tolerances is {}\n",
                fmt_g83(max_tolerance_relative_violation)
            );
        } else if max_allowed_tolerance_relative_violation > 1.0
            && max_tolerance_relative_violation > 1.0
        {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                "Model status is \"Optimal\" since relative violation of tolerances is no more than {}\n",
                fmt_g83(max_tolerance_relative_violation)
            );
        }
    } else if *model_status == HighsModelStatus::kUnknown
        && max_tolerance_relative_violation <= max_allowed_tolerance_relative_violation
    {
        *model_status = HighsModelStatus::kOptimal;
        crate::highs_log_user!(
            log_options,
            HighsLogType::kWarning,
            "Model status changed from \"Unknown\" to \"Optimal\"\n"
        );
    }
    crate::highs_log_user!(log_options, HighsLogType::kInfo, "\n");
}

/// HiGHS: HighsSolution.cpp::reportKktFailures (lines 2045-2133). The
/// `printf_kkt` block is dead (`printf_kkt = false`).
pub fn report_kkt_failures(
    lp: &HighsLp,
    options: &HighsOptions,
    info: &HighsInfo,
    message: &str,
) -> bool {
    let log_options = &options.log_options;
    let mut mip_feasibility_tolerance = options.mip_feasibility_tolerance;
    let mut primal_feasibility_tolerance = options.primal_feasibility_tolerance;
    let mut dual_feasibility_tolerance = options.dual_feasibility_tolerance;
    let mut primal_residual_tolerance = options.primal_residual_tolerance;
    let mut dual_residual_tolerance = options.dual_residual_tolerance;
    let mut optimality_tolerance = options.optimality_tolerance;
    let is_mip = lp.is_mip();
    if is_mip {
        primal_feasibility_tolerance = mip_feasibility_tolerance;
    } else if options.kkt_tolerance != kDefaultKktTolerance {
        mip_feasibility_tolerance = options.kkt_tolerance;
        primal_feasibility_tolerance = options.kkt_tolerance;
        dual_feasibility_tolerance = options.kkt_tolerance;
        primal_residual_tolerance = options.kkt_tolerance;
        dual_residual_tolerance = options.kkt_tolerance;
        optimality_tolerance = options.kkt_tolerance;
    }

    let force_report = options.log_dev_level >= crate::io::highs_io::kHighsLogDevLevelInfo;
    let complementarity_error = !is_mip && info.primal_dual_objective_error > optimality_tolerance;
    let integrality_error = is_mip && info.max_integrality_violation >= mip_feasibility_tolerance;
    let has_kkt_failures = integrality_error
        || info.num_primal_infeasibilities > 0
        || info.num_dual_infeasibilities > 0
        || info.num_primal_residual_errors > 0
        || info.num_dual_residual_errors > 0
        || complementarity_error;
    if !has_kkt_failures && !force_report {
        return has_kkt_failures;
    }

    let log_type = if has_kkt_failures {
        HighsLogType::kWarning
    } else {
        HighsLogType::kInfo
    };

    crate::highs_log_user!(
        log_options,
        log_type,
        "Solution optimality conditions{}{}\n",
        if message.is_empty() { "" } else { ": " },
        message
    );
    if is_mip && info.max_integrality_violation >= 0.0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "    max      {}                                  integrality violations     (tolerance = {})\n",
            fmt_g83(info.max_integrality_violation),
            fmt_tol(mip_feasibility_tolerance)
        );
    }
    if info.num_primal_infeasibilities >= 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "num/max {:6} / {} (relative {:6} / {}) primal infeasibilities     (tolerance = {})\n",
            info.num_primal_infeasibilities,
            fmt_g83(info.max_primal_infeasibility),
            info.num_relative_primal_infeasibilities,
            fmt_g83(info.max_relative_primal_infeasibility),
            fmt_tol(primal_feasibility_tolerance)
        );
    }
    if info.num_dual_infeasibilities >= 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "num/max {:6} / {} (relative {:6} / {})   dual infeasibilities     (tolerance = {})\n",
            info.num_dual_infeasibilities,
            fmt_g83(info.max_dual_infeasibility),
            info.num_relative_dual_infeasibilities,
            fmt_g83(info.max_relative_dual_infeasibility),
            fmt_tol(dual_feasibility_tolerance)
        );
    }
    if info.num_primal_residual_errors >= 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "num/max {:6} / {} (relative {:6} / {}) primal residual errors     (tolerance = {})\n",
            info.num_primal_residual_errors,
            fmt_g83(info.max_primal_residual_error),
            info.num_relative_primal_residual_errors,
            fmt_g83(info.max_relative_primal_residual_error),
            fmt_tol(primal_residual_tolerance)
        );
    }
    if info.num_dual_residual_errors >= 0 {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "num/max {:6} / {} (relative {:6} / {})   dual residual errors     (tolerance = {})\n",
            info.num_dual_residual_errors,
            fmt_g83(info.max_dual_residual_error),
            info.num_relative_dual_residual_errors,
            fmt_g83(info.max_relative_dual_residual_error),
            fmt_tol(dual_residual_tolerance)
        );
    }
    if info.primal_dual_objective_error != kHighsIllegalComplementarityViolation {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "                                         {} / {}  P-D objective error        (tolerance = {})\n",
            if info.primal_dual_objective_error > optimality_tolerance {
                1
            } else {
                0
            },
            fmt_g83(info.primal_dual_objective_error),
            fmt_tol(optimality_tolerance)
        );
    }
    has_kkt_failures
}

#[cfg(test)]
mod h16_obj_cdouble_tests {
    use super::*;

    /// offset 1e16, two columns fixed at 1 with duals [1, 1], no rows: the
    /// bound * dual terms are each 1, which a double sum loses against 1e16.
    fn cancel_case() -> (HighsLp, HighsSolution) {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![1.0, 1.0];
        lp.col_upper_ = vec![1.0, 1.0];
        lp.offset_ = 1e16;
        let solution = HighsSolution {
            value_valid: true,
            dual_valid: true,
            col_value: vec![1.0, 1.0],
            col_dual: vec![1.0, 1.0],
            row_value: vec![],
            row_dual: vec![],
        };
        (lp, solution)
    }

    #[test]
    fn h16_obj_cdouble_dual_objective_on() {
        // HighsSolution.cpp lines 1479-1512 (1.16.0)
        let (lp, solution) = cancel_case();
        let mut value = 0.0;
        assert!(compute_dual_objective_value(
            None, &lp, &solution, &mut value, true
        ));
        assert_eq!(value, 10000000000000002.0);
    }

    #[test]
    fn h16_obj_cdouble_dual_objective_off() {
        let (lp, solution) = cancel_case();
        let mut value = 0.0;
        assert!(compute_dual_objective_value(
            None, &lp, &solution, &mut value, false
        ));
        assert_eq!(value, 1e16);
    }

    /// QP gradient term: cost [1e16, 1, 1], gradient 0, x = 1, duals 0, offset 0.
    /// quad = 1e16 + 2 in HighsCDouble (1e16 in double), dual objective 0.5 * quad.
    #[test]
    fn h16_obj_cdouble_dual_objective_gradient() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![1e16, 1.0, 1.0];
        lp.col_lower_ = vec![1.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        let solution = HighsSolution {
            value_valid: true,
            dual_valid: true,
            col_value: vec![1.0; 3],
            col_dual: vec![0.0; 3],
            row_value: vec![],
            row_dual: vec![],
        };
        let gradient = [0.0; 3];
        let mut on = 0.0;
        assert!(compute_dual_objective_value(
            Some(&gradient),
            &lp,
            &solution,
            &mut on,
            true
        ));
        assert_eq!(on, 5000000000000001.0);
        let mut off = 0.0;
        assert!(compute_dual_objective_value(
            Some(&gradient),
            &lp,
            &solution,
            &mut off,
            false
        ));
        assert_eq!(off, 5e15);
    }
}

/// The 1.16.0 use of `getEffectiveCosts` in `getKktFailures` (#3299, commit
/// cd50bf7f1): HighsSolution.cpp lines 169, 274, 283-284, 397, 433 (1.16.0).
///
/// Every expected value below is hand-computed from the C++ formulas, with the
/// default tolerances (primal/dual feasibility 1e-7, `kkt_tolerance` unset, so
/// the residual tolerances are `kDefaultKktTolerance` = 1e-7 too) and
/// `sense_ = kMinimize`, so `dual` is not flipped. All the values are exact in
/// binary (halves and quarters), so the assertions are exact.
#[cfg(test)]
mod h16_kkt_effective_costs_tests {
    use super::*;

    /// min z s.t. z - 3x = 0, z free, x in [0, 10]. z is a free column
    /// singleton in an equation, so `getEffectiveCosts` (HighsSolution.cpp
    /// 685-781, 1.16.0) gives effective costs [0, 3]: mu = cost_z / 1 = 1, so
    /// effective_costs[x] = 0 - 1 * (-3) = 3 and effective_costs[z] = 0. The
    /// raw costs are [1, 0]. Column-wise matrix, which `getEffectiveCosts`
    /// asserts (HighsSolution.cpp line 705, 1.16.0).
    fn objective_row_lp() -> HighsLp {
        let mut lp = HighsLp {
            num_col_: 2,
            num_row_: 1,
            col_cost_: vec![1.0, 0.0],
            col_lower_: vec![-kHighsInf, 0.0],
            col_upper_: vec![kHighsInf, 10.0],
            row_lower_: vec![0.0],
            row_upper_: vec![0.0],
            ..Default::default()
        };
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, -3.0];
        lp
    }

    /// x = 1 and z = 3 satisfy z - 3x = 0, so all primal infeasibilities and
    /// the primal residual are zero and only the dual measures matter.
    fn solution(col_dual: [f64; 2]) -> HighsSolution {
        HighsSolution {
            value_valid: true,
            dual_valid: true,
            col_value: vec![3.0, 1.0],
            col_dual: col_dual.to_vec(),
            row_value: vec![0.0],
            row_dual: vec![0.0],
        }
    }

    /// getLpKktFailures passes the LP costs as gradient (HighsSolution.cpp
    /// line 76, 1.16.0), and get_residuals = true exercises both denominators
    /// that the cost norm feeds.
    fn kkt(lp: &HighsLp, solution: &HighsSolution, h16_kkt_effective_costs: bool) -> HighsInfo {
        let options = HighsOptions {
            h16_kkt_effective_costs,
            ..Default::default()
        };
        let mut info = HighsInfo::default();
        get_kkt_failures(
            &options,
            false,
            lp,
            &lp.col_cost_,
            solution,
            &mut info,
            true,
        );
        info
    }

    /// The premise: z is a free column singleton in an equation row.
    #[test]
    fn h16_kkt_effective_costs_premise() {
        assert_eq!(
            crate::lp_data::highs_solution::get_effective_costs(
                &objective_row_lp(),
                &HighsOptions::default()
            ),
            vec![0.0, 3.0]
        );
    }

    /// z has dual 0 (so it feeds the cost norm) and x has dual -0.25.
    ///
    /// x: value 1 < (0 + 10) / 2, so getVariableKktFailures gives
    /// dual_infeasibility = max(-dual, 0) = 0.25, and its raw cost is 0, so the
    /// local measure of line 388-397 (1.16.0) is not applied: the relative
    /// value is dual_infeasibility / (1 + cost_norm).
    ///
    /// Cost norm over the columns with dual * dual < 1e-7, i.e. only z:
    /// 1.15.1 max(|cost_z|) = 1, 1.16.0 max(|effective_z|) = 0. Hence
    /// relative dual infeasibility = 0.25 / (1 + 1) = 0.125 off and
    /// 0.25 / (1 + 0) = 0.25 on.
    ///
    /// Dual residual of z: |A^T y - c + d| = |0 - 1 + 0| = 1, so
    /// relative dual residual error (line 433, 1.16.0) = 1 / (1 + 1) = 0.5 off
    /// and 1 / (1 + 0) = 1 on. x contributes a 0.25 residual in both, smaller.
    #[test]
    fn h16_kkt_effective_costs_relative_measures_off_and_on() {
        let lp = objective_row_lp();
        let solution = solution([0.0, -0.25]);
        let off = kkt(&lp, &solution, false);
        let on = kkt(&lp, &solution, true);

        // The absolute measures are identical: only the denominators change.
        for info in [&off, &on] {
            assert_eq!(info.num_dual_infeasibilities, 1);
            assert_eq!(info.max_dual_infeasibility, 0.25);
            assert_eq!(info.sum_dual_infeasibilities, 0.25);
            assert_eq!(info.num_relative_dual_infeasibilities, 1);
            assert_eq!(info.num_dual_residual_errors, 2);
            assert_eq!(info.max_dual_residual_error, 1.0);
            assert_eq!(info.num_relative_dual_residual_errors, 2);
            assert_eq!(info.num_primal_infeasibilities, 0);
            assert_eq!(info.max_primal_residual_error, 0.0);
        }
        assert_eq!(off.max_relative_dual_infeasibility, 0.125);
        assert_eq!(on.max_relative_dual_infeasibility, 0.25);
        assert_eq!(off.max_relative_dual_residual_error, 0.5);
        assert_eq!(on.max_relative_dual_residual_error, 1.0);
    }

    /// Here the dual-infeasible column is z itself (free, so
    /// dual_infeasibility = |dual| = 0.25) and x has dual 0, so x feeds the
    /// cost norm instead: 1.15.1 max(|cost_x|) = 0, 1.16.0 max(|effective_x|) =
    /// 3. z has raw cost 1 and a dual of 0.25, so the local measure of
    /// HighsSolution.cpp line 397 (1.16.0) applies to it:
    /// 1.15.1 max(|cost_z|, 0) = 1 -> 0.25 / (1 + 1) = 0.125;
    /// 1.16.0 max(|effective_z|, 3) = 3 -> 0.25 / (1 + 3) = 0.0625.
    ///
    /// The dual residual of z is |0 - 1 + 0.25| = 0.75, the only one, so the
    /// relative dual residual error is 0.75 / (1 + 0) = 0.75 off and
    /// 0.75 / (1 + 3) = 0.1875 on.
    #[test]
    fn h16_kkt_effective_costs_local_relative_measure_off_and_on() {
        let lp = objective_row_lp();
        let solution = solution([0.25, 0.0]);
        let off = kkt(&lp, &solution, false);
        let on = kkt(&lp, &solution, true);

        for info in [&off, &on] {
            assert_eq!(info.num_dual_infeasibilities, 1);
            assert_eq!(info.max_dual_infeasibility, 0.25);
            assert_eq!(info.num_dual_residual_errors, 1);
            assert_eq!(info.max_dual_residual_error, 0.75);
        }
        assert_eq!(off.max_relative_dual_infeasibility, 0.125);
        assert_eq!(on.max_relative_dual_infeasibility, 0.0625);
        assert_eq!(off.max_relative_dual_residual_error, 0.75);
        assert_eq!(on.max_relative_dual_residual_error, 0.1875);
    }

    /// The option must not touch the matrix: a row-wise LP gets the same
    /// effective costs as a column-wise one (the C++ assert of line 705 rules
    /// the case out; here the port transposes a local copy).
    #[test]
    fn h16_kkt_effective_costs_row_wise_lp() {
        let lp = objective_row_lp();
        let mut row_wise = lp.clone();
        row_wise.ensure_rowwise();
        assert!(!row_wise.a_matrix_.is_colwise());
        assert_eq!(
            kkt_effective_costs(&row_wise, &HighsOptions::default()),
            kkt_effective_costs(&lp, &HighsOptions::default())
        );
        // and the KKT measures of the two formats agree
        let solution = solution([0.0, -0.25]);
        let col_wise = kkt(&lp, &solution, true);
        let row_wise = kkt(&row_wise, &solution, true);
        assert_eq!(
            col_wise.max_relative_dual_infeasibility,
            row_wise.max_relative_dual_infeasibility
        );
        assert_eq!(col_wise.max_relative_dual_infeasibility, 0.25);
    }
}
