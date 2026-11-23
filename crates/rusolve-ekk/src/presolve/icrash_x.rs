//! HiGHS source: presolve/ICrashX.h, presolve/ICrashX.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Drives the Rust IPX port (`crate::ipm::ipx`) crossover from
//! the iCrash point: `ipx::LpSolver::LoadModel`, `LoadIPMStartingPoint`,
//! `CrossoverFromStartingPoint`, `GetBasicSolution`, then
//! `ipxBasicSolutionToHighsBasicSolution` (IpxWrapper, already ported).

use crate::ipm::ipx::lp_solver::LpSolver;
use crate::ipm::ipx::{self, IPX_STATUS_imprecise, IPX_STATUS_optimal, IPX_STATUS_time_limit, Int};
use crate::ipm::ipx_wrapper::{IpxSolution, fill_in_ipx_data, vec_data, vec_data_mut};
use crate::lp_data::highs_solution_ipx::ipx_basic_solution_to_highs_basic_solution;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g;
use crate::util::highs_cdouble::HighsCDouble;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: ICrashX.cpp::callCrossover (lines 17-173)
/// HiGHS 1.16 (#3466, #3185): ICrashX.cpp::callCrossover (lines 17-318)
///
/// PORT NOTE: every libc++ `std::vector::data()` argument goes through
/// [`vec_data`] / [`vec_data_mut`] (empty vector -> `None`, as `nullptr`).
/// The C++ `return`s after `lps.SetCallback(&highs_callback)` are `return`s
/// inside a closure, so the callback moved into `lps` is moved back to the
/// caller on every path.
pub fn call_crossover(
    options: &HighsOptions,
    lp: &HighsLp,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
    model_status: &mut HighsModelStatus,
    highs_info: &mut HighsInfo,
    highs_callback: &mut HighsCallback,
) -> HighsStatus {
    let mut num_col: Int = 0;
    let mut num_row: Int = 0;
    let mut offset: f64 = 0.0;
    let mut objective: Vec<f64> = Vec::new();
    let mut col_lb: Vec<f64> = Vec::new();
    let mut col_ub: Vec<f64> = Vec::new();
    let mut Ap: Vec<Int> = Vec::new();
    let mut Ai: Vec<Int> = Vec::new();
    let mut Av: Vec<f64> = Vec::new();
    let mut rhs: Vec<f64> = Vec::new();
    let mut constraint_type: Vec<u8> = Vec::new();

    fill_in_ipx_data(
        lp,
        &mut num_col,
        &mut num_row,
        &mut offset,
        &mut objective,
        &mut col_lb,
        &mut col_ub,
        &mut Ap,
        &mut Ai,
        &mut Av,
        &mut rhs,
        &mut constraint_type,
        // rusolve (HiGHS 1.16 backport, default off): compile-only plumbing of the
        // fill_in_ipx_data argument; read by the #3319 free-row port
        options.h16_ipx_free_rows,
    );
    // C++ line 30: if (res != IpxStatus::OK) return HighsStatus::kError; (commented out)

    let log_options = &options.log_options;
    // HiGHS 1.16 #3466 #3185 (h16_crossover_x): the one read of this option
    let h16_crossover_x = options.h16_crossover_x;

    let mut parameters = ipx::Parameters::default();
    parameters.run_crossover = 1; // 1 = "on"
    parameters.crash_basis = 1; // 0 = slack basis; 1 = crash basis
    parameters.display = 1;
    if !options.output_flag {
        parameters.display = 0;
    }
    // Modify parameters.debug according to log_dev_level
    parameters.debug = 0;
    if options.log_dev_level == kHighsLogDevLevelDetailed {
        // Default options.log_dev_level setting is kHighsLogDevLevelNone, yielding
        // default setting debug = 0
    } else if options.log_dev_level == kHighsLogDevLevelInfo {
        parameters.debug = 2;
    } else if options.log_dev_level == kHighsLogDevLevelVerbose {
        parameters.debug = 4;
    }
    parameters.highs_logging = true;
    parameters.log_options = Some(options.log_options.clone());
    // rusolve (HiGHS 1.16 backport, default off): IPX flag plumbing, as in
    // solve_lp_ipx; all false (the Parameters defaults) when the options are off
    parameters.h16_ipx_stall = options.h16_ipx_stall;
    parameters.h16_ipx_scaled = options.h16_ipx_scaled;
    parameters.h16_ipm_gap_offset = options.h16_ipm_gap_offset;
    // HiGHS 1.16 #3185 (h16_crossover_x): ICrashX.cpp line 51 (1.16.0). A separate
    // h16_ipx_timeout option was not added, so this option
    // gates the line.
    if h16_crossover_x {
        parameters.time_limit = options.time_limit;
    }

    let mut lps = LpSolver::new();
    lps.set_parameters(parameters);

    // Set pointer to any callback. The callback is moved into `lps` here and
    // moved back below, after the closure has returned on any path.
    lps.set_callback(Some(std::mem::take(highs_callback)));

    let status = (|| -> HighsStatus {
        let load_status = lps.load_model(
            num_col,
            offset,
            vec_data(&objective),
            vec_data(&col_lb),
            vec_data(&col_ub),
            num_row,
            vec_data(&Ap),
            vec_data(&Ai),
            vec_data(&Av),
            vec_data(&rhs),
            vec_data(&constraint_type),
        );
        if load_status != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Error loading ipx model\n"
            );
            return HighsStatus::kError;
        }

        // HiGHS 1.16 #3466 (h16_crossover_x): ICrashX.cpp lines 68-251 (1.16.0)
        let crossover_status: Int = if h16_crossover_x {
            h16_crossover_from_start(
                log_options,
                lp,
                highs_solution,
                num_col,
                num_row,
                &objective,
                &col_lb,
                &col_ub,
                &Ap,
                &Ai,
                &Av,
                &rhs,
                &constraint_type,
                &mut lps,
            )
        } else {
            // Set x values within bounds.
            let mut x: Vec<f64> = highs_solution.col_value.clone();
            for i in 0..num_col {
                x[i as usize] = cpp_max(x[i as usize], col_lb[i as usize]);
                x[i as usize] = cpp_min(x[i as usize], col_ub[i as usize]);
            }

            // Build slack variables from rhs-A*x but subject to sign conditions.
            let mut slack: Vec<f64> = rhs.clone();
            for i in 0..num_col {
                let mut p = Ap[i as usize];
                while p < Ap[(i + 1) as usize] {
                    slack[Ai[p as usize] as usize] -= Av[p as usize] * x[i as usize];
                    p += 1;
                }
            }
            for i in 0..num_row {
                match constraint_type[i as usize] {
                    b'=' => {
                        slack[i as usize] = 0.0;
                    }
                    b'<' => {
                        slack[i as usize] = cpp_max(slack[i as usize], 0.0);
                    }
                    b'>' => {
                        slack[i as usize] = cpp_min(slack[i as usize], 0.0);
                    }
                    _ => {}
                }
            }

            let crossover_status: Int = if highs_solution.dual_valid
                && (highs_solution.col_dual.len() as HighsInt) == num_col
                && (highs_solution.row_dual.len() as HighsInt) == num_row
            {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kInfo,
                    "Calling IPX crossover with primal and dual values\n"
                );
                lps.crossover_from_starting_point(
                    vec_data(&x),
                    vec_data(&slack),
                    vec_data(&highs_solution.row_dual),
                    vec_data(&highs_solution.col_dual),
                )
            } else {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kInfo,
                    "Calling IPX crossover with only primal values\n"
                );
                lps.crossover_from_starting_point(vec_data(&x), vec_data(&slack), None, None)
            };
            crossover_status
        };

        if crossover_status != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "IPX crossover error: flag = {}\n",
                crossover_status as i32
            );
            return HighsStatus::kError;
        }

        let ipx_info = lps.get_info();
        highs_info.crossover_iteration_count += ipx_info.updates_crossover as HighsInt;
        let imprecise_solution = ipx_info.status_crossover == IPX_STATUS_imprecise;
        if ipx_info.status_crossover != IPX_STATUS_optimal
            && ipx_info.status_crossover != IPX_STATUS_imprecise
            && ipx_info.status_crossover != IPX_STATUS_time_limit
        {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "IPX crossover failed: status = {}\n",
                ipx_info.status_crossover
            );
            return HighsStatus::kError;
        }
        if ipx_info.status_crossover == IPX_STATUS_time_limit {
            *model_status = HighsModelStatus::kTimeLimit;
            return HighsStatus::kWarning;
        }

        // Get basis
        let mut ipx_solution = IpxSolution {
            num_col,
            num_row,
            ipx_col_value: vec![0.0; num_col as usize],
            ipx_row_value: vec![0.0; num_row as usize],
            ipx_col_dual: vec![0.0; num_col as usize],
            ipx_row_dual: vec![0.0; num_row as usize],
            ipx_row_status: vec![0; num_row as usize],
            ipx_col_status: vec![0; num_col as usize],
        };

        let errflag = lps.get_basic_solution(
            vec_data_mut(&mut ipx_solution.ipx_col_value),
            vec_data_mut(&mut ipx_solution.ipx_row_value),
            vec_data_mut(&mut ipx_solution.ipx_row_dual),
            vec_data_mut(&mut ipx_solution.ipx_col_dual),
            vec_data_mut(&mut ipx_solution.ipx_row_status),
            vec_data_mut(&mut ipx_solution.ipx_col_status),
        );
        if errflag != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "IPX crossover getting basic solution: flag = {}\n",
                errflag as i32
            );
            return HighsStatus::kError;
        }

        // Convert the IPX basic solution to a HiGHS basic solution
        let status = ipx_basic_solution_to_highs_basic_solution(
            &options.log_options,
            lp,
            &rhs,
            &constraint_type,
            &ipx_solution,
            highs_basis,
            highs_solution,
        );

        if status != HighsStatus::kOk {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Failed to convert IPX basic solution to Highs basic solution\n"
            );
            return HighsStatus::kError;
        }

        highs_info.basis_validity = if highs_basis.valid {
            kBasisValidityValid
        } else {
            kBasisValidityInvalid
        };

        if imprecise_solution {
            *model_status = HighsModelStatus::kUnknown;
            HighsStatus::kWarning
        } else {
            *model_status = HighsModelStatus::kOptimal;
            HighsStatus::kOk
        }
    })();

    // The callback goes back to the caller on every path (see above).
    *highs_callback = lps
        .take_callback()
        .expect("IPX callback was moved in by SetCallback");
    status
}

/// HiGHS 1.16 (#3466): ICrashX.cpp::callCrossover (lines 68-251)
///
/// The 1.16 starting point (LP columns plus the explicit slacks of boxed
/// rows, free rows skipped, maximisation duals negated, made exactly
/// complementary) and its `CrossoverFromStartingPoint` call. Only the
/// `h16_crossover_x` branch of [`call_crossover`] calls it.
#[allow(clippy::too_many_arguments, non_snake_case)]
fn h16_crossover_from_start(
    log_options: &HighsLogOptions,
    lp: &HighsLp,
    highs_solution: &HighsSolution,
    num_col: Int,
    num_row: Int,
    objective: &[f64],
    col_lb: &[f64],
    col_ub: &[f64],
    Ap: &[Int],
    Ai: &[Int],
    Av: &[f64],
    rhs: &[f64],
    constraint_type: &[u8],
    lps: &mut LpSolver,
) -> Int {
    // IPX's columns are the LP's plus one slack per row with two
    // finite, distinct bounds (fillInIpxData; coefficient -1, the row
    // becomes '='), and its rows are the LP's less any free rows.
    let use_duals = highs_solution.dual_valid
        && (highs_solution.col_dual.len() as HighsInt) == lp.num_col_
        && (highs_solution.row_dual.len() as HighsInt) == lp.num_row_;

    // IPX accepts only an exactly complementary start: a column with z
    // > 0 must sit at its lower bound, z < 0 at its upper; a row's
    // slack must be at its bound wherever its dual is nonzero. An
    // approximate (first-order or early-stopped) point is not, and IPX
    // would refuse it (invalid vector).  Per coordinate, either move
    // the primal onto the bound or zero the dual, whichever is the
    // smaller relative change.
    let mut cost_scale: f64 = 1.0;
    for &cost in objective {
        cost_scale = cpp_max(cost_scale, cost.abs());
    }
    let mut num_snapped: HighsInt = 0;
    let mut num_dual_zeroed: HighsInt = 0;
    let mut sum_snapped: f64 = 0.0;
    let mut sum_dual_zeroed: f64 = 0.0;
    // C++ line 90-91: use_primal_tolerance = 0 is declared but never read
    let use_dual_tolerance: f64 = 0.0; // options.dual_feasibility_tolerance;

    let mut x: Vec<f64> = vec![0.0; num_col as usize];
    let mut y: Vec<f64> = Vec::new();
    let mut z: Vec<f64> = Vec::new();

    // Duals of maximizations must be negated
    let sense: HighsInt = if lp.sense_ == ObjSense::kMinimize {
        1
    } else {
        -1
    };

    // Boxed rows l <= a^Tx <= u are converted using explicit slack s to
    // a^Tx - s = 0; l <= s <= u
    //
    // Hence, i \in [lp.num_col_, num_col), Ap[i + 1] == Ap[i] + 1,
    // Av[Ap[i]] = -1 and rhs[Ai[Ap[i]]] = 0
    //
    // Since IPX excludes any free rows, skip them in mapping the LP's
    // row duals onto y
    let mut ipx_num_row: HighsInt = 0;
    for i_row in 0..lp.num_row_ {
        let r = i_row as usize;
        if lp.row_lower_[r] <= -kHighsInf && lp.row_upper_[r] >= kHighsInf {
            // Free row not in IPX: any dual it has should be small, but is ignored
            //
            // PORT NOTE: C++ reads row_dual[iRow] even when use_duals is false,
            // undefined behaviour when row_dual is shorter than num_row_. The
            // two counters are only reported when use_duals is true (then the
            // read is in bounds), so an out-of-range read is skipped here.
            if let Some(&row_dual) = highs_solution.row_dual.get(r) {
                let dual = sense as f64 * row_dual;
                if dual != 0.0 {
                    num_dual_zeroed += 1;
                    sum_dual_zeroed += dual.abs();
                }
            }
        } else {
            if use_duals {
                y.push(sense as f64 * highs_solution.row_dual[r]);
            }
            ipx_num_row += 1;
        }
    }
    // C++: ipx_num_row is set but never read
    let _ = ipx_num_row;

    for i in 0..lp.num_col_ as usize {
        x[i] = cpp_min(cpp_max(highs_solution.col_value[i], col_lb[i]), col_ub[i]);
    }
    if use_duals {
        z = vec![0.0; num_col as usize];
        for i in 0..lp.num_col_ as usize {
            z[i] = sense as f64 * highs_solution.col_dual[i];
        }
        for i in lp.num_col_ as usize..num_col as usize {
            debug_assert!(Av[Ap[i] as usize] == -1.0);
            z[i] = y[Ai[Ap[i] as usize] as usize];
        }
    }

    let mut complement = |value: &mut f64, dual: &mut f64, lower: f64, upper: f64| {
        if *dual > use_dual_tolerance {
            if *value == lower {
                return;
            }
            if lower.is_finite() && (*value - lower) / (1.0 + lower.abs()) <= *dual / cost_scale {
                sum_snapped += (*value - lower).abs();
                num_snapped += 1;
                *value = lower;
            } else {
                sum_dual_zeroed += dual.abs();
                num_dual_zeroed += 1;
                *dual = 0.0;
            }
        } else if *dual < -use_dual_tolerance {
            if *value == upper {
                return;
            }
            if upper.is_finite() && (upper - *value) / (1.0 + upper.abs()) <= -*dual / cost_scale {
                sum_snapped += (*value - upper).abs();
                num_snapped += 1;
                *value = upper;
            } else {
                sum_dual_zeroed += dual.abs();
                num_dual_zeroed += 1;
                *dual = 0.0;
            }
        }
    };
    if use_duals {
        for i in 0..lp.num_col_ as usize {
            complement(&mut x[i], &mut z[i], col_lb[i], col_ub[i]);
        }
    }

    // Compute the row activites for the original LP
    let mut activity: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); num_row as usize];
    for i in 0..lp.num_col_ as usize {
        let mut p = Ap[i];
        while p < Ap[i + 1] {
            activity[Ai[p as usize] as usize] += Av[p as usize] * x[i];
            p += 1;
        }
    }

    if num_col > lp.num_col_ {
        for i in lp.num_col_ as usize..num_col as usize {
            let row = Ai[Ap[i] as usize] as usize;
            let value = f64::from(activity[row]);
            x[i] = cpp_min(cpp_max(value, col_lb[i]), col_ub[i]);
            if use_duals {
                complement(&mut x[i], &mut z[i], col_lb[i], col_ub[i]);
                debug_assert!(rhs[row] == 0.0);
                y[row] = z[i];
            }
        }
    }

    // Build slack variables from rhs-A*x but subject to sign conditions.
    // Complete the calculation of the activities for IPX
    for i in lp.num_col_ as usize..num_col as usize {
        activity[Ai[Ap[i] as usize] as usize] -= x[i];
    }

    let mut slack: Vec<f64> = vec![0.0; num_row as usize];
    for i in 0..num_row as usize {
        slack[i] = f64::from(HighsCDouble::from(rhs[i]) - activity[i]);
        match constraint_type[i] {
            b'=' => {
                slack[i] = 0.0;
            }
            b'<' => {
                // slack in [0, inf), reduced cost -y
                slack[i] = cpp_max(slack[i], 0.0);
                if use_duals {
                    let mut dual = -y[i];
                    complement(&mut slack[i], &mut dual, 0.0, kHighsInf);
                    y[i] = -dual;
                }
            }
            b'>' => {
                // slack in (-inf, 0], reduced cost -y
                slack[i] = cpp_min(slack[i], 0.0);
                if use_duals {
                    let mut dual = -y[i];
                    complement(&mut slack[i], &mut dual, -kHighsInf, 0.0);
                    y[i] = -dual;
                }
            }
            _ => {}
        }
    }

    let crossover_status: HighsInt = if use_duals {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "Calling IPX crossover with primal and dual values"
        );
        if num_snapped != 0 || num_dual_zeroed != 0 {
            crate::highs_log_user!(
                log_options,
                HighsLogType::kInfo,
                ": made complementary with\n"
            );
            if num_snapped != 0 {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kInfo,
                    "   {} primal value{} moved to a bound (total changes = {})\n",
                    num_snapped,
                    if num_snapped == 1 { "" } else { "s" },
                    c_fmt_g(sum_snapped, 6)
                );
            }
            if num_dual_zeroed != 0 {
                crate::highs_log_user!(
                    log_options,
                    HighsLogType::kInfo,
                    "   {} dual{} zeroed (total changes = {})\n",
                    num_dual_zeroed,
                    if num_dual_zeroed == 1 { "" } else { "s" },
                    c_fmt_g(sum_dual_zeroed, 6)
                );
            }
        } else {
            crate::highs_log_user!(log_options, HighsLogType::kInfo, "\n");
        }
        lps.crossover_from_starting_point(
            vec_data(&x),
            vec_data(&slack),
            vec_data(&y),
            vec_data(&z),
        )
    } else {
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "Calling IPX crossover with only primal values\n"
        );
        lps.crossover_from_starting_point(vec_data(&x), vec_data(&slack), None, None)
    };
    crossover_status
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

    fn quiet(lp: HighsLp) -> Highs {
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        highs
    }

    #[test]
    fn icrash_then_crossover_reaches_optimum() {
        // min x0 + 2 x1 s.t. x0 + x1 = 1, x >= 0: optimum (1, 0), objective 1
        let lp = lp(
            &[
                (1.0, 0.0, kHighsInf, &[(0, 1.0)]),
                (2.0, 0.0, kHighsInf, &[(0, 1.0)]),
            ],
            &[(1.0, 1.0)],
        );
        let mut highs = quiet(lp);
        assert_eq!(
            highs.set_option_value_bool("icrash", true),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_int("simplex_strategy", kSimplexStrategyPrimal),
            HighsStatus::kOk
        );
        assert_ne!(highs.run(), HighsStatus::kError);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        assert!((highs.get_info().objective_function_value - 1.0).abs() <= 1e-9);
        assert_eq!(highs.icrash_info_.x_values.len(), 2);
    }

    // -----------------------------------------------------------------------
    // HiGHS 1.16 #3466 (TestCrossover.cpp "issue-3459", PR diff
    // lines 1-127) and #3185 (crossover time limit), behind h16_crossover_x.
    // -----------------------------------------------------------------------

    /// 3459.mps (3466.diff lines 128-412): columns in file order as
    /// (cost, [(row, value)]); every row is `L` with a RANGES entry.
    #[allow(clippy::type_complexity)]
    const ISSUE_3459_COLS: &[(f64, &[(HighsInt, f64)])] = &[
        (
            2.0,
            &[
                (0, 3.0),
                (1, 3.0),
                (2, 3.0),
                (3, 1.0),
                (7, 1.0),
                (8, 2.0),
                (12, 2.0),
            ],
        ),
        (8.0, &[(3, 3.0), (5, 1.0), (11, 1.0), (15, 2.0), (19, 2.0)]),
        (
            7.0,
            &[
                (3, 1.0),
                (4, 3.0),
                (5, 3.0),
                (7, 2.0),
                (12, 3.0),
                (15, 3.0),
                (16, 1.0),
            ],
        ),
        (5.0, &[(0, 1.0), (7, 1.0), (10, 1.0), (15, 1.0)]),
        (4.0, &[(3, 2.0), (6, 3.0), (8, 3.0), (13, 2.0), (19, 1.0)]),
        (8.0, &[(1, 2.0), (7, 3.0), (11, 3.0)]),
        (
            4.0,
            &[
                (4, 2.0),
                (5, 2.0),
                (9, 2.0),
                (13, 2.0),
                (14, 2.0),
                (16, 2.0),
                (18, 1.0),
                (19, 2.0),
            ],
        ),
        (
            3.0,
            &[
                (4, 2.0),
                (8, 2.0),
                (10, 2.0),
                (17, 2.0),
                (18, 2.0),
                (19, 2.0),
            ],
        ),
        (4.0, &[(2, 3.0), (13, 2.0), (15, 1.0)]),
        (
            7.0,
            &[
                (4, 1.0),
                (5, 1.0),
                (8, 1.0),
                (9, 1.0),
                (11, 2.0),
                (15, 1.0),
                (17, 3.0),
                (18, 3.0),
            ],
        ),
        (
            5.0,
            &[
                (5, 3.0),
                (12, 3.0),
                (13, 3.0),
                (14, 1.0),
                (16, 1.0),
                (19, 1.0),
            ],
        ),
        (7.0, &[(3, 3.0), (4, 2.0), (9, 2.0), (18, 3.0), (19, 1.0)]),
        (
            4.0,
            &[
                (0, 2.0),
                (4, 3.0),
                (11, 2.0),
                (12, 2.0),
                (13, 3.0),
                (14, 3.0),
                (19, 3.0),
            ],
        ),
        (5.0, &[(2, 1.0), (5, 1.0), (15, 1.0), (16, 2.0), (17, 1.0)]),
        (
            3.0,
            &[
                (2, 2.0),
                (4, 3.0),
                (5, 1.0),
                (6, 2.0),
                (10, 3.0),
                (11, 2.0),
                (16, 1.0),
                (17, 2.0),
            ],
        ),
        (
            2.0,
            &[
                (2, 3.0),
                (3, 3.0),
                (5, 2.0),
                (7, 2.0),
                (9, 3.0),
                (11, 1.0),
                (15, 2.0),
            ],
        ),
        (
            4.0,
            &[
                (0, 2.0),
                (7, 1.0),
                (12, 2.0),
                (13, 2.0),
                (14, 2.0),
                (18, 1.0),
            ],
        ),
        (4.0, &[(6, 3.0), (15, 3.0), (19, 1.0)]),
        (8.0, &[(8, 3.0), (17, 1.0), (18, 2.0)]),
        (
            3.0,
            &[
                (0, 3.0),
                (1, 2.0),
                (5, 1.0),
                (6, 3.0),
                (12, 2.0),
                (14, 3.0),
                (16, 3.0),
                (19, 1.0),
            ],
        ),
        (
            8.0,
            &[
                (1, 3.0),
                (3, 2.0),
                (4, 2.0),
                (8, 3.0),
                (10, 1.0),
                (13, 3.0),
                (17, 2.0),
            ],
        ),
        (
            8.0,
            &[
                (6, 1.0),
                (7, 2.0),
                (8, 2.0),
                (16, 1.0),
                (17, 1.0),
                (18, 2.0),
            ],
        ),
        (
            2.0,
            &[
                (0, 2.0),
                (2, 1.0),
                (3, 3.0),
                (4, 3.0),
                (6, 1.0),
                (7, 3.0),
                (9, 3.0),
                (11, 1.0),
                (12, 3.0),
                (13, 3.0),
                (14, 1.0),
            ],
        ),
        (6.0, &[(1, 3.0), (4, 3.0), (10, 1.0), (16, 1.0)]),
        (
            6.0,
            &[
                (0, 1.0),
                (2, 2.0),
                (3, 3.0),
                (9, 2.0),
                (10, 1.0),
                (11, 2.0),
                (16, 2.0),
                (18, 2.0),
            ],
        ),
        (
            8.0,
            &[
                (2, 2.0),
                (3, 3.0),
                (4, 1.0),
                (5, 3.0),
                (11, 3.0),
                (18, 3.0),
                (19, 3.0),
            ],
        ),
        (
            8.0,
            &[
                (1, 1.0),
                (3, 2.0),
                (9, 3.0),
                (10, 1.0),
                (12, 1.0),
                (13, 2.0),
                (15, 3.0),
            ],
        ),
        (3.0, &[(1, 2.0), (6, 3.0), (13, 1.0), (18, 1.0)]),
        (
            5.0,
            &[
                (2, 3.0),
                (3, 3.0),
                (6, 3.0),
                (8, 1.0),
                (9, 3.0),
                (10, 2.0),
                (12, 2.0),
                (13, 2.0),
            ],
        ),
        (
            4.0,
            &[
                (2, 1.0),
                (6, 1.0),
                (10, 3.0),
                (16, 1.0),
                (18, 3.0),
                (19, 3.0),
            ],
        ),
    ];
    const ISSUE_3459_RHS: [f64; 20] = [
        6.269489302794252,
        5.4189127335612373,
        5.825998508615795,
        6.3972282288045239,
        4.1131438974023098,
        5.3592446195886989,
        4.0008867988411927,
        5.9301504880066922,
        6.4320026603217606,
        4.3719618402626708,
        4.9683776085791873,
        7.2159672842360543,
        7.3611262393983914,
        5.5509330170814586,
        7.2568949221649071,
        5.1085610121285709,
        6.8244328892539272,
        6.1818264946284582,
        5.7603963156106497,
        6.6257691046951059,
    ];
    const ISSUE_3459_RANGE: [f64; 20] = [
        3.0000000000000004,
        3.0,
        3.0,
        3.0,
        3.0000000000000004,
        3.0000000000000004,
        3.0,
        3.0000000000000004,
        2.9999999999999996,
        3.0000000000000004,
        3.0000000000000004,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0000000000000004,
        3.0000000000000004,
        3.0,
        2.9999999999999996,
    ];

    /// The 3459.mps LP as HMpsFF reads it: `L` row upper = RHS, lower =
    /// upper - |RANGE| (HMpsFF.cpp lines 1556-1558), columns in [0, inf).
    #[allow(clippy::excessive_precision)]
    fn issue_3459_mps_lp() -> HighsLp {
        let cols: Vec<(f64, f64, f64, &[(HighsInt, f64)])> = ISSUE_3459_COLS
            .iter()
            .map(|&(cost, entries)| (cost, 0.0, kHighsInf, entries))
            .collect();
        let rows: Vec<(f64, f64)> = (0..20)
            .map(|i| {
                (
                    ISSUE_3459_RHS[i] - ISSUE_3459_RANGE[i].abs(),
                    ISSUE_3459_RHS[i],
                )
            })
            .collect();
        lp(&cols, &rows)
    }

    /// TestCrossover.cpp issue-3459, k == 1 (3466.diff lines 35-74): the
    /// maximisation with boxed columns and boxed rows.
    fn issue_3459_max_lp() -> HighsLp {
        let mut lp = lp(
            &[
                (-10.0, 20.0, 30.0, &[(0, 1.0), (1, 1.0), (2, 1.0)]),
                (-16.0, 5.0, 7.0, &[(0, 1.25), (1, 2.25), (2, 3.0)]),
                (-24.0, 20.0, 40.0, &[(0, 2.0), (1, 3.0), (2, 4.0)]),
            ],
            &[(70.0, 90.0), (100.0, 120.0), (125.0, 145.0)],
        );
        lp.sense_ = ObjSense::kMaximize;
        for i_col in 0..lp.num_col_ as usize {
            lp.col_cost_[i_col] *= -1.0;
        }
        lp
    }

    /// Whether the #3319 free-row port is in this tree: with
    /// h16_ipx_free_rows on, the slack of a boxed row that follows a free row
    /// is in reduced row 0 (IpxWrapper.cpp line 770, 1.16.0), not row 1.
    fn ipx_free_rows_ported() -> bool {
        let lp = lp(
            &[(1.0, 0.0, kHighsInf, &[(0, 1.0), (1, 1.0)])],
            &[(-kHighsInf, kHighsInf), (1.0, 2.0)],
        );
        let (mut num_col, mut num_row, mut offset) = (0, 0, 0.0);
        let (mut obj, mut col_lb, mut col_ub) = (Vec::new(), Vec::new(), Vec::new());
        let (mut ap, mut ai, mut av) = (Vec::new(), Vec::new(), Vec::new());
        let (mut rhs, mut ctype) = (Vec::new(), Vec::new());
        fill_in_ipx_data(
            &lp,
            &mut num_col,
            &mut num_row,
            &mut offset,
            &mut obj,
            &mut col_lb,
            &mut col_ub,
            &mut ap,
            &mut ai,
            &mut av,
            &mut rhs,
            &mut ctype,
            true,
        );
        ai.last() == Some(&0)
    }

    /// TestCrossover.cpp issue-3459 loop body (3466.diff lines 76-127), for
    /// the solvers built here (kPdlpString, kIpxString; HiPO is `#ifdef HIPO`).
    fn issue_3459_check(
        lp: HighsLp,
        objective_value: f64,
        free_row0: HighsInt,
        free_row1: HighsInt,
    ) {
        let mut h = quiet(lp);
        assert_eq!(
            h.set_option_value_bool("h16_crossover_x", true),
            HighsStatus::kOk
        );
        assert_eq!(
            h.set_option_value_str("presolve", kHighsOffString),
            HighsStatus::kOk
        );
        assert_eq!(
            h.set_option_value_str("run_crossover", kHighsOffString),
            HighsStatus::kOk
        );

        let objective_error = |h: &Highs| {
            (h.get_info().objective_function_value - objective_value).abs()
                / (1.0 + objective_value.abs())
        };
        let test_solvers = |h: &mut Highs, kind: &str| {
            for solver in [kPdlpString, kIpxString] {
                assert_eq!(h.set_option_value_str("solver", solver), HighsStatus::kOk);
                assert_eq!(h.run(), HighsStatus::kOk, "{solver} ({kind}) run");
                assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
                assert!(objective_error(h) < 1e-6, "{solver} ({kind}) run objective");
                let solution = h.get_solution().clone();
                assert_eq!(
                    h.crossover(&solution),
                    HighsStatus::kOk,
                    "{solver} ({kind}) crossover"
                );
                assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
                assert!(
                    objective_error(h) < 1e-6,
                    "{solver} ({kind}) crossover objective"
                );
            }
        };

        test_solvers(&mut h, "original");

        // The free-row half needs #3319 too (1.16 fillInIpxData maps a boxed
        // row's slack through reduced_rowmap; 1.15.1 does not, so IPX gets a
        // slack in the wrong row). That is the #3319 free-row port, a separate
        // backport: run it only once that port is in the tree.
        if !ipx_free_rows_ported() {
            eprintln!("#3319 free-row port not in this tree: free-row half of issue-3459 skipped");
            return;
        }
        assert_eq!(
            h.set_option_value_bool("h16_ipx_free_rows", true),
            HighsStatus::kOk
        );

        // Now free two of the rows to test the ability to handle them in
        // Highs::crossover
        let tol = h.get_options().dual_feasibility_tolerance;
        assert!(h.get_solution().row_dual[free_row0 as usize].abs() < tol);
        assert!(h.get_solution().row_dual[free_row1 as usize].abs() < tol);
        assert_eq!(
            h.change_row_bounds(free_row0, -kHighsInf, kHighsInf),
            HighsStatus::kOk
        );
        assert_eq!(
            h.change_row_bounds(free_row1, -kHighsInf, kHighsInf),
            HighsStatus::kOk
        );

        test_solvers(&mut h, "with free row");
    }

    #[test]
    fn h16_crossover_x_issue_3459_max_lp() {
        // TestCrossover.cpp: objective_value = -980, negated for kMaximize
        issue_3459_check(issue_3459_max_lp(), 980.0, 0, 1);
    }

    #[test]
    fn h16_crossover_x_issue_3459_mps() {
        // TestCrossover.cpp: objective_value = 17.555148585, free rows 10 and 14
        issue_3459_check(issue_3459_mps_lp(), 17.555148585, 10, 14);
    }

    /// min x0 + 2 x1 s.t. x0 + x1 = 1, x >= 0, with an approximate start:
    /// x1 = 1e-7 off its lower bound while its reduced cost is 1.
    fn approximate_start() -> (Highs, HighsSolution) {
        let lp = lp(
            &[
                (1.0, 0.0, kHighsInf, &[(0, 1.0)]),
                (2.0, 0.0, kHighsInf, &[(0, 1.0)]),
            ],
            &[(1.0, 1.0)],
        );
        let solution = HighsSolution {
            value_valid: true,
            dual_valid: true,
            col_value: vec![1.0 - 1e-7, 1e-7],
            row_value: vec![1.0],
            col_dual: vec![0.0, 1.0],
            row_dual: vec![1.0],
        };
        (quiet(lp), solution)
    }

    #[test]
    fn h16_crossover_x_on_snaps_approximate_start() {
        // ICrashX.cpp lines 120-157 (1.16.0): z1 = 1 > 0 and
        // (1e-7 - 0) / (1 + 0) <= 1 / cost_scale (= 2), so x1 moves to 0
        let (mut h, solution) = approximate_start();
        assert_eq!(
            h.set_option_value_bool("h16_crossover_x", true),
            HighsStatus::kOk
        );
        assert_eq!(h.crossover(&solution), HighsStatus::kOk);
        assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
        assert!((h.get_info().objective_function_value - 1.0).abs() <= 1e-9);
        assert!(h.basis_.valid);
    }

    #[test]
    fn h16_crossover_x_off_approximate_start() {
        // 1.15.1 passes the start unchanged
        let (mut h, solution) = approximate_start();
        let status = h.crossover(&solution);
        assert_eq!(status, HighsStatus::kError);
    }

    #[test]
    fn h16_crossover_x_time_limit() {
        // #3185 (ICrashX.cpp line 51, 1.16.0): with the flag on, IPX gets
        // options.time_limit; with it off, IPX keeps its default (-1, none).
        // The rows are one-sided (>=, so that the optimum is not x = 0 and
        // crossover has primal pushes to make): with boxed rows the 1.15.1 path indexes x
        // past col_value (ICrashX.cpp 1.15.1 lines 70-73, undefined
        // behaviour in C++, a panic here).
        for flag in [false, true] {
            let mut lp = issue_3459_mps_lp();
            for i_row in 0..lp.num_row_ as usize {
                lp.row_upper_[i_row] = kHighsInf;
            }
            let mut h = quiet(lp);
            assert_eq!(
                h.set_option_value_str("solver", kIpxString),
                HighsStatus::kOk
            );
            assert_eq!(
                h.set_option_value_str("run_crossover", kHighsOffString),
                HighsStatus::kOk
            );
            assert_eq!(h.run(), HighsStatus::kOk);
            let mut solution = h.get_solution().clone();
            solution.dual_valid = false;
            assert_eq!(
                h.set_option_value_bool("h16_crossover_x", flag),
                HighsStatus::kOk
            );
            assert_eq!(
                h.set_option_value_double("time_limit", 0.0),
                HighsStatus::kOk
            );
            let status = h.crossover(&solution);
            if flag {
                assert_eq!(status, HighsStatus::kWarning);
                assert_eq!(h.get_model_status(), HighsModelStatus::kTimeLimit);
            } else {
                assert_eq!(status, HighsStatus::kOk);
                assert_eq!(h.get_model_status(), HighsModelStatus::kOptimal);
            }
        }
    }
}
