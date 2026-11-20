//! HiGHS source: presolve/ICrash.h, presolve/ICrash.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Only the types are in the crate prelude; the free functions
//! have generic names (`initialize`, `update`) and are reached by path
//! (`crate::presolve::icrash::call_icrash`). ICrash.cpp's `initialize`
//! overloads ICrashUtil.cpp's: reach the latter as
//! `crate::presolve::icrash_util::initialize`.
//! `int` in this file is C `int` = `i32` (same as `HighsInt` here).
//! `callCrossover(lp, options, x_values, solution, basis, callback)` is
//! declared in ICrash.h (line 120) but never defined (ICrashX.cpp defines a
//! different overload), so it is not ported.

use crate::prelude::*;
use crate::presolve::icrash_util;
use crate::util::c_printf;

/// HiGHS: ICrash.cpp::kExitTolerance (line 28)
const kExitTolerance: f64 = 0.00000001;

/// HiGHS: ICrash.h::ICrashStrategy (lines 23-29)
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ICrashStrategy {
    kPenalty = 0,
    kAdmm,
    kICA,
    kUpdatePenalty,
    kUpdateAdmm,
}

/// HiGHS: ICrash.h::ICrashIterationDetails (lines 31-41)
#[derive(Clone, Copy, Debug, Default)]
pub struct ICrashIterationDetails {
    pub num: i32,
    pub weight: f64,
    pub lambda_norm_2: f64,

    pub lp_objective: f64,
    pub quadratic_objective: f64,
    pub residual_norm_2: f64,

    pub time: f64,
}

/// HiGHS: ICrash.h::ICrashInfo (lines 43-57)
///
/// C++ scalars are uninitialised; `Default` gives 0 (every reader runs after
/// `fillICrashInfo`).
#[derive(Clone, Debug, Default)]
pub struct ICrashInfo {
    pub num_iterations: i32,

    pub final_lp_objective: f64,
    pub final_quadratic_objective: f64,
    pub final_residual_norm_2: f64,

    pub starting_weight: f64,
    pub final_weight: f64,

    pub details: Vec<ICrashIterationDetails>,
    pub x_values: Vec<f64>,

    pub total_time: f64,
}

/// HiGHS: ICrash.h::ICrashOptions (lines 59-68)
#[derive(Clone)]
pub struct ICrashOptions {
    pub dualize: bool,
    pub strategy: ICrashStrategy,
    pub starting_weight: f64,
    pub iterations: HighsInt,
    pub approximate_minimization_iterations: HighsInt,
    pub exact: bool,
    /// gets ignored if exact is set to true
    pub breakpoints: bool,
    pub log_options: HighsLogOptions,
}

/// HiGHS: ICrash.h::Quadratic (lines 70-93)
#[derive(Clone)]
pub struct Quadratic {
    pub lp: HighsLp,
    pub options: ICrashOptions,
    pub details: Vec<ICrashIterationDetails>,

    pub xk: HighsSolution,

    pub lp_objective: f64,
    pub quadratic_objective: f64,
    pub residual: Vec<f64>,
    pub residual_norm_2: f64,

    pub mu: f64,
    pub lambda: Vec<f64>,
}

impl Quadratic {
    /// HiGHS: ICrash.h::Quadratic::Quadratic (lines 85-92)
    pub fn new(lp_: HighsLp, options_: ICrashOptions) -> Self {
        Self {
            lp: lp_,
            options: options_,
            details: Vec::new(),
            xk: HighsSolution::default(),
            lp_objective: 0.0,
            quadratic_objective: 0.0,
            residual: Vec::new(),
            residual_norm_2: 0.0,
            mu: 0.0,
            lambda: Vec::new(),
        }
    }
}

/// HiGHS: ICrash.cpp::parseICrashStrategy (lines 30-50)
pub fn parse_icrash_strategy(strategy: &str, icrash_strategy: &mut ICrashStrategy) -> bool {
    let mut lower = strategy.to_string();
    // trim(lower): strips leading/trailing chars in "\t\n\v\f\r "
    lower = lower
        .trim_matches(|ch: char| {
            ch == '\t' || ch == '\n' || ch == '\u{0B}' || ch == '\u{0C}' || ch == '\r' || ch == ' '
        })
        .to_string();
    // std::transform tolower
    lower = lower.to_ascii_lowercase();

    if lower == "penalty" {
        *icrash_strategy = ICrashStrategy::kPenalty;
    } else if lower == "admm" {
        *icrash_strategy = ICrashStrategy::kAdmm;
    } else if lower == "ica" {
        *icrash_strategy = ICrashStrategy::kICA;
    } else if lower == "update_penalty" {
        *icrash_strategy = ICrashStrategy::kUpdatePenalty;
    } else if lower == "update_admm" {
        *icrash_strategy = ICrashStrategy::kUpdateAdmm;
    } else {
        return false;
    }
    true
}

/// HiGHS: ICrash.cpp::checkOptions (lines 52-88)
pub fn check_options(lp: &HighsLp, options: ICrashOptions) -> bool {
    if options.exact {
        // std::cout << "ICrashError: exact subproblem solution not available "
        //              "at the moment." << std::endl;
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "ICrashError: exact subproblem solution not available at the moment.\n"
        );
        return false;
    }

    if options.breakpoints {
        if options.exact {
            // std::cout << "ICrashError: exact strategy not allowed for "
            //              "breakpoints minimization." << std::endl;
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "ICrashError: exact strategy not allowed for breakpoints minimization.\n"
            );
            return false;
        }
        if options.dualize {
            // std::cout << "ICrashError: breakpoints does not support dualize option."
            //           << std::endl;
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "ICrashError: breakpoints does not support dualize option.\n"
            );
            return false;
        }
        // std::cout << "ICrashError: breakpoints not implemented yet." << std::endl;
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kInfo,
            "ICrashError: breakpoints not implemented yet.\n"
        );
        return false;
    }

    true
}

/// HiGHS: ICrash.cpp::parseOptions (lines 90-128)
pub fn parse_options(lp: &HighsLp, options: ICrashOptions) -> Quadratic {
    let mut ilp = lp.clone();
    // HighsLp local_lp;
    // HighsStatus status;
    icrash_util::convert_to_minimization(&mut ilp);
    if icrash_util::is_equality_problem(&ilp) {
        if options.dualize {
            // status = dualizeEqualityProblem(ilp, local_lp);
            // if (status == HighsStatus::kOk) {
            //   ilp = local_lp;
            // } else {
            //   printf("Cannot dualize equality problem\n");
            // }
        }
    } else {
        // not equality problem.
        // assert(!options.breakpoints);  // remove when implementing breakpoints and
        //                                // add if else.
        // status = transformIntoEqualityProblem(ilp, local_lp);
        // if (status == HighsStatus::kOk) {
        //   ilp = local_lp;
        // } else {
        //   printf("Cannot transform into equality problem\n");
        // }
        // if (options.dualize) {
        //   // Add slacks & dualize.
        //   // dualizeEqualityProblem returns a minimization equality problem.
        //   // status = dualizeEqualityProblem(ilp, local_lp);
        //   // if (status == HighsStatus::kOk) {
        //   //   ilp = local_lp;
        //   // } else {
        //   //   printf("Cannot dualize equality problem\n");
        //   // }
        // }
    }

    Quadratic::new(ilp, options)
}

/// HiGHS: ICrash.cpp::getQuadraticObjective (lines 130-142)
pub fn get_quadratic_objective(idata: &Quadratic) -> f64 {
    // c'x
    let mut quadratic = icrash_util::vector_product(&idata.lp.col_cost_, &idata.xk.col_value);

    // lambda'x
    quadratic += icrash_util::vector_product(&idata.lambda, &idata.residual);

    // 1/2mu r'r
    let rtr = icrash_util::vector_product(&idata.residual, &idata.lambda);
    quadratic += rtr / (2.0 * idata.mu);

    quadratic
}

/// HiGHS: ICrash.cpp::initialize (lines 144-151)
pub fn initialize(idata: &mut Quadratic, options: &ICrashOptions) -> bool {
    if !icrash_util::initialize(&idata.lp, &mut idata.xk, &mut idata.lambda) {
        return false;
    }

    idata.mu = options.starting_weight;

    // Maybe other values for x0.
    true
}

/// HiGHS: ICrash.cpp::update (lines 153-167)
pub fn update(idata: &mut Quadratic) {
    // lp_objective
    idata.lp_objective = icrash_util::vector_product(&idata.lp.col_cost_, &idata.xk.col_value);

    // residual & residual_norm_2
    let _ =
        crate::lp_data::presolve_support::calculate_row_values_quad(&idata.lp, &mut idata.xk, -1);
    icrash_util::update_residual(
        idata.options.breakpoints,
        &idata.lp,
        &idata.xk,
        &mut idata.residual,
    );
    idata.residual_norm_2 = crate::util::highs_utils::get_norm2(&idata.residual);

    // quadratic_objective
    idata.quadratic_objective = idata.lp_objective;
    idata.quadratic_objective += icrash_util::vector_product(&idata.lambda, &idata.residual);
    idata.quadratic_objective +=
        icrash_util::vector_product(&idata.residual, &idata.residual) / (2.0 * idata.mu);
}

/// HiGHS: ICrash.cpp::fillDetails (lines 169-178)
pub fn fill_details(num: i32, idata: &Quadratic) -> ICrashIterationDetails {
    let lambda_norm_2 = crate::util::highs_utils::get_norm2(&idata.lambda);
    ICrashIterationDetails {
        num,
        weight: idata.mu,
        lambda_norm_2,
        lp_objective: idata.lp_objective,
        quadratic_objective: idata.quadratic_objective,
        residual_norm_2: idata.residual_norm_2,
        time: 0.0,
    }
}

/// HiGHS: ICrash.cpp::fillICrashInfo (lines 180-191)
pub fn fill_icrash_info(n_iterations: i32, result: &mut ICrashInfo) {
    debug_assert!(result.details.len() as i32 == n_iterations + 1);
    result.num_iterations = n_iterations;

    result.final_lp_objective = result.details[n_iterations as usize].lp_objective;
    result.final_quadratic_objective = result.details[n_iterations as usize].quadratic_objective;
    result.final_residual_norm_2 = result.details[n_iterations as usize].residual_norm_2;

    result.starting_weight = result.details[0].weight;
    result.final_weight = result.details[n_iterations as usize].weight;
}

/// HiGHS: ICrash.cpp::updateParameters (lines 193-242)
pub fn update_parameters(idata: &mut Quadratic, options: &ICrashOptions, iteration: i32) {
    if iteration == 1 {
        return;
    }

    // The other strategies are WIP.
    match options.strategy {
        ICrashStrategy::kPenalty => {
            idata.mu = 0.1 * idata.mu;
        }
        ICrashStrategy::kAdmm => {
            // std::cout << "ICrash Error: ADMM parameter update not implemented." << std::endl;
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "ICrashError: ADMM parameter update not implemented yet."
            );
        }
        ICrashStrategy::kICA => {
            // Update mu every third iteration, otherwise update lambda.
            if iteration % 3 == 0 {
                idata.mu = 0.1 * idata.mu;
            } else {
                let mut residual_ica = vec![0.0; idata.lp.num_row_ as usize];
                icrash_util::update_residual_ica(&idata.lp, &idata.xk, &mut residual_ica);
                for row in 0..idata.lp.num_row_ {
                    idata.lambda[row as usize] = idata.mu * residual_ica[row as usize];
                }
            }
        }
        ICrashStrategy::kUpdatePenalty => {
            // Update mu every third iteration, otherwise do nothing.
            if iteration % 3 == 0 {
                idata.mu = 0.1 * idata.mu;
            }
        }
        ICrashStrategy::kUpdateAdmm => {
            // Update mu every third iteration, otherwise update lambda.
            if iteration % 3 == 0 {
                idata.mu = 0.1 * idata.mu;
            } else {
                let _ = crate::lp_data::presolve_support::calculate_row_values_quad(
                    &idata.lp,
                    &mut idata.xk,
                    -1,
                );
                let mut residual = vec![0.0; idata.lp.num_row_ as usize];
                icrash_util::update_residual_fast(&idata.lp, &idata.xk, &mut residual);
                for row in 0..idata.lp.num_row_ {
                    // todo: double check clp.
                    idata.lambda[row as usize] =
                        idata.lambda[row as usize] + idata.mu * residual[row as usize];
                }
            }
        }
    }
}

/// HiGHS: ICrash.cpp::solveSubproblemICA (lines 244-284)
fn solve_subproblem_ica(idata: &mut Quadratic, options: &ICrashOptions) {
    let minor_iteration_details = false;

    let mut residual_ica = vec![0.0; idata.lp.num_row_ as usize];
    icrash_util::update_residual_ica(&idata.lp, &idata.xk, &mut residual_ica);
    let mut objective_ica = 0.0;

    for k in 0..options.approximate_minimization_iterations {
        for col in 0..idata.lp.num_col_ {
            // determine whether to minimize for col.
            // if empty skip.
            if idata.lp.a_matrix_.start_[col as usize]
                == idata.lp.a_matrix_.start_[(col + 1) as usize]
            {
                continue;
            }

            let old_value = idata.xk.col_value[col as usize];
            icrash_util::minimize_component_ica(
                col,
                idata.mu,
                &idata.lambda,
                &idata.lp,
                &mut objective_ica,
                &mut residual_ica,
                &mut idata.xk,
            );

            let new_value = idata.xk.col_value[col as usize];
            let delta_x = new_value - old_value;
            if minor_iteration_details {
                let quadratic_objective = get_quadratic_objective(idata);
                icrash_util::print_minor_iteration_details(
                    k as f64,
                    col as f64,
                    idata.xk.col_value[col as usize] - delta_x,
                    delta_x,
                    objective_ica,
                    &residual_ica,
                    quadratic_objective,
                    options.log_options.clone(),
                );
            }

            debug_assert!(
                (objective_ica
                    - icrash_util::vector_product(&idata.lp.col_cost_, &idata.xk.col_value))
                .abs()
                    < 1e08
            );
        }

        // code below just for checking. Can comment out later if speed up is
        // needed.
        let mut residual_ica_check = vec![0.0; idata.lp.num_row_ as usize];
        icrash_util::update_residual_ica(&idata.lp, &idata.xk, &mut residual_ica_check);
        let difference = crate::util::highs_utils::get_norm2(&residual_ica)
            - crate::util::highs_utils::get_norm2(&residual_ica_check);
        debug_assert!(difference.abs() < 1e08);
        let _ = difference;
    }
}

/// HiGHS: ICrash.cpp::solveSubproblemQP (lines 286-316)
fn solve_subproblem_qp(idata: &mut Quadratic, options: &ICrashOptions) {
    let minor_iteration_details = false;

    let _ =
        crate::lp_data::presolve_support::calculate_row_values_quad(&idata.lp, &mut idata.xk, -1);
    let mut residual = vec![0.0; idata.lp.num_row_ as usize];
    icrash_util::update_residual_fast(&idata.lp, &idata.xk, &mut residual);
    let mut objective = 0.0;

    // todo: Ax = rv
    for k in 0..options.approximate_minimization_iterations {
        for col in 0..idata.lp.num_col_ {
            // determine whether to minimize for col.
            // if empty skip.
            if idata.lp.a_matrix_.start_[col as usize]
                == idata.lp.a_matrix_.start_[(col + 1) as usize]
            {
                continue;
            }

            let old_value = idata.xk.col_value[col as usize];
            icrash_util::minimize_component_qp(
                col,
                idata.mu,
                &idata.lp,
                &mut objective,
                &mut residual,
                &mut idata.xk,
            );

            let new_value = idata.xk.col_value[col as usize];
            let delta_x = new_value - old_value;
            if minor_iteration_details {
                let quadratic_objective = get_quadratic_objective(idata);
                icrash_util::print_minor_iteration_details(
                    k as f64,
                    col as f64,
                    idata.xk.col_value[col as usize] - delta_x,
                    delta_x,
                    objective,
                    &residual,
                    quadratic_objective,
                    options.log_options.clone(),
                );
            }
        }
    }
}

/// HiGHS: ICrash.cpp::solveSubproblem (lines 318-341)
pub fn solve_subproblem(idata: &mut Quadratic, options: &ICrashOptions) -> bool {
    match options.strategy {
        ICrashStrategy::kUpdatePenalty | ICrashStrategy::kUpdateAdmm | ICrashStrategy::kICA => {
            debug_assert!(!options.exact);
            solve_subproblem_ica(idata, options);
        }
        ICrashStrategy::kPenalty => {
            debug_assert!(!options.exact);
            solve_subproblem_qp(idata, options);
        }
        _ => {
            // std::cout <<               "ICrashError: Not implemented yet." << std::endl;
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "ICrashError: not implemented yet.\n"
            );
            return false;
        }
    }
    true
}

/// HiGHS: ICrash.cpp::reportSubproblem (lines 343-360)
pub fn report_subproblem(options: ICrashOptions, idata: &Quadratic, iteration: i32) {
    let mut ss = String::new();
    // Report outcome.
    if iteration == 0 {
        // ss << "Iteration " << std::setw(3) << 0 << ": objective " << std::setw(3)
        //    << std::fixed << std::setprecision(2) << idata.lp_objective
        //    << " residual " << std::setw(5) << std::scientific
        //    << idata.residual_norm_2 << std::endl;
        ss.push_str("Iteration ");
        ss.push_str(&c_printf::c_pad_right("0", 3));
        ss.push_str(": objective ");
        ss.push_str(&c_printf::c_fmt_f_width(idata.lp_objective, 2, 3));
        ss.push_str(" residual ");
        // setprecision(2) is still in effect: `%5.2e`.
        ss.push_str(&c_printf::c_pad_right(
            &c_printf::c_fmt_e(idata.residual_norm_2, 2),
            5,
        ));
        ss.push('\n');
    } else {
        // ss << "Iter " << std::setw(3) << iteration << ", mu " << idata.mu
        //    << std::scientific << ", c'x " << std::setprecision(5)
        //    << idata.lp_objective << ", res " << idata.residual_norm_2
        //    << ", quad_obj " << idata.quadratic_objective << std::endl;
        ss.push_str("Iter ");
        ss.push_str(&c_printf::c_pad_right(&iteration.to_string(), 3));
        ss.push_str(", mu ");
        ss.push_str(&c_printf::c_fmt_g(idata.mu, 6));
        ss.push_str(", c'x ");
        ss.push_str(&c_printf::c_fmt_e(idata.lp_objective, 5));
        ss.push_str(", res ");
        ss.push_str(&c_printf::c_fmt_e(idata.residual_norm_2, 5));
        ss.push_str(", quad_obj ");
        ss.push_str(&c_printf::c_fmt_e(idata.quadratic_objective, 5));
        ss.push('\n');
    }
    // std::cout << ss.str();
    crate::highs_log_user!(&options.log_options, HighsLogType::kInfo, "{}", ss);
}

/// HiGHS: ICrash.cpp::ICrashtrategyToString (lines 362-376)
pub fn icrashtrategy_to_string(strategy: ICrashStrategy) -> String {
    match strategy {
        ICrashStrategy::kPenalty => "Penalty".to_string(),
        ICrashStrategy::kAdmm => "ADMM".to_string(),
        ICrashStrategy::kICA => "ICA".to_string(),
        ICrashStrategy::kUpdatePenalty => "UpdatePenalty".to_string(),
        ICrashStrategy::kUpdateAdmm => "UpdateAdmm".to_string(),
    }
}

/// HiGHS: ICrash.cpp::reportOptions (lines 378-397)
pub fn report_options(options: &ICrashOptions) {
    let mut ss = String::new();
    // Report outcome.
    ss.push_str("ICrashOptions \n");
    ss.push_str("dualize: ");
    ss.push_str(if options.dualize { "true" } else { "false" });
    ss.push('\n');
    ss.push_str("strategy: ");
    ss.push_str(&icrashtrategy_to_string(options.strategy));
    ss.push('\n');
    ss.push_str("starting_weight: ");
    ss.push_str(&c_printf::c_fmt_e(options.starting_weight, 6));
    ss.push('\n');
    ss.push_str("iterations: ");
    ss.push_str(&options.iterations.to_string());
    ss.push('\n');
    if !options.exact {
        ss.push_str("approximate_minimization_iterations: ");
        ss.push_str(&options.approximate_minimization_iterations.to_string());
        ss.push('\n');
        ss.push_str("breakpoints: ");
        ss.push_str(if options.breakpoints { "true" } else { "false" });
        ss.push('\n');
    } else {
        ss.push_str("exact: true\n");
    }
    ss.push('\n');
    // std::cout << ss.str();
    crate::highs_log_user!(&options.log_options, HighsLogType::kInfo, "{}", ss);
}

/// HiGHS: ICrash.cpp::callICrash (lines 399-474)
pub fn call_icrash(lp: &HighsLp, options: &ICrashOptions, result: &mut ICrashInfo) -> HighsStatus {
    if !check_options(lp, options.clone()) {
        return HighsStatus::kError;
    }
    debug_assert!(lp.a_matrix_.is_colwise());

    // Initialize data structures and initial values.
    let mut idata = parse_options(lp, options.clone());
    report_options(options);
    initialize(&mut idata, options);
    update(&mut idata);
    report_subproblem(options.clone(), &idata, 0);
    idata.details.push(fill_details(0, &idata));

    let residual_0 = idata.residual_norm_2;
    // Initialize clocks.
    let start = std::time::Instant::now();

    // Main loop.
    let mut iteration: i32;
    iteration = 1;
    while iteration <= options.iterations {
        update_parameters(&mut idata, options, iteration);

        // Solve subproblem.
        let start_iteration = std::time::Instant::now();
        let success = solve_subproblem(&mut idata, options);
        if !success {
            return HighsStatus::kError;
        }
        let end_iteration = std::time::Instant::now();
        let elapsed_seconds = (end_iteration - start_iteration).as_secs_f64();

        update(&mut idata);

        report_subproblem(options.clone(), &idata, iteration);
        idata.details.push(fill_details(iteration, &idata));
        debug_assert!(iteration + 1 == idata.details.len() as i32);
        idata.details[iteration as usize].time = elapsed_seconds;

        // Exit if feasible.
        if idata.residual_norm_2 < kExitTolerance {
            // std::cout << "Solution feasible within exit tolerance: " <<
            // kExitTolerance << std::endl;
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "Solution feasible within exit tolerance: {}\n",
                c_printf::c_fmt_g(kExitTolerance, 6)
            );
            iteration += 1;
            break;
        }

        // Exit if residual much larger
        if idata.residual_norm_2 > (5.0 * residual_0) {
            crate::highs_log_user!(
                &options.log_options,
                HighsLogType::kInfo,
                "Residual growing too large: exit iCrash\n"
            );
            return HighsStatus::kError;
        }

        iteration += 1;
    }
    // Fill in return values.
    iteration -= 1;
    result.details = std::mem::take(&mut idata.details);
    // reportICrashIterationDetails(result.details);
    fill_icrash_info(iteration, result);
    result.x_values = idata.xk.col_value.clone();

    let end = std::time::Instant::now();
    let elapsed_seconds = (end - start).as_secs_f64();
    result.total_time = elapsed_seconds;

    // std::cout << "\nICrash finished successfully after: " << result.total_time
    // << "sec." << std::endl;
    crate::highs_log_user!(
        &options.log_options,
        HighsLogType::kInfo,
        "\nICrash finished successfully after: {} sec.\n",
        c_printf::c_fmt_g(result.total_time, 3)
    );

    HighsStatus::kOk
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

    #[test]
    fn strategy_parse_and_print() {
        let mut s = ICrashStrategy::kPenalty;
        assert!(parse_icrash_strategy(" ICA ", &mut s));
        assert_eq!(s, ICrashStrategy::kICA);
        assert!(parse_icrash_strategy("update_admm", &mut s));
        assert_eq!(s, ICrashStrategy::kUpdateAdmm);
        assert!(!parse_icrash_strategy("foo", &mut s));
        assert_eq!(s, ICrashStrategy::kUpdateAdmm);
        assert_eq!(
            icrashtrategy_to_string(ICrashStrategy::kUpdatePenalty),
            "UpdatePenalty"
        );
    }

    #[test]
    fn exact_and_breakpoints_are_rejected() {
        let lp = HighsLp::default();
        let base = ICrashOptions {
            dualize: false,
            strategy: ICrashStrategy::kICA,
            starting_weight: 1e-3,
            iterations: 30,
            approximate_minimization_iterations: 50,
            exact: false,
            breakpoints: false,
            log_options: HighsLogOptions::default(),
        };
        assert!(check_options(&lp, base.clone()));
        let mut o = base.clone();
        o.exact = true;
        assert!(!check_options(&lp, o));
        let mut o = base;
        o.breakpoints = true;
        assert!(!check_options(&lp, o));
    }

    fn quiet_options(strategy: ICrashStrategy) -> ICrashOptions {
        let mut log_options = HighsLogOptions::default();
        log_options.output_flag = false;
        ICrashOptions {
            dualize: false,
            strategy,
            starting_weight: 1e-3,
            iterations: 30,
            approximate_minimization_iterations: 50,
            exact: false,
            breakpoints: false,
            log_options,
        }
    }

    /// min x0 + x1 s.t. x0 + x1 = 1, 0 <= x <= 10 (an equality problem).
    fn tiny_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.row_lower_ = vec![1.0];
        lp.row_upper_ = vec![1.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        lp
    }

    /// ICrash.cpp:195-226: no update on iteration 1; kPenalty scales mu by 0.1
    /// every later iteration, kUpdatePenalty only when iteration % 3 == 0, kAdmm
    /// never.
    #[test]
    fn update_parameters_weight_schedule() {
        let mut idata = Quadratic::new(tiny_lp(), quiet_options(ICrashStrategy::kPenalty));
        idata.mu = 1.0;
        update_parameters(&mut idata, &quiet_options(ICrashStrategy::kPenalty), 1);
        assert_eq!(idata.mu, 1.0);
        update_parameters(&mut idata, &quiet_options(ICrashStrategy::kPenalty), 2);
        assert_eq!(idata.mu, 0.1 * 1.0);
        idata.mu = 1.0;
        update_parameters(
            &mut idata,
            &quiet_options(ICrashStrategy::kUpdatePenalty),
            4,
        );
        assert_eq!(idata.mu, 1.0);
        update_parameters(
            &mut idata,
            &quiet_options(ICrashStrategy::kUpdatePenalty),
            3,
        );
        assert_eq!(idata.mu, 0.1);
        idata.mu = 1.0;
        update_parameters(&mut idata, &quiet_options(ICrashStrategy::kAdmm), 3);
        assert_eq!(idata.mu, 1.0);
    }

    /// ICrash.cpp:130-142: the "1/2mu r'r" term uses residual'lambda (kept as in C++).
    #[test]
    fn quadratic_objective_uses_residual_lambda() {
        let mut idata = Quadratic::new(tiny_lp(), quiet_options(ICrashStrategy::kPenalty));
        idata.xk.col_value = vec![3.0, 4.0];
        idata.lambda = vec![1.0];
        idata.residual = vec![2.0];
        idata.mu = 0.5;
        // 7 + 2 + 2 / 1
        assert_eq!(get_quadratic_objective(&idata), 11.0);
    }

    /// ICrash.cpp:332-338: kAdmm has no subproblem solver.
    #[test]
    fn admm_subproblem_is_rejected() {
        let options = quiet_options(ICrashStrategy::kAdmm);
        let mut idata = Quadratic::new(tiny_lp(), options.clone());
        assert!(!solve_subproblem(&mut idata, &options));
        let mut result = ICrashInfo::default();
        assert_eq!(
            call_icrash(&tiny_lp(), &options, &mut result),
            HighsStatus::kError
        );
    }

    /// ICrash.cpp:180-191 and 399-474: details hold iterations 0..=num_iterations
    /// and the weights follow the kPenalty schedule.
    #[test]
    fn call_icrash_penalty_details() {
        let options = quiet_options(ICrashStrategy::kPenalty);
        let mut result = ICrashInfo::default();
        let status = call_icrash(&tiny_lp(), &options, &mut result);
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(result.details.len() as i32, result.num_iterations + 1);
        assert_eq!(result.x_values.len(), 2);
        let mut mu = options.starting_weight;
        for (i, d) in result.details.iter().enumerate() {
            assert_eq!(d.num, i as i32);
            if i >= 2 {
                mu = 0.1 * mu;
            }
            assert_eq!(d.weight, mu);
        }
        assert_eq!(result.starting_weight, options.starting_weight);
        let last = &result.details[result.num_iterations as usize];
        assert_eq!(result.final_weight, last.weight);
        assert_eq!(result.final_residual_norm_2, last.residual_norm_2);
        assert_eq!(result.final_lp_objective, last.lp_objective);
    }
}
