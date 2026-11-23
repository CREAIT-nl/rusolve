//! HiGHS source: lp_data/HighsSolution.cpp (ipxSolutionToHighsSolution,
//! ipxBasicSolutionToHighsBasicSolution)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The two IPX -> HiGHS solution copies, kept apart from `highs_solution.rs`
//! so the IPX port owns them. `constraint_type` is C++ `std::vector<char>`.

use crate::ipm::ipx::{IPX_basic, IPX_nonbasic_lb, IPX_nonbasic_ub, IPX_superbasic};
use crate::ipm::ipx_wrapper::IpxSolution;
use crate::prelude::*;
use crate::util::c_printf::c_fmt_g_width;

/// HiGHS: HighsSolution.cpp::ipxSolutionToHighsSolution (lines 1500-1588)
/// HiGHS 1.16 (#3319): HighsSolution.cpp::ipxSolutionToHighsSolution (lines 1625-1713)
///
/// `options` is unused by the C++ function too.
#[allow(clippy::needless_range_loop)]
pub fn ipx_solution_to_highs_solution(
    options: &HighsOptions,
    lp: &HighsLp,
    rhs: &[f64],
    constraint_type: &[u8],
    ipx_num_col: HighsInt,
    ipx_num_row: HighsInt,
    ipx_x: &[f64],
    ipx_slack_vars: &[f64],
    ipx_y: &[f64],
    ipx_zl: &[f64],
    ipx_zu: &[f64],
    highs_solution: &mut HighsSolution,
) -> HighsStatus {
    let _ = options;
    let num_col = lp.num_col_ as usize;
    let num_row = lp.num_row_ as usize;

    // Resize the HighsSolution
    highs_solution.col_value.resize(num_col, 0.0);
    highs_solution.row_value.resize(num_row, 0.0);
    highs_solution.col_dual.resize(num_col, 0.0);
    highs_solution.row_dual.resize(num_row, 0.0);

    let ipx_col_value = ipx_x;
    let ipx_row_value = ipx_slack_vars;

    // Row activities are needed to set activity values of free rows -
    // which are ignored by IPX
    let get_row_activities = ipx_num_row < lp.num_row_;
    let mut row_activity: Vec<f64> = if get_row_activities {
        vec![0.0; num_row]
    } else {
        Vec::new()
    };
    // HiGHS 1.16 #3319 (h16_ipx_free_rows): HighsSolution.cpp line 1647 (1.16.0), the
    // assert relaxed to `<=` on the flag path (rule 9: the off path keeps `==`).
    if !options.h16_ipx_free_rows {
        debug_assert!(ipx_num_row == lp.num_row_);
    } else {
        debug_assert!(ipx_num_row <= lp.num_row_);
    }
    for i_col in 0..num_col {
        let value = ipx_col_value[i_col];
        if get_row_activities {
            // Accumulate row activities to assign value to free rows
            let start = lp.a_matrix_.start_[i_col] as usize;
            let end = lp.a_matrix_.start_[i_col + 1] as usize;
            for el in start..end {
                let row = lp.a_matrix_.index_[el] as usize;
                row_activity[row] += value * lp.a_matrix_.value_[el];
            }
        }
        let dual = ipx_zl[i_col] - ipx_zu[i_col];
        highs_solution.col_value[i_col] = value;
        highs_solution.col_dual[i_col] = dual;
    }
    let mut ipx_row: HighsInt = 0;
    let mut ipx_slack: HighsInt = lp.num_col_;
    for i_row in 0..num_row {
        let lower = lp.row_lower_[i_row];
        let upper = lp.row_upper_[i_row];
        if lower <= -kHighsInf && upper >= kHighsInf {
            // Free row - removed by IPX so set it to its row activity
            highs_solution.row_value[i_row] = row_activity[i_row];
            highs_solution.row_dual[i_row] = 0.0;
            continue;
        }
        // Non-free row, so IPX will have it
        let (value, dual) = if (lower > -kHighsInf && upper < kHighsInf) && (lower < upper) {
            debug_assert!(constraint_type[ipx_row as usize] == b'=');
            // Boxed row - look at its slack
            let value = ipx_col_value[ipx_slack as usize];
            let dual = ipx_zl[ipx_slack as usize] - ipx_zu[ipx_slack as usize];
            // Update the slack to be used for boxed rows
            ipx_slack += 1;
            (value, dual)
        } else {
            let value = rhs[ipx_row as usize] - ipx_row_value[ipx_row as usize];
            let dual = ipx_y[ipx_row as usize];
            (value, dual)
        };
        highs_solution.row_value[i_row] = value;
        highs_solution.row_dual[i_row] = dual;
        // Update the IPX row index
        ipx_row += 1;
    }
    debug_assert!(ipx_row == ipx_num_row);
    debug_assert!(ipx_slack == ipx_num_col);
    if lp.sense_ == ObjSense::kMaximize {
        // Flip dual values since original LP is maximization
        for i_col in 0..num_col {
            highs_solution.col_dual[i_col] *= -1.0;
        }
        for i_row in 0..num_row {
            highs_solution.row_dual[i_row] *= -1.0;
        }
    }

    // Indicate that the primal and dual solution are known
    highs_solution.value_valid = true;
    highs_solution.dual_valid = true;
    HighsStatus::kOk
}

/// HiGHS: HighsSolution.cpp::ipxBasicSolutionToHighsBasicSolution (lines 1590-1841)
#[allow(clippy::needless_range_loop)]
pub fn ipx_basic_solution_to_highs_basic_solution(
    log_options: &HighsLogOptions,
    lp: &HighsLp,
    rhs: &[f64],
    constraint_type: &[u8],
    ipx_solution: &IpxSolution,
    highs_basis: &mut HighsBasis,
    highs_solution: &mut HighsSolution,
) -> HighsStatus {
    let num_col = lp.num_col_ as usize;
    let num_row = lp.num_row_ as usize;

    // Resize the HighsSolution and HighsBasis
    highs_solution.col_value.resize(num_col, 0.0);
    highs_solution.row_value.resize(num_row, 0.0);
    highs_solution.col_dual.resize(num_col, 0.0);
    highs_solution.row_dual.resize(num_row, 0.0);
    highs_basis
        .col_status
        .resize(num_col, HighsBasisStatus::kLower);
    highs_basis
        .row_status
        .resize(num_row, HighsBasisStatus::kLower);

    let ipx_col_value = &ipx_solution.ipx_col_value;
    let ipx_row_value = &ipx_solution.ipx_row_value;
    let ipx_col_dual = &ipx_solution.ipx_col_dual;
    let ipx_row_dual = &ipx_solution.ipx_row_dual;
    let ipx_col_status = &ipx_solution.ipx_col_status;
    let ipx_row_status = &ipx_solution.ipx_row_status;

    // IPX's status constants (ipx::Int values). The C++ locals are named
    // ipx_basic, ipx_nonbasic_at_lb, ipx_nonbasic_at_ub and ipx_superbasic.
    // The IPX constants have the same values, so they are used directly.

    // Row activities are needed to set activity values of free rows -
    // which are ignored by IPX
    let get_row_activities = ipx_solution.num_row < lp.num_row_;
    let mut row_activity: Vec<f64> = if get_row_activities {
        vec![0.0; num_row]
    } else {
        Vec::new()
    };
    let mut num_basic_variables: HighsInt = 0;
    for col in 0..num_col {
        let mut unrecognised = false;
        if ipx_col_status[col] == IPX_basic {
            // Column is basic
            highs_basis.col_status[col] = HighsBasisStatus::kBasic;
            highs_solution.col_value[col] = ipx_col_value[col];
            highs_solution.col_dual[col] = 0.0;
        } else if ipx_col_status[col] == IPX_nonbasic_lb {
            // Column is at lower bound
            highs_basis.col_status[col] = HighsBasisStatus::kLower;
            highs_solution.col_value[col] = ipx_col_value[col];
            highs_solution.col_dual[col] = ipx_col_dual[col];
        } else if ipx_col_status[col] == IPX_nonbasic_ub {
            // Column is at upper bound
            highs_basis.col_status[col] = HighsBasisStatus::kUpper;
            highs_solution.col_value[col] = ipx_col_value[col];
            highs_solution.col_dual[col] = ipx_col_dual[col];
        } else if ipx_col_status[col] == IPX_superbasic {
            // Column is superbasic
            highs_basis.col_status[col] = HighsBasisStatus::kZero;
            highs_solution.col_value[col] = ipx_col_value[col];
            highs_solution.col_dual[col] = ipx_col_dual[col];
        } else {
            unrecognised = true;
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kError,
                "\nError in IPX conversion: Unrecognised value ipx_col_status[{:>2}] = {}\n",
                col,
                ipx_col_status[col]
            );
        }
        if unrecognised {
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kError,
                "Bounds [{}, {}]\n",
                c_fmt_g_width(lp.col_lower_[col], 4, 11),
                c_fmt_g_width(lp.col_upper_[col], 4, 11)
            );
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kError,
                "Col {:>2} ipx_col_status[{:>2}] = {:>2}; x[{:>2}] = {}; z[{:>2}] = {}\n",
                col,
                col,
                ipx_col_status[col],
                col,
                c_fmt_g_width(ipx_col_value[col], 4, 11),
                col,
                c_fmt_g_width(ipx_col_dual[col], 4, 11)
            );
            debug_assert!(!unrecognised);
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Unrecognised ipx_col_status value from IPX\n"
            );
            return HighsStatus::kError;
        }
        if get_row_activities {
            // Accumulate row activities to assign value to free rows
            let start = lp.a_matrix_.start_[col] as usize;
            let end = lp.a_matrix_.start_[col + 1] as usize;
            for el in start..end {
                let row = lp.a_matrix_.index_[el] as usize;
                row_activity[row] += highs_solution.col_value[col] * lp.a_matrix_.value_[el];
            }
        }
        if highs_basis.col_status[col] == HighsBasisStatus::kBasic {
            num_basic_variables += 1;
        }
    }
    let mut ipx_row: HighsInt = 0;
    let mut ipx_slack: HighsInt = lp.num_col_;
    let mut num_boxed_rows: HighsInt = 0;
    let mut num_boxed_rows_basic: HighsInt = 0;
    let mut num_boxed_row_slacks_basic: HighsInt = 0;
    for row in 0..num_row {
        let mut unrecognised = false;
        let lower = lp.row_lower_[row];
        let upper = lp.row_upper_[row];
        let this_ipx_row = ipx_row;
        if lower <= -kHighsInf && upper >= kHighsInf {
            // Free row - removed by IPX so make it basic at its row activity
            highs_basis.row_status[row] = HighsBasisStatus::kBasic;
            highs_solution.row_value[row] = row_activity[row];
            highs_solution.row_dual[row] = 0.0;
        } else {
            // Non-free row, so IPX will have it
            if (lower > -kHighsInf && upper < kHighsInf) && (lower < upper) {
                // Boxed row - look at its slack
                num_boxed_rows += 1;
                let slack_value = ipx_col_value[ipx_slack as usize];
                let slack_dual = ipx_col_dual[ipx_slack as usize];
                let value = slack_value;
                // @FlipRowDual -slack_dual became slack_dual
                let dual = slack_dual;
                if ipx_row_status[ipx_row as usize] == IPX_basic {
                    // Row is basic
                    num_boxed_rows_basic += 1;
                    highs_basis.row_status[row] = HighsBasisStatus::kBasic;
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = 0.0;
                } else if ipx_col_status[ipx_slack as usize] == IPX_basic {
                    // Slack is basic
                    num_boxed_row_slacks_basic += 1;
                    highs_basis.row_status[row] = HighsBasisStatus::kBasic;
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = 0.0;
                } else if ipx_col_status[ipx_slack as usize] == IPX_nonbasic_lb {
                    // Slack at lower bound
                    highs_basis.row_status[row] = HighsBasisStatus::kLower;
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = dual;
                } else if ipx_col_status[ipx_slack as usize] == IPX_nonbasic_ub {
                    // Slack is at its upper bound
                    debug_assert!(ipx_col_status[ipx_slack as usize] == IPX_nonbasic_ub);
                    highs_basis.row_status[row] = HighsBasisStatus::kUpper;
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = dual;
                } else {
                    unrecognised = true;
                    crate::highs_log_dev!(
                        log_options,
                        HighsLogType::kError,
                        "Error in IPX conversion: Row {:>2} (IPX row {:>2}) has unrecognised value ipx_col_status[{:>2}] = {}\n",
                        row,
                        ipx_row,
                        ipx_slack,
                        ipx_col_status[ipx_slack as usize]
                    );
                }
                // Update the slack to be used for boxed rows
                ipx_slack += 1;
            } else if ipx_row_status[ipx_row as usize] == IPX_basic {
                // Row is basic
                highs_basis.row_status[row] = HighsBasisStatus::kBasic;
                highs_solution.row_value[row] =
                    rhs[ipx_row as usize] - ipx_row_value[ipx_row as usize];
                highs_solution.row_dual[row] = 0.0;
            } else {
                // Nonbasic row at fixed value, lower bound or upper bound
                // (ipx_nonbasic_row = -1)
                debug_assert!(ipx_row_status[ipx_row as usize] == -1);
                let value = rhs[ipx_row as usize] - ipx_row_value[ipx_row as usize];
                // @FlipRowDual -ipx_row_dual[ipx_row]; became ipx_row_dual[ipx_row];
                let dual = ipx_row_dual[ipx_row as usize];
                // IPX quirk: the C++ prints constraint_type with %d, and char is
                // signed on the Apple arm64 target, so the byte is sign-extended.
                if constraint_type[ipx_row as usize] == b'>' {
                    // Row is at its lower bound
                    highs_basis.row_status[row] = HighsBasisStatus::kLower;
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = dual;
                } else if constraint_type[ipx_row as usize] == b'<' {
                    // Row is at its upper bound
                    highs_basis.row_status[row] = HighsBasisStatus::kUpper;
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = dual;
                } else if constraint_type[ipx_row as usize] == b'=' {
                    // Row is at its fixed value: set HighsBasisStatus according
                    // to sign of dual.
                    //
                    // Don't worry about maximization problems. IPX solves them
                    // as minimizations with negated costs, so a negative dual
                    // yields HighsBasisStatus::kUpper here, and dual signs are
                    // then flipped below, so HighsBasisStatus::kUpper will have
                    // corresponding positive dual.
                    highs_basis.row_status[row] = if dual >= 0.0 {
                        HighsBasisStatus::kLower
                    } else {
                        HighsBasisStatus::kUpper
                    };
                    highs_solution.row_value[row] = value;
                    highs_solution.row_dual[row] = dual;
                } else {
                    unrecognised = true;
                    crate::highs_log_dev!(
                        log_options,
                        HighsLogType::kError,
                        "Error in IPX conversion: Row {:>2}: cannot handle constraint_type[{:>2}] = {}\n",
                        row,
                        ipx_row,
                        constraint_type[ipx_row as usize] as i8 as i32
                    );
                }
            }
            // Update the IPX row index
            ipx_row += 1;
        }
        if unrecognised {
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kError,
                "Bounds [{}, {}]\n",
                c_fmt_g_width(lp.row_lower_[row], 4, 11),
                c_fmt_g_width(lp.row_upper_[row], 4, 11)
            );
            crate::highs_log_dev!(
                log_options,
                HighsLogType::kError,
                "Row {:>2} ipx_row_status[{:>2}] = {:>2}; s[{:>2}] = {}; y[{:>2}] = {}\n",
                row,
                this_ipx_row,
                ipx_row_status[this_ipx_row as usize],
                this_ipx_row,
                c_fmt_g_width(ipx_row_value[this_ipx_row as usize], 4, 11),
                this_ipx_row,
                c_fmt_g_width(ipx_row_dual[this_ipx_row as usize], 4, 11)
            );
            debug_assert!(!unrecognised);
            crate::highs_log_user!(
                log_options,
                HighsLogType::kError,
                "Unrecognised ipx_row_status value from IPX\n"
            );
            return HighsStatus::kError;
        }
        if highs_basis.row_status[row] == HighsBasisStatus::kBasic {
            num_basic_variables += 1;
        }
    }
    debug_assert!(num_basic_variables == lp.num_row_);
    debug_assert!(ipx_row == ipx_solution.num_row);
    debug_assert!(ipx_slack == ipx_solution.num_col);

    if lp.sense_ == ObjSense::kMaximize {
        // Flip dual values since original LP is maximization
        for i_col in 0..num_col {
            highs_solution.col_dual[i_col] *= -1.0;
        }
        for i_row in 0..num_row {
            highs_solution.row_dual[i_row] *= -1.0;
        }
    }

    if num_boxed_rows != 0 {
        crate::highs_log_dev!(
            log_options,
            HighsLogType::kInfo,
            "Of {} boxed rows: {} are basic and {} have basic slacks\n",
            num_boxed_rows,
            num_boxed_rows_basic,
            num_boxed_row_slacks_basic
        );
    }
    // Indicate that the primal solution, dual solution and basis are valid
    highs_solution.value_valid = true;
    highs_solution.dual_valid = true;
    highs_basis.valid = true;
    highs_basis.useful = true;
    HighsStatus::kOk
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two structural columns, four rows: row 0 free, row 1 boxed [1, 3],
    /// row 2 '<' with upper bound 5, row 3 '=' fixed at 2.
    fn tiny_lp(sense: ObjSense) -> HighsLp {
        let mut a = HighsSparseMatrix::default();
        a.num_col_ = 2;
        a.num_row_ = 4;
        a.start_ = vec![0, 2, 4];
        a.index_ = vec![0, 1, 0, 3];
        a.value_ = vec![1.0, 2.0, -1.0, 1.0];
        HighsLp {
            num_col_: 2,
            num_row_: 4,
            col_lower_: vec![0.0, 0.0],
            col_upper_: vec![10.0, 10.0],
            row_lower_: vec![-kHighsInf, 1.0, -kHighsInf, 2.0],
            row_upper_: vec![kHighsInf, 3.0, 5.0, 2.0],
            a_matrix_: a,
            sense_: sense,
            ..Default::default()
        }
    }

    /// IPX sees 2 structural columns plus 1 slack (for the boxed row) and
    /// 3 non-free rows (boxed, '<', '=').
    const RHS: [f64; 3] = [4.0, 3.0, 2.0];
    const CONSTRAINT_TYPE: [u8; 3] = *b"=<=";

    #[test]
    fn ipx_solution_maps_free_boxed_and_inequality_rows() {
        let x = [1.0, 2.0, 2.5];
        let slack = [0.5, 1.0, 0.25];
        let y = [0.0, 0.75, -0.5];
        let zl = [0.5, 0.0, 0.75];
        let zu = [0.0, 0.25, 0.5];

        let run = |sense: ObjSense| {
            let lp = tiny_lp(sense);
            let mut sol = HighsSolution::default();
            // IPX drops free row 0 (3 IPX rows for 4 LP rows): the #3319 case, so the
            // flag path's relaxed `ipx_num_row <= num_row` assert applies.
            let options = HighsOptions {
                h16_ipx_free_rows: true,
                ..HighsOptions::default()
            };
            let status = ipx_solution_to_highs_solution(
                &options,
                &lp,
                &RHS,
                &CONSTRAINT_TYPE,
                3,
                3,
                &x,
                &slack,
                &y,
                &zl,
                &zu,
                &mut sol,
            );
            assert_eq!(status, HighsStatus::kOk);
            assert!(sol.value_valid && sol.dual_valid);
            sol
        };

        // Minimisation: duals are taken as IPX reports them.
        let sol = run(ObjSense::kMinimize);
        assert_eq!(sol.col_value, vec![1.0, 2.0]);
        assert_eq!(sol.col_dual, vec![0.5, -0.25]);
        // Free row 0 gets its activity 1*1 + 2*(-1) = -1.
        // Boxed row 1 takes its slack column (value 2.5, dual zl-zu = 0.25).
        // '<' row 2 uses IPX row 1 (3 - 1 = 2, dual 0.75).
        // '=' row 3 uses IPX row 2 (2 - 0.25 = 1.75, dual -0.5).
        assert_eq!(sol.row_value, vec![-1.0, 2.5, 2.0, 1.75]);
        assert_eq!(sol.row_dual, vec![0.0, 0.25, 0.75, -0.5]);

        // Maximisation: both dual vectors are negated.
        let sol = run(ObjSense::kMaximize);
        assert_eq!(sol.col_dual, vec![-0.5, 0.25]);
        assert_eq!(sol.row_dual, vec![-0.0, -0.25, -0.75, 0.5]);
        assert_eq!(sol.row_value, vec![-1.0, 2.5, 2.0, 1.75]);
    }

    #[test]
    fn ipx_basic_solution_sets_statuses_and_flips_max_duals() {
        // Columns 0 and 1 basic, boxed row's slack basic, rows nonbasic.
        let ipx = IpxSolution {
            num_col: 3,
            num_row: 3,
            ipx_col_value: vec![1.0, 2.0, 2.5],
            ipx_row_value: vec![0.5, 1.0, 0.25],
            ipx_col_dual: vec![0.5, -0.25, 0.25],
            ipx_row_dual: vec![0.0, 0.75, 0.5],
            ipx_col_status: vec![IPX_basic, IPX_basic, IPX_basic],
            ipx_row_status: vec![-1, -1, -1],
        };

        let run = |sense: ObjSense| {
            let lp = tiny_lp(sense);
            let mut basis = HighsBasis::default();
            let mut sol = HighsSolution::default();
            let status = ipx_basic_solution_to_highs_basic_solution(
                &HighsLogOptions::default(),
                &lp,
                &RHS,
                &CONSTRAINT_TYPE,
                &ipx,
                &mut basis,
                &mut sol,
            );
            assert_eq!(status, HighsStatus::kOk);
            (basis, sol)
        };

        let (basis, sol) = run(ObjSense::kMinimize);
        assert!(basis.valid && basis.useful);
        assert_eq!(
            basis.col_status,
            vec![HighsBasisStatus::kBasic, HighsBasisStatus::kBasic]
        );
        // Free row (basic via activity), boxed row (basic slack), '<' row
        // (nonbasic, kUpper), '=' row (nonbasic with dual 0.5 >= 0, kLower).
        // Basic count = 2 columns + 2 rows = num_row.
        assert_eq!(
            basis.row_status,
            vec![
                HighsBasisStatus::kBasic,
                HighsBasisStatus::kBasic,
                HighsBasisStatus::kUpper,
                HighsBasisStatus::kLower,
            ]
        );
        assert_eq!(sol.col_value, vec![1.0, 2.0]);
        assert_eq!(sol.col_dual, vec![0.0, 0.0]);
        assert_eq!(sol.row_value, vec![-1.0, 2.5, 2.0, 1.75]);
        assert_eq!(sol.row_dual, vec![0.0, 0.0, 0.75, 0.5]);

        let (_, sol) = run(ObjSense::kMaximize);
        assert_eq!(sol.row_dual, vec![-0.0, -0.0, -0.75, -0.5]);
    }
}
