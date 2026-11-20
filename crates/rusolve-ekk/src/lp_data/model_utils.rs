//! HiGHS source: lp_data/HighsModelUtils.cpp (utilModelStatusToString,
//! writeLpObjective, writeObjectiveValue, writePrimalSolution),
//! simplex/HSimplexReport.cpp (reportSimplexPhaseIterations)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use crate::prelude::*;
use crate::simplex::ekk::HEkk;

/// HiGHS: HighsModelUtils.cpp::utilModelStatusToString (lines 1415-1481, 67 lines)
pub fn util_model_status_to_string(model_status: HighsModelStatus) -> String {
    let s: &str = match model_status {
        HighsModelStatus::kNotset => "Not Set",
        HighsModelStatus::kLoadError => "Load error",
        HighsModelStatus::kModelError => "Model error",
        HighsModelStatus::kPresolveError => "Presolve error",
        HighsModelStatus::kSolveError => "Solve error",
        HighsModelStatus::kPostsolveError => "Postsolve error",
        HighsModelStatus::kModelEmpty => "Empty",
        HighsModelStatus::kMemoryLimit => "Memory limit reached",
        HighsModelStatus::kOptimal => "Optimal",
        HighsModelStatus::kInfeasible => "Infeasible",
        HighsModelStatus::kUnboundedOrInfeasible => "Primal infeasible or unbounded",
        HighsModelStatus::kUnbounded => "Unbounded",
        HighsModelStatus::kObjectiveBound => "Bound on objective reached",
        HighsModelStatus::kObjectiveTarget => "Target for objective reached",
        HighsModelStatus::kTimeLimit => "Time limit reached",
        HighsModelStatus::kIterationLimit => "Iteration limit reached",
        HighsModelStatus::kSolutionLimit => "Solution limit reached",
        HighsModelStatus::kInterrupt => "Interrupted by user",
        HighsModelStatus::kHighsInterrupt => "Interrupted by HiGHS",
        HighsModelStatus::kUnknown => "Unknown",
    };
    s.to_string()
}

/// HiGHS: HSimplexReport.cpp::reportSimplexPhaseIterations (lines 16-77, 62 lines)
pub fn report_simplex_phase_iterations(
    log_options: &HighsLogOptions,
    iteration_count: HighsInt,
    info: &mut HighsSimplexInfo,
    initialise: bool,
) {
    // Out of scope (reporting only): no-op.
}

/// HiGHS: HighsModelUtils.cpp::utilBasisStatusToString (lines 1382-1403, 22 lines)
pub fn util_basis_status_to_string(basis_status: HighsBasisStatus) -> String {
    let s = match basis_status {
        HighsBasisStatus::kLower => "At lower/fixed bound",
        HighsBasisStatus::kBasic => "Basic",
        HighsBasisStatus::kUpper => "At upper bound",
        HighsBasisStatus::kZero => "Free at zero",
        HighsBasisStatus::kNonbasic => "Nonbasic",
    };
    s.to_string()
}

/// HiGHS: HighsModelUtils.cpp::utilPresolveRuleTypeToString (lines 1483-1527, 45 lines)
pub fn util_presolve_rule_type_to_string(rule_type: HighsInt) -> String {
    let s = if rule_type == kPresolveRuleEmptyRow {
        "Empty row"
    } else if rule_type == kPresolveRuleSingletonRow {
        "Singleton row"
    } else if rule_type == kPresolveRuleRedundantRow {
        "Redundant row"
    } else if rule_type == kPresolveRuleEmptyCol {
        "Empty column"
    } else if rule_type == kPresolveRuleFixedCol {
        "Fixed column"
    } else if rule_type == kPresolveRuleDominatedCol {
        "Dominated col"
    } else if rule_type == kPresolveRuleForcingRow {
        "Forcing row"
    } else if rule_type == kPresolveRuleForcingCol {
        "Forcing col"
    } else if rule_type == kPresolveRuleFreeColSubstitution {
        "Free col substitution"
    } else if rule_type == kPresolveRuleDoubletonEquation {
        "Doubleton equation"
    } else if rule_type == kPresolveRuleDependentEquations {
        "Dependent equations"
    } else if rule_type == kPresolveRuleDependentFreeCols {
        "Dependent free columns"
    } else if rule_type == kPresolveRuleAggregator {
        "Aggregator"
    } else if rule_type == kPresolveRuleParallelRowsAndCols {
        "Parallel rows and columns"
    } else if rule_type == kPresolveRuleSparsify {
        "Sparsify"
    } else if rule_type == kPresolveRuleProbing {
        "Probing"
    } else if rule_type == kPresolveRuleEnumeration {
        "Enumeration"
    } else if rule_type == kPresolveRuleDualFixing {
        "Dual fixing"
    } else if rule_type == kPresolveRuleColStuffing {
        "Col stuffing"
    } else if rule_type == kPresolveRuleInitialSweep {
        "Initial sweep"
    } else {
        if cfg!(debug_assertions) {
            panic!("utilPresolveRuleTypeToString: illegal rule type");
        }
        "????"
    };
    s.to_string()
}

/// HiGHS: HighsModelUtils.cpp::writeLpObjective (lines 189-194).
/// `FILE* file` -> `&mut dyn std::io::Write`: the only
/// caller writes the MIP improving-solution file, never `stdout`, so the
/// `stdout` branch of `highsFprintfString` (a `highsLogUser`) is not reached
/// and `log_options` is unused. Write errors are ignored, as `fprintf`'s are.
pub fn write_lp_objective(
    file: &mut dyn std::io::Write,
    log_options: &HighsLogOptions,
    lp: &HighsLp,
    primal_solution: &[f64],
) {
    let objective_value = lp.objective_c_double_value(primal_solution);
    write_objective_value(file, log_options, f64::from(objective_value));
}

/// HiGHS: HighsModelUtils.cpp::writeObjectiveValue (lines 196-202)
pub fn write_objective_value(
    file: &mut dyn std::io::Write,
    log_options: &HighsLogOptions,
    objective_value: f64,
) {
    let obj_str = highs_double_to_string(objective_value, kHighsSolutionValueToStringTolerance);
    let s = format!("Objective {}\n", obj_str);
    let _ = file.write_all(s.as_bytes());
}

/// HiGHS: HighsModelUtils.cpp::writePrimalSolution (lines 204-244).
/// `highsFormatToString("%-s %s", ..)` is `name value`; `" %d"` appends the
/// column index when `sparse`. Ends with `fflush(file)` -> `flush`.
pub fn write_primal_solution(
    file: &mut dyn std::io::Write,
    log_options: &HighsLogOptions,
    lp: &HighsLp,
    primal_solution: &[f64],
    sparse: bool,
) {
    // Use when writing out the solution file (when names can be assumed
    // to exist) and the improving solution in the MIP solver (when
    // names cannot be assumed to exist)
    let mut num_nonzero_primal_value: HighsInt = 0;
    let have_col_names = !lp.col_names_.is_empty();
    if have_col_names {
        debug_assert!(lp.col_names_.len() == lp.num_col_ as usize);
    }
    if sparse {
        // Determine the number of nonzero primal solution values
        for i_col in 0..lp.num_col_ as usize {
            if primal_solution[i_col] != 0.0 {
                num_nonzero_primal_value += 1;
            }
        }
    }
    // Indicate the number of column values to be written out, depending
    // on whether format is sparse: either lp.num_col_ if not sparse, or
    // the negation of the number of nonzero values, if sparse
    let num_col_field = if sparse {
        -num_nonzero_primal_value
    } else {
        lp.num_col_
    };
    let _ = file.write_all(format!("# Columns {}\n", num_col_field).as_bytes());
    for ix in 0..lp.num_col_ as usize {
        if sparse && primal_solution[ix] == 0.0 {
            continue;
        }
        let val_str =
            highs_double_to_string(primal_solution[ix], kHighsSolutionValueToStringTolerance);
        // Don't invent names locallty: if none exist, then indicate this
        // - so that the (sparse) solution line format remains "name value
        // (index)"
        let name: &str = if have_col_names {
            &lp.col_names_[ix]
        } else {
            "NoName"
        };
        let mut ss = format!("{} {}", name, val_str);
        if sparse {
            ss += &format!(" {}", ix);
        }
        ss.push('\n');
        let _ = file.write_all(ss.as_bytes());
    }
    let _ = file.flush();
}

#[cfg(test)]
mod g1_writer_tests {
    //! The improving-solution writers.
    use super::*;

    fn lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.col_cost_ = vec![1.0, 2.0, -0.5];
        lp.offset_ = 0.25;
        lp
    }

    #[test]
    fn objective_and_dense_solution_text() {
        let lp = lp();
        let log_options = HighsLogOptions::default();
        let mut out: Vec<u8> = Vec::new();
        let x = [1.0, 0.0, 3.0];
        write_lp_objective(&mut out, &log_options, &lp, &x);
        write_primal_solution(&mut out, &log_options, &lp, &x, false);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "Objective -0.25\n# Columns 3\nNoName 1\nNoName 0\nNoName 3\n"
        );
    }

    #[test]
    fn sparse_solution_text_uses_names_and_indices() {
        let mut lp = lp();
        lp.col_names_ = vec!["a".into(), "b".into(), "c".into()];
        let log_options = HighsLogOptions::default();
        let mut out: Vec<u8> = Vec::new();
        write_primal_solution(&mut out, &log_options, &lp, &[0.0, 2.5, 0.0], true);
        assert_eq!(String::from_utf8(out).unwrap(), "# Columns -1\nb 2.5 1\n");
    }
}
