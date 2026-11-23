//! HiGHS source: lp_data/HighsSolution.cpp (formSimplexLpBasisAndFactorReturn,
//! formSimplexLpBasisAndFactor, accommodateAlienBasis)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! API slice. Reached from two places:
//! * `solveLpSimplex` (simplex/happ.rs, HApp.h:179-190) when there is no
//!   simplex basis and the HiGHS basis is useful but not valid. That is the
//!   state `deleteRows` leaves when a NONBASIC row is deleted
//!   (`deleteBasisRows` sets `valid = false`, `useful` stays true).
//! * `Highs::setBasis` with an alien basis (lp_data/highs_api.rs).
//!
//! Ownership (crate rules (a), (e)): the Rust `HighsLpSolverObject` has no
//! `ekk_instance_`, so the `HEkk` is a separate `ekk` parameter; every C++
//! `solver_object.ekk_instance_` / `ekk_instance` is `ekk`. `ekk.move_lp`
//! takes the LP out of `solver_object.lp_`; the `Return` helper moves it
//! back with `solver_object.lp_.move_back_lp_and_unapply_scaling(&mut ekk.lp_)`.

use crate::prelude::*;

/// HiGHS: HighsSolution.cpp::formSimplexLpBasisAndFactorReturn (lines 1842-1848, 7 lines)
pub fn form_simplex_lp_basis_and_factor_return(
    return_status: HighsStatus,
    ekk: &mut HEkk,
    solver_object: &mut HighsLpSolverObject,
) -> HighsStatus {
    // C++ `HighsLp& lp = solver_object.lp_;` and
    // `HighsLp& ekk_lp = solver_object.ekk_instance_.lp_;` are the two
    // disjoint LPs; `ekk` is the crate-rule (e) parameter standing for
    // `solver_object.ekk_instance_`.
    if solver_object.lp_.is_moved_ {
        solver_object
            .lp_
            .move_back_lp_and_unapply_scaling(&mut ekk.lp_);
    }
    return_status
}

/// HiGHS: HighsSolution.cpp::formSimplexLpBasisAndFactor (lines 1850-1912, 63 lines)
/// C++ default arguments: only_from_known_basis = false
///
/// * `considerScaling(options, lp)` is `crate::lp_data::lp_scale::consider_scaling`.
/// * The `check_basis` branch never touches `ekk`: it sets `basis.alien =
///   true`, calls [`accommodate_alien_basis`], resets `alien`, and unapplies
///   any scaling it applied (`lp.unapply_scale()`).
/// * `ekk_instance.setBasis(basis)` is `ekk.set_basis_from(&*solver_object.basis_)`
///   (simplex/ekk/solve.rs:242); `initialiseSimplexLpBasisAndFactor` is
///   `ekk.initialise_simplex_lp_basis_and_factor(only_from_known_basis)`
///   (solve.rs:529).
pub fn form_simplex_lp_basis_and_factor(
    ekk: &mut HEkk,
    solver_object: &mut HighsLpSolverObject,
    only_from_known_basis: bool,
) -> HighsStatus {
    // Ideally, forms a SimplexBasis from the HighsBasis in the
    // HighsLpSolverObject
    //
    // If only_from_known_basis is true and
    // initialiseSimplexLpBasisAndFactor finds that there is no simplex
    // basis, then its error return is passed down
    //
    // If only_from_known_basis is false, then the basis is completed
    // with logicals if it is rank deficient (from singularity or being
    // incomplete)
    //
    let mut return_status = HighsStatus::kOk;
    let mut call_status: HighsStatus;
    // C++ references: `lp` = solver_object.lp_, `basis` = solver_object.basis_,
    // `options` = solver_object.options_, `ekk_instance` = ekk (crate rule (e)),
    // `ekk_status` = ekk.status_ (a reference in C++; read at each use, same struct).
    solver_object.lp_.ensure_colwise();
    let passed_scaled = solver_object.lp_.is_scaled_;
    // Consider scaling the LP
    if !passed_scaled {
        crate::lp_data::lp_scale::consider_scaling(solver_object.options_, solver_object.lp_);
    }
    let check_basis =
        solver_object.basis_.alien || (!solver_object.basis_.valid && solver_object.basis_.useful);
    if check_basis {
        // The basis needs to be checked for rank deficiency, and possibly
        // completed if it is rectangular
        //
        // If it's not valid but useful, but not alien,
        // accommodateAlienBasis will assert, so make the basis alien
        solver_object.basis_.alien = true;
        debug_assert!(!only_from_known_basis);
        accommodate_alien_basis(solver_object);
        solver_object.basis_.alien = false;
        // Unapply any scaling used only for factorization to check and
        // complete the basis
        if !passed_scaled {
            solver_object.lp_.unapply_scale();
        }
        // Check that any scaling the LP arrived with has not been removed
        debug_assert!(solver_object.lp_.is_scaled_ == passed_scaled);
        return HighsStatus::kOk;
    }
    // Move the HighsLpSolverObject's LP to EKK
    ekk.move_lp(solver_object);
    if !ekk.status_.has_basis {
        // The Ekk instance has no simplex basis, so pass the HiGHS basis
        let call_status = ekk.set_basis_from(&*solver_object.basis_);
        return_status = interpret_call_status(
            &solver_object.options_.log_options,
            call_status,
            return_status,
            "setBasis",
        );
        if return_status == HighsStatus::kError {
            return form_simplex_lp_basis_and_factor_return(return_status, ekk, solver_object);
        }
    }
    // Now form the invert
    debug_assert!(ekk.status_.has_basis);
    call_status = ekk.initialise_simplex_lp_basis_and_factor(only_from_known_basis);
    // If the current basis cannot be inverted, return an error
    if call_status != HighsStatus::kOk {
        return form_simplex_lp_basis_and_factor_return(HighsStatus::kError, ekk, solver_object);
    }
    // Once the invert is formed, move back the LP and remove any scaling.
    form_simplex_lp_basis_and_factor_return(HighsStatus::kOk, ekk, solver_object)
}

/// HiGHS: HighsSolution.cpp::accommodateAlienBasis (lines 1914-1974, 61 lines)
///
/// A local `HFactor` (util/hfactor): `factor.setup_general(&lp.a_matrix_,
/// num_basic_variables, kDefaultPivotThreshold, kDefaultPivotTolerance,
/// kHighsDebugLevelMin, Some(&options.log_options))`, then
/// `factor.build(&mut basic_index, None)` (the Rust `HFactor` takes
/// `basic_index` per call, crate rule (c)). The rank deficiency indexes
/// `factor.row_with_no_pivot[rank_deficiency + k]`, exactly as C++.
pub fn accommodate_alien_basis(solver_object: &mut HighsLpSolverObject) {
    // C++ references: `lp` = solver_object.lp_, `basis` = solver_object.basis_,
    // `options` = solver_object.options_ (disjoint fields).
    let lp = &*solver_object.lp_;
    let basis = &mut *solver_object.basis_;
    let options = &*solver_object.options_;
    debug_assert!(basis.alien);
    let num_row = lp.num_row_;
    let num_col = lp.num_col_;
    debug_assert!(basis.col_status.len() as HighsInt >= num_col);
    debug_assert!(basis.row_status.len() as HighsInt >= num_row);
    let mut basic_index: Vec<HighsInt> = Vec::new();
    for i_col in 0..num_col {
        if basis.col_status[i_col as usize] == HighsBasisStatus::kBasic {
            basic_index.push(i_col);
        }
    }
    for i_row in 0..num_row {
        if basis.row_status[i_row as usize] == HighsBasisStatus::kBasic {
            basic_index.push(num_col + i_row);
        }
    }
    let mut num_basic_variables: HighsInt = basic_index.len() as HighsInt;
    let mut factor = HFactor::default();
    // HiGHS 1.16 #3354 (h16_hfactor_pivot): the C++ local `HFactor factor;` runs the 1.16.0
    // kernel, so the flag that HSimplexNla::setup sets for the simplex factor is set
    // here too (roll3000).
    factor.h16_hfactor_pivot = options.h16_hfactor_pivot;
    factor.setup_general(
        &lp.a_matrix_,
        num_basic_variables,
        kDefaultPivotThreshold,
        kDefaultPivotTolerance,
        kHighsDebugLevelMin,
        Some(&options.log_options),
    );
    // C++ `HighsInt rank_deficiency = factor.build();` — the Rust HFactor
    // takes `basic_index` per call (crate rule (c)).
    let rank_deficiency = factor.build(&mut basic_index, None);
    // Must not have timed out
    debug_assert!(rank_deficiency >= 0);
    // Deduce the basis from basic_index
    //
    // Set all basic variables to nonbasic
    for i_col in 0..num_col {
        if basis.col_status[i_col as usize] == HighsBasisStatus::kBasic {
            basis.col_status[i_col as usize] = HighsBasisStatus::kNonbasic;
        }
    }
    for i_row in 0..num_row {
        if basis.row_status[i_row as usize] == HighsBasisStatus::kBasic {
            basis.row_status[i_row as usize] = HighsBasisStatus::kNonbasic;
        }
    }
    // Set at most the first num_row variables in basic_index to basic
    let use_basic_variables = cpp_min(num_row, num_basic_variables);
    // num_basic_variables is no longer needed, so can be used as a check
    num_basic_variables = 0;
    for i_row in 0..use_basic_variables {
        let i_var = basic_index[i_row as usize];
        if i_var < num_col {
            basis.col_status[i_var as usize] = HighsBasisStatus::kBasic;
        } else {
            basis.row_status[(i_var - num_col) as usize] = HighsBasisStatus::kBasic;
        }
        num_basic_variables += 1;
    }
    // Complete the assignment of basic variables using the logicals of
    // non-pivotal rows
    let num_missing = num_row - num_basic_variables;
    for k in 0..num_missing {
        let i_row = factor.row_with_no_pivot[(rank_deficiency + k) as usize];
        basis.row_status[i_row as usize] = HighsBasisStatus::kBasic;
        num_basic_variables += 1;
    }
    debug_assert!(num_basic_variables == num_row);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `Highs::setBasis` with an alien basis on the near-dependent LP of
    /// check/TestFactor.cpp (Factor-zero-col-pivot-tolerance): all four columns and row 0
    /// basic. The `accommodateAlienBasis` HFactor runs the #3354 kernel only when it gets
    /// h16_hfactor_pivot (roll3000).
    fn alien_basis_run(flag: bool) -> (HighsStatus, Vec<HighsBasisStatus>, Vec<HighsBasisStatus>) {
        let mut lp = HighsLp {
            num_col_: 4,
            num_row_: 5,
            col_cost_: vec![0.0; 4],
            col_lower_: vec![0.0, 0.0, 0.0, -kHighsInf],
            col_upper_: vec![kHighsInf; 4],
            row_lower_: vec![-kHighsInf; 5],
            row_upper_: vec![0.0, 0.0, 0.0, 0.0, -1.0],
            ..Default::default()
        };
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 5;
        lp.a_matrix_.start_ = vec![0, 4, 8, 11, 14];
        lp.a_matrix_.index_ = vec![0, 2, 3, 4, 0, 2, 3, 4, 0, 1, 4, 1, 2, 4];
        lp.a_matrix_.value_ = vec![
            2.0, -1.0, 3e6, -1.0, -1e13, 1.0, -3000001.0, 1.0, 1.0, 1.0, -2.0, 1e7, -1e7, 1.0,
        ];
        let mut highs = Highs::default();
        assert_eq!(
            highs.set_option_value_bool("output_flag", false),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.set_option_value_bool("h16_hfactor_pivot", flag),
            HighsStatus::kOk
        );
        assert_eq!(highs.pass_model(lp), HighsStatus::kOk);
        let basis = HighsBasis {
            col_status: vec![HighsBasisStatus::kBasic; 4],
            row_status: vec![
                HighsBasisStatus::kBasic,
                HighsBasisStatus::kUpper,
                HighsBasisStatus::kUpper,
                HighsBasisStatus::kUpper,
                HighsBasisStatus::kUpper,
            ],
            valid: true,
            alien: true,
            ..Default::default()
        };
        let status = highs.set_basis(&basis, "");
        let out = highs.get_basis();
        (status, out.col_status.clone(), out.row_status.clone())
    }

    /// Golden values from the highspy 1.16.0 wheel (`Highs.setBasis` on the same LP and
    /// basis): the basis is kept as given.
    #[test]
    fn h16_hfactor_pivot_alien_basis() {
        use HighsBasisStatus::*;
        assert_eq!(
            alien_basis_run(true),
            (
                HighsStatus::kOk,
                vec![kBasic, kBasic, kBasic, kBasic],
                vec![kBasic, kUpper, kUpper, kUpper, kUpper]
            )
        );
    }

    /// Golden values from the highspy 1.15.1 wheel (a release build): the 1.15.1 kernel
    /// finds column 1 rank deficient and replaces it with the logical of row 1. A debug
    /// build fails the 1.15.1 `assert(abs_value < pivot_tolerance)` that #3354 fixes, as
    /// in util/hfactor/kernel.rs.
    #[test]
    #[cfg(not(debug_assertions))]
    fn h16_hfactor_pivot_alien_basis_off() {
        use HighsBasisStatus::*;
        assert_eq!(
            alien_basis_run(false),
            (
                HighsStatus::kOk,
                vec![kBasic, kNonbasic, kBasic, kBasic],
                vec![kBasic, kBasic, kUpper, kUpper, kUpper]
            )
        );
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "abs_value < self.pivot_tolerance")]
    fn h16_hfactor_pivot_alien_basis_off_trips_cpp_assert() {
        alien_basis_run(false);
    }
}
