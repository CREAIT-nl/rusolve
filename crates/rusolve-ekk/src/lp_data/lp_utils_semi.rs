//! HiGHS source: lp_data/HighsLpUtils.cpp (assessSemiVariables,
//! relaxSemiVariables), lp_data/HighsLp.cpp (HighsLp::unapplyMods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The semi-variable modifications `calledOptimizeModel` makes to
//! `model_.lp_` before solving, recorded in `HighsLp::mods_`, and their undo.
//!
//! * `assess_semi_variables` (unless `solve_relaxation`): semi-variables with
//!   zero lower bound become continuous/integer, inconsistent ones are fixed
//!   at zero, and upper bounds above `kMaxSemiVariableUpper` are tightened.
//! * `relax_semi_variables` (with `solve_relaxation`): lower bounds to zero.
//! * `HighsLp::unapply_mods` (from `returnFromOptimizeModel` when
//!   `undo_mods`): restores all of these, then `mods_.clear()`.
//!
//! `highs_run.rs` calls `assess_semi_variables` and `relax_semi_variables` on
//! every path where C++ calls them (no semi-variable pre-check).
//!

use crate::prelude::*;

/// HiGHS: HighsLpUtils.cpp::assessSemiVariables (lines 472-662).
/// `const double kLowerBoundMu = 10.0` is a local, as in C++.
pub fn assess_semi_variables(
    lp: &mut HighsLp,
    options: &HighsOptions,
    made_semi_variable_mods: &mut bool,
) -> HighsStatus {
    *made_semi_variable_mods = false;
    let mut return_status = HighsStatus::kOk;
    if lp.integrality_.is_empty() {
        return return_status;
    }
    debug_assert!(lp.integrality_.len() as HighsInt == lp.num_col_);
    let mut num_illegal_lower: HighsInt = 0;
    let mut num_illegal_upper: HighsInt = 0;
    let mut num_tightened_upper: HighsInt = 0;
    let mut num_inconsistent_semi: HighsInt = 0;
    let mut num_non_semi: HighsInt = 0;
    let mut num_non_continuous_variables: HighsInt = 0;
    const k_lower_bound_mu: f64 = 10.0;

    // lp.mods_.save_non_semi_variable_index
    // lp.mods_.save_inconsistent_semi_variable_index
    // lp.mods_.save_inconsistent_semi_variable_lower_bound_value
    // lp.mods_.save_inconsistent_semi_variable_upper_bound_value
    // lp.mods_.save_inconsistent_semi_variable_type
    // lp.mods_.save_tightened_semi_variable_upper_bound_index
    // lp.mods_.save_tightened_semi_variable_upper_bound_value
    debug_assert!(lp.mods_.save_inconsistent_semi_variable_index.is_empty());
    for i_col in 0..lp.num_col_ {
        if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
            || lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger
        {
            if lp.col_lower_[i_col as usize] > lp.col_upper_[i_col as usize] {
                // Semi-variables with inconsistent bounds become continuous
                // and fixed at zero
                num_inconsistent_semi += 1;
                lp.mods_.save_inconsistent_semi_variable_index.push(i_col);
                lp.mods_
                    .save_inconsistent_semi_variable_lower_bound_value
                    .push(lp.col_lower_[i_col as usize]);
                lp.mods_
                    .save_inconsistent_semi_variable_upper_bound_value
                    .push(lp.col_upper_[i_col as usize]);
                lp.mods_
                    .save_inconsistent_semi_variable_type
                    .push(lp.integrality_[i_col as usize]);
                continue;
            }
            // Semi-variables with zero lower bound are not semi
            if lp.col_lower_[i_col as usize] == 0.0 {
                num_non_semi += 1;
                lp.mods_.save_non_semi_variable_index.push(i_col);
                // Semi-integer become integer so still have a non-continuous variable
                if lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger {
                    num_non_continuous_variables += 1;
                }
                continue;
            }
            if lp.col_lower_[i_col as usize] < 0.0 {
                // Semi-variables must have a positive lower bound
                num_illegal_lower += 1;
            } else if lp.col_upper_[i_col as usize] > kMaxSemiVariableUpper {
                // Semi-variables must have upper bound that's not too large,
                // so see whether the limiting value is sufficiently larger than the
                // lower bound
                if k_lower_bound_mu * lp.col_lower_[i_col as usize] > kMaxSemiVariableUpper {
                    num_illegal_upper += 1;
                } else {
                    // Record the upper bound change to be made later
                    lp.mods_
                        .save_tightened_semi_variable_upper_bound_index
                        .push(i_col);
                    lp.mods_
                        .save_tightened_semi_variable_upper_bound_value
                        .push(kMaxSemiVariableUpper);
                    num_tightened_upper += 1;
                }
            }
            num_non_continuous_variables += 1;
        } else if lp.integrality_[i_col as usize] == HighsVarType::kInteger {
            num_non_continuous_variables += 1;
        }
    }
    if num_inconsistent_semi != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "{} semi-continuous/integer variable(s) have inconsistent bounds so are fixed at zero\n",
            num_inconsistent_semi
        );
        return_status = HighsStatus::kWarning;
    }
    if num_non_semi != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "{} semi-continuous/integer variable(s) have zero lower bound so are continuous/integer\n",
            num_non_semi
        );
        return_status = HighsStatus::kWarning;
    }
    if num_non_continuous_variables == 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "No semi-integer/integer variables in model with non-empty integrality\n"
        );
        return_status = HighsStatus::kWarning;
    }
    let has_illegal_bounds = num_illegal_lower != 0 || num_illegal_upper != 0;
    if num_tightened_upper != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kWarning,
            "{} semi-continuous/integer variable(s) have upper bounds exceeding {} that can be tightened to {} > {}*lower)\n",
            num_tightened_upper,
            crate::util::c_printf::c_fmt_g(kMaxSemiVariableUpper, 6),
            crate::util::c_printf::c_fmt_g(kMaxSemiVariableUpper, 6),
            crate::util::c_printf::c_fmt_g(k_lower_bound_mu, 6)
        );
        return_status = HighsStatus::kWarning;
        if has_illegal_bounds {
            // Don't apply upper bound tightenings if there are illegal bounds
            debug_assert!(num_illegal_lower != 0 || num_illegal_upper != 0);
            lp.mods_
                .save_tightened_semi_variable_upper_bound_index
                .clear();
            lp.mods_
                .save_tightened_semi_variable_upper_bound_value
                .clear();
            num_tightened_upper = 0;
        } else {
            // Apply the upper bound tightenings, saving the over-written
            // values
            for k in 0..num_tightened_upper {
                let use_upper_bound =
                    lp.mods_.save_tightened_semi_variable_upper_bound_value[k as usize];
                let i_col = lp.mods_.save_tightened_semi_variable_upper_bound_index[k as usize];
                lp.mods_.save_tightened_semi_variable_upper_bound_value[k as usize] =
                    lp.col_upper_[i_col as usize];
                lp.col_upper_[i_col as usize] = use_upper_bound;
            }
        }
    }
    if num_inconsistent_semi != 0 {
        if has_illegal_bounds {
            // Don't apply bound changes if there are illegal bounds
            lp.mods_.save_inconsistent_semi_variable_index.clear();
            lp.mods_
                .save_inconsistent_semi_variable_lower_bound_value
                .clear();
            lp.mods_
                .save_inconsistent_semi_variable_upper_bound_value
                .clear();
            lp.mods_.save_inconsistent_semi_variable_type.clear();
            num_inconsistent_semi = 0;
        } else {
            for k in 0..num_inconsistent_semi {
                let i_col = lp.mods_.save_inconsistent_semi_variable_index[k as usize];
                lp.col_lower_[i_col as usize] = 0.0;
                lp.col_upper_[i_col as usize] = 0.0;
                lp.integrality_[i_col as usize] = HighsVarType::kContinuous;
            }
        }
    }
    if num_non_semi != 0 {
        if has_illegal_bounds {
            // Don't apply type changes if there are illegal bounds
            lp.mods_.save_non_semi_variable_index.clear();
        } else {
            for k in 0..num_non_semi {
                let i_col = lp.mods_.save_non_semi_variable_index[k as usize];
                if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous {
                    // Semi-continuous become continuous
                    lp.integrality_[i_col as usize] = HighsVarType::kContinuous;
                } else {
                    // Semi-integer become integer
                    lp.integrality_[i_col as usize] = HighsVarType::kInteger;
                }
            }
        }
    }
    if num_illegal_lower != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "{} semi-continuous/integer variable(s) have negative lower bounds\n",
            num_illegal_lower
        );
        return_status = HighsStatus::kError;
    }
    if num_illegal_upper != 0 {
        crate::highs_log_user!(
            &options.log_options,
            HighsLogType::kError,
            "{} semi-continuous/integer variables have upper bounds exceeding {} that cannot be modified due to large lower bounds\n",
            num_illegal_upper,
            crate::util::c_printf::c_fmt_g(kMaxSemiVariableUpper, 6)
        );
        return_status = HighsStatus::kError;
    }
    *made_semi_variable_mods =
        num_non_semi > 0 || num_inconsistent_semi > 0 || num_tightened_upper > 0;
    debug_assert!(num_non_semi as usize <= lp.mods_.save_non_semi_variable_index.len());
    debug_assert!(
        num_inconsistent_semi as usize <= lp.mods_.save_inconsistent_semi_variable_index.len()
    );
    debug_assert!(
        num_tightened_upper as usize
            <= lp
                .mods_
                .save_tightened_semi_variable_upper_bound_index
                .len()
    );
    // C++ lines 658-660: commented-out alternative for made_semi_variable_mods
    return_status
}

/// HiGHS: HighsLpUtils.cpp::relaxSemiVariables (lines 664-685).
pub fn relax_semi_variables(lp: &mut HighsLp, made_semi_variable_mods: &mut bool) {
    *made_semi_variable_mods = false;
    if lp.integrality_.is_empty() {
        return;
    }
    debug_assert!(lp.integrality_.len() as HighsInt == lp.num_col_);
    // lp.mods_.save_relaxed_semi_variable_lower_bound_index
    // lp.mods_.save_relaxed_semi_variable_lower_bound_value
    debug_assert!(
        lp.mods_
            .save_relaxed_semi_variable_lower_bound_index
            .is_empty()
    );
    for i_col in 0..lp.num_col_ {
        if lp.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
            || lp.integrality_[i_col as usize] == HighsVarType::kSemiInteger
        {
            lp.mods_
                .save_relaxed_semi_variable_lower_bound_index
                .push(i_col);
            lp.mods_
                .save_relaxed_semi_variable_lower_bound_value
                .push(lp.col_lower_[i_col as usize]);
            lp.col_lower_[i_col as usize] = 0.0;
        }
    }
    *made_semi_variable_mods = lp.mods_.save_relaxed_semi_variable_lower_bound_index.len() > 0;
}

impl HighsLp {
    /// HiGHS: HighsLp.cpp::HighsLp::unapplyMods (lines 428-493).
    pub fn unapply_mods(&mut self) {
        // Restore any non-semi types
        let num_non_semi = self.mods_.save_non_semi_variable_index.len() as HighsInt;
        for k in 0..num_non_semi {
            let i_col = self.mods_.save_non_semi_variable_index[k as usize];
            debug_assert!(
                self.integrality_[i_col as usize] == HighsVarType::kContinuous
                    || self.integrality_[i_col as usize] == HighsVarType::kInteger
            );
            if self.integrality_[i_col as usize] == HighsVarType::kContinuous {
                self.integrality_[i_col as usize] = HighsVarType::kSemiContinuous;
            } else {
                self.integrality_[i_col as usize] = HighsVarType::kSemiInteger;
            }
        }
        // Restore any inconsistent semi variables
        let num_inconsistent_semi =
            self.mods_.save_inconsistent_semi_variable_index.len() as HighsInt;
        if num_inconsistent_semi == 0 {
            debug_assert!(
                self.mods_
                    .save_inconsistent_semi_variable_lower_bound_value
                    .is_empty()
            );
            debug_assert!(
                self.mods_
                    .save_inconsistent_semi_variable_upper_bound_value
                    .is_empty()
            );
            debug_assert!(self.mods_.save_inconsistent_semi_variable_type.is_empty());
        }
        for k in 0..num_inconsistent_semi {
            let i_col = self.mods_.save_inconsistent_semi_variable_index[k as usize];
            self.col_lower_[i_col as usize] =
                self.mods_.save_inconsistent_semi_variable_lower_bound_value[k as usize];
            self.col_upper_[i_col as usize] =
                self.mods_.save_inconsistent_semi_variable_upper_bound_value[k as usize];
            self.integrality_[i_col as usize] =
                self.mods_.save_inconsistent_semi_variable_type[k as usize];
        }
        // Restore any relaxed lower bounds
        // relaxed_semi_variable_lower_index = this->mods_.save_relaxed_semi_variable_lower_bound_index
        // relaxed_semi_variable_lower_value = this->mods_.save_relaxed_semi_variable_lower_bound_value
        let num_lower_bound = self
            .mods_
            .save_relaxed_semi_variable_lower_bound_index
            .len() as HighsInt;
        if num_lower_bound == 0 {
            debug_assert!(
                self.mods_
                    .save_relaxed_semi_variable_lower_bound_value
                    .is_empty()
            );
        }
        for k in 0..num_lower_bound {
            let i_col = self.mods_.save_relaxed_semi_variable_lower_bound_index[k as usize];
            debug_assert!(
                self.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
                    || self.integrality_[i_col as usize] == HighsVarType::kSemiInteger
            );
            self.col_lower_[i_col as usize] =
                self.mods_.save_relaxed_semi_variable_lower_bound_value[k as usize];
        }
        // Restore any tightened upper bounds
        // tightened_semi_variable_upper_bound_index = this->mods_.save_tightened_semi_variable_upper_bound_index
        // tightened_semi_variable_upper_bound_value = this->mods_.save_tightened_semi_variable_upper_bound_value
        let num_upper_bound = self
            .mods_
            .save_tightened_semi_variable_upper_bound_index
            .len() as HighsInt;
        if num_upper_bound == 0 {
            debug_assert!(
                self.mods_
                    .save_tightened_semi_variable_upper_bound_value
                    .is_empty()
            );
        }
        for k in 0..num_upper_bound {
            let i_col = self.mods_.save_tightened_semi_variable_upper_bound_index[k as usize];
            debug_assert!(
                self.integrality_[i_col as usize] == HighsVarType::kSemiContinuous
                    || self.integrality_[i_col as usize] == HighsVarType::kSemiInteger
            );
            self.col_upper_[i_col as usize] =
                self.mods_.save_tightened_semi_variable_upper_bound_value[k as usize];
        }

        self.mods_.clear();
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    /// Five columns: sc with lb 0, sc with lb 2 ub 1 (inconsistent),
    /// si with lb 1 ub 1e6 (tightened), integer, continuous.
    fn semi_lp() -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 5;
        lp.col_cost_ = vec![1.0; 5];
        lp.col_lower_ = vec![0.0, 2.0, 1.0, 0.0, 0.0];
        lp.col_upper_ = vec![5.0, 1.0, 1e6, 3.0, 4.0];
        lp.integrality_ = vec![
            HighsVarType::kSemiContinuous,
            HighsVarType::kSemiContinuous,
            HighsVarType::kSemiInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        lp.a_matrix_.start_ = vec![0; 6];
        lp
    }

    #[test]
    fn assess_then_unapply_round_trips() {
        let mut lp = semi_lp();
        let before = lp.clone();
        let options = HighsOptions::default();
        let mut made = false;
        let status = assess_semi_variables(&mut lp, &options, &mut made);
        assert_eq!(status, HighsStatus::kWarning);
        assert!(made);
        assert_eq!(lp.integrality_[0], HighsVarType::kContinuous);
        assert_eq!((lp.col_lower_[1], lp.col_upper_[1]), (0.0, 0.0));
        assert_eq!(lp.integrality_[1], HighsVarType::kContinuous);
        assert_eq!(lp.col_upper_[2], kMaxSemiVariableUpper);
        lp.unapply_mods();
        assert_eq!(lp.col_lower_, before.col_lower_);
        assert_eq!(lp.col_upper_, before.col_upper_);
        assert_eq!(lp.integrality_, before.integrality_);
        assert!(lp.mods_.is_clear());
    }

    #[test]
    fn negative_semi_lower_bound_is_an_error() {
        let mut lp = semi_lp();
        lp.col_lower_[0] = -1.0;
        let mut made = true;
        let status = assess_semi_variables(&mut lp, &HighsOptions::default(), &mut made);
        assert_eq!(status, HighsStatus::kError);
        // has_illegal_bounds: no modification is made
        assert!(!made);
    }

    #[test]
    fn relax_zeroes_semi_lower_bounds() {
        let mut lp = semi_lp();
        let mut made = false;
        relax_semi_variables(&mut lp, &mut made);
        assert!(made);
        assert_eq!(lp.col_lower_, vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(
            lp.mods_.save_relaxed_semi_variable_lower_bound_index,
            vec![0, 1, 2]
        );
        lp.unapply_mods();
        assert_eq!(lp.col_lower_, semi_lp().col_lower_);
    }

    #[test]
    fn large_semi_lower_bound_with_large_upper_is_an_error() {
        // C++ 534-535: kLowerBoundMu * lower > kMaxSemiVariableUpper
        // Column 0 gets a nonzero lower bound: C++ 616-620 clears
        // save_non_semi_variable_index without zeroing num_non_semi, so a
        // non-semi column would leave made_semi_variable_mods true (and trip
        // the C++ assert at 652).
        let mut lp = semi_lp();
        lp.col_lower_[0] = 1.0;
        lp.col_lower_[2] = 2e4;
        let mut made = true;
        let status = assess_semi_variables(&mut lp, &HighsOptions::default(), &mut made);
        assert_eq!(status, HighsStatus::kError);
        // C++ 581-620: illegal bounds discard every recorded modification
        assert!(!made);
        assert!(lp.mods_.save_non_semi_variable_index.is_empty());
        assert!(lp.mods_.save_inconsistent_semi_variable_index.is_empty());
        assert!(
            lp.mods_
                .save_tightened_semi_variable_upper_bound_index
                .is_empty()
        );
        assert_eq!(lp.col_upper_, semi_lp().col_upper_);
        assert_eq!(lp.integrality_, semi_lp().integrality_);
        assert_eq!(lp.col_lower_[1], 2.0);
    }

    #[test]
    fn tightening_is_dropped_when_another_bound_is_illegal() {
        // Column 2 would be tightened, column 0 has a negative lower bound
        let mut lp = semi_lp();
        lp.col_lower_[0] = -1.0;
        let mut made = true;
        let status = assess_semi_variables(&mut lp, &HighsOptions::default(), &mut made);
        assert_eq!(status, HighsStatus::kError);
        assert!(!made);
        assert_eq!(lp.col_upper_[2], 1e6);
        assert!(
            lp.mods_
                .save_tightened_semi_variable_upper_bound_value
                .is_empty()
        );
    }

    #[test]
    fn continuous_only_integrality_warns_without_mods() {
        // C++ 566-571: no semi-integer/integer variables
        let mut lp = semi_lp();
        lp.integrality_ = vec![HighsVarType::kContinuous; 5];
        let mut made = true;
        let status = assess_semi_variables(&mut lp, &HighsOptions::default(), &mut made);
        assert_eq!(status, HighsStatus::kWarning);
        assert!(!made);
        assert!(lp.mods_.is_clear());
    }

    #[test]
    fn semi_integer_zero_lower_becomes_integer_and_restores() {
        let mut lp = semi_lp();
        lp.col_lower_[2] = 0.0;
        let mut made = false;
        let status = assess_semi_variables(&mut lp, &HighsOptions::default(), &mut made);
        assert_eq!(status, HighsStatus::kWarning);
        assert!(made);
        assert_eq!(lp.integrality_[2], HighsVarType::kInteger);
        // Zero lower bound: not tightened (C++ 519-525 continue)
        assert_eq!(lp.col_upper_[2], 1e6);
        lp.unapply_mods();
        assert_eq!(lp.integrality_[2], HighsVarType::kSemiInteger);
        assert!(lp.mods_.is_clear());
    }

    #[test]
    fn lp_without_integrality_is_untouched() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        let mut made = true;
        assert_eq!(
            assess_semi_variables(&mut lp, &HighsOptions::default(), &mut made),
            HighsStatus::kOk
        );
        assert!(!made);
        made = true;
        relax_semi_variables(&mut lp, &mut made);
        assert!(!made);
        lp.mods_.save_inf_cost_variable_index.push(0);
        lp.unapply_mods();
        assert!(lp.mods_.save_inf_cost_variable_index.is_empty());
    }
}
