//! HiGHS source: presolve/HPresolve.cpp (probing: lines 1456-2067)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! MIP-only methods. Read `presolve/hpresolve/mod.rs` first
//! (ownership, "MIP mode", containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

/// HiGHS: HighsIO.cpp::highsTimeToString, NDEBUG branch
/// (`std::to_string(static_cast<int>(time))`). PORT NOTE: the release build is the oracle.
fn highs_time_to_string(time: f64) -> String {
    (time as HighsInt).to_string()
}

/// HiGHS: HighsIO.cpp::highsTimeSecondToString
fn highs_time_second_to_string(time: f64) -> String {
    format!("{}s", highs_time_to_string(time))
}

/// HiGHS: HPresolve.cpp::HPresolve::computeProbingScore (lines 1589-1599).
/// A free function so that `runProbing` can call it inside `with_mip_model!`, where
/// `self` cannot be borrowed.
fn probing_score(
    cliquetable: &HighsCliqueTable,
    num_probes: &[u16],
    col: HighsInt,
) -> (i64, HighsInt) {
    let implics_up = cliquetable.get_num_implications_2(col, true);
    let implics_down = cliquetable.get_num_implications_2(col, false);
    (
        cpp_min(5000_i64, (implics_up as i64) * (implics_down as i64))
            / (1_i64 + num_probes[col as usize] as i64),
        cpp_min(100, implics_up + implics_down),
    )
}

/// HiGHS: HPresolve.cpp::HPresolve::runProbing, `calcNumDel` lambda (lines 1645-1649).
fn probing_calc_num_del(
    probing_num_del_col: HighsInt,
    num_del_start: HighsInt,
    num_implic_substitutions: usize,
    num_clique_substitutions: usize,
) -> HighsInt {
    probing_num_del_col - num_del_start
        + (num_implic_substitutions + num_clique_substitutions) as HighsInt
}

/// HiGHS: HPresolve.cpp::HPresolve::runProbing, `probingLog` lambda (lines 1694-1726).
/// The C++ lambda captures the loop state by reference; here it is passed explicitly.
#[allow(clippy::too_many_arguments)]
fn probing_log_step(
    log_options: &HighsLogOptions,
    silent: bool,
    timeless_log: bool,
    tt: f64,
    tt0: f64,
    i_bin: HighsInt,
    i_bin_probed: HighsInt,
    num_binary: HighsInt,
    log_tt: &mut f64,
    log_i_bin_probed: &mut HighsInt,
) {
    if silent || timeless_log {
        return;
    }
    // Ensure that enough time has elapsed, and at least on binary has been probed
    let log_tt_interval = 5.0;
    if tt > *log_tt + log_tt_interval && i_bin_probed > *log_i_bin_probed {
        // Get the average rate from the start of probing
        debug_assert!(i_bin_probed > 0);
        let rate0 = (tt - tt0) / (i_bin_probed as f64);
        // Get the average rate since last logging
        let dl_i_bin_probed = i_bin_probed - *log_i_bin_probed;
        debug_assert!(dl_i_bin_probed > 0);
        let rate1 = (tt - *log_tt) / (dl_i_bin_probed as f64);
        // Assess the time for probing based on the greater rate
        let rate = cpp_max(rate0, rate1);
        let rate_str = format!(" (rate {}/ms", highs_time_to_string(1e3 * rate));
        let expected_probing_finish_time = tt + rate * ((num_binary - i_bin_probed) as f64);
        let expected_probing_finish_time_str = format!(
            " => expected probing finish time {})",
            highs_time_second_to_string(expected_probing_finish_time)
        );
        let time_str = highs_time_second_to_string(tt);
        crate::highs_log_user!(
            log_options,
            HighsLogType::kInfo,
            "   Considered {} / {} binaries; {} probed {}{} {}\n",
            i_bin,
            num_binary,
            i_bin_probed,
            rate_str,
            expected_probing_finish_time_str,
            time_str
        );
        *log_tt = tt;
        *log_i_bin_probed = i_bin_probed;
    }
}

/// HiGHS: HPresolve.cpp::liftingForProbing, `computeOverallScore` lambda (lines 1932-1938).
fn compute_overall_score(score: f64, numelms: HighsInt, numfillin: HighsInt, best: f64) -> f64 {
    let weight = 0.5;
    weight * (score / best) + (1.0 - weight) * ((numelms - numfillin) as f64) / (numelms as f64)
}

impl HPresolve {
    /// HiGHS latest (#3286, 843e9fdb6): HPresolve.cpp::isBinary (lines 241-244)
    pub fn is_binary(&self, col: HighsInt) -> bool {
        let col_index = col as usize;
        self.model.integrality_[col_index] == HighsVarType::kInteger
            && self.model.col_lower_[col_index] == 0.0
            && self.model.col_upper_[col_index] == 1.0
    }

    /// HiGHS 1.16.0 (#3318): HPresolve.cpp::HPresolve::matrixNonZeroChanged (lines 926-962)
    ///
    /// Only called under `mip_coef_change_impl_bounds`, from
    /// `add_to_matrix` and `unlink`; the 1.15.1 path keeps the unconditional resets inline
    /// there. The `implications.rowModified` call is `x1_row_modified`, as at the 1.15.1 sites.
    pub fn matrix_non_zero_changed(
        &mut self,
        row: HighsInt,
        col: HighsInt,
        old_coef: f64,
        new_coef: f64,
        row_is_deleted: bool,
        col_is_deleted: bool,
    ) {
        let (ri, ci) = (row as usize, col as usize);
        if !col_is_deleted {
            // update dual activity; relaxed activity bound invalidates derived implied
            // bounds
            if !self
                .implied_dual_row_bounds
                .impl_bnds_valid_after_coef_change(
                    &row_dual_bound_arrays!(self),
                    col,
                    row,
                    old_coef,
                    new_coef,
                )
            {
                self.reset_row_dual_implied_bounds_derived_from_col(col);
            }
            // not relaxed; only reset if sourced from this column
            else if self.row_dual_lower_source[ri] == col || self.row_dual_upper_source[ri] == col
            {
                self.reset_row_dual_implied_bounds(row, col);
            }
        } else {
            // col is deleted; unconditionally reset derived implied bounds
            self.reset_row_dual_implied_bounds_derived_from_col(col);
        }

        if !row_is_deleted {
            // update primal activity; relaxed activity bound invalidates derived
            // implied bounds
            if !self.implied_row_bounds.impl_bnds_valid_after_coef_change(
                &col_bound_arrays!(self),
                row,
                col,
                old_coef,
                new_coef,
            ) {
                self.reset_col_implied_bounds_derived_from_row(row);
            }
            // not relaxed; only reset if sourced from this row
            else if self.col_lower_source[ci] == row || self.col_upper_source[ci] == row {
                self.reset_col_implied_bounds(col, row);
            }
        } else {
            // row is deleted; unconditionally reset derived implied bounds
            self.reset_col_implied_bounds_derived_from_row(row);
        }

        // modifications to row invalidate lifting opportunities
        self.clear_lifting_opportunities(row);

        // inform implication table about row modification
        self.x1_row_modified(row);
    }

    /// HiGHS 1.16.0 (#3332): HPresolve.cpp::HPresolve::isStrictlyImpliedForDual (lines 257-267
    /// at v1.16.0)
    pub fn is_strictly_implied_for_dual(&self, col: HighsInt, is_lower_bound: bool) -> bool {
        // checks whether the lower (upper) bound of a column is strictly implied by
        // the rows; then the bound is never active, the reduced cost of the column
        // cannot be positive (negative), and the dual constraint of the column has
        // the side sum_i a_ij y_i >= (<=) c_j. for singleton columns, a bound that
        // is implied within the primal feasibility tolerance suffices.
        let implied_margin = if self.colsize[col as usize] != 1 {
            self.primal_feastol
        } else {
            -self.primal_feastol
        };
        if is_lower_bound {
            self.is_lower_strictly_implied(col, Some(&implied_margin))
        } else {
            self.is_upper_strictly_implied(col, Some(&implied_margin))
        }
    }

    /// HiGHS 1.16.0 (#3332): HPresolve.cpp::HPresolve::resetRowDualImpliedBoundsIfDualSideLost
    /// (lines 899-924)
    ///
    /// Rust-only gate (`mip_dual_side_lost`): a no-op when the
    /// flag is off, so the four call sites in `matrix.rs` (changeColUpper/Lower,
    /// changeImplColUpper/Lower) keep the 1.15.1 path. The 1.16.0
    /// checkUpdateRowDualImpliedBounds refactor onto `is_strictly_implied_for_dual` is
    /// behaviour-preserving and not ported.
    pub fn x1_reset_row_dual_implied_bounds_if_dual_side_lost(
        &mut self,
        col: HighsInt,
        is_lower_bound: bool,
    ) {
        if !self.options.mip_dual_side_lost {
            return;
        }
        // while the column's lower (upper) bound is strictly implied, the dual
        // constraint of the column has the side sum_i a_ij y_i >= (<=) c_j, which
        // yields lower (upper) bounds on y_i for positive coefficients and upper
        // (lower) bounds for negative coefficients. the side >= (<=) is lost when
        // the lower (upper) bound is no longer strictly implied, i.e. when it is
        // tightened such that it may become active, or when the implied lower (upper)
        // bound from the rows is relaxed. then the implied row dual bounds derived
        // from the lost side are reset.
        if self.is_strictly_implied_for_dual(col, is_lower_bound) {
            return;
        }
        if self.impl_row_dual_source_by_col[col as usize].is_empty() {
            return;
        }
        // walk the column rather than implRowDualSourceByCol[col], since the sign of
        // the coefficient is needed and looking it up with findNonzero would splay
        // the row, which invalidates a traversal of the row by a caller
        let mut nonzeros = Vec::new();
        let mut it = self.get_column_vector(col);
        while let Some(nonzero) = it.next_nz(self) {
            nonzeros.push((nonzero.index(), nonzero.value()));
        }
        for (row, value) in nonzeros {
            let r = row as usize;
            if (value > 0.0) == is_lower_bound {
                if self.row_dual_lower_source[r] == col {
                    self.change_impl_row_dual_lower(row, -kHighsInf, -1);
                }
            } else if self.row_dual_upper_source[r] == col {
                self.change_impl_row_dual_upper(row, kHighsInf, -1);
            }
        }
    }

    /// HiGHS 1.16.0 (#3286 + #3318): HPresolveProbingUtil.cpp::HPresolve::normaliseCliqueRows
    /// (lines 16-176)
    ///
    /// Only reached under `mip_clique_normalise` (`prepare_probing`). The #3318 form: skip
    /// redundant rows and rows that already are set packing/partitioning, keep only the
    /// row sides that are not implied by the activity bounds, transform the row into the
    /// complemented space (`transformRow`) and round its sides there, and modify the stored
    /// row only once it is known to be a clique. The row sides change through
    /// `h16_change_row_upper` / `h16_change_row_lower` (`row_bounds_h16.rs`); a redundant lhs
    /// is kept.
    pub fn normalise_clique_rows(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        struct NonZero {
            index: HighsInt,
            value: f64,
            complementation: i8,
            position: HighsInt,
        }

        // `transformRow` lambda (lines 25-42)
        fn transform_row(
            p: &HPresolve,
            nzs: &mut Vec<NonZero>,
            scale: f64,
            update: &mut HighsCDouble,
        ) {
            nzs.clear();
            *update = HighsCDouble::from(0.0);

            for &rowiter in &p.rowpositions {
                let col = p.a_col[rowiter as usize];
                let val = scale * p.a_value[rowiter as usize];
                let cu = col as usize;

                if val < 0.0 {
                    nzs.push(NonZero {
                        index: col,
                        value: -val,
                        complementation: -1,
                        position: rowiter,
                    });
                    *update += HighsCDouble::from(val) * p.model.col_upper_[cu];
                } else {
                    nzs.push(NonZero {
                        index: col,
                        value: val,
                        complementation: 1,
                        position: rowiter,
                    });
                    *update += HighsCDouble::from(val) * p.model.col_lower_[cu];
                }
            }
        }

        // maximum dynamism
        const MAX_DYNAMISM: f64 = 1e5;

        let mut nzs: Vec<NonZero> = Vec::new();
        let mut perm: Vec<HighsInt> = Vec::new();
        let mut row_coefs_int: Vec<f64> = Vec::new();

        for row in 0..self.model.num_row_ {
            let r = row as usize;

            // skip deleted and non-all-integer rows
            if self.row_deleted[r] != 0
                || self.rowsize[r] <= 1
                || self.rowsize_integer[r] != self.rowsize[r]
            {
                continue;
            }

            // skip redundant rows
            if self.is_redundant(row) {
                continue;
            }

            // store row
            self.store_row(row);

            // skip rows that are not all-binary
            let mut all_binary = true;
            let mut is_set_ppc = true;
            row_coefs_int.clear();
            let mut n_comp: HighsInt = 0;
            {
                let stored = self.get_stored_row();
                for nz in stored.iter(self) {
                    all_binary = self.is_binary(nz.index());
                    if !all_binary {
                        break;
                    }
                    is_set_ppc = is_set_ppc && nz.value().abs() == 1.0;
                    row_coefs_int.push(nz.value());
                    if nz.value() < 0.0 {
                        n_comp += 1;
                    }
                }
            }
            if !all_binary {
                continue;
            }

            // skip row if it is already setppc
            if is_set_ppc && self.model.row_upper_[r] == 1.0 - n_comp as f64 {
                continue;
            }

            // skip row if dynamism is too large
            if self.compute_dynamism(self.get_stored_row()) > MAX_DYNAMISM {
                continue;
            }

            // compute scaling factor that makes all coefficients integral
            let mut int_scale = HighsIntegers::integral_scale_vec(
                &row_coefs_int,
                self.options.small_matrix_value,
                self.options.small_matrix_value,
            );
            if int_scale == 0.0 || int_scale > 1e3 {
                continue;
            }

            // check if lhs / rhs are redundant
            let mut lhs = -kHighsInf;
            let mut rhs = kHighsInf;
            if self.implied_row_bounds.get_sum_lower(row, 0.0)
                < self.model.row_lower_[r] - self.primal_feastol
            {
                lhs = self.model.row_lower_[r];
            }
            if self.implied_row_bounds.get_sum_upper(row, 0.0)
                > self.model.row_upper_[r] + self.primal_feastol
            {
                rhs = self.model.row_upper_[r];
            }

            // transform row
            let direction: HighsInt = if rhs < kHighsInf { 1 } else { -1 };
            let mut update = HighsCDouble::from(0.0);
            transform_row(self, &mut nzs, direction as f64 * int_scale, &mut update);

            // update lhs / rhs
            lhs *= direction as f64 * int_scale;
            rhs *= direction as f64 * int_scale;
            if direction < 0 {
                std::mem::swap(&mut lhs, &mut rhs);
            }
            if lhs > -kHighsInf {
                lhs = f64::from((lhs - update - self.primal_feastol).ceil());
            }
            if rhs < kHighsInf {
                rhs = f64::from((rhs - update + self.primal_feastol).floor());
            }

            // equation?
            let equation = rhs == lhs;

            // after rounding, a ranged row may have become an equation; skip if still
            // ranged and not an equation
            if lhs > -kHighsInf && rhs < kHighsInf && !equation {
                continue;
            }

            // sort by descending coefficient value
            let num_bin = nzs.len() as HighsInt;
            perm.clear();
            perm.extend(0..num_bin);
            pdqsort_by(&mut perm, |a, b| {
                nzs[*a as usize].value > nzs[*b as usize].value
            });

            // trivial fixings: variables with coefficient > rhs must be zero
            // in the complemented space
            let mut start: HighsInt = 0;
            for j in 0..num_bin {
                if nzs[perm[j as usize] as usize].value <= rhs + self.primal_feastol {
                    break;
                }
                let nz = &nzs[perm[j as usize] as usize];
                let col = nz.index;
                if nz.complementation == -1 {
                    hpresolve_checked_call!(self.fix_col_to_upper(postsolve_stack, col));
                } else {
                    hpresolve_checked_call!(self.fix_col_to_lower(postsolve_stack, col));
                }
                start += 1;
            }

            // skip row if there are less than two remaining variables or two smallest
            // remaining coefficients do not form a clique
            if num_bin - start < 2
                || nzs[perm[(num_bin - 2) as usize] as usize].value
                    + nzs[perm[(num_bin - 1) as usize] as usize].value
                    <= rhs + self.primal_feastol
            {
                continue;
            }

            // for equations, normalisation to x1 + ... + xn = 1 is only valid
            // for set partitioning rows where all coefficients equal rhs
            if equation
                && (nzs[perm[start as usize] as usize].value
                    != nzs[perm[(num_bin - 1) as usize] as usize].value
                    || nzs[perm[start as usize] as usize].value != rhs)
            {
                continue;
            }

            // re-store row since trivial fixings may have removed nonzeros
            self.store_row(row);

            // scale so that we have a <= inequality or equation
            int_scale *= direction as f64;
            self.scale_stored_row(row, int_scale, true);

            // normalize remaining coefficients to ±1
            let mut num_comp: HighsInt = 0;
            for j in start..num_bin {
                let nz = &nzs[perm[j as usize] as usize];
                let pos = nz.position as usize;
                let delta = nz.complementation as f64 - self.a_value[pos];
                if delta != 0.0 {
                    self.add_to_matrix(row, nz.index, delta);
                }
                if nz.complementation == -1 {
                    num_comp += 1;
                }
            }

            // update row bounds
            self.h16_change_row_upper(row, 1.0 - num_comp as f64, false);
            if equation {
                self.h16_change_row_lower(row, 1.0 - num_comp as f64, false);
            }
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::prepareProbing (lines 1456-1521, 66 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn prepare_probing(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        first_call: &mut bool,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::PROBING_INIT);
        let x1_clique_normalise = self.options.mip_clique_normalise;

        self.shrink_problem(postsolve_stack);

        // HiGHS #3286 (mip_clique_normalise): HPresolve.cpp lines 1754-1758 at 843e9fdb6
        if x1_clique_normalise {
            *first_call = !self.mip().mipdata_.cliques_extracted;
            if *first_call {
                hpresolve_checked_call!(self.normalise_clique_rows(postsolve_stack));
            }
        }

        let mut a_value = std::mem::take(&mut self.model.a_matrix_.value_);
        let mut a_index = std::mem::take(&mut self.model.a_matrix_.index_);
        let mut a_start = std::mem::take(&mut self.model.a_matrix_.start_);
        self.to_csc(&mut a_value, &mut a_index, &mut a_start);
        self.ok_from_csc(&a_value, &a_index, &a_start);
        self.model.a_matrix_.value_ = a_value;
        self.model.a_matrix_.index_ = a_index;
        self.model.a_matrix_.start_ = a_start;

        let nnz = self.num_nonzeros();
        self.mip_mut().mipdata_.cliquetable.set_max_entries(nnz);

        // first tighten all bounds if they have an implied bound that is tighter
        // than their column bound before probing this is not done for continuous
        // columns since it may allow stronger dual presolve and more aggregations
        let huge_bound = self.primal_feastol / kHighsTiny;
        for i in 0..self.model.num_col_ {
            let iu = i as usize;
            if self.impl_col_lower[iu].abs() <= huge_bound
                && self.impl_col_lower[iu] > self.model.col_lower_[iu]
            {
                let v = self.impl_col_lower[iu];
                hpresolve_checked_call!(self.change_col_lower(i, v));
            }

            if self.impl_col_upper[iu].abs() <= huge_bound
                && self.impl_col_upper[iu] < self.model.col_upper_[iu]
            {
                let v = self.impl_col_upper[iu];
                hpresolve_checked_call!(self.change_col_upper(i, v));
            }
        }

        // prepare for domain propagation
        // `setupDomainPropagation` builds `HighsPseudocost(mipsolver)`, which reads
        // `mipdata_->postSolveStack` (HighsPseudocost.cpp:50): the live stack.
        with_mip_model_pss!(self, postsolve_stack, ms => ms.setup_domain_propagation());

        // first call?
        if !x1_clique_normalise {
            *first_call = !self.mip().mipdata_.cliques_extracted;
        }

        self.probing_propagate_domain();
        if self.mip().mipdata_.domains[0].infeasible() {
            return HPresolveResult::kPrimalInfeasible;
        }

        // extract cliques that are part of the formulation every time before probing
        // after the first call we only add cliques that directly correspond to set
        // packing constraints so that the clique merging step can extend/delete them
        if *first_call {
            self.mip_mut().mipdata_.cliques_extracted = true;

            self.probing_extract_cliques(postsolve_stack, true);
            if self.mip().mipdata_.domains[0].infeasible() {
                return HPresolveResult::kPrimalInfeasible;
            }

            // during presolve we keep the objective upper bound without the current
            // offset so we need to update it
            if self.mip().mipdata_.upper_limit != kHighsInf {
                let tmp_limit = self.mip().mipdata_.upper_limit;
                self.mip_mut().mipdata_.upper_limit = tmp_limit - self.model.offset_;
                self.probing_extract_obj_cliques();
                self.mip_mut().mipdata_.upper_limit = tmp_limit;
                if self.mip().mipdata_.domains[0].infeasible() {
                    return HPresolveResult::kPrimalInfeasible;
                }
            }

            self.probing_propagate_domain();
            if self.mip().mipdata_.domains[0].infeasible() {
                return HPresolveResult::kPrimalInfeasible;
            }
        }

        self.probing_cleanup_fixed();
        if self.mip().mipdata_.domains[0].infeasible() {
            return HPresolveResult::kPrimalInfeasible;
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::finaliseProbing (lines 1523-1587, 65 lines)
    /// HiGHS latest (#3264, f2ae2b320): HPresolve.cpp::HPresolve::finaliseProbing (lines 1945-2016)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::finaliseProbing (line 1998)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn finalise_probing(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        first_call: bool,
        num_vars_fixed: &mut HighsInt,
        num_bnds_tightened: &mut HighsInt,
        num_vars_substituted: &mut HighsInt,
        lifted_non_zeros: &mut HighsInt,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::PROBING_FINAL);
        // HiGHS #3220 (mip_fix_last_active): HPresolve.cpp lines 1714-1717 at b7d4e98c5
        let x1b_fix_last = self.options.mip_fix_last_active;

        self.probing_cleanup_fixed();

        if !first_call {
            self.probing_extract_cliques(postsolve_stack, false);
        }
        self.probing_run_clique_merging();

        // apply changes from probing

        // first delete redundant clique inequalities
        let deleted_rows = self
            .mip_mut()
            .mipdata_
            .cliquetable
            .get_deleted_rows()
            .clone();
        for delrow in deleted_rows {
            if self.row_deleted[delrow as usize] == 0 {
                self.remove_row(delrow);
            }
        }
        self.mip_mut()
            .mipdata_
            .cliquetable
            .get_deleted_rows()
            .clear();

        // add nonzeros from clique lifting before removing fixed variables, since
        // this might lead to stronger constraint sides
        let extensionvars =
            std::mem::take(self.mip_mut().mipdata_.cliquetable.get_clique_extensions());
        *lifted_non_zeros += extensionvars.len() as HighsInt;
        for (row, cliqueextension) in extensionvars {
            if self.row_deleted[row as usize] != 0 {
                *lifted_non_zeros -= 1;
                continue;
            }
            let mut val = 1.0;
            if cliqueextension.val() == 0 && self.options.mip_coef_change_impl_bounds {
                // HiGHS #3318 (mip_coef_change_impl_bounds): HPresolveProbingUtil.cpp lines
                // 279-280 at v1.16.0 (#3270's addToRowLower/Upper). Needed
                // once add_to_matrix stops resetting the implied column bounds of the row:
                // lowering the row lower bound may loosen it.
                self.h16_add_to_row_lower(row, HighsCDouble::from(-1.0));
                self.h16_add_to_row_upper(row, HighsCDouble::from(-1.0));
                val = -1.0;
            } else if cliqueextension.val() == 0 {
                self.model.row_lower_[row as usize] -= 1.0;
                self.model.row_upper_[row as usize] -= 1.0;
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 1975-1976 at f2ae2b320
                self.x1_row_modified(row);
                val = -1.0;
            }
            self.add_to_matrix(row, cliqueextension.col() as HighsInt, val);
        }

        // now remove fixed columns and tighten domains
        for i in 0..self.model.num_col_ {
            let iu = i as usize;
            if self.col_deleted[iu] != 0 {
                continue;
            }
            let dom_lower = self.mip().mipdata_.domains[0].col_lower_[iu];
            let dom_upper = self.mip().mipdata_.domains[0].col_upper_[iu];
            let new_lower_bnd = self.model.col_lower_[iu] < dom_lower;
            let new_upper_bnd = self.model.col_upper_[iu] > dom_upper;
            if new_lower_bnd {
                hpresolve_checked_call!(self.change_col_lower(i, dom_lower));
            }
            if new_upper_bnd {
                hpresolve_checked_call!(self.change_col_upper(i, dom_upper));
            }
            if self.mip().mipdata_.domains[0].is_fixed(i) {
                *num_vars_fixed += 1;
                let fixed_value = self.model.col_lower_[iu];
                postsolve_stack.removed_fixed_col(i, fixed_value, 0.0, HighsEmptySlice);
                // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 1998 at 7b890c255
                if self.options.mip_presolve_clique_table {
                    hpresolve_checked_call!(self.remove_fixed_col_x1b(i));
                } else {
                    self.remove_fixed_col(i);
                }
            } else {
                if new_lower_bnd {
                    *num_bnds_tightened += 1;
                }
                if new_upper_bnd {
                    *num_bnds_tightened += 1;
                }
            }
            // HiGHS #3220 (mip_fix_last_active): HPresolve.cpp lines 1714-1717 at b7d4e98c5
            if !x1b_fix_last {
                hpresolve_checked_call!(self.check_limits(postsolve_stack));
            }
        }

        // finally apply substitutions
        hpresolve_checked_call!(
            self.apply_conflict_graph_substitutions(postsolve_stack, num_vars_substituted)
        );

        // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 2012 at f2ae2b320
        let x1_aggregate = self.options.mip_aggregate_implications;
        if x1_aggregate {
            self.aggregate_var_bounds();
        }

        self.check_limits(postsolve_stack)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::computeProbingScore (lines 1589-1599, 11 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn compute_probing_score(&self, col: HighsInt) -> (i64, HighsInt) {
        probing_score(&self.mip().mipdata_.cliquetable, &self.num_probes, col)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::runProbing (lines 1601-1873, 273 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::runProbing (lines 2039-2312)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    ///
    /// PORT NOTE: the C++ `return`s inside the probing loop become `outcome` values
    /// (the loop runs inside one `with_mip_model!`). The `kStopped` return skips
    /// `profiling_->stop`, as in C++.
    pub fn run_probing(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::PROBING);
        self.mip_mut()
            .profiling_
            .start(kMipClockProbingPresolve, false);
        self.probing_early_abort = false;

        let old_num_probed = self.num_probed;

        // prepare probing
        let mut first_call = false;
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 2039-2043 at 7b890c255
        let x1b = self.options.mip_presolve_clique_table;
        if x1b {
            self.mip_mut().mipdata_.cliquetable.set_presolve_flag(true);
        }
        let prepare_result = self.prepare_probing(postsolve_stack, &mut first_call);
        if prepare_result != HPresolveResult::kOk {
            self.mip_mut().profiling_.stop(kMipClockProbingPresolve);
            if x1b {
                self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
            }
            return prepare_result;
        }

        let num_col = self.model.num_col_;
        let num_row = self.model.num_row_;

        // store binary variables in vector with their number of implications on other binaries
        let mut binaries = self.probing_collect_binaries(num_col);

        if !binaries.is_empty() {
            // sort variables with many implications on other binaries first
            pdqsort(&mut binaries[..]);

            let mut num_changed_cols: usize = 0;
            while self.mip().mipdata_.domains[0].get_changed_cols().len() != num_changed_cols {
                let col = self.mip().mipdata_.domains[0].get_changed_cols()[num_changed_cols];
                num_changed_cols += 1;
                if self.mip().mipdata_.domains[0].is_fixed(col) {
                    self.probing_num_del_col += 1;
                }
            }

            let num_cliques_start = self.mip().mipdata_.cliquetable.num_cliques();
            let num_implics_start = self.mip().mipdata_.implications.get_num_implications();
            let num_del_start = self.probing_num_del_col;

            let mut num_del = probing_calc_num_del(
                self.probing_num_del_col,
                num_del_start,
                self.mip().mipdata_.implications.substitutions.len(),
                self.mip_mut()
                    .mipdata_
                    .cliquetable
                    .get_substitutions()
                    .len(),
            );

            let nnz = self.num_nonzeros();
            let submip = self.mip().submip;
            let mut splay_contingent: i64 =
                self.mip().mipdata_.cliquetable.num_neighbourhood_queries
                    + cpp_max(if submip { 0 } else { 100000 }, 10 * nnz) as i64;
            let mut num_fail: HighsInt = 0;

            // collect up to 10 lifting opportunities per row
            let max_num_lift_opps: usize = cpp_max(100000_usize, 10_usize * (num_row as usize));

            // only search for lifting opportunities if at least 2 percent of the
            // variables in the problem are continuous
            let mut num_lift_opps: usize = 0;
            if self.mip().options_mip_.mip_lifting_for_probing != -1
                && self.probing_model_has_percentage_cont_vars(2)
            {
                // store lifting opportunities (drained from the buffer after each probe)
                self.mip_mut()
                    .mipdata_
                    .implications
                    .store_lifting_opportunity = Some(Vec::new());
            }

            // Set up logging for probing
            let silent = self.probing_silent_log();
            let mut i_bin: HighsInt = -1;
            let mut i_bin_probed: HighsInt = -1;
            let num_binary = binaries.len() as HighsInt;
            let tt0 = self.timer.as_ref().unwrap().read();
            let mut tt: f64;
            let mut log_tt = tt0;
            let mut log_i_bin_probed = i_bin_probed;

            // HiGHS #3244 (mip_dual_fix_probing): HPresolve.cpp lines 2217-2218 at 8a5f9c1d8
            let dual_fix_probing_option = self.options.mip_dual_fix_probing;
            let dual_fix_probing_enabled = dual_fix_probing_option && !submip;

            let outcome: Option<(HPresolveResult, bool)> = with_mip_model!(self, ms => {
                let mut outcome: Option<(HPresolveResult, bool)> = None;

                // HiGHS #3244 (mip_dual_fix_probing): HPresolve.cpp lines 2219-2221 at 8a5f9c1d8
                if dual_fix_probing_enabled {
                    crate::detach!(ms.mipdata_.domains[0], |d| d.dfp_recompute_locks(&ms));
                }

                let mut probing_useless: HighsInt = 0;
                for binvar in binaries.iter() {
                    // Count the binaries considered
                    i_bin += 1;
                    let i = binvar.3;
                    if ms.mipdata_.cliquetable.get_substitution(i).is_some()
                        || !ms.mipdata_.domains[0].is_binary(&ms, i)
                    {
                        continue;
                    }

                    // Count the binaries probed
                    i_bin_probed += 1;

                    // Check for timeout
                    tt = self.timer.as_ref().unwrap().read();
                    if tt > self.options.time_limit {
                        crate::highs_log_user!(
                            &self.options.log_options,
                            HighsLogType::kInfo,
                            "Time limit reached in probing: consider not using probing by setting option presolve_rule_off to 2^{} = {}\n",
                            kPresolveRuleProbing,
                            2_i32.pow(kPresolveRuleProbing as u32)
                        );
                        outcome = Some((HPresolveResult::kStopped, false));
                        break;
                    }

                    // Possibly log probing
                    probing_log_step(
                        &self.options.log_options,
                        silent,
                        self.options.timeless_log,
                        tt,
                        tt0,
                        i_bin,
                        i_bin_probed,
                        num_binary,
                        &mut log_tt,
                        &mut log_i_bin_probed,
                    );

                    let tighten_limits = (self.num_probed - old_num_probed) >= 2500;

                    // when a large percentage of columns have been deleted, stop this round
                    // of probing
                    if !tighten_limits {
                        self.probing_early_abort =
                            num_del > cpp_max(1000, (num_row + num_col) / 20);
                    } else {
                        self.probing_early_abort =
                            num_del > cpp_min(1000, (num_row + num_col) / 20);
                    }
                    if self.probing_early_abort {
                        break;
                    }

                    // break in case of too many new implications to not spent ages in probing
                    if ms.mipdata_.cliquetable.is_full()
                        || ms.mipdata_.cliquetable.num_cliques() - num_cliques_start
                            > cpp_max(1000000, 2 * nnz)
                        || ms.mipdata_.implications.get_num_implications() - num_implics_start
                            > cpp_max(1000000, 2 * nnz)
                    {
                        break;
                    }

                    if ms.mipdata_.cliquetable.num_neighbourhood_queries > splay_contingent {
                        break;
                    }

                    if self.probing_contingent - (self.num_probed as i64) < 0 {
                        break;
                    }

                    // rusolve mip_presolve_probing_useless
                    if self.options.mip_presolve_probing_useless > 0
                        && probing_useless >= self.options.mip_presolve_probing_useless
                    {
                        self.probing_contingent = self.num_probed as i64;
                        break;
                    }

                    // rusolve mip_presolve_probing_cap
                    if self.options.mip_presolve_probing_cap && self.num_probed >= 10000 {
                        // exhaust the contingent so the outer loop does not re-enter probing
                        self.probing_contingent = self.num_probed as i64;
                        break;
                    }

                    let mut num_bound_chgs: HighsInt = 0;
                    let mut num_new_cliques: HighsInt = -ms.mipdata_.cliquetable.num_cliques();

                    // HiGHS #3244 (mip_dual_fix_probing): HPresolve.cpp line 2288 at 8a5f9c1d8
                    // (rule F: the C++ calls are unconditional; when not enabled they write
                    // false to a field that is already false, so they sit inside the gate)
                    if dual_fix_probing_enabled {
                        ms.mipdata_.domains[0].set_dual_fix_probing_active(dual_fix_probing_enabled);
                    }

                    let probing_result = crate::detach!(ms.mipdata_.implications, |imp| imp
                        .run_probing(&mut ms, i, &mut num_bound_chgs));

                    // HiGHS #3244 (mip_dual_fix_probing): HPresolve.cpp line 2290 at 8a5f9c1d8
                    if dual_fix_probing_enabled {
                        ms.mipdata_.domains[0].set_dual_fix_probing_active(false);
                    }

                    // drain the lifting opportunities stored by this probe (C++ storeLiftingOpportunity)
                    if let Some(buffer) = ms.mipdata_.implications.store_lifting_opportunity.as_mut()
                    {
                        for (row, col, val, coef) in buffer.drain(..) {
                            // HiGHS #3244 (mip_dual_fix_probing): HPresolve.cpp line 2170 at 8a5f9c1d8
                            // (unconditional in C++, so gated on the option alone, not on
                            // dual_fix_probing_enabled, which also excludes submips)
                            if dual_fix_probing_option && coef == 0.0 {
                                continue;
                            }
                            let (_, inserted) = self
                                .lifting_opportunities
                                .entry(row)
                                .or_default()
                                .insert_or_get(HighsHashTableEntry::new(
                                    (col as HighsUInt, val as HighsUInt),
                                    coef,
                                ));
                            debug_assert!(inserted);
                            num_lift_opps += 1;
                        }
                    }

                    if !probing_result {
                        continue;
                    }
                    // rusolve mip_presolve_probing_useless
                    if self.options.mip_presolve_probing_useless > 0 {
                        let mut progress = num_bound_chgs > 0;
                        let new_cliques = ms.mipdata_.cliquetable.num_cliques() + num_new_cliques;
                        progress |= new_cliques > 0;
                        if progress {
                            probing_useless = 0;
                        } else {
                            probing_useless += 1;
                        }
                    }
                    self.probing_contingent += num_bound_chgs as i64;
                    num_new_cliques += ms.mipdata_.cliquetable.num_cliques();
                    num_new_cliques = cpp_max(num_new_cliques, 0);
                    while ms.mipdata_.domains[0].get_changed_cols().len() != num_changed_cols {
                        let col = ms.mipdata_.domains[0].get_changed_cols()[num_changed_cols];
                        num_changed_cols += 1;
                        if ms.mipdata_.domains[0].is_fixed(col) {
                            self.probing_num_del_col += 1;
                        }
                    }
                    let new_num_del = probing_calc_num_del(
                        self.probing_num_del_col,
                        num_del_start,
                        ms.mipdata_.implications.substitutions.len(),
                        ms.mipdata_.cliquetable.get_substitutions().len(),
                    );

                    if new_num_del > num_del {
                        self.probing_contingent += num_del as i64;
                        if !submip {
                            splay_contingent += (100 * (new_num_del + num_del_start)) as i64;
                            splay_contingent += (1000 * num_new_cliques) as i64;
                        }
                        num_del = new_num_del;
                        num_fail = 0;
                    } else if submip || num_new_cliques == 0 {
                        splay_contingent -=
                            ((if tighten_limits { 250 } else { 100 }) * num_fail) as i64;
                        num_fail += 1;
                    } else {
                        splay_contingent += (1000 * num_new_cliques) as i64;
                        num_fail = 0;
                    }

                    self.num_probed += 1;
                    let np = &mut self.num_probes[i as usize];
                    *np = np.wrapping_add(1);

                    // Stop collecting lifting opportunities if maximum is reached
                    if num_lift_opps >= max_num_lift_opps {
                        ms.mipdata_.implications.store_lifting_opportunity = None;
                    }

                    if ms.mipdata_.domains[0].infeasible() {
                        outcome = Some((HPresolveResult::kPrimalInfeasible, true));
                        break;
                    }
                }
                outcome
            });
            if let Some((result, stop_profiling)) = outcome {
                if stop_profiling {
                    self.mip_mut().profiling_.stop(kMipClockProbingPresolve);
                    // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 2271 at 7b890c255
                    if x1b {
                        self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
                    }
                }
                return result;
            }

            // finalise probing
            let mut num_vars_fixed: HighsInt = 0;
            let mut num_bnds_tightened: HighsInt = 0;
            let mut num_vars_substituted: HighsInt = 0;
            let mut liftednonzeros: HighsInt = 0;
            let finalise_result = self.finalise_probing(
                postsolve_stack,
                first_call,
                &mut num_vars_fixed,
                &mut num_bnds_tightened,
                &mut num_vars_substituted,
                &mut liftednonzeros,
            );
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 2281-2288 at 7b890c255
            if x1b {
                if finalise_result != HPresolveResult::kOk {
                    self.mip_mut().profiling_.stop(kMipClockProbingPresolve);
                    self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
                    return finalise_result;
                }
            } else {
                hpresolve_checked_call!(finalise_result);
            }
            self.probing_num_del_col += num_vars_substituted;

            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "{} probing evaluations: {} deleted rows, {} deleted columns, {} lifted nonzeros\n",
                self.num_probed - old_num_probed,
                self.num_deleted_rows,
                self.num_deleted_cols,
                liftednonzeros
            );

            // lifting for probing
            if self.mip().options_mip_.mip_lifting_for_probing != -1 {
                // only perform lifting if probing did not modify the problem so far
                if self.num_deleted_rows == 0 && self.num_deleted_cols == 0 && liftednonzeros == 0 {
                    hpresolve_checked_call!(self.lifting_for_probing(postsolve_stack));
                }
                // clear lifting opportunities
                self.lifting_opportunities.clear();
                self.mip_mut()
                    .mipdata_
                    .implications
                    .store_lifting_opportunity = None;
            }
        }

        self.mip_mut().profiling_.stop(kMipClockProbingPresolve);
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 2312 at 7b890c255
        if x1b {
            self.mip_mut().mipdata_.cliquetable.set_presolve_flag(false);
        }
        self.check_limits(postsolve_stack)
    }

    /// HiGHS: HPresolve.cpp::HPresolve::liftingForProbing (lines 1875-2067, 193 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    ///
    /// PORT NOTE: `liftingOpportunities` is a `BTreeMap` (row order) here, not the libc++
    /// `unordered_map` order. The order only decides the `redundantRow` postsolve order
    /// of rows removed in the loop below; it is non-default (`mip_lifting_for_probing`
    /// defaults to -1), so the emulation of the libc++ bucket order is left out.
    /// A `redundant_row` error path leaves the map taken out of `self` (the C++ keeps it),
    /// which is harmless because a non-OK result ends presolve.
    pub fn lifting_for_probing(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        // this method implements lifting for probing as described by Achterberg et
        // al. (2019) Presolve Reductions in Mixed Integer Programming. INFORMS
        // Journal on Computing 32(2):473-506.

        // collect best lifting opportunity for each row in a vector
        // (row, best clique as (CliqueVar, coefficient), best score, best number of fill-ins)
        let mut liftingtable: Vec<(
            HighsInt,
            Vec<(HighsCliqueTable_CliqueVar, f64)>,
            f64,
            HighsInt,
        )> = Vec::new();

        // remember overall best score
        let mut bestscoretotal = -kHighsInf;

        // is lifting allowed to add non-zeros?
        let fillallowed = self.mip().options_mip_.mip_lifting_for_probing > 0;

        // store clique variables and coefficients in a map, ordered by (col, val) as the
        // C++ comparator `c1.col < c2.col || (c1.col == c2.col && c1.val < c2.val)`
        let mut coefficients: BTreeMap<(HighsUInt, HighsUInt), (f64, HighsInt)> = BTreeMap::new();

        // consider lifting opportunities
        let mut numrowsremoved: usize = 0;
        let lifting_opportunities = std::mem::take(&mut self.lifting_opportunities);
        for (&row, htree) in lifting_opportunities.iter() {
            // skip deleted rows
            if self.row_deleted[row as usize] != 0 {
                continue;
            }

            // do not add non-zeros to dense rows
            let dense = self.rowsize[row as usize]
                > cpp_max(1000, (self.model.num_col_ - self.num_deleted_cols) / 20);

            // iterate over elements in hash tree
            let mut entries: Vec<((HighsUInt, HighsUInt), f64)> = Vec::new();
            htree.for_each(&mut |data, coef| entries.push((*data, *coef)));
            let mut isredundant = false;
            for ((col_u, val_u), coef) in entries {
                let col = col_u as HighsInt;
                let val = val_u as HighsInt;
                let pos = self.find_nonzero(row, col);
                isredundant = isredundant || htree.contains(&(col_u, (1 - val) as HighsUInt));
                if !dense
                    && (fillallowed || pos != -1)
                    && self.col_deleted[col as usize] == 0
                    && !self.mip().mipdata_.domains[0].is_fixed(col)
                {
                    coefficients.insert((col_u, val_u), (coef, pos));
                }
            }

            // remove redundant rows
            if isredundant {
                numrowsremoved += 1;
                postsolve_stack.redundant_row(row);
                self.remove_row(row);
                coefficients.clear();
                hpresolve_checked_call!(self.check_limits(postsolve_stack));
                continue;
            }

            // skip rows with empty coefficient map
            if coefficients.is_empty() {
                continue;
            }

            // vector to hold best clique
            let mut bestclique: Vec<(HighsCliqueTable_CliqueVar, f64)> = Vec::new();
            let mut bestscore = -kHighsInf;
            let mut bestnfill: HighsInt = 0;

            // store candidates in a vector
            let mut candidates: Vec<HighsCliqueTable_CliqueVar> =
                Vec::with_capacity(coefficients.len());
            for (&(col, val), &(coef, pos)) in coefficients.iter() {
                let cliquevar = HighsCliqueTable_CliqueVar::from_parts(col, val);
                candidates.push(cliquevar);
                // initialize best clique
                let score = coeff_diff(coef, pos, &self.a_value);
                if score > bestscore {
                    bestscore = score;
                    bestnfill = if pos == -1 { 1 } else { 0 };
                    bestclique = vec![(cliquevar, coef)];
                }
            }

            if candidates.len() > 1 {
                // compute cliques
                let feastol = self.primal_feastol;
                let cliques = self
                    .mip_mut()
                    .mipdata_
                    .cliquetable
                    .compute_maximal_cliques(&candidates, feastol);

                // identify clique with highest score
                for clique in cliques.iter() {
                    let mut score = HighsCDouble::from(0.0);
                    let mut nfill: HighsInt = 0;
                    for cliquevar in clique.iter() {
                        let (coef, pos) = *coefficients
                            .entry((cliquevar.col(), cliquevar.val()))
                            .or_insert((0.0, 0));
                        score += coeff_diff(coef, pos, &self.a_value);
                        if pos == -1 {
                            nfill += 1;
                        }
                    }
                    if score > bestscore {
                        bestscore = f64::from(score);
                        bestnfill = nfill;
                        bestclique.clear();
                        for cliquevar in clique.iter() {
                            let (coef, _) = *coefficients
                                .entry((cliquevar.col(), cliquevar.val()))
                                .or_insert((0.0, 0));
                            bestclique.push((*cliquevar, coef));
                        }
                    }
                }
            }

            // store best clique
            liftingtable.push((row, bestclique, bestscore, bestnfill));
            bestscoretotal = cpp_max(bestscoretotal, bestscore);
            coefficients.clear();
        }

        // sort according to score
        pdqsort_by(&mut liftingtable[..], |opp1, opp2| {
            let score1 =
                compute_overall_score(opp1.2, opp1.1.len() as HighsInt, opp1.3, bestscoretotal);
            let score2 =
                compute_overall_score(opp2.2, opp2.1.len() as HighsInt, opp2.3, bestscoretotal);
            if score1 == score2 {
                opp1.0 < opp2.0
            } else {
                score1 > score2
            }
        });

        // perform actual lifting
        let mut nfill: usize = 0;
        let mut nmod: usize = 0;
        let mut numrowsmodified: usize = 0;
        let maxnfill: usize = cpp_max(
            10 * liftingtable.len(),
            (self.num_nonzeros() as usize) / 100,
        );
        for lifting in liftingtable.iter() {
            // get clique
            let row = lifting.0;
            let bestclique = &lifting.1;

            // check against max. fill-in
            let newfill = lifting.3 as usize;
            if nfill + newfill > maxnfill {
                break;
            }
            nfill += newfill;
            nmod += bestclique.len() - newfill;

            // update matrix
            let mut update = HighsCDouble::from(0.0);
            for (cliquevar, coeff) in bestclique.iter() {
                // add non-zero to matrix
                self.add_to_matrix(row, cliquevar.col() as HighsInt, *coeff);
                // compute term to update left-hand / right-hand side
                if cliquevar.val() == 0 {
                    update += *coeff;
                }
            }

            // update left-hand / right-hand sides
            numrowsmodified += 1;
            if self.model.row_lower_[row as usize] != -kHighsInf {
                self.model.row_lower_[row as usize] += f64::from(update);
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 2482 at f2ae2b320 (1.15.1 2054)
                self.x1_row_modified(row);
            }
            if self.model.row_upper_[row as usize] != kHighsInf {
                self.model.row_upper_[row as usize] += f64::from(update);
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 2483 at f2ae2b320 (1.15.1 2056)
                self.x1_row_modified(row);
            }
        }

        crate::highs_log_dev!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "Lifting for probing removed {} and modified {} row(s), added {} new and modified {} existing nonzero(s)\n",
            numrowsremoved,
            numrowsmodified,
            nfill,
            nmod
        );

        HPresolveResult::kOk
    }

    // ----- module-private helpers (C++ lambdas and the MIP `with_mip_model!` bodies) -----

    /// HiGHS: HPresolve.cpp::HPresolve::silentLog (lines 6326-6328), inlined:
    /// `mipsolver && mipsolver->mipdata_->numRestarts > 0`.
    fn probing_silent_log(&self) -> bool {
        self.mipsolver.is_some() && self.mip().mipdata_.num_restarts > 0
    }

    /// HiGHS: HPresolve.cpp::HPresolve::runProbing, `modelHasPercentageContVars` lambda (lines 1667-1676).
    fn probing_model_has_percentage_cont_vars(&self, percentage: usize) -> bool {
        let mut num_cols: usize = 0;
        let mut num_cont_cols: usize = 0;
        for col in 0..self.colsize.len() {
            if self.col_deleted[col] != 0 {
                continue;
            }
            num_cols += 1;
            if self.model.integrality_[col] == HighsVarType::kContinuous {
                num_cont_cols += 1;
            }
        }
        100 * num_cont_cols >= percentage * num_cols
    }

    /// HiGHS: HPresolve.cpp::HPresolve::runProbing, binaries collection (lines 1629-1643).
    fn probing_collect_binaries(
        &mut self,
        num_col: HighsInt,
    ) -> Vec<(i64, HighsInt, HighsInt, HighsInt)> {
        with_mip_model!(self, ms => {
            let mut binaries: Vec<(i64, HighsInt, HighsInt, HighsInt)> = Vec::new();
            if !ms.mipdata_.cliquetable.is_full() {
                let mut random = HighsRandom::new(self.options.random_seed as HighsUInt);
                for i in 0..num_col {
                    if ms.mipdata_.domains[0].is_binary(&ms, i) {
                        let probing_score =
                            probing_score(&ms.mipdata_.cliquetable, &self.num_probes, i);
                        binaries.push((-probing_score.0, -probing_score.1, random.integer(), i));
                    }
                }
            }
            binaries
        })
    }

    /// `domain.propagate()` on the global domain (`mipsolver->mipdata_->getDomain()`).
    fn probing_propagate_domain(&mut self) {
        let _prof = super::presolve_prof::scope(super::presolve_prof::DOMAIN_PROPAGATE);
        with_mip_model!(self, ms => crate::detach!(ms.mipdata_.domains[0], |d| {
            let _ = d.propagate(&mut ms);
        }));
    }

    /// `cliquetable.extractCliques(*mipsolver, transformRows)`.
    ///
    /// PORT NOTE: takes `postsolve_stack` because `extractCliques` reads
    /// `mipdata_->postSolveStack.getOrigRowIndex` (HighsCliqueTable.cpp:1286), which is
    /// the stack presolve runs on (`with_mip_model_pss!`).
    fn probing_extract_cliques(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        transform_rows: bool,
    ) {
        let _prof = super::presolve_prof::scope(super::presolve_prof::CLIQUE_EXTRACT);
        with_mip_model_pss!(self, postsolve_stack, ms => crate::detach!(ms.mipdata_.cliquetable, |ct| {
            ct.extract_cliques_2(&mut ms, transform_rows);
        }));
    }

    /// `cliquetable.extractObjCliques(*mipsolver)`.
    fn probing_extract_obj_cliques(&mut self) {
        with_mip_model!(self, ms => crate::detach!(ms.mipdata_.cliquetable, |ct| {
            ct.extract_obj_cliques(&mut ms);
        }));
    }

    /// `cliquetable.cleanupFixed(domain)`.
    fn probing_cleanup_fixed(&mut self) {
        with_mip_model!(self, ms => crate::detach!(ms.mipdata_.cliquetable, |ct| {
            crate::detach!(ms.mipdata_.domains[0], |d| ct.cleanup_fixed(&mut ms, d));
        }));
    }

    /// `cliquetable.runCliqueMerging(domain)` (the one-argument overload).
    fn probing_run_clique_merging(&mut self) {
        let _prof = super::presolve_prof::scope(super::presolve_prof::CLIQUE_MERGE);
        with_mip_model!(self, ms => crate::detach!(ms.mipdata_.cliquetable, |ct| {
            crate::detach!(ms.mipdata_.domains[0], |d| ct.run_clique_merging_2(&mut ms, d));
        }));
    }
}

/// HiGHS: HPresolve.cpp::liftingForProbing, `computeCoeffDiff` lambda (lines 1905-1907):
/// `std::fabs(newvalue - (nzpos == -1 ? 0 : Avalue[nzpos]))`.
fn coeff_diff(new_value: f64, nzpos: HighsInt, a_value: &[f64]) -> f64 {
    (new_value
        - if nzpos == -1 {
            0.0
        } else {
            a_value[nzpos as usize]
        })
    .abs()
}

#[cfg(test)]
mod x1_clqn_b_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::util::highs_timer::HighsTimer;

    fn make_presolve_with_lp(lp: &HighsLp, flag: bool) -> HPresolve {
        let mut presolve = HPresolve::default();
        let mut options = HighsOptions::default();
        options.mip_clique_normalise = flag;
        options.small_matrix_value = 1e-12;
        assert!(presolve.ok_set_input(lp.clone(), &options, 0, Some(HighsTimer::new())));
        // LP mode makes every column continuous before the CSC scan, so rowsize_integer
        // would be 0; restore the declared integrality and rebuild the triplets.
        presolve.model.integrality_ = lp.integrality_.clone();
        assert!(presolve.ok_from_csc(
            &lp.a_matrix_.value_,
            &lp.a_matrix_.index_,
            &lp.a_matrix_.start_
        ));
        presolve
    }

    fn make_lp_4col_4row(
        row_lower: Vec<f64>,
        row_upper: Vec<f64>,
        a_start: Vec<HighsInt>,
        a_index: Vec<HighsInt>,
        a_value: Vec<f64>,
    ) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 4;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![1.0, 1.0, 1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 1.0, 1.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = row_lower;
        lp.row_upper_ = row_upper;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = a_start;
        lp.a_matrix_.index_ = a_index;
        lp.a_matrix_.value_ = a_value;
        lp
    }

    fn make_lp_3col_1row(
        row_lower: Vec<f64>,
        row_upper: Vec<f64>,
        a_start: Vec<HighsInt>,
        a_index: Vec<HighsInt>,
        a_value: Vec<f64>,
    ) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 1;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![1.0, 1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 1.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = row_lower;
        lp.row_upper_ = row_upper;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = a_start;
        lp.a_matrix_.index_ = a_index;
        lp.a_matrix_.value_ = a_value;
        lp
    }

    fn make_lp_2col_1row(
        row_lower: Vec<f64>,
        row_upper: Vec<f64>,
        a_start: Vec<HighsInt>,
        a_index: Vec<HighsInt>,
        a_value: Vec<f64>,
    ) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.sense_ = ObjSense::kMinimize;
        lp.col_cost_ = vec![1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0];
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        lp.row_lower_ = row_lower;
        lp.row_upper_ = row_upper;
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = a_start;
        lp.a_matrix_.index_ = a_index;
        lp.a_matrix_.value_ = a_value;
        lp
    }

    #[test]
    fn x1_clqn_b_clique_extract_origin() {
        // row 0: x0 + x1 + x2 <= 1 (set packing, 3-clique)
        // row 1: -x0 + x3 >= 0
        // row 2: -x1 + x3 >= 0
        // row 3: -x2 + x3 >= 0
        // C++ checks mipsolver.numRow() <= 1 after full presolve; here we check
        // the normalisation result on row 0 directly.
        let lp = make_lp_4col_4row(
            vec![-kHighsInf, 0.0, 0.0, 0.0],
            vec![1.0, kHighsInf, kHighsInf, kHighsInf],
            vec![0, 2, 4, 6, 9],
            vec![0, 1, 0, 2, 0, 3, 1, 2, 3],
            vec![1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(4, 4);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        // row 0 should remain as a clique: x0 + x1 + x2 <= 1
        assert_eq!(presolve.model.row_upper_[0], 1.0);
        assert_eq!(presolve.model.row_lower_[0], -kHighsInf);
    }

    #[test]
    fn x1_clqn_b_normalise_unequal_coeff() {
        // row 0: x0 + 2*x1 + 2*x2 <= 2
        // row 1: -x0 + x3 >= 0
        // row 2: -x1 + x3 >= 0
        // row 3: -x2 + x3 >= 0
        // C++ checks mipsolver.numRow() <= 1 after full presolve; here we check
        // the normalisation result on row 0 directly.
        let lp = make_lp_4col_4row(
            vec![-kHighsInf, 0.0, 0.0, 0.0],
            vec![2.0, kHighsInf, kHighsInf, kHighsInf],
            vec![0, 2, 4, 6, 9],
            vec![0, 1, 0, 2, 0, 3, 1, 2, 3],
            vec![1.0, -1.0, 2.0, -1.0, 2.0, -1.0, 1.0, 1.0, 1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(4, 4);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        // row 0 normalised to x0 + x1 + x2 <= 1
        assert_eq!(presolve.model.row_upper_[0], 1.0);
        assert_eq!(presolve.model.row_lower_[0], -kHighsInf);
    }

    #[test]
    fn x1_clqn_b_normalise_equation() {
        // 2*x0 + 2*x1 + 2*x2 = 2 (set partitioning)
        // C++ checks cliquetable.haveCommonClique; here we check the normalised
        // row bounds which imply the clique structure.
        let lp = make_lp_3col_1row(
            vec![2.0],
            vec![2.0],
            vec![0, 1, 2, 3],
            vec![0, 0, 0],
            vec![2.0, 2.0, 2.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, 3);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        // normalised to x0 + x1 + x2 = 1
        assert_eq!(presolve.model.row_upper_[0], 1.0);
        assert_eq!(presolve.model.row_lower_[0], 1.0);
    }

    #[test]
    fn x1_clqn_b_trivial_fixing() {
        // 3*x0 + x1 + x2 <= 2: x0 has coefficient 3 > 2, fixed to lower bound
        let lp = make_lp_3col_1row(
            vec![-kHighsInf],
            vec![2.0],
            vec![0, 1, 2, 3],
            vec![0, 0, 0],
            vec![3.0, 1.0, 1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, 3);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        // x0 fixed to lower bound -> 1 reduction on postsolve stack
        assert_eq!(postsolve.num_reductions(), 1);
    }

    #[test]
    fn x1_clqn_b_complemented() {
        // 3*x0 - x1 - x2 <= -1: after complementing x1,x2, x0 coeff 3 > 1, fixed to 0.
        // Remaining complemented vars normalised to -1, row_upper = 1 - 2 = -1.
        let lp = make_lp_3col_1row(
            vec![-kHighsInf],
            vec![-1.0],
            vec![0, 1, 2, 3],
            vec![0, 0, 0],
            vec![3.0, -1.0, -1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, 3);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        assert_eq!(postsolve.num_reductions(), 1);
        assert_eq!(presolve.model.row_upper_[0], -1.0);
        assert_eq!(presolve.model.row_lower_[0], -kHighsInf);
    }

    #[test]
    fn x1_clqn_b_fix_to_upper() {
        // -3*x0 - x1 - x2 <= -4: all complemented, x0 fixed to upper bound.
        // Remaining normalised to -1, row_upper = 1 - 2 = -1.
        let lp = make_lp_3col_1row(
            vec![-kHighsInf],
            vec![-4.0],
            vec![0, 1, 2, 3],
            vec![0, 0, 0],
            vec![-3.0, -1.0, -1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, 3);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        assert_eq!(postsolve.num_reductions(), 1);
        assert_eq!(presolve.model.row_upper_[0], -1.0);
        assert_eq!(presolve.model.row_lower_[0], -kHighsInf);
    }

    #[test]
    fn x1_clqn_b_non_integral_rhs() {
        // 2*x0 + 2*x1 <= 3.7: integralScale=0.5 gives x0+x1<=1.85, floor=1.
        let lp = make_lp_2col_1row(
            vec![-kHighsInf],
            vec![3.7],
            vec![0, 1, 2],
            vec![0, 0],
            vec![2.0, 2.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, 2);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        assert_eq!(presolve.model.row_upper_[0], 1.0);
        assert_eq!(presolve.model.row_lower_[0], -kHighsInf);
    }

    #[test]
    fn x1_clqn_b_no_delete_ranged_row() {
        // Ranged rows must not be deleted by clique merging.
        // row 0: x0 + x1 + x2 <= 1
        // row 1: 0 <= -x0 + x3 <= 1 (ranged)
        // row 2: 0 <= -x1 + x3 <= 1 (ranged)
        // row 3: 0 <= -x2 + x3 <= 1 (ranged)
        // C++ checks cliquetable.getDeletedRows().empty(); here we check that
        // ranged rows are not normalised (their bounds stay ranged).
        let lp = make_lp_4col_4row(
            vec![-kHighsInf, 0.0, 0.0, 0.0],
            vec![1.0, 1.0, 1.0, 1.0],
            vec![0, 2, 4, 6, 9],
            vec![0, 1, 0, 2, 0, 3, 1, 2, 3],
            vec![1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(4, 4);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        // Row 0 is a valid clique and gets normalised.
        assert_eq!(presolve.model.row_upper_[0], 1.0);
        // 1.16.0 (#3318): the upper side of rows 1-3 is implied (activity <= 1), so the
        // row is the clique x0' + x3 <= 1 in the complemented space: (0, 1) becomes
        // (-1, 0), its lhs kept (the 1.16.0 oracle does the same to row c09_ranged_lhs of
        // tools/highs-oracle/lps/g1/clqn_rows.mps).
        assert_eq!(presolve.model.row_upper_[1], 0.0);
        assert_eq!(presolve.model.row_lower_[1], -1.0);
    }

    #[test]
    fn x1_clqn_b_implied_equality() {
        // 0.5 <= x0 + x1 <= 1: after floor/ceil rounding, becomes equality.
        // C++ checks domain propagation forcing x1=1 when x0=0; here we check
        // the normalised row bounds which make it an equality.
        let lp = make_lp_2col_1row(
            vec![0.5],
            vec![1.0],
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 1.0],
        );
        let mut presolve = make_presolve_with_lp(&lp, true);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, 2);
        let result = presolve.normalise_clique_rows(&mut postsolve);
        assert_eq!(result, HPresolveResult::kOk);
        // 1.16.0 (#3318): x0 + x1 <= 1 already is set packing, so the row is skipped
        // and keeps its bounds (as the 1.16.0 oracle does for row c07_ranged_setppc of
        // tools/highs-oracle/lps/g1/clqn_rows.mps); the C++ test's propagation sees the
        // ranged row directly.
        assert_eq!(presolve.model.row_upper_[0], 1.0);
        assert_eq!(presolve.model.row_lower_[0], 0.5);
    }

    /// HiGHS 1.16.0 (#3332) check/TestMipSolver.cpp "implied-row-dual-bound-invalidation":
    /// the fuzzer model whose optimum (-4) a stale implied row dual bound cut off. Solved
    /// with presolve on and off under every combination of the `mip_*` clique/implied-bound flags: without
    /// `mip_dual_side_lost`, presolve declares it infeasible, as HiGHS 1.15.1 does.
    #[test]
    fn x1_clqn_d_implied_row_dual_bound_invalidation() {
        let mut lp = HighsLp::default();
        lp.num_col_ = 6;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![0.0, 0.0, 0.0, 0.0, 4.0, -2.0];
        lp.col_lower_ = vec![-kHighsInf, -5.0, 0.0, -5.0, 0.0, 0.0];
        lp.col_upper_ = vec![-10.0, kHighsInf, kHighsInf, 5.0, 7.0, 2.0];
        lp.row_lower_ = vec![-13.0, 11.0];
        lp.row_upper_ = vec![kHighsInf, 11.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.start_ = vec![0, 2, 4, 5, 6, 8, 10];
        lp.a_matrix_.index_ = vec![0, 1, 0, 1, 1, 1, 0, 1, 0, 1];
        lp.a_matrix_.value_ = vec![-3.0, 1.0, -1.0, -1.0, -2.0, -1.0, -8.0, 1.0, -1.0, 10.0];

        for flags in 0..8 {
            for presolve in ["on", "off"] {
                let mut h = Highs::default();
                h.set_option_value_bool("output_flag", false);
                h.set_option_value_str("presolve", presolve);
                h.set_option_value_bool("mip_clique_normalise", flags & 1 != 0);
                h.set_option_value_bool("mip_coef_change_impl_bounds", flags & 2 != 0);
                h.set_option_value_bool("mip_dual_side_lost", flags & 4 != 0);
                h.pass_model(lp.clone());
                h.run();
                let status = h.get_model_status();
                let obj = h.get_info().objective_function_value;
                if flags & 4 == 0 && presolve == "on" {
                    // 1.15.1 parity: HiGHS 1.15.1 (pip highspy) also reports infeasible
                    assert_eq!(status, HighsModelStatus::kInfeasible, "flags {flags}");
                } else {
                    assert_eq!(
                        status,
                        HighsModelStatus::kOptimal,
                        "flags {flags} {presolve}"
                    );
                    assert!((obj + 4.0).abs() < 1e-6, "flags {flags} {presolve}: {obj}");
                }
            }
        }
    }

    #[test]
    fn x1_clqn_b_off_unchanged() {
        // With flag off, scale_stored_row uses the old round() path.
        let mut presolve = HPresolve::default();
        presolve.options.mip_clique_normalise = false;
        presolve.model.num_row_ = 1;
        presolve.model.num_col_ = 2;
        presolve.model.row_upper_ = vec![3.7];
        presolve.model.row_lower_ = vec![-kHighsInf];
        presolve.row_dual_lower = vec![0.0];
        presolve.row_dual_upper = vec![0.0];
        presolve.impl_row_dual_lower = vec![0.0];
        presolve.impl_row_dual_upper = vec![0.0];
        presolve.row_dual_lower_source = vec![0];
        presolve.row_dual_upper_source = vec![0];
        presolve.rowpositions = vec![];
        presolve.a_value = vec![];
        presolve.implied_row_bounds.set_num_sums(1);
        presolve.primal_feastol = 1e-7;
        presolve.options.small_matrix_value = 1e-12;

        presolve.scale_stored_row(0, 0.5, true);
        // Off path: round(3.7 * 0.5) = round(1.85) = 2.0
        assert_eq!(presolve.model.row_upper_[0], 2.0);
    }
}

#[cfg(test)]
mod x1b_clqp_e_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;

    /// Pure binary MIP (r0: x0 + x1 + x2 <= 1, r1: x0 + x2 <= 1) in MIP-mode presolve, with
    /// the clique-table presolve flag set to `flag_before`; runs runProbing and returns
    /// (result, presolve flag afterwards, cliques extracted).
    fn probe(flag: bool, flag_before: bool) -> (HPresolveResult, bool, bool) {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 2;
        lp.col_cost_ = vec![-1.0, -2.0, -1.0];
        lp.col_lower_ = vec![0.0; 3];
        lp.col_upper_ = vec![1.0; 3];
        lp.integrality_ = vec![HighsVarType::kInteger; 3];
        lp.row_lower_ = vec![-kHighsInf; 2];
        lp.row_upper_ = vec![1.0; 2];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 2;
        lp.a_matrix_.start_ = vec![0, 2, 3, 5];
        lp.a_matrix_.index_ = vec![0, 1, 0, 0, 1];
        lp.a_matrix_.value_ = vec![1.0; 5];
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.mip_presolve_clique_table = flag;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.timer_.start();
        ms.new_data();
        ms.init();
        ms.setup_domain_propagation();
        let mut postsolve_stack = std::mem::take(&mut ms.mipdata_.post_solve_stack);
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), -1));
        let silent = presolve.silent_log();
        presolve
            .analysis_
            .setup(&presolve.model, &presolve.options, silent);
        presolve
            .mip_mut()
            .mipdata_
            .cliquetable
            .set_presolve_flag(flag_before);
        let result = presolve.run_probing(&mut postsolve_stack);
        let ms = presolve.mip();
        (
            result,
            ms.mipdata_.cliquetable.get_presolve_flag(),
            ms.mipdata_.cliques_extracted,
        )
    }

    /// On (HPresolve.cpp 2039-2312 at 7b890c255): runProbing sets the flag for its own
    /// scope and clears it on exit, whatever it was before. Off (1.15.1): runProbing leaves
    /// the flag that `presolve()` set.
    #[test]
    fn x1b_clqp_e_probing_flag_scoped() {
        let (result, after, extracted) = probe(true, true);
        assert_eq!(result, HPresolveResult::kOk);
        assert!(!after);
        assert!(extracted);
        assert!(!probe(true, false).1);

        let (result, after, extracted) = probe(false, true);
        assert_eq!(result, HPresolveResult::kOk);
        assert!(after);
        assert!(extracted);
        assert!(!probe(false, false).1);
    }
}
