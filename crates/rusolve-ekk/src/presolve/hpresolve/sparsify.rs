//! HiGHS source: presolve/HPresolve.cpp (lines 7962-8591)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Equality-row addition, variable bounds, sparsify, final checks. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

/// HiGHS 1.16: HighsImplications.h::HighsImplications::VarBound::minValue (lines 58-63)
///
/// Inf-safe form used by `aggregate_var_bounds` (#3562); the crate-wide
/// `VarBound::min_value` is the 1.15.1 form without the infinity guard.
fn var_bound_min_value_116(vb: &HighsImplications_VarBound) -> f64 {
    let m = cpp_min(vb.coef, 0.0);
    if vb.constant.abs() >= kHighsInf || m.abs() >= kHighsInf {
        return vb.constant + m;
    }
    f64::from(HighsCDouble::from(vb.constant) + m)
}

/// HiGHS 1.16: HighsImplications.h::HighsImplications::VarBound::maxValue (lines 64-69)
fn var_bound_max_value_116(vb: &HighsImplications_VarBound) -> f64 {
    let m = cpp_max(vb.coef, 0.0);
    if vb.constant.abs() >= kHighsInf || m.abs() >= kHighsInf {
        return vb.constant + m;
    }
    f64::from(HighsCDouble::from(vb.constant) + m)
}

impl HPresolve {
    /// HiGHS latest (#3264, f2ae2b320): HPresolve.cpp::HPresolve::aggregateVarBounds (lines 9787-9960)
    /// HiGHS 1.16 (#3562, #3596): HPresolveVarBndAggregation.cpp::aggregateVarBounds (lines 118-327)
    ///
    /// Only reached from `mip_probing.rs` under `mip_aggregate_implications`, so the 1.16
    /// statements are the on path and need no separate gate.
    pub fn aggregate_var_bounds(&mut self) {
        if self.mipsolver.is_none() {
            return;
        }

        // data structures
        let mut vlbs: HighsHashTree<
            HighsInt,
            (HighsImplications_VarBound, HighsImplications_VarBound),
        > = HighsHashTree::default();
        let mut vubs: HighsHashTree<
            HighsInt,
            (HighsImplications_VarBound, HighsImplications_VarBound),
        > = HighsHashTree::default();
        let mut vlbs_clique: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
        let mut vubs_clique: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
        let mut vlbs_cover: Vec<Vec<HighsCliqueTable_CliqueVar>> = Vec::new();
        let mut vubs_cover: Vec<Vec<HighsCliqueTable_CliqueVar>> = Vec::new();
        let mut row_consumed = vec![false; self.model.num_row_ as usize];

        // initialise counters
        let mut num_rows_removed: HighsInt = 0;
        let mut num_rows_modified: HighsInt = 0;
        let mut num_vars_lifted: HighsInt = 0;

        for col in 0..self.model.num_col_ {
            // skip deleted and binary columns (aggregate VLBs/VUBs on non-binary
            // columns)
            if self.col_deleted[col as usize] != 0 || self.is_binary(col) {
                continue;
            }

            // clear vectors
            vlbs.clear();
            vubs.clear();
            vlbs_clique.clear();
            vubs_clique.clear();

            // get lower bound and upper bound
            let mut lb = self.model.col_lower_[col as usize];
            let mut ub = self.model.col_upper_[col as usize];
            let feastol = self.mip().mipdata_.feastol;

            // the minimum (maximum) of every VLB (VUB) is a valid lower (upper) bound,
            // so tighten lb (ub) such that the standardized forms below imply the VLBs
            // (VUBs) exactly; the column bounds in the model are not changed
            self.mip().mipdata_.implications.get_vlbs(col).for_each(
                &mut |binary_col: &HighsInt, vlb: &HighsImplications_VarBound| {
                    if self.col_deleted[*binary_col as usize] == 0 {
                        lb = cpp_max(lb, var_bound_min_value_116(vlb));
                    }
                },
            );
            self.mip().mipdata_.implications.get_vubs(col).for_each(
                &mut |binary_col: &HighsInt, vub: &HighsImplications_VarBound| {
                    if self.col_deleted[*binary_col as usize] == 0 {
                        ub = cpp_min(ub, var_bound_max_value_116(vub));
                    }
                },
            );

            if self.model.integrality_[col as usize] != HighsVarType::kContinuous {
                lb = (lb - feastol).ceil();
                ub = (ub + feastol).floor();
            }

            // nothing to aggregate if the column has no VLBs and VUBs
            if lb <= -kHighsInf && ub >= kHighsInf {
                continue;
            }

            // nothing to aggregate if the tightened bounds fix the column or cross; the
            // latter indicates infeasibility, which is left to other reductions
            if lb >= ub {
                continue;
            }

            // compute range
            let mut range = kHighsInf;
            if lb > -kHighsInf && ub < kHighsInf {
                range = f64::from(HighsCDouble::from(ub) - lb);
            }

            // number of VLBs (VUBs) whose origin row could be rewritten
            let mut num_vlb_origins: HighsInt = 0;
            let mut num_vub_origins: HighsInt = 0;

            // collect VLBs: y >= coef * x + constant
            // standardized form: y >= a * x + lb, where a > 0
            if lb > -kHighsInf {
                let mut vlb_entries: Vec<(HighsInt, HighsImplications_VarBound)> = Vec::new();
                self.mip().mipdata_.implications.get_vlbs(col).for_each(
                    &mut |binary_col: &HighsInt, vlb: &HighsImplications_VarBound| {
                        vlb_entries.push((*binary_col, *vlb));
                    },
                );
                for (binary_col, vlb) in vlb_entries {
                    // skip deleted cols
                    if self.col_deleted[binary_col as usize] != 0 {
                        continue;
                    }
                    // skip if the VLB is dominated by the global lower bound
                    if self
                        .mip()
                        .mipdata_
                        .implications
                        .redundant_vlb(self.mip(), &vlb, lb)
                    {
                        continue;
                    }
                    // tighten so that minValue() >= lb
                    let mut v = vlb;
                    self.mip()
                        .mipdata_
                        .implications
                        .tighten_vlb(self.mip(), &mut v, lb);
                    // standardize: a = maxValue() - lb
                    let mut new_coef = v.constant - HighsCDouble::from(lb);
                    if v.coef > 0.0 {
                        new_coef += v.coef;
                    }
                    // skip if the standardized coefficient exceeds the variable range
                    if new_coef > range + feastol {
                        continue;
                    }
                    let standard_bound = HighsImplications_VarBound {
                        coef: f64::from(new_coef),
                        constant: lb,
                        origin: v.origin,
                    };
                    vlbs.insert(HighsHashTableEntry::new(binary_col, (v, standard_bound)));
                    if v.origin >= 0 {
                        num_vlb_origins += 1;
                    }
                }
            }
            // collect VUBs: y <= coef * x + constant
            // standardized form: y <= ub - a * x, where a > 0
            if ub < kHighsInf {
                let mut vub_entries: Vec<(HighsInt, HighsImplications_VarBound)> = Vec::new();
                self.mip().mipdata_.implications.get_vubs(col).for_each(
                    &mut |binary_col: &HighsInt, vub: &HighsImplications_VarBound| {
                        vub_entries.push((*binary_col, *vub));
                    },
                );
                for (binary_col, vub) in vub_entries {
                    // skip deleted cols
                    if self.col_deleted[binary_col as usize] != 0 {
                        continue;
                    }
                    // skip if the VUB is dominated by the global upper bound
                    if self
                        .mip()
                        .mipdata_
                        .implications
                        .redundant_vub(self.mip(), &vub, ub)
                    {
                        continue;
                    }
                    // tighten so that maxValue() <= ub
                    let mut v = vub;
                    self.mip()
                        .mipdata_
                        .implications
                        .tighten_vub(self.mip(), &mut v, ub);
                    // standardize: a = ub - minValue()
                    let mut new_coef = HighsCDouble::from(ub) - v.constant;
                    if v.coef < 0.0 {
                        new_coef -= v.coef;
                    }
                    // skip if the standardized coefficient exceeds the variable range
                    if new_coef > range + feastol {
                        continue;
                    }
                    let standard_bound = HighsImplications_VarBound {
                        coef: f64::from(new_coef),
                        constant: ub,
                        origin: v.origin,
                    };
                    vubs.insert(HighsHashTableEntry::new(binary_col, (v, standard_bound)));
                    if v.origin >= 0 {
                        num_vub_origins += 1;
                    }
                }
            }

            // only bounds with an origin row can be merged, so there is nothing to do
            // if no VLB and no VUB has one
            if num_vlb_origins == 0 && num_vub_origins == 0 {
                continue;
            }

            // set up cliques for the sides with an origin row
            if num_vlb_origins > 0 {
                vlbs.for_each(&mut |binary_col: &HighsInt,
                                    bounds: &(
                    HighsImplications_VarBound,
                    HighsImplications_VarBound,
                )| {
                    let val: HighsInt = if bounds.0.coef > 0.0 { 1 } else { 0 };
                    vlbs_clique.push(HighsCliqueTable_CliqueVar::from_parts(
                        *binary_col as HighsUInt,
                        val as HighsUInt,
                    ));
                });
            }
            if num_vub_origins > 0 {
                vubs.for_each(&mut |binary_col: &HighsInt,
                                    bounds: &(
                    HighsImplications_VarBound,
                    HighsImplications_VarBound,
                )| {
                    let val: HighsInt = if bounds.0.coef < 0.0 { 1 } else { 0 };
                    vubs_clique.push(HighsCliqueTable_CliqueVar::from_parts(
                        *binary_col as HighsUInt,
                        val as HighsUInt,
                    ));
                });
            }

            // clique cover (clears the covers, which stay empty without cliques)
            self.mip_mut()
                .mipdata_
                .cliquetable
                .clique_cover(&mut vlbs_clique, &mut vlbs_cover);
            self.mip_mut()
                .mipdata_
                .cliquetable
                .clique_cover(&mut vubs_clique, &mut vubs_cover);

            // mergeCliques for VLBs
            if lb > -kHighsInf {
                for clique in vlbs_cover.iter() {
                    if clique.len() < 2 {
                        continue;
                    }

                    // find an unconsumed row to reuse and remove the rest
                    let mut row: HighsInt = -1;
                    for var in clique.iter() {
                        let bounds = vlbs.find(&(var.col() as HighsInt));
                        let currentrow = bounds.unwrap().0.origin;
                        if currentrow < 0 {
                            continue;
                        }
                        if !row_consumed[currentrow as usize] {
                            row_consumed[currentrow as usize] = true;
                            if row == -1 {
                                row = currentrow;
                            } else {
                                self.remove_row(currentrow);
                                num_rows_removed += 1;
                            }
                        }
                    }
                    if row == -1 {
                        continue;
                    }

                    // origin is invalidated when row bounds change, so the
                    // row must still be one-sided
                    debug_assert!(!self.is_ranged(row));

                    // rewrite the reused row with the aggregated constraint
                    self.unlink_row(row);
                    self.add_to_matrix(row, col, 1.0);
                    num_rows_modified += 1;

                    let mut row_bound = HighsCDouble::from(lb);
                    for var in clique.iter() {
                        let bounds = vlbs.find(&(var.col() as HighsInt));
                        let a = bounds.unwrap().1.coef;
                        num_vars_lifted += 1;
                        let direction: HighsInt = 1;
                        if var.val() == 1 {
                            self.add_to_matrix(row, var.col() as HighsInt, -(direction as f64) * a);
                        } else {
                            self.add_to_matrix(row, var.col() as HighsInt, (direction as f64) * a);
                            row_bound += (direction as f64) * a;
                        }
                    }

                    // direction > 0
                    self.h16_change_row_lower(row, f64::from(row_bound), false);
                    self.h16_change_row_upper(row, kHighsInf, false);
                }
            }

            // mergeCliques for VUBs
            if ub < kHighsInf {
                for clique in vubs_cover.iter() {
                    if clique.len() < 2 {
                        continue;
                    }

                    // find an unconsumed row to reuse and remove the rest
                    let mut row: HighsInt = -1;
                    for var in clique.iter() {
                        let bounds = vubs.find(&(var.col() as HighsInt));
                        let currentrow = bounds.unwrap().0.origin;
                        if currentrow < 0 {
                            continue;
                        }
                        if !row_consumed[currentrow as usize] {
                            row_consumed[currentrow as usize] = true;
                            if row == -1 {
                                row = currentrow;
                            } else {
                                self.remove_row(currentrow);
                                num_rows_removed += 1;
                            }
                        }
                    }
                    if row == -1 {
                        continue;
                    }

                    // origin is invalidated when row bounds change, so the
                    // row must still be one-sided
                    debug_assert!(!self.is_ranged(row));

                    // rewrite the reused row with the aggregated constraint
                    self.unlink_row(row);
                    self.add_to_matrix(row, col, 1.0);
                    num_rows_modified += 1;

                    let mut row_bound = HighsCDouble::from(ub);
                    for var in clique.iter() {
                        let bounds = vubs.find(&(var.col() as HighsInt));
                        let a = bounds.unwrap().1.coef;
                        num_vars_lifted += 1;
                        let direction: HighsInt = -1;
                        if var.val() == 1 {
                            self.add_to_matrix(row, var.col() as HighsInt, -(direction as f64) * a);
                        } else {
                            self.add_to_matrix(row, var.col() as HighsInt, (direction as f64) * a);
                            row_bound += (direction as f64) * a;
                        }
                    }

                    // direction < 0
                    self.h16_change_row_lower(row, -kHighsInf, false);
                    self.h16_change_row_upper(row, f64::from(row_bound), false);
                }
            }
        }

        if num_rows_removed > 0 || num_rows_modified > 0 {
            crate::highs_log_dev!(
                &self.options.log_options,
                HighsLogType::kInfo,
                "Implied variable bound aggregation: {} rows removed, {} rows modified, {} vars lifted\n",
                num_rows_removed,
                num_rows_modified,
                num_vars_lifted
            );
        }
    }
    /// HiGHS: HPresolve.cpp::HPresolve::equalityRowAddition (lines 7962-8000, 39 lines)
    ///
    /// C++ template over `HighsMatrixSlice<RowStorageFormat>`.
    /// `postsolve_stack.equalityRowAddition(removerow, stayrow, scale, rowvector)` ->
    /// `postsolve_stack.equality_row_addition(removerow, stayrow, scale, rowvector.iter(self))`.
    /// `static_cast<HighsCDouble>(scale) * x` -> `HighsCDouble::from(scale) * x`.
    pub fn equality_row_addition<S: HighsTripletSlice>(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        stayrow: HighsInt,
        removerow: HighsInt,
        scale: f64,
        rowvector: S,
    ) -> HPresolveResult {
        // extract non-zero positions
        let mut stay_rowpositions: Vec<HighsInt> = Vec::new();
        self.get_row_positions(stayrow, &mut stay_rowpositions);

        // update postsolve information
        postsolve_stack.equality_row_addition(removerow, stayrow, scale, rowvector.iter(self));

        // iterate over non-zero positions instead of iterating over the
        // HighsMatrixSlice because the latter contains pointers to Acol and Avalue
        // that may be invalidated if these vectors are reallocated
        // (see std::vector::push_back performed in HPresolve::addToMatrix).
        for &rowiter in stay_rowpositions.iter() {
            let stay_col = self.a_col[rowiter as usize];
            let pos = self.find_nonzero(removerow, stay_col);
            if pos != -1 {
                // all common nonzeros are cancelled, as the rows are parallel
                self.unlink(pos);
            } else {
                // might introduce a singleton
                let val = scale * self.a_value[rowiter as usize];
                self.add_to_matrix(removerow, stay_col, val);
            }
        }

        if self.options.h16_row_duals {
            // HiGHS 1.16 #3270 (h16_row_duals): HPresolveParallelRowsAndCols.cpp lines 749-752 (1.16.0)
            // NOTE: the C++ reads `row_upper_[stayrow]` for both sides (kept as in HiGHS).
            let shift = HighsCDouble::from(scale) * self.model.row_upper_[stayrow as usize];
            self.h16_add_to_row_lower(removerow, shift);
            let shift = HighsCDouble::from(scale) * self.model.row_upper_[stayrow as usize];
            self.h16_add_to_row_upper(removerow, shift);
            // HiGHS 1.16 #3367 (h16_row_duals): HPresolveParallelRowsAndCols.cpp lines 753-759 (1.16.0)
            // equalityRowAddition changes the size of removerow, so it must be re-keyed in
            // the size-ordered equations set
            self.reinsert_equation(removerow);
        } else {
            if self.model.row_upper_[removerow as usize] != kHighsInf {
                let v = self.model.row_upper_[removerow as usize]
                    + HighsCDouble::from(scale) * self.model.row_upper_[stayrow as usize];
                self.model.row_upper_[removerow as usize] = f64::from(v);
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 9672 at f2ae2b320 (1.15.1 7987)
                self.x1_row_modified(removerow);
            }
            if self.model.row_lower_[removerow as usize] != -kHighsInf {
                // NOTE: the C++ reads `row_upper_[stayrow]` here too (kept as in HiGHS).
                let v = self.model.row_lower_[removerow as usize]
                    + HighsCDouble::from(scale) * self.model.row_upper_[stayrow as usize];
                self.model.row_lower_[removerow as usize] = f64::from(v);
                // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 9670 at f2ae2b320 (1.15.1 7991)
                self.x1_row_modified(removerow);
            }
        }

        // row is now a singleton row, doubleton equation, or a row
        // that contains only singletons and we let the normal row presolve
        // handle the cases
        hpresolve_checked_call!(self.row_presolve(postsolve_stack, removerow));
        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::setRelaxedImpliedBounds (lines 8002-8037, 36 lines)
    pub fn set_relaxed_implied_bounds(&mut self) {
        let huge_bound = self.primal_feastol / kHighsTiny;
        for i in 0..self.model.num_col_ {
            let iu = i as usize;
            if self.model.col_lower_[iu] >= self.impl_col_lower[iu]
                && self.model.col_upper_[iu] <= self.impl_col_upper[iu]
            {
                continue;
            }

            if self.impl_col_lower[iu].abs() <= huge_bound {
                // if the bound is derived from a small nonzero value
                // then we want to increase the margin so that we make sure
                // the row it was derived from is violated if the column sits
                // at this relaxed bound in the final solution.
                let source_row = self.col_lower_source[iu];
                let nz_pos = self.find_nonzero(source_row, i);

                let bound_relax = cpp_max(1000.0, self.impl_col_lower[iu].abs())
                    * self.primal_feastol
                    / cpp_min(1.0, self.a_value[nz_pos as usize].abs());

                let new_lb = self.impl_col_lower[iu] - bound_relax;
                if new_lb > self.model.col_lower_[iu] + bound_relax {
                    self.model.col_lower_[iu] = new_lb;
                }
            }

            if self.impl_col_upper[iu].abs() <= huge_bound {
                let source_row = self.col_upper_source[iu];
                let nz_pos = self.find_nonzero(source_row, i);

                let bound_relax = cpp_max(1000.0, self.impl_col_upper[iu].abs())
                    * self.primal_feastol
                    / cpp_min(1.0, self.a_value[nz_pos as usize].abs());

                let new_ub = self.impl_col_upper[iu] + bound_relax;
                if new_ub < self.model.col_upper_[iu] - bound_relax {
                    self.model.col_upper_[iu] = new_ub;
                }
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::extractVarBounds (lines 8039-8137, 99 lines)
    ///
    /// Returns at once on LP (`mipsolver == nullptr || ...`); translate the rest anyway.
    pub fn extract_var_bounds(&mut self, row: HighsInt) {
        // return if row is empty or a singleton or contains no integer variables
        if self.mipsolver.is_none()
            || self.rowsize[row as usize] <= 1
            || self.rowsize_integer[row as usize] == 0
        {
            return;
        }

        // check if variable bounds can be derived from row
        let num_inf_sum_lower = self.implied_row_bounds.get_num_inf_sum_lower(row);
        let num_inf_sum_upper = self.implied_row_bounds.get_num_inf_sum_upper(row);
        let mut use_lhs =
            self.model.row_lower_[row as usize] != -kHighsInf && num_inf_sum_upper <= 1;
        let mut use_rhs =
            self.model.row_upper_[row as usize] != kHighsInf && num_inf_sum_lower <= 1;
        if !use_lhs && !use_rhs {
            return;
        }

        let flag = self.options.mip_aggregate_implications;

        // check if row contains a single binary variable
        let mut bin_col: HighsInt = -1;
        let mut bin_coef: f64 = 0.0;
        for nonzero in self.get_row_vector(row).iter(self) {
            // skip fixed variables
            if self.model.col_lower_[nonzero.index() as usize]
                == self.model.col_upper_[nonzero.index() as usize]
            {
                continue;
            }

            // find binary variable
            if self.model.integrality_[nonzero.index() as usize] == HighsVarType::kInteger
                && self.model.col_lower_[nonzero.index() as usize] == 0.0
                && self.model.col_upper_[nonzero.index() as usize] == 1.0
            {
                // return if there is more than one binary variable
                if bin_col != -1 {
                    return;
                }
                bin_col = nonzero.index();
                bin_coef = nonzero.value();
            }
        }
        // return if there is no binary variable
        if bin_col == -1 {
            return;
        }

        let mut it = self.get_row_vector(row);
        while let Some(nonzero) = it.next_nz(self) {
            // skip fixed variables
            if self.model.col_lower_[nonzero.index() as usize]
                == self.model.col_upper_[nonzero.index() as usize]
            {
                continue;
            }

            // skip binary variable
            if nonzero.index() == bin_col {
                continue;
            }

            // compute VLB constant
            let mut vlb_constant = -kHighsInf;
            if use_lhs {
                let residual = self.implied_row_bounds.get_residual_sum_upper(
                    &col_bound_arrays!(self),
                    row,
                    nonzero.index(),
                    nonzero.value(),
                    bin_col,
                    bin_coef,
                    0.0,
                );
                if residual != kHighsInf {
                    vlb_constant = f64::from(
                        (HighsCDouble::from(self.model.row_lower_[row as usize]) - residual)
                            / nonzero.value().abs(),
                    );
                    use_lhs = num_inf_sum_upper == 0;
                }
            }

            // compute VUB constant
            let mut vub_constant = kHighsInf;
            if use_rhs {
                let residual = self.implied_row_bounds.get_residual_sum_lower(
                    &col_bound_arrays!(self),
                    row,
                    nonzero.index(),
                    nonzero.value(),
                    bin_col,
                    bin_coef,
                    0.0,
                );
                if residual != -kHighsInf {
                    vub_constant = f64::from(
                        (HighsCDouble::from(self.model.row_upper_[row as usize]) - residual)
                            / nonzero.value().abs(),
                    );
                    use_rhs = num_inf_sum_lower == 0;
                }
            }

            // switch sign if continuous variable has a negative coefficient
            if nonzero.value() < 0.0 {
                vlb_constant *= -1.0;
                vub_constant *= -1.0;
                std::mem::swap(&mut vlb_constant, &mut vub_constant);
            }

            // compute coefficient for binary variable
            let vb_coef = -bin_coef / nonzero.value();

            // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 9763-9764, 9768-9782 at f2ae2b320
            let origin = if flag && self.rowsize[row as usize] == 2 && !self.is_ranged(row) {
                row
            } else {
                -1
            };

            // add VLB
            if vlb_constant != -kHighsInf {
                let col = nonzero.index();
                let col_lower = self.model.col_lower_[col as usize];
                let integral = self.model.integrality_[col as usize] != HighsVarType::kContinuous;
                if flag {
                    with_mip_solver!(self, ms => crate::detach!(ms.mipdata_.implications, |imp| imp.add_vlb_2_o(
                        ms, col, bin_col, vb_coef, vlb_constant, col_lower, integral, origin
                    )));
                } else {
                    with_mip_solver!(self, ms => crate::detach!(ms.mipdata_.implications, |imp| imp.add_vlb_2(
                        ms, col, bin_col, vb_coef, vlb_constant, col_lower, integral
                    )));
                }
            }

            // add VUB
            if vub_constant != kHighsInf {
                let col = nonzero.index();
                let col_upper = self.model.col_upper_[col as usize];
                let integral = self.model.integrality_[col as usize] != HighsVarType::kContinuous;
                if flag {
                    with_mip_solver!(self, ms => crate::detach!(ms.mipdata_.implications, |imp| imp.add_vub_2_o(
                        ms, col, bin_col, vb_coef, vub_constant, col_upper, integral, origin
                    )));
                } else {
                    with_mip_solver!(self, ms => crate::detach!(ms.mipdata_.implications, |imp| imp.add_vub_2(
                        ms, col, bin_col, vb_coef, vub_constant, col_upper, integral
                    )));
                }
            }

            // stop if no additional variable bounds can be found
            if !use_lhs && !use_rhs {
                break;
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::sparsify (lines 8301-8546, 246 lines)
    ///
    /// `std::map<double, HighsInt> possibleScales` -> `BTreeMap<OrderedF64, HighsInt>`;
    /// `possibleScales.lower_bound(x)` -> `possible_scales.range(OrderedF64(x)..)` (first element).
    pub fn sparsify(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::SPARSIFY);
        debug_assert!(self.analysis_.allow_rule_[kPresolveRuleSparsify as usize]);
        let mut sparsify_rows: Vec<Nonzero> = Vec::new();
        let logging_on = self.analysis_.logging_on_;
        if logging_on {
            self.analysis_start_presolve_rule_log(kPresolveRuleSparsify);
        }
        hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
        hpresolve_checked_call!(self.remove_doubleton_equations(postsolve_stack));
        let mut tmp_equations: Vec<HighsInt> = Vec::with_capacity(self.equations.len());

        let min_nonzero_val = self.primal_feastol.sqrt();

        for eq in self.equations.iter() {
            tmp_equations.push(eq.1);
        }
        for &eqrow in tmp_equations.iter() {
            if self.row_deleted[eqrow as usize] != 0 {
                continue;
            }

            debug_assert!(self.row_deleted[eqrow as usize] == 0);
            debug_assert!(self.is_equation(eqrow));

            self.store_row(eqrow);

            let mut sparsest_col: HighsInt = -1;
            let mut second_sparsest_col: HighsInt = -1;
            let mut sparsest_col_len: HighsInt = kHighsIInf;
            let mut second_sparsest_col_len: HighsInt = kHighsIInf;
            for &nz_pos in self.rowpositions.iter() {
                let col = self.a_col[nz_pos as usize];
                if self.colsize[col as usize] < sparsest_col_len {
                    second_sparsest_col = sparsest_col;
                    second_sparsest_col_len = sparsest_col_len;
                    sparsest_col = col;
                    sparsest_col_len = self.colsize[col as usize];
                } else if self.colsize[col as usize] < second_sparsest_col_len {
                    second_sparsest_col = col;
                    second_sparsest_col_len = self.colsize[col as usize];
                }
            }

            debug_assert!(sparsest_col != -1 && second_sparsest_col != -1);
            debug_assert!(
                self.colsize[sparsest_col as usize] <= self.colsize[second_sparsest_col as usize]
            );

            let mut possible_scales: BTreeMap<OrderedF64, HighsInt> = BTreeMap::new();
            sparsify_rows.clear();

            // cursor over the column; the body may call `find_nonzero` (`&mut self`)
            let mut col_cursor = self.get_column_vector(sparsest_col);
            while let Some(col_nz) = col_cursor.next_nz(self) {
                let cand_row = col_nz.index();
                if cand_row == eqrow {
                    continue;
                }

                possible_scales.clear();

                let mut misses: HighsInt = 0;
                // allow no fillin if a completely continuous row is used to cancel a row
                // that has integers as there are instances where this leads to a huge
                // deterioration of cut performance
                let mut max_misses: HighsInt = 1;
                if self.rowsize_integer[eqrow as usize] == 0
                    && self.rowsize_integer[cand_row as usize] != 0
                {
                    max_misses -= 1;
                }
                let mut stored = self.get_stored_row();
                while let Some(nonzero) = stored.next_nz(self) {
                    #[allow(clippy::needless_late_init)]
                    let cand_row_val: f64;
                    if nonzero.index() == sparsest_col {
                        cand_row_val = col_nz.value();
                    } else {
                        let nz_pos = self.find_nonzero(cand_row, nonzero.index());
                        if nz_pos == -1 {
                            if self.model.integrality_[nonzero.index() as usize]
                                == HighsVarType::kInteger
                                && self.model.col_upper_[nonzero.index() as usize]
                                    - self.model.col_lower_[nonzero.index() as usize]
                                    > 1.5
                            {
                                // do not allow fillin of general integers
                                misses = 2;
                                break;
                            }
                            misses += 1;
                            if misses > max_misses {
                                break;
                            }
                            continue;
                        }
                        cand_row_val = self.a_value[nz_pos as usize];
                    }

                    let scale = -cand_row_val / nonzero.value();
                    if scale.abs() > 1e3 {
                        continue;
                    }

                    let scale_tolerance = min_nonzero_val / nonzero.value().abs();
                    let lower = possible_scales
                        .range(OrderedF64(scale - scale_tolerance)..)
                        .next()
                        .map(|(k, _)| *k);
                    match lower {
                        Some(key) if (key.0 - scale).abs() <= scale_tolerance => {
                            // there already is a scale that is very close and could produce
                            // a matrix value for this nonzero that is below the allowed
                            // threshold. Therefore we check if the matrix value is small enough
                            // for this nonzero to be deleted, in which case the number of
                            // deleted nonzeros for the other scale is increased. If it is not
                            // small enough we do not use this scale or the other one because
                            // such small matrix values may lead to numerical troubles.
                            let count = possible_scales.get_mut(&key).unwrap();

                            // scale is already marked to be numerically bad
                            if *count == -1 {
                                continue;
                            }

                            if (key.0 * nonzero.value() + cand_row_val).abs()
                                <= self.options.small_matrix_value
                            {
                                *count += 1;
                            } else {
                                *count = -1;
                            }
                        }
                        _ => {
                            possible_scales.entry(OrderedF64(scale)).or_insert(1);
                        }
                    }
                }

                if misses > max_misses || possible_scales.is_empty() {
                    continue;
                }

                let mut num_cancel: HighsInt = 0;
                let mut scale: f64 = 0.0;

                for (s_key, s_val) in possible_scales.iter() {
                    if *s_val <= misses {
                        continue;
                    }

                    if *s_val > num_cancel || (*s_val == num_cancel && s_key.0.abs() < scale.abs())
                    {
                        scale = s_key.0;
                        num_cancel = *s_val;
                    }
                }

                debug_assert!(scale != 0.0 || num_cancel == 0);

                // cancels at least one nonzero if the scale cancels more than there is
                // fillin
                if num_cancel > misses {
                    sparsify_rows.push(Nonzero::new(cand_row, scale));
                }
            }

            if self.model.integrality_[sparsest_col as usize] != HighsVarType::kInteger
                || (self.model.col_upper_[sparsest_col as usize]
                    - self.model.col_lower_[sparsest_col as usize])
                    < 1.5
            {
                // now check for rows which do not contain the sparsest column but all
                // other columns by scanning the second sparsest column
                let mut col_cursor = self.get_column_vector(second_sparsest_col);
                while let Some(col_nz) = col_cursor.next_nz(self) {
                    let cand_row = col_nz.index();
                    if cand_row == eqrow {
                        continue;
                    }

                    if self.rowsize_integer[eqrow as usize] == 0
                        && self.rowsize_integer[cand_row as usize] != 0
                    {
                        continue;
                    }

                    let sparsest_col_pos = self.find_nonzero(cand_row, sparsest_col);

                    // if the row has a nonzero for the sparsest column we have already
                    // checked it
                    if sparsest_col_pos != -1 {
                        continue;
                    }

                    possible_scales.clear();
                    let mut skip = false;
                    let mut stored = self.get_stored_row();
                    while let Some(nonzero) = stored.next_nz(self) {
                        #[allow(clippy::needless_late_init)]
                        let cand_row_val: f64;
                        if nonzero.index() == second_sparsest_col {
                            cand_row_val = col_nz.value();
                        } else {
                            let nz_pos = self.find_nonzero(cand_row, nonzero.index());
                            // we already have a miss for the sparsest column, so with another
                            // one we want to skip the row
                            skip = nz_pos == -1;
                            if skip {
                                break;
                            }

                            cand_row_val = self.a_value[nz_pos as usize];
                        }

                        let scale = -cand_row_val / nonzero.value();
                        if scale.abs() > 1e3 {
                            continue;
                        }

                        let scale_tolerance = min_nonzero_val / nonzero.value().abs();
                        let lower = possible_scales
                            .range(OrderedF64(scale - scale_tolerance)..)
                            .next()
                            .map(|(k, _)| *k);
                        match lower {
                            Some(key) if (key.0 - scale).abs() <= scale_tolerance => {
                                // there already is a scale that is very close (see above)
                                let count = possible_scales.get_mut(&key).unwrap();

                                // scale is already marked to be numerically bad
                                if *count == -1 {
                                    continue;
                                }

                                if (key.0 * nonzero.value() + cand_row_val).abs()
                                    <= self.options.small_matrix_value
                                {
                                    *count += 1;
                                } else {
                                    // mark scale to be numerically bad
                                    *count = -1;
                                    continue;
                                }
                            }
                            _ => {
                                possible_scales.entry(OrderedF64(scale)).or_insert(1);
                            }
                        }
                    }

                    if skip || possible_scales.is_empty() {
                        continue;
                    }

                    let mut num_cancel: HighsInt = 0;
                    let mut scale: f64 = 0.0;

                    for (s_key, s_val) in possible_scales.iter() {
                        if *s_val <= 1 {
                            continue;
                        }
                        if *s_val > num_cancel
                            || (*s_val == num_cancel && s_key.0.abs() < scale.abs())
                        {
                            scale = s_key.0;
                            num_cancel = *s_val;
                        }
                    }

                    debug_assert!(scale != 0.0 || num_cancel == 0);

                    // cancels at least one nonzero if the scale cancels more than there is
                    // fillin
                    if num_cancel > 1 {
                        sparsify_rows.push(Nonzero::new(cand_row, scale));
                    }
                }
            }

            if sparsify_rows.is_empty() {
                continue;
            }

            let stored_row_vec = self.get_stored_row().iter(self);
            // HiGHS 1.16 #3324 (h16_eq_row_adds): HPresolveSparsify.cpp line 233 (1.16.0)
            let h16_eq_row_adds = self.options.h16_eq_row_adds;
            postsolve_stack.equality_row_additions(
                eqrow,
                stored_row_vec,
                &sparsify_rows,
                h16_eq_row_adds,
            );
            let rhs = self.model.row_lower_[eqrow as usize];
            let h16_row_duals = self.options.h16_row_duals;
            for sparsify_row in sparsify_rows.iter() {
                let row = sparsify_row.index;
                let scale = sparsify_row.value;

                if h16_row_duals {
                    // HiGHS 1.16 #3270 (h16_row_duals): HPresolveSparsify.cpp lines 239-240 (1.16.0)
                    self.h16_add_to_row_lower(row, HighsCDouble::from(scale) * rhs);
                    self.h16_add_to_row_upper(row, HighsCDouble::from(scale) * rhs);
                } else {
                    if self.model.row_lower_[row as usize] != -kHighsInf {
                        self.model.row_lower_[row as usize] += scale * rhs;
                        // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 10187 at f2ae2b320 (1.15.1 8526)
                        self.x1_row_modified(row);
                    }

                    if self.model.row_upper_[row as usize] != kHighsInf {
                        self.model.row_upper_[row as usize] += scale * rhs;
                        // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 10188 at f2ae2b320 (1.15.1 8529)
                        self.x1_row_modified(row);
                    }
                }

                // C++ range-for over `rowpositions`: the length is fixed at loop start
                let n = self.rowpositions.len();
                for k in 0..n {
                    let pos = self.rowpositions[k] as usize;
                    let col = self.a_col[pos];
                    let val = scale * self.a_value[pos];
                    self.add_to_matrix(row, col, val);
                }

                self.reinsert_equation(row);
            }

            hpresolve_checked_call!(self.check_limits(postsolve_stack));
            hpresolve_checked_call!(self.remove_row_singletons(postsolve_stack));
            hpresolve_checked_call!(self.remove_doubleton_equations(postsolve_stack));
        }

        self.analysis_.logging_on_ = logging_on;
        if logging_on {
            self.analysis_stop_presolve_rule_log(kPresolveRuleSparsify);
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::zeroRowActivityFeasible (lines 8548-8557, 10 lines)
    pub fn zero_row_activity_feasible(&self) -> bool {
        // Check that zero row activity is feasible - called when reduced model
        // has no columns to assess whether the HighsModelStatus returned is
        // kOptimal or kInfeasible (as was required for 2326)
        for i_row in 0..self.model.num_row_ {
            let r = i_row as usize;
            if self.model.row_lower_[r] > self.primal_feastol
                || self.model.row_upper_[r] < -self.primal_feastol
            {
                return false;
            }
        }
        true
    }

    /// HiGHS: HPresolve.cpp::HPresolve::debugGetCheckCol (lines 8559-8574, 16 lines)
    ///
    /// `check_col_name` is the constant `""`, so the C++ returns `-1` before the name lookup;
    /// the dead lookup (which needs `HighsLp::col_hash_`) is not ported.
    pub fn debug_get_check_col(&self) -> HighsInt {
        let check_col: HighsInt = -1;
        // `if (check_col_name == "") return check_col;` always taken.
        check_col
    }

    /// HiGHS: HPresolve.cpp::HPresolve::debugGetCheckRow (lines 8576-8591, 16 lines)
    ///
    /// As `debug_get_check_col`.
    pub fn debug_get_check_row(&self) -> HighsInt {
        let check_row: HighsInt = -1;
        // `if (check_row_name == "") return check_row;` always taken.
        check_row
    }
}

#[cfg(test)]
mod x1_agg_d_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;

    /// C++ TestPresolveRules "test-implied-bound-aggregation" (VLB, min z) and
    /// "test-implied-bound-aggregation-vub" (VUB, max z), 4 columns (z, x0, x1, x2).
    fn build_lp(vubs: bool) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 4;
        lp.sense_ = if vubs {
            ObjSense::kMaximize
        } else {
            ObjSense::kMinimize
        };
        lp.col_cost_ = vec![1.0, 0.0, 0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0];
        lp.col_upper_ = if vubs {
            vec![30.0, 1.0, 1.0, 1.0]
        } else {
            vec![kHighsInf, 1.0, 1.0, 1.0]
        };
        lp.integrality_ = vec![
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = if vubs {
            vec![2.0, -kHighsInf, -kHighsInf, -kHighsInf]
        } else {
            vec![1.0, -kHighsInf, -kHighsInf, -kHighsInf]
        };
        lp.row_upper_ = if vubs {
            vec![2.0, 20.0, 10.0, 0.0]
        } else {
            vec![1.0, 0.0, 0.0, 0.0]
        };
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 5, 7, 9];
        lp.a_matrix_.index_ = vec![1, 2, 3, 0, 1, 0, 2, 0, 3];
        lp.a_matrix_.value_ = if vubs {
            vec![1.0, 1.0, 1.0, 1.0, -10.0, 1.0, -20.0, 1.0, -30.0]
        } else {
            vec![-1.0, -1.0, -1.0, 1.0, 10.0, 1.0, 20.0, 1.0, 30.0]
        };
        lp
    }

    /// Sets up MIP presolve as the C++ test does (mipdata_ init, setupDomainPropagation,
    /// addVLB / addVUB with origins 1, 2, 3, setInput), then runs what `presolve()` runs with
    /// `presolve_rule_test = kPresolveRuleProbing`: its prefix up to `presolveRuleTest`, then
    /// `presolveRuleTestProbing` (latest HPresolveTest.cpp 53-60, i.e. `runProbing`) and the
    /// `presolveReturn` lambda (`scaleMIP`). The 1.15.1 `presolve_rule_test` port has no
    /// probing branch, so the steps are called here. Returns, for every row other than the
    /// original row 0 that is left, its dense coefficients and bounds (no row has been
    /// shrunk away, so row indices are the original ones).
    fn presolve_rows(flag: bool, vubs: bool) -> Vec<(Vec<f64>, f64, f64)> {
        let lp = build_lp(vubs);
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.presolve_rule_test = kPresolveRuleProbing;
        options.mip_aggregate_implications = flag;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.timer_.start();
        ms.new_data();
        ms.init();
        ms.setup_domain_propagation();
        if vubs {
            crate::detach!(ms.mipdata_.implications, |imp| {
                imp.add_vub_2_o(&mut ms, 0, 1, 10.0, 20.0, 30.0, false, 1);
                imp.add_vub_2_o(&mut ms, 0, 2, 20.0, 10.0, 30.0, false, 2);
                imp.add_vub_2_o(&mut ms, 0, 3, 30.0, 0.0, 30.0, false, 3);
            });
        } else {
            crate::detach!(ms.mipdata_.implications, |imp| {
                imp.add_vlb_o(&mut ms, 0, 1, 10.0, 0.0, 1);
                imp.add_vlb_o(&mut ms, 0, 2, 20.0, 0.0, 2);
                imp.add_vlb_o(&mut ms, 0, 3, 30.0, 0.0, 3);
            });
        }
        let mut postsolve_stack = std::mem::take(&mut ms.mipdata_.post_solve_stack);

        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), -1));

        // HPresolve::presolve prefix (presolve_loop.rs) up to presolveRuleTest.
        if presolve.model.sense_ == ObjSense::kMaximize {
            for i in 0..presolve.model.num_col_ {
                presolve.model.col_cost_[i as usize] = -presolve.model.col_cost_[i as usize];
            }
            presolve.model.offset_ = -presolve.model.offset_;
            presolve.model.sense_ = ObjSense::kMinimize;
        }
        let silent = presolve.silent_log();
        presolve
            .analysis_
            .setup(&presolve.model, &presolve.options, silent);
        presolve
            .mip_mut()
            .mipdata_
            .cliquetable
            .set_presolve_flag(true);

        // presolveRuleTestProbing, then presolveReturn.
        assert!(presolve.run_probing(&mut postsolve_stack) == HPresolveResult::kOk);
        assert!(presolve.scale_mip(&mut postsolve_stack) == HPresolveResult::kOk);

        let mut rows = Vec::new();
        for row in 1..presolve.model.num_row_ {
            if presolve.row_deleted[row as usize] != 0 {
                continue;
            }
            let mut coefs = vec![0.0; presolve.model.num_col_ as usize];
            let mut it = presolve.get_row_vector(row);
            while let Some(nz) = it.next_nz(&presolve) {
                coefs[nz.index() as usize] = nz.value();
            }
            rows.push((
                coefs,
                presolve.model.row_lower_[row as usize],
                presolve.model.row_upper_[row as usize],
            ));
        }
        rows
    }

    #[test]
    fn x1_agg_d_vlb_aggregation() {
        // z - 10 x0 - 20 x1 - 30 x2 >= 0 is the only row left besides row 0.
        let rows = presolve_rows(true, false);
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert_eq!(rows[0].0, vec![1.0, -10.0, -20.0, -30.0]);
        assert_eq!(rows[0].1, 0.0);
        assert_eq!(rows[0].2, kHighsInf);
    }

    #[test]
    fn x1_agg_d_vub_aggregation() {
        // z - 10 x0 - 20 x1 - 30 x2 <= -30 is the only row left besides row 0.
        let rows = presolve_rows(true, true);
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert_eq!(rows[0].0, vec![1.0, -10.0, -20.0, -30.0]);
        assert_eq!(rows[0].1, -kHighsInf);
        assert_eq!(rows[0].2, -30.0);
    }

    #[test]
    fn x1_agg_d_off_no_aggregation() {
        // Option off: the three variable-bound rows stay as they are.
        for vubs in [false, true] {
            let rows = presolve_rows(false, vubs);
            assert_eq!(rows.len(), 3, "{rows:?}");
        }
    }
}

#[cfg(test)]
mod h16_agg_fix_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;

    /// C++ TestPresolveRules "test-implied-bound-aggregation-vlb-floor" (`vubs == false`,
    /// min z, lb of z raised to 5 by the VLB minima) and "...-vub-ceiling" (`vubs == true`,
    /// max w, ub of w lowered to -5 by the VUB maxima). Columns are (z, x0, x1, x2).
    fn build_lp(vubs: bool) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 4;
        lp.sense_ = if vubs {
            ObjSense::kMaximize
        } else {
            ObjSense::kMinimize
        };
        lp.col_cost_ = vec![1.0, 0.0, 0.0, 0.0];
        lp.col_lower_ = if vubs {
            vec![-1000.0, 0.0, 0.0, 0.0]
        } else {
            vec![0.0, 0.0, 0.0, 0.0]
        };
        lp.col_upper_ = if vubs {
            vec![0.0, 1.0, 1.0, 1.0]
        } else {
            vec![1000.0, 1.0, 1.0, 1.0]
        };
        lp.integrality_ = vec![
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = vec![-kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf];
        lp.row_upper_ = vec![1.0, 0.0, -5.0, 0.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 5, 7, 9];
        lp.a_matrix_.index_ = vec![1, 2, 3, 0, 1, 0, 2, 0, 3];
        lp.a_matrix_.value_ = if vubs {
            vec![1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 15.0, 1.0, 30.0]
        } else {
            vec![-1.0, -1.0, -1.0, 1.0, 10.0, 1.0, 15.0, 1.0, 30.0]
        };
        lp
    }

    /// Same set-up and presolve steps as `x1_agg_d_tests::presolve_rows`, with the VLB/VUB
    /// bounds of the C++ tests. Returns the dense coefficients and bounds of every row other
    /// than row 0 that is left.
    fn presolve_rows(flag: bool, vubs: bool) -> Vec<(Vec<f64>, f64, f64)> {
        let lp = build_lp(vubs);
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.presolve_rule_test = kPresolveRuleProbing;
        options.mip_aggregate_implications = flag;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.timer_.start();
        ms.new_data();
        ms.init();
        ms.setup_domain_propagation();
        if vubs {
            crate::detach!(ms.mipdata_.implications, |imp| {
                imp.add_vub_o(&mut ms, 0, 1, -10.0, 0.0, 1);
                imp.add_vub_o(&mut ms, 0, 2, -15.0, -5.0, 2);
                imp.add_vub_o(&mut ms, 0, 3, -30.0, 0.0, 3);
            });
        } else {
            crate::detach!(ms.mipdata_.implications, |imp| {
                imp.add_vlb_o(&mut ms, 0, 1, 10.0, 0.0, 1);
                imp.add_vlb_o(&mut ms, 0, 2, 15.0, 5.0, 2);
                imp.add_vlb_o(&mut ms, 0, 3, 30.0, 0.0, 3);
            });
        }
        let mut postsolve_stack = std::mem::take(&mut ms.mipdata_.post_solve_stack);

        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input_mip(Box::new(ms), -1));

        // HPresolve::presolve prefix (presolve_loop.rs) up to presolveRuleTest.
        if presolve.model.sense_ == ObjSense::kMaximize {
            for i in 0..presolve.model.num_col_ {
                presolve.model.col_cost_[i as usize] = -presolve.model.col_cost_[i as usize];
            }
            presolve.model.offset_ = -presolve.model.offset_;
            presolve.model.sense_ = ObjSense::kMinimize;
        }
        let silent = presolve.silent_log();
        presolve
            .analysis_
            .setup(&presolve.model, &presolve.options, silent);
        presolve
            .mip_mut()
            .mipdata_
            .cliquetable
            .set_presolve_flag(true);

        // presolveRuleTestProbing, then presolveReturn.
        assert!(presolve.run_probing(&mut postsolve_stack) == HPresolveResult::kOk);
        assert!(presolve.scale_mip(&mut postsolve_stack) == HPresolveResult::kOk);

        let mut rows = Vec::new();
        for row in 1..presolve.model.num_row_ {
            if presolve.row_deleted[row as usize] != 0 {
                continue;
            }
            let mut coefs = vec![0.0; presolve.model.num_col_ as usize];
            let mut it = presolve.get_row_vector(row);
            while let Some(nz) = it.next_nz(&presolve) {
                coefs[nz.index() as usize] = nz.value();
            }
            rows.push((
                coefs,
                presolve.model.row_lower_[row as usize],
                presolve.model.row_upper_[row as usize],
            ));
        }
        rows
    }

    #[test]
    fn h16_agg_fix_vlb_floor_lb_tightened() {
        // z >= 5 (the VLB minimum of 15*x1 + 5 raises lb from 0), so the three VLBs merge into
        // z - 5 x0 - 15 x1 - 25 x2 >= 5, the only row left besides row 0.
        let rows = presolve_rows(true, false);
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert_eq!(rows[0].0, vec![1.0, -5.0, -15.0, -25.0]);
        assert_eq!(rows[0].1, 5.0);
        assert_eq!(rows[0].2, kHighsInf);
    }

    #[test]
    fn h16_agg_fix_vub_ceiling_ub_tightened() {
        // w <= -5 (the VUB maximum of -15 x1 - 5 lowers ub from 0), so the three VUBs merge into
        // w + 5 x0 + 15 x1 + 25 x2 <= -5.
        let rows = presolve_rows(true, true);
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert_eq!(rows[0].0, vec![1.0, 5.0, 15.0, 25.0]);
        assert_eq!(rows[0].1, -kHighsInf);
        assert_eq!(rows[0].2, -5.0);
    }

    #[test]
    fn h16_agg_fix_vlb_only_column_aggregates() {
        // z has VLBs and a finite upper bound but no VUBs: only the VLB side is built and merged.
        // The VUB side has no origin rows and no clique work.
        let rows = presolve_rows(true, false);
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert_eq!(rows[0].2, kHighsInf);
    }

    #[test]
    fn h16_agg_fix_var_bound_value_inf_safe() {
        // 1.16 VarBound::minValue/maxValue return constant + m without CDouble for infinite parts.
        let vb = HighsImplications_VarBound {
            coef: 15.0,
            constant: 5.0,
            origin: 2,
        };
        assert_eq!(var_bound_min_value_116(&vb), 5.0);
        assert_eq!(var_bound_max_value_116(&vb), 20.0);
        let vb_inf = HighsImplications_VarBound {
            coef: -2.0,
            constant: -kHighsInf,
            origin: -1,
        };
        assert_eq!(var_bound_min_value_116(&vb_inf), -kHighsInf);
        assert_eq!(var_bound_max_value_116(&vb_inf), -kHighsInf);
    }

    #[test]
    fn h16_agg_fix_off_no_aggregation() {
        // Option off: the aggregation is not run and the three variable-bound rows stay.
        for vubs in [false, true] {
            let rows = presolve_rows(false, vubs);
            assert_eq!(rows.len(), 3, "{rows:?}");
        }
    }
}
