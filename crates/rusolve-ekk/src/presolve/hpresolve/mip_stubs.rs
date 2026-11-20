//! HiGHS source: presolve/HPresolve.cpp (okSetInput (MIP) 151-172, scaleMIP 2441-2486, applyConflictGraphSubstitutions 2488-2541, strengthenInequalities 7023-7259, detectImpliedIntegers 7261-7265)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! MIP-only methods. Read `presolve/hpresolve/mod.rs` first
//! (ownership, "MIP mode", containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::okSetInput (lines 151-172, 22 lines)
    ///
    /// C++ overload `okSetInput(HighsMipSolver&, HighsInt)`.
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::setInput (lines 111-114)
    ///
    /// PORT NOTE: takes the solver by value (`mod.rs` "MIP mode"); the caller
    /// gets it back with [`HPresolve::take_mip_solver`] whatever this returns. The C++
    /// passes `presolvedModel` by reference; here it is moved into `self.model`.
    pub fn ok_set_input_mip(
        &mut self,
        mut mipsolver: Box<HighsMipSolver>,
        presolve_reduction_limit: HighsInt,
    ) -> bool {
        // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 114 at 7b890c255
        if mipsolver.options_mip_.mip_presolve_clique_table {
            self.presolve_clique_table
                .rebuild(&mipsolver.mipdata_.cliquetable);
        }
        self.probing_contingent = 1000;
        self.probing_num_del_col = 0;
        self.num_probed = 0;
        self.num_probes = vec![0; mipsolver.num_col() as usize];

        if mipsolver.model_ != HighsMipModel::Presolved {
            mipsolver.mipdata_.presolved_model = crate::mip_model!(mipsolver).clone();
            mipsolver.model_ = HighsMipModel::Presolved;
        } else {
            let col_lower = mipsolver.mipdata_.domains[0].col_lower_.clone();
            let col_upper = mipsolver.mipdata_.domains[0].col_upper_.clone();
            mipsolver.mipdata_.presolved_model.col_lower_ = col_lower;
            mipsolver.mipdata_.presolved_model.col_upper_ = col_upper;
        }

        // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 130-131 (setInput) at f2ae2b320
        {
            let flag = mipsolver.options_mip_.mip_aggregate_implications;
            if flag {
                let num_row = mipsolver.num_row();
                mipsolver.mipdata_.implications.resize_rows(num_row);
            }
        }

        let model = std::mem::take(&mut mipsolver.mipdata_.presolved_model);
        let options = mipsolver.options_mip_.clone();
        let timer = Some(mipsolver.timer_.clone());
        // `this->mipsolver = &mipsolver` comes first in C++: okSetInput tests it.
        self.mipsolver = Some(mipsolver);
        self.ok_set_input(model, &options, presolve_reduction_limit, timer)
    }

    /// Rust-only: move the solver back out after `run` (or after a failed
    /// [`HPresolve::ok_set_input_mip`]), with the presolved model put back into
    /// `mipdata_.presolved_model`, where the C++ modified it in place.
    pub fn take_mip_solver(&mut self) -> Box<HighsMipSolver> {
        let mut ms = self
            .mipsolver
            .take()
            .expect("take_mip_solver: not in MIP mode");
        ms.mipdata_.presolved_model = self.take_model();
        ms
    }

    /// C++ `*mipsolver` (read). The presolved model is in `self.model`, not here.
    pub fn mip(&self) -> &HighsMipSolver {
        self.mipsolver.as_deref().expect("MIP mode")
    }

    /// C++ `*mipsolver` (write). The presolved model is in `self.model`, not here.
    pub fn mip_mut(&mut self) -> &mut HighsMipSolver {
        self.mipsolver.as_deref_mut().expect("MIP mode")
    }

    /// HiGHS: HPresolve.cpp::HPresolve::scaleMIP (lines 2441-2486, 46 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn scale_mip(&mut self, postsolve_stack: &mut HighsPostsolveStack) -> HPresolveResult {
        for i in 0..self.model.num_row_ {
            let r = i as usize;
            if self.row_deleted[r] != 0
                || self.rowsize[r] < 1
                || self.rowsize_integer[r] + self.rowsize_impl_int[r] == self.rowsize[r]
            {
                continue;
            }

            self.store_row(i);

            let mut max_abs_val: f64 = 0.0;

            for j in 0..self.rowpositions.len() {
                let nz_pos = self.rowpositions[j] as usize;
                if self.model.integrality_[self.a_col[nz_pos] as usize] != HighsVarType::kContinuous
                {
                    continue;
                }

                max_abs_val = cpp_max(self.a_value[nz_pos].abs(), max_abs_val);
            }

            // C++ `assert(maxAbsVal != 0.0)` is compiled out in the NDEBUG oracle build.

            let mut scale = (-(max_abs_val.log2())).round().exp2();
            if scale == 1.0 {
                continue;
            }

            if self.model.row_upper_[r] == kHighsInf {
                scale = -scale;
            }

            self.scale_stored_row(i, scale, false);
        }

        for i in 0..self.model.num_col_ {
            let col = i as usize;
            if self.col_deleted[col] != 0
                || self.colsize[col] < 1
                || self.model.integrality_[col] != HighsVarType::kContinuous
            {
                continue;
            }

            let mut max_abs_val: f64 = 0.0;

            for nonz in self.get_column_vector(i).iter(self) {
                max_abs_val = cpp_max(nonz.value().abs(), max_abs_val);
            }

            let scale = (-(max_abs_val.log2())).round().exp2();
            if scale == 1.0 {
                continue;
            }

            hpresolve_checked_call!(self.transform_column(postsolve_stack, i, scale, 0.0));
        }

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::applyConflictGraphSubstitutions (lines 2488-2541, 54 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::applyConflictGraphSubstitutions (lines 3034-3064)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn apply_conflict_graph_substitutions(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        num_del_col: &mut HighsInt,
    ) -> HPresolveResult {
        // C++ holds references into mipdata_. The loop bodies never touch these
        // vectors, so iterate over copies and clear the originals after each loop,
        // as the C++ does (an early return leaves them uncleared in both versions).
        let implication_subs: Vec<HighsSubstitution> =
            self.mip().mipdata_.implications.substitutions.clone();
        for substitution in implication_subs {
            let substcol = substitution.substcol as usize;
            let staycol = substitution.staycol as usize;
            if self.col_deleted[substcol] != 0 || self.col_deleted[staycol] != 0 {
                continue;
            }

            *num_del_col += 1;

            postsolve_stack.doubleton_equation(
                -1,
                substitution.substcol,
                substitution.staycol,
                1.0,
                -substitution.scale,
                substitution.offset,
                self.model.col_lower_[substcol],
                self.model.col_upper_[substcol],
                0.0,
                false,
                false,
                RowType::kEq,
                HighsEmptySlice,
            );
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 3034-3036 at 7b890c255
            if self.options.mip_presolve_clique_table {
                hpresolve_checked_call!(self.substitute_x1b(
                    substitution.substcol,
                    substitution.staycol,
                    substitution.offset,
                    substitution.scale,
                    -1,
                ));
            } else {
                self.mark_col_deleted(substitution.substcol);
                self.substitute_col(
                    substitution.substcol,
                    substitution.staycol,
                    substitution.offset,
                    substitution.scale,
                );
            }
            hpresolve_checked_call!(self.check_limits(postsolve_stack));
        }

        self.mip_mut().mipdata_.implications.substitutions.clear();

        let clique_subs: Vec<HighsCliqueTable_Substitution> = self
            .mip()
            .mipdata_
            .cliquetable
            .get_substitutions_2()
            .to_vec();
        for subst in clique_subs {
            let substcol = subst.substcol as usize;
            let replace_col = subst.replace.col() as HighsInt;
            if self.col_deleted[substcol] != 0 || self.col_deleted[replace_col as usize] != 0 {
                continue;
            }

            *num_del_col += 1;

            let (scale, offset): (f64, f64) = if subst.replace.val() == 0 {
                (-1.0, 1.0)
            } else {
                (1.0, 0.0)
            };

            postsolve_stack.doubleton_equation(
                -1,
                subst.substcol,
                replace_col,
                1.0,
                -scale,
                offset,
                self.model.col_lower_[substcol],
                self.model.col_upper_[substcol],
                0.0,
                false,
                false,
                RowType::kEq,
                HighsEmptySlice,
            );
            // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp lines 3063-3064 at 7b890c255
            if self.options.mip_presolve_clique_table {
                hpresolve_checked_call!(self.substitute_x1b(
                    subst.substcol,
                    replace_col,
                    offset,
                    scale,
                    -1
                ));
            } else {
                self.mark_col_deleted(subst.substcol);
                self.substitute_col(subst.substcol, replace_col, offset, scale);
            }
            hpresolve_checked_call!(self.check_limits(postsolve_stack));
        }

        self.mip_mut()
            .mipdata_
            .cliquetable
            .get_substitutions()
            .clear();

        HPresolveResult::kOk
    }

    /// HiGHS: HPresolve.cpp::HPresolve::strengthenInequalities (lines 7023-7259, 237 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    pub fn strengthen_inequalities(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        num_strengthened: &mut HighsInt,
    ) -> HPresolveResult {
        let _prof = super::presolve_prof::scope(super::presolve_prof::STRENGTHEN);
        // HiGHS #3468 (mip_strengthen_small_contrib): gate site; the latest body is
        // strengthen_inequalities_x1_small_contrib, the 1.15.1 body below is unchanged.
        if self.options.mip_strengthen_small_contrib {
            return self
                .strengthen_inequalities_x1_small_contrib(postsolve_stack, num_strengthened);
        }

        let mut complementation: Vec<i8> = Vec::new();
        let mut reducedcost: Vec<f64> = Vec::new();
        let mut upper: Vec<f64> = Vec::new();
        let mut indices: Vec<HighsInt> = Vec::new();
        let mut positions: Vec<HighsInt> = Vec::new();
        let mut stack: Vec<HighsInt> = Vec::new();
        let mut coefs: Vec<f64> = Vec::new();
        let mut cover: Vec<HighsInt> = Vec::new();

        *num_strengthened = 0;
        // Check for timeout according to this frequency
        let check_time_frequency: HighsInt = 100;

        for row in 0..self.model.num_row_ {
            let r = row as usize;
            if self.rowsize[r] <= 1 {
                continue;
            }
            if self.is_ranged(row) {
                continue;
            }

            // do not run on very dense rows as this could get expensive
            let rowsize_limit: HighsInt =
                cpp_max(1000, (self.model.num_col_ - self.num_deleted_cols) / 20);
            if self.rowsize[r] > rowsize_limit {
                continue;
            }

            let mut maxviolation: HighsCDouble;
            let mut continuouscontribution = HighsCDouble::from(0.0);
            let scale: f64;

            if self.model.row_lower_[r] != -kHighsInf {
                // ax >= lb, treat as -ax <= -lb --> -ax + lb <= 0
                maxviolation = HighsCDouble::from(self.model.row_lower_[r]);
                scale = -1.0;
            } else {
                // ax <= ub --> ax - ub <= 0
                maxviolation = HighsCDouble::from(-self.model.row_upper_[r]);
                scale = 1.0;
            }

            complementation.clear();
            reducedcost.clear();
            upper.clear();
            indices.clear();
            positions.clear();
            stack.clear();
            stack.push(self.rowroot[r]);

            let mut skiprow = false;

            while let Some(pos) = stack.pop() {
                let pos_u = pos as usize;
                // add non-zeros to stack
                if self.a_r_right[pos_u] != -1 {
                    stack.push(self.a_r_right[pos_u]);
                }
                if self.a_r_left[pos_u] != -1 {
                    stack.push(self.a_r_left[pos_u]);
                }

                // get column index
                let col = self.a_col[pos_u];
                let cu = col as usize;

                // skip row if a column bound is not finite
                skiprow = self.model.col_lower_[cu] == -kHighsInf
                    || self.model.col_upper_[cu] == kHighsInf;
                if skiprow {
                    break;
                }

                // compute maximum violation
                // scale =  1:  ax <=  ub --> violation =  ax - ub > 0
                // scale = -1: -ax <= -lb --> violation = -ax + lb > 0
                // this means that for scale = 1 we sum up an upper bound on constraint
                // activity, and for scale = -1 we sum up a lower bound on constraint
                // activity.
                let comp: i8;
                let mut weight = self.a_value[pos_u] * scale;
                let ub = self.model.col_upper_[cu] - self.model.col_lower_[cu];
                if weight > 0.0 {
                    comp = 1;
                    maxviolation += self.model.col_upper_[cu] * weight;
                } else {
                    comp = -1;
                    maxviolation += self.model.col_lower_[cu] * weight;
                    weight = -weight;
                }

                if ub <= self.primal_feastol || weight <= self.primal_feastol {
                    continue;
                }

                if self.model.integrality_[cu] == HighsVarType::kContinuous {
                    continuouscontribution += weight * ub;
                    continue;
                }

                indices.push(reducedcost.len() as HighsInt);
                positions.push(pos);
                reducedcost.push(weight);
                complementation.push(comp);
                upper.push(ub);
            }

            // Check for timeout according to frequency, unless a particularly
            // dense row has just been analysed
            if (row & check_time_frequency) == 0 || 10 * self.rowsize[r] > rowsize_limit {
                hpresolve_checked_call!(self.check_time_limit());
            }

            if skiprow {
                stack.clear();
                continue;
            }

            // maxviolation <= 0 implies that the constraint is redundant:
            // scale =  1: upper bound on activity <= model->row_upper_[row]
            // scale = -1: model->row_lower_[row]  <= lower bound on activity
            if maxviolation <= self.primal_feastol {
                hpresolve_checked_call!(self.row_presolve(postsolve_stack, row));
                continue;
            }

            let small_val = cpp_max(
                100.0 * self.primal_feastol,
                self.primal_feastol * f64::from(maxviolation),
            );
            loop {
                if (maxviolation - continuouscontribution) <= small_val || indices.is_empty() {
                    break;
                }

                // C++ pair comparison `make_pair(rc[i1], i1) > make_pair(rc[i2], i2)`,
                // written as the libc++ `operator>` it expands to (`p2 < p1`), so NaN
                // reduced costs order exactly as in C++.
                pdqsort_by(&mut indices[..], |&i1, &i2| {
                    let (rc1, rc2) = (reducedcost[i1 as usize], reducedcost[i2 as usize]);
                    rc2 < rc1 || (!(rc1 < rc2) && i2 < i1)
                });

                let mut lambda = maxviolation - continuouscontribution;

                cover.clear();
                for k in (0..indices.len()).rev() {
                    let index = indices[k];
                    let iu = index as usize;
                    let delta = upper[iu] * reducedcost[iu];

                    if upper[iu] <= 1000.0
                        && reducedcost[iu] > small_val
                        && (lambda - delta) <= small_val
                    {
                        cover.push(index);
                    } else {
                        lambda -= delta;
                    }
                }

                if cover.is_empty() || lambda <= small_val {
                    break;
                }

                // std::min_element with the C++ comparator: libc++ keeps the first
                // element, then replaces it whenever comp(*it, *smallest) holds.
                let cover_less = |i1: HighsInt, i2: HighsInt| -> bool {
                    let (r1, r2) = (reducedcost[i1 as usize], reducedcost[i2 as usize]);
                    if r1 <= 1e-3 || r2 <= 1e-3 {
                        return r1 > r2;
                    }
                    r1 < r2
                };
                let mut alpos = cover[0];
                for &cv in cover.iter().skip(1) {
                    if cover_less(cv, alpos) {
                        alpos = cv;
                    }
                }

                let al = reducedcost[alpos as usize];
                coefs.resize(cover.len(), 0.0);
                let coverrhs = cpp_max(f64::from(lambda / al - self.primal_feastol).ceil(), 1.0);
                let mut slackupper = HighsCDouble::from(-coverrhs);

                let mut step = kHighsInf;
                for i in 0..cover.len() {
                    let ci = cover[i] as usize;
                    coefs[i] = (cpp_min(reducedcost[ci], f64::from(lambda)) / al
                        - self.options.small_matrix_value)
                        .ceil();
                    slackupper += upper[ci] * coefs[i];
                    step = cpp_min(step, reducedcost[ci] / coefs[i]);
                }
                step = cpp_min(step, f64::from(maxviolation / coverrhs));
                maxviolation -= step * coverrhs;

                let slackind = reducedcost.len() as HighsInt;
                reducedcost.push(step);
                upper.push(f64::from(slackupper));

                for i in 0..cover.len() {
                    let ci = cover[i] as usize;
                    reducedcost[ci] -= step * coefs[i];
                }

                let primal_feastol = self.primal_feastol;
                indices.retain(|&i| !(reducedcost[i as usize] <= primal_feastol));
                indices.push(slackind);
            }

            let threshold = f64::from(maxviolation + self.primal_feastol);

            indices.retain(|&i| {
                !((i as usize) >= positions.len() || reducedcost[i as usize].abs() <= threshold)
            });
            if indices.is_empty() {
                continue;
            }

            // update / add non-zeros
            if scale < 0.0 {
                let mut lhs = HighsCDouble::from(self.model.row_lower_[r]);
                self.strengthen_update_non_zeros(
                    row,
                    &mut lhs,
                    -1,
                    &indices,
                    &positions,
                    &reducedcost,
                    &complementation,
                    maxviolation,
                );
                if self.options.h16_row_duals {
                    // HiGHS 1.16 #3270 (h16_row_duals): HPresolveStrengthenIneq.cpp lines 262 (1.16.0)
                    self.h16_change_row_lower(row, f64::from(lhs), false);
                } else {
                    self.model.row_lower_[r] = f64::from(lhs);
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8918-8922 at f2ae2b320 (1.15.1 7248-7252)
                    self.x1_row_modified(row);
                }
            } else {
                let mut rhs = HighsCDouble::from(self.model.row_upper_[r]);
                self.strengthen_update_non_zeros(
                    row,
                    &mut rhs,
                    1,
                    &indices,
                    &positions,
                    &reducedcost,
                    &complementation,
                    maxviolation,
                );
                if self.options.h16_row_duals {
                    // HiGHS 1.16 #3270 (h16_row_duals): HPresolveStrengthenIneq.cpp lines 266 (1.16.0)
                    self.h16_change_row_upper(row, f64::from(rhs), false);
                } else {
                    self.model.row_upper_[r] = f64::from(rhs);
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8918-8922 at f2ae2b320 (1.15.1 7248-7252)
                    self.x1_row_modified(row);
                }
            }

            *num_strengthened += indices.len() as HighsInt;
        }

        HPresolveResult::kOk
    }

    /// HiGHS latest (#3468, 3af3ff719): HPresolveStrengthenIneq.cpp::HPresolve::strengthenInequalities (lines 14-275)
    ///
    /// Called only when `mip_strengthen_small_contrib` is on (see
    /// [`HPresolve::strengthen_inequalities`]); every intermediate is `HighsCDouble`, as
    /// in the merged C++ (`std::vector<HighsCDouble> reducedcost`, `upper`).
    pub fn strengthen_inequalities_x1_small_contrib(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
        num_strengthened: &mut HighsInt,
    ) -> HPresolveResult {
        let mut complementation: Vec<i8> = Vec::new();
        let mut reducedcost: Vec<HighsCDouble> = Vec::new();
        let mut upper: Vec<HighsCDouble> = Vec::new();
        let mut indices: Vec<HighsInt> = Vec::new();
        let mut positions: Vec<HighsInt> = Vec::new();
        let mut stack: Vec<HighsInt> = Vec::new();
        let mut coefs: Vec<f64> = Vec::new();
        let mut cover: Vec<HighsInt> = Vec::new();

        *num_strengthened = 0;
        // Check for timeout according to this frequency
        let check_time_frequency: HighsInt = 100;

        for row in 0..self.model.num_row_ {
            let r = row as usize;
            if self.rowsize[r] <= 1 {
                continue;
            }
            if self.is_ranged(row) {
                continue;
            }

            // do not run on very dense rows as this could get expensive
            let rowsize_limit: HighsInt =
                cpp_max(1000, (self.model.num_col_ - self.num_deleted_cols) / 20);
            if self.rowsize[r] > rowsize_limit {
                continue;
            }

            let mut maxviolation: HighsCDouble;
            let mut continuouscontribution = HighsCDouble::from(0.0);
            // maximum contribution of negligible terms not used in covers; a term is
            // negligible if its maximum contribution is at most primal_feastol
            let mut smallcontribution = HighsCDouble::from(0.0);
            let scale: f64;

            if self.model.row_lower_[r] != -kHighsInf {
                // ax >= lb, treat as -ax <= -lb --> -ax + lb <= 0
                maxviolation = HighsCDouble::from(self.model.row_lower_[r]);
                scale = -1.0;
            } else {
                // ax <= ub --> ax - ub <= 0
                maxviolation = HighsCDouble::from(-self.model.row_upper_[r]);
                scale = 1.0;
            }

            complementation.clear();
            reducedcost.clear();
            upper.clear();
            indices.clear();
            positions.clear();
            stack.clear();
            stack.push(self.rowroot[r]);

            let mut skiprow = false;

            while let Some(pos) = stack.pop() {
                let pos_u = pos as usize;
                // add non-zeros to stack
                if self.a_r_right[pos_u] != -1 {
                    stack.push(self.a_r_right[pos_u]);
                }
                if self.a_r_left[pos_u] != -1 {
                    stack.push(self.a_r_left[pos_u]);
                }

                // get column index
                let col = self.a_col[pos_u];
                let cu = col as usize;

                // skip row if a column bound is not finite
                skiprow = self.model.col_lower_[cu] == -kHighsInf
                    || self.model.col_upper_[cu] == kHighsInf;
                if skiprow {
                    break;
                }

                // compute maximum violation
                // scale =  1:  ax <=  ub --> violation =  ax - ub > 0
                // scale = -1: -ax <= -lb --> violation = -ax + lb > 0
                // this means that for scale = 1 we sum up an upper bound on constraint
                // activity, and for scale = -1 we sum up a lower bound on constraint
                // activity.
                let comp: i8;
                let mut weight = self.a_value[pos_u] * scale;
                let ub = HighsCDouble::from(self.model.col_upper_[cu])
                    - HighsCDouble::from(self.model.col_lower_[cu]);
                if weight > 0.0 {
                    comp = 1;
                    maxviolation += HighsCDouble::from(self.model.col_upper_[cu]) * weight;
                } else {
                    comp = -1;
                    maxviolation += HighsCDouble::from(self.model.col_lower_[cu]) * weight;
                    weight = -weight;
                }

                if self.model.integrality_[cu] == HighsVarType::kContinuous {
                    continuouscontribution += weight * ub;
                    continue;
                }

                if weight * ub <= self.primal_feastol {
                    smallcontribution += weight * ub;
                    continue;
                }

                indices.push(reducedcost.len() as HighsInt);
                positions.push(pos);
                reducedcost.push(HighsCDouble::from(weight));
                complementation.push(comp);
                upper.push(ub);
            }

            // Check for timeout according to frequency, unless a particularly
            // dense row has just been analysed
            if (row & check_time_frequency) == 0 || 10 * self.rowsize[r] > rowsize_limit {
                hpresolve_checked_call!(self.check_time_limit());
            }

            if skiprow {
                stack.clear();
                continue;
            }

            // maxviolation <= 0 implies that the constraint is redundant:
            // scale =  1: upper bound on activity <= model->row_upper_[row]
            // scale = -1: model->row_lower_[row]  <= lower bound on activity
            if maxviolation <= self.primal_feastol {
                hpresolve_checked_call!(self.row_presolve(postsolve_stack, row));
                continue;
            }

            let small_val = cpp_max(
                100.0 * self.primal_feastol,
                self.primal_feastol * f64::from(maxviolation),
            );
            loop {
                if (maxviolation - continuouscontribution - smallcontribution) <= small_val
                    || indices.is_empty()
                {
                    break;
                }

                // C++ pair comparison `make_pair(rc[i1], i1) > make_pair(rc[i2], i2)`,
                // written as the libc++ `operator>` it expands to (`p2 < p1`), so NaN
                // reduced costs order exactly as in C++.
                pdqsort_by(&mut indices[..], |&i1, &i2| {
                    let (rc1, rc2) = (reducedcost[i1 as usize], reducedcost[i2 as usize]);
                    rc2 < rc1 || (!(rc1 < rc2) && i2 < i1)
                });

                let mut lambda = maxviolation - continuouscontribution - smallcontribution;

                cover.clear();
                for k in (0..indices.len()).rev() {
                    let index = indices[k];
                    let iu = index as usize;
                    let delta = upper[iu] * reducedcost[iu];

                    if upper[iu] <= 1000.0
                        && reducedcost[iu] > small_val
                        && (lambda - delta) <= small_val
                    {
                        cover.push(index);
                    } else {
                        lambda -= delta;
                    }
                }

                if cover.is_empty() || lambda <= small_val {
                    break;
                }

                // std::min_element with the C++ comparator: libc++ keeps the first
                // element, then replaces it whenever comp(*it, *smallest) holds.
                let cover_less = |i1: HighsInt, i2: HighsInt| -> bool {
                    let (r1, r2) = (reducedcost[i1 as usize], reducedcost[i2 as usize]);
                    if r1 <= 1e-3 || r2 <= 1e-3 {
                        return r1 > r2;
                    }
                    r1 < r2
                };
                let mut alpos = cover[0];
                for &cv in cover.iter().skip(1) {
                    if cover_less(cv, alpos) {
                        alpos = cv;
                    }
                }

                let al = reducedcost[alpos as usize];
                coefs.resize(cover.len(), 0.0);
                // `max` / `min` on HighsCDouble are the HighsCDouble.h friends (not
                // std::max / std::min): they keep the second argument on a tie.
                let coverrhs = cdouble_max(
                    (lambda / al - self.primal_feastol).ceil(),
                    HighsCDouble::from(1.0),
                );
                let mut slackupper = -coverrhs;

                let mut step = HighsCDouble::from(kHighsInf);
                for i in 0..cover.len() {
                    let ci = cover[i] as usize;
                    coefs[i] = f64::from(
                        (cdouble_min(reducedcost[ci], lambda) / al
                            - self.options.small_matrix_value)
                            .ceil(),
                    );
                    slackupper += upper[ci] * coefs[i];
                    step = cdouble_min(step, reducedcost[ci] / coefs[i]);
                }
                debug_assert!(f64::from(slackupper) >= 0.0);
                // rusolve guard (mip_strengthen_cover_guard): the
                // C++ asserts slackupper >= 0 (compiled out); past it the loop has a fixed point.
                if self.options.mip_strengthen_cover_guard && f64::from(slackupper) < 0.0 {
                    break;
                }
                step = cdouble_min(step, maxviolation / coverrhs);
                maxviolation -= step * coverrhs;

                // only the cover members' reduced costs change, so only they can become
                // negligible; move their maximum contribution to smallcontribution and
                // set their reduced cost to zero, so that they are removed below
                for i in 0..cover.len() {
                    let ci = cover[i] as usize;
                    reducedcost[ci] -= step * coefs[i];
                    if reducedcost[ci] * upper[ci] <= self.primal_feastol {
                        if reducedcost[ci] > 0.0 {
                            smallcontribution += reducedcost[ci] * upper[ci];
                        }
                        reducedcost[ci] = HighsCDouble::from(0.0);
                    }
                }

                // terms with non-positive reduced cost contribute nothing
                indices.retain(|&i| !(reducedcost[i as usize] <= 0.0));
                if step * slackupper <= self.primal_feastol {
                    smallcontribution += step * slackupper;
                } else {
                    indices.push(reducedcost.len() as HighsInt);
                    reducedcost.push(step);
                    upper.push(slackupper);
                }
            }

            let threshold = f64::from(maxviolation + self.primal_feastol);

            indices.retain(|&i| {
                !((i as usize) >= positions.len() || reducedcost[i as usize].abs() <= threshold)
            });
            if indices.is_empty() {
                continue;
            }

            // update / add non-zeros
            // The latest C++ calls changeRowLower / changeRowUpper here; their row-dual
            // update is #3270, ported behind h16_row_duals.
            if scale < 0.0 {
                let mut lhs = HighsCDouble::from(self.model.row_lower_[r]);
                self.strengthen_update_non_zeros_x1_small_contrib(
                    row,
                    &mut lhs,
                    -1,
                    &indices,
                    &positions,
                    &reducedcost,
                    &complementation,
                    maxviolation,
                );
                if self.options.h16_row_duals {
                    // HiGHS 1.16 #3270 (h16_row_duals): HPresolveStrengthenIneq.cpp lines 262 (1.16.0)
                    self.h16_change_row_lower(row, f64::from(lhs), false);
                } else {
                    self.model.row_lower_[r] = f64::from(lhs);
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8918-8922 at f2ae2b320 (1.15.1 7248-7252)
                    self.x1_row_modified(row);
                }
            } else {
                let mut rhs = HighsCDouble::from(self.model.row_upper_[r]);
                self.strengthen_update_non_zeros_x1_small_contrib(
                    row,
                    &mut rhs,
                    1,
                    &indices,
                    &positions,
                    &reducedcost,
                    &complementation,
                    maxviolation,
                );
                if self.options.h16_row_duals {
                    // HiGHS 1.16 #3270 (h16_row_duals): HPresolveStrengthenIneq.cpp lines 266 (1.16.0)
                    self.h16_change_row_upper(row, f64::from(rhs), false);
                } else {
                    self.model.row_upper_[r] = f64::from(rhs);
                    // HiGHS #3264 (mip_aggregate_implications): HPresolve.cpp lines 8918-8922 at f2ae2b320 (1.15.1 7248-7252)
                    self.x1_row_modified(row);
                }
            }

            *num_strengthened += indices.len() as HighsInt;
        }

        HPresolveResult::kOk
    }

    /// HiGHS latest (#3468, 3af3ff719): HPresolveStrengthenIneq.cpp, the `updateNonZeros`
    /// lambda of `HPresolve::strengthenInequalities` (lines 232-248), with `HighsCDouble`
    /// `coefdelta`. Module-private, as [`HPresolve::strengthen_update_non_zeros`].
    #[allow(clippy::too_many_arguments)]
    fn strengthen_update_non_zeros_x1_small_contrib(
        &mut self,
        row: HighsInt,
        rhs: &mut HighsCDouble,
        direction: HighsInt,
        indices: &[HighsInt],
        positions: &[HighsInt],
        reducedcost: &[HighsCDouble],
        complementation: &[i8],
        maxviolation: HighsCDouble,
    ) {
        for &i in indices {
            let iu = i as usize;
            // C++ `assert(Arow[positions[i]] == row)` is compiled out in the NDEBUG oracle.
            // C++ `HighsInt * HighsCDouble`: the friend operator*(double, HighsCDouble).
            let coefdelta: HighsCDouble = direction as f64 * (reducedcost[iu] - maxviolation);
            let col = self.a_col[positions[iu] as usize];
            let cu = col as usize;

            if complementation[iu] == -1 {
                *rhs += coefdelta * self.model.col_lower_[cu];
                self.add_to_matrix(row, col, f64::from(coefdelta));
            } else {
                *rhs -= coefdelta * self.model.col_upper_[cu];
                self.add_to_matrix(row, col, f64::from(-coefdelta));
            }
        }
    }

    /// HiGHS: HPresolve.cpp::HPresolve::detectImpliedIntegers (lines 7261-7265, 5 lines)
    ///
    /// MIP-only (see `mod.rs` "MIP mode").
    ///
    /// Called only when `tryProbing` (initialised to `mipsolver != nullptr`).
    pub fn detect_implied_integers(&mut self) -> HPresolveResult {
        for col in 0..self.model.num_col_ {
            hpresolve_checked_call!(self.convert_implied_integer(col, -1, false).result());
        }
        HPresolveResult::kOk
    }

    /// C++ lambda `updateNonZeros` inside `HPresolve::strengthenInequalities`
    /// (HPresolve.cpp 7196-7214). Captures `this`, `indices`, `positions`,
    /// `reducedcost`, `complementation` and `maxviolation` by reference.
    /// Module-private: it mirrors the lambda, not a separate C++ method.
    #[allow(clippy::too_many_arguments)]
    fn strengthen_update_non_zeros(
        &mut self,
        row: HighsInt,
        rhs: &mut HighsCDouble,
        direction: HighsInt,
        indices: &[HighsInt],
        positions: &[HighsInt],
        reducedcost: &[f64],
        complementation: &[i8],
        maxviolation: HighsCDouble,
    ) {
        for &i in indices {
            let iu = i as usize;
            // C++ `assert(Arow[positions[i]] == row)` is compiled out in the NDEBUG oracle.
            let coefdelta: f64 = direction as f64 * f64::from(reducedcost[iu] - maxviolation);
            let col = self.a_col[positions[iu] as usize];
            let cu = col as usize;

            if complementation[iu] == -1 {
                *rhs += coefdelta * self.model.col_lower_[cu];
                self.add_to_matrix(row, col, coefdelta);
            } else {
                *rhs -= coefdelta * self.model.col_upper_[cu];
                self.add_to_matrix(row, col, -coefdelta);
            }
        }
    }
}

/// HiGHS latest (#3468, 3af3ff719): HighsCDouble.h friend min (lines 322-324)
///
/// `a < b ? a : b`, which keeps `b` on a tie of the double values (std::min keeps `a`).
/// Used only by [`HPresolve::strengthen_inequalities_x1_small_contrib`].
fn cdouble_min(a: HighsCDouble, b: HighsCDouble) -> HighsCDouble {
    if a < b { a } else { b }
}

/// HiGHS latest (#3468, 3af3ff719): HighsCDouble.h friend max (lines 326-328)
///
/// `a > b ? a : b`, which keeps `b` on a tie of the double values (std::max keeps `a`).
/// Used only by [`HPresolve::strengthen_inequalities_x1_small_contrib`].
fn cdouble_max(a: HighsCDouble, b: HighsCDouble) -> HighsCDouble {
    if a > b { a } else { b }
}

#[cfg(test)]
mod x1_stri_tests {
    use super::*;

    /// One row `sum coefs[j] x_j <= rhs` over binary columns, presolved by `ok_set_input`.
    ///
    /// `small_contrib` selects the 1.16.0 body (#3468) or the 1.15.1 body; `guard` is the
    /// rusolve-only `mip_strengthen_cover_guard`.
    fn one_row_presolve(coefs: &[f64], rhs: f64, small_contrib: bool, guard: bool) -> HPresolve {
        let n = coefs.len();
        let mut lp = HighsLp::default();
        lp.num_col_ = n as HighsInt;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0; n];
        lp.col_lower_ = vec![0.0; n];
        lp.col_upper_ = vec![1.0; n];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![rhs];
        lp.integrality_ = vec![HighsVarType::kInteger; n];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = n as HighsInt;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = (0..=n as HighsInt).collect();
        lp.a_matrix_.index_ = vec![0; n];
        lp.a_matrix_.value_ = coefs.to_vec();
        let mut options = HighsOptions::default();
        options.mip_strengthen_small_contrib = small_contrib;
        options.mip_strengthen_cover_guard = guard;
        let mut presolve = HPresolve::default();
        assert!(presolve.ok_set_input(lp, &options, 0, Some(HighsTimer::new())));
        // LP mode (no mipsolver) makes every column continuous; the MIP path keeps
        // the declared integrality, which strengthenInequalities depends on.
        presolve.model.integrality_ = vec![HighsVarType::kInteger; n];
        presolve
    }

    /// Coefficients of row 0 in the current presolve matrix, as (column, value).
    fn row_zero_coefs(presolve: &HPresolve) -> Vec<(HighsInt, f64)> {
        let mut out = Vec::new();
        let mut stack = vec![presolve.rowroot[0]];
        while let Some(pos) = stack.pop() {
            if pos == -1 {
                continue;
            }
            let pu = pos as usize;
            out.push((presolve.a_col[pu], presolve.a_value[pu]));
            stack.push(presolve.a_r_left[pu]);
            stack.push(presolve.a_r_right[pu]);
        }
        out.sort_by_key(|&(col, _)| col);
        out
    }

    /// Runs `strengthen_inequalities` on one row and returns everything the guard
    /// on/off comparison compares: the return value, `num_strengthened`, the row
    /// coefficients, and the row upper bound.
    #[allow(clippy::type_complexity)]
    fn run_one_row(
        coefs: &[f64],
        rhs: f64,
        small_contrib: bool,
        guard: bool,
    ) -> (HPresolveResult, HighsInt, Vec<(HighsInt, f64)>, f64) {
        let mut presolve = one_row_presolve(coefs, rhs, small_contrib, guard);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, coefs.len() as HighsInt);
        let mut num_strengthened: HighsInt = 0;
        let result = presolve.strengthen_inequalities(&mut postsolve, &mut num_strengthened);
        (
            result,
            num_strengthened,
            row_zero_coefs(&presolve),
            presolve.model.row_upper_[0],
        )
    }

    /// Feasibility of `sum coefs[j] x_j <= rhs` at binary point `x` (bit j of `mask`).
    /// Returns the slack; points with |slack| below 1e-6 are not compared.
    fn slack(coefs: &[(HighsInt, f64)], rhs: f64, mask: u32) -> f64 {
        let act: f64 = coefs
            .iter()
            .map(|&(col, v)| if (mask >> col) & 1 == 1 { v } else { 0.0 })
            .sum();
        rhs - act
    }

    /// Every binary point keeps its feasibility status after strengthening.
    fn assert_binary_feasibility_preserved(coefs: &[f64], rhs: f64, small_contrib: bool) {
        let original: Vec<(HighsInt, f64)> = coefs
            .iter()
            .enumerate()
            .map(|(j, &v)| (j as HighsInt, v))
            .collect();
        let mut presolve = one_row_presolve(coefs, rhs, small_contrib, false);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, coefs.len() as HighsInt);
        let mut num_strengthened: HighsInt = 0;
        presolve.strengthen_inequalities(&mut postsolve, &mut num_strengthened);
        let new_coefs = row_zero_coefs(&presolve);
        let new_rhs = presolve.model.row_upper_[0];
        for mask in 0..(1u32 << coefs.len()) {
            let s_old = slack(&original, rhs, mask);
            let s_new = slack(&new_coefs, new_rhs, mask);
            if s_old.abs() < 1e-6 || s_new.abs() < 1e-6 {
                continue;
            }
            assert_eq!(
                s_old >= 0.0,
                s_new >= 0.0,
                "binary point {mask:b}: original slack {s_old}, strengthened slack {s_new}"
            );
        }
    }

    /// Row with the issue 3461 coefficient 1.0000010735453797 and a negligible term
    /// (5e-8, `weight * ub <= primal_feastol`). The latest path strengthens the 2 and
    /// keeps the 1.0000010735453797 term; the 1.15.1 path leaves the row unchanged.
    #[test]
    fn x1_stri_small_contribution_kept() {
        let coefs = [2.0, 1.0000010735453797, 3.0, 4.0, 5e-8];
        let rhs = 5.0;
        assert_binary_feasibility_preserved(&coefs, rhs, true);
        let mut presolve = one_row_presolve(&coefs, rhs, true, false);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, coefs.len() as HighsInt);
        let mut num_strengthened: HighsInt = 0;
        presolve.strengthen_inequalities(&mut postsolve, &mut num_strengthened);
        assert_eq!(num_strengthened, 1);
        let new_coefs = row_zero_coefs(&presolve);
        assert_eq!(new_coefs[1], (1, 1.0000010735453797));
        // a row with coefficients 1 and 1.0000010735453797 (no negligible term)
        assert_binary_feasibility_preserved(&[1.0, 1.0000010735453797, 3.0, 4.0], 4.5, true);
    }

    /// Option off: the same row is left exactly as the 1.15.1 path leaves it.
    #[test]
    fn x1_stri_small_contribution_off_unchanged() {
        let coefs = [2.0, 1.0000010735453797, 3.0, 4.0, 5e-8];
        let mut presolve = one_row_presolve(&coefs, 5.0, false, false);
        let mut postsolve = HighsPostsolveStack::default();
        postsolve.initialize_index_maps(1, coefs.len() as HighsInt);
        let mut num_strengthened: HighsInt = 0;
        presolve.strengthen_inequalities(&mut postsolve, &mut num_strengthened);
        assert_eq!(num_strengthened, 0);
        assert_eq!(presolve.model.row_upper_[0], 5.0);
        let expected: Vec<(HighsInt, f64)> = coefs
            .iter()
            .enumerate()
            .map(|(j, &v)| (j as HighsInt, v))
            .collect();
        assert_eq!(row_zero_coefs(&presolve), expected);
    }

    /// rusolve guard: `slackupper >= 0` holds for every
    /// cover these rows produce, so the guard never fires and on is identical to off — the
    /// same return value, `num_strengthened`, row coefficients and row bound. With the option
    /// off the code path is unchanged, which is the property that matters for parity.
    #[test]
    fn x1_stri_cover_guard_on_matches_off() {
        let rows: [(&[f64], f64); 2] = [
            (&[2.0, 1.0000010735453797, 3.0, 4.0, 5e-8][..], 5.0),
            (&[1.0, 1.0000010735453797, 3.0, 4.0][..], 4.5),
        ];
        for (coefs, rhs) in rows {
            for small_contrib in [true, false] {
                let off = run_one_row(coefs, rhs, small_contrib, false);
                let on = run_one_row(coefs, rhs, small_contrib, true);
                assert_eq!(
                    on, off,
                    "guard changed the outcome for coefs {coefs:?} rhs {rhs} \
                     small_contrib {small_contrib}"
                );
            }
        }
    }
}
