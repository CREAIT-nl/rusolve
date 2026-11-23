//! HiGHS source: presolve/HPresolve.cpp (lines 827-1043)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! shrinkProblem. Read `presolve/hpresolve/mod.rs` first (ownership, MIP
//! mode, containers, slices, bound arrays, lambdas).

use super::*;
use crate::prelude::*;

impl HPresolve {
    /// HiGHS: HPresolve.cpp::HPresolve::shrinkProblem (lines 827-1043, 217 lines)
    /// HiGHS latest (#3312, 7b890c255): HPresolve.cpp::HPresolve::shrinkProblem (line 1258)
    ///
    /// `std::set<HighsInt> newSet` -> `BTreeSet<HighsInt>`.
    /// `eqiters`/`equations` rebuild (C++ 1008-1011): see the module doc.
    /// `impliedRowBounds.shrink(...)` takes no bound arrays.
    pub fn shrink_problem(&mut self, postsolve_stack: &mut HighsPostsolveStack) {
        let old_num_col = self.model.num_col_;
        self.model.num_col_ = 0;
        let mut new_col_index: Vec<HighsInt> = vec![0; old_num_col as usize];
        let have_col_names = self.model.col_names_.len() > 0;
        debug_assert!(!have_col_names || self.model.col_names_.len() == old_num_col as usize);
        for i in 0..old_num_col {
            let iu = i as usize;
            if self.col_deleted[iu] != 0 {
                new_col_index[iu] = -1;
            } else {
                new_col_index[iu] = self.model.num_col_;
                self.model.num_col_ += 1;
                let ni = new_col_index[iu];
                if ni < i {
                    let nu = ni as usize;
                    self.model.col_cost_[nu] = self.model.col_cost_[iu];
                    self.model.col_lower_[nu] = self.model.col_lower_[iu];
                    self.model.col_upper_[nu] = self.model.col_upper_[iu];
                    debug_assert!(!self.model.col_lower_[nu].is_nan());
                    debug_assert!(!self.model.col_upper_[nu].is_nan());
                    self.model.integrality_[nu] = self.model.integrality_[iu];
                    self.impl_col_lower[nu] = self.impl_col_lower[iu];
                    self.impl_col_upper[nu] = self.impl_col_upper[iu];
                    self.col_lower_source[nu] = self.col_lower_source[iu];
                    self.col_upper_source[nu] = self.col_upper_source[iu];
                    self.impl_row_dual_source_by_col[nu] =
                        self.impl_row_dual_source_by_col[iu].clone();
                    self.colhead[nu] = self.colhead[iu];
                    self.colsize[nu] = self.colsize[iu];
                    if have_col_names {
                        self.model.col_names_[nu] = std::mem::take(&mut self.model.col_names_[iu]);
                    }
                    self.changed_col_flag[nu] = self.changed_col_flag[iu];
                }
            }
        }
        // HiGHS 1.16 (#3255): HPresolve.cpp lines 1097-1098; -1 unless FM reformulated.
        if self.model.fme_obj_col_ >= 0 {
            self.model.fme_obj_col_ = new_col_index[self.model.fme_obj_col_ as usize];
        }
        self.col_deleted = vec![0; self.model.num_col_ as usize];
        let num_col = self.model.num_col_ as usize;
        self.model.col_cost_.resize(num_col, 0.0);
        self.model.col_lower_.resize(num_col, 0.0);
        self.model.col_upper_.resize(num_col, 0.0);
        self.model
            .integrality_
            .resize(num_col, HighsVarType::default());
        self.impl_col_lower.resize(num_col, 0.0);
        self.impl_col_upper.resize(num_col, 0.0);
        self.col_lower_source.resize(num_col, 0);
        self.col_upper_source.resize(num_col, 0);
        self.impl_row_dual_source_by_col
            .resize(num_col, BTreeSet::new());
        self.colhead.resize(num_col, 0);
        self.colsize.resize(num_col, 0);
        if have_col_names {
            self.model.col_names_.resize(num_col, String::new());
        }
        self.changed_col_flag.resize(num_col, 0);
        self.num_deleted_cols = 0;

        let old_num_row = self.model.num_row_;
        let have_row_names = self.model.row_names_.len() > 0;
        debug_assert!(!have_row_names || self.model.row_names_.len() == old_num_row as usize);
        self.model.num_row_ = 0;
        let mut new_row_index: Vec<HighsInt> = vec![0; old_num_row as usize];
        for i in 0..old_num_row {
            let iu = i as usize;
            if self.row_deleted[iu] != 0 {
                new_row_index[iu] = -1;
            } else {
                new_row_index[iu] = self.model.num_row_;
                self.model.num_row_ += 1;
                let ni = new_row_index[iu];
                if ni < i {
                    let nu = ni as usize;
                    self.model.row_lower_[nu] = self.model.row_lower_[iu];
                    self.model.row_upper_[nu] = self.model.row_upper_[iu];
                    debug_assert!(!self.model.row_lower_[nu].is_nan());
                    debug_assert!(!self.model.row_upper_[nu].is_nan());
                    self.row_dual_lower[nu] = self.row_dual_lower[iu];
                    self.row_dual_upper[nu] = self.row_dual_upper[iu];
                    self.impl_row_dual_lower[nu] = self.impl_row_dual_lower[iu];
                    self.impl_row_dual_upper[nu] = self.impl_row_dual_upper[iu];
                    self.row_dual_lower_source[nu] = self.row_dual_lower_source[iu];
                    self.row_dual_upper_source[nu] = self.row_dual_upper_source[iu];
                    self.col_impl_source_by_row[nu] = self.col_impl_source_by_row[iu].clone();
                    self.rowroot[nu] = self.rowroot[iu];
                    self.rowsize[nu] = self.rowsize[iu];
                    self.rowsize_integer[nu] = self.rowsize_integer[iu];
                    self.rowsize_impl_int[nu] = self.rowsize_impl_int[iu];
                    // HiGHS 1.16 #3367: HPresolve.cpp line 1136 (1.16.0)
                    self.rowsize_singleton[nu] = self.rowsize_singleton[iu];
                    if have_row_names {
                        self.model.row_names_[nu] = std::mem::take(&mut self.model.row_names_[iu]);
                    }
                    self.changed_row_flag[nu] = self.changed_row_flag[iu];
                    self.single_equation_checked[nu] = self.single_equation_checked[iu];
                }
            }
        }
        for i in 0..self.model.num_col_ as usize {
            if self.col_lower_source[i] != -1 {
                self.col_lower_source[i] = new_row_index[self.col_lower_source[i] as usize];
            }
            if self.col_upper_source[i] != -1 {
                self.col_upper_source[i] = new_row_index[self.col_upper_source[i] as usize];
            }
        }

        for i in 0..self.model.num_row_ as usize {
            if self.row_dual_lower_source[i] != -1 {
                self.row_dual_lower_source[i] =
                    new_col_index[self.row_dual_lower_source[i] as usize];
            }
            if self.row_dual_upper_source[i] != -1 {
                self.row_dual_upper_source[i] =
                    new_col_index[self.row_dual_upper_source[i] as usize];
            }
        }

        for i in 0..self.model.num_col_ as usize {
            let mut new_set: BTreeSet<HighsInt> = BTreeSet::new();
            for &row in self.impl_row_dual_source_by_col[i].iter() {
                if new_row_index[row as usize] != -1 {
                    new_set.insert(new_row_index[row as usize]);
                }
            }
            self.impl_row_dual_source_by_col[i] = new_set;
        }

        for i in 0..self.model.num_row_ as usize {
            let mut new_set: BTreeSet<HighsInt> = BTreeSet::new();
            for &col in self.col_impl_source_by_row[i].iter() {
                if new_col_index[col as usize] != -1 {
                    new_set.insert(new_col_index[col as usize]);
                }
            }
            self.col_impl_source_by_row[i] = new_set;
        }
        let num_row = self.model.num_row_ as usize;
        self.row_deleted = vec![0; num_row];
        self.model.row_lower_.resize(num_row, 0.0);
        self.model.row_upper_.resize(num_row, 0.0);
        self.row_dual_lower.resize(num_row, 0.0);
        self.row_dual_upper.resize(num_row, 0.0);
        self.impl_row_dual_lower.resize(num_row, 0.0);
        self.impl_row_dual_upper.resize(num_row, 0.0);
        self.row_dual_lower_source.resize(num_row, 0);
        self.row_dual_upper_source.resize(num_row, 0);
        self.col_impl_source_by_row.resize(num_row, BTreeSet::new());
        self.rowroot.resize(num_row, 0);
        self.rowsize.resize(num_row, 0);
        self.rowsize_integer.resize(num_row, 0);
        self.rowsize_impl_int.resize(num_row, 0);
        self.rowsize_singleton.resize(num_row, 0);
        if have_row_names {
            self.model.row_names_.resize(num_row, String::new());
        }
        self.changed_row_flag.resize(num_row, 0);
        self.single_equation_checked.resize(num_row, 0);

        self.num_deleted_rows = 0;
        postsolve_stack.compress_index_maps(&new_row_index, &new_col_index);
        self.implied_row_bounds
            .shrink(&new_row_index, self.model.num_row_);
        self.implied_dual_row_bounds
            .shrink(&new_col_index, self.model.num_col_);

        for i in 0..self.a_value.len() {
            if self.a_value[i] == 0.0 {
                continue;
            }
            debug_assert!(new_col_index[self.a_col[i] as usize] != -1);
            debug_assert!(new_row_index[self.a_row[i] as usize] != -1);
            self.a_col[i] = new_col_index[self.a_col[i] as usize];
            self.a_row[i] = new_row_index[self.a_row[i] as usize];
        }

        // update index sets
        for sing_col in self.singleton_columns.iter_mut() {
            *sing_col = new_col_index[*sing_col as usize];
        }
        self.singleton_columns.retain(|&x| x != -1);

        for chg_col in self.changed_col_indices.iter_mut() {
            *chg_col = new_col_index[*chg_col as usize];
        }
        self.changed_col_indices.retain(|&x| x != -1);

        for sing_row in self.singleton_rows.iter_mut() {
            *sing_row = new_row_index[*sing_row as usize];
        }
        self.singleton_rows.retain(|&x| x != -1);

        for chg_row in self.changed_row_indices.iter_mut() {
            *chg_row = new_row_index[*chg_row as usize];
        }
        self.changed_row_indices.retain(|&x| x != -1);

        for row_col_pair in self.substitution_opportunities.iter_mut() {
            // skip deleted elements
            if row_col_pair.0 == -1 {
                continue;
            }
            row_col_pair.0 = new_row_index[row_col_pair.0 as usize];
            row_col_pair.1 = new_col_index[row_col_pair.1 as usize];
        }
        self.substitution_opportunities
            .retain(|p| !(p.0 == -1 || p.1 == -1));

        // todo remove equation set and replace with a vector of doubleton eqs
        self.equations.clear();
        self.eqiters = vec![None; self.model.num_row_ as usize];
        for i in 0..self.model.num_row_ {
            if self.is_equation(i) {
                let key = (self.rowsize[i as usize], i);
                self.equations.insert(key);
                self.eqiters[i as usize] = Some(key);
            }
        }

        // HiGHS: `if (mipsolver != nullptr) { ... }` (HPresolve.cpp 1014-1034)
        if self.mipsolver.is_some() {
            // C++ `model->num_col_` is the already shrunk column count.
            let ncols = self.model.num_col_;
            with_mip_model!(self, ms => {
                ms.mipdata_.row_matrix_set = false;
                let objective_function = HighsObjectiveFunction::new(&ms);
                *ms.mipdata_.objective_function = objective_function;
                // C++ `getDomain() = HighsDomain(*mipsolver)` is HighsDomain::operator=
                // (HighsDomain.h:388) from a temporary without propagators: the deque
                // assignments destroy (unregister) the old cut/conflict propagators, and
                // `infeasible_pos` etc. keep their values. `HPresolve::run` re-adds the pools.
                let dom_tmp = HighsDomain::new(&ms);
                crate::detach!(ms.mipdata_.domains[0], |d| d.op_assign(&mut ms, &dom_tmp));
                crate::detach!(ms.mipdata_.cliquetable, |ct| {
                    crate::detach!(ms.mipdata_.domains[0], |d| ct.rebuild(
                        &mut ms,
                        ncols,
                        postsolve_stack,
                        d,
                        &new_col_index,
                        &new_row_index,
                    ))
                });
                // HiGHS #3312 (mip_presolve_clique_table): HPresolve.cpp line 1258 at 7b890c255
                if self.options.mip_presolve_clique_table {
                    self.presolve_clique_table
                        .rebuild(&ms.mipdata_.cliquetable);
                }
                // C++ `implications.rebuild` reads `mipdata_->postSolveStack`, which is
                // `postsolve_stack` (see `with_mip_model_pss!`): swap it in for the call.
                std::mem::swap(&mut *postsolve_stack, &mut ms.mipdata_.post_solve_stack);
                crate::detach!(ms.mipdata_.implications, |imp| imp.rebuild(
                    &mut ms,
                    ncols,
                    &new_col_index,
                    &new_row_index,
                ));
                std::mem::swap(&mut *postsolve_stack, &mut ms.mipdata_.post_solve_stack);
                let agelim = ms.options_mip_.mip_pool_age_limit;
                let softlimit = ms.options_mip_.mip_pool_soft_limit;
                // HiGHS 1.16 #3574 (h16_cutpool_age): HPresolve.cpp shrinkProblem (1.16.0), flag from the MIP options
                let h16_cutpool_age = ms.options_mip_.h16_cutpool_age;
                *ms.mipdata_.cutpools[0] = HighsCutPool::new_with_age_flag(
                    ncols,
                    agelim,
                    softlimit,
                    0,
                    h16_cutpool_age,
                );
                *ms.mipdata_.conflictpools[0] =
                    HighsConflictPool::new(5 * agelim, softlimit);

                for i in 0..old_num_col as usize {
                    if new_col_index[i] != -1 {
                        self.num_probes[new_col_index[i] as usize] = self.num_probes[i];
                    }
                }

                ms.mipdata_.debug_solution.shrink(&new_col_index);
                self.num_probes.resize(ncols as usize, 0);
                // Need to set the constraint matrix dimensions
                ms.mipdata_.presolved_model.set_matrix_dimensions();
            });
        }
        // Need to set the constraint matrix dimensions
        self.model.set_matrix_dimensions();
        // Need to reset current number of deleted rows and columns in logging
        self.analysis_.reset_num_deleted();
    }
}
