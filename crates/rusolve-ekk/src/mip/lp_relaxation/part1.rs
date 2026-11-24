//! HiGHS source: mip/HighsLpRelaxation.cpp (lines 1-926)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

impl HighsLpRelaxation {
    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::setProfiling (lines 23-26)
    /// C++: `void setProfiling(HighsProfiling * profiling)`
    // PORT NOTE: the profiler is `ms.profiling_` (a no-op stub); `Highs::set_profiling` takes
    // no argument, so neither does this.
    pub fn set_profiling(&mut self) {
        self.lpsolver.set_profiling();
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::getCutPool (lines 28-87)
    /// C++: `void getCutPool(HighsInt & num_col, HighsInt & num_cut, std::vector<double> & cut_lower, std::vector<double> & cut_upper, HighsSparseMatrix & cut_matrix) const`
    pub fn get_cut_pool(
        &self,
        ms: &HighsMipSolver,
        num_col: &mut HighsInt,
        num_cut: &mut HighsInt,
        cut_lower: &mut Vec<f64>,
        cut_upper: &mut Vec<f64>,
        cut_matrix: &mut HighsSparseMatrix,
    ) {
        let lp = self.lpsolver.get_lp();
        *num_col = lp.num_col_;
        let num_lp_row = lp.num_row_;
        let num_model_row = ms.num_row();
        let n_cut = num_lp_row - num_model_row;
        *num_cut = n_cut;
        cut_lower.resize(n_cut as usize, 0.0);
        cut_upper.resize(n_cut as usize, 0.0);
        // Get a map from row index to cut row index
        let mut cut_row_index: Vec<HighsInt> = vec![-1; num_lp_row as usize];
        let mut cut_num: HighsInt = 0;
        for i_row in 0..lp.num_row_ as usize {
            if self.lprows[i_row].origin != HighsLpRelaxation_LpRow::kCutPool {
                continue;
            }
            cut_row_index[i_row] = cut_num;
            cut_lower[cut_num as usize] = lp.row_lower_[i_row];
            cut_upper[cut_num as usize] = lp.row_upper_[i_row];
            cut_num += 1;
        }
        debug_assert!(cut_num == n_cut);

        cut_matrix.num_col_ = lp.num_col_;
        cut_matrix.num_row_ = n_cut;
        cut_matrix.format_ = MatrixFormat::kRowwise;

        let mut cut_matrix_length: Vec<HighsInt> = vec![0; n_cut as usize];
        for i_col in 0..lp.num_col_ as usize {
            for i_el in lp.a_matrix_.start_[i_col]..lp.a_matrix_.start_[i_col + 1] {
                let i_cut = cut_row_index[lp.a_matrix_.index_[i_el as usize] as usize];
                if i_cut >= 0 {
                    cut_matrix_length[i_cut as usize] += 1;
                }
            }
        }
        cut_matrix.start_.resize(n_cut as usize + 1, 0);
        cut_matrix.start_[0] = 0;
        let mut num_cut_nz: HighsInt = 0;
        for i_cut in 0..n_cut as usize {
            let length = cut_matrix_length[i_cut];
            cut_matrix_length[i_cut] = cut_matrix.start_[i_cut];
            num_cut_nz += length;
            cut_matrix.start_[i_cut + 1] = num_cut_nz;
        }
        cut_matrix.index_.resize(num_cut_nz as usize, 0);
        cut_matrix.value_.resize(num_cut_nz as usize, 0.0);
        for i_col in 0..lp.num_col_ as usize {
            for i_el in lp.a_matrix_.start_[i_col]..lp.a_matrix_.start_[i_col + 1] {
                let i_cut = cut_row_index[lp.a_matrix_.index_[i_el as usize] as usize];
                if i_cut >= 0 {
                    let pos = cut_matrix_length[i_cut as usize] as usize;
                    cut_matrix.index_[pos] = i_col as HighsInt;
                    cut_matrix.value_[pos] = lp.a_matrix_.value_[i_el as usize];
                    cut_matrix_length[i_cut as usize] += 1;
                }
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::slackLower (lines 141-156)
    /// C++: `double slackLower(HighsInt row, const HighsDomain & globaldom) const`
    // PORT NOTE: `globaldom` is an index into `ms.mipdata_.domains` (the worker's global domain
    // or `domains[0]`). `HighsDomain::get_min_cut_activity` replays the pool event log first
    // (`&mut`, domain/stack.rs), so the domain is detached for that call; the replay is what
    // the C++ eager pool -> domain notification had already done, so the value is identical.
    pub fn slack_lower(&self, ms: &mut HighsMipSolver, row: HighsInt, globaldom: usize) -> f64 {
        let lprow = self.lprows[row as usize];
        match lprow.origin {
            HighsLpRelaxation_LpRow::kCutPool => {
                let cutpool = lprow.cutpoolindex as usize;
                // Nothing to replay: read the domain in place (no detach round trip).
                if !ms.mipdata_.domains[globaldom].has_pending_pool_events(ms) {
                    return ms.mipdata_.domains[globaldom]
                        .get_min_cut_activity_replayed(cutpool, lprow.index);
                }
                detach!(ms.mipdata_.domains[globaldom], |d| d.get_min_cut_activity(
                    ms,
                    cutpool,
                    lprow.index
                ))
            }
            _ => {
                let rowlower = self.row_lower(row);
                if rowlower != -kHighsInf {
                    return rowlower;
                }
                ms.mipdata_.domains[globaldom].get_min_activity(lprow.index)
            }
        }
    }

    /// PORT NOTE (perf): [`HighsLpRelaxation::slack_lower`] for a caller that knows no pool event
    /// is pending for `globaldom` (`has_pending_pool_events` is false). Then the detach-and-replay
    /// branch of `slack_lower` is a no-op and this reads the domain in place. Same value bit for
    /// bit. `HighsTransformedLp::transform` uses it after one check per call.
    // PORT NOTE (perf): `#[inline]` only.
    #[inline]
    pub fn slack_lower_no_pending(
        &self,
        ms: &HighsMipSolver,
        row: HighsInt,
        globaldom: usize,
    ) -> f64 {
        let lprow = self.lprows[row as usize];
        match lprow.origin {
            HighsLpRelaxation_LpRow::kCutPool => ms.mipdata_.domains[globaldom]
                .get_min_cut_activity_replayed(lprow.cutpoolindex as usize, lprow.index),
            _ => {
                let rowlower = self.row_lower(row);
                if rowlower != -kHighsInf {
                    return rowlower;
                }
                ms.mipdata_.domains[globaldom].get_min_activity(lprow.index)
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::slackUpper (lines 158-171)
    /// C++: `double slackUpper(HighsInt row, const HighsDomain & globaldom) const`
    // PORT NOTE: `globaldom` is an index into `ms.mipdata_.domains` (see `slack_lower`).
    // PORT NOTE (perf): `#[inline]` only.
    #[inline]
    pub fn slack_upper(&self, ms: &HighsMipSolver, row: HighsInt, globaldom: usize) -> f64 {
        let rowupper = self.row_upper(row);
        let lprow = self.lprows[row as usize];
        match lprow.origin {
            HighsLpRelaxation_LpRow::kCutPool => rowupper,
            _ => {
                if rowupper != kHighsInf {
                    return rowupper;
                }
                ms.mipdata_.domains[globaldom].get_max_activity(lprow.index)
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::HighsLpRelaxation (lines 179-206)
    /// C++: `HighsLpRelaxation(const HighsMipSolver & mipsolver)`
    pub fn new(ms: &HighsMipSolver) -> Self {
        let mut this = HighsLpRelaxation {
            worker_: None,
            ..Default::default()
        };
        this.lpsolver.set_option_value_bool("output_flag", false);
        // MIP relaxations keep the faithful orientation
        this.lpsolver
            .set_option_value_bool("lp_dualize_auto", false);
        // h16_sub_solver_flags: C++ `Highs lpsolver` has default options, which in 1.16.0
        // run 1.16.0 code; the port carries the parent's 1.16 flags
        if ms.options_mip_.h16_sub_solver_flags {
            this.lpsolver
                .options_
                .copy_version_flags_from(&ms.options_mip_);
        }
        this.lpsolver
            .set_option_value_int("random_seed", ms.options_mip_.random_seed);
        // Set primal feasibility tolerance for LP solves according to
        // mip_feasibility_tolerance, and smaller tolerance for dual
        // feasibility
        let mip_primal_feasibility_tolerance = ms.options_mip_.mip_feasibility_tolerance;
        let mip_dual_feasibility_tolerance = ms.options_mip_.mip_feasibility_tolerance * 0.1;
        this.lpsolver.set_option_value_double(
            "primal_feasibility_tolerance",
            mip_primal_feasibility_tolerance,
        );
        this.lpsolver
            .set_option_value_double("dual_feasibility_tolerance", mip_dual_feasibility_tolerance);
        // rusolve lever (default off, not in HiGHS): most LP solves of a MIP search are short
        // node and sub-MIP re-solves, and each optimality rebuild runs one FTRAN and one BTRAN
        // (HEkk::factorSolveError) to test the updated factor before keeping it. A tolerance
        // <= 0 is the HiGHS option value that skips the test and keeps the factor.
        if ms.options_mip_.mip_lp_skip_refactor_check {
            this.lpsolver
                .set_option_value_double("rebuild_refactor_solution_error_tolerance", 0.0);
        }
        // rusolve lever (default off, not in HiGHS): keep the dual steepest-edge weights by
        // variable across the basis restores and cut row changes of the search, see
        // simplex/ekk/dse_keep.rs.
        if ms.options_mip_.mip_lp_keep_dse_weights > 0 {
            this.lpsolver.set_option_value_int(
                "mip_lp_keep_dse_weights",
                ms.options_mip_.mip_lp_keep_dse_weights,
            );
        }
        // rusolve lever (default off, not in HiGHS): `storeBasis` records the cut rows of each
        // stored basis, see node_basis.rs.
        if ms.options_mip_.mip_node_install_basis {
            this.lpsolver
                .set_option_value_bool("mip_node_install_basis", true);
        }
        this.status = HighsLpRelaxation_Status::kNotSet;
        this.numlpiters = 0;
        this.avg_solve_iters = 0.0;
        this.num_solved = 0;
        this.epochs = 0;
        this.max_num_fractional = 0;
        this.last_age_call = 0;
        this.objective = -kHighsInf;
        this.currentbasisstored = false;
        this.adjust_sym_branching_col = true;
        this.solved_first_lp = true;
        this.row_ep.size = 0;
        this
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::HighsLpRelaxation (lines 211-236)
    /// C++: `HighsLpRelaxation(const HighsLpRelaxation & other)`
    pub fn copy_from(ms: &HighsMipSolver, other: &HighsLpRelaxation) -> Self {
        Self::copy_from_num_col(other, ms.num_col() as usize)
    }

    /// Body of the C++ copy constructor; `num_col` is `mipsolver.numCol()` (see the `Clone`
    /// impl in mod.rs).
    pub(crate) fn copy_from_num_col(other: &HighsLpRelaxation, num_col: usize) -> Self {
        let mut this = HighsLpRelaxation {
            lprows: other.lprows.clone(),
            fractionalints: other.fractionalints.clone(),
            objective: other.objective,
            basischeckpoint: other.basischeckpoint.clone(),
            currentbasisstored: other.currentbasisstored,
            adjust_sym_branching_col: other.adjust_sym_branching_col,
            worker_: None,
            ..Default::default()
        };
        this.lpsolver.set_option_value_bool("output_flag", false);
        // MIP relaxations keep the faithful orientation
        this.lpsolver
            .set_option_value_bool("lp_dualize_auto", false);
        this.lpsolver.pass_options(other.lpsolver.get_options());
        this.lpsolver.pass_model(other.lpsolver.get_lp().clone());
        // rusolve mip_node_bigm: the copy serves another domain, so it starts from the global rows.
        if other.nlocal_rows == 0 {
            Self::undo_mods(&mut this.lpsolver, &other.local_mods);
        } else {
            let n = this.lpsolver.get_num_row();
            this.lpsolver
                .delete_rows_interval(n - other.nlocal_rows as HighsInt, n - 1);
            let keep = this.lprows.len() - other.nlocal_rows;
            this.lprows.truncate(keep);
        }
        // C++ `setBasis(basis)` uses the default origin "".
        this.lpsolver.set_basis(other.lpsolver.get_basis(), "");
        this.col_lb_buffer.resize(num_col, 0.0);
        this.col_ub_buffer.resize(num_col, 0.0);
        this.status = HighsLpRelaxation_Status::kNotSet;
        this.numlpiters = 0;
        this.avg_solve_iters = 0.0;
        this.num_solved = 0;
        this.epochs = 0;
        this.max_num_fractional = 0;
        this.last_age_call = 0;
        this.objective = -kHighsInf;
        this.solved_first_lp = true;
        this.row_ep.size = 0;
        this
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::loadModel (lines 238-256)
    /// C++: `void loadModel()`
    // PORT NOTE: `worker_->getGlobalDomain()` is `domains[workers[w].globaldom_]`.
    pub fn load_model(&mut self, ms: &HighsMipSolver) {
        let mut lpmodel: HighsLp = mip_model!(ms).clone();
        let globaldom = match self.worker_ {
            Some(w) => ms.mipdata_.workers[w].globaldom_,
            None => 0,
        };
        lpmodel.col_lower_ = ms.mipdata_.domains[globaldom].col_lower_.clone();
        lpmodel.col_upper_ = ms.mipdata_.domains[globaldom].col_upper_.clone();
        lpmodel.offset_ = 0.0;
        self.local_mods.clear();
        self.nlocal_rows = 0;
        self.lprows.clear();
        self.lprows.reserve(lpmodel.num_row_ as usize);
        for i in 0..lpmodel.num_row_ {
            self.lprows.push(HighsLpRelaxation_LpRow::model(i));
        }
        lpmodel.integrality_.clear();
        let num_col = lpmodel.num_col_ as usize;
        self.lpsolver.clear_solver();
        self.lpsolver.clear_model();
        self.lpsolver.pass_model(lpmodel);
        self.col_lb_buffer.resize(num_col, 0.0);
        self.col_ub_buffer.resize(num_col, 0.0);
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::resetToGlobalDomain (lines 258-262)
    /// C++: `void resetToGlobalDomain(const HighsDomain & globaldom)`
    pub fn reset_to_global_domain(&mut self, ms: &HighsMipSolver, globaldom: &HighsDomain) {
        self.lpsolver.change_cols_bounds_interval(
            0,
            ms.num_col() - 1,
            &globaldom.col_lower_,
            &globaldom.col_upper_,
        );
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::computeBasicDegenerateDuals (lines 264-426)
    /// C++: `void computeBasicDegenerateDuals(double threshold, HighsDomain & localdom, HighsDomain & globaldom, HighsConflictPool & conflictpool, HighsPseudocost & pseudocost, bool getdualproof)`
    pub fn compute_basic_degenerate_duals(
        &mut self,
        ms: &mut HighsMipSolver,
        threshold: f64,
        localdom: &mut HighsDomain,
        globaldom: usize,
        conflictpool: usize,
        pseudocost: &mut HighsPseudocost,
        getdualproof: bool,
    ) {
        // PORT NOTE: `globaldom` / `conflictpool` are indices into `ms.mipdata_.domains` /
        // `conflictpools` (the HighsDomain::conflictAnalyzeReconvergence convention); the
        // `localdom` is detached by the caller.
        if !self.lpsolver.has_invert() || self.has_local_mods() {
            return;
        }

        let feastol = ms.mipdata_.feastol;
        let eps = ms.mipdata_.epsilon;

        let mut k: HighsInt = 0;
        for idx in 0..ms.mipdata_.integral_cols.len() {
            let col = ms.mipdata_.integral_cols[idx];
            let ci = col as usize;
            if self.lpsolver.get_basis().col_status[ci] != HighsBasisStatus::kBasic {
                continue;
            }
            let lb = self.lpsolver.get_lp().col_lower_[ci];
            let ub = self.lpsolver.get_lp().col_upper_[ci];
            if ub - lb < feastol {
                continue;
            }

            let value = self.lpsolver.solution_.col_value[ci];
            if value - lb < ub - value {
                if value > lb + feastol {
                    continue;
                }

                self.lpsolver.solution_.col_dual[ci] = 1.0;
                k += 1;
            } else {
                if value < ub - feastol {
                    continue;
                }

                self.lpsolver.solution_.col_dual[ci] = -1.0;
                k += 1;
            }
        }

        if k == 0 {
            return;
        }

        let num_row = self.lpsolver.get_lp().num_row_;
        let num_col = self.lpsolver.get_lp().num_col_;

        if self.row_ep.size < num_row {
            self.row_ep.setup(num_row);

            if (self.row_ap.values.len() as HighsInt) < num_col {
                self.row_ap.set_dimension(num_col);
                self.dualproofvals.reserve(num_col as usize);
                self.dualproofinds.reserve(num_col as usize);
            }
        }

        // HiGHS quirk: the C++ `for (row = 0; k > 0; row++)` runs `row++` on `continue` too, so
        // the increment is taken at the top of the body.
        let mut row: HighsInt = 0;
        while k > 0 {
            let cur_row = row;
            row += 1;
            // `basicIndex` is re-read per iteration: BTRAN needs `&mut self.lpsolver`.
            let var = self.lpsolver.get_basic_variables_array()[cur_row as usize];
            if var >= num_col {
                continue;
            }
            if self.lpsolver.solution_.col_dual[var as usize].abs() != 1.0 {
                continue;
            }

            k -= 1;
            let _ = self
                .lpsolver
                .get_basis_inverse_row_sparse(cur_row, &mut self.row_ep);

            let sign = self.lpsolver.solution_.col_dual[var as usize];
            self.lpsolver.solution_.col_dual[var as usize] = 0.0;
            let mut degenerate_col_dual = kHighsInf;
            for i in 0..self.row_ep.count {
                let i_row = self.row_ep.index[i as usize];
                let lb = self.lpsolver.get_lp().row_lower_[i_row as usize];
                let ub = self.lpsolver.get_lp().row_upper_[i_row as usize];
                if lb == ub {
                    continue;
                }
                let dual = self.lpsolver.solution_.row_dual[i_row as usize];

                let val = -sign * self.row_ep.array[i_row as usize];
                if val > 0.0 {
                    if self.lpsolver.solution_.row_value[i_row as usize] - lb > feastol {
                        degenerate_col_dual = cpp_min(degenerate_col_dual, -dual / val);
                        if degenerate_col_dual < threshold {
                            break;
                        }
                    }
                } else if ub - self.lpsolver.solution_.row_value[i_row as usize] > feastol {
                    degenerate_col_dual = cpp_min(degenerate_col_dual, -dual / val);
                    if degenerate_col_dual < threshold {
                        break;
                    }
                }
            }

            if degenerate_col_dual < threshold {
                continue;
            }

            self.row_ap.clear();
            for i in 0..self.row_ep.count {
                let i_row = self.row_ep.index[i as usize];

                let mut len: HighsInt = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                self.get_row(&*ms, i_row, &mut len, &mut inds, &mut vals);

                for j in 0..len as usize {
                    self.row_ap
                        .add(inds[j], self.row_ep.array[i_row as usize] * vals[j]);
                }
            }

            self.row_ap.cleanup(|_, val| val.abs() <= eps);

            for idx in 0..self.row_ap.nonzeroinds.len() {
                let i_col = self.row_ap.nonzeroinds[idx];
                if i_col == var {
                    continue;
                }

                let lb = self.lpsolver.get_lp().col_lower_[i_col as usize];
                let ub = self.lpsolver.get_lp().col_upper_[i_col as usize];
                if lb == ub {
                    continue;
                }
                let dual = self.lpsolver.solution_.col_dual[i_col as usize];

                let val = sign * self.row_ap.get_value(i_col);
                if val > eps {
                    if self.lpsolver.solution_.col_value[i_col as usize] - lb > feastol {
                        degenerate_col_dual = cpp_min(degenerate_col_dual, -dual / val);
                        if degenerate_col_dual < threshold {
                            break;
                        }
                    }
                } else if val < -eps {
                    if ub - self.lpsolver.solution_.col_value[i_col as usize] > feastol {
                        degenerate_col_dual = cpp_min(degenerate_col_dual, -dual / val);
                        if degenerate_col_dual < threshold {
                            break;
                        }
                    }
                }
            }

            if degenerate_col_dual < threshold {
                continue;
            }

            if degenerate_col_dual == kHighsInf && getdualproof {
                let mut rhs = HighsCDouble::from(0.0);
                for i in 0..self.row_ep.count {
                    let i_row = self.row_ep.index[i as usize];
                    let lb = self.lpsolver.get_lp().row_lower_[i_row as usize];
                    let ub = self.lpsolver.get_lp().row_upper_[i_row as usize];

                    let val = sign * self.row_ep.array[i_row as usize];
                    if ub == lb || val > eps {
                        rhs += val * ub;
                    } else if val < -eps {
                        rhs += val * lb;
                    } else {
                        rhs += val * self.lpsolver.solution_.row_value[i_row as usize];
                    }
                }

                let nnz = self.row_ap.nonzeroinds.len();
                self.dualproofvals.resize(nnz, 0.0);
                for i in 0..nnz {
                    let ci = self.row_ap.nonzeroinds[i];
                    self.dualproofvals[i] = sign * self.row_ap.get_value(ci);
                }

                let mut domchg = HighsDomainChange::default();
                domchg.column = var;
                if sign == 1.0 {
                    domchg.boundtype = HighsBoundType::kUpper;
                    domchg.boundval = self.lpsolver.get_lp().col_lower_[var as usize];
                } else {
                    domchg.boundtype = HighsBoundType::kLower;
                    domchg.boundval = self.lpsolver.get_lp().col_upper_[var as usize];
                }

                localdom.conflict_analyze_reconvergence(
                    ms,
                    &domchg,
                    &self.row_ap.nonzeroinds,
                    &self.dualproofvals,
                    nnz as HighsInt,
                    f64::from(rhs),
                    conflictpool,
                    globaldom,
                    pseudocost,
                );

                continue;
            }

            self.lpsolver.solution_.col_dual[var as usize] = sign * degenerate_col_dual;
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::computeBestEstimate (lines 428-450)
    /// C++: `double computeBestEstimate(const HighsPseudocost & ps) const`
    pub fn compute_best_estimate(&self, ms: &HighsMipSolver, ps: &HighsPseudocost) -> f64 {
        let mut estimate = HighsCDouble::from(self.objective);

        if !self.fractionalints.is_empty() {
            // because the pseudocost may be zero, we add an offset to the pseudocost so
            // that we always have some part of the estimate depending on the
            // fractionality.

            let mut increase = HighsCDouble::from(0.0);
            let offset = ms.mipdata_.feastol * cpp_max(self.objective.abs(), 1.0)
                / (ms.mipdata_.integral_cols.len() as f64);

            for &(col, frac) in self.fractionalints.iter() {
                increase += cpp_min(
                    ps.get_pseudocost_up(col, frac, offset),
                    ps.get_pseudocost_down(col, frac, offset),
                );
            }

            estimate += f64::from(increase);
        }

        f64::from(estimate)
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::computeLPDegneracy (lines 452-507)
    /// C++: `double computeLPDegneracy(const HighsDomain & localdomain) const`
    pub fn compute_lp_degneracy(&self, localdomain: &HighsDomain) -> f64 {
        if !self.lpsolver.get_solution().dual_valid || !self.lpsolver.get_basis().valid {
            return 1.0;
        }

        let dual_feas_tol = self.lpsolver.get_info().max_dual_infeasibility;

        let basis = self.lpsolver.get_basis();
        let sol = self.lpsolver.get_solution();
        let lp = self.lpsolver.get_lp();

        let mut num_fixed_rows: HighsInt = 0;
        let mut num_inequalities: HighsInt = 0;
        let mut num_basic_equalities: HighsInt = 0;

        for i in 0..self.num_rows() as usize {
            if lp.row_lower_[i] != lp.row_upper_[i] {
                num_inequalities += 1;

                if basis.row_status[i] != HighsBasisStatus::kBasic
                    && sol.row_dual[i].abs() > dual_feas_tol
                {
                    num_fixed_rows += 1;
                }
            } else {
                num_basic_equalities +=
                    (basis.row_status[i] == HighsBasisStatus::kBasic) as HighsInt;
            }
        }

        let mut num_already_fixed_cols: HighsInt = 0;
        let mut num_fixed_cols: HighsInt = 0;
        for i in 0..self.num_cols() as usize {
            if basis.col_status[i] != HighsBasisStatus::kBasic {
                if sol.col_dual[i].abs() > dual_feas_tol {
                    num_fixed_cols += 1;
                } else if localdomain.col_lower_[i] == localdomain.col_upper_[i] {
                    num_already_fixed_cols += 1;
                }
            }
        }

        let base: HighsInt =
            self.num_cols() - num_already_fixed_cols + num_inequalities + num_basic_equalities
                - self.num_rows();

        let degenerate_column_share: f64 = if base > 0 {
            1.0 - ((num_fixed_cols + num_fixed_rows) as f64) / (base as f64)
        } else {
            0.0
        };

        let var_cons_ratio: f64 = if self.num_rows() > 0 {
            ((self.num_cols() + num_inequalities + num_basic_equalities
                - num_fixed_cols
                - num_fixed_rows
                - num_already_fixed_cols) as f64)
                / (self.num_rows() as f64)
        } else {
            1.0
        };
        let fac1: f64 = if degenerate_column_share < 0.8 {
            1.0
        } else {
            10.0_f64.powf(10.0 * (degenerate_column_share - 0.7))
        };
        let fac2: f64 = if var_cons_ratio < 2.0 {
            1.0
        } else {
            10.0 * var_cons_ratio
        };

        fac1 * fac2
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::addCuts (lines 509-534)
    /// C++: `void addCuts(HighsCutSet & cutset)`
    pub fn add_cuts(&mut self, cutset: &mut HighsCutSet) {
        let numcuts = cutset.num_cuts();
        if numcuts > 0 {
            self.restore_local_mods();
        }
        debug_assert!(
            self.lpsolver.get_lp().num_row_ == self.lpsolver.get_lp().row_lower_.len() as HighsInt
        );
        debug_assert!(self.lpsolver.get_lp().num_row_ == self.lprows.len() as HighsInt);
        if numcuts > 0 {
            self.status = HighsLpRelaxation_Status::kNotSet;
            self.currentbasisstored = false;
            self.basischeckpoint = None;

            self.lprows.reserve(numcuts as usize);
            for i in 0..numcuts as usize {
                self.lprows.push(HighsLpRelaxation_LpRow::cut(
                    cutset.cutindices[i],
                    cutset.cutpools[i],
                ));
            }

            let success = self.lpsolver.add_rows(
                numcuts,
                &cutset.lower_,
                &cutset.upper_,
                cutset.a_rvalue_.len() as HighsInt,
                &cutset.a_rstart_,
                &cutset.a_rindex_,
                &cutset.a_rvalue_,
            ) == HighsStatus::kOk;
            debug_assert!(success);
            let _ = success;
            debug_assert!(
                self.lpsolver.get_lp().num_row_
                    == self.lpsolver.get_lp().row_lower_.len() as HighsInt
            );
            cutset.clear();
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::removeObsoleteRows (lines 536-556)
    /// C++: `void removeObsoleteRows(bool notifyPool = true)`
    /// C++ default arguments: `notify_pool = true`
    pub fn remove_obsolete_rows(&mut self, ms: &mut HighsMipSolver, notify_pool: bool) {
        let nlprows = self.num_rows();
        let nummodelrows = self.get_num_model_rows(ms);
        let mut deletemask: Vec<HighsInt> = Vec::new();

        let mut ndelcuts: HighsInt = 0;
        for i in nummodelrows..nlprows {
            let iu = i as usize;
            if self.lprows[iu].origin == HighsLpRelaxation_LpRow::kLocal {
                continue;
            }
            debug_assert!(self.lprows[iu].origin == HighsLpRelaxation_LpRow::kCutPool);
            if self.lpsolver.get_basis().row_status[iu] == HighsBasisStatus::kBasic {
                if ndelcuts == 0 {
                    deletemask.resize(nlprows as usize, 0);
                }
                ndelcuts += 1;
                deletemask[iu] = 1;
                if notify_pool {
                    let cutpoolindex = self.lprows[iu].cutpoolindex as usize;
                    let index = self.lprows[iu].index;
                    let thread_safe = ms.parallel_lock_active();
                    ms.mipdata_.cutpools[cutpoolindex].lp_cut_removed(index, thread_safe);
                }
            }
        }

        self.remove_cuts(ms, ndelcuts, &mut deletemask);
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::removeWorkerSpecificRows (lines 558-580)
    /// C++: `void removeWorkerSpecificRows()`
    pub fn remove_worker_specific_rows(&mut self, ms: &mut HighsMipSolver) {
        let nlprows = self.num_rows();
        let nummodelrows = self.get_num_model_rows(ms);
        let mut deletemask: Vec<HighsInt> = Vec::new();

        let mut ndelcuts: HighsInt = 0;
        for i in nummodelrows..nlprows {
            let iu = i as usize;
            if self.lprows[iu].origin == HighsLpRelaxation_LpRow::kLocal {
                continue;
            }
            debug_assert!(self.lprows[iu].origin == HighsLpRelaxation_LpRow::kCutPool);
            if self.lprows[iu].cutpoolindex > 0 {
                if ndelcuts == 0 {
                    deletemask.resize(nlprows as usize, 0);
                }
                ndelcuts += 1;
                deletemask[iu] = 1;
            }
        }

        if ndelcuts > 0 {
            let mut root_basis: HighsBasis = ms.mipdata_.firstrootbasis.clone();
            root_basis
                .row_status
                .resize(self.num_rows() as usize, HighsBasisStatus::kBasic);
            self.lpsolver.set_basis(&root_basis, "");
        }

        self.remove_cuts(ms, ndelcuts, &mut deletemask);
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::removeCuts (lines 582-609)
    /// C++: `void removeCuts(HighsInt ndelcuts, std::vector<HighsInt> & deletemask)`
    pub fn remove_cuts(
        &mut self,
        ms: &mut HighsMipSolver,
        ndelcuts: HighsInt,
        deletemask: &mut Vec<HighsInt>,
    ) {
        debug_assert!(
            self.lpsolver.get_lp().num_row_ == self.lpsolver.get_lp().row_lower_.len() as HighsInt
        );
        if ndelcuts > 0 {
            let mut basis: HighsBasis = self.lpsolver.get_basis().clone();
            let nlprows = self.lpsolver.get_num_row();
            self.lpsolver.delete_rows_mask(&mut deletemask[..]);
            for i in ms.num_row()..nlprows {
                let iu = i as usize;
                if deletemask[iu] >= 0 {
                    let dst = deletemask[iu] as usize;
                    self.lprows[dst] = self.lprows[iu];
                    basis.row_status[dst] = basis.row_status[iu];
                }
            }

            debug_assert!(
                self.lpsolver.get_lp().num_row_
                    == self.lpsolver.get_lp().row_lower_.len() as HighsInt
            );

            let new_basis_len = basis.row_status.len() - ndelcuts as usize;
            basis.row_status.truncate(new_basis_len);
            let new_lprows_len = self.lprows.len() - ndelcuts as usize;
            self.lprows.truncate(new_lprows_len);

            debug_assert!(self.lpsolver.get_lp().num_row_ == self.lprows.len() as HighsInt);
            basis.debug_origin_name = "HighsLpRelaxation::removeCuts".to_string();
            self.lpsolver.set_basis(&basis, "");
            ms.profiling_.solve_call("LP0", ms.submip);
            let _ = self.lpsolver.optimize_lp();
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::removeCuts (lines 611-627)
    /// C++: `void removeCuts()`
    pub fn remove_cuts_2(&mut self, ms: &mut HighsMipSolver) {
        debug_assert!(
            self.lpsolver.get_lp().num_row_ == self.lpsolver.get_lp().row_lower_.len() as HighsInt
        );
        let nlprows = self.lpsolver.get_num_row();
        let modelrows = ms.num_row();

        self.lpsolver.delete_rows_interval(modelrows, nlprows - 1);
        for i in modelrows..nlprows {
            let iu = i as usize;
            if self.lprows[iu].origin == HighsLpRelaxation_LpRow::kCutPool {
                let cutpoolindex = self.lprows[iu].cutpoolindex as usize;
                let index = self.lprows[iu].index;
                let thread_safe = ms.parallel_lock_active();
                ms.mipdata_.cutpools[cutpoolindex].lp_cut_removed(index, thread_safe);
            }
        }
        self.lprows.truncate(modelrows as usize);

        debug_assert!(
            self.lpsolver.get_lp().num_row_ == self.lpsolver.get_lp().row_lower_.len() as HighsInt
        );
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::performAging (lines 629-677)
    /// C++: `void performAging(bool deleteRows = false)`
    /// C++ default arguments: `delete_rows = false`
    pub fn perform_aging(&mut self, ms: &mut HighsMipSolver, delete_rows: bool) {
        debug_assert!(
            self.lpsolver.get_lp().num_row_ == self.lpsolver.get_lp().row_lower_.len() as HighsInt
        );
        if self.lpsolver.get_info().basis_validity == kBasisValidityInvalid
            || self.lpsolver.get_info().max_dual_infeasibility > ms.mipdata_.feastol
            || !self.lpsolver.get_solution().dual_valid
        {
            return;
        }

        let agelimit: HighsInt = if delete_rows {
            let mut limit = ms.options_mip_.mip_lp_age_limit;

            self.epochs += 1;
            // HiGHS quirk: `epochs` is size_t, so the modulus is unsigned; the divisor is
            // the HighsInt `max(agelimit >> 1, 2)` converted to size_t.
            if !(self.epochs as u64).is_multiple_of(cpp_max(limit >> 1, 2) as u64) {
                limit = kHighsIInf;
            } else if (self.epochs as HighsInt) < limit {
                limit = self.epochs as HighsInt;
            }
            limit
        } else {
            if self.last_age_call == self.numlpiters {
                return;
            }
            kHighsIInf
        };

        self.last_age_call = self.numlpiters;

        let nlprows = self.num_rows();
        let nummodelrows = self.get_num_model_rows(ms);
        let mut deletemask: Vec<HighsInt> = Vec::new();

        let mut ndelcuts: HighsInt = 0;
        for i in nummodelrows..nlprows {
            let iu = i as usize;
            if self.lprows[iu].origin == HighsLpRelaxation_LpRow::kLocal {
                continue;
            }
            debug_assert!(self.lprows[iu].origin == HighsLpRelaxation_LpRow::kCutPool);
            if self.lpsolver.get_basis().row_status[iu] == HighsBasisStatus::kBasic {
                // HiGHS quirk: `age += (deleteRows || age != 0)` adds a bool.
                let bump = (delete_rows || self.lprows[iu].age != 0) as HighsInt;
                self.lprows[iu].age += bump;
                if self.lprows[iu].age > agelimit {
                    if ndelcuts == 0 {
                        deletemask.resize(nlprows as usize, 0);
                    }
                    ndelcuts += 1;
                    deletemask[iu] = 1;
                    let cutpoolindex = self.lprows[iu].cutpoolindex as usize;
                    let index = self.lprows[iu].index;
                    let thread_safe = ms.parallel_lock_active();
                    ms.mipdata_.cutpools[cutpoolindex].lp_cut_removed(index, thread_safe);
                }
            } else if self.lpsolver.get_solution().row_dual[iu].abs()
                > self.lpsolver.get_options().dual_feasibility_tolerance
            {
                self.lprows[iu].age = 0;
            }
        }

        self.remove_cuts(ms, ndelcuts, &mut deletemask);
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::resetAges (lines 679-698)
    /// C++: `void resetAges()`
    pub fn reset_ages(&mut self, ms: &mut HighsMipSolver) {
        debug_assert!(
            self.lpsolver.get_lp().num_row_ == self.lpsolver.get_lp().row_lower_.len() as HighsInt
        );

        if self.lpsolver.get_info().basis_validity == kBasisValidityInvalid
            || self.lpsolver.get_info().max_dual_infeasibility > ms.mipdata_.feastol
            || !self.lpsolver.get_solution().dual_valid
        {
            return;
        }

        let nlprows = self.num_rows();
        let nummodelrows = self.get_num_model_rows(ms);

        for i in nummodelrows..nlprows {
            let iu = i as usize;
            if self.lprows[iu].origin == HighsLpRelaxation_LpRow::kLocal {
                continue;
            }
            debug_assert!(self.lprows[iu].origin == HighsLpRelaxation_LpRow::kCutPool);
            if self.lpsolver.get_basis().row_status[iu] != HighsBasisStatus::kBasic
                && self.lpsolver.get_solution().row_dual[iu].abs()
                    > self.lpsolver.get_options().dual_feasibility_tolerance
            {
                self.lprows[iu].age = 0;
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::notifyCutPoolsLpCopied (lines 700-709)
    /// C++: `void notifyCutPoolsLpCopied(HighsInt n)`
    pub fn notify_cut_pools_lp_copied(&mut self, ms: &mut HighsMipSolver, n: HighsInt) {
        let nlprows = self.num_rows();
        let modelrows = ms.num_row();
        for i in modelrows..nlprows {
            let iu = i as usize;
            if self.lprows[iu].origin == HighsLpRelaxation_LpRow::kCutPool {
                let cutpoolindex = self.lprows[iu].cutpoolindex as usize;
                let index = self.lprows[iu].index;
                ms.mipdata_.cutpools[cutpoolindex].increase_num_lps(index, n);
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::flushDomain (lines 711-730)
    /// C++: `void flushDomain(HighsDomain & domain, bool continuous = false)`
    /// C++ default arguments: `continuous = false`
    // PORT NOTE: the C++ `&domain == &mipsolver.mipdata_->domain` identity test is
    // `ms.mipdata_.domains[0].holds(domain)` (the caller detaches the global domain to pass
    // it here) or `std::ptr::eq` with the attached `domains[0]`.
    pub fn flush_domain(
        &mut self,
        ms: &HighsMipSolver,
        domain: &mut HighsDomain,
        continuous: bool,
    ) {
        let mut continuous = continuous;
        if !domain.get_changed_cols().is_empty() {
            // C++ `&domain == &mipsolver.mipdata_->getDomain()`; the global domain is slot 0.
            if ms.mipdata_.domains[0].holds(domain as *const HighsDomain) {
                continuous = true;
            }
            self.currentbasisstored = false;
            if !continuous {
                domain.remove_continuous_changed_cols(ms);
            }
            let num_chg_cols = domain.get_changed_cols().len() as HighsInt;
            if num_chg_cols == 0 {
                return;
            }
            let chg_cols: Vec<HighsInt> = domain.get_changed_cols().to_vec();
            for i in 0..num_chg_cols as usize {
                let col = chg_cols[i] as usize;
                self.col_lb_buffer[i] = domain.col_lower_[col];
                self.col_ub_buffer[i] = domain.col_upper_[col];
            }

            self.lpsolver.change_cols_bounds_set(
                num_chg_cols,
                &chg_cols,
                &self.col_lb_buffer,
                &self.col_ub_buffer,
            );

            domain.clear_changed_cols();
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::computeDualProof (lines 732-926)
    /// C++: `bool computeDualProof(const HighsDomain & globaldomain, double upperbound, std::vector<HighsInt> & inds, std::vector<double> & vals, double & rhs, bool extractCliques = true) const`
    /// C++ default arguments: `extract_cliques = true`
    // PORT NOTE: `globaldomain` is an index into `ms.mipdata_.domains`. The body reads the
    // attached domain and clique table, then calls
    // `detach!(ms.mipdata_.cliquetable, |c| c.extract_cliques_from_cut(ms, ..))`.
    /// HiGHS latest (#3223, 7db1c8180): HighsLpRelaxation.cpp::HighsLpRelaxation::computeDualProof (lines 820-916)
    #[allow(clippy::too_many_arguments)]
    pub fn compute_dual_proof(
        &self,
        ms: &mut HighsMipSolver,
        globaldomain: usize,
        upperbound: f64,
        inds: &mut Vec<HighsInt>,
        vals: &mut Vec<f64>,
        rhs: &mut f64,
        extract_cliques: bool,
    ) -> bool {
        // HiGHS quirk: this is the `#else` branch of HighsLpRelaxation.cpp; the `#if 0` branch is
        // not compiled in C++ and is not ported.
        // HiGHS #3223 (mip_cdouble_proofs): HighsLpRelaxation.cpp lines 820-916 at 7db1c8180
        let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;

        if self.has_local_mods() {
            return false;
        }
        let mut row_dual: Vec<f64> = self.lpsolver.get_solution().row_dual.clone();

        let lp = self.lpsolver.get_lp();

        debug_assert!(upperbound.is_finite());
        let mut upper = HighsCDouble::from(upperbound);

        for i in 0..lp.num_row_ as usize {
            // @FlipRowDual row_dual[i] < 0 became row_dual[i] > 0
            if row_dual[i] > 0.0 {
                if lp.row_lower_[i] != -kHighsInf {
                    // @FlipRowDual += became -=
                    if x1b_cdouble {
                        upper -= HighsCDouble::from(row_dual[i]) * lp.row_lower_[i];
                    } else {
                        upper -= row_dual[i] * lp.row_lower_[i];
                    }
                } else {
                    row_dual[i] = 0.0;
                }
                // @FlipRowDual row_dual[i] > 0 became row_dual[i] < 0
            } else if row_dual[i] < 0.0 {
                if lp.row_upper_[i] != kHighsInf {
                    // @FlipRowDual += became -=
                    if x1b_cdouble {
                        upper -= HighsCDouble::from(row_dual[i]) * lp.row_upper_[i];
                    } else {
                        upper -= row_dual[i] * lp.row_upper_[i];
                    }
                } else {
                    row_dual[i] = 0.0;
                }
            }
        }

        inds.clear();
        vals.clear();
        inds.reserve(lp.num_col_ as usize);
        vals.reserve(lp.num_col_ as usize);

        let gd: &HighsDomain = &ms.mipdata_.domains[globaldomain];
        for i in 0..lp.num_col_ as usize {
            let start = lp.a_matrix_.start_[i];
            let end = lp.a_matrix_.start_[i + 1];

            let mut sum = HighsCDouble::from(lp.col_cost_[i]);

            for j in start..end {
                let ju = j as usize;
                let row = lp.a_matrix_.index_[ju] as usize;
                if row_dual[row] == 0.0 {
                    continue;
                }
                // @FlipRowDual += became -=
                if x1b_cdouble {
                    sum -= HighsCDouble::from(lp.a_matrix_.value_[ju]) * row_dual[row];
                } else {
                    sum -= lp.a_matrix_.value_[ju] * row_dual[row];
                }
            }

            let exact_val = sum;
            let val = f64::from(sum);

            if x1b_cdouble {
                if HighsDomain::term_is_zero(
                    val,
                    gd.col_lower_[i],
                    gd.col_upper_[i],
                    ms.mipdata_.epsilon,
                ) {
                    continue;
                }
            } else if val.abs() <= ms.options_mip_.small_matrix_value {
                continue;
            }

            let mut remove_value = val.abs() <= ms.mipdata_.feastol;

            if !remove_value
                && (gd.col_lower_[i] == gd.col_upper_[i] || ms.is_col_continuous(i as HighsInt))
            {
                if val > 0.0 {
                    remove_value = self.lpsolver.get_solution().col_value[i] - gd.col_lower_[i]
                        <= ms.mipdata_.feastol;
                } else {
                    remove_value = gd.col_upper_[i] - self.lpsolver.get_solution().col_value[i]
                        <= ms.mipdata_.feastol;
                }
            }

            if remove_value {
                if val < 0.0 {
                    if gd.col_upper_[i] == kHighsInf {
                        return false;
                    }
                    if x1b_cdouble {
                        upper -= exact_val * gd.col_upper_[i];
                    } else {
                        upper -= val * gd.col_upper_[i];
                    }
                } else {
                    if gd.col_lower_[i] == -kHighsInf {
                        return false;
                    }

                    if x1b_cdouble {
                        upper -= exact_val * gd.col_lower_[i];
                    } else {
                        upper -= val * gd.col_lower_[i];
                    }
                }

                continue;
            }

            vals.push(val);
            inds.push(i as HighsInt);
        }

        *rhs = f64::from(upper);
        debug_assert!(rhs.is_finite());
        let len = inds.len() as HighsInt;
        gd.tighten_coefficients(&*ms, &mut inds[..], &mut vals[..], len, rhs);

        ms.mipdata_
            .debug_solution
            .check_cut(&inds[..], &vals[..], inds.len() as HighsInt, *rhs);

        if extract_cliques && !ms.parallel_lock_active() {
            let len = inds.len() as HighsInt;
            let rhs_value = *rhs;
            detach!(ms.mipdata_.cliquetable, |cliques| cliques
                .extract_cliques_from_cut(
                    ms,
                    &inds[..],
                    &vals[..],
                    len,
                    rhs_value
                ));
        }

        true
    }
}

impl HighsLpRelaxation_LpRow {
    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::LpRow::get (lines 89-99)
    /// C++: `void get(const HighsMipSolver & mipsolver, HighsInt & len, const HighsInt *& inds, const double *& vals) const`
    pub fn get<'s>(
        &self,
        ms: &'s HighsMipSolver,
        len: &mut HighsInt,
        inds: &mut &'s [HighsInt],
        vals: &mut &'s [f64],
    ) {
        match self.origin {
            HighsLpRelaxation_LpRow::kCutPool => ms.mipdata_.cutpools[self.cutpoolindex as usize]
                .get_cut(self.index, len, inds, vals),
            _ => ms.get_row(self.index, len, inds, vals),
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::LpRow::getRowLen (lines 101-113)
    /// C++: `HighsInt getRowLen(const HighsMipSolver & mipsolver) const`
    pub fn get_row_len(&self, ms: &HighsMipSolver) -> HighsInt {
        match self.origin {
            HighsLpRelaxation_LpRow::kCutPool => {
                ms.mipdata_.cutpools[self.cutpoolindex as usize].get_row_length(self.index)
            }
            _ => {
                ms.mipdata_.a_rstart_[(self.index + 1) as usize]
                    - ms.mipdata_.a_rstart_[self.index as usize]
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::LpRow::isIntegral (lines 115-126)
    /// C++: `bool isIntegral(const HighsMipSolver & mipsolver) const`
    pub fn is_integral(&self, ms: &HighsMipSolver) -> bool {
        match self.origin {
            HighsLpRelaxation_LpRow::kCutPool => {
                ms.mipdata_.cutpools[self.cutpoolindex as usize].cut_is_integral(self.index)
            }
            HighsLpRelaxation_LpRow::kLocal => false,
            _ => ms.mipdata_.rowintegral[self.index as usize] != 0,
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::LpRow::getMaxAbsVal (lines 128-139)
    /// C++: `double getMaxAbsVal(const HighsMipSolver & mipsolver) const`
    pub fn get_max_abs_val(&self, ms: &HighsMipSolver) -> f64 {
        match self.origin {
            HighsLpRelaxation_LpRow::kCutPool => {
                ms.mipdata_.cutpools[self.cutpoolindex as usize].get_max_abs_cut_coef(self.index)
            }
            _ => ms.mipdata_.max_abs_row_coef[self.index as usize],
        }
    }
}
