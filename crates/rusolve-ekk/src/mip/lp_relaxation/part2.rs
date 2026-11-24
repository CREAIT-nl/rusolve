//! HiGHS source: mip/HighsLpRelaxation.cpp (lines 927-1663)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// C++ `worker_` / `parallelLockActive()` test used by `storeDualInfProof`, `resolveLp` and
/// the `trySolution` / `addIncumbent` branches: `worker_ && parallelLockActive()`.
fn active_worker(lp: &HighsLpRelaxation, ms: &HighsMipSolver) -> Option<usize> {
    match lp.worker_ {
        Some(w) if ms.parallel_lock_active() => Some(w),
        _ => None,
    }
}

/// C++ `(worker_ && parallelLockActive()) ? worker_->getGlobalDomain() : getDomain()` as an index
/// into `ms.mipdata_.domains` (`HighsMipWorker::globaldom_`, else `domains[0]`).
fn global_domain_index(lp: &HighsLpRelaxation, ms: &HighsMipSolver) -> usize {
    match active_worker(lp, ms) {
        Some(w) => ms.mipdata_.workers[w].globaldom_,
        None => 0,
    }
}

/// C++ `globaldomain.col_lower_[col]` / `col_upper_[col]` for the global domain `g`. `domain` is
/// the caller's domain, which may be the detached `domains[0]` (the caller detaches it for
/// `resolveLp`). Then the identity test (`holds`) answers the read, as C++ reads the same object
/// through `getDomain()`.
fn global_col_bounds(
    ms: &HighsMipSolver,
    domain: Option<&HighsDomain>,
    g: usize,
    col: usize,
) -> (f64, f64) {
    let glb: &HighsDomain = match domain {
        Some(d) if ms.mipdata_.domains[g].holds(d as *const HighsDomain) => d,
        _ => &ms.mipdata_.domains[g],
    };
    (glb.col_lower_[col], glb.col_upper_[col])
}

/// C++ `domain` pointer passed to a call that re-enters through `ms` (`addIncumbent`,
/// `trySolution`, `checkSolution`, `run`): the detached `domain` is lent back to `domains` for
/// the call. A no-op for a domain outside `domains`.
fn with_lent_domain<R>(
    ms: &mut HighsMipSolver,
    domain: Option<&mut HighsDomain>,
    body: impl FnOnce(&mut HighsMipSolver) -> R,
) -> R {
    match domain {
        Some(d) => lend_in!(ms.mipdata_.domains, d, body(&mut *ms)),
        None => body(ms),
    }
}

impl HighsLpRelaxation {
    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::storeDualInfProof (lines 928-1045)
    /// C++: `void storeDualInfProof()`
    // PORT NOTE: the global domain (`worker_ && parallelLockActive()` ? the worker's :
    // `domains[0]`) is read attached; `extractCliquesFromCut` runs as
    // `detach!(ms.mipdata_.cliquetable, |c| c.extract_cliques_from_cut(ms, ..))`. `self` must
    // be detached by the caller (it is: `run` is `&mut self`).
    /// HiGHS latest (#3223, 7db1c8180): HighsLpRelaxation.cpp::HighsLpRelaxation::storeDualInfProof (lines 940-1039)
    pub fn store_dual_inf_proof(&mut self, ms: &mut HighsMipSolver) {
        // HiGHS #3223 (mip_cdouble_proofs): HighsLpRelaxation.cpp lines 940-1039 at 7db1c8180
        let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;

        debug_assert!(matches!(
            self.lpsolver.get_model_status(),
            HighsModelStatus::kInfeasible
        ));
        self.hasdualproof = false;
        if self.lpsolver.get_info().basis_validity == kBasisValidityInvalid {
            return;
        }
        let num_row = self.lpsolver.get_num_row();
        let num_col = self.lpsolver.get_num_col();

        if self.row_ep.size < num_row {
            self.row_ep.setup(num_row);

            if (self.row_ap.values.len() as HighsInt) < num_col {
                self.row_ap.set_dimension(num_col);
                self.dualproofvals.reserve(num_col as usize);
                self.dualproofinds.reserve(num_col as usize);
            }
        }

        self.lpsolver
            .get_dual_ray_sparse(&mut self.hasdualproof, &mut self.row_ep);

        if !self.hasdualproof {
            crate::highs_log_dev!(
                &ms.options_mip_.log_options,
                HighsLogType::kVerbose,
                "no dual ray stored\n"
            );
            return;
        }

        self.dualproofinds.clear();
        self.dualproofvals.clear();
        self.dualproofrhs = kHighsInf;

        debug_assert!(self.hasdualproof);

        let mut upper = HighsCDouble::from(0.0);
        self.row_ap.clear();

        for i in 0..self.row_ep.count as usize {
            let i_row = self.row_ep.index[i];
            let weight = -self.row_ep.array[i_row as usize];
            if weight.abs() * self.get_max_abs_row_val(ms, i_row) <= ms.mipdata_.epsilon {
                continue;
            } else if weight > 0.0 {
                if self.lpsolver.get_lp().row_upper_[i_row as usize] == kHighsInf {
                    continue;
                }
                if x1b_cdouble {
                    upper += HighsCDouble::from(weight)
                        * self.lpsolver.get_lp().row_upper_[i_row as usize];
                } else {
                    upper += weight * self.lpsolver.get_lp().row_upper_[i_row as usize];
                }
            } else {
                if self.lpsolver.get_lp().row_lower_[i_row as usize] == -kHighsInf {
                    continue;
                }
                if x1b_cdouble {
                    upper += HighsCDouble::from(weight)
                        * self.lpsolver.get_lp().row_lower_[i_row as usize];
                } else {
                    upper += weight * self.lpsolver.get_lp().row_lower_[i_row as usize];
                }
            }

            let mut len: HighsInt = 0;
            let mut inds: &[HighsInt] = &[];
            let mut vals: &[f64] = &[];
            self.get_row(ms, i_row, &mut len, &mut inds, &mut vals);

            for j in 0..len as usize {
                if x1b_cdouble {
                    self.row_ap
                        .add_cdouble(inds[j], HighsCDouble::from(weight) * vals[j]);
                } else {
                    self.row_ap.add(inds[j], weight * vals[j]);
                }
            }
        }

        let g = global_domain_index(self, ms);
        let globaldomain: &HighsDomain = &ms.mipdata_.domains[g];

        for &i in self.row_ap.get_nonzeros() {
            let iu = i as usize;
            let val;
            let exact_val;
            if x1b_cdouble {
                exact_val = self.row_ap.get_value_cd(i);
                val = f64::from(exact_val);
            } else {
                val = self.row_ap.get_value(i);
                exact_val = HighsCDouble::from(val);
            }

            if x1b_cdouble {
                if HighsDomain::term_is_zero(
                    val,
                    globaldomain.col_lower_[iu],
                    globaldomain.col_upper_[iu],
                    ms.mipdata_.epsilon,
                ) {
                    continue;
                }
            } else if val.abs() <= ms.mipdata_.epsilon {
                continue;
            }

            let mut remove_value = val.abs() <= ms.mipdata_.feastol;

            if !remove_value
                && (globaldomain.col_lower_[iu] == globaldomain.col_upper_[iu]
                    || ms.is_col_continuous(i))
            {
                // remove continuous entries and globally fixed entries whenever the
                // local LP's bound is not tighter than the global bound
                if val > 0.0 {
                    remove_value = self.lpsolver.get_lp().col_lower_[iu]
                        - globaldomain.col_lower_[iu]
                        <= ms.mipdata_.feastol;
                } else {
                    remove_value = globaldomain.col_upper_[iu]
                        - self.lpsolver.get_lp().col_upper_[iu]
                        <= ms.mipdata_.feastol;
                }
            }

            if remove_value {
                if val < 0.0 {
                    if globaldomain.col_upper_[iu] == kHighsInf {
                        self.hasdualproof = false;
                        return;
                    }
                    if x1b_cdouble {
                        upper -= exact_val * globaldomain.col_upper_[iu];
                    } else {
                        upper -= val * globaldomain.col_upper_[iu];
                    }
                } else {
                    if globaldomain.col_lower_[iu] == -kHighsInf {
                        self.hasdualproof = false;
                        return;
                    }
                    if x1b_cdouble {
                        upper -= exact_val * globaldomain.col_lower_[iu];
                    } else {
                        upper -= val * globaldomain.col_lower_[iu];
                    }
                }

                continue;
            }

            self.dualproofvals.push(val);
            self.dualproofinds.push(i);
        }

        self.dualproofrhs = f64::from(upper);
        let len = self.dualproofinds.len() as HighsInt;
        globaldomain.tighten_coefficients(
            ms,
            &mut self.dualproofinds[..],
            &mut self.dualproofvals[..],
            len,
            &mut self.dualproofrhs,
        );

        // HiGHS quirk: `debugSolution.checkCut(..)` is a HIGHS_DEBUGSOL check, dropped.

        if !ms.parallel_lock_active() {
            let len = self.dualproofinds.len() as HighsInt;
            let rhs = self.dualproofrhs;
            detach!(ms.mipdata_.cliquetable, |ct| ct.extract_cliques_from_cut(
                ms,
                &self.dualproofinds[..],
                &self.dualproofvals[..],
                len,
                rhs
            ));
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::storeDualUBProof (lines 1047-1066)
    /// C++: `void storeDualUBProof()`
    pub fn store_dual_ub_proof(&mut self, ms: &mut HighsMipSolver) {
        debug_assert!(matches!(
            self.lpsolver.get_model_status(),
            HighsModelStatus::kObjectiveBound
        ));

        self.dualproofinds.clear();
        self.dualproofvals.clear();

        if self.lpsolver.get_solution().dual_valid {
            // C++ `use_worker_info ? worker_->getGlobalDomain() : getDomain()` and
            // `use_worker_info ? worker_->upper_limit : upper_limit`.
            let g = global_domain_index(self, ms);
            let upper_limit = match active_worker(self, ms) {
                Some(w) => ms.mipdata_.workers[w].upper_limit,
                None => ms.mipdata_.upper_limit,
            };
            // `computeDualProof` is `const` and writes the proof vectors: they are taken out
            // of `self` for the call, then put back.
            let mut inds = std::mem::take(&mut self.dualproofinds);
            let mut vals = std::mem::take(&mut self.dualproofvals);
            let mut rhs = self.dualproofrhs;
            self.hasdualproof =
                self.compute_dual_proof(ms, g, upper_limit, &mut inds, &mut vals, &mut rhs, true);
            self.dualproofinds = inds;
            self.dualproofvals = vals;
            self.dualproofrhs = rhs;
        } else {
            self.hasdualproof = false;
        }

        if !self.hasdualproof {
            self.dualproofrhs = kHighsInf;
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::checkDualProof (lines 1068-1091)
    /// C++: `bool checkDualProof() const`
    pub fn check_dual_proof(&self, ms: &HighsMipSolver) -> bool {
        if !self.hasdualproof {
            return true;
        }
        if self.dualproofrhs == kHighsInf {
            return false;
        }

        let len = self.dualproofinds.len();

        let mut viol = HighsCDouble::from(-self.dualproofrhs);

        let lp = self.lpsolver.get_lp();

        for i in 0..len {
            let col = self.dualproofinds[i] as usize;
            if self.dualproofvals[i] > 0.0 {
                if lp.col_lower_[col] == -kHighsInf {
                    return false;
                }
                viol += self.dualproofvals[i] * lp.col_lower_[col];
            } else {
                debug_assert!(self.dualproofvals[i] < 0.0);
                if lp.col_upper_[col] == kHighsInf {
                    return false;
                }
                viol += self.dualproofvals[i] * lp.col_upper_[col];
            }
        }

        viol > ms.mipdata_.feastol
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::computeDualInfProof (lines 1093-1105)
    /// C++: `bool computeDualInfProof(const HighsDomain & globaldomain, std::vector<HighsInt> & inds, std::vector<double> & vals, double & rhs) const`
    pub fn compute_dual_inf_proof(
        &self,
        globaldomain: &HighsDomain,
        inds: &mut Vec<HighsInt>,
        vals: &mut Vec<f64>,
        rhs: &mut f64,
    ) -> bool {
        // rusolve mip_node_bigm: proofs from a node-tightened LP are not globally valid
        if !self.hasdualproof || self.has_local_mods() {
            return false;
        }

        debug_assert!(self.dualproofrhs.is_finite());

        inds.clone_from(&self.dualproofinds);
        vals.clone_from(&self.dualproofvals);
        *rhs = self.dualproofrhs;
        true
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::recoverBasis (lines 1107-1112)
    /// C++: `void recoverBasis()`
    pub fn recover_basis(&mut self) {
        if let Some(basis) = self.basischeckpoint.clone() {
            let status = self
                .lpsolver
                .set_basis(&basis, "HighsLpRelaxation::recoverBasis");
            // rusolve lever `mip_lp_keep_dse_weights` = 2 (not in HiGHS): the weights stored
            // with the basis replace the ones kept from the last solve (simplex/ekk/dse_keep.rs).
            if let Some(weight) = &basis.dse_weight {
                if status == HighsStatus::kOk && self.lpsolver.options_.mip_lp_keep_dse_weights == 2
                {
                    self.lpsolver
                        .ekk_instance_
                        .dse_keep_set(weight, basis.col_status.len());
                }
            }
            self.currentbasisstored = true;
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::setObjectiveLimit (lines 1114-1123)
    /// C++: `void setObjectiveLimit(double objlim = kHighsInf)`
    /// C++ default arguments: `objlim = kHighsInf`
    pub fn set_objective_limit(&mut self, ms: &HighsMipSolver, objlim: f64) {
        let offset = if ms.mipdata_.objective_function.is_integral() {
            0.5 / ms.mipdata_.objective_function.integral_scale()
        } else {
            cpp_max(1000.0 * ms.mipdata_.feastol, objlim.abs() * kHighsTiny)
        };

        self.lpsolver
            .set_option_value_double("objective_bound", objlim + offset);
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::run (lines 1125-1439)
    /// C++: `Status run(bool resolve_on_error = true)`
    /// C++ default arguments: `resolve_on_error = true`
    // PORT NOTE: `ipm.setProfiling` is not ported; the `ipm_logging` / `dump_ipm_lp`
    // blocks are dead (`const bool ... = false`).
    // PORT NOTE: the `profiling_->mip_` "MIP-Timing" logs are dropped: `HighsProfiling` is a
    // stub without that flag.
    pub fn run(
        &mut self,
        ms: &mut HighsMipSolver,
        resolve_on_error: bool,
    ) -> HighsLpRelaxation_Status {
        let this_time_limit = cpp_max(
            self.lpsolver.get_run_time() + ms.options_mip_.time_limit - ms.timer_.read(),
            0.0,
        );
        self.lpsolver
            .set_option_value_double("time_limit", this_time_limit);
        let mut valid_basis = self.lpsolver.get_basis().valid;
        // rusolve mip_root_lp_concurrent (`root_race.rs`): race dual simplex
        // against IPM + crossover for the first LP; the winner's basis makes the solve below warm.
        if !valid_basis
            && ms.options_mip_.mip_root_lp_concurrent
            && ms.options_mip_.threads >= 2
            && !ms.submip
            && !self.solved_first_lp
            && self.race_root_lp(ms, this_time_limit)
        {
            valid_basis = true;
        }

        // Determine the solver (cpp 1158-1185).
        let solver = self
            .lpsolver
            .get_option_value_str("solver")
            .unwrap_or_default();
        let use_solver: String = if valid_basis {
            kSimplexString.to_string()
        } else {
            let mip_lp_solver = ms.options_mip_.mip_lp_solver.clone();
            if use_ipm(&mip_lp_solver) {
                let use_hipo = mip_lp_solver == kHipoString;
                if use_hipo {
                    kHipoString.to_string()
                } else {
                    kIpxString.to_string()
                }
            } else {
                kSimplexString.to_string()
            }
        };
        let mut callstatus = HighsStatus::kOk;
        // Now allowing the use of IPM at the root node
        self.lpsolver.set_option_value_str("solver", &use_solver);
        let use_ipm = use_ipm(&use_solver);
        let mut use_simplex = !use_ipm;
        if use_ipm {
            debug_assert!(!valid_basis);
            ms.profiling_.solve_call("LP1", ms.submip);
            callstatus = self.lpsolver.optimize_lp();
            if callstatus == HighsStatus::kError {
                crate::highs_log_dev!(
                    &ms.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "HighsLpRelaxation::run HiPO has failed : status = {} Try IPX\n",
                    self.lpsolver
                        .model_status_to_string(self.lpsolver.get_model_status())
                );
                self.lpsolver.set_option_value_str("solver", kSimplexString);
                use_simplex = true;
            }
        }
        if use_simplex {
            // The profiling sub-MIP bookkeeping (cpp 1203-1211) is a stub.
            ms.profiling_.solve_call("LP2", ms.submip);
            callstatus = self.lpsolver.optimize_lp();
        }

        // Revert the value of lpsolver.options_.solver
        self.lpsolver.set_option_value_str("solver", &solver);

        self.solved_first_lp = true;
        let info = self.lpsolver.get_info().clone();
        let mut itercount: HighsInt = -1;
        if use_simplex {
            itercount = cpp_max(0, info.simplex_iteration_count);
            self.numlpiters += itercount as i64;
        }
        // reporting only
        crate::work::add_lp_solve(cpp_max(0, itercount) as u64);

        if callstatus == HighsStatus::kError {
            self.lpsolver.clear_solver();
            // HiGHS quirk: the `#if 0` retry with primal simplex from the last basis (cpp
            // 1300-1313) is not compiled in HiGHS and is not ported.

            if resolve_on_error {
                // still an error: now try to solve with presolve from scratch
                self.lpsolver
                    .set_option_value_int("simplex_strategy", kSimplexStrategyDual);
                self.lpsolver
                    .set_option_value_str("presolve", kHighsOnString);
                let retval = self.run(ms, false);
                self.lpsolver
                    .set_option_value_str("presolve", kHighsOffString);

                return retval;
            }

            self.recover_basis();

            return HighsLpRelaxation_Status::kError;
        }

        let model_status = self.lpsolver.get_model_status();
        match model_status {
            HighsModelStatus::kObjectiveBound => {
                self.num_solved += 1;
                if itercount >= 0 {
                    self.avg_solve_iters +=
                        (itercount as f64 - self.avg_solve_iters) / (self.num_solved as f64);
                }

                self.store_dual_ub_proof(ms);
                HighsLpRelaxation_Status::kInfeasible
            }
            HighsModelStatus::kInfeasible => {
                self.num_solved += 1;
                if itercount >= 0 {
                    self.avg_solve_iters +=
                        (itercount as f64 - self.avg_solve_iters) / (self.num_solved as f64);
                }

                self.store_dual_inf_proof(ms);
                // HiGHS quirk: `if (true || checkDualProof()) return Status::kInfeasible;`
                // short-circuits, so the proof is never checked and the rest of this case
                // (C++ 1383-1438, `hasdualproof = false` and the commented-out re-solves) is
                // unreachable. It is not ported.
                HighsLpRelaxation_Status::kInfeasible
            }
            HighsModelStatus::kUnbounded => {
                // If unboundedness is detected in the presolved LP, then postsolve cannot be run,
                // so there is no basis. Returning Status::kError as a result yielded #1962, where
                // the root node is unbounded.
                if info.basis_validity == kBasisValidityInvalid {
                    crate::highs_log_user!(
                        &ms.options_mip_.log_options,
                        HighsLogType::kWarning,
                        "HighsLpRelaxation::run LP is unbounded with no basis, but not returning Status::kError\n"
                    );
                }
                if info.primal_solution_status == kSolutionStatusFeasible {
                    let sol = self.lpsolver.get_solution().col_value.clone();
                    let worker = active_worker(self, ms);
                    lend_in!(
                        ms.mipdata_.lps,
                        self,
                        match worker {
                            None => {
                                ms.try_solution(&sol, kSolutionSourceUnbounded);
                            }
                            Some(w) => {
                                detach!(ms.mipdata_.workers[w], |wk| wk.try_solution(
                                    ms,
                                    &sol,
                                    kSolutionSourceUnbounded
                                ));
                            }
                        }
                    );
                }

                HighsLpRelaxation_Status::kUnbounded
            }
            HighsModelStatus::kUnknown | HighsModelStatus::kOptimal => {
                // HiGHS: `case kUnknown` returns kError with no basis, then falls through.
                if model_status == HighsModelStatus::kUnknown
                    && info.basis_validity == kBasisValidityInvalid
                {
                    return HighsLpRelaxation_Status::kError;
                }
                debug_assert!(info.max_primal_infeasibility >= 0.0);
                debug_assert!(info.max_dual_infeasibility >= 0.0);
                self.num_solved += 1;
                if itercount >= 0 {
                    self.avg_solve_iters +=
                        (itercount as f64 - self.avg_solve_iters) / (self.num_solved as f64);
                }
                let feastol = ms.mipdata_.feastol;
                if info.max_primal_infeasibility <= feastol
                    && info.max_dual_infeasibility <= feastol
                {
                    return HighsLpRelaxation_Status::kOptimal;
                }

                if info.max_primal_infeasibility <= feastol {
                    return HighsLpRelaxation_Status::kUnscaledPrimalFeasible;
                }

                if info.max_dual_infeasibility <= feastol {
                    return HighsLpRelaxation_Status::kUnscaledDualFeasible;
                }

                if model_status == HighsModelStatus::kOptimal {
                    return HighsLpRelaxation_Status::kUnscaledInfeasible;
                }

                HighsLpRelaxation_Status::kError
            }
            HighsModelStatus::kIterationLimit => {
                if !ms.submip && resolve_on_error {
                    // Highs instantiation
                    let mut ipm = Highs::default();
                    ipm.set_option_value_bool("output_flag", false);
                    // h16_sub_solver_flags: the C++ `Highs ipm` runs 1.16.0 code
                    if ms.options_mip_.h16_sub_solver_flags {
                        ipm.options_.copy_version_flags_from(&ms.options_mip_);
                    }
                    // check if only root presolve is allowed
                    let use_presolve = !ms.options_mip_.mip_root_presolve_only;
                    let presolve = if use_presolve {
                        kHighsChooseString
                    } else {
                        kHighsOffString
                    };
                    ipm.set_option_value_str("presolve", presolve);
                    // Determine the solver
                    let mip_ipm_solver = ms.options_mip_.mip_ipm_solver.clone();
                    // Currently use HiPO by default and take action on failure
                    // here. Later pass mip_ipm_solver and take action on failure in
                    // solveLp
                    let use_hipo = mip_ipm_solver == kHipoString;
                    // Later still, pass mip_ipm_solver and take action on failure in
                    // solveLp
                    let ipm_solver = if use_hipo { kHipoString } else { kIpxString };
                    ipm.set_option_value_str("solver", ipm_solver);
                    ipm.set_option_value_int("ipm_iteration_limit", 200);
                    ipm.pass_model(self.lpsolver.get_lp().clone());
                    // todo @ Julian : If you remove this you can see the looping on
                    // istanbul-no-cutoff
                    ipm.set_option_value_int(
                        "simplex_iteration_limit",
                        info.simplex_iteration_count,
                    );
                    ms.profiling_.solve_call("LP3", ms.submip);
                    ipm.optimize_lp();
                    if use_hipo && !ipm.get_basis().valid {
                        // HiPO has failed to get a solution, so try IPX
                        crate::highs_log_dev!(
                            &ms.options_mip_.log_options,
                            HighsLogType::kInfo,
                            "HighsLpRelaxation::run HiPO has failed to get a valid basis: status = {} Try IPX\n",
                            ipm.model_status_to_string(ipm.get_model_status())
                        );
                        ipm.set_option_value_str("solver", kIpxString);
                        ms.profiling_.solve_call("LP4", ms.submip);
                        ipm.optimize_lp();
                    }
                    // Port-only: compare with the oracle `lp_ipm` record when
                    // RUSOLVE_MIP_AC_FIXTURE is set (mip/ipm_replay.rs).
                    ipm_replay_check(
                        kMipIpmRecordLpIpm,
                        ms.submip_level,
                        false,
                        ipm.get_model_status(),
                        ipm.get_solution(),
                        ipm.get_basis(),
                    );
                    self.lpsolver
                        .set_basis(ipm.get_basis(), "HighsLpRelaxation::run IPM basis");
                    return self.run(ms, false);
                }

                // printf("error: lpsolver reached iteration limit\n");
                HighsLpRelaxation_Status::kError
            }
            HighsModelStatus::kTimeLimit => HighsLpRelaxation_Status::kError,
            _ => {
                crate::highs_log_user!(
                    &ms.options_mip_.log_options,
                    HighsLogType::kWarning,
                    "LP solved to unexpected status: {}\n",
                    self.lpsolver.model_status_to_string(model_status)
                );
                HighsLpRelaxation_Status::kError
            }
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation::resolveLp (lines 1441-1662)
    /// C++: `Status resolveLp(HighsDomain * domain = nullptr)`
    /// C++ default arguments: `domain = nullptr`
    // PORT NOTE: `domain` is detached by the caller (the node domain, or `domains[0]` detached
    // for the call). `addIncumbent` / `trySolution` re-enter through `ms` (they read every LP
    // and domain), so the body wraps them as
    // `lend_in!(ms.mipdata_.lps, self, lend_in!(ms.mipdata_.domains, d, ms.add_incumbent(..)))`
    // after cloning the solution it passes (`self` is `Default` while lent). Under
    // `parallelLockActive()` with `worker_ = Some(w)` the C++ calls the worker's
    // `addIncumbent` / `trySolution` instead; same lending.
    // PORT NOTE: the C++ `const HighsDomain & globaldom` is read per access through
    // `global_col_bounds`, because the detached `domain` may be that same domain.
    pub fn resolve_lp(
        &mut self,
        ms: &mut HighsMipSolver,
        domain: Option<&mut HighsDomain>,
    ) -> HighsLpRelaxation_Status {
        let mut domain = domain;
        self.fractionalints.clear();
        // HiGHS 1.16 #3149 (h16_cdouble): HighsLpRelaxation.cpp::resolveLp objsum loops (1.16.0)
        let h16_cdouble_obj = ms.options_mip_.h16_cdouble;

        let mut solveagain = true;
        while solveagain {
            solveagain = false;
            if let Some(d) = domain.as_deref_mut() {
                self.flush_domain(ms, d, false);
            }
            let status = with_lent_domain(ms, domain.as_deref_mut(), |ms| self.run(ms, true));
            self.status = status;

            match self.status {
                HighsLpRelaxation_Status::kUnscaledInfeasible
                | HighsLpRelaxation_Status::kUnscaledDualFeasible
                | HighsLpRelaxation_Status::kUnscaledPrimalFeasible
                | HighsLpRelaxation_Status::kOptimal => {
                    let mut fracints: HighsHashTable<HighsInt, (f64, i32)> =
                        HighsHashTable::with_min_capacity(self.max_num_fractional as u64);

                    let mut objsum = HighsCDouble::from(0.0);
                    let mut roundable = true;
                    let g = global_domain_index(self, ms);

                    for k in 0..ms.mipdata_.integral_cols.len() {
                        let i = ms.mipdata_.integral_cols[k];
                        let iu = i as usize;

                        // for the fractionality we assume that LP bounds are not violated
                        // bounds that are violated by the unscaled LP are indicated by the
                        // return status already
                        let mut val = cpp_max(
                            cpp_min(
                                self.lpsolver.get_solution().col_value[iu],
                                self.lpsolver.get_lp().col_upper_[iu],
                            ),
                            self.lpsolver.get_lp().col_lower_[iu],
                        );

                        if fractionality(val, None) > ms.mipdata_.feastol {
                            let mut col = i;
                            roundable = roundable
                                && (ms.mipdata_.uplocks[iu] == 0 || ms.mipdata_.downlocks[iu] == 0);

                            let mut subst = ms.mipdata_.cliquetable.get_substitution(col).cloned();
                            while let Some(sub) = subst {
                                let replace_col = sub.replace.col() as usize;
                                if self.lpsolver.get_lp().col_lower_[replace_col]
                                    == self.lpsolver.get_lp().col_upper_[replace_col]
                                {
                                    match domain.as_deref_mut() {
                                        Some(d) => {
                                            let weight = sub
                                                .replace
                                                .weight(&self.lpsolver.get_lp().col_lower_);
                                            d.fix_col(
                                                ms,
                                                col,
                                                weight,
                                                HighsDomain_Reason::unspecified(),
                                            );
                                            if d.infeasible() {
                                                self.objective = kHighsInf;
                                                self.status = HighsLpRelaxation_Status::kInfeasible;
                                                return self.status;
                                            }
                                        }
                                        None => break,
                                    }
                                }

                                // Spot cases where a global domain change during the search has
                                // led to a clique and the local domain has branched into a state
                                // where the clique subst. would be invalid for the local domain.
                                // TODO: Turn into assert when search becomes parallel
                                if let Some(d) = domain.as_deref() {
                                    let replace_val = if sub.replace.val() == 0 {
                                        1.0 - val
                                    } else {
                                        val
                                    };
                                    let replace_lb = d.col_lower_[replace_col];
                                    let replace_ub = d.col_upper_[replace_col];
                                    let feastol = ms.mipdata_.feastol;
                                    if replace_val < replace_lb - feastol
                                        || replace_val > replace_ub + feastol
                                    {
                                        break;
                                    }
                                }

                                col = sub.replace.col() as HighsInt;
                                if sub.replace.val() == 0 {
                                    val = 1.0 - val;
                                }

                                subst = ms.mipdata_.cliquetable.get_substitution(col).cloned();
                            }

                            if self.adjust_sym_branching_col {
                                col = ms.mipdata_.symmetries.get_branching_column(
                                    &self.lpsolver.get_lp().col_lower_,
                                    &self.lpsolver.get_lp().col_upper_,
                                    col,
                                );
                            }

                            let pair = fracints.index_or_insert(&col);
                            pair.0 += val;
                            pair.1 += 1;
                        } else {
                            if self.lpsolver.get_basis().col_status[iu] == HighsBasisStatus::kBasic
                            {
                                continue;
                            }

                            let (glb, gub) = global_col_bounds(ms, domain.as_deref(), g, iu);
                            let sol_val = self.lpsolver.get_solution().col_value[iu];
                            if cpp_min(gub - sol_val, sol_val - glb) <= ms.mipdata_.feastol {
                                continue;
                            }

                            let matrix = &self.lpsolver.get_lp().a_matrix_;
                            let model_start = &mip_model!(ms).a_matrix_.start_;
                            let col_start =
                                matrix.start_[iu] + (model_start[iu + 1] - model_start[iu]);
                            let col_end = matrix.start_[iu + 1];

                            // skip further checks if the column has no entry in any cut
                            if col_start == col_end {
                                continue;
                            }

                            for j in col_start..col_end {
                                let row = matrix.index_[j as usize];
                                debug_assert!(row >= ms.num_row());

                                // age is already zero, so irrelevant whether we reset it
                                if self.lprows[row as usize].age == 0 {
                                    continue;
                                }

                                // check that row has no slack in the current solution
                                if self.lpsolver.get_solution().row_value[row as usize]
                                    < self.lpsolver.get_lp().row_upper_[row as usize]
                                        - ms.mipdata_.feastol
                                {
                                    continue;
                                }

                                self.lprows[row as usize].age = 0;
                            }
                        }
                    }

                    self.max_num_fractional =
                        cpp_max(fracints.size() as HighsInt, self.max_num_fractional);

                    if domain
                        .as_deref()
                        .is_some_and(|d| !d.get_changed_cols().is_empty())
                    {
                        // resolving due to fixings of substituted columns
                        solveagain = true;
                        continue;
                    }

                    for entry in fracints.iter() {
                        self.fractionalints
                            .push((*entry.key(), entry.value().0 / (entry.value().1 as f64)));
                    }

                    if roundable && !self.fractionalints.is_empty() {
                        let mut roundsol: Vec<f64> = self.lpsolver.get_solution().col_value.clone();

                        for k in 0..self.fractionalints.len() {
                            // get column index
                            let (col, second) = self.fractionalints[k];
                            let cu = col as usize;
                            let lp = self.lpsolver.get_lp();
                            // round based on locks and sign of objective coefficient
                            if ms.mipdata_.uplocks[cu] == 0
                                && (ms.col_cost_2(col) < 0.0 || ms.mipdata_.downlocks[cu] != 0)
                            {
                                // round up
                                roundsol[cu] = cpp_min(
                                    (second - ms.mipdata_.feastol).ceil(),
                                    if lp.col_upper_[cu] == kHighsInf {
                                        kHighsInf
                                    } else {
                                        (lp.col_upper_[cu] + ms.mipdata_.feastol).floor()
                                    },
                                );
                            } else {
                                // round down
                                roundsol[cu] = cpp_max(
                                    (second + ms.mipdata_.feastol).floor(),
                                    if lp.col_lower_[cu] == -kHighsInf {
                                        -kHighsInf
                                    } else {
                                        (lp.col_lower_[cu] - ms.mipdata_.feastol).ceil()
                                    },
                                );
                            }
                        }

                        // C++ `getSubstitutions()` is the non-const accessor: copy the list out.
                        let cliquesubst: Vec<HighsCliqueTable_Substitution> =
                            detach!(ms.mipdata_.cliquetable, |ct| ct.get_substitutions().clone());
                        for k in (0..cliquesubst.len()).rev() {
                            let sub = &cliquesubst[k];
                            if sub.replace.val() == 0 {
                                roundsol[sub.substcol as usize] =
                                    1.0 - roundsol[sub.replace.col() as usize];
                            } else {
                                roundsol[sub.substcol as usize] =
                                    roundsol[sub.replace.col() as usize];
                            }
                        }

                        for i in 0..ms.num_col() {
                            if h16_cdouble_obj {
                                objsum +=
                                    HighsCDouble::from(roundsol[i as usize]) * ms.col_cost_2(i);
                            } else {
                                objsum += roundsol[i as usize] * ms.col_cost_2(i);
                            }
                        }

                        let solobj = f64::from(objsum);
                        let worker = active_worker(self, ms);
                        lend_in!(
                            ms.mipdata_.lps,
                            self,
                            match worker {
                                None => {
                                    with_lent_domain(ms, domain.as_deref_mut(), |ms| {
                                        ms.add_incumbent(
                                            &roundsol,
                                            solobj,
                                            kSolutionSourceSolveLp,
                                            true,
                                            false,
                                        )
                                    });
                                }
                                Some(w) => {
                                    detach!(ms.mipdata_.workers[w], |wk| wk.add_incumbent(
                                        ms,
                                        &roundsol,
                                        solobj,
                                        kSolutionSourceSolveLp
                                    ));
                                }
                            }
                        );
                        objsum = HighsCDouble::from(0.0);
                    }

                    for i in 0..ms.num_col() {
                        if h16_cdouble_obj {
                            objsum += HighsCDouble::from(
                                self.lpsolver.get_solution().col_value[i as usize],
                            ) * ms.col_cost_2(i);
                        } else {
                            objsum += self.lpsolver.get_solution().col_value[i as usize]
                                * ms.col_cost_2(i);
                        }
                    }

                    if self.fractionalints.is_empty()
                        && !HighsLpRelaxation::unscaled_primal_feasible(self.status)
                    {
                        let mut fix_sol: Vec<f64> = self.lpsolver.get_solution().col_value.clone();
                        for i in 0..ms.num_col() {
                            let iu = i as usize;
                            let lp = self.lpsolver.get_lp();
                            if fix_sol[iu] < lp.col_lower_[iu] {
                                fix_sol[iu] = lp.col_lower_[iu];
                            } else if fix_sol[iu] > lp.col_upper_[iu] {
                                fix_sol[iu] = lp.col_upper_[iu];
                            } else if ms.is_col_integral(i) {
                                fix_sol[iu] = fix_sol[iu].round();
                            }
                        }

                        let feasible = lend_in!(
                            ms.mipdata_.lps,
                            self,
                            with_lent_domain(ms, domain.as_deref_mut(), |ms| {
                                ms.check_solution(&fix_sol)
                            })
                        );
                        if feasible {
                            // C++ `const_cast<std::vector<double>&>(sol.col_value) = std::move(fixSol)`
                            self.lpsolver.solution_.col_value = fix_sol;
                            if HighsLpRelaxation::unscaled_dual_feasible(self.status) {
                                self.status = HighsLpRelaxation_Status::kOptimal;
                            } else {
                                self.status = HighsLpRelaxation_Status::kUnscaledPrimalFeasible;
                            }
                        }
                    }

                    self.objective = f64::from(objsum);
                }
                HighsLpRelaxation_Status::kInfeasible => {
                    self.objective = kHighsInf;
                }
                _ => {}
            }
        }

        self.status
    }
}
