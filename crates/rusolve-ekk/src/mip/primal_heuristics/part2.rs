//! HiGHS source: mip/HighsPrimalHeuristics.cpp (lines 971-1862)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use std::collections::{HashMap, HashSet};

use crate::lp_data::presolve_support::calculate_row_values_quad_vec;
use crate::prelude::*;

impl HighsPrimalHeuristics {
    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::tryRoundedPoint (lines 971-1059)
    /// C++: `bool tryRoundedPoint(HighsMipWorker & worker, const std::vector<double> & point, const int solution_source)`
    pub fn try_rounded_point(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        point: &[f64],
        solution_source: i32,
    ) -> bool {
        // C++ `HighsDomain localdom = worker.getGlobalDomain();`, released on every exit.
        let mut localdom = local_copy_of_global_domain(ms, worker);
        let result = self.try_rounded_point_body(ms, worker, &mut localdom, point, solution_source);
        localdom.release(ms);
        result
    }

    fn try_rounded_point_body(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        localdom: &mut HighsDomain,
        point: &[f64],
        solution_source: i32,
    ) -> bool {
        let mut integer_feasible = true;

        let numintcols: HighsInt = self.intcols.len() as HighsInt;
        for i in 0..numintcols as usize {
            let col = self.intcols[i];
            let mut intval = point[col as usize];
            let mut rounded: f64 = 0.0;
            // check if solution value of integer-constrained variable is actually integral
            let feasible = fractionality(intval, Some(&mut rounded)) <= ms.mipdata_.feastol;
            integer_feasible = integer_feasible && feasible;
            if !feasible {
                continue;
            }
            // use rounded solution value and check against bounds
            intval = rounded;
            intval = cpp_min(localdom.col_upper_[col as usize], intval);
            intval = cpp_max(localdom.col_lower_[col as usize], intval);

            localdom.fix_col(ms, col, intval, HighsDomain_Reason::branching());
            if localdom.infeasible() {
                conflict_analysis_for_worker(ms, worker, localdom);
                return false;
            }
            let _ = localdom.propagate(ms);
            if localdom.infeasible() {
                conflict_analysis_for_worker(ms, worker, localdom);
                return false;
            }
        }

        if numintcols != ms.num_col() {
            // LP relaxation instantiation
            let mut lprelax = HighsLpRelaxation::new(ms);
            lprelax.set_mip_worker(worker);
            lprelax.set_profiling();
            lprelax.load_model(ms);
            lprelax.set_iteration_limit(
                cpp_max(10000i64, 2 * ms.mipdata_.firstrootlpiters) as HighsInt
            );
            let _ = lprelax.get_lp_solver().change_cols_bounds_interval(
                0,
                ms.num_col() - 1,
                &localdom.col_lower_,
                &localdom.col_upper_,
            );

            // check if only root presolve is allowed
            if ms.options_mip_.mip_root_presolve_only {
                let _ = lprelax
                    .get_lp_solver()
                    .set_option_value_str("presolve", kHighsOffString);
            }
            if !ms.options_mip_.mip_root_presolve_only && (5 * numintcols) / ms.num_col() >= 1 {
                let _ = lprelax
                    .get_lp_solver()
                    .set_option_value_str("presolve", kHighsOnString);
            } else {
                let _ = lprelax.get_lp_solver().set_basis(
                    &ms.mipdata_.firstrootbasis,
                    "HighsPrimalHeuristics::tryRoundedPoint",
                );
            }

            let st = lprelax.resolve_lp(ms, None);

            if st == HighsLpRelaxation_Status::kInfeasible {
                let mut inds: Vec<HighsInt> = Vec::new();
                let mut vals: Vec<f64> = Vec::new();
                let mut rhs: f64 = 0.0;
                let gd = ms.mipdata_.workers[worker].get_global_domain();
                if lprelax.compute_dual_inf_proof(
                    &ms.mipdata_.domains[gd],
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                ) {
                    let cutpool = ms.mipdata_.workers[worker].get_cut_pool();
                    let mut cut_gen = HighsCutGeneration::new(ms, &lprelax, cutpool);
                    let _ = cut_gen.generate_conflict(
                        ms, &lprelax, &*localdom, gd, &mut inds, &mut vals, &mut rhs,
                    );
                }
                return false;
            } else if HighsLpRelaxation::unscaled_primal_feasible(st) {
                let lpsol: Vec<f64> = lprelax.get_lp_solver().get_solution().col_value.clone();
                if !integer_feasible {
                    // there may be fractional integer variables -> try ziRound heuristic
                    self.zi_round(ms, worker, &lpsol);
                    return self.try_solution(ms, &lpsol, solution_source, worker);
                } else {
                    // all integer variables are fixed -> add incumbent
                    let objective = lprelax.get_objective();
                    self.add_incumbent(ms, &lpsol, objective, solution_source, worker);
                    return true;
                }
            }
        }

        self.try_solution(ms, &localdom.col_lower_, solution_source, worker)
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::linesearchRounding (lines 1061-1115)
    /// C++: `bool linesearchRounding(HighsMipWorker & worker, const std::vector<double> & point1, const std::vector<double> & point2, const int solution_source)`
    pub fn linesearch_rounding(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        point1: &[f64],
        point2: &[f64],
        solution_source: i32,
    ) -> bool {
        let numintcols: HighsInt = self.intcols.len() as HighsInt;
        let mut roundedpoint: Vec<f64> = vec![0.0; ms.num_col() as usize];

        let mut alpha: f64 = 0.0;
        debug_assert!(ms.mipdata_.uplocks.len() as HighsInt == ms.num_col());
        debug_assert!(point1.len() as HighsInt == ms.num_col());
        debug_assert!(point2.len() as HighsInt == ms.num_col());

        while alpha < 1.0 {
            let mut nextalpha: f64 = 1.0;
            let mut reachedpoint2 = true;
            for i in 0..numintcols as usize {
                let col = self.intcols[i] as usize;
                debug_assert!(col < ms.num_col() as usize);
                let feastol = ms.mipdata_.feastol;
                if ms.mipdata_.uplocks[col] == 0 {
                    roundedpoint[col] = (cpp_max(point1[col], point2[col]) - feastol).ceil();
                    continue;
                }

                if ms.mipdata_.downlocks[col] == 0 {
                    roundedpoint[col] = (cpp_min(point1[col], point2[col]) + feastol).floor();
                    continue;
                }

                let convexcomb = (1.0 - alpha) * point1[col] + alpha * point2[col];
                let intpoint2 = (point2[col] + 0.5).floor();
                roundedpoint[col] = (convexcomb + 0.5).floor();

                if roundedpoint[col] == intpoint2 {
                    continue;
                }

                reachedpoint2 = false;
                let tmpalpha = (roundedpoint[col] + 0.5 + feastol - point1[col])
                    / (point2[col] - point1[col]).abs();
                if tmpalpha < nextalpha && tmpalpha > alpha + 1e-2 {
                    nextalpha = tmpalpha;
                }
            }

            if self.try_rounded_point(ms, worker, &roundedpoint, solution_source) {
                return true;
            }

            if reachedpoint2 {
                return false;
            }

            alpha = nextalpha;
        }

        false
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::randomizedRounding (lines 1117-1202)
    /// C++: `void randomizedRounding(HighsMipWorker & worker, const std::vector<double> & relaxationsol)`
    pub fn randomized_rounding(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        relaxationsol: &[f64],
    ) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_ROUNDING);
        if relaxationsol.len() != ms.num_col() as usize {
            return;
        }

        let locked = ms.parallel_lock_active();
        let mut randgen = self.take_randgen(ms, worker, locked);
        // C++ `HighsDomain localdom = worker.getGlobalDomain();`, released on every exit.
        let mut localdom = local_copy_of_global_domain(ms, worker);
        self.randomized_rounding_body(ms, worker, relaxationsol, &mut localdom, &mut randgen);
        localdom.release(ms);
        self.put_randgen(ms, worker, locked, randgen);
    }

    fn randomized_rounding_body(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        relaxationsol: &[f64],
        localdom: &mut HighsDomain,
        randgen: &mut HighsRandom,
    ) {
        for &i in self.intcols.iter() {
            let iu = i as usize;
            let feastol = ms.mipdata_.feastol;
            let mut intval: f64;
            if ms.mipdata_.uplocks[iu] == 0 {
                intval = (relaxationsol[iu] - feastol).ceil();
            } else if ms.mipdata_.downlocks[iu] == 0 {
                intval = (relaxationsol[iu] + feastol).floor();
            } else {
                intval = (relaxationsol[iu] + randgen.real(0.1, 0.9)).floor();
            }

            intval = cpp_min(localdom.col_upper_[iu], intval);
            intval = cpp_max(localdom.col_lower_[iu], intval);

            localdom.fix_col(ms, i, intval, HighsDomain_Reason::branching());
            if localdom.infeasible() {
                conflict_analysis_for_worker(ms, worker, localdom);
                return;
            }
            let _ = localdom.propagate(ms);
            if localdom.infeasible() {
                conflict_analysis_for_worker(ms, worker, localdom);
                return;
            }
        }

        if ms.mipdata_.integer_cols.len() != ms.num_col() as usize {
            // LP relaxation instantiation
            let mut lprelax = HighsLpRelaxation::new(ms);
            lprelax.set_mip_worker(worker);
            lprelax.set_profiling();
            lprelax.load_model(ms);
            lprelax.set_iteration_limit(
                cpp_max(10000i64, 2 * ms.mipdata_.firstrootlpiters) as HighsInt
            );
            let _ = lprelax.get_lp_solver().change_cols_bounds_interval(
                0,
                ms.num_col() - 1,
                &localdom.col_lower_,
                &localdom.col_upper_,
            );

            // check if only root presolve is allowed
            if ms.options_mip_.mip_root_presolve_only {
                let _ = lprelax
                    .get_lp_solver()
                    .set_option_value_str("presolve", kHighsOffString);
            }

            if !ms.options_mip_.mip_root_presolve_only
                && (5 * self.intcols.len()) / (ms.num_col() as usize) >= 1
            {
                // LP to solve is very much smaller, so use presolve rather than the root basis
                let _ = lprelax
                    .get_lp_solver()
                    .set_option_value_str("presolve", kHighsOnString);
            } else {
                let _ = lprelax.get_lp_solver().set_basis(
                    &ms.mipdata_.firstrootbasis,
                    "HighsPrimalHeuristics::randomizedRounding",
                );
            }

            let st = lprelax.resolve_lp(ms, None);

            if st == HighsLpRelaxation_Status::kInfeasible {
                let mut inds: Vec<HighsInt> = Vec::new();
                let mut vals: Vec<f64> = Vec::new();
                let mut rhs: f64 = 0.0;
                let gd = ms.mipdata_.workers[worker].get_global_domain();
                if lprelax.compute_dual_inf_proof(
                    &ms.mipdata_.domains[gd],
                    &mut inds,
                    &mut vals,
                    &mut rhs,
                ) {
                    let cutpool = ms.mipdata_.workers[worker].get_cut_pool();
                    let mut cut_gen = HighsCutGeneration::new(ms, &lprelax, cutpool);
                    let _ = cut_gen.generate_conflict(
                        ms, &lprelax, &*localdom, gd, &mut inds, &mut vals, &mut rhs,
                    );
                }
            } else if HighsLpRelaxation::unscaled_primal_feasible(st) {
                let lpsol: Vec<f64> = lprelax.get_lp_solver().get_solution().col_value.clone();
                let objective = lprelax.get_objective();
                self.add_incumbent(
                    ms,
                    &lpsol,
                    objective,
                    kSolutionSourceRandomizedRounding,
                    worker,
                );
            }
        } else {
            let sol = localdom.col_lower_.clone();
            self.try_solution(ms, &sol, kSolutionSourceRandomizedRounding, worker);
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::shifting (lines 1204-1463)
    /// C++: `void shifting(HighsMipWorker & worker, const std::vector<double> & relaxationsol)`
    pub fn shifting(&mut self, ms: &mut HighsMipSolver, worker: usize, relaxationsol: &[f64]) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_ROUNDING);
        if relaxationsol.len() != ms.num_col() as usize {
            return;
        }

        let locked = ms.parallel_lock_active();
        let mut randgen = self.take_randgen(ms, worker, locked);
        self.shifting_body(ms, worker, relaxationsol, &mut randgen);
        self.put_randgen(ms, worker, locked, randgen);
    }

    fn shifting_body(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        relaxationsol: &[f64],
        randgen: &mut HighsRandom,
    ) {
        let mut current_relax_solution: Vec<f64> = relaxationsol.to_vec();
        let mut t: HighsInt = 0;
        let current_lp: &HighsLp = mip_model!(ms);
        // LP relaxation instantiation
        let lp_index = ms.mipdata_.workers[worker].lp_;
        let mut lprelax = HighsLpRelaxation::copy_from(ms, &ms.mipdata_.lps[lp_index]);
        lprelax.set_mip_worker(worker);
        lprelax.set_profiling();
        let mut current_fractional_integers: Vec<(HighsInt, f64)> =
            lprelax.get_fractional_integers().to_vec();
        let mut current_infeasible_rows: Vec<(HighsInt, HighsInt, f64)> =
            ms.get_infeasible_rows(&current_relax_solution);
        let mut previous_infeasible_rows_size: usize = current_infeasible_rows.len();
        let mut has_infeasible_constraints: bool = !current_infeasible_rows.is_empty();
        let mut iterations_without_reductions: HighsInt = 0;
        let max_iterations_without_reductions: HighsInt = 5;
        let mut shift_iterations_set: HashMap<HighsInt, Vec<HighsInt>> = HashMap::new();

        while (!current_fractional_integers.is_empty() || has_infeasible_constraints)
            && iterations_without_reductions <= max_iterations_without_reductions
            && t <= ms.mipdata_.integer_cols.len() as HighsInt
        {
            t += 1;
            // rusolve: every pass rescans the infeasible rows, so a large model overruns the
            // time limit; honour it.
            if t % 16 == 0 && ms.timer_.read() >= ms.options_mip_.time_limit {
                break;
            }
            let mut fractional_integers_reduced = false;
            iterations_without_reductions += 1;
            if has_infeasible_constraints {
                // find an infeasible row that has a non-zero coefficient a(row,col) where
                // col in currentFractInt
                let mut fractional_integer_found = false;
                let mut r_index: HighsInt = 0;
                while !fractional_integer_found
                    && r_index != current_infeasible_rows.len() as HighsInt
                {
                    let r = current_infeasible_rows[r_index as usize].0 as usize;
                    let start = ms.mipdata_.a_rstart_[r];
                    let end = ms.mipdata_.a_rstart_[r + 1];
                    let mut jind = start;
                    while jind != end {
                        let col = ms.mipdata_.a_rindex_[jind as usize];
                        fractional_integer_found =
                            current_fractional_integers.iter().any(|p| p.0 == col);
                        if fractional_integer_found {
                            break;
                        }
                        jind += 1;
                    }
                    r_index += 1;
                }
                let mut row_index: HighsInt = r_index - 1;
                if !fractional_integer_found {
                    // otherwise select a random infeasible row
                    row_index = randgen.integer_sup(current_infeasible_rows.len() as HighsInt);
                }

                let row = current_infeasible_rows[row_index as usize].0;
                let row_sense = current_infeasible_rows[row_index as usize].1;
                let infeasibility = current_infeasible_rows[row_index as usize].2;
                let mut score_min: f64 = kHighsInf;
                let mut j_min: HighsInt = HighsInt::MAX;
                let mut x_j_min: f64 = kHighsInf;
                let mut aij_min: f64 = 0.0;
                let mut move_value_up = false;
                let start = ms.mipdata_.a_rstart_[row as usize];
                let end = ms.mipdata_.a_rstart_[row as usize + 1];
                let is_max = ms.orig_model_.sense_ == ObjSense::kMaximize;
                for jind in start as usize..end as usize {
                    let j = ms.mipdata_.a_rindex_[jind];
                    let ju = j as usize;

                    // skip fixed variables
                    if current_lp.col_lower_[ju] == current_lp.col_upper_[ju] {
                        continue;
                    }

                    // lambda for finding best shift
                    let mut repair = |col: HighsInt,
                                      direction: HighsInt,
                                      row_sense: HighsInt,
                                      num_locks: HighsInt,
                                      coef: f64,
                                      cost: f64,
                                      is_maximization: bool,
                                      is_integer: bool,
                                      is_at_bound: bool| {
                        let dir = direction as f64;
                        if (row_sense < 0 || dir * coef > 0.0)
                            && (row_sense > 0 || dir * coef < 0.0)
                        {
                            return;
                        }

                        // skip variables at bounds
                        if is_at_bound {
                            return;
                        }

                        // search for column
                        let found = current_fractional_integers.iter().any(|p| p.0 == col);

                        // add data
                        let mut score: f64;
                        if found {
                            score = -1.0 + 1.0 / ((num_locks + 1) as f64);
                        } else {
                            match shift_iterations_set.get(&col) {
                                Some(shifts) if !shifts.is_empty() => {
                                    score = 0.0;
                                    for &shift in shifts {
                                        let ds = dir * (shift as f64);
                                        if ds > 0.0 {
                                            score += 1.1f64.powf(ds - (t as f64));
                                        }
                                    }
                                }
                                _ => {
                                    score = dir * (if is_maximization { -cost } else { cost });
                                }
                            }
                            if is_integer {
                                score += 1.0;
                            }
                        }

                        if score < score_min {
                            score_min = score;
                            j_min = col;
                            aij_min = coef;
                            x_j_min = current_relax_solution[col as usize];
                            move_value_up = direction > 0;
                        }
                    };

                    // repair with up-rounding
                    repair(
                        j,
                        1,
                        row_sense,
                        ms.mipdata_.downlocks[ju],
                        ms.mipdata_.a_rvalue_[jind],
                        current_lp.col_cost_[ju],
                        is_max,
                        current_lp.integrality_[ju] == HighsVarType::kInteger,
                        (current_lp.col_upper_[ju] - current_relax_solution[ju]).abs()
                            <= ms.mipdata_.feastol,
                    );

                    // repair with down-rounding
                    repair(
                        j,
                        -1,
                        row_sense,
                        ms.mipdata_.uplocks[ju],
                        ms.mipdata_.a_rvalue_[jind],
                        current_lp.col_cost_[ju],
                        is_max,
                        current_lp.integrality_[ju] == HighsVarType::kInteger,
                        (current_relax_solution[ju] - current_lp.col_lower_[ju]).abs()
                            <= ms.mipdata_.feastol,
                    );
                }

                if j_min != HighsInt::MAX {
                    // Update current_fractional_integers
                    if let Some(pos) = current_fractional_integers
                        .iter()
                        .position(|p| p.0 == j_min)
                    {
                        current_fractional_integers.remove(pos);
                        fractional_integers_reduced = true;
                    }
                    // Update current_relax_solution and shift_iterations_set (for not fractional integers)
                    let jm = j_min as usize;
                    let feastol = ms.mipdata_.feastol;
                    if move_value_up {
                        if fractional_integers_reduced {
                            current_relax_solution[jm] = (x_j_min - feastol).ceil();
                        } else {
                            if current_lp.integrality_[jm] == HighsVarType::kInteger {
                                // variable is integer and not at the upper bound, so increment by 1.
                                current_relax_solution[jm] = x_j_min + 1.0;
                            } else {
                                current_relax_solution[jm] = cpp_min(
                                    x_j_min + infeasibility / aij_min.abs(),
                                    current_lp.col_upper_[jm] + feastol,
                                );
                            }
                            shift_iterations_set.entry(j_min).or_default().push(t);
                        }
                    } else if fractional_integers_reduced {
                        current_relax_solution[jm] = (x_j_min + feastol).floor();
                    } else {
                        if current_lp.integrality_[jm] == HighsVarType::kInteger {
                            // variable is integer and not at the lower bound, so decrement by 1.
                            current_relax_solution[jm] = x_j_min - 1.0;
                        } else {
                            current_relax_solution[jm] = cpp_max(
                                x_j_min - infeasibility / aij_min.abs(),
                                current_lp.col_lower_[jm] - feastol,
                            );
                        }

                        shift_iterations_set.entry(j_min).or_default().push(-t);
                    }
                }
            } else {
                let mut xi_max: f64 = -1.0;
                let mut delta_c_min: f64 = kHighsInf;
                let mut pind_j_min: HighsInt = HighsInt::MAX;
                let mut j_min: HighsInt = HighsInt::MAX;
                let mut x_j_min: f64 = kHighsInf;
                let mut sigma: HighsInt = 0;
                for i in 0..current_fractional_integers.len() {
                    let it = current_fractional_integers[i];
                    let col = it.0;
                    debug_assert!(col >= 0);
                    debug_assert!(col < ms.num_col());

                    let mut is_better =
                        |col: HighsInt, xi: f64, rounded_val: f64, direction: HighsInt| {
                            let c_min = current_lp.col_cost_[col as usize] * (rounded_val - it.1);
                            if xi > xi_max || (xi == xi_max && c_min < delta_c_min) {
                                xi_max = xi;
                                delta_c_min = c_min;
                                pind_j_min = i as HighsInt;
                                j_min = col;
                                x_j_min = rounded_val;
                                sigma = direction;
                            }
                        };

                    let feastol = ms.mipdata_.feastol;
                    is_better(
                        col,
                        ms.mipdata_.uplocks[col as usize] as f64,
                        (it.1 + feastol).floor(),
                        -1,
                    );
                    is_better(
                        col,
                        ms.mipdata_.downlocks[col as usize] as f64,
                        (it.1 - feastol).ceil(),
                        1,
                    );
                }
                if sigma != 0 {
                    current_relax_solution[j_min as usize] = x_j_min;
                }
                if pind_j_min != HighsInt::MAX {
                    current_fractional_integers.remove(pind_j_min as usize);
                    fractional_integers_reduced = true;
                }
            }
            current_infeasible_rows = ms.get_infeasible_rows(&current_relax_solution);
            has_infeasible_constraints = !current_infeasible_rows.is_empty();
            if current_infeasible_rows.len() < previous_infeasible_rows_size
                || fractional_integers_reduced
            {
                iterations_without_reductions = 0;
            }
            previous_infeasible_rows_size = current_infeasible_rows.len();
        }
        // re-check for feasibility and add incumbent
        if has_infeasible_constraints {
            self.try_rounded_point(ms, worker, &current_relax_solution, kSolutionSourceShifting);
        } else if !current_fractional_integers.is_empty() {
            self.zi_round(ms, worker, &current_relax_solution);
        } else {
            self.try_solution(ms, &current_relax_solution, kSolutionSourceShifting, worker);
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::ziRound (lines 1465-1576)
    /// C++: `void ziRound(HighsMipWorker & worker, const std::vector<double> & relaxationsol)`
    pub fn zi_round(&mut self, ms: &mut HighsMipSolver, worker: usize, relaxationsol: &[f64]) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_ROUNDING);
        // if (mipsolver.submip) return;
        if relaxationsol.len() != ms.num_col() as usize {
            return;
        }

        let mut current_relax_solution: Vec<f64> = relaxationsol.to_vec();

        let feastol = ms.mipdata_.feastol;
        let zi = |x: f64| -> f64 { cpp_min((x - feastol).ceil() - x, x - (x + feastol).floor()) };

        let mut zi_total: HighsCDouble = HighsCDouble::from(0.0);
        for &i in self.intcols.iter() {
            zi_total += zi(current_relax_solution[i as usize]);
        }

        if zi_total <= feastol {
            return;
        }

        let current_lp: &HighsLp = mip_model!(ms);

        let num_row = current_lp.num_row_ as usize;
        let mut row_activities: Vec<f64> = vec![0.0; num_row];
        let mut xrow_lower: Vec<f64> = vec![0.0; num_row];
        let mut xrow_upper: Vec<f64> = vec![0.0; num_row];

        let mut loop_count: HighsInt = 0;
        let max_loop_count: HighsInt = 5;
        let mut previous_zi_total: HighsCDouble;
        let mut improvement_in_feasibility: HighsCDouble = HighsCDouble::from(kHighsInf);

        'zi_loops: while zi_total > feastol
            && improvement_in_feasibility > feastol
            && loop_count <= max_loop_count
        {
            previous_zi_total = zi_total;
            loop_count += 1;

            // getLpRowBounds(currentLp, 0, num_row - 1, ...), inlined: no-op when there are no rows.
            if current_lp.num_row_ > 0 {
                xrow_lower.copy_from_slice(&current_lp.row_lower_[..num_row]);
                xrow_upper.copy_from_slice(&current_lp.row_upper_[..num_row]);
            }

            let mut n_frac_cols: HighsInt = 0;
            for &j in self.intcols.iter() {
                let ju = j as usize;
                let relax_solution = current_relax_solution[ju];
                if fractionality(relax_solution, None) <= feastol {
                    continue;
                }
                // rusolve: each fractional column recomputes every row activity, so the loop
                // is quadratic on large models; honour the time limit.
                n_frac_cols += 1;
                if n_frac_cols % 16 == 0 && ms.timer_.read() >= ms.options_mip_.time_limit {
                    break 'zi_loops;
                }

                row_activities.clear();
                row_activities.resize(num_row, 0.0);
                calculate_row_values_quad_vec(
                    current_lp,
                    &current_relax_solution,
                    &mut row_activities,
                    -1,
                );

                let mut min_row_ratio_for_upper = kHighsInf;
                let mut min_row_ratio_for_lower = kHighsInf;

                for el in current_lp.a_matrix_.start_[ju]..current_lp.a_matrix_.start_[ju + 1] {
                    let i = current_lp.a_matrix_.index_[el as usize] as usize;
                    let aij = current_lp.a_matrix_.value_[el as usize];

                    let slack_upper = xrow_upper[i] - row_activities[i];
                    let slack_lower = row_activities[i] - xrow_lower[i];
                    min_row_ratio_for_upper = cpp_min(
                        min_row_ratio_for_upper,
                        (if aij > 0.0 { slack_upper } else { -slack_lower }) / aij,
                    );
                    min_row_ratio_for_lower = cpp_min(
                        min_row_ratio_for_lower,
                        (if aij > 0.0 { slack_lower } else { -slack_upper }) / aij,
                    );
                }

                let upper_bound = cpp_min(
                    current_lp.col_upper_[ju] - relax_solution,
                    min_row_ratio_for_upper,
                );
                let lower_bound = cpp_min(
                    relax_solution - current_lp.col_lower_[ju],
                    min_row_ratio_for_lower,
                );

                let mut perform_updates = |col: HighsInt, change: f64| {
                    let old_relax_solution = current_relax_solution[col as usize];
                    current_relax_solution[col as usize] += change;
                    zi_total = zi_total - zi(old_relax_solution)
                        + zi(current_relax_solution[col as usize]);
                };

                if (zi(relax_solution + upper_bound) - zi(relax_solution - lower_bound)).abs()
                    <= feastol
                    && zi(relax_solution + upper_bound) < zi(relax_solution)
                {
                    let xcol_cost = current_lp.col_cost_[ju];
                    let ub_obj_change_smaller = xcol_cost * (relax_solution + upper_bound)
                        <= xcol_cost * (relax_solution - lower_bound);
                    let is_minimization = current_lp.sense_ == ObjSense::kMinimize;
                    if (is_minimization && ub_obj_change_smaller)
                        || (!is_minimization && !ub_obj_change_smaller)
                    {
                        perform_updates(j, upper_bound);
                    } else {
                        perform_updates(j, -lower_bound);
                    }
                } else if zi(relax_solution + upper_bound) < zi(relax_solution - lower_bound)
                    && zi(relax_solution + upper_bound) < zi(relax_solution)
                {
                    perform_updates(j, upper_bound);
                } else if zi(relax_solution + upper_bound) > zi(relax_solution - lower_bound)
                    && zi(relax_solution - lower_bound) < zi(relax_solution)
                {
                    perform_updates(j, -lower_bound);
                }
            }
            improvement_in_feasibility = previous_zi_total - zi_total;
        }
        // re-check for feasibility and add incumbent
        self.try_solution(ms, &current_relax_solution, kSolutionSourceZiRound, worker);
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::feasibilityPump (lines 1578-1695)
    /// C++: `void feasibilityPump(HighsMipWorker & worker)`
    pub fn feasibility_pump(&mut self, ms: &mut HighsMipSolver, worker: usize) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_FEAS_PUMP);
        let locked = ms.parallel_lock_active();
        let mut randgen = self.take_randgen(ms, worker, locked);
        self.feasibility_pump_body(ms, worker, &mut randgen);
        self.put_randgen(ms, worker, locked, randgen);
    }

    fn feasibility_pump_body(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        randgen: &mut HighsRandom,
    ) {
        // LP relaxation instantiation
        let lp_index = ms.mipdata_.workers[worker].lp_;
        let mut lprelax = HighsLpRelaxation::copy_from(ms, &ms.mipdata_.lps[lp_index]);
        lprelax.set_mip_worker(worker);
        lprelax.set_profiling();
        let mut referencepoints: HashSet<Vec<HighsInt>> = HashSet::new();
        let mut status = lprelax.resolve_lp(ms, None);
        *ms.mipdata_.workers[worker].get_heur_lp_iterations() += lprelax.get_num_lp_iterations();

        let numcol = ms.num_col() as usize;
        let mask: Vec<HighsInt> = vec![1; numcol];
        let mut cost: Vec<f64> = vec![0.0; numcol];

        let _ = lprelax
            .get_lp_solver()
            .set_option_value_int("simplex_strategy", kSimplexStrategyPrimal);
        lprelax.set_objective_limit(ms, kHighsInf);
        let _ = lprelax
            .get_lp_solver()
            .set_option_value_double("primal_simplex_bound_perturbation_multiplier", 0.0);

        lprelax.set_iteration_limit((5.0 * ms.mipdata_.avgrootlpiters) as HighsInt);

        while !lprelax.get_fractional_integers().is_empty() {
            let lpsol: Vec<f64> = lprelax.get_lp_solver().get_solution().col_value.clone();
            let mut roundedsol: Vec<f64> = lprelax.get_lp_solver().get_solution().col_value.clone();

            let mut referencepoint: Vec<HighsInt> =
                Vec::with_capacity(ms.mipdata_.integer_cols.len());

            // HighsDomain localdom = worker.getGlobalDomain(); local to this iteration,
            // released right after its last use (before any exit below).
            let mut localdom = local_copy_of_global_domain(ms, worker);
            for ii in 0..ms.mipdata_.integer_cols.len() {
                let i = ms.mipdata_.integer_cols[ii];
                let iu = i as usize;
                let mut intval = (roundedsol[iu] + randgen.real(0.4, 0.6)).floor();
                intval = cpp_max(intval, localdom.col_lower_[iu]);
                intval = cpp_min(intval, localdom.col_upper_[iu]);
                roundedsol[iu] = intval;
                referencepoint.push(intval as HighsInt);
                if !localdom.infeasible() {
                    localdom.fix_col(ms, i, intval, HighsDomain_Reason::branching());
                    if localdom.infeasible() {
                        conflict_analysis_for_worker(ms, worker, &mut localdom);
                        continue;
                    }
                    let _ = localdom.propagate(ms);
                    if localdom.infeasible() {
                        conflict_analysis_for_worker(ms, worker, &mut localdom);
                        continue;
                    }
                }
            }
            localdom.release(ms);

            let mut havecycle = !referencepoints.insert(referencepoint.clone());
            let mut k: HighsInt = 0;
            while havecycle && k < 2 {
                for _ in 0..10 {
                    let flippos =
                        randgen.integer_sup(ms.mipdata_.integer_cols.len() as HighsInt) as usize;
                    let col = ms.mipdata_.integer_cols[flippos] as usize;
                    let gd = ms.mipdata_.workers[worker].get_global_domain();
                    if roundedsol[col] > lpsol[col] {
                        roundedsol[col] = (lpsol[col].floor() as HighsInt) as f64;
                    } else if roundedsol[col] < lpsol[col] {
                        roundedsol[col] = (lpsol[col].ceil() as HighsInt) as f64;
                    } else if roundedsol[col] < ms.mipdata_.domains[gd].col_upper_[col] {
                        roundedsol[col] = ms.mipdata_.domains[gd].col_upper_[col];
                    } else {
                        roundedsol[col] = ms.mipdata_.domains[gd].col_lower_[col];
                    }

                    referencepoint[flippos] = roundedsol[col] as HighsInt;
                }
                havecycle = !referencepoints.insert(referencepoint.clone());
                k += 1;
            }

            if havecycle {
                return;
            }

            if self.linesearch_rounding(
                ms,
                worker,
                &lpsol,
                &roundedsol,
                kSolutionSourceFeasibilityPump,
            ) {
                return;
            }

            if (lprelax.get_num_lp_iterations() as f64) >= 1000.0 + ms.mipdata_.avgrootlpiters * 5.0
            {
                break;
            }

            let feastol = ms.mipdata_.feastol;
            for ii in 0..ms.mipdata_.integer_cols.len() {
                let i = ms.mipdata_.integer_cols[ii];
                let iu = i as usize;
                if ms.mipdata_.uplocks[iu] == 0 || ms.mipdata_.downlocks[iu] == 0 {
                    cost[iu] = 0.0;
                } else if lpsol[iu] > roundedsol[iu] - feastol {
                    cost[iu] = -1.0 + randgen.real(-1e-4, 1e-4);
                } else {
                    cost[iu] = 1.0 + randgen.real(-1e-4, 1e-4);
                }
            }

            let _ = lprelax.get_lp_solver().change_cols_cost_mask(&mask, &cost);
            let mut niters: i64 = -lprelax.get_num_lp_iterations();
            status = lprelax.resolve_lp(ms, None);
            niters += lprelax.get_num_lp_iterations();
            if niters == 0 {
                break;
            }
            *ms.mipdata_.workers[worker].get_heur_lp_iterations() += niters;
        }

        if lprelax.get_fractional_integers().is_empty()
            && HighsLpRelaxation::unscaled_primal_feasible(status)
        {
            let sol: Vec<f64> = lprelax.get_lp_solver().get_solution().col_value.clone();
            let objective = lprelax.get_objective();
            self.add_incumbent(ms, &sol, objective, kSolutionSourceFeasibilityPump, worker);
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::centralRounding (lines 1697-1714)
    /// C++: `void centralRounding(HighsMipWorker & worker)`
    pub fn central_rounding(&mut self, ms: &mut HighsMipSolver, worker: usize) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_CENTRAL);
        if ms.mipdata_.analytic_center.len() != ms.num_col() as usize {
            return;
        }

        // Clone the mipdata_ members: `linesearch_rounding` takes `ms` mutably.
        let analytic_center = ms.mipdata_.analytic_center.clone();
        if !ms.mipdata_.firstlpsol.is_empty() {
            let firstlpsol = ms.mipdata_.firstlpsol.clone();
            self.linesearch_rounding(
                ms,
                worker,
                &firstlpsol,
                &analytic_center,
                kSolutionSourceCentralRounding,
            );
        } else if !ms.mipdata_.rootlpsol.is_empty() {
            let rootlpsol = ms.mipdata_.rootlpsol.clone();
            self.linesearch_rounding(
                ms,
                worker,
                &rootlpsol,
                &analytic_center,
                kSolutionSourceCentralRounding,
            );
        } else {
            self.linesearch_rounding(
                ms,
                worker,
                &analytic_center,
                &analytic_center,
                kSolutionSourceCentralRounding,
            );
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::addIncumbent (lines 1786-1795)
    /// C++: `bool addIncumbent(const std::vector<double> & sol, double solobj, const int solution_source, HighsMipWorker & worker)`
    pub fn add_incumbent(
        &mut self,
        ms: &mut HighsMipSolver,
        sol: &[f64],
        solobj: f64,
        solution_source: i32,
        worker: usize,
    ) -> bool {
        if ms.parallel_lock_active() {
            detach!(ms.mipdata_.workers[worker], |wk| wk.add_incumbent(
                ms,
                sol,
                solobj,
                solution_source
            ))
        } else {
            ms.add_incumbent(sol, solobj, solution_source, true, false)
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::trySolution (lines 1797-1805)
    /// C++: `bool trySolution(const std::vector<double> & solution, const int solution_source, HighsMipWorker & worker)`
    pub fn try_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        solution: &[f64],
        solution_source: i32,
        worker: usize,
    ) -> bool {
        if ms.parallel_lock_active() {
            detach!(ms.mipdata_.workers[worker], |wk| wk.try_solution(
                ms,
                solution,
                solution_source
            ))
        } else {
            ms.try_solution(solution, solution_source)
        }
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::getNumSuccessObservations (lines 1807-1810)
    /// C++: `HighsInt getNumSuccessObservations(HighsMipWorker & worker) const`
    // PORT NOTE: the worker is passed by reference (`&ms.mipdata_.workers[w]`); read-only.
    pub fn get_num_success_observations(&self, worker: &HighsMipWorker) -> HighsInt {
        self.num_success_observations + worker.heur_stats.num_success_observations
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::getNumInfeasObservations (lines 1812-1815)
    /// C++: `HighsInt getNumInfeasObservations(HighsMipWorker & worker) const`
    // PORT NOTE: the worker is passed by reference (`&ms.mipdata_.workers[w]`); read-only.
    pub fn get_num_infeas_observations(&self, worker: &HighsMipWorker) -> HighsInt {
        self.num_infeas_observations + worker.heur_stats.num_infeas_observations
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::getSuccessObservations (lines 1817-1820)
    /// C++: `double getSuccessObservations(HighsMipWorker & worker) const`
    // PORT NOTE: the worker is passed by reference (`&ms.mipdata_.workers[w]`); read-only.
    pub fn get_success_observations(&self, worker: &HighsMipWorker) -> f64 {
        self.success_observations + worker.heur_stats.success_observations
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::getInfeasObservations (lines 1822-1825)
    /// C++: `double getInfeasObservations(HighsMipWorker & worker) const`
    // PORT NOTE: the worker is passed by reference (`&ms.mipdata_.workers[w]`); read-only.
    pub fn get_infeas_observations(&self, worker: &HighsMipWorker) -> f64 {
        self.infeas_observations + worker.heur_stats.infeas_observations
    }

    /// HiGHS: HighsPrimalHeuristics.cpp::HighsPrimalHeuristics::flushStatistics (lines 1827-1861)
    /// C++: `void flushStatistics(HighsMipSolver & mipsolver, HighsMipWorker & worker)`
    // PORT NOTE: every C++ call passes the heuristics' own `mipsolver`: one `ms` parameter.
    pub fn flush_statistics(&mut self, ms: &mut HighsMipSolver, worker: usize) {
        let mut total_repair_lp: i64 = 0;
        let mut total_repair_lp_feasible: i64 = 0;
        let mut total_repair_lp_iterations: i64 = 0;
        let mut lp_iterations: i64 = 0;
        let mut success_observations: f64 = 0.0;
        let mut num_success_observations: HighsInt = 0;
        let mut infeas_observations: f64 = 0.0;
        let mut num_infeas_observations: HighsInt = 0;
        let mut max_submip_level: HighsInt = 0;
        let mut termination_status = HighsModelStatus::kNotset;
        ms.mipdata_.workers[worker].get_heur_stats_values(
            &mut total_repair_lp,
            &mut total_repair_lp_feasible,
            &mut total_repair_lp_iterations,
            &mut lp_iterations,
            &mut success_observations,
            &mut num_success_observations,
            &mut infeas_observations,
            &mut num_infeas_observations,
            &mut max_submip_level,
            &mut termination_status,
        );

        ms.mipdata_.total_repair_lp += total_repair_lp;
        ms.mipdata_.total_repair_lp_feasible += total_repair_lp_feasible;
        ms.mipdata_.total_repair_lp_iterations += total_repair_lp_iterations;
        ms.mipdata_.heuristic_lp_iterations += lp_iterations;
        ms.mipdata_.total_lp_iterations += lp_iterations;
        ms.max_submip_level = cpp_max(ms.max_submip_level, max_submip_level);
        if termination_status != HighsModelStatus::kNotset
            && ms.termination_status_ == HighsModelStatus::kNotset
        {
            ms.termination_status_ = termination_status;
        }
        self.success_observations += success_observations;
        self.num_success_observations += num_success_observations;
        self.infeas_observations += infeas_observations;
        self.num_infeas_observations += num_infeas_observations;
        ms.mipdata_.workers[worker].reset_heur_stats();
    }

    /// C++ `HighsRandom& randgen = parallelLockActive() ? worker.randgen : this->randgen;`
    /// The chosen generator is moved out for the body; `put_randgen` restores it.
    fn take_randgen(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        locked: bool,
    ) -> HighsRandom {
        if locked {
            std::mem::take(&mut ms.mipdata_.workers[worker].randgen)
        } else {
            std::mem::take(&mut self.randgen)
        }
    }

    /// Puts back the generator taken by `take_randgen` (called on every exit of the caller).
    fn put_randgen(
        &mut self,
        ms: &mut HighsMipSolver,
        worker: usize,
        locked: bool,
        randgen: HighsRandom,
    ) {
        if locked {
            ms.mipdata_.workers[worker].randgen = randgen;
        } else {
            self.randgen = randgen;
        }
    }
}

/// C++ `HighsDomain localdom = worker.getGlobalDomain();` (the copy constructor).
fn local_copy_of_global_domain(ms: &mut HighsMipSolver, worker: usize) -> HighsDomain {
    let gd = ms.mipdata_.workers[worker].get_global_domain();
    detach!(ms.mipdata_.domains[gd], |g| HighsDomain::copy_from(ms, g))
}

/// C++ `localdom.conflictAnalysis(worker.getConflictPool(), worker.getGlobalDomain(), worker.getPseudocost())`.
fn conflict_analysis_for_worker(
    ms: &mut HighsMipSolver,
    worker: usize,
    localdom: &mut HighsDomain,
) {
    let _src = crate::work::conflict_src(crate::work::ConflictSrc::Heur);
    let cp = ms.mipdata_.workers[worker].get_conflict_pool();
    let gd = ms.mipdata_.workers[worker].get_global_domain();
    let pk = ms.mipdata_.workers[worker].get_pseudocost();
    detach!(ms.mipdata_.pseudocosts[pk], |pc| localdom
        .conflict_analysis(ms, cp, gd, pc));
}

#[cfg(test)]
mod time_limit_tests {
    use crate::lp_data::hconst::{HighsVarType, ObjSense};
    use crate::lp_data::highs_run::Highs;
    use crate::prelude::*;
    use std::time::Instant;

    /// Disjoint 5-cycles, maximise the node count: the root LP is 0.5 everywhere (all columns
    /// fractional) and the shifting and ZI-round passes are quadratic in the model size. With a
    /// short time limit the solve must stop near it, not run those passes to completion. This
    /// size is a smoke bound: the model does not separate the fixed from the unfixed loops (the
    /// overrun reproduces on gfd-schedulen180f7d50m30k18).
    ///
    /// Timing test calibrated for an optimised build: a debug build needs about 10 times longer
    /// per pass and overshoots the smoke bound (33.6 s seen), so it is ignored there. Run it with
    /// `cargo test --release -p rusolve-ekk --lib shifting_and_zi` (release: 3 s limit, 4.5 to
    /// 7.7 s wall with 3000 and 6000 cycles, so no limit check is missing).
    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "timing test: a debug build overshoots the smoke bound; run with --release"
    )]
    fn shifting_and_zi_round_honour_the_time_limit() {
        const CYCLES: usize = 3000;
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_int("threads", 1);
        h.set_option_value_double("time_limit", 3.0);
        h.set_option_value_bool("mip_heuristic_run_shifting", true);
        h.set_option_value_bool("mip_heuristic_run_zi_round", true);
        for j in 0..5 * CYCLES {
            h.add_vars(1, &[0.0], &[1.0]);
            h.change_col_cost(j as HighsInt, 1.0);
            h.change_col_integrality(j as HighsInt, HighsVarType::kInteger);
        }
        for cyc in 0..CYCLES {
            for kk in 0..5 {
                let (a, b) = (
                    (5 * cyc + kk) as HighsInt,
                    (5 * cyc + (kk + 1) % 5) as HighsInt,
                );
                h.add_row(-kHighsInf, 1.0, 2, &[a, b], &[1.0, 1.0]);
            }
        }
        h.change_objective_sense(ObjSense::kMaximize);
        let t = Instant::now();
        let _ = h.run();
        let secs = t.elapsed().as_secs_f64();
        assert!(secs < 30.0, "solve ran {secs:.1} s against a 3 s limit");
    }
}
