//! rusolve `mip_cdcl_dive`: an LP-free conflict-driven depth-first search with restarts,
//! run as a primal heuristic on the propagated global domain.
//!
//! Source: the technique is the CDCL / nogood-learning search of
//! constraint and SAT solvers applied to MIP feasibility: Achterberg, "Conflict analysis in
//! mixed integer programming", Discrete Optimization 4 (2007) 4-20; Achterberg, "SCIP: solving
//! constraint integer programs", Math. Prog. Comp. 1 (2009), ch. 6 (constraint integer
//! programming searches without an LP); Moskewicz et al., "Chaff", DAC 2001 (activity-based
//! decision order, VSIDS); Luby, Sinclair, Zuckerman, "Optimal speedup of Las Vegas
//! algorithms", IPL 1993 (restart schedule). The learned nogoods are the port's own
//! `HighsDomain::conflictAnalysis` conflicts (plus resolution when on), which the
//! domain's conflict-pool propagation then enforces.
//!
//! No LP is solved: a node costs one propagation, so a run visits orders of magnitude more
//! nodes than the LP tree. Pure-integer models only (a leaf is a complete assignment).

use crate::mip::domain::HighsDomain_Reason;
use crate::prelude::*;

/// Luby sequence value (1-based): 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8 ...
fn luby(i: u64) -> u64 {
    let mut k = 1u32;
    while (1u64 << k) - 1 < i {
        k += 1;
    }
    let mut i = i;
    loop {
        if i == (1u64 << k) - 1 {
            return 1u64 << (k - 1);
        }
        i -= (1u64 << (k - 1)) - 1;
        k = 1;
        while (1u64 << k) - 1 < i {
            k += 1;
        }
    }
}

impl HighsMipSolver {
    /// Run the conflict-driven DFS with at most `max_conflicts` conflicts. Returns true when
    /// it submitted a feasible solution.
    pub fn cdcl_dive(&mut self, max_conflicts: u64) -> bool {
        let n = self.num_col() as usize;
        if self.submip || self.mipdata_.integral_cols.len() != n || n == 0 {
            return false;
        }
        let gd = self.mipdata_.workers[0].get_global_domain();
        if self.options_mip_.mip_cdcl_uip {
            return self.cdcl_uip_dive(gd, 25.0);
        }
        let cp = self.mipdata_.workers[0].get_conflict_pool();
        let pc = self.mipdata_.workers[0].get_pseudocost();
        let mut localdom: HighsDomain =
            detach!(self.mipdata_.domains[gd], |g| HighsDomain::copy_from(
                self, g
            ));
        let cols: Vec<HighsInt> = self.mipdata_.integral_cols.clone();
        let mut activity = vec![0.0f64; n];
        // phase saving: start with the cost-minimising bound
        let mut phase_up: Vec<bool> = (0..n as HighsInt)
            .map(|j| self.col_cost_2(j) < 0.0)
            .collect();
        let mut bump = 1.0f64;
        let mut flipped: Vec<bool> = Vec::new();
        let mut conflicts = 0u64;
        let mut nodes = 0u64;
        let mut restarts = 1u64;
        let mut next_restart = 100 * luby(restarts);
        let mut found = false;
        let t0 = std::time::Instant::now();
        let tmax: f64 = 10.0; // seconds per run

        'search: loop {
            localdom.propagate(self);
            if localdom.infeasible() {
                conflicts += 1;
                // decision-level activity bump (the branching columns on the trail)
                for &p in localdom.branch_pos_.iter() {
                    let j = localdom.domchgstack_[p as usize].column as usize;
                    activity[j] += bump;
                }
                bump *= 1.0 / 0.95;
                if bump > 1e100 {
                    activity.iter_mut().for_each(|a| *a *= 1e-100);
                    bump *= 1e-100;
                }
                let _src = crate::work::conflict_src(crate::work::ConflictSrc::Heur);
                detach!(self.mipdata_.pseudocosts[pc], |p| localdom
                    .conflict_analysis(self, cp, gd, p));
                if conflicts >= max_conflicts {
                    break 'search;
                }
                if conflicts >= next_restart {
                    localdom.backtrack_to_global(self);
                    flipped.clear();
                    restarts += 1;
                    next_restart = conflicts + 100 * luby(restarts);
                    continue;
                }
                // chronological backtrack: flip the deepest unflipped decision
                loop {
                    if flipped.is_empty() {
                        break 'search; // search space exhausted under the cutoff
                    }
                    let undone = localdom.backtrack(self);
                    let was_flipped = flipped.pop().unwrap();
                    if undone.column < 0 {
                        break 'search;
                    }
                    if !was_flipped {
                        let sib = if undone.boundtype == HighsBoundType::kUpper {
                            HighsDomainChange {
                                boundval: undone.boundval + 1.0,
                                column: undone.column,
                                boundtype: HighsBoundType::kLower,
                            }
                        } else {
                            HighsDomainChange {
                                boundval: undone.boundval - 1.0,
                                column: undone.column,
                                boundtype: HighsBoundType::kUpper,
                            }
                        };
                        flipped.push(true);
                        localdom.change_bound(self, sib, HighsDomain_Reason::branching());
                        break;
                    }
                }
                continue;
            }
            nodes += 1;
            if nodes & 63 == 0 && (self.terminate() || t0.elapsed().as_secs_f64() > tmax) {
                break;
            }
            // pick the unfixed column of maximal activity
            let mut best: Option<usize> = None;
            let mut best_act = -1.0f64;
            for &j in cols.iter() {
                let cu = j as usize;
                if localdom.col_lower_[cu] < localdom.col_upper_[cu] && activity[cu] > best_act {
                    best_act = activity[cu];
                    best = Some(cu);
                }
            }
            let Some(col) = best else {
                // all columns fixed: a complete assignment
                let sol: Vec<f64> = localdom.col_lower_.clone();
                found = self.try_solution(&sol, kSolutionSourceHeuristic);
                break;
            };
            let (lb, ub) = (localdom.col_lower_[col], localdom.col_upper_[col]);
            let (chg, up) = if phase_up[col] {
                let v = if ub.is_finite() { ub } else { lb.max(0.0) };
                (
                    HighsDomainChange {
                        boundval: v,
                        column: col as HighsInt,
                        boundtype: HighsBoundType::kLower,
                    },
                    true,
                )
            } else {
                let v = if lb.is_finite() { lb } else { ub.min(0.0) };
                (
                    HighsDomainChange {
                        boundval: v,
                        column: col as HighsInt,
                        boundtype: HighsBoundType::kUpper,
                    },
                    false,
                )
            };
            phase_up[col] = up;
            flipped.push(false);
            localdom.change_bound(self, chg, HighsDomain_Reason::branching());
        }
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "CDCL dive: {} nodes, {} conflicts, {} restarts, {}\n",
            nodes,
            conflicts,
            restarts,
            if found {
                "solution found"
            } else {
                "no solution"
            }
        );
        localdom.release(self);
        found
    }
}

impl HighsMipSolver {
    /// `mip_zero_obj_cdcl`: with an all-zero objective every feasible point is optimal, so run
    /// the CDCL engine as a root feasibility search (10 % of the time limit, 5-60 s). Returns
    /// true when it submitted a solution.
    pub fn cdcl_zero_obj(&mut self) -> bool {
        if !self.options_mip_.mip_zero_obj_cdcl
            || self.submip
            || self.mipdata_.num_restarts != 0
            || self.mipdata_.upper_limit != kHighsInf
            || self.num_col() == 0
            || self.col_cost().iter().any(|&cc| cc != 0.0)
        {
            return false;
        }
        let gd = self.mipdata_.workers[0].get_global_domain();
        let budget = (0.1 * self.options_mip_.time_limit).clamp(5.0, 60.0);
        self.cdcl_uip_dive(gd, budget)
    }

    /// `mip_cdcl_uip`: the first-UIP / two-watched-literal engine of `cdcl_uip.rs` on the
    /// global domain and the model rows (binaries are decided, other finite-bound columns are
    /// propagated).
    fn cdcl_uip_dive(&mut self, gd: usize, budget: f64) -> bool {
        use crate::mip::cdcl_uip::{Outcome, UipEngine};
        let n = self.num_col() as usize;
        let (lb, ub): (Vec<f64>, Vec<f64>) = detach!(self.mipdata_.domains[gd], |g| (
            g.col_lower_.clone(),
            g.col_upper_.clone()
        ));
        let mut rows: Vec<(Vec<(u32, f64)>, f64)> = Vec::new();
        let (rl, ru) = (self.row_lower().to_vec(), self.row_upper().to_vec());
        for r in 0..self.num_row() as usize {
            let (s, e) = (
                self.mipdata_.a_rstart_[r] as usize,
                self.mipdata_.a_rstart_[r + 1] as usize,
            );
            let ent: Vec<(u32, f64)> = (s..e)
                .map(|k| {
                    (
                        self.mipdata_.a_rindex_[k] as u32,
                        self.mipdata_.a_rvalue_[k],
                    )
                })
                .collect();
            if ru[r] < kHighsInf {
                rows.push((ent.clone(), ru[r]));
            }
            if rl[r] > -kHighsInf {
                rows.push((ent.iter().map(|&(j, a)| (j, -a)).collect(), -rl[r]));
            }
        }
        let phase: Vec<bool> = (0..n as HighsInt)
            .map(|j| self.col_cost_2(j) < 0.0)
            .collect();
        let is_int: Vec<bool> = self
            .variable_type()
            .iter()
            .map(|t| *t != HighsVarType::kContinuous)
            .collect();
        let Some(mut eng) = UipEngine::new(n, &rows, &lb, &ub, &is_int, &phase) else {
            return false;
        };
        let t0 = std::time::Instant::now();
        let mut stop = || self.terminate() || t0.elapsed().as_secs_f64() > budget;
        let out = eng.solve(u64::MAX, &mut stop);
        let mut found = false;
        let what = match out {
            Outcome::Solution(sol) => {
                found = self.try_solution(&sol, kSolutionSourceHeuristic);
                if found {
                    "solution found"
                } else {
                    "solution rejected"
                }
            }
            Outcome::Exhausted => "search space exhausted",
            Outcome::Stopped => "no solution",
        };
        let secs = t0.elapsed().as_secs_f64();
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "CDCL UIP: {} decisions, {} conflicts ({:.0}/s), {} restarts, learnt lits {} -> {}, {}\n",
            eng.decisions,
            eng.conflicts,
            eng.conflicts as f64 / secs.max(1e-9),
            eng.restarts,
            eng.learnt_lits,
            eng.minimized_lits,
            what
        );
        crate::highs_log_user!(
            &self.options_mip_.log_options,
            HighsLogType::kInfo,
            "CDCL UIP stats: {:.2}s total, propagate {:.2}s, analyze {:.2}s, rows scanned {} ({} entries), watch visits {}\n",
            secs,
            eng.st_prop_ns as f64 * 1e-9,
            eng.st_analyze_ns as f64 * 1e-9,
            eng.st_rows,
            eng.st_row_entries,
            eng.st_watch
        );
        found
    }
}

#[cfg(test)]
mod tests {
    use super::luby;
    #[test]
    fn luby_prefix() {
        let v: Vec<u64> = (1..=15).map(luby).collect();
        assert_eq!(v, vec![1, 1, 2, 1, 1, 2, 4, 1, 1, 2, 1, 1, 2, 4, 8]);
    }
}
