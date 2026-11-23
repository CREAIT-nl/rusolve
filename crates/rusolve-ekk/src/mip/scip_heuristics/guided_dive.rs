//! Translated from SCIP (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! `mip_guided_dive`: SCIP's guided diving, an LP dive toward the incumbent
//!.
//!
//! * SCIP `src/scip/heur_guideddiving.c`: `divesetGetScoreGuideddiving` (round each fractional
//!   variable toward its incumbent value; score = -fractionality of that direction, with the
//!   0.01 penalty x0.1 and the non-binary x0.1 factors) and the defaults `maxlpiterquot` 0.05,
//!   `maxlpiterofs` 1000, backtrack = one level.
//! * SCIP `src/scip/heur_dive.c` / `SCIPperformGenericDivingAlgorithm`: one bound change per LP
//!   resolve, propagation, one backtrack to the other child on infeasibility.
//!
//! Divergences: every candidate is scored as "cannot be rounded" (the port keeps no variable
//! locks here), so the score is the directional fractionality alone; the dive runs on the
//! port's `HighsSearch` heuristic search and LP copy as RINS does; the LP-iteration budget is
//! `0.05 * node LP iterations + 1000` summed over all calls.

use super::*;
use crate::mip::mip_solver::HighsMipSolver;

impl HighsPrimalHeuristics {
    /// Guided dive from the current node LP solution toward `ms.mipdata_.incumbent`.
    pub fn guided_dive(&mut self, ms: &mut HighsMipSolver, worker: usize) {
        let _lp_src_guard = crate::work::lp_heur(crate::work::LP_H_GUIDED_DIVE);
        let gd = ms.mipdata_.workers[worker].get_global_domain();
        if ms.submip || ms.mipdata_.incumbent.is_empty() || ms.mipdata_.domains[gd].infeasible() {
            return;
        }
        let node_lp = ms.mipdata_.total_lp_iterations
            - ms.mipdata_.heuristic_lp_iterations
            - ms.mipdata_.sb_lp_iterations;
        let budget = 1000 + node_lp / 20;
        if self.guided_lp_iterations > budget {
            return;
        }

        let pscost_k = ms.mipdata_.workers[worker].get_pseudocost();
        let mut pscost: HighsPseudocost = detach!(ms.mipdata_.pseudocosts[pscost_k], |p| p.clone());
        let mut heur = HighsSearch::new(ms, gd);
        heur.set_heuristic(true);
        let lp_k = ms.mipdata_.workers[worker].get_lp_relaxation();
        let mut heurlp = HighsLpRelaxation::copy_from(ms, &ms.mipdata_.lps[lp_k]);
        heurlp.set_mip_worker(worker);
        heurlp.set_profiling();
        let worker_upper_limit = ms.mipdata_.workers[worker].upper_limit;
        heurlp.set_objective_limit(ms, worker_upper_limit);
        heurlp.set_adjust_symmetric_branching_col(false);
        let num_col_minus_one = ms.num_col() - 1;
        let _ = heurlp.get_lp_solver().change_cols_bounds_interval(
            0,
            num_col_minus_one,
            &heur.localdom.col_lower_,
            &heur.localdom.col_upper_,
        );
        heur.localdom.clear_changed_cols();
        heur.create_new_node();

        let sols_before = ms.mipdata_.num_improving_sols;
        let mut backtracked = false;
        let iters_start = heur.get_local_lp_iterations();

        loop {
            {
                let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                heur.evaluate_node(&mut ctx);
            }
            if ms.mipdata_.num_improving_sols != sols_before {
                break;
            }
            if self.guided_lp_iterations + heur.get_local_lp_iterations() - iters_start > budget {
                break;
            }
            if heur.current_node_pruned() {
                if backtracked || ms.mipdata_.domains[gd].infeasible() {
                    break;
                }
                backtracked = true;
                let ok = {
                    let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                    heur.backtrack(&mut ctx, true)
                };
                if !ok || !heur.has_node() {
                    break;
                }
                continue;
            }
            // pick the candidate: best score = smallest fractionality toward the incumbent
            let mut best: Option<(f64, HighsInt, f64, bool)> = None;
            for &(col, val) in heurlp.get_fractional_integers() {
                let ci = col as usize;
                let bestval = ms.mipdata_.incumbent[ci];
                let roundup = val < bestval;
                let mut frac = val - val.floor();
                if roundup {
                    frac = 1.0 - frac;
                }
                if frac < 0.01 {
                    frac *= 0.1;
                }
                if !heur.localdom.is_binary(ms, col) {
                    frac *= 0.1;
                }
                let score = -frac;
                if best.is_none_or(|b| score > b.0) {
                    best = Some((score, col, val, roundup));
                }
            }
            let Some((_, col, val, roundup)) = best else {
                break;
            };
            {
                let mut ctx = SearchCtx::new(&mut *ms, worker, &mut heurlp, &mut pscost);
                if roundup {
                    heur.branch_upwards(&mut ctx, col, val.ceil(), val);
                } else {
                    heur.branch_downwards(&mut ctx, col, val.floor(), val);
                }
            }
            heur.localdom.propagate(ms);
            if heur.localdom.infeasible() {
                // the node evaluation prunes it and the pruned branch backtracks once
                continue;
            }
            heurlp.flush_domain(ms, &mut heur.localdom, false);
        }
        let used = heur.get_local_lp_iterations();
        self.guided_lp_iterations += used - iters_start;
        *ms.mipdata_.workers[worker].get_heur_lp_iterations() += used;
        heur.release(ms);
    }
}
