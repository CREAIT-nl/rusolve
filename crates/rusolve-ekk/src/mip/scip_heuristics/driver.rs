//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.

use super::*;
use crate::mip::mip_solver::HighsMipSolver;

/// SCIP 11.0.0 (32dd425): scip/heur.c::SCIPheurExec (lines 1264-1340)
pub(crate) fn s6_run_tree_heuristics(ms: &mut HighsMipSolver) {
    if ms.options_mip_.mip_heuristic_run_oneopt {
        s6_run_oneopt(ms);
    }
    // mip_alns runs crossover as a bandit arm, not twice.
    if ms.options_mip_.mip_heuristic_run_crossover && !ms.submip && !ms.options_mip_.mip_alns {
        s6_run_crossover(ms);
    }
}

/// SCIP 11.0.0 (32dd425): scip/heur_oneopt.c::heurInitsolOneopt (lines 422-437)
pub(crate) fn s6_run_root_oneopt(ms: &mut HighsMipSolver) {
    if ms.options_mip_.mip_heuristic_run_oneopt {
        s6_run_oneopt(ms);
    }
}

fn s6_run_oneopt(ms: &mut HighsMipSolver) {
    if let Some(s6) = ms.mipdata_.s6.as_mut() {
        s6.cur_origin = S6Origin::Oneopt;
    }
    let _ = crate::mip::scip_heuristics::oneopt::s6_oneopt_exec(ms);
    if let Some(s6) = ms.mipdata_.s6.as_mut() {
        s6.cur_origin = S6Origin::Highs(crate::mip::mip_solver_data::kSolutionSourceNone);
    }
}

fn s6_run_crossover(ms: &mut HighsMipSolver) {
    if let Some(s6) = ms.mipdata_.s6.as_mut() {
        s6.cur_origin = S6Origin::Crossover;
    }
    let _ = crate::mip::scip_heuristics::crossover::s6_crossover_exec(ms);
    if let Some(s6) = ms.mipdata_.s6.as_mut() {
        s6.cur_origin = S6Origin::Highs(crate::mip::mip_solver_data::kSolutionSourceNone);
    }
}
