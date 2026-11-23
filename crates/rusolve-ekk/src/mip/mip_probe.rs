//! A read-only snapshot of the counters a finished
//! `HighsMipSolver` keeps, for `rusolve probe`.
//!
//! Not a HiGHS translation: every field is a value the port already records (the ones
//! `HighsMipSolverData::printDisplayLine` and `HighsMipSolver::solvingReport` print, plus the
//! per-separator `numCalls` / `numCutsFound` of HighsSeparator.h:34-58). The snapshot is taken
//! after `HighsMipSolver::run` returned and only when `Highs::capture_mip_probe_` is set, so it
//! adds no work to the search and nothing at all to a default solve.

use crate::mip::HighsSeparatorBox;
use crate::mip::mip_solver::HighsMipSolver;
use crate::prelude::*;

/// Counters of one separator (HighsSeparator.h numCalls / numCutsFound).
#[derive(Clone, Debug, Default)]
pub struct MipProbeSeparator {
    pub name: &'static str,
    pub num_calls: i64,
    pub num_cuts_found: i64,
}

/// What a finished MIP solve left in `HighsMipSolverData`. Objective values are in the solver's
/// internal minimisation sense *without* offset; `sense` and `offset` convert them.
#[derive(Clone, Debug, Default)]
pub struct MipProbeSnapshot {
    /// Dimensions of the model the MIP solver ended on (presolved, after the last restart).
    pub model_num_row: i64,
    pub model_num_col: i64,
    pub model_num_nz: i64,
    pub model_num_integer: i64,
    pub sense: f64,
    pub offset: f64,
    pub num_restarts: i64,
    pub num_nodes: i64,
    pub num_leaves: i64,
    pub total_lp_iterations: i64,
    pub heuristic_lp_iterations: i64,
    pub sepa_lp_iterations: i64,
    pub sb_lp_iterations: i64,
    /// `firstrootlpiters` (HighsMipSolverData.h): LP iterations of the first root LP solve.
    pub first_root_lp_iterations: i64,
    pub num_improving_sols: i64,
    pub lower_bound: f64,
    pub upper_bound: f64,
    /// Cuts alive in the global cut pool (printDisplayLine "Cuts" column).
    pub cutpool_cuts: i64,
    /// Open nodes left in the queue and the deepest of them (a lower bound on the depth reached).
    pub open_nodes: i64,
    pub open_nodes_max_depth: Option<i64>,
    pub separators: Vec<MipProbeSeparator>,
}

impl HighsMipSolver {
    /// Read-only copy of the statistics counters; see the module doc.
    pub fn probe_snapshot(&self) -> MipProbeSnapshot {
        let model = crate::mip_model!(self);
        let d = &self.mipdata_;
        let num_col = model.num_col_ as usize;
        let num_nz = model.a_matrix_.start_.get(num_col).copied().unwrap_or(0);
        let mut separators = Vec::new();
        for w in d.workers.iter() {
            if let Some(sepa) = w.sepa_ptr_.as_ref() {
                for s in sepa.separators.iter() {
                    let name = match s {
                        HighsSeparatorBox::Tableau(_) => "tableau",
                        HighsSeparatorBox::Path(_) => "path",
                        HighsSeparatorBox::Modk(_) => "modk",
                        HighsSeparatorBox::MachineSched(_) => "machine_sched",
                        HighsSeparatorBox::KnapsackCover(_) => "knapsack_cover",
                        HighsSeparatorBox::TiPrecedence(_) => "ti_precedence",
                        HighsSeparatorBox::Netdesign(_) => "netdesign",
                        HighsSeparatorBox::Dff(_) => "dff",
                    };
                    let b = s.base();
                    separators.push(MipProbeSeparator {
                        name,
                        num_calls: i64::from(b.num_calls),
                        num_cuts_found: i64::from(b.num_cuts_found),
                    });
                }
                break;
            }
        }
        let open_nodes_max_depth = d.nodequeue.nodes.iter().map(|n| i64::from(n.depth)).max();
        MipProbeSnapshot {
            model_num_row: i64::from(model.num_row_),
            model_num_col: i64::from(model.num_col_),
            model_num_nz: i64::from(num_nz),
            model_num_integer: i64::from(d.numintegercols),
            sense: f64::from(self.orig_model_.sense_ as i32),
            offset: model.offset_,
            num_restarts: i64::from(d.num_restarts),
            num_nodes: d.num_nodes,
            num_leaves: d.num_leaves,
            total_lp_iterations: d.total_lp_iterations,
            heuristic_lp_iterations: d.heuristic_lp_iterations,
            sepa_lp_iterations: d.sepa_lp_iterations,
            sb_lp_iterations: d.sb_lp_iterations,
            first_root_lp_iterations: d.firstrootlpiters,
            num_improving_sols: i64::from(d.num_improving_sols),
            lower_bound: d.lower_bound,
            upper_bound: d.upper_bound,
            cutpool_cuts: i64::from(d.cutpools[0].get_num_cuts()),
            open_nodes: d.nodequeue.num_active_nodes(),
            open_nodes_max_depth,
            separators,
        }
    }
}
