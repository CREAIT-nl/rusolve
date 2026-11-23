//! Node-LP engine selection: [`MipEngine`] is what `bnb::engine` hands the
//! tree search. The default is [`RusolveEngine`] (unchanged behaviour);
//! `RUSOLVE_MIP_EKK=1` (any value but `0`) selects
//! [`crate::engine::ekk::EkkEngine`], the faithful HiGHS simplex driven
//! through the ported `Highs` API (crates/rusolve-ekk).
//!
//! An enum rather than `Box<dyn LpEngine>` keeps `Search<E>` monomorphic
//! and the default path's dispatch a single predictable branch.

use crate::engine::ekk::EkkEngine;
use crate::engine::{
    CutRow, LpDuals, LpEngine, LpError, LpOutcome, LpStats, RusolveEngine, TableauAgg, WarmStart,
};

/// Whether `RUSOLVE_MIP_EKK` asks for the ekk engine.
pub fn ekk_requested() -> bool {
    std::env::var("RUSOLVE_MIP_EKK").is_ok_and(|v| v != "0")
}

// One value per search, moved into `Search` once: the variant size gap
// costs nothing, and boxing `RusolveEngine` would add a pointer chase to
// every default-path engine call.
#[allow(clippy::large_enum_variant)]
pub enum MipEngine {
    Rusolve(RusolveEngine),
    Ekk(Box<EkkEngine>),
}

macro_rules! delegate {
    ($self:ident, $e:ident => $call:expr) => {
        match $self {
            MipEngine::Rusolve($e) => $call,
            MipEngine::Ekk($e) => $call,
        }
    };
}

impl LpEngine for MipEngine {
    fn solve(&mut self, col_lower: &[f64], col_upper: &[f64]) -> Result<LpOutcome, LpError> {
        delegate!(self, e => e.solve(col_lower, col_upper))
    }

    fn solve_from(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
    ) -> Result<LpOutcome, LpError> {
        delegate!(self, e => e.solve_from(col_lower, col_upper, start, cutoff))
    }

    fn probe(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
        iter_limit: u64,
    ) -> Result<LpOutcome, LpError> {
        delegate!(self, e => e.probe(col_lower, col_upper, start, cutoff, iter_limit))
    }

    fn set_iteration_limit(&mut self, limit: Option<u64>) {
        delegate!(self, e => e.set_iteration_limit(limit))
    }

    fn set_root_lp(&mut self, dual: bool, retry: bool) {
        delegate!(self, e => e.set_root_lp(dual, retry))
    }

    fn set_deadline(&mut self, deadline: Option<std::time::Instant>) {
        delegate!(self, e => e.set_deadline(deadline))
    }

    fn add_cuts(&mut self, cuts: &[CutRow]) -> bool {
        delegate!(self, e => e.add_cuts(cuts))
    }

    fn retain_cuts(&mut self, keep: &[bool]) {
        delegate!(self, e => e.retain_cuts(keep))
    }

    fn num_lp_rows(&self) -> usize {
        delegate!(self, e => e.num_lp_rows())
    }

    fn stats(&self) -> LpStats {
        delegate!(self, e => e.stats())
    }

    fn duals(&mut self, lp: &LpOutcome) -> Option<LpDuals> {
        delegate!(self, e => e.duals(lp))
    }

    fn cut_rows_basic(&mut self, lp: &LpOutcome) -> Option<Vec<bool>> {
        delegate!(self, e => e.cut_rows_basic(lp))
    }

    fn infeasibility_proof(&self) -> Option<(Vec<f64>, f64)> {
        delegate!(self, e => e.infeasibility_proof())
    }

    fn cutoff_proof(&self) -> Option<(Vec<f64>, f64)> {
        delegate!(self, e => e.cutoff_proof())
    }

    fn tableau_weights(
        &mut self,
        lp: &LpOutcome,
        lower: &[f64],
        upper: &[f64],
        cols: &[usize],
    ) -> Option<Vec<TableauAgg>> {
        delegate!(self, e => e.tableau_weights(lp, lower, upper, cols))
    }
}
