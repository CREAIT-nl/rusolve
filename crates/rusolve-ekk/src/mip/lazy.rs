//! rusolve extension: lazy constraints in the Rust callback API. Design: docs/design/hooks.md,
//! section "Lazy constraints".
//!
//! This file is rusolve's own (not a HiGHS translation); HiGHS 1.16 defines
//! `kCallbackMipDefineLazyConstraints` but never fires it.
//!
//! A lazy constraint is a row the model needs but that is too many to list (subtour elimination
//! is the standard case: Dantzig, Fulkerson, Johnson, Operations Research 2 (1954) 393-410;
//! lazy constraints as a solver callback: the SCIP constraint-handler `enforce` callback,
//! `SCIP_CONSHDLR`). The user oracle receives every candidate optimal solution of
//! the current model; it returns the rows that solution violates, and an empty list accepts it.
//!
//! Mechanism: an outer loop around `Highs::run` (solve, ask the oracle, add the violated rows to
//! the model, re-solve). Because it works on the model, not inside the tree, it is independent of
//! presolve (rows and indices are in the original space), the search is untouched when it is not
//! used, and the final solution is exactly optimal for the model plus every row the oracle added.
//! The price against in-tree lazy constraints is that each round re-solves the MIP.
//!
//! Gates: a returned row is dropped when it is malformed (length mismatch, empty, non-finite,
//! index out of range, duplicate column) or when the candidate solution does not violate it
//! (it could not make progress). Counts are in [`LazyReport`].

use crate::lp_data::highs_run::Highs;
use crate::mip::hooks::UserCut;
use crate::prelude::*;
use std::collections::HashSet;

/// The lazy-constraint oracle. `solution` is the candidate (original column space).
pub trait LazyConstraints {
    /// Rows `a x <= rhs` that every feasible solution satisfies; return the ones `solution`
    /// violates (rows it satisfies are ignored). An empty list accepts the solution.
    fn lazy_constraints(&mut self, solution: &[f64]) -> Vec<UserCut>;
}

impl<F: FnMut(&[f64]) -> Vec<UserCut>> LazyConstraints for F {
    fn lazy_constraints(&mut self, solution: &[f64]) -> Vec<UserCut> {
        self(solution)
    }
}

/// What the loop did.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LazyReport {
    /// MIP solves performed.
    pub rounds: usize,
    /// Rows added to the model.
    pub rows_added: usize,
    /// Rows dropped because malformed.
    pub rows_rejected_malformed: usize,
    /// Rows the candidate solution already satisfied (ignored).
    pub rows_not_violated: usize,
    /// True when the oracle accepted the last solution; false at the round limit or when a
    /// solve did not end optimal.
    pub accepted: bool,
}

impl Highs {
    /// Solve the loaded MIP with lazy constraints: `run`, offer the solution to `lazy`, add the
    /// violated rows, repeat until the oracle accepts, a solve is not optimal, or `max_rounds`
    /// solves were made. Returns the last `run` status and the report. On acceptance the model
    /// status is optimal and `get_solution` is the accepted solution.
    pub fn run_with_lazy_constraints(
        &mut self,
        lazy: &mut dyn LazyConstraints,
        max_rounds: usize,
    ) -> (HighsStatus, LazyReport) {
        let mut rep = LazyReport::default();
        let mut status = HighsStatus::kOk;
        while rep.rounds < max_rounds {
            status = self.run();
            rep.rounds += 1;
            if self.get_model_status() != HighsModelStatus::kOptimal {
                return (status, rep);
            }
            let sol = self.get_solution().col_value.clone();
            let feastol = self.options_.mip_feasibility_tolerance;
            let num_col = sol.len() as HighsInt;
            let mut added = 0usize;
            let mut malformed_round = 0usize;
            for cut in lazy.lazy_constraints(&sol) {
                let mut seen = HashSet::new();
                let malformed = cut.indices.len() != cut.values.len()
                    || cut.indices.is_empty()
                    || !cut.rhs.is_finite()
                    || cut.values.iter().any(|v| !v.is_finite())
                    || cut
                        .indices
                        .iter()
                        .any(|&col| col < 0 || col >= num_col || !seen.insert(col));
                if malformed {
                    rep.rows_rejected_malformed += 1;
                    malformed_round += 1;
                    continue;
                }
                let act: f64 = cut
                    .indices
                    .iter()
                    .zip(&cut.values)
                    .map(|(&col, &v)| v * sol[col as usize])
                    .sum();
                if act <= cut.rhs + feastol {
                    rep.rows_not_violated += 1;
                    continue;
                }
                self.add_row(
                    -kHighsInf,
                    cut.rhs,
                    cut.indices.len() as HighsInt,
                    &cut.indices,
                    &cut.values,
                );
                added += 1;
            }
            rep.rows_added += added;
            if added == 0 {
                // Only an empty answer, or one whose rows are all well-formed and already
                // satisfied, accepts. A non-empty answer made only of malformed rows is an oracle
                // bug: re-solving would return the same solution, so stop with accepted = false.
                rep.accepted = malformed_round == 0;
                return (status, rep);
            }
        }
        (status, rep)
    }
}
