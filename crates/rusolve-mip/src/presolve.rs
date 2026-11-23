//! MIP presolve + primal postsolve.
//!
//! Runs once on the [`MipProblem`] before branch-and-bound. Reductions, each
//! valid for the *integer* feasible set (not only the LP relaxation):
//!
//! * **Column bounds**: integer bounds rounded inward, fixed columns removed
//!   (their row contribution moves into the row sides, their cost into the
//!   objective offset), empty columns fixed at their cost-optimal bound.
//! * **Dual fixing**: a column whose cost does not oppose moving it towards a
//!   finite bound, and that no row locks in that direction, is fixed there
//!   (Achterberg et al. 2020 §4.4; HiGHS `HPresolve::dualFixing`,
//!   HPresolve.cpp:4670).
//! * **Rows**: empty rows checked and removed, singleton rows turned into
//!   column bounds (HiGHS `singletonRow`, HPresolve.cpp:3205), infeasibility
//!   and redundancy from activity bounds (HiGHS `rowPresolve`
//!   `checkRowInfeasible`/`checkRowRedundant`, HPresolve.cpp:3446–3466); a
//!   side implied by the activity bounds is dropped, which turns ranged rows
//!   into inequalities that coefficient strengthening can work on.
//! * **Row gcd / integral scaling** on rows whose columns are all integer:
//!   scale the row so its coefficients are coprime integers, then round the
//!   sides inward (`⌊hi·s + ε⌋`, `⌈lo·s − ε⌉`). An equation whose scaled
//!   rhs is fractional — e.g. all-even coefficients with an odd rhs — is
//!   infeasible. HiGHS `rowPresolve` inequality branch
//!   (HPresolve.cpp:3866–4054: `checkScaleRow`, `scaleRowIntVals` and its
//!   `maxVal > 1000 && intScale > 100` guard) and the common-denominator
//!   idea of `HighsIntegers::integralScale` (highs/util/HighsIntegers.h:127);
//!   reimplemented with a continued-fraction denominator per coefficient and
//!   an exact lcm/gcd, and only for rows that scale to integers *exactly*
//!   (no HiGHS-style `coefDelta` relaxation of near-integral coefficients).
//! * **Coefficient strengthening** of one-sided inequalities with integer
//!   columns: for `aᵀx ≤ b` with max activity `M` and `d = M − b > 0`, an
//!   integer column with `a_j > d` gets coefficient `d` and `b += (d − a_j)·u_j`
//!   (mirrored for `a_j < −d` with `l_j`). Same integer points, tighter LP.
//!   Achterberg et al. 2020 §3.2; HiGHS `strengthenCoefs` lambda in
//!   `rowPresolve` (HPresolve.cpp:4217–4283).
//! * **Doubleton equations** `a_e x_e + a_t x_t = b`: substitute
//!   `x_e = (b − a_t x_t)/a_e` everywhere, moving `x_e`'s bounds onto `x_t`.
//!   Integer-safe choice of the substituted column (continuous first; two
//!   integers only when `a_t/a_e` is integral, then `b/a_e` must be too or the
//!   row is infeasible). HiGHS `doubletonEq` (HPresolve.cpp:3036–3203,
//!   `colAtPos1Better` substitution order and integrality conditions).
//!
//! * **Activity bound tightening** (Savelsbergh 1994; Achterberg et al. 2020
//!   §3.1), which feeds coefficient strengthening of big-M / VUB rows.
//! * **Pair reductions**, once per outer round: parallel rows merged into
//!   one (ranged) row and identical columns merged into their sum (HiGHS
//!   `detectParallelRowsAndCols`, HPresolve.cpp:7267), dominated binary
//!   columns fixed (Achterberg et al. 2020 §6.4; Gamrath et al. 2015; HiGHS
//!   `dominatedColumns`, HPresolve.cpp:1045), set-packing conflicts fixed
//!   (Balas & Padberg 1976).
//!
//! Postsolve is primal only (B&B needs no duals): a stack of fixed values and
//! substitutions replayed in reverse, after rounding integer columns of the
//! reduced solution. The postsolve design follows HiGHS
//! `HighsPostsolveStack` (highs/presolve/HighsPostsolveStack.h) reduced to its
//! primal part.
//!
//! Provenance: T. Achterberg, R. E. Bixby, Z. Gu, E. Rothberg, D. Weninger,
//! "Presolve reductions in mixed integer programming", INFORMS J. Comput.
//! 32(2):473–506 (2020). HiGHS (MIT License, Copyright (c) 2016–2025 the
//! HiGHS developers) `highs/presolve/HPresolve.cpp` read as a reference text;
//! algorithms and constants reimplemented, not translated.
//!
//! No `unsafe`.
// TODO(perf): row → column index lists are plain `Vec`s with linear-search
// removal; HiGHS keeps the matrix in linked row/column splay trees
// (HPresolve.cpp `ARleft/ARright`) for O(log) updates on MIPLIB-scale models.

use crate::bnb::{MipConfig, MipError, MipSolution, MipStatus};
use crate::problem::MipProblem;
use rusolve_core::{CscMatrix, Triplet};
use std::time::Instant;

/// Same acceptance tolerance as the tree's incumbent check (bnb.rs).
const CANDIDATE_TOL: f64 = 1e-6;

/// Coefficients of this magnitude are dropped / treated as equal (HiGHS
/// `small_matrix_value` default).
const SMALL: f64 = 1e-9;
/// Activity-based redundancy margin: a side is dropped only when implied
/// to well within the checker's tolerance.
const REDUNDANT_TOL: f64 = 1e-9;
/// Largest denominator searched per coefficient (HiGHS `integralScale` uses
/// 1000 for its continued-fraction step).
const MAX_DENOM: i64 = 1000;
/// Upper limit on the common denominator of a row.
const MAX_SCALE: i64 = 1_000_000;
/// Doubleton substitution is skipped when `|a_t/a_e|` exceeds this (error
/// amplification in postsolve).
const MAX_SUBST_RATIO: f64 = 1e3;
/// Fill-in guard: the substituted column may appear in at most this many rows.
const MAX_SUBST_COL_LEN: usize = 50;
const MAX_ROUNDS: usize = 50;
/// Fast-pass fixpoints interleaved with cover strengthening.
const MAX_OUTER: usize = 5;

/// Reduction counts, for logging and A/B measurements.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PresolveStats {
    pub rounds: usize,
    pub rows_removed: usize,
    pub cols_removed: usize,
    pub cols_fixed: usize,
    pub dual_fixed: usize,
    pub rows_scaled: usize,
    pub coefs_strengthened: usize,
    /// Coefficients reduced by the cover-based rule.
    pub cover_strengthened: usize,
    pub bounds_tightened: usize,
    pub doubletons: usize,
    pub sides_dropped: usize,
    /// Rows merged into a parallel row.
    pub parallel_rows: usize,
    /// Columns merged into an identical column.
    pub parallel_cols: usize,
    /// Binary columns fixed by column domination.
    pub dominated_fixed: usize,
    /// Columns fixed by set-partitioning row domination.
    pub setppc_fixed: usize,
}

impl std::fmt::Display for PresolveStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "presolve: {} rounds, removed {} rows {} cols (fixed {}, dual-fixed {}, \
             doubleton {}), {} rows gcd-scaled, {} coefs strengthened ({} by covers), \
             {} bounds tightened, {} sides dropped, {} parallel rows, \
             {} parallel cols, {} dominated cols, {} setppc fixings",
            self.rounds,
            self.rows_removed,
            self.cols_removed,
            self.cols_fixed,
            self.dual_fixed,
            self.doubletons,
            self.rows_scaled,
            self.coefs_strengthened,
            self.cover_strengthened,
            self.bounds_tightened,
            self.sides_dropped,
            self.parallel_rows,
            self.parallel_cols,
            self.dominated_fixed,
            self.setppc_fixed
        )
    }
}

/// Result of [`presolve`].
#[derive(Debug, Clone)]
pub enum Presolve {
    /// Presolve proved the problem has no integer-feasible point.
    Infeasible,
    Reduced(Box<Presolved>),
}

/// A reduced problem plus what is needed to map its solutions back.
#[derive(Debug, Clone)]
pub struct Presolved {
    pub problem: MipProblem,
    /// Internal objective = reduced objective + `offset`.
    pub offset: f64,
    pub stats: PresolveStats,
    orig_num_col: usize,
    /// Reduced column → original column.
    col_map: Vec<usize>,
    stack: Vec<Step>,
}

#[derive(Debug, Clone)]
enum Step {
    Fixed {
        col: usize,
        value: f64,
    },
    /// `x_col = (rhs − other_coef·x_other) / coef`.
    Substituted {
        col: usize,
        rhs: f64,
        coef: f64,
        other: usize,
        other_coef: f64,
    },
    /// Identical columns `col` and `into` merged into `z = x_col + x_into`
    /// (stored in `into`); the bounds are the ones at merge time.
    Merged {
        col: usize,
        into: usize,
        lo_col: f64,
        up_col: f64,
        lo_into: f64,
        up_into: f64,
    },
}

impl Presolved {
    /// Maps a reduced-space point back to the original columns. Integer
    /// columns are rounded first so substituted values are computed from
    /// exact integers.
    pub fn postsolve(&self, x: &[f64]) -> Vec<f64> {
        self.postsolve_with(x, true)
    }

    /// [`Self::postsolve`] with optional integer rounding.
    pub fn postsolve_with(&self, x: &[f64], round_integers: bool) -> Vec<f64> {
        let mut out = vec![0.0; self.orig_num_col];
        for (k, &j) in self.col_map.iter().enumerate() {
            out[j] = if round_integers && self.problem.integer[k] {
                x[k].round()
            } else {
                x[k]
            };
        }
        for step in self.stack.iter().rev() {
            match *step {
                Step::Fixed { col, value } => out[col] = value,
                Step::Substituted {
                    col,
                    rhs,
                    coef,
                    other,
                    other_coef,
                } => out[col] = (rhs - other_coef * out[other]) / coef,
                Step::Merged {
                    col,
                    into,
                    lo_col,
                    up_col,
                    lo_into,
                    up_into,
                } => {
                    // Put x_col at a finite bound if it can be, the rest
                    // on x_into; z ∈ [lo_col + lo_into, up_col + up_into]
                    // keeps both in range.
                    let z = out[into];
                    let anchor = if lo_col.is_finite() {
                        lo_col
                    } else if up_col.is_finite() {
                        up_col
                    } else {
                        z - 0.0f64.clamp(lo_into, up_into)
                    };
                    let x_into = (z - anchor).clamp(lo_into, up_into);
                    out[into] = x_into;
                    out[col] = z - x_into;
                }
            }
        }
        out
    }

    /// Reduced column → original column.
    pub fn col_map(&self) -> &[usize] {
        &self.col_map
    }

    /// Maps an original-space point to the reduced space (the inverse of
    /// [`Self::postsolve`] for a *point*, the direction a root restart needs
    /// to carry its incumbent into the new reduced problem; HiGHS keeps
    /// incumbents in the original space and re-reduces them the same way,
    /// `Highs::getReducedPrimalSolution`, MIT).
    ///
    /// Each `Merged { col, into }` step folds the removed column's value
    /// into the survivor (`z = x_col + x_into`) in stack order: a column is
    /// only merged while it still exists, so one forward pass over the stack
    /// sees every merge before its argument disappears. `Fixed` and
    /// `Substituted` steps remove columns whose value postsolve recomputes
    /// from the rows, so the surviving columns simply keep their values.
    ///
    /// Panics if `x.len()` differs from the original column count.
    pub fn forward_map(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.orig_num_col);
        let mut cur = x.to_vec();
        for step in &self.stack {
            if let Step::Merged { col, into, .. } = *step {
                cur[into] += cur[col];
            }
        }
        self.col_map.iter().map(|&j| cur[j]).collect()
    }
}

#[derive(Debug)]
struct Infeasible;

type Res = Result<(), Infeasible>;

/// Presolves `p`. `feastol` is the primal feasibility tolerance used for
/// rounding and infeasibility decisions (HiGHS `primal_feastol`).
pub fn presolve(p: &MipProblem, feastol: f64) -> Presolve {
    presolve_until(p, feastol, None)
}

/// [`presolve`] that stops reducing once `deadline` passes and returns the
/// reductions applied so far (each one is complete, so the partial result
/// is a valid presolve). Polled between passes and, inside the cover-based
/// strengthening, every [`DEADLINE_WORK`] work units — HiGHS checks its
/// presolve time limit the same way, between rules and every few rows of
/// `strengthenInequalities` (HPresolve.cpp:7036, :7130, MIT).
/// Measured: s100's strengthening of its 350k-column
/// rows ran past three minutes under a 60 s limit.
pub fn presolve_until(p: &MipProblem, feastol: f64, deadline: Option<Instant>) -> Presolve {
    let mut w = Work::new(p, feastol);
    w.deadline = deadline;
    match w.run() {
        Ok(()) => Presolve::Reduced(Box::new(w.finish(p))),
        Err(Infeasible) => Presolve::Infeasible,
    }
}

/// Presolve → branch-and-bound on the reduced problem → postsolve. The
/// returned point is verified against the *original* rows and bounds; if
/// postsolve ever produced a violating point the original problem is solved
/// without presolve instead (and the event is reported on stderr), so a
/// presolve bug costs time, never a wrong answer.
///
/// Root restarts: when the root search
/// asks for one, the loop builds the restart problem (root box + active cut
/// rows), presolves it **again on the reduced problem** (HiGHS
/// `performRestart`, HighsMipSolverData.cpp:1317, MIT), carries the
/// incumbent forward through the new stack, and re-enters the search. The
/// postsolve stacks chain; the final point postsolves through every level
/// and is verified against the original exactly as for one level.
pub fn solve_presolved(p: &MipProblem, config: &MipConfig) -> Result<MipSolution, MipError> {
    use crate::restart;
    let start = Instant::now();
    let plain = MipConfig {
        presolve: false,
        ..*config
    };
    let deadline = config
        .time_limit
        .map(|t| start + std::time::Duration::from_secs_f64(t));
    let pre = match presolve_until(p, config.feastol, deadline) {
        Presolve::Infeasible => {
            return Ok(MipSolution {
                status: MipStatus::Infeasible,
                objective: f64::NAN,
                bound: p.caller_objective(f64::INFINITY),
                x: vec![],
                nodes: 0,
                lp_iterations: 0,
                first_incumbent: None,
                incumbents: vec![],
                cuts: Default::default(),
                lp_stats: Default::default(),
                redcost: Default::default(),
                submip: Default::default(),
                conflict: Default::default(),
                restarts: 0,
                restart: None,
                cliques: crate::cuts::CliqueFixStats::default(),
                probing: crate::probing::ProbingStats::default(),
            });
        }
        Presolve::Reduced(pre) => pre,
    };
    if config.audit {
        // Audit: presolve sizes and the
        // objective map `caller = sense·(internal + offset)`.
        eprintln!(
            "audit presolve orig_rows={} orig_cols={} orig_nnz={} rows={} cols={} nnz={} \
             offset={:e} maximize={} t={:.6}",
            p.num_row,
            p.num_col,
            p.matrix.nnz(),
            pre.problem.num_row,
            pre.problem.num_col,
            pre.problem.matrix.nnz(),
            pre.offset,
            p.maximize as u8,
            start.elapsed().as_secs_f64()
        );
    }
    if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
        eprintln!(
            "{} -> {} rows {} cols ({:.3}s)",
            pre.stats,
            pre.problem.num_row,
            pre.problem.num_col,
            start.elapsed().as_secs_f64()
        );
    }
    // Presolve levels: chain[0] reduces `p`, chain[k] reduces the restart
    // problem of level k-1. A restart round needs a single search (the race
    // shares one problem across racers), so `threads > 1` never restarts.
    let mut chain: Vec<Box<Presolved>> = vec![pre];
    let mut carry: Option<restart::Carry> = None;
    let mut restarts: u32 = 0;
    let mut lp_stats_total = crate::engine::LpStats::default();
    let sol;
    loop {
        let q = &chain[chain.len() - 1].problem;
        let mut round = if q.num_col == 0 {
            let c = carry.take().unwrap_or_default();
            MipSolution {
                status: MipStatus::Optimal,
                objective: q.caller_objective(0.0),
                bound: q.caller_objective(0.0),
                x: vec![],
                nodes: c.nodes,
                lp_iterations: c.lp_iterations,
                first_incumbent: c.first_incumbent,
                incumbents: c.incumbents,
                cuts: Default::default(),
                lp_stats: Default::default(),
                redcost: Default::default(),
                submip: Default::default(),
                conflict: Default::default(),
                restarts,
                restart: None,
                cliques: crate::cuts::CliqueFixStats::default(),
                probing: crate::probing::ProbingStats::default(),
            }
        } else {
            // `objective_cutoff` is caller-sense for the problem a config is
            // applied to: re-reference it into this level's reduced space
            // (internal_p = internal_q + Σ offsets, as `lift` below).
            let cum_offset: f64 = chain.iter().map(|pre| pre.offset).sum();
            let inner = MipConfig {
                time_limit: config
                    .time_limit
                    .map(|t| (t - start.elapsed().as_secs_f64()).max(0.0)),
                objective_cutoff: config
                    .objective_cutoff
                    .map(|c| q.caller_objective(p.caller_objective(c) - cum_offset)),
                ..plain
            };
            if config.audit {
                eprintln!(
                    "audit round level={} offset={:e} rows={} cols={} t={:.6}",
                    chain.len() - 1,
                    cum_offset,
                    q.num_row,
                    q.num_col,
                    start.elapsed().as_secs_f64()
                );
            }
            if inner.threads > 1 {
                // Racing shares one problem across racers: no restart.
                debug_assert!(carry.is_none());
                crate::bnb::solve(q, &inner)?
            } else {
                crate::bnb::solve_tree(q, &inner, None, carry.take())?
            }
        };
        lp_stats_total.absorb(&round.lp_stats);
        let Some(req) = round.restart.take() else {
            sol = round;
            break;
        };
        // The restart problem: current problem + root box + active cut rows,
        // then presolve it all over again (HiGHS re-presolves the reduced
        // LP, not the original model, HighsMipSolverData.cpp:1386–1409).
        let base = restart::restart_problem(q, &req);
        match presolve_until(&base, config.feastol, deadline) {
            Presolve::Infeasible => {
                // The root box plus the valid cuts hold no integer point.
                // Reduced-cost fixing only removes values that cannot be
                // part of a point better than the incumbent, so: with an
                // incumbent it is optimal; without one the problem is
                // infeasible.
                round.status = if req.incumbent.is_some() {
                    MipStatus::Optimal
                } else {
                    MipStatus::Infeasible
                };
                match req.incumbent {
                    Some((z, x)) => {
                        round.objective = q.caller_objective(z);
                        round.bound = q.caller_objective(z);
                        round.x = x;
                    }
                    None => {
                        round.objective = f64::NAN;
                        round.bound = q.caller_objective(f64::INFINITY);
                        round.x = vec![];
                    }
                }
                sol = round;
                break;
            }
            Presolve::Reduced(pre) => {
                let mut carry_next = req.carry;
                carry_next.incumbent = req
                    .incumbent
                    .map(|(z, x)| (z - pre.offset, pre.forward_map(&x)));
                // Re-reference the ladder into the new space (HiGHS
                // re-references the bounds and limits across the offset
                // change, HighsMipSolverData.cpp:1367–1376).
                carry_next
                    .incumbents
                    .iter_mut()
                    .for_each(|f| f.value -= pre.offset);
                if let Some(f) = &mut carry_next.first_incumbent {
                    f.value -= pre.offset;
                }
                if std::env::var_os("RUSOLVE_MIP_LOG").is_some() {
                    eprintln!(
                        "root restart {}: {} rows {} cols, offset {} ({:.3}s)",
                        carry_next.restarts_done,
                        pre.problem.num_row,
                        pre.problem.num_col,
                        pre.offset,
                        start.elapsed().as_secs_f64()
                    );
                }
                restarts = carry_next.restarts_done;
                carry = Some(carry_next);
                chain.push(pre);
            }
        }
    }
    let mut sol = sol;
    sol.restarts = restarts;
    let q_last = &chain[chain.len() - 1].problem;
    let has_point = sol.status != MipStatus::Infeasible
        && (!sol.x.is_empty() || (q_last.num_col == 0 && sol.status == MipStatus::Optimal));
    let postsolve_chain = |x: &[f64], round_integers: bool| -> Vec<f64> {
        let mut x = x.to_vec();
        for pre in chain.iter().rev() {
            x = pre.postsolve_with(&x, round_integers);
        }
        x
    };
    let x = if has_point {
        let violation = |x: &[f64]| {
            let v = p.violation(x);
            v.row.max(v.bound).max(v.integrality)
        };
        let rounded = postsolve_chain(&sol.x, true);
        let x = if violation(&rounded) <= CANDIDATE_TOL {
            rounded
        } else {
            postsolve_chain(&sol.x, false)
        };
        if violation(&x) > CANDIDATE_TOL {
            eprintln!(
                "mip presolve: postsolved point violates the original by {:e}; \
                 re-solving without presolve",
                violation(&x)
            );
            return crate::bnb::solve(p, &plain);
        }
        x
    } else {
        vec![]
    };
    // Reduced internal values + offsets = original internal values: the
    // chain telescopes `p_internal = q_last_internal + Σ offsets` (each
    // level's `internal = reduced + offset`, and a restart problem reuses
    // the parent's cost with no new offset). `q_last` is the space `sol`'s
    // values live in; with one level this is the original single-offset lift.
    let total_offset: f64 = chain.iter().map(|pre| pre.offset).sum();
    let lift = |v: f64| p.caller_objective(q_last.caller_objective(v) + total_offset);
    let objective = match sol.status {
        MipStatus::Unbounded => sol.objective,
        _ if has_point => p.caller_objective(p.objective(&x)),
        _ => f64::NAN,
    };
    // The root-bound stats are caller-sense values of the reduced problem.
    let cuts = crate::cuts::CutStats {
        root_bound_before: lift(q_last.caller_objective(sol.cuts.root_bound_before)),
        root_bound_after: lift(q_last.caller_objective(sol.cuts.root_bound_after)),
        ..sol.cuts
    };
    Ok(MipSolution {
        status: sol.status,
        objective,
        bound: lift(sol.bound),
        x,
        nodes: sol.nodes,
        lp_iterations: sol.lp_iterations,
        first_incumbent: sol.first_incumbent.map(|f| crate::bnb::FirstIncumbent {
            value: lift(q_last.caller_objective(f.value)),
            ..f
        }),
        incumbents: sol
            .incumbents
            .iter()
            .map(|f| crate::bnb::FirstIncumbent {
                value: lift(q_last.caller_objective(f.value)),
                ..*f
            })
            .collect(),
        cuts,
        lp_stats: lp_stats_total,
        redcost: sol.redcost,
        submip: sol.submip,
        conflict: sol.conflict,
        cliques: sol.cliques,
        probing: sol.probing,
        restarts: sol.restarts,
        restart: None,
    })
}

struct Work {
    feastol: f64,
    /// Row-wise matrix (original row/col indices), the primary storage.
    rows: Vec<Vec<(usize, f64)>>,
    /// Rows containing each column (kept exactly in sync with `rows`).
    cols: Vec<Vec<usize>>,
    row_lo: Vec<f64>,
    row_up: Vec<f64>,
    row_alive: Vec<bool>,
    lo: Vec<f64>,
    up: Vec<f64>,
    cost: Vec<f64>,
    integer: Vec<bool>,
    col_alive: Vec<bool>,
    offset: f64,
    stack: Vec<Step>,
    stats: PresolveStats,
    changed: bool,
    /// Run the pair reductions (parallel rows/columns, dominated columns,
    /// setppc conflicts); off only in tests isolating single-row rules.
    pairs: bool,
    /// Wall-clock stop ([`presolve_until`]); `expired` latches once seen.
    deadline: Option<Instant>,
    expired: bool,
    /// Work units since the last deadline poll.
    work: usize,
}

/// Work units (row entries scanned or sorted) between two deadline polls.
const DEADLINE_WORK: usize = 1 << 20;

impl Work {
    fn new(p: &MipProblem, feastol: f64) -> Self {
        let mut rows = vec![Vec::new(); p.num_row];
        let mut cols = vec![Vec::new(); p.num_col];
        let (ptr, idx, val) = (p.matrix.col_ptr(), p.matrix.row_idx(), p.matrix.values());
        for j in 0..p.num_col {
            for k in ptr[j]..ptr[j + 1] {
                if val[k] != 0.0 {
                    let i = idx[k] as usize;
                    rows[i].push((j, val[k]));
                    cols[j].push(i);
                }
            }
        }
        Self {
            feastol,
            rows,
            cols,
            row_lo: p.row_lower.clone(),
            row_up: p.row_upper.clone(),
            row_alive: vec![true; p.num_row],
            lo: p.col_lower.clone(),
            up: p.col_upper.clone(),
            cost: p.cost.clone(),
            integer: p.integer.clone(),
            col_alive: vec![true; p.num_col],
            offset: 0.0,
            stack: Vec::new(),
            stats: PresolveStats::default(),
            changed: false,
            pairs: true,
            deadline: None,
            expired: false,
            work: 0,
        }
    }

    /// Charges `units` of work and reports whether the deadline has passed,
    /// reading the clock only every [`DEADLINE_WORK`] units.
    fn out_of_time(&mut self, units: usize) -> bool {
        if self.expired {
            return true;
        }
        let Some(d) = self.deadline else {
            return false;
        };
        self.work += units;
        if self.work >= DEADLINE_WORK {
            self.work = 0;
            self.expired = Instant::now() >= d;
        }
        self.expired
    }

    /// Unconditional poll, between passes.
    fn past_deadline(&mut self) -> bool {
        self.out_of_time(DEADLINE_WORK)
    }

    /// Fast passes to a fixpoint, then the expensive cover-based
    /// strengthening, then the fast passes again (the order of HiGHS's
    /// `presolve` loop around `strengthenInequalities`, HPresolve.cpp:5957).
    fn run(&mut self) -> Res {
        for _ in 0..MAX_OUTER {
            for _ in 0..MAX_ROUNDS {
                self.changed = false;
                self.stats.rounds += 1;
                self.col_pass()?;
                self.row_pass()?;
                // TODO(perf): tightened continuous bounds cost jobshop_bigM
                // 1.08x nodes (sgm seeds 1-12, 381 -> 413, mip_perf A/B
                // 2026-10-22) while helping lot_sizing/knapsack_cover; revisit
                // tightening only integer columns if node LPs stay sensitive.
                self.bound_pass()?;
                self.doubleton_pass()?;
                if !self.changed || self.past_deadline() {
                    break;
                }
            }
            if self.past_deadline() {
                break;
            }
            self.changed = false;
            if self.pairs {
                self.setppc_conflicts();
                self.dominated_columns();
                self.parallel_rows()?;
                self.parallel_cols();
            }
            if self.past_deadline() {
                break;
            }
            self.strengthen_inequalities();
            if !self.changed || self.past_deadline() {
                break;
            }
        }
        Ok(())
    }

    // ---------------------------------------------------------------- matrix

    fn coef(&self, i: usize, j: usize) -> f64 {
        self.rows[i].iter().find(|e| e.0 == j).map_or(0.0, |e| e.1)
    }

    /// Sets `A[i, j] = v`, dropping it when `|v| ≤ SMALL`.
    fn set_coef(&mut self, i: usize, j: usize, v: f64) {
        let pos = self.rows[i].iter().position(|e| e.0 == j);
        match (pos, v.abs() > SMALL) {
            (Some(k), true) => self.rows[i][k].1 = v,
            (Some(k), false) => {
                self.rows[i].swap_remove(k);
                self.cols[j].retain(|&r| r != i);
            }
            (None, true) => {
                self.rows[i].push((j, v));
                self.cols[j].push(i);
            }
            (None, false) => {}
        }
    }

    fn remove_row(&mut self, i: usize) {
        for &(j, _) in &self.rows[i] {
            self.cols[j].retain(|&r| r != i);
        }
        self.rows[i].clear();
        self.row_alive[i] = false;
        self.stats.rows_removed += 1;
        self.changed = true;
    }

    /// Removes column `j` at value `v`, moving its contribution into the
    /// row sides and the objective offset.
    fn fix_col(&mut self, j: usize, v: f64) {
        let rows = std::mem::take(&mut self.cols[j]);
        for i in rows {
            let k = self.rows[i].iter().position(|e| e.0 == j).unwrap();
            let a = self.rows[i].swap_remove(k).1;
            self.row_lo[i] -= a * v;
            self.row_up[i] -= a * v;
        }
        self.offset += self.cost[j] * v;
        self.lo[j] = v;
        self.up[j] = v;
        self.col_alive[j] = false;
        self.stack.push(Step::Fixed { col: j, value: v });
        self.stats.cols_removed += 1;
        self.changed = true;
    }

    /// `(min, max)` activity of row `i` over the column box (±∞ when a
    /// contribution is unbounded).
    fn activity(&self, i: usize) -> (f64, f64) {
        let (mut min, mut max) = (0.0, 0.0);
        for &(j, a) in &self.rows[i] {
            if a > 0.0 {
                min += a * self.lo[j];
                max += a * self.up[j];
            } else {
                min += a * self.up[j];
                max += a * self.lo[j];
            }
        }
        // ∞ − ∞ cannot occur: each sum only collects one sign of infinity.
        (min, max)
    }

    // ---------------------------------------------------------------- bounds

    fn tighten_lower(&mut self, j: usize, v: f64) -> Res {
        let v = if self.integer[j] {
            (v - self.feastol).ceil()
        } else {
            v
        };
        // Any tightening is kept (the bound may replace a removed row); only
        // one that moves by more than feastol counts as progress.
        if v > self.lo[j] {
            if v > self.lo[j] + self.feastol {
                self.stats.bounds_tightened += 1;
                self.changed = true;
            }
            self.lo[j] = v;
        }
        self.check_col(j)
    }

    fn tighten_upper(&mut self, j: usize, v: f64) -> Res {
        let v = if self.integer[j] {
            (v + self.feastol).floor()
        } else {
            v
        };
        if v < self.up[j] {
            if v < self.up[j] - self.feastol {
                self.stats.bounds_tightened += 1;
                self.changed = true;
            }
            self.up[j] = v;
        }
        self.check_col(j)
    }

    fn check_col(&mut self, j: usize) -> Res {
        if self.lo[j] > self.up[j] + self.feastol {
            return Err(Infeasible);
        }
        if self.lo[j] > self.up[j] {
            // Crossed within tolerance (continuous only: integer bounds
            // are integral here).
            let mid = 0.5 * (self.lo[j] + self.up[j]);
            self.lo[j] = mid;
            self.up[j] = mid;
        }
        Ok(())
    }

    // ---------------------------------------------------------------- columns

    fn col_pass(&mut self) -> Res {
        // Locks from one row sweep. Fixing columns during the pass only
        // removes locks, so stale entries are conservative.
        let n = self.lo.len();
        let (mut down_locked, mut up_locked) = (vec![false; n], vec![false; n]);
        for i in 0..self.rows.len() {
            if !self.row_alive[i] {
                continue;
            }
            let (lo_fin, up_fin) = (self.row_lo[i].is_finite(), self.row_up[i].is_finite());
            for &(j, a) in &self.rows[i] {
                let (down, up) = if a > 0.0 {
                    (lo_fin, up_fin)
                } else {
                    (up_fin, lo_fin)
                };
                down_locked[j] |= down;
                up_locked[j] |= up;
            }
        }
        for j in 0..n {
            if !self.col_alive[j] {
                continue;
            }
            if self.integer[j] {
                let (l, u) = (
                    (self.lo[j] - self.feastol).ceil(),
                    (self.up[j] + self.feastol).floor(),
                );
                if l != self.lo[j] || u != self.up[j] {
                    self.lo[j] = l;
                    self.up[j] = u;
                    self.changed = true;
                }
            }
            self.check_col(j)?;
            if self.lo[j] == self.up[j] {
                self.stats.cols_fixed += 1;
                self.fix_col(j, self.lo[j]);
                continue;
            }
            self.dual_fix(j, down_locked[j], up_locked[j]);
        }
        Ok(())
    }

    /// Empty columns and dual fixing (Achterberg et al. 2020 §4.4): if
    /// `c_j ≥ 0` and no row is violated by decreasing `x_j`, every
    /// solution stays feasible and no worse with `x_j = l_j` (mirrored for
    /// the upper bound). An infinite target bound is left alone so the tree
    /// still sees unboundedness.
    fn dual_fix(&mut self, j: usize, down_locked: bool, up_locked: bool) {
        let c = self.cost[j];
        let target = if c >= 0.0 && !down_locked && self.lo[j].is_finite() {
            Some(self.lo[j])
        } else if c <= 0.0 && !up_locked && self.up[j].is_finite() {
            Some(self.up[j])
        } else if c == 0.0 && !down_locked && !up_locked {
            // Free, cost-free column in no binding row: any value works.
            Some(0.0f64.clamp(self.lo[j], self.up[j]))
        } else {
            None
        };
        if let Some(v) = target {
            self.stats.dual_fixed += 1;
            self.fix_col(j, v);
        }
    }

    // ---------------------------------------------------------------- rows

    fn row_pass(&mut self) -> Res {
        for i in 0..self.rows.len() {
            if self.row_alive[i] {
                self.row_presolve(i)?;
            }
        }
        Ok(())
    }

    fn row_presolve(&mut self, i: usize) -> Res {
        let tol = self.feastol;
        if self.rows[i].is_empty() {
            if self.row_lo[i] > tol || self.row_up[i] < -tol {
                return Err(Infeasible);
            }
            self.remove_row(i);
            return Ok(());
        }
        let (min, max) = self.activity(i);
        if min > self.row_up[i] + tol || max < self.row_lo[i] - tol {
            return Err(Infeasible);
        }
        if self.rows[i].len() == 1 {
            return self.singleton_row(i);
        }
        if self.drop_redundant_sides(i, min, max) {
            return Ok(());
        }
        if self.rows[i].iter().all(|&(j, _)| self.integer[j]) {
            self.scale_integer_row(i)?;
        }
        self.strengthen_coefs(i);
        Ok(())
    }

    /// Activity-based bound tightening on every row (Savelsbergh 1994;
    /// Achterberg et al. 2020 §3.1): with `minact_{−j}` the minimum activity
    /// of the other entries, `a_j x_j ≤ hi − minact_{−j}` (and the mirror
    /// with `lo` and the maximum). Implied bounds on continuous columns are
    /// what lets coefficient strengthening tighten big-M / VUB rows
    /// (`f − 27y ≤ 0` → `f − 7y ≤ 0` once flow conservation bounds `f ≤ 7`).
    /// Candidates with `|b|·1e-14 >
    /// feastol` are ignored (HighsDomain.cpp:1349 constant, see
    /// propagate.rs).
    fn bound_pass(&mut self) -> Res {
        for i in 0..self.rows.len() {
            if self.row_alive[i] && self.rows[i].len() >= 2 {
                self.tighten_from_row(i)?;
            }
        }
        Ok(())
    }

    fn tighten_from_row(&mut self, i: usize) -> Res {
        // Finite parts and counts of infinite contributions.
        let (mut min_fin, mut min_inf, mut max_fin, mut max_inf) = (0.0, 0usize, 0.0, 0usize);
        for &(j, a) in &self.rows[i] {
            let (cmin, cmax) = contribution(a, self.lo[j], self.up[j]);
            if cmin.is_finite() {
                min_fin += cmin;
            } else {
                min_inf += 1;
            }
            if cmax.is_finite() {
                max_fin += cmax;
            } else {
                max_inf += 1;
            }
        }
        let (row_lo, row_up) = (self.row_lo[i], self.row_up[i]);
        let usable_up = row_up.is_finite() && min_inf <= 1;
        let usable_lo = row_lo.is_finite() && max_inf <= 1;
        if !usable_up && !usable_lo {
            return Ok(());
        }
        let entries = self.rows[i].clone();
        for (j, a) in entries {
            let (cmin, cmax) = contribution(a, self.lo[j], self.up[j]);
            // a_j x_j ≤ row_up − (min activity of the others).
            if usable_up {
                let rest = match (min_inf, cmin.is_finite()) {
                    (0, _) => Some(min_fin - cmin),
                    (1, false) => Some(min_fin),
                    _ => None,
                };
                if let Some(rest) = rest {
                    let b = (row_up - rest) / a;
                    if b.abs() * 1e-14 <= self.feastol {
                        if a > 0.0 {
                            self.tighten_upper(j, b)?;
                        } else {
                            self.tighten_lower(j, b)?;
                        }
                    }
                }
            }
            // a_j x_j ≥ row_lo − (max activity of the others).
            if usable_lo {
                let rest = match (max_inf, cmax.is_finite()) {
                    (0, _) => Some(max_fin - cmax),
                    (1, false) => Some(max_fin),
                    _ => None,
                };
                if let Some(rest) = rest {
                    let b = (row_lo - rest) / a;
                    if b.abs() * 1e-14 <= self.feastol {
                        if a > 0.0 {
                            self.tighten_lower(j, b)?;
                        } else {
                            self.tighten_upper(j, b)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `lo ≤ a x_j ≤ hi` ⇒ bounds on `x_j`; the row is then removed.
    fn singleton_row(&mut self, i: usize) -> Res {
        let (j, a) = self.rows[i][0];
        // Coefficient writes go through `set_coef` (drops |a| ≤ SMALL) or
        // stay above feastol (`strengthen_coefs`); a zero here would turn
        // the sides into ±∞/NaN bounds.
        debug_assert!(a.abs() > SMALL, "row {i} keeps a zero entry for col {j}");
        let (lo, hi) = if a > 0.0 {
            (self.row_lo[i] / a, self.row_up[i] / a)
        } else {
            (self.row_up[i] / a, self.row_lo[i] / a)
        };
        if lo > f64::NEG_INFINITY {
            self.tighten_lower(j, lo)?;
        }
        if hi < f64::INFINITY {
            self.tighten_upper(j, hi)?;
        }
        self.remove_row(i);
        Ok(())
    }

    /// Drops sides implied by the activity bounds; removes the row when both
    /// are. Returns `true` iff the row was removed.
    fn drop_redundant_sides(&mut self, i: usize, min: f64, max: f64) -> bool {
        let lo_red = self.row_lo[i] == f64::NEG_INFINITY
            || min >= self.row_lo[i] - REDUNDANT_TOL * self.row_lo[i].abs().max(1.0);
        let up_red = self.row_up[i] == f64::INFINITY
            || max <= self.row_up[i] + REDUNDANT_TOL * self.row_up[i].abs().max(1.0);
        if lo_red && up_red {
            self.remove_row(i);
            return true;
        }
        if lo_red && self.row_lo[i] > f64::NEG_INFINITY {
            self.row_lo[i] = f64::NEG_INFINITY;
            self.stats.sides_dropped += 1;
            self.changed = true;
        }
        if up_red && self.row_up[i] < f64::INFINITY {
            self.row_up[i] = f64::INFINITY;
            self.stats.sides_dropped += 1;
            self.changed = true;
        }
        false
    }

    /// Row gcd reasoning on an all-integer row: scale to coprime integer
    /// coefficients and round the sides inward (HiGHS HPresolve.cpp:3866).
    fn scale_integer_row(&mut self, i: usize) -> Res {
        let vals: Vec<f64> = self.rows[i].iter().map(|e| e.1).collect();
        let Some((scale, ints)) = integral_scale(&vals) else {
            return Ok(());
        };
        let tol = self.feastol;
        let lo = if self.row_lo[i].is_finite() {
            (self.row_lo[i] * scale - tol).ceil()
        } else {
            f64::NEG_INFINITY
        };
        let hi = if self.row_up[i].is_finite() {
            (self.row_up[i] * scale + tol).floor()
        } else {
            f64::INFINITY
        };
        if lo > hi {
            return Err(Infeasible);
        }
        let max_int = ints.iter().fold(0i64, |m, &v| m.max(v.abs()));
        if max_int > 1000 && scale > 100.0 {
            // Large scaled values: keep the coefficients, only tighten the
            // sides (scaled back) when that gains more than feastol.
            if lo / scale > self.row_lo[i] + tol {
                self.row_lo[i] = lo / scale;
                self.stats.rows_scaled += 1;
                self.changed = true;
            }
            if hi / scale < self.row_up[i] - tol {
                self.row_up[i] = hi / scale;
                self.stats.rows_scaled += 1;
                self.changed = true;
            }
            return Ok(());
        }
        let mut differs = lo != self.row_lo[i] || hi != self.row_up[i];
        for (e, &v) in self.rows[i].iter_mut().zip(&ints) {
            let v = v as f64;
            differs |= e.1 != v;
            e.1 = v;
        }
        self.row_lo[i] = lo;
        self.row_up[i] = hi;
        if differs {
            self.stats.rows_scaled += 1;
            self.changed = true;
        }
        Ok(())
    }

    /// Coefficient strengthening on a one-sided row (Achterberg et al. 2020
    /// §3.2; HiGHS `strengthenCoefs`, HPresolve.cpp:4217). With `dir = +1`
    /// for `aᵀx ≤ b` and `−1` for `aᵀx ≥ b` (read as `−aᵀx ≤ −b`), and
    /// `d = maxact − b > 0`: an integer column with `dir·a_j > d` takes
    /// coefficient `d` and the rhs drops by `(dir·a_j − d)·u_j`; with
    /// `dir·a_j < −d` it takes `−d` and the rhs grows by `(−d − dir·a_j)·l_j`.
    /// At the bound the row is unchanged; one unit away it is redundant both
    /// before and after.
    fn strengthen_coefs(&mut self, i: usize) {
        let dir = match (self.row_lo[i].is_finite(), self.row_up[i].is_finite()) {
            (false, true) => 1.0,
            (true, false) => -1.0,
            _ => return,
        };
        if !self.rows[i].iter().any(|&(j, _)| self.integer[j]) {
            return;
        }
        let (min, max) = self.activity(i);
        let (maxact, mut rhs) = if dir > 0.0 {
            (max, self.row_up[i])
        } else {
            (-min, -self.row_lo[i])
        };
        if !maxact.is_finite() {
            return;
        }
        let d = maxact - rhs;
        if d <= self.feastol {
            return;
        }
        let tol = self.feastol;
        let mut count = 0;
        for k in 0..self.rows[i].len() {
            let (j, a) = self.rows[i][k];
            if !self.integer[j] {
                continue;
            }
            let v = dir * a;
            if v > d + tol {
                rhs += (d - v) * self.up[j];
                self.rows[i][k].1 = dir * d;
                count += 1;
            } else if v < -d - tol {
                rhs += (-d - v) * self.lo[j];
                self.rows[i][k].1 = -dir * d;
                count += 1;
            }
        }
        if count > 0 {
            if dir > 0.0 {
                self.row_up[i] = rhs;
            } else {
                self.row_lo[i] = -rhs;
            }
            self.stats.coefs_strengthened += count;
            self.changed = true;
        }
    }

    /// Cover-based coefficient strengthening (HiGHS
    /// `HPresolve::strengthenInequalities`, HPresolve.cpp:7023–7254), a
    /// generalization of [`Self::strengthen_coefs`].
    ///
    /// A one-sided row over finitely bounded columns is read as
    /// `Σ w_j z_j + (continuous) ≥ V`, with `z_j ≥ 0` the distance of `x_j`
    /// from the bound that maximizes its contribution, `w_j = |a_j|` and `V`
    /// the max activity minus the rhs. Repeatedly, the integer items that can
    /// each cover the remaining need `λ` (after all other items are taken at
    /// full range) form the valid integer inequality
    /// `Σ_K c_i z_i ≥ r` with `c_i = ⌈min(w_i, λ)/α⌉`, `r = ⌈λ/α⌉` (α the
    /// smallest weight in the cover); `step` times it is split off the row
    /// as a new non-negative integer slack. What is left of `V` bounds every
    /// remaining weight: `w_j > V` ⇒ `w_j := V` exactly as in the simple
    /// rule, but on the reduced weights. The split is an identity and every
    /// slack is non-negative on integer points, so the strengthened row has
    /// the same integer solutions.
    fn strengthen_inequalities(&mut self) {
        struct Item {
            rc: f64,
            upper: f64,
            /// Position in the row (`None` for a cover slack).
            pos: Option<usize>,
        }
        let tol = self.feastol;
        for i in 0..self.rows.len() {
            if !self.row_alive[i] || self.rows[i].len() <= 1 {
                continue;
            }
            let scale = match (self.row_lo[i].is_finite(), self.row_up[i].is_finite()) {
                (true, false) => -1.0,
                (false, true) => 1.0,
                _ => continue,
            };
            let mut maxviol = if scale < 0.0 {
                self.row_lo[i]
            } else {
                -self.row_up[i]
            };
            let mut cont = 0.0;
            let mut items: Vec<Item> = Vec::new();
            let mut skip = false;
            for (k, &(j, a)) in self.rows[i].iter().enumerate() {
                if !self.lo[j].is_finite() || !self.up[j].is_finite() {
                    skip = true;
                    break;
                }
                let mut w = a * scale;
                let ub = self.up[j] - self.lo[j];
                if w > 0.0 {
                    maxviol += self.up[j] * w;
                } else {
                    maxviol += self.lo[j] * w;
                    w = -w;
                }
                if ub <= tol || w <= tol {
                    continue;
                }
                if !self.integer[j] {
                    cont += w * ub;
                    continue;
                }
                items.push(Item {
                    rc: w,
                    upper: ub,
                    pos: Some(k),
                });
            }
            if skip || maxviol <= tol {
                continue; // a redundant row is the row pass's job
            }
            let small = (100.0 * tol).max(tol * maxviol);
            let mut idx: Vec<usize> = (0..items.len()).collect();
            let mut cover: Vec<usize> = Vec::new();
            for _ in 0..(4 * items.len() + 8) {
                if maxviol - cont <= small || idx.is_empty() {
                    break;
                }
                if self.out_of_time(idx.len()) {
                    // Abandon this row before any update: rows already
                    // strengthened stay valid on their own.
                    return;
                }
                idx.sort_by(|&a, &b| items[b].rc.total_cmp(&items[a].rc).then(b.cmp(&a)));
                let mut lambda = maxviol - cont;
                cover.clear();
                for &t in idx.iter().rev() {
                    let delta = items[t].upper * items[t].rc;
                    if items[t].upper <= 1000.0 && items[t].rc > small && lambda - delta <= small {
                        cover.push(t);
                    } else {
                        lambda -= delta;
                    }
                }
                if cover.is_empty() || lambda <= small {
                    break;
                }
                // Smallest weight above 1e-3, else the largest (the order
                // HiGHS's min_element comparator produces).
                let al = cover
                    .iter()
                    .map(|&t| items[t].rc)
                    .filter(|&rc| rc > 1e-3)
                    .reduce(f64::min)
                    .unwrap_or_else(|| cover.iter().map(|&t| items[t].rc).fold(0.0, f64::max));
                let cover_rhs = (lambda / al - tol).ceil().max(1.0);
                let mut slack_upper = -cover_rhs;
                let mut step = f64::INFINITY;
                let coefs: Vec<f64> = cover
                    .iter()
                    .map(|&t| {
                        let c = (items[t].rc.min(lambda) / al - SMALL).ceil();
                        slack_upper += items[t].upper * c;
                        step = step.min(items[t].rc / c);
                        c
                    })
                    .collect();
                step = step.min(maxviol / cover_rhs);
                maxviol -= step * cover_rhs;
                for (&t, &c) in cover.iter().zip(&coefs) {
                    items[t].rc -= step * c;
                }
                idx.retain(|&t| items[t].rc > tol);
                idx.push(items.len());
                items.push(Item {
                    rc: step,
                    upper: slack_upper,
                    pos: None,
                });
            }
            let threshold = maxviol + tol;
            let mut side = if scale < 0.0 {
                -self.row_lo[i]
            } else {
                self.row_up[i]
            };
            let mut updates = Vec::new();
            for &t in &idx {
                let (Some(k), rc) = (items[t].pos, items[t].rc) else {
                    continue;
                };
                if rc.abs() <= threshold {
                    continue;
                }
                let delta = rc - maxviol;
                let (j, a) = self.rows[i][k];
                let a_scaled = a * scale;
                let new_scaled = if a_scaled > 0.0 {
                    side -= delta * self.up[j];
                    a_scaled - delta
                } else {
                    side += delta * self.lo[j];
                    a_scaled + delta
                };
                updates.push((j, new_scaled * scale));
            }
            // A coefficient reduced to zero leaves the row (HiGHS
            // `addToMatrix` unlinks |a| ≤ small_matrix_value,
            // HPresolve.cpp:2106–2109). Kept as an explicit 0 it reached
            // `singleton_row` as `x/0` bounds and a false Infeasible
            // (block_angular s46). Applied after the loop: `set_coef`
            // reorders the row, invalidating the stored positions.
            let count = updates.len();
            for (j, v) in updates {
                self.set_coef(i, j, v);
            }
            if count > 0 {
                if scale < 0.0 {
                    self.row_lo[i] = -side;
                } else {
                    self.row_up[i] = side;
                }
                self.stats.cover_strengthened += count;
                self.changed = true;
            }
        }
    }

    // ---------------------------------------------------------------- pairs

    /// Parallel rows (HiGHS `detectParallelRowsAndCols`, HPresolve.cpp:7267,
    /// the case without extra column singletons; Achterberg et al. 2020
    /// "parallel rows"): if row `i = s · row r` entry by entry, `i`'s sides
    /// are moved into `r`'s space (swapped when `s < 0`), intersected with
    /// `r`'s, and `i` is deleted. Pairs that would turn two opposite
    /// one-sided rows into a new ranged row are left alone (see the
    /// TODO(perf) below): on disjunctive big-M pairs (`dA: S_a − S_b + M z ≥
    /// p_b`, `dB: S_b − S_a − M z ≥ p_a − M`) HiGHS/SCIP form the ranged row
    /// `[p_b, M − p_a]` (jobshop_bigM).
    fn parallel_rows(&mut self) -> Res {
        use std::collections::HashMap;
        let mut seen: HashMap<Vec<(usize, i64)>, usize> = HashMap::new();
        for i in 0..self.rows.len() {
            if !self.row_alive[i] || self.rows[i].len() < 2 {
                continue;
            }
            let mut entries = self.rows[i].clone();
            entries.sort_unstable_by_key(|e| e.0);
            let pivot = entries[0].1;
            // Hash key only; a match is verified with tolerances below.
            let key: Vec<(usize, i64)> = entries
                .iter()
                .map(|&(j, a)| (j, (a / pivot * 1e8).round() as i64))
                .collect();
            let Some(&r) = seen.get(&key) else {
                seen.insert(key, i);
                continue;
            };
            let mut other = self.rows[r].clone();
            other.sort_unstable_by_key(|e| e.0);
            let s = pivot / other[0].1;
            let parallel = other.len() == entries.len()
                && entries.iter().zip(&other).all(|(&(j, a), &(jr, b))| {
                    j == jr && (a - s * b).abs() <= SMALL * a.abs().max(1.0)
                });
            if !parallel {
                continue;
            }
            let (lo_i, up_i) = if s > 0.0 {
                (self.row_lo[i] / s, self.row_up[i] / s)
            } else {
                (self.row_up[i] / s, self.row_lo[i] / s)
            };
            let lo = self.row_lo[r].max(lo_i);
            let up = self.row_up[r].min(up_i);
            let one_sided = |l: f64, u: f64| l.is_finite() != u.is_finite();
            if lo.is_finite()
                && up.is_finite()
                && lo < up
                && one_sided(self.row_lo[r], self.row_up[r])
                && one_sided(lo_i, up_i)
            {
                // Two opposite one-sided rows would become a new ranged row.
                // TODO(perf): that is HiGHS/SCIP behaviour, but it costs our
                // tree 1.68x nodes (sgm, jobshop_bigM seeds 1-12, mip_perf
                // A/B 2026-10-22) on an unchanged LP relaxation; revisit
                // once node LPs / branching treat ranged rows on par.
                continue;
            }
            if lo > up {
                // A violation of `gap` in r's space is `gap·|s|` in i's.
                if lo - up > self.feastol * (1.0f64).max(1.0 / s.abs()) {
                    return Err(Infeasible);
                }
                let mid = 0.5 * (lo + up);
                self.row_lo[r] = mid;
                self.row_up[r] = mid;
            } else {
                self.row_lo[r] = lo;
                self.row_up[r] = up;
            }
            self.remove_row(i);
            self.stats.parallel_rows += 1;
        }
        Ok(())
    }

    /// Parallel columns with equal scale (HiGHS `detectParallelRowsAndCols`
    /// column part, HPresolve.cpp:7267, restricted to identical columns):
    /// two columns of the same type with identical coefficients and cost
    /// only ever appear as `x_j + x_k`, so they are replaced by one column
    /// `z` with bounds `[l_j + l_k, u_j + u_k]`. For integers every integer
    /// `z` in that range splits into integers within the original bounds,
    /// so the integer feasible set maps onto itself. Removes the symmetry of
    /// duplicated items (knapsack_cover: HiGHS keeps 44% of the
    /// columns).
    fn parallel_cols(&mut self) {
        use std::collections::HashMap;
        let n = self.lo.len();
        let mut vecs: Vec<Vec<(usize, u64)>> = vec![Vec::new(); n];
        for (i, row) in self.rows.iter().enumerate() {
            if self.row_alive[i] {
                for &(j, a) in row {
                    vecs[j].push((i, (a + 0.0).to_bits()));
                }
            }
        }
        type Key = (Vec<(usize, u64)>, u64, bool);
        let mut seen: HashMap<Key, usize> = HashMap::new();
        for (j, mut v) in vecs.into_iter().enumerate() {
            if !self.col_alive[j] || v.is_empty() {
                continue;
            }
            v.sort_unstable();
            let key = (v, (self.cost[j] + 0.0).to_bits(), self.integer[j]);
            let Some(&k) = seen.get(&key) else {
                seen.insert(key, j);
                continue;
            };
            self.stack.push(Step::Merged {
                col: j,
                into: k,
                lo_col: self.lo[j],
                up_col: self.up[j],
                lo_into: self.lo[k],
                up_into: self.up[k],
            });
            self.lo[k] += self.lo[j];
            self.up[k] += self.up[j];
            for i in std::mem::take(&mut self.cols[j]) {
                self.rows[i].retain(|e| e.0 != j);
            }
            self.col_alive[j] = false;
            self.stats.parallel_cols += 1;
            self.stats.cols_removed += 1;
            self.changed = true;
        }
    }

    /// Set-packing conflicts (Balas & Padberg 1976, "Set partitioning: a
    /// survey", SIAM Review 18, row/column reductions; SCIP setppc fixings).
    /// A *packing* row is
    /// `Σ x ≤ 1` and a *covering* row `Σ x ≥ 1` (an equation is both), all
    /// coefficients 1 over nonnegative integer columns. `x_j = 1` zeroes
    /// every column sharing a packing row with `j`; if that empties some
    /// covering row without `j`, then `x_j = 0`. This subsumes the classic
    /// row domination `S ⊂ T` (`Σ_S x = 1`, `Σ_T x ≤ 1` ⇒ `x_{T∖S} = 0`);
    /// the dominating row then duplicates `S` and goes in
    /// [`Self::parallel_rows`].
    fn setppc_conflicts(&mut self) {
        let unit = |w: &Self, i: usize| {
            w.row_alive[i]
                && !w.rows[i].is_empty()
                && w.rows[i]
                    .iter()
                    .all(|&(j, a)| a == 1.0 && w.integer[j] && w.lo[j] >= 0.0)
        };
        let m = self.rows.len();
        let packing: Vec<bool> = (0..m)
            .map(|i| self.row_up[i] == 1.0 && unit(self, i))
            .collect();
        let covering: Vec<usize> = (0..m)
            .filter(|&i| self.row_lo[i] == 1.0 && unit(self, i))
            .collect();
        if covering.is_empty() || !packing.contains(&true) {
            return;
        }
        // stamp[c] == j + 1 ⇔ x_j = 1 forces x_c = 0.
        let mut stamp = vec![0usize; self.lo.len()];
        let mut budget: usize = 4_000_000;
        for j in 0..self.lo.len() {
            if !self.col_alive[j] || !self.integer[j] || self.lo[j] != 0.0 {
                continue;
            }
            let mut any = false;
            for &i in &self.cols[j] {
                if packing[i] && self.row_alive[i] {
                    any = true;
                    budget = budget.saturating_sub(self.rows[i].len());
                    for &(c, _) in &self.rows[i] {
                        if c != j {
                            stamp[c] = j + 1;
                        }
                    }
                }
            }
            if !any {
                continue;
            }
            let emptied = covering.iter().any(|&r| {
                self.row_alive[r]
                    && !self.rows[r].is_empty()
                    && self.rows[r].iter().all(|&(c, _)| stamp[c] == j + 1)
            });
            budget = budget.saturating_sub(covering.len());
            if emptied {
                self.stats.setppc_fixed += 1;
                self.fix_col(j, 0.0);
            }
            if budget == 0 {
                return;
            }
        }
    }

    /// Dominated binary columns (Achterberg et al. 2020 §6.4; Gamrath et
    /// al., "Progress in presolving for mixed integer programming", Math.
    /// Prog. Comp. 7 (2015), Thm. 3 worst-case bounds; HiGHS
    /// `dominatedColumns`, HPresolve.cpp:1045, its coefficient test
    /// `checkDominationNonZero`). For binaries `j ≠ k`, write every
    /// one-sided row as `ãᵀx ≤ b`. `j` dominates `k` when `c_j ≤ c_k`,
    /// `ã_j ≤ ã_k` in every one-sided row and `a_j = a_k` in every two-sided
    /// row. Then some optimum has `x_k = 0` if either
    ///
    /// * (A) `x_j = x_k = 1` violates some row for every value of the other
    ///   columns, or
    /// * (B) `c_k ≥ 0`, `k` is in no two-sided row, and every one-sided row
    ///   that lowering `x_k` can violate (`ã_k < 0`) holds for *every* point
    ///   with `x_j = 1, x_k = 0`.
    ///
    /// Proof: take an optimum with `x_k = 1`. If `x_j = 0`, moving to
    /// `x_j = 1, x_k = 0` changes each row activity by `ã_j − ã_k ≤ 0` (by
    /// `0` in two-sided rows) and the cost by `c_j − c_k ≤ 0`. Otherwise
    /// `x_j = 1`, impossible under (A); under (B) setting `x_k = 0` only
    /// relaxes rows with `ã_k ≥ 0`, keeps the others by the worst-case test,
    /// and changes the cost by `−c_k ≤ 0`. Each fixing is valid for the
    /// current problem (and `k` leaves it), so successive fixings compose.
    fn dominated_columns(&mut self) {
        let binary =
            |w: &Self, j: usize| w.col_alive[j] && w.integer[j] && w.lo[j] == 0.0 && w.up[j] == 1.0;
        // Column vectors (sorted by row) and row activities, built once.
        // Fixing a binary at 0 removes its contribution, so stale minimum
        // activities are too low and stale maxima too high: both tests
        // below stay conservative.
        let n = self.lo.len();
        let mut colvec: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut act = vec![(0.0, 0.0); self.rows.len()];
        for (i, row) in self.rows.iter().enumerate() {
            if self.row_alive[i] {
                for &(j, a) in row {
                    colvec[j].push((i, a));
                }
                act[i] = self.activity(i);
            }
        }
        let mut budget: usize = 4_000_000;
        for k in 0..n {
            if !binary(self, k) {
                continue;
            }
            // Rows where lowering x_k can hurt. The dominating column must be
            // in each (ã_j ≤ ã_k < 0, or a_j = a_k ≠ 0), so candidates come
            // from the shortest one.
            let Some(best) = colvec[k]
                .iter()
                .filter(|&&(i, a)| self.norm(i, a).is_none_or(|v| v < 0.0))
                .map(|e| e.0)
                .min_by_key(|&i| self.rows[i].len())
            else {
                continue; // no down-lock: dual fixing's job
            };
            budget = budget.saturating_sub(self.rows[best].len());
            let cands: Vec<usize> = self.rows[best].iter().map(|e| e.0).collect();
            for j in cands {
                if j == k || !binary(self, j) || self.cost[j] > self.cost[k] + SMALL {
                    continue;
                }
                let work = colvec[j].len() + colvec[k].len();
                if budget < work {
                    return;
                }
                budget -= work;
                if self.dominates(j, k, &colvec[j], &colvec[k], &act) {
                    self.stats.dominated_fixed += 1;
                    self.fix_col(k, 0.0);
                    break;
                }
            }
        }
    }

    /// `ã` for a one-sided row, `None` for a two-sided one.
    fn norm(&self, i: usize, a: f64) -> Option<f64> {
        match (self.row_lo[i].is_finite(), self.row_up[i].is_finite()) {
            (false, true) => Some(a),
            (true, false) => Some(-a),
            _ => None,
        }
    }

    /// The row tests of [`Self::dominated_columns`], merging the two
    /// row-sorted column vectors.
    fn dominates(
        &self,
        j: usize,
        k: usize,
        vj: &[(usize, f64)],
        vk: &[(usize, f64)],
        act: &[(f64, f64)],
    ) -> bool {
        let (mut exclusive, mut worst_ok) = (false, self.cost[k] >= 0.0);
        let (mut p, mut q) = (0, 0);
        while p < vj.len() || q < vk.len() {
            let ij = vj.get(p).map_or(usize::MAX, |e| e.0);
            let ik = vk.get(q).map_or(usize::MAX, |e| e.0);
            let i = ij.min(ik);
            let aj = if ij == i { vj[p].1 } else { 0.0 };
            let ak = if ik == i { vk[q].1 } else { 0.0 };
            p += usize::from(ij == i);
            q += usize::from(ik == i);
            debug_assert!(j != k);
            let (min, max) = act[i];
            match (self.norm(i, aj), self.norm(i, ak)) {
                (Some(nj), Some(nk)) => {
                    if nj > nk + SMALL {
                        return false;
                    }
                    let (minact, maxact, rhs) = if self.row_up[i].is_finite() {
                        (min, max, self.row_up[i])
                    } else {
                        (-max, -min, -self.row_lo[i])
                    };
                    let tol = REDUNDANT_TOL * rhs.abs().max(1.0);
                    // Binary contributions: min(ã, 0) and max(ã, 0).
                    let both = minact - nj.min(0.0) - nk.min(0.0) + nj + nk;
                    exclusive |= both > rhs + self.feastol.max(tol);
                    // +∞ max activity compares greater, as wanted.
                    if nk < 0.0 && maxact - nj.max(0.0) - nk.max(0.0) + nj > rhs + tol {
                        worst_ok = false;
                    }
                }
                _ => {
                    if (aj - ak).abs() > SMALL {
                        return false;
                    }
                    if ak != 0.0 {
                        worst_ok = false;
                    }
                    // Two-sided: x_j = x_k = 1 exceeds the upper side, or
                    // falls below the lower side, for every other value.
                    let lo_both = min - aj.min(0.0) - ak.min(0.0) + aj + ak;
                    let up_both = max - aj.max(0.0) - ak.max(0.0) + aj + ak;
                    exclusive |= lo_both > self.row_up[i] + self.feastol
                        || up_both < self.row_lo[i] - self.feastol;
                }
            }
        }
        exclusive || worst_ok
    }

    // ---------------------------------------------------------------- doubletons

    fn doubleton_pass(&mut self) -> Res {
        for i in 0..self.rows.len() {
            if self.row_alive[i]
                && self.rows[i].len() == 2
                && self.row_lo[i] == self.row_up[i]
                && self.row_lo[i].is_finite()
            {
                self.doubleton_eq(i)?;
            }
        }
        Ok(())
    }

    /// Which entry of a doubleton equation to substitute (HiGHS
    /// `colAtPos1Better`, HPresolve.cpp:3051): continuous before integer; two
    /// integers → smaller |coefficient|, then sparser column; two
    /// continuous → sparser column unless the coefficients differ by more
    /// than 2×, then the larger |coefficient|.
    fn first_is_better(&self, (j1, a1): (usize, f64), (j2, a2): (usize, f64)) -> bool {
        let (n1, n2) = (self.cols[j1].len(), self.cols[j2].len());
        match (self.integer[j1], self.integer[j2]) {
            (true, true) => {
                if a1.abs() < a2.abs() - SMALL {
                    true
                } else if a2.abs() < a1.abs() - SMALL {
                    false
                } else {
                    n1 < n2
                }
            }
            (true, false) => false,
            (false, true) => true,
            (false, false) => {
                if n1 == 1 {
                    true
                } else if n2 == 1 {
                    false
                } else if n1 != n2 && a1.abs().max(a2.abs()) <= 2.0 * a1.abs().min(a2.abs()) {
                    n1 < n2
                } else {
                    a1.abs() > a2.abs()
                }
            }
        }
    }

    fn doubleton_eq(&mut self, i: usize) -> Res {
        let (e1, e2) = (self.rows[i][0], self.rows[i][1]);
        let ((e, ae), (t, mut at)) = if self.first_is_better(e1, e2) {
            (e1, e2)
        } else {
            (e2, e1)
        };
        let mut b = self.row_up[i];
        if self.integer[e] {
            if !self.integer[t] {
                return Ok(()); // continuous x_t would make x_e fractional
            }
            // x_e = b/a_e − (a_t/a_e)·x_t is integral for every integral x_t
            // iff both ratios are integers (HPresolve.cpp:3122).
            let r = (at / ae).round();
            if (r * ae - at).abs() > SMALL {
                return Ok(());
            }
            at = r * ae;
            let rb = (b / ae).round() * ae;
            if (rb - b).abs() > self.feastol {
                return Err(Infeasible);
            }
            b = rb;
        }
        if (at / ae).abs() > MAX_SUBST_RATIO || self.cols[e].len() > MAX_SUBST_COL_LEN {
            return Ok(());
        }
        // Bounds of x_e carried over to x_t = (b − a_e x_e)/a_t.
        let at_lo = |v: f64| {
            if v.is_finite() {
                (b - ae * v) / at
            } else {
                // x_e → ±∞ drives x_t to ∓sign(a_e/a_t)·∞.
                -(ae / at).signum() * v
            }
        };
        let (t1, t2) = (at_lo(self.lo[e]), at_lo(self.up[e]));
        let (tl, tu) = (t1.min(t2), t1.max(t2));
        if tl > f64::NEG_INFINITY {
            self.tighten_lower(t, tl)?;
        }
        if tu < f64::INFINITY {
            self.tighten_upper(t, tu)?;
        }
        self.stack.push(Step::Substituted {
            col: e,
            rhs: b,
            coef: ae,
            other: t,
            other_coef: at,
        });
        self.remove_row(i);
        // Substitute x_e = b/a_e − (a_t/a_e)·x_t into every other row.
        let (shift, slope) = (b / ae, at / ae);
        for r in std::mem::take(&mut self.cols[e]) {
            let k = self.rows[r].iter().position(|en| en.0 == e).unwrap();
            let a = self.rows[r].swap_remove(k).1;
            self.row_lo[r] -= a * shift;
            self.row_up[r] -= a * shift;
            let new = self.coef(r, t) - a * slope;
            self.set_coef(r, t, new);
        }
        self.offset += self.cost[e] * shift;
        self.cost[t] -= self.cost[e] * slope;
        self.cost[e] = 0.0;
        self.col_alive[e] = false;
        self.stats.cols_removed += 1;
        self.stats.doubletons += 1;
        self.changed = true;
        Ok(())
    }

    // ---------------------------------------------------------------- output

    fn finish(self, p: &MipProblem) -> Presolved {
        let col_map: Vec<usize> = (0..p.num_col).filter(|&j| self.col_alive[j]).collect();
        let mut new_index = vec![usize::MAX; p.num_col];
        for (k, &j) in col_map.iter().enumerate() {
            new_index[j] = k;
        }
        let row_map: Vec<usize> = (0..p.num_row).filter(|&i| self.row_alive[i]).collect();
        let mut trip: Vec<Triplet> = Vec::new();
        for (r, &i) in row_map.iter().enumerate() {
            for &(j, a) in &self.rows[i] {
                trip.push((r, new_index[j], a));
            }
        }
        let pick = |v: &[f64]| col_map.iter().map(|&j| v[j]).collect::<Vec<f64>>();
        let problem = MipProblem {
            num_col: col_map.len(),
            num_row: row_map.len(),
            matrix: CscMatrix::from_triplets(row_map.len(), col_map.len(), &trip),
            row_lower: row_map.iter().map(|&i| self.row_lo[i]).collect(),
            row_upper: row_map.iter().map(|&i| self.row_up[i]).collect(),
            cost: pick(&self.cost),
            maximize: p.maximize,
            col_lower: pick(&self.lo),
            col_upper: pick(&self.up),
            integer: col_map.iter().map(|&j| self.integer[j]).collect(),
        };
        Presolved {
            problem,
            offset: self.offset,
            stats: self.stats,
            orig_num_col: p.num_col,
            col_map,
            stack: self.stack,
        }
    }
}

/// `(min, max)` of `a·x` over `x ∈ [lo, up]` (`a ≠ 0`).
fn contribution(a: f64, lo: f64, up: f64) -> (f64, f64) {
    if a > 0.0 {
        (a * lo, a * up)
    } else {
        (a * up, a * lo)
    }
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

/// Smallest denominator `q ≤ MAX_DENOM` with `x·q` integral (to `SMALL`
/// relative), via the convergents of the continued fraction of `|x|`.
fn denominator(x: f64) -> Option<i64> {
    let x = x.abs();
    let integral = |q: i64| {
        let v = x * q as f64;
        (v - v.round()).abs() <= SMALL * v.abs().max(1.0)
    };
    let (mut h_prev, mut k_prev) = (1i64, 0i64);
    let mut a = x.floor();
    let (mut h, mut k) = (a as i64, 1i64);
    let mut frac = x - a;
    loop {
        if integral(k) {
            return Some(k);
        }
        if frac <= 0.0 {
            return None;
        }
        let inv = 1.0 / frac;
        a = inv.floor();
        frac = inv - a;
        let ai = a as i64;
        let (h_next, k_next) = (ai * h + h_prev, ai * k + k_prev);
        if k_next > MAX_DENOM || ai <= 0 {
            return None;
        }
        (h_prev, k_prev, h, k) = (h, k, h_next, k_next);
    }
}

/// Positive `s` such that `s·vals` are coprime integers (returned), or `None`
/// when no common denominator `≤ MAX_SCALE` makes them integral. Same goal as
/// HiGHS `HighsIntegers::integralScale` (HighsIntegers.h:127).
fn integral_scale(vals: &[f64]) -> Option<(f64, Vec<i64>)> {
    let mut lcm: i64 = 1;
    for &v in vals {
        let q = denominator(v)?;
        lcm = lcm / gcd(lcm, q) * q;
        if lcm > MAX_SCALE {
            return None;
        }
    }
    let l = lcm as f64;
    let mut ints = Vec::with_capacity(vals.len());
    let mut g = 0i64;
    for &v in vals {
        let s = v * l;
        let r = s.round();
        if (s - r).abs() > SMALL * s.abs().max(1.0) || r.abs() > 2f64.powi(52) {
            return None;
        }
        let r = r as i64;
        g = gcd(g, r);
        ints.push(r);
    }
    if g == 0 {
        return None;
    }
    for r in &mut ints {
        *r /= g;
    }
    Some((l / g as f64, ints))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bnb::{MipConfig, MipStatus, solve};
    use rusolve_lp::rng::XorShift64;

    const INF: f64 = f64::INFINITY;

    /// Builds a problem from dense rows.
    fn problem(
        a: &[&[f64]],
        row: &[(f64, f64)],
        cost: &[f64],
        bounds: &[(f64, f64)],
        integer: &[bool],
    ) -> MipProblem {
        let (m, n) = (a.len(), cost.len());
        let mut trip = Vec::new();
        for (i, r) in a.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                if v != 0.0 {
                    trip.push((i, j, v));
                }
            }
        }
        MipProblem {
            num_col: n,
            num_row: m,
            matrix: CscMatrix::from_triplets(m, n, &trip),
            row_lower: row.iter().map(|r| r.0).collect(),
            row_upper: row.iter().map(|r| r.1).collect(),
            cost: cost.to_vec(),
            maximize: false,
            col_lower: bounds.iter().map(|b| b.0).collect(),
            col_upper: bounds.iter().map(|b| b.1).collect(),
            integer: integer.to_vec(),
        }
    }

    fn reduced(p: &MipProblem) -> Presolved {
        match presolve(p, 1e-6) {
            Presolve::Reduced(r) => *r,
            Presolve::Infeasible => panic!("unexpected infeasible"),
        }
    }

    /// Presolve without the pair reductions, for tests of single-row rules
    /// whose instances have interchangeable columns.
    fn reduced_single(p: &MipProblem) -> Presolved {
        let mut w = Work::new(p, 1e-6);
        w.pairs = false;
        w.run().expect("unexpected infeasible");
        w.finish(p)
    }

    #[test]
    fn integral_scale_finds_coprime_integers() {
        assert_eq!(integral_scale(&[2.0, 4.0, 6.0]), Some((0.5, vec![1, 2, 3])));
        let (s, ints) = integral_scale(&[0.5, 1.5, -0.25]).unwrap();
        assert_eq!((s, ints), (4.0, vec![2, 6, -1]));
        let (s, ints) = integral_scale(&[1.0 / 3.0, 0.5]).unwrap();
        assert!((s - 6.0).abs() < 1e-12);
        assert_eq!(ints, vec![2, 3]);
        assert_eq!(integral_scale(&[1.0, std::f64::consts::PI]), None);
        let (s, ints) = integral_scale(&[0.1, 0.3]).unwrap();
        assert!((s - 10.0).abs() < 1e-12);
        assert_eq!(ints, vec![1, 3]);
    }

    #[test]
    fn even_equation_with_odd_rhs_is_infeasible() {
        // 2x + 4y + 6z = 7 over unbounded integers.
        let p = problem(
            &[&[2.0, 4.0, 6.0]],
            &[(7.0, 7.0)],
            &[1.0, 1.0, 1.0],
            &[(-INF, INF); 3],
            &[true; 3],
        );
        assert!(matches!(presolve(&p, 1e-6), Presolve::Infeasible));
    }

    #[test]
    fn gcd_rounds_inequality_sides() {
        // 2x + 4y − 2z ≤ 7, 2x − 2y + 4z ≥ −3 ⇒ x + 2y − z ≤ 3,
        // x − y + 2z ≥ −1.
        let p = problem(
            &[&[2.0, 4.0, -2.0], &[2.0, -2.0, 4.0]],
            &[(-INF, 7.0), (-3.0, INF)],
            &[-1.0, -1.0, 1.0],
            &[(-10.0, 10.0); 3],
            &[true; 3],
        );
        let r = reduced(&p);
        assert!(r.stats.rows_scaled >= 2, "{:?}", r.stats);
        let q = &r.problem;
        assert_eq!(q.num_row, 2);
        let row = |i| (0..3).map(|j| q.matrix.get(i, j)).collect::<Vec<_>>();
        assert_eq!(row(0), vec![1.0, 2.0, -1.0]);
        assert_eq!(row(1), vec![1.0, -1.0, 2.0]);
        assert_eq!(q.row_upper[0], 3.0);
        assert_eq!(q.row_lower[1], -1.0);
        check_against_original(&p);
    }

    #[test]
    fn coefficient_strengthening_binary_knapsack() {
        // x + y + 5z ≤ 6 with binaries: max act 7, d = 1 ⇒ 5z → z, rhs 2.
        let p = problem(
            &[&[1.0, 1.0, 5.0]],
            &[(-INF, 6.0)],
            &[-1.0, -1.0, -1.0],
            &[(0.0, 1.0); 3],
            &[true; 3],
        );
        let r = reduced_single(&p);
        assert_eq!(r.stats.coefs_strengthened, 1);
        let q = &r.problem;
        assert_eq!(q.matrix.get(0, 2), 1.0);
        assert_eq!(q.row_upper[0], 2.0);
    }

    #[test]
    fn cover_strengthening_beyond_the_simple_rule() {
        // 3x + 3y + 2z ≤ 4 (binaries): max act 8, d = 4 ⇒ the simple rule
        // does nothing; the cover {x, y} splits off and z's weight drops to
        // 1: 3x + 3y + z ≤ 3.
        let p = problem(
            &[&[3.0, 3.0, 2.0]],
            &[(-INF, 4.0)],
            &[-1.0, -1.0, -1.0],
            &[(0.0, 1.0); 3],
            &[true; 3],
        );
        let r = reduced_single(&p);
        assert_eq!(r.stats.coefs_strengthened, 0);
        assert!(r.stats.cover_strengthened >= 1, "{:?}", r.stats);
        let q = &r.problem;
        assert_eq!(q.matrix.get(0, 2), 1.0);
        assert_eq!(q.row_upper[0], 3.0);
        assert_same_integer_points(&p, &r);
        // The same row with z complemented, written as a ≥ row:
        // −3x − 3y + 2z' ≥ −2.
        let p = problem(
            &[&[-3.0, -3.0, 2.0]],
            &[(-2.0, INF)],
            &[-1.0, -1.0, 1.0],
            &[(0.0, 1.0); 3],
            &[true; 3],
        );
        let r = reduced_single(&p);
        assert!(r.stats.cover_strengthened >= 1, "{:?}", r.stats);
        assert_same_integer_points(&p, &r);
    }

    #[test]
    fn cover_strengthening_drops_coefficients_reduced_to_zero() {
        // Minimised from block_angular s46: four jobs, each assigned to
        // machine A (x) or B (y), with a knapsack capacity per machine.
        // The doubletons substitute y_j = 1 − x_j into B's row, and cover
        // strengthening then reduces one weight to exactly 0. Kept as an
        // explicit 0, it later became a singleton row giving `x/0` bounds,
        // and presolve declared this feasible problem Infeasible.
        let p = problem(
            &[
                &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                &[0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                &[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                &[7.0, 33.0, 68.0, 4.0, 0.0, 0.0, 0.0, 0.0],
                &[0.0, 0.0, 0.0, 0.0, 16.0, 86.0, 33.0, 73.0],
            ],
            &[
                (1.0, 1.0),
                (1.0, 1.0),
                (1.0, 1.0),
                (1.0, 1.0),
                (-INF, 44.0),
                (-INF, 122.0),
            ],
            &[101.0, 76.0, 50.0, 115.0, 95.0, 23.0, 78.0, 32.0],
            &[(0.0, 1.0); 8],
            &[true; 8],
        );
        let r = reduced(&p);
        assert!(
            r.problem.matrix.values().iter().all(|v| v.abs() > SMALL),
            "zero coefficient left in the reduced matrix"
        );
        assert_same_integer_points(&p, &r);
        check_against_original(&p);
    }

    #[test]
    fn implied_flow_bound_tightens_big_m_row() {
        // f = g + h (conservation), g ≤ 3, h ≤ 4 ⇒ f ≤ 7; f − 27y ≤ 0 with y
        // binary becomes f − 7y ≤ 0 (vub32).
        let p = problem(
            &[&[1.0, -1.0, -1.0, 0.0], &[1.0, 0.0, 0.0, -27.0]],
            &[(0.0, 0.0), (-INF, 0.0)],
            &[-1.0, 0.5, 0.5, 10.0],
            &[(0.0, INF), (0.0, 3.0), (0.0, 4.0), (0.0, 1.0)],
            &[false, false, false, true],
        );
        let r = reduced(&p);
        let q = &r.problem;
        let (kf, ky) = (
            r.col_map().iter().position(|&j| j == 0).unwrap(),
            r.col_map().iter().position(|&j| j == 3).unwrap(),
        );
        assert_eq!(q.col_upper[kf], 7.0);
        let vub = (0..q.num_row)
            .find(|&i| q.matrix.get(i, ky) != 0.0)
            .unwrap();
        assert_eq!(q.matrix.get(vub, ky), -7.0, "{:?}", r.stats);
        assert_eq!(q.row_upper[vub], 0.0);
        check_against_original(&p);
    }

    #[test]
    fn parallel_rows_merge_unless_they_form_a_new_ranged_row() {
        // dA: Sa − Sb + 10z ≥ 3; its scaled negation −Sa + Sb − 10z ≤ −5
        // (dA ≥ 5) merges into dA. dB: Sb − Sa − 10z ≥ −6 (dA ≤ 6) would
        // make a new ranged row and stays separate (TODO(perf) in
        // parallel_rows).
        let p = problem(
            &[&[1.0, -1.0, 10.0], &[-2.0, 2.0, -20.0], &[-1.0, 1.0, -10.0]],
            &[(3.0, INF), (-INF, -10.0), (-6.0, INF)],
            &[1.0, 1.0, 0.0],
            &[(0.0, INF), (0.0, INF), (0.0, 1.0)],
            &[false, false, true],
        );
        let r = reduced(&p);
        assert_eq!(r.stats.parallel_rows, 1, "{:?}", r.stats);
        let q = &r.problem;
        assert_eq!(q.num_row, 2);
        assert_eq!((q.row_lower[0], q.row_upper[0]), (5.0, INF));
        check_against_original(&p);
    }

    #[test]
    fn identical_columns_merge_and_split_back() {
        // 3x + 3y + 5z ≥ 7 with x, y ∈ {0..3} identical: z' = x + y ∈
        // [0, 6]. (General integers: binaries would go to domination.)
        let p = problem(
            &[&[3.0, 3.0, 5.0]],
            &[(7.0, INF)],
            &[2.0, 2.0, 3.0],
            &[(0.0, 3.0), (0.0, 3.0), (0.0, 1.0)],
            &[true; 3],
        );
        let r = reduced(&p);
        assert_eq!(r.stats.parallel_cols, 1, "{:?}", r.stats);
        assert_eq!(r.problem.num_col, 2);
        let k = r.col_map().iter().position(|&j| j == 0).unwrap();
        assert_eq!(r.problem.col_upper[k], 6.0);
        // z' = 4 splits into integers in range summing to 4.
        let mut x = vec![0.0; 2];
        x[k] = 4.0;
        let y = r.postsolve(&x);
        assert_eq!(y[0] + y[1], 4.0);
        assert!(
            y[..2]
                .iter()
                .all(|&v| (0.0..=3.0).contains(&v) && v == v.round())
        );
        assert_same_integer_points(&p, &r);
        check_against_original(&p);
    }

    #[test]
    fn forward_map_carries_a_point_into_the_reduced_space() {
        // The merged-columns model plus a singleton row fixing w = 2: the
        // forward map must fold the merge (x + y lands on the survivor) and
        // skip the fixed column, so postsolve(forward_map(x)) lands back on
        // an equivalent integer point.
        let p = problem(
            &[&[3.0, 3.0, 5.0, 1.0], &[0.0, 0.0, 0.0, 1.0]],
            &[(7.0, INF), (2.0, 2.0)],
            &[2.0, 2.0, 3.0, 1.0],
            &[(0.0, 3.0), (0.0, 3.0), (0.0, 1.0), (0.0, 5.0)],
            &[true; 4],
        );
        let r = reduced(&p);
        assert!(r.stats.parallel_cols >= 1, "{:?}", r.stats);
        for x in [
            vec![1.0, 3.0, 0.0, 2.0],
            vec![3.0, 3.0, 1.0, 2.0],
            vec![0.0, 0.0, 1.0, 2.0],
        ] {
            assert!(
                p.violation(&x).row.max(p.violation(&x).bound) <= 1e-9,
                "{x:?}"
            );
            let z = r.forward_map(&x);
            assert_eq!(z.len(), r.problem.num_col);
            assert!(
                r.problem
                    .violation(&z)
                    .row
                    .max(r.problem.violation(&z).bound)
                    <= 1e-9,
                "{x:?} -> {z:?}"
            );
            let back = r.postsolve(&z);
            assert!(
                p.violation(&back).row.max(p.violation(&back).bound) <= 1e-9,
                "{x:?} -> {z:?} -> {back:?}"
            );
            assert!((p.objective(&back) - p.objective(&x)).abs() <= 1e-9);
        }
    }

    #[test]
    fn dominated_set_cover_columns_are_fixed() {
        // x0 covers {r0, r1} ⊇ x1's {r0} at lower cost; x2 covers
        // {r1, r2} ⊇ x3's {r2} at equal cost.
        let p = problem(
            &[
                &[1.0, 1.0, 0.0, 0.0],
                &[1.0, 0.0, 1.0, 0.0],
                &[0.0, 0.0, 1.0, 1.0],
            ],
            &[(1.0, INF); 3],
            &[2.0, 3.0, 2.0, 2.0],
            &[(0.0, 1.0); 4],
            &[true; 4],
        );
        let r = reduced(&p);
        assert!(r.stats.dominated_fixed >= 2, "{:?}", r.stats);
        assert_same_integer_points(&p, &r);
        check_against_original(&p);
    }

    #[test]
    fn setppc_row_domination_fixes_extra_columns() {
        // x0 + x1 + x2 = 1 and x0 + x1 + x2 + x3 ≤ 1 ⇒ x3 = 0, although
        // x3 is the most attractive column.
        let p = problem(
            &[&[1.0, 1.0, 1.0, 0.0], &[1.0, 1.0, 1.0, 1.0]],
            &[(1.0, 1.0), (-INF, 1.0)],
            &[1.0, 2.0, 3.0, -5.0],
            &[(0.0, 1.0); 4],
            &[true; 4],
        );
        let r = reduced(&p);
        assert!(r.stats.setppc_fixed >= 1, "{:?}", r.stats);
        assert!(!r.col_map().contains(&3));
        assert_same_integer_points(&p, &r);
        check_against_original(&p);
    }

    #[test]
    fn coefficient_strengthening_ge_row_negative_coef() {
        // −4x + y + z ≥ −3, x ∈ {0,1,2}, y,z ∈ {0,1}: min act −8; at x = 2
        // only −8 + 2 ≥ −3 fails ⇒ x ≤ ... The reduction must keep the
        // integer feasible set; compare by enumeration.
        let p = problem(
            &[&[-4.0, 1.0, 1.0]],
            &[(-3.0, INF)],
            &[-1.0, 0.5, 0.25],
            &[(0.0, 2.0), (0.0, 1.0), (0.0, 1.0)],
            &[true; 3],
        );
        let r = reduced(&p);
        assert_same_integer_points(&p, &r);
    }

    #[test]
    fn singleton_rows_become_bounds_and_fixed_columns_vanish() {
        // 2x ≤ 5 (x int) ⇒ x ≤ 2; y fixed at 3 enters row 1's sides.
        let p = problem(
            &[&[2.0, 0.0, 0.0], &[1.0, 1.0, 1.0]],
            &[(-INF, 5.0), (4.0, 9.0)],
            &[-1.0, 2.0, 1.0],
            &[(0.0, 10.0), (3.0, 3.0), (0.0, 10.0)],
            &[true, false, false],
        );
        let r = reduced(&p);
        assert_eq!(r.offset, 6.0);
        let x = r.postsolve(&vec![0.0; r.problem.num_col]);
        assert_eq!(x[1], 3.0);
        let k = r.col_map().iter().position(|&j| j == 0);
        if let Some(k) = k {
            assert_eq!(r.problem.col_upper[k], 2.0);
        }
        let sol = solve(&r.problem, &MipConfig::default()).unwrap();
        let full = solve(
            &p,
            &MipConfig {
                presolve: false,
                ..MipConfig::default()
            },
        )
        .unwrap();
        let x = r.postsolve(&sol.x);
        assert!((p.objective(&x) - full.objective).abs() < 1e-9);
    }

    #[test]
    fn root_bound_stats_include_the_presolve_offset() {
        // A 3-row multi-knapsack over 15 binaries plus y fixed at 3 with
        // cost 2, which leaves offset 6 in presolve. Presolve only tightens,
        // so its root bound may not fall below the unpresolved run's;
        // without the offset it would sit up to 6 lower.
        let n = 15;
        let a: Vec<Vec<f64>> = (0..3)
            .map(|r| {
                let mut row: Vec<f64> = (0..n).map(|i| ((i * 7 + r * 3) % 11 + 1) as f64).collect();
                row.push(0.0);
                row
            })
            .collect();
        let rows: Vec<&[f64]> = a.iter().map(|r| r.as_slice()).collect();
        let sides: Vec<(f64, f64)> = a
            .iter()
            .map(|r| (-INF, (r.iter().sum::<f64>() / 2.0).floor() + 0.5))
            .collect();
        let mut cost: Vec<f64> = (0..n).map(|i| -(((i * 5) % 9 + 1) as f64)).collect();
        cost.push(2.0);
        let mut bounds = vec![(0.0, 1.0); n];
        bounds.push((3.0, 3.0));
        let mut integer = vec![true; n];
        integer.push(false);
        let p = problem(&rows, &sides, &cost, &bounds, &integer);
        let on = solve_presolved(&p, &MipConfig::default()).unwrap();
        let off = solve(
            &p,
            &MipConfig {
                presolve: false,
                ..MipConfig::default()
            },
        )
        .unwrap();
        assert!(on.cuts.rounds > 0 && off.cuts.rounds > 0);
        assert!(on.cuts.root_bound_before >= off.cuts.root_bound_before - 1e-9);
        assert!(on.cuts.root_bound_after <= on.objective + 1e-9);
    }

    #[test]
    fn doubleton_equation_substitutes_continuous_column() {
        // x int, y cont: x + 2y = 3, y ∈ [0, 1] ⇒ x ∈ {1, 3} ... x ∈ [1,3];
        // second row x + y ≤ 2.5.
        let p = problem(
            &[&[1.0, 2.0], &[1.0, 1.0]],
            &[(3.0, 3.0), (-INF, 2.5)],
            &[-1.0, -3.0],
            &[(0.0, 5.0), (0.0, 1.0)],
            &[true, false],
        );
        let r = reduced(&p);
        assert!(r.stats.doubletons >= 1, "{:?}", r.stats);
        assert!(!r.col_map().contains(&1));
        check_against_original(&p);
    }

    #[test]
    fn doubleton_with_two_integers_requires_integral_ratio() {
        // 2x + 3y = 6: ratio 1.5 ⇒ no integer-safe substitution.
        let p = problem(
            &[&[2.0, 3.0], &[1.0, 1.0]],
            &[(6.0, 6.0), (-INF, 5.0)],
            &[1.0, 1.0],
            &[(-5.0, 5.0), (-5.0, 5.0)],
            &[true, true],
        );
        let r = reduced(&p);
        assert_eq!(r.stats.doubletons, 0);
        check_against_original(&p);
        // 2x + 4y = 6 ⇒ x + 2y = 3 ⇒ x = 3 − 2y.
        let p = problem(
            &[&[2.0, 4.0], &[1.0, 1.0]],
            &[(6.0, 6.0), (-INF, 5.0)],
            &[1.0, -1.0],
            &[(-5.0, 5.0), (-5.0, 5.0)],
            &[true, true],
        );
        let r = reduced(&p);
        assert_eq!(r.stats.doubletons, 1);
        check_against_original(&p);
    }

    #[test]
    fn dual_fixing_and_empty_columns() {
        // y appears only in a ≥ row with positive coefficient and has negative
        // cost ⇒ fixed at its upper bound; z is empty with positive cost.
        let p = problem(
            &[&[1.0, 1.0, 0.0]],
            &[(1.0, INF)],
            &[1.0, -1.0, 2.0],
            &[(0.0, 4.0), (0.0, 3.0), (-1.0, 5.0)],
            &[true, true, false],
        );
        let r = reduced(&p);
        assert!(r.stats.dual_fixed >= 2, "{:?}", r.stats);
        check_against_original(&p);
    }

    #[test]
    fn all_columns_removed_and_empty_problem_solves() {
        let p = problem(
            &[&[1.0, 0.0]],
            &[(1.0, 1.0)],
            &[3.0, -1.0],
            &[(0.0, 2.0), (0.0, 1.0)],
            &[true, false],
        );
        let r = reduced(&p);
        assert_eq!(r.problem.num_col, 0);
        let sol = solve(&p, &MipConfig::default()).unwrap();
        assert_eq!(sol.status, MipStatus::Optimal);
        assert!((sol.objective - 2.0).abs() < 1e-12);
        assert_eq!(sol.x, vec![1.0, 1.0]);
    }

    /// Every integer point in the (finite, small) box: feasible in the
    /// original iff its reduced image is feasible in the reduced problem —
    /// checked here by enumerating reduced points, postsolving and verifying,
    /// and by comparing the optima.
    fn assert_same_integer_points(p: &MipProblem, r: &Presolved) {
        assert_eq!(enumerate(p), enumerate_reduced(r, p));
    }

    fn enumerate(p: &MipProblem) -> Option<f64> {
        if (0..p.num_col).any(|j| p.col_lower[j] > p.col_upper[j]) {
            return None;
        }
        let mut x: Vec<f64> = p.col_lower.clone();
        let mut best: Option<f64> = None;
        loop {
            if p.violation(&x).row <= 1e-9 {
                let z = p.objective(&x);
                best = Some(best.map_or(z, |b: f64| b.min(z)));
            }
            let mut j = 0;
            loop {
                if j == p.num_col {
                    return best;
                }
                if x[j] < p.col_upper[j] {
                    x[j] += 1.0;
                    break;
                }
                x[j] = p.col_lower[j];
                j += 1;
            }
        }
    }

    /// Optimum over reduced integer points whose postsolved image is
    /// original-feasible (and every reduced-feasible point must be).
    fn enumerate_reduced(r: &Presolved, p: &MipProblem) -> Option<f64> {
        let q = &r.problem;
        let mut best = enumerate_with(q, &mut |x| {
            let y = r.postsolve(x);
            let v = p.violation(&y);
            assert!(
                v.row <= 1e-9 && v.bound <= 1e-9 && v.integrality <= 1e-9,
                "postsolved point infeasible: {v:?}"
            );
            let z = p.objective(&y);
            assert!((z - (q.objective(x) + r.offset)).abs() <= 1e-9);
            z
        });
        if q.num_col == 0 {
            best = Some(p.objective(&r.postsolve(&[])));
        }
        best
    }

    fn enumerate_with(q: &MipProblem, f: &mut dyn FnMut(&[f64]) -> f64) -> Option<f64> {
        if (0..q.num_col).any(|j| q.col_lower[j] > q.col_upper[j]) {
            return None;
        }
        let mut x: Vec<f64> = q.col_lower.clone();
        let mut best: Option<f64> = None;
        loop {
            if q.violation(&x).row <= 1e-9 {
                let z = f(&x);
                best = Some(best.map_or(z, |b: f64| b.min(z)));
            }
            let mut j = 0;
            loop {
                if j == q.num_col {
                    return best;
                }
                if x[j] < q.col_upper[j] {
                    x[j] += 1.0;
                    break;
                }
                x[j] = q.col_lower[j];
                j += 1;
            }
        }
    }

    /// Solves with and without presolve and compares status / optimum; the
    /// presolved answer is checked against the original rows and bounds.
    fn check_against_original(p: &MipProblem) {
        let on = solve(p, &MipConfig::default()).unwrap();
        let off = solve(
            p,
            &MipConfig {
                presolve: false,
                ..MipConfig::default()
            },
        )
        .unwrap();
        assert_eq!(on.status, off.status);
        if on.status == MipStatus::Optimal {
            assert!(
                (on.objective - off.objective).abs() <= 1e-6 * off.objective.abs().max(1.0),
                "{} vs {}",
                on.objective,
                off.objective
            );
            let v = p.violation(&on.x);
            assert!(v.row <= 1e-6 && v.bound <= 1e-6 && v.integrality <= 1e-6);
            assert!((p.caller_objective(p.objective(&on.x)) - on.objective).abs() <= 1e-9);
        }
    }

    /// Random small MIPs shaped to trigger every reduction: even
    /// coefficients, singleton and doubleton rows, fixed and empty columns,
    /// ranged rows, mixed integrality.
    fn random_mip(rng: &mut XorShift64, n: usize, m: usize, pure: bool) -> MipProblem {
        let pick =
            |rng: &mut XorShift64, xs: &[f64]| xs[(rng.next_u64() % xs.len() as u64) as usize];
        let mut trip = Vec::new();
        for i in 0..m {
            let shape = rng.next_u64() % 4;
            let even = rng.next_u64().is_multiple_of(3);
            for j in 0..n {
                let keep = match shape {
                    0 => j == (i % n),                       // singleton
                    1 => j == (i % n) || j == ((i + 1) % n), // doubleton
                    _ => !rng.next_u64().is_multiple_of(3),
                };
                if !keep {
                    continue;
                }
                let mut a = pick(rng, &[1.0, 2.0, -1.0, 3.0, -2.0, 0.5, 5.0, -4.0]);
                if even {
                    a = 2.0 * a.round().max(1.0) * a.signum();
                }
                trip.push((i, j, a));
            }
        }
        let (mut row_lower, mut row_upper) = (Vec::new(), Vec::new());
        for _ in 0..m {
            let b = pick(rng, &[-2.0, 0.0, 1.0, 3.0, 4.5, 6.0, 7.0]);
            let (lo, hi) = match rng.next_u64() % 4 {
                0 => (-INF, b),
                1 => (b, INF),
                2 => (b, b),
                _ => (b - 2.0, b + 1.5),
            };
            row_lower.push(lo);
            row_upper.push(hi);
        }
        let col_lower: Vec<f64> = (0..n).map(|_| pick(rng, &[-2.0, -1.0, 0.0, 0.0])).collect();
        let col_upper = col_lower
            .iter()
            .map(|&l| l + pick(rng, &[0.0, 1.0, 2.0, 3.0, 4.0]))
            .collect();
        let integer = (0..n)
            .map(|_| pure || !rng.next_u64().is_multiple_of(3))
            .collect();
        MipProblem {
            num_col: n,
            num_row: m,
            matrix: CscMatrix::from_triplets(m, n, &trip),
            row_lower,
            row_upper,
            cost: (0..n)
                .map(|_| pick(rng, &[-3.0, -1.0, 0.0, 1.0, 2.0, 1.5]))
                .collect(),
            maximize: rng.next_u64().is_multiple_of(4),
            col_lower,
            col_upper,
            integer,
        }
    }

    #[test]
    fn presolve_preserves_integer_optimum_by_enumeration() {
        let mut rng = XorShift64::new(6_2026_1022);
        let (mut feasible, mut infeasible, mut reduced_any) = (0, 0, 0);
        let mut total = PresolveStats::default();
        let mut presolve_infeasible = 0;
        for case in 0..3000 {
            let n = 2 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 5) as usize;
            let p = random_mip(&mut rng, n, m, true);
            let truth = enumerate(&p);
            match presolve(&p, 1e-6) {
                Presolve::Infeasible => {
                    assert_eq!(truth, None, "case {case}: presolve said infeasible");
                    infeasible += 1;
                    presolve_infeasible += 1;
                }
                Presolve::Reduced(r) => {
                    if r.problem.num_col < p.num_col || r.problem.num_row < p.num_row {
                        reduced_any += 1;
                    }
                    let s = r.stats;
                    total.dual_fixed += s.dual_fixed;
                    total.cols_fixed += s.cols_fixed;
                    total.rows_scaled += s.rows_scaled;
                    total.coefs_strengthened += s.coefs_strengthened;
                    total.cover_strengthened += s.cover_strengthened;
                    total.doubletons += s.doubletons;
                    total.sides_dropped += s.sides_dropped;
                    total.bounds_tightened += s.bounds_tightened;
                    let got = enumerate_reduced(&r, &p);
                    match (truth, got) {
                        (None, None) => infeasible += 1,
                        (Some(a), Some(b)) => {
                            feasible += 1;
                            assert!((a - b).abs() <= 1e-9, "case {case}: {a} vs {b}");
                        }
                        _ => panic!("case {case}: {truth:?} vs {got:?} {:?}", r.stats),
                    }
                }
            }
            // End to end through B&B as well.
            let sol = solve(&p, &MipConfig::default()).unwrap();
            match truth {
                Some(z) => {
                    assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                    assert!((p.caller_objective(sol.objective) - z).abs() <= 1e-6);
                    let v = p.violation(&sol.x);
                    assert!(v.row <= 1e-6 && v.bound <= 1e-6 && v.integrality <= 1e-6);
                }
                None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
            }
        }
        assert!(
            feasible > 500 && infeasible > 250 && reduced_any > 750,
            "{feasible}/{infeasible}/{reduced_any}"
        );
        // Every reduction must actually be exercised by this property test.
        assert!(presolve_infeasible > 100, "{presolve_infeasible}");
        for (name, n) in [
            ("dual_fixed", total.dual_fixed),
            ("cols_fixed", total.cols_fixed),
            ("rows_scaled", total.rows_scaled),
            ("coefs_strengthened", total.coefs_strengthened),
            ("doubletons", total.doubletons),
            ("sides_dropped", total.sides_dropped),
            ("bounds_tightened", total.bounds_tightened),
        ] {
            assert!(n > 40, "{name} fired only {n} times: {total:?}");
        }
    }

    /// A presolve stopped by its deadline (already passed: it polls
    /// after the first round of passes) applies fewer reductions, and what
    /// it applied is still a valid presolve — same integer optimum.
    #[test]
    fn presolve_stopped_by_its_deadline_stays_valid() {
        let mut rng = XorShift64::new(48_2026_1029);
        let mut cut_short = 0;
        for case in 0..600 {
            let n = 2 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 5) as usize;
            let p = random_mip(&mut rng, n, m, true);
            let truth = enumerate(&p);
            let full = match presolve(&p, 1e-6) {
                Presolve::Reduced(r) => r.stats.rounds,
                Presolve::Infeasible => continue,
            };
            match presolve_until(&p, 1e-6, Some(Instant::now())) {
                Presolve::Infeasible => assert_eq!(truth, None, "case {case}"),
                Presolve::Reduced(r) => {
                    assert_eq!(r.stats.rounds, 1, "case {case}");
                    cut_short += usize::from(full > 1);
                    let got = enumerate_reduced(&r, &p);
                    match (truth, got) {
                        (None, None) => {}
                        (Some(a), Some(b)) => assert!((a - b).abs() <= 1e-9, "case {case}"),
                        _ => panic!("case {case}: {truth:?} vs {got:?}"),
                    }
                }
            }
        }
        assert!(cut_short > 50, "{cut_short}");
    }

    /// Random binary / small-integer knapsacks: the shape where covers
    /// exist. Integer points must be preserved exactly.
    #[test]
    fn cover_strengthening_preserves_knapsack_points() {
        let mut rng = XorShift64::new(7_2026_1022);
        let mut total = PresolveStats::default();
        for case in 0..800 {
            let n = 3 + (rng.next_u64() % 5) as usize;
            let m = 1 + (rng.next_u64() % 2) as usize;
            let mut trip = Vec::new();
            let mut row_lower = Vec::new();
            let mut row_upper = Vec::new();
            let col_upper: Vec<f64> = (0..n)
                .map(|_| {
                    if rng.next_u64().is_multiple_of(4) {
                        2.0
                    } else {
                        1.0
                    }
                })
                .collect();
            for i in 0..m {
                let mut maxact = 0.0;
                for (j, &u) in col_upper.iter().enumerate() {
                    let mut a = 1.0 + (rng.next_u64() % 9) as f64;
                    if rng.next_u64().is_multiple_of(5) {
                        a = -a;
                    }
                    if a > 0.0 {
                        maxact += a * u;
                    }
                    trip.push((i, j, a));
                }
                let b = (maxact * (0.3 + 0.4 * rng.next_f64())).floor();
                if rng.next_u64().is_multiple_of(2) {
                    row_lower.push(-INF);
                    row_upper.push(b);
                } else {
                    // Same row negated as a ≥ row.
                    for t in trip.iter_mut().filter(|t| t.0 == i) {
                        t.2 = -t.2;
                    }
                    row_lower.push(-b);
                    row_upper.push(INF);
                }
            }
            let p = MipProblem {
                num_col: n,
                num_row: m,
                matrix: CscMatrix::from_triplets(m, n, &trip),
                row_lower,
                row_upper,
                cost: (0..n).map(|_| -1.0 - (rng.next_u64() % 7) as f64).collect(),
                maximize: false,
                col_lower: vec![0.0; n],
                col_upper,
                integer: vec![true; n],
            };
            let truth = enumerate(&p);
            match presolve(&p, 1e-6) {
                Presolve::Infeasible => assert_eq!(truth, None, "case {case}"),
                Presolve::Reduced(r) => {
                    total.cover_strengthened += r.stats.cover_strengthened;
                    total.coefs_strengthened += r.stats.coefs_strengthened;
                    let got = enumerate_reduced(&r, &p);
                    match (truth, got) {
                        (None, None) => {}
                        (Some(a), Some(b)) => assert!((a - b).abs() <= 1e-9, "case {case}"),
                        _ => panic!("case {case}: {truth:?} vs {got:?}"),
                    }
                }
            }
        }
        assert!(
            total.cover_strengthened > 100 && total.coefs_strengthened > 100,
            "{total:?}"
        );
    }

    /// Random set covering / packing / partitioning structure with
    /// duplicated columns, parallel (scaled, negated) rows and small knapsack
    /// rows: the shapes the pair reductions (parallel rows and columns,
    /// dominated columns, setppc conflicts) act on.
    fn random_setppc(rng: &mut XorShift64, n: usize, m: usize) -> MipProblem {
        let pick =
            |rng: &mut XorShift64, xs: &[f64]| xs[(rng.next_u64() % xs.len() as u64) as usize];
        let mut supports: Vec<Vec<usize>> = Vec::new();
        let mut cost: Vec<f64> = Vec::new();
        for j in 0..n {
            if j > 0 && rng.next_u64().is_multiple_of(3) {
                let k = (rng.next_u64() % j as u64) as usize;
                supports.push(supports[k].clone());
                let c = if rng.next_u64().is_multiple_of(2) {
                    cost[k]
                } else {
                    pick(rng, &[-1.0, 0.0, 1.0, 2.0, 3.0])
                };
                cost.push(c);
            } else {
                supports.push(
                    (0..m)
                        .filter(|_| rng.next_u64().is_multiple_of(2))
                        .collect(),
                );
                cost.push(pick(rng, &[-1.0, 0.0, 1.0, 2.0, 3.0]));
            }
        }
        let mut trip = Vec::new();
        let (mut row_lower, mut row_upper) = (Vec::new(), Vec::new());
        let mut rows = 0;
        for i in 0..m {
            let knap = rng.next_u64().is_multiple_of(5);
            let coefs: Vec<(usize, f64)> = (0..n)
                .filter(|&j| supports[j].contains(&i))
                .map(|j| {
                    (
                        j,
                        if knap {
                            pick(rng, &[1.0, 2.0, 3.0])
                        } else {
                            1.0
                        },
                    )
                })
                .collect();
            let (lo, up) = if knap {
                (-INF, pick(rng, &[2.0, 3.0, 4.0]))
            } else {
                match rng.next_u64() % 3 {
                    0 => (1.0, INF),
                    1 => (-INF, 1.0),
                    _ => (1.0, 1.0),
                }
            };
            let negate = rng.next_u64().is_multiple_of(5);
            let copies = if rng.next_u64().is_multiple_of(4) {
                2
            } else {
                1
            };
            for c in 0..copies {
                let mut s = if negate { -1.0 } else { 1.0 };
                let (mut lo, mut up) = (lo, up);
                if c == 1 {
                    s *= pick(rng, &[-2.0, 0.5, 3.0, -1.0]);
                    // A looser or tighter copy of the same row.
                    lo -= pick(rng, &[0.0, 1.0]);
                    up += pick(rng, &[0.0, 1.0, -1.0]);
                }
                for &(j, a) in &coefs {
                    trip.push((rows, j, s * a));
                }
                let (l, u) = if s > 0.0 {
                    (s * lo, s * up)
                } else {
                    (s * up, s * lo)
                };
                row_lower.push(l);
                row_upper.push(u);
                rows += 1;
            }
        }
        MipProblem {
            num_col: n,
            num_row: rows,
            matrix: CscMatrix::from_triplets(rows, n, &trip),
            row_lower,
            row_upper,
            cost,
            maximize: rng.next_u64().is_multiple_of(5),
            col_lower: vec![0.0; n],
            col_upper: (0..n)
                .map(|_| {
                    if rng.next_u64().is_multiple_of(6) {
                        2.0
                    } else {
                        1.0
                    }
                })
                .collect(),
            integer: vec![true; n],
        }
    }

    #[test]
    fn pair_reductions_preserve_integer_optimum_by_enumeration() {
        let mut rng = XorShift64::new(8_2026_1022);
        let mut total = PresolveStats::default();
        let (mut feasible, mut infeasible) = (0, 0);
        for case in 0..2000 {
            let n = 4 + (rng.next_u64() % 6) as usize;
            let m = 2 + (rng.next_u64() % 4) as usize;
            let p = random_setppc(&mut rng, n, m);
            let truth = enumerate(&p);
            match presolve(&p, 1e-6) {
                Presolve::Infeasible => assert_eq!(truth, None, "case {case}"),
                Presolve::Reduced(r) => {
                    let s = r.stats;
                    total.parallel_rows += s.parallel_rows;
                    total.parallel_cols += s.parallel_cols;
                    total.dominated_fixed += s.dominated_fixed;
                    total.setppc_fixed += s.setppc_fixed;
                    let got = enumerate_reduced(&r, &p);
                    match (truth, got) {
                        (None, None) => {}
                        (Some(a), Some(b)) => {
                            assert!((a - b).abs() <= 1e-9, "case {case}: {a} vs {b} {s:?}")
                        }
                        _ => panic!("case {case}: {truth:?} vs {got:?} {s:?}"),
                    }
                }
            }
            if truth.is_some() {
                feasible += 1;
            } else {
                infeasible += 1;
            }
            if case % 4 == 0 {
                let sol = solve(&p, &MipConfig::default()).unwrap();
                match truth {
                    Some(z) => {
                        assert_eq!(sol.status, MipStatus::Optimal, "case {case}");
                        assert!((p.caller_objective(sol.objective) - z).abs() <= 1e-6);
                    }
                    None => assert_eq!(sol.status, MipStatus::Infeasible, "case {case}"),
                }
            }
        }
        assert!(
            feasible > 300 && infeasible > 100,
            "{feasible}/{infeasible}"
        );
        for (name, n) in [
            ("parallel_rows", total.parallel_rows),
            ("parallel_cols", total.parallel_cols),
            ("dominated_fixed", total.dominated_fixed),
            ("setppc_fixed", total.setppc_fixed),
        ] {
            assert!(n > 40, "{name} fired only {n} times: {total:?}");
        }
    }

    #[test]
    fn presolve_matches_plain_bnb_on_random_mixed_integer() {
        let mut rng = XorShift64::new(20261022);
        for _ in 0..300 {
            let n = 2 + (rng.next_u64() % 6) as usize;
            let m = 1 + (rng.next_u64() % 6) as usize;
            check_against_original(&random_mip(&mut rng, n, m, false));
        }
    }
}
