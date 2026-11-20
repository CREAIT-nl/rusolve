//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: conflict_resolution.c (conflictRowCopy 98-136, conflictRowReplace 140-180,
//!   conflictRowClear 1257-1287, conflictRowRemoveZeroVar 782-799, conflictRowRemoveZeroVars 803-816,
//!   computeSlack 1291-1368); def.h 162-193; set.c 113-186, 6824, 6906, 6946, 7025, 7037, 7113, 7206.

pub mod mir;
pub mod output;
pub mod queue;
pub mod reason;
pub mod reduce;
pub mod resolve;
pub mod row;
pub mod step;

use std::collections::BinaryHeap;

use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;

// ---------------------------------------------------------------------------
// D7 constants
// ---------------------------------------------------------------------------

/// SCIP set.c:179
pub const S2_MAXVARSFRACRES: f64 = 0.20;
/// SCIP set.c:113
pub const S2_MINMAXVARS: HighsInt = 0;
/// SCIP set.c:180
pub const S2_RESFUIPLEVELS: HighsInt = 1;
/// SCIP set.c:182
pub const S2_FIXANDCONTINUE: bool = false;
/// SCIP set.c:183
pub const S2_MAXCOEFQUOT: f64 = 1e6;
/// SCIP set.c:186
pub const S2_MBREDUCTION: bool = true;
/// Not SCIP, Q6
pub const S2_RESOLVE_LP_PROOFS: bool = true;

// MIR settings (conflict_resolution.c:87-94)

/// conflict_resolution.c:87
pub const S2_BOUNDSWITCH: f64 = 0.51;
/// conflict_resolution.c:88
pub const S2_POSTPROCESS: bool = false;
/// conflict_resolution.c:89
pub const S2_USEVBDS: bool = false;
/// conflict_resolution.c:90
pub const S2_FIXINTEGRALRHS: bool = false;
/// conflict_resolution.c:91
pub const S2_MINFRAC: f64 = 0.05;
/// conflict_resolution.c:92
pub const S2_MAXFRAC: f64 = 0.999;
/// cuts.c:4090
pub const S2_MAXCMIRSCALE: f64 = 1e6;
/// set.c:460-461, norm 'e' at set.c:463
pub const S2_MINEFFICACY: f64 = 1e-4;

// Tolerances (def.h)

/// def.h:163
pub const S2_EPSILON: f64 = 1e-9;
/// def.h:164
pub const S2_SUMEPSILON: f64 = 1e-6;
/// def.h:165
pub const S2_FEASTOL: f64 = 1e-6;

// ---------------------------------------------------------------------------
// Tolerance helpers (def.h:181-193, set.c)
// ---------------------------------------------------------------------------

/// def.h:187, set.c:6946
pub fn s2_is_zero(x: f64) -> bool {
    x.abs() <= S2_EPSILON
}

/// def.h:185, set.c:6906
pub fn s2_is_gt(x: f64, y: f64) -> bool {
    x - y > S2_EPSILON
}

/// def.h:186
pub fn s2_is_ge(x: f64, y: f64) -> bool {
    x - y >= -S2_EPSILON
}

/// def.h:183
pub fn s2_is_lt(x: f64, y: f64) -> bool {
    x - y < -S2_EPSILON
}

/// def.h:184
pub fn s2_is_le(x: f64, y: f64) -> bool {
    x - y <= S2_EPSILON
}

/// def.h:182
pub fn s2_is_eq(x: f64, y: f64) -> bool {
    (x - y).abs() <= S2_EPSILON
}

/// def.h:188
pub fn s2_is_positive(x: f64) -> bool {
    x > S2_EPSILON
}

/// def.h:189
pub fn s2_is_negative(x: f64) -> bool {
    x < -S2_EPSILON
}

/// def.h:162, set.c:6824
pub fn s2_is_infinity(x: f64) -> bool {
    x >= kHighsInf
}

/// def.h:190, set.c:7025
pub fn s2_floor(x: f64) -> f64 {
    (x + S2_EPSILON).floor()
}

/// def.h:191, set.c:7037
pub fn s2_ceil(x: f64) -> f64 {
    (x - S2_EPSILON).ceil()
}

/// set.c:7113
pub fn s2_is_sum_le(x: f64, y: f64) -> bool {
    x - y <= S2_SUMEPSILON
}

/// set.c:7206
pub fn s2_sum_floor(x: f64) -> f64 {
    (x + S2_SUMEPSILON).floor()
}

// ---------------------------------------------------------------------------
// S2ConflictRow
// ---------------------------------------------------------------------------

/// SCIP_CONFLICTROW (conflict_resolution.c)
#[derive(Clone, Default, Debug)]
pub struct S2ConflictRow {
    pub vals: Vec<f64>,
    pub inds: Vec<HighsInt>,
    pub lhs: f64,
    pub slack: f64,
    pub coefquotient: f64,
}

impl S2ConflictRow {
    /// nnz = number of indices
    pub fn nnz(&self) -> usize {
        self.inds.len()
    }

    /// Ensure vals has length n; if it changes, zero all entries.
    pub fn ensure_dim(&mut self, n: usize) {
        if self.vals.len() != n {
            self.vals = vec![0.0; n];
        }
    }

    /// SCIP conflictRowClear (lines 1257-1287)
    pub fn clear(&mut self, nvars: usize) {
        self.ensure_dim(nvars);
        for v in self.vals.iter_mut() {
            *v = 0.0;
        }
        self.inds.clear();
        self.lhs = 0.0;
        self.slack = 0.0;
        self.coefquotient = 0.0;
    }

    /// SCIP conflictRowReplace (lines 140-180)
    pub fn replace_from(&mut self, src: &S2ConflictRow) {
        let nvars = src.vals.len();

        // copy all data from source to target (nnz is inds.len(), so inds is set to the source's)
        self.inds.clear();
        self.inds.extend_from_slice(&src.inds);
        self.ensure_dim(nvars);
        self.vals.copy_from_slice(&src.vals);

        self.lhs = src.lhs;
        self.slack = src.slack;
        self.coefquotient = src.coefquotient;
    }
}

// ---------------------------------------------------------------------------
// S2Change
// ---------------------------------------------------------------------------

/// A bound change from the domain stack.
#[derive(Clone, Copy, Debug)]
pub struct S2Change {
    pub col: HighsInt,
    pub is_upper: bool,
    pub newbound: f64,
    pub oldbound: f64,
    pub depth: HighsInt,
    pub branching: bool,
    pub resolvable: bool,
}

// ---------------------------------------------------------------------------
// S2Context trait
// ---------------------------------------------------------------------------

/// Bound access trait for conflict resolution.
pub trait S2Context {
    fn num_col(&self) -> usize;
    fn num_integral(&self) -> usize;
    fn is_integral(&self, col: HighsInt) -> bool;
    fn global_lb(&self, col: HighsInt) -> f64;
    fn global_ub(&self, col: HighsInt) -> f64;
    fn lb_upto(&self, col: HighsInt, upto: HighsInt) -> (f64, HighsInt);
    fn ub_upto(&self, col: HighsInt, upto: HighsInt) -> (f64, HighsInt);
    fn change(&self, pos: HighsInt) -> S2Change;
    fn reason_row(&mut self, pos: HighsInt, out: &mut S2ConflictRow) -> bool;
}

// ---------------------------------------------------------------------------
// Free functions: bound-at-index, binary check
// ---------------------------------------------------------------------------

/// Get the lower bound of col at or before pos (after=true) or before pos (after=false).
/// pos == kHighsIInf means current bounds.
pub fn s2_lb_at<C: S2Context>(ctx: &C, col: HighsInt, pos: HighsInt, after: bool) -> f64 {
    let upto = if pos == kHighsIInf {
        kHighsIInf
    } else if after {
        pos
    } else {
        pos - 1
    };
    ctx.lb_upto(col, upto).0
}

/// Get the upper bound of col at or before pos (after=true) or before pos (after=false).
/// pos == kHighsIInf means current bounds.
pub fn s2_ub_at<C: S2Context>(ctx: &C, col: HighsInt, pos: HighsInt, after: bool) -> f64 {
    let upto = if pos == kHighsIInf {
        kHighsIInf
    } else if after {
        pos
    } else {
        pos - 1
    };
    ctx.ub_upto(col, upto).0
}

/// SCIPvarIsBinary: integral, global lb == 0 and global ub == 1.
pub fn s2_is_binary<C: S2Context>(ctx: &C, col: HighsInt) -> bool {
    ctx.is_integral(col) && ctx.global_lb(col) == 0.0 && ctx.global_ub(col) == 1.0
}

// ---------------------------------------------------------------------------
// compute_slack
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::computeSlack (lines 1291-1368)
///
/// Calculates the slack (maxact - rhs) for a generalized resolution row.
/// Documented deviation: if any used bound is infinite, returns kHighsInf
/// (SCIP's +1e20-dominated result; HighsCDouble two-sum over infinity gives NaN).
pub fn compute_slack<C: S2Context>(
    ctx: &C,
    row: &mut S2ConflictRow,
    pos: HighsInt,
    fix: Option<(&[f64], &[i8])>,
) {
    let mut slack = HighsCDouble::from(0.0);

    for i in 0..row.nnz() {
        let v = row.inds[i];
        let coef = row.vals[v as usize];
        let bound: f64;

        if coef > 0.0 {
            if let Some((fixbounds, fixsides)) = fix {
                if fixsides[v as usize] == 1 {
                    bound = fixbounds[v as usize];
                } else {
                    bound = s2_ub_at(ctx, v, pos, true);
                }
            } else {
                bound = s2_ub_at(ctx, v, pos, true);
            }
        } else {
            if let Some((fixbounds, fixsides)) = fix {
                if fixsides[v as usize] == -1 {
                    bound = fixbounds[v as usize];
                } else {
                    bound = s2_lb_at(ctx, v, pos, true);
                }
            } else {
                bound = s2_lb_at(ctx, v, pos, true);
            }
        }

        // Infinite-bound early exit (documented deviation)
        if bound.abs() >= kHighsInf {
            row.slack = kHighsInf;
            return;
        }

        // SCIPquadprecProdDD(delta, coef, bound)
        let delta = HighsCDouble::from(coef) * bound;
        // SCIPquadprecSumQQ(slack, slack, delta)
        slack += delta;
    }

    // SCIPquadprecSumQD(slack, slack, -row->lhs)
    slack += -row.lhs;
    row.slack = f64::from(slack);
}

// ---------------------------------------------------------------------------
// remove_zero_var / remove_zero_vars
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictRowRemoveZeroVar (lines 782-799)
pub fn remove_zero_var(row: &mut S2ConflictRow, pos: usize) {
    let idx = row.inds[pos] as usize;
    row.vals[idx] = 0.0;
    row.inds.swap_remove(pos);
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictRowRemoveZeroVars (lines 803-816)
pub fn remove_zero_vars(row: &mut S2ConflictRow) {
    let nnz = row.nnz();
    for i in (0..nnz).rev() {
        let idx = row.inds[i] as usize;
        if s2_is_zero(row.vals[idx]) {
            remove_zero_var(row, i);
        }
    }
}

// ---------------------------------------------------------------------------
// Diagnostic: conflict-row validity against the debug solution
// ---------------------------------------------------------------------------

/// Diagnostic only (never changes the solve): `RUSOLVE_S2_TRACE=1` reports the first
/// resolution step at which the conflict row stops being satisfied by the debug solution
/// (`mip_debug_solution_file`, see `mip::debug_sol`). Every transformation in
/// `conflict_resolution.c` is a non-negative combination of valid rows, a valid weakening,
/// or a valid MIR/CMIR reduction, so the row stays valid at every step: a step that makes
/// the debug solution violate the row is the first bad step.
pub(crate) fn s2_trace_row(tag: &str, row: &S2ConflictRow, dbg: &S2DebugSol) {
    if !s2_trace_env() {
        return;
    }
    let step = S2_TRACE_STEP.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let Some(x) = dbg.x.as_deref() else {
        return;
    };
    if x.len() != row.vals.len() {
        return;
    }
    // The row is sum_j vals[j] x_j >= lhs; its slack w.r.t. x is activity - lhs.
    let mut act = HighsCDouble::from(0.0);
    for &col in &row.inds {
        act += x[col as usize] * row.vals[col as usize];
    }
    let viol = f64::from(HighsCDouble::from(row.lhs) - act);
    if viol <= S2_FEASTOL {
        return;
    }
    eprintln!(
        "s2-trace step={step} {tag}: nnz={} lhs={:e} slack={:e} VIOL={viol:e}",
        row.nnz(),
        row.lhs,
        row.slack,
    );
    if S2_TRACE_DUMP.load(std::sync::atomic::Ordering::Relaxed) != 0 {
        let mut terms = String::new();
        for &col in &row.inds {
            terms += &format!(" {}:{:e}", col, row.vals[col as usize]);
        }
        eprintln!("s2-trace step={step} row:{terms}");
    }
}

/// `RUSOLVE_S2_TRACE_DUMP=1` also dumps the terms of the first invalid row.
pub(crate) fn s2_trace_dump_enable(on: bool) {
    S2_TRACE_DUMP.store(on as i64, std::sync::atomic::Ordering::Relaxed);
}

static S2_TRACE_STEP: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static S2_TRACE_DUMP: std::sync::atomic::AtomicI64 = std::sync::atomic::AtomicI64::new(0);

fn s2_trace_env() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("RUSOLVE_S2_TRACE").is_ok_and(|v| v == "1"))
}

/// `RUSOLVE_S2_TRACE=1`? Exported for the diagnostic wrappers in the other files.
pub(crate) fn s2_trace_enabled() -> bool {
    s2_trace_env()
}

/// Violation of a ">=" row `sum terms.1 * x[terms.0] >= lhs` by the cached debug solution:
/// positive when the solution does not satisfy the row. Diagnostic only.
pub(crate) fn s2_debugsol_viol_ge(dbg: &S2DebugSol, terms: &[(HighsInt, f64)], lhs: f64) -> f64 {
    let Some(x) = dbg.x.as_deref() else {
        return 0.0;
    };
    let mut act = HighsCDouble::from(0.0);
    for &(col, coef) in terms {
        let Some(v) = x.get(col as usize) else {
            return 0.0;
        };
        act += *v * coef;
    }
    f64::from(HighsCDouble::from(lhs) - act)
}

// ---------------------------------------------------------------------------
// Data-only structs (methods added by later jobs)
// ---------------------------------------------------------------------------

/// SCIP pqueue for bound-change indices (max-heap over positions).
#[derive(Clone, Default, Debug)]
pub struct S2BdchgQueue {
    pub heap: BinaryHeap<HighsInt>,
    pub conflictvarslbs: Vec<f64>,
    pub conflictvarsubs: Vec<f64>,
}

/// Statistics for conflict resolution.
#[derive(Clone, Default, Debug)]
pub struct S2Stats {
    pub calls: i64,
    pub initial_fail: i64,
    pub nressteps: i64,
    pub nreslongconfs: i64,
    pub rows_found: i64,
    pub rows_added: i64,
    pub empty_infeasible: i64,
    pub debugsol_checked: i64,
    pub debugsol_violations: i64,
}

impl S2Stats {
    fn add(&mut self, o: &S2Stats) {
        self.calls += o.calls;
        self.initial_fail += o.initial_fail;
        self.nressteps += o.nressteps;
        self.nreslongconfs += o.nreslongconfs;
        self.rows_found += o.rows_found;
        self.rows_added += o.rows_added;
        self.empty_infeasible += o.empty_infeasible;
        self.debugsol_checked += o.debugsol_checked;
        self.debugsol_violations += o.debugsol_violations;
    }
}

/// Debug solution state: a copy of the `HighsDebugSol` solution mapped into the presolved space
/// (`mip::debug_sol`), refreshed from it on every `s2_conflict_analysis` call, and empty whenever
/// no debug solution is active. `announced` prints the column mapping once per solution.
#[derive(Clone, Default, Debug)]
pub struct S2DebugSol {
    pub x: Option<Vec<f64>>,
    pub announced: bool,
}

/// The conflict resolution workspace. Inert on the off path.
#[derive(Clone, Default, Debug)]
pub struct HighsConflictResolution {
    pub conflictrow: S2ConflictRow,
    pub resolvedconflictrow: S2ConflictRow,
    pub reasonrow: S2ConflictRow,
    pub reducedreasonrow: S2ConflictRow,
    pub conflictrows: Vec<S2ConflictRow>,
    pub resqueue: S2BdchgQueue,
    pub contqueue: BinaryHeap<HighsInt>,
    pub fixsides: Vec<i8>,
    pub fixbounds: Vec<f64>,
    pub nressteps: HighsInt,
    pub nfuips: HighsInt,
    pub stats: S2Stats,
    pub debug_sol: S2DebugSol,
    pub stats_env: Option<bool>,
}

/// One [`HighsConflictResolution`] per conflict pool (port glue, no SCIP counterpart).
///
/// SCIP runs one conflict analysis at a time. In the port, the MIP workers of
/// `HighsMipSolver::run_task` analyse conflicts concurrently, each into its own conflict
/// pool, so the workspace is kept per pool: [`Self::take`] moves pool `i`'s workspace out
/// and [`Self::put`] returns it, each under a short lock, and the analysis itself runs
/// unlocked. A single shared workspace taken with `mem::take` from two workers at once
/// gave both of them the same buffers (double free at threads=8, parallel=on).
#[derive(Default, Debug)]
pub struct S2Workspaces(std::sync::Mutex<Vec<HighsConflictResolution>>);

impl Clone for S2Workspaces {
    fn clone(&self) -> Self {
        Self(std::sync::Mutex::new(self.lock().clone()))
    }
}

impl S2Workspaces {
    fn lock(&self) -> std::sync::MutexGuard<'_, Vec<HighsConflictResolution>> {
        self.0.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Move conflict pool `pool`'s workspace out (empty on first use).
    pub fn take(&self, pool: usize) -> HighsConflictResolution {
        let mut v = self.lock();
        if v.len() <= pool {
            v.resize_with(pool + 1, HighsConflictResolution::default);
        }
        std::mem::take(&mut v[pool])
    }

    /// Return conflict pool `pool`'s workspace.
    pub fn put(&self, pool: usize, ws: HighsConflictResolution) {
        let mut v = self.lock();
        if v.len() <= pool {
            v.resize_with(pool + 1, HighsConflictResolution::default);
        }
        v[pool] = ws;
    }

    /// Statistics summed over all pools.
    pub fn stats(&self) -> S2Stats {
        let v = self.lock();
        let mut total = S2Stats::default();
        for ws in v.iter() {
            total.add(&ws.stats);
        }
        total
    }
}

impl HighsConflictResolution {
    /// Resize all internal structures to n columns. On path only.
    pub fn ensure_dim(&mut self, n: usize) {
        self.conflictrow.ensure_dim(n);
        self.resolvedconflictrow.ensure_dim(n);
        self.reasonrow.ensure_dim(n);
        self.reducedreasonrow.ensure_dim(n);
        for row in self.conflictrows.iter_mut() {
            row.ensure_dim(n);
        }
        self.fixsides.resize(n, 0);
        self.fixbounds.resize(n, 0.0);
    }
}

// ---------------------------------------------------------------------------
// Test context
// ---------------------------------------------------------------------------

#[cfg(test)]
pub(crate) struct S2TestChg {
    pub col: HighsInt,
    pub is_upper: bool,
    pub newbound: f64,
    pub depth: HighsInt,
    pub branching: bool,
    pub reason: Option<(Vec<(HighsInt, f64)>, f64)>,
}

#[cfg(test)]
pub(crate) struct S2TestCtx {
    pub glb: Vec<f64>,
    pub gub: Vec<f64>,
    pub integral: Vec<bool>,
    pub num_integral: usize,
    pub stack: Vec<S2TestChg>,
}

#[cfg(test)]
impl S2Context for S2TestCtx {
    fn num_col(&self) -> usize {
        self.glb.len()
    }

    fn num_integral(&self) -> usize {
        self.num_integral
    }

    fn is_integral(&self, col: HighsInt) -> bool {
        self.integral[col as usize]
    }

    fn global_lb(&self, col: HighsInt) -> f64 {
        self.glb[col as usize]
    }

    fn global_ub(&self, col: HighsInt) -> f64 {
        self.gub[col as usize]
    }

    fn lb_upto(&self, col: HighsInt, upto: HighsInt) -> (f64, HighsInt) {
        let mut best_bound = self.glb[col as usize];
        let mut best_pos: HighsInt = -1;
        let limit = if upto == kHighsIInf {
            self.stack.len() as HighsInt
        } else {
            std::cmp::min(upto + 1, self.stack.len() as HighsInt)
        };
        for p in 0..limit {
            let chg = &self.stack[p as usize];
            // Only a change that moved the bound counts (manifest the port).
            if chg.col == col && !chg.is_upper && chg.newbound != best_bound {
                best_bound = chg.newbound;
                best_pos = p;
            }
        }
        (best_bound, best_pos)
    }

    fn ub_upto(&self, col: HighsInt, upto: HighsInt) -> (f64, HighsInt) {
        let mut best_bound = self.gub[col as usize];
        let mut best_pos: HighsInt = -1;
        let limit = if upto == kHighsIInf {
            self.stack.len() as HighsInt
        } else {
            std::cmp::min(upto + 1, self.stack.len() as HighsInt)
        };
        for p in 0..limit {
            let chg = &self.stack[p as usize];
            // Only a change that moved the bound counts (manifest the port).
            if chg.col == col && chg.is_upper && chg.newbound != best_bound {
                best_bound = chg.newbound;
                best_pos = p;
            }
        }
        (best_bound, best_pos)
    }

    fn change(&self, pos: HighsInt) -> S2Change {
        let chg = &self.stack[pos as usize];
        let col = chg.col;
        let is_upper = chg.is_upper;

        // Compute oldbound: the bound of that side before pos
        let oldbound = if is_upper {
            let mut b = self.gub[col as usize];
            for p in 0..pos {
                let chg2 = &self.stack[p as usize];
                if chg2.col == col && chg2.is_upper {
                    b = chg2.newbound;
                }
            }
            b
        } else {
            let mut b = self.glb[col as usize];
            for p in 0..pos {
                let chg2 = &self.stack[p as usize];
                if chg2.col == col && !chg2.is_upper {
                    b = chg2.newbound;
                }
            }
            b
        };

        S2Change {
            col,
            is_upper,
            newbound: chg.newbound,
            oldbound,
            depth: chg.depth,
            branching: chg.branching,
            resolvable: chg.reason.is_some() && !chg.branching,
        }
    }

    fn reason_row(&mut self, pos: HighsInt, out: &mut S2ConflictRow) -> bool {
        let chg = &self.stack[pos as usize];
        if let Some((ref terms, rhs)) = chg.reason {
            let nvars = self.num_col();
            out.clear(nvars);
            for &(col2, v) in terms {
                out.vals[col2 as usize] = v;
                out.inds.push(col2);
            }
            out.lhs = rhs;
            true
        } else {
            false
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::mip_solver_data::HighsMipSolverData;
    use crate::prelude::*;

    #[test]
    fn s2_row_slack_basic() {
        // 2x0 - 3x1 + x2 >= 1, lb [0,0,0], ub [1,2,4], no changes
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![1.0, 2.0, 4.0],
            integral: vec![true, true, true],
            num_integral: 3,
            stack: vec![],
        };
        let mut row = S2ConflictRow {
            vals: vec![2.0, -3.0, 1.0],
            inds: vec![0, 1, 2],
            lhs: 1.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        compute_slack(&ctx, &mut row, kHighsIInf, None);
        // maxact = 2*1 + (-3)*0 + 1*4 = 2 + 0 + 4 = 6; slack = 6 - 1 = 5
        assert_eq!(row.slack, 5.0);
    }

    #[test]
    fn s2_row_slack_fixsides() {
        // Same row, fixsides[1] = -1, fixbounds[1] = 1.0
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![1.0, 2.0, 4.0],
            integral: vec![true, true, true],
            num_integral: 3,
            stack: vec![],
        };
        let mut row = S2ConflictRow {
            vals: vec![2.0, -3.0, 1.0],
            inds: vec![0, 1, 2],
            lhs: 1.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        let fixsides: Vec<i8> = vec![0, -1, 0];
        let fixbounds: Vec<f64> = vec![0.0, 1.0, 0.0];
        compute_slack(&ctx, &mut row, kHighsIInf, Some((&fixbounds, &fixsides)));
        // coef[1] = -3 < 0, fixsides[1] == -1 -> bound = fixbounds[1] = 1.0
        // maxact = 2*1 + (-3)*1 + 1*4 = 2 - 3 + 4 = 3; slack = 3 - 1 = 2
        assert_eq!(row.slack, 2.0);
    }

    #[test]
    fn s2_row_slack_infinite() {
        // Same row but gub[0] = kHighsInf
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![kHighsInf, 2.0, 4.0],
            integral: vec![true, true, true],
            num_integral: 3,
            stack: vec![],
        };
        let mut row = S2ConflictRow {
            vals: vec![2.0, -3.0, 1.0],
            inds: vec![0, 1, 2],
            lhs: 1.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        compute_slack(&ctx, &mut row, kHighsIInf, None);
        assert_eq!(row.slack, kHighsInf);
    }

    #[test]
    fn s2_row_remove_zero_vars() {
        // inds [0,1,2], vals [0, 1, 0]
        let mut row = S2ConflictRow {
            vals: vec![0.0, 1.0, 0.0],
            inds: vec![0, 1, 2],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        remove_zero_vars(&mut row);
        // Backward loop: i=2, idx=2, vals[2]=0 -> remove (swap_remove(2) removes last, inds=[0,1])
        // i=1, idx=1, vals[1]=1 -> keep
        // i=0, idx=0, vals[0]=0 -> remove (swap_remove(0) puts last(1) at 0, inds=[1])
        assert_eq!(row.inds, vec![1]);
    }

    #[test]
    fn s2_row_at_index() {
        // glb [0.0], gub [5.0], stack: change 0 (col 0, lb 0->1), change 1 (col 0, lb 1->2)
        let ctx = S2TestCtx {
            glb: vec![0.0],
            gub: vec![5.0],
            integral: vec![true],
            num_integral: 1,
            stack: vec![
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 1.0,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 2.0,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
            ],
        };
        assert_eq!(s2_lb_at(&ctx, 0, 1, true), 2.0);
        assert_eq!(s2_lb_at(&ctx, 0, 1, false), 1.0);
        assert_eq!(s2_lb_at(&ctx, 0, 0, false), 0.0);
    }

    #[test]
    fn s2_row_default_inert() {
        let hms = HighsMipSolverData::default();
        assert!(hms.conflict_resolution.take(0).conflictrows.is_empty());
        assert_eq!(hms.conflict_resolution.stats().calls, 0);
    }
}
