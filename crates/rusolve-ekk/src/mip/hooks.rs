//! rusolve extension hooks (design: docs/design/hooks.md): optional, default-off hook points
//! in the MIP search, and the validity gates on what they add.
//!
//! This file is rusolve's own (not a HiGHS translation). The HiGHS callback system itself
//! (`HighsCallback`, `kCallbackMip*`) is ported in `lp_data/highs_callback.rs` and fires at the
//! HiGHS call sites; this file adds the hook points HiGHS does not have. Each site sits beside
//! the HiGHS code named in `docs/design/hooks.md`.
//!
//! With no hooks installed (`HighsCallback::hooks == None`) every call site reduces to one
//! `Option::is_some` test: no RNG draw, no allocation, no extra work, so the search is identical.
//!
//! Space: hooks run inside the MIP solver, i.e. in the presolved column space. Set `presolve=off`
//! (or `mip_presolve`-equivalent) when a hook needs original column indices. User solutions
//! offered through the HiGHS `kCallbackMipUserSolution` callback are in the original space.
//!
//! Validity gates: every cut a `separate` hook returns is checked before it reaches the cut
//! pool: malformed rows are rejected, and when the debug-solution harness is active
//! (`mip_debug_solution_file`, `mip/debug_sol.rs`) a cut that cuts off the debug solution is
//! rejected and reported on stderr with the hook name. Rejections are counted in [`HookStats`].

use crate::lp_data::highs_run::Highs;
use crate::mip::node_queue::HighsNodeQueue_OpenNode;
use crate::prelude::*;
use std::collections::HashSet;
use std::sync::{Arc, Mutex, MutexGuard};

/// One fractional integer column at the current node.
#[derive(Clone, Copy, Debug)]
pub struct BranchingCandidate {
    pub col: HighsInt,
    pub value: f64,
}

/// Read-only view handed to [`MipHooks::branching_score`].
pub struct BranchingInfo<'a> {
    pub depth: HighsInt,
    pub node_lower_bound: f64,
    pub upper_limit: f64,
    pub candidates: &'a [BranchingCandidate],
    /// Node LP solution (column values).
    pub solution: &'a [f64],
    /// Local (node) column bounds.
    pub col_lower: &'a [f64],
    pub col_upper: &'a [f64],
}

/// One open node of the global queue, handed to [`MipHooks::select_node`].
#[derive(Clone, Copy, Debug)]
pub struct OpenNodeInfo {
    pub lower_bound: f64,
    pub estimate: f64,
    pub depth: HighsInt,
    pub num_domain_changes: usize,
}

/// A violated pool cut about to be scored, handed to [`MipHooks::cut_score`].
#[derive(Clone, Copy, Debug)]
pub struct CutScoreInfo {
    pub pool_row: HighsInt,
    pub violation: f64,
    pub num_nonzeros: HighsInt,
    pub num_active_nonzeros: HighsInt,
    /// Euclidean norm over the active nonzeros (the HiGHS efficacy denominator part).
    pub active_norm: f64,
    pub default_score: f64,
}

/// Read-only view handed to [`MipHooks::separate`].
pub struct SeparationInfo<'a> {
    pub lp_objective: f64,
    pub upper_limit: f64,
    pub feastol: f64,
    pub num_nodes: i64,
    pub solution: &'a [f64],
    pub col_lower: &'a [f64],
    pub col_upper: &'a [f64],
}

/// A cut `sum_j values[j] * x[indices[j]] <= rhs`, assumed globally valid.
#[derive(Clone, Debug)]
pub struct UserCut {
    pub indices: Vec<HighsInt>,
    pub values: Vec<f64>,
    pub rhs: f64,
}

/// rusolve hook points. Every method has a no-op default; `None` / an empty list means
/// "no opinion, use the HiGHS rule". The hooks can only reorder or add rows; whatever they add is
/// gated in this file.
pub trait MipHooks: Send {
    /// Name used in validity-gate reports.
    fn name(&self) -> &str {
        "user"
    }

    /// Sits beside `HighsSearch::selectBranchingCandidate` (HighsSearch.cpp:247). One score per
    /// candidate; the highest wins (ties: first). Returning `Some` skips the HiGHS reliability
    /// strong-branching selection at that node. `None` keeps the HiGHS rule.
    fn branching_score(&mut self, _info: &BranchingInfo<'_>) -> Option<Vec<f64>> {
        None
    }

    /// Sits beside `HighsNodeQueue::popBestNode` (HighsMipSolver.cpp, `installNodes` and
    /// `prepareNodes`). Return the index into `open` of the node to dive into next.
    fn select_node(&mut self, _open: &[OpenNodeInfo]) -> Option<usize> {
        None
    }

    /// Sits beside the efficacy score in `HighsCutPool::separate` (HighsCutPool.cpp:229). The
    /// returned score replaces the default one (larger = better); `None` keeps it.
    fn cut_score(&mut self, _info: &CutScoreInfo) -> Option<f64> {
        None
    }

    /// Sits after the built-in separators in `HighsSeparation::separationRound`
    /// (HighsSeparation.cpp:43). Cuts returned are checked by the validity gate and then added to the
    /// cut pool; only violated ones reach the LP.
    fn separate(&mut self, _info: &SeparationInfo<'_>) -> Vec<UserCut> {
        Vec::new()
    }
}

/// Counters of what the hooks did and what the gates rejected (D3: "counted and reported").
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HookStats {
    pub branching_calls: u64,
    pub branching_overrides: u64,
    pub node_select_calls: u64,
    pub node_select_overrides: u64,
    pub cut_score_calls: u64,
    pub separator_calls: u64,
    pub user_cuts_offered: u64,
    pub user_cuts_accepted: u64,
    pub user_cuts_duplicate: u64,
    pub user_cuts_rejected_malformed: u64,
    pub user_cuts_rejected_debug_solution: u64,
    /// `kCallbackMipUserSolution` solutions that passed the original-model feasibility check.
    pub user_solutions_checked_feasible: u64,
    /// `kCallbackMipUserSolution` solutions rejected by the original-model feasibility check.
    pub user_solutions_rejected: u64,
}

/// Shared, thread-safe counters.
pub type HookStatsHandle = Arc<Mutex<HookStats>>;

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

struct HooksShared {
    hooks: Mutex<Box<dyn MipHooks>>,
    name: String,
    stats: HookStatsHandle,
}

/// The installed hooks (`Option<MipHooksHandle>` on [`HighsCallback`]). Cloning shares them.
#[derive(Clone)]
pub struct MipHooksHandle(Arc<HooksShared>);

impl std::fmt::Debug for MipHooksHandle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MipHooksHandle({})", self.0.name)
    }
}

impl MipHooksHandle {
    pub fn new(hooks: Box<dyn MipHooks>, stats: HookStatsHandle) -> Self {
        let name = hooks.name().to_string();
        MipHooksHandle(Arc::new(HooksShared {
            hooks: Mutex::new(hooks),
            name,
            stats,
        }))
    }

    pub fn name(&self) -> &str {
        &self.0.name
    }

    /// Index of the best-scored candidate, or `None` to keep the HiGHS rule.
    pub fn branching_choice(&self, info: &BranchingInfo<'_>) -> Option<usize> {
        let scores = lock(&self.0.hooks).branching_score(info);
        let mut st = lock(&self.0.stats);
        st.branching_calls += 1;
        let scores = scores?;
        if scores.len() != info.candidates.len() || scores.is_empty() {
            return None;
        }
        let mut best = 0usize;
        for (i, &s) in scores.iter().enumerate() {
            if s.is_nan() {
                return None;
            }
            if s > scores[best] {
                best = i;
            }
        }
        st.branching_overrides += 1;
        Some(best)
    }

    pub fn node_choice(&self, open: &[OpenNodeInfo]) -> Option<usize> {
        let pick = lock(&self.0.hooks).select_node(open);
        let mut st = lock(&self.0.stats);
        st.node_select_calls += 1;
        let pick = pick.filter(|&i| i < open.len())?;
        st.node_select_overrides += 1;
        Some(pick)
    }

    /// The replacement score for one violated cut (the default when the hook has no opinion or
    /// returns a non-finite value).
    pub fn cut_score(&self, info: &CutScoreInfo) -> f64 {
        let s = lock(&self.0.hooks).cut_score(info);
        lock(&self.0.stats).cut_score_calls += 1;
        match s {
            Some(v) if v.is_finite() => v,
            _ => info.default_score,
        }
    }

    fn separate(&self, info: &SeparationInfo<'_>) -> Vec<UserCut> {
        let cuts = lock(&self.0.hooks).separate(info);
        let mut st = lock(&self.0.stats);
        st.separator_calls += 1;
        st.user_cuts_offered += cuts.len() as u64;
        cuts
    }
}

impl Highs {
    /// Install rusolve MIP hooks. Replaces any previous hooks and resets the counters.
    /// Independent of the HiGHS `set_callback` system.
    pub fn set_mip_hooks(&mut self, hooks: Box<dyn MipHooks>) {
        let stats: HookStatsHandle = Arc::new(Mutex::new(HookStats::default()));
        self.callback_.hook_stats = Arc::clone(&stats);
        self.callback_.hooks = Some(MipHooksHandle::new(hooks, stats));
    }

    /// Remove the hooks; the search is the unmodified HiGHS one again.
    pub fn clear_mip_hooks(&mut self) {
        self.callback_.hooks = None;
    }

    /// Counters of the last run (hook activity, gate rejections, rejected user solutions).
    pub fn mip_hook_stats(&self) -> HookStats {
        lock(&self.callback_.hook_stats).clone()
    }
}

impl HighsMipSolver {
    /// The installed hooks, or `None` in sub-MIPs (they inherit the parent's `callback_`, but
    /// hooks steer the main search only, as the HiGHS callbacks do: `!submip`).
    pub(crate) fn user_hooks(&self) -> Option<MipHooksHandle> {
        if self.submip {
            None
        } else {
            self.callback_.hooks.clone()
        }
    }

    /// `HighsNodeQueue::popBestNode`, unless a `select_node` hook picks another open node.
    /// With no hooks this is exactly `nodequeue.pop_best_node()`.
    pub(crate) fn pop_best_node_hooked(&mut self) -> HighsNodeQueue_OpenNode {
        if let Some(h) = self.user_hooks() {
            let q = &self.mipdata_.nodequeue;
            let free: HashSet<i64> = q.freeslots.iter().map(|r| r.0).collect();
            let mut ids = Vec::new();
            let mut open = Vec::new();
            for (pos, n) in q.nodes.iter().enumerate() {
                // active (not freed, not suboptimal) nodes only
                if free.contains(&(pos as i64)) || n.estimate == kHighsInf {
                    continue;
                }
                ids.push(pos as i64);
                open.push(OpenNodeInfo {
                    lower_bound: n.lower_bound,
                    estimate: n.estimate,
                    depth: n.depth,
                    num_domain_changes: n.domchgstack.len(),
                });
            }
            if let Some(k) = h.node_choice(&open) {
                // same steps as HighsNodeQueue::popBestNode with another slot
                let pos = ids[k];
                self.mipdata_.nodequeue.unlink(pos);
                return std::mem::take(&mut self.mipdata_.nodequeue.nodes[pos as usize]);
            }
        }
        self.mipdata_.nodequeue.pop_best_node()
    }
}

/// Run the user separator and pass its cuts through the validity gate into the cut pool `cutpool`.
/// Returns the number of cuts added to the pool.
pub(crate) fn run_user_separator(
    ms: &mut HighsMipSolver,
    lp: &HighsLpRelaxation,
    dom: &HighsDomain,
    cutpool: usize,
) -> HighsInt {
    let Some(h) = ms.user_hooks() else {
        return 0;
    };
    let cuts = {
        let info = SeparationInfo {
            lp_objective: lp.get_objective(),
            upper_limit: ms.mipdata_.upper_limit,
            feastol: ms.mipdata_.feastol,
            num_nodes: ms.mipdata_.num_nodes,
            solution: &lp.get_solution().col_value,
            col_lower: &dom.col_lower_,
            col_upper: &dom.col_upper_,
        };
        h.separate(&info)
    };
    let num_col = dom.col_lower_.len() as HighsInt;
    let mut added = 0;
    for cut in cuts {
        let UserCut {
            mut indices,
            mut values,
            rhs,
        } = cut;
        // validity gate 1: well-formed
        let mut seen = HashSet::new();
        let malformed = indices.len() != values.len()
            || indices.is_empty()
            || !rhs.is_finite()
            || values.iter().any(|v| !v.is_finite())
            || indices
                .iter()
                .any(|&col| col < 0 || col >= num_col || !seen.insert(col));
        if malformed {
            lock(&h.0.stats).user_cuts_rejected_malformed += 1;
            eprintln!(
                "HOOK VALIDITY GATE: hook '{}' separator returned a malformed cut (rejected; num_col {})",
                h.name(),
                num_col
            );
            continue;
        }
        // validity gate 2: the debug solution (mip_debug_solution_file), when active
        let dbg = &ms.mipdata_.debug_solution;
        if dbg.debug_sol_active && dbg.debug_solution.len() == num_col as usize {
            let mut viol = HighsCDouble::from(-rhs);
            for (&col, &v) in indices.iter().zip(&values) {
                viol += HighsCDouble::from(dbg.debug_solution[col as usize]) * v;
            }
            let viol = f64::from(viol);
            if viol > dbg.feastol {
                lock(&h.0.stats).user_cuts_rejected_debug_solution += 1;
                let msg = format!(
                    "HOOK VALIDITY GATE: hook '{}' separator cut of length {} rhs {} cuts off the debug solution (violation {}); rejected",
                    h.name(),
                    indices.len(),
                    rhs,
                    viol
                );
                if std::env::var("RUSOLVE_DEBUG_SOL_ABORT").is_ok_and(|v| v == "1") {
                    panic!("{msg}");
                }
                eprintln!("{msg}");
                continue;
            }
        }
        let len = indices.len() as HighsInt;
        let r = detach!(ms.mipdata_.cutpools[cutpool], |cp| cp.add_cut(
            ms,
            &mut indices,
            &mut values,
            len,
            rhs,
            false,
            false,
            false,
            false
        ));
        let mut st = lock(&h.0.stats);
        if r < 0 {
            st.user_cuts_duplicate += 1;
        } else {
            st.user_cuts_accepted += 1;
            added += 1;
        }
    }
    added
}
