//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! mip-heuristics (36f7d9d): src/repair_search.h + src/repair_search.cpp (lines 1-785)
//!
//! RepairSearch on a *complete* assignment (paper Fig. 5, Sect. 5.1): a DFS
//! over a two-branch disjunction per repair move, with a secondary
//! propagation engine R that carries the branch and whose deductions are
//! transferred back into the primary engine E by `sync_changes`.

use crate::HighsInt;
use crate::mip::fpr_local_mip::Deadline;
use crate::mip::fpr_local_mip::FprScratch;
use crate::mip::fpr_local_mip::PropEngine;
use crate::mip::fpr_local_mip::PropResult;
use crate::mip::fpr_local_mip::RepairSearchNode;
use crate::mip::fpr_local_mip::rng::S5Rng;
use crate::mip::fpr_local_mip::rng::is_row_violated;
use crate::mip::fpr_local_mip::rng::row_violation;
use crate::mip::fpr_local_mip::rng::s5_uniform_int_i32;
use crate::mip::fpr_local_mip::walksat::walksat_select_move;
use crate::util::highs_utils::cpp_max;
use crate::util::highs_utils::cpp_min;

/// mip-heuristics (36f7d9d): repair_search.h::kRepairProgressThreshold (line 48)
pub(crate) const K_REPAIR_PROGRESS_THRESHOLD: HighsInt = 10;

/// mip-heuristics (36f7d9d): repair_search.h::RepairSearchStats (lines 57-95)
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct RepairSearchStats {
    /// Times the stall gate fired and jumped to the best open node.
    pub(crate) best_open_jumps: usize,
    /// DFS nodes popped from Q.
    pub(crate) nodes_visited: usize,
    /// Open nodes discarded by `backtrack_best_open`'s jumps (issue #158).
    pub(crate) nodes_abandoned_by_jump: usize,
    /// Nodes popped with at least one undo mark above the live stack size.
    pub(crate) mark_overshoots: usize,
    /// Nodes whose bound branch actually moved the incumbent point (#131).
    pub(crate) bound_branch_moves: usize,
}

/// mip-heuristics (36f7d9d): repair_search.h::RepairBranch (lines 101-106)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct RepairBranch {
    pub(crate) var: HighsInt,
    pub(crate) val: f64,
    pub(crate) is_fix: bool,
    pub(crate) is_lb: bool,
}

/// mip-heuristics (36f7d9d): repair_search.cpp::sync_changes (lines 58-138)
pub(crate) fn sync_changes(e: &mut PropEngine, r: &PropEngine) -> bool {
    let mut any_seeded = false;
    let feastol = e.feastol();
    for j in 0..e.ncol() {
        let rv = r.var(j);
        // Copy E's state for j before any mutable borrow of E below: the
        // C++ reads `ev` through a const ref while mutating E.
        let ev = *e.var(j);
        // R's effective domain for j: a single point if R has fixed it
        // (regardless of whether R's own vs_[j].lb/ub were narrowed --
        // R.fix() does not narrow them, so the raw bounds can understate
        // what R has actually decided).
        let r_lb = if rv.fixed { rv.val } else { rv.lb };
        let r_ub = if rv.fixed { rv.val } else { rv.ub };

        if ev.fixed {
            if ev.val >= r_lb - feastol && ev.val <= r_ub + feastol {
                continue; // domains agree -- leave as is
            }
            let target = if ev.val < r_lb { r_lb } else { r_ub };
            if !e.refix(j, target) {
                return false;
            }
            e.seed_worklist(j);
            any_seeded = true;
            continue;
        }

        // E unfixed.  Dr ⊇ D means R hasn't restricted this column at
        // all relative to what E already knows -- nothing to sync.
        if r_lb <= ev.lb + feastol && r_ub >= ev.ub - feastol {
            continue;
        }

        let new_lb = cpp_max(ev.lb, r_lb);
        let new_ub = cpp_min(ev.ub, r_ub);
        if new_lb <= new_ub + feastol {
            // Case 1: non-empty intersection -- tighten D to D ∩ Dr.
            // `new_lb`/`new_ub` are within E's current [ev.lb, ev.ub] by
            // construction, so `fix()`'s current-domain validation is the
            // right (and sufficient) check here.
            if new_ub - new_lb < feastol {
                if !e.fix(j, (new_lb + new_ub) * 0.5) {
                    return false;
                }
                e.seed_worklist(j);
            } else {
                if new_lb > ev.lb + feastol && !e.tighten_lb(j, new_lb) {
                    return false;
                }
                if new_ub < ev.ub - feastol && !e.tighten_ub(j, new_ub) {
                    return false;
                }
            }
            any_seeded = true;
        } else {
            // Case 2: disjoint intervals -- fix to the endpoint of Dr
            // closer to D.  The target is outside E's current domain by
            // definition, hence `refix`.
            let target = if r_lb > ev.ub { r_lb } else { r_ub };
            if !e.refix(j, target) {
                return false;
            }
            e.seed_worklist(j);
            any_seeded = true;
        }
    }
    if any_seeded {
        // Budget exhaustion (issue #127) and deadline expiry (#151) are
        // both sound-but-incomplete, not verdicts; only a proven
        // inconsistency fails the sync.  An expiry deliberately does *not*
        // propagate outward through this `bool`: the only thing a `false`
        // can mean to `repair_search` is "prune this node", which is the
        // exact misreading #127 and #151 exist to prevent.  Promptness is
        // already covered one level up -- the RepairSearch loop polls the
        // same deadline on every node, so an expiry costs at most the
        // remainder of the node it fired in.
        if e.propagate(-1) == PropResult::kInfeasible {
            return false;
        }
    }
    true
}

/// mip-heuristics (36f7d9d): repair_search.cpp::move_to_disjunction (lines 208-243)
pub(crate) fn move_to_disjunction(
    e: &PropEngine,
    r: &PropEngine,
    var: HighsInt,
    cur_val: f64,
    move_val: f64,
) -> (RepairBranch, RepairBranch) {
    let ev = e.var(var);

    // Binary: fix to move_val vs fix to 1-move_val
    if ev.lb == 0.0 && ev.ub == 1.0 && e.is_int(var) {
        let alt = if move_val < 0.5 { 1.0 } else { 0.0 };
        return (
            RepairBranch {
                var,
                val: move_val,
                is_fix: true,
                is_lb: false,
            },
            RepairBranch {
                var,
                val: alt,
                is_fix: true,
                is_lb: false,
            },
        );
    }

    // D, the node's current interval for this column, and Dr = [c, d].
    let d_lb = if ev.fixed { ev.val } else { ev.lb };
    let d_ub = if ev.fixed { ev.val } else { ev.ub };
    let rv = r.var(var);
    let c = if rv.fixed { rv.val } else { rv.lb };
    let d = if rv.fixed { rv.val } else { rv.ub };

    // [a, b]: D translated by the repair move's shift, clipped into
    // [c, d].  Clipping is monotone, so a <= b survives it.
    let shift = move_val - cur_val;
    let a = cpp_max(c, cpp_min(d, d_lb + shift));
    let b = cpp_max(c, cpp_min(d, d_ub + shift));

    let left_gap = a - c;
    let right_gap = d - b;

    if left_gap <= right_gap {
        // "if l <= r, we impose x_j <= b \/ x_j >= b"
        return (
            RepairBranch {
                var,
                val: b,
                is_fix: false,
                is_lb: false,
            }, // tighten_ub(b) — holds [a, b]
            RepairBranch {
                var,
                val: b,
                is_fix: false,
                is_lb: true,
            }, // tighten_lb(b) — the larger gap
        );
    }
    // "otherwise x_j <= a \/ x_j >= a"
    (
        RepairBranch {
            var,
            val: a,
            is_fix: false,
            is_lb: true,
        }, // tighten_lb(a) — holds [a, b]
        RepairBranch {
            var,
            val: a,
            is_fix: false,
            is_lb: false,
        }, // tighten_ub(a) — the larger gap
    )
}

/// mip-heuristics (36f7d9d): repair_search.cpp::backtrack_best_open (lines 251-317)
pub(crate) fn backtrack_best_open(q: &mut Vec<RepairSearchNode>) -> usize {
    if q.is_empty() {
        return 0;
    }
    let mut best = 0usize;
    for i in 1..q.len() {
        if q[i].violation < q[best].violation {
            best = i;
        }
    }
    let p = best;
    // By value, not by reference: the compaction loop below writes into
    // `Q` (never at index `p`, since `out` starts past it), and a
    // reference into a container being rewritten is the kind of aliasing
    // that survives review and not a refactor.
    let ref_node = q[p];

    // A node is *deeper* than the promoted node B (`ref`) iff any one of
    // its seven undo marks exceeds B's.  Its state was recorded on a
    // trail that a backtrack to B's marks unwinds away, so after the jump
    // its marks name a prefix that no longer exists: popping it either
    // grows a stack through `resize` or replays an unrelated one.
    // `e_pq_mark` is -1 on every node when E has no domain PQ, so it
    // never discriminates in that case.
    let deeper_than_best = |n: &RepairSearchNode| -> bool {
        n.e_vs_mark > ref_node.e_vs_mark
            || n.e_sol_mark > ref_node.e_sol_mark
            || n.e_pq_mark > ref_node.e_pq_mark
            || n.r_vs_mark > ref_node.r_vs_mark
            || n.r_sol_mark > ref_node.r_sol_mark
            || n.sol_undo_mark > ref_node.sol_undo_mark
            || n.lhs_undo_mark > ref_node.lhs_undo_mark
    };

    // The scan covers the **whole** of Q, not the suffix after `p`, and
    // that is load-bearing rather than defensive.  Q is sorted while the
    // invariant holds, which would make the prefix vacuous -- but this
    // function is the only thing that breaks the sort, and it runs many
    // times in one search.  A previous jump swapped its own B to the
    // back and left that B's deeper former neighbours to its left; once
    // B is popped and its children pushed at the new, lower live marks,
    // those older nodes sit *before* shallower ones.  The next jump can
    // then promote a node with deeper nodes on its left, and a
    // suffix-only scan leaves exactly the node that later pops into the
    // `resize`-grows path this issue is about.  Scanning all of Q is what
    // makes the post-condition below hold rather than merely hold on the
    // fixtures.
    let mut out = 0usize;
    let mut new_p = 0usize;
    for i in 0..q.len() {
        if i != p && deeper_than_best(&q[i]) {
            continue;
        }
        if i == p {
            new_p = out;
        }
        q[out] = q[i];
        out += 1;
    }
    let dropped = q.len() - out;
    q.truncate(out);

    // Safe after the truncation, and it preserves sortedness.  Q is
    // sorted on entry (this function is the only thing that can break
    // that, and the full scan above is what makes it restore it, so the
    // property is inductive), so a survivor sitting *after* B is both
    // >= B by sortedness and <= B by the filter -- it carries B's own
    // marks.  Moving B past that tail is therefore mark-neutral, and the
    // survivors before B are shallower and stay put.
    let last = q.len() - 1;
    q.swap(new_p, last);
    dropped
}

/// mip-heuristics (36f7d9d): repair_search.cpp::apply_branch_to_r (lines 328-352)
fn apply_branch_to_r(r: &mut PropEngine, node: &RepairSearchNode) -> bool {
    if node.var < 0 {
        return true; // root node — no branch
    }
    let mut ok;
    let mut pr = PropResult::kFixpoint;
    if node.is_fix {
        ok = r.fix(node.var, node.val);
        if ok {
            pr = r.propagate(node.var);
        }
    } else if node.is_lb {
        ok = r.tighten_lb(node.var, node.val);
        if ok {
            pr = r.propagate(-1);
        }
    } else {
        ok = r.tighten_ub(node.var, node.val);
        if ok {
            pr = r.propagate(-1);
        }
    }
    ok = ok && pr != PropResult::kInfeasible;
    ok
}

/// mip-heuristics (36f7d9d): repair_search.cpp::repair_search (lines 362-785)
pub(crate) fn repair_search<'p>(
    e: &mut PropEngine<'p>,
    solution: &mut Vec<f64>,
    lhs_cache: &mut Vec<f64>,
    col_lb: &'p [f64],
    col_ub: &'p [f64],
    row_lo: &'p [f64],
    row_hi: &'p [f64],
    repair_iterations: HighsInt,
    repair_noise: f64,
    repair_track_best: bool,
    rng: &mut S5Rng,
    effort_out: &mut usize,
    scratch: &mut FprScratch<'p>,
    deadline: &Deadline,
    stats: Option<&mut RepairSearchStats>,
    progress_threshold: HighsInt,
) -> bool {
    // Counters go somewhere unconditionally; `local_stats` is the sink when
    // the caller does not want them, which keeps the two increment sites
    // below free of null checks.
    let mut local_stats = RepairSearchStats::default();
    let st: &mut RepairSearchStats = match stats {
        Some(s) => s,
        None => &mut local_stats,
    };
    let ncol = e.ncol();
    let nrow = e.nrow();
    let feastol = e.feastol();
    let csc_start = e.csc_start();
    let csc_row = e.csc_row();
    let csc_val = e.csc_val();

    // --- Violated set (reuses WalkSatScratch buffers; never runs concurrently
    // with walksat_repair — the two functions are alternative Phase 3 paths
    // at the same call site in fpr_attempt). ---
    let violated = &mut scratch.walksat.violated;
    let violated_pos = &mut scratch.walksat.violated_pos;
    violated.clear();
    if violated_pos.len() < nrow as usize {
        violated_pos.clear();
        violated_pos.resize(nrow as usize, -1);
    } else {
        for slot in violated_pos.iter_mut().take(nrow as usize) {
            *slot = -1;
        }
    }
    if violated.capacity() < nrow as usize {
        violated.reserve(nrow as usize);
    }

    // --- Initialize total violation and violated set ---
    let mut total_viol = 0.0f64;
    for i in 0..nrow {
        total_viol += row_violation(
            lhs_cache[i as usize],
            row_lo[i as usize],
            row_hi[i as usize],
        );
    }

    // rebuild_violated (C++ lambda, lines 420-430)
    {
        let violated = &mut scratch.walksat.violated;
        let violated_pos = &mut scratch.walksat.violated_pos;
        for &vi in violated.iter() {
            violated_pos[vi as usize] = -1;
        }
        violated.clear();
        for i in 0..nrow {
            let lhs = lhs_cache[i as usize];
            if lhs > row_hi[i as usize] + feastol || lhs < row_lo[i as usize] - feastol {
                if violated_pos[i as usize] != -1 {
                    continue;
                }
                violated_pos[i as usize] = violated.len() as HighsInt;
                violated.push(i);
            }
        }
    }
    if scratch.walksat.violated.is_empty() {
        *effort_out = 0;
        return true;
    }

    // --- Solution/LHS undo stacks (reuse WalkSatScratch buffers) ---
    let sol_undo = &mut scratch.walksat.sol_undo;
    let lhs_undo = &mut scratch.walksat.lhs_undo;
    sol_undo.clear();
    lhs_undo.clear();

    let mut best_viol = total_viol;
    let best_solution = &mut scratch.repair_best_solution;
    let best_lhs = &mut scratch.repair_best_lhs;
    if repair_track_best {
        best_solution.clear();
        best_solution.extend(solution.iter().copied());
        best_lhs.clear();
        best_lhs.extend(lhs_cache.iter().copied());
    } else {
        best_solution.clear();
        best_lhs.clear();
    }

    // --- Secondary engine R from global bounds ---
    // Reused across calls via FprScratch to eliminate the per-call allocation
    // (vs_/solution_/prop_in_wl_/undo reserves) that mirrors the primary
    // engine's reuse pattern.  Guard: R is always constructed from E's
    // problem-data pointers, so any cached R built against a prior E is
    // safe to reuse iff every pointer E exposes still matches.  Since E
    // is itself pooled with a pointer-identity guard in fpr_attempt, the
    // only way to invalidate R mid-solve is to swap E's underlying
    // problem data — which implies fpr_attempt already re-emplaced E;
    // in that case R's guard below also mismatches and R is re-emplaced.
    let engine_r_opt = &mut scratch.repair_prop_engine_r;
    let r_valid = match engine_r_opt.as_ref() {
        Some(r) => {
            r.ncol() == ncol
                && r.nrow() == nrow
                && std::ptr::eq(r.ar_start(), e.ar_start())
                && std::ptr::eq(r.ar_index(), e.ar_index())
                && std::ptr::eq(r.ar_value(), e.ar_value())
                && std::ptr::eq(r.csc_start(), e.csc_start())
                && std::ptr::eq(r.csc_row(), e.csc_row())
                && std::ptr::eq(r.csc_val(), e.csc_val())
                && std::ptr::eq(r.col_lb(), col_lb)
                && std::ptr::eq(r.col_ub(), col_ub)
                && std::ptr::eq(r.row_lo(), row_lo)
                && std::ptr::eq(r.row_hi(), row_hi)
                && std::ptr::eq(r.integrality(), e.integrality())
                && r.feastol() == feastol
        }
        None => false,
    };
    if r_valid {
        engine_r_opt.as_mut().unwrap().reset();
    } else {
        *engine_r_opt = Some(PropEngine::new(
            ncol,
            nrow,
            e.ar_start(),
            e.ar_index(),
            e.ar_value(),
            e.csc_start(),
            e.csc_row(),
            e.csc_val(),
            col_lb,
            col_ub,
            row_lo,
            row_hi,
            e.integrality(),
            feastol,
        ));
    }
    // NOLINT rationale: `R` is the paper's own symbol for this object
    // (Fig. 5), used under that name in the surrounding prose in
    // fpr_core.h, repair_search.h and fpr_strategies.h, and — for the
    // primary engine — as the parameter name of repair_search() itself.
    // Lower-casing only the locals would split one symbol across two
    // spellings; renaming the documentation too would cost the mapping
    // back to the paper, which the standards rank above naming.
    let r = engine_r_opt.as_mut().unwrap();

    // Arm both engines' internal fixpoint polls (issue #151) with this
    // call's deadline, on the reuse path as much as on the construction
    // path.  E is armed too, not only R: `sync_changes` propagates on E,
    // and a standalone caller (a unit test, or any future one that does
    // not come through `fpr_attempt_begin`'s `acquire_engine`) would
    // otherwise leave that half of every node unbounded.
    e.set_deadline(deadline);
    r.set_deadline(deadline);

    let mut total_effort: usize = 0;
    let e_effort_baseline = e.effort();
    let r_effort_baseline = r.effort();
    let mut nodes_without_progress: HighsInt = 0;

    // --- DFS stack (paper Fig. 5, lines 3-4).  Reused across calls via
    // scratch to avoid per-call allocations. ---
    // NOLINT rationale: `Q` is the paper's own symbol for this object
    // (Fig. 5), used under that name in the surrounding prose in
    // fpr_core.h, repair_search.h and fpr_strategies.h, and — for the
    // primary engine — as the parameter name of repair_search() itself.
    // Lower-casing only the locals would split one symbol across two
    // spellings; renaming the documentation too would cost the mapping
    // back to the paper, which the standards rank above naming.
    let q = &mut scratch.repair_dfs_stack;
    q.clear();
    let root_e_pq = if e.pq_initialized() { e.pq_mark() } else { -1 };
    q.push(RepairSearchNode {
        var: -1,
        val: 0.0,
        is_fix: true,
        is_lb: false,
        e_vs_mark: e.vs_mark(),
        e_sol_mark: e.sol_mark(),
        e_pq_mark: root_e_pq,
        r_vs_mark: r.vs_mark(),
        r_sol_mark: r.sol_mark(),
        sol_undo_mark: 0,
        lhs_undo_mark: 0,
        violation: total_viol,
    });

    let mut nodes_visited: HighsInt = 0;

    // What this call charges its caller: the WalkSAT move counters plus
    // both engines' propagation effort, which is where the bulk of a node
    // goes (one AC-3 fixpoint on `R` in apply_branch_to_r, one on `E` in
    // sync_changes).  Reported only -- it is no longer also a gate (issue
    // #156): the node limit is the effort governor here, since each of a
    // node's two fixpoints already answers to `kPropagateBudgetPerNnz`.
    // effort_spent (C++ lambda, lines 570-572) is inlined at its two call
    // sites below.

    // The wall clock joins the node gate (issue #117).  Polled on every
    // node rather than on a cadence: `repair_iterations` is 50, so fifty
    // clock reads is the whole cost, while one node is two propagation
    // fixpoints — the reason this loop can outlive a time limit at all.
    while !q.is_empty() && nodes_visited < repair_iterations && !deadline.expired() {
        let node = q[q.len() - 1];
        q.pop();
        nodes_visited += 1;
        st.nodes_visited += 1;

        // The mark-ordering invariant, observed rather than asserted
        // (issue #158): a node may only ever restore *downward*, so
        // every one of its seven marks must sit at or below the
        // corresponding live stack.  A node that fails this restores to
        // a state the trail no longer holds, and both backtracks answer
        // an over-large target with a `resize` that grows the stack into
        // value-initialized entries -- see `RepairSearchStats::
        // mark_overshoots`.  Checked in Release too, because that is
        // where the tests run; the cost is seven integer compares
        // against a node worth two propagation fixpoints.
        if node.e_vs_mark > e.vs_mark()
            || node.e_sol_mark > e.sol_mark()
            || (e.pq_initialized() && node.e_pq_mark > e.pq_mark())
            || node.r_vs_mark > r.vs_mark()
            || node.r_sol_mark > r.sol_mark()
            || (node.sol_undo_mark as i64) > scratch.walksat.sol_undo.len() as i64
            || (node.lhs_undo_mark as i64) > scratch.walksat.lhs_undo.len() as i64
        {
            st.mark_overshoots += 1;
        }

        // Restore parent state (paper lines 7-8).  Pass `node.e_pq_mark`
        // explicitly: when E was `init_domain_pq`'d in Phase 2 (any
        // dynamic-var strategy), omitting the PQ undo target here leaves
        // E's heap inconsistent with vs_ across backtracks — a subsequent
        // `pq_notify` triggered by the new branch's E.fix/tighten then
        // erases a var that's no longer in the heap.  R has no PQ active
        // (repair_search never calls init_domain_pq on R), so the default
        // -1 there is harmless.
        // backtrack_sol_lhs (C++ lambda, lines 475-484)
        {
            let sol_undo = &mut scratch.walksat.sol_undo;
            let lhs_undo = &mut scratch.walksat.lhs_undo;
            let s_mark = node.sol_undo_mark;
            let l_mark = node.lhs_undo_mark;
            let mut u = sol_undo.len() as HighsInt - 1;
            while u >= s_mark {
                solution[sol_undo[u as usize].idx as usize] = sol_undo[u as usize].old_val;
                u -= 1;
            }
            // C++ `resize`: grows with value-initialized entries on a mark
            // overshoot (see `RepairSearchStats::mark_overshoots`).
            sol_undo.resize(
                s_mark as usize,
                crate::mip::fpr_local_mip::WalkSatUndoEntry::default(),
            );
            let mut u = lhs_undo.len() as HighsInt - 1;
            while u >= l_mark {
                lhs_cache[lhs_undo[u as usize].idx as usize] = lhs_undo[u as usize].old_val;
                u -= 1;
            }
            lhs_undo.resize(
                l_mark as usize,
                crate::mip::fpr_local_mip::WalkSatUndoEntry::default(),
            );
        }
        e.backtrack_to(node.e_vs_mark, node.e_sol_mark, -1, node.e_pq_mark);
        r.backtrack_to(node.r_vs_mark, node.r_sol_mark, -1, -1);
        total_viol = node.violation;
        // rebuild_violated (C++ lambda, lines 420-430)
        {
            let violated = &mut scratch.walksat.violated;
            let violated_pos = &mut scratch.walksat.violated_pos;
            for &vi in violated.iter() {
                violated_pos[vi as usize] = -1;
            }
            violated.clear();
            for i in 0..nrow {
                let lhs = lhs_cache[i as usize];
                if lhs > row_hi[i as usize] + feastol || lhs < row_lo[i as usize] - feastol {
                    if violated_pos[i as usize] != -1 {
                        continue;
                    }
                    violated_pos[i as usize] = violated.len() as HighsInt;
                    violated.push(i);
                }
            }
        }

        // Apply branch to R, propagate (paper lines 8-10)
        if !apply_branch_to_r(r, &node) {
            continue; // infeasible — prune (paper lines 11-12)
        }

        // SyncChanges R→E (paper line 13)
        if !sync_changes(e, r) {
            continue; // E infeasible after sync
        }

        // Apply branch to solution/lhs (our extension for complete-assignment)
        // apply_move updates total_viol and violated set incrementally.
        if node.var >= 0 {
            if node.is_fix {
                // apply_move (C++ lambda, lines 451-473)
                {
                    let var = node.var;
                    let new_val = node.val;
                    let old_val = solution[var as usize];
                    scratch
                        .walksat
                        .sol_undo
                        .push(crate::mip::fpr_local_mip::WalkSatUndoEntry { idx: var, old_val });
                    solution[var as usize] = new_val;
                    let delta = new_val - old_val;
                    total_effort +=
                        (csc_start[var as usize + 1] - csc_start[var as usize]) as usize;
                    let mut p = csc_start[var as usize];
                    while p < csc_start[var as usize + 1] {
                        let row = csc_row[p as usize];
                        let old_v = row_violation(
                            lhs_cache[row as usize],
                            row_lo[row as usize],
                            row_hi[row as usize],
                        );
                        scratch.walksat.lhs_undo.push(
                            crate::mip::fpr_local_mip::WalkSatUndoEntry {
                                idx: row,
                                old_val: lhs_cache[row as usize],
                            },
                        );
                        lhs_cache[row as usize] += csc_val[p as usize] * delta;
                        let new_v = row_violation(
                            lhs_cache[row as usize],
                            row_lo[row as usize],
                            row_hi[row as usize],
                        );
                        total_viol += new_v - old_v;

                        let was = scratch.walksat.violated_pos[row as usize] != -1;
                        let now = is_row_violated(
                            lhs_cache[row as usize],
                            row_lo[row as usize],
                            row_hi[row as usize],
                            feastol,
                        );
                        if was && !now {
                            // remove_violated (C++ lambda, lines 409-419)
                            let pos = scratch.walksat.violated_pos[row as usize];
                            if pos != -1 {
                                let last =
                                    scratch.walksat.violated[scratch.walksat.violated.len() - 1];
                                scratch.walksat.violated[pos as usize] = last;
                                scratch.walksat.violated_pos[last as usize] = pos;
                                scratch.walksat.violated.pop();
                                scratch.walksat.violated_pos[row as usize] = -1;
                            }
                        } else if !was && now {
                            // add_violated (C++ lambda, lines 402-408)
                            if scratch.walksat.violated_pos[row as usize] == -1 {
                                scratch.walksat.violated_pos[row as usize] =
                                    scratch.walksat.violated.len() as HighsInt;
                                scratch.walksat.violated.push(row);
                            }
                        }
                        p += 1;
                    }
                }
            } else {
                // A bound branch moves the point only through what
                // `sync_changes` transferred into E, so the clamp has to
                // read the interval E actually holds: `fix()` leaves
                // `[lb, ub]` wide on a column it fixed, so a column the
                // sync fixed inside its own bounds would otherwise clamp
                // against a domain the search no longer believes and the
                // move would be dropped (issue #131 -- the same reading
                // `move_to_disjunction` and `PropEngine::shift_domain`
                // take of a fixed column's interval).
                let bv = *e.var(node.var);
                let cur = solution[node.var as usize];
                let new_lb = if bv.fixed { bv.val } else { bv.lb };
                let new_ub = if bv.fixed { bv.val } else { bv.ub };
                let clamped = cpp_max(new_lb, cpp_min(new_ub, cur));
                if (clamped - cur).abs() > feastol {
                    // apply_move (C++ lambda, lines 451-473)
                    {
                        let var = node.var;
                        let new_val = clamped;
                        let old_val = solution[var as usize];
                        scratch.walksat.sol_undo.push(
                            crate::mip::fpr_local_mip::WalkSatUndoEntry { idx: var, old_val },
                        );
                        solution[var as usize] = new_val;
                        let delta = new_val - old_val;
                        total_effort +=
                            (csc_start[var as usize + 1] - csc_start[var as usize]) as usize;
                        let mut p = csc_start[var as usize];
                        while p < csc_start[var as usize + 1] {
                            let row = csc_row[p as usize];
                            let old_v = row_violation(
                                lhs_cache[row as usize],
                                row_lo[row as usize],
                                row_hi[row as usize],
                            );
                            scratch.walksat.lhs_undo.push(
                                crate::mip::fpr_local_mip::WalkSatUndoEntry {
                                    idx: row,
                                    old_val: lhs_cache[row as usize],
                                },
                            );
                            lhs_cache[row as usize] += csc_val[p as usize] * delta;
                            let new_v = row_violation(
                                lhs_cache[row as usize],
                                row_lo[row as usize],
                                row_hi[row as usize],
                            );
                            total_viol += new_v - old_v;

                            let was = scratch.walksat.violated_pos[row as usize] != -1;
                            let now = is_row_violated(
                                lhs_cache[row as usize],
                                row_lo[row as usize],
                                row_hi[row as usize],
                                feastol,
                            );
                            if was && !now {
                                // remove_violated (C++ lambda, lines 409-419)
                                let pos = scratch.walksat.violated_pos[row as usize];
                                if pos != -1 {
                                    let last = scratch.walksat.violated
                                        [scratch.walksat.violated.len() - 1];
                                    scratch.walksat.violated[pos as usize] = last;
                                    scratch.walksat.violated_pos[last as usize] = pos;
                                    scratch.walksat.violated.pop();
                                    scratch.walksat.violated_pos[row as usize] = -1;
                                }
                            } else if !was && now {
                                // add_violated (C++ lambda, lines 402-408)
                                if scratch.walksat.violated_pos[row as usize] == -1 {
                                    scratch.walksat.violated_pos[row as usize] =
                                        scratch.walksat.violated.len() as HighsInt;
                                    scratch.walksat.violated.push(row);
                                }
                            }
                            p += 1;
                        }
                    }
                    st.bound_branch_moves += 1;
                }
            }
        }

        if scratch.walksat.violated.is_empty() {
            // Feasible! (paper lines 15-16)
            // effort_spent (C++ lambda, lines 570-572)
            *effort_out =
                total_effort + (e.effort() - e_effort_baseline) + (r.effort() - r_effort_baseline);
            return true;
        }

        // Update best state (paper line 17)
        if total_viol < best_viol - feastol {
            best_viol = total_viol;
            if repair_track_best {
                scratch.repair_best_solution.clear();
                scratch
                    .repair_best_solution
                    .extend(solution.iter().copied());
                scratch.repair_best_lhs.clear();
                scratch.repair_best_lhs.extend(lhs_cache.iter().copied());
            }
            nodes_without_progress = 0;
        } else {
            nodes_without_progress += 1;
        }

        // FindRepairMove: WalkSAT on current solution (paper line 20)
        let pick = s5_uniform_int_i32(rng, 0, scratch.walksat.violated.len() as HighsInt - 1);
        let row = scratch.walksat.violated[pick as usize];

        let mv = walksat_select_move(
            row,
            solution,
            lhs_cache,
            col_lb,
            col_ub,
            e,
            repair_noise,
            rng,
            &mut total_effort,
            &mut scratch.walksat,
        );
        if mv.var < 0 {
            continue; // no valid move (paper lines 21-22)
        }

        // MoveToDisjunction (paper lines 24-26)
        let (preferred, alternative) =
            move_to_disjunction(e, r, mv.var, solution[mv.var as usize], mv.val);

        // Save current state marks
        let cur_e_vs = e.vs_mark();
        let cur_e_sol = e.sol_mark();
        let cur_r_vs = r.vs_mark();
        let cur_r_sol = r.sol_mark();
        let cur_sol = scratch.walksat.sol_undo.len() as HighsInt;
        let cur_lhs = scratch.walksat.lhs_undo.len() as HighsInt;

        // Push alternative first (explored second), then preferred (explored first)
        let cur_e_pq = if e.pq_initialized() { e.pq_mark() } else { -1 };
        q.push(RepairSearchNode {
            var: alternative.var,
            val: alternative.val,
            is_fix: alternative.is_fix,
            is_lb: alternative.is_lb,
            e_vs_mark: cur_e_vs,
            e_sol_mark: cur_e_sol,
            e_pq_mark: cur_e_pq,
            r_vs_mark: cur_r_vs,
            r_sol_mark: cur_r_sol,
            sol_undo_mark: cur_sol,
            lhs_undo_mark: cur_lhs,
            violation: total_viol,
        });
        q.push(RepairSearchNode {
            var: preferred.var,
            val: preferred.val,
            is_fix: preferred.is_fix,
            is_lb: preferred.is_lb,
            e_vs_mark: cur_e_vs,
            e_sol_mark: cur_e_sol,
            e_pq_mark: cur_e_pq,
            r_vs_mark: cur_r_vs,
            r_sol_mark: cur_r_sol,
            sol_undo_mark: cur_sol,
            lhs_undo_mark: cur_lhs,
            violation: total_viol,
        });

        // Check progress -- jump to the best open node if stuck (paper
        // Fig. 5 lines 18-19, Sect. 5.1: "if we detect that we are not
        // making enough progress in the current subtree, we backtrack
        // directly to the most promising open node").  This gate is the
        // only thing that changes the shape of Q (issue #130).  Until
        // #130 the same `BacktrackBestOpen` also ran unconditionally, so
        // the node popped next was the lowest-violation open node on
        // every step and the search was best-first rather than the
        // paper's DFS-with-occasional-jumps; the gate could not change
        // anything.
        //
        // **A jump gives the abandoned subtree up; it does not reorder
        // it** (issue #158).  Every open node restores its parent state
        // by replaying an undo trail down to a mark, so the marks along
        // Q have to stay non-decreasing front-to-back for a pop to
        // unwind downward -- and a permutation that merely moved the
        // best node to the back left the strictly-deeper open nodes
        // stranded above a trail the search then unwound beneath them,
        // corrupting engine state when one was later popped.  So
        // `backtrack_best_open` discards them, which is also what
        // Sect. 5.1 describes when it prices the jump: "we backtrack
        // directly to the most promising open node, at the cost of
        // giving up on completeness".  A permutation gives up no
        // completeness at all.  `st.nodes_abandoned_by_jump` counts what
        // this costs; `st.mark_overshoots`, checked at the pop site
        // above, is the invariant it buys.
        //
        // **The position is deliberate and is not the figure's line
        // order.**  Q is a LIFO stack, so a promotion made where lines
        // 18-19 sit -- before this node's two children are pushed at
        // lines 24-26 -- is immediately buried under them and resurfaces
        // only once their whole subtree has been explored: precisely the
        // "inability to escape a wrong subtree without exploring it all"
        // Sect. 5.1 gives as the jump's reason to exist, and why it says
        // to backtrack *directly*.  Promoting after the pushes is what
        // makes the promoted node the one popped next.
        //
        // Consequence worth naming: a node that reaches neither push --
        // refuted by R, refuted by the sync, or offering no repair move
        // -- never jumps.  The two refutation paths `continue` before
        // `nodes_without_progress` is even touched; a node with no move
        // increments it and leaves the jump to the next node that does
        // push.
        if nodes_without_progress >= progress_threshold {
            st.nodes_abandoned_by_jump += backtrack_best_open(q);
            nodes_without_progress = 0;
            st.best_open_jumps += 1;
        }
    }

    // Fig. 5 line 27 ends the search by backtracking to the best open node
    // and returning the state E associated with it.  This restores the
    // best state *visited* instead -- the (solution, lhs_cache) snapshot
    // line 17's improvement test took -- and that is a **choice, not a
    // necessity**: every `RepairSearchNode` carries a complete set of
    // restore marks (`e_vs_mark`, `e_sol_mark`, `e_pq_mark`, `r_vs_mark`,
    // `r_sol_mark`, `sol_undo_mark`, `lhs_undo_mark`), which is exactly
    // what the head of the loop above replays on every iteration, so line
    // 27 is implementable here verbatim.
    //
    // That last claim holds *because* the mark-ordering invariant is
    // maintained (issue #158), and not otherwise: replaying an open
    // node's marks is only meaningful while they sit at or below the
    // live stacks, which is precisely what `backtrack_best_open`'s
    // discard restores after a jump.  Under the pre-#158 permutation an
    // open node's marks could exceed them, and line 27 would then
    // restore a state the trail no longer holds.
    //
    // The reason for the choice is that the best visited state is at
    // least as good, under the paper's own line-17 measure, and needs no
    // re-walk: a node's recorded `violation` is the total violation of
    // the state that pushed it, i.e. some visited node's post-move total,
    // so every open node's violation is >= `best_viol` by construction.
    // With `repair_track_best` false no snapshot was ever taken, so
    // nothing is restored and the caller keeps whatever assignment the
    // last expanded node left behind.
    if repair_track_best && best_viol < total_viol {
        solution.clear();
        solution.extend(scratch.repair_best_solution.iter().copied());
        lhs_cache.clear();
        lhs_cache.extend(scratch.repair_best_lhs.iter().copied());
        // rebuild_violated (C++ lambda, lines 420-430)
        {
            let violated = &mut scratch.walksat.violated;
            let violated_pos = &mut scratch.walksat.violated_pos;
            for &vi in violated.iter() {
                violated_pos[vi as usize] = -1;
            }
            violated.clear();
            for i in 0..nrow {
                let lhs = lhs_cache[i as usize];
                if lhs > row_hi[i as usize] + feastol || lhs < row_lo[i as usize] - feastol {
                    if violated_pos[i as usize] != -1 {
                        continue;
                    }
                    violated_pos[i as usize] = violated.len() as HighsInt;
                    violated.push(i);
                }
            }
        }
    }

    // effort_spent (C++ lambda, lines 570-572)
    *effort_out =
        total_effort + (e.effort() - e_effort_baseline) + (r.effort() - r_effort_baseline);
    scratch.walksat.violated.is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;
    use crate::lp_data::hconst::kHighsInf;
    use crate::mip::fpr_local_mip::FprScratch;
    use crate::mip::fpr_local_mip::PropEngine;
    use crate::mip::fpr_local_mip::PropResult;
    use crate::mip::fpr_local_mip::RepairSearchNode;
    use crate::mip::fpr_local_mip::rng::S5Rng;
    use crate::mip::fpr_local_mip::rng::build_csc;
    // approx is not a dependency of this crate: the same check as
    // approx::assert_relative_eq with its defaults (epsilon = max_relative = EPSILON).
    macro_rules! assert_relative_eq {
        ($a:expr, $b:expr $(,)?) => {{
            let (a, b): (f64, f64) = ($a, $b);
            let diff = (a - b).abs();
            assert!(
                diff <= f64::EPSILON || diff <= f64::EPSILON * a.abs().max(b.abs()),
                "assertion failed: {} != {}",
                a,
                b
            );
        }};
    }

    // ---------------------------------------------------------------------
    // sync_changes / RepairSearch tests (issue #125).
    //
    // Salvagnin, Roberti, Fischetti, MPC 17:111-139, 2025, Sect. 5.1: two
    // domain-synchronization cases (`SyncChanges`, Fig. 5 line 13) that the
    // pre-#125 code got wrong.
    //
    //   1. A column already fixed in the primary engine E, but whose value
    //      the secondary engine R has since ruled out, must be re-fixed
    //      ("flipped", for a binary) to R's domain.
    //   2. A column still unfixed in E whose domain is disjoint from R's
    //      must be fixed to the endpoint of R's domain closer to E's
    //      (paper worked example: D=[1,3], Dr=[4,5] -> 4).
    //
    // `E`/`R` throughout are the paper's own symbols for the
    // primary/secondary propagation engines (Fig. 5).
    // ---------------------------------------------------------------------

    // 2 binary variables, 1 clique constraint: x0 + x1 <= 1.
    struct CliqueModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    impl CliqueModel {
        fn new() -> Self {
            let ar_start = vec![0, 2];
            let ar_index = vec![0, 1];
            let ar_value = vec![1.0, 1.0];
            let csc = build_csc(2, 1, &ar_start, &ar_index, &ar_value);
            CliqueModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0, 0.0],
                col_ub: vec![1.0, 1.0],
                row_lo: vec![-kHighsInf],
                row_hi: vec![1.0],
                integrality: vec![HighsVarType::kInteger, HighsVarType::kInteger],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                2,
                1,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                feastol,
            )
        }
    }

    #[test]
    fn s5_rsearch_sync_changes_flips_an_already_fixed_binary_via_clique_propagation() {
        const FEASTOL: f64 = 1e-6;
        let m = CliqueModel::new();

        // E: the primary engine, as a real DFS would leave it after fixing
        // x0=1 and propagating.  Clique propagation auto-fixes x1 to 0 with
        // its bounds genuinely narrowed to [0,0].
        let mut e = m.make_engine(FEASTOL);
        assert!(e.fix(0, 1.0));
        assert_eq!(e.propagate(0), PropResult::kFixpoint);
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 1.0);
        assert!(e.var(1).fixed);
        assert_relative_eq!(e.var(1).val, 0.0);
        assert_relative_eq!(e.var(1).lb, 0.0);
        assert_relative_eq!(e.var(1).ub, 0.0);

        // R: the secondary engine, after RepairSearch applies the opposite
        // branch (fix x1=1) and propagates it.  Clique propagation forces
        // x0 down to 0 in R -- the swap RepairSearch is supposed to detect.
        let mut r = m.make_engine(FEASTOL);
        assert!(r.fix(1, 1.0));
        assert_eq!(r.propagate(1), PropResult::kFixpoint);
        assert!(r.var(0).fixed);
        assert_relative_eq!(r.var(0).val, 0.0);
        assert!(r.var(1).fixed);
        assert_relative_eq!(r.var(1).val, 1.0);

        assert!(sync_changes(&mut e, &r));

        // The swap: x0 flips 1 -> 0, x1 flips 0 -> 1, both while already
        // fixed in E.
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 0.0);
        assert!(e.var(1).fixed);
        assert_relative_eq!(e.var(1).val, 1.0);
    }

    #[test]
    fn s5_rsearch_sync_changes_agreeing_domains_are_left_alone() {
        const FEASTOL: f64 = 1e-6;
        let m = CliqueModel::new();

        let mut e = m.make_engine(FEASTOL);
        assert!(e.fix(0, 1.0));
        assert_eq!(e.propagate(0), PropResult::kFixpoint);

        let mut r = m.make_engine(FEASTOL);
        assert!(r.fix(0, 1.0));
        assert_eq!(r.propagate(0), PropResult::kFixpoint);

        assert!(sync_changes(&mut e, &r));
        assert_relative_eq!(e.var(0).val, 1.0);
        assert_relative_eq!(e.var(1).val, 0.0);
    }

    // 1 general-integer variable; the row is a slack constraint that never
    // binds -- only tighten_lb/tighten_ub, called directly, drive the two
    // domains here.
    struct SlackModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    impl SlackModel {
        fn new() -> Self {
            let ar_start = vec![0, 1];
            let ar_index = vec![0];
            let ar_value = vec![1.0];
            let csc = build_csc(1, 1, &ar_start, &ar_index, &ar_value);
            SlackModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0],
                col_ub: vec![10.0],
                row_lo: vec![-kHighsInf],
                row_hi: vec![100.0],
                integrality: vec![HighsVarType::kInteger],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                1,
                1,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                feastol,
            )
        }
    }

    #[test]
    fn s5_rsearch_sync_changes_disjoint_domains_fix_to_the_nearest_endpoint_of_r_s_domain() {
        const FEASTOL: f64 = 1e-6;
        let m = SlackModel::new();

        // E: domain D = [1, 3], still unfixed.
        let mut e = m.make_engine(FEASTOL);
        assert!(e.tighten_lb(0, 1.0));
        assert!(e.tighten_ub(0, 3.0));
        assert!(!e.var(0).fixed);
        assert_relative_eq!(e.var(0).lb, 1.0);
        assert_relative_eq!(e.var(0).ub, 3.0);

        // R: domain Dr = [4, 5] -- the paper's own worked example.
        let mut r = m.make_engine(FEASTOL);
        assert!(r.tighten_lb(0, 4.0));
        assert!(r.tighten_ub(0, 5.0));
        assert!(!r.var(0).fixed);
        assert_relative_eq!(r.var(0).lb, 4.0);
        assert_relative_eq!(r.var(0).ub, 5.0);

        // The fix: sync_changes must succeed and fix to 4, the endpoint of
        // Dr closer to D.
        assert!(sync_changes(&mut e, &r));
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 4.0);
    }

    // Two continuous columns, two rows: x - 2y >= 0 and y - 2x >= 0.
    // Feasible at x = y = 0, but each row visit only halves the other
    // column's upper bound, so a fixpoint takes ~log2(ub / feastol) row
    // visits.  (Same construction as the halving model in
    // tests/test_prop_engine.cpp.)
    struct HalvingModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    impl HalvingModel {
        fn new() -> Self {
            let ar_start = vec![0, 2, 4];
            let ar_index = vec![0, 1, 1, 0];
            let ar_value = vec![1.0, -2.0, 1.0, -2.0];
            let csc = build_csc(2, 2, &ar_start, &ar_index, &ar_value);
            HalvingModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0, 0.0],
                col_ub: vec![1e300, 1e300],
                row_lo: vec![0.0, 0.0],
                row_hi: vec![kHighsInf, kHighsInf],
                integrality: vec![HighsVarType::kContinuous, HighsVarType::kContinuous],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                2,
                2,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                feastol,
            )
        }
    }

    #[test]
    fn s5_rsearch_sync_changes_an_expired_deadline_truncates_propagation_without_pruning() {
        const FEASTOL: f64 = 1e-6;
        let m = HalvingModel::new();

        // R restricts column 0, which is what gives `sync_changes` something
        // to transfer into E and therefore something to propagate.
        let mut r = m.make_engine(FEASTOL);
        assert!(r.tighten_ub(0, 1e299));

        // Control: with no deadline the sync succeeds having spent the whole
        // per-call work budget (100 * nnz = 400 counted accesses).
        let mut e_free = m.make_engine(FEASTOL);
        assert!(sync_changes(&mut e_free, &r));
        assert!(e_free.effort() > 400);

        // Armed with a limit already behind the timer's origin: no dependence
        // on real elapsed time, and expired on the fixpoint's first poll.
        // (S5: `Deadline` is ported with `timer` always `None`, so
        // `expired()` is always false -- Q4.  A deadline armed with a
        // negative limit is therefore indistinguishable from a never-expiring
        // one here, and the C++'s effort assertions below cannot be
        // reproduced: the sync runs to its work budget exactly as the control
        // does.  The faithful check is that the sync still succeeds and the
        // resulting domain is a legitimate partial fixpoint.)
        let mut e = m.make_engine(FEASTOL);
        e.set_deadline(&Deadline::default());

        assert!(sync_changes(&mut e, &r));
        assert!(e.effort() > 0);
        // Sound: the domain E was left with is a legitimate partial fixpoint.
        assert!(e.var(0).lb <= e.var(0).ub + FEASTOL);
        assert!(e.var(1).lb <= e.var(1).ub + FEASTOL);
    }

    // Same clique model as the sync_changes swap test, plus a second row
    // that only x1 can satisfy (x1 >= 1), so the starting assignment
    // (x0=1, x1=0) violates exactly one row and the only useful repair move
    // is to raise x1.
    struct CliqueSwapModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    impl CliqueSwapModel {
        fn new() -> Self {
            let ar_start = vec![0, 2, 3];
            let ar_index = vec![0, 1, 1];
            let ar_value = vec![1.0, 1.0, 1.0];
            let csc = build_csc(2, 2, &ar_start, &ar_index, &ar_value);
            CliqueSwapModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0, 0.0],
                col_ub: vec![1.0, 1.0],
                row_lo: vec![-kHighsInf, 1.0],
                row_hi: vec![1.0, kHighsInf],
                integrality: vec![HighsVarType::kInteger, HighsVarType::kInteger],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                2,
                2,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                feastol,
            )
        }
    }

    #[test]
    fn s5_rsearch_repair_search_sync_changes_flip_is_visible_in_e_after_a_full_search() {
        const FEASTOL: f64 = 1e-6;
        let m = CliqueSwapModel::new();

        // E: x0 decision-fixed to 1 (wide bounds, as a DFS `fix()` call
        // leaves them); x1 left unfixed.  Deliberately *not* propagated
        // here, so x1 stays open for RepairSearch's own move selection.
        let mut e = m.make_engine(FEASTOL);
        assert!(e.fix(0, 1.0));

        let mut solution = vec![1.0, 0.0];
        let mut lhs_cache = vec![1.0, 0.0]; // row0: 1*1+1*0=1; row1: 1*0=0

        let mut scratch = FprScratch::default();
        let mut rng = S5Rng::new(42);
        let deadline = Deadline::default(); // never expires
        let mut effort_out: usize = 0;

        let feasible = repair_search(
            &mut e,
            &mut solution,
            &mut lhs_cache,
            &m.col_lb,
            &m.col_ub,
            &m.row_lo,
            &m.row_hi,
            50,
            0.75,
            true,
            &mut rng,
            &mut effort_out,
            &mut scratch,
            &deadline,
            None,
            K_REPAIR_PROGRESS_THRESHOLD,
        );

        assert!(feasible);
        // Both rows genuinely satisfied.
        assert!(solution[0] + solution[1] <= 1.0 + FEASTOL);
        assert!(solution[1] >= 1.0 - FEASTOL);
        // The swap happened inside E: x0, decided at 1, was flipped to 0 by
        // sync_changes once the clique forced it.
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 0.0);
    }

    #[test]
    fn s5_rsearch_repair_search_e_and_r_are_armed_with_the_call_s_deadline() {
        const FEASTOL: f64 = 1e-6;
        let m = CliqueSwapModel::new();

        let mut e = m.make_engine(FEASTOL);
        assert!(e.fix(0, 1.0));
        // Un-armed going in, which is the state a production E is in only if
        // `acquire_engine` failed to arm it.
        assert!(e.deadline().timer.is_none());

        let mut solution = vec![1.0, 0.0];
        let mut lhs_cache = vec![1.0, 0.0];

        let mut scratch = FprScratch::default();
        let mut rng = S5Rng::new(42);
        // `Deadline` is ported with `timer` always `None`, so the
        // C++'s `timer`-identity assertions cannot be expressed; the `limit`
        // ones can, with a timer-less deadline carrying the C++'s 3600.0.
        let deadline = Deadline {
            timer: None,
            limit: 3600.0,
        };
        let mut effort_out: usize = 0;

        let _ = repair_search(
            &mut e,
            &mut solution,
            &mut lhs_cache,
            &m.col_lb,
            &m.col_ub,
            &m.row_lo,
            &m.row_hi,
            50,
            0.75,
            true,
            &mut rng,
            &mut effort_out,
            &mut scratch,
            &deadline,
            None,
            K_REPAIR_PROGRESS_THRESHOLD,
        );

        assert_eq!(e.deadline().limit, 3600.0);
        let r_engine = scratch.repair_prop_engine_r.as_ref();
        assert!(r_engine.is_some());
        assert_eq!(r_engine.unwrap().deadline().limit, 3600.0);
    }

    // ---------------------------------------------------------------------
    // The stall gate is the node loop's only steering (issue #130).
    //
    // 15 binaries, 7 rows of +/-1 coefficients, each row bounded on one
    // side.  At RNG seed 52: threshold 1 jumps out of the dead subtree and
    // reaches a feasible assignment; threshold 10^6 is a pure DFS that
    // spends all 50 nodes in the dead subtree.  The shipped default of 10
    // sits on the escaping side too.
    // ---------------------------------------------------------------------

    struct StallModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    impl StallModel {
        const NCOL: HighsInt = 15;
        const NROW: HighsInt = 7;

        fn new() -> Self {
            let ar_start = vec![0, 7, 16, 26, 33, 39, 47, 56];
            let ar_index = vec![
                4, 6, 7, 9, 11, 13, 14, 0, 1, 2, 3, 4, 8, 11, 12, 13, 0, 1, 3, 4, 5, 9, 10, 11, 12,
                13, 0, 1, 2, 4, 5, 13, 14, 4, 6, 8, 10, 12, 13, 0, 1, 5, 6, 7, 8, 10, 14, 0, 3, 5,
                6, 7, 8, 10, 11, 12,
            ];
            let ar_value = vec![
                1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0,
                -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0,
                1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0,
                1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0,
            ];
            let csc = build_csc(Self::NCOL, Self::NROW, &ar_start, &ar_index, &ar_value);
            StallModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0; Self::NCOL as usize],
                col_ub: vec![1.0; Self::NCOL as usize],
                row_lo: vec![1.0, 0.0, -kHighsInf, -kHighsInf, 2.0, -kHighsInf, 0.0],
                row_hi: vec![kHighsInf, kHighsInf, -2.0, -2.0, kHighsInf, -3.0, kHighsInf],
                integrality: vec![HighsVarType::kInteger; Self::NCOL as usize],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                Self::NCOL,
                Self::NROW,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                feastol,
            )
        }
    }

    struct StallRun {
        feasible: bool,
        stats: RepairSearchStats,
        solution: Vec<f64>,
        lhs_cache: Vec<f64>,
    }

    // One `repair_search` run on the model above, from the all-zero
    // assignment, with a fresh engine and a fresh RNG so the two runs differ
    // in the threshold and in nothing else.
    fn run_stall_model(m: &StallModel, progress_threshold: HighsInt, seed: u64) -> StallRun {
        const FEASTOL: f64 = 1e-6;
        let mut e = m.make_engine(FEASTOL);
        let mut solution = vec![0.0; StallModel::NCOL as usize];
        let mut lhs_cache = vec![0.0; StallModel::NROW as usize];
        let mut scratch = FprScratch::default();
        let mut rng = S5Rng::new(seed);
        let deadline = Deadline::default(); // never expires
        let mut effort_out: usize = 0;
        let mut stats = RepairSearchStats::default();
        let feasible = repair_search(
            &mut e,
            &mut solution,
            &mut lhs_cache,
            &m.col_lb,
            &m.col_ub,
            &m.row_lo,
            &m.row_hi,
            50,
            0.75,
            true,
            &mut rng,
            &mut effort_out,
            &mut scratch,
            &deadline,
            Some(&mut stats),
            progress_threshold,
        );
        StallRun {
            feasible,
            stats,
            solution,
            lhs_cache,
        }
    }

    #[test]
    fn s5_rsearch_repair_search_the_progress_threshold_decides_the_search() {
        const SEED: u64 = 52;
        const FEASTOL: f64 = 1e-6;
        let m = StallModel::new();

        // Out of reach: 10^6 exceeds `repair_iterations`, so the gate cannot
        // fire even once and the node loop is a pure DFS.
        let dfs = run_stall_model(&m, 1000000, SEED);
        assert_eq!(dfs.stats.best_open_jumps, 0);
        assert_eq!(dfs.stats.nodes_visited, 50); // the whole node budget, in one dead subtree
        assert!(!dfs.feasible);

        // Threshold 1: abandon the subtree at the first node that fails to
        // improve the best violation.
        let jumping = run_stall_model(&m, 1, SEED);
        assert!(jumping.stats.best_open_jumps > 0);
        assert!(jumping.stats.nodes_visited < dfs.stats.nodes_visited);
        assert!(jumping.feasible);

        // The value itself, which mip-heuristics' PARAMETERS.md documents and nothing
        // else checks.
        assert_eq!(K_REPAIR_PROGRESS_THRESHOLD, 10);

        // The production call site names no threshold -- it takes the
        // parameter's default -- so this is what `fpr_core.cpp` actually
        // runs, and it must be `kRepairProgressThreshold`'s search and no
        // other.
        let production = run_stall_model(&m, K_REPAIR_PROGRESS_THRESHOLD, SEED);
        let named = run_stall_model(&m, K_REPAIR_PROGRESS_THRESHOLD, SEED);
        assert_eq!(production.feasible, named.feasible);
        assert_eq!(production.stats.nodes_visited, named.stats.nodes_visited);
        assert_eq!(
            production.stats.best_open_jumps,
            named.stats.best_open_jumps
        );
        assert_eq!(production.solution, named.solution);
        // Between the two ends, and strictly inside the node budget, so the
        // pin above is on a search that actually escapes.
        assert!(named.stats.best_open_jumps > 0);
        assert!(named.stats.nodes_visited < 50);
        assert!(named.feasible);

        // ... and what it returns really is feasible, checked against the
        // rows rather than trusting the return value.
        for i in 0..StallModel::NROW {
            let mut lhs = 0.0;
            let mut k = m.ar_start[i as usize];
            while k < m.ar_start[i as usize + 1] {
                lhs += m.ar_value[k as usize] * jumping.solution[m.ar_index[k as usize] as usize];
                k += 1;
            }
            assert!(lhs >= m.row_lo[i as usize] - FEASTOL, "row {}", i);
            assert!(lhs <= m.row_hi[i as usize] + FEASTOL, "row {}", i);
        }
    }

    // ---------------------------------------------------------------------
    // BacktrackBestOpen gives its subtree up (issue #158).
    // ---------------------------------------------------------------------

    // A node carrying the same value in all seven undo marks, which is the
    // shape the invariant is about; `label` is the name the cases identify
    // it by and plays no part in the drop rule.
    fn marked_node(label: HighsInt, mark: HighsInt, violation: f64) -> RepairSearchNode {
        RepairSearchNode {
            var: label,
            val: 0.0,
            is_fix: true,
            is_lb: false,
            e_vs_mark: mark,
            e_sol_mark: mark,
            e_pq_mark: mark,
            r_vs_mark: mark,
            r_sol_mark: mark,
            sol_undo_mark: mark,
            lhs_undo_mark: mark,
            violation,
        }
    }

    #[test]
    fn s5_rsearch_backtrack_best_open_a_jump_drops_the_open_nodes_beneath_it() {
        // The issue's own trace: Q = [A(0), C(m1), E(m2)] with m1 < m2 and A
        // the lowest-violation node.  After #158 the two nodes beneath A are
        // given up instead.
        let mut q = vec![
            marked_node(10, 0, 1.0),
            marked_node(11, 4, 5.0),
            marked_node(12, 9, 7.0),
        ];

        let dropped = backtrack_best_open(&mut q);

        assert_eq!(dropped, 2);
        assert_eq!(q.len(), 1);
        assert_eq!(q[0].var, 10);
        assert_relative_eq!(q[0].violation, 1.0);
    }

    #[test]
    fn s5_rsearch_backtrack_best_open_the_alt_pref_pair_at_the_best_node_s_own_marks_survives() {
        // The two children of one parent are pushed together at the same
        // state, so their marks are *equal*, not deeper.  `min_element`
        // returns the *first* minimum, so with the pair tied on violation
        // the alternative is chosen and the swap sends it to the back.
        let mut q = vec![
            marked_node(20, 0, 5.0),
            marked_node(21, 3, 1.0),
            marked_node(22, 3, 1.0),
        ];

        let dropped = backtrack_best_open(&mut q);

        assert_eq!(dropped, 0);
        assert_eq!(q.len(), 3);
        assert_eq!(q[0].var, 20);
        assert_eq!(q[1].var, 22); // pref, left where it was
        assert_eq!(q[2].var, 21); // alt, the first minimum, swapped to the back
    }

    #[test]
    fn s5_rsearch_backtrack_best_open_any_one_deeper_mark_is_enough_to_drop_a_node() {
        // The drop rule is a *componentwise* comparison over all seven
        // marks; each is checked on its own here.
        for k in 0..7usize {
            let mut best = marked_node(30, 2, 1.0);
            let mut other = marked_node(31, 2, 9.0);
            match k {
                0 => other.e_vs_mark = 3,
                1 => other.e_sol_mark = 3,
                2 => other.e_pq_mark = 3,
                3 => other.r_vs_mark = 3,
                4 => other.r_sol_mark = 3,
                5 => other.sol_undo_mark = 3,
                _ => other.lhs_undo_mark = 3,
            }

            let mut q = vec![best, other];
            let dropped = backtrack_best_open(&mut q);

            assert_eq!(dropped, 1, "mark component {}", k);
            assert_eq!(q.len(), 1, "mark component {}", k);
            assert_eq!(q[0].var, 30, "mark component {}", k);
        }
    }

    #[test]
    fn s5_rsearch_backtrack_best_open_a_deeper_node_before_the_best_one_is_dropped_too() {
        // A suffix-only scan would strand the deep leftover at 9.  Q here is
        // that state: a deep leftover at 9, the promoted node at 3 in the
        // middle, and a deeper node after it.  Both deep nodes go.
        let mut q = vec![
            marked_node(40, 9, 7.0),
            marked_node(41, 3, 1.0),
            marked_node(42, 5, 6.0),
        ];

        let dropped = backtrack_best_open(&mut q);

        assert_eq!(dropped, 2);
        assert_eq!(q.len(), 1);
        assert_eq!(q[0].var, 41);

        // And the shallower prefix survives, in place, with the promoted
        // node still swapped to the back where the LIFO pop takes it next.
        let mut q2 = vec![
            marked_node(50, 1, 7.0),
            marked_node(51, 3, 1.0),
            marked_node(52, 5, 6.0),
        ];

        let dropped2 = backtrack_best_open(&mut q2);

        assert_eq!(dropped2, 1);
        assert_eq!(q2.len(), 2);
        assert_eq!(q2[0].var, 50);
        assert_eq!(q2[1].var, 51);
    }

    #[test]
    fn s5_rsearch_repair_search_a_jump_abandons_its_subtree_instead_of_corrupting_the_trail() {
        // The whole-search half of #158, on the same `StallModel` and at the
        // fixture's own seed 42, where the defect is reachable: at
        // `progress_threshold` 1 the gate fires often enough to strand a
        // deeper open node above a rewritten trail.
        let m = StallModel::new();
        let run = run_stall_model(&m, 1, 42);

        assert!(run.stats.best_open_jumps > 0);
        assert!(run.stats.nodes_abandoned_by_jump > 0);
        assert_eq!(run.stats.mark_overshoots, 0);
    }

    // ---------------------------------------------------------------------
    // MoveToDisjunction (issue #131).
    // ---------------------------------------------------------------------

    // x0 integer in [0, 10]; one row `row_lo0 <= x0 <= row_hi0`.  The caller
    // sets up E and the incumbent point, which is what distinguishes the
    // cases.
    struct OneIntModel {
        ar_start: Vec<HighsInt>,
        ar_index: Vec<HighsInt>,
        ar_value: Vec<f64>,
        col_lb: Vec<f64>,
        col_ub: Vec<f64>,
        row_lo: Vec<f64>,
        row_hi: Vec<f64>,
        integrality: Vec<HighsVarType>,
        csc_start: Vec<HighsInt>,
        csc_row: Vec<HighsInt>,
        csc_val: Vec<f64>,
    }

    impl OneIntModel {
        fn new(lo: f64, hi: f64) -> Self {
            let ar_start = vec![0, 1];
            let ar_index = vec![0];
            let ar_value = vec![1.0];
            let csc = build_csc(1, 1, &ar_start, &ar_index, &ar_value);
            OneIntModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0],
                col_ub: vec![10.0],
                row_lo: vec![lo],
                row_hi: vec![hi],
                integrality: vec![HighsVarType::kInteger],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                1,
                1,
                &self.ar_start,
                &self.ar_index,
                &self.ar_value,
                &self.csc_start,
                &self.csc_row,
                &self.csc_val,
                &self.col_lb,
                &self.col_ub,
                &self.row_lo,
                &self.row_hi,
                &self.integrality,
                feastol,
            )
        }
    }

    #[test]
    fn s5_rsearch_move_to_disjunction_a_decision_fixed_binary_still_gets_the_flip_pair() {
        const FEASTOL: f64 = 1e-6;
        let mut m = OneIntModel::new(0.0, kHighsInf);
        m.col_ub = vec![1.0];
        let mut e = m.make_engine(FEASTOL);
        let r = m.make_engine(FEASTOL);
        assert!(e.fix(0, 0.0));

        let (preferred, alternative) = move_to_disjunction(&e, &r, 0, 0.0, 1.0);
        assert!(preferred.is_fix);
        assert_relative_eq!(preferred.val, 1.0);
        assert!(alternative.is_fix);
        assert_relative_eq!(alternative.val, 0.0);
    }

    #[test]
    fn s5_rsearch_move_to_disjunction_the_split_is_a_point_split_at_one_endpoint_of_a_b() {
        // D = [2, 4] unfixed, point 3, move to 7: s = +4, so the shifted
        // interval is [6, 8] -- width 2 -- against R's root domain [0, 10].
        // l = 6 > r = 2, so the paper splits at a = 6.
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(7.0, kHighsInf);
        let mut e = m.make_engine(FEASTOL);
        let r = m.make_engine(FEASTOL);
        assert!(e.tighten_lb(0, 2.0));
        assert!(e.tighten_ub(0, 4.0));

        let (preferred, alternative) = move_to_disjunction(&e, &r, 0, 3.0, 7.0);
        // A point split: same value, opposite senses, neither a fix.
        assert_relative_eq!(preferred.val, alternative.val);
        assert!(!preferred.is_fix);
        assert!(!alternative.is_fix);
        assert_ne!(preferred.is_lb, alternative.is_lb);
        // Split at a = D.lb + s = 6, which is neither D's own 2 nor the move
        // value 7.
        assert_relative_eq!(preferred.val, 6.0);
        // l > r, so the preferred child is the one holding [a, b]: x0 >= a.
        assert!(preferred.is_lb);
        assert!(!alternative.is_lb);
    }

    #[test]
    fn s5_rsearch_move_to_disjunction_l_le_r_splits_at_b_instead_with_the_senses_swapped() {
        // Same column, mirrored: D = [6, 8] unfixed, point 7, move to 3, so
        // s = -4 and the shifted interval is [2, 4].  Now l = 2 <= r = 6, so
        // the split is at b = 4 and the child holding [a, b] is x0 <= b.
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(0.0, 3.0);
        let mut e = m.make_engine(FEASTOL);
        let r = m.make_engine(FEASTOL);
        assert!(e.tighten_lb(0, 6.0));
        assert!(e.tighten_ub(0, 8.0));

        let (preferred, alternative) = move_to_disjunction(&e, &r, 0, 7.0, 3.0);
        assert_relative_eq!(preferred.val, alternative.val);
        assert_relative_eq!(preferred.val, 4.0);
        assert!(!preferred.is_lb);
        assert!(alternative.is_lb);
    }

    #[test]
    fn s5_rsearch_move_to_disjunction_the_shifted_interval_is_clipped_into_r_s_domain() {
        // D is the singleton {8} of a decision-fixed column and the point
        // has drifted to 0, so a move to 3 is a shift of +3 and the shifted
        // interval is {11}, outside R's [0, 10].  Clipped, the split is at
        // 10 and both children are real bound changes.
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(3.0, kHighsInf);
        let mut e = m.make_engine(FEASTOL);
        let r = m.make_engine(FEASTOL);
        assert!(e.fix(0, 8.0));

        let (preferred, alternative) = move_to_disjunction(&e, &r, 0, 0.0, 3.0);
        assert_relative_eq!(preferred.val, 10.0);
        assert_relative_eq!(alternative.val, 10.0);
        assert!(preferred.is_lb);
    }

    struct OneIntRun {
        feasible: bool,
        stats: RepairSearchStats,
        solution: Vec<f64>,
        lhs_cache: Vec<f64>,
    }

    // One `repair_search` run on `m` from the given E state and point.
    fn run_one_int<'a>(m: &'a OneIntModel, e: &mut PropEngine<'a>, point: f64) -> OneIntRun {
        let mut solution = vec![point];
        let mut lhs_cache = vec![point];
        let mut scratch = FprScratch::default();
        let mut rng = S5Rng::new(42);
        let deadline = Deadline::default(); // never expires
        let mut effort_out: usize = 0;
        let mut stats = RepairSearchStats::default();
        let feasible = repair_search(
            e,
            &mut solution,
            &mut lhs_cache,
            &m.col_lb,
            &m.col_ub,
            &m.row_lo,
            &m.row_hi,
            50,
            0.75,
            true,
            &mut rng,
            &mut effort_out,
            &mut scratch,
            &deadline,
            Some(&mut stats),
            K_REPAIR_PROGRESS_THRESHOLD,
        );
        OneIntRun {
            feasible,
            stats,
            solution,
            lhs_cache,
        }
    }

    #[test]
    fn s5_rsearch_repair_search_a_fixed_general_integer_is_moved_by_a_bound_branch() {
        // x0 decision-fixed to 2 (so `[lb, ub]` stays the wide [0, 10]) and
        // the row wants x0 >= 7.  The one move is the shift +5, the shifted
        // interval is the singleton [7, 7], l = 7 > r = 3, so the split is
        // at 7 and the preferred child `x0 >= 7` drives R to [7, 10];
        // `sync_changes` then re-fixes E, which is what moves the point.
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(7.0, kHighsInf);
        let mut e = m.make_engine(FEASTOL);
        assert!(e.fix(0, 2.0));
        assert_relative_eq!(e.var(0).lb, 0.0); // wide, as `fix()` leaves them
        assert_relative_eq!(e.var(0).ub, 10.0); // -- which is why D is the value

        let r = run_one_int(&m, &mut e, 2.0);
        assert!(r.feasible);
        assert_eq!(r.stats.bound_branch_moves, 1);
        assert_relative_eq!(r.solution[0], 7.0);
        assert_relative_eq!(r.lhs_cache[0], 7.0);
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 7.0);
    }

    #[test]
    fn s5_rsearch_repair_search_an_unfixed_general_integer_is_moved_by_a_bound_branch() {
        // The disjunction of the point-split test above, run: D = [2, 4],
        // point 3, row x0 >= 7, split at 6, preferred `x0 >= 6`.  R
        // propagates that to [7, 10], which is disjoint from D, so
        // `sync_changes` takes its disjoint case and fixes x0 to R's nearer
        // endpoint.
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(7.0, kHighsInf);
        let mut e = m.make_engine(FEASTOL);
        assert!(e.tighten_lb(0, 2.0));
        assert!(e.tighten_ub(0, 4.0));
        assert!(!e.var(0).fixed);

        let r = run_one_int(&m, &mut e, 3.0);
        assert!(r.feasible);
        assert_eq!(r.stats.bound_branch_moves, 1);
        assert_relative_eq!(r.solution[0], 7.0);
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 7.0);
    }

    #[test]
    fn s5_rsearch_repair_search_a_sync_that_fixes_inside_e_s_own_bounds_still_moves_the_point() {
        // x0 unfixed on the full [0, 10], point 0, and the row is the
        // equality x0 = 8.  The move is +8, the shifted interval [8, 18]
        // clips to [8, 10], l = 8 > r = 0, so the preferred child is
        // `x0 >= 8`; R propagates the equality to [8, 8], and that
        // intersects D in the single point 8, which is `sync_changes`'
        // intersection case -- committed through `PropEngine::fix`, which
        // leaves `[lb, ub]` at [0, 10].
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(8.0, 8.0);
        let mut e = m.make_engine(FEASTOL);
        assert!(!e.var(0).fixed);

        let r = run_one_int(&m, &mut e, 0.0);
        assert!(r.feasible);
        assert_eq!(r.stats.bound_branch_moves, 1);
        assert_relative_eq!(r.solution[0], 8.0);
        assert!(e.var(0).fixed);
        assert_relative_eq!(e.var(0).val, 8.0);
        assert_relative_eq!(e.var(0).lb, 0.0); // `fix()` left them wide ...
        assert_relative_eq!(e.var(0).ub, 10.0); // ... which is the point of this case
    }

    #[test]
    fn s5_rsearch_repair_search_the_clip_decides_which_node_lands_the_move() {
        // The clip case above, run.  x0 decision-fixed to 8, point drifted
        // to 0, row x0 >= 3: the shifted interval is {11}, outside R's
        // [0, 10].  Clipped, the split is at 10, the preferred child
        // `x0 >= 10` drives R to [10, 10], the sync re-fixes E, and the
        // second node ends the search at 10.
        const FEASTOL: f64 = 1e-6;
        let m = OneIntModel::new(3.0, kHighsInf);
        let mut e = m.make_engine(FEASTOL);
        assert!(e.fix(0, 8.0));

        let r = run_one_int(&m, &mut e, 0.0);
        assert!(r.feasible);
        assert_eq!(r.stats.bound_branch_moves, 1);
        assert_eq!(r.stats.nodes_visited, 2);
        assert_relative_eq!(r.solution[0], 10.0);
    }
}
