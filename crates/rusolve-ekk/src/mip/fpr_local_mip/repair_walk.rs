//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! mip-heuristics (36f7d9d): src/repair_walk.h + src/repair_walk.cpp (lines 1-478)
//!
//! RepairWalk on a *partial* assignment (paper Fig. 1 line 8, Sect. 5).
//! The repair the paper's framework is built around: at an infeasible node of
//! the fix-and-propagate DFS most columns are not fixed, so "the partial
//! assignment at the current node is encoded by the current domain".  A row's
//! violation is the distance between its activity range `[m_i, M_i]` and its
//! bounds; a repair move is a *shift* of a column's whole current domain,
//! clipped to the column's structural bounds, keeping its width (no domain
//! enlargement).  Every change goes through `PropEngine`, so the engine's own
//! undo stacks are the best-state snapshot.  `repair_walk` deliberately does
//! **not** propagate (Sect. 5.1: "it is not possible to do constraint
//! propagation from an infeasible state").

use crate::HighsInt;
use crate::lp_data::hconst::kHighsInf;
use crate::mip::fpr_local_mip::Deadline;
use crate::mip::fpr_local_mip::PropEngine;
use crate::mip::fpr_local_mip::RepairWalkCandidate;
use crate::mip::fpr_local_mip::RepairWalkScratch;
use crate::mip::fpr_local_mip::VarState;
use crate::mip::fpr_local_mip::kRepairWalkBudgetPerNnz;
use crate::mip::fpr_local_mip::kSoftRestartPeriod;
use crate::mip::fpr_local_mip::kTabuLength;
use crate::mip::fpr_local_mip::rng::S5Rng;
use crate::mip::fpr_local_mip::rng::s5_uniform_int_i32;
use crate::mip::fpr_local_mip::rng::s5_uniform_real;
use crate::util::highs_utils::cpp_max;
use crate::util::highs_utils::cpp_min;

/// mip-heuristics (36f7d9d): src/repair_walk.cpp::range_violation (lines 60-68)
///
/// A row's violation on a *partial* assignment (paper Sect. 5): the distance
/// between the activity range `[m, mx]` the current domain admits and the
/// row's bounds.  Zero unless the row is unsatisfiable by *every* completion.
/// Always finite despite the infinities in play.
fn range_violation(min_act: f64, max_act: f64, lo: f64, hi: f64) -> f64 {
    if max_act < lo {
        return lo - max_act;
    }
    if min_act > hi {
        return min_act - hi;
    }
    0.0
}

/// mip-heuristics (36f7d9d): src/repair_walk.cpp::range_violated (lines 70-72)
fn range_violated(min_act: f64, max_act: f64, lo: f64, hi: f64, feastol: f64) -> bool {
    max_act < lo - feastol || min_act > hi + feastol
}

/// mip-heuristics (36f7d9d): src/deadline.h::Deadline::expired (Q4: `timer`
/// is always `None` in this port, so `expired()` is always false; the poll
/// stays in the code as a literal port).
fn deadline_expired(d: &Deadline) -> bool {
    match d.timer {
        Some(never) => match never {},
        None => false,
    }
}

/// mip-heuristics (36f7d9d): src/repair_walk.cpp::repair_walk `rescan` lambda (lines 201-214)
///
/// Rebuild the violated set and the total violation from the engine's
/// activities; returns the total (the C++ lambda writes the captured
/// `total_viol`).  `add_violated` (lines 157-163) is inlined.
fn rescan(
    violated: &mut Vec<HighsInt>,
    violated_pos: &mut [HighsInt],
    min_act: &[f64],
    max_act: &[f64],
    row_lo: &[f64],
    row_hi: &[f64],
    feastol: f64,
    nrow: HighsInt,
    effort: &mut usize,
) -> f64 {
    for vi in violated.iter() {
        violated_pos[*vi as usize] = -1;
    }
    violated.clear();
    let mut total_viol = 0.0f64;
    *effort += nrow as usize;
    let mut i: HighsInt = 0;
    while i < nrow {
        total_viol += range_violation(
            min_act[i as usize],
            max_act[i as usize],
            row_lo[i as usize],
            row_hi[i as usize],
        );
        if range_violated(
            min_act[i as usize],
            max_act[i as usize],
            row_lo[i as usize],
            row_hi[i as usize],
            feastol,
        ) {
            // add_violated(i)
            if violated_pos[i as usize] == -1 {
                violated_pos[i as usize] = violated.len() as HighsInt;
                violated.push(i);
            }
        }
        i += 1;
    }
    total_viol
}

/// mip-heuristics (36f7d9d): src/repair_walk.cpp::any_violated_row_in_column (lines 76-98)
///
/// Whether applying a fixing to column `j` left one of the rows `j` appears
/// in unsatisfiable -- the infeasibility half of `Apply(fixing, P)`, paper
/// Fig. 1 line 4.  Only `j`'s own rows are scanned, because only their
/// activities moved.  `E` must have had `init_activities()` called; returns
/// false otherwise.  `effort` is incremented by the coefficient accesses
/// consumed.
pub(crate) fn any_violated_row_in_column(e: &PropEngine, j: HighsInt, effort: &mut usize) -> bool {
    if !e.activities_initialized() {
        return false;
    }
    let csc_start = e.csc_start();
    let csc_row = e.csc_row();
    let row_lo = e.row_lo();
    let row_hi = e.row_hi();
    let min_act = e.min_activity_data();
    let max_act = e.max_activity_data();
    let feastol = e.feastol();

    let cbeg = csc_start[j as usize];
    let cend = csc_start[(j + 1) as usize];
    *effort += (cend - cbeg) as usize;
    let mut p = cbeg;
    while p < cend {
        let i = csc_row[p as usize];
        if range_violated(
            min_act[i as usize],
            max_act[i as usize],
            row_lo[i as usize],
            row_hi[i as usize],
            feastol,
        ) {
            return true;
        }
        p += 1;
    }
    false
}

/// mip-heuristics (36f7d9d): src/repair_walk.cpp::repair_walk (lines 106-478)
///
/// Cognitive complexity 66 (threshold 25).  Kept whole: RepairWalk from
/// Sect. 5 in full -- the activity-range violation bookkeeping, the
/// shift/clip/damage candidate rule and the WalkSAT selection are one
/// algorithm and share the per-step state.
///
/// Returns true when no row is violated any more, i.e. the node was repaired.
/// `effort_out` is the coefficient accesses consumed.  `max_steps` is the
/// paper's per-call iteration limit; the effort valve is *internal*
/// (`kRepairWalkBudgetPerNnz`).  A deadline expiry is **not** signalled
/// outward: the only thing this `bool` can mean to its caller is "the node is
/// still infeasible".
pub(crate) fn repair_walk(
    e: &mut PropEngine,
    max_steps: HighsInt,
    noise: f64,
    rng: &mut S5Rng,
    effort_out: &mut usize,
    scratch: &mut RepairWalkScratch,
    deadline: &Deadline,
) -> bool {
    // Sect. 5 leans on activities explicitly, so an engine without them
    // cannot measure a partial assignment's violation at all.
    debug_assert!(
        e.activities_initialized(),
        "repair_walk requires PropEngine::init_activities() (fpr_attempt_begin arms it)"
    );
    if !e.activities_initialized() {
        // Release-build backstop for the assert above: without activities
        // there is nothing to read violations from.  "Still infeasible" is
        // the honest answer.
        *effort_out = 0;
        return false;
    }

    let nrow = e.nrow();
    let feastol = e.feastol();
    let row_lo = e.row_lo();
    let row_hi = e.row_hi();
    let ar_start = e.ar_start();
    let ar_index = e.ar_index();
    let ar_value = e.ar_value();
    let csc_start = e.csc_start();
    let csc_row = e.csc_row();
    let csc_val = e.csc_val();
    let col_lb = e.col_lb();
    let col_ub = e.col_ub();

    let max_effort =
        (kRepairWalkBudgetPerNnz * (ar_start[nrow as usize] as usize)) + (nrow as usize);
    let mut effort: usize = 0;

    let violated = &mut scratch.violated;
    let violated_pos = &mut scratch.violated_pos;
    if violated_pos.len() < nrow as usize {
        violated_pos.clear();
        violated_pos.resize(nrow as usize, -1);
    } else {
        for slot in violated_pos[..nrow as usize].iter_mut() {
            *slot = -1;
        }
    }
    violated.clear();
    if violated.capacity() < nrow as usize {
        violated.reserve(nrow as usize);
    }

    // add_violated / remove_violated / refresh_violated are inlined below
    // (Rust closures cannot capture `violated` and `violated_pos` together
    // while also being called from a loop that reads `min_act`/`max_act`).

    // Rebuild the violated set and the total violation from the engine's
    // activities.  Nothing upstream maintains a violated-row set, so the walk
    // has to find its starting point this way; within the walk the set is
    // then maintained incrementally, one column at a time.  O(nrow), at every
    // refuted node; a soft restart pays it again.
    //
    // The C++ `rescan` lambda reads `min_act`/`max_act` through pointers that
    // `shift_domain`/`backtrack_to` mutate in place; in Rust the slices are
    // re-fetched from `e` at each call site (identical values, no held borrow
    // across a mutation).
    let mut total_viol = rescan(
        violated,
        violated_pos,
        e.min_activity_data(),
        e.max_activity_data(),
        row_lo,
        row_hi,
        feastol,
        nrow,
        &mut effort,
    );
    if violated.is_empty() {
        *effort_out = effort;
        return true;
    }

    // Every repair move goes through `E`, so the engine's own undo stacks are
    // the best-state snapshot: a mark is four integers, and restoring it is
    // the same `backtrack_to` the DFS uses.
    #[derive(Clone, Copy)]
    struct Mark {
        vs: HighsInt,
        sol: HighsInt,
        act: HighsInt,
        pq: HighsInt,
    }
    // mark_now() inlined:
    let mark_now = |e: &PropEngine| -> Mark {
        Mark {
            vs: e.vs_mark(),
            sol: e.sol_mark(),
            act: e.act_mark(),
            pq: if e.pq_initialized() { e.pq_mark() } else { -1 },
        }
    };

    let mut best_mark = mark_now(e);
    let mut best_viol = total_viol;

    let cand = &mut scratch.cand;
    let best_indices = &mut scratch.best_indices;
    let tabu = &mut scratch.tabu;
    tabu.clear();

    let mut since_restart: HighsInt = 0;
    let mut step: HighsInt = 0;
    // stats.rs throughput counters (never read by the walk).
    let (mut ls_evals, mut ls_commits) = (0u64, 0u64);
    while step < max_steps && !violated.is_empty() {
        if effort >= max_effort {
            break;
        }
        // Polled per step, not on a cadence: a step is O(row degree x column
        // degree) and this runs inside the DFS's own node.  The expiry stops
        // the walk and nothing more.
        if deadline_expired(deadline) {
            break;
        }

        // Re-read the engine's activity slices at the top of the step (the
        // previous step's `shift_domain`/`backtrack_to` may have mutated them
        // in place; C++ reads through the same pointers).
        let min_act = e.min_activity_data();
        let max_act = e.max_activity_data();

        // Pick a violated row uniformly at random (paper Sect. 5).
        let row = violated[s5_uniform_int_i32(rng, 0, violated.len() as HighsInt - 1) as usize];

        let lo_r = row_lo[row as usize];
        let hi_r = row_hi[row as usize];
        let m_r = min_act[row as usize];
        let mx_r = max_act[row as usize];
        let cur_viol = range_violation(m_r, mx_r, lo_r, hi_r);
        // Which side the activity range sits on decides both the target bound
        // and which endpoint of the range has to reach it.
        let overshoot_above = m_r > hi_r + feastol;
        let target = if overshoot_above { hi_r } else { lo_r };
        let anchor = if overshoot_above { m_r } else { mx_r };

        cand.clear();
        let mut best_damage = f64::INFINITY;
        let kbeg = ar_start[row as usize];
        let kend = ar_start[(row + 1) as usize];
        effort += (kend - kbeg) as usize;
        let mut k = kbeg;
        while k < kend {
            let j = ar_index[k as usize];
            let a = ar_value[k as usize];
            if a.abs() < 1e-15 {
                k += 1;
                continue;
            }
            if tabu.contains(&j) {
                k += 1;
                continue;
            }

            let vj: &VarState = e.var(j);
            let cur_lb = if vj.fixed { vj.val } else { vj.lb };
            let cur_ub = if vj.fixed { vj.val } else { vj.ub };
            // An unbounded interval cannot be translated, and its activity
            // contribution would not move if it were.
            if cur_lb <= -kHighsInf || cur_ub >= kHighsInf {
                k += 1;
                continue;
            }

            // Shifting the interval by `s` translates *both* activity
            // endpoints of every row containing `j` by `a * s`, whichever sign
            // `a` has.  The minimal shift that would satisfy this row, rounded
            // (paper Sect. 5: "compute the minimal shift that would make the
            // constraint satisfied, round it if the variable is general
            // integer, and then clip it").
            let mut s = (target - anchor) / a;
            let integer = e.is_int(j);
            if integer {
                let floor_it = overshoot_above == (a > 0.0);
                s = if floor_it {
                    (s + feastol).floor()
                } else {
                    (s - feastol).ceil()
                };
            }
            // Clip so the *translated interval* still lies inside the column's
            // structural bounds -- the paper's "clip it so that the variable
            // is still within its global bounds".  Both ends of the clip window
            // bracket zero, so this can only shorten a shift, never reverse
            // one.  "No domain enlargement" is enforced a line earlier, by the
            // interval keeping its width; the clip is a separate rule.
            s = cpp_max(s, col_lb[j as usize] - cur_lb);
            s = cpp_min(s, col_ub[j as usize] - cur_ub);
            if integer {
                // A fractional structural bound could have made the clip
                // fractional; round back toward zero, which can only shrink
                // the shift and so cannot re-break the clip.
                s = if s > 0.0 { s.floor() } else { s.ceil() };
            }
            if s.abs() < 1e-15 {
                k += 1;
                continue;
            }

            // Skip the variables that do not reduce this row's violation
            // (paper Sect. 5).
            let shifted = a * s;
            let new_viol = range_violation(m_r + shifted, mx_r + shifted, lo_r, hi_r);
            if new_viol >= cur_viol - feastol {
                k += 1;
                continue;
            }

            // Damage: the *increases* in violation over the other rows `j`
            // appears in.  As in the original WalkSAT, improvements are
            // ignored.  Measured as summed magnitude, not as a count of newly
            // violated rows (see the C++ comment: the paper's partial-assignment
            // reading is a magnitude, and this matches `walksat_select_move`).
            let cbeg = csc_start[j as usize];
            let cend = csc_start[(j + 1) as usize];
            effort += (cend - cbeg) as usize;
            let mut damage = 0.0f64;
            let mut p = cbeg;
            while p < cend {
                let i2 = csc_row[p as usize];
                if i2 == row {
                    p += 1;
                    continue;
                }
                let d = csc_val[p as usize] * s;
                let dv = range_violation(
                    min_act[i2 as usize] + d,
                    max_act[i2 as usize] + d,
                    row_lo[i2 as usize],
                    row_hi[i2 as usize],
                ) - range_violation(
                    min_act[i2 as usize],
                    max_act[i2 as usize],
                    row_lo[i2 as usize],
                    row_hi[i2 as usize],
                );
                if dv > 0.0 {
                    damage += dv;
                }
                p += 1;
            }
            best_damage = cpp_min(damage, best_damage);
            ls_evals += 1;
            cand.push(RepairWalkCandidate {
                var: j,
                shift: s,
                damage,
            });
            k += 1;
        }

        if cand.is_empty() {
            step += 1;
            continue;
        }

        // WalkSAT selection (paper Sect. 5): a zero-damage shift is taken
        // greedily; otherwise a random candidate with probability `noise` and
        // a minimum-damage one with probability 1 - `noise`.
        let pick: HighsInt = if best_damage > feastol && s5_uniform_real(rng, 0.0, 1.0) < noise {
            s5_uniform_int_i32(rng, 0, cand.len() as HighsInt - 1)
        } else {
            best_indices.clear();
            let threshold = best_damage + feastol;
            let mut ci: usize = 0;
            while ci < cand.len() {
                if cand[ci].damage <= threshold {
                    best_indices.push(ci as HighsInt);
                }
                ci += 1;
            }
            best_indices[s5_uniform_int_i32(rng, 0, best_indices.len() as HighsInt - 1) as usize]
        };
        let var = cand[pick as usize].var;
        let shift = cand[pick as usize].shift;

        // Apply, then re-read the engine's activities rather than trusting the
        // predicted translation: `shift_domain` may clamp, and a fixed-but-wide
        // column collapses to a singleton, so the applied move is not always
        // the predicted one.
        let cbeg = csc_start[var as usize];
        let cend = csc_start[(var + 1) as usize];
        effort += 2 * (cend - cbeg) as usize;
        let mut removed = 0.0f64;
        {
            let mut p = cbeg;
            while p < cend {
                let i2 = csc_row[p as usize];
                removed += range_violation(
                    min_act[i2 as usize],
                    max_act[i2 as usize],
                    row_lo[i2 as usize],
                    row_hi[i2 as usize],
                );
                p += 1;
            }
        }
        if !e.shift_domain(var, shift) {
            // Unreachable given the clip above; a guard, not a path.
            step += 1;
            continue;
        }
        ls_commits += 1;
        // Re-read the engine's activity slices after the shift (the borrow of
        // `e` above ended at `shift_domain`).
        let min_act = e.min_activity_data();
        let max_act = e.max_activity_data();
        let mut added = 0.0f64;
        {
            let mut p = cbeg;
            while p < cend {
                let i2 = csc_row[p as usize];
                added += range_violation(
                    min_act[i2 as usize],
                    max_act[i2 as usize],
                    row_lo[i2 as usize],
                    row_hi[i2 as usize],
                );
                // refresh_violated(i2)
                let was = violated_pos[i2 as usize] != -1;
                let now = range_violated(
                    min_act[i2 as usize],
                    max_act[i2 as usize],
                    row_lo[i2 as usize],
                    row_hi[i2 as usize],
                    feastol,
                );
                if was && !now {
                    // remove_violated(i2)
                    let pos = violated_pos[i2 as usize];
                    if pos != -1 {
                        let last = *violated.last().unwrap();
                        violated[pos as usize] = last;
                        violated_pos[last as usize] = pos;
                        violated.pop();
                        violated_pos[i2 as usize] = -1;
                    }
                } else if !was && now {
                    // add_violated(i2)
                    if violated_pos[i2 as usize] == -1 {
                        violated_pos[i2 as usize] = violated.len() as HighsInt;
                        violated.push(i2);
                    }
                }
                p += 1;
            }
        }
        total_viol += added - removed;

        tabu.push(var);
        if tabu.len() > kTabuLength as usize {
            tabu.remove(0);
        }

        if total_viol < best_viol - feastol {
            best_viol = total_viol;
            best_mark = mark_now(e);
        }

        since_restart += 1;
        if since_restart >= kSoftRestartPeriod {
            since_restart = 0;
            if total_viol > best_viol + feastol {
                // restore_best(): backtrack_to reverts exactly the shifts made
                // since best_mark, and the rescan then reads the restored
                // activities, so the set and the total are correct by
                // construction.
                e.backtrack_to(best_mark.vs, best_mark.sol, best_mark.act, best_mark.pq);
                total_viol = rescan(
                    violated,
                    violated_pos,
                    e.min_activity_data(),
                    e.max_activity_data(),
                    row_lo,
                    row_hi,
                    feastol,
                    nrow,
                    &mut effort,
                );
            }
        }
        // The tabu list deliberately survives a soft restart.  The paper
        // introduces the two tricks side by side and says nothing about their
        // interaction, so this is our reading, not its text: the list is a
        // rolling window over the last three shifts, full stop.

        step += 1;
    }

    // Leave the node in the least-violated state the walk saw.  **Not a paper
    // citation**: Sect. 5 names exactly two additions to the basic walk, the
    // tabu list and the soft restart, and no end-of-walk restore.  This is
    // ours, and it is the same property `FprConfig::repair_track_best` gives
    // the leaf-time walk.  Unlike at the leaf it is not optional: the soft
    // restart above needs the same best-state bookkeeping, and that one *is*
    // the paper's.
    if !violated.is_empty() && total_viol > best_viol + feastol {
        // restore_best()
        e.backtrack_to(best_mark.vs, best_mark.sol, best_mark.act, best_mark.pq);
        rescan(
            violated,
            violated_pos,
            e.min_activity_data(),
            e.max_activity_data(),
            row_lo,
            row_hi,
            feastol,
            nrow,
            &mut effort,
        );
    }

    super::stats::FPR.add_moves(ls_evals, ls_commits);
    *effort_out = effort;
    violated.is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::HighsVarType;
    use crate::mip::fpr_local_mip::PropEngine;
    use crate::mip::fpr_local_mip::RepairWalkScratch;
    use crate::mip::fpr_local_mip::rng::S5Rng;
    use crate::mip::fpr_local_mip::rng::build_csc;

    // x0 + x1 >= 1, 3*x0 + x1 <= 2, both binary.
    //
    // Fixing x0 = 1 makes the second row unsatisfiable by any completion, and
    // the repair that recovers it is available on exactly one column: x0's
    // singleton domain {1} can slide down to {0} inside its structural [0, 1],
    // while x1's *unfixed* [0, 1] cannot slide at all.  So the walk is
    // deterministic here: one candidate, zero damage, chosen greedily.
    struct SwapModel {
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

    impl SwapModel {
        fn new() -> Self {
            let ar_start = vec![0, 2, 4];
            let ar_index = vec![0, 1, 0, 1];
            let ar_value = vec![1.0, 1.0, 3.0, 1.0];
            let csc = build_csc(2, 2, &ar_start, &ar_index, &ar_value);
            SwapModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![0.0, 0.0],
                col_ub: vec![1.0, 1.0],
                row_lo: vec![1.0, -kHighsInf],
                row_hi: vec![kHighsInf, 2.0],
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

    // One general integer column and one row `coeff * y {<=,>=} rhs`, so the
    // only repair available is a shift of `y`'s own interval and a test can
    // assert the *whole* resulting domain rather than a value.
    struct OneRowModel {
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

    impl OneRowModel {
        fn new(coeff: f64, lo: f64, hi: f64, structural_lb: f64) -> Self {
            let ar_start = vec![0, 1];
            let ar_index = vec![0];
            let ar_value = vec![coeff];
            let csc = build_csc(1, 1, &ar_start, &ar_index, &ar_value);
            OneRowModel {
                ar_start,
                ar_index,
                ar_value,
                col_lb: vec![structural_lb],
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

    // Narrow column 0 to `[lo, hi]` the way propagation would have.
    fn narrow(engine: &mut PropEngine, lo: f64, hi: f64) {
        assert!(engine.tighten_lb(0, lo));
        assert!(engine.tighten_ub(0, hi));
    }

    fn walk(
        engine: &mut PropEngine,
        max_steps: HighsInt,
        rng: &mut S5Rng,
        effort: &mut usize,
        scratch: &mut RepairWalkScratch,
        noise: f64,
    ) -> bool {
        repair_walk(
            engine,
            max_steps,
            noise,
            rng,
            effort,
            scratch,
            &Deadline::default(),
        )
    }

    #[test]
    fn s5_walk_repair_walk_partial_assignment() {
        // Port of test_repair_walk.cpp "repairs a partial assignment the DFS
        // propagated to infeasibility" (:123).
        let m = SwapModel::new();
        let mut e = m.make_engine(1e-6);
        e.init_activities();

        // The node: x0 decided to 1, x1 still *unfixed*.  This is the state
        // Fig. 1 line 7 hands to repair, and it is a partial assignment.
        assert!(e.fix(0, 1.0));
        assert_eq!(
            e.propagate(0),
            crate::mip::fpr_local_mip::PropResult::kInfeasible
        );
        assert!(!e.var(1).fixed);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(11);
        let mut effort: usize = 0;
        assert!(walk(&mut e, 200, &mut rng, &mut effort, &mut scratch, 0.75));
        assert!(effort > 0);

        // x0's singleton domain slid to {0}: the decision was overridden by
        // the repair, which is the point of repairing rather than pruning.
        assert!(e.var(0).fixed);
        assert!((e.var(0).val - 0.0).abs() < 1e-6);
        // x1 was never a legal candidate (its [0, 1] cannot slide down inside
        // [0, 1]) and must therefore be untouched -- not fixed to something
        // convenient, and not widened.
        assert!(!e.var(1).fixed);
        assert!((e.var(1).lb - 0.0).abs() < 1e-6);
        assert!((e.var(1).ub - 1.0).abs() < 1e-6);
    }

    #[test]
    fn s5_walk_repair_walk_shifts_domain_not_global_bounds() {
        // Port of test_repair_walk.cpp "shifts the node's current domain,
        // never the global bounds" (:156).
        // y in [0, 10] structurally, narrowed to [5, 7] at this node.  The row
        // `y <= 2` wants the activity down to 2, which is a shift of -3.
        let m = OneRowModel::new(1.0, -kHighsInf, 2.0, 0.0);
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        narrow(&mut e, 5.0, 7.0);
        assert!(!e.var(0).fixed);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(3);
        let mut effort: usize = 0;
        assert!(walk(&mut e, 200, &mut rng, &mut effort, &mut scratch, 0.75));

        // The paper's shift: the interval *translates*, keeping its width.  A
        // repair that operated on the global bounds instead -- picking a value
        // anywhere in [0, 10], or reopening the domain to it -- would not land
        // on [2, 4].
        assert!(!e.var(0).fixed);
        assert!((e.var(0).lb - 2.0).abs() < 1e-6);
        assert!((e.var(0).ub - 4.0).abs() < 1e-6);
        assert!(((e.var(0).ub - e.var(0).lb) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn s5_walk_repair_walk_shift_clipped_by_structural_bounds() {
        // Port of test_repair_walk.cpp "a shift is clipped by the structural
        // bounds, not by the row's demand" (:184).
        // Same node domain [5, 7] and same demand for a -3 shift, but the
        // structural lower bound is now 3, so only -2 of it is available.  The
        // walk must take the clipped shift and report the node still
        // infeasible, rather than either overshooting the structural bound or
        // refusing the partial improvement.
        let m = OneRowModel::new(1.0, -kHighsInf, 2.0, 3.0);
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        narrow(&mut e, 5.0, 7.0);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(5);
        let mut effort: usize = 0;
        assert!(!walk(
            &mut e,
            200,
            &mut rng,
            &mut effort,
            &mut scratch,
            0.75
        ));

        assert!((e.var(0).lb - 3.0).abs() < 1e-6);
        assert!((e.var(0).ub - 5.0).abs() < 1e-6);
    }

    #[test]
    fn s5_walk_repair_walk_undershoot_rounds_other_way() {
        // Port of test_repair_walk.cpp "an undershooting row rounds the shift
        // the other way" (:207).
        // `3y >= 8` with y narrowed to [1, 2]: the activity range [3, 6] sits
        // *below* the row, so the minimal shift is 2/3 of a unit and must
        // round **up** to 1.  Rounding it the way an overshooting row rounds
        // (down, to 0) yields no move at all and the node stays infeasible.
        let m = OneRowModel::new(3.0, 8.0, kHighsInf, 0.0);
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        narrow(&mut e, 1.0, 2.0);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(13);
        let mut effort: usize = 0;
        assert!(walk(&mut e, 200, &mut rng, &mut effort, &mut scratch, 0.75));

        assert!((e.var(0).lb - 2.0).abs() < 1e-6);
        assert!((e.var(0).ub - 3.0).abs() < 1e-6);
    }

    #[test]
    fn s5_walk_repair_walk_negative_coefficient_flips_direction() {
        // Port of test_repair_walk.cpp "a negative coefficient flips which way
        // the activity moves" (:230).
        // `-2y <= -9` with y narrowed to [1, 3]: the activity range is
        // [-6, -2], which overshoots the row's upper bound of -9 from above,
        // and the shift that fixes it is *positive* because the coefficient is
        // negative.  1.5 units of it, rounding up to 2.
        let m = OneRowModel::new(-2.0, -kHighsInf, -9.0, 0.0);
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        narrow(&mut e, 1.0, 3.0);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(17);
        let mut effort: usize = 0;
        assert!(walk(&mut e, 200, &mut rng, &mut effort, &mut scratch, 0.75));

        assert!((e.var(0).lb - 3.0).abs() < 1e-6);
        assert!((e.var(0).ub - 5.0).abs() < 1e-6);
    }

    #[test]
    fn s5_walk_repair_walk_zero_steps_is_noop() {
        // Port of test_repair_walk.cpp "a zero step limit leaves the node
        // exactly as it found it" (:251).
        let m = SwapModel::new();
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        assert!(e.fix(0, 1.0));
        assert_eq!(
            e.propagate(0),
            crate::mip::fpr_local_mip::PropResult::kInfeasible
        );

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(7);
        let mut effort: usize = 0;
        assert!(!walk(&mut e, 0, &mut rng, &mut effort, &mut scratch, 0.75));
        assert!(e.var(0).fixed);
        assert!((e.var(0).val - 1.0).abs() < 1e-6);
    }

    // A general model: `ncol` integer columns, rows given row-wise.
    struct GenModel {
        ncol: HighsInt,
        nrow: HighsInt,
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

    impl GenModel {
        #[allow(clippy::too_many_arguments)]
        fn new(
            ncol: HighsInt,
            ar_start: Vec<HighsInt>,
            ar_index: Vec<HighsInt>,
            ar_value: Vec<f64>,
            col_lb: Vec<f64>,
            col_ub: Vec<f64>,
            row_lo: Vec<f64>,
            row_hi: Vec<f64>,
        ) -> Self {
            let nrow = (ar_start.len() - 1) as HighsInt;
            let csc = build_csc(ncol, nrow, &ar_start, &ar_index, &ar_value);
            GenModel {
                ncol,
                nrow,
                ar_start,
                ar_index,
                ar_value,
                col_lb,
                col_ub,
                row_lo,
                row_hi,
                integrality: vec![HighsVarType::kInteger; ncol as usize],
                csc_start: csc.col_start,
                csc_row: csc.col_row,
                csc_val: csc.col_val,
            }
        }

        fn make_engine(&self, feastol: f64) -> PropEngine<'_> {
            PropEngine::new(
                self.ncol,
                self.nrow,
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
    fn s5_walk_repair_walk_damage_picks_move() {
        // Port of test_repair_walk.cpp "damage picks the move, not arrival
        // order" (:301).  DamageModel: p + q <= 2, q >= 1, p >= 4, both
        // general integers on [0, 10], narrowed to [3, 4].
        for seed in [23u64, 5, 7, 11, 13, 17, 19, 29] {
            let m = GenModel::new(
                2,
                vec![0, 2, 3, 4],
                vec![0, 1, 1, 0],
                vec![1.0, 1.0, 1.0, 1.0],
                vec![0.0, 0.0],
                vec![10.0, 10.0],
                vec![-kHighsInf, 1.0, 4.0],
                vec![2.0, kHighsInf, kHighsInf],
            );
            let mut e = m.make_engine(1e-6);
            e.init_activities();
            for j in 0..2 {
                assert!(e.tighten_lb(j, 3.0));
                assert!(e.tighten_ub(j, 4.0));
            }

            let mut scratch = RepairWalkScratch::default();
            let mut rng = S5Rng::new(seed);
            let mut effort: usize = 0;
            // Exactly one step, so the assertion is about *this* move.
            let _ = walk(&mut e, 1, &mut rng, &mut effort, &mut scratch, 0.75);

            assert!((e.var(1).lb - 0.0).abs() < 1e-6, "seed {seed}");
            assert!((e.var(1).ub - 1.0).abs() < 1e-6, "seed {seed}");
            assert!((e.var(0).lb - 3.0).abs() < 1e-6, "seed {seed}");
            assert!((e.var(0).ub - 4.0).abs() < 1e-6, "seed {seed}");
        }
    }

    #[test]
    fn s5_walk_repair_walk_tabu_stops_two_move_cycle() {
        // Port of test_repair_walk.cpp "the tabu list stops a two-move cycle
        // after one move" (:366).  OscillateModel: `y <= 2` and `y >= 5` over
        // one column narrowed to [5, 7].
        let m = GenModel::new(
            1,
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 1.0],
            vec![0.0],
            vec![10.0],
            vec![-kHighsInf, 5.0],
            vec![2.0, kHighsInf],
        );
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        narrow(&mut e, 5.0, 7.0);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(29);
        let mut effort: usize = 0;
        assert!(!walk(
            &mut e,
            200,
            &mut rng,
            &mut effort,
            &mut scratch,
            0.75
        ));

        assert!((e.var(0).lb - 2.0).abs() < 1e-6);
        assert!((e.var(0).ub - 4.0).abs() < 1e-6);
        assert!(effort < 500);
    }

    #[test]
    fn s5_walk_repair_walk_failed_walk_leaves_best_state() {
        // Port of test_repair_walk.cpp "a failed walk leaves the node in the
        // best state it saw" (:419).  WorseningModel: `y <= 2` and
        // `10y >= 50` over one column narrowed to [5, 7].
        let m = GenModel::new(
            1,
            vec![0, 1, 2],
            vec![0, 0],
            vec![1.0, 10.0],
            vec![0.0],
            vec![10.0],
            vec![-kHighsInf, 50.0],
            vec![2.0, kHighsInf],
        );
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        narrow(&mut e, 5.0, 7.0);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(31);
        let mut effort: usize = 0;
        assert!(!walk(
            &mut e,
            200,
            &mut rng,
            &mut effort,
            &mut scratch,
            0.75
        ));

        assert!((e.var(0).lb - 5.0).abs() < 1e-6);
        assert!((e.var(0).ub - 7.0).abs() < 1e-6);
    }

    #[test]
    fn s5_walk_repair_walk_backtrack_restores_every_array() {
        // Port of test_repair_walk.cpp "a backtrack past a repair restores
        // every engine array" (:442).
        let m = SwapModel::new();
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        e.init_domain_pq();
        assert!(e.fix(0, 1.0));
        assert_eq!(
            e.propagate(0),
            crate::mip::fpr_local_mip::PropResult::kInfeasible
        );

        let vs_m = e.vs_mark();
        let sol_m = e.sol_mark();
        let act_m = e.act_mark();
        let pq_m = e.pq_mark();
        let min_act_before = e.row_min_activity(1);
        assert!((min_act_before - 3.0).abs() < 1e-6);
        assert!((e.sol(0) - 1.0).abs() < 1e-6);
        assert_eq!(e.pq_top(), 1);

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(37);
        let mut effort: usize = 0;
        assert!(walk(&mut e, 200, &mut rng, &mut effort, &mut scratch, 0.75));
        assert!((e.var(0).val - 0.0).abs() < 1e-6);
        assert!((e.sol(0) - 0.0).abs() < 1e-6);
        assert!((e.row_min_activity(1) - 0.0).abs() < 1e-6);

        e.backtrack_to(vs_m, sol_m, act_m, pq_m);

        assert!(e.var(0).fixed);
        assert!((e.var(0).val - 1.0).abs() < 1e-6);
        assert!((e.var(0).lb - 0.0).abs() < 1e-6);
        assert!((e.var(0).ub - 1.0).abs() < 1e-6);
        assert!((e.sol(0) - 1.0).abs() < 1e-6);
        assert!((e.row_min_activity(1) - min_act_before).abs() < 1e-6);
        assert!(!e.var(1).fixed);
        assert_eq!(e.pq_top(), 1);
    }

    #[test]
    fn s5_walk_repair_walk_long_drift_restored() {
        // Port of test_repair_walk.cpp "a long drifting walk is restored to
        // the state it started in" (:1107).  DriftModel: twelve binaries
        // fixed to 1, each `x_j <= 0`, plus `sum 3 x_j >= 36`.
        const NCOL: HighsInt = 12;
        let mut ar_start = vec![0];
        let mut ar_index = Vec::new();
        let mut ar_value = Vec::new();
        let mut row_lo = Vec::new();
        let mut row_hi = Vec::new();
        for j in 0..NCOL {
            ar_index.push(j);
            ar_value.push(1.0);
            ar_start.push(ar_index.len() as HighsInt);
            row_lo.push(-kHighsInf);
            row_hi.push(0.0);
        }
        for j in 0..NCOL {
            ar_index.push(j);
            ar_value.push(3.0);
        }
        ar_start.push(ar_index.len() as HighsInt);
        row_lo.push(3.0 * NCOL as f64);
        row_hi.push(kHighsInf);
        let m = GenModel::new(
            NCOL,
            ar_start,
            ar_index,
            ar_value,
            vec![0.0; NCOL as usize],
            vec![1.0; NCOL as usize],
            row_lo,
            row_hi,
        );
        let mut e = m.make_engine(1e-6);
        e.init_activities();
        for j in 0..NCOL {
            assert!(e.fix(j, 1.0));
        }

        let mut scratch = RepairWalkScratch::default();
        let mut rng = S5Rng::new(53);
        let mut effort: usize = 0;
        assert!(!walk(
            &mut e,
            200,
            &mut rng,
            &mut effort,
            &mut scratch,
            0.75
        ));

        for j in 0..NCOL {
            assert!(e.var(j).fixed, "column {j}");
            assert!((e.var(j).val - 1.0).abs() < 1e-6, "column {j}");
        }
    }
}
