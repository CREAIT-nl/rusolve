//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP 11.0.0 (32dd425): conflict_resolution.c::existsResolvablebdchginfo (1041-1060),
//!   isBdchgConflictRelevant (1064-1091), conflictFirstCand (1095-1137),
//!   conflictRemoveCand (1140-1157), updateBdchgQueue (1678-1726);
//!   conflict_graphanalysis.c::betterBoundInResolutionQueue (2341-2380);
//!   conflict_general.c::conflictBdchginfoComp (233-252).

use super::{S2BdchgQueue, S2Change, S2ConflictRow, S2Context, s2_is_gt, s2_is_lt};
use crate::prelude::*;

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::isBdchgConflictRelevant (lines 1064-1091)
pub fn is_bdchg_conflict_relevant(row: &S2ConflictRow, chg: &S2Change, initial: bool) -> bool {
    // the initial bound change is always relevant
    if initial {
        return true;
    }
    // if the conflict row is empty, we have a global infeasibility
    if row.nnz() == 0 {
        return false;
    }

    let idx = chg.col as usize;

    if row.vals[idx] > 0.0 && chg.is_upper {
        return true;
    }
    if row.vals[idx] < 0.0 && !chg.is_upper {
        return true;
    }
    false
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::updateBdchgQueue (lines 1678-1726)
pub fn update_bdchg_queue<C: S2Context>(
    ctx: &C,
    q: &mut S2BdchgQueue,
    row: &S2ConflictRow,
    upto: HighsInt,
) {
    for i in 0..row.nnz() {
        let v = row.inds[i];
        let coef = row.vals[v as usize];

        if coef > 0.0 {
            let (b, p) = ctx.ub_upto(v, upto);
            if p != -1 && s2_is_lt(b, ctx.global_ub(v)) {
                q.add_bound(ctx, p);
            }
        } else {
            let (b, p) = ctx.lb_upto(v, upto);
            if p != -1 && s2_is_gt(b, ctx.global_lb(v)) {
                q.add_bound(ctx, p);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// impl S2BdchgQueue
// ---------------------------------------------------------------------------

impl S2BdchgQueue {
    /// Reset the queue for a new conflict analysis with `ncol` columns.
    pub fn reset(&mut self, ncol: usize) {
        self.heap.clear();
        self.conflictvarslbs = vec![f64::MIN; ncol];
        self.conflictvarsubs = vec![f64::MAX; ncol];
    }

    /// SCIP 11.0.0 (32dd425): conflict_graphanalysis.c::betterBoundInResolutionQueue (lines 2341-2380)
    /// Insert a bound change into the resolution queue only if it is strictly tighter
    /// than what is already recorded for that column.
    pub fn add_bound<C: S2Context>(&mut self, ctx: &C, pos: HighsInt) {
        let chg = ctx.change(pos);
        let col = chg.col as usize;

        if !chg.is_upper {
            // lower bound
            if self.conflictvarslbs[col] < chg.newbound {
                self.conflictvarslbs[col] = chg.newbound;
                self.heap.push(pos);
            }
        } else {
            // upper bound
            if self.conflictvarsubs[col] > chg.newbound {
                self.conflictvarsubs[col] = chg.newbound;
                self.heap.push(pos);
            }
        }
    }

    /// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictFirstCand (lines 1095-1137)
    /// Iterative form of the recursion: peek the top of the heap; while the candidate
    /// is not relevant, pop it and reset the corresponding bound array entry.
    pub fn first_cand<C: S2Context>(
        &mut self,
        ctx: &C,
        row: &S2ConflictRow,
        initial: bool,
    ) -> Option<HighsInt> {
        loop {
            let top = match self.heap.peek() {
                Some(&p) => p,
                None => return None,
            };
            let chg = ctx.change(top);
            if is_bdchg_conflict_relevant(row, &chg, initial) {
                return Some(top);
            }
            // not relevant: pop and reset
            self.heap.pop();
            let col = chg.col as usize;
            if !chg.is_upper {
                self.conflictvarslbs[col] = f64::MIN;
            } else {
                self.conflictvarsubs[col] = f64::MAX;
            }
        }
    }

    /// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictRemoveCand (lines 1140-1157)
    pub fn remove_cand(&mut self) -> Option<HighsInt> {
        self.heap.pop()
    }

    /// SCIP 11.0.0 (32dd425): conflict_resolution.c::existsResolvablebdchginfo (lines 1041-1060)
    pub fn exists_resolvable<C: S2Context>(&self, ctx: &C) -> bool {
        self.heap.iter().any(|&pos| ctx.change(pos).resolvable)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::{S2TestChg, S2TestCtx};
    use super::*;

    #[test]
    fn s2_queue_dedup() {
        // 1 column x0, glb [0], gub [1], not integral
        let ctx = S2TestCtx {
            glb: vec![0.0],
            gub: vec![1.0],
            integral: vec![false],
            num_integral: 0,
            stack: vec![
                // pos 0: dummy lb change on x0 (unused)
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 1: x0 ub 1 -> 0.5
                S2TestChg {
                    col: 0,
                    is_upper: true,
                    newbound: 0.5,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
                // pos 2: dummy lb change on x0 (unused)
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 3: x0 ub 0.5 -> 0.0
                S2TestChg {
                    col: 0,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 2,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let mut q = S2BdchgQueue::default();
        q.reset(1);

        // add pos 3 first (ub 0.0), then pos 1 (ub 0.5)
        q.add_bound(&ctx, 3);
        q.add_bound(&ctx, 1);

        // heap should hold only pos 3 (the tighter bound)
        assert_eq!(q.heap.len(), 1);
        assert_eq!(q.heap.peek(), Some(&3));
        assert_eq!(q.conflictvarsubs[0], 0.0);
    }

    #[test]
    fn s2_queue_order() {
        // 3 columns, changes at positions 1, 4, 2 on different cols
        // Fill other positions with dummy changes on col 0 lb with newbound equal to glb
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![1.0, 1.0, 1.0],
            integral: vec![false, false, false],
            num_integral: 0,
            stack: vec![
                // pos 0: dummy on col 0
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 1: col 0 ub change
                S2TestChg {
                    col: 0,
                    is_upper: true,
                    newbound: 0.5,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
                // pos 2: col 1 ub change
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.5,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
                // pos 3: dummy on col 0
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 4: col 2 ub change
                S2TestChg {
                    col: 2,
                    is_upper: true,
                    newbound: 0.5,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let mut q = S2BdchgQueue::default();
        q.reset(3);

        q.add_bound(&ctx, 1);
        q.add_bound(&ctx, 4);
        q.add_bound(&ctx, 2);

        // BinaryHeap is a max-heap, so pops give 4, 2, 1
        assert_eq!(q.remove_cand(), Some(4));
        assert_eq!(q.remove_cand(), Some(2));
        assert_eq!(q.remove_cand(), Some(1));
        assert_eq!(q.remove_cand(), None);
    }

    #[test]
    fn s2_queue_first_cand_irrelevant() {
        // 2 columns x0, x1; row: nnz 2, inds [0,1], vals [1.0, 1.0], lhs 0
        // pos 2 = x1 ub change (newbound 0.5, resolvable via reason Some)
        // pos 5 = x0 lb change newbound 0.5 (not relevant because vals[0]>0 but is_upper=false)
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![1.0, 1.0],
            integral: vec![false, false],
            num_integral: 0,
            stack: vec![
                // pos 0: dummy
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 1: dummy
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 2: x1 ub change, resolvable
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.5,
                    depth: 1,
                    branching: false,
                    reason: Some((vec![], 0.0)),
                },
                // pos 3: dummy
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 4: dummy
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                // pos 5: x0 lb change newbound 0.5
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 0.5,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let row = S2ConflictRow {
            vals: vec![1.0, 1.0],
            inds: vec![0, 1],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 1.0,
        };

        let mut q = S2BdchgQueue::default();
        q.reset(2);

        // add pos 5 (x0 lb 0.5) and pos 2 (x1 ub 0.5)
        // For pos 5: is_upper=false, conflictvarslbs[0] starts at f64::MIN < 0.5, so add succeeds
        // For pos 2: is_upper=true, conflictvarsubs[1] starts at f64::MAX > 0.5, so add succeeds
        q.add_bound(&ctx, 5);
        q.add_bound(&ctx, 2);

        // first_cand: top is 5 (max-heap). chg at pos 5: col 0, is_upper=false.
        // is_bdchg_conflict_relevant: vals[0]=1.0 > 0, but is_upper=false -> not relevant.
        // Pop 5, reset conflictvarslbs[0] = f64::MIN.
        // Next top is 2. chg at pos 2: col 1, is_upper=true.
        // is_bdchg_conflict_relevant: vals[1]=1.0 > 0, is_upper=true -> relevant. Return Some(2).
        let result = q.first_cand(&ctx, &row, false);
        assert_eq!(result, Some(2));
        assert_eq!(q.conflictvarslbs[0], f64::MIN);
        assert_eq!(q.heap.peek(), Some(&2));
    }

    #[test]
    fn s2_queue_update() {
        // 2 columns, row x0 - x1 >= 0: inds [0,1], vals [1.0, -1.0]
        // binary both: integral true, glb [0,0], gub [1,1]
        // Stack: pos 0 x1 lb 0->1, pos 1 x0 ub 1->0
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![1.0, 1.0],
            integral: vec![true, true],
            num_integral: 2,
            stack: vec![
                // pos 0: x1 lb 0 -> 1
                S2TestChg {
                    col: 1,
                    is_upper: false,
                    newbound: 1.0,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
                // pos 1: x0 ub 1 -> 0
                S2TestChg {
                    col: 0,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 1,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let row = S2ConflictRow {
            vals: vec![1.0, -1.0],
            inds: vec![0, 1],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 1.0,
        };

        // update_bdchg_queue with upto=1
        let mut q = S2BdchgQueue::default();
        q.reset(2);
        update_bdchg_queue(&ctx, &mut q, &row, 1);
        // heap should contain {0, 1}
        assert_eq!(q.heap.len(), 2);

        // fresh q with upto=0
        let mut q2 = S2BdchgQueue::default();
        q2.reset(2);
        update_bdchg_queue(&ctx, &mut q2, &row, 0);
        // heap should contain {0} only
        assert_eq!(q2.heap.len(), 1);
        assert_eq!(q2.heap.peek(), Some(&0));
    }
}
