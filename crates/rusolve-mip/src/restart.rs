//! Root restart.
//!
//! When the root cut loop and reduced-cost fixing have fixed a significant
//! share of the integer columns, the search stops, MIP presolve runs again
//! on the tightened problem (with the active cut rows appended to the
//! model), and the root is redone on the smaller problem. The incumbent and
//! the cumulative counters survive the restart; the postsolve stacks chain.
//!
//! Provenance: ported from HiGHS (MIT, read as a reference text):
//! `HighsMipSolverData::evaluateRootNode` restart checks
//! (reference/HiGHS/highs/mip/HighsMipSolverData.cpp:2127–2136 pre-separation
//! — not ported — and :2486–2505 post-separation, the trigger thresholds),
//! `performRestart` (:1317: append active cuts to the model, re-presolve the
//! reduced problem, keep the incumbent, re-reference the objective limits
//! across the offset change) and `percentageInactiveIntegers` (:1310).
//! The handshake split between the search (`bnb`) and the presolve chain
//! (`presolve::solve_presolved`) is ours, forced by running presolve outside
//! the search.
//!
//! No tree-size-estimate restarts yet (HighsMipSolver.cpp `checkRestart`
//! :475–540; Anderson, Hendel, Le Bodic & Viernickel 2019).

use crate::bnb::FirstIncumbent;
use crate::cuts;
use crate::engine::CutRow;
use crate::problem::MipProblem;
use rusolve_core::{CscMatrix, Triplet};

/// Hard cap on restarts per solve, guarding against a presolve that keeps
/// finding the same fixings (HiGHS prevents the cycle with
/// `restart_presolve_reduction_limit >= 1`, HighsMipSolverData.cpp:1395–1402;
/// a count cap is our cruder equivalent).
pub const MAX_RESTARTS: u32 = 10;

/// Our design, no HiGHS counterpart: a restart is only worth its cost when
/// the root LP bound is still far from the incumbent, i.e. the search has
/// real work left for a tightened presolve to reshape. HiGHS can restart at
/// any gap because its root cut loop usually closes the root and the restart
/// rides the presolve it needed anyway; our weaker root closure lets the
/// faithful trigger fire where branching would finish immediately. Measured
/// on the round-1 A/B: the small families fire the rate trigger and inflate
/// (set_packing 50.8% fixed at root gap 0.013, set_cover 59.6%/0.035,
/// budget_max_cover 26.4%/0.020, fcn 9.3%/0.047), while the instance the
/// restart exists for — mikety_mandi_1, 41.1 s → 29.1 s and 198,666 → 16,466
/// nodes — checks in at gap 0.163. 0.10 separates the two populations with
/// margin on both sides. A missing incumbent is an infinite gap: the gate
/// never blocks a cold start, where a tightened re-presolve helps most.
pub const MIN_ROOT_GAP: f64 = 0.10;

/// Stall restart: minimum tree size before an improvement-stalled restart
/// may fire (nodes processed since the stall window opened). The stall restart
/// covers the instances HiGHS restarts through clique-table substitutions we
/// do not have — they reach the root check with 0 % fixings, a wide gap and
/// cut mass (pg: 0 % / gap 0.31 / 2 473 tight cuts), a shape no root-only
/// signal separates from int_production s1 (0 % / 0.61 / 170) or jobshop s1
/// (0 % / 0.61 / 106): a census of every root input on our benchmark set
/// measured them identical. The separating evidence is
/// the tree itself, as in HiGHS's tree-size-estimate restart
/// (HighsMipSolver.cpp `checkRestart` :475–540; Anderson, Hendel, Le Bodic &
/// Viernickel 2019, *Tree-size estimation*): a tree that has processed this
/// many nodes without the incumbent moving is demonstrably not clean, so the
/// re-presolve + root redo amortizes. Measured on pg (the stall restart's
/// target shape; window opens at its last significant improvement, node 237;
/// baseline 24.4 s / 2 863 nodes): bar 1 000 fires at node 1 261 and only
/// ties on time (26.9 s, nodes 0.76×), 500 fires at 795 and ties (24.9 s,
/// 0.49×), 300 fires at 561 and wins (18.8 s, 0.77× time / 0.41× nodes) —
/// the second tree is ~600 nodes at every bar, so the earlier the fire the
/// less tree it discards against the fixed ~5 s re-presolve of the cut-heavy
/// restart problem. Every census family whose whole tree fits under the bar
/// (jobshop 90, facility_mixed 92) stays untouched, and the ones over it are
/// held back by [`STALL_SECS`] (market_split: 10 344 nodes in 0.06 s,
/// block_angular 663 in 0.15 s).
pub const STALL_NODES: u64 = 300;

/// Stall restart: the stall window must also span this many seconds since the
/// window opened, so a tree that runs hundreds of nodes in milliseconds
/// (market_split s1: 10 344 nodes in 0.06 s) never pays a restart whose
/// re-presolve costs more than the whole solve. Wall-time amortization in the
/// style of the sub-MIP cost-aware trigger (`submip::SUBMIP_TIME_SHARE`,
/// ours; HiGHS's analog effort accounting is `moreHeuristicsAllowed`,
/// HighsMipSolverData.cpp:672–731, used to break the root loop at
/// :2392–2400).
pub const STALL_SECS: f64 = 2.0;

/// Stall restart: tight cut rows per model row the root must leave for the
/// re-presolve to have something to strip (arm-time proxy for HiGHS's
/// `restart_presolve_reduction_limit ≥ 1`, HighsMipSolverData.cpp:1390–1402:
/// a restart whose re-presolve reduces nothing cycles and never pays).
/// Measured 2026-10-24: with any-cut arming, enlight_hard (106 tight cuts on
/// 201 rows = 0.5/row, nothing for presolve to deduce) fired at node 19 401
/// and *lost* 2.1× — Optimal 37 in 46.4 s / 478 645 nodes vs the no-restart
/// 22.1 s / 302 714, the restart problem presolving to 206×200 (no
/// reduction). pg, the shape the stall restart exists for, sits at 19.8 cuts/row
/// (2 473 on 125) — its nodes are expensive *because* of those rows, and a
/// presolve that strips them pays. The census danger families all sit ≤ 1.1
/// cuts/row (int_production s1 0.51, jobshop 0.44), so the bar excludes them
/// too. 2.0 keeps pg alone in the measured set.
pub const STALL_RESTART_CUT_DENSITY: f64 = 2.0;

/// Stall restart: an incumbent improvement smaller than this relative share
/// does not re-open the stall window — a tree whose incumbent only crawls
/// (pg: −8661.77 → −8668.54 = 0.008 % over 1 900 nodes) is stalled for
/// restart purposes exactly like one that stopped moving. The ladder that
/// motivates it (pg, branch point): improvements at nodes 51 / 173 / 1056 /
/// 2107 / 2633 — the first two significant (29.9 %, 0.28 %), the last three
/// 0.04–0.08 % — so with a strict any-improvement re-arm the window never
/// spans 1 000 nodes (largest gap 883) and the stall restart cannot fire on the one
/// instance shaped for it. 0.1 % keeps every significant move (a fresh
/// incumbent plan) re-arming while treating the crawl as the stall it is
/// (HiGHS breaks its root loop on the same logic from the effort side,
/// `moreHeuristicsAllowed`, HighsMipSolverData.cpp:2392–2400).
pub const STALL_IMPROVEMENT: f64 = 1e-3;

/// Integers fixed to a point by the root (global box equals at a value the
/// model's own bounds did not already force) — our stand-in for HiGHS
/// `percentageInactiveIntegers` (HighsMipSolverData.cpp:1310), which counts
/// clique-table substitutions too (we have none: bound fixings only).
pub fn inactive_integers(
    global_lower: &[f64],
    global_upper: &[f64],
    col_lower: &[f64],
    col_upper: &[f64],
    integer: &[bool],
) -> usize {
    (0..integer.len())
        .filter(|&j| {
            integer[j] && global_lower[j] == global_upper[j] && col_lower[j] != col_upper[j]
        })
        .count()
}

/// Counters and the incumbent carried from one search round into the next
/// across a restart. Per-subsystem stats (cuts, redcost, conflict) are NOT
/// merged across rounds: the final solution reports the last round's.
#[derive(Debug, Clone, Default)]
pub struct Carry {
    /// Incumbent already mapped into the next round's reduced space:
    /// (internal objective, point).
    pub incumbent: Option<(f64, Vec<f64>)>,
    pub nodes: u64,
    pub lp_iterations: u64,
    pub heur_lp_iterations: u64,
    /// Restarts already performed (this round is `restarts_done + 1`-th).
    pub restarts_done: u32,
    pub first_incumbent: Option<FirstIncumbent>,
    pub incumbents: Vec<FirstIncumbent>,
}

/// What the root search hands the presolve loop when it asks for a restart:
/// all in the *current* (pre-restart) reduced space.
#[derive(Debug, Clone)]
pub struct RestartRequest {
    /// Root box after propagation and reduced-cost fixing: the new column
    /// bounds (HiGHS applies the root domain to the LP before restarting).
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
    /// Cut rows tight at the final root LP point, to append to the model
    /// (`postSolveStack::appendCutsToModel`, HighsMipSolverData.cpp:1336).
    pub cuts: Vec<CutRow>,
    /// The incumbent to carry, in the current space (the presolve loop
    /// forward-maps it into the new reduced space).
    pub incumbent: Option<(f64, Vec<f64>)>,
    /// The final root LP objective: a valid dual bound to keep while the
    /// new search runs.
    pub bound: f64,
    /// Cumulative counters and history.
    pub carry: Carry,
}

/// The cuts of the root loop's log that are tight at the final root LP
/// point `x` (`a·x ≥ rhs − feastol`, the same activity test
/// `cuts::root_separation` uses for `drop_slack_at_end`, cuts/mod.rs:761–768).
/// Cuts aged out of the LP mid-loop but tight at the
/// end are included: they are globally valid, so re-adding can only cut.
pub fn active_cuts(log: &[(cuts::CutClass, CutRow)], x: &[f64], feastol: f64) -> Vec<CutRow> {
    log.iter()
        .filter_map(|(_, cut)| {
            let activity: f64 = cut
                .inds
                .iter()
                .zip(&cut.vals)
                .map(|(&j, &v)| v * x.get(j).copied().unwrap_or(0.0))
                .sum();
            (activity >= cut.rhs - feastol).then(|| cut.clone())
        })
        .collect()
}

/// The problem a restart re-presolves: the current problem with its column
/// bounds tightened to the root box and the active cut rows appended as
/// model rows (HiGHS `performRestart`: `appendCutsToModel` then the LP with
/// its root-domain bounds becomes the new model, :1331–1341).
pub fn restart_problem(p: &MipProblem, req: &RestartRequest) -> MipProblem {
    let ncut = req.cuts.len();
    let mut triplets: Vec<Triplet> =
        Vec::with_capacity(p.matrix.nnz() + req.cuts.iter().map(|c| c.inds.len()).sum::<usize>());
    for c in 0..p.matrix.ncol() {
        for k in p.matrix.col_ptr()[c]..p.matrix.col_ptr()[c + 1] {
            triplets.push((p.matrix.row_idx()[k] as usize, c, p.matrix.values()[k]));
        }
    }
    for (i, cut) in req.cuts.iter().enumerate() {
        triplets.extend(
            cut.inds
                .iter()
                .zip(&cut.vals)
                .map(|(&j, &v)| (p.num_row + i, j, v)),
        );
    }
    let matrix = CscMatrix::from_triplets(p.num_row + ncut, p.num_col, &triplets);
    let mut row_lower = p.row_lower.clone();
    let mut row_upper = p.row_upper.clone();
    row_lower.resize(p.num_row + ncut, f64::NEG_INFINITY);
    row_upper.extend(req.cuts.iter().map(|c| c.rhs));
    MipProblem {
        num_col: p.num_col,
        num_row: p.num_row + ncut,
        matrix,
        row_lower,
        row_upper,
        cost: p.cost.clone(),
        maximize: p.maximize,
        col_lower: req.lower.clone(),
        col_upper: req.upper.clone(),
        integer: p.integer.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cut(inds: &[usize], vals: &[f64], rhs: f64) -> CutRow {
        CutRow {
            inds: inds.to_vec(),
            vals: vals.to_vec(),
            rhs,
        }
    }

    #[test]
    fn inactive_integers_counts_only_root_fixings() {
        // Model bounds: [0,1] bin, [0,2] general int, [1,1] born fixed,
        // continuous ignored. Root fixes the bin to 1 and the continuous to
        // 0 — only the general integer and the bin count; the born-fixed
        // column must not (it would make every model trigger), the
        // continuous must not (not an integer).
        let col_lower = [0.0, 0.0, 1.0, 0.0];
        let col_upper = [1.0, 2.0, 1.0, 5.0];
        let integer = [true, true, true, false];
        let g_lower = [1.0, 1.0, 1.0, 0.0];
        let g_upper = [1.0, 1.0, 1.0, 0.0];
        assert_eq!(
            inactive_integers(&g_lower, &g_upper, &col_lower, &col_upper, &integer),
            2
        );
        // Nothing fixed beyond the model: zero even though a column is a
        // point by its own bounds.
        assert_eq!(
            inactive_integers(&col_lower, &col_upper, &col_lower, &col_upper, &integer),
            0
        );
    }

    #[test]
    fn active_cuts_keeps_only_rows_tight_at_the_point() {
        let log = vec![
            (cuts::CutClass::Clique, cut(&[0, 1], &[1.0, 1.0], 1.0)),
            (cuts::CutClass::Cmir, cut(&[0], &[2.0], 4.0)),
            (cuts::CutClass::Gomory, cut(&[1], &[1.0], 0.5)),
        ];
        // x = (1, 0): cut 0 tight (1 = 1), cut 1 slack (2 < 4), cut 2 slack
        // (0 < 0.5 − feastol).
        let kept = active_cuts(&log, &[1.0, 0.0], 1e-6);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].rhs, 1.0);
        // x = (2, −1.5): only cut 1 is tight (4 = 4); cut 0 falls under its
        // rhs (0.5 < 1) and cut 2 far under (−1.5 < 0.5).
        let kept = active_cuts(&log, &[2.0, -1.5], 1e-6);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].rhs, 4.0);
    }

    fn tiny_problem() -> MipProblem {
        // min x0 + x2  s.t.  x0 + x1 ≤ 1.5,  both binaries.
        MipProblem {
            num_col: 2,
            num_row: 1,
            matrix: CscMatrix::from_triplets(1, 2, &[(0, 0, 1.0), (0, 1, 1.0)]),
            row_lower: vec![f64::NEG_INFINITY],
            row_upper: vec![1.5],
            cost: vec![1.0, 1.0],
            maximize: false,
            col_lower: vec![0.0, 0.0],
            col_upper: vec![1.0, 1.0],
            integer: vec![true, true],
        }
    }

    /// End-to-end guard for the restart handshake (bnb → presolve chain →
    /// re-entered search → chained postsolve): a model whose root fixes an
    /// integer by reduced-cost fixing, so the first restart fires on the
    /// `fixingRate > 0` rule (verified with `RUSOLVE_MIP_LOG=1 cargo test
    /// -- --nocapture restart_end_to_end` at this commit). The restart-on
    /// answer must equal the restart-off answer on the original problem.
    #[test]
    fn restart_end_to_end_matches_no_restart() {
        // min x1+x2+x3 + 100(y1+y2)  s.t.  x1+x2+x3 >= 1,  x1 <= y1+y2,
        // all binary. LP optimum 34 1/3 at x=(1/3,1/3,1/3), y=(1/3,0);
        // the rounding incumbent is 101, so the internal objective is
        // integral and the cutoff row is y1+y2 <= 1 + x-overflow —
        // reduced-cost fixing can pin y2 and the first restart fires while
        // the root LP stays fractional.
        let p = MipProblem {
            num_col: 5,
            num_row: 2,
            matrix: CscMatrix::from_triplets(
                2,
                5,
                &[
                    (0, 0, 1.0),
                    (0, 1, 1.0),
                    (0, 2, 1.0),
                    (1, 0, 1.0),
                    (1, 3, -1.0),
                    (1, 4, -1.0),
                ],
            ),
            row_lower: vec![1.0, f64::NEG_INFINITY],
            row_upper: vec![f64::INFINITY, 0.0],
            cost: vec![1.0, 1.0, 1.0, 100.0, 100.0],
            maximize: false,
            col_lower: vec![0.0; 5],
            col_upper: vec![1.0; 5],
            integer: [true; 5].to_vec(),
        };
        let off = crate::presolve::solve_presolved(
            &p,
            &crate::bnb::MipConfig {
                restart: false,
                ..Default::default()
            },
        )
        .unwrap();
        let on = crate::presolve::solve_presolved(
            &p,
            &crate::bnb::MipConfig {
                restart: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(on.status, off.status);
        assert_eq!(on.status, crate::bnb::MipStatus::Optimal);
        assert!(
            (on.objective - off.objective).abs() <= 1e-6,
            "{} vs {}",
            on.objective,
            off.objective
        );
        assert_eq!(on.objective, 1.0); // x2 = 1, x1 = x3 = y = 0
        let v = p.violation(&on.x);
        assert!(
            v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6,
            "{v:?} {:?}",
            on.x
        );
    }

    /// Stall restart end-to-end: a
    /// market-split-shaped instance (harness `market_split` 2×6, seed 1)
    /// whose root fixes nothing (0 % inactive integers) but leaves cut mass
    /// (19 tight cuts on 2 rows) at a wide-open gap (1.0) — the pg/timtab1
    /// shape the fixing-rate rule cannot serve. With the stall bars at
    /// their floor the stall restart fires in-tree; `restarts >= 1` with zero
    /// root fixings pins the stall restart (the HiGHS rate rule needs fixings to fire),
    /// and the answer must equal the restart-off answer. Arming verified
    /// with RUSOLVE_MIP_LOG ("stall restart armed: 19 cuts tight,
    /// root gap 1.0000").
    #[test]
    fn stall_restart_restart_end_to_end_matches_no_restart() {
        let mps = "\
NAME market_split_2x6
OBJSENSE
    MIN
ROWS
 N  OBJ
 E  R0
 E  R1
COLUMNS
    M0  'MARKER'  'INTORG'
    X0  R0  17.0
    X0  R1  63.0
    X1  R0  72.0
    X1  R1  97.0
    X2  R0  97.0
    X2  R1  57.0
    X3  R0  8.0
    X3  R1  60.0
    X4  R0  32.0
    X4  R1  83.0
    X5  R0  15.0
    X5  R1  48.0
    M6  'MARKER'  'INTEND'
    X6  OBJ  1.0
    X6  R0  1.0
    X7  OBJ  1.0
    X7  R1  1.0
    X8  OBJ  1.0
    X8  R0  -1.0
    X9  OBJ  1.0
    X9  R1  -1.0
RHS
    RHS  R0  120.0
    RHS  R1  204.0
BOUNDS
 BV BND X0
 BV BND X1
 BV BND X2
 BV BND X3
 BV BND X4
 BV BND X5
ENDATA
";
        let p = MipProblem::from_model(&rusolve_io::read_mps(mps).unwrap());
        let off = crate::presolve::solve_presolved(
            &p,
            &crate::bnb::MipConfig {
                restart: false,
                ..Default::default()
            },
        )
        .unwrap();
        let on = crate::presolve::solve_presolved(
            &p,
            &crate::bnb::MipConfig {
                restart: true,
                restart_stall_nodes: 1,
                restart_stall_secs: 0.0,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(on.status, off.status);
        assert_eq!(on.status, crate::bnb::MipStatus::Optimal);
        assert!(
            (on.objective - off.objective).abs() <= 1e-6,
            "{} vs {}",
            on.objective,
            off.objective
        );
        assert!(on.restarts >= 1, "stall restart did not fire");
        let v = p.violation(&on.x);
        assert!(
            v.bound <= 1e-6 && v.row <= 1e-6 && v.integrality <= 1e-6,
            "{v:?} {:?}",
            on.x
        );
    }

    #[test]
    fn restart_problem_appends_cuts_and_tightens_bounds() {
        let p = tiny_problem();
        let req = RestartRequest {
            lower: vec![1.0, 0.0],
            upper: vec![1.0, 1.0],
            cuts: vec![cut(&[0, 1], &[1.0, -1.0], 0.0)],
            incumbent: None,
            bound: 0.0,
            carry: Carry::default(),
        };
        let q = restart_problem(&p, &req);
        assert_eq!(q.num_row, p.num_row + 1);
        assert_eq!(q.col_lower, vec![1.0, 0.0]);
        assert_eq!(q.col_upper, vec![1.0, 1.0]);
        assert_eq!(q.row_lower[p.num_row], f64::NEG_INFINITY);
        assert_eq!(q.row_upper[p.num_row], 0.0);
        // Both row sets are live: x = (1, 0) breaks the cut (1 > 0) but not
        // the model row; x = (1, 1) breaks the model row (2 > 1.5) only.
        assert!(q.violation(&[1.0, 0.0]).row >= 0.9);
        let v = q.violation(&[1.0, 1.0]).row;
        assert!((v - 0.5).abs() < 1e-9, "{v}");
    }
}
