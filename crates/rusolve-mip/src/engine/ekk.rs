//! `EkkEngine`: the [`LpEngine`] node-relaxation contract served by
//! `rusolve-ekk`, the faithful port of the HiGHS serial dual simplex, driven
//! through the ported public `Highs` API exactly the way HiGHS's own MIP
//! solver drives its `lpsolver` (`HighsLpRelaxation`,
//! reference/HiGHS/highs/mip/HighsLpRelaxation.cpp / .h, MIT; the API slice is
//! in crates/rusolve-ekk).
//!
//! Opt-in: `bnb::engine` builds it only when `RUSOLVE_MIP_EKK` is set and not
//! `0` ([`crate::engine_select::MipEngine`]); the default MIP path is
//! [`super::RusolveEngine`], unchanged.
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # Rows
//! The model rows are passed as HiGHS ranged rows, all `p.num_row` of them,
//! in model order (no L/G split, free rows kept), so LP row `i < m0` IS model
//! row `i` and the separator row index of cut `k` is `m0 + k`. Cuts are rows
//! `-inf <= a.x <= rhs` (`HighsLpRelaxation::addCuts`, :509-534).
//! `num_lp_rows = m0 + number of cuts`.
//!
//! # Hot start (the HiGHS rule, not the RusolveEngine one)
//! HiGHS never cold-starts a node LP. While diving, a child is solved from
//! whatever state the LP solver is in, after `flushDomain` (:710-730) changed
//! only the column bounds (`changeColsBounds` -> `updateStatus(kNewBounds)`:
//! invert, basis and DSE weights all survive). On a backtrack the search
//! calls `recoverBasis` -> `setBasis(stored node basis)` (HighsSearch.cpp:
//! 787-793, :1503-1517, HighsLpRelaxation.cpp:1107-1112), which keeps only the
//! statuses (refactorisation, DSE weights reset). A stored basis whose row
//! count no longer fits is ignored and the LP solver keeps its state
//! (HighsSearch.cpp:790). Mapped onto the trait:
//! * `current` is the snapshot the HiGHS state descends from: set to the
//!   new snapshot after every Optimal solve, cleared (`None`) after any other
//!   outcome, left alone by `add_cuts` / `retain_cuts` (cuts carry their ids
//!   and [`EkkEngine::adapt`] maps across them).
//! * `start == None`, or `start` pointer-equal to `current`: keep the HiGHS
//!   state, flush bounds only (dive).
//! * Otherwise: [`EkkEngine::adapt`] the snapshot to the current rows and
//!   `highs.set_basis(..)` (backtrack). When it does not fit (a removed cut
//!   row was nonbasic) or `set_basis` returns kError, keep the HiGHS state
//!   (the :790 analogue) and count `stats.fallback_unfit`.
//!
//! # Strong-branching probes (the Playground, HighsLpRelaxation.h:110-153)
//! The first probe after a node solve does `putIterate` then `flushDomain`;
//! every later probe from the same node does `flushDomain` then `getIterate`;
//! each runs `run(false)`. Closing the playground (C++ destructor:
//! `getIterate(); run()`) happens lazily in the next `solve_from` /
//! `add_cuts` / `retain_cuts`: `getIterate` only, because the caller's next
//! solve runs the LP anyway (one fewer LP per node; the HiGHS state after it
//! is the one the destructor's `run` would start from).
//!
//! # Options
//! Written straight to `highs.options_` (`setOptionValue` is not ported; its
//! only side effect is IIS bookkeeping). Constructor (:179-205): output off,
//! `random_seed`, primal tolerance `mip_feasibility_tolerance` (1e-6), dual
//! tolerance `0.1 *` that. Presolve: the HiGHS default for the first solve
//! (the root LP), `"off"` afterwards (HighsMipSolverData.cpp:2053-2055).

use super::{
    CutRow, LpDuals, LpEngine, LpError, LpOutcome, LpStats, LpStatus, TableauAgg, WarmStart,
};
use crate::problem::MipProblem;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::{
    HVector, HighsBasis, HighsBasisStatus, HighsCDouble, HighsInt, HighsLp, HighsModelStatus,
    HighsSparseMatrix, HighsSparseVectorSum, HighsStatus, MatrixFormat, ObjSense, cpp_max,
    kBasisValidityInvalid, kHighsIInf, kHighsInf, kHighsOffString, kHighsOnString,
    kSimplexStrategyDual, kSolutionStatusFeasible,
};
use rusolve_lp::{Basis, BasisStatus, SimplexError};
use std::rc::Rc;
use std::time::Instant;

/// `mipdata_->feastol` (= `mip_feasibility_tolerance`, HighsMipSolverData.cpp:755).
const FEASTOL: f64 = 1e-6;
/// `mipdata_->epsilon` (= `small_matrix_value`, HighsMipSolverData.cpp:756).
const EPSILON: f64 = 1e-9;

/// `HighsLpRelaxation::Status` (HighsLpRelaxation.h:40-50), the outcome of
/// [`EkkEngine::run`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RunStatus {
    Optimal,
    Infeasible,
    UnscaledDualFeasible,
    UnscaledPrimalFeasible,
    UnscaledInfeasible,
    Unbounded,
    Error,
}

pub struct EkkEngine {
    highs: Highs,
    /// Model rows (`m0`): LP rows `0..m0`; cut rows follow in `cuts` order.
    num_model_rows: usize,
    /// Row-wise copy of the model rows (`getRow` for the dual proof, :971).
    row_start: Vec<usize>,
    row_index: Vec<usize>,
    row_value: Vec<f64>,
    /// `getMaxAbsRowVal` per LP row (model rows, then cuts), maintained by
    /// `add_cuts` / `retain_cuts`.
    row_max_abs: Vec<f64>,
    /// The cuts currently held, in LP row order, and their ids (ids only
    /// grow, so `cut_ids` is sorted).
    cuts: Vec<CutRow>,
    cut_ids: Rc<[u64]>,
    next_cut_id: u64,
    /// The model box and integrality: the "global domain" of the proof's
    /// removal step (see [`Self::store_dual_inf_proof`]).
    global_lower: Vec<f64>,
    global_upper: Vec<f64>,
    integer: Vec<bool>,
    /// The snapshot the HiGHS state descends from (module docs).
    current: Option<Rc<WarmStart>>,
    /// Playground: an iterate is stored (`putIterate`) and `iterate_owner`
    /// is the `current` it was taken at.
    iterate_stored: bool,
    iterate_owner: Option<Rc<WarmStart>>,
    /// The first solve (the root LP) keeps the default presolve.
    solved_first: bool,
    iteration_limit: Option<u64>,
    deadline: Option<Instant>,
    stats: LpStats,
    /// `(coefs, rhs)` of the latest Infeasible solve's dual proof
    /// (`storeDualInfProof`); cleared by every solve.
    proof: Option<(Vec<f64>, f64)>,
}

/// HiGHS -> rusolve-lp basis status. Lossless on a solved basis: kZero (a
/// free nonbasic column at 0) has no rusolve-lp twin and travels as Fixed,
/// which [`to_highs_status`] maps back; kNonbasic never occurs after a solve
/// (`refineBasis` resolves it) and is mapped like kLower.
fn to_lp_status(s: HighsBasisStatus) -> BasisStatus {
    match s {
        HighsBasisStatus::kLower | HighsBasisStatus::kNonbasic => BasisStatus::AtLower,
        HighsBasisStatus::kBasic => BasisStatus::Basic,
        HighsBasisStatus::kUpper => BasisStatus::AtUpper,
        HighsBasisStatus::kZero => BasisStatus::Fixed,
    }
}

/// Inverse of [`to_lp_status`].
fn to_highs_status(s: BasisStatus) -> HighsBasisStatus {
    match s {
        BasisStatus::AtLower => HighsBasisStatus::kLower,
        BasisStatus::Basic => HighsBasisStatus::kBasic,
        BasisStatus::AtUpper => HighsBasisStatus::kUpper,
        BasisStatus::Fixed => HighsBasisStatus::kZero,
    }
}

/// The relaxation of `p` as a `HighsLp`: minimise `p.cost` (already in the
/// internal minimisation sense; `sense_ = kMinimize`, `offset_ = 0`),
/// columns in model order with `p.col_lower/col_upper`, rows
/// `p.row_lower <= A x <= p.row_upper` as given (ranged, none dropped), the
/// CSC matrix copied entry for entry into `a_matrix_` (`format_ =
/// kColwise`; `start_` as HighsInt, `num_col + 1` entries; `passModel`
/// drops tiny values itself). No integrality: the LP relaxation.
fn highs_lp_from_problem(p: &MipProblem) -> HighsLp {
    let num_col = p.num_col as HighsInt;
    let num_row = p.num_row as HighsInt;
    let a = &p.matrix;
    HighsLp {
        num_col_: num_col,
        num_row_: num_row,
        col_cost_: p.cost.clone(),
        col_lower_: p.col_lower.clone(),
        col_upper_: p.col_upper.clone(),
        row_lower_: p.row_lower.clone(),
        row_upper_: p.row_upper.clone(),
        sense_: ObjSense::kMinimize,
        offset_: 0.0,
        a_matrix_: HighsSparseMatrix {
            format_: MatrixFormat::kColwise,
            num_col_: num_col,
            num_row_: num_row,
            start_: a.col_ptr().iter().map(|&k| k as HighsInt).collect(),
            p_end_: Vec::new(),
            index_: a.row_idx().to_vec(),
            value_: a.values().to_vec(),
        },
        ..Default::default()
    }
}

impl EkkEngine {
    /// HiGHS: HighsLpRelaxation.cpp::HighsLpRelaxation (lines 179-205) +
    /// `loadModel` (`passModel`).
    ///
    /// `Highs::default()`, options as in the module docs (`output_flag =
    /// false`, then `options_.sync_log_options()`; `random_seed = 0`, the
    /// HiGHS MIP default), `highs.pass_model(highs_lp_from_problem(p))`
    /// (must return kOk; panics otherwise: the model came from a valid
    /// MipProblem). Builds the row-wise copy and `row_max_abs` (max |a_ij|
    /// per model row; 0 for an empty row), copies the box and `p.integer`.
    pub fn new(p: &MipProblem) -> Self {
        let mut highs = Highs::default();
        highs.options_.output_flag = false;
        highs.options_.sync_log_options();
        highs.options_.random_seed = 0;
        highs.options_.primal_feasibility_tolerance = FEASTOL;
        highs.options_.dual_feasibility_tolerance = FEASTOL * 0.1;
        let status = highs.pass_model(highs_lp_from_problem(p));
        assert!(
            status == HighsStatus::kOk,
            "EkkEngine: passModel failed: {status:?}"
        );

        let m0 = p.num_row;
        let a = &p.matrix;
        let (ptr, idx, val) = (a.col_ptr(), a.row_idx(), a.values());
        // Row-wise copy: prefix counts per row, then fill by ascending column.
        let mut row_start = vec![0usize; m0 + 1];
        for &i in idx {
            row_start[i as usize + 1] += 1;
        }
        for i in 0..m0 {
            row_start[i + 1] += row_start[i];
        }
        let mut next = row_start.clone();
        let nnz = val.len();
        let mut row_index = vec![0usize; nnz];
        let mut row_value = vec![0.0f64; nnz];
        for j in 0..p.num_col {
            for k in ptr[j]..ptr[j + 1] {
                let i = idx[k] as usize;
                row_index[next[i]] = j;
                row_value[next[i]] = val[k];
                next[i] += 1;
            }
        }
        let row_max_abs: Vec<f64> = (0..m0)
            .map(|i| {
                row_value[row_start[i]..row_start[i + 1]]
                    .iter()
                    .fold(0.0, |m: f64, v| cpp_max(m, v.abs()))
            })
            .collect();
        Self {
            highs,
            num_model_rows: m0,
            row_start,
            row_index,
            row_value,
            row_max_abs,
            cuts: Vec::new(),
            cut_ids: Rc::from(Vec::<u64>::new()),
            next_cut_id: 0,
            global_lower: p.col_lower.clone(),
            global_upper: p.col_upper.clone(),
            integer: p.integer.clone(),
            current: None,
            iterate_stored: false,
            iterate_owner: None,
            solved_first: false,
            iteration_limit: None,
            deadline: None,
            stats: LpStats::default(),
            proof: None,
        }
    }

    /// The HiGHS basis after an Optimal solve as a snapshot over the current
    /// rows: `highs.get_basis()` statuses through [`to_lp_status`], `dse:
    /// None` (faithful: `setBasis` drops the DSE weights), `cut_ids:
    /// self.cut_ids.clone()`.
    fn snapshot(&self) -> Rc<WarmStart> {
        let basis = self.highs.get_basis();
        Rc::new(WarmStart {
            basis: Basis {
                col_status: basis.col_status.iter().map(|&s| to_lp_status(s)).collect(),
                row_status: basis.row_status.iter().map(|&s| to_lp_status(s)).collect(),
            },
            dse: None,
            cut_ids: self.cut_ids.clone(),
        })
    }

    /// `ws` re-expressed over the current rows as a HiGHS basis (`valid =
    /// true`, `alien = false`, `useful = true`), or `None` when it does not
    /// fit. Same row rule as `RusolveEngine::adapt` (engine.rs:783-819): the
    /// column count must match and the model-row prefix is copied; cuts kept
    /// since the snapshot keep their status, cuts added after it are kBasic,
    /// and a removed cut that was not Basic makes it unfit. No column
    /// repair: `set_basis` checks `isBasisConsistent` and the bound flush
    /// that follows fixes nonbasic moves (`setNonbasicStatusInterface`).
    fn adapt(&self, ws: &WarmStart) -> Option<HighsBasis> {
        let m0 = self.num_model_rows;
        let snap = &ws.basis;
        let num_col = self.highs.get_num_col() as usize;
        if snap.col_status.len() != num_col || snap.row_status.len() != m0 + ws.cut_ids.len() {
            return None;
        }
        let row_status: Vec<BasisStatus> = if ws.cut_ids == self.cut_ids {
            snap.row_status.clone()
        } else {
            let old = &ws.cut_ids;
            let mut rows = Vec::with_capacity(m0 + self.cut_ids.len());
            rows.extend_from_slice(&snap.row_status[..m0]);
            let mut k = 0;
            for &id in self.cut_ids.iter() {
                while k < old.len() && old[k] < id {
                    if snap.row_status[m0 + k] != BasisStatus::Basic {
                        return None;
                    }
                    k += 1;
                }
                if k < old.len() && old[k] == id {
                    rows.push(snap.row_status[m0 + k]);
                    k += 1;
                } else {
                    rows.push(BasisStatus::Basic);
                }
            }
            if snap.row_status[m0 + k..]
                .iter()
                .any(|&s| s != BasisStatus::Basic)
            {
                return None;
            }
            rows
        };
        Some(HighsBasis {
            valid: true,
            alien: false,
            useful: true,
            col_status: snap
                .col_status
                .iter()
                .map(|&s| to_highs_status(s))
                .collect(),
            row_status: row_status.iter().map(|&s| to_highs_status(s)).collect(),
            ..Default::default()
        })
    }

    /// Positions the HiGHS state for a solve from `start` (module docs,
    /// "Hot start"): first closes an open playground (`highs.get_iterate()`,
    /// `iterate_stored = false`, `current = iterate_owner.take()`); then
    /// keep the state if `start` is `None` or `Rc::ptr_eq`-equal to
    /// `current` (compare `std::ptr::eq(start, &**current)`); otherwise
    /// `adapt` + `highs.set_basis(&b, "EkkEngine")`. Returns whether the
    /// solve is warm (any existing basis: always, except before the first
    /// solve) for the stats.
    fn sync_start(&mut self, start: Option<&WarmStart>) -> bool {
        if self.iterate_stored {
            self.highs.get_iterate();
            self.iterate_stored = false;
            self.current = self.iterate_owner.take();
        }
        let warm = self.solved_first;
        if let Some(s) = start {
            let same = matches!(&self.current, Some(c) if std::ptr::eq(s, &**c));
            if !same {
                let fits = match self.adapt(s) {
                    Some(b) => self.highs.set_basis(&b, "EkkEngine") != HighsStatus::kError,
                    None => false,
                };
                if !fits {
                    self.stats.fallback_unfit += 1;
                }
            }
        }
        warm
    }

    /// HiGHS: HighsLpRelaxation.cpp::flushDomain (lines 710-730)
    ///
    /// The set of columns whose bound in `highs.get_lp()` differs from
    /// `col_lower` / `col_upper` (`!=` on f64, ascending), passed in one
    /// `highs.change_cols_bounds_set(..)` with the new values. Nothing to
    /// change: no call (C++ `if (!domain.getChangedCols().empty())`).
    fn flush_domain(&mut self, col_lower: &[f64], col_upper: &[f64]) {
        let mut set: Vec<HighsInt> = Vec::new();
        let mut lower: Vec<f64> = Vec::new();
        let mut upper: Vec<f64> = Vec::new();
        {
            let lp = self.highs.get_lp();
            for (j, (&lo, &up)) in col_lower.iter().zip(col_upper).enumerate() {
                if lp.col_lower_[j] != lo || lp.col_upper_[j] != up {
                    set.push(j as HighsInt);
                    lower.push(lo);
                    upper.push(up);
                }
            }
        }
        if !set.is_empty() {
            self.highs
                .change_cols_bounds_set(set.len() as HighsInt, &set, &lower, &upper);
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::run (lines 1125-1441)
    ///
    /// * `options_.time_limit = max(get_run_time() + remaining, 0)` with
    ///   `remaining` the seconds to `deadline` (`kHighsInf` without one; the
    ///   C++ `time_limit - mipsolver.timer_.read()`);
    ///   `options_.simplex_iteration_limit` from `iter_cap` (`kHighsIInf`
    ///   when `None`, else `min(cap, HighsInt::MAX)`). Profiling, the IPM
    ///   `use_solver` choice and the timing log are dropped (simplex only).
    /// * `optimize_lp()`, then `solved_first = true`; once set, every later
    ///   run has `presolve = "off"` (written by `new` after the first run;
    ///   HighsMipSolverData.cpp:2053-2055).
    /// * Iterations: `max(0, info_.simplex_iteration_count)` (:1228).
    /// * `kError` (:1232-1256): `clear_solver()`; with `resolve_on_error`:
    ///   `simplex_strategy = kSimplexStrategyDual`, `presolve = "on"`,
    ///   `run(false, iter_cap)`, `presolve = "off"`, return its result (its
    ///   iterations added). Without: `recoverBasis` (no checkpoint here:
    ///   `current = None`) and Error. The primal-retry block is `#if 0` in
    ///   C++ and is not ported.
    /// * Model status (:1258-1440): kObjectiveBound -> Infeasible
    ///   (`storeDualUBProof` out of scope; [`Self::outcome`] tells Cutoff
    ///   apart by `get_model_status()`); kInfeasible ->
    ///   [`Self::store_dual_inf_proof`], Infeasible; kUnbounded -> Unbounded
    ///   (`trySolution` dropped); kUnknown with `info_.basis_validity ==
    ///   kBasisValidityInvalid` -> Error, else falls through to kOptimal:
    ///   Optimal iff `max_primal_infeasibility <= FEASTOL &&
    ///   max_dual_infeasibility <= FEASTOL`, else UnscaledPrimalFeasible,
    ///   UnscaledDualFeasible, UnscaledInfeasible (kOptimal only) or Error,
    ///   in that order; kIterationLimit -> Error (the IPM fallback under
    ///   `resolve_on_error` is out of scope); kTimeLimit and the rest ->
    ///   Error. [`Self::outcome`] reads the model status to map Error.
    fn run(&mut self, resolve_on_error: bool, iter_cap: Option<u64>) -> (RunStatus, u64) {
        let remaining = match self.deadline {
            Some(d) => {
                let now = Instant::now();
                if d >= now {
                    (d - now).as_secs_f64()
                } else {
                    -(now - d).as_secs_f64()
                }
            }
            None => kHighsInf,
        };
        self.highs.options_.time_limit = cpp_max(self.highs.get_run_time() + remaining, 0.0);
        self.highs.options_.simplex_iteration_limit = match iter_cap {
            None => kHighsIInf,
            Some(cap) => cap.min(HighsInt::MAX as u64) as HighsInt,
        };
        let call_status = self.highs.optimize_lp();
        self.solved_first = true;
        self.highs.options_.presolve = kHighsOffString.to_string();
        let iterations = self.highs.get_info().simplex_iteration_count.max(0) as u64;
        if call_status == HighsStatus::kError {
            self.highs.clear_solver();
            if resolve_on_error {
                self.highs.options_.simplex_strategy = kSimplexStrategyDual;
                self.highs.options_.presolve = kHighsOnString.to_string();
                let (status, more) = self.run(false, iter_cap);
                self.highs.options_.presolve = kHighsOffString.to_string();
                return (status, iterations + more);
            }
            return (RunStatus::Error, iterations);
        }
        let model_status = self.highs.get_model_status();
        match model_status {
            HighsModelStatus::kObjectiveBound => (RunStatus::Infeasible, iterations),
            HighsModelStatus::kInfeasible => {
                self.store_dual_inf_proof();
                (RunStatus::Infeasible, iterations)
            }
            HighsModelStatus::kUnbounded => (RunStatus::Unbounded, iterations),
            HighsModelStatus::kUnknown | HighsModelStatus::kOptimal => {
                let info = self.highs.get_info();
                if model_status == HighsModelStatus::kUnknown
                    && info.basis_validity == kBasisValidityInvalid
                {
                    return (RunStatus::Error, iterations);
                }
                let (primal, dual) = (info.max_primal_infeasibility, info.max_dual_infeasibility);
                let status = if primal <= FEASTOL && dual <= FEASTOL {
                    RunStatus::Optimal
                } else if primal <= FEASTOL {
                    RunStatus::UnscaledPrimalFeasible
                } else if dual <= FEASTOL {
                    RunStatus::UnscaledDualFeasible
                } else if model_status == HighsModelStatus::kOptimal {
                    RunStatus::UnscaledInfeasible
                } else {
                    RunStatus::Error
                };
                (status, iterations)
            }
            _ => (RunStatus::Error, iterations),
        }
    }

    /// HiGHS: HighsLpRelaxation.cpp::storeDualInfProof (lines 928-1046)
    ///
    /// Sets `self.proof`. `highs.get_dual_ray_sparse(&mut has, &mut row_ep)`
    /// (row_ep `setup(num_row)`); no ray -> `None`. Then verbatim: `weight =
    /// -row_ep.array[iRow]` over `row_ep.index[..count]`, skip `|weight| *
    /// row_max_abs[iRow] <= EPSILON`, the bound by the sign of weight (skip
    /// an infinite one), `upper: HighsCDouble += weight * bound`, row
    /// entries (model rows from the row-wise copy, cut `k = iRow - m0` from
    /// `cuts[k]`) into a `HighsSparseVectorSum` `row_ap` over `num_col`.
    /// The removal loop runs as in C++ with the global domain = the model
    /// box (`global_lower/upper`) and `isColContinuous(i) = !integer[i]`
    /// (valid in every node: the box only shrinks). Kept entries go into a
    /// dense `coefs[num_col]`; `rhs = f64::from(upper)`.
    /// `tightenCoefficients`, `checkCut` and `extractCliquesFromCut` are out
    /// of scope.
    fn store_dual_inf_proof(&mut self) {
        self.proof = None;
        if self.highs.get_info().basis_validity == kBasisValidityInvalid {
            return;
        }
        let num_row = self.highs.get_num_row();
        let num_col = self.highs.get_num_col();
        let mut row_ep = HVector::default();
        row_ep.setup(num_row);
        let mut has_ray = false;
        self.highs.get_dual_ray_sparse(&mut has_ray, &mut row_ep);
        if !has_ray {
            return;
        }
        let m0 = self.num_model_rows;
        let mut row_ap = HighsSparseVectorSum::new(num_col);
        let mut upper = HighsCDouble::from(0.0);
        for i in 0..row_ep.count as usize {
            let i_row = row_ep.index[i] as usize;
            let weight = -row_ep.array[i_row];
            if weight.abs() * self.row_max_abs[i_row] <= EPSILON {
                continue;
            } else if weight > 0.0 {
                let row_upper = self.highs.get_lp().row_upper_[i_row];
                if row_upper == kHighsInf {
                    continue;
                }
                upper += weight * row_upper;
            } else {
                let row_lower = self.highs.get_lp().row_lower_[i_row];
                if row_lower == -kHighsInf {
                    continue;
                }
                upper += weight * row_lower;
            }
            if i_row < m0 {
                for k in self.row_start[i_row]..self.row_start[i_row + 1] {
                    row_ap.add(self.row_index[k] as HighsInt, weight * self.row_value[k]);
                }
            } else {
                let cut = &self.cuts[i_row - m0];
                for (&j, &v) in cut.inds.iter().zip(&cut.vals) {
                    row_ap.add(j as HighsInt, weight * v);
                }
            }
        }
        let lp = self.highs.get_lp();
        let mut coefs = vec![0.0; num_col as usize];
        for &i in row_ap.get_nonzeros() {
            let val = row_ap.get_value(i);
            if val.abs() <= EPSILON {
                continue;
            }
            let col = i as usize;
            let (glo, gup) = (self.global_lower[col], self.global_upper[col]);
            let mut remove_value = val.abs() <= FEASTOL;
            if !remove_value && (glo == gup || !self.integer[col]) {
                remove_value = if val > 0.0 {
                    lp.col_lower_[col] - glo <= FEASTOL
                } else {
                    gup - lp.col_upper_[col] <= FEASTOL
                };
            }
            if remove_value {
                if val < 0.0 {
                    if gup == kHighsInf {
                        self.proof = None;
                        return;
                    }
                    upper -= val * gup;
                } else {
                    if glo == -kHighsInf {
                        self.proof = None;
                        return;
                    }
                    upper -= val * glo;
                }
                continue;
            }
            coefs[col] = val;
        }
        let rhs = f64::from(upper);
        self.proof = Some((coefs, rhs));
    }

    /// The shared tail of `solve_from` / `probe`: maps the HiGHS outcome to
    /// the trait. kObjectiveBound -> `Cutoff` (objective =
    /// `info_.objective_function_value`, a valid bound; `stats.cutoffs`);
    /// Infeasible -> `Infeasible` (`stats.infeasible_trusted` when a proof
    /// was stored); Optimal / UnscaledPrimalFeasible -> `Optimal` with `x =
    /// solution_.col_value`, objective `info_.objective_function_value`,
    /// `basis = Some(snapshot())`, and `current` = that snapshot unless
    /// `probe`; Unbounded -> `Unbounded`; kIterationLimit / kTimeLimit ->
    /// `Ok(IterationLimit)` for a probe (objective: the dual bound, i.e.
    /// `objective_function_value` when `dual_solution_status` is feasible,
    /// else `-inf`),
    /// `Err(LpError::IterationLimit)` otherwise; any other Error /
    /// UnscaledDualFeasible / UnscaledInfeasible ->
    /// `Err(LpError::Simplex(SimplexError::Unsupported(<model status>)))`.
    /// Every non-Optimal outcome of a non-probe sets `current = None`.
    fn outcome(
        &mut self,
        status: RunStatus,
        iterations: u64,
        probe: bool,
    ) -> Result<LpOutcome, LpError> {
        let model_status = self.highs.get_model_status();
        let out = match status {
            RunStatus::Infeasible if model_status == HighsModelStatus::kObjectiveBound => {
                self.stats.cutoffs += 1;
                Ok(LpOutcome {
                    status: LpStatus::Cutoff,
                    objective: self.highs.get_info().objective_function_value,
                    x: Vec::new(),
                    iterations,
                    basis: None,
                })
            }
            RunStatus::Infeasible => {
                if self.proof.is_some() {
                    self.stats.infeasible_trusted += 1;
                }
                Ok(LpOutcome {
                    status: LpStatus::Infeasible,
                    objective: f64::INFINITY,
                    x: Vec::new(),
                    iterations,
                    basis: None,
                })
            }
            RunStatus::Optimal | RunStatus::UnscaledPrimalFeasible => {
                let snapshot = self.snapshot();
                if !probe {
                    self.current = Some(Rc::clone(&snapshot));
                }
                Ok(LpOutcome {
                    status: LpStatus::Optimal,
                    objective: self.highs.get_info().objective_function_value,
                    x: self.highs.get_solution().col_value.clone(),
                    iterations,
                    basis: Some(snapshot),
                })
            }
            RunStatus::Unbounded => Ok(LpOutcome {
                status: LpStatus::Unbounded,
                objective: f64::NEG_INFINITY,
                x: Vec::new(),
                iterations,
                basis: None,
            }),
            RunStatus::Error
                if matches!(
                    model_status,
                    HighsModelStatus::kIterationLimit | HighsModelStatus::kTimeLimit
                ) =>
            {
                if probe {
                    // The trait wants the dual bound reached. After a
                    // capped run `returnFromSolve` has removed the cost
                    // perturbations and recomputed the duals, so a
                    // dual-feasible basic solution's objective (= its dual
                    // objective) is a valid lower bound; otherwise none.
                    let info = self.highs.get_info();
                    let objective = if info.dual_solution_status == kSolutionStatusFeasible
                        && info.objective_function_value.is_finite()
                    {
                        info.objective_function_value
                    } else {
                        f64::NEG_INFINITY
                    };
                    Ok(LpOutcome {
                        status: LpStatus::IterationLimit,
                        objective,
                        x: Vec::new(),
                        iterations,
                        basis: None,
                    })
                } else {
                    Err(LpError::IterationLimit)
                }
            }
            _ => Err(LpError::Simplex(SimplexError::Unsupported(format!(
                "{model_status:?}"
            )))),
        };
        if !probe && !matches!(&out, Ok(o) if o.status == LpStatus::Optimal) {
            self.current = None;
        }
        out
    }
}

impl LpEngine for EkkEngine {
    fn solve(&mut self, col_lower: &[f64], col_upper: &[f64]) -> Result<LpOutcome, LpError> {
        self.solve_from(col_lower, col_upper, None, f64::INFINITY)
    }

    /// `sync_start(start)`; `options_.objective_bound = cutoff`
    /// (`setObjectiveLimit`, :1114-1123: the search's cutoff already holds
    /// its integrality offset, so no second offset); `flush_domain`;
    /// `run(true, self.iteration_limit)` (`resolveLp` / `evaluateNode` call
    /// `run()` with its default `resolve_on_error = true`); `outcome(..,
    /// false)`. Clears `proof` first. Stats: `warm`/`cold` and
    /// `warm_iterations`/`cold_iterations` from `sync_start`'s answer.
    fn solve_from(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
    ) -> Result<LpOutcome, LpError> {
        self.proof = None;
        let warm = self.sync_start(start);
        self.highs.options_.objective_bound = cutoff;
        self.flush_domain(col_lower, col_upper);
        let (status, iterations) = self.run(true, self.iteration_limit);
        if warm {
            self.stats.warm += 1;
            self.stats.warm_iterations += iterations;
        } else {
            self.stats.cold += 1;
            self.stats.cold_iterations += iterations;
        }
        self.outcome(status, iterations, false)
    }

    /// Playground (module docs): with an iterate stored and `start` the
    /// owner (`None`, or pointer-equal to `iterate_owner`): `flush_domain`
    /// then `highs.get_iterate()`. Otherwise `sync_start(start)`, then
    /// `highs.put_iterate()` (kOk only with an invert; without one run the
    /// probe outside a playground and set `current = None` afterwards),
    /// `iterate_owner = current.clone()`, `iterate_stored = true`, then
    /// `flush_domain`. Then `objective_bound = cutoff`,
    /// `run(false, Some(iter_limit))`, `outcome(.., true)`. Stats: `probes`,
    /// `probe_iterations`, `probe_iter_limits`.
    fn probe(
        &mut self,
        col_lower: &[f64],
        col_upper: &[f64],
        start: Option<&WarmStart>,
        cutoff: f64,
        iter_limit: u64,
    ) -> Result<LpOutcome, LpError> {
        self.proof = None;
        let reuse = self.iterate_stored
            && match start {
                None => true,
                Some(s) => matches!(&self.iterate_owner, Some(o) if std::ptr::eq(s, &**o)),
            };
        let (warm, stored) = if reuse {
            self.flush_domain(col_lower, col_upper);
            self.highs.get_iterate();
            (true, true)
        } else {
            let warm = self.sync_start(start);
            let stored = self.highs.put_iterate() == HighsStatus::kOk;
            if stored {
                self.iterate_owner = self.current.clone();
                self.iterate_stored = true;
            }
            self.flush_domain(col_lower, col_upper);
            (warm, stored)
        };
        self.highs.options_.objective_bound = cutoff;
        let (status, iterations) = self.run(false, Some(iter_limit));
        if !stored {
            self.current = None;
        }
        if warm {
            self.stats.warm += 1;
            self.stats.warm_iterations += iterations;
        } else {
            self.stats.cold += 1;
            self.stats.cold_iterations += iterations;
        }
        let out = self.outcome(status, iterations, true);
        if let Ok(o) = &out {
            self.stats.probes += 1;
            self.stats.probe_iterations += o.iterations;
            if o.status == LpStatus::IterationLimit {
                self.stats.probe_iter_limits += 1;
            }
        }
        out
    }

    /// `setIterationLimit` (HighsLpRelaxation.h:370): applied in `run`.
    fn set_iteration_limit(&mut self, limit: Option<u64>) {
        self.iteration_limit = limit;
    }

    /// Applied in `run` as the HiGHS `time_limit` (:1127-1131).
    fn set_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.deadline = deadline;
    }

    /// HiGHS: HighsLpRelaxation.cpp::addCuts (lines 509-534)
    ///
    /// Closes an open playground first (`get_iterate`, as `sync_start`).
    /// One `highs.add_rows(n, lower = -kHighsInf each, upper = rhs, nnz,
    /// starts, inds, vals)` for all cuts (CSR; entries in the cut's own
    /// order, as `addCuts` copies the cut pool rows); ids from
    /// `next_cut_id`; `row_max_abs` extended. Empty `cuts`: return true.
    fn add_cuts(&mut self, cuts: &[CutRow]) -> bool {
        if cuts.is_empty() {
            return true;
        }
        if self.iterate_stored {
            self.highs.get_iterate();
            self.iterate_stored = false;
            self.current = self.iterate_owner.take();
        }
        let mut starts: Vec<HighsInt> = Vec::with_capacity(cuts.len());
        let mut indices: Vec<HighsInt> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        let mut lower: Vec<f64> = Vec::with_capacity(cuts.len());
        let mut upper: Vec<f64> = Vec::with_capacity(cuts.len());
        let mut ids = self.cut_ids.to_vec();
        for cut in cuts {
            starts.push(indices.len() as HighsInt);
            for (&j, &v) in cut.inds.iter().zip(&cut.vals) {
                indices.push(j as HighsInt);
                values.push(v);
            }
            lower.push(-kHighsInf);
            upper.push(cut.rhs);
            self.row_max_abs
                .push(cut.vals.iter().fold(0.0, |m: f64, v| cpp_max(m, v.abs())));
            self.cuts.push(cut.clone());
            ids.push(self.next_cut_id);
            self.next_cut_id += 1;
        }
        let added = self.highs.add_rows(
            cuts.len() as HighsInt,
            &lower,
            &upper,
            indices.len() as HighsInt,
            &starts,
            &indices,
            &values,
        );
        debug_assert!(added == HighsStatus::kOk);
        self.cut_ids = Rc::from(ids);
        true
    }

    /// HiGHS: HighsLpRelaxation.cpp::removeCuts (lines 581-608)
    ///
    /// Closes an open playground first. Nothing dropped: return. Else
    /// `basis = highs.get_basis().clone()`; mask over all LP rows (model
    /// rows 0, dropped cuts 1); `highs.delete_rows_mask(&mut mask)`;
    /// compact `basis.row_status` by the returned mask (`mask[i] != -1`
    /// keeps row `i` at `mask[i]`); `highs.set_basis(&basis, "")`; then
    /// `highs.optimize_lp()` immediately, exactly as C++ (its iterations go
    /// to `stats.warm_iterations`). Drops the cuts, ids and
    /// `row_max_abs` entries (`v.retain(|x| !dropped(x))`, the negation is
    /// mandatory). `current` unchanged.
    fn retain_cuts(&mut self, keep: &[bool]) {
        if self.iterate_stored {
            self.highs.get_iterate();
            self.iterate_stored = false;
            self.current = self.iterate_owner.take();
        }
        if keep.iter().all(|&k| k) {
            return;
        }
        let m0 = self.num_model_rows;
        debug_assert_eq!(keep.len(), self.cuts.len());
        let nlprows = self.highs.get_num_row() as usize;
        let ndel = keep.iter().filter(|&&k| !k).count();
        let mut basis = self.highs.get_basis().clone();
        let mut mask: Vec<HighsInt> = vec![0; nlprows];
        for (k, &kept) in keep.iter().enumerate() {
            if !kept {
                mask[m0 + k] = 1;
            }
        }
        self.highs.delete_rows_mask(&mut mask);
        // A basis that does not span the LP's rows (no valid basis) is left
        // as it is; set_basis then refuses it and optimize_lp starts cold.
        if basis.row_status.len() == nlprows {
            for (i, &m) in mask.iter().enumerate() {
                if m >= 0 {
                    basis.row_status[m as usize] = basis.row_status[i];
                }
            }
            let new_len = nlprows - ndel;
            basis.row_status.truncate(new_len);
        }
        self.highs
            .set_basis(&basis, "HighsLpRelaxation::removeCuts");
        self.highs.optimize_lp();
        self.stats.warm_iterations += self.highs.get_info().simplex_iteration_count.max(0) as u64;
        let mut k = 0;
        self.cuts.retain(|_| {
            let kept = keep[k];
            k += 1;
            kept
        });
        let mut i = 0;
        self.row_max_abs.retain(|_| {
            let kept = i < m0 || keep[i - m0];
            i += 1;
            kept
        });
        let ids: Vec<u64> = self
            .cut_ids
            .iter()
            .zip(keep)
            .filter(|(_, kept)| **kept)
            .map(|(&id, _)| id)
            .collect();
        self.cut_ids = Rc::from(ids);
    }

    fn num_lp_rows(&self) -> usize {
        self.highs.get_num_row() as usize
    }

    fn stats(&self) -> LpStats {
        self.stats
    }

    /// `None` unless `lp` is Optimal, no playground is open, its basis is
    /// pointer-equal to `current` and `solution_.dual_valid` (the HiGHS
    /// solution is that solve's). `y =
    /// solution_.row_dual` (HiGHS: `d = c - A^T y`, min sense, `y >= 0` at a
    /// lower row bound); clip a wrong-signed entry to 0 (`y > 0` needs a
    /// finite `row_lower`, `y < 0` a finite `row_upper`). Built field by
    /// field as `LpDuals::from_rows` (engine.rs:225) but over ranged rows:
    /// `d_j = c_j - sum_i y_i a_ij` and `scale_j = |c_j| + sum |y_i a_ij|`
    /// over the model columns plus the cut rows, `row_term = sum (y_i > 0 ?
    /// y_i * row_lower_i : y_i * row_upper_i)`, `row_scale` the sum of the
    /// absolute terms. `stats.duals_lp`.
    fn duals(&mut self, lp: &LpOutcome) -> Option<LpDuals> {
        // An open playground holds a probe's solution while `current` still
        // names the node's basis: serving it would hand out the probe's duals.
        if lp.status != LpStatus::Optimal || self.iterate_stored {
            return None;
        }
        let ws = lp.basis.as_ref()?;
        if !matches!(&self.current, Some(c) if Rc::ptr_eq(c, ws)) {
            return None;
        }
        let nrows = self.highs.get_num_row() as usize;
        let ncol = self.highs.get_num_col() as usize;
        let model = self.highs.get_lp();
        let sol = self.highs.get_solution();
        if !sol.dual_valid || sol.row_dual.len() != nrows {
            return None;
        }
        let mut y: Vec<f64> = sol.row_dual.clone();
        for (i, yi) in y.iter_mut().enumerate() {
            if (*yi > 0.0 && !model.row_lower_[i].is_finite())
                || (*yi < 0.0 && !model.row_upper_[i].is_finite())
            {
                *yi = 0.0;
            }
        }
        let cost = &model.col_cost_;
        let (start, index, value) = (
            &model.a_matrix_.start_,
            &model.a_matrix_.index_,
            &model.a_matrix_.value_,
        );
        let mut reduced_costs = Vec::with_capacity(ncol);
        let mut scale = Vec::with_capacity(ncol);
        for (j, &cj) in cost.iter().enumerate() {
            let (mut dot, mut mag) = (0.0, 0.0);
            for k in start[j] as usize..start[j + 1] as usize {
                let t = value[k] * y[index[k] as usize];
                dot += t;
                mag += t.abs();
            }
            reduced_costs.push(cj - dot);
            scale.push(cj.abs() + mag);
        }
        let (mut row_term, mut row_scale) = (0.0, 0.0);
        for (i, &yi) in y.iter().enumerate() {
            // A zero multiplier contributes nothing (0 * infinity would be NaN).
            let term = if yi > 0.0 {
                yi * model.row_lower_[i]
            } else if yi < 0.0 {
                yi * model.row_upper_[i]
            } else {
                0.0
            };
            row_term += term;
            row_scale += term.abs();
        }
        self.stats.duals_lp += 1;
        Some(LpDuals {
            reduced_costs,
            row_duals: y,
            scale,
            row_term,
            row_scale,
        })
    }

    /// The cut-row statuses of `adapt(lp.basis)` (`None` when unfit or not
    /// Optimal): `row_status[m0..] == kBasic` (HighsLpRelaxation.cpp:540,
    /// :659).
    fn cut_rows_basic(&mut self, lp: &LpOutcome) -> Option<Vec<bool>> {
        if lp.status != LpStatus::Optimal {
            return None;
        }
        let ws = lp.basis.as_ref()?;
        let basis = self.adapt(ws)?;
        let m0 = self.num_model_rows;
        Some(
            basis.row_status[m0..]
                .iter()
                .map(|&s| s == HighsBasisStatus::kBasic)
                .collect(),
        )
    }

    fn infeasibility_proof(&self) -> Option<(Vec<f64>, f64)> {
        self.proof.clone()
    }

    /// HiGHS: HighsTableauSeparator.cpp:131-201 (the BTRAN half).
    ///
    /// `None` unless `lp` is Optimal, `lp.basis` is pointer-equal to
    /// `current`, no iterate is stored, and `highs.has_invert()`. For each
    /// `j` in `cols` with a position `pos` in `get_basic_variables_array()`
    /// (`== j`; build a column -> position map once): `row_ep.setup(num_row)`,
    /// `highs.get_basis_inverse_row_sparse(pos, &mut row_ep)`, weights =
    /// `(i, row_ep.array[i])` over `row_ep.index[..count]` with nonzero
    /// value, sorted by `i` (LP row `i` is separator row `i`: model rows by
    /// index, cut `k` at `m0 + k`). Nonbasic `j` omitted. The sign
    /// convention must pass the RusolveEngine test
    /// `tableau_weights_are_the_tableau_row_of_each_candidate`
    /// (engine.rs:2053) run against this engine; HiGHS logicals are `+e_i`
    /// columns, like the RusolveEngine slack the weights were defined for.
    fn tableau_weights(
        &mut self,
        lp: &LpOutcome,
        lower: &[f64],
        upper: &[f64],
        cols: &[usize],
    ) -> Option<Vec<TableauAgg>> {
        if lp.status != LpStatus::Optimal || self.iterate_stored {
            return None;
        }
        let ws = lp.basis.as_ref()?;
        if !matches!(&self.current, Some(c) if Rc::ptr_eq(c, ws)) {
            return None;
        }
        if !self.highs.has_invert() {
            return None;
        }
        // The basis inverse rows depend on the basis only, not on the bounds.
        let _ = (lower, upper);
        let num_row = self.highs.get_num_row();
        let num_col = self.highs.get_num_col() as usize;
        let total = num_col + num_row as usize;
        // Column -> position in the basic variables (logicals sit at >= num_col).
        let mut pos_of: Vec<HighsInt> = vec![-1; total];
        for (pos, &j) in self.highs.get_basic_variables_array().iter().enumerate() {
            if let Some(slot) = pos_of.get_mut(j as usize) {
                *slot = pos as HighsInt;
            }
        }
        let mut row_ep = HVector::default();
        row_ep.setup(num_row);
        let mut out = Vec::with_capacity(cols.len());
        for &j in cols {
            let pos = pos_of.get(j).copied().unwrap_or(-1);
            if pos < 0 {
                continue;
            }
            self.highs.get_basis_inverse_row_sparse(pos, &mut row_ep);
            let mut weights: Vec<(usize, f64)> = Vec::with_capacity(row_ep.count as usize);
            for &i in &row_ep.index[..row_ep.count as usize] {
                let i = i as usize;
                let v = row_ep.array[i];
                if v != 0.0 {
                    weights.push((i, v));
                }
            }
            weights.sort_unstable_by_key(|&(i, _)| i);
            out.push(TableauAgg { col: j, weights });
        }
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_io::read_mps;

    /// The LADDER LP of engine.rs's tests: root optimum fractional.
    const LADDER: &str = "\
NAME LADDER
OBJSENSE
 MAX
ROWS
 N  OBJ
 L  R1
 L  R2
 L  R3
COLUMNS
    X0  OBJ  3  R1  1
    X0  R2  4  R3  2
    X1  OBJ  5  R1  1
    X1  R2  3  R3  5
    X2  OBJ  4  R1  2
    X2  R2  1  R3  3
    X3  OBJ  2  R1  1
    X3  R2  2  R3  1
RHS
    RHS  R1  12  R2  18
    RHS  R3  15
BOUNDS
 UP BND X0 6
 UP BND X1 6
 UP BND X2 6
 UP BND X3 6
ENDATA
";

    fn problem() -> MipProblem {
        MipProblem::from_model(&read_mps(LADDER).unwrap())
    }

    /// Bounds of the two children of the first fractional column of `x`.
    fn children(p: &MipProblem, x: &[f64]) -> [(Vec<f64>, Vec<f64>); 2] {
        let j = (0..p.num_col)
            .find(|&j| (x[j] - x[j].round()).abs() > 1e-6)
            .expect("fractional root");
        let mut down = (p.col_lower.clone(), p.col_upper.clone());
        down.1[j] = x[j].floor();
        let mut up = (p.col_lower.clone(), p.col_upper.clone());
        up.0[j] = x[j].ceil();
        [down, up]
    }

    /// A probe capped before it can finish reports a finite dual bound no
    /// larger than the child's optimum (a bound change keeps the parent basis
    /// dual feasible), not `-inf`: strong branching scores on it.
    #[test]
    fn capped_probe_reports_a_finite_valid_bound() {
        let p = problem();
        let mut e = EkkEngine::new(&p);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        assert_eq!(root.status, LpStatus::Optimal);
        let ws = root.basis.clone().unwrap();
        let mut capped_seen = 0;
        for (lo, hi) in children(&p, &root.x) {
            let exact = EkkEngine::new(&p).solve(&lo, &hi).unwrap();
            let capped = e.probe(&lo, &hi, Some(&ws), f64::INFINITY, 0).unwrap();
            match capped.status {
                LpStatus::IterationLimit => {
                    capped_seen += 1;
                    assert!(capped.objective.is_finite(), "{capped:?}");
                    assert!(capped.objective >= root.objective - 1e-9);
                    if exact.status == LpStatus::Optimal {
                        assert!(capped.objective <= exact.objective + 1e-9);
                    }
                }
                LpStatus::Optimal => {
                    assert!((capped.objective - exact.objective).abs() < 1e-9)
                }
                s => panic!("{s:?}"),
            }
        }
        assert!(capped_seen >= 1);
    }

    /// While a probe's playground is open the HiGHS solution is the probe's,
    /// so the node's duals must not be served from it.
    #[test]
    fn duals_are_refused_after_a_probe_and_valid_after_a_solve() {
        let p = problem();
        let mut e = EkkEngine::new(&p);
        let root = e.solve(&p.col_lower, &p.col_upper).unwrap();
        let ws = root.basis.clone().unwrap();
        let d = e.duals(&root).expect("root duals");
        let bound = d.bound(&p.col_lower, &p.col_upper).expect("finite box");
        assert!(
            (bound - root.objective).abs() < 1e-8,
            "{bound} vs {}",
            root.objective
        );
        let [(lo, hi), _] = children(&p, &root.x);
        let _ = e.probe(&lo, &hi, Some(&ws), f64::INFINITY, 1000).unwrap();
        assert!(e.duals(&root).is_none());
        let again = e
            .solve_from(&p.col_lower, &p.col_upper, Some(&ws), f64::INFINITY)
            .unwrap();
        let d2 = e.duals(&again).expect("duals after re-solve");
        let b2 = d2.bound(&p.col_lower, &p.col_upper).unwrap();
        assert!((b2 - again.objective).abs() < 1e-8);
    }
}
