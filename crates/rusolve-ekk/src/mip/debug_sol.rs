//! HiGHS source: mip/HighsDebugSol.h, mip/HighsDebugSol.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The `HIGHS_DEBUGSOL` build of HiGHS, switched at run time: everything here is inert unless
//! the `mip_debug_solution_file` option names a solution file (HiGHS reads the same option in
//! `activate`). The default (empty) option leaves every method a no-op, as the HiGHS build
//! without `HIGHS_DEBUGSOL`.
//!
//! PORT NOTE: HiGHS `assert`s on a violated check. The port reports each violation on stderr
//! (`HighsDebugSol: ...`, the first ones with a backtrace) and continues, so that one run
//! localises the first step that cuts the debug solution off; `RUSOLVE_DEBUG_SOL_ABORT=1`
//! panics instead, as the assert.
//!
//! PORT NOTE: HiGHS keys `conflictingBounds` by `&domain`. Port domains move (`detach!`,
//! worker vectors), so a pointer key is not stable. The multiset is non-empty exactly when
//! the domain's current bounds exclude the debug solution (bounds only tighten on the stack,
//! and backtracking pops the offending change with its entry), so the port evaluates that
//! condition directly on the domain bounds instead of keeping the multiset. Unlike HiGHS,
//! unregistered domains (heuristic copies) are checked too; each report names the reason
//! type so such reports can be told apart.
//!
//! PORT NOTE: deriving the multiset condition from the domain bounds covers every change
//! that is *already on the stack*, but not the one that is currently being applied. HiGHS
//! exploits that with the `branching == true` call of
//! `HighsDomain::setDomainChangeStack` (HighsDomain.cpp:2467-2468 at 1.16.0, 2225-2226 at
//! 1.15.1): it inserts the replayed branching change into the multiset *before*
//! `changeBound`, so the `boundChangeAdded` that `changeBound` then makes for the same
//! change — with `branching == false`, since it is replayed with `Reason::unspecified()` —
//! finds a non-empty multiset and takes the silent path instead of the assert. The port
//! cannot see that through its bounds, so [`HighsDebugSol::begin_branching_replay`] holds
//! the entry for the duration of that one `change_bound` (the window is exactly the span in
//! which the C++ entry exists before the change has taken effect). Two divergences remain,
//! both in the silent direction and only reachable when the pre-recorded change never takes
//! effect: HiGHS keeps its entry until the change is popped by `backtrackToGlobal`, so if
//! `changeBound` returns without tightening a bound (a no-op or an infeasible replay) HiGHS
//! stays silent for the following changes of that domain while the port resumes reporting.
//! Neither changes what a solve does: everything here is inert unless
//! `mip_debug_solution_file` names a file.

use crate::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

static NUM_REPORTS: AtomicUsize = AtomicUsize::new(0);

/// HiGHS: HighsDebugSol.h::HighsDebugSol (lines 31-80)
#[derive(Clone, Debug)]
pub struct HighsDebugSol {
    /// C++: `double debugSolObjective`
    pub debug_sol_objective: f64,
    /// C++: `std::vector<double> debugOrigSolution`
    pub debug_orig_solution: Vec<f64>,
    /// C++: `std::vector<double> debugSolution`
    pub debug_solution: Vec<f64>,
    /// C++: `bool debugSolActive`
    pub debug_sol_active: bool,
    /// The `bool debugSolActive` local of HighsMipSolver::run / HighsMipSolverData::performRestart
    /// that the C++ swaps the flag into while presolve runs.
    pub saved_active: bool,
    /// `mipsolver->mipdata_->feastol`, captured when the solution is mapped (constant in a solve).
    pub feastol: f64,
    /// PORT NOTE: stands in for the one `conflictingBounds[&domain]` entry that HiGHS
    /// inserts before `changeBound` in `setDomainChangeStack` (see the module PORT NOTE).
    /// Set only for the span of that single `change_bound` call. A `Cell` rather than an
    /// atomic: `HighsDebugSol` is `Clone`, and the window is only ever opened and closed on
    /// the thread that is running the replay.
    branching_replay_window: std::cell::Cell<bool>,
}

impl Default for HighsDebugSol {
    fn default() -> Self {
        Self {
            debug_sol_objective: -kHighsInf,
            debug_orig_solution: Vec::new(),
            debug_solution: Vec::new(),
            debug_sol_active: false,
            saved_active: false,
            feastol: 1e-6,
            branching_replay_window: std::cell::Cell::new(false),
        }
    }
}

impl HighsDebugSol {
    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::HighsDebugSol (lines 20-25)
    /// C++: `HighsDebugSol(HighsMipSolver & )`
    pub fn new(_arg0: &mut HighsMipSolver) -> Self {
        Default::default()
    }

    /// Report a violated check (C++ `assert(false)`).
    fn report(&self, what: std::fmt::Arguments<'_>) {
        let n = NUM_REPORTS.fetch_add(1, Ordering::Relaxed);
        if std::env::var("RUSOLVE_DEBUG_SOL_ABORT").is_ok_and(|v| v == "1") {
            panic!("HighsDebugSol: {what}");
        }
        eprintln!("HighsDebugSol: {what}");
        if n < 2 {
            eprintln!("{}", std::backtrace::Backtrace::force_capture());
        }
    }

    /// True when the domain's current bounds contain the debug solution (the C++
    /// `conflictingBounds[&domain].empty()`, see the module doc).
    pub(crate) fn domain_contains(&self, domain: &HighsDomain) -> bool {
        let n = self.debug_solution.len();
        if domain.col_lower_.len() != n {
            return false;
        }
        (0..n).all(|i| {
            domain.col_lower_[i] <= self.debug_solution[i] + self.feastol
                && domain.col_upper_[i] >= self.debug_solution[i] - self.feastol
        })
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::activate (lines 27-89)
    /// C++: `void activate()`
    /// PORT NOTE: a solver-level function (`mipsolver->` reads); the `debug_mip.mps` dump of the
    /// could-not-open branch is reduced to its warning.
    pub fn activate(ms: &mut HighsMipSolver) {
        if ms.submip
            || ms.mipdata_.debug_solution.debug_sol_objective > ms.mipdata_.upper_limit
            || ms.options_mip_.mip_debug_solution_file.is_empty()
        {
            return;
        }
        let path = ms.options_mip_.mip_debug_solution_file.clone();
        let Ok(content) = std::fs::read_to_string(&path) else {
            eprintln!("debug solution: could not open file '{path}'");
            return;
        };
        let num_col = ms.orig_model_.num_col_ as usize;
        let mut nametoidx: HashMap<&str, usize> = HashMap::new();
        for i in 0..num_col {
            if let Some(name) = ms.orig_model_.col_names_.get(i) {
                nametoidx.insert(name.as_str(), i);
            }
        }
        // PORT NOTE: a model passed without column names (the CLI does not pass them) reads
        // the `# Columns` section by position (a HiGHS solution file lists columns in order).
        let by_position = nametoidx.len() != num_col;
        let mut position = 0usize;
        let mut sol = vec![0.0; num_col];
        let mut incolsection = false;
        for line in content.lines() {
            if line.contains("# Columns") {
                incolsection = true;
                continue;
            }
            if line.contains("# Rows") {
                break;
            }
            if !incolsection {
                continue;
            }
            let mut it = line.split_whitespace();
            if let Some(name) = it.next() {
                let idx = if by_position {
                    position += 1;
                    (position <= num_col).then_some(position - 1)
                } else {
                    nametoidx.get(name).copied()
                };
                if let Some(idx) = idx {
                    if let Some(v) = it.next().and_then(|v| v.parse::<f64>().ok()) {
                        sol[idx] = v;
                    }
                }
            }
        }
        let mut debugsolobj = HighsCDouble::from(0.0);
        for i in 0..num_col {
            debugsolobj += HighsCDouble::from(sol[i]) * ms.orig_model_.col_cost_[i];
        }
        let dbg = &mut ms.mipdata_.debug_solution;
        dbg.debug_sol_objective = f64::from(debugsolobj + ms.orig_model_.offset_);
        dbg.debug_orig_solution = sol.clone();
        dbg.debug_solution = sol;
        dbg.debug_sol_active = true;
        dbg.feastol = ms.mipdata_.feastol;
        eprintln!(
            "debug sol active (objective {})",
            ms.mipdata_.debug_solution.debug_sol_objective
        );
    }

    /// The `std::swap(debugSolution.debugSolActive, debugSolActive)` before presolve
    /// (HighsMipSolver.cpp:120, HighsMipSolverData.cpp:1342).
    pub fn suspend(&mut self) {
        self.saved_active = self.debug_sol_active;
        self.debug_sol_active = false;
    }

    /// The `debugSolution.debugSolActive = debugSolActive` after presolve
    /// (HighsMipSolver.cpp:155, HighsMipSolverData.cpp:1443).
    pub fn resume(&mut self) {
        self.debug_sol_active = self.saved_active;
    }

    /// HiGHS: HighsMipSolverData.cpp::runSetup `#ifdef HIGHS_DEBUGSOL` block (lines 1132-1145
    /// at 1.15.1, 1001-1014 at 1.16.0)
    pub fn run_setup_block(ms: &mut HighsMipSolver) {
        if !ms.mipdata_.debug_solution.debug_sol_active {
            return;
        }
        let red = ms
            .mipdata_
            .post_solve_stack
            .get_reduced_primal_solution(&ms.mipdata_.debug_solution.debug_orig_solution);
        let mut debugsolobj = HighsCDouble::from(0.0);
        for i in 0..ms.num_col() {
            debugsolobj += HighsCDouble::from(ms.col_cost()[i as usize]) * red[i as usize];
        }
        // Not `+ model_->offset_`: C++ stops at `static_cast<double>(debugsolobj)` because
        // `upper_bound`/`upper_limit`, which `newIncumbentFound` compares this with, are kept
        // without the offset too (`solobj = ... - model_->offset_`, HighsMipSolverData.cpp:889
        // at 1.16.0). pg5_34's -13137.35 is therefore the right number, not the -14339.35 model
        // value, and `debug sol inactive` once the optimum is found is the intended exit.
        let ok = ms.check_solution(&red);
        let dbg = &mut ms.mipdata_.debug_solution;
        dbg.debug_solution = red;
        dbg.debug_sol_objective = f64::from(debugsolobj);
        dbg.feastol = ms.mipdata_.feastol;
        eprintln!(
            "debug sol mapped to the presolved space: {} columns, objective {}",
            dbg.debug_solution.len(),
            dbg.debug_sol_objective
        );
        // C++ registerDomain(getDomain())
        ms.mipdata_
            .debug_solution
            .register_domain(ms.get_domain_2());
        if !ok {
            let dbg = &ms.mipdata_.debug_solution;
            let model = mip_model!(ms);
            let x = &dbg.debug_solution;
            let mut details = String::new();
            for i in 0..x.len() {
                if x[i] < model.col_lower_[i] - dbg.feastol
                    || x[i] > model.col_upper_[i] + dbg.feastol
                {
                    details += &format!(
                        " col {i}={} not in [{}, {}];",
                        x[i], model.col_lower_[i], model.col_upper_[i]
                    );
                }
            }
            let mut row_value = vec![HighsCDouble::from(0.0); model.num_row_ as usize];
            if model.a_matrix_.is_colwise() {
                for j in 0..model.num_col_ as usize {
                    for k in model.a_matrix_.start_[j]..model.a_matrix_.start_[j + 1] {
                        let k = k as usize;
                        row_value[model.a_matrix_.index_[k] as usize] +=
                            HighsCDouble::from(x[j]) * model.a_matrix_.value_[k];
                    }
                }
                for r in 0..model.num_row_ as usize {
                    let v = f64::from(row_value[r]);
                    if v < model.row_lower_[r] - dbg.feastol
                        || v > model.row_upper_[r] + dbg.feastol
                    {
                        details += &format!(
                            " row {r}={v} not in [{}, {}];",
                            model.row_lower_[r], model.row_upper_[r]
                        );
                    }
                }
            }
            dbg.report(format_args!(
                "checkSolution failed after presolve:{details}"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::newIncumbentFound (lines 124-129 at 1.15.1,
    /// 127-132 at 1.16.0)
    /// C++: `void newIncumbentFound()`
    pub fn new_incumbent_found(&mut self, upper_limit: f64) {
        if self.debug_sol_active && self.debug_sol_objective > upper_limit {
            eprintln!("debug sol inactive");
            self.debug_sol_active = false;
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::shrink (lines 94-108)
    /// C++: `void shrink(const std::vector<HighsInt> & )`
    /// PORT NOTE: `&self` as the call site; it runs only during a restart presolve, where the
    /// flag is suspended (as in C++), so the C++ body never executes.
    pub fn shrink(&self, _new_col_index: &[HighsInt]) {}

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::registerDomain (lines 105-115)
    /// C++: `void registerDomain(const HighsDomain & )`
    pub fn register_domain(&self, domain: &HighsDomain) {
        if !self.debug_sol_active {
            return;
        }
        if !self.domain_contains(domain) {
            self.report(format_args!(
                "registerDomain: the domain does not contain the debug solution"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::boundChangeAdded (lines 125-153)
    /// C++: `void boundChangeAdded(const HighsDomain & , const HighsDomainChange & , bool  = false)`
    /// PORT NOTE: `reason` is extra (reported only).
    pub fn bound_change_added(
        &self,
        domain: &HighsDomain,
        domchg: &HighsDomainChange,
        branching: bool,
        reason: Option<&HighsDomain_Reason>,
    ) {
        if !self.debug_sol_active {
            return;
        }
        let col = domchg.column as usize;
        if col >= self.debug_solution.len() || domain.col_lower_.len() != self.debug_solution.len()
        {
            return;
        }
        // C++ lines 133-145: keep only the changes that exclude the debug solution.
        if domchg.boundtype == HighsBoundType::kLower {
            if domchg.boundval <= self.debug_solution[col] + self.feastol {
                return;
            }
        } else if domchg.boundval >= self.debug_solution[col] - self.feastol {
            return;
        }
        // C++ lines 148-151 at 1.15.1, 151-154 at 1.16.0: `if (branching ||
        // !conflictingBounds[&domain].empty()) { conflictingBounds[&domain].insert(domchg);
        // return; }`, so only the *first* non-branching change that cuts the debug solution
        // off reaches the assert. The
        // multiset test is `!domain_contains` (module PORT NOTE), plus the window that
        // stands in for the entry inserted by `setDomainChangeStack` (see
        // [`HighsDebugSol::begin_branching_replay`]).
        if branching || self.branching_replay_window.get() || !self.domain_contains(domain) {
            return;
        }
        self.report(format_args!(
            "boundChangeAdded cuts off the debug solution: col {col} {:?} {} (debug value {}), reason {:?}, depth {}",
            domchg.boundtype,
            domchg.boundval,
            self.debug_solution[col],
            reason,
            domain.branch_pos_.len()
        ));
    }

    /// PORT NOTE: opens the window in which C++'s `conflictingBounds[&domain]` holds the
    /// entry that `HighsDomain::setDomainChangeStack` inserted for the change it is about to
    /// replay (HighsDomain.cpp:2467-2468 at 1.16.0, 2225-2226 at 1.15.1); see the module
    /// PORT NOTE. Call [`HighsDebugSol::end_branching_replay`] right after the
    /// `change_bound` that the C++ call precedes.
    pub fn begin_branching_replay(&self) {
        if !self.debug_sol_active {
            return;
        }
        self.branching_replay_window.set(true);
    }

    /// PORT NOTE: closes [`HighsDebugSol::begin_branching_replay`]. From here on the domain
    /// bounds carry the condition, as the replayed change has excluded the debug solution.
    pub fn end_branching_replay(&self) {
        self.branching_replay_window.set(false);
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::boundChangeRemoved (lines 155-164)
    /// PORT NOTE: nothing to do without the multiset (module doc).
    pub fn bound_change_removed(&self, _domain: &HighsDomain, _domchg: &HighsDomainChange) {}

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::resetDomain (lines 211-215)
    /// PORT NOTE: nothing to do without the multiset (module doc).
    pub fn reset_domain(&self, _domain: &HighsDomain) {}

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::nodePruned (lines 217-223)
    /// C++: `void nodePruned(const HighsDomain & )`
    pub fn node_pruned(&self, localdomain: &HighsDomain) {
        if !self.debug_sol_active {
            return;
        }
        if self.domain_contains(localdomain) {
            self.report(format_args!(
                "nodePruned: a node containing the debug solution was pruned (depth {}, infeasible {})",
                localdomain.branch_pos_.len(),
                localdomain.infeasible_
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkCut (lines 166-176)
    /// C++: `void checkCut(const HighsInt * , const double * , HighsInt , double )`
    pub fn check_cut(&self, rindex: &[HighsInt], rvalue: &[f64], rlen: HighsInt, rhs: f64) {
        if !self.debug_sol_active {
            return;
        }
        let mut violation = HighsCDouble::from(-rhs);
        for i in 0..rlen as usize {
            violation += HighsCDouble::from(self.debug_solution[rindex[i] as usize]) * rvalue[i];
        }
        let v = f64::from(violation);
        if v > self.feastol {
            self.report(format_args!(
                "checkCut: cut of length {rlen} rhs {rhs} cuts off the debug solution (violation {v})"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkRow (lines 199-209)
    /// C++: `void checkRow(const HighsInt * , const double * , HighsInt , double , double )`
    pub fn check_row(
        &self,
        rindex: &[HighsInt],
        rvalue: &[f64],
        rlen: HighsInt,
        rlower: f64,
        rupper: f64,
    ) {
        if !self.debug_sol_active {
            return;
        }
        let mut activity = HighsCDouble::from(0.0);
        for i in 0..rlen as usize {
            activity += HighsCDouble::from(self.debug_solution[rindex[i] as usize]) * rvalue[i];
        }
        let a = f64::from(activity);
        if a - self.feastol > rupper || a + self.feastol < rlower {
            self.report(format_args!(
                "checkRow: activity {a} not in [{rlower}, {rupper}]"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkRowAggregation (lines 178-197)
    /// C++: `void checkRowAggregation(const HighsLp & , const HighsInt * , const double * , HighsInt )`
    pub fn check_row_aggregation(
        &self,
        lp: &HighsLp,
        rindex: &[HighsInt],
        rvalue: &[f64],
        rlen: HighsInt,
    ) {
        if !self.debug_sol_active || !lp.a_matrix_.is_colwise() {
            return;
        }
        let num_col = lp.num_col_ as usize;
        // calculateRowValuesQuad
        let mut row_value = vec![HighsCDouble::from(0.0); lp.num_row_ as usize];
        for j in 0..num_col {
            for k in lp.a_matrix_.start_[j]..lp.a_matrix_.start_[j + 1] {
                let k = k as usize;
                row_value[lp.a_matrix_.index_[k] as usize] +=
                    HighsCDouble::from(self.debug_solution[j]) * lp.a_matrix_.value_[k];
            }
        }
        let mut violation = HighsCDouble::from(0.0);
        for i in 0..rlen as usize {
            let r = rindex[i] as usize;
            if r < num_col {
                violation += HighsCDouble::from(self.debug_solution[r]) * rvalue[i];
            } else {
                violation += HighsCDouble::from(f64::from(row_value[r - num_col])) * rvalue[i];
            }
        }
        let viol = f64::from(violation).abs();
        if viol > self.feastol {
            self.report(format_args!(
                "checkRowAggregation: aggregation of length {rlen} violated by {viol}"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkClique (lines 225-234)
    /// C++: `void checkClique(const HighsCliqueTable::CliqueVar * , HighsInt )`
    // PORT NOTE: signature changed because C++ takes a `CliqueVar*` plus length; a slice carries both.
    pub fn check_clique(&self, clq: &[HighsCliqueTable_CliqueVar]) {
        if !self.debug_sol_active {
            return;
        }
        let mut violation: HighsInt = -1;
        for v in clq {
            violation += (v.weight(&self.debug_solution) + 0.5) as HighsInt;
        }
        if violation > 0 {
            self.report(format_args!(
                "checkClique: clique of length {} cuts off the debug solution",
                clq.len()
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkVub (lines 236-242)
    pub fn check_vub(&self, col: HighsInt, vubcol: HighsInt, vubcoef: f64, vubconstant: f64) {
        if !self.debug_sol_active || vubcoef.abs() == kHighsInf {
            return;
        }
        let (ci, vi) = (col as usize, vubcol as usize);
        if self.debug_solution[ci] > self.debug_solution[vi] * vubcoef + vubconstant + self.feastol
        {
            self.report(format_args!(
                "checkVub: col {col} <= {vubcoef} col {vubcol} + {vubconstant} cuts off the debug solution"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkVlb (lines 244-250)
    pub fn check_vlb(&self, col: HighsInt, vlbcol: HighsInt, vlbcoef: f64, vlbconstant: f64) {
        if !self.debug_sol_active || vlbcoef.abs() == kHighsInf {
            return;
        }
        let (ci, vi) = (col as usize, vlbcol as usize);
        if self.debug_solution[ci] < self.debug_solution[vi] * vlbcoef + vlbconstant - self.feastol
        {
            self.report(format_args!(
                "checkVlb: col {col} >= {vlbcoef} col {vlbcol} + {vlbconstant} cuts off the debug solution"
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkConflictReasonFrontier (lines 252-275)
    pub fn check_conflict_reason_frontier(
        &self,
        reason_side_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        _domchgstack: &[HighsDomainChange],
    ) {
        if !self.debug_sol_active {
            return;
        }
        let mut num_active: usize = 0;
        for d in reason_side_frontier {
            let col = d.domchg.column as usize;
            if col >= self.debug_solution.len() {
                return;
            }
            if d.domchg.boundtype == HighsBoundType::kLower {
                if self.debug_solution[col] >= d.domchg.boundval - self.feastol {
                    num_active += 1;
                }
            } else if self.debug_solution[col] <= d.domchg.boundval + self.feastol {
                num_active += 1;
            }
        }
        if num_active >= reason_side_frontier.len() {
            self.report(format_args!(
                "checkConflictReasonFrontier: conflict of length {} cuts off the debug solution",
                reason_side_frontier.len()
            ));
        }
    }

    /// HiGHS: HighsDebugSol.cpp::HighsDebugSol::checkConflictReconvergenceFrontier (lines 277-311)
    pub fn check_conflict_reconvergence_frontier(
        &self,
        reconvergence_frontier: &BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
        reconv_domchg: &HighsDomain_ConflictSet_LocalDomChg,
        _domchgstack: &[HighsDomainChange],
    ) {
        if !self.debug_sol_active {
            return;
        }
        let mut num_active: usize = 0;
        let x = &self.debug_solution;
        for d in reconvergence_frontier {
            let col = d.domchg.column as usize;
            if col >= x.len() {
                return;
            }
            if d.domchg.boundtype == HighsBoundType::kLower {
                if x[col] >= d.domchg.boundval {
                    num_active += 1;
                }
            } else if x[col] <= d.domchg.boundval {
                num_active += 1;
            }
        }
        // C++ `getDomain().flip(reconvDomchg.domchg)`
        let rc = &reconv_domchg.domchg;
        let col = rc.column as usize;
        if col >= x.len() {
            return;
        }
        // PORT NOTE: HighsDomain::flip needs the solver (`is_col_integral`); the port flips as
        // for an integral column, which over-reports only on a continuous reconvergence column.
        if rc.boundtype == HighsBoundType::kLower {
            if x[col] <= (rc.boundval - self.feastol).floor() {
                num_active += 1;
            }
        } else if x[col] >= (rc.boundval + self.feastol).ceil() {
            num_active += 1;
        }
        if num_active > reconvergence_frontier.len() {
            self.report(format_args!(
                "checkConflictReconvergenceFrontier: reconvergence cut of length {} cuts off the debug solution",
                reconvergence_frontier.len()
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    //! Covers the offset-free objective of the
    //! `runSetup` mapping (faithful, not a bug), and the `branching == true` suppression of
    //! `boundChangeAdded`.

    use super::*;
    use std::sync::Mutex;

    /// `NUM_REPORTS` is process wide and `cargo test` runs tests on many threads, so the
    /// tests that read it hold this lock.
    static REPORT_LOCK: Mutex<()> = Mutex::new(());

    /// Test hook: the number of violated checks reported so far in this process.
    fn reports() -> usize {
        NUM_REPORTS.load(Ordering::Relaxed)
    }

    fn lower(column: HighsInt, boundval: f64) -> HighsDomainChange {
        HighsDomainChange {
            boundval,
            column,
            boundtype: HighsBoundType::kLower,
        }
    }

    /// A two-column integer model with no rows, so that `change_bound`, `check_solution` and
    /// the postsolve stack all find a model to read with no presolve or LP state.
    fn two_col_model(offset: f64) -> HighsLp {
        let mut lp = HighsLp::default();
        lp.num_col_ = 2;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![3.0, 5.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.offset_ = offset;
        lp.integrality_ = vec![HighsVarType::kInteger, HighsVarType::kInteger];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0, 0, 0];
        lp
    }

    /// A solver in the state `activate` leaves behind before `runSetup`: the debug solution
    /// is loaded and active, and the global domain holds it.
    fn solver_with_debug_sol(debug_sol: &[f64]) -> HighsMipSolver {
        let mut options = HighsOptions::default();
        options.output_flag = false;
        options.log_options.output_flag = false;
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms =
            HighsMipSolver::new(&mut cb, &options, &two_col_model(0.0), &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.setup_domain_propagation();
        ms.mipdata_
            .post_solve_stack
            .initialize_index_maps(0, ms.num_col());
        let feastol = ms.mipdata_.feastol;
        let dbg = &mut ms.mipdata_.debug_solution;
        dbg.debug_orig_solution = debug_sol.to_vec();
        dbg.debug_solution = debug_sol.to_vec();
        dbg.debug_sol_active = true;
        dbg.feastol = feastol;
        ms
    }

    /// `runSetup` maps the debug solution into the presolved space and, as in C++, reports
    /// its objective without the offset: the space `upper_limit` lives in.
    #[test]
    fn run_setup_block_objective_excludes_the_model_offset() {
        // cost [3, 5] at x [2, 3] is 21, whatever the offset.
        let mut ms = solver_with_debug_sol(&[2.0, 3.0]);
        ms.orig_model_.offset_ = -1202.0;
        HighsDebugSol::run_setup_block(&mut ms);
        assert_eq!(ms.mipdata_.debug_solution.debug_sol_objective, 21.0);

        let mut ms = solver_with_debug_sol(&[2.0, 3.0]);
        let mut presolved = two_col_model(-100.0);
        presolved.col_cost_ = vec![1.0, 1.0];
        ms.mipdata_.presolved_model = presolved;
        ms.model_ = HighsMipModel::Presolved;
        HighsDebugSol::run_setup_block(&mut ms);
        assert_eq!(ms.mipdata_.debug_solution.debug_sol_objective, 5.0);
    }

    #[test]
    fn branching_bound_changes_are_recorded_not_reported() {
        let _serialise = REPORT_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let base = reports();
        let mut ms = solver_with_debug_sol(&[5.0, 5.0]);

        detach!(ms.mipdata_.domains[0], |dom| dom.change_bound(
            &mut ms,
            lower(0, 7.0),
            HighsDomain_Reason::branching()
        ));
        assert_eq!(
            reports(),
            base,
            "a branching change that excludes the debug solution must not be reported"
        );

        // the propagated change under it is silent too: in C++ because conflictingBounds is
        // non-empty, in the port because the domain bounds now exclude the debug solution
        detach!(ms.mipdata_.domains[0], |dom| dom.change_bound(
            &mut ms,
            lower(1, 6.0),
            HighsDomain_Reason::unspecified()
        ));
        assert_eq!(
            reports(),
            base,
            "a change propagated under a branching change must not be reported"
        );
    }

    /// The ported call: `HighsDomain::setDomainChangeStack` pre-records the trailing change
    /// it replays (HighsDomain.cpp:2467-2468 at 1.16.0, 2225-2226 at 1.15.1), so the
    /// `boundChangeAdded` that `changeBound` makes for that same change stays quiet even
    /// though it is the first change to cut off the debug solution.
    #[test]
    fn replayed_node_changes_are_recorded_not_reported() {
        let _serialise = REPORT_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let base = reports();
        let mut ms = solver_with_debug_sol(&[5.0, 5.0]);
        // col 0 branches to 3, which still holds the debug solution at 5; the trailing
        // change on col 1 is the one that cuts it off.
        let stack = vec![lower(0, 3.0), lower(1, 6.0)];
        detach!(ms.mipdata_.domains[0], |dom| dom.set_domain_change_stack_2(
            &mut ms,
            &stack,
            &[0]
        ));
        assert_eq!(
            reports(),
            base,
            "a replayed node bound change that cuts off the debug solution must not be reported"
        );
        assert!(
            !ms.mipdata_.debug_solution.branching_replay_window.get(),
            "the suppression window must not outlive the change it was opened for"
        );
    }

    /// The suppression must not swallow the real signal: the first non-branching change at
    /// a node that still holds the debug solution reaches C++'s assert, so it is reported.
    #[test]
    fn a_non_branching_change_that_cuts_off_the_debug_solution_is_reported() {
        let _serialise = REPORT_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let base = reports();
        let mut ms = solver_with_debug_sol(&[5.0, 5.0]);
        detach!(ms.mipdata_.domains[0], |dom| dom.change_bound(
            &mut ms,
            lower(0, 7.0),
            HighsDomain_Reason::unspecified()
        ));
        assert_eq!(
            reports(),
            base + 1,
            "the first non-branching change that cuts off the debug solution must be reported"
        );
    }
}
