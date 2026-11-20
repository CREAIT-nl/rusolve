//! Translated from SCIP (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! `mip_shift_and_propagate`: the Shift-and-Propagate primal heuristic, translated from
//! SCIP `src/scip/heur_shiftandpropagate.c` (scipopt/scip commit 63417ae) and `sorttpl.c`.
//! Berthold and Hendel, "Shift-and-propagate", J. Heuristics 21 (2015) 73-106.
//!
//! LP-free: the model is transformed so every discrete variable starts at its lower bound
//! (zero in the transformed space), then the variables are visited in a fixed order; each is
//! shifted to the value that minimises the weighted number of violated rows, fixed, and the
//! fixing is propagated (`HighsDomain`, standing in for `SCIPpropagateProbing`). A cutoff
//! backtracks the fixing, tightens the bound by one, and postpones the variable. Continuous
//! columns are relaxed out of the rows and are left to the LP of `try_rounded_point`
//! (SCIP: `SCIPsolveProbingLP`).
//!
//! Defaults of SCIP's parameters are constants here. Not ported, because they are off by
//! default: `binlocksfirst`, `selectbest`, `nozerofixing`, `updateweights`, the random sort
//! key `r`. Differences from SCIP: `HighsDomain::propagate` runs to a fixpoint (SCIP stops
//! after `nproprounds` = 10 rounds); the probing depth is the count of branching positions of
//! the local domain; time is checked every 10 variables against the MIP time limit.

use super::sort::sorttpl_sort;
use super::*;
use crate::mip::HighsMipSolver;
use crate::mip::mip_solver_data::kSolutionSourceShiftProp;

/// heur_shiftandpropagate.c (lines 76-77)
const WEIGHT_INEQUALITY: i32 = 1;
const WEIGHT_EQUALITY: i32 = 3;
/// heur_shiftandpropagate.c (line 83): fixed maximum number of propagations
const PROPBREAKER: i64 = 65000;
/// heur_shiftandpropagate.c (line 84): allowed cutoffs before the heuristic stops
const CUTOFFBREAKER: i32 = 15;
/// heur_shiftandpropagate.c (line 88-89): sort key `v`, sort variables
const SORTKEY: char = 'v';
const SORTVARS: bool = true;
/// heur_shiftandpropagate.c (line 91-92)
const PREFERBINARIES: bool = true;
/// heur_shiftandpropagate.c (line 94): `maxcutoffquot`
const MAXCUTOFFQUOT: f64 = 0.0;
/// heur_shiftandpropagate.c (line 99): `fixbinlocks`
const FIXBINLOCKS: bool = true;
/// heur_shiftandpropagate.c (line 101): `normalize`
const NORMALIZE: bool = true;
/// SCIP `numerics/boundstreps` default (set.h SCIP_DEFAULT_NUM_BOUNDSTREPS): minimum relative
/// range improvement for a continuous bound change in propagation.
const SCIP_BOUNDSTREPS: f64 = 0.05;

/// Check the time limit every this many variables (SCIP: `c % 10 == 0 && SCIPisStopped`).
const STOP_CHECK_EVERY: usize = 10;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Transform {
    None,
    /// shifted by the lower bound: x - lb
    Lb,
    /// negated by the upper bound: ub - x
    Neg,
    /// free, no shift
    Free,
}

#[inline]
fn is_inf(v: f64) -> bool {
    v >= S6_INFINITY
}

/// Bounds of the port are +-inf, SCIP's are +-1e20.
#[inline]
fn sb(v: f64) -> f64 {
    v.clamp(-S6_INFINITY, S6_INFINITY)
}

/// SCIPisLE (absolute epsilon)
#[inline]
fn is_le(a: f64, b: f64) -> bool {
    a - b <= S6_EPSILON
}

/// SCIPisGT (absolute epsilon)
#[inline]
fn is_gt(a: f64, b: f64) -> bool {
    a - b > S6_EPSILON
}

/// SCIPisFeasEQ
#[inline]
fn feas_eq(a: f64, b: f64, ft: f64) -> bool {
    let d = s6_rel_diff(a, b);
    d <= ft && d >= -ft
}

/// Information about the matrix after its heuristic transformation (`CONSTRAINTMATRIX`).
struct Matrix {
    /// discrete columns by row, `nrows + 1` entries
    rowbegin: Vec<usize>,
    rowind: Vec<usize>,
    colbegin: Vec<usize>,
    colvals: Vec<f64>,
    colind: Vec<usize>,
    violrows: Vec<i32>,
    transform: Vec<Transform>,
    lhs: Vec<f64>,
    rhs: Vec<f64>,
    colnorms: Vec<f64>,
    upperbounds: Vec<f64>,
    shiftvals: Vec<f64>,
    nrows: usize,
    ndisc: usize,
    feastol: f64,
}

/// Violated rows, their positions, and the row weights.
struct Viol {
    rows: Vec<usize>,
    pos: Vec<i64>,
    weights: Vec<i32>,
}

impl Matrix {
    fn col(&self, cj: usize) -> (&[f64], &[usize]) {
        let (s, e) = (self.colbegin[cj], self.colbegin[cj + 1]);
        (&self.colvals[s..e], &self.colind[s..e])
    }

    /// heur_shiftandpropagate.c::transformVariable
    fn transform_variable(&mut self, colpos: usize, lb: f64, ub: f64) {
        let ft = self.feastol;
        let mut negatecoeffs = false;
        let deltashift;
        if is_inf(-lb) && is_inf(ub) {
            if self.transform[colpos] == Transform::Neg {
                negatecoeffs = true;
            }
            deltashift = self.shiftvals[colpos];
            self.shiftvals[colpos] = 0.0;
            self.transform[colpos] = Transform::Free;
        } else if is_le(lb.abs(), ub.abs()) {
            self.transform[colpos] = Transform::Lb;
            deltashift = lb;
            self.shiftvals[colpos] = lb;
        } else {
            if self.transform[colpos] != Transform::Neg {
                negatecoeffs = true;
            }
            self.transform[colpos] = Transform::Neg;
            deltashift = ub;
            self.shiftvals[colpos] = ub;
        }

        // upper bound in the transformation (the lower bound is implicitly zero)
        if !is_inf(ub) && !is_inf(lb) {
            self.upperbounds[colpos] = (ub - lb).min(S6_INFINITY);
        } else {
            self.upperbounds[colpos] = S6_INFINITY;
        }

        if !s6_feas_zero(deltashift, ft) || negatecoeffs {
            let (s, e) = (self.colbegin[colpos], self.colbegin[colpos + 1]);
            for i in s..e {
                let rowpos = self.colind[i];
                let v = self.colvals[i];
                if !is_inf(-self.lhs[rowpos]) {
                    self.lhs[rowpos] -= v * deltashift;
                }
                if !is_inf(self.rhs[rowpos]) {
                    self.rhs[rowpos] -= v * deltashift;
                }
                if negatecoeffs {
                    self.colvals[i] = -v;
                }
            }
        }
    }

    /// heur_shiftandpropagate.c::checkRowViolation (updateweights is off)
    fn check_row_violation(&mut self, viol: &mut Viol, rowindex: usize) {
        let ft = self.feastol;
        let mut violadd = 0i32;
        if viol.pos[rowindex] == -1
            && (s6_feas_gt(self.lhs[rowindex], 0.0, ft) || s6_feas_lt(self.rhs[rowindex], 0.0, ft))
        {
            // row is now violated: enqueue it
            viol.pos[rowindex] = viol.rows.len() as i64;
            viol.rows.push(rowindex);
            violadd = 1;
        } else if viol.pos[rowindex] >= 0
            && s6_feas_le(self.lhs[rowindex], 0.0, ft)
            && s6_feas_ge(self.rhs[rowindex], 0.0, ft)
        {
            // row is now feasible: swap it with the last violated row and drop it
            let n = viol.rows.len();
            let p = viol.pos[rowindex] as usize;
            if p != n - 1 {
                let last = viol.rows[n - 1];
                viol.rows[p] = last;
                viol.pos[last] = p as i64;
            }
            viol.pos[rowindex] = -1;
            viol.rows.pop();
            violadd = -1;
        }
        if violadd != 0 {
            let (s, e) = (self.rowbegin[rowindex], self.rowbegin[rowindex + 1]);
            for k in s..e {
                self.violrows[self.rowind[k]] += violadd;
            }
        }
    }

    /// heur_shiftandpropagate.c::checkViolations; `colidx` is `None` for all rows.
    fn check_violations(&mut self, viol: &mut Viol, colidx: Option<usize>) {
        match colidx {
            Some(cj) => {
                let (s, e) = (self.colbegin[cj], self.colbegin[cj + 1]);
                for i in s..e {
                    let row = self.colind[i];
                    self.check_row_violation(viol, row);
                }
            }
            None => {
                viol.rows.clear();
                for p in viol.pos.iter_mut() {
                    *p = -1;
                }
                for v in self.violrows.iter_mut() {
                    *v = 0;
                }
                for row in 0..self.nrows {
                    self.check_row_violation(viol, row);
                }
            }
        }
    }

    /// heur_shiftandpropagate.c::retransformVariable
    fn retransform(&self, varindex: usize, solvalue: f64) -> f64 {
        match self.transform[varindex] {
            Transform::Lb => solvalue + self.shiftvals[varindex],
            Transform::Neg => self.shiftvals[varindex] - solvalue,
            _ => solvalue,
        }
    }

    /// heur_shiftandpropagate.c::updateTransformation
    fn update_transformation(&mut self, viol: &mut Viol, varindex: usize, lb: f64, ub: f64) {
        let ft = self.feastol;
        let mut deltashift = 0.0;
        let mut checkviolations = false;
        match self.transform[varindex] {
            Transform::Lb => {
                if is_inf(-lb) {
                    self.transform_variable(varindex, lb, ub);
                    checkviolations = true;
                } else {
                    deltashift = lb - self.shiftvals[varindex];
                    self.shiftvals[varindex] = lb;
                    self.upperbounds[varindex] = if !is_inf(ub) { ub - lb } else { S6_INFINITY };
                }
            }
            Transform::Neg => {
                if is_inf(ub) {
                    self.transform_variable(varindex, lb, ub);
                    checkviolations = true;
                } else {
                    deltashift = self.shiftvals[varindex] - ub;
                    self.shiftvals[varindex] = ub;
                    self.upperbounds[varindex] = if !is_inf(-lb) {
                        (ub - lb).min(S6_INFINITY)
                    } else {
                        S6_INFINITY
                    };
                }
            }
            Transform::Free => {
                // a bound became finite: transform to a lower-bound or negated variable
                if !is_inf(-lb) || !is_inf(ub) {
                    self.transform_variable(varindex, lb, ub);
                    checkviolations = true;
                }
            }
            Transform::None => {}
        }
        if !s6_feas_zero(deltashift, ft) {
            let (s, e) = (self.colbegin[varindex], self.colbegin[varindex + 1]);
            for i in s..e {
                let r = self.colind[i];
                let v = self.colvals[i];
                if !is_inf(-self.lhs[r]) {
                    self.lhs[r] -= v * deltashift;
                }
                if !is_inf(self.rhs[r]) {
                    self.rhs[r] -= v * deltashift;
                }
            }
            checkviolations = true;
        }
        if checkviolations {
            self.check_violations(viol, Some(varindex));
        }
    }

    /// heur_shiftandpropagate.c::getOptimalShiftingValue. Returns `(beststep, rowviolations)`;
    /// `rowviolations` is 0 when the variable cannot affect any row (SCIP leaves the caller's 0).
    fn optimal_shifting_value(
        &self,
        varindex: usize,
        direction: i32,
        weights: &[i32],
        steps: &mut Vec<f64>,
        violationchange: &mut Vec<i32>,
    ) -> (f64, i32) {
        let ft = self.feastol;
        let upperbound = self.upperbounds[varindex];
        let (vals, rows) = self.col(varindex);
        let nrows = rows.len();
        steps.clear();
        violationchange.clear();
        steps.resize(nrows, 0.0);
        violationchange.resize(nrows, 0);

        let mut allzero = true;
        let mut slacksurplus = 0.0;
        for i in 0..nrows {
            let rowpos = rows[i];
            let lhs = self.lhs[rowpos];
            let rhs = self.rhs[rowpos];
            let rowweight = weights[rowpos];
            let val = direction as f64 * vals[i];

            let rowisviolated = s6_feas_lt(rhs, 0.0, ft) || s6_feas_lt(-lhs, 0.0, ft);

            if !rowisviolated {
                // minimum integer shift that violates this feasible row
                let mut maxfeasshift = S6_INFINITY;
                if val < 0.0 {
                    if !is_inf(-lhs) {
                        maxfeasshift = s6_feas_floor(lhs.min(0.0) / val, ft);
                    } else {
                        slacksurplus -= val;
                    }
                } else if !is_inf(rhs) {
                    maxfeasshift = s6_feas_floor(rhs.max(0.0) / val, ft);
                } else {
                    slacksurplus += val;
                }

                if !is_inf(maxfeasshift) && s6_feas_le(maxfeasshift + 1.0, upperbound, ft) {
                    steps[i] = maxfeasshift + 1.0;
                    violationchange[i] = rowweight;
                    allzero = false;
                } else {
                    steps[i] = upperbound;
                    violationchange[i] = 0;
                }
            } else {
                // minimum integer shift that makes this violated row feasible
                let mut minfeasshift = S6_INFINITY;
                if val < 0.0 {
                    if s6_feas_lt(rhs, 0.0, ft) {
                        minfeasshift = s6_feas_ceil(rhs.min(val) / val, ft);
                    }
                } else if s6_feas_lt(-lhs, 0.0, ft) {
                    minfeasshift = s6_feas_ceil(lhs.max(val) / val, ft);
                }

                if !is_inf(minfeasshift) && s6_feas_le(minfeasshift, upperbound, ft) {
                    steps[i] = minfeasshift;
                    violationchange[i] = -rowweight;
                    allzero = false;
                } else {
                    steps[i] = upperbound;
                    violationchange[i] = 0;
                }
            }
        }

        // the variable cannot affect any row, but may add slack to feasible rows
        if allzero {
            if !is_inf(upperbound) && is_gt(slacksurplus, 0.0) {
                return (direction as f64 * upperbound, 0);
            }
            return (0.0, 0);
        }

        sorttpl_sort(steps, violationchange, &|x: f64, y: f64| x - y < 0.0);

        let mut beststep = 0.0;
        let mut rowviolations = 0i32;
        let mut sum = 0i32;
        let mut i = 0;
        while i < nrows && !is_inf(steps[i]) {
            sum += violationchange[i];
            if (i == nrows - 1 || steps[i + 1] > steps[i]) && sum < rowviolations {
                rowviolations = sum;
                beststep = direction as f64 * steps[i];
            }
            i += 1;
        }
        (beststep, rowviolations)
    }
}

/// Counters of one run (reporting only).
#[derive(Default)]
struct RunStats {
    probings: i64,
    cutoffs: i64,
    violations_left: i64,
}

/// Where the local domain stands, as the heuristic's probing path.
struct Probing {
    /// probing depth (SCIPgetProbingDepth): the root probing node plus one per node
    depth: i32,
}

/// Domain bounds of a discrete column in SCIP's +-1e20 convention.
#[inline]
fn bounds(dom: &HighsDomain, col: usize) -> (f64, f64) {
    (sb(dom.col_lower_[col]), sb(dom.col_upper_[col]))
}

/// Propagates the changes the local domain made since stack position `from` into the matrix
/// (SCIP: the bound-change event handler `eventExecShiftandpropagate`).
fn sync_changes_cols(
    m: &mut Matrix,
    viol: &mut Viol,
    dom: &HighsDomain,
    colpos_of: &[usize],
    disc_cols: &[usize],
    from: usize,
) {
    let stack = dom.get_domain_change_stack();
    if from >= stack.len() {
        return;
    }
    let mut touched: Vec<usize> = stack[from..]
        .iter()
        .map(|chg| colpos_of[chg.column as usize])
        .filter(|&p| p != usize::MAX)
        .collect();
    touched.sort_unstable();
    touched.dedup();
    for p in touched {
        let (lb, ub) = bounds(dom, disc_cols[p]);
        m.update_transformation(viol, p, lb, ub);
    }
}

/// Undo the probing nodes above `branching_len` branching positions (SCIPbacktrackProbing),
/// keeping the matrix in step with the restored bounds.
fn backtrack_to(
    ms: &mut HighsMipSolver,
    dom: &mut HighsDomain,
    m: &mut Matrix,
    viol: &mut Viol,
    colpos_of: &[usize],
    disc_cols: &[usize],
    branching_len: usize,
) {
    while dom.get_branching_positions().len() > branching_len {
        let from = *dom.get_branching_positions().last().unwrap() as usize;
        let touched: Vec<usize> = dom.get_domain_change_stack()[from..]
            .iter()
            .map(|chg| colpos_of[chg.column as usize])
            .filter(|&p| p != usize::MAX)
            .collect();
        dom.backtrack(ms);
        let mut touched = touched;
        touched.sort_unstable();
        touched.dedup();
        for p in touched {
            let (lb, ub) = bounds(dom, disc_cols[p]);
            m.update_transformation(viol, p, lb, ub);
        }
    }
}

/// heur_shiftandpropagate.c::heurExecShiftandpropagate. Returns whether a new incumbent was
/// found. The caller guarantees the heuristic is wanted (option on, not a sub-MIP).
pub(crate) fn shift_and_propagate(ms: &mut HighsMipSolver) -> bool {
    let t0 = ms.timer_.read();
    ms.stats_.shiftprop_calls += 1;
    let mut stats = RunStats::default();
    let found = shift_and_propagate_body(ms, &mut stats);
    ms.stats_.shiftprop_found += found as i64;
    ms.stats_.shiftprop_probings += stats.probings;
    ms.stats_.shiftprop_cutoffs += stats.cutoffs;
    ms.stats_.shiftprop_violations_left += stats.violations_left;
    ms.stats_.shiftprop_s += ms.timer_.read() - t0;
    found
}

fn shift_and_propagate_body(ms: &mut HighsMipSolver, stats: &mut RunStats) -> bool {
    let ncol = ms.num_col().max(0) as usize;
    let nrows = ms.num_row().max(0) as usize;
    let ft = ms.mipdata_.feastol;

    // heuristic is obsolete without discrete variables; "onlywithoutsol"; "we need an LP"
    if ms.mipdata_.integer_cols.is_empty() || ms.mipdata_.upper_bound < kHighsInf {
        return false;
    }
    if !ms.mipdata_.incumbent.is_empty() || nrows == 0 || ncol == 0 {
        return false;
    }

    // local copy of the global domain: the probing tree of the heuristic
    let gd = ms.mipdata_.workers[0].get_global_domain();
    let mut dom = detach!(ms.mipdata_.domains[gd], |g| HighsDomain::copy_from(ms, g));
    // SCIP accepts a continuous bound tightening that shrinks the range by `boundstreps`;
    // HiGHS wants 0.3, which hides the pairwise-disjunction infeasibility on fhnw.
    dom.cont_min_improve = SCIP_BOUNDSTREPS;
    let found = run(ms, &mut dom, stats, ncol, nrows, ft);
    dom.release(ms);
    found
}

#[allow(clippy::too_many_arguments)]
fn run(
    ms: &mut HighsMipSolver,
    dom: &mut HighsDomain,
    stats: &mut RunStats,
    ncol: usize,
    nrows: usize,
    ft: f64,
) -> bool {
    // ---- columns sorted by type: binary, integer, implied integer, continuous ----
    // (SCIPsortPtr with heurSortColsShiftandpropagate)
    let mut is_disc = vec![false; ncol];
    let mut is_bin = vec![false; ncol];
    for &cc in &ms.mipdata_.integer_cols {
        let cu = cc as usize;
        is_disc[cu] = true;
        is_bin[cu] = dom.col_lower_[cu] == 0.0 && dom.col_upper_[cu] == 1.0;
    }
    let mut type_key: Vec<i32> = (0..ncol)
        .map(|cc| {
            if is_bin[cc] {
                0
            } else if is_disc[cc] {
                1
            } else if ms.variable_type_2(cc as HighsInt) == HighsVarType::kImplicitInteger {
                2
            } else {
                3
            }
        })
        .collect();
    let mut lpcols: Vec<usize> = (0..ncol).collect();
    sorttpl_sort(&mut type_key, &mut lpcols, &|x: i32, y: i32| x - y < 0);
    let ndisc = ms.mipdata_.integer_cols.len();
    let nbin = type_key.iter().filter(|&&k| k == 0).count();
    let disc_cols: Vec<usize> = lpcols[..ndisc].to_vec();
    let mut colpos_of = vec![usize::MAX; ncol];
    for (p, &cc) in disc_cols.iter().enumerate() {
        debug_assert!(is_disc[cc]);
        colpos_of[cc] = p;
    }

    // ---- initMatrix ----
    let (a_rstart, a_rindex, a_rvalue) = (
        &ms.mipdata_.a_rstart_,
        &ms.mipdata_.a_rindex_,
        &ms.mipdata_.a_rvalue_,
    );
    let mut colcount = vec![0usize; ndisc];
    let mut nnonzs = 0usize;
    for k in 0..a_rindex.len() {
        let p = colpos_of[a_rindex[k] as usize];
        if p != usize::MAX {
            colcount[p] += 1;
            nnonzs += 1;
        }
    }
    if nnonzs == 0 {
        return false;
    }
    let mut m = Matrix {
        rowbegin: Vec::with_capacity(nrows + 1),
        rowind: Vec::with_capacity(nnonzs),
        colbegin: vec![0; ndisc + 1],
        colvals: vec![0.0; nnonzs],
        colind: vec![0; nnonzs],
        violrows: vec![0; ndisc],
        transform: vec![Transform::None; ndisc],
        lhs: vec![0.0; nrows],
        rhs: vec![0.0; nrows],
        colnorms: vec![0.0; ndisc],
        upperbounds: vec![0.0; ndisc],
        shiftvals: vec![0.0; ndisc],
        nrows,
        ndisc,
        feastol: ft,
    };
    for j in 0..ndisc {
        m.colbegin[j + 1] = m.colbegin[j] + colcount[j];
    }
    let mut fill = m.colbegin.clone();
    let mut rowmax = vec![0.0f64; nrows];
    for i in 0..nrows {
        let (s, e) = (a_rstart[i] as usize, a_rstart[i + 1] as usize);
        m.rowbegin.push(m.rowind.len());
        let l = ms.row_lower_2(i as HighsInt);
        let u = ms.row_upper_2(i as HighsInt);
        m.lhs[i] = if l <= -S6_INFINITY { -S6_INFINITY } else { l };
        m.rhs[i] = if u >= S6_INFINITY { S6_INFINITY } else { u };
        // empty row with an unsatisfiable side: infeasible
        if s == e && (s6_feas_gt(m.lhs[i], 0.0, ft) || s6_feas_lt(m.rhs[i], 0.0, ft)) {
            return false;
        }
        for k in s..e {
            rowmax[i] = rowmax[i].max(a_rvalue[k].abs());
            let p = colpos_of[a_rindex[k] as usize];
            if p == usize::MAX {
                continue;
            }
            m.rowind.push(p);
            m.colvals[fill[p]] = a_rvalue[k];
            m.colind[fill[p]] = i;
            fill[p] += 1;
        }
    }
    m.rowbegin.push(m.rowind.len());
    let mut nmaxrows = 0usize;
    for j in 0..ndisc {
        let ncolnonz = colcount[j];
        m.colnorms[j] = ncolnonz as f64;
        nmaxrows = nmaxrows.max(ncolnonz);
        for i in m.colbegin[j]..m.colbegin[j + 1] {
            let mut normval = m.colvals[i];
            if NORMALIZE {
                normval /= rowmax[m.colind[i]];
            }
            m.colnorms[j] += normval.abs();
        }
    }
    // discrete variables are transformed, continuous ones relaxed out of their rows
    for j in 0..ndisc {
        m.shiftvals[j] = 0.0;
        let (lb, ub) = bounds(dom, disc_cols[j]);
        m.transform_variable(j, lb, ub);
    }
    for i in 0..nrows {
        let (s, e) = (a_rstart[i] as usize, a_rstart[i + 1] as usize);
        for k in s..e {
            let col = a_rindex[k] as usize;
            if colpos_of[col] != usize::MAX {
                continue;
            }
            let colval = a_rvalue[k];
            let (lb, ub) = bounds(dom, col);
            let (lhsvarbound, rhsvarbound) = if colval > ft {
                (ub, lb)
            } else if colval < -ft {
                (lb, ub)
            } else {
                continue;
            };
            if !is_inf(-m.lhs[i]) && !is_inf(lhsvarbound.abs()) {
                m.lhs[i] -= colval * lhsvarbound;
            } else {
                m.lhs[i] = -S6_INFINITY;
            }
            if !is_inf(m.rhs[i]) && !is_inf(rhsvarbound.abs()) {
                m.rhs[i] -= colval * rhsvarbound;
            } else {
                m.rhs[i] = S6_INFINITY;
            }
        }
    }

    // ---- row weights, initial violations ----
    let mut viol = Viol {
        rows: Vec::new(),
        pos: vec![-1; nrows],
        weights: (0..nrows)
            .map(|r| {
                if !is_inf(-m.lhs[r]) && !is_inf(m.rhs[r]) {
                    WEIGHT_EQUALITY
                } else {
                    WEIGHT_INEQUALITY
                }
            })
            .collect(),
    };
    m.check_violations(&mut viol, None);

    // ---- variable order ----
    let mut permutation: Vec<usize> = (0..ndisc).collect();
    let violatedvarrows: Vec<i32> = m.violrows.clone();
    if SORTVARS {
        // SCIPsortDownRealInt / SCIPsortRealInt on the column norms,
        // SCIPsortDownIntInt / SCIPsortIntInt on the violated-row counts
        let (head, tail) = if PREFERBINARIES {
            (nbin, ndisc)
        } else {
            (ndisc, ndisc)
        };
        let ranges: Vec<(usize, usize)> = if PREFERBINARIES {
            vec![(0, head), (head, tail)]
        } else {
            vec![(0, ndisc)]
        };
        for (a, b) in ranges {
            if a == b {
                continue;
            }
            match SORTKEY {
                'n' => {
                    let mut keys = m.colnorms[a..b].to_vec();
                    sorttpl_sort(&mut keys, &mut permutation[a..b], &|x: f64, y: f64| {
                        y - x < 0.0
                    });
                }
                'u' => {
                    let mut keys = m.colnorms[a..b].to_vec();
                    sorttpl_sort(&mut keys, &mut permutation[a..b], &|x: f64, y: f64| {
                        x - y < 0.0
                    });
                }
                't' => {
                    let mut keys = violatedvarrows[a..b].to_vec();
                    sorttpl_sort(&mut keys, &mut permutation[a..b], &|x: i32, y: i32| {
                        x - y < 0
                    });
                }
                _ => {
                    let mut keys = violatedvarrows[a..b].to_vec();
                    sorttpl_sort(&mut keys, &mut permutation[a..b], &|x: i32, y: i32| {
                        y - x < 0
                    });
                }
            }
        }
    }

    // ---- main loop ----
    let mut sol = vec![0.0f64; ncol];
    let mut steps: Vec<f64> = Vec::with_capacity(nmaxrows);
    let mut violationchange: Vec<i32> = Vec::with_capacity(nmaxrows);
    let mut probing = Probing { depth: 1 };
    let mut cutoff = false;
    let mut lastindexofsusp: i64 = -1;
    let mut do_probing = true;
    let mut nprobings: i64 = 0;
    let mut ncutoffs: i32 = 0;
    let time_limit = ms.options_mip_.time_limit;

    for ci in 0..ndisc {
        let permutedvarindex = permutation[ci];
        let col = disc_cols[permutedvarindex];
        let (lb, ub) = bounds(dom, col);
        let mut optimalshiftvalue = 0.0;
        let mut nviolations = 0;

        if ci % STOP_CHECK_EVERY == 0
            && (ms.terminate() || (time_limit < kHighsInf && ms.timer_.read() >= time_limit))
        {
            stats.probings += nprobings;
            stats.cutoffs += ncutoffs as i64;
            stats.violations_left = viol.rows.len() as i64;
            return false;
        }

        // propagation may have changed the bounds: update the transformed upper bound
        m.update_transformation(&mut viol, permutedvarindex, lb, ub);

        // ignore the variable if propagation fixed it
        if s6_feas_zero(m.upperbounds[permutedvarindex], ft) {
            sol[col] = ub;
            continue;
        }

        let mut marksuspicious = false;
        let uplocks = ms.mipdata_.uplocks[col];
        let downlocks = ms.mipdata_.downlocks[col];
        #[allow(clippy::needless_late_init)] // mirrors the C control flow
        let origsolval;
        if FIXBINLOCKS && is_bin[col] && (uplocks == 0 || downlocks == 0) {
            origsolval = if uplocks == 0 { ub } else { lb };
        } else {
            // the expensive best shift selection only while a row is violated
            if !viol.rows.is_empty() {
                let (s, v) = m.optimal_shifting_value(
                    permutedvarindex,
                    1,
                    &viol.weights,
                    &mut steps,
                    &mut violationchange,
                );
                optimalshiftvalue = s;
                nviolations = v;
                // free variables are handled twice, once per direction
                if m.transform[permutedvarindex] == Transform::Free {
                    let (ds, dv) = m.optimal_shifting_value(
                        permutedvarindex,
                        -1,
                        &viol.weights,
                        &mut steps,
                        &mut violationchange,
                    );
                    if dv < nviolations {
                        optimalshiftvalue = ds;
                    }
                }
            }
            origsolval = m.retransform(permutedvarindex, optimalshiftvalue);
        }

        if nprobings > PROPBREAKER {
            do_probing = false;
        }

        // fix the variable to the shifted value and propagate the fixing
        if !marksuspicious && do_probing {
            let node_start = dom.get_branching_positions().len();
            let from = dom.get_domain_change_stack().len();
            probing.depth += 1;
            dom.fix_col(
                ms,
                col as HighsInt,
                origsolval,
                HighsDomain_Reason::branching(),
            );
            if !dom.infeasible() {
                let _ = dom.propagate(ms);
            }
            nprobings += 1;
            cutoff = dom.infeasible();
            if !cutoff {
                sync_changes_cols(&mut m, &mut viol, dom, &colpos_of, &disc_cols, from);
            }

            // propagation led to an empty domain: backtrack and postpone the variable
            if cutoff {
                ncutoffs += 1;
                // only continue if the number of cutoffs is reasonably small
                if CUTOFFBREAKER >= 0
                    && ncutoffs as f64
                        >= (MAXCUTOFFQUOT * probing.depth as f64) + CUTOFFBREAKER as f64
                {
                    break;
                }
                cutoff = false;
                backtrack_to(
                    ms, dom, &mut m, &mut viol, &colpos_of, &disc_cols, node_start,
                );
                probing.depth -= 1;

                let (clb, cub) = bounds(dom, col);
                if feas_eq(cub, origsolval, ft) && feas_eq(clb, origsolval, ft) {
                    // trapped at an infeasible node
                    cutoff = true;
                    break;
                } else if feas_eq(clb, origsolval, ft) && origsolval.abs() < 1.0 / S6_EPSILON {
                    // at the lower bound: repropagate with the bound tightened by one
                    let from2 = dom.get_domain_change_stack().len();
                    probing.depth += 1;
                    dom.change_bound_2(
                        ms,
                        HighsBoundType::kLower,
                        col as HighsInt,
                        origsolval + 1.0,
                        HighsDomain_Reason::branching(),
                    );
                    if !dom.infeasible() {
                        let _ = dom.propagate(ms);
                    }
                    cutoff = dom.infeasible();
                    if !cutoff {
                        sync_changes_cols(&mut m, &mut viol, dom, &colpos_of, &disc_cols, from2);
                    }
                } else if feas_eq(cub, origsolval, ft) && origsolval.abs() < 1.0 / S6_EPSILON {
                    let from2 = dom.get_domain_change_stack().len();
                    probing.depth += 1;
                    dom.change_bound_2(
                        ms,
                        HighsBoundType::kUpper,
                        col as HighsInt,
                        origsolval - 1.0,
                        HighsDomain_Reason::branching(),
                    );
                    if !dom.infeasible() {
                        let _ = dom.propagate(ms);
                    }
                    cutoff = dom.infeasible();
                    if !cutoff {
                        sync_changes_cols(&mut m, &mut viol, dom, &colpos_of, &disc_cols, from2);
                    }
                }

                // both subproblems infeasible: stop; else the variable led to a cutoff in one direction
                if cutoff {
                    break;
                }
                marksuspicious = true;
            }
        }

        if marksuspicious {
            // postpone the variable to the front of the permutation
            lastindexofsusp += 1;
            let l = lastindexofsusp as usize;
            debug_assert!(l <= ci);
            permutation[ci] = permutation[l];
            permutation[l] = permutedvarindex;
        } else {
            sol[col] = origsolval;
            if !do_probing {
                // only to keep the local domain consistent (SCIPfixVarProbing)
                dom.fix_col(
                    ms,
                    col as HighsInt,
                    origsolval,
                    HighsDomain_Reason::branching(),
                );
            }
        }
        let _ = nviolations;
    }

    stats.probings += nprobings;
    stats.cutoffs += ncutoffs as i64;
    stats.violations_left = viol.rows.len() as i64;

    if !(viol.rows.is_empty() && !cutoff) {
        return false;
    }

    // set the postponed variables: their transformed value is zero
    for v in 0..=lastindexofsusp {
        if v < 0 {
            break;
        }
        let p = permutation[v as usize];
        let col = disc_cols[p];
        let (lb, ub) = bounds(dom, col);
        m.update_transformation(&mut viol, p, lb, ub);
        let origsolval = m.retransform(p, 0.0);
        sol[col] = origsolval;
        let from = dom.get_domain_change_stack().len();
        dom.fix_col(
            ms,
            col as HighsInt,
            origsolval,
            HighsDomain_Reason::branching(),
        );
        if dom.infeasible() {
            stats.violations_left = viol.rows.len() as i64;
            return false;
        }
        sync_changes_cols(&mut m, &mut viol, dom, &colpos_of, &disc_cols, from);
    }
    stats.violations_left = viol.rows.len() as i64;
    // fixing the remaining variables led to infeasibility
    if !viol.rows.is_empty() {
        return false;
    }

    // ---- solution: pure integer, or the LP over the continuous columns ----
    let ub_before = ms.mipdata_.upper_bound;
    if ndisc == ncol {
        let _ = ms.try_solution(&sol, kSolutionSourceShiftProp);
    } else {
        // SCIPsolveProbingLP with all discrete variables fixed
        let r = detach!(ms.mipdata_.heuristics, |h| {
            h.try_rounded_point(ms, 0, &sol, kSolutionSourceShiftProp)
        });
    }
    ms.mipdata_.upper_bound < ub_before
}
