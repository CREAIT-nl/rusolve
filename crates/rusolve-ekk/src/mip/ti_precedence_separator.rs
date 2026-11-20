//! Time-indexed scheduling separator (rusolve lever `mip_ti_precedence_sepa`, default off).
//!
//! Not a HiGHS file. Provenance (literature, no solver source was copied):
//! * step (disaggregated) precedence inequalities of the time-indexed formulation, Sousa and
//!   Wolsey, "A time indexed formulation of non-preemptive single machine scheduling
//!   problems", Math. Programming 54 (1992);
//! * machine work inequalities for the makespan: the one-machine relaxation with release
//!   dates and tails, Carlier, "The one-machine sequencing problem", EJOR 11 (1982); in the
//!   time-indexed space Dyer and Wolsey, "Formulating the single machine sequencing problem
//!   with release dates as a mixed integer program", Discrete Appl. Math. 26 (1990).
//!
//! Detected structure (no names are read): "assignment sets" are binaries with an equality row
//! `sum x = 1` (one start time per operation). A row `g1 + g2 >= p`, `p > 0`, with
//! `g_s = sum g_i x_i` (a member missing from the row has `g_i = 0`) is an aggregated precedence
//! `sum t x_b - sum t x_a >= p`; a row `C + g(S) >= p` with one continuous `C` is the makespan link
//! of a last operation (`g = -t`); `sum x <= 1` over unit coefficients on two or more sets is a
//! machine time slot.
//!
//! Precedence cut. With `f1 = g1`, `f2 = -g2` every integer point satisfies, for every `theta`,
//! `sum_{f1_i <= theta} x_i - sum_{f2_i <= theta - p} x_i <= 0` (set 1 at `f1 <= theta` forces
//! `f2 <= f1 - p <= theta - p`). Thresholds `theta in {f1_i}` suffice.
//!
//! Machine cut. Ops sharing time slots form a machine with processing times `p_o` and start
//! times `t_i`. All work of the machine that lies at or after `tau` needs one time unit per
//! unit of work, so the last operation completes at `tau + W(tau)`, `W(tau) = sum x_i *
//! min(p_o, max(0, t_i + p_o - tau))`, and the makespan covers that completion plus the
//! smallest tail `q` (work that must follow an operation on its job): `C >= tau + q + W(tau)`.
//! The tail comes from the precedence chain to the makespan rows; if any chain is unknown `q = 0`.

use crate::prelude::*;
use std::collections::HashMap;

/// Cuts added per precedence row and call.
const CUTS_PER_ROW: usize = 3;
/// Machine cuts added per machine and call.
const CUTS_PER_MACHINE: usize = 2;
/// Cuts added per call.
const MAX_CUTS_PER_CALL: usize = 4_000;
/// Minimum violation of a cut to add it.
const MIN_VIOLATION: f64 = 1e-3;
/// Evaluations without a new best before the walk restarts.
const STALL_RESTART: usize = 1500;
/// Column-scan budget of one list scheduling call.
const HILL_CLIMB_WORK: f64 = 2.0e8;

/// One detected precedence row. `s1` is sorted by `f1`, `s2` by `f2`.
#[derive(Clone, Debug, Default)]
struct TiPair {
    s1: Vec<(f64, HighsInt)>,
    s2: Vec<(f64, HighsInt)>,
    p: f64,
}

/// One machine: operations (processing time, members as (start time, column)).
#[derive(Clone, Debug, Default)]
struct TiMachine {
    ops: Vec<(f64, Vec<(f64, HighsInt)>)>,
    tail: f64,
}

/// One operation for list scheduling: processing time, start times (time, column), job
/// predecessor (its end is a lower bound on the start), machine id (-1: unconstrained).
#[derive(Clone, Debug, Default)]
struct TiOp {
    p: f64,
    mem: Vec<(f64, HighsInt)>,
    pred: i32,
    machine: i32,
    last: bool,
    release: f64,
}

#[derive(Clone, Default, Debug)]
pub struct HighsTiPrecedenceSeparator {
    pub base: HighsSeparator,
    pairs: Vec<TiPair>,
    machines: Vec<TiMachine>,
    makespan: HighsInt,
    detected: bool,
    ops: Vec<TiOp>,
    /// slot clique rows (columns) and, per column, the ids of the cliques containing it
    cliques: Vec<Vec<HighsInt>>,
    col_cliques: Vec<Vec<u32>>,
    n_comp: usize,
    heur_calls: u32,
    rng: u64,
    /// see `HighsNetdesignSeparator::pending`
    pub pending: Vec<(u64, Vec<f64>)>,
}

fn find(uf: &mut Vec<usize>, x: usize) -> usize {
    let mut r = x;
    while uf[r] != r {
        r = uf[r];
    }
    let mut y = x;
    while uf[y] != r {
        let n = uf[y];
        uf[y] = r;
        y = n;
    }
    r
}

impl HighsTiPrecedenceSeparator {
    pub fn new(ms: &HighsMipSolver) -> Self {
        HighsTiPrecedenceSeparator {
            base: HighsSeparator::new(ms, "Separation: TI precedence"),
            pairs: Vec::new(),
            machines: Vec::new(),
            makespan: -1,
            detected: false,
            ops: Vec::new(),
            cliques: Vec::new(),
            col_cliques: Vec::new(),
            n_comp: 0,
            heur_calls: 0,
            rng: 0x9E37_79B9_7F4A_7C15,
            pending: Vec::new(),
        }
    }

    /// Finds the assignment sets, precedence rows, makespan links and machines (once).
    fn detect(&mut self, ms: &HighsMipSolver, lp: &HighsLpRelaxation, g: usize) {
        self.detected = true;
        let nrow = ms.num_row();
        let ncol = ms.num_col() as usize;
        let lower = &ms.mipdata_.domains[g].col_lower_;
        let upper = &ms.mipdata_.domains[g].col_upper_;
        let model = lp.get_lp();
        let is_bin = |col: usize| {
            !matches!(
                ms.variable_type_2(col as HighsInt),
                HighsVarType::kContinuous
            ) && lower[col] == 0.0
                && upper[col] == 1.0
        };
        let get = |r: HighsInt| -> (&[HighsInt], &[f64]) {
            let (mut len, mut inds, mut vals): (HighsInt, &[HighsInt], &[f64]) = (0, &[], &[]);
            lp.get_row(ms, r, &mut len, &mut inds, &mut vals);
            (inds, vals)
        };
        // assignment sets: rows `sum x = 1` over binaries with unit coefficients
        let mut set_of: Vec<i32> = vec![-1; ncol];
        let mut sets: Vec<Vec<HighsInt>> = Vec::new();
        for r in 0..nrow {
            let (rl, ru) = (model.row_lower_[r as usize], model.row_upper_[r as usize]);
            if rl != 1.0 || ru != 1.0 {
                continue;
            }
            let (inds, vals) = get(r);
            if inds.len() < 2
                || vals.iter().any(|&v| v != 1.0)
                || inds
                    .iter()
                    .any(|&cc| !is_bin(cc as usize) || set_of[cc as usize] != -1)
            {
                continue;
            }
            let id = sets.len() as i32;
            for &cc in inds {
                set_of[cc as usize] = id;
            }
            sets.push(inds.to_vec());
        }
        if sets.is_empty() {
            return;
        }
        let nsets = sets.len();
        let mut times: Vec<HashMap<HighsInt, f64>> = vec![HashMap::new(); nsets];
        let mut p_of: Vec<f64> = vec![f64::NAN; nsets];
        let mut succ: Vec<i32> = vec![-1; nsets];
        let mut last: Vec<bool> = vec![false; nsets];
        let mut release: Vec<f64> = vec![f64::NEG_INFINITY; nsets];
        let mut clq_count: Vec<u32> = vec![0; ncol];
        let mut clique_sets: Vec<Vec<i32>> = Vec::new();
        let mut all_cliques: Vec<Vec<HighsInt>> = Vec::new();
        let mut pair_set: std::collections::HashSet<(HighsInt, HighsInt)> =
            std::collections::HashSet::new();
        let mut makespan: HighsInt = -1;
        let mut makespan_ok = true;
        for r in 0..nrow {
            let (rl, ru) = (model.row_lower_[r as usize], model.row_upper_[r as usize]);
            let (inds, vals) = get(r);
            if inds.len() < 2 {
                continue;
            }
            // machine time slot: sum x <= 1, unit coefficients over two or more sets
            if rl <= -kHighsInf
                && ru == 1.0
                && vals.iter().all(|&v| v == 1.0)
                && inds.iter().all(|&cc| set_of[cc as usize] >= 0)
            {
                let mut ss: Vec<i32> = inds.iter().map(|&cc| set_of[cc as usize]).collect();
                ss.sort_unstable();
                ss.dedup();
                if ss.len() >= 2 {
                    for &cc in inds {
                        clq_count[cc as usize] += 1;
                    }
                    for (x, &c1) in inds.iter().enumerate() {
                        for &c2 in &inds[x + 1..] {
                            pair_set.insert((c1.min(c2), c1.max(c2)));
                        }
                    }
                    clique_sets.push(ss);
                    all_cliques.push(inds.to_vec());
                }
                continue;
            }
            let (sign, p) = if rl.is_finite() && ru >= kHighsInf {
                (1.0, rl)
            } else if ru.is_finite() && rl <= -kHighsInf {
                (-1.0, -ru)
            } else {
                continue;
            };
            if p <= 0.0 {
                continue;
            }
            // makespan link: one continuous column with coefficient 1 and one set
            let cont: Vec<usize> = (0..inds.len())
                .filter(|&k| set_of[inds[k] as usize] < 0)
                .collect();
            if cont.len() == 1 {
                let k = cont[0];
                let col = inds[k];
                if sign * vals[k] != 1.0
                    || matches!(ms.variable_type_2(col), HighsVarType::kInteger)
                {
                    continue;
                }
                let sid = set_of[inds[(k + 1) % inds.len()] as usize];
                if sid < 0
                    || inds
                        .iter()
                        .enumerate()
                        .any(|(q, &cc)| q != k && set_of[cc as usize] != sid)
                {
                    continue;
                }
                if makespan != -1 && makespan != col {
                    makespan_ok = false;
                }
                makespan = col;
                last[sid as usize] = true;
                p_of[sid as usize] = p;
                for (q, &cc) in inds.iter().enumerate() {
                    if q != k {
                        times[sid as usize].insert(cc, -sign * vals[q]);
                    }
                }
                continue;
            }
            if cont.is_empty() {
                let mut ids: [i32; 2] = [-1, -1];
                let mut ok = true;
                for &cc in inds {
                    let s = set_of[cc as usize];
                    if s < 0 {
                        ok = false;
                        break;
                    }
                    if ids[0] == -1 || ids[0] == s {
                        ids[0] = s;
                    } else if ids[1] == -1 || ids[1] == s {
                        ids[1] = s;
                    } else {
                        ok = false;
                        break;
                    }
                }
                if ok && ids[1] == -1 && ids[0] >= 0 && sign > 0.0 && vals.iter().all(|&v| v > 0.0)
                {
                    // presolve fixed a predecessor: `sum t x >= r` is a release time
                    release[ids[0] as usize] = release[ids[0] as usize].max(rl);
                    continue;
                }
                if !ok || ids[1] == -1 {
                    continue;
                }
                let coef = |set: i32, mult: f64| -> Vec<(f64, HighsInt)> {
                    let mut v: Vec<(f64, HighsInt)> =
                        sets[set as usize].iter().map(|&cc| (0.0, cc)).collect();
                    for (k, &cc) in inds.iter().enumerate() {
                        if set_of[cc as usize] == set {
                            let pos = v.iter().position(|e| e.1 == cc).unwrap();
                            v[pos].0 = mult * sign * vals[k];
                        }
                    }
                    v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
                    v
                };
                let s1 = coef(ids[0], 1.0);
                let s2 = coef(ids[1], -1.0);
                // the set with positive coefficients in the `>= p` row is the successor b
                let sum1: f64 = s1.iter().map(|e| e.0).sum();
                let (a_set, b_set) = if sum1 > 0.0 {
                    (ids[1], ids[0])
                } else {
                    (ids[0], ids[1])
                };
                for (k, &cc) in inds.iter().enumerate() {
                    let t = sign * vals[k];
                    let sidx = set_of[cc as usize];
                    times[sidx as usize].insert(cc, if sidx == a_set { -t } else { t });
                }
                p_of[a_set as usize] = p;
                succ[a_set as usize] = b_set;
                self.pairs.push(TiPair { s1, s2, p });
            }
        }
        // tails along the precedence chains
        let mut tail: Vec<f64> = vec![f64::NAN; nsets];
        for s in 0..nsets {
            let mut chain = vec![s];
            let mut cur = s;
            let mut guard = 0;
            let mut t_end = f64::NAN;
            loop {
                if last[cur] {
                    t_end = 0.0;
                    break;
                }
                let nx = succ[cur];
                guard += 1;
                if nx < 0 || guard > nsets {
                    break;
                }
                cur = nx as usize;
                chain.push(cur);
            }
            if t_end.is_nan() {
                continue;
            }
            // accumulate from the chain end back to s: tail(prev) = p(cur) + tail(cur)
            let mut acc = 0.0;
            tail[*chain.last().unwrap()] = 0.0;
            for w in (0..chain.len() - 1).rev() {
                acc += p_of[chain[w + 1]];
                tail[chain[w]] = acc;
            }
        }
        // machines: union of sets that share a time slot
        let mut uf: Vec<usize> = (0..nsets).collect();
        for ss in &clique_sets {
            for &x in &ss[1..] {
                let (ra, rb) = (find(&mut uf, ss[0] as usize), find(&mut uf, x as usize));
                if ra != rb {
                    uf[ra] = rb;
                }
            }
        }
        let in_clique: Vec<bool> = {
            let mut v = vec![false; nsets];
            for ss in &clique_sets {
                for &x in ss {
                    v[x as usize] = true;
                }
            }
            v
        };
        let mut machine_of: Vec<i32> = vec![-1; nsets];
        let mut n_machines: i32 = 0;
        let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
        for s in 0..nsets {
            if in_clique[s] {
                let r = find(&mut uf, s);
                by_root.entry(r).or_default().push(s);
            }
        }
        if makespan_ok && makespan >= 0 {
            let mut roots: Vec<usize> = by_root.keys().copied().collect();
            roots.sort_unstable();
            for r in roots {
                let members = &by_root[&r];
                // usable ops of this component
                let mut usable: Vec<(usize, f64, Vec<(f64, HighsInt)>)> = Vec::new();
                for &s in members {
                    let p = p_of[s];
                    if p.is_nan() || p <= 0.0 || times[s].len() < 2 {
                        continue;
                    }
                    // columns fixed to 0 in the global domain cannot be chosen: leave them out
                    let mem: Vec<(f64, HighsInt)> = sets[s]
                        .iter()
                        .filter(|&&cc| upper[cc as usize] > 0.5)
                        .map(|&cc| (times[s].get(&cc).copied().unwrap_or(0.0), cc))
                        .collect();
                    if mem.len() < 2 {
                        continue;
                    }
                    usable.push((s, p, mem));
                }
                // exact no-overlap check: two ops share a machine only if every pair of
                // overlapping starts lies in one slot row
                let no_overlap = |a: &(usize, f64, Vec<(f64, HighsInt)>),
                                  b: &(usize, f64, Vec<(f64, HighsInt)>)|
                 -> bool {
                    for &(ta, ca) in &a.2 {
                        for &(tb, cb) in &b.2 {
                            if ta < tb + b.1
                                && tb < ta + a.1
                                && !pair_set.contains(&(ca.min(cb), ca.max(cb)))
                            {
                                return false;
                            }
                        }
                    }
                    true
                };
                let mut groups: Vec<Vec<usize>> = Vec::new();
                for i in 0..usable.len() {
                    let mut placed = false;
                    for gr in groups.iter_mut() {
                        if gr.iter().all(|&j| no_overlap(&usable[i], &usable[j])) {
                            gr.push(i);
                            placed = true;
                            break;
                        }
                    }
                    if !placed {
                        groups.push(vec![i]);
                    }
                }
                for gr in groups {
                    if gr.len() < 2 {
                        continue;
                    }
                    for &i in &gr {
                        machine_of[usable[i].0] = n_machines;
                    }
                    n_machines += 1;
                    let mut q = f64::INFINITY;
                    for &i in &gr {
                        let s = usable[i].0;
                        q = q.min(if tail[s].is_nan() { 0.0 } else { tail[s] });
                    }
                    // ops of the component that are not in the group may still share nothing
                    // with it; the tail minimum over the group only is valid (ops in the cut)
                    let ops: Vec<(f64, Vec<(f64, HighsInt)>)> = gr
                        .iter()
                        .map(|&i| (usable[i].1, usable[i].2.clone()))
                        .collect();
                    self.machines.push(TiMachine { ops, tail: q });
                }
            }
            self.makespan = makespan;
            // list scheduling data: every set needs a known processing time
            let mut pred_of: Vec<i32> = vec![-1; nsets];
            for a in 0..nsets {
                if succ[a] >= 0 {
                    pred_of[succ[a] as usize] = a as i32;
                }
            }
            if (0..nsets).all(|q| !p_of[q].is_nan() && p_of[q] > 0.0) {
                let mut cc_list: Vec<Vec<u32>> = vec![Vec::new(); ncol];
                for (id, cl) in all_cliques.iter().enumerate() {
                    for &cc in cl {
                        cc_list[cc as usize].push(id as u32);
                    }
                }
                self.cliques = std::mem::take(&mut all_cliques);
                self.col_cliques = cc_list;
                self.ops = (0..nsets)
                    .map(|q| TiOp {
                        p: p_of[q],
                        mem: sets[q]
                            .iter()
                            .filter(|&&cc| upper[cc as usize] > 0.5)
                            .map(|&cc| (times[q].get(&cc).copied().unwrap_or(0.0), cc))
                            .collect(),
                        pred: pred_of[q],
                        // list scheduling uses the slot-clique components as machines
                        machine: if in_clique[q] {
                            find(&mut uf, q) as i32
                        } else {
                            -1
                        },
                        last: last[q],
                        release: release[q],
                    })
                    .collect();
                self.n_comp = nsets;
                for o in self.ops.iter_mut() {
                    o.mem
                        .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
                }
            }
        }
        // Assignment rows alone are common outside scheduling (drayage, atlanta-ip). Without both
        // precedence rows and machine slots the schedule builder would only perturb the search
        //, so it stays empty.
        if self.pairs.is_empty() || self.machines.is_empty() {
            self.ops.clear();
        }
        if !self.pairs.is_empty() || !self.machines.is_empty() {
            crate::highs_log_user!(
                &ms.options_mip_.log_options,
                HighsLogType::kInfo,
                "TI separator: {} assignment sets, {} precedence rows, {} machines\n",
                sets.len(),
                self.pairs.len(),
                self.machines.len()
            );
        }
    }

    /// Serial schedule generation on the model's own columns: the eligible op with the smallest
    /// key takes its earliest start time (column) that is after its job predecessor, not before
    /// its release time and not in a slot clique already used (Kolisch 1996). An op with no such
    /// column is counted as a failure. Returns (failures, tie-break, chosen column per op).
    fn sgs(&self, key: &[f64]) -> (usize, f64, Vec<HighsInt>) {
        let n = self.ops.len();
        let mut start: Vec<f64> = vec![f64::NAN; n];
        let mut chosen: Vec<HighsInt> = vec![-1; n];
        let mut used: Vec<bool> = vec![false; self.cliques.len()];
        let mut fails = 0;
        // intervals of placed ops per machine component, to measure how far past the model's
        // start times a failed op would have to go (a gradient for the search)
        let mut busy: Vec<Vec<(f64, f64)>> = vec![Vec::new(); self.n_comp];
        let mut late = 0.0;
        for _ in 0..n {
            let mut pick = usize::MAX;
            for q in 0..n {
                if start[q].is_nan()
                    && (self.ops[q].pred < 0 || !start[self.ops[q].pred as usize].is_nan())
                    && (pick == usize::MAX || key[q] < key[pick])
                {
                    pick = q;
                }
            }
            let o = &self.ops[pick];
            let ready = if o.pred >= 0 {
                start[o.pred as usize] + self.ops[o.pred as usize].p
            } else {
                f64::NEG_INFINITY
            };
            let mut found = false;
            for &(t, cc) in &o.mem {
                if t + 1e-9 < ready || t + 1e-9 < o.release {
                    continue;
                }
                if self.col_cliques[cc as usize]
                    .iter()
                    .any(|&k| used[k as usize])
                {
                    continue;
                }
                start[pick] = t;
                chosen[pick] = cc;
                for &k in &self.col_cliques[cc as usize] {
                    used[k as usize] = true;
                }
                if o.machine >= 0 {
                    busy[o.machine as usize].push((t, t + o.p));
                }
                found = true;
                break;
            }
            if !found {
                fails += 1;
                // virtual start: earliest gap on its component after the ready time
                let mut t = ready.max(o.release).max(0.0);
                if o.machine >= 0 {
                    let b = &mut busy[o.machine as usize];
                    let mut cands: Vec<f64> = vec![t];
                    cands.extend(b.iter().map(|x| x.1).filter(|&e| e > t));
                    cands.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    for c0 in cands {
                        if !b
                            .iter()
                            .any(|&(b0, b1)| c0 < b1 - 1e-9 && b0 < c0 + o.p - 1e-9)
                        {
                            t = c0;
                            break;
                        }
                    }
                    b.push((t, t + o.p));
                }
                let latest = o.mem.last().map(|e| e.0).unwrap_or(0.0);
                late += (t - latest).max(0.0) + 1.0;
                start[pick] = t;
            }
        }
        let _ = &start;
        (fails, late, chosen)
    }

    /// LP-guided list scheduling (Kolisch 1996 serial SGS; alpha points: Hall, Schulz, Shmoys,
    /// Wein 1997) with a hill climb over the priority keys. Returns the full column vector of a
    /// schedule that fits the model, if one is found.
    fn list_schedule(&mut self, ncol: usize, sol: &[f64]) -> Option<Vec<f64>> {
        let n = self.ops.len();
        if n == 0 || self.makespan < 0 {
            return None;
        }
        let alpha_point = |o: &TiOp, alpha: f64| -> f64 {
            let mut cum = 0.0;
            for &(t, cc) in &o.mem {
                cum += sol[cc as usize].max(0.0);
                if cum >= alpha - 1e-9 {
                    return t;
                }
            }
            o.mem.last().map(|e| e.0).unwrap_or(0.0)
        };
        let span = self
            .ops
            .iter()
            .flat_map(|o| o.mem.iter().map(|e| e.0))
            .fold(0.0, f64::max)
            .max(1.0);
        let rnd = |rng: &mut u64| -> f64 {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*rng >> 33) as f64) / ((1u64 << 31) as f64)
        };
        let mut rng = self.rng;
        let mut best: Option<(usize, f64, Vec<HighsInt>, Vec<f64>)> = None;
        // work budget: columns scanned
        let cost_per_eval = (self.ops.iter().map(|o| o.mem.len()).sum::<usize>() as f64).max(1.0);
        let evals = ((HILL_CLIMB_WORK / cost_per_eval) as usize).clamp(20, 30000);
        for alpha in [0.5, 0.0, 1.0] {
            let key: Vec<f64> = self.ops.iter().map(|o| alpha_point(o, alpha)).collect();
            let (f, tie, ch) = self.sgs(&key);
            if best.as_ref().is_none_or(|b| (f, tie) < (b.0, b.1)) {
                best = Some((f, tie, ch, key));
            }
        }
        let (mut bf, mut bt, mut bc, bk) = best?;
        // current point of the walk (restarts jump to a fresh alpha-point order)
        let (mut cf, mut ct, mut ck) = (bf, bt, bk);
        let mut stall = 0;
        for _ in 0..evals {
            if bf == 0 {
                break;
            }
            let mut key = ck.clone();
            let flips = 1 + (rnd(&mut rng) * 3.0) as usize;
            for _ in 0..flips {
                let q = (rnd(&mut rng) * n as f64) as usize % n;
                key[q] += (rnd(&mut rng) - 0.5) * span * 0.3;
            }
            let (f, tie, ch) = self.sgs(&key);
            if (f, tie) < (bf, bt) {
                bf = f;
                bt = tie;
                bc = ch.clone();
                stall = 0;
            } else {
                stall += 1;
            }
            if (f, tie) <= (cf, ct) {
                cf = f;
                ct = tie;
                ck = key;
            }
            if stall > STALL_RESTART {
                // restart: alpha-point order for a random alpha, plus noise
                let alpha = rnd(&mut rng);
                ck = self
                    .ops
                    .iter()
                    .map(|o| alpha_point(o, alpha) + (rnd(&mut rng) - 0.5) * span * 0.1)
                    .collect();
                let (f, tie, _) = self.sgs(&ck);
                cf = f;
                ct = tie;
                stall = 0;
            }
        }
        self.rng = rng;
        if bf > 0 {
            return None;
        }
        let mut x = vec![0.0; ncol];
        let mut cmax = f64::NEG_INFINITY;
        for (q, o) in self.ops.iter().enumerate() {
            x[bc[q] as usize] = 1.0;
            if o.last {
                let t = o
                    .mem
                    .iter()
                    .find(|e| e.1 == bc[q])
                    .map(|e| e.0)
                    .unwrap_or(0.0);
                cmax = cmax.max(t + o.p);
            }
        }
        x[self.makespan as usize] = cmax;
        Some(x)
    }

    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        _lp_aggregator: &mut HighsLpAggregator,
        trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        // sub-MIPs fix and remove variables, which shifts row sides: the structure is read from
        // the main problem only
        if ms.submip {
            return;
        }
        if !self.detected {
            self.detect(ms, lp_relaxation, trans_lp.globaldom_);
        }
        let sol: Vec<f64> = lp_relaxation.get_solution().col_value.clone();
        if ms.options_mip_.mip_ti_list_sched && self.heur_calls < 20 && !self.ops.is_empty() {
            self.heur_calls += 1;
            if let Some(x) = self.list_schedule(sol.len(), &sol) {
                crate::highs_log_user!(
                    &ms.options_mip_.log_options,
                    HighsLogType::kInfo,
                    "TI list scheduling: makespan {}\n",
                    x[self.makespan as usize]
                );
                self.pending
                    .push((crate::mip::separation::next_pending_seq(), x));
            }
        }
        if !ms.options_mip_.mip_ti_precedence_sepa {
            return;
        }
        let mut cuts: Vec<(Vec<HighsInt>, Vec<f64>, f64)> = Vec::new();
        {
            for pair in &self.pairs {
                let mut cand: Vec<(f64, usize)> = Vec::new();
                let (mut sum1, mut sum2) = (0.0, 0.0);
                let (mut i, mut j) = (0usize, 0usize);
                while i < pair.s1.len() {
                    let theta = pair.s1[i].0;
                    while i < pair.s1.len() && pair.s1[i].0 <= theta {
                        sum1 += sol[pair.s1[i].1 as usize];
                        i += 1;
                    }
                    while j < pair.s2.len() && pair.s2[j].0 <= theta - pair.p {
                        sum2 += sol[pair.s2[j].1 as usize];
                        j += 1;
                    }
                    if sum1 - sum2 > MIN_VIOLATION {
                        cand.push((sum1 - sum2, i));
                    }
                }
                cand.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
                for &(_, end) in cand.iter().take(CUTS_PER_ROW) {
                    let theta = pair.s1[end - 1].0;
                    let mut inds: Vec<HighsInt> = Vec::new();
                    let mut vals: Vec<f64> = Vec::new();
                    for &(_, cc) in &pair.s1[..end] {
                        inds.push(cc);
                        vals.push(1.0);
                    }
                    for &(_, cc) in pair.s2.iter().take_while(|e| e.0 <= theta - pair.p) {
                        inds.push(cc);
                        vals.push(-1.0);
                    }
                    cuts.push((inds, vals, 0.0));
                }
                if cuts.len() >= MAX_CUTS_PER_CALL {
                    break;
                }
            }
        }
        if self.makespan >= 0 {
            let cval = sol[self.makespan as usize];
            for m in &self.machines {
                let tmin = m
                    .ops
                    .iter()
                    .flat_map(|o| o.1.iter().map(|e| e.0))
                    .fold(f64::INFINITY, f64::min);
                let tmax = m
                    .ops
                    .iter()
                    .flat_map(|o| o.1.iter().map(|e| e.0 + o.0))
                    .fold(0.0, f64::max);
                let mut cand: Vec<(f64, f64)> = Vec::new();
                let mut tau = tmin.floor();
                while tau <= tmax {
                    let mut w = 0.0;
                    for (p, mem) in &m.ops {
                        for &(t, cc) in mem {
                            let k = (t + p - tau).clamp(0.0, *p);
                            if k > 0.0 {
                                w += k * sol[cc as usize];
                            }
                        }
                    }
                    let viol = tau + m.tail + w - cval;
                    if viol > MIN_VIOLATION {
                        cand.push((viol, tau));
                    }
                    tau += 1.0;
                }
                cand.sort_by(|a, b| {
                    b.0.partial_cmp(&a.0)
                        .unwrap()
                        .then(a.1.partial_cmp(&b.1).unwrap())
                });
                for &(_, tau) in cand.iter().take(CUTS_PER_MACHINE) {
                    // C - sum w x >= tau + q   <=>   -C + sum w x <= -(tau + q)
                    let mut inds: Vec<HighsInt> = vec![self.makespan];
                    let mut vals: Vec<f64> = vec![-1.0];
                    for (p, mem) in &m.ops {
                        for &(t, cc) in mem {
                            let k = (t + p - tau).clamp(0.0, *p);
                            if k > 0.0 {
                                inds.push(cc);
                                vals.push(k);
                            }
                        }
                    }
                    cuts.push((inds, vals, -(tau + m.tail)));
                }
            }
        }
        for (mut inds, mut vals, rhs) in cuts {
            let len = inds.len() as HighsInt;
            detach!(ms.mipdata_.cutpools[cutpool], |p| p.add_cut(
                ms,
                &mut inds[..],
                &mut vals[..],
                len,
                rhs,
                true,
                true,
                false,
                false,
            ));
        }
    }
}

/// Synthetic time-indexed job shop (`jobs` x `machines`, horizon `horizon`) for the criterion
/// bench of the list scheduling hot loop (`benches/ti_list_sched.rs`).
pub struct TiBenchInstance {
    sep: HighsTiPrecedenceSeparator,
    ncol: usize,
}

impl TiBenchInstance {
    pub fn new(jobs: usize, machines: usize, horizon: usize) -> Self {
        let mut sep = HighsTiPrecedenceSeparator::default();
        sep.rng = 0x9E37_79B9_7F4A_7C15;
        let mut ncol = 0usize;
        let mut slot_cliques: HashMap<(usize, usize), Vec<HighsInt>> = HashMap::new();
        for j in 0..jobs {
            for k in 0..machines {
                let m = (j + k) % machines;
                let p = 1 + (j * 7 + k * 3) % 9;
                let mut mem = Vec::new();
                for t in 0..horizon.saturating_sub(p) {
                    let col = ncol as HighsInt;
                    ncol += 1;
                    mem.push((t as f64, col));
                    for slot in t..t + p {
                        slot_cliques.entry((m, slot)).or_default().push(col);
                    }
                }
                sep.ops.push(TiOp {
                    p: p as f64,
                    mem,
                    pred: if k == 0 {
                        -1
                    } else {
                        (j * machines + k - 1) as i32
                    },
                    machine: m as i32,
                    last: k + 1 == machines,
                    release: f64::NEG_INFINITY,
                });
            }
        }
        sep.makespan = ncol as HighsInt;
        ncol += 1;
        sep.n_comp = machines;
        let mut keys: Vec<_> = slot_cliques.into_iter().collect();
        keys.sort_by_key(|e| e.0);
        sep.col_cliques = vec![Vec::new(); ncol];
        for (id, (_, cols)) in keys.into_iter().enumerate() {
            for &cc in &cols {
                sep.col_cliques[cc as usize].push(id as u32);
            }
            sep.cliques.push(cols);
        }
        TiBenchInstance { sep, ncol }
    }

    /// One list scheduling call against an all-zero LP point; true if a schedule was found.
    pub fn run(&mut self) -> bool {
        let sol = vec![0.0; self.ncol];
        self.sep.list_schedule(self.ncol, &sol).is_some()
    }
}
