//! rusolve `mip_cdcl_uip` / `mip_zero_obj_cdcl`: a CDCL engine for feasibility.
//!
//! Sources: Marques-Silva and Sakallah, "GRASP",
//! IEEE Trans. Comput. 1999 (first-UIP conflict analysis); Moskewicz et al., "Chaff", DAC 2001
//! (two watched literals, VSIDS); Een and Soerensson, "An extensible SAT-solver", SAT 2003 (trail,
//! heap order, phase saving, activity decay); Soerensson and Biere, "Minimizing learned clauses",
//! SAT 2009 (recursive minimisation); Achterberg, "Conflict analysis in mixed integer
//! programming", Discrete Optimization 4 (2007) (explaining a linear row by the bound changes that
//! made its minimal activity exceed the side); Achterberg, Berthold, "Improving the feasibility
//! pump" and Berthold-Hendel, "Shift-and-Propagate" 2015 (propagation-driven feasibility search);
//! Luby, Sinclair, Zuckerman 1993 (restarts).
//!
//! Model: rows `sum a_j x_j <= rhs` (a range row is two). Binary columns are decision variables
//! with literal `x = 0/1`; clauses are over binaries only and are propagated with two watched
//! literals. Other columns (general integer, continuous, finite bounds) are *numeric*: they are
//! never decided; a row tightens their bounds by minimal activity, each tightening is a trail
//! entry with a lazily computed reason, and conflict analysis resolves through those entries so a
//! learned clause contains only binary literals. When every binary is assigned the numeric
//! propagation is run to a tight fixpoint and every numeric column is set to its lower bound (a
//! solution for difference-type rows). No LP, no objective.

/// Literal `2 * var + value`: "x_var = value".
type Lit = u32;
const EPS: f64 = 1e-7;
const HIGH: u32 = 1 << 31;
const NONE: u32 = u32::MAX;

#[derive(Clone, Copy, PartialEq)]
enum Reason {
    Decision,
    Clause(u32),
    Row(u32),
}

enum Conflict {
    Clause(u32),
    Row(u32),
    /// a row whose rounded implication empties a numeric domain
    Round(u32),
}

/// A node of the implication graph: binary literal, or numeric bound entry (`HIGH | index`).
type Node = u32;

struct Entry {
    var: u32,
    upper: bool,
    old: f64,
    new: f64,
    row: u32,
    pos: u32,
    level: u32,
    prev: u32,
}

pub struct UipEngine {
    n: usize,
    is_num: Vec<bool>,
    is_int: Vec<bool>,
    lb0: Vec<f64>,
    ub0: Vec<f64>,
    nlb: Vec<f64>,
    nub: Vec<f64>,
    entries: Vec<Entry>,
    last_entry: Vec<[u32; 2]>,
    eseen: Vec<u8>,
    exact: bool,
    exact_count: usize,
    // rows
    cstart: Vec<usize>,
    cvar: Vec<u32>,
    ccoef: Vec<f64>,
    crhs: Vec<f64>,
    cmaxabs: Vec<f64>,
    cbase: Vec<f64>,
    minact: Vec<f64>,
    queued: Vec<bool>,
    dirty: Vec<u32>,
    // occurrences per var
    ostart: Vec<usize>,
    ocon: Vec<u32>,
    ocoef: Vec<f64>,
    // assignment of binaries
    assign: Vec<i8>,
    level: Vec<u32>,
    pos: Vec<u32>,
    reason: Vec<Reason>,
    trail: Vec<Node>,
    trail_lim: Vec<usize>,
    qhead: usize,
    // clauses
    clauses: Vec<Vec<Lit>>,
    learnt: Vec<u32>,
    /// per literal: (clause, blocker literal); a true blocker skips the clause (Chu, Stuckey and
    /// Harwood 2008; Glucose 3 / MiniSat 2.2 `Watcher`)
    watches: Vec<Vec<(u32, Lit)>>,
    // order
    activity: Vec<f64>,
    bump: f64,
    heap: Vec<u32>,
    heap_pos: Vec<i32>,
    phase: Vec<i8>,
    seen: Vec<u8>,
    pub conflicts: u64,
    pub decisions: u64,
    pub restarts: u64,
    pub learnt_lits: u64,
    pub minimized_lits: u64,
    pub num_entries: u64,
    pub st_rows: u64,
    pub st_row_entries: u64,
    pub st_watch: u64,
    pub st_analyze_ns: u64,
    pub st_prop_ns: u64,
}

pub enum Outcome {
    Solution(Vec<f64>),
    Exhausted,
    Stopped,
}

impl UipEngine {
    /// `rows`: (entries, rhs) meaning `sum a x <= rhs`. A column is binary when it is integer
    /// with bounds in {0,1}; every other column needs finite bounds. `phase_up`: initial saved
    /// phase. Returns None when a non-binary column has an infinite bound.
    pub fn new(
        n: usize,
        rows: &[(Vec<(u32, f64)>, f64)],
        lb: &[f64],
        ub: &[f64],
        is_int: &[bool],
        phase_up: &[bool],
    ) -> Option<Self> {
        let is_num: Vec<bool> = (0..n)
            .map(|j| {
                !(is_int[j]
                    && (lb[j] == 0.0 || lb[j] == 1.0)
                    && (ub[j] == 0.0 || ub[j] == 1.0)
                    && lb[j] <= ub[j])
            })
            .collect();
        for j in 0..n {
            if is_num[j] && !(lb[j].is_finite() && ub[j].is_finite()) {
                return None;
            }
        }
        let mut cstart = vec![0usize];
        let (mut cvar, mut ccoef) = (Vec::new(), Vec::new());
        let (mut crhs, mut cmaxabs, mut minact) = (Vec::new(), Vec::new(), Vec::new());
        let mut cnt = vec![0usize; n + 1];
        for (ent, rhs) in rows {
            let mut m = 0.0f64;
            let mut mn = 0.0;
            for &(j, a) in ent {
                let ju = j as usize;
                cvar.push(j);
                ccoef.push(a);
                let range = if is_num[ju] { ub[ju] - lb[ju] } else { 1.0 };
                m = m.max(a.abs() * range);
                if is_num[ju] {
                    mn += if a > 0.0 { a * lb[ju] } else { a * ub[ju] };
                } else if a < 0.0 {
                    mn += a;
                }
                cnt[ju + 1] += 1;
            }
            cstart.push(cvar.len());
            crhs.push(*rhs);
            cmaxabs.push(m);
            minact.push(mn);
        }
        for j in 0..n {
            cnt[j + 1] += cnt[j];
        }
        let ostart = cnt.clone();
        let mut fill = cnt;
        let mut ocon = vec![0u32; cvar.len()];
        let mut ocoef = vec![0.0; cvar.len()];
        for c in 0..crhs.len() {
            for k in cstart[c]..cstart[c + 1] {
                let j = cvar[k] as usize;
                ocon[fill[j]] = c as u32;
                ocoef[fill[j]] = ccoef[k];
                fill[j] += 1;
            }
        }
        let nc = crhs.len();
        let mut e = UipEngine {
            n,
            is_num,
            is_int: is_int.to_vec(),
            lb0: lb.to_vec(),
            ub0: ub.to_vec(),
            nlb: lb.to_vec(),
            nub: ub.to_vec(),
            entries: Vec::new(),
            last_entry: vec![[NONE; 2]; n],
            eseen: Vec::new(),
            exact: false,
            exact_count: 0,
            cstart,
            cvar,
            ccoef,
            crhs,
            cmaxabs,
            cbase: minact.clone(),
            minact,
            queued: vec![false; nc],
            dirty: Vec::new(),
            ostart,
            ocon,
            ocoef,
            assign: vec![-1; n],
            level: vec![0; n],
            pos: vec![0; n],
            reason: vec![Reason::Decision; n],
            trail: Vec::new(),
            trail_lim: Vec::new(),
            qhead: 0,
            clauses: Vec::new(),
            learnt: Vec::new(),
            watches: vec![Vec::new(); 2 * n],
            activity: vec![0.0; n],
            bump: 1.0,
            heap: Vec::new(),
            heap_pos: vec![-1; n],
            phase: phase_up.iter().map(|&p| p as i8).collect(),
            seen: vec![0; n],
            conflicts: 0,
            decisions: 0,
            restarts: 0,
            learnt_lits: 0,
            minimized_lits: 0,
            num_entries: 0,
            st_rows: 0,
            st_row_entries: 0,
            st_watch: 0,
            st_analyze_ns: 0,
            st_prop_ns: 0,
        };
        for c in 0..nc as u32 {
            e.mark(c);
        }
        for j in 0..n {
            if !e.is_num[j] {
                e.heap_insert(j as u32);
            }
        }
        for j in 0..n {
            if !e.is_num[j] && lb[j] == ub[j] {
                e.enqueue((2 * j) as u32 + lb[j] as u32, Reason::Decision);
            }
        }
        Some(e)
    }

    fn mark(&mut self, c: u32) {
        if !self.queued[c as usize] {
            self.queued[c as usize] = true;
            self.dirty.push(c);
        }
    }

    fn lit_true(&self, l: Lit) -> bool {
        self.assign[(l >> 1) as usize] == (l & 1) as i8
    }
    fn lit_false(&self, l: Lit) -> bool {
        let a = self.assign[(l >> 1) as usize];
        a >= 0 && a != (l & 1) as i8
    }
    fn dl(&self) -> u32 {
        self.trail_lim.len() as u32
    }

    fn enqueue(&mut self, l: Lit, r: Reason) {
        let v = (l >> 1) as usize;
        let val = (l & 1) as i8;
        self.assign[v] = val;
        self.level[v] = self.dl();
        self.pos[v] = self.trail.len() as u32;
        self.reason[v] = r;
        self.trail.push(l);
        for k in self.ostart[v]..self.ostart[v + 1] {
            let a = self.ocoef[k];
            let c = self.ocon[k] as usize;
            if val == 1 && a > 0.0 {
                self.minact[c] += a;
            } else if val == 0 && a < 0.0 {
                self.minact[c] -= a;
            } else {
                continue;
            }
            self.mark(c as u32);
        }
    }

    /// Tighten a numeric bound (entry on the trail, reason = row `row`).
    fn set_bound(&mut self, v: usize, upper: bool, new: f64, row: u32) {
        let old = if upper { self.nub[v] } else { self.nlb[v] };
        let side = upper as usize;
        let e = Entry {
            var: v as u32,
            upper,
            old,
            new,
            row,
            pos: self.trail.len() as u32,
            level: self.dl(),
            prev: self.last_entry[v][side],
        };
        self.last_entry[v][side] = self.entries.len() as u32;
        self.trail.push(HIGH | self.entries.len() as u32);
        self.entries.push(e);
        self.eseen.push(0);
        self.num_entries += 1;
        if self.exact {
            self.exact_count += 1;
        }
        if upper {
            self.nub[v] = new;
        } else {
            self.nlb[v] = new;
        }
        for k in self.ostart[v]..self.ostart[v + 1] {
            let a = self.ocoef[k];
            let c = self.ocon[k] as usize;
            if (upper && a < 0.0) || (!upper && a > 0.0) {
                self.minact[c] += a * (new - old);
                self.mark(c as u32);
            }
        }
    }

    fn cancel_until(&mut self, lvl: u32) {
        if self.dl() <= lvl {
            return;
        }
        self.exact = false;
        self.exact_count = 0;
        let lim = self.trail_lim[lvl as usize];
        for i in (lim..self.trail.len()).rev() {
            let t = self.trail[i];
            if t & HIGH != 0 {
                let e = self.entries.pop().unwrap();
                self.eseen.pop();
                let v = e.var as usize;
                self.last_entry[v][e.upper as usize] = e.prev;
                if e.upper {
                    self.nub[v] = e.old;
                } else {
                    self.nlb[v] = e.old;
                }
                for k in self.ostart[v]..self.ostart[v + 1] {
                    let a = self.ocoef[k];
                    let c = self.ocon[k] as usize;
                    if (e.upper && a < 0.0) || (!e.upper && a > 0.0) {
                        self.minact[c] -= a * (e.new - e.old);
                    }
                }
                continue;
            }
            let l = t;
            let v = (l >> 1) as usize;
            let val = (l & 1) as i8;
            for k in self.ostart[v]..self.ostart[v + 1] {
                let a = self.ocoef[k];
                let c = self.ocon[k] as usize;
                if val == 1 && a > 0.0 {
                    self.minact[c] -= a;
                } else if val == 0 && a < 0.0 {
                    self.minact[c] += a;
                }
            }
            self.assign[v] = -1;
            self.phase[v] = val;
            self.heap_insert(v as u32);
        }
        self.trail.truncate(lim);
        self.trail_lim.truncate(lvl as usize);
        self.qhead = lim;
        for c in std::mem::take(&mut self.dirty) {
            self.queued[c as usize] = false;
        }
    }

    // ---- heap (max activity) ----
    fn heap_insert(&mut self, v: u32) {
        if self.heap_pos[v as usize] >= 0 {
            return;
        }
        self.heap_pos[v as usize] = self.heap.len() as i32;
        self.heap.push(v);
        self.sift_up(self.heap.len() - 1);
    }
    fn sift_up(&mut self, mut i: usize) {
        let x = self.heap[i];
        while i > 0 {
            let p = (i - 1) / 2;
            if self.activity[x as usize] <= self.activity[self.heap[p] as usize] {
                break;
            }
            self.heap[i] = self.heap[p];
            self.heap_pos[self.heap[i] as usize] = i as i32;
            i = p;
        }
        self.heap[i] = x;
        self.heap_pos[x as usize] = i as i32;
    }
    fn sift_down(&mut self, mut i: usize) {
        let x = self.heap[i];
        let len = self.heap.len();
        loop {
            let mut c = 2 * i + 1;
            if c >= len {
                break;
            }
            if c + 1 < len
                && self.activity[self.heap[c + 1] as usize] > self.activity[self.heap[c] as usize]
            {
                c += 1;
            }
            if self.activity[self.heap[c] as usize] <= self.activity[x as usize] {
                break;
            }
            self.heap[i] = self.heap[c];
            self.heap_pos[self.heap[i] as usize] = i as i32;
            i = c;
        }
        self.heap[i] = x;
        self.heap_pos[x as usize] = i as i32;
    }
    fn heap_pop(&mut self) -> Option<u32> {
        if self.heap.is_empty() {
            return None;
        }
        let top = self.heap[0];
        let last = self.heap.pop().unwrap();
        self.heap_pos[top as usize] = -1;
        if !self.heap.is_empty() {
            self.heap[0] = last;
            self.heap_pos[last as usize] = 0;
            self.sift_down(0);
        }
        Some(top)
    }
    fn bump_var(&mut self, v: usize) {
        self.activity[v] += self.bump;
        if self.activity[v] > 1e100 {
            for a in self.activity.iter_mut() {
                *a *= 1e-100;
            }
            self.bump *= 1e-100;
        }
        if self.heap_pos[v] >= 0 {
            self.sift_up(self.heap_pos[v] as usize);
        }
    }

    // ---- propagation ----
    fn propagate(&mut self) -> Option<Conflict> {
        loop {
            if self.qhead < self.trail.len() {
                let t = self.trail[self.qhead];
                self.qhead += 1;
                if t & HIGH != 0 {
                    continue;
                }
                // clauses watching the literal that just became false
                let fl = t ^ 1;
                let mut ws = std::mem::take(&mut self.watches[fl as usize]);
                let mut i = 0;
                let mut j = 0;
                let mut confl = None;
                while i < ws.len() {
                    let (ci, blocker) = ws[i];
                    i += 1;
                    self.st_watch += 1;
                    if self.lit_true(blocker) {
                        ws[j] = (ci, blocker);
                        j += 1;
                        continue;
                    }
                    {
                        let c = &mut self.clauses[ci as usize];
                        if c.is_empty() {
                            continue; // deleted
                        }
                        if c[0] == fl {
                            c.swap(0, 1);
                        }
                    }
                    let first = self.clauses[ci as usize][0];
                    if first != blocker && self.lit_true(first) {
                        ws[j] = (ci, first);
                        j += 1;
                        continue;
                    }
                    let len = self.clauses[ci as usize].len();
                    let mut found = false;
                    for k in 2..len {
                        let lk = self.clauses[ci as usize][k];
                        if !self.lit_false(lk) {
                            self.clauses[ci as usize].swap(1, k);
                            self.watches[lk as usize].push((ci, first));
                            found = true;
                            break;
                        }
                    }
                    if found {
                        continue;
                    }
                    ws[j] = (ci, first);
                    j += 1;
                    if self.lit_false(first) {
                        confl = Some(Conflict::Clause(ci));
                        while i < ws.len() {
                            ws[j] = ws[i];
                            j += 1;
                            i += 1;
                        }
                    } else {
                        self.enqueue(first, Reason::Clause(ci));
                    }
                }
                ws.truncate(j);
                let added = std::mem::take(&mut self.watches[fl as usize]);
                ws.extend(added);
                self.watches[fl as usize] = ws;
                if confl.is_some() {
                    return confl;
                }
                continue;
            }
            let c = self.dirty.pop()?;
            let cu = c as usize;
            self.queued[cu] = false;
            let slack = self.crhs[cu] - self.minact[cu];
            if slack < -EPS {
                return Some(Conflict::Row(c));
            }
            if slack + EPS >= self.cmaxabs[cu] {
                continue;
            }
            if self.exact && self.exact_count > 200_000 {
                self.dirty.clear();
                return None;
            }
            self.st_rows += 1;
            self.st_row_entries += (self.cstart[cu + 1] - self.cstart[cu]) as u64;
            for k in self.cstart[cu]..self.cstart[cu + 1] {
                let v = self.cvar[k] as usize;
                let a = self.ccoef[k];
                if self.is_num[v] {
                    let (lb, ub) = (self.nlb[v], self.nub[v]);
                    if lb >= ub {
                        continue;
                    }
                    let slack = self.crhs[cu] - self.minact[cu];
                    let width = ub - lb;
                    let thr = if self.exact {
                        1e-9
                    } else if self.is_int[v] {
                        0.5
                    } else {
                        (0.05 * width).max(1e-6)
                    };
                    if a > 0.0 {
                        let mut nu = lb + slack / a;
                        if self.is_int[v] {
                            nu = (nu + 1e-6).floor();
                        }
                        if nu < lb - EPS {
                            return Some(Conflict::Round(c));
                        }
                        if nu < ub - thr {
                            self.set_bound(v, true, nu.max(lb), c);
                        }
                    } else if a < 0.0 {
                        let mut nl = ub + slack / a;
                        if self.is_int[v] {
                            nl = (nl - 1e-6).ceil();
                        }
                        if nl > ub + EPS {
                            return Some(Conflict::Round(c));
                        }
                        if nl > lb + thr {
                            self.set_bound(v, false, nl.min(ub), c);
                        }
                    }
                    continue;
                }
                if self.assign[v] >= 0 {
                    continue;
                }
                let slack = self.crhs[cu] - self.minact[cu];
                if a > 0.0 && a > slack + EPS {
                    self.enqueue((2 * v) as u32, Reason::Row(c));
                } else if a < 0.0 && -a > slack + EPS {
                    self.enqueue((2 * v) as u32 + 1, Reason::Row(c));
                }
            }
        }
    }

    fn node_level(&self, nd: Node) -> u32 {
        if nd & HIGH != 0 {
            self.entries[(nd & !HIGH) as usize].level
        } else {
            self.level[(nd >> 1) as usize]
        }
    }

    /// Nodes that explain an implication (`v` = the implied column, None for a conflict) of row
    /// `c`, among the trail items before position `before`. A binary implication is explained by
    /// the earliest items that already push the minimal activity past the side (Achterberg 2007);
    /// numeric implications and conflicts by every item.
    fn row_reason(&self, c: u32, before: u32, v: Option<usize>, all: bool, out: &mut Vec<Node>) {
        let cu = c as usize;
        let mut cand: Vec<(u32, Node, f64)> = Vec::new();
        let mut implied_abs = 0.0;
        let all = all || v.is_none_or(|u| self.is_num[u]);
        for k in self.cstart[cu]..self.cstart[cu + 1] {
            let u = self.cvar[k] as usize;
            let a = self.ccoef[k];
            if self.is_num[u] {
                let side = (a < 0.0) as usize;
                let mut e = self.last_entry[u][side];
                while e != NONE && self.entries[e as usize].pos >= before {
                    e = self.entries[e as usize].prev;
                }
                if e == NONE {
                    continue;
                }
                let en = &self.entries[e as usize];
                let root = if side == 1 { self.ub0[u] } else { self.lb0[u] };
                cand.push((en.pos, HIGH | e, a.abs() * (en.new - root).abs()));
                continue;
            }
            if Some(u) == v {
                implied_abs = a.abs();
                continue;
            }
            if self.assign[u] < 0 || self.pos[u] >= before {
                continue;
            }
            let val = self.assign[u];
            if (val == 1 && a > 0.0) || (val == 0 && a < 0.0) {
                cand.push((self.pos[u], (2 * u) as u32 + val as u32, a.abs()));
            }
        }
        cand.sort_unstable_by_key(|x| x.0);
        let thr = self.crhs[cu] - self.cbase[cu] - implied_abs + EPS;
        let mut sum = 0.0;
        for (_, l, w) in cand {
            out.push(l);
            sum += w;
            if !all && sum > thr {
                break;
            }
        }
    }

    fn reason_nodes(&self, nd: Node, out: &mut Vec<Node>) {
        if nd & HIGH != 0 {
            let e = &self.entries[(nd & !HIGH) as usize];
            self.row_reason(e.row, e.pos, Some(e.var as usize), false, out);
            return;
        }
        let v = (nd >> 1) as usize;
        match self.reason[v] {
            Reason::Decision => {}
            Reason::Clause(ci) => {
                for &l in &self.clauses[ci as usize] {
                    if (l >> 1) as usize != v {
                        out.push(l ^ 1);
                    }
                }
            }
            Reason::Row(c) => self.row_reason(c, self.pos[v], Some(v), false, out),
        }
    }

    fn is_decision(&self, nd: Node) -> bool {
        nd & HIGH == 0 && self.reason[(nd >> 1) as usize] == Reason::Decision
    }

    /// True if every binary the node depends on is in the clause (seen 1), at level 0, or itself
    /// redundant.
    fn redundant(
        &mut self,
        nd: Node,
        depth: u32,
        cb: &mut Vec<usize>,
        ce: &mut Vec<usize>,
    ) -> bool {
        if depth > 48 {
            return false;
        }
        let mut nodes = Vec::new();
        self.reason_nodes(nd, &mut nodes);
        for &t in &nodes {
            if self.node_level(t) == 0 {
                continue;
            }
            let st = if t & HIGH != 0 {
                self.eseen[(t & !HIGH) as usize]
            } else {
                self.seen[(t >> 1) as usize]
            };
            if st == 1 || st == 2 {
                continue;
            }
            if st == 3 || self.is_decision(t) {
                self.mark_seen(t, 3, cb, ce);
                return false;
            }
            if self.redundant(t, depth + 1, cb, ce) {
                self.mark_seen(t, 2, cb, ce);
            } else {
                self.mark_seen(t, 3, cb, ce);
                return false;
            }
        }
        true
    }

    fn mark_seen(&mut self, t: Node, st: u8, cb: &mut Vec<usize>, ce: &mut Vec<usize>) {
        if t & HIGH != 0 {
            let e = (t & !HIGH) as usize;
            if self.eseen[e] == 0 {
                ce.push(e);
            }
            self.eseen[e] = st;
        } else {
            let v = (t >> 1) as usize;
            if self.seen[v] == 0 {
                cb.push(v);
            }
            self.seen[v] = st;
        }
    }

    /// First-UIP analysis. Returns (learnt clause with the asserting literal first, backjump
    /// level), or None when the conflict does not depend on any decision.
    fn analyze(&mut self, confl: Conflict) -> Option<(Vec<Lit>, u32)> {
        let mut learnt: Vec<Lit> = vec![0];
        let mut path = 0i32;
        let mut idx = self.trail.len();
        let mut work: Vec<Node> = Vec::new();
        match confl {
            Conflict::Clause(ci) => {
                for &l in &self.clauses[ci as usize] {
                    work.push(l ^ 1);
                }
            }
            Conflict::Row(c) => self.row_reason(c, u32::MAX, None, true, &mut work),
            Conflict::Round(c) => self.row_reason(c, u32::MAX, None, true, &mut work),
        }
        let cur = self.dl();
        let mut uip: Option<Lit> = None;
        let mut tb: Vec<usize> = Vec::new();
        let mut te: Vec<usize> = Vec::new();
        loop {
            while let Some(nd) = work.pop() {
                if nd & HIGH != 0 {
                    let e = (nd & !HIGH) as usize;
                    if self.eseen[e] != 0 || self.entries[e].level == 0 {
                        continue;
                    }
                    self.eseen[e] = 1;
                    te.push(e);
                    if self.entries[e].level == cur {
                        path += 1;
                    } else {
                        self.reason_nodes(nd, &mut work);
                    }
                } else {
                    let v = (nd >> 1) as usize;
                    if self.seen[v] == 0 && self.level[v] > 0 {
                        self.bump_var(v);
                        self.seen[v] = 1;
                        tb.push(v);
                        if self.level[v] == cur {
                            path += 1;
                        } else {
                            learnt.push(nd ^ 1);
                        }
                    }
                }
            }
            if path <= 0 {
                break;
            }
            let item = loop {
                idx -= 1;
                let t = self.trail[idx];
                let s = if t & HIGH != 0 {
                    self.eseen[(t & !HIGH) as usize]
                } else {
                    self.seen[(t >> 1) as usize]
                };
                if s == 1 {
                    break t;
                }
            };
            if path == 1 && item & HIGH == 0 {
                self.seen[(item >> 1) as usize] = 0;
                uip = Some(item);
                break;
            }
            if item & HIGH != 0 {
                self.eseen[(item & !HIGH) as usize] = 0;
            } else {
                self.seen[(item >> 1) as usize] = 0;
            }
            path -= 1;
            self.reason_nodes(item, &mut work);
        }
        let result = (|| {
            match uip {
                Some(p) => learnt[0] = p ^ 1,
                None => {
                    if learnt.len() == 1 {
                        return None;
                    }
                    // no literal at the current level: the highest-level one asserts
                    let mut at = 1;
                    for i in 2..learnt.len() {
                        if self.level[(learnt[i] >> 1) as usize]
                            > self.level[(learnt[at] >> 1) as usize]
                        {
                            at = i;
                        }
                    }
                    learnt.swap(0, at);
                    let last = learnt.len() - 1;
                    learnt.swap(at, last);
                    learnt.pop();
                    // the placeholder slot was replaced by the old last literal
                }
            }
            Some(())
        })();
        if result.is_none() {
            for v in tb {
                self.seen[v] = 0;
            }
            for e in te {
                self.eseen[e] = 0;
            }
            return None;
        }
        // recursive minimisation of learnt[1..]
        let before = learnt.len();
        let (mut cb, mut ce): (Vec<usize>, Vec<usize>) = (Vec::new(), Vec::new());
        let mut keep = vec![learnt[0]];
        for i in 1..learnt.len() {
            let nd = learnt[i] ^ 1;
            if self.is_decision(nd) || !self.redundant(nd, 0, &mut cb, &mut ce) {
                keep.push(learnt[i]);
            }
        }
        for v in cb.into_iter().chain(tb) {
            self.seen[v] = 0;
        }
        for e in ce.into_iter().chain(te) {
            self.eseen[e] = 0;
        }
        self.learnt_lits += before as u64;
        self.minimized_lits += keep.len() as u64;
        let mut learnt = keep;
        let mut bt = 0;
        let mut at = 1;
        for i in 1..learnt.len() {
            let l = self.level[(learnt[i] >> 1) as usize];
            if l > bt {
                bt = l;
                at = i;
            }
        }
        if learnt.len() > 1 {
            learnt.swap(1, at);
        }
        Some((learnt, bt))
    }

    fn add_learnt(&mut self, c: Vec<Lit>) -> u32 {
        let ci = self.clauses.len() as u32;
        self.watches[c[0] as usize].push((ci, c[1]));
        self.watches[c[1] as usize].push((ci, c[0]));
        self.clauses.push(c);
        self.learnt.push(ci);
        ci
    }

    fn luby(i: u64) -> u64 {
        let mut k = 1u32;
        while (1u64 << k) - 1 < i {
            k += 1;
        }
        let mut i = i;
        loop {
            if i == (1u64 << k) - 1 {
                return 1u64 << (k - 1);
            }
            i -= (1u64 << (k - 1)) - 1;
            k = 1;
            while (1u64 << k) - 1 < i {
                k += 1;
            }
        }
    }

    fn reduce_db(&mut self) {
        // at level 0 only: drop the older half of the long learnt clauses
        let mut long: Vec<u32> = self
            .learnt
            .iter()
            .copied()
            .filter(|&c| self.clauses[c as usize].len() > 2)
            .collect();
        long.sort_unstable();
        let drop = long.len() / 2;
        for &c in &long[..drop] {
            self.clauses[c as usize] = Vec::new();
        }
        let clauses = &self.clauses;
        self.learnt.retain(|&c| !clauses[c as usize].is_empty());
        for w in &mut self.watches {
            w.retain(|&(c, _)| !clauses[c as usize].is_empty());
        }
    }

    /// Search until a solution, exhaustion, `stop()` or `max_conflicts`.
    pub fn solve(&mut self, max_conflicts: u64, stop: &mut dyn FnMut() -> bool) -> Outcome {
        let mut next_restart = 100 * Self::luby(1);
        let mut since = 0u64;
        let mut max_learnts = 20000usize.max(self.n);
        self.restarts = 1;
        loop {
            let tp = std::time::Instant::now();
            let pr = self.propagate();
            self.st_prop_ns += tp.elapsed().as_nanos() as u64;
            if let Some(confl) = pr {
                self.conflicts += 1;
                since += 1;
                if self.dl() == 0 {
                    return Outcome::Exhausted;
                }
                let ta = std::time::Instant::now();
                let an = self.analyze(confl);
                self.st_analyze_ns += ta.elapsed().as_nanos() as u64;
                let Some((learnt, bt)) = an else {
                    return Outcome::Exhausted;
                };
                self.bump *= 1.0 / 0.95;
                self.cancel_until(bt);
                if learnt.len() == 1 {
                    self.enqueue(learnt[0], Reason::Decision);
                } else {
                    let l0 = learnt[0];
                    let ci = self.add_learnt(learnt);
                    self.enqueue(l0, Reason::Clause(ci));
                }
                if self.conflicts >= max_conflicts || (self.conflicts & 255 == 0 && stop()) {
                    return Outcome::Stopped;
                }
                continue;
            }
            if self.exact && self.exact_count > 200_000 {
                // the tight numeric fixpoint did not converge: count it as a failed leaf
                self.conflicts += 1;
                since = next_restart;
                self.cancel_until(0);
            }
            if since >= next_restart {
                self.restarts += 1;
                since = 0;
                next_restart = 100 * Self::luby(self.restarts);
                self.cancel_until(0);
                if self.learnt.len() > max_learnts {
                    self.reduce_db();
                    max_learnts += max_learnts / 10;
                }
                continue;
            }
            self.decisions += 1;
            if self.decisions & 1023 == 0 && stop() {
                return Outcome::Stopped;
            }
            let mut next = None;
            while let Some(v) = self.heap_pop() {
                if self.assign[v as usize] < 0 {
                    next = Some(v);
                    break;
                }
            }
            let Some(v) = next else {
                if !self.exact && self.is_num.iter().any(|&b| b) {
                    // all binaries set: tighten the numeric bounds to a fixpoint
                    self.exact = true;
                    self.exact_count = 0;
                    for c in 0..self.crhs.len() as u32 {
                        self.mark(c);
                    }
                    continue;
                }
                let sol = (0..self.n)
                    .map(|j| {
                        if self.is_num[j] {
                            self.nlb[j]
                        } else {
                            self.assign[j] as f64
                        }
                    })
                    .collect();
                return Outcome::Solution(sol);
            };
            self.trail_lim.push(self.trail.len());
            let ph = self.phase[v as usize] as u32;
            self.enqueue(2 * v + ph, Reason::Decision);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(n: usize, rows: Vec<(Vec<(u32, f64)>, f64)>) -> Outcome {
        let mut e = UipEngine::new(
            n,
            &rows,
            &vec![0.0; n],
            &vec![1.0; n],
            &vec![true; n],
            &vec![false; n],
        )
        .unwrap();
        e.solve(1_000_000, &mut || false)
    }

    #[test]
    fn pigeonhole_is_exhausted() {
        // 5 pigeons, 4 holes: x[p][h]; each pigeon in a hole, each hole at most one pigeon
        let (p, h) = (5usize, 4usize);
        let id = |a: usize, b: usize| (a * h + b) as u32;
        let mut rows = Vec::new();
        for a in 0..p {
            rows.push(((0..h).map(|b| (id(a, b), -1.0)).collect(), -1.0)); // sum >= 1
        }
        for b in 0..h {
            rows.push(((0..p).map(|a| (id(a, b), 1.0)).collect(), 1.0));
        }
        assert!(matches!(run(p * h, rows), Outcome::Exhausted));
    }

    #[test]
    fn assignment_is_solved_and_feasible() {
        // 6 pigeons, 6 holes with a knapsack-ish side row
        let (p, h) = (6usize, 6usize);
        let id = |a: usize, b: usize| (a * h + b) as u32;
        let mut rows = Vec::new();
        for a in 0..p {
            rows.push(((0..h).map(|b| (id(a, b), -1.0)).collect(), -1.0));
            rows.push(((0..h).map(|b| (id(a, b), 1.0)).collect(), 1.0));
        }
        for b in 0..h {
            rows.push(((0..p).map(|a| (id(a, b), 1.0)).collect(), 1.0));
        }
        match run(p * h, rows) {
            Outcome::Solution(s) => {
                for a in 0..p {
                    assert_eq!((0..h).map(|b| s[a * h + b]).sum::<f64>(), 1.0);
                }
            }
            _ => panic!("expected a solution"),
        }
    }

    #[test]
    fn numeric_columns_propagate_and_complete() {
        // binary y: s2 - s1 >= 3 - 10 (1 - y)... a disjunctive pair: s1 + 3 <= s2 or s2 + 3 <= s1
        // (big-M 10) over s1, s2 in [0, 5], with y fixing the order. Row form `<=`.
        let (y, s1, s2) = (0u32, 1u32, 2u32);
        let rows = vec![
            (vec![(s1, 1.0), (s2, -1.0), (y, 10.0)], 7.0), // s1 - s2 + 3 <= 10 (1 - y)
            (vec![(s2, 1.0), (s1, -1.0), (y, -10.0)], -3.0), // s2 - s1 + 3 <= 10 y
        ];
        let mut e = UipEngine::new(
            3,
            &rows,
            &[0.0, 0.0, 0.0],
            &[1.0, 5.0, 5.0],
            &[true, false, false],
            &[true, false, false],
        )
        .unwrap();
        match e.solve(1000, &mut || false) {
            Outcome::Solution(s) => {
                let (y, a, b) = (s[0], s[1], s[2]);
                assert!(a - b + 3.0 <= 10.0 * (1.0 - y) + 1e-6);
                assert!(b - a + 3.0 <= 10.0 * y + 1e-6);
            }
            _ => panic!("expected a solution"),
        }
    }
}
