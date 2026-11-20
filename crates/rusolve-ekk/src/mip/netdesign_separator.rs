//! Network design separator and primal heuristic (rusolve levers `mip_netdesign_sepa` and
//! `mip_netdesign_heur`, default off).
//!
//! Not a HiGHS file. Provenance (literature, no solver source was copied):
//! * Gendron, Crainic, Frangioni, "Multicommodity capacitated network design" (1999): the strong
//!   linking inequalities `x^k_a <= min(d^k, u_a) y_a` and their gap to the weak capacity rows;
//! * Chouman, Crainic, Gendron, "Commodity representations and cut-set-based inequalities for
//!   multicommodity capacitated fixed-charge network design", Transportation Science 51 (2017);
//! * Crainic, Gendron, Hernu, "A slope scaling/Lagrangean perturbation heuristic with long-term
//!   memory for multicommodity capacitated fixed-charge network design", J. Heuristics 10 (2004).
//! * Bienstock, Chopra, Guenluek, Tsai, "Minimum cost capacity installation for multicommodity
//!   network flows", Math. Programming 81 (1998); Raack, Koster, Orlowski, Waessaly, "On cut-based
//!   inequalities for capacitated network design polyhedra", Networks 57 (2011); Atamtuerk, "On
//!   capacitated network design cut-set polyhedra", Math. Programming 92 (2002): cut-set
//!   inequalities and their MIR strengthening (lever `mip_netdesign_cutset`).
//!
//! Detected structure (no names are read):
//! * a commodity is a connected set of continuous columns tied together by equality rows with
//!   all coefficients `+-1` (flow conservation); its supply `d` is half the sum of `|rhs|`;
//! * a capacity row is `sum x_k - u y <= 0` with one binary `y` (`u > 0`) and continuous `x_k`
//!   with coefficient `1`.
//!
//! Validity guard (checked per commodity, the strong row is skipped when it fails). An integer
//! feasible point with a flow cycle can be replaced by one without the cycle (cancel the
//! smallest flow on it), and an acyclic flow has `x_a <= d` on every arc. The replacement keeps
//! feasibility and does not raise the cost when, for every column of the commodity:
//! * the lower bound is 0 and the cost is non-negative (minimisation);
//! * the column appears in exactly two equality rows of the commodity, with coefficients +1 and
//!   -1 (a network matrix, so a cycle exists as a circulation);
//! * every other row holding the column is one-sided with the sign that a smaller `x` keeps
//!   satisfied (`<=` with a positive, `>=` with a negative coefficient).
//!
//! Then `x <= min(d, u) y` removes no optimal solution (cost of `y` is untouched by the swap).
//!
//! Heuristic: slope scaling. The LP with the fixed costs replaced by a linear cost per unit of
//! flow, rescaled each pass by the current arc flow (`f_a / flow_a`), with the design `y = (flow > 0)`
//! read off every pass; the best design is then re-solved with the true costs.

use crate::prelude::*;

/// Cuts added per call.
const MAX_CUTS_PER_CALL: usize = 3_000;
/// Minimum violation of a cut to add it.
const MIN_VIOLATION: f64 = 1e-3;
/// Cut-set cuts added per call.
const MAX_CUTSET_CUTS: usize = 300;

/// Slope scaling passes.
const SLOPE_PASSES: usize = 25;
/// Share of the time limit the heuristic may use.
const HEUR_TIME_SHARE: f64 = 0.2;
/// Warm-start slope scaling from the root basis only once the root LP took this many iterations:
/// below it a cold pass is cheap, and the cold optimum is the one the canad family was tuned on
/// (canad_c46/c50 incumbents are worse from the warm optimum).
const WARM_MIN_ROOT_ITERS: i64 = 20_000;

/// One capacity row: binary `y`, capacity `u`, fixed cost `f`, flow columns.
#[derive(Clone, Debug, Default)]
struct NdArc {
    y: HighsInt,
    u: f64,
    xs: Vec<HighsInt>,
}

#[derive(Clone, Default, Debug)]
pub struct HighsNetdesignSeparator {
    pub base: HighsSeparator,
    detected: bool,
    heur_done: bool,
    /// (flow column, binary column, bound `min(d, u)`)
    links: Vec<(HighsInt, HighsInt, f64)>,
    arcs: Vec<NdArc>,
    cs: CutsetGraph,
    /// solutions built by the heuristic, taken by `HighsSeparation::take_pending_solutions`
    /// (per separator, so a worker never writes shared solver state)
    pub pending: Vec<(u64, Vec<f64>)>,
}

/// Node-and-arc view of the detected network for cut-set separation.
#[derive(Clone, Default, Debug)]
struct CutsetGraph {
    /// (binary column, capacity, tail node, head node); the tail row has coefficient +1
    arcs: Vec<(HighsInt, f64, usize, usize)>,
    /// per node: (commodity, summed rhs of its conservation rows)
    node_rhs: Vec<Vec<(usize, f64)>>,
    out_adj: Vec<Vec<usize>>,
    in_adj: Vec<Vec<usize>>,
    ncomm: usize,
}

/// Cut-set inequality sum_{a in cut} u_a y_a >= D, strengthened by MIR (Bienstock, Chopra,
/// Guenluek, Tsai 1998; Atamtuerk 2002). `caps` are the arc capacities, `y` the LP values.
/// Returns (coefficients g with `sum g_j y_j <= rhs`, violation at `y`) for the best divisor.
fn cutset_mir(caps: &[f64], y: &[f64], d: f64) -> Option<(Vec<f64>, f64, f64)> {
    let trunc: Vec<f64> = caps.iter().map(|&u| u.min(d)).collect();
    let mut divisors: Vec<f64> = trunc.clone();
    divisors.sort_by(|a, b| a.partial_cmp(b).unwrap());
    divisors.dedup_by(|a, b| (*a - *b).abs() < 1e-9 * b.abs().max(1.0));
    if divisors.len() > 8 {
        let step = divisors.len() as f64 / 8.0;
        divisors = (0..8)
            .map(|i| divisors[(i as f64 * step) as usize])
            .collect();
    }
    let mut best: Option<(Vec<f64>, f64, f64)> = None;
    for &dv in &divisors {
        if dv <= 1e-9 {
            continue;
        }
        let beta = -d / dv;
        let fl = (beta + 1e-9).floor();
        let mut f = beta - fl;
        if f < 1e-6 {
            f = 0.0;
        }
        let g: Vec<f64> = trunc
            .iter()
            .map(|&u| {
                let alpha = -u / dv;
                let af = (alpha + 1e-9).floor();
                let fj = (alpha - af).max(0.0);
                af + if fj > f { (fj - f) / (1.0 - f) } else { 0.0 }
            })
            .collect();
        let viol: f64 = g.iter().zip(y).map(|(a, b)| a * b).sum::<f64>() - fl;
        // scale so the coefficients are O(1)
        if best.as_ref().is_none_or(|b| viol > b.1) {
            best = Some((g, viol, fl));
        }
    }
    best.filter(|b| b.1 > MIN_VIOLATION)
}

fn find(uf: &mut [usize], x: usize) -> usize {
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

/// Validity guard per commodity (union-find root -> true when every flow column passes).
#[allow(clippy::too_many_arguments)]
fn component_guard(
    uf: &mut [usize],
    in_flow: &[bool],
    lower: &[f64],
    cost: &[f64],
    sgn: f64,
    flow_pos: &[u8],
    flow_neg: &[u8],
    other_bad: &[bool],
) -> Vec<bool> {
    let mut comp_ok = vec![true; in_flow.len()];
    for cc in 0..in_flow.len() {
        if !in_flow[cc] {
            continue;
        }
        let ok = lower[cc] == 0.0
            && sgn * cost[cc] >= 0.0
            && flow_pos[cc] == 1
            && flow_neg[cc] == 1
            && !other_bad[cc];
        if !ok {
            let r = find(uf, cc);
            comp_ok[r] = false;
        }
    }
    comp_ok
}

impl HighsNetdesignSeparator {
    pub fn new(ms: &HighsMipSolver) -> Self {
        HighsNetdesignSeparator {
            base: HighsSeparator::new(ms, "Separation: network design linking"),
            ..Default::default()
        }
    }

    /// Finds the commodities and the capacity rows (once).
    fn detect(&mut self, ms: &HighsMipSolver, lp: &HighsLpRelaxation, g: usize) {
        self.detected = true;
        let nrow = ms.num_row();
        let ncol = ms.num_col() as usize;
        let lower = &ms.mipdata_.domains[g].col_lower_;
        let upper = &ms.mipdata_.domains[g].col_upper_;
        let model = lp.get_lp();
        let is_cont = |col: usize| {
            matches!(
                ms.variable_type_2(col as HighsInt),
                HighsVarType::kContinuous
            )
        };
        let is_bin = |col: usize| !is_cont(col) && lower[col] == 0.0 && upper[col] == 1.0;
        let get = |r: HighsInt| -> (&[HighsInt], &[f64]) {
            let (mut len, mut inds, mut vals): (HighsInt, &[HighsInt], &[f64]) = (0, &[], &[]);
            lp.get_row(ms, r, &mut len, &mut inds, &mut vals);
            (inds, vals)
        };
        let sgn = if model.sense_ == ObjSense::kMaximize {
            -1.0
        } else {
            1.0
        };
        let mut uf: Vec<usize> = (0..ncol).collect();
        let mut in_flow = vec![false; ncol];
        let mut flow_pos = vec![0u8; ncol];
        let mut flow_neg = vec![0u8; ncol];
        let mut other_bad = vec![false; ncol];
        let mut eq_rows: Vec<(usize, f64)> = Vec::new();
        let mut flow_rows: Vec<(usize, f64, usize)> = Vec::new();
        let mut col_prow = vec![usize::MAX; ncol];
        let mut col_nrow = vec![usize::MAX; ncol];
        for r in 0..nrow {
            let (rl, ru) = (model.row_lower_[r as usize], model.row_upper_[r as usize]);
            let (inds, vals) = get(r);
            let is_flow = rl == ru
                && inds.len() >= 2
                && vals.iter().all(|v| v.abs() == 1.0)
                && inds.iter().all(|&cc| is_cont(cc as usize));
            if is_flow {
                let a = inds[0] as usize;
                for (k, &cc) in inds.iter().enumerate() {
                    in_flow[cc as usize] = true;
                    if vals[k] > 0.0 {
                        col_prow[cc as usize] = r as usize;
                        flow_pos[cc as usize] = flow_pos[cc as usize].saturating_add(1);
                    } else {
                        col_nrow[cc as usize] = r as usize;
                        flow_neg[cc as usize] = flow_neg[cc as usize].saturating_add(1);
                    }
                    let (ra, rb) = (find(&mut uf, a), find(&mut uf, cc as usize));
                    if ra != rb {
                        uf[ra] = rb;
                    }
                }
                eq_rows.push((a, rl.abs()));
                flow_rows.push((r as usize, rl, a));
            } else {
                // a smaller flow must keep the row satisfied
                for (k, &cc) in inds.iter().enumerate() {
                    if !is_cont(cc as usize) {
                        continue;
                    }
                    let keeps = if rl <= -kHighsInf && ru < kHighsInf {
                        vals[k] > 0.0
                    } else if ru >= kHighsInf && rl > -kHighsInf {
                        vals[k] < 0.0
                    } else {
                        false
                    };
                    if !keeps {
                        other_bad[cc as usize] = true;
                    }
                }
            }
        }
        let mut supply = vec![0.0f64; ncol];
        for (a, d) in eq_rows {
            let ra = find(&mut uf, a);
            supply[ra] += 0.5 * d;
        }
        let comp_ok = component_guard(
            &mut uf,
            &in_flow,
            &lower[..ncol],
            &model.col_cost_,
            sgn,
            &flow_pos,
            &flow_neg,
            &other_bad,
        );
        let mut skipped = 0usize;
        for r in 0..nrow {
            let (rl, ru) = (model.row_lower_[r as usize], model.row_upper_[r as usize]);
            if rl > -kHighsInf || ru != 0.0 {
                continue;
            }
            let (inds, vals) = get(r);
            if inds.len() < 2 {
                continue;
            }
            let mut y: Option<(HighsInt, f64)> = None;
            let mut ok = true;
            for (k, &cc) in inds.iter().enumerate() {
                if is_bin(cc as usize) && vals[k] < 0.0 && y.is_none() {
                    y = Some((cc, -vals[k]));
                } else if !(is_cont(cc as usize) && vals[k] == 1.0 && in_flow[cc as usize]) {
                    ok = false;
                    break;
                }
            }
            let Some((yc, u)) = y else { continue };
            if !ok {
                continue;
            }
            let mut arc = NdArc {
                y: yc,
                u,
                xs: Vec::new(),
            };
            for (k, &cc) in inds.iter().enumerate() {
                if cc == yc {
                    continue;
                }
                arc.xs.push(cc);
                let root = find(&mut uf, cc as usize);
                if !comp_ok[root] {
                    skipped += 1;
                    continue;
                }
                let b = supply[root].min(upper[cc as usize].min(u));
                if b > 0.0 && b < u * (1.0 - 1e-9) {
                    self.links.push((cc, yc, b));
                }
            }
            self.arcs.push(arc);
        }
        if ms.options_mip_.mip_netdesign_cutset && !self.arcs.is_empty() {
            self.build_cutset_graph(
                nrow as usize,
                &mut uf,
                &flow_rows,
                &col_prow,
                &col_nrow,
                &in_flow,
                &flow_pos,
                &flow_neg,
                &lower[..ncol],
            );
        }
        if !self.arcs.is_empty() {
            crate::highs_log_user!(
                &ms.options_mip_.log_options,
                HighsLogType::kInfo,
                "Network design: {} capacity rows, {} strong linking candidates, {} flow columns skipped by the validity guard\n",
                self.arcs.len(),
                self.links.len(),
                skipped
            );
        }
    }

    /// Builds the node graph. A commodity is used only when all its columns are network columns
    /// inside a capacity row, with lower bound 0 and its conservation rows balance (sum of rhs is 0, which fails when
    /// presolve dropped a row). Ignoring a commodity keeps the cuts valid: flows are
    /// non-negative, so the rest still has to cross the cut. Arcs carrying no used commodity are
    /// left out for the same reason.
    #[allow(clippy::too_many_arguments)]
    fn build_cutset_graph(
        &mut self,
        nrow: usize,
        uf: &mut [usize],
        flow_rows: &[(usize, f64, usize)],
        col_prow: &[usize],
        col_nrow: &[usize],
        in_flow: &[bool],
        flow_pos: &[u8],
        flow_neg: &[u8],
        lower: &[f64],
    ) {
        let ncol = in_flow.len();
        let mut comp_bad = vec![false; ncol];
        // a flow column without a capacity row can cross any cut for free
        let mut covered = vec![false; ncol];
        for arc in &self.arcs {
            for &q in &arc.xs {
                covered[q as usize] = true;
            }
        }
        for cc in 0..ncol {
            if in_flow[cc]
                && (!covered[cc] || flow_pos[cc] != 1 || flow_neg[cc] != 1 || lower[cc] != 0.0)
            {
                let r = find(uf, cc);
                comp_bad[r] = true;
            }
        }
        let mut bal = vec![0.0f64; ncol];
        for &(_, rl, a) in flow_rows {
            let r = find(uf, a);
            bal[r] += rl;
        }
        for r in 0..ncol {
            if bal[r].abs() > 1e-7 {
                comp_bad[r] = true;
            }
        }
        let good = |uf: &mut [usize], q: usize| in_flow[q] && !comp_bad[find(uf, q)];
        // y columns must appear in one capacity row only
        let mut y_seen = std::collections::HashSet::new();
        if !self.arcs.iter().all(|a| y_seen.insert(a.y)) {
            return;
        }
        let mut ur: Vec<usize> = (0..nrow).collect();
        let mut arc_rows: Vec<Option<(usize, usize)>> = Vec::new();
        for arc in &self.arcs {
            let mut first: Option<(usize, usize)> = None;
            for &q in &arc.xs {
                let q = q as usize;
                if !good(uf, q) {
                    continue;
                }
                let (p, n) = (col_prow[q], col_nrow[q]);
                match first {
                    None => first = Some((p, n)),
                    Some((p0, n0)) => {
                        let (a, b) = (find(&mut ur, p0), find(&mut ur, p));
                        ur[a] = b;
                        let (a, b) = (find(&mut ur, n0), find(&mut ur, n));
                        ur[a] = b;
                    }
                }
            }
            arc_rows.push(first);
        }
        let mut node_id = std::collections::HashMap::new();
        let mut comm_id = std::collections::HashMap::new();
        let mut node_rhs: Vec<std::collections::BTreeMap<usize, f64>> = Vec::new();
        for &(r, rl, a) in flow_rows {
            if !good(uf, a) {
                continue;
            }
            let k = {
                let root = find(uf, a);
                let n = comm_id.len();
                *comm_id.entry(root).or_insert(n)
            };
            let root = find(&mut ur, r);
            let n = node_id.len();
            let v = *node_id.entry(root).or_insert(n);
            if v == node_rhs.len() {
                node_rhs.push(Default::default());
            }
            *node_rhs[v].entry(k).or_insert(0.0) += rl;
        }
        let nn = node_rhs.len();
        let mut g = CutsetGraph {
            node_rhs: node_rhs
                .into_iter()
                .map(|m| m.into_iter().filter(|e| e.1 != 0.0).collect())
                .collect(),
            out_adj: vec![Vec::new(); nn],
            in_adj: vec![Vec::new(); nn],
            ncomm: comm_id.len(),
            ..Default::default()
        };
        for (arc, rows) in self.arcs.iter().zip(&arc_rows) {
            let Some((p, n)) = *rows else { continue };
            let (Some(&t), Some(&h)) = (
                node_id.get(&find(&mut ur, p)),
                node_id.get(&find(&mut ur, n)),
            ) else {
                continue;
            };
            if t == h {
                continue;
            }
            let ai = g.arcs.len();
            g.arcs.push((arc.y, arc.u, t, h));
            g.out_adj[t].push(ai);
            g.in_adj[h].push(ai);
        }
        self.cs = g;
    }

    /// Cut-set separation at the LP point. For dir = +1 the arcs leaving S must carry the net
    /// outflow of S of every commodity; dir = -1 handles the complement (arcs entering S).
    fn separate_cutsets(&self, sol: &[f64]) -> Vec<(Vec<HighsInt>, Vec<f64>, f64, f64)> {
        let g = &self.cs;
        let nn = g.node_rhs.len();
        let mut out: Vec<(Vec<HighsInt>, Vec<f64>, f64, f64)> = Vec::new();
        if nn < 2 || g.arcs.is_empty() {
            return out;
        }
        let max_size = (nn / 2).clamp(1, 12);
        let mut seen = std::collections::HashSet::new();
        for dir in [1.0f64, -1.0] {
            for v0 in 0..nn {
                let mut in_s = vec![false; nn];
                let mut f = vec![0.0f64; g.ncomm];
                let mut members: Vec<usize> = Vec::new();
                let mut d = 0.0f64;
                let mut cap = 0.0f64;
                // adds node v and returns nothing; maintains d, cap, f
                let add = |v: usize,
                           in_s: &mut Vec<bool>,
                           f: &mut Vec<f64>,
                           d: &mut f64,
                           cap: &mut f64| {
                    for &(k, r) in &g.node_rhs[v] {
                        *d -= (dir * f[k]).max(0.0);
                        f[k] += r;
                        *d += (dir * f[k]).max(0.0);
                    }
                    let (leave, enter) = if dir > 0.0 {
                        (&g.out_adj[v], &g.in_adj[v])
                    } else {
                        (&g.in_adj[v], &g.out_adj[v])
                    };
                    for &a in enter {
                        let (y, u, t, h) = g.arcs[a];
                        let other = if dir > 0.0 { t } else { h };
                        if in_s[other] {
                            *cap -= u * sol[y as usize];
                        }
                    }
                    for &a in leave {
                        let (y, u, t, h) = g.arcs[a];
                        let other = if dir > 0.0 { h } else { t };
                        if !in_s[other] {
                            *cap += u * sol[y as usize];
                        }
                    }
                    in_s[v] = true;
                };
                add(v0, &mut in_s, &mut f, &mut d, &mut cap);
                members.push(v0);
                let mut best = (cap - d, 1usize);
                while members.len() < max_size {
                    // best neighbour by the score after adding it
                    let mut cand: Option<(f64, usize)> = None;
                    let mut tried = std::collections::HashSet::new();
                    for &m in &members {
                        for &a in g.out_adj[m].iter().chain(&g.in_adj[m]) {
                            let (_, _, t, h) = g.arcs[a];
                            let u = if in_s[t] { h } else { t };
                            if in_s[u] || !tried.insert(u) {
                                continue;
                            }
                            let (mut i2, mut f2, mut d2, mut c2) =
                                (in_s.clone(), f.clone(), d, cap);
                            add(u, &mut i2, &mut f2, &mut d2, &mut c2);
                            let sc = c2 - d2;
                            if cand.is_none_or(|q| sc < q.0 - 1e-12) {
                                cand = Some((sc, u));
                            }
                        }
                    }
                    let Some((_, u)) = cand else { break };
                    add(u, &mut in_s, &mut f, &mut d, &mut cap);
                    members.push(u);
                    if cap - d < best.0 - 1e-12 {
                        best = (cap - d, members.len());
                    }
                }
                let _ = best;
                for plen in 1..=members.len() {
                    // the cut of the prefix S = members[..plen]
                    let mut s_set = vec![false; nn];
                    let mut fk = vec![0.0f64; g.ncomm];
                    for &m in &members[..plen] {
                        s_set[m] = true;
                        for &(k, r) in &g.node_rhs[m] {
                            fk[k] += r;
                        }
                    }
                    let dem: f64 = fk.iter().map(|x| (dir * x).max(0.0)).sum();
                    if dem <= 1e-9 {
                        continue;
                    }
                    let cut: Vec<usize> = (0..g.arcs.len())
                        .filter(|&a| {
                            let (_, _, t, h) = g.arcs[a];
                            if dir > 0.0 {
                                s_set[t] && !s_set[h]
                            } else {
                                !s_set[t] && s_set[h]
                            }
                        })
                        .collect();
                    if cut.is_empty() || !seen.insert(cut.clone()) {
                        continue;
                    }
                    let caps: Vec<f64> = cut.iter().map(|&a| g.arcs[a].1).collect();
                    let ys: Vec<f64> = cut.iter().map(|&a| sol[g.arcs[a].0 as usize]).collect();
                    if let Some((coef, viol, rhs)) = cutset_mir(&caps, &ys, dem) {
                        let inds: Vec<HighsInt> = cut.iter().map(|&a| g.arcs[a].0).collect();
                        out.push((inds, coef, rhs, viol));
                    }
                }
            }
        }
        out.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
        out.truncate(MAX_CUTSET_CUTS);
        out
    }

    /// Slope scaling (Crainic, Gendron, Hernu 2004) on a copy of the root LP.
    fn slope_scaling(
        &mut self,
        ms: &HighsMipSolver,
        lp: &HighsLpRelaxation,
        nrow: HighsInt,
    ) -> Option<Vec<f64>> {
        let budget = (HEUR_TIME_SHARE * ms.options_mip_.time_limit).max(2.0);
        let t0 = std::time::Instant::now();
        let model = lp.get_lp().clone();
        let ncol = model.num_col_ as usize;
        let cost = model.col_cost_.clone();
        let offset = model.offset_;
        let int_cols: Vec<usize> = (0..ncol)
            .filter(|&c0| {
                !matches!(
                    ms.variable_type_2(c0 as HighsInt),
                    HighsVarType::kContinuous
                )
            })
            .collect();
        let mut arc_of_y = vec![usize::MAX; ncol];
        for (a, arc) in self.arcs.iter().enumerate() {
            if cost[arc.y as usize] < 0.0 {
                return None;
            }
            arc_of_y[arc.y as usize] = a;
        }
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_int("threads", 1);
        h.set_option_value_str("presolve", kHighsOffString);
        h.pass_model(model);
        // Start from the root LP basis: a cold start costs a pass the size of the root solve
        // (dano3_3: 48k iterations, 19 s, twice). Rows past nrow are
        // cuts; deleting them keeps the basis when it can.
        if ms.mipdata_.total_lp_iterations >= WARM_MIN_ROOT_ITERS && lp.lpsolver.get_basis().valid {
            let basis = lp.lpsolver.get_basis().clone();
            h.set_basis(&basis, "netdesign slope scaling");
        }
        if h.get_lp().num_row_ > nrow {
            h.delete_rows_interval(nrow, h.get_lp().num_row_ - 1);
        }
        let mut slope: Vec<f64> = self.arcs.iter().map(|a| cost[a.y as usize] / a.u).collect();
        for arc in &self.arcs {
            h.change_col_cost(arc.y, 0.0);
        }
        // true cost and integrality repair of an LP point: y = (flow > 0)
        let finish = |val: &[f64], arcs: &[NdArc]| -> Option<(Vec<f64>, f64, Vec<bool>)> {
            let mut x = val.to_vec();
            let mut open = vec![false; arcs.len()];
            for (a, arc) in arcs.iter().enumerate() {
                let flow: f64 = arc.xs.iter().map(|&c0| val[c0 as usize]).sum();
                open[a] = flow > 1e-9;
                x[arc.y as usize] = if open[a] { 1.0 } else { 0.0 };
            }
            for &c0 in &int_cols {
                // integer columns outside the capacity rows: nearest value, the solution check
                // rejects the point if that breaks a row
                if arc_of_y[c0] == usize::MAX {
                    x[c0] = x[c0].round();
                }
            }
            let obj: f64 = (0..x.len()).map(|c0| cost[c0] * x[c0]).sum();
            Some((x, obj, open))
        };
        let mut best: Option<(Vec<f64>, f64, Vec<bool>)> = None;
        let mut stale = 0;
        for _ in 0..SLOPE_PASSES {
            let left = budget - t0.elapsed().as_secs_f64();
            if left <= 0.0 {
                break;
            }
            h.set_option_value_double("time_limit", left);
            for (a, arc) in self.arcs.iter().enumerate() {
                for &c0 in &arc.xs {
                    h.change_col_cost(c0, cost[c0 as usize] + slope[a]);
                }
            }
            h.optimize_lp();
            if h.get_model_status() != HighsModelStatus::kOptimal {
                break;
            }
            let val = h.get_solution().col_value.clone();
            let Some(cand) = finish(&val, &self.arcs) else {
                break;
            };
            if best
                .as_ref()
                .is_none_or(|b| cand.1 < b.1 - 1e-9 * b.1.abs())
            {
                best = Some(cand);
                stale = 0;
            } else {
                stale += 1;
                if stale >= 5 {
                    break;
                }
            }
            for (a, arc) in self.arcs.iter().enumerate() {
                let flow: f64 = arc.xs.iter().map(|&c0| val[c0 as usize]).sum();
                if flow > 1e-9 {
                    slope[a] = cost[arc.y as usize] / flow;
                }
            }
        }
        let (mut x, mut obj, open) = best?;
        // true-cost flows on the chosen design
        let left = budget - t0.elapsed().as_secs_f64();
        if left > 0.0 {
            h.set_option_value_double("time_limit", left);
            for (a, arc) in self.arcs.iter().enumerate() {
                h.change_col_cost(arc.y, cost[arc.y as usize]);
                h.change_col_bounds(
                    arc.y,
                    if open[a] { 1.0 } else { 0.0 },
                    if open[a] { 1.0 } else { 0.0 },
                );
                for &c0 in &arc.xs {
                    h.change_col_cost(c0, cost[c0 as usize]);
                }
            }
            h.optimize_lp();
            if h.get_model_status() == HighsModelStatus::kOptimal {
                let val = h.get_solution().col_value.clone();
                if let Some(cand) = finish(&val, &self.arcs) {
                    if cand.1 < obj {
                        x = cand.0;
                        obj = cand.1;
                    }
                }
            }
        }
        crate::highs_log_user!(
            &ms.options_mip_.log_options,
            HighsLogType::kInfo,
            "Network design slope scaling: objective {} in {:.1}s\n",
            obj + offset,
            t0.elapsed().as_secs_f64()
        );
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
        if ms.submip {
            return;
        }
        if !self.detected {
            self.detect(ms, lp_relaxation, trans_lp.globaldom_);
        }
        if ms.options_mip_.mip_netdesign_heur && !self.heur_done && !self.arcs.is_empty() {
            self.heur_done = true;
            let nrow = ms.num_row();
            if let Some(x) = self.slope_scaling(ms, lp_relaxation, nrow) {
                self.pending
                    .push((crate::mip::separation::next_pending_seq(), x));
            }
        }
        let sol = &lp_relaxation.get_solution().col_value;
        if ms.options_mip_.mip_netdesign_cutset && !self.cs.arcs.is_empty() {
            let cuts = self.separate_cutsets(sol);
            for (mut inds, mut vals, rhs, _) in cuts {
                let len = inds.len() as HighsInt;
                detach!(ms.mipdata_.cutpools[cutpool], |p| p.add_cut(
                    ms,
                    &mut inds[..],
                    &mut vals[..],
                    len,
                    rhs,
                    false,
                    true,
                    false,
                    false,
                ));
            }
        }
        if !ms.options_mip_.mip_netdesign_sepa || self.links.is_empty() {
            return;
        }
        let mut cand: Vec<(f64, usize)> = Vec::new();
        for (i, &(x, y, b)) in self.links.iter().enumerate() {
            let viol = sol[x as usize] - b * sol[y as usize];
            if viol > MIN_VIOLATION {
                cand.push((viol, i));
            }
        }
        cand.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
        cand.truncate(MAX_CUTS_PER_CALL);
        for (_, i) in cand {
            let (x, y, b) = self.links[i];
            let mut inds = [x, y];
            let mut vals = [1.0, -b];
            detach!(ms.mipdata_.cutpools[cutpool], |p| p.add_cut(
                ms,
                &mut inds[..],
                &mut vals[..],
                2,
                0.0,
                false,
                true,
                false,
                false,
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn negative_cost_commodity_is_skipped() {
        // two commodities of two columns each; column 2 has a negative cost
        let mut uf = vec![0, 0, 2, 2];
        let in_flow = vec![true; 4];
        let lower = vec![0.0; 4];
        let cost = vec![1.0, 2.0, -1.0, 3.0];
        let one = vec![1u8; 4];
        let bad = vec![false; 4];
        let ok = component_guard(&mut uf, &in_flow, &lower, &cost, 1.0, &one, &one, &bad);
        assert!(ok[find(&mut uf, 0)]);
        assert!(!ok[find(&mut uf, 2)]);
        // maximisation flips the sign: the same costs now fail commodity 0 and pass commodity 1 only if all <= 0
        let ok = component_guard(&mut uf, &in_flow, &lower, &cost, -1.0, &one, &one, &bad);
        assert!(!ok[find(&mut uf, 0)]);
    }

    #[test]
    fn cutset_mir_rounds_up_equal_capacities() {
        // three arcs of capacity 10 must carry demand 25: sum y >= 3; at y = 0.8 each the
        // point (sum 2.4) violates it by 0.6
        let (g, viol, rhs) = cutset_mir(&[10.0; 3], &[0.8; 3], 25.0).unwrap();
        assert_eq!(g, vec![-1.0; 3]);
        assert_eq!(rhs, -3.0);
        assert!((viol - 0.6).abs() < 1e-9);
        assert!(cutset_mir(&[10.0; 3], &[1.0; 3], 25.0).is_none());
    }
}
