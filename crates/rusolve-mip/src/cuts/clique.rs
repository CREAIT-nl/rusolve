//! Clique cuts and root fixings from a conflict graph on binary literals.
//!
//! The conflict graph is stored as a list of cliques: a row (after
//! complementing negative binaries and shifting every other column to its
//! minimum) whose two largest binary coefficients exceed the right-hand side
//! yields cliques of the form "largest prefix ∪ {k}" (Atamtürk, Nemhauser &
//! Savelsbergh, "Conflict graphs in solving integer programming problems",
//! EJOR 121:40–55 (2000), §3). Separation enumerates maximal cliques of the
//! LP-positive literals with Bron–Kerbosch (pivoting) restricted to the most
//! violated ones.
//!
//! Row extraction, the set-packing shortcut and the separation limits follow
//! HiGHS `HighsCliqueTable` (reference/HiGHS/highs/mip/HighsCliqueTable.cpp,
//! MIT; read as a reference text and reimplemented): `extractCliques`
//! :839–977 and :1268–1378, `bronKerboschRecurse` :196–332,
//! `separateCliques` :1605–1714; limits from HighsCliqueTable.h:121–144
//! (minimum weight 1.05, ≤ 10000 calls, ≤ 100 cliques, neighbourhood-query
//! budget 1e6 + 100·nnz). Not ported: clique merging/subsumption, implication
//! and probing-derived cliques. HiGHS's zero-weight extension of the
//! separated cliques sits behind `CutConfig::clique_extend`.
//!
//! The incumbent → fixings link over the same table (root-only,
//! `bnb::Search`): [`CliqueTable::extract_objective_cliques`] re-extracts
//! cliques from the objective cutoff row `cᵀx ≤ cutoff` on every improving
//! incumbent (HiGHS `extractObjCliques` :1380–1461, called from
//! `HighsMipSolverData::addSolutionImprovingSolution`,
//! HighsMipSolverData.cpp:1576–1591, after `propagateRootRedcost`) and
//! [`CliqueTable::propagate_fixings`] turns them into bound tightenings
//! (`addClique`'s `fixAllVarsInClique` :663–773 and the
//! `processInfeasibleVertices` cascade :1463–1567, simplified to a whole-table
//! fixpoint — the root needs no event stacks).

use super::relax::Relaxation;
use std::collections::HashSet;

/// Literal `2·col + 1` is `x_col`, `2·col` is `1 − x_col`.
type Lit = u32;

fn lit(col: usize, positive: bool) -> Lit {
    (2 * col + usize::from(positive)) as Lit
}

fn lit_col(l: Lit) -> usize {
    (l / 2) as usize
}

fn lit_positive(l: Lit) -> bool {
    l & 1 == 1
}

/// Hard cap on stored clique entries (memory guard for long knapsack rows).
const MAX_ENTRIES: usize = 2_000_000;
const MIN_WEIGHT: f64 = 1.05;
const MAX_CALLS: u32 = 10_000;
const MAX_CLIQUES: usize = 100;

/// One xorshift64 step (Marsaglia 2003) on the separator's seed.
fn xorshift(s: &mut u64) -> u64 {
    let mut x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    x
}

/// Outcome of an objective-clique extraction over a box.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ObjCliques {
    /// `n` new cliques stored (possibly 0: no pair conflicts, or all dups).
    Added(u32),
    /// The box's minimum objective activity already exceeds the cutoff: no
    /// improving point exists (the search is over).
    CutoffInfeasible,
}

pub(crate) struct CliqueTable {
    cliques: Vec<Vec<Lit>>,
    lit_cliques: Vec<Vec<u32>>,
    /// Canonical (sorted) member sets, to skip re-adding the same clique when
    /// incumbents tighten the cutoff slack (HiGHS dedups on insert through
    /// its inverted hash lists).
    keys: HashSet<Vec<Lit>>,
    /// Stored entries across all cliques, against [`MAX_ENTRIES`].
    entries: usize,
    mark: Vec<u32>,
    stamp: u32,
    max_queries: u64,
}

struct Bk<'a> {
    weight: &'a [f64],
    p: Vec<Lit>,
    r: Vec<Lit>,
    cliques: Vec<Vec<Lit>>,
    w_r: f64,
    min_w: f64,
    feastol: f64,
    calls: u32,
    queries: u64,
}

impl CliqueTable {
    /// Cliques implied by the model rows under the bounds of `relax`.
    pub fn extract(relax: &Relaxation, feastol: f64, epsilon: f64) -> Self {
        let mut t = CliqueTable {
            cliques: Vec::new(),
            lit_cliques: vec![Vec::new(); 2 * relax.num_col],
            keys: HashSet::new(),
            entries: 0,
            mark: vec![0; 2 * relax.num_col],
            stamp: 0,
            max_queries: 0,
        };
        t.add_rows(&relax.rows[..relax.num_model_rows], relax, feastol, epsilon);
        t.max_queries = 1_000_000 + 100 * t.entries as u64;
        t
    }

    /// Cliques implied by `rows` under the bounds of `relax`, added to this
    /// table. The same pass the model-row build uses; [`Self::extract`] runs
    /// it on the model rows and the root loop runs it again on its accepted
    /// cut rows (HiGHS `extractCliquesFromCut`, HighsCliqueTable.cpp:
    /// 1078–1266: set-packing shortcuts and prefix cliques over separated
    /// rows — cliques from cuts are valid because the cuts are).
    pub(crate) fn add_rows(
        &mut self,
        rows: &[super::relax::LpRow],
        relax: &Relaxation,
        feastol: f64,
        epsilon: f64,
    ) {
        let mut nnz = 0u64;
        for r in rows {
            nnz += r.inds.len() as u64;
            for (rhs, dir) in [(r.row_upper, 1.0), (r.row_lower, -1.0)] {
                if !rhs.is_finite() || self.entries > MAX_ENTRIES {
                    continue;
                }
                // Shift every column to the bound minimizing its term.
                let mut rhs = dir * rhs;
                let mut bins: Vec<(f64, Lit)> = Vec::new();
                let mut ok = true;
                for (&j, &a) in r.inds.iter().zip(&r.vals) {
                    let v = dir * a;
                    if v.abs() < epsilon {
                        continue;
                    }
                    let bound = if v < 0.0 {
                        relax.upper[j]
                    } else {
                        relax.lower[j]
                    };
                    if !bound.is_finite() {
                        ok = false;
                        break;
                    }
                    rhs -= v * bound;
                    if relax.is_binary(j) {
                        bins.push((v.abs(), lit(j, v > 0.0)));
                    }
                }
                if ok && bins.len() > 1 {
                    self.add_row_cliques(&mut bins, rhs, feastol);
                }
            }
        }
        self.max_queries = (self.max_queries).max(1_000_000 + 100 * nnz);
    }

    pub fn is_empty(&self) -> bool {
        self.cliques.is_empty()
    }

    /// HighsCliqueTable::extractCliques on a row `Σ w·lit ≤ rhs` (+ terms
    /// with zero minimum). Returns the number of stored entries.
    fn add_row_cliques(&mut self, bins: &mut [(f64, Lit)], rhs: f64, feastol: f64) -> usize {
        bins.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(b.1.cmp(&a.1)));
        let n = bins.len();
        if bins[0].0 + bins[1].0 <= rhs + feastol {
            return 0;
        }
        if (bins[0].0 - bins[n - 1].0).abs() <= feastol && rhs < 2.0 * bins[n - 1].0 - feastol {
            let c: Vec<Lit> = bins.iter().map(|b| b.1).collect();
            return self.add_clique(c);
        }
        let mut added = 0;
        for k in (1..n).rev() {
            let min_val = rhs - bins[k].0 + feastol;
            let end = bins[..k].partition_point(|b| b.0 > min_val);
            if end == 0 {
                continue;
            }
            let mut c: Vec<Lit> = bins[..end].iter().map(|b| b.1).collect();
            c.push(bins[k].1);
            added += self.add_clique(c);
            // Further cliques are subsets of this one.
            if end == k {
                break;
            }
        }
        added
    }

    fn add_clique(&mut self, mut c: Vec<Lit>) -> usize {
        if c.len() < 2 || self.entries >= MAX_ENTRIES {
            return 0;
        }
        c.sort_unstable();
        if !self.keys.insert(c.clone()) {
            return 0;
        }
        let id = self.cliques.len() as u32;
        for &l in &c {
            self.lit_cliques[l as usize].push(id);
        }
        let len = c.len();
        self.entries += len;
        self.cliques.push(c);
        len
    }

    pub fn num_cliques(&self) -> usize {
        self.cliques.len()
    }

    /// A 2-clique from root probing — the literals `(col_a, pos_a)`
    /// and `(col_b, pos_b)` cannot both hold (probing fixed `col_b` away
    /// from `pos_b` under `col_a = pos_a`; HiGHS adds these as 2-cliques
    /// from `computeImplications`, HighsImplications.cpp:101–112, MIT).
    /// Returns whether a new clique was stored.
    pub(crate) fn add_probing_clique(&mut self, a: (usize, bool), b: (usize, bool)) -> bool {
        if a.0 == b.0 {
            return false;
        }
        self.add_clique(vec![lit(a.0, a.1), lit(b.0, b.1)]) > 0
    }

    pub fn num_entries(&self) -> usize {
        self.entries
    }

    /// The objective cutoff row as a clique source (HiGHS
    /// `HighsCliqueTable::extractObjCliques`, HighsCliqueTable.cpp:1380–1461,
    /// MIT; called on every improving incumbent,
    /// HighsMipSolverData.cpp:1576–1591).
    ///
    /// Every improving point satisfies `cᵀx ≤ cutoff`, and raising a binary
    /// from its minimum-contribution value to the other one adds exactly
    /// `|c_j|` to the row's minimum activity over the box, so two objective
    /// binaries with `|c_i| + |c_j| > cutoff − minact + feastol` can never
    /// both sit at their objective-raising value: a clique over the literals
    /// `(j, c_j > 0)`. Only free binaries enter (HiGHS partitions on
    /// `vals != 0 && !isFixed`); the minimum activity runs over *every*
    /// objective column at the current box, which keeps the clique valid for
    /// all points in the (monotonically tightening) global domain.
    pub fn extract_objective_cliques(
        &mut self,
        cost: &[f64],
        lower: &[f64],
        upper: &[f64],
        integer: &[bool],
        cutoff: f64,
        feastol: f64,
    ) -> ObjCliques {
        // Minimum activity of the objective row: -inf (an unbounded-below
        // contribution) means the cutoff row cannot bound anything.
        let mut minact = 0.0;
        for (j, &c) in cost.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            let (lo, up) = (c * lower[j], c * upper[j]);
            if lo == f64::NEG_INFINITY || up == f64::NEG_INFINITY || lo.min(up) == f64::NEG_INFINITY
            {
                return ObjCliques::Added(0);
            }
            minact += lo.min(up);
        }
        let slack = cutoff - minact;
        if slack < -feastol {
            return ObjCliques::CutoffInfeasible;
        }
        // Free binary objective columns, largest |c| first (HiGHS sorts
        // |vals| descending; ties by position for determinism).
        let mut bins: Vec<(f64, Lit)> = cost
            .iter()
            .enumerate()
            .filter(|&(j, &c)| {
                c != 0.0 && integer[j] && lower[j] >= 0.0 && upper[j] <= 1.0 && lower[j] < upper[j]
            })
            .map(|(j, &c)| (c.abs(), lit(j, c > 0.0)))
            .collect();
        if bins.len() <= 1 {
            return ObjCliques::Added(0);
        }
        bins.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        if bins[0].0 + bins[1].0 <= slack + feastol {
            return ObjCliques::Added(0);
        }
        let mut added = 0;
        for k in (1..bins.len()).rev() {
            // Members whose |c| conflicts with k's: |c_j| + |c_k| > slack.
            let min_val = slack - bins[k].0 + feastol;
            let end = bins[..k].partition_point(|b| b.0 > min_val);
            if end == 0 {
                continue;
            }
            let mut c: Vec<Lit> = bins[..end].iter().map(|b| b.1).collect();
            c.push(bins[k].1);
            added += u32::from(self.add_clique(c) > 0);
            // Further cliques are subsets of this one.
            if end == k {
                break;
            }
        }
        ObjCliques::Added(added)
    }

    /// Root fixings fixpoint over the whole table (HiGHS `addClique`'s
    /// `fixAllVarsInClique`, HighsCliqueTable.cpp:663–773, plus the
    /// `processInfeasibleVertices` cascade :1463–1567): a clique member
    /// binary at its literal value (`x = 1` for `x`, `x = 0` for `1 − x`)
    /// fixes every other member to the opposite value; two members at their
    /// literal value prove the box empty. Event-driven like HiGHS's vertex
    /// stack: a fixing re-queues only the cliques the column belongs to.
    /// Appends tightened columns to `changed`; `false` = no improving point.
    pub fn propagate_fixings(
        &self,
        lower: &mut [f64],
        upper: &mut [f64],
        feastol: f64,
        changed: &mut Vec<usize>,
    ) -> bool {
        let mut queue: Vec<u32> = (0..self.cliques.len() as u32).collect();
        let mut queued = vec![true; self.cliques.len()];
        while let Some(id) = queue.pop() {
            queued[id as usize] = false;
            let members = &self.cliques[id as usize];
            // The (at most one) member sitting at its literal value.
            let mut at_lit: Option<Lit> = None;
            for &l in members {
                let (col, pos) = (lit_col(l), lit_positive(l));
                let at = if pos {
                    lower[col] >= 1.0 - feastol
                } else {
                    upper[col] <= feastol
                };
                if at {
                    if at_lit.is_some() {
                        return false;
                    }
                    at_lit = Some(l);
                }
            }
            let Some(fixed) = at_lit else { continue };
            for &l in members {
                if l == fixed {
                    continue;
                }
                // Fix to the literal's opposite value: a positive literal
                // `x` forces x = 0, a negative `1 − x` forces x = 1. Binary
                // bounds stay integral, so the tolerances compare to 0/1.
                let (col, pos) = (lit_col(l), lit_positive(l));
                let already = if pos {
                    upper[col] <= feastol
                } else {
                    lower[col] >= 1.0 - feastol
                };
                if already {
                    continue;
                }
                let crosses = if pos {
                    lower[col] > feastol
                } else {
                    upper[col] < 1.0 - feastol
                };
                if crosses {
                    // The opposite bound contradicts the box: no improving
                    // point exists.
                    return false;
                }
                if pos {
                    upper[col] = 0.0;
                } else {
                    lower[col] = 1.0;
                }
                changed.push(col);
                // The column moved: re-check every clique it can fix in.
                for l2 in [lit(col, true), lit(col, false)] {
                    for &c2 in &self.lit_cliques[l2 as usize] {
                        if !queued[c2 as usize] {
                            queued[c2 as usize] = true;
                            queue.push(c2);
                        }
                    }
                }
            }
        }
        true
    }

    /// Marks the neighbours of `v` (literals sharing a clique, `v` excluded).
    fn mark_neighbours(&mut self, v: Lit, queries: &mut u64) {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.mark.fill(0);
            self.stamp = 1;
        }
        for &c in &self.lit_cliques[v as usize] {
            let members = &self.cliques[c as usize];
            *queries += members.len() as u64;
            for &l in members {
                self.mark[l as usize] = self.stamp;
            }
        }
        self.mark[v as usize] = 0;
    }

    /// Moves the neighbours of `v` in `list` to its front; returns their
    /// count.
    fn partition_neighbours(&mut self, v: Lit, list: &mut [Lit], queries: &mut u64) -> usize {
        self.mark_neighbours(v, queries);
        let mut k = 0;
        for i in 0..list.len() {
            if self.mark[list[i] as usize] == self.stamp {
                list.swap(k, i);
                k += 1;
            }
        }
        k
    }

    /// Most violated cliques at `x` as cuts `Σ vals·x[inds] ≤ rhs`.
    ///
    /// With `extend` (`CutConfig::clique_extend`) each separated clique
    /// is lifted greedily with the zero-weight literals adjacent to all of
    /// its members, visited in a random order drawn from `rng`: HiGHS
    /// `separateCliques` under its compiled-in `ADD_ZERO_WEIGHT_VARS`
    /// (HighsCliqueTable.cpp:25, :1621–1682). The extension leaves the
    /// violation at `x` unchanged but makes the cut a (near-)maximal clique,
    /// so later LP points cannot move the weight onto the omitted literals.
    pub fn separate(
        &mut self,
        x: &[f64],
        feastol: f64,
        extend: bool,
        rng: &mut u64,
    ) -> Vec<(Vec<usize>, Vec<f64>, f64)> {
        let num_lit = self.lit_cliques.len();
        let weight: Vec<f64> = (0..num_lit)
            .map(|l| {
                let v = x[lit_col(l as Lit)];
                if lit_positive(l as Lit) { v } else { 1.0 - v }
            })
            .collect();
        let p: Vec<Lit> = (0..num_lit as Lit)
            .filter(|&l| !self.lit_cliques[l as usize].is_empty() && weight[l as usize] > feastol)
            .collect();
        let mut bk = Bk {
            weight: &weight,
            p,
            r: Vec::new(),
            cliques: Vec::new(),
            w_r: 0.0,
            min_w: MIN_WEIGHT,
            feastol,
            calls: 0,
            queries: 0,
        };
        let plen = bk.p.len();
        self.recurse(&mut bk, plen, &[]);
        let mut queries = bk.queries;
        let mut cliques = bk.cliques;
        if extend {
            // Zero-weight literals that lie in some clique
            // (HighsCliqueTable.cpp:1621–1633).
            let mut z: Vec<Lit> = (0..num_lit as Lit)
                .filter(|&l| {
                    !self.lit_cliques[l as usize].is_empty() && weight[l as usize] <= feastol
                })
                .collect();
            for c in &mut cliques {
                // Common neighbours of the clique (:1660–1666).
                let mut end = z.len();
                for &v in c.iter() {
                    end = self.partition_neighbours(v, &mut z[..end], &mut queries);
                    if end == 0 {
                        break;
                    }
                }
                if end == 0 {
                    continue;
                }
                // Random greedy maximal extension (:1668–1680).
                for i in (1..end).rev() {
                    let j = (xorshift(rng) % (i as u64 + 1)) as usize;
                    z.swap(i, j);
                }
                let mut i = 0;
                while i < end {
                    let k = i + 1;
                    let (head, tail) = z.split_at_mut(k);
                    end =
                        k + self.partition_neighbours(head[i], &mut tail[..end - k], &mut queries);
                    i += 1;
                }
                c.extend_from_slice(&z[..end]);
            }
        }
        let mut out = Vec::with_capacity(cliques.len());
        for c in cliques {
            let mut rhs = 1.0;
            let mut inds = Vec::with_capacity(c.len());
            let mut vals = Vec::with_capacity(c.len());
            for l in c {
                inds.push(lit_col(l));
                if lit_positive(l) {
                    vals.push(1.0);
                } else {
                    vals.push(-1.0);
                    rhs -= 1.0;
                }
            }
            out.push((inds, vals, rhs));
        }
        out
    }

    fn stop(&self, bk: &Bk) -> bool {
        bk.calls >= MAX_CALLS || bk.cliques.len() >= MAX_CLIQUES || bk.queries > self.max_queries
    }

    fn recurse(&mut self, bk: &mut Bk, mut plen: usize, x: &[Lit]) {
        let weight = bk.weight;
        let wt = |l: Lit| weight[l as usize];
        let mut w = bk.w_r + bk.p[..plen].iter().map(|&l| wt(l)).sum::<f64>();
        if w < bk.min_w - bk.feastol {
            return;
        }
        if plen == 0 && x.is_empty() {
            if bk.min_w < w - bk.feastol {
                bk.cliques.clear();
                bk.min_w = w;
            }
            bk.cliques.push(bk.r.clone());
            return;
        }
        bk.calls += 1;
        if self.stop(bk) {
            return;
        }
        // Pivot: heaviest vertex of X ∪ P.
        let mut pivot = None;
        let mut pivw = -1.0;
        for &l in x.iter().chain(&bk.p[..plen]) {
            if wt(l) > pivw {
                pivw = wt(l);
                pivot = Some(l);
                if pivw >= 1.0 - bk.feastol {
                    break;
                }
            }
        }
        let Some(pivot) = pivot else {
            return;
        };
        self.mark_neighbours(pivot, &mut bk.queries);
        let mut branch: Vec<Lit> = bk.p[..plen]
            .iter()
            .copied()
            .filter(|&l| self.mark[l as usize] != self.stamp)
            .collect();
        branch.sort_unstable_by(|&a, &b| wt(b).total_cmp(&wt(a)).then(b.cmp(&a)));
        let mut local_x: Vec<Lit> = x.to_vec();
        for v in branch {
            let new_plen = self.partition_neighbours(v, &mut bk.p[..plen], &mut bk.queries);
            let new_xlen = self.partition_neighbours(v, &mut local_x, &mut bk.queries);
            let wv = wt(v);
            bk.r.push(v);
            bk.w_r += wv;
            let sub_x: Vec<Lit> = local_x[..new_xlen].to_vec();
            self.recurse(bk, new_plen, &sub_x);
            if self.stop(bk) {
                return;
            }
            bk.r.pop();
            bk.w_r -= wv;
            w -= wv;
            if w < bk.min_w {
                return;
            }
            let Some(vpos) = (new_plen..plen).find(|&i| bk.p[i] == v) else {
                return;
            };
            plen -= 1;
            bk.p.swap(vpos, plen);
            local_x.push(v);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FEASTOL: f64 = 1e-6;

    fn table(num_col: usize) -> CliqueTable {
        CliqueTable {
            cliques: Vec::new(),
            lit_cliques: vec![Vec::new(); 2 * num_col],
            keys: HashSet::new(),
            entries: 0,
            mark: vec![0; 2 * num_col],
            stamp: 0,
            max_queries: 0,
        }
    }

    /// Edges {0,1},{0,2},{1,2} (the violated triangle at x = ½) plus the
    /// zero-weight column 3 adjacent to all of them and column 4 adjacent to
    /// 0 and 1 only: extension appends x3 (and never x4), the violation is
    /// unchanged, and without `extend` the cut is the bare triangle.
    #[test]
    fn clique_extension_adds_common_zero_weight_neighbours() {
        let mut t = table(5);
        t.max_queries = u64::MAX;
        for (a, b) in [
            (0, 1),
            (0, 2),
            (1, 2),
            (0, 3),
            (1, 3),
            (2, 3),
            (0, 4),
            (1, 4),
        ] {
            t.add_probing_clique((a, true), (b, true));
        }
        let x = [0.5, 0.5, 0.5, 0.0, 0.0];
        let mut rng = 0x1234_5678_9abc_def1;
        for (extend, want) in [(false, vec![0, 1, 2]), (true, vec![0, 1, 2, 3])] {
            let cuts = t.separate(&x, FEASTOL, extend, &mut rng);
            assert_eq!(cuts.len(), 1);
            let (mut inds, vals, rhs) = cuts[0].clone();
            inds.sort_unstable();
            assert_eq!(inds, want, "extend {extend}");
            assert!(vals.iter().all(|&v| v == 1.0) && rhs == 1.0);
        }
    }

    /// Three positive-cost binaries: every pair busts a tight cutoff, so the
    /// extraction yields one 3-clique; a looser cutoff drops the weakest pair
    /// and yields the two cliques that remain pairwise valid.
    #[test]
    fn objective_cliques_follow_the_cutoff_slack() {
        let cost = [10.0, 6.0, 3.0];
        let lower = [0.0; 3];
        let upper = [1.0; 3];
        let integer = [true; 3];
        let mut t = table(3);
        assert_eq!(
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, 8.0, FEASTOL),
            ObjCliques::Added(1)
        );
        assert_eq!(t.num_cliques(), 1);
        // 10 + 3 > 10 and 10 + 6 > 10, but 6 + 3 does not.
        let mut t = table(3);
        assert_eq!(
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, 10.0, FEASTOL),
            ObjCliques::Added(2)
        );
        // A re-run at the same cutoff adds nothing (dedup).
        assert_eq!(
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, 10.0, FEASTOL),
            ObjCliques::Added(0)
        );
        // No pair conflicts at all.
        let mut t = table(3);
        assert_eq!(
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, 20.0, FEASTOL),
            ObjCliques::Added(0)
        );
    }

    /// Negative costs complement the literals: both binaries at 0 bust the
    /// cutoff, so the clique is over `1 − x` literals.
    #[test]
    fn objective_cliques_complement_negative_costs() {
        let cost = [-10.0, -6.0];
        let lower = [0.0; 2];
        let upper = [1.0; 2];
        let integer = [true; 2];
        let mut t = table(2);
        assert_eq!(
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, -5.0, FEASTOL),
            ObjCliques::Added(1)
        );
        // Fixing x0 = 0 (the literal's value) must fix x1 = 1.
        let mut lo = [0.0, 0.0];
        let mut up = [0.0, 1.0];
        let mut changed = Vec::new();
        assert!(t.propagate_fixings(&mut lo, &mut up, FEASTOL, &mut changed));
        assert_eq!(lo, [0.0, 1.0]);
        assert_eq!(changed, [1]);
    }

    /// The box's minimum objective already above the cutoff: the search is
    /// over (HiGHS reports this through the domain's objective bound).
    #[test]
    fn objective_cliques_report_cutoff_infeasible() {
        let cost = [10.0, 0.0];
        let lower = [0.0; 2];
        let upper = [1.0; 2];
        let integer = [true; 2];
        let mut t = table(2);
        assert_eq!(
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, -1.0, FEASTOL),
            ObjCliques::CutoffInfeasible
        );
    }

    /// A member at its literal value fixes the other members to the opposite
    /// value, cascading through overlapping cliques, and two members at their
    /// literal value empty the box.
    #[test]
    fn propagate_fixings_cascade_and_conflict() {
        // x0 + x1 <= 1 and x2 <= x1, the latter as (1-x1) + x2 <= 1.
        let mut t = table(3);
        t.add_clique(vec![lit(0, true), lit(1, true)]);
        t.add_clique(vec![lit(1, false), lit(2, true)]);
        // x0 = 1: the first clique fixes x1 = 0, which makes ~x1 true in the
        // second and fixes x2 = 0.
        let mut lo = [1.0, 0.0, 0.0];
        let mut up = [1.0, 1.0, 1.0];
        let mut changed = Vec::new();
        assert!(t.propagate_fixings(&mut lo, &mut up, FEASTOL, &mut changed));
        assert_eq!(up, [1.0, 0.0, 0.0]);
        assert_eq!(changed, [1, 2]);
        // Two members at their literal value: infeasible.
        let mut lo = [1.0, 1.0, 0.0];
        let mut up = [1.0, 1.0, 1.0];
        assert!(!t.propagate_fixings(&mut lo, &mut up, FEASTOL, &mut Vec::new()));
        // Every literal falsified by the box: inert.
        let mut lo = [0.0, 1.0, 0.0];
        let mut up = [0.0, 1.0, 1.0];
        let mut changed = Vec::new();
        assert!(t.propagate_fixings(&mut lo, &mut up, FEASTOL, &mut changed));
        assert!(changed.is_empty());
    }

    /// Objective cliques must never cut off an improving point, whatever the
    /// box: enumerate every binary point of a small knapsack and verify no
    /// point with objective below the cutoff violates a clique, and that the
    /// extraction actually fires on it.
    #[test]
    fn objective_cliques_exclude_no_improving_point() {
        let cost = [30.0, 25.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0];
        let w = [5.0, 5.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0];
        let cap = 9.0;
        let n = cost.len();
        let feasible =
            |x: &[f64]| w.iter().zip(x).map(|(a, b)| a * b).sum::<f64>() <= cap + FEASTOL;
        let obj = |x: &[f64]| cost.iter().zip(x).map(|(c, v)| c * v).sum::<f64>();
        // The optimum takes both negatives plus one unit item; a cutoff three
        // above it leaves improving-but-suboptimal points in the check.
        let opt = (0u32..(1 << n))
            .map(|bits| {
                let x: Vec<f64> = (0..n).map(|j| f64::from((bits >> j) & 1)).collect();
                if feasible(&x) { obj(&x) } else { f64::INFINITY }
            })
            .fold(f64::INFINITY, f64::min);
        let cutoff = opt + 3.0;
        let mut t = table(n);
        let lower = vec![0.0; n];
        let upper = vec![1.0; n];
        let integer = vec![true; n];
        let ObjCliques::Added(k) =
            t.extract_objective_cliques(&cost, &lower, &upper, &integer, cutoff, FEASTOL)
        else {
            panic!("expected cliques");
        };
        assert!(k > 0);
        for bits in 0u32..(1 << n) {
            let x: Vec<f64> = (0..n).map(|j| f64::from((bits >> j) & 1)).collect();
            if !feasible(&x) || obj(&x) >= cutoff {
                continue;
            }
            for cl in &t.cliques {
                let at_value = cl
                    .iter()
                    .filter(|&&l| {
                        let (col, pos) = (lit_col(l), lit_positive(l));
                        x[col] == if pos { 1.0 } else { 0.0 }
                    })
                    .count();
                assert!(
                    at_value <= 1,
                    "clique {cl:?} cuts off improving point {bits:b} (obj {})",
                    obj(&x)
                );
            }
        }
    }
}
