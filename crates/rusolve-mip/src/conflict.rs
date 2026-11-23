//! Conflict analysis from LP dual proofs, and the conflict pool.
//! A port of the dual-proof half of HiGHS
//! conflict analysis; HiGHS is MIT-licensed and was read as a reference text,
//! then reimplemented (not translated). All line numbers refer to
//! `reference/HiGHS/highs/mip/`.
//!
//! * A node LP that is infeasible comes with a Farkas ray, i.e. a globally
//!   valid row `Σ aⱼxⱼ ≤ ρ` whose minimum activity over the node box exceeds
//!   `ρ` ([`clean_proof`]; HighsLpRelaxation.cpp:928–1045 `storeDualInfProof`).
//! * [`explain_leq`] picks the bound changes (node bound tighter than the
//!   root bound) that push the proof's activity past `ρ`, strongest first,
//!   then relaxes them as far as the slack allows
//!   (HighsDomain.cpp:3423–3467 `explainInfeasibilityLeq`, :3094–3203
//!   `resolveLinearLeq`). The result is a *conflict*: a set of bound literals
//!   that cannot all hold in a solution.
//! * [`ConflictPool`] stores conflicts with ageing
//!   (HighsConflictPool.cpp:12–82 `addConflictCut`, :107–128 `removeConflict`,
//!   :130–158 `performAging`) and propagates them
//!   (HighsDomain.cpp:305–372 `ConflictPoolPropagation::propagateConflict`):
//!   with every literal but one holding, the last one is negated; with all
//!   holding, the node is infeasible.
//!
//! Adapted, not ported: HiGHS keeps a stack of domain
//! changes with a reason each and resolves the proof's literals depth by
//! depth towards the branching decisions (`resolveDepth`/`computeCuts`,
//! HighsDomain.cpp:3702–3825) and adds reconvergence (UIP) cuts. Our nodes
//! carry their box and no reasons, so a conflict stays in the literals the
//! proof itself uses — the case HiGHS falls back on when resolution finds
//! nothing (HighsDomain.cpp:3869–3873). The duplicate filter and one `Vec`
//! per conflict instead of a flat array with a free-space set are ours.
//!
//! Propagation drives itself from per-column watch lists like HiGHS's
//! `ConflictPoolPropagation` (HighsDomain.cpp:88–123 link/unlink, :179
//! `conflictAdded`, :269–303 `updateActivityLb/UbChange`, :305
//! `propagateConflict`), with one adaptation: HiGHS owns a domain it
//! backtracks and gets bound-change events, while our pool is shared across
//! unrelated node boxes. So [`Revisit::Columns`] (a node extending its
//! parent's box: bounds only tighten along a lineage, and a false watched
//! literal only becomes true through a move of its own column) serves
//! plunge children, and [`Revisit::All`] — today's full scan — runs
//! whenever the queue jumps the box to another lineage.

use std::collections::HashSet;

/// `x_col ≥ val` (`upper == false`) or `x_col ≤ val` (`upper == true`): one
/// element of a conflict (HiGHS `HighsDomainChange`, HighsDomainChange.h).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundLit {
    pub col: u32,
    pub upper: bool,
    pub val: f64,
}

impl BoundLit {
    /// The literal holds in the box (HiGHS `HighsDomain::isActive`,
    /// HighsDomain.h:474–478).
    #[inline]
    fn holds(&self, lower: &[f64], upper: &[f64]) -> bool {
        let j = self.col as usize;
        if self.upper {
            self.val >= upper[j]
        } else {
            self.val <= lower[j]
        }
    }

    /// The negation as a bound: `x ≥ v` becomes `x ≤ v − feastol`, rounded
    /// down for integers (HighsDomain.cpp:2736–2750 `HighsDomain::flip`).
    fn negated(&self, integer: bool, feastol: f64) -> BoundLit {
        let mut val = if self.upper {
            self.val + feastol
        } else {
            self.val - feastol
        };
        if integer {
            val = if self.upper { val.ceil() } else { val.floor() };
        }
        BoundLit {
            col: self.col,
            upper: !self.upper,
            val,
        }
    }
}

/// Counters of one search.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ConflictStats {
    /// Infeasible node LPs that gave a usable Farkas proof.
    pub proofs: u64,
    /// Of [`Self::proofs`], those from cutoff nodes (dual-bound proofs).
    pub cutoff_proofs: u64,
    /// Infeasible node LPs without a ray, or whose proof was unusable
    /// (infinite root bound in it, no bound change explains it).
    pub no_proof: u64,
    /// Conflicts stored.
    pub added: u64,
    /// Explained conflicts equal to a stored one, or too long to keep.
    pub rejected: u64,
    /// Bound tightenings made by propagating stored conflicts.
    pub tightenings: u64,
    /// Nodes found infeasible by a stored conflict.
    pub pruned: u64,
    /// Conflicts aged out of the pool.
    pub aged_out: u64,
}

/// Inputs of [`explain_leq`] that stay fixed across a search.
pub struct Explain<'a> {
    pub global_lower: &'a [f64],
    pub global_upper: &'a [f64],
    pub integer: &'a [bool],
    pub feastol: f64,
    /// Open/observed nodes branched up and down per column: the tie-break
    /// weight of a candidate (HiGHS `nodequeue.numNodesUp/Down`).
    pub up_count: &'a dyn Fn(usize) -> u32,
    pub down_count: &'a dyn Fn(usize) -> u32,
}

/// `Σ vals·x ≤ rhs`, valid for every solution in the root box.
#[derive(Debug, Clone, PartialEq)]
pub struct Proof {
    pub inds: Vec<usize>,
    pub vals: Vec<f64>,
    pub rhs: f64,
}

/// Entries at most this large (after scaling the row to unit maximum) are
/// moved to the right-hand side through the root bound.
const TINY_COEF: f64 = 1e-6;
/// HiGHS `epsilon`: a relaxed bound this close to the root bound is dropped.
const EPSILON: f64 = 1e-9;

/// Turns the dense aggregated row `Σ coefs·x ≤ rhs` of a Farkas ray into a
/// [`Proof`]: scale to unit maximum, move negligible entries into the
/// right-hand side through the root bound (HighsLpRelaxation.cpp:1004–1033;
/// there the bound is also used for continuous columns the node did not
/// tighten, which [`explain_leq`] covers implicitly), then strengthen
/// integer coefficients (HighsDomain.cpp:2668–2706 `tightenCoefficients`).
/// `None` if a dropped entry meets an infinite root bound.
pub fn clean_proof(
    coefs: &[f64],
    rhs: f64,
    global_lower: &[f64],
    global_upper: &[f64],
    integer: &[bool],
    feastol: f64,
) -> Option<Proof> {
    let max = coefs.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    if !(max.is_finite() && max > 0.0) || !rhs.is_finite() {
        return None;
    }
    let scale = 1.0 / max;
    let mut rhs = rhs * scale;
    let (mut inds, mut vals) = (Vec::new(), Vec::new());
    for (j, &c) in coefs.iter().enumerate() {
        let v = c * scale;
        if v == 0.0 {
            continue;
        }
        if v.abs() <= TINY_COEF {
            // v·x ≥ v·(root bound that minimises it), so the entry can go to
            // the right-hand side.
            let b = if v > 0.0 {
                global_lower[j]
            } else {
                global_upper[j]
            };
            if !b.is_finite() {
                return None;
            }
            rhs -= v * b;
            continue;
        }
        inds.push(j);
        vals.push(v);
    }
    if inds.is_empty() {
        return None;
    }
    tighten_coefficients(
        &inds,
        &mut vals,
        &mut rhs,
        global_lower,
        global_upper,
        integer,
        feastol,
    );
    Some(Proof { inds, vals, rhs })
}

/// A coefficient of an integer column larger than the slack
/// `maxact − rhs` of the row can be cut down to it
/// (HighsDomain.cpp:2668–2706).
fn tighten_coefficients(
    inds: &[usize],
    vals: &mut [f64],
    rhs: &mut f64,
    glo: &[f64],
    ghi: &[f64],
    integer: &[bool],
    feastol: f64,
) {
    let mut maxact = 0.0;
    for (&j, &v) in inds.iter().zip(vals.iter()) {
        let b = if v > 0.0 { ghi[j] } else { glo[j] };
        if !b.is_finite() {
            return;
        }
        maxact += v * b;
    }
    let slack = maxact - *rhs;
    if slack <= feastol {
        return;
    }
    let mut upper = *rhs;
    for (&j, v) in inds.iter().zip(vals.iter_mut()) {
        if !integer[j] {
            continue;
        }
        if *v > slack {
            upper -= (*v - slack) * ghi[j];
            *v = slack;
        } else if *v < -slack {
            upper += (-*v - slack) * glo[j];
            *v = -slack;
        }
    }
    *rhs = upper;
}

struct Candidate {
    /// Position in the proof.
    pos: usize,
    delta: f64,
    base: f64,
    prio: f64,
}

/// The conflict explaining `proof` in the box `lower`/`upper` (a node's):
/// bound literals that, together with the root box, make the proof's
/// minimum activity exceed its right-hand side. `None` when the root box has
/// an infinite contribution to the activity or the node's bound changes do
/// not cover the gap; `Some(vec![])` when the root box alone does.
pub fn explain_leq(
    proof: &Proof,
    lower: &[f64],
    upper: &[f64],
    ex: &Explain<'_>,
) -> Option<Vec<BoundLit>> {
    let (glo, ghi) = (ex.global_lower, ex.global_upper);
    let mut act = 0.0;
    for (&j, &v) in proof.inds.iter().zip(&proof.vals) {
        let b = if v > 0.0 { glo[j] } else { ghi[j] };
        if !b.is_finite() {
            return None;
        }
        act += v * b;
    }
    // Candidates: every proof column whose node bound is tighter than its
    // root bound, with the activity it adds and a priority of
    // |a·Δbound|·(1 + nodes branching that way)
    // (HighsDomain.cpp:3423–3452, :2765–2768 `computePrio`).
    let mut cands: Vec<Candidate> = Vec::new();
    for (pos, (&j, &v)) in proof.inds.iter().zip(&proof.vals).enumerate() {
        let (bound, base, count) = if v > 0.0 {
            (lower[j], glo[j], (ex.up_count)(j))
        } else {
            (upper[j], ghi[j], (ex.down_count)(j))
        };
        let moved = if v > 0.0 { bound > base } else { bound < base };
        if !moved {
            continue;
        }
        let delta = v * (bound - base);
        cands.push(Candidate {
            pos,
            delta,
            base,
            prio: (delta * (1.0 + f64::from(count))).abs(),
        });
    }
    cands.sort_by(|a, b| {
        b.prio
            .total_cmp(&a.prio)
            .then(proof.inds[a.pos].cmp(&proof.inds[b.pos]))
    });
    let target = proof.rhs + proof.rhs.abs().max(10.0) * ex.feastol;
    let mut m = act;
    let mut used = 0;
    if m < target {
        for c in &cands {
            m += c.delta;
            used += 1;
            if m >= target {
                break;
            }
        }
        if m < target {
            return None;
        }
    }
    // The literals as the node's bounds, then relaxed from the weakest
    // candidate up while the activity stays above the target
    // (HighsDomain.cpp:3114–3199).
    let mut lits: Vec<Option<BoundLit>> = cands[..used]
        .iter()
        .map(|c| {
            let j = proof.inds[c.pos];
            let up = proof.vals[c.pos] < 0.0;
            Some(BoundLit {
                col: j as u32,
                upper: up,
                val: if up { upper[j] } else { lower[j] },
            })
        })
        .collect();
    if m - target > ex.feastol {
        for k in (0..used).rev() {
            let c = &cands[k];
            let j = proof.inds[c.pos];
            let v = proof.vals[c.pos];
            let lit = lits[k].as_mut().expect("kept literal");
            let mut relaxed = (target - (m - c.delta)) / v + c.base;
            if v > 0.0 {
                if ex.integer[j] {
                    relaxed = relaxed.ceil();
                }
                if relaxed - lit.val >= -ex.feastol {
                    continue;
                }
                let old = lit.val;
                lit.val = relaxed;
                if relaxed - c.base <= EPSILON {
                    lits[k] = None;
                    m -= c.delta;
                } else {
                    m += v * (relaxed - old);
                }
            } else {
                if ex.integer[j] {
                    relaxed = relaxed.floor();
                }
                if relaxed - lit.val <= ex.feastol {
                    continue;
                }
                let old = lit.val;
                lit.val = relaxed;
                if relaxed - c.base >= -EPSILON {
                    lits[k] = None;
                    m -= c.delta;
                } else {
                    m += v * (relaxed - old);
                }
            }
            if m - target <= ex.feastol {
                break;
            }
        }
    }
    Some(lits.into_iter().flatten().collect())
}

struct Stored {
    lits: Vec<BoundLit>,
    /// Indices of the two watched literals (equal for a unit conflict):
    /// while both are false the conflict cannot propagate. Both are kept
    /// linked into the per-column watch lists of the side their bound
    /// change tightens ([`ConflictPool::link`]).
    watch: [u32; 2],
    age: i32,
    key: u64,
}

/// Empty linked-list slot (HiGHS `-1`, HighsDomain.cpp:70–71).
const NO_WATCH: u32 = u32::MAX;

/// One node of a per-column watch list; watch `2·slot + k` of conflict
/// `slot`, linked by the watched literal's column and side
/// (HiGHS `WatchedLiteral`, HighsDomain.h, linked/unlinked by
/// `linkWatchedLiteral`/`unlinkWatchedLiteral`, HighsDomain.cpp:88–123).
#[derive(Clone, Copy)]
struct WatchLink {
    prev: u32,
    next: u32,
}

/// Which conflicts [`ConflictPool::propagate`] considers.
#[derive(Clone, Copy)]
pub enum Revisit<'a> {
    /// Every live conflict (HiGHS propagates within one domain it
    /// backtracks; our pool is shared across unrelated boxes, so a node
    /// taken from the queue — where the box jumps to a different lineage —
    /// re-checks everything).
    All,
    /// Only conflicts watching one of these columns, plus the marked ones.
    /// Sound for a node whose box extends its parent's (bounds only
    /// tighten along a lineage): a skipped conflict has two false watched
    /// literals, and a false literal can only become true through a bound
    /// move of its own column.
    Columns(&'a [usize]),
}

/// Conflict store with HiGHS-style ageing (HighsConflictPool.h/.cpp) and
/// per-column watch lists so propagation revisits only conflicts whose
/// watched bounds moved (`ConflictPoolPropagation::updateActivityLbChange`
/// / `updateActivityUbChange`, HighsDomain.cpp:269–303, MIT).
pub struct ConflictPool {
    slots: Vec<Option<Stored>>,
    free: Vec<usize>,
    keys: HashSet<u64>,
    live: usize,
    /// Number of live conflicts of each age (`ageDistribution_`).
    age_distribution: Vec<usize>,
    age_limit: i32,
    soft_limit: usize,
    max_len: usize,
    /// Watch-list heads per column and side (`colLowerWatched_` /
    /// `colUpperWatched_`, HighsDomain.cpp:129–133).
    col_lower_watch: Vec<u32>,
    col_upper_watch: Vec<u32>,
    /// List links, two per slot (watch `2·slot + k`).
    watch_links: Vec<WatchLink>,
    /// Conflicts added since the last propagate: like HiGHS's
    /// `conflictAdded` → `markPropagateConflict` (HighsDomain.cpp:179–256)
    /// a new conflict is propagated at the next node whatever moved.
    marked: Vec<u32>,
    /// Work list kept across calls to avoid reallocating per node.
    work: Vec<u32>,
}

impl ConflictPool {
    /// `age_limit = 5·mip_pool_age_limit = 150`, `soft_limit = 10000`
    /// (HighsMipSolverData.cpp:77; options default 30 and 10000,
    /// HighsOptions.h:1184–1193); conflicts longer than
    /// `(1000 + 3·num_integer)/10` are not kept
    /// (HighsDomain.cpp:3841–3843).
    pub fn new(num_col: usize, num_integer: usize) -> Self {
        let (age_limit, soft_limit) = (150, 10_000);
        Self {
            slots: Vec::new(),
            free: Vec::new(),
            keys: HashSet::new(),
            live: 0,
            age_distribution: vec![0; age_limit as usize + 1],
            age_limit,
            soft_limit,
            max_len: (1000 + 3 * num_integer) / 10,
            col_lower_watch: vec![NO_WATCH; num_col],
            col_upper_watch: vec![NO_WATCH; num_col],
            watch_links: Vec::new(),
            marked: Vec::new(),
            work: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.live
    }

    pub fn is_empty(&self) -> bool {
        self.live == 0
    }

    pub fn over_soft_limit(&self) -> bool {
        self.live > self.soft_limit
    }

    /// Stores a conflict; continuous bounds are moved by `feastol` away from
    /// the node's value, so the literal only holds when the bound is clearly
    /// met (HighsConflictPool.cpp:62–75). Returns `false` for a duplicate,
    /// an empty or an over-long conflict.
    pub fn add(&mut self, mut lits: Vec<BoundLit>, integer: &[bool], feastol: f64) -> bool {
        if lits.is_empty() || lits.len() > self.max_len {
            return false;
        }
        for l in &mut lits {
            if !integer[l.col as usize] {
                l.val += if l.upper { -feastol } else { feastol };
            }
        }
        lits.sort_by(|a, b| {
            (a.col, a.upper)
                .cmp(&(b.col, b.upper))
                .then(a.val.total_cmp(&b.val))
        });
        let mut key = 0xcbf2_9ce4_8422_2325_u64;
        for l in &lits {
            for w in [
                (u64::from(l.col) << 1) | u64::from(l.upper),
                l.val.to_bits(),
            ] {
                key = (key ^ w).wrapping_mul(0x0100_0000_01b3);
                key ^= key >> 29;
            }
        }
        if !self.keys.insert(key) {
            return false;
        }
        let last = (lits.len() - 1) as u32;
        let stored = Stored {
            watch: [0, last.min(1)],
            lits,
            age: 0,
            key,
        };
        let i = match self.free.pop() {
            Some(i) => {
                self.slots[i] = Some(stored);
                i
            }
            None => {
                self.slots.push(Some(stored));
                self.watch_links.extend_from_slice(
                    &[WatchLink {
                        prev: NO_WATCH,
                        next: NO_WATCH,
                    }; 2],
                );
                self.slots.len() - 1
            }
        };
        self.link(2 * i as u32);
        self.link(2 * i as u32 + 1);
        self.age_distribution[0] += 1;
        self.live += 1;
        self.marked.push(i as u32);
        true
    }

    /// Links watch `w` into its literal's column-and-side list
    /// (HiGHS `linkWatchedLiteral`, HighsDomain.cpp:88–102).
    fn link(&mut self, w: u32) {
        let c = self.slots[(w >> 1) as usize]
            .as_ref()
            .expect("live slot on link");
        let lit = c.lits[c.watch[(w & 1) as usize] as usize];
        let (col, upper) = (lit.col as usize, lit.upper);
        let head = if upper {
            &mut self.col_upper_watch[col]
        } else {
            &mut self.col_lower_watch[col]
        };
        self.watch_links[w as usize] = WatchLink {
            prev: NO_WATCH,
            next: *head,
        };
        if *head != NO_WATCH {
            self.watch_links[*head as usize].prev = w;
        }
        *head = w;
    }

    /// Removes watch `w` from its list (HiGHS `unlinkWatchedLiteral`,
    /// HighsDomain.cpp:104–123). Both watches of a live slot are always
    /// linked, so this never needs a membership test.
    fn unlink(&mut self, w: u32) {
        let Some(c) = self.slots.get((w >> 1) as usize).and_then(|s| s.as_ref()) else {
            return;
        };
        let lit = c.lits[c.watch[(w & 1) as usize] as usize];
        let col = lit.col as usize;
        let (prev, next) = (
            self.watch_links[w as usize].prev,
            self.watch_links[w as usize].next,
        );
        let head = if lit.upper {
            &mut self.col_upper_watch[col]
        } else {
            &mut self.col_lower_watch[col]
        };
        if prev != NO_WATCH {
            self.watch_links[prev as usize].next = next;
        } else {
            *head = next;
        }
        if next != NO_WATCH {
            self.watch_links[next as usize].prev = prev;
        }
        self.watch_links[w as usize] = WatchLink {
            prev: NO_WATCH,
            next: NO_WATCH,
        };
    }

    /// Moves a conflict's watches, relinking the lists
    /// (HiGHS `propagateConflict` case 2, HighsDomain.cpp:355–367).
    fn relink(&mut self, slot: usize, watch: [u32; 2]) {
        let Some(c) = &self.slots[slot] else {
            return;
        };
        if c.watch == watch {
            return;
        }
        self.unlink(2 * slot as u32);
        self.unlink(2 * slot as u32 + 1);
        self.slots[slot].as_mut().expect("live slot").watch = watch;
        self.link(2 * slot as u32);
        self.link(2 * slot as u32 + 1);
    }

    /// Pushes the conflicts watching a bound of column `j` (one side for a
    /// one-sided move, both when the bounds collapsed): the conflicts whose
    /// watched literal this move can turn true (HiGHS
    /// `updateActivityLbChange`/`updateActivityUbChange`,
    /// HighsDomain.cpp:269–303).
    fn push_watchers(
        &mut self,
        j: usize,
        lower_moved: bool,
        upper_moved: bool,
        work: &mut Vec<u32>,
    ) {
        if lower_moved {
            let mut w = self.col_lower_watch[j];
            while w != NO_WATCH {
                work.push(w >> 1);
                w = self.watch_links[w as usize].next;
            }
        }
        if upper_moved {
            let mut w = self.col_upper_watch[j];
            while w != NO_WATCH {
                work.push(w >> 1);
                w = self.watch_links[w as usize].next;
            }
        }
    }

    fn reset_age(&mut self, i: usize) {
        if let Some(c) = &mut self.slots[i]
            && c.age > 0
        {
            self.age_distribution[c.age as usize] -= 1;
            self.age_distribution[0] += 1;
            c.age = 0;
        }
    }

    /// Ages every conflict by one and drops those past the limit; while the
    /// pool is over its soft limit the limit shrinks (down to 5) by the
    /// oldest age classes (HighsConflictPool.cpp:130–158). Returns how many
    /// were dropped.
    pub fn age(&mut self) -> usize {
        let mut limit = self.age_limit;
        let mut active = self.live;
        while limit > 5 && active > self.soft_limit {
            active -= self.age_distribution[limit as usize];
            limit -= 1;
        }
        let mut dropped = 0;
        for i in 0..self.slots.len() {
            let Some(c) = &mut self.slots[i] else {
                continue;
            };
            self.age_distribution[c.age as usize] -= 1;
            c.age += 1;
            if c.age <= limit {
                self.age_distribution[c.age as usize] += 1;
                continue;
            }
            self.keys.remove(&c.key);
            // HiGHS `conflictDeleted` unlinks both watches
            // (HighsDomain.cpp:173–177).
            self.unlink(2 * i as u32);
            self.unlink(2 * i as u32 + 1);
            self.slots[i] = None;
            self.free.push(i);
            self.live -= 1;
            dropped += 1;
        }
        dropped
    }

    /// Propagates the conflicts `revisit` selects on the node box
    /// (HiGHS `ConflictPoolPropagation::propagateConflict`,
    /// HighsDomain.cpp:305–372, driven by the watch lists :269–303).
    /// Returns `Err(())` if one has all its literals holding (the box is
    /// infeasible); otherwise the columns whose bound it tightened are
    /// appended to `changed`.
    #[allow(clippy::result_unit_err)] // Err = the box is infeasible, no payload
    #[allow(clippy::too_many_arguments)] // box, selectors and sinks of one pass
    pub fn propagate(
        &mut self,
        lower: &mut [f64],
        upper: &mut [f64],
        integer: &[bool],
        feastol: f64,
        revisit: Revisit<'_>,
        changed: &mut Vec<usize>,
        stats: &mut ConflictStats,
    ) -> Result<(), ()> {
        let mut work = std::mem::take(&mut self.work);
        work.clear();
        match revisit {
            Revisit::All => {
                work.extend(
                    (0..self.slots.len() as u32).filter(|&i| self.slots[i as usize].is_some()),
                );
                self.marked.clear();
            }
            Revisit::Columns(cols) => {
                work.append(&mut self.marked);
                for &j in cols {
                    if j < self.col_lower_watch.len() {
                        let mut w = self.col_lower_watch[j];
                        while w != NO_WATCH {
                            work.push(w >> 1);
                            w = self.watch_links[w as usize].next;
                        }
                        let mut w = self.col_upper_watch[j];
                        while w != NO_WATCH {
                            work.push(w >> 1);
                            w = self.watch_links[w as usize].next;
                        }
                    }
                }
                work.sort_unstable();
                work.dedup();
            }
        }
        let mut pos = 0;
        while pos < work.len() {
            let i = work[pos] as usize;
            pos += 1;
            // Scan phase: decide from the box, then apply outside the
            // borrow of the slot.
            enum Action {
                Skip,
                Prune,
                Tighten(BoundLit, [u32; 2]),
                Rewatch([u32; 2]),
            }
            let action = {
                let Some(c) = self.slots.get(i).and_then(|s| s.as_ref()) else {
                    continue;
                };
                let [w0, w1] = c.watch;
                if w0 != w1
                    && !c.lits[w0 as usize].holds(lower, upper)
                    && !c.lits[w1 as usize].holds(lower, upper)
                {
                    Action::Skip
                } else {
                    // Look for two literals that do not hold.
                    let mut open = [u32::MAX; 2];
                    let mut n = 0;
                    for (k, l) in c.lits.iter().enumerate() {
                        if !l.holds(lower, upper) {
                            open[n] = k as u32;
                            n += 1;
                            if n == 2 {
                                break;
                            }
                        }
                    }
                    match n {
                        0 => Action::Prune,
                        1 => {
                            let k = open[0];
                            Action::Tighten(
                                c.lits[k as usize],
                                [k, if k == 0 && c.lits.len() > 1 { 1 } else { 0 }],
                            )
                        }
                        _ => Action::Rewatch(open),
                    }
                }
            };
            match action {
                Action::Skip => {}
                Action::Prune => {
                    stats.pruned += 1;
                    self.reset_age(i);
                    self.work = work;
                    return Err(());
                }
                Action::Tighten(lit, watch) => {
                    self.relink(i, watch);
                    let j = lit.col as usize;
                    let neg = lit.negated(integer[j], feastol);
                    let mut moved = (false, false);
                    if neg.upper {
                        if neg.val < upper[j] {
                            upper[j] = neg.val;
                            moved.1 = true;
                        }
                    } else if neg.val > lower[j] {
                        lower[j] = neg.val;
                        moved.0 = true;
                    }
                    if moved.0 || moved.1 {
                        stats.tightenings += 1;
                        changed.push(j);
                        self.reset_age(i);
                        if lower[j] > upper[j] {
                            if lower[j] - upper[j] > feastol || integer[j] {
                                stats.pruned += 1;
                                self.work = work;
                                return Err(());
                            }
                            // Continuous crossing within tolerance: collapse,
                            // which moves both bounds.
                            let mid = 0.5 * (lower[j] + upper[j]);
                            lower[j] = mid;
                            upper[j] = mid;
                            moved.1 = true;
                        }
                        // Conflicts watching the moved bound may now have a
                        // watched literal holding (HiGHS drives this from
                        // its bound-change events, HighsDomain.cpp:269–303).
                        self.push_watchers(j, moved.0, moved.1, &mut work);
                    }
                }
                Action::Rewatch(watch) => self.relink(i, watch),
            }
        }
        self.work = work;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ex<'a>(
        glo: &'a [f64],
        ghi: &'a [f64],
        integer: &'a [bool],
        z: &'a dyn Fn(usize) -> u32,
    ) -> Explain<'a> {
        Explain {
            global_lower: glo,
            global_upper: ghi,
            integer,
            feastol: 1e-6,
            up_count: z,
            down_count: z,
        }
    }

    #[test]
    fn explanation_keeps_only_the_bounds_the_proof_needs() {
        // x0 + x1 + x2 ≤ 1 over binaries; the node fixes x0 = x1 = 1 (and x2
        // = 1). Two of them already break the row: the literals are the
        // strongest candidates only and are all `x ≥ 1`.
        let proof = Proof {
            inds: vec![0, 1, 2],
            vals: vec![1.0, 1.0, 1.0],
            rhs: 1.0,
        };
        let (glo, ghi) = (vec![0.0; 3], vec![1.0; 3]);
        let integer = vec![true; 3];
        let z = |_| 0;
        let e = ex(&glo, &ghi, &integer, &z);
        let lits = explain_leq(&proof, &[1.0, 1.0, 1.0], &[1.0; 3], &e).unwrap();
        assert_eq!(lits.len(), 2);
        assert!(lits.iter().all(|l| !l.upper && l.val == 1.0));
        // Only one bound moved: nothing explains the gap.
        assert!(explain_leq(&proof, &[1.0, 0.0, 0.0], &[1.0; 3], &e).is_none());
    }

    #[test]
    fn explanation_relaxes_bounds_with_slack() {
        // 2x0 + x1 ≤ 3 with x in [0,10] integer; node x0 ≥ 5, x1 ≥ 5.
        // x0 ≥ 2 alone gives 4 > 3, so x1's literal is dropped and x0's
        // bound relaxes from 5 to 2.
        let proof = Proof {
            inds: vec![0, 1],
            vals: vec![2.0, 1.0],
            rhs: 3.0,
        };
        let (glo, ghi) = (vec![0.0; 2], vec![10.0; 2]);
        let integer = vec![true; 2];
        let z = |_| 0;
        let e = ex(&glo, &ghi, &integer, &z);
        let lits = explain_leq(&proof, &[5.0, 5.0], &[10.0; 2], &e).unwrap();
        assert_eq!(lits.len(), 1);
        assert_eq!((lits[0].col, lits[0].upper, lits[0].val), (0, false, 2.0));
    }

    #[test]
    fn negative_coefficient_gives_an_upper_literal() {
        // -x0 ≤ -3, i.e. x0 ≥ 3; the node has x0 ≤ 2.
        let proof = Proof {
            inds: vec![0],
            vals: vec![-1.0],
            rhs: -3.0,
        };
        let (glo, ghi) = (vec![0.0], vec![10.0]);
        let integer = vec![true];
        let z = |_| 0;
        let e = ex(&glo, &ghi, &integer, &z);
        let lits = explain_leq(&proof, &[0.0], &[2.0], &e).unwrap();
        assert_eq!(lits.len(), 1);
        assert_eq!((lits[0].col, lits[0].upper, lits[0].val), (0, true, 2.0));
    }

    #[test]
    fn pool_propagates_and_prunes() {
        let integer = vec![true; 3];
        let mut pool = ConflictPool::new(3, 3);
        // not (x0 ≥ 1 and x1 ≥ 1 and x2 ≥ 1)
        let lits = (0..3)
            .map(|j| BoundLit {
                col: j,
                upper: false,
                val: 1.0,
            })
            .collect();
        assert!(pool.add(lits, &integer, 1e-6));
        let mut st = ConflictStats::default();
        let mut changed = Vec::new();
        let (mut lo, mut hi) = (vec![1.0, 1.0, 0.0], vec![1.0; 3]);
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::All,
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!((changed.clone(), hi[2], st.tightenings), (vec![2], 0.0, 1));
        let (mut lo, mut hi) = (vec![1.0; 3], vec![1.0; 3]);
        let r = pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::All,
            &mut changed,
            &mut st,
        );
        assert!(r.is_err() && st.pruned == 1);
        // untouched box: nothing happens
        let (mut lo, mut hi) = (vec![0.0; 3], vec![1.0; 3]);
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::All,
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!(hi, vec![1.0; 3]);
    }

    #[test]
    fn column_revisit_only_touches_watched_conflicts() {
        // One conflict on x0, x1, x2. The first call consumes the mark on an
        // untouched box; afterwards a node whose seeds do not name the
        // watched columns must not revisit it even when its box moved, and
        // a node naming x0 (whose watched literal now holds) must.
        let integer = vec![true; 3];
        let mut pool = ConflictPool::new(3, 3);
        let lits = (0..3)
            .map(|j| BoundLit {
                col: j,
                upper: false,
                val: 1.0,
            })
            .collect();
        assert!(pool.add(lits, &integer, 1e-6));
        let mut st = ConflictStats::default();
        let mut changed = Vec::new();
        let (mut lo, mut hi) = (vec![0.0; 3], vec![1.0; 3]);
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::Columns(&[]),
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!((changed.len(), st.tightenings), (0, 0));
        // The box moves x0 and x1 up without the column set naming them.
        lo[0] = 1.0;
        lo[1] = 1.0;
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::Columns(&[]),
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!((changed.len(), hi[2]), (0, 1.0));
        // Naming the watched column x0 (whose literal now holds) revisits.
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::Columns(&[0]),
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!(
            (changed.as_slice(), hi[2], st.tightenings),
            (&[2][..], 0.0, 1)
        );
    }

    #[test]
    fn a_tightening_revisits_the_conflicts_watching_the_moved_bound() {
        // c2 = not (x1 ≤ 0 and x2 ≥ 1) is skipped first (both watches
        // false); c1 = not (x0 ≥ 1 and x1 ≥ 1) then tightens x1 ≤ 0, which
        // makes c2's watched literal x1 ≤ 0 hold, and the watch list
        // revisits c2 to tighten x2 ≤ 0 (HiGHS updateActivityUbChange,
        // HighsDomain.cpp:287–303).
        let integer = vec![true; 3];
        let mut pool = ConflictPool::new(3, 3);
        let lit = |col, upper, val| BoundLit { col, upper, val };
        assert!(pool.add(vec![lit(1, true, 0.0), lit(2, false, 1.0)], &integer, 1e-6));
        assert!(pool.add(vec![lit(0, false, 1.0), lit(1, false, 1.0)], &integer, 1e-6));
        let mut st = ConflictStats::default();
        let mut changed = Vec::new();
        let (mut lo, mut hi) = (vec![1.0, 0.0, 0.0], vec![1.0; 3]);
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::Columns(&[]),
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert_eq!(
            (changed.as_slice(), hi[1], hi[2], st.tightenings),
            (&[1, 2][..], 0.0, 0.0, 2)
        );
    }

    #[test]
    fn pool_filters_duplicates_and_ages_out() {
        let integer = vec![true; 2];
        let mut pool = ConflictPool::new(2, 2);
        let lit = |c| BoundLit {
            col: c,
            upper: false,
            val: 1.0,
        };
        assert!(pool.add(vec![lit(0), lit(1)], &integer, 1e-6));
        assert!(!pool.add(vec![lit(1), lit(0)], &integer, 1e-6));
        assert_eq!(pool.len(), 1);
        let mut dropped = 0;
        for _ in 0..200 {
            dropped += pool.age();
        }
        assert_eq!((dropped, pool.len()), (1, 0));
        // The aged-out conflict left the watch lists empty: revisiting the
        // columns it watched does nothing and does not dangle.
        let mut st = ConflictStats::default();
        let mut changed = Vec::new();
        let (mut lo, mut hi) = (vec![1.0, 1.0], vec![1.0; 2]);
        pool.propagate(
            &mut lo,
            &mut hi,
            &integer,
            1e-6,
            Revisit::Columns(&[0, 1]),
            &mut changed,
            &mut st,
        )
        .unwrap();
        assert!(changed.is_empty());
        assert!(pool.add(vec![lit(0), lit(1)], &integer, 1e-6));
    }

    #[test]
    fn clean_proof_scales_and_moves_tiny_entries() {
        let (glo, ghi) = (vec![0.0; 3], vec![1.0; 3]);
        let integer = vec![false; 3];
        // 4x0 + 4x1 + 1e-9 x2 ≤ 4 -> x0 + x1 ≤ 1 (x2 ≥ 0 contributes 0)
        let p = clean_proof(&[4.0, 4.0, 1e-9], 4.0, &glo, &ghi, &integer, 1e-6).unwrap();
        assert_eq!(p.inds, vec![0, 1]);
        assert!((p.rhs - 1.0).abs() < 1e-12);
        // tiny entry against an infinite root bound
        let inf = vec![f64::NEG_INFINITY; 3];
        assert!(clean_proof(&[4.0, 4.0, 1e-9], 4.0, &inf, &ghi, &integer, 1e-6).is_none());
    }
}
