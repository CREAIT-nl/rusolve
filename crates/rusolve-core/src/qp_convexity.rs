//! Convexity test for a quadratic form `x'Qx` (convex QP). Objective only; a row `q(x) <= b` is the same test on its own Q.
//!
//! The question is whether `Q` is positive semidefinite (PSD) for a minimisation, or negative
//! semidefinite for a maximisation (the code negates Q and tests PSD). The steps, cheapest first:
//!
//! 1. **Diagonal tests.** A negative diagonal entry, or a zero diagonal entry with a nonzero
//!    off-diagonal entry in its row, makes a 1x1 or 2x2 principal minor negative: not PSD.
//!    (Horn and Johnson, *Matrix Analysis*, 2nd ed., 2013, section 7.1, principal minors.)
//! 2. **Gershgorin accept.** If every row has `Q_ii >= sum_j |Q_ij|` then `Q` is PSD
//!    (Gershgorin circle theorem; Horn and Johnson, section 6.1).
//! 3. **Sparse LDL' with a scaled pivot tolerance.** `Q` is scaled to a unit diagonal,
//!    `D^-1/2 Q D^-1/2`, which keeps PSD-ness; the elimination follows a minimum-degree order
//!    (George and Liu 1989, "The evolution of the minimum degree ordering algorithm", SIAM Review
//!    31). A pivot above the tolerance is eliminated; a pivot within the tolerance of zero is
//!    accepted only if its remaining column is also within tolerance of zero (a PSD matrix with a
//!    zero pivot has a zero column in the Schur complement); a pivot below minus the tolerance
//!    proves `Q` indefinite (Golub and Van Loan, *Matrix Computations*, 4th ed., section 4.2.8,
//!    Cholesky for semidefinite matrices, without pivoting). When the Schur complement turns
//!    dense it is finished with a dense kernel.
//! 4. **Lanczos on failure.** A failed factorization is followed by a Lanczos estimate of the
//!    smallest eigenvalue with full reorthogonalization (Parlett, *The Symmetric Eigenvalue
//!    Problem*, 1980, chapter 13; Golub and Van Loan section 10.1). The Ritz vector names the
//!    variables the negative curvature lives on. The estimate is an upper bound of the true
//!    smallest eigenvalue.
//! 5. **Rounding shift.** If the estimate is within [`ROUNDING_TOL`] of zero relative to the
//!    largest `|Q_ij|`, and `Q + shift*I` factorizes, the form is accepted as convex after the
//!    reported shift (a Tikhonov-style diagonal shift; Nocedal and Wright, *Numerical
//!    Optimization*, 2nd ed., 2006, section 3.4, modification by a multiple of the identity).

use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};

/// A pivot (of the unit-diagonal scaled matrix) with `|d| <= ZERO_PIVOT` is a zero pivot.
pub const ZERO_PIVOT: f64 = 1e-10;
/// A zero pivot's remaining column must have every entry within this (scaled) bound. A PSD
/// matrix has `|a_pi|^2 <= d a_ii <= ZERO_PIVOT`, so `1e-5` plus a factor 2 for rounding.
const ZERO_COLUMN: f64 = 2e-5;
/// A smallest eigenvalue above `-ROUNDING_TOL * max|Q_ij|` is treated as rounding noise.
pub const ROUNDING_TOL: f64 = 1e-8;
/// Work cap of the sparse factorization (multiply-adds); past it the verdict is `Unverified`.
const WORK_CAP: f64 = 4e9;
/// Switch to the dense kernel when the remaining active matrix is this dense and small enough.
const DENSE_SWITCH_DENSITY: f64 = 0.3;
const DENSE_SWITCH_MAX: usize = 4000;

/// The verdict of [`analyze`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// PSD (NSD for a maximisation) within the pivot tolerance.
    Convex,
    /// Indefinite by no more than rounding; convex after adding `shift` to the diagonal.
    ConvexWithShift,
    /// Indefinite beyond rounding: the smallest eigenvalue estimate and witnesses are set.
    Nonconvex,
    /// The factorization hit its work cap: neither proven nor refuted.
    Unverified,
}

/// How the verdict was reached.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// No nonzero entries: the objective is linear.
    Empty,
    /// A negative diagonal or a zero diagonal with off-diagonal entries (step 1).
    Diagonal,
    /// Diagonal dominance (step 2).
    Gershgorin,
    /// Sparse LDL' (step 3).
    Ldlt,
    /// LDL' failed, Lanczos estimated the smallest eigenvalue (steps 4 and 5).
    Lanczos,
}

/// The result of [`analyze`]. Indices are the caller's variable indices.
#[derive(Debug, Clone)]
pub struct Convexity {
    pub verdict: Verdict,
    pub method: Method,
    /// Added to every diagonal entry of the (sign-adjusted) Q to make it PSD; 0 unless
    /// [`Verdict::ConvexWithShift`]. For a maximisation the caller subtracts it.
    pub shift: f64,
    /// An estimate of the smallest eigenvalue of the minimisation-oriented Q (`-Q` for a
    /// maximisation), when a failed test produced one. Exact for a diagonal entry (step 1).
    pub min_eig_estimate: Option<f64>,
    /// Variables the negative curvature lives on, with their weight in the unit eigenvector
    /// estimate (largest first, at most 8). Empty when convex.
    pub witness: Vec<(usize, f64)>,
    /// Nonzeros of the symmetrised Q (both triangles).
    pub nnz: usize,
    /// True if the Hessian is diagonal.
    pub diagonal: bool,
    /// `max |Q_ij|`.
    pub scale: f64,
}

/// Sparse symmetric matrix, both triangles, row-wise, sorted.
struct Sym {
    n: usize,
    rows: Vec<Vec<(usize, f64)>>,
}

impl Sym {
    /// Symmetrise `entries` as `(Q + Q')/2`, summing duplicates, negating when `negate`.
    fn build(n: usize, entries: &[(usize, usize, f64)], negate: bool) -> Sym {
        let mut raw: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for &(i, j, v) in entries {
            if v != 0.0 {
                *raw.entry((i, j)).or_insert(0.0) += if negate { -v } else { v };
            }
        }
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut seen: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for (&(i, j), &v) in &raw {
            *seen.entry((i, j)).or_insert(0.0) += v / 2.0;
            *seen.entry((j, i)).or_insert(0.0) += v / 2.0;
        }
        for (&(i, j), &v) in &seen {
            if v != 0.0 {
                rows[i].push((j, v));
            }
        }
        Sym { n, rows }
    }

    fn nnz(&self) -> usize {
        self.rows.iter().map(Vec::len).sum()
    }

    fn diag(&self, i: usize) -> f64 {
        self.rows[i]
            .iter()
            .find(|&&(j, _)| j == i)
            .map_or(0.0, |&(_, v)| v)
    }

    fn matvec(&self, x: &[f64], y: &mut [f64], shift: f64) {
        for i in 0..self.n {
            let mut s = shift * x[i];
            for &(j, v) in &self.rows[i] {
                s += v * x[j];
            }
            y[i] = s;
        }
    }
}

/// Outcome of a factorization attempt.
enum Fact {
    Psd,
    /// Indefinite; the (scaled) pivot value and the original index of the pivot variable.
    #[allow(dead_code)]
    Fail {
        pivot: f64,
        var: usize,
    },
    CapHit,
}

/// Test `Q` (entries `(row, col, value)` of the full symmetric matrix; `n` variables).
///
/// `maximize`: the objective is maximised, so `-Q` is tested for PSD. Duplicates are summed and an
/// asymmetric input is replaced by `(Q + Q')/2`.
pub fn analyze(n: usize, entries: &[(usize, usize, f64)], maximize: bool) -> Convexity {
    let q = Sym::build(n, entries, maximize);
    let nnz = q.nnz();
    let scale = q
        .rows
        .iter()
        .flatten()
        .fold(0.0_f64, |m, &(_, v)| m.max(v.abs()));
    let diagonal = q
        .rows
        .iter()
        .enumerate()
        .all(|(i, r)| r.iter().all(|&(j, _)| j == i));
    let mut out = Convexity {
        verdict: Verdict::Convex,
        method: Method::Empty,
        shift: 0.0,
        min_eig_estimate: None,
        witness: Vec::new(),
        nnz,
        diagonal,
        scale,
    };
    if nnz == 0 {
        return out;
    }
    // Step 1: principal minors of order 1 and 2.
    let mut worst: Option<(usize, f64)> = None;
    for i in 0..n {
        let d = q.diag(i);
        let off = q.rows[i].iter().any(|&(j, _)| j != i);
        let bad = if d < 0.0 {
            Some(d)
        } else if d == 0.0 && off {
            // [[0, a], [a, d_jj]] has determinant -a^2 < 0, eigenvalue <= -|a| ... use
            // the Gershgorin-free bound: the 2x2 [[0,a],[a,t]] has smallest eigenvalue
            // (t - sqrt(t^2 + 4a^2))/2 <= -|a|/... evaluated below with the largest |a|.
            let (j, a) =
                q.rows[i]
                    .iter()
                    .filter(|&&(j, _)| j != i)
                    .fold(
                        (i, 0.0_f64),
                        |b, &(j, v)| if v.abs() > b.1 { (j, v.abs()) } else { b },
                    );
            let t = q.diag(j);
            Some((t - (t * t + 4.0 * a * a).sqrt()) / 2.0)
        } else {
            None
        };
        if let Some(e) = bad
            && worst.is_none_or(|(_, w)| e < w)
        {
            worst = Some((i, e));
        }
    }
    if let Some((i, e)) = worst {
        // A diagonal entry is exact; a 2x2 minor gives an upper bound on the smallest
        // eigenvalue, so refine it with Lanczos when it is not a plain diagonal entry.
        out.verdict = Verdict::Nonconvex;
        out.method = Method::Diagonal;
        out.min_eig_estimate = Some(e);
        out.witness = vec![(i, 1.0)];
        if q.diag(i) == 0.0 {
            let lz = lanczos_min(&q, 0.0);
            out.min_eig_estimate = Some(lz.0.min(e));
            if !lz.1.is_empty() {
                out.witness = lz.1;
            }
        }
        return out;
    }
    // Step 2: Gershgorin.
    let dominant = (0..n).all(|i| {
        let d = q.diag(i);
        let off: f64 = q.rows[i]
            .iter()
            .filter(|&&(j, _)| j != i)
            .map(|&(_, v)| v.abs())
            .sum();
        d >= off
    });
    if dominant {
        out.method = Method::Gershgorin;
        return out;
    }
    // Step 3: LDL'.
    out.method = Method::Ldlt;
    match factor(&q, 0.0) {
        Fact::Psd => return out,
        Fact::CapHit => {
            out.verdict = Verdict::Unverified;
            return out;
        }
        Fact::Fail { var, .. } => {
            // Step 4: Lanczos.
            out.method = Method::Lanczos;
            let (lam, wit) = lanczos_min(&q, 0.0);
            out.min_eig_estimate = Some(lam);
            out.witness = if wit.is_empty() {
                vec![(var, 1.0)]
            } else {
                wit
            };
            // Step 5: rounding shift, when the eigenvalue estimate says "tiny" and the
            // shifted matrix factorizes (the estimate is only an upper bound; the
            // refactorization is the proof; the failing pivot itself is no gauge, since
            // rounding is amplified by earlier small pivots).
            let tiny = ROUNDING_TOL * scale;
            if lam >= -tiny {
                let shift = (-lam).max(0.0) * 1.01 + ZERO_PIVOT * scale;
                if matches!(factor(&q, shift), Fact::Psd) {
                    out.verdict = Verdict::ConvexWithShift;
                    out.shift = shift;
                    out.witness.clear();
                    return out;
                }
            }
            out.verdict = Verdict::Nonconvex;
        }
    }
    out
}

/// LDL' of `Q + shift*I` scaled to a unit diagonal, in a minimum-degree order.
fn factor(q: &Sym, shift: f64) -> Fact {
    let n = q.n;
    // Active symmetric matrix, off-diagonal entries per node, plus the diagonal.
    let mut diag: Vec<f64> = (0..n).map(|i| q.diag(i) + shift).collect();
    let mut dscale = vec![1.0; n];
    for i in 0..n {
        if diag[i] > 0.0 {
            dscale[i] = 1.0 / diag[i].sqrt();
        }
    }
    let mut adj: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); n];
    for i in 0..n {
        for &(j, v) in &q.rows[i] {
            if j != i {
                adj[i].insert(j, v * dscale[i] * dscale[j]);
            }
        }
        if diag[i] > 0.0 {
            diag[i] = 1.0;
        } else {
            // zero (or, with shift, tiny) diagonal: it has no off-diagonal entries here
            // (step 1 ruled them out), so it is an empty row.
            diag[i] *= dscale[i] * dscale[i];
        }
    }
    let mut heap: BinaryHeap<Reverse<(usize, usize)>> =
        (0..n).map(|i| Reverse((adj[i].len(), i))).collect();
    let mut eliminated = vec![false; n];
    let mut remaining = n;
    let mut edges: usize = adj.iter().map(BTreeMap::len).sum::<usize>();
    let mut work = 0.0_f64;
    while let Some(Reverse((deg, p))) = heap.pop() {
        if eliminated[p] || deg != adj[p].len() {
            continue;
        }
        // Dense switch.
        if remaining <= DENSE_SWITCH_MAX
            && remaining > 1
            && edges as f64 >= DENSE_SWITCH_DENSITY * (remaining as f64) * (remaining as f64)
        {
            let nodes: Vec<usize> = (0..n).filter(|&i| !eliminated[i]).collect();
            return dense_kernel(&nodes, &diag, &adj);
        }
        let d = diag[p];
        let nbr: Vec<(usize, f64)> = adj[p].iter().map(|(&j, &v)| (j, v)).collect();
        eliminated[p] = true;
        remaining -= 1;
        edges -= nbr.len() * 2;
        for &(j, _) in &nbr {
            adj[j].remove(&p);
        }
        adj[p].clear();
        if d < -ZERO_PIVOT {
            return Fact::Fail { pivot: d, var: p };
        }
        if d <= ZERO_PIVOT {
            if nbr.iter().any(|&(_, v)| v.abs() > ZERO_COLUMN) {
                return Fact::Fail {
                    pivot: -ZERO_COLUMN,
                    var: p,
                };
            }
            // tiny column: dropped
            for &(j, _) in &nbr {
                heap.push(Reverse((adj[j].len(), j)));
            }
            continue;
        }
        let k = nbr.len();
        work += (k * k) as f64;
        if work > WORK_CAP {
            return Fact::CapHit;
        }
        for a in 0..k {
            let (i, vi) = nbr[a];
            let f = vi / d;
            diag[i] -= f * vi;
            for &(j, vj) in &nbr[a + 1..] {
                let cur = match adj[i].get(&j) {
                    Some(&c) => c,
                    None => {
                        edges += 2;
                        0.0
                    }
                };
                let upd = cur - f * vj;
                adj[i].insert(j, upd);
                adj[j].insert(i, upd);
            }
        }
        for &(j, _) in &nbr {
            heap.push(Reverse((adj[j].len(), j)));
        }
    }
    Fact::Psd
}

/// Finish the factorization on the remaining nodes with a dense LDL' (no pivoting).
fn dense_kernel(nodes: &[usize], diag: &[f64], adj: &[BTreeMap<usize, f64>]) -> Fact {
    let m = nodes.len();
    let pos: BTreeMap<usize, usize> = nodes.iter().enumerate().map(|(k, &v)| (v, k)).collect();
    let mut a = vec![0.0; m * m];
    for (k, &v) in nodes.iter().enumerate() {
        a[k * m + k] = diag[v];
        for (&j, &x) in &adj[v] {
            a[k * m + pos[&j]] = x;
        }
    }
    for k in 0..m {
        let d = a[k * m + k];
        if d < -ZERO_PIVOT {
            return Fact::Fail {
                pivot: d,
                var: nodes[k],
            };
        }
        if d <= ZERO_PIVOT {
            if (k + 1..m).any(|i| a[i * m + k].abs() > ZERO_COLUMN) {
                return Fact::Fail {
                    pivot: -ZERO_COLUMN,
                    var: nodes[k],
                };
            }
            continue;
        }
        for i in k + 1..m {
            let f = a[i * m + k] / d;
            if f == 0.0 {
                continue;
            }
            for j in k + 1..=i {
                a[i * m + j] -= f * a[j * m + k];
            }
        }
    }
    Fact::Psd
}

/// Smallest eigenvalue estimate of `Q + shift*I` by restarted Lanczos with full
/// reorthogonalization, and the heaviest components of its Ritz vector.
fn lanczos_min(q: &Sym, shift: f64) -> (f64, Vec<(usize, f64)>) {
    let n = q.n;
    let steps = n.min(120);
    let mut v0: Vec<f64> = {
        // deterministic start vector (xorshift64*)
        let mut s = 0x9E37_79B9_7F4A_7C15_u64;
        (0..n)
            .map(|_| {
                s ^= s >> 12;
                s ^= s << 25;
                s ^= s >> 27;
                let r = s.wrapping_mul(0x2545_F491_4F6C_DD1D);
                (r >> 11) as f64 / (1u64 << 53) as f64 - 0.5
            })
            .collect()
    };
    let mut best = (f64::INFINITY, Vec::new());
    for _restart in 0..4 {
        let nrm = norm(&v0);
        if nrm == 0.0 {
            break;
        }
        for x in &mut v0 {
            *x /= nrm;
        }
        let mut basis: Vec<Vec<f64>> = vec![v0.clone()];
        let (mut alpha, mut beta): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
        let mut w = vec![0.0; n];
        for k in 0..steps {
            q.matvec(&basis[k], &mut w, shift);
            let a = dot(&w, &basis[k]);
            alpha.push(a);
            // full reorthogonalization, twice
            for _ in 0..2 {
                for b in &basis {
                    let c = dot(&w, b);
                    for (wi, bi) in w.iter_mut().zip(b) {
                        *wi -= c * bi;
                    }
                }
            }
            let bn = norm(&w);
            if k + 1 == steps || bn < 1e-12 * (a.abs() + 1.0) {
                break;
            }
            beta.push(bn);
            basis.push(w.iter().map(|x| x / bn).collect());
        }
        let m = alpha.len();
        let mut t = vec![0.0; m * m];
        for k in 0..m {
            t[k * m + k] = alpha[k];
            if k + 1 < m {
                t[k * m + k + 1] = beta[k];
                t[(k + 1) * m + k] = beta[k];
            }
        }
        let (evals, evecs) = jacobi_eigen(&mut t, m);
        let (kmin, &lam) = evals
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .expect("m >= 1");
        // Ritz vector
        let mut y = vec![0.0; n];
        for (k, b) in basis.iter().enumerate().take(m) {
            let c = evecs[k * m + kmin];
            for (yi, bi) in y.iter_mut().zip(b) {
                *yi += c * bi;
            }
        }
        let yn = norm(&y);
        if yn > 0.0 {
            for x in &mut y {
                *x /= yn;
            }
        }
        if lam < best.0 {
            best = (lam, y.clone());
        }
        // residual |Q y - lam y|
        q.matvec(&y, &mut w, shift);
        let res = w
            .iter()
            .zip(&y)
            .map(|(a, b)| (a - lam * b).powi(2))
            .sum::<f64>()
            .sqrt();
        if res <= 1e-10 * (lam.abs() + 1.0) || m == n {
            break;
        }
        v0 = y;
    }
    let (lam, y) = best;
    let mut idx: Vec<(usize, f64)> = y
        .iter()
        .enumerate()
        .map(|(i, &c)| (i, c.abs()))
        .filter(|&(_, c)| c > 1e-3)
        .collect();
    idx.sort_by(|a, b| b.1.total_cmp(&a.1));
    idx.truncate(8);
    (lam, idx)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Cyclic Jacobi on a dense symmetric `m x m` matrix: eigenvalues and eigenvectors (columns).
fn jacobi_eigen(a: &mut [f64], m: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0; m * m];
    for i in 0..m {
        v[i * m + i] = 1.0;
    }
    for _sweep in 0..60 {
        let off: f64 = (0..m)
            .flat_map(|i| (0..i).map(move |j| (i, j)))
            .map(|(i, j)| a[i * m + j].powi(2))
            .sum();
        let tot: f64 = a.iter().map(|x| x * x).sum();
        if off <= 1e-30 * tot.max(1e-300) {
            break;
        }
        for p in 0..m {
            for q in p + 1..m {
                let apq = a[p * m + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q * m + q] - a[p * m + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..m {
                    let (akp, akq) = (a[k * m + p], a[k * m + q]);
                    a[k * m + p] = c * akp - s * akq;
                    a[k * m + q] = s * akp + c * akq;
                }
                for k in 0..m {
                    let (apk, aqk) = (a[p * m + k], a[q * m + k]);
                    a[p * m + k] = c * apk - s * aqk;
                    a[q * m + k] = s * apk + c * aqk;
                }
                for k in 0..m {
                    let (vkp, vkq) = (v[k * m + p], v[k * m + q]);
                    v[k * m + p] = c * vkp - s * vkq;
                    v[k * m + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..m).map(|i| a[i * m + i]).collect(), v)
}

/// QPLIB-style class code (Furini et al. 2019, Math. Prog. Comp. 11): objective, variables,
/// constraints. Objective: `L` linear, `D` convex diagonal, `C` convex, `Q` nonconvex (also
/// `U` when the verdict is `Unverified`). Variables: `C` continuous, `B` binary, `I` integer,
/// `M` mixed continuous and binary, `G` general mixed. Constraints: `N` none, `B` bounds only,
/// `L` linear.
pub fn class_code(
    c: &Convexity,
    cont: usize,
    binary: usize,
    other_int: usize,
    rows: usize,
    any_finite_bound: bool,
) -> String {
    let o = match (c.method == Method::Empty, c.verdict, c.diagonal) {
        (true, _, _) => 'L',
        (_, Verdict::Nonconvex, _) => 'Q',
        (_, Verdict::Unverified, _) => 'U',
        (_, _, true) => 'D',
        _ => 'C',
    };
    let v = match (cont > 0, binary > 0, other_int > 0) {
        (true, false, false) => 'C',
        (false, true, false) => 'B',
        (false, false, true) => 'I',
        (true, true, false) => 'M',
        _ => 'G',
    };
    let k = if rows > 0 {
        'L'
    } else if any_finite_bound {
        'B'
    } else {
        'N'
    };
    format!("{o}{v}{k}")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense(a: &[&[f64]]) -> (usize, Vec<(usize, usize, f64)>) {
        let mut e = Vec::new();
        for (i, r) in a.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                if v != 0.0 {
                    e.push((i, j, v));
                }
            }
        }
        (a.len(), e)
    }

    /// `F'F` for a k x n pseudo-random `F`, as full entries (PSD, rank k).
    fn gram(k: usize, n: usize, seed: u64) -> Vec<Vec<f64>> {
        let mut s = seed;
        let mut next = || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            (s.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let f: Vec<Vec<f64>> = (0..k).map(|_| (0..n).map(|_| next()).collect()).collect();
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| (0..k).map(|r| f[r][i] * f[r][j]).sum())
                    .collect()
            })
            .collect()
    }

    fn entries(m: &[Vec<f64>]) -> Vec<(usize, usize, f64)> {
        let mut e = Vec::new();
        for (i, r) in m.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                e.push((i, j, v));
            }
        }
        e
    }

    #[test]
    fn linear_diagonal_and_dominant() {
        assert_eq!(analyze(3, &[], false).method, Method::Empty);
        let c = analyze(2, &[(0, 0, 2.0), (1, 1, 0.0)], false);
        assert_eq!((c.verdict, c.diagonal), (Verdict::Convex, true));
        let (n, e) = dense(&[&[3.0, 1.0], &[1.0, 3.0]]);
        let c = analyze(n, &e, false);
        assert_eq!((c.verdict, c.method), (Verdict::Convex, Method::Gershgorin));
    }

    #[test]
    fn negative_diagonal_and_zero_diagonal_minor() {
        let c = analyze(2, &[(0, 0, 1.0), (1, 1, -2.0)], false);
        assert_eq!(c.verdict, Verdict::Nonconvex);
        assert_eq!(c.witness[0].0, 1);
        assert_eq!(c.min_eig_estimate, Some(-2.0));
        // [[0, 1], [1, 0]]: eigenvalues +-1
        let c = analyze(2, &[(0, 1, 1.0), (1, 0, 1.0)], false);
        assert_eq!(
            (c.verdict, c.method),
            (Verdict::Nonconvex, Method::Diagonal)
        );
        assert!((c.min_eig_estimate.unwrap() + 1.0).abs() < 1e-9);
    }

    #[test]
    fn maximise_flips_the_sign() {
        let e = [(0, 0, -1.0), (1, 1, -4.0)];
        assert_eq!(analyze(2, &e, true).verdict, Verdict::Convex);
        assert_eq!(analyze(2, &e, false).verdict, Verdict::Nonconvex);
    }

    #[test]
    fn indefinite_names_the_variables() {
        // [[1, 2], [2, 1]] plus an unrelated PSD variable 2: eigenvalues 3, -1, 5
        let (n, e) = dense(&[&[1.0, 2.0, 0.0], &[2.0, 1.0, 0.0], &[0.0, 0.0, 5.0]]);
        let c = analyze(n, &e, false);
        assert_eq!(c.verdict, Verdict::Nonconvex);
        assert!((c.min_eig_estimate.unwrap() + 1.0).abs() < 1e-8, "{c:?}");
        let vars: Vec<usize> = c.witness.iter().map(|w| w.0).collect();
        assert!(
            vars.contains(&0) && vars.contains(&1) && !vars.contains(&2),
            "{vars:?}"
        );
    }

    #[test]
    fn low_rank_psd_is_convex_without_shift() {
        let g = gram(3, 40, 7);
        let c = analyze(40, &entries(&g), false);
        assert_eq!(c.verdict, Verdict::Convex, "{c:?}");
        assert_eq!(c.shift, 0.0);
    }

    #[test]
    fn rounding_level_indefinite_gets_a_shift() {
        let mut g = gram(3, 40, 11);
        // a perturbation far below ROUNDING_TOL * scale on one diagonal, below the rank
        for (i, row) in g.iter_mut().enumerate() {
            row[i] -= 1e-12;
        }
        let e = entries(&g);
        let c = analyze(40, &e, false);
        assert_eq!(c.verdict, Verdict::ConvexWithShift, "{c:?}");
        assert!(c.shift > 0.0 && c.shift < 1e-9, "{}", c.shift);
    }

    #[test]
    fn clearly_indefinite_dense_is_refused_with_the_eigenvalue() {
        let mut g = gram(60, 60, 3);
        g[5][5] -= 3.0;
        g[9][9] -= 3.0;
        let e = entries(&g);
        let c = analyze(60, &e, false);
        assert_eq!(c.verdict, Verdict::Nonconvex);
        // compare with the exact smallest eigenvalue by dense Jacobi
        let mut a: Vec<f64> = g.iter().flatten().copied().collect();
        let (ev, _) = jacobi_eigen(&mut a, 60);
        let exact = ev.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(exact < -1.0);
        assert!(
            (c.min_eig_estimate.unwrap() - exact).abs() < 1e-6 * exact.abs(),
            "{} vs {exact}",
            c.min_eig_estimate.unwrap()
        );
    }

    #[test]
    fn sparse_tridiagonal_laplacian_psd_and_perturbed() {
        // 1-D Laplacian with Neumann ends: PSD, singular (constant vector)
        let n = 500;
        let mut e = Vec::new();
        for i in 0..n {
            let d = if i == 0 || i == n - 1 { 1.0 } else { 2.0 };
            e.push((i, i, d));
            if i + 1 < n {
                e.push((i, i + 1, -1.0));
                e.push((i + 1, i, -1.0));
            }
        }
        let c = analyze(n, &e, false);
        assert!(
            matches!(c.verdict, Verdict::Convex | Verdict::ConvexWithShift),
            "{c:?}"
        );
        e.push((250, 250, -0.5));
        let c = analyze(n, &e, false);
        assert_eq!(c.verdict, Verdict::Nonconvex, "{c:?}");
    }

    #[test]
    fn class_codes() {
        let c = analyze(2, &[(0, 0, 1.0), (1, 1, 1.0)], false);
        assert_eq!(class_code(&c, 2, 0, 0, 1, true), "DCL");
        let c = analyze(
            2,
            &[(0, 0, 3.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)],
            false,
        );
        assert_eq!(class_code(&c, 0, 2, 0, 0, true), "CBB");
        let c = analyze(2, &[], false);
        assert_eq!(class_code(&c, 1, 1, 0, 3, true), "LML");
    }
}
