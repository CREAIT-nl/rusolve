//! rusolve: solve a small pure-integer MIP in reduced lattice coordinates.
//!
//! Not in HiGHS. Substitute `x = x0 + U y` with `U` unimodular, found by LLL reduction of the
//! columns of `[A; I]`, so the constraint columns `A U` are short and nearly orthogonal; `y` is a free
//! integer vector, and the variable bounds of `x` become `n` extra rows. Branching on `y`
//! is branching along lattice directions rather than coordinate axes, which is what shrinks the tree
//! of general-integer systems with equality rows (held-out Aardal-Lenstra style systems: 10 of 17 solved only with the lever). Inequality-only models are declined: the gen-ip win (gen-ip054: 26 k nodes against more than 1.27 M)
//! did not survive held-out validation.
//! Sources: Aardal, Hurkens and Lenstra, "Solving a system of linear Diophantine equations with lower
//! and upper bounds on the variables" (Math. Oper. Res. 2000); Krishnamoorthy and Pataki, "Column
//! basis reduction and decomposition for integer programs" (2006); Lenstra, Lenstra and Lovasz (1982) for
//! the reduction. Applies only to models with at most `MAX_COLS` integer columns, all with a finite
//! integral lower bound, at least one of them a general integer, and at least one equality row. `mip_lattice`.

use crate::lp_data::highs_solve::highs_status_from_highs_model_status;
use crate::lp_data::presolve_support::get_kkt_failures;
use crate::prelude::*;

const MAX_COLS: usize = 200;
const MAX_ROWS: usize = 120;
const MAX_SWAPS: usize = 200_000;
/// Wall-clock backstop for the reduction (the swap budget is the deterministic limit).
const MAX_SECONDS: f64 = 5.0;
/// Incumbent seeding: a plain branch-and-bound probe of this many nodes on the original model, whose solution
/// (a lattice point, since integrality is preserved) is mapped into the lattice coordinates as the child's
/// starting solution. Lattice-coordinate heuristics find poor incumbents (market split m = 6: 7555 against 15
/// in the original coordinates); the original-coordinate heuristics find good ones.
/// Measured negative: 1000 nodes seeded no incumbent the child did not find on the
/// equality systems and cost 1-3 s on the ones it solves; on market splits m = 6 the incumbent stays worse than
/// without the lever (69 against 15); so the probe is off (0 nodes) and the n - m cap stays at 40.
const PROBE_NODES: HighsInt = 0;
/// Apply the reformulation only when LLL shortens the columns of `[A; I]` by at least this much, as the
/// mean over columns of `ln(original norm / reduced norm)` (0.5 = columns 1.65x shorter in the geometric
/// mean). Tuned on inequality-only gen-ip, where the rule is now superseded by the
/// equality requirement; the held-out equality systems have gains of 6.6 and above. Measured: gen-ip016 2.04 and gen-ip054 1.18 win, gen-ip036 0.76,
/// gen-ip002 0.46 and gen-ip021 0.12 are neutral either way. This is a cost-benefit guard in the spirit of
/// the orthogonality defect of a basis (Lenstra, Lenstra and Lovasz 1982; Aardal and Lenstra 2004).
const MIN_GAIN: f64 = 0.5;
/// Kernel-lattice path limits (Aardal-Hurkens-Lenstra form): rows, columns, and dimension `n - m` of the
/// reformulated problem; and the weight on the equality block in the reduced lattice.
const KERNEL_MAX_ROWS: usize = 12;
const KERNEL_MAX_COLS: usize = 100;
/// Held out: m = 6 market splits (n - m = 44) end TL 120 s with incumbents of 7555 and
/// 6197 against 15 and 22 without the lever; m = 5 (n - m = 35) is fine.
const KERNEL_MAX_DIM: usize = 40;
const KERNEL_WEIGHT: f64 = 1000.0;

/// Deep-insertion depth (Schnorr-Euchner 1994): 1 is plain LLL. Deeper insertions did not pay on the held-in
/// systems, so this stays at 1.
const DEEP_DEPTH: usize = 1;

/// Gram-Schmidt of the columns `b`: orthogonal vectors, coefficients `mu[i][j]` (j < i), squared norms.
fn gram_schmidt(b: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>) {
    let n = b.len();
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, w)| p * w).sum::<f64>();
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut mu = vec![vec![0.0; n]; n];
    let mut bb = vec![0.0; n];
    for i in 0..n {
        let mut v = b[i].clone();
        for j in 0..i {
            mu[i][j] = if bb[j] > 1e-300 {
                dot(&b[i], &q[j]) / bb[j]
            } else {
                0.0
            };
            for (vv, qq) in v.iter_mut().zip(&q[j]) {
                *vv -= mu[i][j] * qq;
            }
        }
        bb[i] = dot(&v, &v);
        q.push(v);
    }
    (q, mu, bb)
}

/// LLL (delta 0.99) on the columns of `m` (`rows x n`, column `j` is `m[j]`); returns the unimodular
/// transform `u` (column `j` of `u` is `u[j]`), or `None` when the swap budget, the time cap or the integer
/// range is exceeded. Floating-point LLL with incremental Gram-Schmidt updates on a swap (Cohen, "A Course in
/// Computational Algebraic Number Theory", Alg. 2.6.3), so a swap costs O(n^2 + rows) instead of a full
/// re-orthogonalisation.
fn lll_columns(m: &[Vec<f64>]) -> Option<(Vec<Vec<i64>>, f64)> {
    let n = m.len();
    let t0 = std::time::Instant::now();
    let mut b: Vec<Vec<f64>> = m.to_vec();
    let mut u: Vec<Vec<i64>> = (0..n)
        .map(|j| (0..n).map(|k| i64::from(j == k)).collect())
        .collect();
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, w)| p * w).sum::<f64>();
    let (_, mut mu, mut bb) = gram_schmidt(&b);
    let depth = DEEP_DEPTH;
    let (mut k, mut swaps) = (1usize, 0usize);
    while k < n {
        for j in (0..k).rev() {
            let r = mu[k][j].round();
            if r != 0.0 {
                let ri = r as i64;
                let (bj, uj) = (b[j].clone(), u[j].clone());
                for (x, y) in b[k].iter_mut().zip(&bj) {
                    *x -= r * y;
                }
                for (x, y) in u[k].iter_mut().zip(&uj) {
                    *x = x.checked_sub(ri.checked_mul(*y)?)?;
                    if x.abs() > 1_000_000 {
                        return None;
                    }
                }
                for l in 0..j {
                    let t = mu[j][l];
                    mu[k][l] -= r * t;
                }
                mu[k][j] -= r;
            }
        }
        // Deep insertion (Schnorr-Euchner 1994): the first position i, at most `depth` back, where the projection
        // of b_k orthogonal to b_0..b_{i-1} is shorter than 0.99 |b*_i|; depth 1 is plain LLL.
        let mut cproj = bb[k] + (0..k).map(|j| mu[k][j] * mu[k][j] * bb[j]).sum::<f64>();
        let mut ins = None;
        for i in 0..k {
            if i + depth >= k && cproj < 0.99 * bb[i] {
                ins = Some(i);
                break;
            }
            cproj -= mu[k][i] * mu[k][i] * bb[i];
        }
        match ins {
            None => k += 1,
            Some(i) => {
                swaps += 1;
                if swaps > MAX_SWAPS
                    || (swaps % 64 == 0 && t0.elapsed().as_secs_f64() > MAX_SECONDS)
                {
                    return None;
                }
                if i + 1 == k {
                    b.swap(k - 1, k);
                    u.swap(k - 1, k);
                    let mk = mu[k][k - 1];
                    let bnew = bb[k] + mk * mk * bb[k - 1];
                    if bnew <= 1e-300 {
                        return None;
                    }
                    mu[k][k - 1] = mk * bb[k - 1] / bnew;
                    bb[k] = bb[k - 1] * bb[k] / bnew;
                    bb[k - 1] = bnew;
                    for l in 0..k - 1 {
                        let t = mu[k - 1][l];
                        mu[k - 1][l] = mu[k][l];
                        mu[k][l] = t;
                    }
                    for r in k + 1..n {
                        let t = mu[r][k];
                        mu[r][k] = mu[r][k - 1] - mk * t;
                        mu[r][k - 1] = t + mu[k][k - 1] * mu[r][k];
                    }
                } else {
                    b[i..=k].rotate_right(1);
                    u[i..=k].rotate_right(1);
                    (_, mu, bb) = gram_schmidt(&b);
                }
                k = i.max(1);
            }
        }
    }
    let ln_norm_sum =
        |cols: &[Vec<f64>]| -> f64 { cols.iter().map(|v| dot(v, v).sqrt().ln()).sum() };
    let gain = (ln_norm_sum(m) - ln_norm_sum(&b)) / n as f64;
    Some((u, gain))
}

/// A change of variables `x = x0 + sum_k basis[k] * y_k` with integral `basis` vectors (n entries each).
struct LatticePlan {
    x0: Vec<f64>,
    basis: Vec<Vec<i64>>,
    /// the original rows are identically satisfied (equality system solved by the kernel form): drop them
    drop_rows: bool,
    /// columns left as they are (continuous columns of a mixed-integer model): no substitution, no bound row
    passthrough: Vec<usize>,
    note: String,
}

impl Highs {
    /// Reduced-basis plan for small pure-integer models (see the module header); `None` when the model
    /// does not qualify or LLL does not shorten the columns enough.
    fn lattice_plan_reduction(&self) -> Option<LatticePlan> {
        let lp = &self.model_lp_;
        let (n, m) = (lp.num_col_ as usize, lp.num_row_ as usize);
        if lp.integrality_.len() != n || m > MAX_ROWS {
            return None;
        }
        let ints: Vec<usize> = (0..n)
            .filter(|&j| lp.integrality_[j] == HighsVarType::kInteger)
            .collect();
        if !(2..=MAX_COLS).contains(&ints.len()) {
            return None;
        }
        let mut general = false;
        for &j in &ints {
            let l = lp.col_lower_[j];
            if !l.is_finite() || l.abs() > 1e9 || l != l.round() {
                return None;
            }
            general |= lp.col_upper_[j] - l > 1.5;
        }
        if !general {
            return None;
        }
        let a = &lp.a_matrix_;
        // Rows with a continuous column are left out of the reduction (the continuous column absorbs them);
        // the rows kept are the pure-integer ones, all rows for a pure-integer model.
        let mut pure = vec![true; m];
        let mut is_int = vec![false; n];
        for &j in &ints {
            is_int[j] = true;
        }
        for j in (0..n).filter(|&j| !is_int[j]) {
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                pure[a.index_[p] as usize] = false;
            }
        }
        // Inequality-only models are declined: the substitution turns every bound into a dense row and loses
        // the cuts, reduced-cost fixing and propagation of the original, and the column-norm gain does not
        // predict a win there. Held-out: 0 wins and 17 of 17 losses (0.1-80 s to the 120 s limit) on random
        // multi-knapsacks and mixed-sign dense systems, integer or real coefficients, against 14 of 14
        // wins on equality systems. Equality rows define a lattice (thin fibre).
        if !(0..m).any(|i| pure[i] && lp.row_lower_[i] == lp.row_upper_[i]) {
            return None;
        }
        let mut rpos = vec![usize::MAX; m];
        let mut np = 0;
        for i in (0..m).filter(|&i| pure[i]) {
            rpos[i] = np;
            np += 1;
        }
        let ni = ints.len();
        // dense columns of [A_pure; I]
        let mut cols: Vec<Vec<f64>> = vec![vec![0.0; np + ni]; ni];
        for (ci, &j) in ints.iter().enumerate() {
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                let i = a.index_[p] as usize;
                if pure[i] {
                    cols[ci][rpos[i]] = a.value_[p];
                }
            }
            cols[ci][np + ci] = 1.0;
        }
        let Some((u, gain)) = lll_columns(&cols) else {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Lattice reduction of {ni} columns stopped by its budget: declined\n"
            );
            return None;
        };
        if gain < MIN_GAIN {
            crate::highs_log_user!(
                &self.options_.log_options,
                HighsLogType::kInfo,
                "Lattice reduction declined: column-norm gain {gain:.2} below {MIN_GAIN}\n"
            );
            return None;
        }
        // new column k = sum_r u[k][r] * old column r, so x_r = x0_r + sum_k u[k][r] y_k; expand to all n columns
        let basis: Vec<Vec<i64>> = u
            .iter()
            .map(|uk| {
                let mut v = vec![0i64; n];
                for (ci, &j) in ints.iter().enumerate() {
                    v[j] = uk[ci];
                }
                v
            })
            .collect();
        let mut x0 = vec![0.0; n];
        for &j in &ints {
            x0[j] = lp.col_lower_[j];
        }
        Some(LatticePlan {
            x0,
            basis,
            drop_rows: false,
            passthrough: (0..n).filter(|&j| !is_int[j]).collect(),
            note: format!("column-norm gain {gain:.2}"),
        })
    }

    /// Aardal-Hurkens-Lenstra kernel-lattice plan for an equality system `A x = d` with integral data: see the
    /// module header. Columns are integer with a finite integral lower bound, or continuous slacks that are the
    /// only continuous column of a row and have one entry +-1 (the row then forces them integral).
    fn lattice_plan_kernel(&self) -> Option<LatticePlan> {
        let lp = &self.model_lp_;
        let (n, m) = (lp.num_col_ as usize, lp.num_row_ as usize);
        if m == 0
            || m > KERNEL_MAX_ROWS
            || n > KERNEL_MAX_COLS
            || n <= m + 1
            || n - m > KERNEL_MAX_DIM
        {
            return None;
        }
        let a = &lp.a_matrix_;
        let integral = |v: f64| v.is_finite() && v.abs() < 1e9 && v == v.round();
        for i in 0..m {
            if lp.row_lower_[i] != lp.row_upper_[i] || !integral(lp.row_lower_[i]) {
                return None;
            }
        }
        let is_int =
            |j: usize| lp.integrality_.len() == n && lp.integrality_[j] == HighsVarType::kInteger;
        let mut cont_in_row = vec![0usize; m];
        for j in 0..n {
            let (s, e) = (a.start_[j] as usize, a.start_[j + 1] as usize);
            if !integral(lp.col_lower_[j])
                || !(lp.col_upper_[j] == kHighsInf || integral(lp.col_upper_[j]))
            {
                return None;
            }
            if (s..e).any(|p| !integral(a.value_[p])) {
                return None;
            }
            if !is_int(j) {
                // continuous: one entry, +-1
                if e - s != 1 || a.value_[s].abs() != 1.0 {
                    return None;
                }
                cont_in_row[a.index_[s] as usize] += 1;
            }
        }
        if cont_in_row.iter().any(|&cnt| cnt > 1) {
            return None;
        }
        // kernel lattice of [A | -d]: columns (e_j ; N a_j) and (e_{n} ; -N d)
        let nn = n + 1;
        let weight = KERNEL_WEIGHT;
        let mut cols: Vec<Vec<f64>> = vec![vec![0.0; nn + m]; nn];
        for j in 0..n {
            cols[j][j] = 1.0;
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                cols[j][nn + a.index_[p] as usize] = weight * a.value_[p];
            }
        }
        cols[n][n] = 1.0;
        for i in 0..m {
            cols[n][nn + i] = -weight * lp.row_lower_[i];
        }
        let (u, _) = lll_columns(&cols)?;
        // reduced column k has a zero lower block iff it lies in the kernel of [A | -d]
        let mut ker: Vec<Vec<i64>> = Vec::new();
        for k in 0..nn {
            let mut low = vec![0.0; m];
            for (r, &coef) in u[k].iter().enumerate() {
                if coef != 0 {
                    for (i, l) in low.iter_mut().enumerate() {
                        *l += coef as f64 * cols[r][nn + i];
                    }
                }
            }
            if low.iter().all(|v| v.abs() < 1e-6) {
                ker.push(u[k].clone());
            }
        }
        if ker.len() != nn - m {
            return None; // rank deficient rows or a numerically failed reduction
        }
        // make the last coordinate (y) vanish on all but one vector, by gcd steps on the y column
        loop {
            let mut nz: Vec<usize> = (0..ker.len()).filter(|&k| ker[k][n] != 0).collect();
            if nz.len() <= 1 {
                break;
            }
            nz.sort_by_key(|&k| ker[k][n].abs());
            let p = ker[nz[0]].clone();
            for &k in &nz[1..] {
                let q = ker[k][n] / p[n];
                for (x, pv) in ker[k].iter_mut().zip(&p) {
                    *x = x.checked_sub(q.checked_mul(*pv)?)?;
                }
            }
        }
        let pi = (0..ker.len()).find(|&k| ker[k][n] != 0)?;
        if ker[pi][n].abs() != 1 {
            return None; // no integral solution of the equality system
        }
        let sign = ker[pi][n];
        let x0: Vec<f64> = ker[pi][..n].iter().map(|&v| (v * sign) as f64).collect();
        let basis: Vec<Vec<i64>> = (0..ker.len())
            .filter(|&k| k != pi)
            .map(|k| ker[k][..n].to_vec())
            .collect();
        if basis.iter().flatten().any(|v| v.abs() > 1_000_000) {
            return None;
        }
        // exact check: A x0 = d and A q = 0
        let mut ax = vec![0.0; m];
        for j in 0..n {
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                ax[a.index_[p] as usize] += a.value_[p] * x0[j];
            }
        }
        if (0..m).any(|i| ax[i] != lp.row_lower_[i]) {
            return None;
        }
        Some(LatticePlan {
            x0,
            basis,
            drop_rows: true,
            passthrough: Vec::new(),
            note: format!(
                "kernel lattice of {m} equality rows, dimension {}",
                nn - m - 1
            ),
        })
    }

    /// Plain probe on the original model; its solution mapped into lattice coordinates (`None` when the probe
    /// has no solution, or the mapping is not exact).
    fn lattice_probe_start(
        &self,
        opts: &HighsOptions,
        x0: &[f64],
        u: &[Vec<i64>],
        passthrough: &[usize],
    ) -> Option<Vec<f64>> {
        if PROBE_NODES == 0 {
            return None;
        }
        let n = self.model_lp_.num_col_ as usize;
        let k = u.len();
        let mut po = opts.clone();
        po.output_flag = false;
        po.log_to_console = false;
        po.mip_max_nodes = PROBE_NODES;
        if po.time_limit.is_finite() {
            po.time_limit = (0.1 * (po.time_limit - self.timer_.read())).max(0.0);
        }
        let mut probe = Highs::default();
        probe.pass_options(&po);
        probe.pass_model(self.model_lp_.clone());
        probe.run();
        if !probe.solution_.value_valid || probe.solution_.col_value.len() != n {
            return None;
        }
        let x = probe.solution_.col_value.clone();
        let mut is_pass = vec![false; n];
        for &j in passthrough {
            is_pass[j] = true;
        }
        let free: Vec<usize> = (0..n).filter(|&r| !is_pass[r]).collect();
        // solve (U U^T) y = U (x - x0) restricted to the substituted columns, k x k Gauss with pivoting
        let mut g = vec![vec![0.0; k + 1]; k];
        for a in 0..k {
            for b in 0..k {
                g[a][b] = free.iter().map(|&r| u[a][r] as f64 * u[b][r] as f64).sum();
            }
            g[a][k] = free
                .iter()
                .map(|&r| u[a][r] as f64 * (x[r].round() - x0[r]))
                .sum();
        }
        for col in 0..k {
            let piv = (col..k).max_by(|&i, &j| g[i][col].abs().total_cmp(&g[j][col].abs()))?;
            if g[piv][col].abs() < 1e-9 {
                return None;
            }
            g.swap(col, piv);
            for r in 0..k {
                if r != col {
                    let f = g[r][col] / g[col][col];
                    if f != 0.0 {
                        for c2 in col..=k {
                            let t = g[col][c2];
                            g[r][c2] -= f * t;
                        }
                    }
                }
            }
        }
        let y: Vec<f64> = (0..k).map(|a| (g[a][k] / g[a][a]).round()).collect();
        for &r in &free {
            let v = x0[r] + (0..k).map(|a| u[a][r] as f64 * y[a]).sum::<f64>();
            if (v - x[r].round()).abs() > 1e-6 {
                return None; // not a lattice point (e.g. a slack outside the kernel form)
            }
        }
        let mut start = y;
        start.extend(passthrough.iter().map(|&j| x[j]));
        Some(start)
    }

    /// Solve the lattice-reformulated MIP and map the solution back. `None` when neither lever is on or
    /// the model does not qualify (the caller then runs the usual path).
    pub(crate) fn call_solve_mip_lattice(&mut self) -> Option<HighsStatus> {
        let lp = &self.model_lp_;
        if (!self.options_.mip_lattice && !self.options_.mip_lattice_kernel)
            || self.solution_.value_valid
            || lp.has_semi_variables()
            || lp.a_matrix_.format_ != MatrixFormat::kColwise
        {
            return None;
        }
        let plan = self
            .options_
            .mip_lattice_kernel
            .then(|| self.lattice_plan_kernel())
            .flatten()
            .or_else(|| {
                self.options_
                    .mip_lattice
                    .then(|| self.lattice_plan_reduction())
                    .flatten()
            })?;
        let LatticePlan {
            x0,
            basis: u,
            drop_rows,
            passthrough,
            note,
        } = plan;
        let lp = &self.model_lp_;
        let (n, m) = (lp.num_col_ as usize, lp.num_row_ as usize);
        let k = u.len();
        let mrows = if drop_rows { 0 } else { m };
        let a = &lp.a_matrix_;
        let mut is_pass = vec![false; n];
        for &j in &passthrough {
            is_pass[j] = true;
        }
        let brows: Vec<usize> = (0..n).filter(|&r| !is_pass[r]).collect();
        let nb = brows.len();
        let np = passthrough.len();
        let mut sub = HighsLp::default();
        sub.num_col_ = (k + np) as HighsInt;
        sub.num_row_ = (mrows + nb) as HighsInt;
        sub.sense_ = lp.sense_;
        sub.offset_ = lp.offset_ + (0..n).map(|j| lp.col_cost_[j] * x0[j]).sum::<f64>();
        sub.col_cost_ = (0..k)
            .map(|t| (0..n).map(|r| lp.col_cost_[r] * u[t][r] as f64).sum())
            .chain(passthrough.iter().map(|&j| lp.col_cost_[j]))
            .collect();
        sub.col_lower_ = vec![-kHighsInf; k];
        sub.col_upper_ = vec![kHighsInf; k];
        sub.integrality_ = vec![HighsVarType::kInteger; k];
        for &j in &passthrough {
            sub.col_lower_.push(lp.col_lower_[j]);
            sub.col_upper_.push(lp.col_upper_[j]);
            sub.integrality_.push(HighsVarType::kContinuous);
        }
        // row activity at x0
        let mut ax0 = vec![0.0; m];
        for r in 0..n {
            for p in a.start_[r] as usize..a.start_[r + 1] as usize {
                ax0[a.index_[p] as usize] += a.value_[p] * x0[r];
            }
        }
        sub.row_lower_ = (0..mrows)
            .map(|i| lp.row_lower_[i] - ax0[i])
            .chain(brows.iter().map(|&r| lp.col_lower_[r] - x0[r]))
            .collect();
        sub.row_upper_ = (0..mrows)
            .map(|i| lp.row_upper_[i] - ax0[i])
            .chain(brows.iter().map(|&r| lp.col_upper_[r] - x0[r]))
            .collect();
        sub.a_matrix_.format_ = MatrixFormat::kColwise;
        sub.a_matrix_.num_col_ = (k + np) as HighsInt;
        sub.a_matrix_.num_row_ = (mrows + nb) as HighsInt;
        sub.a_matrix_.start_ = vec![0];
        for t in 0..k {
            if !drop_rows {
                // rows 0..m: (A U)[i][t] = sum_r A[i][r] u[t][r]
                let mut au = vec![0.0; m];
                for r in 0..n {
                    if u[t][r] == 0 {
                        continue;
                    }
                    for p in a.start_[r] as usize..a.start_[r + 1] as usize {
                        au[a.index_[p] as usize] += a.value_[p] * u[t][r] as f64;
                    }
                }
                for (i, v) in au.iter().enumerate() {
                    if v.abs() > 1e-13 {
                        sub.a_matrix_.index_.push(i as HighsInt);
                        sub.a_matrix_.value_.push(*v);
                    }
                }
            }
            // rows mrows..mrows+n: bound row of x_r
            for (bi, &r) in brows.iter().enumerate() {
                if u[t][r] != 0 {
                    sub.a_matrix_.index_.push((mrows + bi) as HighsInt);
                    sub.a_matrix_.value_.push(u[t][r] as f64);
                }
            }
            sub.a_matrix_
                .start_
                .push(sub.a_matrix_.index_.len() as HighsInt);
        }
        for &j in &passthrough {
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                sub.a_matrix_.index_.push(a.index_[p]);
                sub.a_matrix_.value_.push(a.value_[p]);
            }
            sub.a_matrix_
                .start_
                .push(sub.a_matrix_.index_.len() as HighsInt);
        }
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "MIP solved in lattice coordinates ({} columns, {} rows, {})\n",
            k,
            mrows + nb,
            note
        );
        let mut opts = self.options_.clone();
        opts.mip_lattice = false;
        opts.mip_lattice_kernel = false;
        if opts.time_limit.is_finite() {
            opts.time_limit = (opts.time_limit - self.timer_.read()).max(0.0);
        }
        let mut child = Highs::default();
        child.pass_options(&opts);
        child.pass_model(sub);
        if let Some(start) = self.lattice_probe_start(&opts, &x0, &u, &passthrough) {
            child.set_solution(&HighsSolution {
                col_value: start,
                ..HighsSolution::default()
            });
        }
        if opts.time_limit.is_finite() {
            // the probe used part of the budget
            child.options_.time_limit = (self.options_.time_limit - self.timer_.read()).max(0.0);
        }
        child.run();
        let st = child.get_model_status();
        let info = child.get_info();
        let (nodes, iters) = (
            info.mip_node_count,
            i64::from(info.simplex_iteration_count.max(0)),
        );
        let bound = info.mip_dual_bound;
        let has_sol = child.solution_.value_valid && info.objective_function_value.is_finite();
        let mut obj = kHighsInf;
        if has_sol {
            let y: Vec<f64> = child.solution_.col_value[..k]
                .iter()
                .map(|v| v.round())
                .collect();
            let mut pi = 0;
            let x: Vec<f64> = (0..n)
                .map(|r| {
                    if is_pass[r] {
                        pi += 1;
                        child.solution_.col_value[k + pi - 1]
                    } else {
                        x0[r] + (0..k).map(|t| u[t][r] as f64 * y[t]).sum::<f64>()
                    }
                })
                .collect();
            obj = self.model_lp_.offset_
                + (0..n)
                    .map(|r| self.model_lp_.col_cost_[r] * x[r])
                    .sum::<f64>();
            self.solution_.col_value = x;
            self.model_lp_.a_matrix_.product_quad(
                &mut self.solution_.row_value,
                &self.solution_.col_value,
                kDebugReportOff,
            );
            self.solution_.value_valid = true;
        }
        self.model_status_ = st;
        let mut return_status = highs_status_from_highs_model_status(self.model_status_);
        self.info_.objective_function_value = obj;
        let ptol = self.options_.primal_feasibility_tolerance;
        self.options_.primal_feasibility_tolerance = self.options_.mip_feasibility_tolerance;
        get_kkt_failures(
            &self.options_,
            false,
            &self.model_lp_,
            &self.model_lp_.col_cost_,
            &self.solution_,
            &mut self.info_,
            false,
        );
        self.options_.primal_feasibility_tolerance = ptol;
        self.info_.mip_node_count = nodes;
        self.info_.mip_dual_bound = bound;
        self.info_.mip_gap = if has_sol {
            (obj - bound).abs() / obj.abs().max(1.0)
        } else {
            kHighsInf
        };
        self.info_.simplex_iteration_count = iters.min(i64::from(kHighsIInf)) as HighsInt;
        self.info_.valid = true;
        if self.model_status_ == HighsModelStatus::kOptimal {
            return_status = self.check_optimality("MIP");
        }
        Some(return_status)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(lattice: bool) -> Highs {
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_bool("mip_lattice", lattice);
        for cost in [7.0, 5.0, 9.0] {
            h.add_col(cost, 0.0, kHighsInf, 0, &[], &[]);
        }
        for j in 0..3 {
            h.change_col_integrality(j, HighsVarType::kInteger);
        }
        // nearly parallel columns: the reduced basis is much shorter
        h.add_row(1234.5, kHighsInf, 3, &[0, 1, 2], &[100.0, 99.0, 3.0]);
        h.add_row(4.2, kHighsInf, 3, &[0, 1, 2], &[1.0, 1.0, 1.3]);
        // the reduction path needs an equality row (inequality-only models are declined)
        h.add_row(1396.0, 1396.0, 3, &[0, 1, 2], &[99.0, 98.0, 4.0]);
        h
    }

    #[test]
    fn mixed_integer_model_keeps_continuous_columns() {
        let build = |lattice: bool| {
            let mut h = model(lattice);
            // a continuous column in a fourth row only (that row is outside the reduction)
            h.add_col(1.0, 0.0, 50.0, 0, &[], &[]);
            h.add_row(3.5, kHighsInf, 2, &[0, 3], &[1.0, 1.0]);
            h
        };
        let mut plain = build(false);
        plain.run();
        let mut lat = build(true);
        assert!(
            lat.lattice_plan_reduction()
                .is_some_and(|p| p.passthrough == vec![3])
        );
        lat.run();
        assert_eq!(lat.get_model_status(), HighsModelStatus::kOptimal);
        let (a, b) = (
            plain.get_info().objective_function_value,
            lat.get_info().objective_function_value,
        );
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
        assert!(
            lat.solution_.col_value[..3]
                .iter()
                .all(|v| (v - v.round()).abs() < 1e-9)
        );
    }

    #[test]
    fn inequality_only_models_are_declined() {
        let mut h = model(true);
        assert!(h.lattice_plan_reduction().is_some());
        h.change_row_bounds(2, 1396.0, kHighsInf);
        assert!(h.lattice_plan_reduction().is_none());
    }

    #[test]
    fn lattice_coordinates_give_the_same_optimum() {
        let mut plain = model(false);
        plain.run();
        let mut lat = model(true);
        lat.run();
        assert_eq!(plain.get_model_status(), HighsModelStatus::kOptimal);
        assert_eq!(lat.get_model_status(), HighsModelStatus::kOptimal);
        let (a, b) = (
            plain.get_info().objective_function_value,
            lat.get_info().objective_function_value,
        );
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
        let x = &lat.solution_.col_value;
        assert!(
            x.iter()
                .all(|v| *v >= -1e-9 && (v - v.round()).abs() < 1e-9)
        );
    }

    #[test]
    fn lll_transform_is_unimodular_and_reduces() {
        let cols = vec![vec![1.0, 100.0, 1.0, 0.0], vec![0.0, 99.0, 0.0, 1.0]];
        let (u, gain) = lll_columns(&cols).unwrap();
        assert_eq!((u[0][0] * u[1][1] - u[0][1] * u[1][0]).abs(), 1);
        assert!(gain > MIN_GAIN, "{gain}");
        // already reduced columns: no gain, the guard declines
        let id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(lll_columns(&id).unwrap().1 < MIN_GAIN);
    }

    /// Market-split style equality knapsack with one-sided slacks: the kernel form gives the same optimum.
    fn market(kernel: bool) -> Highs {
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_bool("mip_lattice_kernel", kernel);
        let a = [
            [7.0, 12.0, 5.0, 9.0, 3.0, 11.0],
            [4.0, 8.0, 13.0, 2.0, 6.0, 10.0],
        ];
        for _ in 0..6 {
            h.add_col(0.0, 0.0, 1.0, 0, &[], &[]);
        }
        for _ in 0..2 {
            h.add_col(1.0, 0.0, kHighsInf, 0, &[], &[]);
        }
        for j in 0..6 {
            h.change_col_integrality(j, HighsVarType::kInteger);
        }
        for (r, row) in a.iter().enumerate() {
            let mut idx: Vec<HighsInt> = (0..6).collect();
            let mut val: Vec<f64> = row.to_vec();
            idx.push(6 + r as HighsInt);
            val.push(1.0);
            let d = if r == 0 { 23.0 } else { 20.0 };
            h.add_row(d, d, 7, &idx, &val);
        }
        h
    }

    #[test]
    fn kernel_lattice_matches_the_plain_optimum() {
        let mut plain = market(false);
        plain.run();
        let mut ker = market(true);
        ker.model_lp_.a_matrix_.ensure_colwise();
        assert!(
            ker.lattice_plan_kernel().is_some(),
            "guard should admit the model"
        );
        ker.run();
        assert_eq!(plain.get_model_status(), HighsModelStatus::kOptimal);
        assert_eq!(ker.get_model_status(), HighsModelStatus::kOptimal);
        let (a, b) = (
            plain.get_info().objective_function_value,
            ker.get_info().objective_function_value,
        );
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
        let x = &ker.solution_.col_value;
        assert!(
            x[..6]
                .iter()
                .all(|v| (v - v.round()).abs() < 1e-9 && *v > -1e-9 && *v < 1.0 + 1e-9)
        );
    }
}
