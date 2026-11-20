//! rusolve `mip_node_bigm`: node-local coefficient tightening of big-M rows in the node LP.
//!
//! Source: Savelsbergh, "Preprocessing and probing techniques for mixed integer programming
//! problems", ORSA J. Computing 6 (1994) 445-454, coefficient reduction; applied to the local
//! node bounds instead of the global ones (Achterberg, thesis 2007, section 10.2, "local
//! domain"). Measured on our benchmark set (big-M job-shop models).
//!
//! For a one-sided model row `c x <= r` and an unfixed binary z, let `maxact` be the maximum
//! activity over the local bounds.
//! - `c_z > 0`, `d = r - (maxact - c_z) > 0`: replace `c_z -> c_z - d` and `r -> r - d`.
//! - `c_z < 0`, `d = r - c_z - maxact > 0`: replace `c_z -> c_z + d`.
//!
//! The tightened row has the same integer solutions inside the node and a strictly smaller LP
//! relaxation, so the node LP bound rises. The rows are valid only below the node, so the
//! modification is undone before the next node is evaluated, and every consumer that turns the
//! LP matrix into globally valid information (dual proofs, cuts) is gated off while modified.

use crate::prelude::*;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HighsLpRelaxation_LocalMod {
    pub row: HighsInt,
    pub col: HighsInt,
    pub orig_coef: f64,
    pub orig_lower: f64,
    pub orig_upper: f64,
    pub new_coef: f64,
    pub new_lower: f64,
    pub new_upper: f64,
}

impl HighsLpRelaxation {
    /// True while some model row of the LP differs from the global model row.
    pub fn has_local_mods(&self) -> bool {
        !self.local_mods.is_empty()
    }

    /// Undo every local modification in the LP and forget it.
    pub fn restore_local_mods(&mut self) {
        let mods = std::mem::take(&mut self.local_mods);
        if self.nlocal_rows == 0 {
            Self::undo_mods(&mut self.lpsolver, &mods);
        } else {
            let n = self.lpsolver.get_num_row();
            self.lpsolver
                .delete_rows_interval(n - self.nlocal_rows as HighsInt, n - 1);
            let keep = self.lprows.len() - self.nlocal_rows;
            self.lprows.truncate(keep);
            self.nlocal_rows = 0;
            self.status = HighsLpRelaxation_Status::kNotSet;
            self.currentbasisstored = false;
        }
    }

    /// Variant (b): add one tightened copy of each modified row as a node-local LP row, so the
    /// basis survives (a new row extends it) and the row is dropped by `restore_local_mods`.
    fn add_local_rows(&mut self, ms: &HighsMipSolver, mods: &[HighsLpRelaxation_LocalMod]) {
        let mut lower: Vec<f64> = Vec::new();
        let mut upper: Vec<f64> = Vec::new();
        let mut start: Vec<HighsInt> = Vec::new();
        let mut index: Vec<HighsInt> = Vec::new();
        let mut value: Vec<f64> = Vec::new();
        let mut rows: Vec<HighsInt> = Vec::new();
        let mut i = 0;
        while i < mods.len() {
            let r = mods[i].row;
            let mut j = i;
            while j < mods.len() && mods[j].row == r {
                j += 1;
            }
            let s = ms.mipdata_.a_rstart_[r as usize] as usize;
            let e = ms.mipdata_.a_rstart_[r as usize + 1] as usize;
            start.push(index.len() as HighsInt);
            rows.push(r);
            for k in s..e {
                let col = ms.mipdata_.a_rindex_[k];
                let mut v = ms.mipdata_.a_rvalue_[k];
                if let Some(m) = mods[i..j].iter().find(|m| m.col == col) {
                    v = m.new_coef;
                }
                index.push(col);
                value.push(v);
            }
            // the last entry of a row carries its final bounds
            lower.push(mods[j - 1].new_lower);
            upper.push(mods[j - 1].new_upper);
            i = j;
        }
        let n = lower.len();
        let _ = self.lpsolver.add_rows(
            n as HighsInt,
            &lower,
            &upper,
            value.len() as HighsInt,
            &start,
            &index,
            &value,
        );
        for r in rows {
            self.lprows.push(HighsLpRelaxation_LpRow::local(r));
        }
        self.nlocal_rows = n;
        self.status = HighsLpRelaxation_Status::kNotSet;
        self.currentbasisstored = false;
    }

    pub(crate) fn undo_mods(lpsolver: &mut Highs, mods: &[HighsLpRelaxation_LocalMod]) {
        for m in mods.iter().rev() {
            let _ = lpsolver.change_coeff(m.row, m.col, m.orig_coef);
            let _ = lpsolver.change_row_bounds(m.row, m.orig_lower, m.orig_upper);
        }
    }

    /// Tighten the big-M rows against `dom` (the local domain of the node about to be solved).
    /// Returns the number of modified coefficients.
    pub fn apply_local_bigm(&mut self, ms: &HighsMipSolver, dom: &HighsDomain) -> usize {
        let model = mip_model!(ms);
        let nrow = ms.num_row() as usize;
        let mut desired: Vec<HighsLpRelaxation_LocalMod> = Vec::new();
        let tol = ms.mipdata_.feastol;
        let min_ratio = ms.options_mip_.mip_node_bigm_min_ratio;
        let active = ms.options_mip_.mip_node_bigm_active;
        let rows_mode = ms.options_mip_.mip_node_bigm_rows;
        let sol = self.lpsolver.get_solution();
        let have_sol = active > 0.0 && sol.value_valid && sol.row_value.len() >= nrow;
        for r in 0..nrow {
            let lo = model.row_lower_[r];
            let up = model.row_upper_[r];
            // one-sided rows only
            let (sign, rhs) = if up < kHighsInf && lo <= -kHighsInf {
                (1.0, up)
            } else if lo > -kHighsInf && up >= kHighsInf {
                (-1.0, -lo)
            } else {
                continue;
            };
            let s = ms.mipdata_.a_rstart_[r] as usize;
            let e = ms.mipdata_.a_rstart_[r + 1] as usize;
            // maximum activity of `sign * a x` over the local bounds
            let mut maxact = 0.0f64;
            let mut finite = true;
            let mut has_bin = false;
            for k in s..e {
                let cc = sign * ms.mipdata_.a_rvalue_[k];
                let j = ms.mipdata_.a_rindex_[k] as usize;
                let b = if cc > 0.0 {
                    dom.col_upper_[j]
                } else {
                    dom.col_lower_[j]
                };
                if b.abs() >= kHighsInf {
                    finite = false;
                    break;
                }
                maxact += cc * b;
                if !has_bin
                    && ms.is_col_integral(j as HighsInt)
                    && dom.col_lower_[j] == 0.0
                    && dom.col_upper_[j] == 1.0
                {
                    has_bin = true;
                }
            }
            if !finite || !has_bin || maxact <= rhs + tol {
                continue;
            }
            let slack = if have_sol {
                rhs - sign * sol.row_value[r]
            } else {
                0.0
            };
            let mut r_cur = rhs;
            for k in s..e {
                let j = ms.mipdata_.a_rindex_[k] as usize;
                if !(ms.is_col_integral(j as HighsInt)
                    && dom.col_lower_[j] == 0.0
                    && dom.col_upper_[j] == 1.0)
                {
                    continue;
                }
                let cc = sign * ms.mipdata_.a_rvalue_[k];
                let (new_c, d) = if cc > 0.0 {
                    let d = r_cur - (maxact - cc);
                    (cc - d, d)
                } else {
                    let d = r_cur - cc - maxact;
                    (cc + d, d)
                };
                if d <= min_ratio * cc.abs() + 10.0 * tol {
                    continue;
                }
                // variant (a): only rows the parent LP solution has (nearly) active
                if have_sol && slack > active * cc.abs() {
                    continue;
                }
                if cc > 0.0 {
                    maxact -= d;
                    r_cur -= d;
                }
                let (nl, nu) = if sign > 0.0 {
                    (lo, r_cur)
                } else {
                    (-r_cur, up)
                };
                desired.push(HighsLpRelaxation_LocalMod {
                    row: r as HighsInt,
                    col: j as HighsInt,
                    orig_coef: ms.mipdata_.a_rvalue_[k],
                    orig_lower: lo,
                    orig_upper: up,
                    new_coef: sign * new_c,
                    new_lower: nl,
                    new_upper: nu,
                });
            }
        }
        if desired == self.local_mods {
            return desired.len();
        }
        self.restore_local_mods();
        if desired.is_empty() {
            return 0;
        }
        if rows_mode {
            self.add_local_rows(ms, &desired);
            let n = desired.len();
            self.local_mods = desired;
            return n;
        }
        // later entries of the same row carry that row's final right-hand side
        for m in &desired {
            let _ = self.lpsolver.change_coeff(m.row, m.col, m.new_coef);
        }
        for m in &desired {
            let _ = self
                .lpsolver
                .change_row_bounds(m.row, m.new_lower, m.new_upper);
        }
        self.status = HighsLpRelaxation_Status::kNotSet;
        self.currentbasisstored = false;
        self.basischeckpoint = None;
        let n = desired.len();
        self.local_mods = desired;
        n
    }
}
