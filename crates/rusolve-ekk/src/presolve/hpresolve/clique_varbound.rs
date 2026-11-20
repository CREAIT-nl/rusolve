//! rusolve (not a HiGHS rule): `mip_presolve_clique_varbound`.
//!
//! Clique-based variable-bound lifting. Variable-bound rows `z + a_j x_j >= L` (or `<= U`) that
//! share the continuous column `z` and the side `L` and whose binaries `x_j` lie in one clique
//! (`sum x <= 1`) are merged into the single row `z + sum_j a_j x_j >= L`. Valid because at most
//! one `x_j` is one: with `x_k = 1` the row reads `z >= L - a_k`, with none it reads `z >= L`,
//! exactly what the originals say together. For `a_j <= 0` (lower side; `a_j >= 0` on the upper
//! side) the merged row implies every original at `x >= 0`, so the originals are dropped. The
//! LP bound improves because the merged row is much tighter than the disaggregated rows.
//!
//! Sources: Achterberg, Bixby, Gu, Rothberg, Weninger (2020), "Presolve reductions in mixed
//! integer programming", INFORMS J. Comput. 32(2) (clique merging of variable bounds); SCIP 10
//! (Apache-2.0) `cons_varbound.c` / clique table handling.
//! Re-implemented from the description, no SCIP code translated.
//!
//! MIP only: the merged row changes row duals, so no dual or basis postsolve is recorded for it.

use super::*;
use crate::prelude::*;
use std::collections::HashMap;

impl HPresolve {
    pub fn clique_varbound_presolve(
        &mut self,
        postsolve_stack: &mut HighsPostsolveStack,
    ) -> HPresolveResult {
        if !self.options.mip_presolve_clique_varbound || self.mipsolver.is_none() {
            return HPresolveResult::kOk;
        }
        let small = self.options.small_matrix_value;
        let is_bin = |h: &Self, cu: usize| {
            h.model.integrality_[cu] != HighsVarType::kContinuous
                && h.model.col_lower_[cu] == 0.0
                && h.model.col_upper_[cu] == 1.0
        };
        // clique rows: unit positive binaries, `sum <= 1` (or `= 1`)
        let mut cliques: Vec<Vec<HighsInt>> = Vec::new();
        // groups of variable-bound rows keyed by (z, side, bound bits)
        type Key = (HighsInt, bool, u64);
        let mut groups: HashMap<Key, Vec<(HighsInt, HighsInt, f64, f64)>> = HashMap::new();
        for row in 0..self.model.num_row_ {
            let ru = row as usize;
            if self.row_deleted[ru] != 0 {
                continue;
            }
            let rs = self.rowsize[ru];
            let (rl, rup) = (self.model.row_lower_[ru], self.model.row_upper_[ru]);
            if rs >= 2 && rup == 1.0 && (rl == 1.0 || rl <= -kHighsInf) {
                let ent: Vec<(HighsInt, f64)> = self
                    .get_row_vector(row)
                    .iter(self)
                    .map(|nz| (nz.index(), nz.value()))
                    .collect();
                if ent
                    .iter()
                    .all(|&(cj, v)| v == 1.0 && is_bin(self, cj as usize))
                {
                    cliques.push(ent.iter().map(|&(cj, _)| cj).collect());
                }
                continue;
            }
            if rs != 2 || (rl.is_finite() == rup.is_finite()) {
                continue;
            }
            let ent: Vec<(HighsInt, f64)> = self
                .get_row_vector(row)
                .iter(self)
                .map(|nz| (nz.index(), nz.value()))
                .collect();
            for (&(zc, cz), &(xc, ax)) in [(&ent[0], &ent[1]), (&ent[1], &ent[0])] {
                let (zu, xu) = (zc as usize, xc as usize);
                if !is_bin(self, xu)
                    || is_bin(self, zu)
                    || self.model.integrality_[zu] != HighsVarType::kContinuous
                    || cz.abs() <= small
                {
                    continue;
                }
                // normalise to z + a x {>=,<=} bound with a unit z coefficient
                let a = ax / cz;
                let (lo, up) = if cz > 0.0 {
                    (rl / cz, rup / cz)
                } else {
                    (rup / cz, rl / cz)
                };
                let (lower_side, bound) = if lo.is_finite() {
                    (true, lo)
                } else {
                    (false, up)
                };
                // lower side z >= L - a x needs a <= 0, upper side z <= U - a x needs a >= 0
                if (lower_side && a > 0.0) || (!lower_side && a < 0.0) {
                    continue;
                }
                groups
                    .entry((zc, lower_side, bound.to_bits()))
                    .or_default()
                    .push((row, xc, a, cz));
            }
        }
        if cliques.is_empty() || groups.is_empty() {
            return HPresolveResult::kOk;
        }
        // clique membership per column
        let mut cliques_of: HashMap<HighsInt, Vec<usize>> = HashMap::new();
        for (k, cq) in cliques.iter().enumerate() {
            for &cj in cq {
                cliques_of.entry(cj).or_default().push(k);
            }
        }
        let mut merged_rows = 0usize;
        let mut dropped_rows = 0usize;
        let mut keys: Vec<Key> = groups.keys().copied().collect();
        keys.sort_unstable();
        for key in keys {
            let members = &groups[&key];
            if members.len() < 2 {
                continue;
            }
            // count members per clique, greedy cover largest first
            let mut cnt: HashMap<usize, usize> = HashMap::new();
            for &(_, xc, _, _) in members {
                if let Some(ks) = cliques_of.get(&xc) {
                    for &k in ks {
                        *cnt.entry(k).or_default() += 1;
                    }
                }
            }
            let mut order: Vec<(usize, usize)> = cnt.into_iter().filter(|&(_, n)| n >= 2).collect();
            order.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            let mut used = vec![false; members.len()];
            let pos_of: HashMap<HighsInt, usize> = members
                .iter()
                .enumerate()
                .map(|(i, &(_, xc, _, _))| (xc, i))
                .collect();
            for (k, _) in order {
                let mut sel: Vec<usize> = cliques[k]
                    .iter()
                    .filter_map(|cj| pos_of.get(cj).copied())
                    .filter(|&i| !used[i])
                    .collect();
                if sel.len() < 2 {
                    continue;
                }
                sel.sort_unstable();
                let (keep_row, _, _, keep_cz) = members[sel[0]];
                if self.row_deleted[keep_row as usize] != 0 {
                    continue;
                }
                for &i in &sel {
                    used[i] = true;
                }
                for &i in &sel[1..] {
                    let (r, xc, a, _) = members[i];
                    if self.row_deleted[r as usize] != 0 || self.col_deleted[xc as usize] != 0 {
                        continue;
                    }
                    // the coefficient lives in the surviving row, then the original goes
                    self.add_to_matrix(keep_row, xc, a * keep_cz);
                    postsolve_stack.redundant_row(r);
                    self.remove_row(r);
                    dropped_rows += 1;
                }
                merged_rows += 1;
            }
        }
        crate::highs_log_user!(
            &self.options.log_options,
            HighsLogType::kInfo,
            "Clique variable bounds: {} merged rows, {} rows dropped\n",
            merged_rows,
            dropped_rows
        );
        HPresolveResult::kOk
    }
}
