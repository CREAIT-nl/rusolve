//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! clique covers for the paper's binary-variable rankings
//! (Salvagnin, Roberti, Fischetti, "A fix-propagate-repair heuristic for mixed
//! integer programming", MPC 17:111-139, 2025, Sect. 4.1 and Figs. 2-3).
//!
//! Two constructions:
//!   * `build_clique_cover` is the five-step greedy of Sect. 4.1.
//!   * `cliques2_order` is Fig. 3.
//!
//! Both are pure functions over a clique-table snapshot.

use crate::HighsInt;
use crate::lp_data::hconst::kHighsIInf;
use crate::mip::clique_table::{HighsCliqueTable_Clique, HighsCliqueTable_CliqueVar};

/// mip-heuristics (36f7d9d): clique_cover.h::kCliqueTightnessTol (line 45)
pub(crate) const kCliqueTightnessTol: f64 = 1e-6;

/// mip-heuristics (36f7d9d): clique_cover.h::Cover (lines 49-87)
pub(crate) struct Cover {
    /// Covered binaries, concatenated group by group.
    pub(crate) members: Vec<HighsInt>,
    /// `member_pos[i]` is 1 iff `members[i]` entered its group as a positive literal.
    pub(crate) member_pos: Vec<u8>,
    /// Half-open group bounds into `members`; size is `num_groups + 1`.
    pub(crate) group_start: Vec<HighsInt>,
    /// Binaries no group covers, in formulation order.
    pub(crate) uncovered: Vec<HighsInt>,
    /// The first `num_equality_groups` groups are the equality cliques of step 1.
    pub(crate) num_equality_groups: HighsInt,
}

impl Default for Cover {
    fn default() -> Self {
        Self {
            members: Vec::new(),
            member_pos: Vec::new(),
            group_start: Vec::new(),
            uncovered: Vec::new(),
            num_equality_groups: 0,
        }
    }
}

impl Cover {
    /// mip-heuristics (36f7d9d): clique_cover.h::Cover::num_groups (lines 84-86)
    pub(crate) fn num_groups(&self) -> HighsInt {
        if self.group_start.is_empty() {
            0
        } else {
            self.group_start.len() as HighsInt - 1
        }
    }
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::is_live (lines 13-15)
fn is_live(cl: &HighsCliqueTable_Clique) -> bool {
    cl.start != -1
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::binary_mask (lines 17-23)
fn binary_mask(bin: &[HighsInt], ncol: HighsInt) -> Vec<u8> {
    let mut is_bin = vec![0u8; ncol as usize];
    for &j in bin {
        is_bin[j as usize] = 1;
    }
    is_bin
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::add_equality_groups (lines 42-92)
fn add_equality_groups(
    cliques: &[HighsCliqueTable_Clique],
    entries: &[HighsCliqueTable_CliqueVar],
    is_bin: &[u8],
    covered: &mut [u8],
    stamp: &mut [HighsInt],
    cover: &mut Cover,
) -> HighsInt {
    let mut num_groups: HighsInt = 0;
    let mut picked: Vec<HighsInt> = Vec::new();
    let num_cliques = cliques.len() as HighsInt;

    for c in 0..num_cliques {
        let cl = &cliques[c as usize];
        if !is_live(cl) || !cl.equality {
            continue;
        }

        picked.clear();
        let mut disjoint = true;
        let mut k = cl.start;
        while k < cl.end {
            let j = entries[k as usize].col() as HighsInt;
            if is_bin[j as usize] == 0 {
                k += 1;
                continue;
            }
            if covered[j as usize] != 0 {
                disjoint = false;
                break;
            }
            if stamp[j as usize] == c {
                k += 1;
                continue;
            }
            stamp[j as usize] = c;
            picked.push(k);
            k += 1;
        }
        if !disjoint || picked.is_empty() {
            continue;
        }

        cover.group_start.push(cover.members.len() as HighsInt);
        for &k in &picked {
            let j = entries[k as usize].col() as HighsInt;
            covered[j as usize] = 1;
            cover.members.push(j);
            cover.member_pos.push(entries[k as usize].val() as u8);
        }
        num_groups += 1;
    }
    num_groups
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::ColCliqueIndex (lines 98-102)
struct ColCliqueIndex {
    start: Vec<HighsInt>,
    clq: Vec<HighsInt>,
    size: Vec<HighsInt>,
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::build_incidence (lines 104-146)
fn build_incidence(
    cliques: &[HighsCliqueTable_Clique],
    entries: &[HighsCliqueTable_CliqueVar],
    is_bin: &[u8],
    covered: &[u8],
    stamp: &mut [HighsInt],
    ncol: HighsInt,
) -> ColCliqueIndex {
    let num_cliques = cliques.len() as HighsInt;
    let mut idx = ColCliqueIndex {
        start: vec![0 as HighsInt; ncol as usize + 1],
        clq: Vec::new(),
        size: vec![0 as HighsInt; num_cliques as usize],
    };

    // First scan: count
    stamp.fill(-1);
    for c in 0..num_cliques {
        let cl = &cliques[c as usize];
        if !is_live(cl) {
            continue;
        }
        let mut k = cl.start;
        while k < cl.end {
            let j = entries[k as usize].col() as HighsInt;
            if is_bin[j as usize] == 0 || covered[j as usize] != 0 || stamp[j as usize] == c {
                k += 1;
                continue;
            }
            stamp[j as usize] = c;
            idx.start[j as usize + 1] += 1;
            idx.size[c as usize] += 1;
            k += 1;
        }
    }
    for j in 0..ncol {
        idx.start[j as usize + 1] += idx.start[j as usize];
    }
    idx.clq = vec![0 as HighsInt; idx.start[ncol as usize] as usize];

    // Second scan: fill
    let mut fill: Vec<HighsInt> = idx.start[..ncol as usize].to_vec();
    stamp.fill(-1);
    for c in 0..num_cliques {
        let cl = &cliques[c as usize];
        if !is_live(cl) {
            continue;
        }
        let mut k = cl.start;
        while k < cl.end {
            let j = entries[k as usize].col() as HighsInt;
            if is_bin[j as usize] == 0 || covered[j as usize] != 0 || stamp[j as usize] == c {
                k += 1;
                continue;
            }
            stamp[j as usize] = c;
            idx.clq[fill[j as usize] as usize] = c;
            fill[j as usize] += 1;
            k += 1;
        }
    }
    idx
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::assign_to_cliques (lines 153-198)
fn assign_to_cliques(
    bin: &[HighsInt],
    covered: &[u8],
    idx: &ColCliqueIndex,
    assign: &mut [HighsInt],
    assigned: &mut [HighsInt],
) {
    for &j in bin {
        if covered[j as usize] != 0 {
            continue;
        }
        let mut best: HighsInt = -1;
        let mut best_size: HighsInt = 0;
        let mut p = idx.start[j as usize];
        while p < idx.start[j as usize + 1] {
            let c = idx.clq[p as usize];
            if idx.size[c as usize] > best_size {
                best_size = idx.size[c as usize];
                best = c;
            }
            p += 1;
        }
        assign[j as usize] = best;
        if best >= 0 {
            assigned[best as usize] += 1;
        }
    }

    for &j in bin {
        let cur = assign[j as usize];
        if covered[j as usize] != 0 || cur < 0 {
            continue;
        }
        let mut best = cur;
        let mut best_size = assigned[cur as usize];
        let mut p = idx.start[j as usize];
        while p < idx.start[j as usize + 1] {
            let c = idx.clq[p as usize];
            if assigned[c as usize] > best_size {
                best_size = assigned[c as usize];
                best = c;
            }
            p += 1;
        }
        if best != cur {
            assigned[cur as usize] -= 1;
            assigned[best as usize] += 1;
            assign[j as usize] = best;
        }
    }
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::emit_selected_groups (lines 204-231)
fn emit_selected_groups(
    cliques: &[HighsCliqueTable_Clique],
    entries: &[HighsCliqueTable_CliqueVar],
    assign: &[HighsInt],
    assigned: &[HighsInt],
    covered: &mut [u8],
    cover: &mut Cover,
) {
    let mut selected: Vec<HighsInt> = Vec::new();
    let num_cliques = assigned.len() as HighsInt;
    for c in 0..num_cliques {
        if assigned[c as usize] > 0 {
            selected.push(c);
        }
    }
    selected.sort_by(|&a, &b| assigned[b as usize].cmp(&assigned[a as usize]));

    for &c in &selected {
        let cl = &cliques[c as usize];
        cover.group_start.push(cover.members.len() as HighsInt);
        let mut k = cl.start;
        while k < cl.end {
            let j = entries[k as usize].col() as HighsInt;
            if assign[j as usize] != c || covered[j as usize] != 0 {
                k += 1;
                continue;
            }
            covered[j as usize] = 1;
            cover.members.push(j);
            cover.member_pos.push(entries[k as usize].val() as u8);
            k += 1;
        }
    }
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::CliqueScan (lines 236-241)
struct CliqueScan {
    sum: f64,
    best_var: HighsInt,
    best_value: f64,
    skip: bool,
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::scan_clique (lines 243-269)
fn scan_clique(
    entries: &[HighsCliqueTable_CliqueVar],
    cl: &HighsCliqueTable_Clique,
    lp_ref: &[f64],
    dom_lower: &[f64],
    dom_upper: &[f64],
) -> CliqueScan {
    let mut out = CliqueScan {
        sum: 0.0,
        best_var: -1,
        best_value: 0.0,
        skip: false,
    };
    let mut k = cl.start;
    while k < cl.end {
        let lit = entries[k as usize];
        let j = lit.col() as HighsInt;
        let positive = lit.val() != 0;
        if positive {
            if dom_lower[j as usize] == 1.0 {
                out.skip = true;
                return out;
            }
        } else {
            if dom_upper[j as usize] == 0.0 {
                out.skip = true;
                return out;
            }
        }
        let v = if positive {
            lp_ref[j as usize]
        } else {
            1.0 - lp_ref[j as usize]
        };
        out.sum += v;
        let free_for_true = if positive {
            dom_upper[j as usize] == 1.0
        } else {
            dom_lower[j as usize] == 0.0
        };
        if v > out.best_value && free_for_true {
            out.best_var = j;
            out.best_value = v;
        }
        k += 1;
    }
    out
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::build_clique_cover (lines 273-302)
pub(crate) fn build_clique_cover(
    cliques: &[HighsCliqueTable_Clique],
    entries: &[HighsCliqueTable_CliqueVar],
    bin: &[HighsInt],
    ncol: HighsInt,
) -> Cover {
    let mut cover = Cover::default();
    if bin.is_empty() {
        cover.group_start.push(0);
        return cover;
    }

    let is_bin = binary_mask(bin, ncol);
    let mut covered = vec![0u8; ncol as usize];
    let mut stamp = vec![-1 as HighsInt; ncol as usize];

    cover.num_equality_groups = add_equality_groups(
        cliques,
        entries,
        &is_bin,
        &mut covered,
        &mut stamp,
        &mut cover,
    );

    let idx = build_incidence(cliques, entries, &is_bin, &covered, &mut stamp, ncol);
    let mut assign = vec![-1 as HighsInt; ncol as usize];
    let mut assigned = vec![0 as HighsInt; cliques.len()];
    assign_to_cliques(bin, &covered, &idx, &mut assign, &mut assigned);
    emit_selected_groups(
        cliques,
        entries,
        &assign,
        &assigned,
        &mut covered,
        &mut cover,
    );

    cover.group_start.push(cover.members.len() as HighsInt);

    for &j in bin {
        if covered[j as usize] == 0 {
            cover.uncovered.push(j);
        }
    }
    cover
}

/// mip-heuristics (36f7d9d): clique_cover.cpp::cliques2_order (lines 304-366)
pub(crate) fn cliques2_order(
    cliques: &[HighsCliqueTable_Clique],
    entries: &[HighsCliqueTable_CliqueVar],
    bin: &[HighsInt],
    ncol: HighsInt,
    lp_ref: &[f64],
    dom_lower: &[f64],
    dom_upper: &[f64],
) -> Vec<HighsInt> {
    let mut order: Vec<HighsInt> = Vec::new();
    if bin.is_empty() {
        return order;
    }
    order.reserve(bin.len());

    let is_bin = binary_mask(bin, ncol);
    let mut covered = vec![0u8; ncol as usize];

    let num_cliques = cliques.len() as HighsInt;
    for c in 0..num_cliques {
        let cl = &cliques[c as usize];
        if !is_live(cl) {
            continue;
        }

        let scan = scan_clique(entries, cl, lp_ref, dom_lower, dom_upper);

        if scan.skip || scan.best_var < 0 || (scan.sum - 1.0).abs() > kCliqueTightnessTol {
            continue;
        }

        // append best_var
        {
            let j = scan.best_var;
            if is_bin[j as usize] != 0 && covered[j as usize] == 0 {
                covered[j as usize] = 1;
                order.push(j);
            }
        }
        let mut k = cl.start;
        while k < cl.end {
            let j = entries[k as usize].col() as HighsInt;
            if j != scan.best_var {
                if is_bin[j as usize] != 0 && covered[j as usize] == 0 {
                    covered[j as usize] = 1;
                    order.push(j);
                }
            }
            k += 1;
        }
    }

    for &j in bin {
        if covered[j as usize] == 0 {
            order.push(j);
        }
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::HighsUInt;

    struct FakeTable {
        cliques: Vec<HighsCliqueTable_Clique>,
        entries: Vec<HighsCliqueTable_CliqueVar>,
    }

    impl FakeTable {
        fn new() -> Self {
            Self {
                cliques: Vec::new(),
                entries: Vec::new(),
            }
        }

        fn add(&mut self, lits: &[HighsCliqueTable_CliqueVar], equality: bool) {
            let mut cl = HighsCliqueTable_Clique::default();
            cl.start = self.entries.len() as HighsInt;
            cl.end = cl.start + lits.len() as HighsInt;
            cl.origin = kHighsIInf;
            cl.num_zero_fixed = 0;
            cl.equality = equality;
            self.cliques.push(cl);
            self.entries.extend_from_slice(lits);
        }

        fn kill_last(&mut self) {
            self.cliques.last_mut().unwrap().start = -1;
        }
    }

    fn pos(cols: &[HighsInt]) -> Vec<HighsCliqueTable_CliqueVar> {
        cols.iter()
            .map(|&c| HighsCliqueTable_CliqueVar::from_parts(c as HighsUInt, 1))
            .collect()
    }

    fn iota_cols(ncol: HighsInt) -> Vec<HighsInt> {
        (0..ncol).collect()
    }

    fn flatten(cover: &Cover) -> Vec<HighsInt> {
        let mut out = cover.members.clone();
        out.extend_from_slice(&cover.uncovered);
        out
    }

    fn group(cover: &Cover, g: HighsInt) -> Vec<HighsInt> {
        cover.members
            [cover.group_start[g as usize] as usize..cover.group_start[g as usize + 1] as usize]
            .to_vec()
    }

    fn group_signs(cover: &Cover, g: HighsInt) -> Vec<u8> {
        cover.member_pos
            [cover.group_start[g as usize] as usize..cover.group_start[g as usize + 1] as usize]
            .to_vec()
    }

    #[test]
    fn s5_clique_equality_cliques_come_before_inequality_cliques() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2, 3]), false);
        t.add(&pos(&[4, 5]), true);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(6), 6);

        assert_eq!(cover.num_equality_groups, 1);
        assert_eq!(cover.num_groups(), 2);
        assert_eq!(group(&cover, 0), vec![4, 5]);
        assert_eq!(group(&cover, 1), vec![0, 1, 2, 3]);
    }

    #[test]
    fn s5_clique_an_equality_clique_overlapping_an_accepted_one_is_dropped() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1]), true);
        t.add(&pos(&[1, 5]), true);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(6), 6);

        assert_eq!(cover.num_equality_groups, 1);
        assert_eq!(group(&cover, 0), vec![0, 1]);

        assert_eq!(cover.num_groups(), 2);
        assert_eq!(group(&cover, 1), vec![5]);
        assert_eq!(cover.uncovered, vec![2, 3, 4]);
    }

    #[test]
    fn s5_clique_selected_cliques_are_sorted_by_size_descending() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1]), false);
        t.add(&pos(&[2, 3, 4, 5]), false);
        t.add(&pos(&[6, 7, 8]), false);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(9), 9);

        assert_eq!(cover.num_equality_groups, 0);
        assert_eq!(cover.num_groups(), 3);
        assert_eq!(group(&cover, 0), vec![2, 3, 4, 5]);
        assert_eq!(group(&cover, 1), vec![6, 7, 8]);
        assert_eq!(group(&cover, 2), vec![0, 1]);
    }

    #[test]
    fn s5_clique_within_a_clique_the_tables_entry_order_is_kept() {
        let mut t = FakeTable::new();
        t.add(&pos(&[5, 2, 9]), false);
        t.add(&pos(&[7, 1]), true);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(10), 10);

        assert_eq!(cover.num_groups(), 2);
        assert_eq!(group(&cover, 0), vec![7, 1]);
        assert_eq!(group(&cover, 1), vec![5, 2, 9]);
    }

    #[test]
    fn s5_clique_uncovered_binaries_are_appended_in_formulation_order() {
        let mut t = FakeTable::new();
        t.add(&pos(&[3, 1]), false);
        t.add(&pos(&[6, 5]), false);
        t.kill_last();
        t.add(&pos(&[2, 7]), false);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(8), 8);

        assert_eq!(cover.uncovered, vec![0, 4, 5, 6]);
        assert_eq!(flatten(&cover).len(), 8);
    }

    #[test]
    fn s5_clique_every_binary_appears_exactly_once() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2]), true);
        t.add(&pos(&[1, 3, 4]), false);
        t.add(&pos(&[2, 3, 8]), false);
        t.add(&pos(&[5, 5, 6]), false);
        t.add(&pos(&[4, 7]), false);
        t.kill_last();

        let bin: Vec<HighsInt> = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let cover = build_clique_cover(&t.cliques, &t.entries, &bin, 9);

        let mut flat = flatten(&cover);
        flat.sort();
        assert_eq!(flat, bin);
    }

    #[test]
    fn s5_clique_a_shared_binary_goes_to_the_largest_covering_clique() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2, 3]), false);
        t.add(&pos(&[3, 4, 5]), false);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(6), 6);

        assert_eq!(cover.num_groups(), 2);
        assert_eq!(group(&cover, 0), vec![0, 1, 2, 3]);
        assert_eq!(group(&cover, 1), vec![4, 5]);
    }

    #[test]
    fn s5_clique_the_local_adjustment_moves_a_stranded_variable() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2, 3, 4, 5, 6]), false);
        t.add(&pos(&[0, 1, 2, 3, 4, 9]), false);
        t.add(&pos(&[7, 8, 9]), false);

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(10), 10);

        assert_eq!(cover.num_groups(), 2);
        assert_eq!(group(&cover, 0), vec![0, 1, 2, 3, 4, 5, 6]);
        assert_eq!(group(&cover, 1), vec![7, 8, 9]);
        assert!(cover.uncovered.is_empty());
    }

    #[test]
    fn s5_clique_a_group_records_its_literals_signs() {
        let mut t = FakeTable::new();
        t.add(
            &[
                HighsCliqueTable_CliqueVar::from_parts(2, 0),
                HighsCliqueTable_CliqueVar::from_parts(0, 1),
                HighsCliqueTable_CliqueVar::from_parts(1, 0),
            ],
            false,
        );

        let cover = build_clique_cover(&t.cliques, &t.entries, &iota_cols(3), 3);

        assert_eq!(cover.num_groups(), 1);
        assert_eq!(group(&cover, 0), vec![2, 0, 1]);
        assert_eq!(group_signs(&cover, 0), vec![0, 1, 0]);
    }

    #[test]
    fn s5_clique_cliques2_a_clique_with_a_literal_already_fixed_true_is_skipped_entirely() {
        let mut t = FakeTable::new();
        t.add(&pos(&[2, 1, 0]), false);

        let lp: Vec<f64> = vec![0.5, 1.0, 0.5];
        let lb: Vec<f64> = vec![0.0, 1.0, 0.0];
        let ub: Vec<f64> = vec![1.0, 1.0, 1.0];
        let bin: Vec<HighsInt> = vec![0, 2];

        let order = cliques2_order(&t.cliques, &t.entries, &bin, 3, &lp, &lb, &ub);
        assert_eq!(order, vec![0, 2]);
    }

    #[test]
    fn s5_clique_cliques2_the_most_positive_literal_leads_a_tight_clique() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2]), false);

        let lp: Vec<f64> = vec![0.2, 0.1, 0.7];
        let lb: Vec<f64> = vec![0.0, 0.0, 0.0];
        let ub: Vec<f64> = vec![1.0, 1.0, 1.0];

        let order = cliques2_order(&t.cliques, &t.entries, &iota_cols(3), 3, &lp, &lb, &ub);
        assert_eq!(order, vec![2, 0, 1]);
    }

    #[test]
    fn s5_clique_cliques2_a_clique_that_is_not_lp_tight_contributes_nothing() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2]), false);

        let lp: Vec<f64> = vec![0.2, 0.1, 0.2];
        let lb: Vec<f64> = vec![0.0, 0.0, 0.0];
        let ub: Vec<f64> = vec![1.0, 1.0, 1.0];

        let order = cliques2_order(&t.cliques, &t.entries, &iota_cols(3), 3, &lp, &lb, &ub);
        assert_eq!(order, vec![0, 1, 2]);
    }

    #[test]
    fn s5_clique_cliques2_the_tightness_test_is_two_sided() {
        let mut t = FakeTable::new();
        t.add(&pos(&[0, 1, 2]), false);

        let lp: Vec<f64> = vec![0.9, 0.8, 0.7];
        let lb: Vec<f64> = vec![0.0, 0.0, 0.0];
        let ub: Vec<f64> = vec![1.0, 1.0, 1.0];

        let order = cliques2_order(&t.cliques, &t.entries, &iota_cols(3), 3, &lp, &lb, &ub);
        assert_eq!(order, vec![0, 1, 2]);
    }

    #[test]
    fn s5_clique_cliques2_negative_literals_weigh_1_minus_x_and_die_on_u_j_eq_0() {
        let mut t = FakeTable::new();
        t.add(
            &[
                HighsCliqueTable_CliqueVar::from_parts(0, 0),
                HighsCliqueTable_CliqueVar::from_parts(1, 1),
            ],
            false,
        );
        t.add(
            &[
                HighsCliqueTable_CliqueVar::from_parts(2, 0),
                HighsCliqueTable_CliqueVar::from_parts(3, 1),
            ],
            false,
        );

        let lp: Vec<f64> = vec![0.1, 0.1, 0.0, 0.9];
        let lb: Vec<f64> = vec![0.0, 0.0, 0.0, 0.0];
        let ub: Vec<f64> = vec![1.0, 1.0, 0.0, 1.0];
        let bin: Vec<HighsInt> = vec![0, 1, 3];

        let order = cliques2_order(&t.cliques, &t.entries, &bin, 4, &lp, &lb, &ub);
        assert_eq!(order, vec![0, 1, 3]);
    }
}
