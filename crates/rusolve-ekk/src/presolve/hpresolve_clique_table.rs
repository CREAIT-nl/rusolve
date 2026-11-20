//! HiGHS source: presolve/HPresolveCliqueTable.h, presolve/HPresolveCliqueTable.cpp
//! (HiGHS latest, #3312 "clique table through presolve", merge 7b890c255)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Reached only from `mip_presolve_clique_table` gates: with the option off nothing in
//! this file runs and `HPresolve::presolve_clique_table` stays empty.
//!
//! # Design: no stored table pointer
//!
//! C++ keeps `HighsCliqueTable* table`, set by `rebuild`. The port stores no pointer: every
//! method takes the table (`mipsolver->mipdata_->cliquetable`, the only table C++ ever
//! passes) as `&mut HighsCliqueTable`. `HPresolveCliqueTable` is a friend of
//! `HighsCliqueTable` in C++; the Rust fields it touches are `pub`.

use crate::mip::clique_table::{HighsCliqueTable, HighsCliqueTable_CliqueVar};
use crate::prelude::*;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::cpp_max;

/// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::HPresolveCliqueTable::ColState (lines 18-34)
/// C++: `enum Value : int8_t { kActive, kFixedZero, kFixedOne, kEliminated }`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HPresolveCliqueTable_ColState {
    #[default]
    KActive,
    KFixedZero,
    KFixedOne,
    KEliminated,
}

impl HPresolveCliqueTable_ColState {
    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::ColState::isActive (line 25)
    pub fn is_active(self) -> bool {
        self == Self::KActive
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::ColState::isEliminated (line 26)
    pub fn is_eliminated(self) -> bool {
        self == Self::KEliminated
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::ColState::isFixed (line 27)
    pub fn is_fixed(self) -> bool {
        self == Self::KFixedOne || self == Self::KFixedZero
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::ColState::isFixedTo (lines 28-30)
    pub fn is_fixed_to(self, val: bool) -> bool {
        self == if val {
            Self::KFixedOne
        } else {
            Self::KFixedZero
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::ColState::fixed (lines 31-33)
    pub fn fixed(val: bool) -> Self {
        if val {
            Self::KFixedOne
        } else {
            Self::KFixedZero
        }
    }
}

/// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.h::HPresolveCliqueTable (lines 15-56)
#[derive(Clone, Debug, Default)]
pub struct HPresolveCliqueTable {
    /// C++: `std::vector<ColState> colStates`
    pub col_states: Vec<HPresolveCliqueTable_ColState>,
    /// C++: `std::vector<HighsCliqueTable::CliqueVar> fixingQueue`
    pub fixing_queue: Vec<HighsCliqueTable_CliqueVar>,
    /// C++: `std::vector<HighsInt> incidentCliques`
    pub incident_cliques: Vec<HighsInt>,
    /// C++: `std::vector<HighsCliqueTable::CliqueVar> shortenedClique`
    pub shortened_clique: Vec<HighsCliqueTable_CliqueVar>,
}

/// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::substituteCol, struct `Overlap` (lines 202-206)
#[derive(Clone, Copy, Debug)]
struct Overlap {
    clique_id: HighsInt,
    subst_pos: HighsInt,
    replace_pos: HighsInt,
}

/// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::fixCol, lambda `collectIncidentCliques` (lines 65-73)
fn collect_incident_cliques(
    incident_cliques: &mut Vec<HighsInt>,
    table: &HighsCliqueTable,
    v: HighsCliqueTable_CliqueVar,
) {
    incident_cliques.clear();
    table.inverted_hash_list[v.index() as usize].for_each(&mut |clique_id, _| {
        incident_cliques.push(*clique_id);
    });
    table.inverted_hash_list_size_two[v.index() as usize].for_each(&mut |clique_id, _| {
        incident_cliques.push(*clique_id);
    });
}

impl HPresolveCliqueTable {
    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::rebuild (lines 18-21)
    /// C++: `void rebuild(HighsCliqueTable& table)`
    /// PORT NOTE: `setCliqueTable(&table)` has no counterpart (no stored pointer, module doc).
    pub fn rebuild(&mut self, table: &HighsCliqueTable) {
        self.col_states.clear();
        self.col_states.resize(
            table.col_deleted.len(),
            HPresolveCliqueTable_ColState::default(),
        );
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::checkCompactClique (lines 23-43)
    /// C++: `void checkCompactClique(HighsInt cliqueId, HighsInt threshold, HighsInt activeSize, HighsInt actualSize, bool equality, HighsInt origin)`
    fn check_compact_clique(
        &mut self,
        table: &mut HighsCliqueTable,
        clique_id: HighsInt,
        threshold: HighsInt,
        active_size: HighsInt,
        actual_size: HighsInt,
        equality: bool,
        origin: HighsInt,
    ) {
        let clique = table.cliques[clique_id as usize];
        if active_size == 2 || clique.num_zero_fixed >= cpp_max(threshold, actual_size >> 1) {
            self.shortened_clique.clear();
            self.shortened_clique.reserve(active_size as usize);
            for i in clique.start..clique.end {
                let entry = table.cliqueentries[i as usize];
                if table.col_deleted[entry.col() as usize] == 0 {
                    self.shortened_clique.push(entry);
                }
            }
            table.remove_clique(clique_id, false);
            table.do_add_clique(&self.shortened_clique, equality, origin);
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::fixCol (lines 45-134)
    /// C++: `bool fixCol(HighsInt col, bool val, std::vector<HighsCliqueTable::CliqueVar>& impliedFixings)`
    pub fn fix_col(
        &mut self,
        table: &mut HighsCliqueTable,
        col: HighsInt,
        val: bool,
        implied_fixings: &mut Vec<HighsCliqueTable_CliqueVar>,
    ) -> bool {
        let state = self.col_states[col as usize];
        if !state.is_active() {
            if state.is_eliminated() {
                return true;
            }
            return state.is_fixed_to(val);
        }

        let cu = col as usize;
        if table.inverted_hash_list[2 * cu].empty()
            && table.inverted_hash_list[2 * cu + 1].empty()
            && table.inverted_hash_list_size_two[2 * cu].empty()
            && table.inverted_hash_list_size_two[2 * cu + 1].empty()
        {
            self.col_states[cu] = HPresolveCliqueTable_ColState::fixed(val);
            table.col_deleted[cu] = 1;
            return true;
        }

        self.fixing_queue.clear();
        self.fixing_queue
            .push(HighsCliqueTable_CliqueVar::new(col, val as HighsInt));
        let mut next_fixing: usize = 0;

        while next_fixing != self.fixing_queue.len() {
            let v = self.fixing_queue[next_fixing];
            next_fixing += 1;
            let state = self.col_states[v.col() as usize];
            if state.is_eliminated() || state.is_fixed_to(v.val() != 0) {
                continue;
            }
            if state.is_fixed() {
                return false;
            }
            self.col_states[v.col() as usize] = HPresolveCliqueTable_ColState::fixed(v.val() != 0);
            table.col_deleted[v.col() as usize] = 1;
            if !(v.col() as HighsInt == col && v.val() == val as HighsUInt) {
                implied_fixings.push(v);
            }

            // Fix all other literals in incident cliques to be inactive
            collect_incident_cliques(&mut self.incident_cliques, table, v);
            for k in 0..self.incident_cliques.len() {
                let clique_id = self.incident_cliques[k];
                if table.cliques[clique_id as usize].start == -1 {
                    continue;
                }
                let start = table.cliques[clique_id as usize].start;
                let end = table.cliques[clique_id as usize].end;
                for i in start..end {
                    let v2 = table.cliqueentries[i as usize];
                    if v.col() == v2.col() {
                        continue;
                    }
                    let other_state = self.col_states[v2.col() as usize];
                    if other_state.is_fixed_to(v2.val() != 0) {
                        return false;
                    }
                    if other_state.is_active() {
                        self.fixing_queue.push(HighsCliqueTable_CliqueVar::new(
                            v2.col() as HighsInt,
                            1 - v2.val() as HighsInt,
                        ));
                    }
                }
                table.remove_clique(clique_id, true);
            }

            // Remove complement-literal from cliques
            collect_incident_cliques(&mut self.incident_cliques, table, v.complement());
            table.inverted_hash_list[v.complement().index() as usize].clear();
            table.inverted_hash_list_size_two[v.complement().index() as usize].clear();
            for k in 0..self.incident_cliques.len() {
                let clique_id = self.incident_cliques[k];
                let cid = clique_id as usize;
                if table.cliques[cid].start == -1 {
                    continue;
                }
                let equality = table.cliques[cid].equality;
                let origin = table.cliques[cid].origin;
                table.cliques[cid].num_zero_fixed += 1;
                let actual_size = table.cliques[cid].end - table.cliques[cid].start;
                let active_size = table.cliques[cid].num_active();
                if active_size <= 1 {
                    if equality {
                        if active_size == 0 {
                            return false;
                        }
                        for i in table.cliques[cid].start..table.cliques[cid].end {
                            let entry = table.cliqueentries[i as usize];
                            if table.col_deleted[entry.col() as usize] == 0 {
                                self.fixing_queue.push(entry);
                                break;
                            }
                        }
                    }
                    table.remove_clique(clique_id, true);
                    continue;
                }

                if !table.in_presolve {
                    self.check_compact_clique(
                        table,
                        clique_id,
                        10,
                        active_size,
                        actual_size,
                        equality,
                        origin,
                    );
                }
            }
        }
        true
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::eliminateCol (lines 136-177)
    /// C++: `void eliminateCol(HighsInt col)`
    pub fn eliminate_col(&mut self, table: &mut HighsCliqueTable, col: HighsInt) {
        let cu = col as usize;
        if self.col_states[cu].is_eliminated() {
            return;
        }
        self.col_states[cu] = HPresolveCliqueTable_ColState::KEliminated;
        if table.col_deleted[cu] != 0 {
            return;
        }
        table.col_deleted[cu] = 1;

        let incident_cliques = &mut self.incident_cliques;
        incident_cliques.clear();
        table.inverted_hash_list[2 * cu].for_each(&mut |clique_id, _| {
            incident_cliques.push(*clique_id);
        });
        table.inverted_hash_list[2 * cu].clear();
        table.inverted_hash_list_size_two[2 * cu].for_each(&mut |clique_id, _| {
            incident_cliques.push(*clique_id);
        });
        table.inverted_hash_list_size_two[2 * cu].clear();
        table.inverted_hash_list[2 * cu + 1].for_each(&mut |clique_id, _| {
            incident_cliques.push(*clique_id);
        });
        table.inverted_hash_list[2 * cu + 1].clear();
        table.inverted_hash_list_size_two[2 * cu + 1].for_each(&mut |clique_id, _| {
            incident_cliques.push(*clique_id);
        });
        table.inverted_hash_list_size_two[2 * cu + 1].clear();

        // C++ `pdqsort(begin, end)` (std::less); the keys are integers, so the result does not
        // depend on the pdqsort variant
        pdqsort_by(&mut self.incident_cliques, |a: &HighsInt, b: &HighsInt| {
            a < b
        });

        for k in 0..self.incident_cliques.len() {
            let clique_id = self.incident_cliques[k];
            let cid = clique_id as usize;
            if table.cliques[cid].start == -1 {
                continue;
            }
            table.cliques[cid].num_zero_fixed += 1;
            table.cliques[cid].origin = -1;
            table.cliques[cid].equality = false;
            let actual_size = table.cliques[cid].end - table.cliques[cid].start;
            let active_size = table.cliques[cid].num_active();
            if active_size <= 1 {
                table.remove_clique(clique_id, false);
                continue;
            }
            if !table.in_presolve {
                self.check_compact_clique(
                    table,
                    clique_id,
                    10,
                    active_size,
                    actual_size,
                    false,
                    -1,
                );
            }
        }
    }

    /// HiGHS latest (#3312, 7b890c255): HPresolveCliqueTable.cpp::HPresolveCliqueTable::substituteCol (lines 179-308)
    /// C++: `bool substituteCol(HighsInt substCol, HighsCliqueTable::CliqueVar replacement, std::vector<HighsCliqueTable::CliqueVar>& impliedFixings)`
    pub fn substitute_col(
        &mut self,
        table: &mut HighsCliqueTable,
        subst_col: HighsInt,
        replacement: HighsCliqueTable_CliqueVar,
        implied_fixings: &mut Vec<HighsCliqueTable_CliqueVar>,
    ) -> bool {
        if subst_col == replacement.col() as HighsInt {
            // Substituting itself is either redundant or infeasible
            return replacement.val() == 1;
        }
        if !self.col_states[subst_col as usize].is_active() {
            // Don't do anything if col to be substituted is already handled
            return true;
        }
        if self.col_states[replacement.col() as usize].is_eliminated() {
            // Eliminate substituted column if replacement is eliminated
            self.eliminate_col(table, subst_col);
            return true;
        }
        if self.col_states[replacement.col() as usize].is_fixed() {
            // If replacement is fixed then bypass direct substitutions via fixings
            let fix_val =
                self.col_states[replacement.col() as usize].is_fixed_to(replacement.val() != 0);
            implied_fixings.push(HighsCliqueTable_CliqueVar::new(
                subst_col,
                fix_val as HighsInt,
            ));
            if !self.fix_col(table, subst_col, fix_val, implied_fixings) {
                return false;
            }
            return true;
        }
        let mut overlaps: Vec<Overlap> = Vec::new();
        // lambda `collectOverlaps` (lines 208-229)
        let mut collect_overlaps = |table: &HighsCliqueTable, v: HighsCliqueTable_CliqueVar| {
            table.inverted_hash_list[v.index() as usize].for_each(&mut |clique_id, subst_pos| {
                let mut replace_pos =
                    table.inverted_hash_list[replacement.index() as usize].find(clique_id);
                if replace_pos.is_none() {
                    replace_pos = table.inverted_hash_list
                        [replacement.complement().index() as usize]
                        .find(clique_id);
                }
                if let Some(replace_pos) = replace_pos {
                    overlaps.push(Overlap {
                        clique_id: *clique_id,
                        subst_pos: *subst_pos,
                        replace_pos: *replace_pos,
                    });
                }
            });
            table.inverted_hash_list_size_two[v.index() as usize].for_each(&mut |clique_id, _| {
                let mut subst_pos = table.cliques[*clique_id as usize].start;
                let mut replace_pos = subst_pos + 1;
                if table.cliqueentries[replace_pos as usize] == v {
                    std::mem::swap(&mut subst_pos, &mut replace_pos);
                }
                if table.cliqueentries[replace_pos as usize].col() == replacement.col() {
                    overlaps.push(Overlap {
                        clique_id: *clique_id,
                        subst_pos,
                        replace_pos,
                    });
                }
            });
        };
        collect_overlaps(table, HighsCliqueTable_CliqueVar::new(subst_col, 0));
        collect_overlaps(table, HighsCliqueTable_CliqueVar::new(subst_col, 1));
        pdqsort_by(&mut overlaps, |a: &Overlap, b: &Overlap| {
            a.clique_id < b.clique_id
        });

        let mut potential_fixings: Vec<HighsCliqueTable_CliqueVar> = Vec::new();

        for overlap in &overlaps {
            let clique_id = overlap.clique_id;
            let cid = clique_id as usize;
            if table.cliques[cid].start == -1 {
                continue;
            }
            let subst_pos = overlap.subst_pos;
            let replace_pos = overlap.replace_pos;
            table.cliques[cid].origin = -1;
            let equality = table.cliques[cid].equality;

            let subst_var = if table.cliqueentries[subst_pos as usize].val() != 0 {
                replacement
            } else {
                replacement.complement()
            };
            let replacement_var = table.cliqueentries[replace_pos as usize];

            // If both replacement and its complement will exist in new clique
            // then clique is automatically fulfilled. Set all other literals to 0
            if subst_var.val() != replacement_var.val() {
                for i in table.cliques[cid].start..table.cliques[cid].end {
                    if i != subst_pos
                        && i != replace_pos
                        && table.col_deleted[table.cliqueentries[i as usize].col() as usize] == 0
                    {
                        potential_fixings.push(table.cliqueentries[i as usize].complement());
                    }
                }
                table.remove_clique(clique_id, false);
                continue;
            }

            // If replacement will exist twice in new clique then the literal
            // has to take value 0
            potential_fixings.push(subst_var.complement());
            self.shortened_clique.clear();
            self.shortened_clique
                .reserve((table.cliques[cid].end - table.cliques[cid].start) as usize);
            for i in table.cliques[cid].start..table.cliques[cid].end {
                if i != subst_pos
                    && i != replace_pos
                    && table.col_deleted[table.cliqueentries[i as usize].col() as usize] == 0
                {
                    self.shortened_clique.push(table.cliqueentries[i as usize]);
                }
            }
            table.remove_clique(clique_id, false);

            if self.shortened_clique.len() >= 2 {
                table.do_add_clique(&self.shortened_clique, equality, -1);
            } else if equality {
                if self.shortened_clique.is_empty() {
                    return false;
                }
                potential_fixings.push(self.shortened_clique[0]);
            }
        }

        // Now substitute all entries
        table.replace_literal(HighsCliqueTable_CliqueVar::new(subst_col, 1), replacement);
        table.replace_literal(
            HighsCliqueTable_CliqueVar::new(subst_col, 0),
            replacement.complement(),
        );
        self.col_states[subst_col as usize] = HPresolveCliqueTable_ColState::KEliminated;
        table.col_deleted[subst_col as usize] = 1;

        for v in potential_fixings {
            let state = self.col_states[v.col() as usize];
            if state.is_eliminated() || state.is_fixed_to(v.val() != 0) {
                continue;
            }
            if state.is_fixed() {
                return false;
            }
            implied_fixings.push(v);
            if !self.fix_col(table, v.col() as HighsInt, v.val() != 0, implied_fixings) {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cv(col: HighsInt, val: HighsInt) -> HighsCliqueTable_CliqueVar {
        HighsCliqueTable_CliqueVar::new(col, val)
    }

    fn table_with(
        ncols: HighsInt,
        cliques: &[(&[HighsCliqueTable_CliqueVar], bool)],
    ) -> HighsCliqueTable {
        let mut table = HighsCliqueTable::new(ncols);
        for (clq, equality) in cliques {
            table.do_add_clique(clq, *equality, kHighsIInf);
        }
        table
    }

    fn active_cliques(table: &HighsCliqueTable) -> Vec<Vec<(HighsUInt, HighsUInt)>> {
        let mut out = Vec::new();
        for clq in &table.cliques {
            if clq.start == -1 {
                continue;
            }
            let mut entries: Vec<(HighsUInt, HighsUInt)> = (clq.start..clq.end)
                .map(|i| table.cliqueentries[i as usize])
                .filter(|e| table.col_deleted[e.col() as usize] == 0)
                .map(|e| (e.col(), e.val()))
                .collect();
            entries.sort();
            out.push(entries);
        }
        out.sort();
        out
    }

    #[test]
    fn x1b_clqp_a1_fix_col_queue() {
        // clique {x0, x1, x2}: fixCol(x0, 1) implies x1 = 0 and x2 = 0
        let mut table = table_with(3, &[(&[cv(0, 1), cv(1, 1), cv(2, 1)], false)]);
        let mut pct = HPresolveCliqueTable::default();
        pct.rebuild(&table);
        let mut implied = Vec::new();
        assert!(pct.fix_col(&mut table, 0, true, &mut implied));
        let got: Vec<(HighsUInt, HighsUInt)> = implied.iter().map(|v| (v.col(), v.val())).collect();
        assert_eq!(got, vec![(1, 0), (2, 0)]);
        assert_eq!(pct.col_states[0], HPresolveCliqueTable_ColState::KFixedOne);
        assert_eq!(pct.col_states[1], HPresolveCliqueTable_ColState::KFixedZero);
        assert_eq!(pct.col_states[2], HPresolveCliqueTable_ColState::KFixedZero);
        assert_eq!(table.col_deleted, vec![1, 1, 1]);
        assert!(active_cliques(&table).is_empty());
        // a conflicting fixing is infeasible, a repeated one is a no-op
        assert!(!pct.fix_col(&mut table, 1, true, &mut implied));
        assert!(pct.fix_col(&mut table, 1, false, &mut implied));
    }

    #[test]
    fn x1b_clqp_a1_fix_col_equality_last_literal() {
        // equality clique x0 + x1 = 1 (size two): fixCol(x0, 0) implies x1 = 1
        let mut table = table_with(2, &[(&[cv(0, 1), cv(1, 1)], true)]);
        let mut pct = HPresolveCliqueTable::default();
        pct.rebuild(&table);
        let mut implied = Vec::new();
        assert!(pct.fix_col(&mut table, 0, false, &mut implied));
        let got: Vec<(HighsUInt, HighsUInt)> = implied.iter().map(|v| (v.col(), v.val())).collect();
        assert_eq!(got, vec![(1, 1)]);
    }

    #[test]
    fn x1b_clqp_a1_eliminate_col() {
        // x0 is in {x0, x1, x2} and {~x0, x3, x4}; eliminating it removes it from both
        let mut table = table_with(
            5,
            &[
                (&[cv(0, 1), cv(1, 1), cv(2, 1)], true),
                (&[cv(0, 0), cv(3, 1), cv(4, 1)], false),
            ],
        );
        let mut pct = HPresolveCliqueTable::default();
        pct.rebuild(&table);
        pct.eliminate_col(&mut table, 0);
        assert_eq!(
            pct.col_states[0],
            HPresolveCliqueTable_ColState::KEliminated
        );
        assert_eq!(table.col_deleted[0], 1);
        assert!(table.inverted_hash_list[0].empty() && table.inverted_hash_list[1].empty());
        assert_eq!(
            active_cliques(&table),
            vec![vec![(1, 1), (2, 1)], vec![(3, 1), (4, 1)]]
        );
        // the equality status is dropped (eliminateCol sets equality = false)
        for clq in &table.cliques {
            if clq.start != -1 {
                assert!(!clq.equality);
            }
        }
        // a second call is a no-op
        pct.eliminate_col(&mut table, 0);
        assert_eq!(active_cliques(&table).len(), 2);
    }

    #[test]
    fn x1b_clqp_a2_substitute_col_merges() {
        // x0 := 1 - x1 in clique {x0, x2, x3}: x0 = 1 becomes ~x1, so the clique is {~x1, x2, x3}
        let mut table = table_with(4, &[(&[cv(0, 1), cv(2, 1), cv(3, 1)], false)]);
        let mut pct = HPresolveCliqueTable::default();
        pct.rebuild(&table);
        let mut implied = Vec::new();
        assert!(pct.substitute_col(&mut table, 0, cv(1, 0), &mut implied));
        assert!(implied.is_empty());
        assert_eq!(
            pct.col_states[0],
            HPresolveCliqueTable_ColState::KEliminated
        );
        assert_eq!(active_cliques(&table), vec![vec![(1, 0), (2, 1), (3, 1)]]);
        assert!(table.inverted_hash_list[0].empty() && table.inverted_hash_list[1].empty());
    }

    #[test]
    fn x1b_clqp_a2_substitute_col_overlap_fixes() {
        // x0 := x1 in clique {x0, x1, x2}: x1 appears twice, so x1 = 0 (and the shortened
        // clique {x2} is dropped); x0 := ~x1 in {x0, x1, x2} instead fixes x2 = 0
        let mut table = table_with(3, &[(&[cv(0, 1), cv(1, 1), cv(2, 1)], false)]);
        let mut pct = HPresolveCliqueTable::default();
        pct.rebuild(&table);
        let mut implied = Vec::new();
        assert!(pct.substitute_col(&mut table, 0, cv(1, 1), &mut implied));
        let got: Vec<(HighsUInt, HighsUInt)> = implied.iter().map(|v| (v.col(), v.val())).collect();
        assert_eq!(got, vec![(1, 0)]);

        let mut table = table_with(3, &[(&[cv(0, 1), cv(1, 1), cv(2, 1)], false)]);
        let mut pct = HPresolveCliqueTable::default();
        pct.rebuild(&table);
        let mut implied = Vec::new();
        assert!(pct.substitute_col(&mut table, 0, cv(1, 0), &mut implied));
        let got: Vec<(HighsUInt, HighsUInt)> = implied.iter().map(|v| (v.col(), v.val())).collect();
        assert_eq!(got, vec![(2, 0)]);
    }
}
