//! HiGHS source: mip/HighsMachineSchedSeparator.h, mip/HighsMachineSchedSeparator.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Port of HiGHS PR #2990 (merge dc83470a9).
//!
//! Off path: nothing in this file is reachable until HiGHS backport-MSS-B registers the separator,
//! so no option is read here.

use crate::prelude::*;
use crate::util::highs_hash::{fibonacci_muliplier, hash_u32_pair, pair_hash};
use crate::util::highs_hash_table::{HighsHashKey, HighsHashTable, HighsHashTableEntry};
use crate::util::highs_pdqsort::pdqsort_by;

/// HiGHS latest (#2990, dc83470a9): HighsMachineSchedSeparator.cpp (lines 23-26)
/// C++: `enum class ArcType { kImplicationWhenOne, kImplicationWhenZero }`
#[derive(Clone, Copy, Default, PartialEq, Debug)]
enum ArcType {
    #[default]
    kImplicationWhenOne,
    kImplicationWhenZero,
}

// HighsHashKey for (HighsInt, HighsInt) — C++ 8-byte pair key.
impl HighsHashKey for (HighsInt, HighsInt) {
    fn highs_hash(&self) -> u64 {
        hash_u32_pair(self.0 as u32, self.1 as u32)
    }
}

// HighsHashKey for (HighsInt, HighsInt, HighsInt) — C++ `std::tuple<HighsInt, HighsInt,
// HighsInt>`, a 12-byte key: HighsHash.h::HighsHashHelpers::hash, branch
// `9 <= sizeof(T) <= 16` (lines 578-590), a memcpy of the tuple. libstdc++ lays a
// `std::tuple` out last element first, so the bytes are (get<2>, get<1>, get<0>).
impl HighsHashKey for (HighsInt, HighsInt, HighsInt) {
    fn highs_hash(&self) -> u64 {
        // memcpy of the 12 bytes; bytes[3] stays 0 because sizeof(T) < 16
        let bytes: [u32; 4] = [self.2 as u32, self.1 as u32, self.0 as u32, 0];
        (pair_hash(0, bytes[0], bytes[1]) ^ (pair_hash(1, bytes[2], bytes[3]) >> 32))
            .wrapping_mul(fibonacci_muliplier())
    }
}

/// HiGHS latest (#2990, dc83470a9): HighsMachineSchedSeparator.h (lines 45-69)
/// C++: `class HighsMachineSchedSeparator : public HighsSeparator`
#[derive(Clone, Debug)]
pub struct HighsMachineSchedSeparator {
    pub base: HighsSeparator,
    pub randgen: HighsRandom,
    pub has_single_machine_schedule: bool,
    pub already_tried: bool,
}

impl HighsMachineSchedSeparator {
    /// HiGHS latest (#2990, dc83470a9): HighsMachineSchedSeparator.h (lines 63-68)
    /// C++: `HighsMachineSchedSeparator(const HighsMipSolver& mipsolver)`
    pub fn new(ms: &HighsMipSolver) -> Self {
        let mut randgen = HighsRandom::new(0);
        randgen.initialise(ms.options_mip_.random_seed as HighsUInt);
        HighsMachineSchedSeparator {
            base: HighsSeparator::new(ms, kMachineSchedSepaString),
            randgen,
            has_single_machine_schedule: false,
            already_tried: false,
        }
    }

    /// HiGHS latest (#2990, dc83470a9): HighsMachineSchedSeparator.cpp::HighsMachineSchedSeparator::separateLpSolution (lines 286-314)
    pub fn separate_lp_solution(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &mut HighsLpRelaxation,
        _lp_aggregator: &mut HighsLpAggregator,
        _trans_lp: &mut HighsTransformedLp,
        cutpool: usize,
    ) {
        // Only try to separate once
        if self.already_tried {
            return;
        }
        let mut vals: Vec<Vec<f64>> = Vec::new();
        let mut inds: Vec<Vec<HighsInt>> = Vec::new();
        let mut rhss: Vec<f64> = Vec::new();
        self.has_single_machine_schedule =
            self.find_single_machine_schedule_clique(ms, &mut vals, &mut inds, &mut rhss);
        if !self.has_single_machine_schedule {
            self.already_tried = true;
            return;
        }

        // Load the cuts
        let lp_solution = &lp_relaxation.get_solution().col_value;
        for i in 0..inds.len() {
            let mut viol = -rhss[i];
            for j in 0..inds[i].len() {
                viol += lp_solution[inds[i][j] as usize] * vals[i][j];
            }
            if viol >= 10.0 * ms.mipdata_.feastol {
                let len = inds[i].len() as HighsInt;
                let rhs = rhss[i];
                detach!(ms.mipdata_.cutpools[cutpool], |p| p.add_cut(
                    ms,
                    &mut inds[i][..],
                    &mut vals[i][..],
                    len,
                    rhs,
                    false,
                    true,
                    true,
                    false,
                ));
            }
        }
        self.already_tried = true;
    }

    /// HiGHS latest (#2990, dc83470a9): HighsMachineSchedSeparator.cpp::HighsMachineSchedSeparator::findSingleMachineScheduleClique (lines 19-284)
    /// HiGHS 1.16 (#3591): HighsMachineSchedSeparator.cpp::findSingleMachineScheduleClique (lines 19-306)
    /// C++: `bool findSingleMachineScheduleClique(std::vector<std::vector<double>>& vals, std::vector<std::vector<HighsInt>>& inds, std::vector<double>& rhss, const HighsDomain& globaldom, const HighsMipSolver& mipsolver)`
    // PORT NOTE: `globaldom` is `ms.mipdata_.domains[0]`, read fresh at each use; the
    // borrow is not held across `add_clique` (which takes `&mut ms`).
    pub fn find_single_machine_schedule_clique(
        &mut self,
        ms: &mut HighsMipSolver,
        vals: &mut Vec<Vec<f64>>,
        inds: &mut Vec<Vec<HighsInt>>,
        rhss: &mut Vec<f64>,
    ) -> bool {
        let mut largest_degree: HighsInt = 0;
        let mut largest_degree_col: HighsInt = -1;
        let mut num_rows: HighsInt = 0;
        let mut degrees: Vec<HighsInt> = vec![0; ms.num_col() as usize];
        let max_rows = cpp_min(50000 as HighsInt, 2 * ms.num_row());

        // C++: HighsHashTable<std::pair<HighsInt,HighsInt>, std::tuple<double,HighsInt,ArcType>> adjacency(maxRows + 2);
        let mut adjacency: HighsHashTable<(HighsInt, HighsInt), (f64, HighsInt, ArcType)> =
            HighsHashTable::with_min_capacity((max_rows + 2) as u64);
        // C++: HighsHashTable<std::pair<HighsInt,HighsInt>, std::vector<HighsInt>> arcToBin;
        let mut arc_to_bin: HighsHashTable<(HighsInt, HighsInt), Vec<HighsInt>> =
            HighsHashTable::new();
        // C++: HighsHashTable<std::tuple<HighsInt,HighsInt,HighsInt>, uint8_t> jobOrder;
        let mut job_order: HighsHashTable<(HighsInt, HighsInt, HighsInt), u8> =
            HighsHashTable::new();

        // C++ lambda addEntry (lines 46-82)
        let mut add_entry =
            |pos_col: HighsInt,
             neg_col: HighsInt,
             bin_col: HighsInt,
             p: f64,
             t: ArcType,
             adjacency: &mut HighsHashTable<(HighsInt, HighsInt), (f64, HighsInt, ArcType)>,
             degrees: &mut Vec<HighsInt>,
             largest_degree: &mut HighsInt,
             largest_degree_col: &mut HighsInt,
             job_order: &mut HighsHashTable<(HighsInt, HighsInt, HighsInt), u8>,
             arc_to_bin: &mut HighsHashTable<(HighsInt, HighsInt), Vec<HighsInt>>| {
                // My_ji + s_i - s_j >= p_ji
                // y_ji = val -> s_i >= s_j + p_ji - M * val
                // Make an arc from negCol (j) to posCol (i)
                if p < 0.0 {
                    return;
                }
                let key = (neg_col, pos_col);
                if let Some(it) = adjacency.find(&key) {
                    if it.0 < p {
                        // C++: adjacency[{negCol, posCol}] = std::make_tuple(p, binCol, t);
                        // index_or_insert returns &mut V; overwrite.
                        let slot = adjacency.index_or_insert(&key);
                        *slot = (p, bin_col, t);
                    }
                } else {
                    degrees[pos_col as usize] += 1;
                    if degrees[pos_col as usize] > *largest_degree
                        || (degrees[pos_col as usize] == *largest_degree
                            && pos_col < *largest_degree_col)
                    {
                        *largest_degree_col = pos_col;
                        *largest_degree = degrees[pos_col as usize];
                    }
                    adjacency.insert(HighsHashTableEntry::new(key, (p, bin_col, t)));
                }
                // Store binaries to later check if any two binaries
                // imply the same ordering of two jobs.
                let u = cpp_min(neg_col, pos_col);
                let v = cpp_max(neg_col, pos_col);
                let jo_key = (u, v, bin_col);
                if job_order.find(&jo_key).is_none() {
                    // C++: jobOrder[{u, v, binCol}] = 0;
                    let slot = job_order.index_or_insert(&jo_key);
                    *slot = 0;
                }
                if neg_col < pos_col {
                    // C++: arcToBin[{negCol, posCol}].emplace_back(binCol);
                    let slot = arc_to_bin.index_or_insert(&key);
                    slot.push(bin_col);
                    // C++: jobOrder[{negCol, posCol, binCol}] |= t == kImplicationWhenOne ? 2 : 1;
                    let jo_key2 = (neg_col, pos_col, bin_col);
                    let slot2 = job_order.index_or_insert(&jo_key2);
                    *slot2 |= if t == ArcType::kImplicationWhenOne {
                        2
                    } else {
                        1
                    };
                } else {
                    // C++: jobOrder[{posCol, negCol, binCol}] |= t == kImplicationWhenOne ? 8 : 4;
                    let jo_key2 = (pos_col, neg_col, bin_col);
                    let slot2 = job_order.index_or_insert(&jo_key2);
                    *slot2 |= if t == ArcType::kImplicationWhenOne {
                        8
                    } else {
                        4
                    };
                }
            };

        // C++ lambda tryAddArc (lines 84-95)
        let mut try_add_arc =
            |pos_col: HighsInt,
             neg_col: HighsInt,
             bin_col: HighsInt,
             rhs_0: f64,
             rhs_1: f64,
             adjacency: &mut HighsHashTable<(HighsInt, HighsInt), (f64, HighsInt, ArcType)>,
             degrees: &mut Vec<HighsInt>,
             largest_degree: &mut HighsInt,
             largest_degree_col: &mut HighsInt,
             job_order: &mut HighsHashTable<(HighsInt, HighsInt, HighsInt), u8>,
             arc_to_bin: &mut HighsHashTable<(HighsInt, HighsInt), Vec<HighsInt>>,
             num_rows: &mut HighsInt,
             add_entry: &mut dyn FnMut(
                HighsInt,
                HighsInt,
                HighsInt,
                f64,
                ArcType,
                &mut HighsHashTable<(HighsInt, HighsInt), (f64, HighsInt, ArcType)>,
                &mut Vec<HighsInt>,
                &mut HighsInt,
                &mut HighsInt,
                &mut HighsHashTable<(HighsInt, HighsInt, HighsInt), u8>,
                &mut HighsHashTable<(HighsInt, HighsInt), Vec<HighsInt>>,
            )| {
                if rhs_0 > 0.0 || rhs_1 > 0.0 {
                    if rhs_0 > rhs_1 {
                        add_entry(
                            pos_col,
                            neg_col,
                            bin_col,
                            rhs_0,
                            ArcType::kImplicationWhenZero,
                            adjacency,
                            degrees,
                            largest_degree,
                            largest_degree_col,
                            job_order,
                            arc_to_bin,
                        );
                    } else {
                        add_entry(
                            pos_col,
                            neg_col,
                            bin_col,
                            rhs_1,
                            ArcType::kImplicationWhenOne,
                            adjacency,
                            degrees,
                            largest_degree,
                            largest_degree_col,
                            job_order,
                            arc_to_bin,
                        );
                    }
                    *num_rows += 1;
                }
            };

        for row in 0..ms.num_row() {
            let row_lower = mip_model!(ms).row_lower_[row as usize];
            let row_upper = mip_model!(ms).row_upper_[row as usize];
            if row_lower == row_upper {
                continue;
            }
            let start = ms.mipdata_.a_rstart_[row as usize];
            let end = ms.mipdata_.a_rstart_[(row + 1) as usize];
            if end - start != 3 {
                continue;
            }
            let mut machine_sched_row = true;
            let mut pos_cont_col: HighsInt = -1;
            let mut neg_cont_col: HighsInt = -1;
            let mut bin_col: HighsInt = -1;
            let mut bin_coef: f64 = 0.0;
            let mut i = start;
            while i != end {
                let col = ms.mipdata_.a_rindex_[i as usize];
                if ms.mipdata_.domains[0].col_lower_[col as usize] == -kHighsInf {
                    machine_sched_row = false;
                    break;
                }
                if ms.mipdata_.domains[0].is_binary(ms, col) {
                    if bin_col != -1 {
                        machine_sched_row = false;
                        break;
                    }
                    bin_col = col;
                    bin_coef = ms.mipdata_.a_rvalue_[i as usize];
                } else if ms.mipdata_.a_rvalue_[i as usize] == -1.0 {
                    if neg_cont_col != -1 {
                        machine_sched_row = false;
                        break;
                    }
                    neg_cont_col = col;
                } else if ms.mipdata_.a_rvalue_[i as usize] == 1.0 {
                    if pos_cont_col != -1 {
                        machine_sched_row = false;
                        break;
                    }
                    pos_cont_col = col;
                } else {
                    machine_sched_row = false;
                    break;
                }
                i += 1;
            }
            if !machine_sched_row || bin_col == -1 || neg_cont_col == -1 || pos_cont_col == -1 {
                continue;
            }
            // We want to put the row into form:
            // My_ji + s_i - s_j >= p_ji, p_ji >= 0
            // y_ji = 1 -> s_i >= s_j + p_ji - M
            if row_upper != kHighsInf {
                // Given My_ij + s_i - s_j <= d
                // The row becomes (after multiplying by -1):
                // -My_ij + s_j - s_i >= -d
                // Add implication s_j >= s_i + p_ij + M, p_ij + M > 0 when binCol = 1
                // Add implication s_j >= s_i + p_ij, p_ij > 0 when binCol = 0
                try_add_arc(
                    neg_cont_col,
                    pos_cont_col,
                    bin_col,
                    -row_upper,
                    -row_upper + bin_coef,
                    &mut adjacency,
                    &mut degrees,
                    &mut largest_degree,
                    &mut largest_degree_col,
                    &mut job_order,
                    &mut arc_to_bin,
                    &mut num_rows,
                    &mut add_entry,
                );
            }
            if row_lower != -kHighsInf {
                // Given My_ij + s_i - s_j >= d
                // Add implication s_i >= s_j + p_ji - M, p_ji - M > 0 when binCol = 1
                // Add implication s_i >= s_j + p_ji, p_ji > 0 when binCol = 0
                try_add_arc(
                    pos_cont_col,
                    neg_cont_col,
                    bin_col,
                    row_lower,
                    row_lower - bin_coef,
                    &mut adjacency,
                    &mut degrees,
                    &mut largest_degree,
                    &mut largest_degree_col,
                    &mut job_order,
                    &mut arc_to_bin,
                    &mut num_rows,
                    &mut add_entry,
                );
            }
            if num_rows >= max_rows {
                break;
            }
        }

        // Extract binary variables that imply the same job ordering
        if !ms.parallel_lock_active() {
            let mut clique: Vec<HighsCliqueTable_CliqueVar> =
                vec![HighsCliqueTable_CliqueVar::default(); 2];
            for entry in arc_to_bin.iter() {
                let arc = *entry.key();
                let mut base_bin_col: HighsInt = -1;
                let mut base_arc_implication_when_one = false;
                for &bin_col in entry.value() {
                    let mut arc_implication_when_one = false;
                    let mut implies_order = false;
                    let jo_key = (arc.0, arc.1, bin_col);
                    let jo_val = *job_order.index_or_insert(&jo_key);
                    if (jo_val & 9) == 9 {
                        implies_order = true;
                    } else if (jo_val & 6) == 6 {
                        implies_order = true;
                        arc_implication_when_one = true;
                    }
                    if !implies_order || bin_col == base_bin_col {
                        continue;
                    }
                    if base_bin_col == -1 {
                        base_bin_col = bin_col;
                        base_arc_implication_when_one = arc_implication_when_one;
                        continue;
                    }
                    // Add cliques x1 + ~x2 <= 1 and ~x1 + x2 <= 1 which together imply x1
                    // == x2 (depending on signs may also have x1 == ~x2)
                    let stay_clique_var = HighsCliqueTable_CliqueVar::new(
                        base_bin_col,
                        base_arc_implication_when_one as HighsInt,
                    );
                    let subst_clique_var = HighsCliqueTable_CliqueVar::new(
                        bin_col,
                        arc_implication_when_one as HighsInt,
                    );
                    clique[0] = stay_clique_var;
                    clique[1] = subst_clique_var.complement();
                    // C++: addClique(mipsolver, clique.data(), 2), with the defaults
                    // equality = false, origin = kHighsIInf.
                    detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                        ms,
                        &mut clique[..],
                        false,
                        kHighsIInf
                    ));
                    clique[0] = stay_clique_var.complement();
                    clique[1] = subst_clique_var;
                    detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                        ms,
                        &mut clique[..],
                        false,
                        kHighsIInf
                    ));
                    if ms.mipdata_.domains[0].infeasible() {
                        return false;
                    }
                }
            }
        }

        // Skip any clique smaller than size 3
        if num_rows <= 5 {
            return false;
        }

        // Greedily search neighbours of largest degree column for a double-sided
        // clique (corresponds to a single machine schedule)
        let mut potential_neighbours: Vec<HighsInt> = Vec::with_capacity(largest_degree as usize);
        for entry in adjacency.iter() {
            let arc = *entry.key();
            let col = arc.1;
            if col == largest_degree_col {
                potential_neighbours.push(arc.0);
            }
        }
        pdqsort_by(&mut potential_neighbours, |a, b| {
            degrees[*a as usize] > degrees[*b as usize]
        });

        let mut neighbours: Vec<HighsInt> = Vec::with_capacity((largest_degree + 1) as usize);
        neighbours.push(largest_degree_col);
        let mut release_date = ms.mipdata_.domains[0].col_lower_[largest_degree_col as usize];
        let mut processing_times: Vec<f64> = vec![kHighsInf; (largest_degree + 1) as usize];
        // Iterate over potential neighbours and check validity
        // Greedily add neighbours if they're valid
        for &col in &potential_neighbours {
            let mut valid_neighbour = true;
            for &neighbour in &neighbours {
                let from_arc = adjacency.find(&(col, neighbour));
                let to_arc = adjacency.find(&(neighbour, col));
                // Need to verify that the to and fromArc exist and are opposites
                if to_arc.is_none()
                    || from_arc.is_none()
                    || from_arc.unwrap().1 != to_arc.unwrap().1
                    || from_arc.unwrap().2 == to_arc.unwrap().2
                {
                    valid_neighbour = false;
                    break;
                }
            }
            if !valid_neighbour {
                continue;
            }
            let mut min_processing_time = kHighsInf;
            // Extract the processing times from the arcs
            // Turn pair-wise dependent times p_ji into p_j = min{p_ji : i \in N / j}
            for i in 0..neighbours.len() {
                let neighbour = neighbours[i];
                let from_arc = adjacency.find(&(col, neighbour)).unwrap();
                let to_arc = adjacency.find(&(neighbour, col)).unwrap();
                min_processing_time = cpp_min(min_processing_time, from_arc.0);
                processing_times[i] = cpp_min(processing_times[i], to_arc.0);
            }
            release_date = cpp_min(
                ms.mipdata_.domains[0].col_lower_[col as usize],
                release_date,
            );
            processing_times[neighbours.len()] = min_processing_time;
            neighbours.push(col);
        }
        if neighbours.len() < 3 {
            return false;
        }

        // Now populate the actual inequalities
        vals.resize(neighbours.len(), Vec::new());
        inds.resize(neighbours.len(), Vec::new());
        rhss.resize(neighbours.len(), 0.0);
        // HiGHS 1.16 #3591 (no HiGHS 1.16 flag: reached only when mip_machine_sched_sepa is on, see
        // separation.rs:61): HighsMachineSchedSeparator.cpp lines 261-299 (1.16.0)
        // C++: std::vector<HighsInt> binColCutPositions(mipsolver.numCol(), -1);
        let mut bin_col_cut_positions: Vec<HighsInt> = vec![-1; ms.num_col() as usize];
        for i in 0..neighbours.len() {
            vals[i].reserve(neighbours.len());
            inds[i].reserve(neighbours.len());
            rhss[i] -= release_date;
            let col = neighbours[i];
            for j in 0..neighbours.len() {
                if i == j {
                    continue;
                }
                let neighbour = neighbours[j];
                let to_arc = adjacency.find(&(neighbour, col)).unwrap();
                let bin_col = to_arc.1;
                let mut val = processing_times[j];
                if to_arc.2 == ArcType::kImplicationWhenZero {
                    rhss[i] -= val;
                    val = -val;
                }
                let pos = bin_col_cut_positions[bin_col as usize];
                if pos == -1 {
                    bin_col_cut_positions[bin_col as usize] = inds[i].len() as HighsInt;
                    inds[i].push(bin_col);
                    vals[i].push(val);
                } else {
                    vals[i][pos as usize] += val;
                }
            }
            // Drop any potential zeros (if the same binary sets the order of many jobs)
            let mut j: usize = 0;
            while j != inds[i].len() {
                bin_col_cut_positions[inds[i][j] as usize] = -1;
                if vals[i][j] == 0.0 {
                    let last = inds[i].len() - 1;
                    inds[i].swap(j, last);
                    vals[i].swap(j, last);
                    inds[i].pop();
                    vals[i].pop();
                } else {
                    j += 1;
                }
            }
            // Put the job start time on the LHS
            inds[i].push(col);
            vals[i].push(-1.0);
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x1_mss_a_find_clique_toy() {
        // Toy 3-job disjunctive model: continuous columns s0, s1, s2 with bounds [0, 20];
        // binary columns y01, y02, y12 (bounds [0,1], integral).
        // For each pair (i,j) with i<j, with binary y = y_ij, two rows:
        //   row A: s_j - s_i - 10*y >= -7      (y=1 gives s_j >= s_i + 3)
        //   row B: s_i - s_j + 10*y >= 4       (y=0 gives s_i >= s_j + 4)
        // Total 6 rows, each with exactly 3 nonzeros.

        let mut lp = HighsLp::default();
        lp.num_col_ = 6;
        lp.num_row_ = 6;
        lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        lp.col_upper_ = vec![20.0, 20.0, 20.0, 1.0, 1.0, 1.0];
        lp.col_cost_ = vec![0.0; 6];
        lp.integrality_ = vec![
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = vec![-7.0, 4.0, -7.0, 4.0, -7.0, 4.0];
        lp.row_upper_ = vec![kHighsInf; 6];

        // Build CSR for the 6 rows (each has 3 nonzeros).
        // Row 0 (pair 0,1, A): s1 - s0 - 10*y01 >= -7  => cols 1,0,3 vals 1,-1,-10
        // Row 1 (pair 0,1, B): s0 - s1 + 10*y01 >= 4   => cols 0,1,3 vals 1,-1,10
        // Row 2 (pair 0,2, A): s2 - s0 - 10*y02 >= -7  => cols 2,0,4 vals 1,-1,-10
        // Row 3 (pair 0,2, B): s0 - s2 + 10*y02 >= 4   => cols 0,2,4 vals 1,-1,10
        // Row 4 (pair 1,2, A): s2 - s1 - 10*y12 >= -7  => cols 2,1,5 vals 1,-1,-10
        // Row 5 (pair 1,2, B): s1 - s2 + 10*y12 >= 4   => cols 1,2,5 vals 1,-1,10
        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.start_ = vec![0, 3, 6, 9, 12, 15, 18];
        lp.a_matrix_.index_ = vec![
            1, 0, 3, // row 0
            0, 1, 3, // row 1
            2, 0, 4, // row 2
            0, 2, 4, // row 3
            2, 1, 5, // row 4
            1, 2, 5, // row 5
        ];
        lp.a_matrix_.value_ = vec![
            1.0, -1.0, -10.0, // row 0
            1.0, -1.0, 10.0, // row 1
            1.0, -1.0, -10.0, // row 2
            1.0, -1.0, 10.0, // row 3
            1.0, -1.0, -10.0, // row 4
            1.0, -1.0, 10.0, // row 5
        ];

        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);

        // Fill mipdata_ AR arrays (row-wise) from the LP matrix.
        ms.mipdata_.a_rstart_ = lp.a_matrix_.start_.clone();
        ms.mipdata_.a_rindex_ = lp.a_matrix_.index_.clone();
        ms.mipdata_.a_rvalue_ = lp.a_matrix_.value_.clone();

        // Build one domain with the column bounds.
        let mut dom = HighsDomain::default();
        dom.col_lower_ = lp.col_lower_.clone();
        dom.col_upper_ = lp.col_upper_.clone();
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut sep = HighsMachineSchedSeparator::new(&ms);
        let mut vals: Vec<Vec<f64>> = Vec::new();
        let mut inds: Vec<Vec<HighsInt>> = Vec::new();
        let mut rhss: Vec<f64> = Vec::new();
        assert!(sep.find_single_machine_schedule_clique(&mut ms, &mut vals, &mut inds, &mut rhss));
        assert_eq!(vals.len(), 3);
        assert_eq!(inds.len(), 3);
        assert_eq!(rhss.len(), 3);
        // Expected values traced by hand through the C++ (HighsMachineSchedSeparator.cpp
        // lines 97-281). largestDegreeCol = 0; the adjacency table iterates (2,0) before
        // (1,0) (checked against the C++ HighsHashTable), so neighbours = [0, 2, 1] and
        // processingTimes = [3, 4, 3].
        assert_eq!(inds, vec![vec![4, 3, 0], vec![4, 5, 2], vec![3, 5, 1]]);
        assert_eq!(
            vals,
            vec![
                vec![-4.0, -3.0, -1.0],
                vec![3.0, 3.0, -1.0],
                vec![3.0, -4.0, -1.0]
            ]
        );
        assert_eq!(rhss, vec![-7.0, 0.0, -4.0]);
    }

    #[test]
    fn x1_mss_a_hash_keys_match_cpp() {
        // Values from HighsHashHelpers::hash (HighsHash.h at dc83470a9) built with g++ and
        // libstdc++, where std::tuple<int, int, int> is laid out last element first.
        assert_eq!(
            (0 as HighsInt, 1 as HighsInt, 3 as HighsInt).highs_hash(),
            832224563789260942
        );
        assert_eq!(
            (5 as HighsInt, -1 as HighsInt, 7 as HighsInt).highs_hash(),
            10372906951778696660
        );
        assert_eq!(
            (2 as HighsInt, 1 as HighsInt).highs_hash(),
            3891188638457094794
        );
    }

    #[test]
    fn h16_msched_on_merges_duplicate_bin_and_drops_zeros() {
        // HiGHS 1.16 #3591 (no HiGHS 1.16 flag): HighsMachineSchedSeparator.cpp lines 261-299 (1.16.0).
        // Three jobs s0, s1, s2 and one shared binary y (col 3). For each pair i < j:
        //   A: s_j - s_i - 10*y >= -7   (y=1 gives s_j >= s_i + 3)
        //   B: s_i - s_j + 10*y >= 3    (y=0 gives s_i >= s_j + 3)
        // Every arc has processing time 3, and every job's cut sees y on each arc.
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 6;
        lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0];
        lp.col_upper_ = vec![20.0, 20.0, 20.0, 1.0];
        lp.col_cost_ = vec![0.0; 4];
        lp.integrality_ = vec![
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = vec![-7.0, 3.0, -7.0, 3.0, -7.0, 3.0];
        lp.row_upper_ = vec![kHighsInf; 6];

        lp.a_matrix_.format_ = MatrixFormat::kRowwise;
        lp.a_matrix_.start_ = vec![0, 3, 6, 9, 12, 15, 18];
        lp.a_matrix_.index_ = vec![
            1, 0, 3, // row 0: pair (0,1) A
            0, 1, 3, // row 1: pair (0,1) B
            2, 0, 3, // row 2: pair (0,2) A
            0, 2, 3, // row 3: pair (0,2) B
            2, 1, 3, // row 4: pair (1,2) A
            1, 2, 3, // row 5: pair (1,2) B
        ];
        lp.a_matrix_.value_ = vec![
            1.0, -1.0, -10.0, // row 0
            1.0, -1.0, 10.0, // row 1
            1.0, -1.0, -10.0, // row 2
            1.0, -1.0, 10.0, // row 3
            1.0, -1.0, -10.0, // row 4
            1.0, -1.0, 10.0, // row 5
        ];

        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.mipdata_.a_rstart_ = lp.a_matrix_.start_.clone();
        ms.mipdata_.a_rindex_ = lp.a_matrix_.index_.clone();
        ms.mipdata_.a_rvalue_ = lp.a_matrix_.value_.clone();
        let mut dom = HighsDomain::default();
        dom.col_lower_ = lp.col_lower_.clone();
        dom.col_upper_ = lp.col_upper_.clone();
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut sep = HighsMachineSchedSeparator::new(&ms);
        let mut vals: Vec<Vec<f64>> = Vec::new();
        let mut inds: Vec<Vec<HighsInt>> = Vec::new();
        let mut rhss: Vec<f64> = Vec::new();
        assert!(sep.find_single_machine_schedule_clique(&mut ms, &mut vals, &mut inds, &mut rhss));
        assert_eq!(vals.len(), 3);

        // Each row ends with its job's start time (-1). The row order depends on the
        // neighbour order, so look each row up by its job.
        let row_of = |job: HighsInt| -> usize {
            inds.iter().position(|r| *r.last().unwrap() == job).unwrap()
        };
        // Job 0: the two arcs into job 0 both carry y (-3 each), merged to -6.
        let r0 = row_of(0);
        assert_eq!(inds[r0], vec![3, 0]);
        assert_eq!(vals[r0], vec![-6.0, -1.0]);
        assert_eq!(rhss[r0], -6.0);
        // Job 1: y enters with +3 from job 0 and -3 from job 2; the zero is dropped.
        let r1 = row_of(1);
        assert_eq!(inds[r1], vec![1]);
        assert_eq!(vals[r1], vec![-1.0]);
        assert_eq!(rhss[r1], -3.0);
        // Job 2: both arcs into job 2 carry y (+3 each), merged to +6.
        let r2 = row_of(2);
        assert_eq!(inds[r2], vec![3, 2]);
        assert_eq!(vals[r2], vec![6.0, -1.0]);
        assert_eq!(rhss[r2], 0.0);
    }

    #[test]
    fn h16_msched_off_separator_not_registered() {
        use crate::mip::HighsSeparatorBox;
        use crate::mip::separation::HighsSeparation;
        // With mip_machine_sched_sepa off (default) the machine-scheduling separator is
        // never built, so find_single_machine_schedule_clique is never reached and the
        // separation set stays the three 1.15.1 separators.
        let mut ms = HighsMipSolver::default();
        ms.options_mip_.mip_ti_precedence_sepa = false;
        ms.options_mip_.mip_ti_list_sched = false;
        ms.options_mip_.mip_netdesign_sepa = false;
        ms.options_mip_.mip_netdesign_heur = false;
        ms.options_mip_.mip_netdesign_cutset = false;
        assert!(!ms.options_mip_.mip_machine_sched_sepa);
        let sepa = HighsSeparation::new(&ms);
        assert_eq!(sepa.separators.len(), 3);
        assert!(
            !sepa
                .separators
                .iter()
                .any(|s| matches!(s, HighsSeparatorBox::MachineSched(_)))
        );
    }
}
