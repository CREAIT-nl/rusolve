//! HiGHS source: presolve/HighsSymmetry.cpp (lines 1-1000)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! PORT NOTE: the six functions whose C++ reaches `HighsDomain::changeBound`, `propagate`,
//! `variableType` or `isGlobalBinary` take `ms: &mut HighsMipSolver` (see the note on each one).
//! The C++ `StabilizerOrbits::symmetries` back-pointer is replaced by an explicit `symmetries` argument.

use crate::prelude::*;
use std::sync::Arc;

/// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::orbitalFixingForFullOrbitope, lambda `i_fixed`
/// C++ returns `numDynamicRows` when no index is found; `colj0.len()` is `numDynamicRows` here.
fn i_fixed(colj0: &[i8], colj1: &[i8]) -> HighsInt {
    let n = colj0.len();
    for i in 0..n {
        if colj0[i] != -1 && colj1[i] != -1 && colj0[i] != colj1[i] {
            return i as HighsInt;
        }
    }
    n as HighsInt
}

/// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::orbitalFixingForFullOrbitope, lambda `i_discr`
fn i_discr(colj0: &[i8], colj1: &[i8], i_f: HighsInt) -> HighsInt {
    let mut i = i_f;
    while i >= 0 {
        if colj0[i as usize] != 0 && colj1[i as usize] != 1 {
            return i;
        }
        i -= 1;
    }
    -1
}

/// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::orbitalFixingForFullOrbitope, the
/// `colj0[k] += (isFree & colj1[k]) + isFree` update for `k < upto`.
/// PORT NOTE: int8_t wraparound is kept with `wrapping_add`.
fn update_free_entries(colj0: &mut [i8], colj1: &[i8], upto: usize) {
    for k in 0..upto {
        let is_free: i8 = i8::from(colj0[k] == -1);
        colj0[k] = colj0[k].wrapping_add((is_free & colj1[k]).wrapping_add(is_free));
    }
}

impl HighsSymmetryDetection {
    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::removeFixPoints (lines 24-97)
    /// C++: `void removeFixPoints()`
    pub fn remove_fix_points(&mut self) {
        self.gend.resize(self.num_vertices as usize, 0);
        for i in 0..self.num_vertices as usize {
            let start = self.gstart[i] as usize;
            let end = self.gstart[i + 1] as usize;
            // cellSize(vertexToCell[edge.first]) > 1, with disjoint field borrows
            let cpl = &self.current_partition_links;
            let vtc = &self.vertex_to_cell;
            let k = cpp_partition(&mut self.gedge[start..end], |edge| {
                let cell = vtc[edge.0 as usize];
                cpl[cell as usize] - cell > 1
            });
            self.gend[i] = (start + k) as HighsInt;
            debug_assert!(self.gend[i] >= self.gstart[i] && self.gend[i] <= self.gstart[i + 1]);
        }

        let mut unit_cell_index: HighsInt = self.num_vertices;
        {
            let cpl = &self.current_partition_links;
            let vtc = &mut self.vertex_to_cell;
            // std::remove_if: the predicate removes unit cells; retain keeps the others
            self.current_partition.retain(|&vertex| {
                let cell = vtc[vertex as usize];
                if cpl[cell as usize] - cell == 1 {
                    unit_cell_index -= 1;
                    vtc[vertex as usize] = unit_cell_index;
                    false
                } else {
                    true
                }
            });
        }

        for i in 0..self.num_vertices as usize {
            if self.gend[i] == self.gstart[i + 1] {
                continue;
            }
            for j in self.gend[i] as usize..self.gstart[i + 1] as usize {
                let first = self.gedge[j].0;
                let cell = self.vertex_to_cell[first as usize];
                self.gedge[j].0 = cell;
            }
        }

        if (self.current_partition.len() as HighsInt) < self.num_vertices {
            self.num_vertices = self.current_partition.len() as HighsInt;
            if self.num_vertices == 0 {
                self.num_active_cols = 0;
                return;
            }
            let n = self.num_vertices as usize;
            self.current_partition_links.resize(n, 0);
            // cellInRefinementQueue.assign(numVertices, false)
            self.cell_in_refinement_queue = vec![0u8; n];
            debug_assert!(self.refinement_queue.is_empty());
            self.refinement_queue.clear();
            let mut cell_start: HighsInt = 0;
            let mut cell_number: HighsInt = 0;
            for i in 0..self.num_vertices {
                let vertex = self.current_partition[i as usize];
                // if the cell number is different to the current cell number this is the
                // start of a new cell
                if cell_number != self.vertex_to_cell[vertex as usize] {
                    // remember the number of this cell to identify its end
                    cell_number = self.vertex_to_cell[vertex as usize];
                    // set the link of the cell start to point to its end
                    self.current_partition_links[cell_start as usize] = i;
                    // remember start of this cell
                    cell_start = i;
                }
                // correct the vertexToCell array to store the start index of the cell, not its number
                self.update_cell_membership(i, cell_start, false);
            }

            // set the column partition link of the last started cell to point past the end
            debug_assert!(!self.current_partition_links.is_empty());
            self.current_partition_links[cell_start as usize] = self.num_vertices;

            // PORT NOTE: std::partition_point is a binary search; on the partitioned range both
            // return the unique boundary, so slice::partition_point gives the same index.
            let num_col = self.num_col;
            self.num_active_cols =
                self.current_partition.partition_point(|&v| v < num_col) as HighsInt;
        } else {
            self.num_active_cols = self.num_col;
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::initializeGroundSet (lines 756-770)
    /// C++: `void initializeGroundSet()`
    pub fn initialize_ground_set(&mut self) {
        self.vertex_ground_set = self.current_partition.clone();
        pdqsort_by(&mut self.vertex_ground_set, |a: &HighsInt, b: &HighsInt| {
            a < b
        });
        self.vertex_position.resize(self.vertex_to_cell.len(), -1);
        for i in 0..self.num_vertices as usize {
            let v = self.vertex_ground_set[i] as usize;
            self.vertex_position[v] = i as HighsInt;
        }

        let n = self.num_vertices as usize;
        self.orbit_partition.resize(n, 0);
        // std::iota(orbitPartition.begin(), orbitPartition.end(), 0)
        for (i, o) in self.orbit_partition.iter_mut().enumerate() {
            *o = i as HighsInt;
        }
        // orbitSize.assign(numVertices, 1)
        self.orbit_size = vec![1; n];

        self.automorphisms
            .resize((self.num_vertices * 64) as usize, 0);
        self.num_automorphisms = 0;
        self.curr_node_certificate.reserve(n);
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::mergeOrbits (lines 772-789)
    /// C++: `bool mergeOrbits(HighsInt v1, HighsInt v2)`
    pub fn merge_orbits(&mut self, v1: HighsInt, v2: HighsInt) -> bool {
        if v1 == v2 {
            return false;
        }

        let orbit1 = self.get_orbit(v1);
        let orbit2 = self.get_orbit(v2);

        if orbit1 == orbit2 {
            return false;
        }

        if orbit1 < orbit2 {
            self.orbit_partition[orbit2 as usize] = orbit1;
            self.orbit_size[orbit1 as usize] += self.orbit_size[orbit2 as usize];
        } else {
            self.orbit_partition[orbit1 as usize] = orbit2;
            self.orbit_size[orbit2 as usize] += self.orbit_size[orbit1 as usize];
        }

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::getOrbit (lines 791-809)
    /// C++: `HighsInt getOrbit(HighsInt vertex)`
    pub fn get_orbit(&mut self, vertex: HighsInt) -> HighsInt {
        let mut i = self.vertex_position[vertex as usize];
        let mut orbit = self.orbit_partition[i as usize];
        if orbit != self.orbit_partition[orbit as usize] {
            loop {
                self.link_compression_stack.push(i);
                i = orbit;
                orbit = self.orbit_partition[orbit as usize];
                if orbit == self.orbit_partition[orbit as usize] {
                    break;
                }
            }

            loop {
                i = self.link_compression_stack.pop().unwrap();
                self.orbit_partition[i as usize] = orbit;
                if self.link_compression_stack.is_empty() {
                    break;
                }
            }
        }

        orbit
    }

    /// `HighsHashHelpers::sparse_combine32(vertexHash[vertex], cell, value)` (1.15.1), or
    /// `sparse_combine` on the `u64` hash with `mip_sym_hash64` (1.16.0 #3535,
    /// HighsSymmetry.cpp lines 736-737 and 756-757 at v1.16.0).
    fn combine_vertex_hash(&mut self, vertex: HighsInt, cell: HighsInt, value: u64) {
        let hash64 = self.mip_sym_hash64;
        let h = self.vertex_hash.index_or_insert(&vertex);
        if hash64 {
            sparse_combine_value(h, cell, value);
        } else {
            let mut h32 = *h as u32;
            sparse_combine32(&mut h32, cell, value);
            *h = u64::from(h32);
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::initializeHashValues (lines 811-820)
    /// C++: `void initializeHashValues()`
    pub fn initialize_hash_values(&mut self) {
        for i in 0..self.num_vertices as usize {
            let cell = self.vertex_to_cell[i];
            let start = self.gstart[i] as usize;
            let end = self.gend[i] as usize;
            for j in start..end {
                let key = self.gedge[j].0;
                let val = self.gedge[j].1;
                self.combine_vertex_hash(key, cell, u64::from(val));
            }
            self.mark_cell_for_refinement(cell);
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::updateCellMembership (lines 822-845)
    /// C++: `bool updateCellMembership(HighsInt i, HighsInt cell, bool markForRefinement = true)`
    /// C++ default arguments: `mark_for_refinement = true`
    pub fn update_cell_membership(
        &mut self,
        i: HighsInt,
        cell: HighsInt,
        mark_for_refinement: bool,
    ) -> bool {
        let vertex = self.current_partition[i as usize];
        if self.vertex_to_cell[vertex as usize] != cell {
            // set new cell id
            self.vertex_to_cell[vertex as usize] = cell;
            if i != cell {
                self.current_partition_links[i as usize] = cell;
            }

            // update hashes of affected rows
            if mark_for_refinement {
                let start = self.gstart[vertex as usize] as usize;
                let end = self.gend[vertex as usize] as usize;
                for j in start..end {
                    let first = self.gedge[j].0;
                    let second = self.gedge[j].1;
                    let edge_destination_cell = self.vertex_to_cell[first as usize];
                    if self.cell_size(edge_destination_cell) == 1 {
                        continue;
                    }
                    self.combine_vertex_hash(first, cell, u64::from(second));
                    self.mark_cell_for_refinement(edge_destination_cell);
                }
            }

            return true;
        }

        false
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::splitCell (lines 847-885)
    /// C++: `bool splitCell(HighsInt cell, HighsInt splitPoint)`
    pub fn split_cell(&mut self, cell: HighsInt, split_point: HighsInt) -> bool {
        let vertex_split = self.current_partition[split_point as usize];
        // 1.16.0 #3535: `static_cast<u32>(getVertexHash(..))` (HighsSymmetry.cpp lines
        // 769-770); with `mip_sym_hash64` off the hashes are < 2^31 and the cast is exact.
        let h_split = self.get_vertex_hash(vertex_split) as u32;
        let vertex_cell = self.current_partition[cell as usize];
        let h_cell = self.get_vertex_hash(vertex_cell) as u32;

        let certificate_val: u32 = ((pair_hash(0, h_split, h_cell).wrapping_add(pair_hash(
            1,
            cell as u32,
            (self.current_partition_links[cell as usize] - split_point) as u32,
        )))
        .wrapping_add(pair_hash(
            2,
            split_point as u32,
            (split_point - cell) as u32,
        )) >> 32) as u32;

        // employ prefix pruning scheme as in bliss
        if !self.first_leave_certificate.is_empty() {
            let size = self.curr_node_certificate.len() as HighsInt;
            // C++ evaluates both factors of each product, so both reads are unconditional
            let first_eq = self.first_leave_prefix_len == size;
            let first_cert_eq = certificate_val == self.first_leave_certificate[size as usize];
            self.first_leave_prefix_len += i32::from(first_eq) * i32::from(first_cert_eq);

            let best_eq = self.best_leave_prefix_len == size;
            let best_cert_eq = certificate_val == self.best_leave_certificate[size as usize];
            self.best_leave_prefix_len += i32::from(best_eq) * i32::from(best_cert_eq);

            // if the node certificate is not a prefix of the first leave's certificate
            // and it comes lexicographically after the certificate value of the
            // lexicographically smallest leave certificate we prune the node
            if self.first_leave_prefix_len <= size && self.best_leave_prefix_len <= size {
                let diff_val = if self.best_leave_prefix_len == size {
                    certificate_val
                } else {
                    self.curr_node_certificate[self.best_leave_prefix_len as usize]
                };
                if diff_val > self.best_leave_certificate[self.best_leave_prefix_len as usize] {
                    return false;
                }
            }
        }

        self.current_partition_links[split_point as usize] =
            self.current_partition_links[cell as usize];
        self.current_partition_links[cell as usize] = split_point;
        self.cell_creation_stack.push(split_point);
        self.curr_node_certificate.push(certificate_val);

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::markCellForRefinement (lines 887-894)
    /// C++: `void markCellForRefinement(HighsInt cell)`
    pub fn mark_cell_for_refinement(&mut self, cell: HighsInt) {
        if self.cell_size(cell) == 1 || self.cell_in_refinement_queue[cell as usize] != 0 {
            return;
        }

        self.cell_in_refinement_queue[cell as usize] = 1;
        self.refinement_queue.push(cell);
        // std::push_heap(..., std::greater<HighsInt>())
        cpp_push_heap_by(&mut self.refinement_queue, |a: &HighsInt, b: &HighsInt| {
            a > b
        });
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::getVertexHash (lines 896-900)
    /// C++: `u32 getVertexHash(HighsInt v)` (1.15.1), `u64` since 1.16.0 #3535
    pub fn get_vertex_hash(&mut self, v: HighsInt) -> u64 {
        match self.vertex_hash.find(&v) {
            Some(h) => *h,
            None => 0,
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::partitionRefinement (lines 902-1000)
    /// HiGHS 1.16 (#3182): HighsSymmetry.cpp::partitionRefinement (lines 834-940)
    /// C++: `bool partitionRefinement()`
    pub fn partition_refinement(&mut self) -> bool {
        while !self.refinement_queue.is_empty() {
            // std::pop_heap(..., std::greater<HighsInt>())
            cpp_pop_heap_by(&mut self.refinement_queue, |a: &HighsInt, b: &HighsInt| {
                a > b
            });

            let mut cell_start: HighsInt = *self.refinement_queue.last().unwrap();
            let first_cell_start = cell_start;
            self.refinement_queue.pop();
            self.cell_in_refinement_queue[cell_start as usize] = 0;

            if self.cell_size(cell_start) == 1 {
                continue;
            }
            let cell_end = self.current_partition_links[cell_start as usize];
            debug_assert!(cell_end >= cell_start);

            // HiGHS 1.16 #3182 (h16_sym_refine): HighsSymmetry.cpp lines 848-849 (1.16.0)
            let sym_refine = self.h16_sym_refine;
            let start_creation_stack_size = if sym_refine {
                self.cell_creation_stack.len()
            } else {
                0
            };
            let start_certificate_size = if sym_refine {
                self.curr_node_certificate.len()
            } else {
                0
            };

            // first check which vertices do have updated hash values and put them to
            // the end of the partition
            let refine_start = {
                let vertex_hash = &self.vertex_hash;
                let k = cpp_partition(
                    &mut self.current_partition[cell_start as usize..cell_end as usize],
                    |v: &HighsInt| vertex_hash.find(v).is_none(),
                );
                cell_start + k as HighsInt
            };

            // if there are none there is nothing to refine
            if refine_start == cell_end {
                continue;
            }

            // sort the vertices that have updated hash values by their hash values
            {
                let vh = &mut self.vertex_hash;
                pdqsort_by(
                    &mut self.current_partition[refine_start as usize..cell_end as usize],
                    |v1: &HighsInt, v2: &HighsInt| {
                        let h1 = *vh.index_or_insert(v1);
                        let h2 = *vh.index_or_insert(v2);
                        h1 < h2
                    },
                );
            }

            // if not all vertices have updated hash values directly create the first
            // new cell at the start of the range that we want to refine
            if refine_start != cell_start {
                if !self.split_cell(cell_start, refine_start) {
                    // node can be pruned, make sure hash values are cleared and queue is empty
                    for &cell_id in self.refinement_queue.iter() {
                        self.cell_in_refinement_queue[cell_id as usize] = 0;
                    }
                    self.refinement_queue.clear();
                    self.vertex_hash.clear();
                    return false;
                }
                cell_start = refine_start;
                self.update_cell_membership(cell_start, cell_start, true);
            }

            // now update the remaining vertices
            let mut prune = false;
            debug_assert!(
                self.vertex_hash
                    .find(&self.current_partition[cell_start as usize])
                    .is_some()
            );
            // store value of first hash
            let first_vertex = self.current_partition[cell_start as usize];
            let mut last_hash: u64 = *self.vertex_hash.index_or_insert(&first_vertex);
            let mut i: HighsInt = cell_start + 1;
            while i < cell_end {
                let vertex = self.current_partition[i as usize];
                // get this vertex hash value
                let hash: u64 = *self.vertex_hash.index_or_insert(&vertex);

                if hash != last_hash {
                    // hash values do not match -> start of new cell
                    if !self.split_cell(cell_start, i) {
                        // refinement process yielded bad prefix of certificate -> node can be pruned
                        prune = true;
                        break;
                    }
                    cell_start = i;
                    // remember hash value of this new cell under lastHash
                    last_hash = hash;
                }

                // update membership of vertex to new cell
                self.update_cell_membership(i, cell_start, true);
                i += 1;
            }

            if prune {
                // node can be pruned, make sure hash values are cleared and queue is empty
                for &cell_id in self.refinement_queue.iter() {
                    self.cell_in_refinement_queue[cell_id as usize] = 0;
                }
                self.refinement_queue.clear();
                self.vertex_hash.clear();
                self.current_partition_links[first_cell_start as usize] = cell_end;

                // undo possibly incomplete changes done to the cells
                i -= 1;
                while i >= refine_start {
                    self.update_cell_membership(i, first_cell_start, false);
                    i -= 1;
                }

                if sym_refine {
                    // HiGHS 1.16 #3182 (h16_sym_refine): HighsSymmetry.cpp lines 922-930 (1.16.0)
                    self.cell_creation_stack
                        .resize(start_creation_stack_size, 0);
                    self.curr_node_certificate.resize(start_certificate_size, 0);
                    if !self.first_leave_certificate.is_empty() {
                        let start_certificate_len = start_certificate_size as HighsInt;
                        self.first_leave_prefix_len =
                            self.first_leave_prefix_len.min(start_certificate_len);
                        self.best_leave_prefix_len =
                            self.best_leave_prefix_len.min(start_certificate_len);
                    }
                }

                return false;
            }

            debug_assert!(self.current_partition_links[cell_start as usize] == cell_end);
        }

        self.vertex_hash.clear();

        true
    }
}

impl HighsSymmetries {
    /// HiGHS: HighsSymmetry.cpp::HighsSymmetries::clear (lines 99-110)
    /// C++: `void clear()`
    pub fn clear(&mut self) {
        self.permutation_columns.clear();
        self.permutations.clear();
        self.orbit_partition.clear();
        self.orbit_size.clear();
        self.column_position.clear();
        self.link_compression_stack.clear();
        self.column_to_orbitope.clear();
        self.orbitopes.clear();
        self.num_perms = 0;
        self.num_generators = 0;
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetries::mergeOrbits (lines 112-130)
    /// C++: `void mergeOrbits(HighsInt v1, HighsInt v2, std::vector<HighsInt> & orbitPartition, std::vector<HighsInt> & orbitSize, std::vector<HighsInt> & linkCompressionStack)`
    pub fn merge_orbits(
        &self,
        v1: HighsInt,
        v2: HighsInt,
        orbit_partition: &mut Vec<HighsInt>,
        orbit_size: &mut Vec<HighsInt>,
        link_compression_stack: &mut Vec<HighsInt>,
    ) {
        if v1 == v2 {
            return;
        }

        let orbit1 = self.get_orbit(v1, orbit_partition, link_compression_stack);
        let orbit2 = self.get_orbit(v2, orbit_partition, link_compression_stack);

        if orbit1 == orbit2 {
            return;
        }

        if orbit_size[orbit2 as usize] < orbit_size[orbit1 as usize] {
            orbit_partition[orbit2 as usize] = orbit1;
            orbit_size[orbit1 as usize] += orbit_size[orbit2 as usize];
        } else {
            orbit_partition[orbit1 as usize] = orbit2;
            orbit_size[orbit2 as usize] += orbit_size[orbit1 as usize];
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetries::getOrbit (lines 132-154)
    /// C++: `HighsInt getOrbit(HighsInt col, std::vector<HighsInt> & orbitPartition, std::vector<HighsInt> & linkCompressionStack)`
    pub fn get_orbit(
        &self,
        col: HighsInt,
        orbit_partition: &mut Vec<HighsInt>,
        link_compression_stack: &mut Vec<HighsInt>,
    ) -> HighsInt {
        let mut i = self.column_position[col as usize];
        if i == -1 {
            return -1;
        }

        let mut orbit = orbit_partition[i as usize];
        if orbit != orbit_partition[orbit as usize] {
            loop {
                link_compression_stack.push(i);
                i = orbit;
                orbit = orbit_partition[orbit as usize];
                if orbit == orbit_partition[orbit as usize] {
                    break;
                }
            }

            loop {
                i = link_compression_stack.pop().unwrap();
                orbit_partition[i as usize] = orbit;
                if link_compression_stack.is_empty() {
                    break;
                }
            }
        }

        orbit
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetries::propagateOrbitopes (lines 156-176)
    /// C++: `HighsInt propagateOrbitopes(HighsDomain & domain) const`
    /// PORT NOTE: signature changed because the domain needs ms (`changeBound`, `propagate` via orbital fixing).
    pub fn propagate_orbitopes(
        &self,
        ms: &mut HighsMipSolver,
        domain: &mut HighsDomain,
    ) -> HighsInt {
        if self.column_to_orbitope.size() == 0 || domain.get_branch_depth() == 0 {
            return 0;
        }

        // std::set<HighsInt>: ordered, duplicate-free
        let mut propagation_orbitopes: BTreeSet<HighsInt> = BTreeSet::new();
        let domchgstack = domain.get_domain_change_stack();
        for &pos in domain.get_branching_positions() {
            if let Some(&orbitope) = self
                .column_to_orbitope
                .find(&domchgstack[pos as usize].column)
            {
                propagation_orbitopes.insert(orbitope);
            }
        }

        let mut num_fixed: HighsInt = 0;
        for i in propagation_orbitopes {
            num_fixed += self.orbitopes[i as usize].orbital_fixing(ms, domain);
            if domain.infeasible() {
                break;
            }
        }

        num_fixed
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetries::computeStabilizerOrbits (lines 178-269)
    /// C++: `std::shared_ptr<const StabilizerOrbits> computeStabilizerOrbits(const HighsDomain & localdom, StabilizerOrbitWorkspace & workspace)`
    /// PORT NOTE: signature changed because the domain needs ms (`variableType`, `isGlobalBinary`).
    /// PORT NOTE: C++ `stabilizerOrbits.symmetries = this` has no field here (see `StabilizerOrbits`).
    /// The C++ `make_shared` never returns null, so this always returns `Some`.
    pub fn compute_stabilizer_orbits(
        &self,
        ms: &mut HighsMipSolver,
        localdom: &HighsDomain,
        workspace: &mut StabilizerOrbitWorkspace,
    ) -> Option<Arc<StabilizerOrbits>> {
        let domchg_stack = localdom.get_domain_change_stack();
        let branching_pos = localdom.get_branching_positions();

        let mut stabilizer_orbits = StabilizerOrbits::default();
        stabilizer_orbits
            .stabilized_cols
            .reserve(self.permutation_columns.len());
        for &i in branching_pos {
            let col = domchg_stack[i as usize].column;
            if self.column_position[col as usize] == -1 {
                continue;
            }

            debug_assert!(localdom.variable_type(ms, col) != HighsVarType::kContinuous);

            // if we branch a variable upwards it is either binary and branched to one
            // and needs to be stabilized or it is a general integer and needs to be
            // stabilized regardless of branching direction.
            // if we branch downwards we only need to stabilize on this
            // branching column if it is a general integer
            if domchg_stack[i as usize].boundtype == HighsBoundType::kLower
                || !localdom.is_global_binary(ms, col)
            {
                stabilizer_orbits
                    .stabilized_cols
                    .push(self.column_position[col as usize]);
            }
        }

        let perm_length = self.permutation_columns.len();
        // std::iota over the resized partition
        workspace.orbit_partition = (0..perm_length as HighsInt).collect();
        workspace.orbit_size = vec![1; perm_length];
        workspace.link_compression_stack.clear();

        for i in 0..self.num_perms as usize {
            // C++ `permutations.data() + i * permutationColumns.size()`
            let perm_base = i * perm_length;

            let mut perm_respects_branchings = true;
            for &pos in &stabilizer_orbits.stabilized_cols {
                if self.permutation_columns[pos as usize]
                    != self.permutations[perm_base + pos as usize]
                {
                    perm_respects_branchings = false;
                    break;
                }
            }

            if !perm_respects_branchings {
                continue;
            }

            for j in 0..perm_length {
                let v1 = self.permutation_columns[j];
                let v2 = self.permutations[perm_base + j];
                self.merge_orbits(
                    v1,
                    v2,
                    &mut workspace.orbit_partition,
                    &mut workspace.orbit_size,
                    &mut workspace.link_compression_stack,
                );
            }
        }

        stabilizer_orbits.stabilized_cols.clear();

        stabilizer_orbits.orbit_cols.reserve(perm_length);
        for i in 0..perm_length {
            let col = self.permutation_columns[i];
            if localdom.variable_type(ms, col) == HighsVarType::kContinuous {
                continue;
            }
            let orbit = self.get_orbit(
                col,
                &mut workspace.orbit_partition,
                &mut workspace.link_compression_stack,
            );
            if workspace.orbit_size[orbit as usize] == 1 {
                stabilizer_orbits.stabilized_cols.push(col);
            } else if localdom.is_global_binary(ms, col) {
                stabilizer_orbits.orbit_cols.push(col);
            }
        }

        pdqsort(&mut stabilizer_orbits.stabilized_cols);
        if !stabilizer_orbits.orbit_cols.is_empty() {
            pdqsort_by(
                &mut stabilizer_orbits.orbit_cols,
                |col1: &HighsInt, col2: &HighsInt| {
                    let orbit1 = self.get_orbit(
                        *col1,
                        &mut workspace.orbit_partition,
                        &mut workspace.link_compression_stack,
                    );
                    let orbit2 = self.get_orbit(
                        *col2,
                        &mut workspace.orbit_partition,
                        &mut workspace.link_compression_stack,
                    );
                    orbit1 < orbit2
                },
            );
            let num_orbit_cols = stabilizer_orbits.orbit_cols.len();
            stabilizer_orbits.orbit_starts.reserve(num_orbit_cols + 1);
            stabilizer_orbits.orbit_starts.push(0);

            for i in 1..num_orbit_cols {
                let orbit_i = self.get_orbit(
                    stabilizer_orbits.orbit_cols[i],
                    &mut workspace.orbit_partition,
                    &mut workspace.link_compression_stack,
                );
                let orbit_prev = self.get_orbit(
                    stabilizer_orbits.orbit_cols[i - 1],
                    &mut workspace.orbit_partition,
                    &mut workspace.link_compression_stack,
                );
                if orbit_i != orbit_prev {
                    stabilizer_orbits.orbit_starts.push(i as HighsInt);
                }
            }
            stabilizer_orbits
                .orbit_starts
                .push(num_orbit_cols as HighsInt);
        }

        Some(Arc::new(stabilizer_orbits))
    }
}

impl StabilizerOrbits {
    /// HiGHS: HighsSymmetry.cpp::StabilizerOrbits::orbitalFixing (lines 271-319)
    /// C++: `HighsInt orbitalFixing(HighsDomain & domain) const`
    /// PORT NOTE: signature changed because the domain needs ms (`propagate`, `changeBound`).
    /// PORT NOTE: `symmetries` is an explicit argument, standing in for the C++ back-pointer.
    pub fn orbital_fixing(
        &self,
        ms: &mut HighsMipSolver,
        symmetries: &HighsSymmetries,
        domain: &mut HighsDomain,
    ) -> HighsInt {
        let mut num_fixed = symmetries.propagate_orbitopes(ms, domain);
        if domain.infeasible() || self.orbit_cols.is_empty() {
            return num_fixed;
        }

        let num_orbits = self.orbit_starts.len() as HighsInt - 1;
        let mut i: HighsInt = 0;
        while i < num_orbits {
            let start = self.orbit_starts[i as usize] as usize;
            let end = self.orbit_starts[i as usize + 1] as usize;

            let mut fixcol: HighsInt = -1;
            for j in start..end {
                if domain.is_fixed(self.orbit_cols[j]) {
                    fixcol = self.orbit_cols[j];
                    break;
                }
            }

            if fixcol != -1 {
                let old_num_fixed = num_fixed;
                let old_size = domain.get_domain_change_stack().len();
                if domain.col_lower_[fixcol as usize] == 1.0 {
                    for j in start..end {
                        let col = self.orbit_cols[j];
                        if domain.col_lower_[col as usize] == 1.0 {
                            continue;
                        }
                        num_fixed += 1;
                        domain.change_bound_2(
                            ms,
                            HighsBoundType::kLower,
                            col,
                            1.0,
                            HighsDomain_Reason::unspecified(),
                        );
                        if domain.infeasible() {
                            return num_fixed;
                        }
                    }
                } else {
                    for j in start..end {
                        let col = self.orbit_cols[j];
                        if domain.col_upper_[col as usize] == 0.0 {
                            continue;
                        }
                        num_fixed += 1;
                        domain.change_bound_2(
                            ms,
                            HighsBoundType::kUpper,
                            col,
                            0.0,
                            HighsDomain_Reason::unspecified(),
                        );
                        if domain.infeasible() {
                            return num_fixed;
                        }
                    }
                }

                let new_fixed = num_fixed - old_num_fixed;

                if new_fixed != 0 {
                    domain.propagate(ms);
                    if domain.infeasible() {
                        return num_fixed;
                    }
                    // C++ `(HighsInt)(size() - oldSize)`: unsigned difference, then truncated
                    let grown = domain
                        .get_domain_change_stack()
                        .len()
                        .wrapping_sub(old_size) as HighsInt;
                    if grown > new_fixed {
                        i = -1;
                    }
                }
            }
            i += 1;
        }

        num_fixed
    }

    /// HiGHS: HighsSymmetry.cpp::StabilizerOrbits::isStabilized (lines 321-324)
    /// C++: `bool isStabilized(HighsInt col) const`
    pub fn is_stabilized(&self, symmetries: &HighsSymmetries, col: HighsInt) -> bool {
        symmetries.column_position[col as usize] == -1
            || self.stabilized_cols.binary_search(&col).is_ok()
    }
}

impl HighsOrbitopeMatrix {
    /// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::determineOrbitopeType (lines 326-429)
    /// C++: `void determineOrbitopeType(HighsCliqueTable & cliquetable)`
    pub fn determine_orbitope_type(&mut self, cliquetable: &mut HighsCliqueTable) {
        for j in 0..self.row_length {
            for i in 0..self.num_rows {
                let key = *self.entry_2(i, j);
                self.column_to_row.insert(HighsHashTableEntry::new(key, i));
            }
        }

        // rowIsSetPacking.assign(numRows, kRowUndetermined)
        self.row_is_set_packing = vec![Self::kRowUndetermined; self.num_rows as usize];
        self.num_set_packing_rows = 0;

        for j in 1..self.row_length {
            // HighsInt* colj1 = &entry(0, j): column j starts at offset j * numRows
            let colj1 = j as usize * self.num_rows as usize;

            for j0 in 0..j {
                let colj0 = j0 as usize * self.num_rows as usize;

                for i in 0..self.num_rows {
                    if self.row_is_set_packing[i as usize] != Self::kRowUndetermined {
                        continue;
                    }

                    let xj0 = self.matrix[colj0 + i as usize];
                    let xj1 = self.matrix[colj1 + i as usize];

                    let (common_clique, count) = cliquetable.find_common_clique(
                        HighsCliqueTable_CliqueVar::from_parts(xj0 as HighsUInt, 1),
                        HighsCliqueTable_CliqueVar::from_parts(xj1 as HighsUInt, 1),
                    );

                    // PORT NOTE: C++ `first == nullptr` is an empty slice here; a found clique
                    // always has at least two entries, so empty <=> nullptr.
                    if common_clique.is_empty() {
                        self.row_is_set_packing[i as usize] = Self::kRowNotPacking;
                        continue;
                    }

                    let mut overlap: HighsInt = 0;

                    for k in 0..count as usize {
                        if common_clique[k].val() == 0 {
                            continue;
                        }

                        if let Some(&clique_col_row) = self
                            .column_to_row
                            .find(&(common_clique[k].col() as HighsInt))
                        {
                            if clique_col_row == i {
                                overlap += 1;
                            }
                        }
                    }

                    if overlap == self.row_length {
                        self.row_is_set_packing[i as usize] = Self::kRowPacking;
                        self.num_set_packing_rows += 1;
                        if self.num_set_packing_rows == self.num_rows {
                            break;
                        }
                    }
                }

                if self.num_set_packing_rows == self.num_rows {
                    break;
                }
            }

            if self.num_set_packing_rows == self.num_rows {
                break;
            }
        }

        // now for the rows that do not have a set packing structure check
        // if we have such structure when all columns in the row are negated

        for i in 0..self.num_rows as usize {
            if self.row_is_set_packing[i] == Self::kRowNotPacking {
                self.row_is_set_packing[i] = Self::kRowUndetermined;
            }
        }

        for j in 1..self.row_length {
            let colj1 = j as usize * self.num_rows as usize;

            for j0 in 0..j {
                let colj0 = j0 as usize * self.num_rows as usize;

                for i in 0..self.num_rows {
                    if self.row_is_set_packing[i as usize] != Self::kRowUndetermined {
                        continue;
                    }

                    let xj0 = self.matrix[colj0 + i as usize];
                    let xj1 = self.matrix[colj1 + i as usize];

                    // now look for cliques with value 0
                    let (common_clique, count) = cliquetable.find_common_clique(
                        HighsCliqueTable_CliqueVar::from_parts(xj0 as HighsUInt, 0),
                        HighsCliqueTable_CliqueVar::from_parts(xj1 as HighsUInt, 0),
                    );

                    // PORT NOTE: see the nullptr note above.
                    if common_clique.is_empty() {
                        self.row_is_set_packing[i as usize] = Self::kRowNotPacking;
                        continue;
                    }

                    let mut overlap: HighsInt = 0;

                    for k in 0..count as usize {
                        // skip clique variables with values of 1
                        if common_clique[k].val() == 1 {
                            continue;
                        }

                        if let Some(&clique_col_row) = self
                            .column_to_row
                            .find(&(common_clique[k].col() as HighsInt))
                        {
                            if clique_col_row == i {
                                overlap += 1;
                            }
                        }
                    }

                    if overlap == self.row_length {
                        // mark with value 2, for negated set packing row with at most one zero
                        self.row_is_set_packing[i as usize] = Self::kRowPackingNegated;
                        self.num_set_packing_rows += 1;
                        if self.num_set_packing_rows == self.num_rows {
                            break;
                        }
                    }
                }

                if self.num_set_packing_rows == self.num_rows {
                    break;
                }
            }

            if self.num_set_packing_rows == self.num_rows {
                break;
            }
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::getBranchingColumn (lines 431-444)
    /// C++: `HighsInt getBranchingColumn(const std::vector<double> & colLower, const std::vector<double> & colUpper, HighsInt col) const`
    pub fn get_branching_column(
        &self,
        col_lower: &[f64],
        col_upper: &[f64],
        col: HighsInt,
    ) -> HighsInt {
        if let Some(&i) = self.column_to_row.find(&col) {
            if self.row_is_set_packing[i as usize] > Self::kRowNotPacking {
                for j in 0..self.row_length {
                    let branch_col = *self.entry_2(i, j);
                    if branch_col == col {
                        break;
                    }
                    if col_lower[branch_col as usize] != col_upper[branch_col as usize] {
                        return branch_col;
                    }
                }
            }
        }

        col
    }

    /// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::orbitalFixingForPackingOrbitope (lines 446-601)
    /// C++: `HighsInt orbitalFixingForPackingOrbitope(const std::vector<HighsInt> & rows, HighsDomain & domain) const`
    /// PORT NOTE: signature changed because the domain needs ms (`changeBound`, `propagate`).
    pub fn orbital_fixing_for_packing_orbitope(
        &self,
        ms: &mut HighsMipSolver,
        rows: &[HighsInt],
        domain: &mut HighsDomain,
    ) -> HighsInt {
        let num_dynamic_rows = rows.len();

        let mut first_one_in_row: Vec<HighsInt> = vec![-1; num_dynamic_rows];

        for j in 0..self.row_length {
            for i in 0..num_dynamic_rows {
                if first_one_in_row[i] != -1 {
                    continue;
                }
                let r = rows[i];
                let colrj = *self.entry_2(r, j);
                if self.row_is_set_packing[r as usize] == Self::kRowPacking {
                    if domain.col_lower_[colrj as usize] > 0.5 {
                        first_one_in_row[i] = j;
                    }
                } else {
                    debug_assert!(self.row_is_set_packing[r as usize] == Self::kRowPackingNegated);
                    if domain.col_upper_[colrj as usize] < 0.5 {
                        first_one_in_row[i] = j;
                    }
                }
            }
        }

        // we start looping over the rows and keep the index j at the column that
        // is the right most possible location for a 1 entry.
        // For the first row this position is 0.

        let mut j: HighsInt = 0;
        let mut num_fixed: HighsInt = 0;
        for i in 0..num_dynamic_rows {
            // at this position we know that the last possible position
            // of a 1-entry in this row is j. If we have a 1 behind j
            // the face of the orbitope intersecting the set of initial fixings is empty
            if first_one_in_row[i] > j {
                domain.mark_infeasible(HighsDomain_Reason::unspecified());
                return num_fixed;
            }
            let col_ij = *self.entry_2(rows[i], j);

            let negate_i = self.row_is_set_packing[rows[i] as usize] == Self::kRowPackingNegated;

            let not_zero_fixed = if negate_i {
                domain.col_lower_[col_ij as usize] < 0.5
            } else {
                domain.col_upper_[col_ij as usize] > 0.5
            };

            // as long as the entry is fixed to zero
            // the frontier stays at the same column
            // if we encounter an entry that is not fixed to zero
            // we need to proceed with the next column and found a frontier step
            if not_zero_fixed {
                // found a frontier step. Now we first check for the current column
                // if it can be fixed to 1.
                // For this we check if we find an infeasibility in the case where this
                // column was fixed to zero in which case the frontier would have stayed
                // at j for the following rows.
                let mut j0 = j;
                for k in (i + 1)..num_dynamic_rows {
                    if first_one_in_row[k] > j0 {
                        if negate_i {
                            domain.change_bound_2(
                                ms,
                                HighsBoundType::kUpper,
                                col_ij,
                                0.0,
                                HighsDomain_Reason::unspecified(),
                            );
                        } else {
                            domain.change_bound_2(
                                ms,
                                HighsBoundType::kLower,
                                col_ij,
                                1.0,
                                HighsDomain_Reason::unspecified(),
                            );
                        }

                        num_fixed += 1;
                        if domain.infeasible() {
                            return num_fixed;
                        }
                        break;
                    }

                    let negate_k =
                        self.row_is_set_packing[rows[k] as usize] == Self::kRowPackingNegated;
                    let not_zero_fixed_k = if negate_k {
                        domain.col_lower_[*self.entry_2(rows[k], j0) as usize] < 0.5
                    } else {
                        domain.col_upper_[*self.entry_2(rows[k], j0) as usize] > 0.5
                    };

                    if not_zero_fixed_k {
                        j0 += 1;
                        if j0 == self.row_length {
                            break;
                        }
                    }
                }

                j += 1;
                if j == self.row_length {
                    break;
                }

                for k in 0..=i {
                    // we should have checked this row at the frontier position
                    debug_assert!(first_one_in_row[k] < j);

                    let col_kj = *self.entry_2(rows[k], j);
                    let negate_k =
                        self.row_is_set_packing[rows[k] as usize] == Self::kRowPackingNegated;

                    if negate_k {
                        if domain.col_lower_[col_kj as usize] > 0.5 {
                            continue;
                        }

                        domain.change_bound_2(
                            ms,
                            HighsBoundType::kLower,
                            col_kj,
                            1.0,
                            HighsDomain_Reason::unspecified(),
                        );
                    } else {
                        if domain.col_upper_[col_kj as usize] < 0.5 {
                            continue;
                        }

                        domain.change_bound_2(
                            ms,
                            HighsBoundType::kUpper,
                            col_kj,
                            0.0,
                            HighsDomain_Reason::unspecified(),
                        );
                    }

                    num_fixed += 1;
                    if domain.infeasible() {
                        // this can happen due to deductions from earlier fixings
                        // otherwise it would have been caught by the infeasibility
                        // check within the next loop that goes over i
                        return num_fixed;
                    }
                }
            }
        }

        // check if there are more columns that can be completely fixed to zero
        j += 1;
        while j < self.row_length {
            for k in 0..num_dynamic_rows {
                // we should have checked this row at the frontier position
                debug_assert!(first_one_in_row[k] < j);

                let col_kj = *self.entry_2(rows[k], j);
                let negate_k =
                    self.row_is_set_packing[rows[k] as usize] == Self::kRowPackingNegated;

                if negate_k {
                    if domain.col_lower_[col_kj as usize] > 0.5 {
                        continue;
                    }

                    domain.change_bound_2(
                        ms,
                        HighsBoundType::kLower,
                        col_kj,
                        1.0,
                        HighsDomain_Reason::unspecified(),
                    );
                } else {
                    if domain.col_upper_[col_kj as usize] < 0.5 {
                        continue;
                    }

                    domain.change_bound_2(
                        ms,
                        HighsBoundType::kUpper,
                        col_kj,
                        0.0,
                        HighsDomain_Reason::unspecified(),
                    );
                }
                num_fixed += 1;
                if domain.infeasible() {
                    // this can happen due to deductions from earlier fixings
                    // otherwise it would have been caught by the infeasibility
                    // check within the next loop that goes over i
                    return num_fixed;
                }
            }
            j += 1;
        }

        if !domain.infeasible() && num_fixed != 0 {
            domain.propagate(ms);
        }

        num_fixed
    }

    /// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::orbitalFixingForFullOrbitope (lines 603-728)
    /// C++: `HighsInt orbitalFixingForFullOrbitope(const std::vector<HighsInt> & rows, HighsDomain & domain) const`
    /// PORT NOTE: signature changed because the domain needs ms (`changeBound`, `propagate`).
    pub fn orbital_fixing_for_full_orbitope(
        &self,
        ms: &mut HighsMipSolver,
        rows: &[HighsInt],
        domain: &mut HighsDomain,
    ) -> HighsInt {
        let num_dynamic_rows = rows.len();
        let n = num_dynamic_rows;
        let row_length = self.row_length as usize;
        // std::vector<int8_t> Mminimal(numDynamicRows * rowLength, -1); column j at offset j * n
        let mut mminimal: Vec<i8> = vec![-1; n * row_length];

        for j in 0..self.row_length {
            for i in 0..n {
                let r = rows[i];
                let colij = *self.entry_2(r, j);
                if domain.col_lower_[colij as usize] == 1.0 {
                    mminimal[i + j as usize * n] = 1;
                } else if domain.col_upper_[colij as usize] == 0.0 {
                    mminimal[i + j as usize * n] = 0;
                }
            }
        }

        let mut mmaximal = mminimal.clone();

        debug_assert!(row_length >= 1);
        let last_off = n * (row_length - 1);
        for k in 0..n {
            if mminimal[last_off + k] == -1 {
                mminimal[last_off + k] = 0;
            }
            if mmaximal[k] == -1 {
                mmaximal[k] = 1;
            }
        }

        // minimal assignment, computed from the last column backwards
        let mut j: HighsInt = self.row_length - 2;
        while j >= 0 {
            let off0 = j as usize * n;
            let off1 = off0 + n;
            let (left, right) = mminimal.split_at_mut(off1);
            let colj0 = &mut left[off0..off1];
            let colj1 = &right[..n];
            let i_f = i_fixed(colj0, colj1);

            if i_f == n as HighsInt {
                update_free_entries(colj0, colj1, n);
            } else {
                let i_d = i_discr(colj0, colj1, i_f);
                if i_d == -1 {
                    domain.mark_infeasible(HighsDomain_Reason::unspecified());
                    return 0;
                }

                update_free_entries(colj0, colj1, i_d as usize);
                colj0[i_d as usize] = 1;
                for k in (i_d as usize + 1)..n {
                    colj0[k] += i8::from(colj0[k] == -1);
                }
            }
            j -= 1;
        }

        // maximal assignment, computed from the first column forwards
        for j in 1..self.row_length {
            let off0 = j as usize * n;
            let off1 = off0 - n;
            let (left, right) = mmaximal.split_at_mut(off0);
            let colj0 = &mut right[..n];
            let colj1 = &left[off1..off0];
            let i_f = i_fixed(colj1, colj0);

            if i_f == n as HighsInt {
                update_free_entries(colj0, colj1, n);
            } else {
                let i_d = i_discr(colj1, colj0, i_f);
                if i_d == -1 {
                    domain.mark_infeasible(HighsDomain_Reason::unspecified());
                    return 0;
                }
                update_free_entries(colj0, colj1, i_d as usize);
                colj0[i_d as usize] = 0;
                for k in (i_d as usize + 1)..n {
                    colj0[k] += 2 * i8::from(colj0[k] == -1);
                }
            }
        }

        let mut num_fixed: HighsInt = 0;

        for j in 0..self.row_length {
            let off = j as usize * n;
            let col_maximal = &mmaximal[off..off + n];
            let col_minimal = &mminimal[off..off + n];

            for i in 0..n {
                if col_minimal[i] != col_maximal[i] {
                    debug_assert!(col_minimal[i] < col_maximal[i]);
                    break;
                }

                let r = rows[i];
                let colrj = *self.entry_2(r, j);
                if domain.is_fixed(colrj) {
                    continue;
                }

                num_fixed += 1;
                if col_minimal[i] == 1 {
                    domain.change_bound_2(
                        ms,
                        HighsBoundType::kLower,
                        colrj,
                        1.0,
                        HighsDomain_Reason::unspecified(),
                    );
                } else {
                    domain.change_bound_2(
                        ms,
                        HighsBoundType::kUpper,
                        colrj,
                        0.0,
                        HighsDomain_Reason::unspecified(),
                    );
                }
                if domain.infeasible() {
                    break;
                }
            }
            if domain.infeasible() {
                break;
            }
        }

        if !domain.infeasible() {
            domain.propagate(ms);
        }

        num_fixed
    }

    /// HiGHS: HighsSymmetry.cpp::HighsOrbitopeMatrix::orbitalFixing (lines 730-754)
    /// C++: `HighsInt orbitalFixing(HighsDomain & domain) const`
    /// PORT NOTE: signature changed because the domain needs ms (dispatches to the two routines above).
    pub fn orbital_fixing(&self, ms: &mut HighsMipSolver, domain: &mut HighsDomain) -> HighsInt {
        let mut rows: Vec<HighsInt> = Vec::with_capacity(self.num_rows as usize);
        let mut row_used: Vec<bool> = vec![false; self.num_rows as usize];

        let domchgstack = domain.get_domain_change_stack();
        let branchpos = domain.get_branching_positions();

        let mut is_packing = true;
        for &pos in branchpos {
            if let Some(&i) = self.column_to_row.find(&domchgstack[pos as usize].column) {
                if !row_used[i as usize] {
                    row_used[i as usize] = true;
                    is_packing =
                        is_packing && self.row_is_set_packing[i as usize] > Self::kRowNotPacking;
                    rows.push(i);
                }
            }
        }

        if rows.is_empty() {
            return 0;
        }

        if is_packing {
            return self.orbital_fixing_for_packing_orbitope(ms, &rows, domain);
        }

        self.orbital_fixing_for_full_orbitope(ms, &rows, domain)
    }
}

#[cfg(test)]
mod h16_sym_refine_tests {
    use super::super::HighsSymmetryDetection;
    use super::*;

    /// One cell {0,1,2,3}; vertex 0 has no hash, vertices 1,2,3 hash to 1,2,3.
    /// Entry stacks: cell_creation_stack has one entry (9), curr_node_certificate is
    /// empty. The first split (at position 1) succeeds and pushes to both stacks. The
    /// second split (at position 2, certificate prefix length 1) fails the
    /// best-leave comparison, so partition_refinement prunes and returns false.
    fn prune_setup(sym_refine: bool) -> HighsSymmetryDetection {
        let mut sd = HighsSymmetryDetection::default();
        sd.h16_sym_refine = sym_refine;
        sd.gstart = vec![0; 4];
        sd.gend = vec![0; 4];
        sd.current_partition = vec![0, 1, 2, 3];
        sd.current_partition_links = vec![4, 0, 0, 0];
        sd.vertex_to_cell = vec![0, 0, 0, 0];
        sd.cell_in_refinement_queue = vec![1, 0, 0, 0];
        sd.refinement_queue = vec![0];
        *sd.vertex_hash.index_or_insert(&1) = 1;
        *sd.vertex_hash.index_or_insert(&2) = 2;
        *sd.vertex_hash.index_or_insert(&3) = 3;
        sd.cell_creation_stack = vec![9];
        sd.curr_node_certificate = Vec::new();
        sd.first_leave_certificate = vec![u32::MAX; 4];
        sd.best_leave_certificate = vec![u32::MAX, 0, 0, 0];
        sd.first_leave_prefix_len = 1;
        sd.best_leave_prefix_len = 1;
        sd
    }

    #[test]
    fn h16_sym_refine_on_prune_restores_stacks() {
        let mut sd = prune_setup(true);
        assert!(!sd.partition_refinement());
        // stacks restored to their entry sizes (1 and 0)
        assert_eq!(sd.cell_creation_stack, vec![9]);
        assert!(sd.curr_node_certificate.is_empty());
        // prefix lengths clamped to the entry certificate size (0)
        assert_eq!(sd.first_leave_prefix_len, 0);
        assert_eq!(sd.best_leave_prefix_len, 0);
    }

    #[test]
    fn h16_sym_refine_off_prune_keeps_stacks() {
        let mut sd = prune_setup(false);
        assert!(!sd.partition_refinement());
        // 1.15.1 behaviour: the successful first split is not undone
        assert_eq!(sd.cell_creation_stack.len(), 2);
        assert_eq!(sd.curr_node_certificate.len(), 1);
        assert_eq!(sd.first_leave_prefix_len, 1);
        assert_eq!(sd.best_leave_prefix_len, 1);
    }
}

#[cfg(test)]
mod x1b_sym_hash64_tests {
    //! HiGHS #3535 (mip_sym_hash64): vertex hashes are `sparse_combine32` values with the
    //! option off (1.15.1) and 61-bit `sparse_combine` values with it on (1.16.0).
    use super::super::HighsSymmetryDetection;
    use super::*;

    fn combined(hash64: bool) -> u64 {
        let mut sd = HighsSymmetryDetection::default();
        sd.mip_sym_hash64 = hash64;
        sd.combine_vertex_hash(7, 70, 3);
        sd.combine_vertex_hash(7, 130, 5);
        sd.get_vertex_hash(7)
    }

    #[test]
    fn x1b_sym_hash64_off_matches_sparse_combine32() {
        let mut h: u32 = 0;
        sparse_combine32(&mut h, 70, 3);
        sparse_combine32(&mut h, 130, 5);
        assert_eq!(combined(false), u64::from(h));
    }

    #[test]
    fn x1b_sym_hash64_on_matches_sparse_combine() {
        let mut h: u64 = 0;
        sparse_combine_value(&mut h, 70, 3);
        sparse_combine_value(&mut h, 130, 5);
        assert_eq!(combined(true), h);
        assert_ne!(combined(true), combined(false));
    }
}
