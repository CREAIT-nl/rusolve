//! HiGHS source: presolve/HighsSymmetry.cpp (lines 1001-1937)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

impl HighsSymmetryDetection {
    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::selectTargetCell (lines 1002-1013)
    /// C++: `HighsInt selectTargetCell() const`
    pub fn select_target_cell(&self) -> HighsInt {
        let mut i: HighsInt = 0;
        if self.node_stack.len() > 1 {
            i = self.node_stack[self.node_stack.len() - 2].target_cell;
        }

        while i < self.num_vertices {
            if self.cell_size(i) > 1 {
                return i;
            }

            i += 1;
        }

        -1
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::checkStoredAutomorphism (lines 1015-1036)
    /// C++: `bool checkStoredAutomorphism(HighsInt vertex) const`
    pub fn check_stored_automorphism(&self, vertex: HighsInt) -> bool {
        let num_check = cpp_min(self.num_automorphisms, 64);

        for i in 0..num_check {
            // C++: const HighsInt* automorphism = automorphisms.data() + i * numVertices;
            let automorphism = &self.automorphisms[(i * self.num_vertices) as usize..];
            let mut automorphism_useful = true;
            // PORT NOTE: C++ `nodeStack.size() - 2` is size_t; for sizes below 2 the
            // conversion to HighsInt yields the same negative value as this signed form.
            let mut j: HighsInt = self.node_stack.len() as HighsInt - 2;
            while j >= self.first_path_depth {
                let fix_pos =
                    self.vertex_position[self.node_stack[j as usize].last_distiguished as usize];

                if automorphism[fix_pos as usize] != self.vertex_ground_set[fix_pos as usize] {
                    automorphism_useful = false;
                    break;
                }
                j -= 1;
            }

            if !automorphism_useful {
                continue;
            }

            if automorphism[self.vertex_position[vertex as usize] as usize] < vertex {
                return false;
            }
        }

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::determineNextToDistinguish (lines 1038-1075)
    /// C++: `bool determineNextToDistinguish()`
    pub fn determine_next_to_distinguish(&mut self) -> bool {
        // PORT NOTE: C++ `Node& currNode = nodeStack.back()` is only read here, so a copy is used.
        let curr_node = self.node_stack[self.node_stack.len() - 1];
        // PORT NOTE: C++ `distinguishCands` holds `HighsInt*` into `currentPartition`; the Rust
        // field holds the positions in `currentPartition` instead.
        self.distinguish_cands.clear();
        let cell_start = curr_node.target_cell as usize;
        let cell_end = self.current_partition_links[curr_node.target_cell as usize] as usize;

        if curr_node.last_distiguished == -1 {
            // C++ std::min_element(cellStart, cellEnd): first strictly smallest element.
            let mut next_distinguish_pos = cell_start;
            for i in cell_start..cell_end {
                if self.current_partition[i] < self.current_partition[next_distinguish_pos] {
                    next_distinguish_pos = i;
                }
            }
            self.distinguish_cands.push(next_distinguish_pos);
        } else if self.node_stack.len() as HighsInt > self.first_path_depth {
            for i in cell_start..cell_end {
                let v = self.current_partition[i];
                if v > curr_node.last_distiguished && self.check_stored_automorphism(v) {
                    self.distinguish_cands.push(i);
                }
            }
            if self.distinguish_cands.is_empty() {
                return false;
            }
            // C++ std::min_element over the candidate pointers, then swap to the front and resize(1).
            let mut next_distinguish_idx = 0;
            for k in 0..self.distinguish_cands.len() {
                if self.current_partition[self.distinguish_cands[k]]
                    < self.current_partition[self.distinguish_cands[next_distinguish_idx]]
                {
                    next_distinguish_idx = k;
                }
            }
            self.distinguish_cands.swap(0, next_distinguish_idx);
            self.distinguish_cands.truncate(1);
        } else {
            for i in cell_start..cell_end {
                let v = self.current_partition[i];
                if v > curr_node.last_distiguished {
                    let orbit = self.get_orbit(v);
                    if self.vertex_ground_set[orbit as usize] == v {
                        self.distinguish_cands.push(i);
                    }
                }
            }
            if self.distinguish_cands.is_empty() {
                return false;
            }
            // C++ std::min_element over the candidate pointers, then swap to the front and resize(1).
            let mut next_distinguish_idx = 0;
            for k in 0..self.distinguish_cands.len() {
                if self.current_partition[self.distinguish_cands[k]]
                    < self.current_partition[self.distinguish_cands[next_distinguish_idx]]
                {
                    next_distinguish_idx = k;
                }
            }
            self.distinguish_cands.swap(0, next_distinguish_idx);
            self.distinguish_cands.truncate(1);
        }

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::distinguishVertex (lines 1077-1089)
    /// C++: `bool distinguishVertex(HighsInt targetCell)`
    pub fn distinguish_vertex(&mut self, target_cell: HighsInt) -> bool {
        debug_assert_eq!(self.distinguish_cands.len(), 1);
        let target_cell_end = self.current_partition_links[target_cell as usize];
        let new_cell = target_cell_end - 1;
        // C++ std::swap(*distinguishCands[0], currentPartition[newCell])
        let cand = self.distinguish_cands[0];
        self.current_partition.swap(cand, new_cell as usize);
        let distinguished = self.current_partition[new_cell as usize];
        let node_len = self.node_stack.len();
        self.node_stack[node_len - 1].last_distiguished = distinguished;

        if !self.split_cell(target_cell, new_cell) {
            return false;
        }

        self.update_cell_membership(new_cell, new_cell, true);

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::backtrack (lines 1091-1110)
    /// C++: `void backtrack(HighsInt backtrackStackNewEnd, HighsInt backtrackStackEnd)`
    pub fn backtrack(&mut self, backtrack_stack_new_end: HighsInt, backtrack_stack_end: HighsInt) {
        // we assume that we always backtrack from a leaf node, i.e. a discrete
        // partition therefore we do not need to remember the values of the hash
        // contributions as it is the indentity for each position and all new cells
        // are on the cell creation stack.
        let mut stack_pos = backtrack_stack_end - 1;
        while stack_pos >= backtrack_stack_new_end {
            let cell = self.cell_creation_stack[stack_pos as usize];
            // look up the cell start of the preceding cell with link compression
            let new_start = self.get_cell_start(cell - 1);
            // remember the current end
            let curr_end = self.current_partition_links[cell as usize];
            // change the link to point to the start of the preceding cell
            self.current_partition_links[cell as usize] = new_start;
            // change the link of the start pointer of the preceding cell to point to
            // the end of this cell
            self.current_partition_links[new_start as usize] = curr_end;
            stack_pos -= 1;
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::cleanupBacktrack (lines 1112-1129)
    /// C++: `void cleanupBacktrack(HighsInt cellCreationStackPos)`
    pub fn cleanup_backtrack(&mut self, cell_creation_stack_pos: HighsInt) {
        // the links have been updated. Even though they might still not be fully
        // compressed the cell starts will all point to the correct cell end and the
        // lookup with path compression will give the correct start
        let mut stack_pos = self.cell_creation_stack.len() as HighsInt - 1;
        while stack_pos >= cell_creation_stack_pos {
            let cell = self.cell_creation_stack[stack_pos as usize];

            let cell_start = self.get_cell_start(cell);
            let cell_end = self.current_partition_links[cell_start as usize];

            let mut v = cell;
            while v < cell_end
                && self.vertex_to_cell[self.current_partition[v as usize] as usize] == cell
            {
                self.update_cell_membership(v, cell_start, false);
                v += 1;
            }
            stack_pos -= 1;
        }

        self.cell_creation_stack
            .resize(cell_creation_stack_pos as usize, 0);
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::getCellStart (lines 1131-1148)
    /// C++: `HighsInt getCellStart(HighsInt pos)`
    pub fn get_cell_start(&mut self, pos: HighsInt) -> HighsInt {
        let mut pos = pos;
        let mut start_pos = self.current_partition_links[pos as usize];
        if start_pos > pos {
            return pos;
        }
        if self.current_partition_links[start_pos as usize] < start_pos {
            loop {
                self.link_compression_stack.push(pos);
                pos = start_pos;
                start_pos = self.current_partition_links[start_pos as usize];
                if self.current_partition_links[start_pos as usize] >= start_pos {
                    break;
                }
            }

            loop {
                let top = self.link_compression_stack[self.link_compression_stack.len() - 1];
                self.current_partition_links[top as usize] = start_pos;
                self.link_compression_stack.pop();
                if self.link_compression_stack.is_empty() {
                    break;
                }
            }
        }

        start_pos
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::createNode (lines 1150-1156)
    /// C++: `void createNode()`
    pub fn create_node(&mut self) {
        self.node_stack.push(HighsSymmetryDetection_Node {
            stack_start: self.cell_creation_stack.len() as HighsInt,
            certificate_end: self.curr_node_certificate.len() as HighsInt,
            target_cell: -1,
            last_distiguished: -1,
        });
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::loadModelAsGraph (lines 1180-1318)
    /// C++: `void loadModelAsGraph(const HighsLp & model, double epsilon)`
    pub fn load_model_as_graph(&mut self, model: &HighsLp, epsilon: f64) {
        // PORT NOTE: C++ `this->model = &model` is dropped; the model is passed explicitly.
        self.num_col = model.num_col_;
        self.num_row = model.num_row_;
        self.num_vertices = self.num_col + self.num_row;
        let num_vertices = self.num_vertices as usize;

        self.cell_in_refinement_queue.resize(num_vertices, 0);
        self.vertex_to_cell.resize(num_vertices, 0);
        self.refinement_queue.reserve(num_vertices);
        self.curr_node_certificate.reserve(num_vertices);

        let mut column_set: HighsHashTable<MatrixColumn, HighsInt> = HighsHashTable::new();
        let mut row_set: HighsHashTable<MatrixRow, HighsInt> = HighsHashTable::new();
        let mut coloring = HighsMatrixColoring::new(epsilon);
        self.edge_buffer.resize(num_vertices, (0, 0));
        // set up row and column based incidence matrix
        let num_nz = model.a_matrix_.index_.len() as HighsInt;
        self.gedge.resize(2 * num_nz as usize, (0, 0));
        for (k, &row_index) in model.a_matrix_.index_.iter().enumerate() {
            self.gedge[k] = (row_index + self.num_col, 0);
        }

        self.gstart.resize(num_vertices + 1, 0);
        let start_len = model.a_matrix_.start_.len();
        self.gstart[..start_len].copy_from_slice(&model.a_matrix_.start_);

        // set up the column colors and count row sizes
        let mut row_sizes: Vec<HighsInt> = vec![0; self.num_row as usize];
        for i in 0..self.num_col {
            for j in self.gstart[i as usize]..self.gstart[(i + 1) as usize] {
                self.gedge[j as usize].1 =
                    coloring.color(model.a_matrix_.value_[j as usize]) as HighsUInt;
                row_sizes[model.a_matrix_.index_[j as usize] as usize] += 1;
            }
        }

        // next set up the row starts using the computed row sizes
        let mut offset = num_nz;
        for i in 0..self.num_row {
            self.gstart[(self.num_col + i) as usize] = offset;
            offset += row_sizes[i as usize];
        }
        self.gstart[(self.num_col + self.num_row) as usize] = offset;

        self.gend = self.gstart[1..].to_vec();

        // finally add the nonzeros to the row major matrix
        for i in 0..self.num_col {
            for j in self.gstart[i as usize]..self.gstart[(i + 1) as usize] {
                let row = model.a_matrix_.index_[j as usize];
                let ar_pos =
                    self.gstart[(self.num_col + row + 1) as usize] - row_sizes[row as usize];
                row_sizes[row as usize] -= 1;
                self.gedge[ar_pos as usize].0 = i;
                self.gedge[ar_pos as usize].1 = self.gedge[j as usize].1;
            }
        }

        // loop over the columns and assign them a number that is distinct based on
        // their upper/lower bounds, cost, and integrality status. Use the columnSet
        // hash table to look up whether a column with similar properties exists and
        // use the previous number in that case. The number is stored in the
        // colToCell array which is subsequently used to sort an initial column
        // permutation.
        let mut index_offset = self.num_col + 1;
        for i in 0..self.num_col {
            let iu = i as usize;
            let mut matrix_col = MatrixColumn::default();

            matrix_col.cost = coloring.color(model.col_cost_[iu]);
            matrix_col.lb = coloring.color(model.col_lower_[iu]);
            matrix_col.ub = coloring.color(model.col_upper_[iu]);
            matrix_col.integral = model.integrality_[iu] as u32;
            matrix_col.len = (self.gstart[(i + 1) as usize] - self.gstart[iu]) as u32;

            // PORT NOTE: C++ holds `HighsInt* columnCell = &columnSet[matrixCol]` across
            // `columnSet.size()`. Rust reads the value, then writes it back by key lookup.
            // The key already exists at that point, so no insertion or rehash happens.
            let mut column_cell = *column_set.index_or_insert(&matrix_col);
            if column_cell == 0 {
                column_cell = column_set.size() as HighsInt;

                if model.col_lower_[iu] != 0.0
                    || model.col_upper_[iu] != 1.0
                    || model.integrality_[iu] == HighsVarType::kContinuous
                {
                    column_cell += index_offset;
                }
                *column_set.index_or_insert(&matrix_col) = column_cell;
            }

            self.vertex_to_cell[iu] = column_cell;
        }

        index_offset = 2 * self.num_col + 1;
        for i in 0..self.num_row {
            let iu = i as usize;
            let mut matrix_row = MatrixRow::default();

            matrix_row.lb = coloring.color(model.row_lower_[iu]);
            matrix_row.ub = coloring.color(model.row_upper_[iu]);
            matrix_row.len = (self.gstart[(self.num_col + i + 1) as usize]
                - self.gstart[(self.num_col + i) as usize]) as u32;

            // PORT NOTE: same key write-back as the column loop above.
            let mut row_cell = *row_set.index_or_insert(&matrix_row);
            if row_cell == 0 {
                row_cell = row_set.size() as HighsInt;
                *row_set.index_or_insert(&matrix_row) = row_cell;
            }

            self.vertex_to_cell[(self.num_col + i) as usize] = index_offset + row_cell;
        }

        // set up the initial partition array, sort by the colToCell value
        // assigned above
        self.current_partition.resize(num_vertices, 0);
        for (k, v) in self.current_partition.iter_mut().enumerate() {
            // C++ std::iota(begin, end, 0)
            *v = k as HighsInt;
        }
        pdqsort_by(&mut self.current_partition, |v1, v2| {
            self.vertex_to_cell[*v1 as usize] < self.vertex_to_cell[*v2 as usize]
        });

        // now set up partition links and correct the colToCell array to the
        // correct cell index
        self.current_partition_links.resize(num_vertices, 0);
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

            // correct the colToCell array to not store the start index of the
            // cell, not its number
            self.vertex_to_cell[vertex as usize] = cell_start;
            // set the link of the column to the cellStart
            self.current_partition_links[i as usize] = cell_start;
        }

        // set the column partition link of the last started cell to point past the
        // end
        self.current_partition_links[cell_start as usize] = self.num_vertices;
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::dumpCurrentGraph (lines 1320-1334)
    /// C++: `HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> dumpCurrentGraph()`
    pub fn dump_current_graph(&mut self) -> HighsHashTable<(HighsInt, HighsInt, HighsUInt), ()> {
        let mut graph_triplets: HighsHashTable<(HighsInt, HighsInt, HighsUInt), ()> =
            HighsHashTable::new();

        for i in 0..self.num_col {
            let iu = i as usize;
            let col_cell = self.vertex_to_cell[iu];
            let mut j = self.gstart[iu];
            while j != self.gend[iu] {
                let edge = self.gedge[j as usize];
                graph_triplets.insert(HighsHashTableEntry::new(
                    (self.vertex_to_cell[edge.0 as usize], col_cell, edge.1),
                    (),
                ));
                j += 1;
            }
            let mut j = self.gend[iu];
            while j != self.gstart[(i + 1) as usize] {
                let edge = self.gedge[j as usize];
                graph_triplets.insert(HighsHashTableEntry::new((edge.0, col_cell, edge.1), ()));
                j += 1;
            }
        }

        graph_triplets
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::switchToNextNode (lines 1336-1377)
    /// C++: `void switchToNextNode(HighsInt backtrackDepth)`
    pub fn switch_to_next_node(&mut self, backtrack_depth: HighsInt) {
        let mut stack_end = self.cell_creation_stack.len() as HighsInt;
        // we need to backtrack the datastructures
        self.node_stack.resize(
            backtrack_depth as usize,
            HighsSymmetryDetection_Node::default(),
        );
        if backtrack_depth == 0 {
            return;
        }
        // PORT NOTE: C++ `do { ... continue; ... } while (!nodeStack.empty())`. The `continue`
        // statements become `break 'body`, which falls through to the while-condition check.
        'do_loop: loop {
            'body: {
                // PORT NOTE: C++ `Node& currNode` is only read, so a copy is used.
                let curr_node = self.node_stack[self.node_stack.len() - 1];
                self.backtrack(curr_node.stack_start, stack_end);
                stack_end = curr_node.stack_start;
                self.first_path_depth =
                    cpp_min(self.node_stack.len() as HighsInt, self.first_path_depth);
                self.best_path_depth =
                    cpp_min(self.node_stack.len() as HighsInt, self.best_path_depth);
                self.first_leave_prefix_len =
                    cpp_min(curr_node.certificate_end, self.first_leave_prefix_len);
                self.best_leave_prefix_len =
                    cpp_min(curr_node.certificate_end, self.best_leave_prefix_len);
                self.curr_node_certificate
                    .resize(curr_node.certificate_end as usize, 0);
                if !self.determine_next_to_distinguish() {
                    self.node_stack.pop();
                    break 'body;
                }

                // call cleanup backtrack with the final stackEnd
                // so that all hashes are up to date and the link arrays do not contain
                // chains anymore
                self.cleanup_backtrack(stack_end);
                let target_cell = curr_node.target_cell;

                if !self.distinguish_vertex(target_cell) {
                    // if distinguishing the next vertex fails, it means that its certificate
                    // value is lexicographically larger than that of the best leave
                    self.node_stack.pop();
                    break 'body;
                }

                if !self.partition_refinement() {
                    stack_end = self.cell_creation_stack.len() as HighsInt;
                    break 'body;
                }

                self.create_node();
                break 'do_loop;
            }
            if self.node_stack.is_empty() {
                break 'do_loop;
            }
        }
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::compareCurrentGraph (lines 1379-1407)
    /// C++: `bool compareCurrentGraph(const HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> & otherGraph, HighsInt & wrongCell) const`
    pub fn compare_current_graph(
        &self,
        other_graph: &HighsHashTable<(HighsInt, HighsInt, HighsUInt), ()>,
        wrong_cell: &mut HighsInt,
    ) -> bool {
        for i in 0..self.num_col {
            let iu = i as usize;
            let col_cell = self.vertex_to_cell[iu];

            let mut j = self.gstart[iu];
            while j != self.gend[iu] {
                let edge = self.gedge[j as usize];
                if other_graph
                    .find(&(self.vertex_to_cell[edge.0 as usize], col_cell, edge.1))
                    .is_none()
                {
                    // return which cell does not match in its neighbourhood as this should
                    // have been detected with the hashing it can very rarely happen due to
                    // a hash collision. In such a case we want to backtrack to the last
                    // time where we targeted this particular cell. Otherwise we could spent
                    // a long time searching for a matching leave value until every
                    // combination is exhausted and for each leave in this subtree the graph
                    // comparison will fail on this edge.
                    *wrong_cell = col_cell;
                    return false;
                }
                j += 1;
            }
            let mut j = self.gend[iu];
            while j != self.gstart[(i + 1) as usize] {
                let edge = self.gedge[j as usize];
                if other_graph.find(&(edge.0, col_cell, edge.1)).is_none() {
                    *wrong_cell = col_cell;
                    return false;
                }
                j += 1;
            }
        }

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::isFromBinaryColumn (lines 1409-1419)
    /// C++: `bool isFromBinaryColumn(HighsInt pos) const`
    pub fn is_from_binary_column(&self, model: &HighsLp, pos: HighsInt) -> bool {
        if pos >= self.num_active_cols {
            return false;
        }

        let col = self.current_partition[pos as usize] as usize;

        if model.col_lower_[col] != 0.0
            || model.col_upper_[col] != 1.0
            || model.integrality_[col] == HighsVarType::kContinuous
        {
            return false;
        }

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::computeComponentData (lines 1421-1519)
    /// C++: `ComponentData computeComponentData(const HighsSymmetries & symmetries)`
    pub fn compute_component_data(
        &mut self,
        symmetries: &HighsSymmetries,
    ) -> HighsSymmetryDetection_ComponentData {
        let mut component_data = HighsSymmetryDetection_ComponentData::default();

        component_data.components.reset(self.num_active_cols);
        component_data.first_unfixed = vec![-1; symmetries.num_perms as usize];
        component_data.num_unfixed = vec![0; symmetries.num_perms as usize];
        for i in 0..symmetries.num_perms {
            // C++: const HighsInt* perm = symmetries.permutations.data() + i * numActiveCols;
            let perm = (i * self.num_active_cols) as usize;

            for j in 0..self.num_active_cols {
                if symmetries.permutations[perm + j as usize] != self.vertex_ground_set[j as usize]
                {
                    let pos =
                        self.vertex_position[symmetries.permutations[perm + j as usize] as usize];
                    component_data.num_unfixed[i as usize] += 1;
                    if component_data.first_unfixed[i as usize] != -1 {
                        component_data
                            .components
                            .merge(component_data.first_unfixed[i as usize], pos);
                    } else {
                        component_data.first_unfixed[i as usize] = pos;
                    }
                }
            }
        }

        component_data.component_sets =
            self.vertex_ground_set[..self.num_active_cols as usize].to_vec();
        pdqsort_by(&mut component_data.component_sets, |u, v| {
            let u_comp = component_data
                .components
                .get_set(self.vertex_position[*u as usize]);
            let v_comp = component_data
                .components
                .get_set(self.vertex_position[*v as usize]);
            (component_data.components.get_set_size(u_comp) == 1, u_comp)
                < (component_data.components.get_set_size(v_comp) == 1, v_comp)
        });

        // PORT NOTE: C++ `currentComponentStart = -1` is never read before it is set in the loop,
        // so the Rust variable is declared without an initial value.
        let mut current_component_start: HighsInt;
        let mut current_component: HighsInt = -1;
        let mut curr_component_orbits: HighsHashTable<HighsInt, ()> = HighsHashTable::new();
        for i in 0..self.num_active_cols {
            let comp = component_data
                .components
                .get_set(self.vertex_position[component_data.component_sets[i as usize] as usize]);
            if component_data.components.get_set_size(comp) == 1 {
                break;
            }
            if comp != current_component {
                current_component = comp;
                current_component_start = i;
                component_data
                    .component_starts
                    .push(current_component_start);
                component_data.component_number.push(current_component);
                component_data.component_num_orbits.push(0);
                curr_component_orbits.clear();
            }

            let orbit = self.get_orbit(component_data.component_sets[i as usize]);
            if curr_component_orbits.insert(HighsHashTableEntry::new(orbit, ())) {
                *component_data.component_num_orbits.last_mut().unwrap() += 1;
            }
        }

        for i in 0..symmetries.num_perms {
            if component_data.first_unfixed[i as usize] == -1 {
                continue;
            }
            component_data.perm_components.push(i);
        }

        pdqsort_by(&mut component_data.perm_components, |i, j| {
            let seti = component_data
                .components
                .get_set(component_data.first_unfixed[*i as usize]);
            let setj = component_data
                .components
                .get_set(component_data.first_unfixed[*j as usize]);
            (seti, component_data.num_unfixed[*i as usize])
                < (setj, component_data.num_unfixed[*j as usize])
        });

        // PORT NOTE: C++ `currentComponentStart = -1` dropped: it is overwritten in the loop below
        // before any read (see the declaration above).
        current_component = -1;

        let num_used_perms = component_data.perm_components.len() as HighsInt;

        for i in 0..num_used_perms {
            let p = component_data.perm_components[i as usize];
            let comp = component_data
                .components
                .get_set(component_data.first_unfixed[p as usize]);
            if comp != current_component {
                current_component = comp;
                current_component_start = i;
                component_data
                    .perm_component_starts
                    .push(current_component_start);
            }
        }

        debug_assert_eq!(
            component_data.perm_component_starts.len(),
            component_data.component_starts.len()
        );
        component_data.perm_component_starts.push(num_used_perms);

        component_data.component_starts.push(self.num_active_cols);

        component_data
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::isFullOrbitope (lines 1521-1690)
    /// C++: `bool isFullOrbitope(const ComponentData & componentData, HighsInt component, HighsSymmetries & symmetries)`
    pub fn is_full_orbitope(
        &mut self,
        model: &HighsLp,
        component_data: &HighsSymmetryDetection_ComponentData,
        component: HighsInt,
        symmetries: &mut HighsSymmetries,
    ) -> bool {
        let component_size = component_data.component_starts[(component + 1) as usize]
            - component_data.component_starts[component as usize];
        if component_size == 1 {
            return false;
        }

        // check that component acts only on binary variables
        for i in component_data.component_starts[component as usize]
            ..component_data.component_starts[(component + 1) as usize]
        {
            let col = component_data.component_sets[i as usize] as usize;
            if model.integrality_[col] == HighsVarType::kContinuous
                || model.col_lower_[col] != 0.0
                || model.col_upper_[col] != 1.0
            {
                return false;
            }
        }

        // check that the number of unfixed variables in the first permutation is even
        let p0 = component_data.perm_components
            [component_data.perm_component_starts[component as usize] as usize];
        if component_data.num_unfixed[p0 as usize] & 1 != 0 {
            return false;
        }

        // check that the other permutations in the component have the same number of
        // unfixed variables as the first one
        for k in (component_data.perm_component_starts[component as usize] + 1)
            ..component_data.perm_component_starts[(component + 1) as usize]
        {
            let p = component_data.perm_components[k as usize];
            if component_data.num_unfixed[p as usize] != component_data.num_unfixed[p0 as usize] {
                return false;
            }
        }

        // all unfixed variables in the permutation must be part of a two cycle
        // and each two cycle defines one row of the orbitope. Hence the number
        // of unfixed variables in each permutation must be even and the number of
        // rows is the number of unfixed variables divided by two
        let orbitope_num_rows = component_data.num_unfixed[p0 as usize] >> 1;

        // if this component is a full orbitope the component size must be
        // divisible by the number of orbits and the result is the size of each
        // orbit
        let orbitope_orbit_size = component_size / orbitope_num_rows;
        if orbitope_orbit_size * orbitope_num_rows != component_size {
            return false;
        }

        // the number of permutations must be n-1 where n is the size of the
        // orbits in the orbitope
        let component_num_perms = component_data.perm_component_starts[(component + 1) as usize]
            - component_data.perm_component_starts[component as usize];
        if component_num_perms != orbitope_orbit_size - 1 {
            return false;
        }

        // set up the first two columns of the orbitope matrix based on the first
        // permutation.
        let mut orbitope_matrix = HighsOrbitopeMatrix::default();
        orbitope_matrix.matrix.resize(component_size as usize, -1);
        orbitope_matrix.num_rows = orbitope_num_rows;
        orbitope_matrix.row_length = orbitope_orbit_size;
        debug_assert_eq!(
            component_size,
            orbitope_matrix.num_rows * orbitope_matrix.row_length
        );

        // PORT NOTE: C++ `const HighsInt* perm` pointers into `symmetries.permutations` are
        // stored as offsets of the first entry.
        let mut perm = p0 as usize * self.num_active_cols as usize;
        let mut col_set: HighsHashTable<HighsInt, ()> = HighsHashTable::new();
        let mut m: HighsInt = 0;
        for j in 0..self.num_active_cols {
            let j_image_pos =
                self.vertex_position[symmetries.permutations[perm + j as usize] as usize];
            if j_image_pos <= j {
                continue;
            }
            if m == orbitope_num_rows {
                return false;
            }

            // permutation should consist of two cycles
            if symmetries.permutations[perm + j_image_pos as usize]
                != self.vertex_ground_set[j as usize]
            {
                return false;
            }

            orbitope_matrix.matrix[m as usize] = self.vertex_ground_set[j as usize];
            orbitope_matrix.matrix[(orbitope_num_rows + m) as usize] =
                symmetries.permutations[perm + j as usize];

            // Remember set of variables of the orbtiope matrix. Each variable should
            // occur only once, otherwise the permutation do not work out. Since we
            if !col_set.insert(HighsHashTableEntry::new(
                self.vertex_ground_set[j as usize],
                (),
            )) {
                return false;
            }
            if !col_set.insert(HighsHashTableEntry::new(
                symmetries.permutations[perm + j as usize],
                (),
            )) {
                return false;
            }
            m += 1;
        }

        let mut num_cols_added: HighsInt = 2;
        let mut tried_left_extension = false;

        while num_cols_added < orbitope_matrix.row_length {
            if col_set.size() != (num_cols_added * orbitope_matrix.num_rows) as i64 {
                return false;
            }

            // C++: HighsInt* thisCol = &orbitopeMatrix(0, numColsAdded);
            // C++: HighsInt* prevCol = &orbitopeMatrix(0, numColsAdded - 1);
            // PORT NOTE: column pointers become offsets into `orbitope_matrix.matrix`.
            let this_col = (num_cols_added * orbitope_matrix.num_rows) as usize;
            let mut prev_col = ((num_cols_added - 1) * orbitope_matrix.num_rows) as usize;

            let mut found_cand = false;

            loop {
                let move_pos = self.vertex_position[orbitope_matrix.matrix[prev_col] as usize];
                // PORT NOTE: C++ sets `perm = nullptr` here. It is read only after `foundCand`
                // is set, which requires this loop to assign it, so no stale offset is read.

                for k in (component_data.perm_component_starts[component as usize] + 1)
                    ..component_data.perm_component_starts[(component + 1) as usize]
                {
                    let p = component_data.perm_components[k as usize];
                    perm = p as usize * self.num_active_cols as usize;

                    if symmetries.permutations[perm + move_pos as usize]
                        != self.vertex_ground_set[move_pos as usize]
                        && col_set
                            .find(&symmetries.permutations[perm + move_pos as usize])
                            .is_none()
                    {
                        found_cand = true;
                        break;
                    }
                }

                if !found_cand && !tried_left_extension {
                    // if we fail to find a permutation extending the last column directly
                    // after we set up the first two columns we might try to extend column
                    // zero instead of the previous one
                    prev_col = 0;
                    tried_left_extension = true;
                    continue;
                }

                break;
            }

            if !found_cand {
                return false;
            }

            for j in 0..orbitope_matrix.num_rows {
                let next_vertex = symmetries.permutations[perm
                    + self.vertex_position[orbitope_matrix.matrix[prev_col + j as usize] as usize]
                        as usize];

                orbitope_matrix.matrix[this_col + j as usize] = next_vertex;

                // check if this is a two cycle
                if symmetries.permutations
                    [perm + self.vertex_position[next_vertex as usize] as usize]
                    != orbitope_matrix.matrix[prev_col + j as usize]
                {
                    return false;
                }

                if !col_set.insert(HighsHashTableEntry::new(
                    orbitope_matrix.matrix[this_col + j as usize],
                    (),
                )) {
                    return false;
                }
            }

            num_cols_added += 1;
        }

        if col_set.size() != component_size as i64 {
            return false;
        }

        for &col in orbitope_matrix.matrix.iter() {
            symmetries
                .column_to_orbitope
                .insert(HighsHashTableEntry::new(
                    col,
                    symmetries.orbitopes.len() as HighsInt,
                ));
        }

        symmetries.orbitopes.push(orbitope_matrix);

        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::initializeDetection (lines 1692-1698)
    /// C++: `bool initializeDetection()`
    pub fn initialize_detection(&mut self) -> bool {
        self.initialize_hash_values();
        self.partition_refinement();
        self.remove_fix_points();
        if self.num_active_cols == 0 {
            return false;
        }
        true
    }

    /// HiGHS: HighsSymmetry.cpp::HighsSymmetryDetection::run (lines 1700-1936)
    /// C++: `void run(HighsSymmetries & symmetries)`
    /// 1.16.0 #3535: `run(symmetries, timer, timeLimit)`; `time_limit` is `Some` only under
    /// `mip_sym_hash64` (None = the 1.15.1 loop, and the C++ `timer == nullptr` default).
    pub fn run(
        &mut self,
        model: &HighsLp,
        symmetries: &mut HighsSymmetries,
        time_limit: Option<(&HighsTimer, f64)>,
    ) {
        debug_assert!(self.num_active_cols != 0);
        self.initialize_ground_set();
        self.curr_node_certificate.clear();
        self.cell_creation_stack.clear();
        self.create_node();
        let max_perms: HighsInt = 64000000 / self.num_active_cols;
        // PORT NOTE: C++ `HighsSplitDeque* workerDeque` and its `checkInterrupt()` calls are
        // dropped: multithreading is skipped (spawn runs inline, no interrupts).
        while !self.node_stack.is_empty() {
            // HiGHS #3535 (mip_sym_hash64): HighsSymmetry.cpp line 1680 at v1.16.0, stop once
            // the time limit is reached (the checkInterrupt of line 1679 has no counterpart)
            if let Some((timer, limit)) = time_limit {
                if timer.read() >= limit {
                    break;
                }
            }
            let target_cell = self.select_target_cell();
            if target_cell == -1 {
                if self.first_leave_partition.is_empty() {
                    self.first_leave_partition = self.current_partition.clone();
                    self.first_leave_certificate = self.curr_node_certificate.clone();
                    self.best_leave_certificate = self.curr_node_certificate.clone();
                    self.first_leave_graph = self.dump_current_graph();
                    self.first_path_depth = self.node_stack.len() as HighsInt;
                    self.best_path_depth = self.node_stack.len() as HighsInt;
                    self.first_leave_prefix_len = self.curr_node_certificate.len() as HighsInt;
                    self.best_leave_prefix_len = self.curr_node_certificate.len() as HighsInt;

                    let mut backtrack_depth = self.first_path_depth - 1;
                    while backtrack_depth > 0
                        && !self.is_from_binary_column(
                            model,
                            self.node_stack[(backtrack_depth - 1) as usize].target_cell,
                        )
                    {
                        backtrack_depth -= 1;
                    }
                    self.switch_to_next_node(backtrack_depth);
                } else {
                    let mut wrong_cell: HighsInt = -1;
                    let mut backtrack_depth = self.node_stack.len() as HighsInt - 1;
                    debug_assert_eq!(
                        self.curr_node_certificate.len(),
                        self.first_leave_certificate.len()
                    );
                    let cert_len = self.curr_node_certificate.len() as HighsInt;
                    if self.first_leave_prefix_len == cert_len
                        || self.best_leave_prefix_len == cert_len
                    {
                        if self.first_leave_prefix_len == cert_len
                            && self.compare_current_graph(&self.first_leave_graph, &mut wrong_cell)
                        {
                            let k = self.num_automorphisms & 63;
                            self.num_automorphisms += 1;
                            // C++: HighsInt* permutation = automorphisms.data() + k * numVertices;
                            let perm_base = (k * self.num_vertices) as usize;
                            for i in 0..self.num_vertices as usize {
                                let first_leave_col = self.first_leave_partition[i];
                                let pos = self.vertex_position[self.current_partition[i] as usize];
                                self.automorphisms[perm_base + pos as usize] = first_leave_col;
                            }

                            let mut report = false;
                            for i in 0..self.num_vertices {
                                let iu = i as usize;
                                let perm_value = self.automorphisms[perm_base + iu];
                                let ground_value = self.vertex_ground_set[iu];
                                if self.merge_orbits(perm_value, ground_value)
                                    && i < self.num_active_cols
                                {
                                    debug_assert!(
                                        self.automorphisms[perm_base + iu] < self.num_col
                                    );
                                    report = true;
                                }
                            }

                            if report {
                                symmetries.permutations.extend_from_slice(
                                    &self.automorphisms
                                        [perm_base..perm_base + self.num_active_cols as usize],
                                );
                                symmetries.num_perms += 1;
                                if symmetries.num_perms == max_perms {
                                    break;
                                }
                            }
                            backtrack_depth = cpp_min(backtrack_depth, self.first_path_depth);
                        } else if !self.best_leave_partition.is_empty()
                            && self.best_leave_prefix_len == cert_len
                            && self.compare_current_graph(&self.best_leave_graph, &mut wrong_cell)
                        {
                            let k = self.num_automorphisms & 63;
                            self.num_automorphisms += 1;
                            // C++: HighsInt* permutation = automorphisms.data() + k * numVertices;
                            let perm_base = (k * self.num_vertices) as usize;
                            for i in 0..self.num_vertices as usize {
                                let best_leave_col = self.best_leave_partition[i];
                                let pos = self.vertex_position[self.current_partition[i] as usize];
                                self.automorphisms[perm_base + pos as usize] = best_leave_col;
                            }

                            let mut report = false;
                            for i in 0..self.num_vertices {
                                let iu = i as usize;
                                let perm_value = self.automorphisms[perm_base + iu];
                                let ground_value = self.vertex_ground_set[iu];
                                if self.merge_orbits(perm_value, ground_value)
                                    && i < self.num_active_cols
                                {
                                    debug_assert!(
                                        self.automorphisms[perm_base + iu] < self.num_col
                                    );
                                    report = true;
                                }
                            }

                            if report {
                                symmetries.permutations.extend_from_slice(
                                    &self.automorphisms
                                        [perm_base..perm_base + self.num_active_cols as usize],
                                );
                                symmetries.num_perms += 1;
                                if symmetries.num_perms == max_perms {
                                    break;
                                }
                            }

                            backtrack_depth = cpp_min(backtrack_depth, self.best_path_depth);
                        } else if self.best_leave_prefix_len < cert_len
                            && self.curr_node_certificate[self.best_leave_prefix_len as usize]
                                > self.best_leave_certificate[self.best_leave_prefix_len as usize]
                        {
                            // certificate value is lexicographically above the smallest one
                            // seen so far, so we might be able to backtrack to a higher level
                            let mut possible_backtrack_depth = self.first_path_depth - 1;
                            while self.node_stack[possible_backtrack_depth as usize].certificate_end
                                <= self.best_leave_prefix_len
                            {
                                possible_backtrack_depth += 1;
                            }

                            backtrack_depth = cpp_min(possible_backtrack_depth, backtrack_depth);
                        } else {
                            // This case can be caused by a hash collision which was now
                            // detected in the graph comparison call. The graph comparison call
                            // will return the cell where the vertex neighbourhood caused a
                            // mismatch on the edges. This would have been detected by
                            // an exact partition refinement when we targeted that cell the last
                            // time, so that is where we can backtrack to.
                            let mut possible_backtrack_depth = backtrack_depth;
                            while possible_backtrack_depth >= 0 {
                                if self.node_stack[possible_backtrack_depth as usize].target_cell
                                    == wrong_cell
                                {
                                    backtrack_depth = possible_backtrack_depth;
                                    break;
                                }
                                possible_backtrack_depth -= 1;
                            }
                        }
                    } else {
                        // leave must have a lexicographically smaller certificate value
                        // than the current best leave, because its prefix length is smaller
                        // than the best leaves and it would have been already pruned if
                        // it's certificate value was larger unless it is equal to the first
                        // leaf nodes certificate value which is caught by the first case
                        // of the if condition. Hence, having a lexicographically smaller
                        // certificate value than the best leave is the only way to get
                        // here.
                        let bl = self.best_leave_prefix_len as usize;
                        debug_assert!(
                            self.best_leave_certificate[bl] > self.curr_node_certificate[bl]
                                && self.best_leave_certificate[..bl]
                                    == self.curr_node_certificate[..bl]
                        );
                        self.best_leave_certificate = self.curr_node_certificate.clone();
                        self.best_leave_graph = self.dump_current_graph();
                        self.best_leave_partition = self.current_partition.clone();
                        self.best_path_depth = self.node_stack.len() as HighsInt;
                        self.best_leave_prefix_len = self.curr_node_certificate.len() as HighsInt;
                    }

                    self.switch_to_next_node(backtrack_depth);
                }

                // PORT NOTE: C++ workerDeque->checkInterrupt() dropped (no threads).
            } else {
                let node_len = self.node_stack.len();
                self.node_stack[node_len - 1].target_cell = target_cell;
                let success = self.determine_next_to_distinguish();
                debug_assert!(success);
                if !self.distinguish_vertex(target_cell) {
                    let depth = self.node_stack.len() as HighsInt - 1;
                    self.switch_to_next_node(depth);
                    continue;
                }
                if !self.partition_refinement() {
                    let depth = self.node_stack.len() as HighsInt;
                    self.switch_to_next_node(depth);
                    continue;
                }

                self.create_node();
            }
        }

        symmetries.num_generators = symmetries.num_perms;
        if symmetries.num_perms > 0 {
            self.vertex_position.resize(self.num_col as usize, 0);

            let mut component_data = self.compute_component_data(symmetries);
            let num_components = component_data.num_components();

            for i in 0..num_components {
                if component_data.component_size(i) == 1 {
                    continue;
                }

                self.is_full_orbitope(model, &component_data, i, symmetries);
            }

            let mut deleted_perms: HighsHashTable<HighsInt, ()> = HighsHashTable::new();
            for p in 0..symmetries.num_perms {
                let perm = p as usize * self.num_active_cols as usize;
                for i in 0..self.num_active_cols {
                    let value = symmetries.permutations[perm + i as usize];
                    if value != self.vertex_ground_set[i as usize]
                        && symmetries.column_to_orbitope.find(&value).is_some()
                    {
                        deleted_perms.insert(HighsHashTableEntry::new(p, ()));
                        break;
                    }
                }
            }

            // check which columns have non-trivial orbits
            let mut num_fixed: HighsInt = 0;
            for i in 0..self.num_active_cols {
                let iu = i as usize;
                let ground_value = self.vertex_ground_set[iu];
                let orbit = self.get_orbit(ground_value);
                if self.orbit_size[orbit as usize] == 1
                    || symmetries.column_to_orbitope.find(&ground_value).is_some()
                {
                    self.vertex_position[ground_value as usize] = -1;
                    self.vertex_ground_set[iu] = -1;
                    num_fixed += 1;
                }
            }

            if num_fixed != 0 {
                // now compress symmetries and the groundset to only contain the unfixed
                // columns and columns not handled by full orbitopes
                let mut p: HighsInt = 0;
                let mut perms: usize = 0;
                let perm_end: usize = symmetries.num_perms as usize * self.num_active_cols as usize;
                let mut perm_output: usize = 0;
                while perms != perm_end {
                    if deleted_perms.find(&p).is_none() {
                        for i in 0..self.num_active_cols {
                            if self.vertex_ground_set[i as usize] == -1 {
                                continue;
                            }

                            let value = symmetries.permutations[perms + i as usize];
                            symmetries.permutations[perm_output] = value;
                            perm_output += 1;
                        }
                    } else {
                        symmetries.num_perms -= 1;
                    }
                    perms += self.num_active_cols as usize;
                    p += 1;
                }

                let mut out_pos: HighsInt = 0;
                for i in 0..self.num_active_cols {
                    let iu = i as usize;
                    if self.vertex_ground_set[iu] == -1 {
                        continue;
                    }

                    self.vertex_ground_set[out_pos as usize] = self.vertex_ground_set[iu];
                    let moved = self.vertex_ground_set[out_pos as usize];
                    self.vertex_position[moved as usize] = out_pos;
                    out_pos += 1;
                }

                self.num_active_cols -= num_fixed;
                debug_assert_eq!(
                    perm_output,
                    symmetries.num_perms as usize * self.num_active_cols as usize
                );
            }

            self.vertex_ground_set
                .resize(self.num_active_cols as usize, 0);
            symmetries.permutation_columns = std::mem::take(&mut self.vertex_ground_set);
            symmetries.column_position = std::mem::take(&mut self.vertex_position);
            symmetries.permutations.resize(
                symmetries.num_perms as usize * self.num_active_cols as usize,
                0,
            );
        }
    }
}

impl MatrixColumn {
    /// HiGHS: HighsSymmetry.cpp::MatrixColumn::operator== (lines 1165-1167)
    /// C++: `bool operator==(const MatrixColumn & other) const`
    pub fn op_eq(&self, other: &MatrixColumn) -> bool {
        // PORT NOTE: C++ memcmp over five uint32_t fields with no padding is field-wise equality.
        self.cost == other.cost
            && self.lb == other.lb
            && self.ub == other.ub
            && self.integral == other.integral
            && self.len == other.len
    }
}

impl MatrixRow {
    /// HiGHS: HighsSymmetry.cpp::MatrixRow::operator== (lines 1175-1177)
    /// C++: `bool operator==(const MatrixRow & other) const`
    pub fn op_eq(&self, other: &MatrixRow) -> bool {
        // PORT NOTE: C++ memcmp over three uint32_t fields with no padding is field-wise equality.
        self.lb == other.lb && self.ub == other.ub && self.len == other.len
    }
}
