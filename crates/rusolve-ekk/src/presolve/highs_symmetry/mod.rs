//! HiGHS source: presolve/HighsSymmetry.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

mod part1;
pub use part1::*;
mod part2;
pub use part2::*;

/// HiGHS: HighsSymmetry.h::HighsMatrixColoring (lines 27-52)
#[derive(Clone, Default, Debug)]
pub struct HighsMatrixColoring {
    /// C++: `std::map<double, u32> colorMap`
    pub color_map: BTreeMap<OrderedF64, u32>,
    /// C++: `double tolerance`
    pub tolerance: f64,
}

/// HiGHS: HighsSymmetry.h::StabilizerOrbitWorkspace (lines 57-61)
#[derive(Clone, Default, Debug)]
pub struct StabilizerOrbitWorkspace {
    /// C++: `std::vector<HighsInt> orbitPartition`
    pub orbit_partition: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> orbitSize`
    pub orbit_size: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> linkCompressionStack`
    pub link_compression_stack: Vec<HighsInt>,
}

/// HiGHS: HighsSymmetry.h::StabilizerOrbits (lines 63-72)
///
/// Removed C++ members:
/// * C++: `const HighsSymmetries * symmetries`: removed; methods that use it take `symmetries: &HighsSymmetries`.
#[derive(Clone, Default, Debug)]
pub struct StabilizerOrbits {
    /// C++: `std::vector<HighsInt> orbitCols`
    pub orbit_cols: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> orbitStarts`
    pub orbit_starts: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> stabilizedCols`
    pub stabilized_cols: Vec<HighsInt>,
}

/// HiGHS: HighsSymmetry.h::HighsOrbitopeMatrix (lines 74-120)
#[derive(Clone, Default, Debug)]
pub struct HighsOrbitopeMatrix {
    /// C++: `HighsInt rowLength`
    pub row_length: HighsInt,
    /// C++: `HighsInt numRows`
    pub num_rows: HighsInt,
    /// C++: `HighsInt numSetPackingRows`
    pub num_set_packing_rows: HighsInt,
    /// C++: `HighsHashTable<HighsInt, HighsInt> columnToRow`
    pub column_to_row: HighsHashTable<HighsInt, HighsInt>,
    /// C++: `std::vector<RowPackingStatus> rowIsSetPacking`
    pub row_is_set_packing: Vec<HighsOrbitopeMatrix_RowPackingStatus>,
    /// C++: `std::vector<HighsInt> matrix`
    pub matrix: Vec<HighsInt>,
}

/// HiGHS: C++ plain `enum Type` nested in `HighsOrbitopeMatrix` (values are associated consts).
pub type HighsOrbitopeMatrix_Type = HighsInt;
impl HighsOrbitopeMatrix {
    pub const kFull: HighsInt = 0;
    pub const kPacking: HighsInt = 1;
}

/// HiGHS: C++ plain `enum RowPackingStatus` nested in `HighsOrbitopeMatrix` (values are associated consts).
pub type HighsOrbitopeMatrix_RowPackingStatus = i8;
impl HighsOrbitopeMatrix {
    pub const kRowUndetermined: i8 = -1;
    pub const kRowNotPacking: i8 = 0;
    pub const kRowPacking: i8 = 1;
    pub const kRowPackingNegated: i8 = 2;
}

/// HiGHS: HighsSymmetry.h::HighsSymmetries (lines 122-162)
#[derive(Clone, Debug)]
pub struct HighsSymmetries {
    /// C++: `std::vector<HighsInt> permutationColumns`
    pub permutation_columns: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> permutations`
    pub permutations: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> orbitPartition`
    pub orbit_partition: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> orbitSize`
    pub orbit_size: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> columnPosition`
    pub column_position: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> linkCompressionStack`
    pub link_compression_stack: Vec<HighsInt>,
    /// C++: `std::vector<HighsOrbitopeMatrix> orbitopes`
    pub orbitopes: Vec<HighsOrbitopeMatrix>,
    /// C++: `HighsHashTable<HighsInt, HighsInt> columnToOrbitope`
    pub column_to_orbitope: HighsHashTable<HighsInt, HighsInt>,
    /// C++: `HighsInt numPerms = 0`
    pub num_perms: HighsInt,
    /// C++: `HighsInt numGenerators = 0`
    pub num_generators: HighsInt,
}

impl Default for HighsSymmetries {
    fn default() -> Self {
        Self {
            permutation_columns: Default::default(),
            permutations: Default::default(),
            orbit_partition: Default::default(),
            orbit_size: Default::default(),
            column_position: Default::default(),
            link_compression_stack: Default::default(),
            orbitopes: Default::default(),
            column_to_orbitope: Default::default(),
            num_perms: 0,
            num_generators: 0,
        }
    }
}

/// HiGHS: HighsSymmetry.h::HighsSymmetryDetection (lines 164-304)
///
/// Removed C++ members:
/// * C++ `const HighsLp* model`: removed; methods take `model: &HighsLp`.
#[derive(Clone, Default, Debug)]
pub struct HighsSymmetryDetection {
    /// C++: `std::vector<HighsInt> Gstart`
    pub gstart: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> Gend`
    pub gend: Vec<HighsInt>,
    /// C++: `std::vector<std::pair<HighsInt, HighsUInt>> Gedge`
    pub gedge: Vec<(HighsInt, HighsUInt)>,
    /// C++: `std::vector<std::pair<HighsInt, HighsUInt>> edgeBuffer`
    pub edge_buffer: Vec<(HighsInt, HighsUInt)>,
    /// C++: `std::vector<HighsInt> currentPartition`
    pub current_partition: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> currentPartitionLinks`
    pub current_partition_links: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> vertexToCell`
    pub vertex_to_cell: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> vertexPosition`
    pub vertex_position: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> vertexGroundSet`
    pub vertex_ground_set: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> orbitPartition`
    pub orbit_partition: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> orbitSize`
    pub orbit_size: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> cellCreationStack`
    pub cell_creation_stack: Vec<HighsInt>,
    /// C++: `std::vector<std::uint8_t> cellInRefinementQueue`
    pub cell_in_refinement_queue: Vec<u8>,
    /// C++: `std::vector<HighsInt> refinementQueue`
    pub refinement_queue: Vec<HighsInt>,
    /// C++ `std::vector<HighsInt*> distinguishCands`: positions into `currentPartition`.
    pub distinguish_cands: Vec<usize>,
    /// C++: `std::vector<HighsInt> automorphisms`
    pub automorphisms: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> linkCompressionStack`
    pub link_compression_stack: Vec<HighsInt>,
    /// C++: `std::vector<u32> currNodeCertificate`
    pub curr_node_certificate: Vec<u32>,
    /// C++: `std::vector<u32> firstLeaveCertificate`
    pub first_leave_certificate: Vec<u32>,
    /// C++: `std::vector<u32> bestLeaveCertificate`
    pub best_leave_certificate: Vec<u32>,
    /// C++: `std::vector<HighsInt> firstLeavePartition`
    pub first_leave_partition: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> bestLeavePartition`
    pub best_leave_partition: Vec<HighsInt>,
    /// C++: `HighsHashTable<HighsInt, u32> vertexHash` (1.15.1); `u64` since 1.16.0 #3535.
    /// Stored as `u64` on both paths: with `mip_sym_hash64` off every value is a
    /// `sparse_combine32` result (< 2^31), so orders, equality and the `u32` reads in
    /// `split_cell` match the 1.15.1 `u32` table.
    pub vertex_hash: HighsHashTable<HighsInt, u64>,
    /// C++: `HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> firstLeaveGraph`
    pub first_leave_graph: HighsHashTable<(HighsInt, HighsInt, HighsUInt), ()>,
    /// C++: `HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> bestLeaveGraph`
    pub best_leave_graph: HighsHashTable<(HighsInt, HighsInt, HighsUInt), ()>,
    /// C++: `HighsInt firstLeavePrefixLen`
    pub first_leave_prefix_len: HighsInt,
    /// C++: `HighsInt bestLeavePrefixLen`
    pub best_leave_prefix_len: HighsInt,
    /// C++: `HighsInt firstPathDepth`
    pub first_path_depth: HighsInt,
    /// C++: `HighsInt bestPathDepth`
    pub best_path_depth: HighsInt,
    /// C++: `HighsInt numAutomorphisms`
    pub num_automorphisms: HighsInt,
    /// C++: `HighsInt numCol`
    pub num_col: HighsInt,
    /// C++: `HighsInt numRow`
    pub num_row: HighsInt,
    /// C++: `HighsInt numVertices`
    pub num_vertices: HighsInt,
    /// C++: `HighsInt numActiveCols`
    pub num_active_cols: HighsInt,
    /// C++: `std::vector<Node> nodeStack`
    pub node_stack: Vec<HighsSymmetryDetection_Node>,
    /// HiGHS 1.16 #3182 (h16_sym_refine) flag, copied from `HighsOptions::h16_sym_refine` by
    /// `HighsMipSolverData::start_symmetry_detection` before the task is spawned.
    /// Off (`false`, the `Default`) reproduces 1.15.1; read only by `partition_refinement`.
    pub h16_sym_refine: bool,
    /// HiGHS #3535 (mip_sym_hash64) flag, copied from `HighsOptions::mip_sym_hash64` by
    /// `HighsMipSolverData::start_symmetry_detection` before the task is spawned. On: vertex
    /// hashes use the 61-bit `sparse_combine` (1.16.0); off (`false`, the `Default`) keeps the
    /// 1.15.1 `sparse_combine32`.
    pub mip_sym_hash64: bool,
}

/// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::Node (lines 214-219)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsSymmetryDetection_Node {
    /// C++: `HighsInt stackStart`
    pub stack_start: HighsInt,
    /// C++: `HighsInt certificateEnd`
    pub certificate_end: HighsInt,
    /// C++: `HighsInt targetCell`
    pub target_cell: HighsInt,
    /// C++: `HighsInt lastDistiguished`
    pub last_distiguished: HighsInt,
}

/// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::ComponentData (lines 263-291)
#[derive(Clone, Default, Debug)]
pub struct HighsSymmetryDetection_ComponentData {
    /// C++: `HighsDisjointSets<> components`
    pub components: HighsDisjointSets,
    /// C++: `std::vector<HighsInt> componentStarts`
    pub component_starts: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> componentSets`
    pub component_sets: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> componentNumOrbits`
    pub component_num_orbits: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> componentNumber`
    pub component_number: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> permComponentStarts`
    pub perm_component_starts: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> permComponents`
    pub perm_components: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> firstUnfixed`
    pub first_unfixed: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> numUnfixed`
    pub num_unfixed: Vec<HighsInt>,
}

/// HiGHS: HighsSymmetry.cpp::MatrixColumn (lines 1158-1168)
#[derive(Clone, Copy, Default, Debug)]
pub struct MatrixColumn {
    /// C++: `uint32_t cost`
    pub cost: u32,
    /// C++: `uint32_t lb`
    pub lb: u32,
    /// C++: `uint32_t ub`
    pub ub: u32,
    /// C++: `uint32_t integral`
    pub integral: u32,
    /// C++: `uint32_t len`
    pub len: u32,
}

impl PartialEq for MatrixColumn {
    fn eq(&self, other: &Self) -> bool {
        self.op_eq(other)
    }
}

/// HiGHS: HighsSymmetry.cpp::MatrixRow (lines 1170-1178)
#[derive(Clone, Copy, Default, Debug)]
pub struct MatrixRow {
    /// C++: `uint32_t lb`
    pub lb: u32,
    /// C++: `uint32_t ub`
    pub ub: u32,
    /// C++: `uint32_t len`
    pub len: u32,
}

impl PartialEq for MatrixRow {
    fn eq(&self, other: &Self) -> bool {
        self.op_eq(other)
    }
}

impl HighsMatrixColoring {
    /// HiGHS: HighsSymmetry.h::HighsMatrixColoring::HighsMatrixColoring (lines 36-38)
    /// C++: `HighsMatrixColoring(double tolerance)`
    pub fn new(tolerance: f64) -> Self {
        // colorMap({{0.0, 1}, {1.0, 2}, {-kHighsInf, 3}, {kHighsInf, 4}})
        let mut color_map: BTreeMap<OrderedF64, u32> = BTreeMap::new();
        color_map.insert(OrderedF64(0.0), 1);
        color_map.insert(OrderedF64(1.0), 2);
        color_map.insert(OrderedF64(-kHighsInf), 3);
        color_map.insert(OrderedF64(kHighsInf), 4);
        Self {
            color_map,
            tolerance,
        }
    }

    /// HiGHS: HighsSymmetry.h::HighsMatrixColoring::color (lines 40-51)
    /// C++: `u32 color(double value)`
    pub fn color(&mut self, value: f64) -> u32 {
        // iterator points to smallest element in map which fulfills key >= value - tolerance
        let found: Option<(f64, u32)> = self
            .color_map
            .range(OrderedF64(value - self.tolerance)..)
            .next()
            .map(|(k, col_color)| (k.0, *col_color));
        // no such element, or this element's key is > value + tolerance: create a new color
        let (create, color) = match found {
            None => (true, 0),
            Some((key, col_color)) => (key > value + self.tolerance, col_color),
        };
        if create {
            // static_cast<u32>(colorMap.size()) + 1, evaluated before the emplace
            let new_color = self.color_map.len() as u32 + 1;
            self.color_map.insert(OrderedF64(value), new_color);
            return new_color;
        }
        color
    }
}

impl HighsOrbitopeMatrix {
    /// HiGHS: HighsSymmetry.h::HighsOrbitopeMatrix::entry (lines 92-94)
    /// C++: `HighsInt & entry(HighsInt i, HighsInt j)`
    pub fn entry(&mut self, i: HighsInt, j: HighsInt) -> &mut HighsInt {
        // matrix[i + static_cast<size_t>(j) * numRows]
        let idx = (i as usize).wrapping_add((j as usize).wrapping_mul(self.num_rows as usize));
        &mut self.matrix[idx]
    }

    /// HiGHS: HighsSymmetry.h::HighsOrbitopeMatrix::entry (lines 96-98)
    /// C++: `const HighsInt & entry(HighsInt i, HighsInt j) const`
    pub fn entry_2(&self, i: HighsInt, j: HighsInt) -> &HighsInt {
        // matrix[i + static_cast<size_t>(j) * numRows]
        let idx = (i as usize).wrapping_add((j as usize).wrapping_mul(self.num_rows as usize));
        &self.matrix[idx]
    }

    /// HiGHS: HighsSymmetry.h::HighsOrbitopeMatrix::operator() (lines 100-100)
    /// C++: `HighsInt & operator()(HighsInt i, HighsInt j)`
    pub fn op_call(&mut self, i: HighsInt, j: HighsInt) -> &mut HighsInt {
        self.entry(i, j)
    }

    /// HiGHS: HighsSymmetry.h::HighsOrbitopeMatrix::operator() (lines 102-104)
    /// C++: `const HighsInt & operator()(HighsInt i, HighsInt j) const`
    pub fn op_call_2(&self, i: HighsInt, j: HighsInt) -> &HighsInt {
        self.entry_2(i, j)
    }
}

impl HighsSymmetries {
    /// HiGHS: HighsSymmetry.h::HighsSymmetries::mergeOrbits (lines 139-141)
    /// C++: `void mergeOrbits(HighsInt col1, HighsInt col2)`
    pub fn merge_orbits_2(&mut self, col1: HighsInt, col2: HighsInt) {
        // PORT NOTE: the C++ passes its own members by reference; they are moved out for the call (rule (f)).
        let mut orbit_partition = std::mem::take(&mut self.orbit_partition);
        let mut orbit_size = std::mem::take(&mut self.orbit_size);
        let mut link_compression_stack = std::mem::take(&mut self.link_compression_stack);
        self.merge_orbits(
            col1,
            col2,
            &mut orbit_partition,
            &mut orbit_size,
            &mut link_compression_stack,
        );
        self.orbit_partition = orbit_partition;
        self.orbit_size = orbit_size;
        self.link_compression_stack = link_compression_stack;
    }

    /// HiGHS: HighsSymmetry.h::HighsSymmetries::getOrbit (lines 144-146)
    /// C++: `HighsInt getOrbit(HighsInt col)`
    pub fn get_orbit_2(&mut self, col: HighsInt) -> HighsInt {
        // PORT NOTE: members moved out for the call, as in merge_orbits_2 (rule (f)).
        let mut orbit_partition = std::mem::take(&mut self.orbit_partition);
        let mut link_compression_stack = std::mem::take(&mut self.link_compression_stack);
        let orbit = self.get_orbit(col, &mut orbit_partition, &mut link_compression_stack);
        self.orbit_partition = orbit_partition;
        self.link_compression_stack = link_compression_stack;
        orbit
    }

    /// HiGHS: HighsSymmetry.h::HighsSymmetries::getBranchingColumn (lines 150-158)
    /// C++: `HighsInt getBranchingColumn(const std::vector<double> & colLower, const std::vector<double> & colUpper, HighsInt col) const`
    pub fn get_branching_column(
        &self,
        col_lower: &[f64],
        col_upper: &[f64],
        col: HighsInt,
    ) -> HighsInt {
        if self.column_to_orbitope.size() == 0 {
            return col;
        }
        match self.column_to_orbitope.find(&col) {
            Some(&orbitope) if self.orbitopes[orbitope as usize].num_set_packing_rows != 0 => {
                self.orbitopes[orbitope as usize].get_branching_column(col_lower, col_upper, col)
            }
            _ => col,
        }
    }
}

impl HighsSymmetryDetection {
    /// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::cellSize (lines 257-259)
    /// C++: `HighsInt cellSize(HighsInt cell) const`
    pub fn cell_size(&self, cell: HighsInt) -> HighsInt {
        self.current_partition_links[cell as usize] - cell
    }
}

impl HighsSymmetryDetection_ComponentData {
    /// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::ComponentData::getComponentByIndex (lines 274-276)
    /// C++: `HighsInt getComponentByIndex(HighsInt compIndex) const`
    pub fn get_component_by_index(&self, comp_index: HighsInt) -> HighsInt {
        self.component_number[comp_index as usize]
    }

    /// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::ComponentData::numComponents (lines 277-279)
    /// C++: `HighsInt numComponents() const`
    pub fn num_components(&self) -> HighsInt {
        self.component_starts.len() as HighsInt - 1
    }

    /// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::ComponentData::componentSize (lines 280-282)
    /// C++: `HighsInt componentSize(HighsInt component) const`
    pub fn component_size(&self, component: HighsInt) -> HighsInt {
        self.component_starts[(component + 1) as usize] - self.component_starts[component as usize]
    }

    /// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::ComponentData::getVertexComponent (lines 284-286)
    /// C++: `HighsInt getVertexComponent(HighsInt vertexPosition)`
    pub fn get_vertex_component(&mut self, vertex_position: HighsInt) -> HighsInt {
        self.components.get_set(vertex_position)
    }

    /// HiGHS: HighsSymmetry.h::HighsSymmetryDetection::ComponentData::getPermuationComponent (lines 288-290)
    /// C++: `HighsInt getPermuationComponent(HighsInt permIndex)`
    pub fn get_permuation_component(&mut self, perm_index: HighsInt) -> HighsInt {
        let first = self.first_unfixed[perm_index as usize];
        self.components.get_set(first)
    }
}

/// C++ `HighsHashHelpers::hash(std::tuple<HighsInt, HighsInt, HighsUInt>)` (12-byte key).
/// HiGHS: HighsHash.h::HighsHashHelpers::hash, branch `9 <= sizeof(T) <= 16` (lines 582-592).
impl HighsHashKey for (HighsInt, HighsInt, HighsUInt) {
    fn highs_hash(&self) -> u64 {
        // memcpy of the 12 bytes; bytes[3] stays 0 because sizeof(T) < 16
        let bytes: [u32; 4] = [self.0 as u32, self.1 as u32, self.2, 0];
        (pair_hash(0, bytes[0], bytes[1]) ^ (pair_hash(1, bytes[2], bytes[3]) >> 32))
            .wrapping_mul(fibonacci_muliplier())
    }
}

/// C++ `HighsHashHelpers::hash(const MatrixColumn&)`: 20-byte key, HighsHash.h 17-24-byte branch
/// (lines 594-606). The struct's five u32 fields are the bytes, in declaration order.
impl HighsHashKey for MatrixColumn {
    fn highs_hash(&self) -> u64 {
        // memcpy of the 20 bytes; bytes[5] stays 0 because sizeof(T) < 24
        let bytes: [u32; 6] = [self.cost, self.lb, self.ub, self.integral, self.len, 0];
        (pair_hash(0, bytes[0], bytes[1])
            ^ (pair_hash(1, bytes[2], bytes[3]).wrapping_add(pair_hash(2, bytes[4], bytes[5]))
                >> 32))
            .wrapping_mul(fibonacci_muliplier())
    }
}

/// C++ `HighsHashHelpers::hash(const MatrixRow&)`: 12-byte key, HighsHash.h 9-16-byte branch
/// (lines 578-590). The struct's three u32 fields are the bytes, in declaration order.
impl HighsHashKey for MatrixRow {
    fn highs_hash(&self) -> u64 {
        // memcpy of the 12 bytes; bytes[3] stays 0 because sizeof(T) < 16
        let bytes: [u32; 4] = [self.lb, self.ub, self.len, 0];
        (pair_hash(0, bytes[0], bytes[1]) ^ (pair_hash(1, bytes[2], bytes[3]) >> 32))
            .wrapping_mul(fibonacci_muliplier())
    }
}
