//! HiGHS source: ipm/hipo/factorhighs/Symbolic.h, ipm/hipo/factorhighs/Symbolic.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The symbolic factorisation (output of `Analyse::run`, input of
//! `Factorise`). `friend class Analyse` becomes `pub` fields. The one-line
//! getters (Symbolic.cpp lines 16-49) are complete.

use super::super::auxiliary::logger::{Logger, fix, integer, sci, textline};
use super::super::{Int, Int64};

/// HiGHS: Symbolic.h `class Symbolic`
#[derive(Clone, Debug, Default)]
pub struct Symbolic {
    pub parallel_tree_: bool,
    pub parallel_node_: bool,
    /// No initialiser in C++; `Analyse::run` sets it.
    pub block_size_: Int,
    pub n_: Int,
    pub nz_: Int64,
    pub sn_: Int,
    pub fillin_: f64,
    pub flops_: f64,
    pub spops_: f64,
    pub critops_: f64,
    pub artificial_nz_: Int64,
    pub artificial_ops_: f64,
    pub serial_storage_: f64,
    pub largest_front_: Int,
    pub largest_sn_: Int,
    pub sn_size_1_: Int,
    pub sn_size_10_: Int,
    pub sn_size_100_: Int,
    pub iperm_: Vec<Int>,
    pub rows_: Vec<Int>,
    pub ptr_: Vec<Int64>,
    pub sn_parent_: Vec<Int>,
    pub sn_start_: Vec<Int>,
    pub relind_cols_: Vec<Int>,
    pub relind_clique_: Vec<Vec<Int>>,
    pub consecutive_sums_: Vec<Vec<Int>>,
    pub pivot_sign_: Vec<Int>,
    pub clique_block_start_: Vec<Vec<Int64>>,
    pub max_stack_size_: Int64,
    pub tree_depth_: Int,
}

impl Symbolic {
    /// HiGHS: Symbolic.cpp::Symbolic::Symbolic (line 9)
    pub fn new() -> Self {
        Self::default()
    }

    /// HiGHS: Symbolic.cpp::Symbolic::setParallel (lines 11-14)
    pub fn setParallel(&mut self, par_tree: bool, par_node: bool) {
        self.parallel_tree_ = par_tree;
        self.parallel_node_ = par_node;
    }

    /// HiGHS: Symbolic.cpp::Symbolic::nz (line 16)
    pub fn nz(&self) -> Int64 {
        self.nz_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::flops (line 17)
    pub fn flops(&self) -> f64 {
        self.flops_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::spops (line 18)
    pub fn spops(&self) -> f64 {
        self.spops_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::critops (line 19)
    pub fn critops(&self) -> f64 {
        self.critops_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::blockSize (line 20)
    pub fn blockSize(&self) -> Int {
        self.block_size_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::size (line 21)
    pub fn size(&self) -> Int {
        self.n_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::sn (line 22)
    pub fn sn(&self) -> Int {
        self.sn_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::fillin (line 23)
    pub fn fillin(&self) -> f64 {
        self.fillin_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::rows (line 24)
    pub fn rows(&self, i: Int64) -> Int {
        self.rows_[i as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::ptr (line 25)
    pub fn ptr(&self, i: Int) -> Int64 {
        self.ptr_[i as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::snStart (line 26)
    pub fn snStart(&self, i: Int) -> Int {
        self.sn_start_[i as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::snParent (line 27)
    pub fn snParent(&self, i: Int) -> Int {
        self.sn_parent_[i as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::relindCols (line 28)
    pub fn relindCols(&self, i: Int) -> Int {
        self.relind_cols_[i as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::relindClique (line 29)
    pub fn relindClique(&self, i: Int, j: Int) -> Int {
        self.relind_clique_[i as usize][j as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::consecutiveSums (lines 30-32)
    pub fn consecutiveSums(&self, i: Int, j: Int) -> Int {
        self.consecutive_sums_[i as usize][j as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::cliqueBlockStart (lines 33-35)
    pub fn cliqueBlockStart(&self, sn: Int, bl: Int) -> Int64 {
        self.clique_block_start_[sn as usize][bl as usize]
    }
    /// HiGHS: Symbolic.cpp::Symbolic::cliqueSize (lines 36-38)
    pub fn cliqueSize(&self, sn: Int) -> Int64 {
        *self.clique_block_start_[sn as usize].last().unwrap()
    }
    /// HiGHS: Symbolic.cpp::Symbolic::maxStackSize (line 39)
    pub fn maxStackSize(&self) -> Int64 {
        self.max_stack_size_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::parTree (line 40)
    pub fn parTree(&self) -> bool {
        self.parallel_tree_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::parNode (line 41)
    pub fn parNode(&self) -> bool {
        self.parallel_node_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::storage (line 42)
    pub fn storage(&self) -> f64 {
        self.serial_storage_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::depth (line 43)
    pub fn depth(&self) -> Int {
        self.tree_depth_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::ptr (line 45; vector overload)
    pub fn ptr_vec(&self) -> &[Int64] {
        &self.ptr_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::iperm (line 46)
    pub fn iperm(&self) -> &[Int] {
        &self.iperm_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::snParent (line 47; vector overload)
    pub fn snParent_vec(&self) -> &[Int] {
        &self.sn_parent_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::snStart (line 48; vector overload)
    pub fn snStart_vec(&self) -> &[Int] {
        &self.sn_start_
    }
    /// HiGHS: Symbolic.cpp::Symbolic::pivotSign (line 49)
    pub fn pivotSign(&self) -> &[Int] {
        &self.pivot_sign_
    }

    /// HiGHS: Symbolic.cpp::Symbolic::print (lines 66-98)
    /// C++ default arguments: `verbose = false`.
    pub fn print(&self, logger: &Logger, verbose: bool) {
        let mut log_stream = String::new();
        log_stream.push_str("\nFactorisation statistics\n");
        log_stream.push_str(&textline("Size:"));
        log_stream.push_str(&sci(self.n_ as f64, 0, 2));
        log_stream.push('\n');
        log_stream.push_str(&textline("Nnz:"));
        log_stream.push_str(&sci(self.nz_ as f64, 0, 2));
        log_stream.push('\n');
        log_stream.push_str(&textline("Fill-in:"));
        log_stream.push_str(&fix(self.fillin_, 0, 2));
        log_stream.push('\n');
        log_stream.push_str(&textline("Serial memory:"));
        log_stream.push_str(&memoryString(self.serial_storage_));
        log_stream.push('\n');
        log_stream.push_str(&textline("Flops:"));
        log_stream.push_str(&sci(self.flops_, 0, 1));
        log_stream.push('\n');
        if verbose {
            log_stream.push_str(&textline("Sparse ops:"));
            log_stream.push_str(&sci(self.spops_, 0, 1));
            log_stream.push('\n');
            log_stream.push_str(&textline("Critical ops:"));
            log_stream.push_str(&sci(self.critops_, 0, 1));
            log_stream.push('\n');
            log_stream.push_str(&textline("Max tree speedup:"));
            log_stream.push_str(&fix(self.flops_ / self.critops_, 0, 2));
            log_stream.push('\n');
            log_stream.push_str(&textline("Tree depth:"));
            log_stream.push_str(&integer(self.tree_depth_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Artificial nz:"));
            log_stream.push_str(&sci(self.artificial_nz_ as f64, 0, 1));
            log_stream.push('\n');
            log_stream.push_str(&textline("Artificial ops:"));
            log_stream.push_str(&sci(self.artificial_ops_, 0, 1));
            log_stream.push('\n');
            log_stream.push_str(&textline("Largest front:"));
            log_stream.push_str(&integer(self.largest_front_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Largest supernode:"));
            log_stream.push_str(&integer(self.largest_sn_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Supernodes:"));
            log_stream.push_str(&integer(self.sn_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Sn size <= 1:"));
            log_stream.push_str(&integer(self.sn_size_1_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Sn size <= 10:"));
            log_stream.push_str(&integer(self.sn_size_10_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Sn size <= 100:"));
            log_stream.push_str(&integer(self.sn_size_100_, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Sn avg size:"));
            log_stream.push_str(&sci(self.n_ as f64 / self.sn_ as f64, 0, 1));
            log_stream.push('\n');
        }
        logger.print(&log_stream);
    }
}

/// HiGHS: Symbolic.cpp::memoryString (lines 51-64, file static)
pub fn memoryString(mem: f64) -> String {
    let mut ss = String::new();

    if mem < 1024.0 {
        ss.push_str(&sci(mem, 0, 1));
        ss.push_str(" B");
    } else if mem < 1024.0 * 1024.0 {
        ss.push_str(&sci(mem / 1024.0, 0, 1));
        ss.push_str(" KB");
    } else if mem < 1024.0 * 1024.0 * 1024.0 {
        ss.push_str(&sci(mem / 1024.0 / 1024.0, 0, 1));
        ss.push_str(" MB");
    } else {
        ss.push_str(&sci(mem / 1024.0 / 1024.0 / 1024.0, 0, 1));
        ss.push_str(" GB");
    }

    ss
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn getters_and_memory_string() {
        let mut s = Symbolic::new();
        s.clique_block_start_ = vec![vec![0, 3, 7]];
        assert_eq!(s.cliqueSize(0), 7);
        s.setParallel(true, false);
        assert!(s.parTree() && !s.parNode());
        assert!(!memoryString(1024.0).is_empty());
    }
}
