//! HiGHS source: mip/HighsGFkSolve.h, mip/HighsGFkSolve.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsGFkSolve.h::HighsGFk<k>::powk, specialisations included (lines 37-54).
///
/// The C++ recursion runs on the template parameter `k`. Stable Rust has no
/// const-generic recursion, so the same recursion runs on the runtime value. The
/// branch structure is C++'s: `k == 2` -> `a * a`, `k == 3` -> `a * a * a`,
/// even `k` -> `powk<2>(powk<k/2>(a))`, odd `k` -> `powk<k-1>(a) * a`.
fn gfk_powk(k: HighsInt, a: u32) -> u32 {
    if k == 2 {
        a.wrapping_mul(a)
    } else if k == 3 {
        a.wrapping_mul(a).wrapping_mul(a)
    } else if (k & 1) == 0 {
        gfk_powk(2, gfk_powk(k / 2, a))
    } else {
        gfk_powk(k - 1, a).wrapping_mul(a)
    }
}

/// HiGHS: HighsGFkSolve.h::HighsGFk<k>::inverse, specialisations included (lines 40-58).
///
/// Runtime form of the template, for the same reason as `gfk_powk`.
fn gfk_inverse(k: HighsInt, a: u32) -> u32 {
    if k == 2 {
        1
    } else if k == 3 {
        a
    } else {
        gfk_powk(k - 2, a) % (k as u32)
    }
}

/// HiGHS: HighsGFkSolve.h::HighsGFkSolve (lines 61-437)
#[derive(Clone, Default, Debug)]
pub struct HighsGFkSolve {
    /// C++: `HighsInt numCol`
    pub num_col: HighsInt,
    /// C++: `HighsInt numRow`
    pub num_row: HighsInt,
    /// C++: `std::vector<HighsInt> Arow`
    pub arow: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> Acol`
    pub acol: Vec<HighsInt>,
    /// C++: `std::vector<unsigned int> Avalue`
    pub avalue: Vec<u32>,
    /// C++: `std::vector<HighsInt> rowsize`
    pub rowsize: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colsize`
    pub colsize: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colhead`
    pub colhead: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> Anext`
    pub anext: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> Aprev`
    pub aprev: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> rowroot`
    pub rowroot: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ARleft`
    pub a_rleft: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> ARright`
    pub a_rright: Vec<HighsInt>,
    /// C++: `std::vector<unsigned int> rhs`
    pub rhs: Vec<u32>,
    /// C++: `std::vector<HighsInt> factorColPerm`
    pub factor_col_perm: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> factorRowPerm`
    pub factor_row_perm: Vec<HighsInt>,
    /// C++: `std::vector<int8_t> colBasisStatus`
    pub col_basis_status: Vec<i8>,
    /// C++: `std::vector<int8_t> rowUsed`
    pub row_used: Vec<i8>,
    /// C++: `std::vector<HighsInt> iterstack`
    pub iterstack: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> rowpositions`
    pub rowpositions: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> rowposColsizes`
    pub rowpos_colsizes: Vec<HighsInt>,
    /// C++: `std::priority_queue<HighsInt, std::vector<HighsInt>, std::greater<HighsInt>> freeslots`
    pub freeslots: BinaryHeap<Reverse<HighsInt>>,
}

/// HiGHS: HighsGFkSolve.h::HighsGFkSolve::SolutionEntry (lines 115-122)
#[derive(Clone, Copy, Default, Debug)]
pub struct HighsGFkSolve_SolutionEntry {
    /// C++: `HighsInt index`
    pub index: HighsInt,
    /// C++: `HighsUInt weight`
    pub weight: HighsUInt,
}

/// `std::set`/`std::sort` ordering = C++ `operator<` (equivalence = neither is less).
impl PartialEq for HighsGFkSolve_SolutionEntry {
    fn eq(&self, other: &Self) -> bool {
        !self.op_lt(other) && !other.op_lt(self)
    }
}
impl Eq for HighsGFkSolve_SolutionEntry {}
impl PartialOrd for HighsGFkSolve_SolutionEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HighsGFkSolve_SolutionEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.op_lt(other) {
            std::cmp::Ordering::Less
        } else if other.op_lt(self) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

impl HighsGFkSolve {
    /// HiGHS: HighsGFkSolve.cpp::HighsGFkSolve::link (lines 14-28)
    /// C++: `void link(HighsInt pos)`
    pub fn link(&mut self, pos: HighsInt) {
        let p = pos as usize;
        let col = self.acol[p] as usize;
        self.anext[p] = self.colhead[col];
        self.aprev[p] = -1;
        self.colhead[col] = pos;
        if self.anext[p] != -1 {
            let next = self.anext[p] as usize;
            self.aprev[next] = pos;
        }

        self.colsize[col] += 1;

        let row = self.arow[p] as usize;
        highs_splay_link(
            pos,
            &mut self.rowroot[row],
            &mut self.a_rleft,
            &mut self.a_rright,
            &self.acol,
        );
        self.rowsize[row] += 1;
    }

    /// HiGHS: HighsGFkSolve.cpp::HighsGFkSolve::unlink (lines 30-51)
    /// C++: `void unlink(HighsInt pos)`
    pub fn unlink(&mut self, pos: HighsInt) {
        let p = pos as usize;
        let next = self.anext[p];
        let prev = self.aprev[p];

        if next != -1 {
            self.aprev[next as usize] = prev;
        }

        let col = self.acol[p] as usize;
        if prev != -1 {
            self.anext[prev as usize] = next;
        } else {
            self.colhead[col] = next;
        }
        self.colsize[col] -= 1;

        let row = self.arow[p] as usize;
        highs_splay_unlink(
            pos,
            &mut self.rowroot[row],
            &mut self.a_rleft,
            &mut self.a_rright,
            &self.acol,
        );
        self.rowsize[row] -= 1;

        self.avalue[p] = 0;
        self.freeslots.push(Reverse(pos));
    }

    /// HiGHS: HighsGFkSolve.cpp::HighsGFkSolve::storeRowPositions (lines 53-68)
    /// C++: `void storeRowPositions(HighsInt pos)`
    pub fn store_row_positions(&mut self, pos: HighsInt) {
        if pos == -1 {
            return;
        }

        debug_assert!(self.iterstack.is_empty());
        self.iterstack.push(pos);
        // C++ do/while: the stack is non-empty on entry, so while-let is the same loop.
        while let Some(p) = self.iterstack.pop() {
            let pu = p as usize;
            self.rowpositions.push(p);
            self.rowpos_colsizes
                .push(self.colsize[self.acol[pu] as usize]);

            if self.a_rleft[pu] != -1 {
                self.iterstack.push(self.a_rleft[pu]);
            }
            if self.a_rright[pu] != -1 {
                self.iterstack.push(self.a_rright[pu]);
            }
        }
    }

    /// HiGHS: HighsGFkSolve.cpp::HighsGFkSolve::findNonzero (lines 70-82)
    /// C++: `HighsInt findNonzero(HighsInt row, HighsInt col)`
    pub fn find_nonzero(&mut self, row: HighsInt, col: HighsInt) -> HighsInt {
        let r = row as usize;
        if self.rowroot[r] == -1 {
            return -1;
        }

        self.rowroot[r] = highs_splay(
            col,
            self.rowroot[r],
            &mut self.a_rleft,
            &mut self.a_rright,
            &self.acol,
        );

        if self.acol[self.rowroot[r] as usize] == col {
            return self.rowroot[r];
        }

        -1
    }

    /// HiGHS: HighsGFkSolve.cpp::HighsGFkSolve::addNonzero (lines 84-106)
    /// C++: `void addNonzero(HighsInt row, HighsInt col, unsigned int val)`
    pub fn add_nonzero(&mut self, row: HighsInt, col: HighsInt, val: u32) {
        // The assert calls findNonzero, which splays the row tree. It runs only in
        // debug builds, as the C++ assert does with NDEBUG unset.
        debug_assert!(self.find_nonzero(row, col) == -1);
        let pos: HighsInt;
        if self.freeslots.is_empty() {
            pos = self.avalue.len() as HighsInt;
            self.avalue.push(val);
            self.arow.push(row);
            self.acol.push(col);
            self.anext.push(-1);
            self.aprev.push(-1);
            self.a_rleft.push(-1);
            self.a_rright.push(-1);
        } else {
            pos = self.freeslots.peek().unwrap().0;
            self.freeslots.pop();
            let p = pos as usize;
            self.avalue[p] = val;
            self.arow[p] = row;
            self.acol[p] = col;
            self.aprev[p] = -1;
        }

        self.link(pos);
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::dropIfZero (lines 106-108)
    /// C++: `void dropIfZero(HighsInt pos)`
    pub fn drop_if_zero(&mut self, pos: HighsInt) {
        if self.avalue[pos as usize] == 0 {
            self.unlink(pos);
        }
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::getArow (lines 125-125)
    /// C++: `const std::vector<HighsInt> & getArow() const`
    pub fn get_arow(&self) -> &[HighsInt] {
        &self.arow
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::getAcol (lines 126-126)
    /// C++: `const std::vector<HighsInt> & getAcol() const`
    pub fn get_acol(&self) -> &[HighsInt] {
        &self.acol
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::getAvalue (lines 127-127)
    /// C++: `const std::vector<unsigned int> & getAvalue() const`
    pub fn get_avalue(&self) -> &[u32] {
        &self.avalue
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::numNonzeros (lines 128-128)
    /// C++: `HighsInt numNonzeros() const`
    pub fn num_nonzeros(&self) -> HighsInt {
        (self.avalue.len() - self.freeslots.len()) as HighsInt
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::fromCSC (lines 131-175)
    /// C++: `void fromCSC()`
    pub fn from_csc<const K: u32, const K_NUM_RHS: i32>(
        &mut self,
        aval: &[i64],
        aindex: &[HighsInt],
        astart: &[HighsInt],
        num_row: HighsInt,
    ) {
        self.avalue.clear();
        self.acol.clear();
        self.arow.clear();

        self.freeslots = BinaryHeap::new();

        self.num_col = astart.len() as HighsInt - 1;
        self.num_row = num_row;

        let num_col = self.num_col as usize;
        self.colhead = vec![-1; num_col];
        self.colsize = vec![0; num_col];

        self.rhs = vec![0; (K_NUM_RHS * num_row) as usize];
        self.rowroot = vec![-1; num_row as usize];
        self.rowsize = vec![0; num_row as usize];

        self.avalue.reserve(aval.len());
        self.acol.reserve(aval.len());
        self.arow.reserve(aval.len());

        for i in 0..num_col {
            let start = astart[i] as usize;
            let end = astart[i + 1] as usize;
            for j in start..end {
                // assert(Aval[j] == (int64_t)Aval[j]) is trivially true for i64.
                let mut val: i64 = aval[j] % (K as i64);
                if val == 0 {
                    continue;
                }

                if val < 0 {
                    val += K as i64;
                }
                debug_assert!(val >= 0);

                self.avalue.push(val as u32);
                self.acol.push(i as HighsInt);
                self.arow.push(aindex[j]);
            }
        }

        let nnz = self.avalue.len() as HighsInt;
        self.anext.resize(nnz as usize, 0);
        self.aprev.resize(nnz as usize, 0);
        self.a_rleft.resize(nnz as usize, 0);
        self.a_rright.resize(nnz as usize, 0);
        for pos in 0..nnz {
            self.link(pos);
        }
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::setRhs (lines 177-180)
    /// C++: `void setRhs()`
    pub fn set_rhs<const K: u32, const K_NUM_RHS: i32>(
        &mut self,
        row: HighsInt,
        val: i64,
        rhs_index: i32,
    ) {
        // C++ `(unsigned int)std::abs(val)` truncates the int64 to 32 bits.
        self.rhs[(K_NUM_RHS * row + rhs_index) as usize] = (val.wrapping_abs() as u32) % K;
    }

    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::solve (lines 182-436)
    /// C++: `void solve()`
    pub fn solve<const K: u32, const K_NUM_RHS: i32>(
        &mut self,
        report_solution: &mut dyn FnMut(&[HighsGFkSolve_SolutionEntry], i32),
    ) {
        // C++ kNumRhs is an int template parameter.
        let k_num_rhs: HighsInt = K_NUM_RHS;
        let num_rhs = K_NUM_RHS as usize;

        // std::priority_queue ordered on .first only: ties are possible, so the
        // libc++ heap order from cpp_std is required.
        let cmp_prio = |a: &(HighsInt, HighsInt), b: &(HighsInt, HighsInt)| a.0 > b.0;
        let mut pqueue: CppPriorityQueue<(HighsInt, HighsInt)> = CppPriorityQueue::default();

        for i in 0..self.num_col {
            pqueue.push((self.colsize[i as usize], i), cmp_prio);
        }

        let max_pivot = cpp_min(self.num_row, self.num_col);
        self.factor_col_perm.clear();
        self.factor_row_perm.clear();
        self.factor_col_perm.reserve(max_pivot as usize);
        self.factor_row_perm.reserve(max_pivot as usize);
        self.col_basis_status = vec![0; self.num_col as usize];
        self.row_used = vec![0; self.num_row as usize];
        let mut num_pivot: HighsInt = 0;

        while !pqueue.empty() {
            let top = *pqueue.top();
            pqueue.pop(cmp_prio);
            let (old_col_size, pivot_col) = top;

            if self.colsize[pivot_col as usize] == 0 {
                continue;
            }

            debug_assert!(self.col_basis_status[pivot_col as usize] == 0);

            if self.colsize[pivot_col as usize] != old_col_size {
                pqueue.push((self.colsize[pivot_col as usize], pivot_col), cmp_prio);
                continue;
            }

            let mut pivot: HighsInt = -1;
            let mut pivot_row: HighsInt = -1;
            let mut pivot_row_len: HighsInt = kHighsIInf;
            let mut coliter = self.colhead[pivot_col as usize];
            while coliter != -1 {
                let row = self.arow[coliter as usize];
                if self.row_used[row as usize] == 0 && self.rowsize[row as usize] < pivot_row_len {
                    pivot_row_len = self.rowsize[row as usize];
                    pivot_row = row;
                    pivot = coliter;
                }
                coliter = self.anext[coliter as usize];
            }

            debug_assert!(pivot != -1);
            debug_assert!(self.acol[pivot as usize] == pivot_col);
            debug_assert!(self.arow[pivot as usize] == pivot_row);
            debug_assert!(self.avalue[pivot as usize] > 0);
            debug_assert!(self.avalue[pivot as usize] < K);

            let pivot_inverse: u32 = gfk_inverse(K as HighsInt, self.avalue[pivot as usize]);
            debug_assert!(self.avalue[pivot as usize].wrapping_mul(pivot_inverse) % K == 1);

            self.rowpositions.clear();
            self.rowpos_colsizes.clear();
            let root = self.rowroot[pivot_row as usize];
            self.store_row_positions(root);
            debug_assert!(pivot_row_len == self.rowpositions.len() as HighsInt);

            let mut coliter = self.colhead[pivot_col as usize];
            while coliter != -1 {
                let next = self.anext[coliter as usize];
                if coliter != pivot {
                    debug_assert!(self.acol[coliter as usize] == pivot_col);
                    debug_assert!(self.arow[coliter as usize] != pivot_row);
                    debug_assert!(self.avalue[coliter as usize] != 0);

                    let row = self.arow[coliter as usize];
                    if self.row_used[row as usize] == 0 {
                        let pivot_row_scale: u32 =
                            pivot_inverse.wrapping_mul(K - self.avalue[coliter as usize]);

                        for i in 0..k_num_rhs {
                            let r = (k_num_rhs * row + i) as usize;
                            let p = (k_num_rhs * pivot_row + i) as usize;
                            self.rhs[r] = self.rhs[r]
                                .wrapping_add(pivot_row_scale.wrapping_mul(self.rhs[p]))
                                % K;
                        }

                        // rowpositions is not modified inside this loop; index it so
                        // that findNonzero and addNonzero may borrow self.
                        for idx in 0..self.rowpositions.len() {
                            let pivot_row_pos = self.rowpositions[idx];
                            let nonzero_pos = self.find_nonzero(
                                self.arow[coliter as usize],
                                self.acol[pivot_row_pos as usize],
                            );

                            if nonzero_pos == -1 {
                                debug_assert!(self.acol[pivot_row_pos as usize] != pivot_col);
                                let val: u32 = pivot_row_scale
                                    .wrapping_mul(self.avalue[pivot_row_pos as usize])
                                    % K;
                                if val != 0 {
                                    self.add_nonzero(row, self.acol[pivot_row_pos as usize], val);
                                }
                            } else {
                                let np = nonzero_pos as usize;
                                self.avalue[np] = self.avalue[np].wrapping_add(
                                    pivot_row_scale
                                        .wrapping_mul(self.avalue[pivot_row_pos as usize]),
                                ) % K;
                                debug_assert!(
                                    self.acol[pivot_row_pos as usize] != pivot_col
                                        || self.avalue[np] == 0
                                );
                                self.drop_if_zero(nonzero_pos);
                            }
                        }
                    }
                }
                coliter = next;
            }

            num_pivot += 1;
            self.factor_col_perm.push(pivot_col);
            self.factor_row_perm.push(pivot_row);
            self.col_basis_status[pivot_col as usize] = 1;
            self.row_used[pivot_row as usize] = 1;
            if num_pivot == max_pivot {
                break;
            }

            for i in 0..pivot_row_len {
                let iu = i as usize;
                debug_assert!(self.arow[self.rowpositions[iu] as usize] == pivot_row);
                let col = self.acol[self.rowpositions[iu] as usize];
                let oldsize = self.rowpos_colsizes[iu];

                // we only want to count rows that are not used so far, so we need to
                // decrease the counter for all columns in the pivot row by one
                self.colsize[col as usize] -= 1;

                // the column size should never get negative
                debug_assert!(self.colsize[col as usize] >= 0);

                if self.colsize[col as usize] == 0 {
                    continue;
                }

                // the pivot column should occur in zero unused rows
                debug_assert!(col != pivot_col);

                // only reinsert the column if the size is smaller, as it otherwise
                // either is already at the correct position in the queue, or is at a
                // too good position and will be lazily reinserted when it is extracted
                // from the queue and the size does not match
                if self.colsize[col as usize] < oldsize {
                    pqueue.push((self.colsize[col as usize], col), cmp_prio);
                }
            }
        }

        // check if a solution exists by scanning the linearly dependent rows for
        // nonzero right hand sides
        let mut has_solution: Vec<bool> = vec![false; num_rhs];
        let mut num_rhs_with_solution: HighsInt = 0;
        for rhs_index in 0..num_rhs {
            has_solution[rhs_index] = true;
            for i in 0..self.num_row {
                // if the row was used it is linearly independent
                if self.row_used[i as usize] == 1 {
                    continue;
                }

                // if the row is linearly dependent, the right hand side must be zero,
                // otherwise no solution exists
                if self.rhs[(k_num_rhs * i) as usize + rhs_index] != 0 {
                    has_solution[rhs_index] = false;
                    break;
                }
            }

            num_rhs_with_solution += has_solution[rhs_index] as HighsInt;
        }

        if num_rhs_with_solution == 0 {
            return;
        }

        // now iterate a subset of the basic solutions.
        // When a column leaves the basis we do not allow it to enter again so that
        // we iterate at most one solution for each nonbasic column
        let mut solution: Vec<Vec<HighsGFkSolve_SolutionEntry>> = vec![Vec::new(); num_rhs];
        for rhs_index in 0..num_rhs {
            if has_solution[rhs_index] {
                solution[rhs_index].reserve(self.num_col as usize);
            }
        }

        let num_factor_rows = self.factor_row_perm.len() as HighsInt;

        // create vector for swapping different columns into the basis
        // For each column we want to iterate one basic solution where the
        // column is basic
        let mut basis_swaps: Vec<(HighsInt, HighsInt)> = Vec::new();
        debug_assert!(self.iterstack.is_empty());
        for i in (0..num_factor_rows).rev() {
            let row = self.factor_row_perm[i as usize];
            self.iterstack.push(self.rowroot[row as usize]);

            while let Some(rowpos) = self.iterstack.pop() {
                let rp = rowpos as usize;
                debug_assert!(rowpos != -1);

                if self.a_rleft[rp] != -1 {
                    self.iterstack.push(self.a_rleft[rp]);
                }
                if self.a_rright[rp] != -1 {
                    self.iterstack.push(self.a_rright[rp]);
                }

                let col = self.acol[rp];
                if self.col_basis_status[col as usize] != 0 {
                    continue;
                }

                self.col_basis_status[col as usize] = -1;
                basis_swaps.push((i, col));
            }
        }

        let mut basis_swap_pos: HighsInt = 0;

        // solval is written before it is read for every rhs with a solution.
        let mut solval: Vec<u32> = vec![0; num_rhs];

        // C++ do/while loop
        loop {
            let mut performed_basis_swap = false;

            for rhs_index in 0..num_rhs {
                solution[rhs_index].clear();
            }

            for i in (0..num_factor_rows).rev() {
                let row = self.factor_row_perm[i as usize];

                for rhs_index in 0..num_rhs {
                    if !has_solution[rhs_index] {
                        continue;
                    }
                    solval[rhs_index] = 0;

                    for solentry in solution[rhs_index].iter() {
                        let pos = self.find_nonzero(row, solentry.index);
                        if pos != -1 {
                            solval[rhs_index] = solval[rhs_index].wrapping_add(
                                self.avalue[pos as usize].wrapping_mul(solentry.weight),
                            );
                        }
                    }

                    solval[rhs_index] = self.rhs[(k_num_rhs * row) as usize + rhs_index]
                        .wrapping_add(K)
                        .wrapping_sub(solval[rhs_index] % K);
                }

                let col = self.factor_col_perm[i as usize];
                let pos = self.find_nonzero(row, col);
                debug_assert!(pos != -1);
                let col_val_inverse: u32 = gfk_inverse(K as HighsInt, self.avalue[pos as usize]);

                for rhs_index in 0..num_rhs {
                    if !has_solution[rhs_index] {
                        continue;
                    }
                    // assert(solval >= 0) is trivially true for unsigned int.
                    debug_assert!(col_val_inverse != 0);

                    solval[rhs_index] = solval[rhs_index].wrapping_mul(col_val_inverse) % K;

                    debug_assert!(solval[rhs_index] < K);

                    // only record nonzero solution values
                    if solval[rhs_index] != 0 {
                        solution[rhs_index].push(HighsGFkSolve_SolutionEntry {
                            index: col,
                            weight: solval[rhs_index],
                        });
                    }
                }
            }

            for rhs_index in 0..num_rhs {
                if has_solution[rhs_index] {
                    report_solution(solution[rhs_index].as_slice(), rhs_index as i32);
                }
            }

            if basis_swap_pos < basis_swaps.len() as HighsInt {
                let basis_index = basis_swaps[basis_swap_pos as usize].0;
                let entering_col = basis_swaps[basis_swap_pos as usize].1;
                let leaving_col = self.factor_col_perm[basis_index as usize];
                debug_assert!(self.col_basis_status[leaving_col as usize] == 1);
                self.factor_col_perm[basis_index as usize] = entering_col;
                self.col_basis_status[entering_col as usize] = 1;
                self.col_basis_status[leaving_col as usize] = 0;
                performed_basis_swap = true;
                basis_swap_pos += 1;
            }

            if !performed_basis_swap {
                break;
            }
        }
    }
}

impl HighsGFkSolve_SolutionEntry {
    /// HiGHS: HighsGFkSolve.h::HighsGFkSolve::SolutionEntry::operator< (lines 119-121)
    /// C++: `bool operator<(const SolutionEntry & other) const`
    pub fn op_lt(&self, other: &HighsGFkSolve_SolutionEntry) -> bool {
        self.index < other.index
    }
}

/// HiGHS: HighsGFkSolve.h::HighsGFk<k> (lines 25-58): arithmetic in GF(k)
/// (`HighsGFk<2>` / `HighsGFk<3>` specialisations are the `K == 2` / `K == 3`
/// branches of the bodies).
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsGFk<const K: HighsInt>;

impl<const K: HighsInt> HighsGFk<K> {
    /// HiGHS: HighsGFkSolve.h::HighsGFk::powk (lines 39, 45, 51-54)
    pub fn powk(a: u32) -> u32 {
        gfk_powk(K, a)
    }
    /// HiGHS: HighsGFkSolve.h::HighsGFk::inverse (lines 40, 46, 56-58)
    pub fn inverse(a: u32) -> u32 {
        gfk_inverse(K, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The runtime-k `gfk_inverse` is a true inverse for every k the mod-k
    /// separator instantiates (HighsModkSeparator.cpp:246-264: 2, 3, 5, 7), and
    /// `gfk_powk` is plain `a^k` there (no u32 wrap for a < k <= 7).
    #[test]
    fn gfk_inverse_is_inverse_for_modk_primes() {
        for k in [2u32, 3, 5, 7] {
            for a in 1..k {
                let inv = gfk_inverse(k as HighsInt, a);
                assert_eq!(a * inv % k, 1, "k={k} a={a}");
                assert_eq!(gfk_powk(k as HighsInt, a), a.pow(k));
            }
        }
        for a in 1..7u32 {
            assert_eq!(a * HighsGFk::<7>::inverse(a) % 7, 1);
        }
    }
}
