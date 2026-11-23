//! HiGHS source: ipm/hipo/factorhighs/Factorise.h, ipm/hipo/factorhighs/Factorise.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The multifrontal numeric factorisation over the supernodal tree. Supports both
//! serial (`S_.parTree()` false) and parallel (`S_.parTree()` true) paths.
//! In parallel mode, children are spawned on a `TaskGroup`, `flag_stop_` is an
//! `AtomicBool`, and child contributions are read from `schur_contribution_`.
//! In serial mode, `stack_` (CliqueStack) is used and `flag_stop_` is checked
//! as a plain bool load (relaxed atomic).

use std::sync::atomic::{AtomicBool, Ordering};

use super::super::auxiliary::logger::Logger;
use super::super::auxiliary::{
    childrenLinkedList, freeVector, permuteSym, reverseLinkedList, transpose_val,
};
use super::super::{Int, Int64};
use super::clique_stack::CliqueStack;
use super::data_collector::DataCollector;
use super::format_handler::FormatHandlerOps;
use super::hybrid_hybrid_format_handler::HybridHybridFormatHandler;
use super::numeric::Numeric;
use super::symbolic::Symbolic;
use super::*;
use crate::lp_data::hconst::kHighsInf;
use crate::parallel::highs_parallel::task_group;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// Structure-only work of `Factorise::new`, recorded once per solve (rusolve lever
/// `ipm_hipo_cache_perm`): for each output entry the index of the input entry it came from.
#[derive(Clone, Default)]
pub struct PermCache {
    valid: bool,
    n: Int,
    in_rows: Vec<Int>,
    in_ptr: Vec<Int>,
    src: Vec<u32>,
    rows: Vec<Int>,
    ptr: Vec<Int>,
}

impl PermCache {
    fn matches(&self, S: &Symbolic, rows: &[Int], ptr: &[Int]) -> bool {
        self.valid && self.n == S.size() && self.in_ptr == ptr && self.in_rows == rows
    }
}

/// HiGHS: Factorise.h `class Factorise`
pub struct Factorise<'a> {
    pub rowsM_: Vec<Int>,
    pub ptrM_: Vec<Int>,
    pub valM_: Vec<f64>,
    pub n_: Int,
    pub nzM_: Int,
    pub S_: &'a Symbolic,
    pub first_child_: Vec<Int>,
    pub next_child_: Vec<Int>,
    pub first_child_reverse_: Vec<Int>,
    pub next_child_reverse_: Vec<Int>,
    pub schur_contribution_: Vec<Vec<f64>>,
    pub sn_columns_: &'a mut Vec<Vec<f64>>,
    pub swaps_: Vec<Vec<Int>>,
    pub pivot_2x2_: Vec<Vec<f64>>,
    pub max_diag_: f64,
    pub min_diag_: f64,
    pub M_norm1_: f64,
    pub total_reg_: Vec<f64>,
    pub regul_: Regul,
    /// C++ `std::atomic<bool> flag_stop_{false}`
    pub flag_stop_: AtomicBool,
    /// rusolve: deadline polled once per supernode; `timed_out_` records that it fired
    pub deadline_: Option<std::time::Instant>,
    pub timed_out_: AtomicBool,
    pub logger_: Option<&'a Logger>,
    pub data_: &'a DataCollector,
    pub stack_: Option<&'a mut CliqueStack>,
}

/// Raw pointer wrapper to allow `Factorise` to be shared across spawned tasks.
///
/// # Safety
/// The C++ protocol guarantees that spawned tasks operate on disjoint supernode
/// indices (`sn`), and the parent task `sync`s before reading child results.
/// The `Factorise` object itself is not mutated concurrently except through
/// `flag_stop_` (atomic) and `schur_contribution_` (written by child, read by
/// parent after sync). The `stack_` is only used in serial mode.
#[derive(Clone, Copy)]
struct SendFactorise<'a> {
    ptr: *mut Factorise<'a>,
}

// SAFETY: The C++ protocol ensures disjoint sn indices per task, sync before
// parent reads child results, and flag_stop_ is atomic. The Factorise fields
// accessed by a task are either local to that task's sn or protected by the
// task dependency graph.
unsafe impl<'a> Send for SendFactorise<'a> {}

impl<'a> SendFactorise<'a> {
    fn get(&self) -> *mut Factorise<'a> {
        self.ptr
    }
}

impl<'a> Factorise<'a> {
    /// Permute to the supernodal order and sort the columns (the original `Factorise` constructor
    /// code, unchanged): returns `(rows, ptr, val)`.
    fn permute_sorted(
        S: &Symbolic,
        rowsM: &[Int],
        ptrM: &[Int],
        valM: &[f64],
    ) -> (Vec<Int>, Vec<Int>, Vec<f64>) {
        let n_ = ptrM.len() as Int - 1;
        // Make a copy of the matrix to be factorised
        let mut rowsM_ = rowsM.to_vec();
        let mut valM_ = valM.to_vec();
        let mut ptrM_ = ptrM.to_vec();

        // Permute the matrix.
        // This also removes any entry not in the lower triangle.
        permuteSym(S.iperm(), &mut ptrM_, &mut rowsM_, &mut valM_, true);

        let nzM_ = *ptrM_.last().unwrap_or(&0);

        // Double transpose to sort columns
        let mut temp_ptr = vec![0 as Int; n_ as usize + 1];
        let mut temp_rows = vec![0 as Int; nzM_ as usize];
        let mut temp_val = vec![0.0f64; nzM_ as usize];
        transpose_val(
            &ptrM_,
            &rowsM_,
            &valM_,
            &mut temp_ptr,
            &mut temp_rows,
            &mut temp_val,
        );
        transpose_val(
            &temp_ptr,
            &temp_rows,
            &temp_val,
            &mut ptrM_,
            &mut rowsM_,
            &mut valM_,
        );
        (rowsM_, ptrM_, valM_)
    }

    /// HiGHS: Factorise.cpp::Factorise::Factorise (lines 17-91)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        S: &'a Symbolic,
        rowsM: &[Int],
        ptrM: &[Int],
        valM: &[f64],
        regul: &Regul,
        logger: Option<&'a Logger>,
        data: &'a DataCollector,
        sn_columns: &'a mut Vec<Vec<f64>>,
        stack: Option<&'a mut CliqueStack>,
    ) -> Self {
        Self::new_cached(
            S, rowsM, ptrM, valM, regul, logger, data, sn_columns, stack, None,
        )
    }

    /// As [`Self::new`], with an optional cache of the structure-only work (symmetric permutation
    /// and the double transpose that sorts columns). rusolve lever `ipm_hipo_cache_perm`: these
    /// steps only move entries, so replaying a recorded source-index map gives bit-identical
    /// values while skipping the O(nnz) integer work of every factorisation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_cached(
        S: &'a Symbolic,
        rowsM: &[Int],
        ptrM: &[Int],
        valM: &[f64],
        regul: &Regul,
        logger: Option<&'a Logger>,
        data: &'a DataCollector,
        sn_columns: &'a mut Vec<Vec<f64>>,
        stack: Option<&'a mut CliqueStack>,
        cache: Option<&mut PermCache>,
    ) -> Self {
        // Input the symmetric matrix to be factorised in CSC format and the symbolic
        // factorisation coming from Analyse.
        // Only the lower triangular part of the matrix is used.

        let n_ = ptrM.len() as Int - 1;

        if n_ != S.size() {
            if let Some(lg) = logger {
                lg.printInfo(
                    "Matrix provided to Factorise has size incompatible with symbolic object.\n",
                );
            }
            // C++ returns after n_ is set; every other member keeps its
            // default initialiser.
            return Self {
                rowsM_: Vec::new(),
                ptrM_: Vec::new(),
                valM_: Vec::new(),
                n_,
                nzM_: 0,
                S_: S,
                first_child_: Vec::new(),
                next_child_: Vec::new(),
                first_child_reverse_: Vec::new(),
                next_child_reverse_: Vec::new(),
                schur_contribution_: Vec::new(),
                sn_columns_: sn_columns,
                swaps_: Vec::new(),
                pivot_2x2_: Vec::new(),
                max_diag_: 0.0,
                min_diag_: 0.0,
                M_norm1_: 0.0,
                total_reg_: Vec::new(),
                regul_: *regul,
                flag_stop_: AtomicBool::new(false),
                deadline_: None,
                timed_out_: AtomicBool::new(false),
                logger_: logger,
                data_: data,
                stack_: stack,
            };
        }

        let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(15);
        let hit = cache.as_ref().is_some_and(|c| c.matches(S, rowsM, ptrM));
        let (rowsM_, ptrM_, valM_) = if hit {
            let c = cache.as_ref().unwrap();
            let valM_: Vec<f64> = c.src.iter().map(|&s| valM[s as usize]).collect();
            (c.rows.clone(), c.ptr.clone(), valM_)
        } else {
            let (rowsM_, ptrM_, valM_) = Self::permute_sorted(S, rowsM, ptrM, valM);
            if let Some(c) = cache {
                // record the source index of every output entry by permuting indices
                let idx: Vec<f64> = (0..valM.len()).map(|i| i as f64).collect();
                let (_, _, src) = Self::permute_sorted(S, rowsM, ptrM, &idx);
                *c = PermCache {
                    valid: true,
                    n: S.size(),
                    in_rows: rowsM.to_vec(),
                    in_ptr: ptrM.to_vec(),
                    src: src.iter().map(|&x| x as u32).collect(),
                    rows: rowsM_.clone(),
                    ptr: ptrM_.clone(),
                };
            }
            (rowsM_, ptrM_, valM_)
        };
        let nzM_ = *ptrM_.last().unwrap_or(&0);

        // create linked lists of children in supernodal elimination tree
        let mut first_child_ = vec![0 as Int; n_ as usize];
        let mut next_child_ = vec![0 as Int; n_ as usize];
        childrenLinkedList(S.snParent_vec(), &mut first_child_, &mut next_child_);

        // create reverse linked lists of children
        let mut first_child_reverse_ = first_child_.clone();
        let mut next_child_reverse_ = next_child_.clone();
        reverseLinkedList(&mut first_child_reverse_, &mut next_child_reverse_);

        // compute largest diagonal entry in absolute value
        let mut max_diag_ = 0.0f64;
        let mut min_diag_ = kHighsInf;
        for col in 0..n_ {
            let val = valM_[ptrM_[col as usize] as usize].abs();
            max_diag_ = cpp_max(max_diag_, val);
            min_diag_ = cpp_min(min_diag_, val);
        }

        // one norm of columns of M
        let mut one_norm_cols = vec![0.0f64; n_ as usize];
        for col in 0..n_ {
            for el in ptrM_[col as usize]..ptrM_[col as usize + 1] {
                let row = rowsM_[el as usize];
                let val = valM_[el as usize];
                one_norm_cols[col as usize] += val.abs();
                if row != col {
                    one_norm_cols[row as usize] += val.abs();
                }
            }
        }
        // std::max_element returns the first largest element
        let mut largest = 0usize;
        for i in 1..one_norm_cols.len() {
            if one_norm_cols[largest] < one_norm_cols[i] {
                largest = i;
            }
        }
        let M_norm1_ = one_norm_cols[largest];

        data.setNorms(M_norm1_, max_diag_);
        drop(_pf);

        Self {
            rowsM_,
            ptrM_,
            valM_,
            n_,
            nzM_,
            S_: S,
            first_child_,
            next_child_,
            first_child_reverse_,
            next_child_reverse_,
            schur_contribution_: Vec::new(),
            sn_columns_: sn_columns,
            swaps_: Vec::new(),
            pivot_2x2_: Vec::new(),
            max_diag_,
            min_diag_,
            M_norm1_,
            total_reg_: Vec::new(),
            regul_: *regul,
            flag_stop_: AtomicBool::new(false),
            deadline_: None,
            timed_out_: AtomicBool::new(false),
            logger_: logger,
            data_: data,
            stack_: stack,
        }
    }

    /// HiGHS: Factorise.cpp::Factorise::permute (lines 93-95)
    pub fn permute(&mut self, iperm: &[Int]) {
        permuteSym(
            iperm,
            &mut self.ptrM_,
            &mut self.rowsM_,
            &mut self.valM_,
            true,
        );
    }

    /// HiGHS: Factorise.cpp::Factorise::processSupernode (lines 97-287)
    pub fn processSupernode(&mut self, sn: Int) {
        // Assemble frontal matrix for supernode sn, perform partial factorisation and
        // store the result.

        // Children are spawned as closures that call back into this object through
        // `this` (the C++ lambda captures `this`).
        let this = SendFactorise {
            ptr: self as *mut Factorise<'a>,
        };

        // C++ line 101: highs::parallel::TaskGroup tg; it lives for the whole body, so
        // every return below runs ~TaskGroup (cancel, taskWait, checkInterrupt) when the
        // closure returns.
        task_group(|tg| {
            // C++ line 102: HIPO_CLOCK_CREATE; (no-op at timing level 0)

            // C++ line 104: const bool parallel = S_.parTree();
            // C++ line 105: const bool serial = !parallel;
            let parallel = self.S_.parTree();
            let serial = !parallel;

            // C++ line 107: if (flag_stop_.load(std::memory_order_relaxed)) return;
            if self.flag_stop_.load(Ordering::Relaxed) {
                return;
            }
            // rusolve: stop at the HiGHS time limit (no change to the arithmetic before it)
            if let Some(d) = self.deadline_ {
                if std::time::Instant::now() > d {
                    self.timed_out_.store(true, Ordering::Relaxed);
                    self.flag_stop_.store(true, Ordering::Relaxed);
                    return;
                }
            }

            // C++ lines 109-120: spawn children of this supernode in forward order.
            if parallel {
                let mut child_to_spawn = self.first_child_[sn as usize];
                while child_to_spawn != -1 {
                    let child = child_to_spawn;
                    tg.spawn(move || {
                        let w = this;
                        // SAFETY: the parent's task group waits for (`sync` / `task_wait`)
                        // every child before it returns, so the Factorise outlives the
                        // task. Children of one supernode work on disjoint supernode
                        // indices (`sn`), and the parent reads a child's contribution only
                        // after the `sync` that matches the spawn (C++ 179-181).
                        unsafe { (*w.get()).processSupernode(child) };
                    });
                    // C++ line 114
                    child_to_spawn = self.next_child_[child_to_spawn as usize];
                }

                // C++ line 119: wait for first child to finish, before starting the parent
                if self.first_child_[sn as usize] != -1 {
                    tg.sync();
                }
            }

            // ===================================================
            // Supernode information
            // ===================================================
            // C++ line 125: HIPO_CLOCK_START(2) -- no-op at timing level 0.
            // first and last+1 column of the supernodes
            let sn_begin = self.S_.snStart(sn);
            let sn_end = self.S_.snStart(sn + 1);
            let sn_size = sn_end - sn_begin;

            // C++ lines 131-139: when the tree is processed in serial, use CliqueStack to
            // store the cliques. Otherwise, use local storage in FormatHandler.
            let mut clique_ptr: Option<Int64> = None;
            if serial {
                let stack = self
                    .stack_
                    .as_deref_mut()
                    .expect("serial path requires stack");
                let mut reallocation = false;
                clique_ptr = stack.setup(self.S_.cliqueSize(sn), &mut reallocation);
                if reallocation {
                    if let Some(lg) = self.logger_ {
                        lg.printInfo("Reallocation of CliqueStack\n");
                    }
                }
            }

            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(16);
            // initialise the format handler
            // this also allocates space for the frontal matrix and schur complement
            // C++ lines 143-144
            let mut FH = HybridHybridFormatHandler::new(
                self.S_,
                sn,
                &self.regul_,
                self.data_,
                &mut self.sn_columns_[sn as usize],
                clique_ptr,
            );

            drop(_pf);
            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(17);
            let t_asm = crate::ipm::hipo::factorhighs::call_and_time_blas::bucket_on()
                .then(std::time::Instant::now);

            // ===================================================
            // Assemble original matrix M into frontal
            // ===================================================
            // C++ line 151: HIPO_CLOCK_START(2) -- no-op.
            // j is relative column index in the frontal matrix
            for j in 0..sn_size {
                // column index in the original matrix
                let col = sn_begin + j;

                // go through the column
                for el in self.ptrM_[col as usize]..self.ptrM_[col as usize + 1] {
                    // relative row index in the frontal matrix
                    let i = self.S_.relindCols(el);

                    let val = self.valM_[el as usize];
                    FH.assembleFrontal(i, j, val);
                }
            }
            // C++ line 165
            // (assemble children cliques included until the dense phase)

            // ===================================================
            // Assemble frontal matrices of children
            // ===================================================
            let mut child_sn = self.first_child_reverse_[sn as usize];
            while child_sn != -1 {
                // Child contribution is found:
                // - in cliquestack, if we are processing the tree in serial.
                // - in schur_contribution_ if we are processing the tree in parallel.
                // Children are always summed from last to first.

                // C++ lines 177-197
                let child_off: Int64;
                if parallel {
                    // sync with spawned child, apart from the first one
                    if child_sn != self.first_child_reverse_[sn as usize] {
                        tg.sync();
                    }

                    if self.flag_stop_.load(Ordering::Relaxed) {
                        return;
                    }

                    // child_clique = schur_contribution_[child_sn].data(); a C++ null
                    // pointer is an empty vector here.
                    if self.schur_contribution_[child_sn as usize].is_empty() {
                        if let Some(lg) = self.logger_ {
                            lg.printInfo("Missing child supernode contribution\n");
                        }
                        self.flag_stop_.store(true, Ordering::Relaxed);
                        return;
                    }
                    child_off = 0;
                } else {
                    let stack = self
                        .stack_
                        .as_deref_mut()
                        .expect("serial path requires stack");
                    let mut child: Int = 0;
                    child_off = stack.getChild(&mut child);
                    debug_assert!(child == child_sn);
                }

                // determine size of clique of child
                // C++ lines 200-201
                let child_begin = self.S_.snStart(child_sn);
                let child_end = self.S_.snStart(child_sn + 1);

                // number of nodes in child sn
                // C++ line 204
                let child_size = child_end - child_begin;

                // size of clique of child sn
                // C++ line 207
                let nc = (self.S_.ptr(child_sn + 1) - self.S_.ptr(child_sn)) as Int - child_size;

                // ASSEMBLE INTO FRONTAL
                // C++ line 210: HIPO_CLOCK_START(2) -- no-op.
                // go through the columns of the contribution of the child
                let mut col: Int = 0;
                while col < nc {
                    // relative index of column in the frontal matrix
                    let j = self.S_.relindClique(child_sn, col);

                    if j < sn_size {
                        // assemble into frontal

                        // go through the rows of the contribution of the child
                        let mut row: Int = col;
                        while row < nc {
                            // relative index of the entry in the matrix frontal
                            let i = self.S_.relindClique(child_sn, row);

                            // how many entries to sum
                            let mut consecutive = self.S_.consecutiveSums(child_sn, row);

                            let child_slice: &[f64] = if parallel {
                                &self.schur_contribution_[child_sn as usize]
                            } else {
                                let stack =
                                    self.stack_.as_deref().expect("serial path requires stack");
                                &stack.stack_[child_off as usize..]
                            };
                            FH.assembleFrontalMultiple(
                                &mut consecutive,
                                child_slice,
                                nc,
                                child_sn,
                                row,
                                col,
                                i,
                                j,
                            );

                            row += consecutive;
                        }
                    } else {
                        break;
                    }
                    col += 1;
                }
                // C++ line 236: HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleChildrenFrontal) -- no-op.

                // ASSEMBLE INTO CLIQUE
                // C++ line 239: HIPO_CLOCK_START(2) -- no-op.
                if parallel {
                    FH.assembleClique(
                        &mut self.schur_contribution_[child_sn as usize],
                        0,
                        nc,
                        child_sn,
                    );
                } else {
                    let stack = self
                        .stack_
                        .as_deref_mut()
                        .expect("serial path requires stack");
                    FH.assembleClique(&mut stack.stack_, child_off, nc, child_sn);
                }
                // C++ line 241: HIPO_CLOCK_STOP(2, data_, kTimeFactoriseAssembleChildrenClique) -- no-op.

                // Schur contribution of the child is no longer needed
                // C++ lines 244-248
                if parallel {
                    freeVector(&mut self.schur_contribution_[child_sn as usize]);
                } else {
                    let stack = self
                        .stack_
                        .as_deref_mut()
                        .expect("serial path requires stack");
                    stack.popChild();
                }

                // move on to the next child
                // C++ line 251
                child_sn = self.next_child_reverse_[child_sn as usize];
            }

            // C++ line 254
            if self.flag_stop_.load(Ordering::Relaxed) {
                return;
            }

            // ===================================================
            // Partial factorisation
            // ===================================================
            // C++ line 259: HIPO_CLOCK_START(2) -- no-op.
            // threshold for regularisation
            // C++ line 261: const double reg_thresh = max_diag_ * kDynamicDiagCoeff;
            // C++ line 262
            let reg_thresh = self.M_norm1_ * kDynamicDiagCoeff;

            // C++ line 264: if (Int flag = FH->denseFactorise(reg_thresh))
            let asm_ns = t_asm.map_or(0, |t| t.elapsed().as_nanos() as u64);
            drop(_pf);
            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(18);
            let t_dense = crate::ipm::hipo::factorhighs::call_and_time_blas::bucket_on()
                .then(std::time::Instant::now);
            let flag: Int = if parallel {
                // the clique lives in the FormatHandler (clique_ptr_ is None)
                FH.denseFactorise(&mut [], reg_thresh)
            } else {
                let stack = self
                    .stack_
                    .as_deref_mut()
                    .expect("serial path requires stack");
                FH.denseFactorise(&mut stack.stack_, reg_thresh)
            };

            if let Some(t) = t_dense {
                crate::ipm::hipo::factorhighs::call_and_time_blas::bucket_add(
                    sn_size as i64,
                    (self.S_.ptr(sn + 1) - self.S_.ptr(sn)) - sn_size as i64,
                    asm_ns,
                    t.elapsed().as_nanos() as u64,
                );
            }
            if flag != 0 {
                // C++ line 265
                self.flag_stop_.store(true, Ordering::Relaxed);

                // C++ lines 267-270
                if flag == kRetInvalidInput {
                    if let Some(lg) = self.logger_ {
                        lg.printInfo("DenseFact: invalid input\n");
                    }
                } else if flag == kRetInvalidPivot {
                    if let Some(lg) = self.logger_ {
                        lg.printInfo("DenseFact: invalid pivot\n");
                    }
                }

                return;
            }
            drop(_pf);
            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(19);

            // C++ line 276: HIPO_CLOCK_START(2) -- no-op.
            // compute largest elements in factorisation
            // C++ line 278
            FH.extremeEntries();

            // terminate the format handler
            // C++ lines 281-282
            FH.base_mut().terminate(
                &mut self.schur_contribution_[sn as usize],
                &mut self.total_reg_,
                &mut self.swaps_[sn as usize],
                &mut self.pivot_2x2_[sn as usize],
            );

            drop(_pf);
            let _pf = crate::ipm::hipo::factorhighs::call_and_time_blas::prof(20);
            // C++ line 284
            if serial {
                let stack = self
                    .stack_
                    .as_deref_mut()
                    .expect("serial path requires stack");
                stack.pushWork(sn);
            }

            // C++ line 286: HIPO_CLOCK_STOP(2, data_, kTimeFactoriseTerminate) -- no-op.
        });
    }

    /// HiGHS: Factorise.cpp::Factorise::run (lines 289-341)
    /// Returns true on failure (C++ convention). Sets `num.factorised_`.
    pub fn run(&mut self, num: &mut Numeric) -> bool {
        // Root tasks call back into this object through `this` (the C++ lambda
        // captures `this`).
        let this = SendFactorise {
            ptr: self as *mut Factorise<'a>,
        };

        // C++ line 292: highs::parallel::TaskGroup tg; it lives for the whole body
        // (destroyed at function exit).
        task_group(|tg| {
            // C++ line 290: HIPO_CLOCK_CREATE; -- no-op.

            // C++ line 294: total_reg_.assign(n_, 0.0);
            self.total_reg_.clear();
            self.total_reg_.resize(self.n_ as usize, 0.0);

            // allocate space
            // C++ lines 297-299
            self.schur_contribution_
                .resize(self.S_.sn() as usize, Vec::new());
            self.swaps_.resize(self.S_.sn() as usize, Vec::new());
            self.pivot_2x2_.resize(self.S_.sn() as usize, Vec::new());

            // This should actually allocate only the first time, then sn_columns_ reuses
            // the memory of previous factorisations.
            // C++ line 303
            self.sn_columns_.resize(self.S_.sn() as usize, Vec::new());

            if self.S_.parTree() {
                // C++ line 306: Int spawned_roots{}; only counted, never read.
                // spawn tasks for root supernodes
                // C++ lines 308-313
                for sn in 0..self.S_.sn() {
                    if self.S_.snParent(sn) == -1 {
                        tg.spawn(move || {
                            let w = this;
                            // SAFETY: the task_wait below (and ~TaskGroup) keeps the
                            // Factorise alive until every root task has finished; root
                            // subtrees are disjoint sets of supernodes, so the tasks
                            // write disjoint slots (C++ FormatHandler::terminate note).
                            unsafe { (*w.get()).processSupernode(sn) };
                        });
                        // C++ line 311: ++spawned_roots;
                    }
                }

                // sync tasks for root supernodes
                // C++ line 316
                tg.task_wait();
            } else {
                // processing the tree in serial requires a CliqueStack
                // C++ lines 319-320
                if self.stack_.is_none() {
                    return true;
                }
                {
                    let stack = self
                        .stack_
                        .as_deref_mut()
                        .expect("serial path requires stack");
                    if stack.empty() {
                        stack.init(self.S_.maxStackSize());
                    }
                }

                // go through each supernode serially
                // C++ lines 323-325
                for sn in 0..self.S_.sn() {
                    self.processSupernode(sn);
                }
            }

            // C++ line 328
            if self.flag_stop_.load(Ordering::Relaxed) {
                return true;
            }

            // move factorisation to numerical object
            // C++ lines 331-336: num.S_ = &S_; num.sn_columns_ = &sn_columns_; and
            // num.data_ = &data_; have no Rust counterpart (Numeric has no such fields).
            num.factorised_ = true;
            num.total_reg_ = std::mem::take(&mut self.total_reg_);
            num.swaps_ = std::mem::take(&mut self.swaps_);
            num.pivot_2x2_ = std::mem::take(&mut self.pivot_2x2_);

            // C++ line 338: HIPO_CLOCK_STOP(1, data_, kTimeFactorise) -- no-op.

            false
        })
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::super::super::Int;
    use super::super::analyse::Analyse;
    use super::super::clique_stack::CliqueStack;
    use super::super::data_collector::DataCollector;
    use super::super::numeric::Numeric;
    use super::super::symbolic::Symbolic;
    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    /// Lower-triangular CSC "reversed arrowhead" of size `n` (SPD): the last row is dense,
    /// so every column `c < n - 1` has the single child-to-parent link `c -> n - 1` and the
    /// elimination tree is a star with `n - 1` leaves under the root.
    fn arrow(n: Int) -> (Vec<Int>, Vec<Int>, Vec<f64>) {
        let mut rows = Vec::new();
        let mut vals = Vec::new();
        let mut ptr = vec![0 as Int];
        for c in 0..n - 1 {
            rows.push(c);
            vals.push(4.0);
            rows.push(n - 1);
            vals.push(-1.0);
            ptr.push(rows.len() as Int);
        }
        rows.push(n - 1);
        vals.push(n as f64 + 1.0);
        ptr.push(rows.len() as Int);
        (rows, ptr, vals)
    }

    /// Bit patterns of the factor outputs: (total_reg, swaps, pivot_2x2, sn_columns).
    type Factors = (Vec<u64>, Vec<Int>, Vec<u64>, Vec<u64>);

    /// Analyse and factorise the arrowhead matrix, with `par_tree` as `S_.parTree()`.
    /// With `stopped`, `flag_stop_` is already raised when `run` starts (as after a failed
    /// dense factorisation elsewhere in the tree).
    /// Returns (`run` failed, `num.factorised_`, factors).
    fn factorise_arrow(par_tree: bool, stopped: bool) -> (bool, bool, Factors) {
        let n: Int = 40;
        let (rows, ptr, vals) = arrow(n);
        let signs = vec![1 as Int; n as usize];
        let perm: Vec<Int> = (0..n).collect();
        let data = DataCollector::new();
        let mut S = Symbolic::new();
        let mut an = Analyse::new(&rows, &ptr, &signs, 4, None, &data, &perm);
        assert_eq!(an.run(&mut S), 0);
        S.setParallel(par_tree, false);
        // The tree must have a supernode with at least two children, so that the
        // last-to-first child sum is exercised.
        let mut kids = vec![0 as Int; S.sn() as usize];
        for &p in S.snParent_vec() {
            if p != -1 {
                kids[p as usize] += 1;
            }
        }
        assert!(
            kids.iter().any(|&k| k >= 2),
            "no supernode with two children"
        );

        let regul = Regul::default();
        let mut sn_columns: Vec<Vec<f64>> = Vec::new();
        let mut stack = CliqueStack::default();
        let mut num = Numeric {
            factorised_: false,
            swaps_: Vec::new(),
            pivot_2x2_: Vec::new(),
            total_reg_: Vec::new(),
        };
        let failed;
        {
            let mut fact = Factorise::new(
                &S,
                &rows,
                &ptr,
                &vals,
                &regul,
                None,
                &data,
                &mut sn_columns,
                if par_tree { None } else { Some(&mut stack) },
            );
            if stopped {
                fact.flag_stop_.store(true, Ordering::Relaxed);
            }
            failed = fact.run(&mut num);
        }
        let factors = (
            num.total_reg_.iter().map(|x| x.to_bits()).collect(),
            num.swaps_.iter().flatten().copied().collect(),
            num.pivot_2x2_
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect(),
            sn_columns.iter().flatten().map(|x| x.to_bits()).collect(),
        );
        (failed, num.factorised_, factors)
    }

    /// Run `f` on a fresh std::thread with its own `threads`-worker scheduler, shut down
    /// afterwards. `threads` = 1 also hosts the
    /// serial runs, whose `TaskGroup`s need a worker deque.
    fn on_fresh_thread<R: Send + 'static>(
        threads: i32,
        f: impl FnOnce() -> R + Send + 'static,
    ) -> R {
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            let r = f();
            HighsTaskExecutor::shutdown(true);
            r
        })
        .join()
        .expect("factorisation thread panicked")
    }

    #[test]
    fn w5_hipo_tree_par_tree_spawn_and_sum_order() {
        // Serial reference, then the tree in parallel at 1, 2 and 4 threads. Each thread
        // count runs on a fresh std::thread with its own scheduler, and the results must
        // be bit-identical to the serial run.
        let reference = on_fresh_thread(1, || factorise_arrow(false, false));
        assert!(!reference.0);
        assert!(reference.1);
        for threads in [1, 2, 4] {
            let got = on_fresh_thread(threads, || {
                (factorise_arrow(true, false), factorise_arrow(true, false))
            });
            assert_eq!(got.0, reference, "threads = {threads}");
            assert_eq!(got.1, reference, "threads = {threads} (second run)");
        }
    }

    #[test]
    fn w5_hipo_tree_par_tree_flag_stop() {
        // C++ lines 107 and 328: with flag_stop_ raised, every processSupernode (each root
        // task in the parallel tree) returns at once and run reports failure without
        // moving the factor into Numeric, serially and at 1, 2 and 4 threads. (A NaN pivot
        // cannot raise it here: blockBunchKaufman never selects a NaN diagonal and then
        // swaps with index -1, out of bounds in the C++ too.)
        let serial = on_fresh_thread(1, || factorise_arrow(false, true));
        assert!(serial.0);
        assert!(!serial.1);
        for threads in [1, 2, 4] {
            let got = on_fresh_thread(threads, || factorise_arrow(true, true));
            assert!(got.0, "threads = {threads}");
            assert!(!got.1, "threads = {threads}");
        }
    }
}
