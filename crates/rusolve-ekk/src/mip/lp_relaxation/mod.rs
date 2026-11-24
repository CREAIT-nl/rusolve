//! HiGHS source: mip/HighsLpRelaxation.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use std::sync::Arc;

mod local_bigm;
mod node_basis;
pub use local_bigm::*;
mod part1;
pub use part1::*;
mod part2;
pub use part2::*;
mod root_race;

/// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation (lines 23-376)
///
/// Removed C++ members:
/// * C++: `const HighsMipSolver & mipsolver`: removed; methods that use it take `ms`
///   (`&HighsMipSolver` for reads, `&mut HighsMipSolver` when the body reaches a pool,
///   a domain or the clique table mutably).
///
/// Ownership: an LP relaxation that lives in
/// `ms.mipdata_.lps[k]` is detached by its caller for every `&mut self` method that also takes
/// `ms`. Read-only methods that take `&HighsMipSolver` can be called on the attached value
/// (`ms.mipdata_.lps[k].get_row(ms, ..)`). Methods that re-enter through `ms`
/// (`resolve_lp` / `run` -> `addIncumbent` / `trySolution`, which print `getLp().numRows()`)
/// lend `self` back with `lend_in!(ms.mipdata_.lps, self, ..)` around that call.
#[derive(Default)]
pub struct HighsLpRelaxation {
    /// C++: `Highs lpsolver`
    pub lpsolver: Highs,
    /// rusolve mip_node_bigm: model-row coefficients currently tightened for the node subtree.
    pub local_mods: Vec<HighsLpRelaxation_LocalMod>,
    /// rusolve mip_node_bigm_rows: number of node-local rows at the tail of the LP.
    pub nlocal_rows: usize,
    /// C++: `std::vector<LpRow> lprows`
    pub lprows: Vec<HighsLpRelaxation_LpRow>,
    /// C++: `std::vector<std::pair<HighsInt, double>> fractionalints`
    pub fractionalints: Vec<(HighsInt, f64)>,
    /// C++: `std::vector<double> dualproofvals`
    pub dualproofvals: Vec<f64>,
    /// C++: `std::vector<HighsInt> dualproofinds`
    pub dualproofinds: Vec<HighsInt>,
    /// C++: `std::vector<double> dualproofbuffer`
    pub dualproofbuffer: Vec<f64>,
    /// C++: `std::vector<double> colLbBuffer`
    pub col_lb_buffer: Vec<f64>,
    /// C++: `std::vector<double> colUbBuffer`
    pub col_ub_buffer: Vec<f64>,
    /// C++: `HVector row_ep`
    pub row_ep: HVector,
    /// C++: `HighsSparseVectorSum row_ap`
    pub row_ap: HighsSparseVectorSum,
    /// C++: `double dualproofrhs`
    pub dualproofrhs: f64,
    /// C++: `bool hasdualproof`
    pub hasdualproof: bool,
    /// C++: `double objective`
    pub objective: f64,
    /// C++: `std::shared_ptr<const HighsBasis> basischeckpoint`
    pub basischeckpoint: Option<Arc<HighsBasis>>,
    /// C++: `bool currentbasisstored`
    pub currentbasisstored: bool,
    /// C++: `int64_t numlpiters`
    pub numlpiters: i64,
    /// C++: `int64_t lastAgeCall`
    pub last_age_call: i64,
    /// C++: `double avgSolveIters`
    pub avg_solve_iters: f64,
    /// C++: `int64_t numSolved`
    pub num_solved: i64,
    /// C++: `size_t epochs`
    pub epochs: usize,
    /// C++: `HighsInt maxNumFractional`
    pub max_num_fractional: HighsInt,
    /// C++: `Status status`
    pub status: HighsLpRelaxation_Status,
    /// C++: `bool adjustSymBranchingCol`
    pub adjust_sym_branching_col: bool,
    /// C++: `bool solved_first_lp`
    pub solved_first_lp: bool,
    /// C++ `HighsMipWorker* worker_`: index into `ms.mipdata_.workers` (`None` = nullptr).
    pub worker_: Option<usize>,
}

/// C++ copy constructor; `Highs` itself is not `Clone`. Only for `#[derive(Clone)]` on
/// containers (`HighsMipSolverData`); the faithful call is
/// [`HighsLpRelaxation::copy_from`]`(ms, other)`. PORT NOTE: the C++ sizes the bound buffers
/// with `mipsolver.numCol()`; without `ms` this uses `other.col_lb_buffer.len()`, which is the
/// same once `loadModel` (or `copy_from`) ran on `other`.
impl Clone for HighsLpRelaxation {
    fn clone(&self) -> Self {
        Self::copy_from_num_col(self, self.col_lb_buffer.len())
    }
}

impl std::fmt::Debug for HighsLpRelaxation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsLpRelaxation").finish_non_exhaustive()
    }
}

/// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::Status (enum class)
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsLpRelaxation_Status {
    #[default]
    kNotSet = 0,
    kOptimal = 1,
    kInfeasible = 2,
    kUnscaledDualFeasible = 3,
    kUnscaledPrimalFeasible = 4,
    kUnscaledInfeasible = 5,
    kUnbounded = 6,
    kError = 7,
}

/// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::LpRow (lines 37-61)
#[derive(Clone, Copy, Default, Debug)]
pub struct HighsLpRelaxation_LpRow {
    /// C++: `Origin origin`
    pub origin: HighsLpRelaxation_LpRow_Origin,
    /// C++: `HighsInt index`
    pub index: HighsInt,
    /// C++: `HighsInt age`
    pub age: HighsInt,
    /// C++: `HighsInt cutpoolindex`
    pub cutpoolindex: HighsInt,
}

/// HiGHS: C++ plain `enum Origin` nested in `HighsLpRelaxation::LpRow` (values are associated consts).
pub type HighsLpRelaxation_LpRow_Origin = HighsInt;
impl HighsLpRelaxation_LpRow {
    pub const kModel: HighsInt = 0;
    pub const kCutPool: HighsInt = 1;
    /// rusolve mip_node_bigm_rows: node-local row, never in a pool.
    pub const kLocal: HighsInt = 2;
}

/// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::Playground (lines 110-153)
///
/// Removed C++ members:
/// * C++: `HighsLpRelaxation * lp`: removed; methods that use it take `lp: &mut HighsLpRelaxation`
///   (the caller passes the LP the playground came from).
///
/// PORT NOTE: the C++ destructor (re-solve after a stored iterate) is
/// [`HighsLpRelaxation_Playground::release`]; Rust has no destructor with access to `ms`, so
/// the owner must call `release(ms, lp)` exactly where the C++ object goes out of scope.
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsLpRelaxation_Playground {
    /// C++: `bool iterateStored`
    pub iterate_stored: bool,
}

impl HighsLpRelaxation_LpRow {
    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::LpRow::cut (lines 57-59)
    /// C++: `static LpRow cut(HighsInt index, HighsInt cutpoolindex)`
    pub fn cut(index: HighsInt, cutpoolindex: HighsInt) -> HighsLpRelaxation_LpRow {
        HighsLpRelaxation_LpRow {
            origin: HighsLpRelaxation_LpRow::kCutPool,
            index,
            age: 0,
            cutpoolindex,
        }
    }

    /// rusolve mip_node_bigm_rows: a node-local row.
    /// `index` is the model row it was derived from (length, max coefficient stay readable).
    pub fn local(index: HighsInt) -> HighsLpRelaxation_LpRow {
        HighsLpRelaxation_LpRow {
            origin: HighsLpRelaxation_LpRow::kLocal,
            index,
            age: 0,
            cutpoolindex: -1,
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::LpRow::model (lines 60-60)
    /// C++: `static LpRow model(HighsInt index)`
    pub fn model(index: HighsInt) -> HighsLpRelaxation_LpRow {
        HighsLpRelaxation_LpRow {
            origin: HighsLpRelaxation_LpRow::kModel,
            index,
            age: 0,
            cutpoolindex: -1,
        }
    }
}

impl HighsLpRelaxation_Playground {
    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::Playground::Playground (lines 115-115)
    /// C++: `Playground(HighsLpRelaxation * lp)`
    // PORT NOTE: the `lp` pointer is not stored (see the struct doc).
    pub fn new() -> Self {
        HighsLpRelaxation_Playground {
            iterate_stored: false,
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::Playground::operator= (lines 123-127)
    /// C++: `Playground & operator=(Playground && other)`
    // PORT NOTE: swaps `iterateStored` only; the `lp` pointers are the callers' (both
    // playgrounds must belong to the same LP, as at every C++ call site).
    pub fn op_assign(
        &mut self,
        other: &mut HighsLpRelaxation_Playground,
    ) -> &mut HighsLpRelaxation_Playground {
        std::mem::swap(&mut self.iterate_stored, &mut other.iterate_stored);
        self
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::Playground::solveLp (lines 129-141)
    /// C++: `HighsLpRelaxation::Status solveLp(HighsDomain & localdom)`
    pub fn solve_lp(
        &mut self,
        ms: &mut HighsMipSolver,
        lp: &mut HighsLpRelaxation,
        localdom: &mut HighsDomain,
    ) -> HighsLpRelaxation_Status {
        if self.iterate_stored {
            lp.flush_domain(ms, localdom, false);
            lp.get_lp_solver().get_iterate();
        } else {
            debug_assert!(lp.get_lp_solver().get_info().valid);
            lp.get_lp_solver().put_iterate();
            lp.flush_domain(ms, localdom, false);
            self.iterate_stored = true;
        }

        lp.run(ms, false)
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::Playground::~Playground (lines 146-152)
    /// C++: `~Playground()`
    pub fn release(&mut self, ms: &mut HighsMipSolver, lp: &mut HighsLpRelaxation) {
        if self.iterate_stored {
            lp.get_lp_solver().get_iterate();
            lp.run(ms, true);
            // If desired, here is the place to clear the stored iterate
        }
    }
}

impl HighsLpRelaxation {
    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::playground (lines 155-155)
    /// C++: `Playground playground()`
    pub fn playground(&mut self) -> HighsLpRelaxation_Playground {
        HighsLpRelaxation_Playground::new()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getRow (lines 159-166)
    /// C++: `void getRow(HighsInt row, HighsInt & len, const HighsInt *& inds, const double *& vals) const`
    // PORT NOTE: the row slices borrow `ms` (a cut row lives in `ms.mipdata_.cutpools[..]`, a
    // model row in `ms.mipdata_.ARindex_`), so the pool of a cut row must be attached.
    pub fn get_row<'s>(
        &self,
        ms: &'s HighsMipSolver,
        row: HighsInt,
        len: &mut HighsInt,
        inds: &mut &'s [HighsInt],
        vals: &mut &'s [f64],
    ) {
        if row < ms.num_row() {
            debug_assert!(self.lprows[row as usize].origin == HighsLpRelaxation_LpRow::kModel);
        } else {
            debug_assert!(self.lprows[row as usize].origin == HighsLpRelaxation_LpRow::kCutPool);
        }
        self.lprows[row as usize].get(ms, len, inds, vals);
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::isRowIntegral (lines 168-171)
    /// C++: `bool isRowIntegral(HighsInt row) const`
    // PORT NOTE (perf): `#[inline]`, as the C++ `isRowIntegral`.
    #[inline]
    pub fn is_row_integral(&self, ms: &HighsMipSolver, row: HighsInt) -> bool {
        debug_assert!(row < self.lprows.len() as HighsInt);
        self.lprows[row as usize].is_integral(ms)
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::setAdjustSymmetricBranchingCol (lines 173-175)
    /// C++: `void setAdjustSymmetricBranchingCol(bool adjustSymBranchingCol)`
    pub fn set_adjust_symmetric_branching_col(&mut self, adjust_sym_branching_col: bool) {
        self.adjust_sym_branching_col = adjust_sym_branching_col;
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getAvgSolveIters (lines 185-185)
    /// C++: `double getAvgSolveIters()`
    pub fn get_avg_solve_iters(&mut self) -> f64 {
        self.avg_solve_iters
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getRowLen (lines 187-189)
    /// C++: `HighsInt getRowLen(HighsInt row) const`
    pub fn get_row_len(&self, ms: &HighsMipSolver, row: HighsInt) -> HighsInt {
        self.lprows[row as usize].get_row_len(ms)
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getMaxAbsRowVal (lines 191-193)
    /// C++: `double getMaxAbsRowVal(HighsInt row) const`
    pub fn get_max_abs_row_val(&self, ms: &HighsMipSolver, row: HighsInt) -> f64 {
        self.lprows[row as usize].get_max_abs_val(ms)
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getLp (lines 195-195)
    /// C++: `const HighsLp & getLp() const`
    pub fn get_lp(&self) -> &HighsLp {
        self.lpsolver.get_lp()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getSolution (lines 197-197)
    /// C++: `const HighsSolution & getSolution() const`
    pub fn get_solution(&self) -> &HighsSolution {
        self.lpsolver.get_solution()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::rowLower (lines 203-205)
    /// C++: `double rowLower(HighsInt row) const`
    pub fn row_lower(&self, row: HighsInt) -> f64 {
        self.lpsolver.get_lp().row_lower_[row as usize]
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::rowUpper (lines 207-209)
    /// C++: `double rowUpper(HighsInt row) const`
    pub fn row_upper(&self, row: HighsInt) -> f64 {
        self.lpsolver.get_lp().row_upper_[row as usize]
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::colLower (lines 211-215)
    /// C++: `double colLower(HighsInt col, const HighsDomain & globaldom) const`
    // PORT NOTE: `globaldom` is an index into `ms.mipdata_.domains` (see `slack_lower`).
    pub fn col_lower(&self, ms: &mut HighsMipSolver, col: HighsInt, globaldom: usize) -> f64 {
        let num_col = self.lpsolver.get_lp().num_col_;
        if col < num_col {
            self.lpsolver.get_lp().col_lower_[col as usize]
        } else {
            self.slack_lower(ms, col - num_col, globaldom)
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::colUpper (lines 217-221)
    /// C++: `double colUpper(HighsInt col, const HighsDomain & globaldom) const`
    // PORT NOTE: `globaldom` is an index into `ms.mipdata_.domains` (see `slack_upper`).
    pub fn col_upper(&self, ms: &HighsMipSolver, col: HighsInt, globaldom: usize) -> f64 {
        let num_col = self.lpsolver.get_lp().num_col_;
        if col < num_col {
            self.lpsolver.get_lp().col_upper_[col as usize]
        } else {
            self.slack_upper(ms, col - num_col, globaldom)
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::isColIntegral (lines 223-227)
    /// C++: `bool isColIntegral(HighsInt col) const`
    // PORT NOTE (perf): `#[inline]`, as the C++ `isColIntegral` header
    // inline it translates: `HighsTransformedLp::transform` runs it once per nonzero.
    #[inline]
    pub fn is_col_integral(&self, ms: &HighsMipSolver, col: HighsInt) -> bool {
        let num_col = self.lpsolver.get_lp().num_col_;
        if col < num_col {
            ms.is_col_integral(col)
        } else {
            self.is_row_integral(ms, col - num_col)
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::solutionValue (lines 229-233)
    /// C++: `double solutionValue(HighsInt col) const`
    pub fn solution_value(&self, col: HighsInt) -> f64 {
        let num_col = self.lpsolver.get_lp().num_col_;
        if col < num_col {
            self.get_solution().col_value[col as usize]
        } else {
            self.get_solution().row_value[(col - num_col) as usize]
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getStatus (lines 235-235)
    /// C++: `Status getStatus() const`
    pub fn get_status(&self) -> HighsLpRelaxation_Status {
        self.status
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getSolverInfo (lines 237-237)
    /// C++: `const HighsInfo & getSolverInfo() const`
    pub fn get_solver_info(&self) -> &HighsInfo {
        self.lpsolver.get_info()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getNumLpIterations (lines 239-239)
    /// C++: `int64_t getNumLpIterations() const`
    pub fn get_num_lp_iterations(&self) -> i64 {
        self.numlpiters
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::integerFeasible (lines 241-248)
    /// C++: `bool integerFeasible() const`
    pub fn integer_feasible(&self) -> bool {
        if (self.status == HighsLpRelaxation_Status::kOptimal
            || self.status == HighsLpRelaxation_Status::kUnscaledPrimalFeasible)
            && self.fractionalints.is_empty()
        {
            return true;
        }

        false
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::setMipWorker (lines 250-250)
    /// C++: `void setMipWorker(HighsMipWorker & worker)`
    // PORT NOTE: the worker is its index into `ms.mipdata_.workers` (C++ `&worker`).
    pub fn set_mip_worker(&mut self, worker_: usize) {
        self.worker_ = Some(worker_);
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::scaledOptimal (lines 256-266)
    /// C++: `static bool scaledOptimal(Status status)`
    pub fn scaled_optimal(status: HighsLpRelaxation_Status) -> bool {
        matches!(
            status,
            HighsLpRelaxation_Status::kOptimal
                | HighsLpRelaxation_Status::kUnscaledDualFeasible
                | HighsLpRelaxation_Status::kUnscaledPrimalFeasible
                | HighsLpRelaxation_Status::kUnscaledInfeasible
        )
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::unscaledPrimalFeasible (lines 268-276)
    /// C++: `static bool unscaledPrimalFeasible(Status status)`
    pub fn unscaled_primal_feasible(status: HighsLpRelaxation_Status) -> bool {
        matches!(
            status,
            HighsLpRelaxation_Status::kOptimal | HighsLpRelaxation_Status::kUnscaledPrimalFeasible
        )
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::unscaledDualFeasible (lines 278-286)
    /// C++: `static bool unscaledDualFeasible(Status status)`
    pub fn unscaled_dual_feasible(status: HighsLpRelaxation_Status) -> bool {
        matches!(
            status,
            HighsLpRelaxation_Status::kOptimal | HighsLpRelaxation_Status::kUnscaledDualFeasible
        )
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::storeBasis (lines 292-297)
    /// C++: `void storeBasis()`
    pub fn store_basis(&mut self) {
        if !self.currentbasisstored && self.lpsolver.get_basis().valid {
            let mut basis = self.lpsolver.get_basis().clone();
            // rusolve lever `mip_lp_keep_dse_weights` = 2 (not in HiGHS): store the weights of
            // this basis with it, so that `recover_basis` restores them (simplex/ekk/dse_keep.rs).
            if self.lpsolver.options_.mip_lp_keep_dse_weights == 2 {
                basis.dse_weight = self
                    .lpsolver
                    .ekk_instance_
                    .dse_weights_by_variable()
                    .map(Arc::new);
            }
            // rusolve lever `mip_node_install_basis` (not in HiGHS): the cut rows of this
            // basis, so that a node installed later can match them to its LP (node_basis.rs).
            if self.lpsolver.options_.mip_node_install_basis {
                basis.row_keys = Some(Arc::new(self.cut_row_keys()));
            }
            self.basischeckpoint = Some(Arc::new(basis));
            self.currentbasisstored = true;
        }
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getStoredBasis (lines 299-301)
    /// C++: `std::shared_ptr<const HighsBasis> getStoredBasis() const`
    pub fn get_stored_basis(&self) -> Option<Arc<HighsBasis>> {
        self.basischeckpoint.clone()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::setStoredBasis (lines 303-306)
    /// C++: `void setStoredBasis(std::shared_ptr<const HighsBasis> basis)`
    pub fn set_stored_basis(&mut self, basis: Option<Arc<HighsBasis>>) {
        self.basischeckpoint = basis;
        self.currentbasisstored = false;
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getMipSolver (lines 308-308)
    /// C++: `const HighsMipSolver & getMipSolver() const`
    // PORT NOTE: identity on `ms` (M1: the LP stores no solver reference).
    pub fn get_mip_solver<'s>(&self, ms: &'s HighsMipSolver) -> &'s HighsMipSolver {
        ms
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getNumModelRows (lines 310-310)
    /// C++: `HighsInt getNumModelRows() const`
    pub fn get_num_model_rows(&self, ms: &HighsMipSolver) -> HighsInt {
        ms.num_row()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::numRows (lines 312-312)
    /// C++: `HighsInt numRows() const`
    pub fn num_rows(&self) -> HighsInt {
        self.lpsolver.get_num_row()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::numCols (lines 314-314)
    /// C++: `HighsInt numCols() const`
    pub fn num_cols(&self) -> HighsInt {
        self.lpsolver.get_num_col()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::numNonzeros (lines 316-316)
    /// C++: `HighsInt numNonzeros() const`
    pub fn num_nonzeros(&self) -> HighsInt {
        self.lpsolver.get_num_nz()
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getDualProof (lines 336-342)
    /// C++: `void getDualProof(const HighsInt *& inds, const double *& vals, double & rhs, HighsInt & len)`
    pub fn get_dual_proof<'s>(
        &'s mut self,
        inds: &mut &'s [HighsInt],
        vals: &mut &'s [f64],
        rhs: &mut f64,
        len: &mut HighsInt,
    ) {
        *inds = &self.dualproofinds;
        *vals = &self.dualproofvals;
        *rhs = self.dualproofrhs;
        *len = self.dualproofinds.len() as HighsInt;
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getLpSolver (lines 356-356)
    /// C++: `Highs & getLpSolver()`
    pub fn get_lp_solver(&mut self) -> &mut Highs {
        &mut self.lpsolver
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getLpSolver (lines 357-357)
    /// C++: `const Highs & getLpSolver() const`
    pub fn get_lp_solver_2(&self) -> &Highs {
        &self.lpsolver
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getFractionalIntegers (lines 359-362)
    /// C++: `const std::vector<std::pair<HighsInt, double>> & getFractionalIntegers() const`
    pub fn get_fractional_integers(&self) -> &[(HighsInt, f64)] {
        &self.fractionalints
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getFractionalIntegers (lines 364-366)
    /// C++: `std::vector<std::pair<HighsInt, double>> & getFractionalIntegers()`
    pub fn get_fractional_integers_2(&mut self) -> &mut Vec<(HighsInt, f64)> {
        &mut self.fractionalints
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::getObjective (lines 368-368)
    /// C++: `double getObjective() const`
    pub fn get_objective(&self) -> f64 {
        self.objective
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::setIterationLimit (lines 370-372)
    /// C++: `void setIterationLimit(HighsInt limit = kHighsIInf)`
    /// C++ default arguments: `limit = kHighsIInf`
    pub fn set_iteration_limit(&mut self, limit: HighsInt) {
        self.lpsolver
            .set_option_value_int("simplex_iteration_limit", limit);
    }

    /// HiGHS: HighsLpRelaxation.h::HighsLpRelaxation::setSolvedFirstLp (lines 373-375)
    /// C++: `void setSolvedFirstLp(const bool solved_first_lp_)`
    pub fn set_solved_first_lp(&mut self, solved_first_lp_: bool) {
        self.solved_first_lp = solved_first_lp_;
    }
}
