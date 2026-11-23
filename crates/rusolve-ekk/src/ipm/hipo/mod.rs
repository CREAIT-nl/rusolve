//! HiGHS source: ipm/hipo/ (HiPO interior point solver), extern/{metis,amd,rcm} (orderings)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # HiPO: the HiGHS 1.15.1 interior point solver with a direct factorisation
//!
//! HiPO solves LPs (and convex QPs) with a primal-dual interior point method
//! whose Newton systems are factorised by FactorHiGHS, a supernodal
//! multifrontal LDL^T solver. HiGHS reaches it from `solveLp` / `Highs::solve`
//! through `solveLpHipo` / `solveHipo` (`ipm/IpxWrapper.cpp`, Rust
//! [`crate::ipm::hipo_wrapper`]). The availability switch is the non-upstream
//! option field `HighsOptions::hipo_available` (default `false`): with the
//! switch off every dispatch site behaves exactly as without HiPO.
//!
//! | Rust module | HiGHS source | reached with |
//! |-------------|--------------|--------------|
//! | [`auxiliary`] | `ipm/hipo/auxiliary/{Auxiliary,Logger,VectorOperations,IntConfig}` | everything |
//! | [`blas`] | `extern/blas/mycblas.h` + `tools/highs-oracle/refblas/refblas.c` (netlib loop order) | `factorhighs::call_and_time_blas` only |
//! | [`factorhighs`] | `ipm/hipo/factorhighs/*` | `ipm::factor_highs_solver` |
//! | [`ipm`] | `ipm/hipo/ipm/*` (`Solver`, `Model`, `Iterate`, `KktMatrix`, ...) | `hipo_wrapper::solve_hipo` |
//! | [`orderings`] | `extern/HighsExtrasExternalDeps.h` (amd/metis/rcm API), `extern/metis`, `extern/amd`, `extern/rcm` | `FactorHiGHSSolver::chooseOrdering` |
//! | [`oracle_trace`] | none (oracle hooks, `tools/highs-oracle/FORMAT.md`) | Solver / FactorHiGHSSolver / solveHipo |
//!
//! `KrylovMethods.cpp` and `KrylovMethodsIpm.cpp` are not ported: no HiPO
//! call site reaches them in 1.15.1. Parallel paths (`parallel_tree_`,
//! `parallel_node_`, `highs::parallel::for_each`, `spawn`/`sync`) are never
//! taken in a serial run: `setParallel` always resolves to serial because the
//! Rust run has one thread.
//!
//! # Naming
//!
//! HiPO is C++ with camelCase methods; like HiPDLP its Rust keeps the
//! exact C++ names of types, fields, methods and free functions
//! (`Solver::startingPoint`, `Iterate::residual1234`, `Analyse::eTree`).
//! METIS / AMD / RCM keep their C names (`MlevelNestedDissection`,
//! `Highs_amd_2`, `Highs_genrcm`). HiGHS-side wrappers (`IpxWrapper.cpp`)
//! follow the crate snake_case rule (`solve_lp_hipo_so`, `solve_hipo`).
//! Overloads take a suffix named in their doc line (`transpose_val`,
//! `norm2_2`, `Symbolic::ptr_vec`). `Int` = [`HighsInt`] (`IntConfig.h`),
//! `Int64` = `i64`.
//!
//! # Virtual dispatch (the one mapping used throughout)
//!
//! An abstract C++ class that carries data becomes a *base struct* holding the
//! data plus a *trait* for the virtual methods; implementors embed the base
//! and expose it through `base()` / `base_mut()`; `std::unique_ptr<Base>` is
//! `Box<dyn Trait + 'a>`:
//!
//! | C++ | Rust base struct | Rust trait | implementors |
//! |-----|------------------|------------|--------------|
//! | `FormatHandler` | `factorhighs::format_handler::FormatHandler<'a>` | `FormatHandlerOps` | `HybridHybridFormatHandler<'a>` |
//! | `SolveHandler` | `factorhighs::solve_handler::SolveHandler<'a>` | `SolveHandlerOps` | `HybridSolveHandler<'a>` |
//! | `LinearSolver` | (flag `valid_` only, via `valid()`/`set_valid()`) | `ipm::linear_solver::LinearSolver` | `FactorHiGHSSolver` |
//! | `PreprocessAction` | `ipm::pre_process::PreprocessActionBase` | `PreprocessAction` | `PreprocessEmptyRows`, `PreprocessFixedVars`, `PreprocessScaling`, `PreprocessFormulation`, `PreprocessFreeVars` |
//!
//! A pure virtual with a body-less default in C++ (`FormatHandler::extremeEntries() {}`,
//! `LinearSolver::setup() { return 0; }`, `flops()/spops()/nz() { return 0; }`,
//! `getReg(...) {}`) is a provided trait method with the same body.
//!
//! # Ownership (C++ references and pointers)
//!
//! * `Solver` owns everything (`model_`, `LS_`, `kkt_`, `it_`, `regul_`,
//!   `info_`, `control_`, `options_`, `logger_`, `ipx_lps_`). The references
//!   that `FactorHiGHSSolver` and `KktMatrix` cache in C++ (`kkt_`, `model_`,
//!   `regul_`, `info_`, `data_` = `it_->data`, `logger_`, `options_`) are passed
//!   per call: `KktMatrix` methods take `(model, regul, info, logger)`;
//!   `LinearSolver` methods take an [`ipm::linear_solver::LsCtx`].
//! * `Iterate` drops `model` / `regul`: methods that read them take
//!   `model: &Model` / `regul: &Regularisation`; `Rp` / `Rd` (pointers into
//!   `total_reg`) are `Option<usize>` offsets.
//! * `Model::lp_orig_` (a `const HighsLp*` used only by `loadIntoIpx` and the
//!   KKT check in `solveHipo`) is an owned clone `Option<HighsLp>`.
//! * `Control` owns the callback (moved in by `setCallback`, moved back by
//!   `take_callback`) and a clone of the timer; `interruptCheck` takes the
//!   `Options` it reads.
//! * `Logger` holds a clone `Option<HighsLogOptions>`; its printf-style
//!   templates take an already formatted `&str`.
//! * FactorHiGHS: `Factorise` / `FormatHandler` / `SolveHandler` / `Analyse`
//!   borrow with lifetimes (`&'a Symbolic`, `&'a mut Vec<f64>`, `&'a DataCollector`),
//!   which C++ guarantees by construction. Pointers into the serial
//!   `CliqueStack` buffer are `Int64` offsets; methods that read or write the
//!   stack take `stack: &mut [f64]`.
//! * METIS: graphs live in a [`orderings::metis::GraphPool`] arena addressed
//!   by `GraphId`; `coarser` / `finer` are `Option<GraphId>`.
//!
//! Nothing here is in the crate prelude.

pub mod auxiliary;
pub mod blas;
pub mod factorhighs;
pub mod ipm;
pub mod oracle_trace;
pub mod orderings;

use crate::HighsInt;

/// HiGHS: ipm/hipo/auxiliary/IntConfig.h `typedef HighsInt Int`
pub type Int = HighsInt;
/// HiGHS: ipm/hipo/auxiliary/IntConfig.h `typedef int64_t Int64`
pub type Int64 = i64;
