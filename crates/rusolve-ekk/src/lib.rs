//! HiGHS source: (crate root; no single source file)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # rusolve-ekk: a faithful port of the HiGHS simplex engine
//!
//! This crate is a *faithful translation* of the HiGHS dual simplex
//! ("EKK" = `HEkk`, `HEkkDual`, `HEkkDualRow`, `HEkkDualRHS`, `HEkkPrimal`),
//! its numerical linear algebra (`HSimplexNla`, `HFactor`) and the
//! supporting `util/` and `lp_data/` pieces, as a faithful port
//! (see AGENTS.md). Same algorithms, same constants, same control flow,
//! same data layout. It is *not* a re-derivation: when filling a body, read
//! the C++ function named in its `/// HiGHS: <File.cpp>::<Class>::<method>`
//! doc line and translate it statement by statement.
//!
//! # Module layout (mirrors HiGHS)
//!
//! | Rust module            | HiGHS directory / file                          |
//! |------------------------|-------------------------------------------------|
//! | [`io`]                 | `io/HighsIO.{h,cpp}` (minimal logging)          |
//! | [`lp_data`]            | `lp_data/` subset (HConst, HStruct, HighsLp, HighsOptions, HighsInfo, HighsStatus, ...) |
//! | [`util`]               | `util/` (HVector, HighsSparseMatrix, HFactor, HSet, HighsCDouble, HighsRandom, ...) |
//! | [`simplex`]            | `simplex/` (HEkk*, HSimplex*, SimplexConst, SimplexStruct) |
//! | [`ipm`]                | `ipm/` IPX, BASICLU, HiPO (`ipm/hipo`, `extern/{metis,amd,rcm}`) and `IpxWrapper.cpp` (not in the prelude) |
//! | [`pdlp`]               | `pdlp/` cuPDLP-C, HiPDLP and their wrappers (not in the prelude) |
//!
//! A `.cpp` longer than ~700 lines is split into a directory of Rust files,
//! each holding a separate `impl` block for a contiguous C++ line range.
//!
//! # Naming rules
//!
//! * Types keep their HiGHS names (`HEkk`, `HVector`, `HighsSparseMatrix`).
//! * Fields, methods, parameters and locals are `snake_case` conversions of
//!   the HiGHS names; trailing underscores are kept: `baseLower_` ->
//!   `base_lower_`, `col_BFRT` -> `col_bfrt`, `updateFtranBFRT` ->
//!   `update_ftran_bfrt`, `iRow` -> `i_row`, `workCount` -> `work_count`.
//!   (Rule: split `ABc` -> `A_Bc`, then `aB` -> `a_B`, then lowercase.)
//! * Rust keywords become raw identifiers (`HSet::in` -> `HSet::r#in`).
//! * Constants keep HiGHS names verbatim (`kHighsInf`, `kUpdateMethodFt`).
//! * `enum class` -> `#[repr(i32)]` Rust enum with verbatim variant names
//!   (`HighsModelStatus::kOptimal`); cast with `as i32` where C++ casts.
//! * Plain C `enum` -> `pub const kX: HighsInt = N;` (C++ code uses them as
//!   integers).
//! * Overloads get a suffixed name; each one's doc line names the C++
//!   overload it translates (e.g. `HFactor::setup` vs
//!   `HFactor::setup_raw`, `HEkk::set_basis` vs `HEkk::set_basis_from`).
//! * C++ default arguments become explicit parameters; the doc line
//!   `/// C++ default arguments: ...` lists the default values callers pass.
//!
//! # Type mapping
//!
//! | C++                               | Rust                                   |
//! |-----------------------------------|----------------------------------------|
//! | `HighsInt`                        | [`HighsInt`] (= `i32`, the HiGHS default build) |
//! | `double` / `bool` / `int8_t`      | `f64` / `bool` / `i8`                  |
//! | `uint64_t` / `HighsUInt` / `size_t` | `u64` / `u32` / `usize`              |
//! | `std::vector<T>`                  | `Vec<T>`                               |
//! | `const std::vector<T>&` param     | `&[T]`                                 |
//! | `std::vector<T>&` param           | `&mut Vec<T>`                          |
//! | `HVector&` / `HVector*`           | `&mut HVector` (`&HVector` if const)   |
//! | `HighsInt*` out-param (one value) | `&mut HighsInt`                        |
//! | `const double*` / `HighsInt*` array | `&[f64]` / `&mut [HighsInt]`         |
//! | nullable pointer param            | `Option<&T>` / `Option<&mut T>`        |
//! | `std::pair<A,B>`                  | `(A, B)` (`.first` -> `.0`)            |
//! | `std::set<HighsInt>`              | `BTreeSet<HighsInt>` (same ordering)   |
//! | `HighsHashTable<uint64_t>`        | `HashSet<u64>`                         |
//! | `std::string`                     | `String` (`&str` params)               |
//! | `FILE*`                           | `Option<HighsLogStream>` (stdout marker) |
//! | `HighsTimerClock*`                | `Option<&HighsTimerClock>` (always `None` here) |
//! | `HighsCDouble`                    | [`util::highs_cdouble::HighsCDouble`] (double-double, same ops) |
//!
//! Indices stay `HighsInt`; index a `Vec` with `v[i as usize]`.
//!
//! # Ownership pattern (the one rule for C++ references and pointers)
//!
//! HiGHS objects cache references/pointers into each other. Rust cannot
//! store those without lifetimes or `unsafe`, so the port uses **one
//! pattern: the owner keeps the data, and every borrower receives it as a
//! method parameter.** Concretely:
//!
//! **(a) Back-references to `HEkk` become a parameter.** `HEkkDual`,
//! `HEkkDualRow`, `HEkkDualRHS` and `HEkkPrimal` each hold
//! `HEkk& ekk_instance_` in C++. In Rust they do not store it: every method
//! takes `ekk: &mut HEkk` (`&HEkk` for C++ `const` methods) as its first
//! parameter after `self`, and `ekk_instance_.x` translates to `ekk.x`.
//! Cached raw-pointer aliases into `HEkk` arrays (`workDual`, `baseLower`,
//! `jMove`, `workRange`, ...) are **removed**; each struct's doc carries a
//! table mapping every removed alias to the Rust expression that replaces
//! it (e.g. `workDual[i]` -> `ekk.info_.work_dual_[i as usize]`). `HEkk::solve`
//! creates the solver locally (`let mut dual = HEkkDual::new(self);
//! dual.solve(self, force_phase2)`): `HEkkDual` is not a field of `HEkk`,
//! so there is no aliasing. `HEkkDual` owns its `dual_row`/`dual_rhs`
//! members and passes the same `ekk` down (`self.dual_row.choose_final(ekk)`).
//!
//! **(b) `HSimplexNla`'s pointer members become a borrowed context.**
//! C++ `HSimplexNla` stores `lp_`, `scale_`, `basic_index_`, `options_`
//! pointers into its owning `HEkk`. In Rust those are bundled in
//! [`simplex::nla::NlaRefs`], built on demand from disjoint `HEkk` fields
//! and passed to every `HSimplexNla` method that reads them:
//!
//! ```ignore
//! let refs = NlaRefs::new(&self.lp_, &self.basis_.basic_index_, &self.options_);
//! self.simplex_nla_.ftran(&refs, rhs, expected_density, None);
//! ```
//!
//! `NlaRefs::new` reproduces `setLpAndScalePointers` (`scale_` is
//! `Some(&lp.scale_)` iff `lp.scale_.has_scaling && !lp.is_scaled_`).
//! Methods that mutate the basic index (`invert`) take
//! `basic_index: &mut [HighsInt]` instead. Consequently
//! `setLpAndScalePointers`, `setBasicIndexPointers`, and
//! `HEkk::setNlaPointersForLpAndScale`/`ForTrans` are documented no-ops.
//! `timer_`/`analysis_` pointers are dropped (timing/analysis are stubs).
//!
//! **(c) `HFactor` owns copies of the constraint matrix and borrows the
//! basic index.** `a_start`/`a_index`/`a_value` pointers become owned
//! `Vec`s, refreshed by every `setup*`/`setup_matrix*` call (the only places
//! C++ re-points them). `basic_index` is not stored: every `HFactor` method
//! that reads or permutes it takes `basic_index: &mut [HighsInt]`
//! (`&[HighsInt]` when read-only).
//!
//! **(d) `HEkk` owns what its pointers pointed at.** `options_` is an owned
//! clone of the caller's [`lp_data::highs_options::HighsOptions`] (HEkk never
//! writes options); `timer_` and `callback_` are owned stubs
//! ([`util::highs_timer::HighsTimer`], [`lp_data::highs_callback::HighsCallback`]).
//!
//! **(e) `HighsLpSolverObject`** is a struct of `&'a mut` references, as in
//! C++, minus the `HEkk` reference (passed separately by (a)).
//!
//! # Out of scope (stubbed with the same call surface)
//!
//! * `HighsSimplexAnalysis`: a struct with the HiGHS fields and defaults;
//!   every method is a no-op (all `analyse_*` flags stay `false`).
//! * Timers: same call surface, no-op `start`/`stop`; `read` returns
//!   wall-clock seconds so time limits still work.
//! * `*Debug*.cpp`: functions return `HighsDebugStatus::kNotChecked` (or the
//!   value HiGHS returns when `highs_debug_level` is below the threshold).
//! * Logging: [`highs_log_user!`]/[`highs_log_dev!`] print when the log
//!   options enable them, so traces can be compared with HiGHS output.
//! * Threads: [`parallel`] is the port of the HiGHS task scheduler; PAMI / SIP
//!   parallel dual (`HEkkDualMulti.cpp`) runs on it.
//!
//! **(f) A member passed to a method of the same object.** C++ often passes
//! a pointer to one of its own members (`chooseColumn(&row_ep)`,
//! `updateFtranDSE(&col_DSE)`). Rust cannot borrow `self` mutably twice, so
//! move the member out for the call and back afterwards:
//!
//! ```ignore
//! let mut row_ep = std::mem::take(&mut self.row_ep);
//! self.choose_column(ekk, &mut row_ep);
//! self.row_ep = row_ep;
//! ```
//!
//! `HVector` and every struct passed this way implement `Default`, so the
//! take is cheap (no allocation). When C++ passes the *same* member twice,
//! or the callee also reads it through `this`, read the callee first and
//! keep the callee's view consistent with C++ (it sees the moved-out value).
//!
//! # Numerical fidelity: HARD RULE
//!
//! The golden-fixture oracle measured that FMA contraction changes the pivot
//! path on 35 of 64 traces. Therefore, in every body of this crate:
//!
//! * **Never** use `f64::mul_add`, `fma`, or any fused/compound float
//!   intrinsic, and never rely on a crate that does.
//! * Keep the C++ **evaluation order** of every floating-point expression
//!   exactly: same parenthesisation, same left-to-right association, same
//!   operand order (`a * b * c` is `(a * b) * c`). No reassociation, no
//!   algebraic simplification (`x / y` is not `x * (1.0 / y)`, `1 + base *
//!   r` stays as written, `a - b` is not `-(b - a)`), no hoisting that
//!   changes which intermediate is rounded. Accumulate sums in the C++
//!   loop order.
//! * Keep C++ integer/float conversions where C++ performs them
//!   (`1.0 * count / n`, `(double)i`), and `HighsInt` arithmetic in `i32`.
//! * [`util::highs_cdouble::HighsCDouble`] is a faithful non-FMA TwoSum /
//!   Dekker-split TwoProduct, translated operation by operation.
//! * C++ library semantics, not Rust's: `std::max`/`std::min` ->
//!   [`util::highs_utils::cpp_max`]/[`util::highs_utils::cpp_min`] (never
//!   `f64::max`/`f64::min`, which differ on NaN and signed zero);
//!   `pow(x, y)` -> `x.powf(y)` (never `powi`); `log` -> `ln`; `fabs` ->
//!   `abs`; `v.resize(n)` keeps existing entries (`v.resize(n as usize, 0)`,
//!   not clear + resize), while `v.assign(n, x)` replaces them; `assert` ->
//!   `debug_assert!`.
//!
//! Rust never contracts `a * b + c` into an FMA on its own. Apple-Silicon
//! clang does by default (`-ffp-contract=on`), so the bit-exact C++ oracle
//! is built with `-ffp-contract=off`.
//!
//! # Determinism the port must reproduce
//!
//! Doc lines marked `MUST-PORT (determinism, ...)` flag code whose exact
//! behaviour decides the pivot sequence:
//!
//! * `HEkk::random_` is seeded in `setSimplexOptions`, and
//!   `initialiseSimplexLpRandomVectors` runs twice on that same stream
//!   (`initialiseEkk`, then `initialiseForSolve`); the second draw is the
//!   one the solve uses. The same stream then drives the CHUZR random start
//!   (`HEkkDualRHS::chooseNormal`) and `HEkkDual::correctDualInfeasibilities`.
//!   [`util::highs_random::HighsRandom`] must be bit-exact (wrapping `u64`
//!   arithmetic, HiGHS's skipped output index 8).
//! * Reinversion is triggered by synthetic ticks (`HVector::synthetic_tick`,
//!   `HFactor::build_synthetic_tick`, `HEkk::{build,total}_synthetic_tick_`):
//!   every tick field and every tick increment is part of the contract.
//! * `pdqsort` is replaced by [`util::highs_sort::pdqsort`]; every call on the
//!   simplex path sorts distinct keys, so the result is the same permutation.
//! * Threads: the task-parallel branches (`iterateTasks`, `chooseColumnSlice`,
//!   `HEkkDualMulti.cpp`) run on the scheduler in [`parallel`]; with the default
//!   `simplex_strategy = kSimplexStrategyDual` (serial dual plain) they are not
//!   reached.

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(unused_mut)]
#![allow(dead_code)]
#![allow(unreachable_code)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::module_inception)]
#![allow(clippy::new_without_default)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::comparison_chain)]
#![allow(clippy::upper_case_acronyms)]
#![allow(clippy::type_complexity)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::nonminimal_bool)]
// HiGHS's `if (c) return true; return false;` stays as written (clippy 1.100 extends the lint)
#![allow(clippy::needless_bool)]
#![allow(clippy::if_same_then_else)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::manual_range_contains)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::float_cmp)]
#![allow(clippy::len_zero)]
#![allow(clippy::bool_comparison)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::derivable_impls)]
// C++ template names in `HiGHS:` doc lines (`HVectorBase<double>`) are not HTML.
#![allow(rustdoc::invalid_html_tags)]

/// HiGHS: `HighsInt` (HConfig.h; 32-bit default build).
pub type HighsInt = i32;
/// HiGHS: `HighsUInt` (HConfig.h; 32-bit default build).
pub type HighsUInt = u32;

#[macro_use]
pub mod io;
pub mod ipm;
pub mod lp_data;
pub mod mip;
pub mod model;
pub mod parallel;
pub mod pdlp;
pub mod presolve;
pub mod qpsolver;
pub mod simplex;
pub mod util;
pub mod work;

/// Glob re-export of every type, constant and free function in the crate,
/// so each file can start with `use crate::prelude::*;` (C++ sees all of
/// these through its `#include`s). Free functions whose names would clash
/// across modules are reached by path instead.
pub mod prelude {
    pub use crate::HighsInt;
    pub use crate::HighsUInt;
    pub use crate::io::highs_io::*;
    pub use crate::lp_data::hconst::*;
    pub use crate::lp_data::highs_basis_form::*;
    pub use crate::lp_data::highs_callback::*;
    pub use crate::lp_data::highs_iis::{
        HighsIis, HighsIisInfo, IisBoundStatus, IisModelStatus, kIisBoundStatusBoxed,
        kIisBoundStatusDropped, kIisBoundStatusFree, kIisBoundStatusLower, kIisBoundStatusNull,
        kIisBoundStatusUpper, kIisDevReport, kIisModelStatusFeasible, kIisModelStatusIrreducible,
        kIisModelStatusReducible, kIisModelStatusTimeLimit, kIisModelStatusUnknown,
    };
    pub use crate::lp_data::highs_info::*;
    pub use crate::lp_data::highs_lp::*;
    pub use crate::lp_data::highs_lp_modify::*;
    pub use crate::lp_data::highs_options::*;
    pub use crate::lp_data::highs_ranging::{HighsRanging, HighsRangingRecord};
    pub use crate::lp_data::highs_run::Highs;
    pub use crate::lp_data::highs_run_data::*;
    pub use crate::lp_data::highs_solution_debug::*;
    pub use crate::lp_data::highs_status::*;
    pub use crate::lp_data::hstruct::*;
    pub use crate::lp_data::lp_solver_object::*;
    pub use crate::lp_data::lp_utils::*;
    pub use crate::lp_data::model_utils::*;
    pub use crate::mip::*;
    pub use crate::model::highs_hessian::*;
    pub use crate::model::highs_hessian_utils::*;
    pub use crate::model::highs_model::*;
    pub use crate::parallel::shared_mut::SharedMut;
    pub use crate::presolve::highs_postsolve_stack::*;
    pub use crate::presolve::highs_symmetry::*;
    pub use crate::presolve::hpresolve::*;
    pub use crate::presolve::hpresolve_analysis::*;
    pub use crate::presolve::icrash::{
        ICrashInfo, ICrashIterationDetails, ICrashOptions, ICrashStrategy, Quadratic,
    };
    pub use crate::presolve::presolve_component::*;
    pub use crate::presolve::*;
    pub use crate::simplex::analysis::*;
    pub use crate::simplex::dual::*;
    pub use crate::simplex::dual_rhs::*;
    pub use crate::simplex::dual_row::*;
    pub use crate::simplex::ekk::*;
    pub use crate::simplex::hsimplex::*;
    pub use crate::simplex::hsimplex_debug::*;
    pub use crate::simplex::nla::*;
    pub use crate::simplex::primal::*;
    pub use crate::simplex::simplex_const::*;
    pub use crate::simplex::simplex_struct::*;
    pub use crate::simplex::simplex_timer::*;
    pub use crate::util::cpp_std::*;
    pub use crate::util::factor_timer::*;
    pub use crate::util::hfactor::*;
    pub use crate::util::hfactor_const::*;
    pub use crate::util::highs_cdouble::*;
    pub use crate::util::highs_data_stack::*;
    pub use crate::util::highs_disjoint_sets::*;
    pub use crate::util::highs_hash::*;
    pub use crate::util::highs_hash_table::*;
    pub use crate::util::highs_hash_tree::*;
    pub use crate::util::highs_integers::*;
    pub use crate::util::highs_linear_sum_bounds::*;
    pub use crate::util::highs_matrix_slice::*;
    pub use crate::util::highs_pdqsort::*;
    pub use crate::util::highs_random::*;
    pub use crate::util::highs_rb_tree::*;
    pub use crate::util::highs_sort::*;
    pub use crate::util::highs_sort_set::*;
    pub use crate::util::highs_sparse_vector_sum::*;
    pub use crate::util::highs_splay::*;
    pub use crate::util::highs_timer::*;
    pub use crate::util::highs_utils::*;
    pub use crate::util::hset::*;
    pub use crate::util::hvector::*;
    pub use crate::util::sparse_matrix::*;
    pub use std::cmp::Reverse;
    pub use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};
    pub use std::rc::Rc;
}
