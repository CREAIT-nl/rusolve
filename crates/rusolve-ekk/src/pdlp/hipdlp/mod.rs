//! HiGHS source: pdlp/hipdlp/ (defs.hpp, solver_results.hpp, linalg.{hpp,cc}, logger.{hpp,cc}, restart.{hpp,cc}, scaling.{hpp,cc}, pdhg.{hpp,cc})
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! HiPDLP, the HiGHS C++ PDLP (`solver=hipdlp`), CPU paths only.
//!
//! | Rust module          | HiGHS source                                   |
//! |----------------------|------------------------------------------------|
//! | [`defs`]             | `defs.hpp` + `pdhg.cc` 1903-1938 (the `initialise` methods) |
//! | [`solver_results`]   | `solver_results.hpp`                           |
//! | [`linalg`]           | `linalg.{hpp,cc}` (`namespace linalg`)         |
//! | [`logger`]           | `logger.{hpp,cc}` (state kept, printing is a no-op) |
//! | [`restart`]          | `restart.{hpp,cc}`                             |
//! | [`scaling`]          | `scaling.{hpp,cc}`                             |
//! | [`pdhg`]             | `pdhg.{hpp,cc}` (2910 lines, split by range; see [`pdhg`]) |
//!
//! Not ported: `pdhg.cu`, `pdhg_kernels.hpp`, `pdlp_gpu_debug.hpp` (GPU);
//! `defs.hpp` `PdlpIterate`, `pdlp_iterate_ops`, `DetailedTimings` (declared,
//! never used by the CPU build); `logger.cc` `EnumToString` (19-22),
//! `format` / `integer` / `sci` / `fix` (210-225, string builders for
//! printing only).
//!
//! # Back-pointers
//!
//! `Scaling::lp_` / `params_`, `RestartScheme::params_` / `log_options_` and
//! the `debug_pdlp_*` pointers are not stored: every method that reads them
//! takes the pointee as an explicit parameter, and `passLp` / `passParams` /
//! `passLogOptions` / `passDebugPdlp*` are not ported. The single write
//! through a back-pointer, `RestartScheme::updateBeta` setting
//! `params_->omega`, becomes `update_beta(beta, &mut PrimalDualParams)`.
//! `PDLPSolver` keeps `original_lp_` and `highs_timer_p_` as borrowed
//! `Option<&'a _>`.

pub mod defs;
pub mod linalg;
pub mod logger;
pub mod pdhg;
pub mod restart;
pub mod scaling;
pub mod solver_results;

pub use defs::*;
pub use logger::*;
pub use pdhg::*;
pub use restart::*;
pub use scaling::*;
pub use solver_results::*;

/// HiGHS: pdhg.hpp `kUseCupdlpx`
pub const kUseCupdlpx: bool = true;
/// HiGHS: pdhg.hpp `kTempSetting`
pub const kTempSetting: bool = true;
/// HiGHS: pdhg.cc `#define PDHG_CHECK_INTERVAL 40`
pub const PDHG_CHECK_INTERVAL: crate::HighsInt = 40;
/// HiGHS: pdhg.cc `const double kDivergentMovement = 1e10`
pub const kDivergentMovement: f64 = 1e10;
