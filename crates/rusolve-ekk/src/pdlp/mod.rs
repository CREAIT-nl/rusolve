//! HiGHS source: pdlp/ (cupdlp/*.{c,h}, CupdlpWrapper.{h,cpp}, hipdlp/*.{cc,hpp}, HiPdlpWrapper.{h,cpp})
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # The two first-order LP solvers of HiGHS
//!
//! | Rust module          | HiGHS source                                    | reached with      |
//! |----------------------|-------------------------------------------------|-------------------|
//! | [`cupdlp`]           | `pdlp/cupdlp/*.{c,h}` (cuPDLP-C, CPU build)     | `solver=pdlp`     |
//! | [`cupdlp_wrapper`]   | `pdlp/CupdlpWrapper.{h,cpp}`                    | `solver=pdlp`     |
//! | [`hipdlp`]           | `pdlp/hipdlp/*.{cc,hpp}` (CPU paths only)       | `solver=hipdlp`   |
//! | [`hipdlp_wrapper`]   | `pdlp/HiPdlpWrapper.{h,cpp}`                    | `solver=hipdlp`   |
//! | [`oracle_trace`]     | (not upstream) Rust mirror of the oracle hooks | tests / `RUSOLVE_PDLP_TRACE` |
//!
//! `lp_data/highs_solve.rs::solve_lp` dispatches here exactly as
//! `HighsSolve.cpp::solveLp` (lines 68, 93-116): only when
//! `use_pdlp(options.solver)` and not `use_only_ipm`.
//!
//! The build translated is HiGHS compiled **without CUDA**: `CUPDLP_CPU` is
//! defined, `CUPDLP_GPU` is not, `PDLP_DEBUG_LOG == 0`, `PDLP_PROFILE == 0`,
//! `CUPDLP_DEBUG` is undefined. Every `#ifndef CUPDLP_CPU` / `#ifdef
//! CUPDLP_GPU` branch, every `*Gpu` method and every CUDA struct field is
//! dropped; the macro expansions are listed in `cupdlp/defs.rs`.
//!
//! The module is not re-exported from the crate prelude: its C names
//! (`OPTIMAL`, `CPU`, `dot`, `Ax`, ...) would collide.

pub mod cupdlp;
pub mod cupdlp_wrapper;
pub mod hipdlp;
pub mod hipdlp_wrapper;
pub mod oracle_trace;
