//! HiGHS source: ipm/hipo/ipm/ (Solver, Model, Iterate, KktMatrix, LinearSolver, FactorHiGHSSolver, PreProcess, Refine, Control, Options, Info, Status, Parameters, IpmData)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The HiPO primal-dual interior point method (`ipm/hipo/ipm/*`).
//!
//! | Rust module | HiGHS source |
//! |-------------|--------------|
//! | [`status`], [`parameters`], [`options`], [`info`], [`ipm_data`] | `Status.h`, `Parameters.h`, `Options.h`, `Info.h`, `IpmData.{h,cpp}` (complete) |
//! | [`control`] | `Control.{h,cpp}` |
//! | [`kkt_matrix`] | `KktMatrix.{h,cpp}` |
//! | [`linear_solver`] | `LinearSolver.h` (trait + [`linear_solver::LsCtx`]) |
//! | [`factor_highs_solver`] | `FactorHiGHSSolver.{h,cpp}` |
//! | [`model`] | `Model.{h,cpp}` |
//! | [`pre_process`] | `PreProcess.{h,cpp}` |
//! | [`iterate`] | `Iterate.{h,cpp}` |
//! | [`solver`] | `Solver.{h,cpp}` |
//! | [`refine`] | `Refine.cpp` (methods of `Solver`) |

pub mod control;
pub mod factor_highs_solver;
pub mod info;
pub mod ipm_data;
pub mod iterate;
pub mod kkt_matrix;
pub mod linear_solver;
pub mod model;
pub mod options;
pub mod parameters;
pub mod pre_process;
pub mod refine;
pub mod solver;
pub mod status;
