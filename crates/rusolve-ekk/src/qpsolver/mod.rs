//! HiGHS source: qpsolver/ (QUASS active-set QP solver)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! One Rust file per HiGHS file;
//! `basis.cpp` and `quass.cpp` are split into directories.
//!
//! NOT in the crate prelude (the names `Matrix`, `Instance`, `Settings`,
//! `Runtime`, `Basis`, `Statistics`, `scale`, `reduce` are generic). Each
//! qpsolver file does `use crate::prelude::*; use crate::qpsolver::*;`;
//! outside it, reach items as `crate::qpsolver::X`.
//!
//! Ownership (lib.rs rule (a)): every C++ `Runtime&` / `Basis&` /
//! `Gradient&` / `ReducedCosts&` member reference becomes a method
//! parameter; see each module doc. `Runtime` owns its `Statistics`
//! (`runtime.rs`).

pub mod a_asm;
pub mod a_quass;
pub mod basis;
pub mod dantzigpricing;
pub mod devexpricing;
pub mod eventhandler;
pub mod factor;
pub mod feasibility_bounded;
pub mod feasibility_highs;
pub mod gradient;
pub mod instance;
pub mod matrix;
pub mod oracle_trace;
pub mod perturbation;
pub mod pricing;
pub mod qpconst;
pub mod qpvector;
pub mod quass;
pub mod ratiotest;
pub mod reducedcosts;
pub mod reducedgradient;
pub mod runtime;
pub mod scaling;
pub mod settings;
pub mod snippets;
pub mod statistics;
pub mod steepestedgepricing;

pub use a_asm::*;
pub use a_quass::*;
pub use basis::*;
pub use dantzigpricing::*;
pub use devexpricing::*;
pub use eventhandler::*;
pub use factor::*;
pub use feasibility_bounded::*;
pub use feasibility_highs::*;
pub use gradient::*;
pub use instance::*;
pub use matrix::*;
pub use perturbation::*;
pub use pricing::*;
pub use qpconst::*;
pub use qpvector::*;
pub use quass::Quass;
pub use ratiotest::*;
pub use reducedcosts::*;
pub use reducedgradient::*;
pub use runtime::*;
pub use scaling::*;
pub use settings::*;
pub use statistics::*;
pub use steepestedgepricing::*;
