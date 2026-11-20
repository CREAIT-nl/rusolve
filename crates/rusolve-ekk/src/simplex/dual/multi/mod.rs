//! HiGHS source: simplex/HEkkDualMulti.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! PAMI (parallel multiple-pricing) dual, split by C++ line range:
//!
//! | file | C++ lines | part |
//! |---|---|---|
//! | `choose.rs` | 26-297 | PAMI-A |
//! | `minor.rs` | 299-533 | PAMI-B |
//! | `major_update.rs` | 535-763 | PAMI-C |
//! | `major_primal.rs` | 765-1020 | PAMI-D |
//!
//! Reached only with `simplex_strategy = kSimplexStrategyDualMulti` (3); the default
//! serial dual never calls into here.

mod choose;
mod major_primal;
mod major_update;
mod minor;
