//! HiGHS source: pdlp/cupdlp/ (cupdlp.h and the .c files below)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! cuPDLP-C, CPU build. One Rust module per C translation unit; the
//! prototypes of each `cupdlp_X.h` live in the module of `cupdlp_X.c`.
//! A `.c` longer than ~700 lines is a directory of contiguous line ranges.
//!
//! | Rust module      | C source                                   |
//! |------------------|--------------------------------------------|
//! | [`defs`]         | `glbopts.h`, `cupdlp_defs.h` (types, macros, enums, structs) |
//! | [`cs`]           | `cupdlp_cs.{h,c}` (CSparse subset)          |
//! | [`linalg`]       | `cupdlp_linalg.{h,c}` (802 lines: `mod.rs` 17-390, `part2.rs` 457-802) |
//! | [`proj`]         | `cupdlp_proj.{h,c}`                         |
//! | [`restart`]      | `cupdlp_restart.{h,c}`                      |
//! | [`scaling`]      | `cupdlp_scaling.{h,c}`                      |
//! | [`step`]         | `cupdlp_step.{h,c}`                         |
//! | [`solver`]       | `cupdlp_solver.{h,c}` (1498 lines: `mod.rs` 12-311, `part2.rs` 313-591, `check.rs` 598-897, `solve.rs` 899-1215, `post.rs` 1217-1498) |
//! | [`utils`]        | `cupdlp_utils.{h,c}` (1850 lines: `mod.rs` 651-1030, `work.rs` 1032-1158, `matrix.rs` 1160-1525) |
//!
//! Not ported: `cupdlp_utils.c`
//! clear functions 22-417 (`Drop`), print functions 419-558 and 1527-1628,
//! `getUserParam` 560-649 (command line), `writeJson` / `writeSol`
//! 1630-1783 (file writers, deferred), debug log 1785-1850
//! (`PDLP_DEBUG_LOG == 0`); `cupdlp_cs.c::cupdlp_dcs_print` 145-187 and the
//! raw allocator wrappers 5-29, 61-76; the GPU functions of
//! `cupdlp_linalg.c` 394-456; `cupdlp/cuda/*`.

pub mod cs;
pub mod defs;
pub mod linalg;
pub mod proj;
pub mod restart;
pub mod scaling;
pub mod solver;
pub mod step;
pub mod utils;

pub use defs::*;
