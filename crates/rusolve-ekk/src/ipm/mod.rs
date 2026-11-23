//! HiGHS source: ipm/ (IPX interior point solver, BASICLU, IpxWrapper)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! # The HiGHS 1.15.1 IPX solver, faithful port
//!
//! * [`basiclu`]: `ipm/basiclu/*.c`, the sparse LU factorization with
//!   Forrest-Tomlin updates that IPX uses for its basis.
//! * [`ipx`]: `ipm/ipx/*.{h,cc}`, the interior point method, the basis
//!   preconditioner, crossover and the `LpSolver` driver.
//! * [`ipx_wrapper`]: `ipm/IpxWrapper.cpp` (IPX part only) and
//!   `ipm/IpxSolution.h`: HiGHS LP -> IPX model, IPX status -> HiGHS status,
//!   IPX solution -> HiGHS solution.
//!
//! * [`hipo`]: `ipm/hipo/*` plus the vendored orderings `extern/{metis,amd,rcm}`,
//!   the HiPO interior point solver with the FactorHiGHS direct solver.
//! * [`hipo_wrapper`]: `ipm/IpxWrapper.cpp` (HiPO part: `solveLpHipo`,
//!   `solveHipo`, `reportHipo*`, `getHipoNonVertexSolution`). Reached only
//!   when `HighsOptions::hipo_available` is set.
//!
//! None of these modules is in the crate prelude: IPX names (`Model`,
//! `Info`, `Control`, `Basis`, `Vector`, `SparseMatrix`, `Timer`) would clash
//! with HiGHS names. Import them through `crate::ipm::ipx::...`.

pub mod basiclu;
pub mod hipo;
pub mod hipo_wrapper;
pub mod ipx;
pub mod ipx_wrapper;
