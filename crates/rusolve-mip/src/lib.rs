//! Mixed-integer linear programming: branch-and-bound over an LP relaxation.
//!
//! Contents:
//!
//! * [`problem::MipProblem`] — the solver's own MIP form: CSC matrix, two-sided
//!   row bounds, minimization-form costs, column bounds, integrality flags.
//! * [`engine::LpEngine`] — the narrow interface the tree search uses to solve
//!   node relaxations; [`engine::RusolveEngine`] backs it with `rusolve-lp`.
//! * [`bnb::solve`] — LP-based branch-and-bound: best-bound node selection with
//!   depth-first plunging, most-fractional branching, objective-integrality
//!   aware cutoff.
//! * [`presolve::presolve`] — MIP presolve (gcd row scaling, coefficient
//!   strengthening, singleton/doubleton rows, fixed/dual-fixed columns) with a
//!   primal postsolve; [`bnb::solve`] runs it first unless
//!   `MipConfig::presolve` is off.
//! * [`propagate::Propagator`] — activity-based bound propagation on rows and
//!   the objective cutoff, run at every node before its LP.
//! * [`redcost`] — reduced-cost fixing:
//!   node-local after every node LP with an incumbent, and global from the
//!   stored root reduced costs at every new incumbent (after HiGHS
//!   `HighsRedcostFixing`, MIT).
//! * [`conflict`] — conflict analysis from the Farkas proofs of infeasible
//!   node LPs and a conflict pool propagated at every node
//!   (after HiGHS `HighsConflictPool` and
//!   `HighsDomain::ConflictSet`, MIT).
//! * [`cuts`] — the root cutting-plane loop: separators that need only
//!   the LP point and the rows, a cut pool, and HiGHS-style round control.
//! * [`classify::classify`] — MIPLIB 2017 row-type histogram, column /
//!   objective statistics, structure detectors and a fixed-length feature
//!   vector (input to tuning profiles).
//! * [`submip::rins`] / [`submip::rens`] — sub-MIP primal heuristics that
//!   solve a fixed or rounding-boxed restriction with [`bnb::solve`]
//!   (sources in the module docs).
//! * [`probing`] — root probing: fix binaries, propagate, apply the
//!   implied bound tightenings to the global box and feed the binary
//!   implications to the clique table (after HiGHS
//!   `HighsImplications::computeImplications` and `HPresolve::runProbing`,
//!   MIT; the seed the clique-table fixings link was measured to lack).
//!
//! Provenance (per AGENTS.md, every algorithmic idea cites a source):
//!
//! * LP-based branch-and-bound: A. H. Land & A. G. Doig, "An automatic method
//!   of solving discrete programming problems", Econometrica 28(3):497–520
//!   (1960).
//! * Plunging + best-bound node selection, most-fractional branching as the
//!   baseline rule: T. Achterberg, *Constraint Integer Programming*, PhD
//!   thesis, TU Berlin (2007), §5.3 and §6.1.
//! * Unbounded-relaxation handling: R. R. Meyer, "On the existence of optimal
//!   solutions to integer and mixed-integer programming problems", Math.
//!   Prog. 7:223–235 (1974) — with rational data, a feasible MIP whose LP
//!   relaxation is unbounded is itself unbounded.
//! * Integral-objective cutoff: HiGHS `HighsMipSolverData::computeNewUpperLimit`
//!   (reference/HiGHS/highs/mip/HighsMipSolverData.cpp:628, MIT, read as a
//!   reference text; reimplemented, not translated).
//!
//! No `unsafe`.

pub mod alns;
pub mod bnb;
pub mod branching;
pub mod classify;
pub mod conflict;
pub mod cuts;
pub mod engine;
pub mod engine_select;
pub mod heuristics;
pub mod nodequeue;
pub mod presolve;
pub mod probing;
pub mod problem;
pub mod propagate;
pub mod pseudocost;
mod race;
pub mod redcost;
pub mod restart;
pub mod submip;

pub use bnb::{MipConfig, MipError, MipSolution, MipStatus, solve};
pub use conflict::ConflictStats;
pub use cuts::{CutClass, CutConfig, CutStats};
pub use engine::{
    CutRow, LpDuals, LpEngine, LpError, LpOutcome, LpStats, LpStatus, RusolveEngine, WarmStart,
};
pub use probing::ProbingStats;
pub use problem::MipProblem;
pub use redcost::RedcostStats;
pub use submip::SubMipStats;
