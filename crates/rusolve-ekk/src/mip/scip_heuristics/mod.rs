//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! constants, numerics helpers, KISS RNG, S6Origin/S6Sol/S6Result/S6State
//! for the SCIP oneopt + crossover primal heuristics behind
//! `mip_heuristic_run_oneopt` / `mip_heuristic_run_crossover`.

use crate::mip::cut_selection::{S4RandNumGen, s4_initialize_random_seed_from};
use crate::prelude::*;
use crate::util::highs_utils::{cpp_max, cpp_min};

pub(crate) mod alns;
pub(crate) mod crossover;
pub(crate) mod driver;
pub(crate) mod guided_dive;
pub(crate) mod local_branching;
pub(crate) mod oneopt;
pub(crate) mod pool;
pub(crate) mod proximity;
pub(crate) mod shiftprop;
pub(crate) mod sort;
pub(crate) mod sub_mip_s6;
#[cfg(test)]
mod tests_e2e;

use pool::S6SolPool;

// ---------------------------------------------------------------------------
// Constants (manifest section 0 table)
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): heur_oneopt.c (line 73)
pub(crate) const S6_ONEOPT_WEIGHTEDOBJ: bool = true;
/// SCIP 11.0.0 (32dd425): heur_oneopt.c (line 74)
pub(crate) const S6_ONEOPT_DURINGROOT: bool = true;
/// SCIP 11.0.0 (32dd425): heur_oneopt.c (line 75) — path not ported
pub(crate) const S6_ONEOPT_BEFOREPRESOL: bool = false;
/// SCIP 11.0.0 (32dd425): heur_oneopt.c (line 76)
pub(crate) const S6_ONEOPT_FORCELPCONSTRUCTION: bool = false;
/// SCIP 11.0.0 (32dd425): heur_oneopt.c (line 77)
pub(crate) const S6_ONEOPT_USELOOP: bool = true;

/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 72)
pub(crate) const S6_XOVER_MAXNODES: i64 = 5000;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 73)
pub(crate) const S6_XOVER_MINIMPROVE: f64 = 0.01;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 74)
pub(crate) const S6_XOVER_MINNODES: i64 = 50;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 75)
pub(crate) const S6_XOVER_MINFIXINGRATE: f64 = 0.666;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 76)
pub(crate) const S6_XOVER_NODESOFS: i64 = 500;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 77)
pub(crate) const S6_XOVER_NODESQUOT: f64 = 0.1;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 78) — not used
pub(crate) const S6_XOVER_LPLIMFAC: f64 = 2.0;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 79)
pub(crate) const S6_XOVER_NUSEDSOLS: i32 = 3;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 80)
pub(crate) const S6_XOVER_NWAITINGNODES: i64 = 200;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 81)
pub(crate) const S6_XOVER_RANDOMIZATION: bool = true;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 82)
pub(crate) const S6_XOVER_DONTWAITATROOT: bool = false;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 83)
pub(crate) const S6_XOVER_USELPROWS: bool = false;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 85) — Q8
pub(crate) const S6_XOVER_COPYCUTS: bool = true;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 88)
pub(crate) const S6_XOVER_PERMUTE: bool = false;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 89) — capacity hint only
pub(crate) const S6_XOVER_HASHSIZE_SOLS: usize = 500;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 90)
pub(crate) const S6_XOVER_BESTSOLLIMIT: i64 = -1;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 91)
pub(crate) const S6_XOVER_USEUCT: bool = false;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 92)
pub(crate) const S6_XOVER_RANDSEED: u32 = 7;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (lines 527-531) — nfailures <= 25
pub(crate) const S6_XOVER_MAXFAILURES_SHIFT: i32 = 25;
/// SCIP 11.0.0 (32dd425): heur_crossover.c (line 322)
pub(crate) const S6_XOVER_SELECT_TRIES: i32 = 10;

/// SCIP 11.0.0 (32dd425): set.c (line 241) — limits/maxsol, pool capacity
pub(crate) const S6_MAXSOL: usize = 100;
/// SCIP 11.0.0 (32dd425): def.h (line 163)
pub(crate) const S6_EPSILON: f64 = 1e-9;
/// SCIP 11.0.0 (32dd425): def.h (line 164)
pub(crate) const S6_SUMEPSILON: f64 = 1e-6;
/// SCIP 11.0.0 (32dd425): def.h (line 162)
pub(crate) const S6_INFINITY: f64 = 1e20;

/// SCIP 11.0.0 (32dd425): sorttpl.c (line 49)
pub(crate) const S6_SORTTPL_SHELLSORTMAX: usize = 25;

// ---------------------------------------------------------------------------
// Numerics helpers (literal ports, all #[inline])
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsInfinity (lines 6824-6832)
#[inline]
pub(crate) fn s6_is_inf(v: f64) -> bool {
    v >= S6_INFINITY
}

/// SCIP 11.0.0 (32dd425): misc.c::SCIPrelDiff (lines 11164-11178)
#[inline]
pub(crate) fn s6_rel_diff(a: f64, b: f64) -> f64 {
    let absval1 = a.abs();
    let absval2 = b.abs();
    let quot = cpp_max(cpp_max(1.0, absval1), absval2);
    (a - b) / quot
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsFeasLT (lines 7278-7299)
#[inline]
pub(crate) fn s6_feas_lt(a: f64, b: f64, ft: f64) -> bool {
    let diff = s6_rel_diff(a, b);
    diff < -ft
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsFeasLE (lines 7302-7323)
#[inline]
pub(crate) fn s6_feas_le(a: f64, b: f64, ft: f64) -> bool {
    let diff = s6_rel_diff(a, b);
    !(diff > ft)
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsFeasGT (lines 7326-7347)
#[inline]
pub(crate) fn s6_feas_gt(a: f64, b: f64, ft: f64) -> bool {
    let diff = s6_rel_diff(a, b);
    diff > ft
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsFeasGE (lines 7350-7371)
#[inline]
pub(crate) fn s6_feas_ge(a: f64, b: f64, ft: f64) -> bool {
    let diff = s6_rel_diff(a, b);
    !(diff < -ft)
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsFeasZero (lines 7374-7382)
#[inline]
pub(crate) fn s6_feas_zero(v: f64, ft: f64) -> bool {
    v.abs() <= ft
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetFeasFloor (lines 7433-7442), def.h:190
#[inline]
pub(crate) fn s6_feas_floor(v: f64, ft: f64) -> f64 {
    (v + ft).floor()
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetFeasCeil (lines 7445-7454), def.h:191
#[inline]
pub(crate) fn s6_feas_ceil(v: f64, ft: f64) -> f64 {
    (v - ft).ceil()
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsEQ (lines 6846-6862), def.h:182
#[inline]
pub(crate) fn s6_is_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= S6_EPSILON
}

// ---------------------------------------------------------------------------
// S6Origin
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum S6Origin {
    Highs(i32),
    Oneopt,
    Crossover,
}

impl Default for S6Origin {
    fn default() -> Self {
        S6Origin::Highs(crate::mip::mip_solver_data::kSolutionSourceNone)
    }
}

// ---------------------------------------------------------------------------
// S6Sol
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default)]
pub(crate) struct S6Sol {
    pub index: i64,
    pub obj: f64,
    pub x: Vec<f64>,
    pub origin: S6Origin,
    pub nodenum: i64,
}

// ---------------------------------------------------------------------------
// S6Result
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum S6Result {
    DidNotRun,
    Delayed,
    DidNotFind,
    FoundSol,
}

// ---------------------------------------------------------------------------
// s6_get_int (KISS generator and seed helper are S4's, cut_selection/mod.rs)
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): misc.c::SCIPrandomGetInt (lines 10223-10241)
pub(crate) fn s6_get_int(rng: &mut S4RandNumGen, min: HighsInt, max: HighsInt) -> HighsInt {
    let randnumber = (rng.get_rand() as f64) / (u32::MAX as f64 + 1.0);
    let zeromax: i64 = max as i64 - min as i64 + 1;
    ((zeromax as f64 * randnumber) as i64 + min as i64) as HighsInt
}

// ---------------------------------------------------------------------------
// S6State
// ---------------------------------------------------------------------------

#[derive(Clone, Default, Debug)]
pub(crate) struct S6State {
    pub pool: S6SolPool,
    /// monotone SCIPsolGetIndex counter
    pub next_sol_index: i64,
    /// transient tag set around calls
    pub cur_origin: S6Origin,
    // oneopt heurdata (heur_oneopt.c:85-88)
    /// -1 (heurInit :454-472)
    pub oneopt_lastsolindex: i64,
    pub oneopt_ncalls: i64,
    // crossover heurdata (heur_crossover.c:97-133)
    pub xover_rng: Option<S4RandNumGen>,
    /// SOLTUPLE indices
    pub xover_hashtable: std::collections::HashSet<Vec<i32>>,
    /// kept for parity; only written
    pub xover_lasttuple: Option<Vec<i32>>,
    pub xover_usednodes: i64,
    /// solution index
    pub xover_prevbestsol: Option<i64>,
    pub xover_prevlastsol: Option<i64>,
    pub xover_nfailures: i32,
    pub xover_nextnodenumber: i64,
    pub xover_prevnsols: i32,
    /// SCIPheurGetNCalls (heur.c:1263-1340)
    pub xover_ncalls: i64,
    /// SCIPheurGetNBestSolsFound
    pub xover_nbestsolsfound: i64,
}

impl S6State {
    /// SCIP 11.0.0 (32dd425): heur_oneopt.c::heurInitOneopt (lines 454-472)
    /// SCIP 11.0.0 (32dd425): heur_crossover.c::heurInitCrossover (lines 835-864)
    pub(crate) fn new(ms: &crate::mip::mip_solver::HighsMipSolver) -> S6State {
        S6State::from_flags(
            // mip_alns draws its crossover parents with this RNG (alns.rs)
            ms.options_mip_.mip_heuristic_run_crossover || ms.options_mip_.mip_alns,
            ms.options_mip_.random_seed,
        )
    }

    /// SCIP 11.0.0 (32dd425): heur_oneopt.c::heurInitOneopt (lines 454-472)
    /// SCIP 11.0.0 (32dd425): heur_crossover.c::heurInitCrossover (lines 835-864)
    /// (`S6State::new` with the two option values passed explicitly)
    pub(crate) fn from_flags(crossover_on: bool, random_seed: HighsInt) -> S6State {
        S6State {
            pool: S6SolPool::default(),
            next_sol_index: 0,
            cur_origin: S6Origin::default(),
            // oneopt heurInit (heur_oneopt.c:454-472)
            oneopt_lastsolindex: -1,
            oneopt_ncalls: 0,
            // crossover heurInit (heur_crossover.c:835-864)
            xover_rng: if crossover_on {
                Some(S4RandNumGen::new(s4_initialize_random_seed_from(
                    S6_XOVER_RANDSEED,
                    random_seed,
                )))
            } else {
                None
            },
            xover_hashtable: std::collections::HashSet::new(),
            xover_lasttuple: None,
            xover_usednodes: 0,
            xover_prevbestsol: None,
            xover_prevlastsol: None,
            xover_nfailures: 0,
            xover_nextnodenumber: 0,
            xover_prevnsols: 0,
            xover_ncalls: 0,
            xover_nbestsolsfound: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s6_data_feas_helpers_edges() {
        // s6_rel_diff: relative vs absolute (|a| < 1 uses the 1.0 denominator)
        // rel_diff(0.5, 0.5 + 1e-10) = (0.5 - (0.5+1e-10)) / max(1, 0.5, 0.5+1e-10)
        //                            = -1e-10 / 1.0 = -1e-10
        let rd = s6_rel_diff(0.5, 0.5 + 1e-10);
        // exact f64 value in C order: (a - b) / max(1, |a|, |b|)
        assert_eq!(rd, (0.5 - (0.5 + 1e-10)) / 1.0);

        // rel_diff(100.0, 100.0 + 1e-8) = -1e-8 / max(1, 100, 100+1e-8) = -1e-8/100 = -1e-10
        let rd2 = s6_rel_diff(100.0, 100.0 + 1e-8);
        assert_eq!(rd2, (100.0 - (100.0 + 1e-8)) / (100.0 + 1e-8));

        // s6_is_inf
        assert!(s6_is_inf(1e20));
        assert!(!s6_is_inf(1e19));

        // s6_feas_zero boundary
        let ft = 1e-6;
        assert!(s6_feas_zero(ft, ft));
        assert!(s6_feas_zero(-ft, ft));
        assert!(!s6_feas_zero(ft + 1e-15, ft));

        // s6_feas_floor / s6_feas_ceil with ft 1e-6
        // s6_feas_floor(2.0000001, 1e-6) = floor(2.0000001 + 1e-6) = floor(2.0000011) = 2.0
        assert_eq!(s6_feas_floor(2.0000001, 1e-6), 2.0);
        // s6_feas_ceil(1.9999999, 1e-6) = ceil(1.9999999 - 1e-6) = ceil(1.9999989) = 2.0
        assert_eq!(s6_feas_ceil(1.9999999, 1e-6), 2.0);
    }

    #[test]
    fn s6_data_initialize_random_seed() {
        assert_eq!(s4_initialize_random_seed_from(S6_XOVER_RANDSEED, 0), 7);
        assert_eq!(s4_initialize_random_seed_from(S6_XOVER_RANDSEED, 1), 68);
    }

    #[test]
    fn s6_data_get_int_seed7() {
        let mut rng2 = S4RandNumGen::new(7);
        assert_eq!(s6_get_int(&mut rng2, 0, 9), 0);
        assert_eq!(s6_get_int(&mut rng2, 0, 9), 3);
        assert_eq!(s6_get_int(&mut rng2, 0, 9), 1);
    }

    #[test]
    fn s6_data_get_int_range() {
        // s6_get_int(rng, a, a) == a for a = -3 and a = 5 (fresh rng each)
        let mut rng = S4RandNumGen::new(7);
        assert_eq!(s6_get_int(&mut rng, -3, -3), -3);

        let mut rng2 = S4RandNumGen::new(7);
        assert_eq!(s6_get_int(&mut rng2, 5, 5), 5);

        // For seed 7, s6_get_int(rng, 0, 99) five times gives 9, 37, 19, 57, 28
        let mut rng3 = S4RandNumGen::new(7);
        assert_eq!(s6_get_int(&mut rng3, 0, 99), 9);
        assert_eq!(s6_get_int(&mut rng3, 0, 99), 37);
        assert_eq!(s6_get_int(&mut rng3, 0, 99), 19);
        assert_eq!(s6_get_int(&mut rng3, 0, 99), 57);
        assert_eq!(s6_get_int(&mut rng3, 0, 99), 28);
    }

    #[test]
    fn s6_data_state_default_none() {
        let state = S6State::from_flags(false, 0);
        assert!(state.xover_rng.is_none());
        assert!(state.pool.sols.is_empty());

        let state_default = S6State::default();
        assert!(state_default.xover_rng.is_none());
    }

    #[test]
    fn s6_data_state_crossover_rng() {
        let mut state = S6State::from_flags(true, 0);
        assert!(state.xover_rng.is_some());
        let rng = state.xover_rng.as_mut().unwrap();
        assert_eq!(s6_get_int(rng, 0, 9), 0);
    }
}
