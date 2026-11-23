//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: cutsel_hybrid.c 42-74; misc.c 10131-10258; pub_misc.h 568;
//!   scip_randnumgen.c 56-76; set.c 343-344, 460-481, 8109-8117; def.h 163-165.

pub(crate) mod apply;
pub(crate) mod hybrid;
pub(crate) mod pool_s4;
pub(crate) mod row_metrics;
#[cfg(test)]
mod tests_e2e;

use crate::prelude::*;

/// SCIP cutsel_hybrid.c:46
#[allow(dead_code)] // not used yet
pub(crate) const S4_RANDSEED: u32 = 0x5EED;

/// SCIP cutsel_hybrid.c:47
#[allow(dead_code)] // not used yet
pub(crate) const S4_GOODSCORE: f64 = 0.9;

/// SCIP cutsel_hybrid.c:48
#[allow(dead_code)] // not used yet
pub(crate) const S4_BADSCORE: f64 = 0.0;

/// SCIP cutsel_hybrid.c:50
#[allow(dead_code)] // not used yet
pub(crate) const S4_EFFICACYWEIGHT: f64 = 1.0;

/// SCIP cutsel_hybrid.c:51
#[allow(dead_code)] // not used yet
pub(crate) const S4_DIRCUTOFFDISTWEIGHT: f64 = 0.0;

/// SCIP cutsel_hybrid.c:52
#[allow(dead_code)] // not used yet
pub(crate) const S4_OBJPARALWEIGHT: f64 = 0.1;

/// SCIP cutsel_hybrid.c:53
#[allow(dead_code)] // not used yet
pub(crate) const S4_INTSUPPORTWEIGHT: f64 = 0.1;

/// SCIP cutsel_hybrid.c:54
#[allow(dead_code)] // not used yet
pub(crate) const S4_MINORTHO: f64 = 0.90;

/// SCIP cutsel_hybrid.c:55
#[allow(dead_code)] // not used yet
pub(crate) const S4_MINORTHOROOT: f64 = 0.90;

/// SCIP set.c:480 (separating/maxcuts)
#[allow(dead_code)] // not used yet
pub(crate) const S4_MAXCUTS: HighsInt = 100;

/// SCIP set.c:481 (separating/maxcutsroot)
#[allow(dead_code)] // not used yet
pub(crate) const S4_MAXCUTSROOT: HighsInt = 2000;

/// SCIP set.c:463 (separating/efficacynorm, only 'e' ported)
#[allow(dead_code)] // not used yet
pub(crate) const S4_EFFICACYNORM: char = 'e';

/// SCIP set.c:344 (randomization/randseedshiftmultiplier default)
pub(crate) const S4_RANDSEEDMULT: u32 = 10;

/// SCIP set.c:343 (randomization/randseedshift default)
#[allow(dead_code)] // not used yet
pub(crate) const S4_RANDSEEDSHIFT_DEFAULT: HighsInt = 0;

/// SCIP def.h:163
#[allow(dead_code)] // not used yet
pub(crate) const S4_EPSILON: f64 = 1e-9;

/// SCIP def.h:164
#[allow(dead_code)] // not used yet
pub(crate) const S4_SUMEPSILON: f64 = 1e-6;

/// SCIP def.h:162 (SCIP_DEFAULT_INFINITY, numerics/infinity default)
#[allow(dead_code)] // not used yet
pub(crate) const S4_INFINITY: f64 = 1e20;

/// SCIP misc.c:10132
pub(crate) const S4_KISS_DEFAULT_SEED: u32 = 123456789;

/// SCIP misc.c:10133
pub(crate) const S4_KISS_DEFAULT_XOR: u32 = 362436000;

/// SCIP misc.c:10134
pub(crate) const S4_KISS_DEFAULT_MWC: u32 = 521288629;

/// SCIP misc.c:10135
pub(crate) const S4_KISS_DEFAULT_CST: u32 = 7654321;

/// SCIP 11.0.0 (32dd425): pub_misc.h::SCIPhashTwo (line 568)
pub(crate) fn s4_hash_two(a: u32, b: u32) -> u32 {
    (((a as u64).wrapping_add(0xd37e9a1ce2148403))
        .wrapping_mul((b as u64).wrapping_add(0xe5fcc163aef32782))
        >> 32) as u32
}

/// SCIP 11.0.0 (32dd425): misc.c::SCIP_RANDNUMGEN (KISS random number generator)
#[derive(Clone, Debug)]
pub(crate) struct S4RandNumGen {
    pub(crate) seed: u32,
    pub(crate) xor_seed: u32,
    pub(crate) mwc_seed: u32,
    pub(crate) cst_seed: u32,
}

impl S4RandNumGen {
    /// SCIP 11.0.0 (32dd425): misc.c::SCIPrandomCreate (lines 10193-10206)
    pub(crate) fn new(initialseed: u32) -> Self {
        let mut rng = S4RandNumGen {
            seed: 0,
            xor_seed: 0,
            mwc_seed: 0,
            cst_seed: 0,
        };
        rng.set_seed(initialseed);
        rng
    }

    /// SCIP 11.0.0 (32dd425): misc.c::SCIPrandomSetSeed (lines 10139-10155)
    pub(crate) fn set_seed(&mut self, initseed: u32) {
        self.seed = s4_hash_two(S4_KISS_DEFAULT_SEED, initseed).max(1);
        self.xor_seed = s4_hash_two(S4_KISS_DEFAULT_XOR, initseed).max(1);
        self.mwc_seed = s4_hash_two(S4_KISS_DEFAULT_MWC, initseed).max(1);
        self.cst_seed = s4_hash_two(S4_KISS_DEFAULT_CST, initseed);
    }

    /// SCIP 11.0.0 (32dd425): misc.c::randomGetRand (lines 10169-10190)
    pub(crate) fn get_rand(&mut self) -> u32 {
        // linear congruential
        self.seed = (self.seed as u64)
            .wrapping_mul(1103515245)
            .wrapping_add(12345) as u32;

        // Xorshift
        self.xor_seed ^= self.xor_seed.wrapping_shl(13);
        self.xor_seed ^= self.xor_seed >> 17;
        self.xor_seed ^= self.xor_seed.wrapping_shl(5);

        // Multiply-with-carry
        let t: u64 = 698769069u64
            .wrapping_mul(self.mwc_seed as u64)
            .wrapping_add(self.cst_seed as u64);
        self.cst_seed = (t >> 32) as u32;
        self.mwc_seed = t as u32;

        self.seed
            .wrapping_add(self.xor_seed)
            .wrapping_add(self.mwc_seed)
    }

    /// SCIP 11.0.0 (32dd425): misc.c::SCIPrandomGetReal (lines 10245-10259)
    pub(crate) fn get_real(&mut self, minrandval: f64, maxrandval: f64) -> f64 {
        let randnumber = self.get_rand() as f64 / u32::MAX as f64;
        minrandval * (1.0 - randnumber) + maxrandval * randnumber
    }
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetInitializeRandomSeed (lines 8109-8117)
/// and cutsel_hybrid.c:314 (SCIPcreateRandom with useglobalseed TRUE, scip_randnumgen.c:56-76)
pub(crate) fn s4_initialize_random_seed(shift: HighsInt) -> u32 {
    s4_initialize_random_seed_from(S4_RANDSEED, shift)
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetInitializeRandomSeed (lines 8109-8117) for any
/// component's initial seed (SCIPcreateRandom with useglobalseed TRUE); shared with crossover
pub(crate) fn s4_initialize_random_seed_from(initialseed: u32, shift: HighsInt) -> u32 {
    initialseed.wrapping_add((1 + 6 * S4_RANDSEEDMULT).wrapping_mul(shift as u32))
}

/// candidate cut with precomputed metrics per round.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct S4Cand {
    pub(crate) pool: HighsInt,
    pub(crate) row: HighsInt,
    pub(crate) efficacy: f64,
    pub(crate) objpar: f64,
    /// SCIPgetRowNumIntCols, as the C converts it for `w * n / nnz` (cutsel_hybrid.c:112)
    pub(crate) numintcols: f64,
    /// SCIProwGetNNonz as (SCIP_Real)
    pub(crate) nnz: f64,
    pub(crate) cutoffdist: f64,
    pub(crate) in_global: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s4_data_hash_two_fixed() {
        assert_eq!(s4_hash_two(123456789, 0x5EED), 668882365);
        assert_eq!(s4_hash_two(1, 2), 2939976213);
    }

    #[test]
    fn s4_data_set_seed_state() {
        let r = S4RandNumGen::new(0x5EED);
        assert_eq!(r.seed, 668882365);
        assert_eq!(r.xor_seed, 2559441955);
        assert_eq!(r.mwc_seed, 4068887404);
        assert_eq!(r.cst_seed, 2102711353);
    }

    #[test]
    fn s4_data_get_rand_first_three() {
        let mut r = S4RandNumGen::new(0x5EED);
        assert_eq!(r.get_rand(), 2969581358);
        assert_eq!(r.get_rand(), 3308917181);
        assert_eq!(r.get_rand(), 924377526);
    }

    #[test]
    fn s4_data_get_real_first() {
        let mut r = S4RandNumGen::new(0x5EED);
        let x = r.get_real(0.0, 1e-6);
        assert!((x - 6.914095391266536e-07).abs() < 1e-15);
    }

    #[test]
    fn s4_data_seed_shift() {
        assert_eq!(s4_initialize_random_seed(0), 0x5EED);
        assert_eq!(s4_initialize_random_seed(1), 24362);
    }

    #[test]
    fn s4_data_consts() {
        assert_eq!(S4_RANDSEED, 0x5EED);
        assert_eq!(S4_RANDSEEDMULT, 10);
        assert_eq!(S4_MAXCUTS, 100);
        assert_eq!(S4_MAXCUTSROOT, 2000);
    }
}
