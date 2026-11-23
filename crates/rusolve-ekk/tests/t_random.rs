//! Golden test for manifest row 3, `src/util/highs_random.rs`
//! (HighsRandom.h). Fixture: `random.jsonl` (committed), written by
//! `tools/highs-oracle/src/random.cpp`.
//!
//! For each seed the oracle constructs `HighsRandom r(seed)` and draws, in
//! this order: 16 x `integer()`, 64 x `integer(sup)` (sup from the record),
//! 16 x `fraction()`, 16 x `closedFraction()`, 16 x `real(-3, 5)`,
//! 16 x `bit()`, then `shuffle(data, 20)` of `0..20`. The test replays the
//! same draws on one generator, so a wrong draw count shifts every later
//! value: read the FIRST mismatch.

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::prelude::*;

#[test]
fn random_draw_sequences_match_highs() {
    let fx = load_required("random.jsonl");
    let mut n = 0;
    for rec in recs_of(&fx, "random") {
        let r = RandomRec::from_rec(&rec);
        let ctx = format!("{} seed={} (HighsRandom)", rec.ctx, r.seed);
        let mut g = HighsRandom::new(r.seed);
        let mut chk = Checker::new(&ctx);

        let integer: Vec<i32> = (0..r.integer.len()).map(|_| g.integer()).collect();
        chk.ints("integer() x16 [HighsRandom::integer]", &r.integer, &integer);

        let integer_sup: Vec<i32> = r.sup.iter().map(|&s| g.integer_sup(s)).collect();
        chk.ints(
            "integer(sup) x64 [HighsRandom::integer(HighsInt sup)]",
            &r.integer_sup,
            &integer_sup,
        );

        let fraction: Vec<f64> = (0..r.fraction.len()).map(|_| g.fraction()).collect();
        chk.f64s(
            "fraction() x16 [HighsRandom::fraction]",
            &r.fraction,
            &fraction,
        );

        let closed: Vec<f64> = (0..r.closed_fraction.len())
            .map(|_| g.closed_fraction())
            .collect();
        chk.f64s(
            "closedFraction() x16 [HighsRandom::closedFraction]",
            &r.closed_fraction,
            &closed,
        );

        let real: Vec<f64> = (0..r.real_m3_5.len()).map(|_| g.real(-3.0, 5.0)).collect();
        chk.f64s("real(-3,5) x16 [HighsRandom::real]", &r.real_m3_5, &real);

        let bit: Vec<i32> = (0..r.bit.len()).map(|_| g.bit() as i32).collect();
        chk.ints("bit() x16 [HighsRandom::bit]", &r.bit, &bit);

        let mut data: Vec<i32> = (0..20).collect();
        g.shuffle(&mut data, 20);
        chk.ints(
            "shuffle(0..20, 20) [HighsRandom::shuffle]",
            &r.shuffle20,
            &data,
        );

        chk.finish();
        n += 1;
    }
    assert!(n >= 4, "random.jsonl: expected >= 4 seeds, found {n}");
}

/// A generator re-initialised with `initialise(seed)` must restart the
/// stream (HiGHS: `HighsRandom::initialise`, used by `setSimplexOptions`).
#[test]
fn random_initialise_restarts_stream() {
    let fx = load_required("random.jsonl");
    for rec in recs_of(&fx, "random") {
        let r = RandomRec::from_rec(&rec);
        let mut g = HighsRandom::new(r.seed.wrapping_add(1));
        let _ = g.integer();
        g.initialise(r.seed);
        let integer: Vec<i32> = (0..r.integer.len()).map(|_| g.integer()).collect();
        let mut chk = Checker::new(format!(
            "{} seed={} after initialise (HighsRandom::initialise)",
            rec.ctx, r.seed
        ));
        chk.ints("integer() x16", &r.integer, &integer);
        chk.finish();
    }
}
