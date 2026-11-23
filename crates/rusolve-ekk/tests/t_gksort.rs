//! Differential test for the GKlib sorts behind libmetis `isortd` and
//! `ikvsorti` (src/ipm/hipo/orderings/metis/gklib.rs).
//!
//! GKlib's GK_MKQSORT is LGPL-2.1+, so the Rust sort is an independent
//! implementation written from a description of its ordering, not a translation. This
//! test is the evidence that it orders ties exactly like GK_MKQSORT: METIS
//! orderings, and hence HiPO's factorisation, depend on the order of equal
//! keys. The expected outputs were produced by compiling the vendored
//! GK_MKQSORT in tools/highs-oracle/gksort (regenerate with its gen.sh);
//! the fixture holds data only.
//!
//! Lines: `d n in.. | out..` (isortd), `kv n keys.. | vals..` (ikvsorti with
//! val_i = i on input), `kvh n seed range | fnv1a64` (large ikvsorti; keys
//! from the LCG below, hash over the sorted (key, val) pairs as i32 LE).

use rusolve_ekk::ipm::hipo::orderings::metis::gklib::{ikv_t, ikvsorti, isortd};

fn fixture() -> String {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/gksort/gksort.txt"
    );
    std::fs::read_to_string(path).expect("missing tests/fixtures/gksort/gksort.txt")
}

fn ints(s: &str) -> Vec<i32> {
    s.split_whitespace().map(|t| t.parse().unwrap()).collect()
}

fn split(rest: &str) -> (usize, Vec<i32>, &str) {
    let (lhs, rhs) = rest.split_once('|').expect("missing '|'");
    let mut it = lhs.split_whitespace();
    let n: usize = it.next().unwrap().parse().unwrap();
    let inp: Vec<i32> = it.map(|t| t.parse().unwrap()).collect();
    (n, inp, rhs)
}

#[test]
fn isortd_matches_gk_mkqsort() {
    let mut cases = 0;
    for line in fixture().lines() {
        let Some(rest) = line.strip_prefix("d ") else {
            continue;
        };
        let (n, inp, rhs) = split(rest);
        let want = ints(rhs);
        assert_eq!(inp.len(), n);
        let mut got = inp.clone();
        isortd(n, &mut got);
        assert_eq!(got, want, "isortd on {inp:?}");
        cases += 1;
    }
    assert_eq!(cases, 120);
}

#[test]
fn ikvsorti_matches_gk_mkqsort_tie_order() {
    let mut cases = 0;
    for line in fixture().lines() {
        let Some(rest) = line.strip_prefix("kv ") else {
            continue;
        };
        let (n, keys, rhs) = split(rest);
        let want = ints(rhs);
        let mut a: Vec<ikv_t> = keys
            .iter()
            .enumerate()
            .map(|(i, &key)| ikv_t { key, val: i as i32 })
            .collect();
        ikvsorti(n, &mut a);
        let vals: Vec<i32> = a.iter().map(|e| e.val).collect();
        assert_eq!(vals, want, "ikvsorti on keys {keys:?}");
        cases += 1;
    }
    assert_eq!(cases, 200);
}

#[test]
fn ikvsorti_large_matches_gk_mkqsort_hash() {
    let mut cases = 0;
    for line in fixture().lines() {
        let Some(rest) = line.strip_prefix("kvh ") else {
            continue;
        };
        let (lhs, rhs) = rest.split_once('|').unwrap();
        let p: Vec<u64> = lhs.split_whitespace().map(|t| t.parse().unwrap()).collect();
        let (n, seed, range) = (p[0] as usize, p[1] as u32, p[2] as u32);
        let want: u64 = rhs.trim().parse().unwrap();
        let mut lcg = seed;
        let mut a: Vec<ikv_t> = (0..n)
            .map(|i| {
                lcg = lcg.wrapping_mul(1103515245).wrapping_add(12345);
                ikv_t {
                    key: ((lcg >> 8) % range) as i32,
                    val: i as i32,
                }
            })
            .collect();
        ikvsorti(n, &mut a);
        let mut h: u64 = 1469598103934665603;
        for e in &a {
            for b in e.key.to_le_bytes().into_iter().chain(e.val.to_le_bytes()) {
                h ^= u64::from(b);
                h = h.wrapping_mul(1099511628211);
            }
        }
        assert_eq!(h, want, "ikvsorti n={n} seed={seed} range={range}");
        cases += 1;
    }
    assert_eq!(cases, 12);
}
