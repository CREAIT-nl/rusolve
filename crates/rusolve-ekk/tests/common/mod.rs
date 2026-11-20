//! Shared support for the `rusolve-ekk` golden-fixture tests.
//!
//! The fixtures are dumped from real HiGHS 1.15.1 internals by
//! `tools/highs-oracle` (see its README.md and FORMAT.md). This module reads
//! them and compares a Rust port against them **bit for bit**. It does not
//! depend on the `rusolve-ekk` API: tests convert their actual values into
//! the plain views/records defined here.
//!
//! Layout:
//! - [`fixture`]: MANIFEST.tsv, locating committed vs target-only files,
//!   loading JSONL, and [`fixture::Rec`], a context-carrying record accessor.
//! - [`records`]: typed decoders for the composite objects (HVector, matrix,
//!   basis, work arrays, LP, factor results, INVERT arrays, solves).
//! - [`hash`]: the oracle's FNV-1a 64 array hashes.
//! - [`emit`]: a JSON builder producing records in exactly the oracle's
//!   encoding, plus HVector / matrix views.
//! - [`events`]: typed simplex-trace records mirroring `simplex.cpp`
//!   (`from_rec` decodes a fixture record, `to_obj` emits the port's).
//! - [`check`]: bit-exact comparison helpers ([`check::Checker`]) and a
//!   generic, encoding-aware JSON differ.
//! - [`trace`]: simplex-trace summaries and the first-divergence reporter.
//! - [`sha256`]: used to verify fixtures against MANIFEST.tsv.
//!
//! Use from an integration test with `mod common;` then `use common::*;`.
#![allow(dead_code, unused_imports)]

pub mod check;
pub mod emit;
pub mod events;
pub mod fixture;
pub mod hash;
pub mod records;
pub mod sha256;
pub mod trace;

pub use check::{
    Checker, Mismatch, assert_f64_bits, assert_f64_slice, assert_i8_slice, assert_int_slice,
    diff_values, fmt_f64,
};
pub use emit::{HVectorView, MatrixView, Obj, PackView};
pub use fixture::{Fixture, Location, Rec, load, load_required, manifest};
pub use hash::{fnv_f64, fnv_i8, fnv_i32};
pub use records::{
    BasisRec, BuildResultRec, HVectorRec, InvertRec, LpRec, MatrixRec, PackRec, RandomRec, SolveIo,
    SolveRec, WorkArraysRec,
};
pub use trace::{
    DivergenceKind, TraceCompare, TraceDivergence, TraceKind, TraceSummary, assert_trace_matches,
    committed_trace_files, compare_traces, summarize, trace_files, trace_lp_name,
};

/// Decode an oracle f64: exactly 16 lowercase hex digits of the IEEE-754
/// binary64 bit pattern (FORMAT.md, "Scalars").
pub fn hexf(s: &str) -> Result<f64, String> {
    hex64(s).map(f64::from_bits)
}

/// Decode a 16-lowercase-hex-digit u64 (f64 bit pattern or hash).
pub fn hex64(s: &str) -> Result<u64, String> {
    if s.len() != 16
        || !s
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(format!("expected 16 lowercase hex digits, found {s:?}"));
    }
    u64::from_str_radix(s, 16).map_err(|e| format!("bad hex {s:?}: {e}"))
}

/// Encode an f64 the way the oracle does (`f64::to_bits` as `{:016x}`).
pub fn f64_hex(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

/// `kHighsZero` (1e-50), which HiGHS stores in place of exact cancellations.
pub const K_HIGHS_ZERO_BITS: u64 = 0x358dee7a4ad4b81f;
