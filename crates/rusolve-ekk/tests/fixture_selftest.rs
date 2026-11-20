//! Self-test of the fixture support in `tests/common/` (no port code
//! involved): the manifest, the loader, the decoders, the re-encoders and the
//! hash scheme must agree with every committed fixture, so a golden-test
//! failure can only mean the port differs from HiGHS.
//!
//! Run: `cargo test -p rusolve-ekk --test fixture_selftest -- --nocapture`
//! (prints per-file record counts). With target fixtures generated
//! (`tools/highs-oracle/gen_fixtures.sh`), they are checked too.

mod common;

use common::events::{EkkHash, Event};
use common::fixture::{Location, locate};
use common::hash::{FNV_OFFSET, fnv1a_bytes};
use common::sha256::sha256_hex;
use common::*;
use serde_json::Value;

#[test]
fn manifest_pins_highs_commit() {
    let m = manifest();
    // Another oracle's fixture set (RUSOLVE_FIXTURES_DIR) has its own commit.
    if common::fixture::fixtures_root().is_none() {
        assert_eq!(m.highs_commit, common::fixture::PINNED_HIGHS_COMMIT);
    }
    assert_eq!(m.patch_sha256.len(), 64, "patch_sha256 header");
    assert!(!m.committed().is_empty());
}

#[test]
fn hex_encoding_round_trips() {
    for x in [
        0.0,
        -0.0,
        1.0,
        -1.5,
        1e-50,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::MIN_POSITIVE,
    ] {
        let s = f64_hex(x);
        assert_eq!(hexf(&s).unwrap().to_bits(), x.to_bits(), "{s}");
    }
    assert_eq!(1e-50f64.to_bits(), K_HIGHS_ZERO_BITS);
    assert_eq!(f64_hex(1.0), "3ff0000000000000");
    assert!(
        hexf("3FF0000000000000").is_err(),
        "uppercase is not the oracle encoding"
    );
    assert!(hexf("3ff000000000000").is_err(), "15 digits");
}

#[test]
fn fnv_and_sha256_known_vectors() {
    assert_eq!(fnv_f64(&[]), FNV_OFFSET);
    assert_eq!(fnv_i32(&[]), FNV_OFFSET);
    // Standard FNV-1a 64 test vectors.
    assert_eq!(fnv1a_bytes(FNV_OFFSET, b"a"), 0xaf63dc4c8601ec8c);
    assert_eq!(fnv1a_bytes(FNV_OFFSET, b"foobar"), 0x85944171f73967e8);
    // Element encodings: little-endian bytes of the bit pattern.
    assert_eq!(fnv_i32(&[0x64636261]), fnv1a_bytes(FNV_OFFSET, b"abcd"));
    assert_eq!(fnv_i8(&[b'a' as i8]), 0xaf63dc4c8601ec8c);
    assert_eq!(
        fnv_f64(&[1.0]),
        fnv1a_bytes(FNV_OFFSET, &1.0f64.to_bits().to_le_bytes())
    );
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    let long = vec![b'a'; 1000];
    assert_eq!(
        sha256_hex(&long),
        "41edece42d63e8d9bf515a9ba6932e1c20cbc9f5a5d134645adb5db1b9737ea3"
    );
}

/// Every fixture that is present (all committed ones, plus the target-only
/// ones that have been generated), loaded and sha256-verified.
fn present_fixtures() -> Vec<Fixture> {
    let mut out = Vec::new();
    let mut missing = Vec::new();
    for e in &manifest().entries {
        match locate(&e.file) {
            Ok(_) => out.push(load_required(&e.file)),
            Err(_) => missing.push(e.file.clone()),
        }
    }
    if !missing.is_empty() {
        eprintln!(
            "NOTE: {} target-only fixture(s) not generated (run tools/highs-oracle/gen_fixtures.sh); not checked: {}",
            missing.len(),
            missing.join(", ")
        );
    }
    out
}

#[test]
fn every_committed_fixture_loads() {
    let committed = manifest().committed();
    let mut total = 0;
    for e in &committed {
        let fx = load_required(&e.file);
        assert_eq!(fx.location, Location::Committed);
        assert!(!fx.is_empty(), "{} is empty", e.file);
        total += fx.len();
        let counts: Vec<String> = fx
            .type_counts()
            .iter()
            .map(|(t, c)| format!("{t}={c}"))
            .collect();
        eprintln!("{:40} {:6} records  {}", e.file, fx.len(), counts.join(" "));
    }
    eprintln!("{} committed fixtures, {total} records", committed.len());
}

/// Decode every record with the typed decoders and re-encode it with the
/// builders; the result must equal the fixture record exactly.
#[test]
fn decoders_and_builders_round_trip() {
    let mut checked = 0usize;
    for fx in present_fixtures() {
        for r in fx.recs() {
            let rebuilt = rebuild(&r);
            if let Some(v) = rebuilt {
                let d = diff_values(&r.ctx, r.v, &v);
                assert!(
                    d.is_empty(),
                    "round trip failed:\n{}",
                    d.iter()
                        .map(|m| m.to_string())
                        .collect::<Vec<_>>()
                        .join("\n")
                );
                checked += 1;
            }
        }
    }
    eprintln!("{checked} records round-tripped");
}

/// Re-encode one record from its decoded form (typed events for traces;
/// every nested HVector / matrix for the other fixtures, after running the
/// typed decoders over the record so they are exercised too).
fn rebuild(r: &Rec) -> Option<Value> {
    let file = r.ctx.split(':').next().unwrap_or("");
    if file.starts_with("simplex/") {
        return Some(Event::from_rec(r).to_value());
    }
    // BASICLU records share type names with the HFactor ones (`update_step`)
    // but not their shape; they hold no HVector / matrix objects, so only the
    // generic re-encoding applies.
    if file.starts_with("basiclu/") {
        return Some(reencode(r.v, &r.ctx));
    }
    match r.t() {
        "random" => {
            let _ = RandomRec::from_rec(r);
        }
        "build" | "refactor" => {
            let _ = r.build_result("result");
            let _ = r.solves("solves");
        }
        "update_step" => {
            let _ = r.solves("solves");
        }
        "update_sequence_end" => {
            let _ = r.invert("invert");
        }
        _ => {}
    }
    Some(reencode(r.v, &r.ctx))
}

fn reencode(v: &Value, ctx: &str) -> Value {
    match v {
        Value::Object(m) if m.contains_key("nz_i") && m.contains_key("size") => {
            let rec = HVectorRec::from_rec(&Rec::new(v, ctx.to_string()));
            Obj::new().hv_rec("x", &rec).build()["x"].clone()
        }
        // HighsHessian objects carry `dim` instead of `num_col`/`num_row`
        // and fall through to the generic re-encoding.
        Value::Object(m)
            if m.contains_key("format") && m.contains_key("start") && m.contains_key("num_col") =>
        {
            let rec = MatrixRec::from_rec(&Rec::new(v, ctx.to_string()));
            Obj::new().mat("x", &rec.view()).build()["x"].clone()
        }
        Value::Object(m) => Value::Object(
            m.iter()
                .map(|(k, x)| (k.clone(), reencode(x, &format!("{ctx}.{k}"))))
                .collect(),
        ),
        Value::Array(a) => Value::Array(
            a.iter()
                .enumerate()
                .map(|(i, x)| reencode(x, &format!("{ctx}[{i}]")))
                .collect(),
        ),
        _ => v.clone(),
    }
}

/// Hash identities that tie the trace records together (checked on every
/// present trace): they prove the FNV scheme matches the oracle's, and they
/// document two HiGHS facts the port must reproduce.
#[test]
fn trace_hash_identities() {
    for fx in present_fixtures()
        .into_iter()
        .filter(|f| f.file.starts_with("simplex/"))
    {
        // 1. INVERT of the initial (logical) basis: h_basicIndex of the first
        //    `invert` is the FNV of the basis at ekk_solve_initialised.
        let init = fx.first("ekk_solve_initialised");
        let basis = init.basis("basis");
        let inv = fx.first("invert");
        assert_eq!(
            inv.hash("h_basicIndex"),
            fnv_i32(&basis.basic_index),
            "{}: first invert h_basicIndex != fnv(initialised basicIndex)",
            fx.file
        );
        // Identities 2 and 3 hold per HEkk::solve call: a trace has one, or
        // two when HApp re-solves the unscaled LP.
        let mut solve_start = 0;
        for (j, v) in fx.records.iter().enumerate() {
            if v["t"] != "ekk_solve_exit" {
                continue;
            }
            let exit = fx.rec(j);
            let mut this_solve = (solve_start..j).map(|i| fx.rec(i));
            // 2. After an Optimal exit, the last cleanup's hash equals the
            //    exit arrays, except workDual: HEkk::returnFromSolve calls
            //    zeroBasicDuals() before the exit hook. (Not after other
            //    exits: e.g. adlittlecut's primal solve changes workCost after
            //    its cleanup before ending Infeasible.)
            let ex_basis = exit.basis("basis");
            let ex = EkkHash::of(
                &ex_basis.basic_index,
                &ex_basis.nonbasic_move,
                &exit.arrays("arrays"),
            );
            let cleanup = this_solve.rfind(|r| r.t() == "cleanup");
            if let Some(cleanup) = cleanup.filter(|_| exit.int("model_status") == 7) {
                let mut h = EkkHash::from_rec(&cleanup.obj("hash"));
                // (They may also agree, when the basic duals were already 0.)
                h.work_dual = ex.work_dual;
                assert_eq!(
                    h, ex,
                    "{}: last cleanup hash vs ekk_solve_exit arrays",
                    exit.ctx
                );
            }
            // 3. The exit basicIndex is the latest one dumped: by an `iter`
            //    (HEkk::updatePivots) or an `invert` (HFactor::build permutes
            //    it, e.g. kb2 reinverts after its last iteration).
            let latest = (0..j)
                .map(|i| fx.rec(i))
                .rfind(|r| matches!(r.t(), "iter" | "invert") && r.has("h_basicIndex"))
                .map(|r| (r.ctx.clone(), r.hash("h_basicIndex")));
            if let Some((ctx, h)) = latest {
                assert_eq!(
                    h, ex.basic_index,
                    "{ctx}: latest h_basicIndex vs ekk_solve_exit basis"
                );
            }
            solve_start = j + 1;
        }
    }
}

/// The comparison helpers report divergences rather than passing them.
#[test]
fn trace_compare_detects_divergence() {
    let fx = load_required("simplex/afiro.dual.jsonl");
    assert!(compare_traces(&fx, &fx.records, &TraceCompare::default()).is_ok());

    let mut bad = fx.records.clone();
    let i = fx.records.iter().position(|v| v["t"] == "iter").unwrap();
    bad[i]["alpha_col"] = Value::String(f64_hex(12345.0));
    let d = compare_traces(&fx, &bad, &TraceCompare::default()).unwrap_err();
    assert_eq!(d.kind, DivergenceKind::Fields);
    assert_eq!(d.line, Some(i + 1));
    let msg = d.to_string();
    assert!(msg.contains("alpha_col") && msg.contains("ftran"), "{msg}");
    assert!(
        compare_traces(&fx, &bad, &TraceCompare::default().skip("iter/alpha_col")).is_ok(),
        "skip spec"
    );

    let short = &fx.records[..fx.len() - 1];
    let d = compare_traces(&fx, short, &TraceCompare::default()).unwrap_err();
    assert_eq!(d.kind, DivergenceKind::ActualEnded);

    // A rebuild where HiGHS iterated: control-flow divergence.
    let mut swapped = fx.records.clone();
    let rebuild = fx
        .records
        .iter()
        .find(|v| v["t"] == "rebuild")
        .unwrap()
        .clone();
    swapped.insert(i, rebuild);
    let d = compare_traces(&fx, &swapped, &TraceCompare::default()).unwrap_err();
    assert_eq!(d.kind, DivergenceKind::Type);
    eprintln!("{d}");
}

#[test]
fn trace_summaries() {
    for kind in [TraceKind::Dual, TraceKind::Primal] {
        for f in committed_trace_files(kind) {
            let fx = load_required(&f);
            let s = summarize(&fx);
            // The `<lp>cut` variants (tools/highs-oracle/mkvariants.py) are
            // primal infeasible; every other committed trace is optimal.
            let expected = if trace_lp_name(&f).ends_with("cut") {
                8
            } else {
                7
            };
            assert_eq!(
                s.model_status,
                Some(expected),
                "{f}: committed trace model status"
            );
            eprintln!("{s}");
        }
    }
}

#[test]
fn checker_reports_bits() {
    let mut c = Checker::new("ctx");
    c.f64("x", 1.0, 1.0);
    c.f64("y", 0.0, -0.0);
    c.ints("v", &[1, 2, 3], &[1, 2]);
    assert_eq!(c.failures.len(), 2);
    let s = c.failures[0].to_string();
    assert!(
        s.contains("0000000000000000") && s.contains("8000000000000000"),
        "{s}"
    );
}
