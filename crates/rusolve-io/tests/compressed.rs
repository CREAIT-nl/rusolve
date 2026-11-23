//! Compressed model files read exactly like their plain originals.
//!
//! Fixtures in `tests/compressed/` were made from `tests/dialects/` with the
//! reference command-line tools (gzip 1.x `-n`, bzip2 1.0.8, xz 5.8.4,
//! zstd 1.5.7), so they test the decoders on what users actually produce:
//! * one file per wrapper and inner format;
//! * `cplex_bz2_named.mps.gz`: bzip2 content under a `.gz` name (magic
//!   bytes, not the name, choose the decoder);
//! * `pulp_lp_noext` / `xpress_mps_noext`: no extension at all, inner
//!   format sniffed from the decoded text;
//! * `scip_two_members.mps.gz` / `scip_two_frames.mps.zst`: the file split
//!   in two halves compressed separately and concatenated, which `gzip -d`
//!   and `zstd -d` decode to the whole file.

use rusolve_io::{Compression, read_file_bytes, read_model_file};

fn path(dir: &str, f: &str) -> String {
    format!("{}/tests/{dir}/{f}", env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn every_wrapper_decodes_to_the_plain_model() {
    let cases = [
        ("highs.mps.gz", "highs.mps"),
        ("highs.mps.bz2", "highs.mps"),
        ("highs.mps.xz", "highs.mps"),
        ("highs.mps.zst", "highs.mps"),
        ("scip.lp.bz2", "scip.lp"),
        ("gurobi.lp.xz", "gurobi.lp"),
        ("copt.lp.zst", "copt.lp"),
        ("cplex_bz2_named.mps.gz", "cplex.mps"),
        ("pulp_lp_noext", "pulp.lp"),
        ("xpress_mps_noext", "xpress.mps"),
        ("scip_two_members.mps.gz", "scip.mps"),
        ("scip_two_frames.mps.zst", "scip.mps"),
    ];
    for (packed, plain) in cases {
        let want = std::fs::read(path("dialects", plain)).unwrap();
        let got = read_file_bytes(&path("compressed", packed)).unwrap();
        assert_eq!(got, want, "{packed}: decoded bytes");
        let m = read_model_file(&path("compressed", packed)).unwrap();
        let m0 = read_model_file(&path("dialects", plain)).unwrap();
        assert_eq!(m, m0, "{packed}: model");
    }
}

#[test]
fn magic_bytes_pick_the_wrapper() {
    let sniff = |f: &str| Compression::sniff(&std::fs::read(path("compressed", f)).unwrap());
    assert_eq!(sniff("highs.mps.gz"), Compression::Gzip);
    assert_eq!(sniff("highs.mps.bz2"), Compression::Bzip2);
    assert_eq!(sniff("highs.mps.xz"), Compression::Xz);
    assert_eq!(sniff("highs.mps.zst"), Compression::Zstd);
    assert_eq!(sniff("cplex_bz2_named.mps.gz"), Compression::Bzip2);
    assert_eq!(
        Compression::sniff(&std::fs::read(path("dialects", "highs.mps")).unwrap()),
        Compression::None
    );
    // "BZh" must be followed by a block-size digit: a plain file that
    // happens to start with those letters is not bzip2.
    assert_eq!(Compression::sniff(b"BZhello"), Compression::None);
}

#[test]
fn truncated_archives_are_errors_not_short_models() {
    let dir = std::env::temp_dir().join(format!("rusolve-io-trunc-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for f in [
        "highs.mps.gz",
        "highs.mps.bz2",
        "highs.mps.xz",
        "highs.mps.zst",
    ] {
        let bytes = std::fs::read(path("compressed", f)).unwrap();
        let cut = dir.join(f);
        std::fs::write(&cut, &bytes[..bytes.len() * 2 / 3]).unwrap();
        let r = read_file_bytes(&cut.display().to_string());
        assert!(
            r.is_err(),
            "{f}: truncated archive decoded to {:?} bytes",
            r.map(|b| b.len())
        );
    }
}
