//! The contributor map stays true.
//!
//! `tools/gen_codemap.py --check` fails when `docs/agent/codemap.md` is stale (it is generated from
//! the `//! HiGHS source:` headers), when a crate or rusolve-ekk subsystem has no `AGENTS.md`, or
//! when an `AGENTS.md` names a repo path that does not exist. `tools/gen_agent_docs.py --check`
//! runs the same check (see `agent_docs.rs`); this file also proves the check can fail.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn codemap(args: &[&str]) -> Output {
    match Command::new("python3")
        .arg(root().join("tools/gen_codemap.py"))
        .args(args)
        .output()
    {
        Ok(o) => o,
        Err(e) => panic!("codemap: python3 is required to run tools/gen_codemap.py: {e}"),
    }
}

#[test]
fn codemap_is_current_and_agents_files_name_real_paths() {
    let out = codemap(&["--check"]);
    assert!(
        out.status.success(),
        "codemap or AGENTS.md is stale (run `python3 tools/gen_codemap.py`, fix the named path):\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn an_agents_file_naming_a_missing_path_fails_the_check() {
    let dir = std::env::temp_dir().join(format!("rusolve-codemap-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("AGENTS.md");
    let body =
        "# broken\n\nSee `crates/rusolve-ekk/src/simplex/does_not_exist.rs` for the driver.\n";
    // pad to the 30-60 line window so only the path is wrong
    let padded = format!("{body}{}", "line\n".repeat(32));
    std::fs::write(&file, padded).unwrap();
    let out = codemap(&["--check-agents", file.to_str().unwrap()]);
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    let _ = std::fs::remove_dir_all(&dir);
    assert!(!out.status.success(), "a missing path must fail the check");
    assert!(
        err.contains("does_not_exist.rs"),
        "the message names the path: {err}"
    );
}

#[test]
fn a_real_agents_file_passes_the_single_file_check() {
    let file = root().join("crates/rusolve-ekk/src/simplex/AGENTS.md");
    let out = codemap(&["--check-agents", file.to_str().unwrap()]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
