//! Captures the commit the binary
//! was built from, as `RUSOLVE_GIT` for `envelope::Solver::git`.
//!
//! `git rev-parse --short=12 HEAD`, plus `-dirty` when `git status --porcelain`
//! lists a change to a tracked file (untracked files do not change a build, so
//! they are left out with `-uno`). When git or the repository is not available
//! (a source tarball, a vendored build) nothing is emitted and the envelope omits
//! `solver.git`, which the schema allows.
//!
//! The script reruns when HEAD moves: it watches HEAD itself, the ref HEAD points
//! to, `packed-refs` and the index. The paths come from `git rev-parse
//! --git-path`, so a worktree (whose `.git` is a file pointing at a per-worktree
//! directory, with refs in the shared common directory) is handled by git itself.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn watch(git_path: &str) {
    if let Some(p) = git(&["rev-parse", "--git-path", git_path]) {
        println!("cargo:rerun-if-changed={p}");
    }
}

/// The skills under `docs/skills/<name>/SKILL.md`
/// become MCP prompts and the schemas under `docs/schema/*.json` become MCP resources. Both
/// directories are listed here at build time, so a skill or schema added later shows up in the
/// next build without touching Rust (`$OUT_DIR/mcp_docs.rs`, included by `src/mcp_resources.rs`).
/// A variable cargo sets for build scripts. Read through `env::vars()` on purpose: the registry
/// test (`src/env_registry.rs`) lists the variables the solver reads, and cargo's own are not
/// part of that surface.
fn cargo_env(name: &str) -> String {
    std::env::vars()
        .find(|(k, _)| k == name)
        .map(|(_, v)| v)
        .unwrap_or_else(|| panic!("cargo did not set {name}"))
}

fn embed_mcp_docs() {
    use std::fmt::Write;
    let root = std::path::PathBuf::from(cargo_env("CARGO_MANIFEST_DIR")).join("../../docs");
    let list = |sub: &str, want: &dyn Fn(&std::path::Path) -> Option<std::path::PathBuf>| {
        let dir = root.join(sub);
        println!("cargo:rerun-if-changed={}", dir.display());
        let mut found: Vec<(String, std::path::PathBuf)> = std::fs::read_dir(&dir)
            .map(|rd| {
                rd.filter_map(|e| e.ok())
                    .filter_map(|e| {
                        let p = e.path();
                        let file = want(&p)?;
                        let name = p.file_name()?.to_string_lossy().into_owned();
                        Some((name, file))
                    })
                    .collect()
            })
            .unwrap_or_default();
        found.sort();
        found
    };
    // A skill whose front matter says `audience: contributor` (the contrib-* skills) is for
    // people changing the solver, not for agents using it; it is not an MCP prompt.
    let skills = list("skills", &|p| {
        let f = p.join("SKILL.md");
        let text = std::fs::read_to_string(&f).ok()?;
        let fm = text.strip_prefix("---\n")?.split("\n---").next()?;
        (!fm.lines().any(|l| l.trim() == "audience: contributor")).then_some(f)
    });
    let schemas = list("schema", &|p| {
        (p.extension().is_some_and(|e| e == "json")).then(|| p.to_path_buf())
    });
    let mut out = String::new();
    for (table, items) in [("SKILLS", &skills), ("SCHEMAS", &schemas)] {
        writeln!(out, "pub const {table}: &[(&str, &str)] = &[").unwrap();
        for (name, file) in items {
            let file = file.canonicalize().unwrap();
            println!("cargo:rerun-if-changed={}", file.display());
            writeln!(
                out,
                "    ({name:?}, include_str!({:?})),",
                file.display().to_string()
            )
            .unwrap();
        }
        writeln!(out, "];").unwrap();
    }
    let reference = root.join("agent/reference.md");
    println!("cargo:rerun-if-changed={}", reference.display());
    writeln!(
        out,
        "pub const AGENT_REFERENCE: &str = include_str!({:?});",
        reference.canonicalize().unwrap().display().to_string()
    )
    .unwrap();
    let dest = std::path::PathBuf::from(cargo_env("OUT_DIR")).join("mcp_docs.rs");
    std::fs::write(dest, out).unwrap();
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    embed_mcp_docs();
    let Some(hash) = git(&["rev-parse", "--short=12", "HEAD"]).filter(|h| !h.is_empty()) else {
        return;
    };
    watch("HEAD");
    watch("index");
    watch("packed-refs");
    // `HEAD` is a symbolic ref on a branch (watch the ref file), a bare hash when
    // detached (HEAD itself is then the file that moves).
    if let Some(head) = git(&["symbolic-ref", "-q", "HEAD"]) {
        watch(&head);
    }
    let dirty = git(&["status", "--porcelain", "-uno"]).is_some_and(|s| !s.is_empty());
    println!(
        "cargo:rustc-env=RUSOLVE_GIT={hash}{}",
        if dirty { "-dirty" } else { "" }
    );
}
