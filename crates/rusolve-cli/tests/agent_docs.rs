//! The agent docs ship with the solver and stay true.
//!
//! `tools/gen_agent_docs.py --check` regenerates `docs/agent/reference.md` from this very binary
//! (`rusolve options --json`, `env --json`, the MCP `tools/list`, `--help --json`)
//! and fails when the committed file differs, and when `llms.txt`, `README.md` or a
//! `docs/skills/*/SKILL.md` names an MCP tool or subcommand that does not exist or links to a
//! missing file. Fix a failure with `python3 tools/gen_agent_docs.py` (and edit the docs).

use std::path::Path;
use std::process::Command;

#[test]
fn agent_docs_reference_is_current_and_skills_point_at_real_commands() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let out = match Command::new("python3")
        .arg(root.join("tools/gen_agent_docs.py"))
        .arg("--check")
        .env("RSV_CLI", env!("CARGO_BIN_EXE_rusolve-cli"))
        .output()
    {
        Ok(o) => o,
        Err(e) => panic!("agent_docs: python3 is required to run tools/gen_agent_docs.py: {e}"),
    };
    assert!(
        out.status.success(),
        "agent docs are stale or inconsistent:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn agent_docs_llms_txt_and_four_skills_exist() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    assert!(root.join("llms.txt").is_file());
    for skill in ["solve", "debug-infeasible", "tune", "extend"] {
        assert!(
            root.join("docs/skills")
                .join(skill)
                .join("SKILL.md")
                .is_file(),
            "docs/skills/{skill}/SKILL.md"
        );
    }
}
