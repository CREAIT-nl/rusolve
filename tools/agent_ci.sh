#!/usr/bin/env bash
# CI check: the agent docs match the binary, every skill has a graded eval set, and every
# grader accepts its reference answer and rejects its wrong ones. Needs no LLM.
#   tools/agent_ci.sh            # builds rusolve-cli first
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build --release -p rusolve-cli
export RSV_CLI="${RSV_CLI:-target/release/rusolve-cli}"
python3 tools/gen_agent_docs.py --check
python3 evals/graders.py --selfcheck >/dev/null
cargo test --release -p rusolve-cli --test agent_docs
echo "agent_ci: ok"
