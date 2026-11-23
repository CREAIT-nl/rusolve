#!/bin/bash
# Build the abi3 wheel that bundles the CLI.
# Usage: tools/build_wheel.sh [outdir]   (env: PROFILE=dist|release, default release; MATURIN=maturin)
# Steps: build rusolve-cli, copy it to crates/rusolve-py/python/rusolve/_bin/, maturin build.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:-dist-wheels}
PROFILE=${PROFILE:-release}
MATURIN=${MATURIN:-maturin}
cargo build --profile "$PROFILE" -p rusolve-cli
EXE=rusolve-cli; case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*) EXE=rusolve-cli.exe;; esac
mkdir -p crates/rusolve-py/python/rusolve/_bin
cp "target/$PROFILE/$EXE" crates/rusolve-py/python/rusolve/_bin/
chmod +x "crates/rusolve-py/python/rusolve/_bin/$EXE"
"$MATURIN" build --profile "$PROFILE" -m crates/rusolve-py/Cargo.toml -o "$OUT" ${MATURIN_ARGS:-}
