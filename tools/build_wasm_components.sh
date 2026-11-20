#!/bin/bash
# Builds the Rust Wasm hook components into crates/rusolve-wasm/components/*.wasm.
# Needs the wasm32-unknown-unknown target (user-local: `rustup target add wasm32-unknown-unknown`).
# Reproducible: fixed flags, stripped, path-remapped, no timestamps. The built .wasm is committed so
# tests need no target; `tools/build_wasm_components.sh --check` rebuilds and compares byte for byte.
#   tools/build_wasm_components.sh           # rebuild and overwrite the committed .wasm
#   tools/build_wasm_components.sh --check   # rebuild to a temp file and fail if it differs
set -e
export PATH="$HOME/.cargo/bin:$PATH"
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$ROOT/crates/rusolve-wasm/components/rust
OUT=$ROOT/crates/rusolve-wasm/components/hooks_rust.wasm
# Host data and stack go above the first 256 KiB: the host writes its input at address 0.
export RUSTFLAGS="-C link-arg=--global-base=262144 -C link-arg=--no-stack-first -C link-arg=-zstack-size=65536 -C link-arg=--initial-memory=1048576 -C link-arg=--max-memory=1048576 --remap-path-prefix=$SRC=/src --remap-path-prefix=$HOME=/home"
cd "$SRC"
CARGO_TARGET_DIR="$SRC/target" cargo build --release --target wasm32-unknown-unknown
BUILT=$SRC/target/wasm32-unknown-unknown/release/rusolve_wasm_components.wasm
if [ "$1" = "--check" ]; then
  cmp "$BUILT" "$OUT" && echo "hooks_rust.wasm: byte-identical rebuild ($(wc -c < "$OUT") bytes)"
else
  cp "$BUILT" "$OUT"
  echo "wrote $OUT ($(wc -c < "$OUT") bytes)"
fi
