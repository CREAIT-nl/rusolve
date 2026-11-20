#!/usr/bin/env bash
# Build the instrumented HiGHS oracle and regenerate the rusolve-ekk golden
# fixtures. See README.md (what/why) and FORMAT.md (encoding).
#
#   tools/highs-oracle/gen_fixtures.sh                 # build + regenerate
#   tools/highs-oracle/gen_fixtures.sh --determinism   # also: regenerate twice, byte-compare
#   tools/highs-oracle/gen_fixtures.sh --w5            # regenerate only w5/ (the rest of
#                                                      # target/highs-fixtures must match MANIFEST.tsv)
#   tools/highs-oracle/gen_fixtures.sh --fp-contract-check
#                                                      # also: compare against a -ffp-contract=on build
#   tools/highs-oracle/gen_fixtures.sh --mip           # regenerate only mip/ (traces + .ac at
#                                                      # threads=1) into the target fixture dir;
#                                                      # needs benches/data/miplib, no MANIFEST change
#   ORACLE_HIGHS=1.16.0 tools/highs-oracle/gen_fixtures.sh
#                                                      # the HiGHS 1.16.0 oracle, side by side (below)
#
# Environment overrides:
#   ORACLE_HIGHS HiGHS version to build  (default 1.15.1). 1.16.0 uses the patch set in
#               patches-1.16.0/, builds in target/highs-oracle-1.16.0/, writes every fixture
#               (and its MANIFEST.tsv) to target/highs-fixtures-1.16.0/ and never touches the
#               committed fixtures (crates/rusolve-ekk/tests/fixtures/).
#   HIGHS_REF   pristine HiGHS checkout   (default: <repo>/reference/HiGHS, else the main checkout's;
#               1.16.0: reference/HiGHS-1.16.0, else target/ref/HiGHS-1.16.0, in either checkout)
#   LP47_DIR    lp47 MPS directory         (default: <repo>/benches/data/lp47, else the main checkout's)
#   JOBS        parallel build jobs        (default: all cores; keep it low on a RAM-tight box)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
# reference/ and benches/data/lp47 are gitignored and live only in the main
# checkout; a git worktree finds it through the shared git common dir.
MAIN="$(cd "$(git -C "$ROOT" rev-parse --path-format=absolute --git-common-dir)/.." && pwd)"

pick_dir() {  # pick_dir <relative path> -> first existing of $ROOT/<p>, $MAIN/<p>
  if [ -d "$ROOT/$1" ]; then echo "$ROOT/$1"; else echo "$MAIN/$1"; fi
}
ORACLE_HIGHS="${ORACLE_HIGHS:-1.15.1}"
case "$ORACLE_HIGHS" in
  1.15.1)
    HIGHS_REF="${HIGHS_REF:-$(pick_dir reference/HiGHS)}"
    PATCH_DIR="$HERE"
    VSUFFIX="" ;;
  1.16.0)
    if [ -z "${HIGHS_REF:-}" ]; then
      HIGHS_REF="$(pick_dir reference/HiGHS-1.16.0)"
      [ -d "$HIGHS_REF" ] || HIGHS_REF="$(pick_dir target/ref/HiGHS-1.16.0)"
    fi
    PATCH_DIR="$HERE/patches-1.16.0"
    VSUFFIX="-1.16.0" ;;
  *) printf '[gen_fixtures] ORACLE_HIGHS=%s: expected 1.15.1 or 1.16.0\n' "$ORACLE_HIGHS" >&2; exit 64 ;;
esac
LP47_DIR="${LP47_DIR:-$(pick_dir benches/data/lp47)}"
NETLIB_DIR="$ROOT/benches/data/netlib"
# getconf works on macOS and Linux (sysctl hw.ncpu is BSD-only, nproc GNU-only).
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || nproc)}"

WORK="$ROOT/target/highs-oracle$VSUFFIX"
SRC="$WORK/highs-src"
PATCH="$PATCH_DIR/highs-oracle-hooks.patch"
# Second patch (IPX per-iteration hooks, oracle (b)), applied after $PATCH
# with its own stamp so that $PATCH, its hash and the MANIFEST header stay
# unchanged.
IPX_PATCH="$PATCH_DIR/highs-oracle-ipx-hooks.patch"
# Third patch (QUASS per-pass hooks, W1 fixtures), applied after $IPX_PATCH with its own stamp, so that the first two
# patch hashes and their fixtures stay unchanged.
W1_PATCH="$PATCH_DIR/highs-oracle-w1-hooks.patch"
# Fourth patch (cuPDLP-C / HiPDLP trace hooks, W3 fixtures), applied after $W1_PATCH with its own stamp, so that the first three
# patch hashes and their fixtures stay unchanged.
W3_PATCH="$PATCH_DIR/highs-oracle-w3-hooks.patch"
# Fifth patch (HiPO trace hooks, W4 fixtures),
# applied after $W3_PATCH with its own stamp, so that the first four patch
# hashes and their fixtures stay unchanged. The hooks are inert unless
# HIGHS_ORACLE_HIPO_DUMP is set, so the default build's output is unchanged.
W4_PATCH="$PATCH_DIR/highs-oracle-w4-hooks.patch"
# Sixth patch (W5 fixtures): hipo_ord records are
# buffered per ordering and written after chooseOrdering's for_each, in index
# order, so they do not race at threads > 1. Applied after $W4_PATCH with its
# own stamp; at one thread the output is byte-identical to patch 5's.
W5_PATCH="$PATCH_DIR/highs-oracle-w5-hooks.patch"
SMALL="$ROOT/crates/rusolve-ekk/tests/fixtures"
LARGE="$ROOT/target/highs-fixtures$VSUFFIX"
VARIANT_DIR="$WORK/variants"

# ---------------------------------------------------------------- fixture sets
# Unit (HighsSparseMatrix / HVector) and HFactor fixtures.
UNIT_LPS="afiro sc50a kb2 adlittle blend share2b"
FACTOR_LPS="afiro sc50a kb2 adlittle blend share2b"
# BASICLU (IPX LU kernel) fixtures.
BASICLU_LPS="afiro sc50a kb2 adlittle blend share2b"
BASICLU_SYNTH="dense70 rand120 neardep40"
# IPX fixtures: the per-iteration trace, ipx::Model and the IPX
# linear-algebra component dumps. Synthetic LPs in $SYNTH_DIR: ipxmix (free / boxed /
# fixed / upper-only / ranged / badly scaled), ipxtall (dualized by the
# default strategy) and ipxdep (dependent rows, duplicate columns). Configurations are <tag>:<overrides> ("default" = none).
SYNTH_DIR="$HERE/lps"
# ipx-model: every netlib LP + the synthetic ones with the default
# parameters, and forced dualization for a few.
IPX_MODEL_DUAL_LPS="afiro sc50a kb2 ipxmix"
# ipx-linalg: ipx::Parameters overrides; lu1 = ForrestTomlin kernel and the
# sequential maxvolume heuristic, slack = slack starting basis.
IPX_LINALG_LPS="afiro sc50a kb2 adlittle blend share2b ipxmix ipxtall ipxdep"
IPX_LINALG_CFGS="default:default lu1:lu_kernel=1,update_heuristic=0 dual:dualize=1 slack:crash_basis=0"
# ipx (oracle (c), every LP): end-to-end Highs::run with solver=ipx in the
# configurations xon/xoff/ipm/dual (src/ipx.cpp runIpx).
# ipx-trace: HiGHS option overrides (solver=ipx, presolve off).
IPX_TRACE_CFGS="default:default dual:ipx_dualize_strategy=1"
# Simplex traces: every netlib instance shipped in benches/data/netlib, plus
# lp47 LP relaxations that solve in well under 2 s with presolve off.
LP47_LPS="beasleyC3 binkar10_1 mc11 milo-v12-6-r2-40-1"
# Non-optimal variants written by mkvariants.py into $VARIANT_DIR:
# <lp>cut (primal infeasible) and <lp>max (maximised; some unbounded).
VARIANT_LPS="afiro adlittle sc50a sc105 blend kb2 share2b"
# Strategies: 1 = serial dual (kSimplexStrategyDual), 4 = primal.
STRATEGIES="1 4"
# W1 fixtures (all target-only).
# HiGHS check instances are read from $HIGHS_REF/check/instances, the
# synthetic IIS LPs from $HERE/lps/w1 (outside lps/*.mps, which the IPX loop
# reads). QP: w1-qp; ranging and iCrash: netlib; IIS: check + synthetic.
CHECK_DIR="$HIGHS_REF/check/instances"
W1_SYNTH_DIR="$HERE/lps/w1"
W1_QP="qjh qjh_uncon qjh_quadobj qjh_qmatrix 2821 2894 atwood0 primal1"
W1_RANGING="afiro sc50a kb2 adlittle blend share2b sc105 scagr7"
W1_IIS="galenet forest6 woodinfe klein1 iis_rows iis_bnds"
W1_ICRASH="afiro sc50a adlittle"
# W3 fixtures (all target-only): the
# same netlib LPs through cuPDLP-C (w3-cupdlp) and HiPDLP (w3-hipdlp).
W3_LPS="afiro sc50a kb2 adlittle blend share2b sc105 scagr7"
# W4 fixtures (all target-only): HiPO
# on the same netlib LPs, written by the second oracle configuration
# ($WORK/build-hipo: HIPO=ON, BLAS pinned to refblas/refblas.c), once per
# forced ordering (w4-orderings) and once with HiPO's own choices (w4-hipo).
W4_LPS="afiro sc50a kb2 adlittle blend share2b sc105 scagr7"
W4_ORDERINGS="metis amd rcm"
# W5 fixtures (all target-only): the
# W1 netlib LP set (= W3_LPS / W4_LPS) through SIP (simplex_strategy 2) and
# PAMI (3) on the default build, and through HiPO on build-hipo once per
# hipo_parallel_type and Newton system (as / ne at every thread count,
# choose at threads 1 only); the MIP trace (mip_trace, default build) on the
# instances of crates/rusolve-ekk/tests/fixtures/mip, read from
# benches/data/miplib (gitignored; this checkout, else the main one).
W5_LPS="$W4_LPS"
W5_STRATEGIES="2 3"
W5_THREADS="1 2 4 8"
W5_HIPO_TYPES="tree node both"
W5_MIP="air05 enlight_hard pg sp150x300d"
W5_MIP_THREADS="2 4 8"
MIPLIB_DIR="${MIPLIB_DIR:-$(pick_dir benches/data/miplib)}"
# G1 fixtures (all target-only): the
# Highs API around run, one g1-<x> command per directory, written to
# g1/<dir>/<name>.<scenario>.jsonl. Each set is "<instances>|<scenarios>";
# an instance is a netlib LP, a HiGHS check instance or a handcrafted model
# in $G1_SYNTH_DIR (.lp or .mps), looked up in that order (g1_path).
G1_SYNTH_DIR="$HERE/lps/g1"
G1_SETS="user-scale:user_scale:afiro sc50a adlittle blend|o3 b-2 o-4b5 o40
user-scale:user_scale:flugpl|o3b2
inf-cost:inf_cost:inf_cost_min inf_cost_max inf_cost_unbounded inf_cost_infeasible inf_cost_mip|default ic1e15
semi:semi:semi-continuous semi-integer semi_inf_upper semi_inconsistent semi_neg_lower|default relax
set-solution:set_solution:afiro adlittle sc50a|full full-dual sparse cold
set-solution:set_solution:flugpl egout small_mip|mip-full mip-partial mip-fractional
multi-objective:multi_objective:multiobj_lp multiobj_mip afiro|blend lex
presolve-api:presolve_api:afiro adlittle sc50a scagr7|presolve postsolve crossover crossover-dual
presolve-api:presolve_api:egout flugpl|presolve postsolve
presolve-api:presolve_api:semi-continuous|presolve
standard-form:standard_form:afiro adlittle blend sc50a std_form|default
rays:rays:ray_infeasible ray_unbounded ray_unbounded_dual|simplex presolve
rays:rays:flugpl egout small_mip|fixed
ill-conditioning:ill_conditioning:afiro sc50a blend ill_cond|m0c m0v m1c m1v
modify:modify:afiro adlittle sc50a|coeff scale add
mipfile:mipfile:flugpl egout small_mip|dense sparse"
# 1.16 only: normaliseCliqueRows is HiGHS #3286/#3318, after 1.15.1 (clqn.cpp).
if [ "$ORACLE_HIGHS" = 1.16.0 ]; then
  G1_SETS="$G1_SETS
clique-normalise:clique_normalise:clqn_rows p0548 lseu egout bell5 dcmulti|default"
fi
# Files copied into the committed fixture dir (keep the total < ~5 MB).
COMMIT_FILES="random.jsonl
unit/afiro.jsonl unit/sc50a.jsonl unit/kb2.jsonl
factor/afiro.jsonl
basiclu/sc50a.jsonl basiclu/dense70.jsonl basiclu/rand120.jsonl basiclu/neardep40.jsonl
simplex/afiro.dual.jsonl simplex/afiro.primal.jsonl
simplex/sc50a.dual.jsonl simplex/sc50a.primal.jsonl
simplex/kb2.dual.jsonl simplex/kb2.primal.jsonl
simplex/adlittle.dual.jsonl simplex/adlittle.primal.jsonl
simplex/blend.dual.jsonl simplex/blend.primal.jsonl
simplex/sc105.dual.jsonl simplex/sc105.primal.jsonl
simplex/share2b.dual.jsonl simplex/share2b.primal.jsonl
simplex/afirocut.dual.jsonl simplex/afirocut.primal.jsonl
simplex/afiromax.dual.jsonl simplex/afiromax.primal.jsonl
ipx/model/afiro.default.jsonl ipx/model/afiro.dual.jsonl
ipx/model/ipxmix.default.jsonl ipx/model/ipxmix.dual.jsonl ipx/model/ipxtall.default.jsonl
ipx/linalg/afiro.default.jsonl ipx/linalg/ipxmix.default.jsonl ipx/linalg/ipxmix.lu1.jsonl
ipx/linalg/ipxmix.dual.jsonl ipx/linalg/ipxtall.default.jsonl ipx/linalg/ipxdep.default.jsonl
ipx/linalg/ipxdep.dual.jsonl ipx/linalg/ipxdep.slack.jsonl
ipx/trace/afiro.default.jsonl ipx/trace/ipxmix.default.jsonl ipx/trace/ipxtall.default.jsonl
ipx/run/afiro.jsonl ipx/run/ipxmix.jsonl ipx/run/ipxtall.jsonl ipx/run/ipxdep.jsonl"

log() { printf '[gen_fixtures] %s\n' "$*" >&2; }

highs_commit() {
  if git -C "$HIGHS_REF" rev-parse HEAD >/dev/null 2>&1; then
    git -C "$HIGHS_REF" rev-parse HEAD
  else
    echo unknown
  fi
}

# shasum ships with macOS (perl); minimal Linux installs may only have coreutils.
if command -v shasum >/dev/null 2>&1; then
  sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
else
  sha() { sha256sum "$1" | cut -d' ' -f1; }
fi

# ------------------------------------------------- 1. patched copy of HiGHS
prepare_source() {
  [ -f "$HIGHS_REF/highs/Highs.h" ] || { log "HiGHS not found at $HIGHS_REF"; exit 1; }
  local want have="" want_ipx have_ipx="" want_w1 have_w1="" want_w3 have_w3="" want_w4 have_w4=""
  local want_w5 have_w5=""
  want="$(highs_commit) $(sha "$PATCH")"
  want_ipx="$(sha "$IPX_PATCH")"
  want_w1="$(sha "$W1_PATCH")"
  want_w3="$(sha "$W3_PATCH")"
  want_w4="$(sha "$W4_PATCH")"
  want_w5="$(sha "$W5_PATCH")"
  [ -f "$WORK/highs-src.stamp" ] && have="$(cat "$WORK/highs-src.stamp")"
  [ -f "$WORK/highs-src.ipx.stamp" ] && have_ipx="$(cat "$WORK/highs-src.ipx.stamp")"
  [ -f "$WORK/highs-src.w1.stamp" ] && have_w1="$(cat "$WORK/highs-src.w1.stamp")"
  [ -f "$WORK/highs-src.w3.stamp" ] && have_w3="$(cat "$WORK/highs-src.w3.stamp")"
  [ -f "$WORK/highs-src.w4.stamp" ] && have_w4="$(cat "$WORK/highs-src.w4.stamp")"
  [ -f "$WORK/highs-src.w5.stamp" ] && have_w5="$(cat "$WORK/highs-src.w5.stamp")"
  if [ "$want" = "$have" ] && [ "$want_ipx" = "$have_ipx" ] && [ "$want_w1" = "$have_w1" ] \
     && [ "$want_w3" = "$have_w3" ] && [ "$want_w4" = "$have_w4" ] && [ "$want_w5" = "$have_w5" ] \
     && [ -d "$SRC" ]; then
    log "patched HiGHS copy up to date ($want; ipx $want_ipx; w1 $want_w1; w3 $want_w3; w4 $want_w4; w5 $want_w5)"
    return
  fi
  log "copying $HIGHS_REF -> $SRC and applying $(basename "$PATCH") + $(basename "$IPX_PATCH") + $(basename "$W1_PATCH") + $(basename "$W3_PATCH") + $(basename "$W4_PATCH") + $(basename "$W5_PATCH")"
  mkdir -p "$WORK"
  rm -f "$WORK/highs-src.stamp" "$WORK/highs-src.ipx.stamp" "$WORK/highs-src.w1.stamp" \
    "$WORK/highs-src.w3.stamp" "$WORK/highs-src.w4.stamp" "$WORK/highs-src.w5.stamp"
  rsync -a --delete --exclude .git "$HIGHS_REF/" "$SRC/"
  (cd "$SRC" && patch -p1 --forward --quiet < "$PATCH")
  (cd "$SRC" && patch -p1 --forward --quiet < "$IPX_PATCH")
  (cd "$SRC" && patch -p1 --forward --quiet < "$W1_PATCH")
  (cd "$SRC" && patch -p1 --forward --quiet < "$W3_PATCH")
  (cd "$SRC" && patch -p1 --forward --quiet < "$W4_PATCH")
  (cd "$SRC" && patch -p1 --forward --quiet < "$W5_PATCH")
  echo "$want" > "$WORK/highs-src.stamp"
  echo "$want_ipx" > "$WORK/highs-src.ipx.stamp"
  echo "$want_w1" > "$WORK/highs-src.w1.stamp"
  echo "$want_w3" > "$WORK/highs-src.w3.stamp"
  echo "$want_w4" > "$WORK/highs-src.w4.stamp"
  echo "$want_w5" > "$WORK/highs-src.w5.stamp"
}

# ------------------------------------------------------------------ 2. build
build_oracle() {  # build_oracle <build dir> <fp-contract>
  local bdir="$1" fpc="$2"
  if [ ! -f "$bdir/CMakeCache.txt" ]; then
    cmake -S "$HERE" -B "$bdir" -DHIGHS_SRC="$SRC" -DCMAKE_BUILD_TYPE=Release \
      -DORACLE_FP_CONTRACT="$fpc" > "$bdir.configure.log" 2>&1 \
      || { cat "$bdir.configure.log"; exit 1; }
  fi
  cmake --build "$bdir" --target highs-oracle mip_trace api_seq -j "$JOBS" > "$bdir.build.log" 2>&1 \
    || { tail -50 "$bdir.build.log"; exit 1; }
}

# The HIPO=ON configuration: same patched source and flags, plus
# -DORACLE_HIPO=ON (CMakeLists.txt: HiPO linked against refblas). Checks the
# configure log for the pin, then runs the refblas unit checks.
build_oracle_hipo() {  # build_oracle_hipo <build dir>
  local bdir="$1"
  if [ ! -f "$bdir/CMakeCache.txt" ]; then
    cmake -S "$HERE" -B "$bdir" -DHIGHS_SRC="$SRC" -DCMAKE_BUILD_TYPE=Release \
      -DORACLE_FP_CONTRACT=off -DORACLE_HIPO=ON > "$bdir.configure.log" 2>&1 \
      || { cat "$bdir.configure.log"; exit 1; }
  fi
  grep -q "BLAS pinned to refblas" "$bdir.configure.log" \
    || { log "HiPO configure did not pin refblas (see $bdir.configure.log)"; exit 1; }
  if grep -Eq "Looking for blas|Found BLAS|Found OpenBLAS|Fetching OpenBLAS|Linking against|BLAS_ROOT is|OpenBLAS CMake config" "$bdir.configure.log"; then
    log "HiPO configure searched for another BLAS (see $bdir.configure.log)"; exit 1
  fi
  cmake --build "$bdir" --target highs-oracle test_refblas -j "$JOBS" > "$bdir.build.log" 2>&1 \
    || { tail -50 "$bdir.build.log"; exit 1; }
  "$bdir/test_refblas" > "$bdir.test_refblas.log" 2>&1 \
    || { cat "$bdir.test_refblas.log"; log "refblas unit checks FAILED"; exit 1; }
  log "refblas unit checks: $(tail -1 "$bdir.test_refblas.log")"
}

# ---------------------------------------------------------------- 3. generate
generate() {  # generate <oracle binary> <out dir>
  local bin="$1" out="$2" name mps s kind
  rm -rf "$out"
  mkdir -p "$out/unit" "$out/factor" "$out/basiclu" "$out/simplex"
  "$bin" random "$out/random.jsonl"
  for name in $UNIT_LPS; do
    "$bin" unit "$NETLIB_DIR/$name.mps" "$name" "$out/unit/$name.jsonl"
  done
  for name in $FACTOR_LPS; do
    "$bin" factor "$NETLIB_DIR/$name.mps" "$name" "$out/factor/$name.jsonl"
  done
  for name in $BASICLU_LPS; do
    "$bin" basiclu "$NETLIB_DIR/$name.mps" "$name" "$out/basiclu/$name.jsonl"
  done
  for name in $BASICLU_SYNTH; do
    "$bin" basiclu-synth "$name" "$out/basiclu/$name.jsonl"
  done
  for mps in "$NETLIB_DIR"/*.mps; do
    name="$(basename "$mps" .mps)"
    for s in $STRATEGIES; do
      kind=dual; [ "$s" = 4 ] && kind=primal
      "$bin" simplex "$mps" "$name" "$s" "$out/simplex/$name.$kind.jsonl"
    done
  done
  for name in $LP47_LPS; do
    for s in $STRATEGIES; do
      kind=dual; [ "$s" = 4 ] && kind=primal
      "$bin" simplex "$LP47_DIR/$name.mps" "$name" "$s" "$out/simplex/$name.$kind.jsonl"
    done
  done
  generate_ipx "$bin" "$out"
  generate_w1 "$bin" "$out"
  generate_w3 "$bin" "$out"
  generate_g1 "$bin" "$out"
  rm -rf "$VARIANT_DIR"
  python3 -I "$HERE/mkvariants.py" "$NETLIB_DIR" "$VARIANT_DIR"
  for base in $VARIANT_LPS; do
    for name in "${base}cut" "${base}max"; do
      for s in $STRATEGIES; do
        kind=dual; [ "$s" = 4 ] && kind=primal
        "$bin" simplex "$VARIANT_DIR/$name.mps" "$name" "$s" "$out/simplex/$name.$kind.jsonl"
      done
    done
  done
}

lp_path() {  # lp_path <name> -> netlib or synthetic MPS
  if [ -f "$NETLIB_DIR/$1.mps" ]; then echo "$NETLIB_DIR/$1.mps"; else echo "$SYNTH_DIR/$1.mps"; fi
}

generate_ipx() {  # generate_ipx <oracle binary> <out dir>
  local bin="$1" out="$2" mps name c
  mkdir -p "$out/ipx/model" "$out/ipx/linalg" "$out/ipx/trace" "$out/ipx/run"
  for mps in "$NETLIB_DIR"/*.mps "$SYNTH_DIR"/*.mps; do
    name="$(basename "$mps" .mps)"
    "$bin" ipx-model "$mps" "$name" default "$out/ipx/model/$name.default.jsonl"
    "$bin" ipx "$mps" "$name" "$out/ipx/run/$name.jsonl"
    for c in $IPX_TRACE_CFGS; do
      "$bin" ipx-trace "$mps" "$name" "${c#*:}" "$out/ipx/trace/$name.${c%%:*}.jsonl"
    done
  done
  for name in $IPX_MODEL_DUAL_LPS; do
    "$bin" ipx-model "$(lp_path "$name")" "$name" dualize=1 "$out/ipx/model/$name.dual.jsonl"
  done
  for name in $IPX_LINALG_LPS; do
    for c in $IPX_LINALG_CFGS; do
      "$bin" ipx-linalg "$(lp_path "$name")" "$name" "${c#*:}" "$out/ipx/linalg/$name.${c%%:*}.jsonl"
    done
  done
}

w1_path() {  # w1_path <name> -> HiGHS check instance or synthetic W1 MPS
  if [ -f "$CHECK_DIR/$1.mps" ]; then echo "$CHECK_DIR/$1.mps"; else echo "$W1_SYNTH_DIR/$1.mps"; fi
}

generate_w1() {  # generate_w1 <oracle binary> <out dir>
  local bin="$1" out="$2" name
  [ -d "$CHECK_DIR" ] || { log "HiGHS check instances not found at $CHECK_DIR"; exit 1; }
  rm -rf "$out/w1"
  mkdir -p "$out/w1/qp" "$out/w1/ranging" "$out/w1/iis" "$out/w1/icrash"
  for name in $W1_QP; do
    "$bin" w1-qp "$CHECK_DIR/$name.mps" "$name" "$out/w1/qp/$name.jsonl"
  done
  for name in $W1_RANGING; do
    "$bin" w1-ranging "$NETLIB_DIR/$name.mps" "$name" "$out/w1/ranging/$name.jsonl"
  done
  for name in $W1_IIS; do
    "$bin" w1-iis "$(w1_path "$name")" "$name" "$out/w1/iis/$name.jsonl"
  done
  for name in $W1_ICRASH; do
    "$bin" w1-icrash "$NETLIB_DIR/$name.mps" "$name" "$out/w1/icrash/$name.jsonl"
  done
}

generate_w3() {  # generate_w3 <oracle binary> <out dir>
  local bin="$1" out="$2" name
  rm -rf "$out/w3"
  mkdir -p "$out/w3/cupdlp" "$out/w3/hipdlp"
  for name in $W3_LPS; do
    "$bin" w3-cupdlp "$NETLIB_DIR/$name.mps" "$name" "$out/w3/cupdlp/$name.jsonl"
    "$bin" w3-hipdlp "$NETLIB_DIR/$name.mps" "$name" "$out/w3/hipdlp/$name.jsonl"
  done
}

generate_w4() {  # generate_w4 <HIPO=ON oracle binary> <out dir>
  local bin="$1" out="$2" name o
  rm -rf "$out/w4"
  mkdir -p "$out/w4/orderings" "$out/w4/hipo"
  for name in $W4_LPS; do
    for o in $W4_ORDERINGS; do
      "$bin" w4-orderings "$NETLIB_DIR/$name.mps" "$name" "$o" "$out/w4/orderings/$name.$o.jsonl"
    done
    "$bin" w4-hipo "$NETLIB_DIR/$name.mps" "$name" "$out/w4/hipo/$name.jsonl"
  done
}

# SIP / PAMI and HiPO at
# threads 1/2/4/8 with parallel=on, and the MIP trace at threads 2/4/8. No
# HIGHS_ORACLE_AC_DUMP at threads > 1 (its writer is not thread-safe and a
# stolen-then-cancelled analytic centre adds a record).
generate_w5() {  # generate_w5 <oracle binary> <HIPO=ON oracle binary> <mip_trace binary> <out dir>
  local bin="$1" hipo_bin="$2" mip_bin="$3" out="$4" name s n ty sys systems
  rm -rf "$out/w5"
  mkdir -p "$out/w5/simplex" "$out/w5/hipo" "$out/w5/mip"
  for name in $W5_LPS; do
    for s in $W5_STRATEGIES; do
      for n in $W5_THREADS; do
        "$bin" w5-simplex "$NETLIB_DIR/$name.mps" "$name" "$s" "$n" \
          "$out/w5/simplex/${name}_s${s}_t${n}.jsonl"
      done
    done
    for ty in $W5_HIPO_TYPES; do
      for n in $W5_THREADS; do
        systems="as ne"; [ "$n" = 1 ] && systems="choose as ne"
        for sys in $systems; do
          "$hipo_bin" w5-hipo "$NETLIB_DIR/$name.mps" "$name" "$ty" "$(w5_system "$sys")" "$n" \
            "$out/w5/hipo/${name}_${ty}_${sys}_t${n}.jsonl"
        done
      done
    done
  done
  for name in $W5_MIP; do
    [ -f "$MIPLIB_DIR/$name.mps" ] || { log "w5: MIPLIB instance $MIPLIB_DIR/$name.mps not found"; exit 1; }
    for n in $W5_MIP_THREADS; do
      env -u HIGHS_ORACLE_AC_DUMP "$mip_bin" "$MIPLIB_DIR/$name.mps" "threads=$n" parallel=on \
        > "$out/w5/mip/$name.t$n.trace"
      grep -q '^final ' "$out/w5/mip/$name.t$n.trace" \
        || { log "w5: $name.t$n.trace has no final line"; exit 1; }
    done
  done
}

w5_system() {  # w5_system <choose|as|ne> -> hipo_system value
  case "$1" in
    choose) echo choose ;;
    as) echo augmented ;;
    ne) echo normaleq ;;
    *) log "w5: unknown HiPO system $1"; exit 2 ;;
  esac
}

g1_path() {  # g1_path <name> -> netlib, check instance or handcrafted model
  local p
  for p in "$NETLIB_DIR/$1.mps" "$CHECK_DIR/$1.mps" "$G1_SYNTH_DIR/$1.lp" "$G1_SYNTH_DIR/$1.mps"; do
    if [ -f "$p" ]; then echo "$p"; return; fi
  done
  log "g1: no model $1"; exit 2
}

generate_g1() {  # generate_g1 <oracle binary> <out dir>
  local bin="$1" out="$2" line cmd dir names scns name s
  rm -rf "$out/g1"
  while IFS= read -r line; do
    cmd="${line%%:*}"; line="${line#*:}"
    dir="${line%%:*}"; line="${line#*:}"
    names="${line%%|*}"; scns="${line#*|}"
    mkdir -p "$out/g1/$dir"
    for name in $names; do
      for s in $scns; do
        "$bin" "g1-$cmd" "$(g1_path "$name")" "$name" "$s" "$out/g1/$dir/$name.$s.jsonl"
      done
    done
  done <<< "$G1_SETS"
}

# MIP traces (FORMAT.md "MIP IPM records"; the committed mip/ fixtures are the
# 1.15.1 macOS output): mip_trace at its defaults (threads=1, parallel=off)
# with HIGHS_ORACLE_AC_DUMP, on W5_MIP, plus sp150x300d with presolve off.
# Written to <out>/mip/; never copied into the committed dir.
generate_mip() {  # generate_mip <mip_trace binary> <out dir>
  local mip_bin="$1" out="$2" name tag
  mkdir -p "$out/mip"
  for tag in $W5_MIP sp150x300d.presolve_off; do
    name="${tag%.presolve_off}"
    [ -f "$MIPLIB_DIR/$name.mps" ] || { log "mip: MIPLIB instance $MIPLIB_DIR/$name.mps not found"; exit 1; }
    if [ "$tag" = "$name" ]; then
      HIGHS_ORACLE_AC_DUMP="$out/mip/$tag.ac" "$mip_bin" "$MIPLIB_DIR/$name.mps" > "$out/mip/$tag.trace"
    else
      HIGHS_ORACLE_AC_DUMP="$out/mip/$tag.ac" "$mip_bin" "$MIPLIB_DIR/$name.mps" presolve=off \
        > "$out/mip/$tag.trace"
    fi
    grep -q '^final ' "$out/mip/$tag.trace" || { log "mip: $tag.trace has no final line"; exit 1; }
  done
}

# -------------------------------------------------- 4. publish + manifest
# The 1.15.1 oracle copies COMMIT_FILES into the committed dir and writes the
# manifest there (then copies it to $LARGE). Any other ORACLE_HIGHS writes
# its manifest to $LARGE only: the location column still names the files the
# 1.15.1 run commits, so a test run pointed at $LARGE (RUSOLVE_FIXTURES_DIR)
# sees the same layout.
publish() {
  local f committed mdir="$SMALL"
  if [ -n "$VSUFFIX" ]; then
    mdir="$LARGE"
  else
    mkdir -p "$SMALL"
    find "$SMALL" -name '*.jsonl' -delete
    for f in $COMMIT_FILES; do
      mkdir -p "$SMALL/$(dirname "$f")"
      cp "$LARGE/$f" "$SMALL/$f"
    done
  fi
  {
    printf '# rusolve-ekk golden fixtures, generated by tools/highs-oracle/gen_fixtures.sh\n'
    printf '# highs_commit %s\n' "$(highs_commit)"
    printf '# patch_sha256 %s\n' "$(sha "$PATCH")"
    printf '# ipx_patch_sha256 %s\n' "$(sha "$IPX_PATCH")"
    printf '# w1_patch_sha256 %s\n' "$(sha "$W1_PATCH")"
    printf '# w3_patch_sha256 %s\n' "$(sha "$W3_PATCH")"
    printf '# w4_patch_sha256 %s\n' "$(sha "$W4_PATCH")"
    printf '# w5_patch_sha256 %s\n' "$(sha "$W5_PATCH")"
    printf '# location: committed = crates/rusolve-ekk/tests/fixtures/<file>;'
    printf ' target = target/highs-fixtures/<file> (regenerate with gen_fixtures.sh)\n'
    printf '# file\tbytes\tsha256\tlocation\n'
    (cd "$LARGE" && find . -name '*.jsonl' | sed 's|^\./||' | LC_ALL=C sort) | while read -r f; do
      committed=target
      if printf '%s\n' $COMMIT_FILES | grep -qx "$f"; then committed=committed; fi
      printf '%s\t%s\t%s\t%s\n' "$f" "$(wc -c < "$LARGE/$f" | tr -d ' ')" "$(sha "$LARGE/$f")" "$committed"
    done
    # The w5/mip traces are plain text, not JSONL records, so they are listed
    # as comments (tests/common/fixture.rs loads every row as JSONL):
    # "# trace <file> <bytes> <sha256>", target-only.
    (cd "$LARGE" && find . -path './w5/*' -name '*.trace' | sed 's|^\./||' | LC_ALL=C sort) | while read -r f; do
      printf '# trace %s %s %s\n' "$f" "$(wc -c < "$LARGE/$f" | tr -d ' ')" "$(sha "$LARGE/$f")"
    done
  } > "$mdir/MANIFEST.tsv.new"
  # Keep hand-maintained per-platform sha alternatives ("<mac>|<linux>", see
  # tests/common/fixture.rs) when this platform's sha is one of them.
  awk -F'\t' -v OFS='\t' 'NR == FNR { if ($3 ~ /\|/) alt[$1] = $3; next }
    ($1 in alt) && index("|" alt[$1] "|", "|" $3 "|") { $3 = alt[$1] } { print }' \
    "$mdir/MANIFEST.tsv" "$mdir/MANIFEST.tsv.new" > "$mdir/MANIFEST.tsv.tmp" 2>/dev/null \
    || cp "$mdir/MANIFEST.tsv.new" "$mdir/MANIFEST.tsv.tmp"
  mv "$mdir/MANIFEST.tsv.tmp" "$mdir/MANIFEST.tsv"
  rm -f "$mdir/MANIFEST.tsv.new"
  if [ -z "$VSUFFIX" ]; then
    cp "$SMALL/MANIFEST.tsv" "$LARGE/MANIFEST.tsv"
    log "committed fixtures: $(du -sh "$SMALL" | cut -f1) in $SMALL"
  fi
  log "all fixtures:       $(du -sh "$LARGE" | cut -f1) in $LARGE"
}

# --w5 regenerates only w5/, so the other rows of the committed manifest must
# already be in $LARGE with their recorded sha256 (else publish would drop or
# rewrite them).
check_non_w5() {
  local f bytes want got bad=0 manifest="$SMALL/MANIFEST.tsv"
  [ -n "$VSUFFIX" ] && manifest="$LARGE/MANIFEST.tsv"
  [ -f "$manifest" ] || { log "--w5: $manifest missing; run without --w5 first"; exit 1; }
  while IFS=$'\t' read -r f bytes want _; do
    case "$f" in '#'*|w5/*|'') continue ;; esac
    if [ ! -f "$LARGE/$f" ]; then log "--w5: $LARGE/$f missing"; bad=1; continue; fi
    got="$(sha "$LARGE/$f")"
    case "|$want|" in *"|$got|"*) ;; *) log "--w5: $LARGE/$f sha256 differs from MANIFEST.tsv"; bad=1 ;; esac
  done < "$manifest"
  [ $bad = 0 ] || { log "--w5: run without --w5 to regenerate everything"; exit 1; }
}

main() {
  local determinism=0 fpcheck=0 w5only=0 miponly=0 arg
  for arg in "$@"; do
    case "$arg" in
      --determinism) determinism=1 ;;
      --fp-contract-check) fpcheck=1 ;;
      --w5) w5only=1 ;;
      --mip) miponly=1 ;;
      *) log "unknown argument $arg"; exit 64 ;;
    esac
  done
  prepare_source
  build_oracle "$WORK/build" off
  build_oracle_hipo "$WORK/build-hipo"
  local bin="$WORK/build/highs-oracle" hipo_bin="$WORK/build-hipo/highs-oracle"
  local mip_bin="$WORK/build/mip_trace"
  if [ $miponly = 1 ]; then
    log "generating mip fixtures"
    generate_mip "$mip_bin" "$LARGE"
    if [ $determinism = 1 ]; then
      log "determinism (mip): regenerating once more"
      rm -rf "$WORK/det-mip"
      generate_mip "$mip_bin" "$WORK/det-mip"
      if diff -r -q "$LARGE/mip" "$WORK/det-mip/mip"; then
        log "determinism (mip): OK (2 runs byte-identical)"
      else
        log "determinism (mip): FAILED"; exit 1
      fi
    fi
    return
  fi
  if [ $w5only = 1 ]; then
    check_non_w5
    log "generating w5 fixtures"
    generate_w5 "$bin" "$hipo_bin" "$mip_bin" "$LARGE"
    publish
    if [ $determinism = 1 ]; then
      log "determinism (w5): regenerating twice more"
      generate_w5 "$bin" "$hipo_bin" "$mip_bin" "$WORK/det-a"
      generate_w5 "$bin" "$hipo_bin" "$mip_bin" "$WORK/det-b"
      if diff -r -q "$WORK/det-a/w5" "$WORK/det-b/w5" && diff -r -q "$LARGE/w5" "$WORK/det-a/w5"; then
        log "determinism (w5): OK (3 runs byte-identical)"
      else
        log "determinism (w5): FAILED"; exit 1
      fi
    fi
    return
  fi
  log "generating fixtures"
  generate "$bin" "$LARGE"
  generate_w4 "$hipo_bin" "$LARGE"
  generate_w5 "$bin" "$hipo_bin" "$mip_bin" "$LARGE"
  publish
  if [ $determinism = 1 ]; then
    log "determinism: regenerating twice more"
    generate "$bin" "$WORK/det-a"
    generate_w4 "$hipo_bin" "$WORK/det-a"
    generate_w5 "$bin" "$hipo_bin" "$mip_bin" "$WORK/det-a"
    generate "$bin" "$WORK/det-b"
    generate_w4 "$hipo_bin" "$WORK/det-b"
    generate_w5 "$bin" "$hipo_bin" "$mip_bin" "$WORK/det-b"
    if diff -r -q "$WORK/det-a" "$WORK/det-b" && diff -r -q -x MANIFEST.tsv "$LARGE" "$WORK/det-a"; then
      log "determinism: OK (3 runs byte-identical)"
    else
      log "determinism: FAILED"; exit 1
    fi
  fi
  if [ $fpcheck = 1 ]; then
    log "fp-contract check: building with -ffp-contract=on"
    build_oracle "$WORK/build-fpc-on" on
    generate "$WORK/build-fpc-on/highs-oracle" "$WORK/fpc-on"
    # w4/ comes from the separate HIPO=ON configuration (refblas is pinned
    # to -ffp-contract=off whatever this check asks), so it is not compared;
    # neither is w5/ (generate_w5 is not run on this build).
    (cd "$LARGE" && find . -name '*.jsonl' -not -path './w4/*' -not -path './w5/*' | LC_ALL=C sort) | while read -r f; do
      if cmp -s "$LARGE/$f" "$WORK/fpc-on/$f"; then echo "same   $f"; else echo "DIFFER $f"; fi
    done > "$WORK/fpc-on.report"
    log "fp-contract check: $(grep -c DIFFER "$WORK/fpc-on.report" || true) of" \
        "$(wc -l < "$WORK/fpc-on.report" | tr -d ' ') files differ (see $WORK/fpc-on.report)"
  fi
}

main "$@"
