#!/bin/sh
# rusolve installer: the rusolve-cli binary from a GitHub release, and optionally the
# OR-Tools shared libraries that the CP-SAT engine (mip_engine=cpsat) loads at run time.
#
#   curl -fsSL https://raw.githubusercontent.com/CREAIT-nl/rusolve/main/install.sh | sh
#   sh install.sh --help
#
# Design and the choice of where OR-Tools comes from: docs/design/install.md.
# Needs: curl or wget, tar, sha256sum or shasum; for --with-cpsat also unzip or python3.
#
# By default it also registers the MCP server with the agent CLIs it finds (claude, codex,
# gemini) and installs the rusolve skills for Claude Code; --no-agent, --no-mcp and --no-skills
# opt out. A failure there warns and never fails the install. When the bin directory is not on
# PATH it also adds it to the shell profiles (one sourced line each, recorded in the manifest);
# --no-modify-path opts out.
set -eu

# ---- pins: the one place to change when OR-Tools moves --------------------------------------
# OR-Tools (Apache-2.0, https://github.com/google/or-tools) is downloaded from PyPI on request
# and is never redistributed by rusolve. Only the shared libraries inside the wheel are used
# (they do not link libpython), so any CPython tag works: cp312 is used for every platform.
ORTOOLS_VERSION=9.15.6755
ORTOOLS_SIZE_MB=30

ortools_pin() { # platform -> sets ORT_URL and ORT_SHA ("" when none is pinned)
  case $1 in
    linux-x86_64)
      ORT_URL=https://files.pythonhosted.org/packages/49/0f/6d6d722102a0ceccf4a5038e2bc91d023da84a6dba98482a4634df3d27ab/ortools-9.15.6755-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl
      ORT_SHA=033836c0eb33bc72697a299e0caedbb25fc9d1cee0b13832d69cb30405f57b3e ;;
    linux-aarch64)
      ORT_URL=https://files.pythonhosted.org/packages/86/a9/37cb31fc5ffbec2650ebb0d2538a83842b5693a788a0ec6057559dab1169/ortools-9.15.6755-cp312-cp312-manylinux_2_26_aarch64.manylinux_2_28_aarch64.whl
      ORT_SHA=3cd6bec0a2e00e3891a53e3b436f45a1000269f302085572f49e9856b7f8eaf0 ;;
    macos-arm64)
      ORT_URL=https://files.pythonhosted.org/packages/f1/b6/7e6618ef7a88e8eb706a8a876806b4d336f1bef8c574f8a02d2da3e483ef/ortools-9.15.6755-cp312-cp312-macosx_11_0_arm64.whl
      ORT_SHA=1e16686c2b457fa6242c474ab890ee1712347ab53678e0d2fab307ae03e97a4b ;;
    macos-x86_64)
      ORT_URL=https://files.pythonhosted.org/packages/6a/fc/9fa53f1a13710e6183df4d00fe4988c79a55b501e282645d49f1e250437f/ortools-9.15.6755-cp312-cp312-macosx_10_15_x86_64.whl
      ORT_SHA=ae1c6e1fd844b4d756b22eb6c0ed574ea4342ee206d807c4f903039e748228fa ;;
    *) ORT_URL=; ORT_SHA= ;;
  esac
}

# The plugin marketplace (.claude-plugin/marketplace.json in the repository) that carries the
# skills: a GitHub owner/repo, a git URL or a local path. Overridable for forks and tests.
MARKETPLACE_SRC=${RUSOLVE_INSTALL_MARKETPLACE:-CREAIT-nl/rusolve}
MCP_NAME=rusolve
MKT_NAME=rusolve
PLUGIN_ID=rusolve-skills@rusolve

DEFAULT_BASE_URL=https://github.com/CREAIT-nl/rusolve/releases
BASE_URL=${RUSOLVE_INSTALL_BASE_URL:-$DEFAULT_BASE_URL}
BASE_URL=${BASE_URL%/}
# Layout under BASE_URL (the shape of GitHub release downloads):
#   --version vX.Y.Z (a tag):  <base>/download/<tag>/<file>
#   the default (latest):      <base>/latest/download/<file>

VERSION=latest
PREFIX=
CPSAT=ask
UNINSTALL=0
DRY=0
YES=0
AG_MCP=1
AG_SKILLS=1
MODPATH=1
AGOUT=
TMP=

usage() {
  cat <<'EOF'
Install rusolve-cli from a GitHub release.

Usage: install.sh [options]
  --version vX.Y.Z   release to install (default: latest)
  --prefix DIR       install under DIR (default: $HOME/.local); the binary goes to DIR/bin
  --with-cpsat       also install the OR-Tools libraries for mip_engine=cpsat
  --no-cpsat         do not install them (no question asked)
  --uninstall        remove what a previous run installed (reads DIR/lib/rusolve/install-manifest)
  --dry-run          print what would happen, change nothing, download nothing
  --yes              never ask; CP-SAT stays off unless --with-cpsat is given
  --no-mcp           do not register the MCP server with claude, codex or gemini
  --no-skills        do not install the rusolve skills plugin for Claude Code
  --no-agent         both of the above
  --no-modify-path   do not add DIR/bin to PATH in your shell profiles (done by default when it
                     is not on PATH yet: one line sourcing DIR/lib/rusolve/env in ~/.profile (created
                     if missing) and in ~/.bashrc, ~/.bash_profile and ~/.zshenv where they
                     exist; ~/.zshenv is also created when zsh is installed)
  -h, --help         this text

Environment: RUSOLVE_INSTALL_BASE_URL replaces https://github.com/CREAIT-nl/rusolve/releases
(a mirror or a fork; same layout as GitHub). With --version vX.Y.Z the files are fetched from
<base>/download/<tag>/<file>; without it, from <base>/latest/download/<file>.
RUSOLVE_INSTALL_MARKETPLACE replaces CREAIT-nl/rusolve as the source of the skills plugin
marketplace (an owner/repo, a git URL or a local path).
By default the MCP server is registered (user scope, absolute path) with each of claude, codex
and gemini found on PATH, and the skills are installed for Claude Code. A registration that
exists with another command is left alone and reported. --uninstall undoes only what was added,
profile lines included.
With `curl ... | sh` stdin is the script, so nothing is asked: pass --with-cpsat to opt in.
EOF
}

say() { printf '%s\n' "$*"; }
WRITING=0
err() {
  printf 'install.sh: error: %s\n' "$*" >&2
  # a failure after files were written still leaves a manifest, so --uninstall can clean up
  if [ "$WRITING" = 1 ]; then WRITING=0; write_manifest || true; fi
  exit 1
}
warn() { printf 'install.sh: warning: %s\n' "$*" >&2; }

cleanup() { if [ -n "$TMP" ]; then rm -rf "$TMP"; fi; }
trap cleanup EXIT
trap 'exit 130' INT TERM

while [ $# -gt 0 ]; do
  case $1 in
    --version) [ $# -ge 2 ] || err "--version needs a value"; VERSION=$2; shift 2 ;;
    --version=*) VERSION=${1#--version=}; shift ;;
    --prefix) [ $# -ge 2 ] || err "--prefix needs a value"; PREFIX=$2; shift 2 ;;
    --prefix=*) PREFIX=${1#--prefix=}; shift ;;
    --with-cpsat) CPSAT=yes; shift ;;
    --no-cpsat) CPSAT=no; shift ;;
    --uninstall) UNINSTALL=1; shift ;;
    --dry-run) DRY=1; shift ;;
    --yes|-y) YES=1; shift ;;
    --no-mcp) AG_MCP=0; shift ;;
    --no-skills) AG_SKILLS=0; shift ;;
    --no-agent) AG_MCP=0; AG_SKILLS=0; shift ;;
    --no-modify-path) MODPATH=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; err "unknown option: $1" ;;
  esac
done

case $VERSION in
  latest) ;;
  *[!A-Za-z0-9._-]*|'') err "invalid --version: $VERSION" ;;
  v*) ;;
  *) VERSION=v$VERSION ;;
esac

[ -n "${HOME:-}" ] || [ -n "$PREFIX" ] || err "HOME is not set; pass --prefix"
[ -n "$PREFIX" ] || PREFIX=$HOME/.local
case $PREFIX in
  /*) ;;
  *) PREFIX=$PWD/$PREFIX ;;
esac
PREFIX=${PREFIX%/}
BIN_DIR=$PREFIX/bin
LIB_DIR=$PREFIX/lib/rusolve
ORT_DIR=$LIB_DIR/ortools
MANIFEST=$LIB_DIR/install-manifest
BIN=$BIN_DIR/rusolve-cli
ENV_FILE=$LIB_DIR/env
# the one line added to a shell profile; the marker makes it recognisable and removable
PROFILE_LINE=". '$ENV_FILE' # added by rusolve install.sh"


# ---- agent integration ----------------------------------------------------------------------
# Each CLI's syntax was read from its --help: claude mcp add|get|remove, codex mcp add|get|remove,
# gemini mcp add|list|remove, claude plugin marketplace add|list|remove, claude plugin
# install|list|uninstall. Nothing here may fail the install: every function returns 0.
have() { command -v "$1" >/dev/null 2>&1; }
ag() { AGOUT="$AGOUT  $*
"; }
shq() { printf "'%s'" "$1"; }

# run a command, or only show it on --dry-run; its output is left in AG_OUT
ag_run() {
  if [ "$DRY" = 1 ]; then
    ag "would run: $*"
    AG_OUT=
    return 0
  fi
  AG_OUT=$("$@" 2>&1 </dev/null) || return 1
  return 0
}

mcp_state() { # cli -> MCP_STATE: none | same (registered with $BIN) | other (another command)
  MCP_STATE=none
  case $1 in
    claude)
      out=$(claude mcp get "$MCP_NAME" 2>&1 </dev/null) || return 0
      MCP_STATE=other
      if printf '%s\n' "$out" | awk -v w="Command: $BIN" '{ sub(/^[ \t]+/, ""); if ($0 == w) f = 1 } END { exit !f }'; then
        MCP_STATE=same
      fi ;;
    codex)
      out=$(codex mcp get "$MCP_NAME" 2>&1 </dev/null) || return 0
      MCP_STATE=other
      if printf '%s\n' "$out" | awk -v w="command: $BIN" '{ sub(/^[ \t]+/, ""); if ($0 == w) f = 1 } END { exit !f }'; then
        MCP_STATE=same
      fi ;;
    gemini)
      out=$(gemini mcp list 2>&1 </dev/null) || return 0
      if printf '%s\n' "$out" | grep -Fq " $MCP_NAME: "; then
        MCP_STATE=other
        if printf '%s\n' "$out" | grep -Fq " $MCP_NAME: $BIN mcp (stdio)"; then MCP_STATE=same; fi
      fi ;;
  esac
  return 0
}

mcp_add_cmd() { # cli -> the registration command line, quoted for a person
  case $1 in
    claude) printf '%s' "claude mcp add --scope user $MCP_NAME -- $(shq "$BIN") mcp" ;;
    codex) printf '%s' "codex mcp add $MCP_NAME -- $(shq "$BIN") mcp" ;;
    gemini) printf '%s' "gemini mcp add --scope user $MCP_NAME $(shq "$BIN") mcp" ;;
  esac
}

mcp_remove_cmd() { # cli
  case $1 in
    claude) printf '%s' "claude mcp remove --scope user $MCP_NAME" ;;
    codex) printf '%s' "codex mcp remove $MCP_NAME" ;;
    gemini) printf '%s' "gemini mcp remove --scope user $MCP_NAME" ;;
  esac
}

mcp_remove_run() { # cli -> AG_OUT, status of the removal
  case $1 in
    claude) AG_OUT=$(claude mcp remove --scope user "$MCP_NAME" 2>&1 </dev/null) ;;
    codex) AG_OUT=$(codex mcp remove "$MCP_NAME" 2>&1 </dev/null) ;;
    gemini) AG_OUT=$(gemini mcp remove --scope user "$MCP_NAME" 2>&1 </dev/null) ;;
  esac
}

mcp_register() { # cli
  c=$1
  mcp_state "$c"
  case $MCP_STATE in
    same) ag "$c: MCP server '$MCP_NAME' already registered with $BIN"; return 0 ;;
    other)
      ag "$c: an MCP server '$MCP_NAME' is already registered with another command; left alone."
      ag "    To point it at this install: $(mcp_remove_cmd "$c") && $(mcp_add_cmd "$c")"
      return 0 ;;
  esac
  case $c in
    claude) set -- claude mcp add --scope user "$MCP_NAME" -- "$BIN" mcp ;;
    codex) set -- codex mcp add "$MCP_NAME" -- "$BIN" mcp ;;
    gemini) set -- gemini mcp add --scope user "$MCP_NAME" "$BIN" mcp ;;
  esac
  if ag_run "$@"; then
    if [ "$DRY" = 0 ]; then
      ag "$c: MCP server '$MCP_NAME' registered (user scope) -> $BIN"
      printf 'm %s %s\n' "$c" "$MCP_NAME" >>"$NEWMAN"
    fi
  else
    warn "$c: could not register the MCP server: $(printf '%s' "$AG_OUT" | tail -n 1)"
    ag "$c: MCP registration failed; run it yourself: $(mcp_add_cmd "$c")"
  fi
  return 0
}

skills_install() { # Claude Code only
  mkt_cmd="claude plugin marketplace add $(shq "$MARKETPLACE_SRC") --scope user"
  plug_cmd="claude plugin install $PLUGIN_ID --scope user"
  have_mkt=0; have_plug=0
  if [ "$DRY" = 0 ]; then
    out=$(claude plugin marketplace list 2>&1 </dev/null) || out=
    if printf '%s\n' "$out" | grep -Eq "[[:space:]]$MKT_NAME\$"; then have_mkt=1; fi
    out=$(claude plugin list 2>&1 </dev/null) || out=
    if printf '%s\n' "$out" | grep -Fq "$PLUGIN_ID"; then have_plug=1; fi
  fi
  if [ "$have_plug" = 1 ]; then
    ag "claude: skills plugin $PLUGIN_ID already installed"
    return 0
  fi
  if [ "$have_mkt" = 0 ]; then
    if ag_run claude plugin marketplace add "$MARKETPLACE_SRC" --scope user; then
      if [ "$DRY" = 0 ]; then printf 's claude %s\n' "$MKT_NAME" >>"$NEWMAN"; fi
    else
      warn "claude: could not add the plugin marketplace $MARKETPLACE_SRC: $(printf '%s' "$AG_OUT" | tail -n 1)"
      ag "claude: skills not installed; run these yourself: $mkt_cmd && $plug_cmd"
      return 0
    fi
  fi
  if ag_run claude plugin install "$PLUGIN_ID" --scope user; then
    if [ "$DRY" = 0 ]; then
      ag "claude: skills plugin $PLUGIN_ID installed (user scope)"
      printf 'p claude %s\n' "$PLUGIN_ID" >>"$NEWMAN"
    fi
  else
    warn "claude: could not install $PLUGIN_ID: $(printf '%s' "$AG_OUT" | tail -n 1)"
    ag "claude: skills not installed; run it yourself: $plug_cmd"
  fi
  return 0
}

agent_setup() { # fills AGOUT
  if [ "$AG_MCP" = 0 ]; then ag "MCP registration skipped (--no-mcp or --no-agent)"; fi
  if [ "$AG_SKILLS" = 0 ]; then ag "skills skipped (--no-skills or --no-agent)"; fi
  found=0
  for c in claude codex gemini; do
    if ! have "$c"; then continue; fi
    found=1
    if [ "$AG_MCP" = 1 ]; then mcp_register "$c" || true; fi
    if [ "$c" = claude ] && [ "$AG_SKILLS" = 1 ]; then skills_install || true; fi
  done
  if [ "$found" = 0 ] && { [ "$AG_MCP" = 1 ] || [ "$AG_SKILLS" = 1 ]; }; then
    ag "none of claude, codex, gemini is on PATH; nothing registered"
  fi
  return 0
}

agent_manual() { # what to do for every other MCP client
  if [ "$AG_MCP" = 0 ]; then return 0; fi
  say "Other MCP clients (Claude Desktop, Cursor, Windsurf, ...): add this server to their MCP config:"
  say "  {\"mcpServers\": {\"$MCP_NAME\": {\"command\": \"$BIN\", \"args\": [\"mcp\"]}}}"
}

agent_undo() { # manifest text: remove the registrations it lists, only while they are still ours
  # servers first, then the plugin, then its marketplace (removing the marketplace drops the plugin)
  { printf '%s\n' "$1" | grep '^m ' || true; printf '%s\n' "$1" | grep '^p ' || true; printf '%s\n' "$1" | grep '^s ' || true; } | while IFS= read -r line; do
    kind=${line%% *}
    rest=${line#* }
    case $kind in
      m)
        [ "$AG_MCP" = 1 ] || continue
        c=${rest%% *}
        if ! have "$c"; then
          say "  $c is not on PATH; remove the MCP server yourself: $(mcp_remove_cmd "$c")"
          continue
        fi
        mcp_state "$c"
        case $MCP_STATE in
          none) say "  $c: MCP server '$MCP_NAME' is already gone" ;;
          other) say "  $c: MCP server '$MCP_NAME' now has another command; left alone" ;;
          same)
            if [ "$DRY" = 1 ]; then
              say "  would run: $(mcp_remove_cmd "$c")"
            elif mcp_remove_run "$c"; then
              say "  $c: removed MCP server '$MCP_NAME'"
            else
              warn "$c: could not remove the MCP server: $AG_OUT"
              say "    run it yourself: $(mcp_remove_cmd "$c")"
            fi ;;
        esac ;;
      p)
        [ "$AG_SKILLS" = 1 ] || continue
        id=${rest#* }
        if ! have claude; then say "  claude is not on PATH; run: claude plugin uninstall $id --scope user"; continue; fi
        if [ "$DRY" = 1 ]; then say "  would run: claude plugin uninstall $id --scope user"
        elif o=$(claude plugin uninstall "$id" --scope user 2>&1 </dev/null); then say "  claude: removed plugin $id"
        else warn "claude: could not uninstall $id: $o"; fi ;;
      s)
        [ "$AG_SKILLS" = 1 ] || continue
        n=${rest#* }
        if ! have claude; then say "  claude is not on PATH; run: claude plugin marketplace remove $n"; continue; fi
        if [ "$DRY" = 1 ]; then say "  would run: claude plugin marketplace remove $n"
        elif o=$(claude plugin marketplace remove "$n" 2>&1 </dev/null); then say "  claude: removed marketplace $n"
        else warn "claude: could not remove marketplace $n: $o"; fi ;;
    esac
  done
}

# ---- PATH in the shell profiles -------------------------------------------------------------
# Like rustup: an env file (PREFIX/lib/rusolve/env) puts the bin directory on PATH once, and one
# line in each profile sources it. The profiles are ~/.profile (created if missing, read by
# login shells), and ~/.bash_profile, ~/.bashrc and ${ZDOTDIR:-~}/.zshenv where they exist
# (~/.zshenv is created when zsh is installed). Manifest lines: `e FILE` (the line was added to
# an existing file), `E FILE` (the file was created for it). Every function returns 0.
env_file_text() {
  printf '%s\n' "# rusolve: put the install directory on PATH (sourced from your shell profile)"
  printf '%s\n' "case \":\${PATH}:\" in"
  printf '%s\n' "  *:'$BIN_DIR':*) ;;"
  printf '%s\n' "  *) export PATH='$BIN_DIR':\"\$PATH\" ;;"
  printf '%s\n' "esac"
}

path_add_to() { # file create(0|1): append PROFILE_LINE once
  f=$1
  if [ -f "$f" ]; then
    if grep -Fxq -- "$PROFILE_LINE" "$f" 2>/dev/null; then
      ag "PATH: $f already sources $ENV_FILE"
      PT_FOUND=1
      return 0
    fi
    kind=e
  elif [ "$2" = 1 ]; then
    kind=E
  else
    return 0
  fi
  if [ "$DRY" = 1 ]; then
    if [ "$kind" = E ]; then ag "would create $f with the line: $PROFILE_LINE"; else ag "would append to $f the line: $PROFILE_LINE"; fi
    return 0
  fi
  if [ -s "$f" ] && [ -n "$(tail -c 1 "$f" 2>/dev/null)" ]; then printf '\n' >>"$f" 2>/dev/null || true; fi
  if printf '%s\n' "$PROFILE_LINE" >>"$f" 2>/dev/null; then
    printf '%s %s\n' "$kind" "$f" >>"$NEWMAN"
    PT_ADDED="$PT_ADDED $f"
  else
    warn "could not write $f"
  fi
  return 0
}

path_setup() { # sets PT_ADDED (files edited) and PT_FOUND
  PT_ADDED=; PT_FOUND=0
  [ "$MODPATH" = 1 ] || return 0
  case ":${PATH:-}:" in *":$BIN_DIR:"*) return 0 ;; esac
  if [ -z "${HOME:-}" ]; then ag "PATH: HOME is not set; no shell profile edited"; return 0; fi
  case $ENV_FILE in
    *"'"*) warn "the prefix contains a quote; not editing shell profiles"; return 0 ;;
  esac
  if [ "$DRY" = 1 ]; then
    ag "would write $ENV_FILE:"
    env_file_text | while IFS= read -r l; do ag "    $l"; done
  else
    mkdir_tracked "$LIB_DIR"
    env_file_text >"$ENV_FILE" || { warn "could not write $ENV_FILE; PATH not edited"; return 0; }
    printf 'f %s\n' "$ENV_FILE" >>"$NEWMAN"
  fi
  path_add_to "$HOME/.profile" 1
  path_add_to "$HOME/.bash_profile" 0
  path_add_to "$HOME/.bashrc" 0
  if have zsh; then zcreate=1; else zcreate=0; fi
  path_add_to "${ZDOTDIR:-$HOME}/.zshenv" "$zcreate"
  return 0
}

profile_undo() { # manifest text: remove exactly the profile lines this installer added
  printf '%s\n' "$1" | grep -E '^[eE] ' | while IFS= read -r line; do
    kind=${line%% *}
    f=${line#? }
    if [ ! -f "$f" ]; then continue; fi
    if ! grep -Fxq -- "$PROFILE_LINE" "$f" 2>/dev/null; then
      say "  $f: the rusolve line is already gone"
      continue
    fi
    if [ "$DRY" = 1 ]; then say "  would remove the rusolve line from $f"; continue; fi
    t=$(mktemp "${TMPDIR:-/tmp}/rusolve-profile.XXXXXX") || { warn "cannot edit $f"; continue; }
    # grep exits 1 when no line is left (the file held only ours); anything above 1 is a read
    # error, and writing back a partial copy would truncate the user's profile
    rc=0; grep -vxF -- "$PROFILE_LINE" "$f" >"$t" || rc=$?
    if [ "$rc" -gt 1 ]; then
      warn "could not read $f; remove this line yourself: $PROFILE_LINE"; rm -f "$t"; continue
    fi
    if cat "$t" >"$f"; then
      say "  removed the rusolve line from $f"
      if [ "$kind" = E ] && [ ! -s "$f" ]; then rm -f "$f"; say "  removed $f (created by the installer, now empty)"; fi
    else
      warn "could not rewrite $f; remove this line yourself: $PROFILE_LINE"
    fi
    rm -f "$t"
  done
}

# ---- uninstall ------------------------------------------------------------------------------
if [ "$UNINSTALL" = 1 ]; then
  [ -f "$MANIFEST" ] || err "no manifest at $MANIFEST: nothing to uninstall under $PREFIX"
  say "Uninstalling what $MANIFEST lists:"
  # The manifest lives in a directory it lists: hold it in memory, remove it, then the rest.
  MAN=$(cat "$MANIFEST")
  [ "$DRY" = 1 ] || rm -f "$MANIFEST"
  # Files first, then directories deepest first (written parents first, so reverse them);
  # a directory is removed only if it is empty.
  agent_undo "$MAN"
  profile_undo "$MAN"
  printf '%s\n' "$MAN" | while IFS= read -r line; do
    case $line in
      "f "*|"l "*)
        p=${line#? }
        say "  remove $p"
        [ "$DRY" = 1 ] || rm -f "$p" ;;
    esac
  done
  printf '%s\n' "$MAN" | grep '^d ' | sed '1!G;h;$!d' | while IFS= read -r line; do
    p=${line#d }
    say "  rmdir  $p (if empty)"
    [ "$DRY" = 1 ] || rmdir "$p" 2>/dev/null || true
  done
  if [ "$DRY" = 1 ]; then say "(dry run: nothing removed)"; fi
  [ "$DRY" = 1 ] || say "Done."
  exit 0
fi

# ---- platform -------------------------------------------------------------------------------
os=$(uname -s 2>/dev/null || echo unknown)
mach=$(uname -m 2>/dev/null || echo unknown)
case $os in
  Linux) OSN=linux ;;
  Darwin) OSN=macos ;;
  MINGW*|MSYS*|CYGWIN*|Windows_NT)
    err "Windows is not handled by this script. Use 'pip install rusolve' or the rusolve-cli-windows-x86_64.zip from $DEFAULT_BASE_URL" ;;
  *) err "unsupported OS: $os (supported: Linux, macOS). Build from source: cargo build --release -p rusolve-cli" ;;
esac
case $mach in
  x86_64|amd64) ARCH=x86_64 ;;
  aarch64|arm64) ARCH=aarch64 ;;
  *) err "unsupported CPU: $mach (supported: x86_64, aarch64/arm64)" ;;
esac
if [ "$OSN" = macos ]; then
  # macOS reports arm64; the release names it macos-arm64
  if [ "$ARCH" = aarch64 ]; then ARCH=arm64; fi
fi
PLATFORM=$OSN-$ARCH
case $PLATFORM in
  linux-x86_64|linux-aarch64|macos-arm64|macos-x86_64) ;;
  *) err "no release build for $PLATFORM" ;;
esac
if [ "$OSN" = linux ]; then
  # The release binary is built against glibc (manylinux); musl systems (Alpine) cannot run it.
  for m in /lib/ld-musl-*.so.1 /usr/lib/ld-musl-*.so.1; do
    if [ -e "$m" ]; then
      err "this looks like a musl system (Alpine): the release binary needs glibc. Build from source instead: cargo build --release -p rusolve-cli"
    fi
  done
fi
ARCHIVE=rusolve-cli-$PLATFORM.tar.gz

if [ "$VERSION" = latest ]; then
  DL_BASE=$BASE_URL/latest/download
else
  DL_BASE=$BASE_URL/download/$VERSION
fi

# ---- tools ----------------------------------------------------------------------------------
if command -v curl >/dev/null 2>&1; then
  DL=curl
elif command -v wget >/dev/null 2>&1; then
  DL=wget
else
  err "need curl or wget"
fi
command -v tar >/dev/null 2>&1 || err "need tar"
if command -v sha256sum >/dev/null 2>&1; then
  SHA=sha256sum
elif command -v shasum >/dev/null 2>&1; then
  SHA='shasum -a 256'
else
  err "need sha256sum or shasum"
fi

fetch() { # url dest
  if [ "$DL" = curl ]; then
    curl -fsSL --retry 2 -o "$2" "$1" || err "download failed: $1"
  else
    wget -q -O "$2" "$1" || err "download failed: $1"
  fi
}

sha_of() { # file -> hex digest
  # shellcheck disable=SC2086 # $SHA is a command with an argument, split on purpose
  $SHA "$1" | awk '{print $1}'
}

# ---- decide on CP-SAT -----------------------------------------------------------------------
ORT_URL=; ORT_SHA=
if [ "$CPSAT" != no ]; then
  ortools_pin "$PLATFORM"
fi
if [ "$CPSAT" = ask ]; then
  if [ "$YES" = 0 ] && [ -t 0 ] && [ -n "$ORT_URL" ]; then
    say "Optional: CP-SAT, Google OR-Tools $ORTOOLS_VERSION, an extra engine for integer programs"
    say "(mip_engine=cpsat). Licence: Apache-2.0 plus bundled libraries; downloaded from PyPI,"
    say "about $ORTOOLS_SIZE_MB MB to download. Installs to $ORT_DIR."
    printf 'Install CP-SAT? [y/N] '
    ans=
    read -r ans || ans=
    case $ans in y|Y|yes|YES|Yes) CPSAT=yes ;; *) CPSAT=no ;; esac
  else
    CPSAT=no
  fi
fi
if [ "$CPSAT" = yes ]; then
  [ -n "$ORT_URL" ] || err "no OR-Tools build is pinned for $PLATFORM"
  if command -v unzip >/dev/null 2>&1; then
    UNZIP=unzip
  elif command -v python3 >/dev/null 2>&1; then
    UNZIP=python3
  else
    err "--with-cpsat needs unzip or python3 to open the OR-Tools wheel"
  fi
fi

# ---- plan -----------------------------------------------------------------------------------
say "rusolve installer"
say "  platform : $PLATFORM"
say "  release  : $VERSION ($DL_BASE/$ARCHIVE)"
say "  binary   : $BIN_DIR/rusolve-cli"
if [ "$CPSAT" = yes ]; then
  say "  cp-sat   : OR-Tools $ORTOOLS_VERSION from $ORT_URL"
  say "             into $ORT_DIR"
else
  say "  cp-sat   : not installed (pass --with-cpsat to add it)"
fi
if [ "$DRY" = 1 ]; then
  say "  agents   : MCP server and skills (opt out with --no-agent, --no-mcp, --no-skills)"
  say "  PATH     : added to the shell profiles unless already on PATH (opt out with --no-modify-path)"
  NEWMAN=/dev/null
  agent_setup
  path_setup
  printf '%s' "$AGOUT"
  say "(dry run: nothing downloaded or written)"
  exit 0
fi

# ---- install --------------------------------------------------------------------------------
TMP=$(mktemp -d "${TMPDIR:-/tmp}/rusolve-install.XXXXXX") || err "cannot create a temporary directory"
NEWMAN=$TMP/manifest
: >"$NEWMAN"

mkdir_tracked() { # dir: create it (and parents below PREFIX), record each new one, parents first
  if [ -d "$1" ]; then return 0; fi
  mkdir_tracked "$(dirname "$1")"
  mkdir "$1" || err "cannot create $1"
  printf 'd %s\n' "$1" >>"$NEWMAN"
}

copy_tree() { # srcdir dstdir: copy every entry, record it in the manifest, set COPIED
  COPIED=0
  for f in "$1"/*; do
    if [ ! -e "$f" ] && [ ! -L "$f" ]; then continue; fi
    b=$(basename "$f")
    rm -f "${2:?}/$b"
    if [ -L "$f" ]; then
      cp -P "$f" "$2/$b" || err "cannot write $2/$b"
      printf 'l %s\n' "$2/$b" >>"$NEWMAN"
    else
      cp "$f" "$2/$b" || err "cannot write $2/$b"
      printf 'f %s\n' "$2/$b" >>"$NEWMAN"
    fi
    COPIED=$((COPIED + 1))
  done
}

say "Downloading $ARCHIVE ..."
fetch "$DL_BASE/$ARCHIVE" "$TMP/$ARCHIVE"
fetch "$DL_BASE/SHA256SUMS" "$TMP/SHA256SUMS"
want=$(awk -v f="$ARCHIVE" '$2 == f || $2 == "*" f {print $1; exit}' "$TMP/SHA256SUMS")
[ -n "$want" ] || err "$ARCHIVE is not listed in SHA256SUMS; refusing to install"
got=$(sha_of "$TMP/$ARCHIVE")
[ "$want" = "$got" ] || err "sha256 mismatch for $ARCHIVE (expected $want, got $got); refusing to install"
mkdir "$TMP/x"
tar -xzf "$TMP/$ARCHIVE" -C "$TMP/x" rusolve-cli || err "cannot extract rusolve-cli from $ARCHIVE"
[ -f "$TMP/x/rusolve-cli" ] || err "the archive has no rusolve-cli"

write_manifest() { # this run's files and directories plus those of an earlier run, no duplicates
  mkdir_tracked "$LIB_DIR"
  {
    if [ -f "$MANIFEST" ]; then cat "$MANIFEST"; fi
    cat "$NEWMAN"
  } | awk '!seen[$0]++' >"$TMP/manifest.all"
  cp "$TMP/manifest.all" "$MANIFEST"
}

# Everything is downloaded and verified before anything is written under the prefix.
if [ "$CPSAT" = yes ]; then
  WHEEL=$TMP/ortools.whl
  say "Downloading OR-Tools $ORTOOLS_VERSION (about $ORTOOLS_SIZE_MB MB) ..."
  fetch "$ORT_URL" "$WHEEL"
  got=$(sha_of "$WHEEL")
  [ "$ORT_SHA" = "$got" ] || err "sha256 mismatch for the OR-Tools wheel (expected $ORT_SHA, got $got); refusing to install"
  mkdir "$TMP/w"
  if [ "$UNZIP" = unzip ]; then
    # only the Linux wheels have a top-level ortools.libs; asking for it on macOS makes unzip warn
    if [ "$OSN" = linux ]; then
      unzip -q "$WHEEL" 'ortools/.libs/*' 'ortools.libs/*' -d "$TMP/w" || err "cannot open the OR-Tools wheel"
    else
      unzip -q "$WHEEL" 'ortools/.libs/*' -d "$TMP/w" || err "cannot open the OR-Tools wheel"
    fi
  else
    python3 -I -m zipfile -e "$WHEEL" "$TMP/w" || err "cannot open the OR-Tools wheel"
  fi
  [ -d "$TMP/w/ortools/.libs" ] || err "the wheel has no ortools/.libs directory"
fi

WRITING=1
mkdir_tracked "$BIN_DIR"
cp "$TMP/x/rusolve-cli" "$BIN_DIR/.rusolve-cli.new.$$" || err "cannot write to $BIN_DIR"
chmod 755 "$BIN_DIR/.rusolve-cli.new.$$"
mv -f "$BIN_DIR/.rusolve-cli.new.$$" "$BIN_DIR/rusolve-cli"
printf 'f %s\n' "$BIN_DIR/rusolve-cli" >>"$NEWMAN"

if [ "$CPSAT" = yes ]; then
  mkdir_tracked "$ORT_DIR"
  copy_tree "$TMP/w/ortools/.libs" "$ORT_DIR"
  [ "$COPIED" -gt 0 ] || err "the wheel's ortools/.libs is empty"
  n=$COPIED
  # On Linux the wheel's libutf8_validity has RPATH $ORIGIN/../../ortools.libs (auditwheel's
  # renamed Abseil copies). From $ORT_DIR that is $PREFIX/lib/ortools.libs, so it goes there.
  if [ -d "$TMP/w/ortools.libs" ]; then
    mkdir_tracked "$PREFIX/lib/ortools.libs"
    copy_tree "$TMP/w/ortools.libs" "$PREFIX/lib/ortools.libs"
    n=$((n + COPIED))
  fi
  printf '%s\n' "$ORTOOLS_VERSION" >"$ORT_DIR/VERSION"
  printf 'f %s\n' "$ORT_DIR/VERSION" >>"$NEWMAN"
  ORT_LIB=
  for c in "$ORT_DIR"/libortools.so* "$ORT_DIR"/libortools*.dylib; do
    if [ -e "$c" ]; then ORT_LIB=$c; break; fi
  done
  [ -n "$ORT_LIB" ] || err "libortools is missing from the wheel"
  if [ "$OSN" = linux ] && command -v ldd >/dev/null 2>&1; then
    if ldd "$ORT_LIB" 2>&1 | grep -q 'not found'; then
      warn "some libraries OR-Tools needs are missing on this system:"
      ldd "$ORT_LIB" 2>&1 | grep 'not found' >&2
    fi
  fi
fi

agent_setup
path_setup
write_manifest
WRITING=0

# ---- report ---------------------------------------------------------------------------------
say ""
say "Installed:"
say "  $BIN_DIR/rusolve-cli"
if [ "$CPSAT" = yes ]; then
  say "  $ORT_DIR/ (OR-Tools $ORTOOLS_VERSION; $n library files in all)"
  if [ -d "$PREFIX/lib/ortools.libs" ]; then say "  $PREFIX/lib/ortools.libs/ (dependencies of the above, found through its RPATH)"; fi
fi
say "  $MANIFEST (what --uninstall removes)"
case ":${PATH:-}:" in
  *":$BIN_DIR:"*) ;;
  *)
    say ""
    if [ -n "$PT_ADDED" ] || [ "$PT_FOUND" = 1 ]; then
      say "$BIN_DIR is not on PATH in this shell. New login or interactive shells get it from"
      say "the line that sources $ENV_FILE in:"
      if [ -n "$PT_ADDED" ]; then say "  added now:$PT_ADDED"; fi
      say "For this shell now: export PATH=\"$BIN_DIR:\$PATH\"   (or: . '$ENV_FILE')"
    else
      say "$BIN_DIR is not on your PATH. Add this line to your shell profile:"
      say "  export PATH=\"$BIN_DIR:\$PATH\""
    fi ;;
esac
say ""
say "Agents:"
printf '%s' "$AGOUT"
agent_manual
say ""
say "Verify:"
say "  $BIN_DIR/rusolve-cli --version      (prints the version and whether CP-SAT is available)"
say "  $BIN_DIR/rusolve-cli --help"
say "  printf 'Maximize\\n obj: x + y\\nSubject To\\n c1: x + 2 y <= 4\\nBounds\\n 0 <= x <= 3\\n 0 <= y <= 3\\nGeneral\\n x y\\nEnd\\n' > model.lp"
say "  $BIN_DIR/rusolve-cli mip model.lp --format json"
if [ "$CPSAT" = yes ]; then
  say "  $BIN_DIR/rusolve-cli mip model.lp --opt mip_engine=cpsat --format json"
  say "  (the CP-SAT run looks for the libraries in $ORT_DIR, or where RUSOLVE_ORTOOLS_LIB points)"
fi
say "Uninstall: sh install.sh --uninstall --prefix $PREFIX"
