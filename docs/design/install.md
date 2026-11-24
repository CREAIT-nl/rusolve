# Installing rusolve outside Python: `install.sh`

`install.sh` at the repository root installs the `rusolve-cli` binary from a GitHub release and,
on request, the Google OR-Tools shared libraries that the optional CP-SAT engine
(`mip_engine=cpsat`, see [cpsat-engine.md](cpsat-engine.md)) loads at run time. Python users keep
`pip install rusolve`.

```
curl -fsSL https://raw.githubusercontent.com/CREAIT-nl/rusolve/main/install.sh | sh
sh install.sh --with-cpsat --prefix /opt/rusolve
```

## Behaviour

| Flag | Meaning |
|---|---|
| `--version vX.Y.Z` | Release to install (default: latest). |
| `--prefix DIR` | Install under `DIR` (default `~/.local`). The binary goes to `DIR/bin/rusolve-cli`. |
| `--with-cpsat` / `--no-cpsat` | Install or skip the OR-Tools libraries. Neither flag: ask once on a terminal (default No); never ask when stdin is not a terminal, which is the case for `curl ... \| sh`. |
| `--uninstall` | Remove exactly the files and the (then empty) directories listed in `DIR/lib/rusolve/install-manifest`. |
| `--dry-run` | Print the plan; download and write nothing. |
| `--no-mcp` / `--no-skills` / `--no-agent` | Skip the MCP registration, the skills plugin, or both. |
| `--no-modify-path` | Do not add `DIR/bin` to the shell profiles (see "PATH" below). |
| `--yes` | Never ask. CP-SAT stays off unless `--with-cpsat` is given. |

`RUSOLVE_INSTALL_BASE_URL` replaces `https://github.com/CREAIT-nl/rusolve/releases` (a mirror or
a fork). It must have the layout of GitHub release downloads, and the script uses one of two forms
depending on `--version`: with `--version vX.Y.Z` the files are fetched from
`<base>/download/<tag>/<file>`; without it (latest, the default) from `<base>/latest/download/<file>`.
A mirror that serves only one form works for that form only. `<file>` is the archive and
`SHA256SUMS`. Unset it after the install: `rusolve-cli` rejects `RUSOLVE_*`
variables it does not know.

## Agent integration (on by default)

After the files are written, `install.sh` looks for agent CLIs on `PATH`. Syntax was read from
each CLI's `--help` and exercised in containers (Claude Code 2.1.330, codex-cli 0.171.0,
gemini 0.70.0); none of them needs a login for these commands.

| Client | Detect | Register | Undo |
|---|---|---|---|
| claude | `claude mcp get rusolve` (exit 1 when absent; `Command:` line) | `claude mcp add --scope user rusolve -- BIN mcp` | `claude mcp remove --scope user rusolve` |
| codex | `codex mcp get rusolve` (`command:` line) | `codex mcp add rusolve -- BIN mcp` | `codex mcp remove rusolve` |
| gemini | `gemini mcp list` (line `rusolve: BIN mcp (stdio)`) | `gemini mcp add --scope user rusolve BIN mcp` | `gemini mcp remove --scope user rusolve` |
| claude skills | `claude plugin marketplace list`, `claude plugin list` | `claude plugin marketplace add SRC --scope user`, `claude plugin install rusolve-skills@rusolve --scope user` | `claude plugin uninstall rusolve-skills@rusolve --scope user`, `claude plugin marketplace remove rusolve` |

`codex mcp add` and `gemini mcp add` overwrite a server of the same name, so the script checks
first: a registration with the same binary is reported as already done, one with another command
is left alone and the replace command is printed. Clients without an `mcp add` command (Claude
Desktop, Cursor, ...) get the JSON snippet printed; their config files are never edited from sh.

Each addition the script made is a line in the manifest (`m CLIENT NAME`, `p claude PLUGIN`,
`s claude MARKETPLACE`); `--uninstall` removes exactly those, and an `m` line only while the
registration still points at this prefix's binary. An already present registration or marketplace
is never recorded, so it is never removed. A failing agent CLI makes the script warn and print the
manual command; it never fails the install. `RUSOLVE_INSTALL_MARKETPLACE` (default
`CREAIT-nl/rusolve`; an owner/repo, git URL or local path) is the skills marketplace source.

## PATH (on by default)

If `DIR/bin` is not on `PATH` when the installer runs, it makes later shells find it, the way
rustup does. It writes `DIR/lib/rusolve/env`, a POSIX `sh` file that prepends the bin directory
to `PATH` unless it is already there (sourcing it twice is harmless), and appends one line to
each profile:

```
. 'DIR/lib/rusolve/env' # added by rusolve install.sh
```

Profiles: `~/.profile` (created if missing; login `sh` and `bash` read it when there is no
`~/.bash_profile`), and `~/.bash_profile`, `~/.bashrc` and `${ZDOTDIR:-~}/.zshenv` when they
exist (`.zshenv` is also created when `zsh` is on `PATH`; zsh reads it for every invocation). A
line that is already in a file is not added again. Fish and csh are not handled. A
non-interactive `bash -c` reads none of these files.

Manifest lines: `e FILE` (the line was appended to an existing file) and `E FILE` (the file
did not exist; the installer created it). `--uninstall` deletes exactly the line (an exact
whole-line match, so a line the user wrote is never touched) from each listed file and deletes an
`E` file only if it is then empty; `--dry-run` shows the `env` file and each line it would add or
remove. Nothing is edited when the directory is already on `PATH` or with `--no-modify-path`.

Platforms: Linux x86_64 and aarch64 (glibc), macOS arm64 and x86_64. Windows and musl (Alpine)
get a message that says what to use instead. Needs `curl` or `wget`, `tar`, `sha256sum` or
`shasum`; `--with-cpsat` also needs `unzip` or `python3`.

Integrity. The CLI archive is checked against `SHA256SUMS`, which the `github-release` job of
`.github/workflows/release.yml` publishes next to the archives. The OR-Tools wheel is checked
against a SHA-256 pinned in the script. Nothing is written under the prefix until every download
has been verified; a failure later still leaves a manifest, so `--uninstall` cleans up.

Layout:

```
PREFIX/bin/rusolve-cli
PREFIX/lib/rusolve/install-manifest
PREFIX/lib/rusolve/ortools/          every file of the wheel's ortools/.libs, plus VERSION
PREFIX/lib/ortools.libs/             Linux only: the wheel's renamed Abseil copies (see below)
```

The binary looks for `libortools` in `<dir of the executable>/../lib/rusolve/ortools/`, after an
explicit `RUSOLVE_ORTOOLS_LIB`.

## Where OR-Tools comes from

Two sources were measured for 9.15.6755 (x86_64 Linux, 2026-11-23).

| | PyPI wheel (`ortools`, manylinux) | Google C++ archive (`or-tools_amd64_debian-12_cpp_v9.15.6755.tar.gz`) |
|---|---|---|
| Download | 29.8 MB (aarch64 27.6, macOS arm64 21.9, macOS x86_64 23.9) | 55.4 MB (macOS arm64 45.9) |
| Installed, shared libraries only | 67 MB in `.libs` (macOS arm64: 53 MB) | 95 MB in `lib/` (254 entries; 262 MB for the whole archive with headers and CMake files) |
| Highest glibc symbol the library needs | GLIBC_2.26 (x86_64), GLIBC_2.23 (aarch64); GLIBCXX_3.4.22 | GLIBC_2.34, GLIBCXX_3.4.30 |
| One file per platform | one wheel per platform, any CPython tag | one archive per distribution (Debian 11/12, Ubuntu 20.04/22.04/24.04, Alma, Fedora, ...) |
| Exports `SolveCpModelWithParameters` | yes | yes |
| Dependencies | vendored next to it, found by `$ORIGIN` RUNPATH, plus 9 renamed Abseil libraries (500 KB) in `ortools.libs`; no `libpython` in any file | vendored next to it, `$ORIGIN` RUNPATH |
| Checksum source | PyPI JSON API, pinned in the script | none published in the release |

Choice: the PyPI wheel. It is a third smaller, covers Linux x86_64 and aarch64 and macOS arm64
and x86_64 with four known files, runs on older systems (glibc 2.26 floor against 2.34 for the
Debian 12 archive; the other archives are built per distribution, so choosing one for the host
would need distribution detection), and ships hashes through PyPI. The wheel is a zip. Only
`ortools/.libs` and `ortools.libs` are taken, and none of those files links `libpython`
(checked with `readelf -d` for every Linux library and `otool -L` for the macOS library), so
the CPython tag in the file name does not matter and cp312 is used everywhere. The script pins
the four URLs and hashes instead of calling the PyPI JSON API at install time, which avoids
parsing JSON in `sh` and makes the hash independent of the channel that delivered it.

The Linux `libutf8_validity` carries `RPATH $ORIGIN/../../ortools.libs`. From
`PREFIX/lib/rusolve/ortools` that resolves to `PREFIX/lib/ortools.libs`, so the script puts the
wheel's `ortools.libs` there. Without it `ldd` reports nine missing `libabsl_*-<hash>` libraries.
The macOS libraries use `@loader_path`, so the single directory is enough.

OR-Tools is not redistributed: it is downloaded from PyPI on request, see
`THIRD_PARTY_NOTICES.md`. To move to a new OR-Tools, change the pinned block at the top of
`install.sh` (version, four URLs, four hashes, size) and the vendored protos of
`crates/rusolve-cpsat` together.

## Not done

* Using the installed libraries needs a binary that opens `libortools` at run time. Today the
  `cpsat` feature links it at build time (with a run path into the OR-Tools prefix), so a
  released binary cannot select `mip_engine=cpsat` until the runtime-loading change lands. The
  install layout and this document are written for that final state.
* `rusolve-cli` has no `--version`; the script's check line is `rusolve-cli --help`.
* Packaging: apt (OR-Tools is not in Debian or Ubuntu), Homebrew (a `rusolve` formula could
  depend on the existing `or-tools` formula), and a Windows installer are possible follow-ups.
