# Installing rusolve (instructions for an AI agent)

You were pointed at this file to install rusolve on the user's machine. Do the steps in order.
Decide everything yourself from what the user told you and from what the machine has; the one
question you may have to ask is in section 6 (a different `rusolve` MCP server already
registered). Do not use `sudo`; every route below installs for the current user.

rusolve is an LP and MIP solver. One binary, `rusolve-cli`, is both the command-line tool and a
stdio MCP server (`rusolve-cli mcp`). The Python package `rusolve` bundles the same binary.

## 1. Choose the route

Run:

    uname -s; uname -m; ldd --version 2>&1 | head -1
    for t in curl wget tar sha256sum shasum unzip python3 pip uv uvx cargo claude codex gemini; do
      p=$(command -v "$t") && echo "$t: $p" || echo "$t: MISSING"
    done

What the output means:

- `curl` or `wget`, `tar`, and `sha256sum` or `shasum` (either of each pair suffices) are what
  route A needs. `unzip` or `python3` (either) is needed only with `--with-cpsat`.
- `pip`, `uv` or `uvx` (any one) is what route B needs; `cargo` is what route C needs.
- `claude`, `codex`, `gemini` decide what route A registers: each one that is present gets the
  MCP server, and only `claude` also gets the skills plugin. With none of them present, nothing
  is registered and you give the user a JSON snippet (section 2).

| If | Route |
|---|---|
| `Linux` with glibc (the `ldd` line says GNU libc / glibc), `x86_64` or `aarch64`; or `Darwin` with `arm64` or `x86_64` | **A. `install.sh`** (binary, MCP registration, skills). |
| Windows (`uname` missing or `MINGW*`/`MSYS*`), or the user works in Python and wants `import rusolve` | **B. Python**: `pip install rusolve` |
| Linux with musl (`ldd --version` says musl, Alpine), another CPU, or A failed because no release is published | **C. From source**: needs a Rust toolchain |

What each route gives the user:

| | CLI (`rusolve-cli`) | MCP server (all tools, including CP objects through `build_model`'s `cp` section) | `import rusolve` (the Python `Model` API, `add_interval`, ...) | Registers the agents, installs the skills, CP-SAT libraries |
|---|---|---|---|---|
| A | yes | yes | no | yes (`--with-cpsat` for CP-SAT) |
| B | yes, as `rusolve` | yes (`rusolve mcp`) | yes | no: you register by hand; CP-SAT is `pip install 'rusolve[cpsat]'` |
| C | yes | yes | no | no: by hand |

Route A is enough for an agent that builds and solves models through MCP or the CLI, scheduling
and other constraint-programming (CP) models included: the `cp` section of `build_model` is the
equivalent of `Model.add_interval` and friends (scheduling skill, `docs/skills/scheduling`).
Choose B when the user writes Python (scripts, notebooks) with `import rusolve`, or is on
Windows. Install both (A and B) when the user wants both: an agent set up by the installer (MCP
registered, skills installed) and `import rusolve` in their own code. They do not conflict;
the commands are `rusolve-cli` (A) and `rusolve` (B).

## 2. Route A: `install.sh`

Decide one thing first: whether to install CP-SAT, an optional second engine for integer
programs (`mip_engine=cpsat`, Google OR-Tools, licence Apache-2.0; a 30 MB download from PyPI
into the install prefix; never needed for LP, MIP, QP or conic models).

The rule: if the user's stated use involves scheduling, timetabling, sequencing, rostering or
other constraint programming (job shop, nurse rostering, no-overlap, resource-constrained
projects), pass `--with-cpsat`; do not ask. Otherwise pass `--no-cpsat`, and in your report to the
user say that CP-SAT can be added later by running the same command again with `--with-cpsat`.

Run (installs to `~/.local/bin`; shown without CP-SAT):

    curl -fsSL https://raw.githubusercontent.com/CREAIT-nl/rusolve/main/install.sh | sh -s -- --yes --no-cpsat

or, from a clone of this repository, `sh install.sh --yes --no-cpsat`. With CP-SAT, replace
`--no-cpsat` by `--with-cpsat`. To see first what it would do, add `--dry-run`: it prints the
commands and the profile lines, and changes nothing.

Flags (all optional):

| Flag | Effect |
|---|---|
| `--yes` | Never ask. Always pass it. |
| `--with-cpsat` / `--no-cpsat` | Install or skip the OR-Tools libraries. |
| `--prefix DIR` | Install under `DIR` instead of `~/.local`. |
| `--version vX.Y.Z` | A specific release (default: latest). |
| `--no-mcp` | Do not register the MCP server with agent CLIs. |
| `--no-skills` | Do not install the skills plugin for Claude Code. |
| `--no-agent` | Both of the above. Use it if the user said they do not want their agent configuration touched. |
| `--no-modify-path` | Do not add the bin directory to the shell profiles (see below). Use it if the user said they do not want their shell configuration touched. |

What the defaults do, after the binary is installed:

- For each of `claude`, `codex`, `gemini` found on `PATH`, registers an MCP server named
  `rusolve` at user scope with the absolute path of the installed binary. An existing server
  named `rusolve` is left alone and reported (never overwritten); re-running is safe.
- If `claude` is found, adds the plugin marketplace `CREAIT-nl/rusolve` and installs the plugin
  `rusolve-skills@rusolve` (the playbook skills for formulating, solving and tuning models).
- Prints the JSON snippet for MCP clients that have no command line (Claude Desktop, Cursor,
  Windsurf, ...). You cannot configure those from a shell; give the snippet to the user.
- A failure in any of this prints a warning and the manual command; the install itself still
  succeeds. Read the "Agents:" block in the output and run any manual command it prints.

- If the bin directory (`~/.local/bin`) is not already on `PATH`, makes it so for later shells,
  the way rustup does: it writes `~/.local/lib/rusolve/env` and appends one line,
  `. '<prefix>/lib/rusolve/env' # added by rusolve install.sh`, to `~/.profile` (created if
  missing) and to `~/.bashrc`, `~/.bash_profile` and `~/.zshenv` where those exist (`~/.zshenv` is
  also created when `zsh` is installed). Running the installer again adds nothing twice. Fish and
  csh are not handled.

Everything it registers or edits (profile lines included) is recorded in
`~/.local/lib/rusolve/install-manifest`.

`RUSOLVE_INSTALL_BASE_URL` (a release mirror) and `RUSOLVE_INSTALL_MARKETPLACE` (the skills
marketplace: an `owner/repo`, a git URL or a path) are for forks and tests. Leave them unset.

**PATH.** The profile edit takes effect in shells started afterwards (new terminals, new login
shells), not in the shell you are running in. For your own session either run
`export PATH="$HOME/.local/bin:$PATH"` or call the binary by its full path
(`~/.local/bin/rusolve-cli`); the commands below assume you did one of these. Tell the user:
"I added `~/.local/bin` to your PATH through `~/.profile` (and the other shell startup files
listed in the output); open a new terminal to pick it up. To undo it, run the installer with
`--uninstall`." With `--no-modify-path`, or when the directory is already on `PATH`, nothing is
edited; say that instead. A non-interactive `bash -c` does not read the profiles, so a tool that
launches `rusolve-cli` that way still needs the full path or its own PATH setting (an MCP
registration already uses the full path).

**Skills source.** The skills come from a Claude Code plugin marketplace. With the default source
(`CREAIT-nl/rusolve` on GitHub), Claude Code (checked with 2.1.330) clones it into
`~/.claude/plugins/marketplaces/rusolve` and copies the plugin into
`~/.claude/plugins/cache/rusolve/rusolve-skills/<version>`; `claude plugin list` then shows no
`Read from:` line, and deleting any checkout you used changes nothing. A marketplace added from a
local path (`RUSOLVE_INSTALL_MARKETPLACE=/some/path`, shown as `Source: Folder (/some/path)` by
`claude plugin marketplace list`) is different: Claude Code reads it live from that path
(`claude plugin list` shows `Read from: /some/path`). If the path is deleted or moved, the plugin
shows `Status: ✘ failed to load` with `Marketplace rusolve failed to load: cache-miss`, and the
skills are gone. Never point it at a temporary directory.

A newly registered MCP server or skills plugin is loaded by an agent session when it starts:
tell the user to restart their agent (or run `/mcp` / `/reload-plugins` in Claude Code) to
see them in the session that is already open.

## 3. Route B: Python

    pip install rusolve                      # or: uv pip install rusolve
    pip install 'rusolve[cpsat]'             # CP-SAT: same rule as in section 2

The wheel provides the `rusolve` command (same commands as `rusolve-cli`). Register the MCP
server with the agents you have, for example:

    claude mcp add --scope user rusolve -- uvx --from rusolve rusolve mcp
    codex mcp add rusolve -- uvx --from rusolve rusolve mcp
    gemini mcp add --scope user rusolve uvx --from rusolve rusolve mcp

(Skip the ones whose CLI is absent. Without `uvx`, use the absolute path of the `rusolve`
command from `command -v rusolve` with the argument `mcp`.) For the skills, run
`claude plugin marketplace add CREAIT-nl/rusolve --scope user` and
`claude plugin install rusolve-skills@rusolve --scope user`.

In the commands below, replace `rusolve-cli` by `rusolve` on this route.

## 4. Route C: from source

    git clone https://github.com/CREAIT-nl/rusolve && cd rusolve
    cargo build --release -p rusolve-cli
    install -D target/release/rusolve-cli ~/.local/bin/rusolve-cli

Then register the MCP server and skills as in route B, with the path
`$HOME/.local/bin/rusolve-cli` and the argument `mcp`.

## 5. Verify

Use the binary's full path if it is not on `PATH`. Do all of these.

1. Version and capabilities:

       rusolve-cli --version

   Expected: a line `rusolve <version> ...` and a line `cpsat: available, OR-Tools ...` or
   `cpsat: not available ...`. For scripts: `rusolve-cli version --format json` has
   `cpsat.available` (true/false). `cpsat.available` must be true if and only if you installed
   CP-SAT. If `--version` is not a command, the binary predates it: reinstall the latest
   release.

2. A solve:

       printf 'Maximize\n obj: x + y\nSubject To\n c1: x + 2 y <= 4\nBounds\n 0 <= x <= 3\n 0 <= y <= 3\nGeneral\n x y\nEnd\n' > /tmp/rusolve-check.lp
       rusolve-cli mip /tmp/rusolve-check.lp --format json

   Expected: one JSON object on stdout with top-level keys `"status": "Optimal"`,
   `"objective": 3` (x=3, y=0, since both variables are integer) and `"engine": "highs-port"`,
   exit code 0. Anything else is a failure; read the `diagnostics` array (see section 6).

3. CP-SAT (only if you installed it):

       rusolve-cli mip /tmp/rusolve-check.lp --opt mip_engine=cpsat --format json

   Expected: top-level `"status": "Optimal"`, `"objective": 3` and `"engine": "cp-sat"` (the
   top-level key `.engine`; check that one). Do not use `.solver.engine`: it names the underlying
   solver code (`highs-port`) on every run, CP-SAT or not.

4. MCP registration (Claude Code):

       claude mcp get rusolve

   Expected: `Scope: User config`, `Status: ✔ Connected` and `Command:` with the absolute path.
   Other agents: `codex mcp get rusolve`, `gemini mcp list`. The handshake can also be checked
   without any agent:

       printf '%s\n' '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-11-25","capabilities":{},"clientInfo":{"name":"check","version":"0"}}}' '{"jsonrpc":"2.0","method":"notifications/initialized"}' '{"jsonrpc":"2.0","id":2,"method":"tools/list"}' | rusolve-cli mcp

   Expected: two JSON lines; the second lists `load_model`, `build_model`, `solve`,
   `get_solution`, and the other tools.

   If the user does scheduling or other CP work, check the CP objects through the same
   handshake. This builds two jobs on one machine and solves them (the server numbers handles
   `m1`, `m2`, ... from 1 in each process, so the first model built is `m1`):

       printf '%s\n' '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-11-25","capabilities":{},"clientInfo":{"name":"check","version":"0"}}}' '{"jsonrpc":"2.0","method":"notifications/initialized"}' '{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"build_model","arguments":{"variables":[{"name":"s1","type":"integer","lower":0,"upper":9},{"name":"e1","type":"integer","lower":0,"upper":9},{"name":"s2","type":"integer","lower":0,"upper":9},{"name":"e2","type":"integer","lower":0,"upper":9},{"name":"mk","type":"integer","lower":0,"upper":9,"objective":1}],"rows":[{"name":"mk1","terms":{"mk":1,"e1":-1},"sense":">=","rhs":0},{"name":"mk2","terms":{"mk":1,"e2":-1},"sense":">=","rhs":0}],"cp":{"intervals":[{"name":"a","start":"s1","size":2,"end":"e1"},{"name":"b","start":"s2","size":3,"end":"e2"}],"no_overlap":[{"name":"one","intervals":["a","b"]}]}}}}' '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"solve","arguments":{"model_id":"m1"}}}' | rusolve-cli mcp

   Expected: three JSON lines (long ones). In the second, `structuredContent.model_id` is `m1`
   (the same JSON is also in `content[0].text`, as an escaped string). In the third,
   `structuredContent` has `"status":"Optimal"`, `"objective":5.0`, `"engine":"cp-linearization"`
   and, inside `cp`, `"verified":true`. Add `"options":{"mip_engine":"cpsat"}` to the solve arguments to run it on
   CP-SAT (needs `--with-cpsat`).

5. Skills (Claude Code). Neither check needs a login; both read what is installed on disk:

       claude plugin list
       claude plugin details rusolve-skills@rusolve

   Expected: `rusolve-skills@rusolve` with `Status: ✔ enabled`, and in `details` a
   `Skills (N)` line that includes `rusolve-start-here` and `rusolve-scheduling`. If the plugin is
   `failed to load`, see "Skills source" in section 2. These commands show the install record
   only. Whether a running agent session has the skills can be checked from inside that session
   (its skill list), so tell the user to restart their agent.

6. PATH in a new shell (if the installer edited the profiles): `bash -lc 'command -v rusolve-cli'`
   prints the path of the binary. (This starts a login shell, which reads `~/.profile`.)

Report to the user what you installed, the `--version` line, and anything you could not do.

## 6. Errors you can meet

| Message | Cause and fix |
|---|---|
| `install.sh: error: no release build for ...` or `download failed: ...releases/latest/download/...` | No release is published for this platform (or yet). Use route B or C. |
| `install.sh: error: this looks like a musl system` | The release binary needs glibc. Route C. |
| `install.sh: error: sha256 mismatch ...` | The download is corrupt or tampered with. Retry once; if it persists, stop and tell the user. |
| `install.sh: warning: claude: could not ...` | Agent registration failed; the install is complete. Run the manual command printed in the "Agents:" block. |
| `an MCP server 'rusolve' is already registered with another command; left alone` | The user has a different rusolve. Ask whether to replace it; if yes, run the `remove && add` line printed. |
| `RS0005` unknown environment variable | A `RUSOLVE_*` variable in the environment is misspelled or retired (for example `ORTOOLS_PREFIX` is replaced by `RUSOLVE_ORTOOLS_LIB`). `rusolve-cli env` lists the valid ones; unset the bad one. |
| `RS0012` with `mip_engine=cpsat` | CP-SAT is not installed or not found. `rusolve-cli version --format json` shows `cpsat.error` and the places searched. Install it: `sh install.sh --yes --with-cpsat` (or `pip install 'rusolve[cpsat]'`), or set `RUSOLVE_ORTOOLS_LIB` to the library file. The model must also have only integer columns with finite bounds. |
| `Status: ✘ Failed to connect` in `claude mcp get` | The binary path is wrong or the binary cannot start. Run the binary by hand (`rusolve-cli --version`) and read its error. |

`rusolve-cli --explain RS0005` explains any code.

## 7. Uninstall

    sh install.sh --uninstall            # from a clone, or: curl -fsSL <install.sh URL> | sh -s -- --uninstall
    sh install.sh --uninstall --dry-run  # shows the commands first

It removes the binary, the CP-SAT libraries, the `env` file, exactly the profile lines the
installer added (every other line of those files stays; a profile file the installer created is
deleted only if nothing else is in it), and exactly the MCP registrations, plugin and marketplace
that the installer added (it does not touch a registration it did not add, or one that now points
elsewhere). Pass `--prefix DIR` if you installed with one. Route B: `pip
uninstall rusolve` and the `mcp remove` command of each agent.

## 8. What to read next

- [`llms.txt`](llms.txt): the index of everything an agent needs to use rusolve.
- [`docs/skills/start-here`](docs/skills/start-here/SKILL.md): classify a problem and pick the
  playbook.
- [`docs/design/install.md`](docs/design/install.md): how the installer works.
