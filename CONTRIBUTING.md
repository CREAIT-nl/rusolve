# Contributing to rusolve

Thank you for helping. Bug reports, models that rusolve solves badly, documentation fixes,
ported HiGHS functions, new techniques and new skills for agents are all welcome.

## Before you start

- **Bugs:** open an issue with the model file (or the smallest model that shows the problem),
  the exact command, and the output. Models that crash rusolve, give a wrong answer or are
  much slower than HiGHS are the most useful reports there are.
- **Larger changes:** open an issue first and describe the plan, so we can agree on the
  approach before you invest in it.
- **Small fixes** (typos, docs, an obvious bug): a pull request directly is fine.

## AI-assisted contributions

Contributions written with AI coding agents are welcome, on the same terms as any other:

- **You are responsible for what you submit.** Read and understand the change before you open
  the pull request, and be ready to answer review questions about it.
- **Describe it well.** The pull request explains what changed, why, which source the change
  follows, and how it was verified (the commands you ran and what they printed).
- **Follow the patterns of the repository.** Agents should start at [`AGENTS.md`](AGENTS.md),
  which points to the per-directory instructions, the code map and the contribution recipes.
- **Say that AI was used.** The pull request template has a checkbox for it.
- **Keep it focused.** One change per pull request. Large unrequested rewrites, bulk
  reformatting and changes without a clear purpose will be closed.

## Development setup

- A Rust toolchain (stable, via [rustup](https://rustup.rs)).
- For the Python package: Python 3.11 or newer and `pip install maturin`.
- For the differential tests against HiGHS (`tools/diff`): `pip install highspy pyscipopt`.

Build:

    cargo build --release -p rusolve-cli

## Checks

These are the checks CI runs. Run them before you open a pull request:

    cargo fmt --all --check
    cargo clippy -p rusolve-ekk -p rusolve-cli -p rusolve-wasm --release --tests --benches -- -D warnings
    cargo test --release -p rusolve-ekk --lib
    cargo test --release -p rusolve-cli --no-fail-fast -- --skip sp150x300d
    tools/agent_ci.sh                     # generated docs are current, every skill's eval graders pass

If you changed the Python package:

    python -m venv .venv && . .venv/bin/activate
    pip install maturin pytest numpy scipy jsonschema mypy
    maturin develop --release -m crates/rusolve-py/Cargo.toml
    mypy crates/rusolve-py/python/rusolve
    pytest -q crates/rusolve-py/tests

Some tests need benchmark instances that are downloaded rather than tracked
(`python tools/bench/fetch_corpus.py`, see `benches/data/ATTRIBUTION.md`); they are skipped
when the files are missing.

## Parity check

The default search is deterministic: the same model, options and seed give the same pivots,
nodes and result. A change that is not meant to change results (a refactor, performance work,
a new option that is off by default) must leave that untouched. To check, build `main` and your
branch, then compare them on 20 MIPLIB instances:

    # one-time: download the 20 instances (see the header of tools/bench/identity.sh)
    git worktree add ../rusolve-main main
    (cd ../rusolve-main && cargo build --release -p rusolve-cli)
    cargo build --release -p rusolve-cli
    tools/bench/identity.sh ../rusolve-main/target/release/rusolve-cli target/release/rusolve-cli

The instances are read from `$RUSOLVE_MIPLIB_DIR` (default
`~/.cache/rusolve-bench/miplib-benchmark`).

It must print `IDENTICAL 20 of 20`. If your change is meant to alter the search, say so in the
pull request and include before-and-after numbers on the instances it affects. The
[`contrib-identity-abtest`](docs/skills/contrib-identity-abtest/SKILL.md) skill walks through
both cases.

## Rules for solver code

- **Ported code stays faithful.** Code translated from HiGHS keeps the same algorithm,
  constants and control flow, and a header naming the HiGHS source file and its MIT notice. A
  difference in behaviour from HiGHS is a bug. Tests for ported code use golden fixtures dumped
  from HiGHS itself (`tools/highs-oracle`). The
  [`contrib-port-function`](docs/skills/contrib-port-function/SKILL.md) skill shows how.
- **Improvements go behind an option** first, off by default, and become the default only with
  benchmark evidence.
- **Cite the source** of every algorithmic idea taken from elsewhere: a paper, a textbook, public
  documentation, or the source of an open-source solver whose license allows reuse. Put the
  citation in a comment next to the code. An original idea needs no citation: explain how it
  works in the comment, and give the benchmark evidence in the pull request.
- **Only permissively licensed code** (MIT, BSD, Apache-2.0, or similar) may be translated or
  copied in, with its notice in the file header and in
  [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). Never copy GPL code or code from
  proprietary software.
- **Performance is measured, never claimed.** A change to a hot path comes with a benchmark
  (criterion, or a before-and-after run on named instances).
- `cargo fmt`, `cargo clippy -- -D warnings` and the tests are clean. `unsafe` needs a written
  justification in the pull request.
- **Generated files** (`docs/agent/reference.md`, `docs/agent/diagnostics.md`, `docs/agent/codemap.md`) are
  regenerated with `python3 tools/gen_agent_docs.py` and `python3 tools/gen_codemap.py`, never
  edited by hand.

## Commits and pull requests

- Commit messages say *why* the change is made; the diff already shows what.
- Pull requests target `main`, pass CI, and are reviewed by a maintainer before they are
  merged. They are squash-merged.
- Fill in the pull request template.

## License

rusolve is distributed under the [Business Source License 1.1](LICENSE). By submitting a contribution
you confirm that you have the right to submit it, and you license it to the maintainers under the
[Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0), so that it can be
distributed under the project's license, including any license the project adopts later.

## Conduct

Everyone taking part follows the [code of conduct](CODE_OF_CONDUCT.md).
