//! The registry of every environment variable the rusolve crates read
//! (`docs/agent/reference.md` is generated from it).
//!
//! One entry per variable, one source of truth. Three things read this table:
//!
//! * `rusolve env [--json]` ([`table_text`], [`table_json`]) — discovery: an
//!   agent asks the solver what it honours instead of grepping for it;
//! * [`unknown_env`] — `main()` refuses to start on a `RUSOLVE_*`
//!   variable that is not here, so a typo is an error rather than a silent
//!   no-op or a silent engine change;
//! * [`legacy_engine_lever`] — `use_ekk_engine` switches `rusolve solve` to the
//!   legacy simplex only for variables the legacy engine alone honours
//!   (`EnvScope::LegacyLp` with `EnvKind::Lever`), replacing the old inverse
//!   allowlist that flagged every other `RUSOLVE_*` name.
//!
//! # Kind
//!
//! * `Lever` — setting it changes the search, or changes which engine has to
//!   run for the variable to mean anything. The second half is why some
//!   counter-only variables are `Lever` and not `Debug`: a knob only the legacy
//!   engine reads (RUSOLVE_DUAL_STATS, RUSOLVE_DEBUG) does nothing on ekk, so
//!   honouring it *is* an engine choice, and census scripts rely on
//!   that switch. `Debug` is output the path already running can produce.
//! * `Profile` — a named bundle of options (RUSOLVE_PROFILE).
//! * `Debug` — extra output: traces, counters, dumps, debug-solution checks.
//!   Never changes the search.
//! * `Harness` — test and benchmark plumbing: fixture paths, seed counts,
//!   extra bench instances. Never read by a shipped solve.
//!
//! # Scope
//!
//! Which path honours the variable: `LegacyLp` (the legacy `rusolve solve`
//! simplex and IPM in crates/rusolve-lp), `Ekk` (the HiGHS-port LP), `PortMip`
//! (the `RUSOLVE_MIP_HIGHS=1` path), `NativeMip` (crates/rusolve-mip), `Any`
//! (more than one path, or the engine choice itself).
//!
//! # Keeping it complete
//!
//! `env_registry_covers_every_read` walks `crates/` and fails if any source
//! reads a variable this table lacks, or lists an entry no source reads. Add
//! the entry next to the read and the gate passes.

use serde::Serialize;

/// What kind of thing an environment variable is (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EnvKind {
    /// Changes the search, or the engine that must run to honour it.
    Lever,
    /// A named option bundle.
    Profile,
    /// Extra output only; never changes the search.
    Debug,
    /// Test or benchmark plumbing.
    Harness,
}

impl EnvKind {
    /// The label `rusolve env` prints.
    pub fn label(self) -> &'static str {
        match self {
            EnvKind::Lever => "Lever",
            EnvKind::Profile => "Profile",
            EnvKind::Debug => "Debug",
            EnvKind::Harness => "Harness",
        }
    }
}

/// Which solve path honours a variable (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum EnvScope {
    /// The legacy `rusolve solve` simplex and IPM (crates/rusolve-lp).
    #[serde(rename = "legacy_lp")]
    LegacyLp,
    /// The rusolve-ekk HiGHS-port LP.
    #[serde(rename = "ekk")]
    Ekk,
    /// The port MIP path (`RUSOLVE_MIP_HIGHS=1`).
    #[serde(rename = "port_mip")]
    PortMip,
    /// The native MIP (crates/rusolve-mip).
    #[serde(rename = "native_mip")]
    NativeMip,
    /// More than one path, or the engine choice itself.
    #[serde(rename = "any")]
    Any,
}

impl EnvScope {
    /// The label `rusolve env` prints.
    pub fn label(self) -> &'static str {
        match self {
            EnvScope::LegacyLp => "LegacyLp",
            EnvScope::Ekk => "Ekk",
            EnvScope::PortMip => "PortMip",
            EnvScope::NativeMip => "NativeMip",
            EnvScope::Any => "Any",
        }
    }
}

/// One environment variable: who reads it, what it is, and what it does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct EnvVar {
    /// The variable's exact name.
    pub name: &'static str,
    /// Every source file that reads it, repo-relative.
    pub read_by: &'static [&'static str],
    /// [`EnvKind`]: lever, profile, debug output or harness plumbing.
    pub kind: EnvKind,
    /// [`EnvScope`]: which solve path honours it.
    pub scope: EnvScope,
    /// What the values mean, in a phrase.
    pub values: &'static str,
    /// One sentence saying what the code that reads it does with it.
    pub doc: &'static str,
}

impl EnvVar {
    /// `true` when setting this variable has to move `rusolve solve` off the
    /// ekk default: only the legacy engine honours it, and it is not plumbing.
    pub fn is_legacy_engine_lever(&self) -> bool {
        self.kind == EnvKind::Lever && self.scope == EnvScope::LegacyLp
    }
}

/// Every environment variable any crate under `crates/` reads.
pub const ENV_VARS: &[EnvVar] = &[
    EnvVar {
        name: "P42_COLD",
        read_by: &["crates/rusolve-lp/src/lu.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Times the cold `factor_sparse` path in the `factor_ws_timing` probe instead of the ",
            "recycled-workspace path."
        ),
    },
    EnvVar {
        name: "P42_REPS",
        read_by: &["crates/rusolve-lp/src/lu.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a repetition count (default 50)",
        doc: "Repetitions of the ignored `factor_ws_timing` timing probe in rusolve-lp's LU tests.",
    },
    EnvVar {
        name: "PD_DETAIL",
        read_by: &["crates/rusolve-cli/src/mip_highs.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::PortMip,
        values: "set = on",
        doc: "`rusolve presolve-dump`: appends the rows and columns each recorded reduction touched.",
    },
    EnvVar {
        name: "RUSOLVE_AGG_FUZZ_FROM",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a seed number (default 1)",
        doc: concat!(
            "Seed the aggregator differential fuzz test starts from; it checks presolve-ON against ",
            "presolve-OFF on random equality-heavy LPs. rusolve-lp unit test only; never read by a ",
            "shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_AGG_FUZZ_SEEDS",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a seed count (default 400)",
        doc: concat!(
            "Last seed (so seed count) of the aggregator differential fuzz test against presolve off. ",
            "rusolve-lp unit test only; never read by a shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_AUDIT",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::NativeMip,
        values: "set = on, 0 = off",
        doc: concat!(
            "Native `rusolve mip`: prints `audit ...` stderr lines and per-LP timings for the gap ",
            "audit, without changing the search."
        ),
    },
    EnvVar {
        name: "RUSOLVE_AUTOPSY_EVERY",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an iteration period (default 64)",
        doc: concat!(
            "With RUSOLVE_DEBUG=1, run the dual autopsy (ground-truth recomputation of the dual ",
            "state) every N iterations instead of 64. Read only by the legacy engine, so setting it ",
            "also selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_BENCH_MPS",
        read_by: &["crates/rusolve-lp/benches/ipm_chol.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "comma-separated MPS paths",
        doc: concat!(
            "Extra instances the ipm_chol criterion benchmark adds to its matrix set; unreadable ",
            "paths are skipped."
        ),
    },
    EnvVar {
        name: "RUSOLVE_BIN",
        read_by: &["tools/diff/dual_diff.py", "tools/diff/mip_diff.py"],
        kind: EnvKind::Harness,
        scope: EnvScope::Any,
        values: "a path to a rusolve-cli binary",
        doc: concat!(
            "Which binary the tools/diff differential scripts run. The solver itself ignores it; it is ",
            "registered because the scripts' children inherit it."
        ),
    },
    EnvVar {
        name: "RUSOLVE_BLAND_AFTER",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a pivot count (default 25)",
        doc: concat!(
            "Stall threshold in pivots after which the Bland anti-cycling rule engages; a huge value ",
            "keeps Bland out of the solve entirely."
        ),
    },
    EnvVar {
        name: "RUSOLVE_BLAND_EPISODE",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a pivot count (default 25)",
        doc: concat!(
            "Budget of the bounded Bland episode: Bland is released early on non-degenerate pivots ",
            "and re-entry needs doubled thresholds; 0 restores sticky Bland."
        ),
    },
    EnvVar {
        name: "RUSOLVE_BRANCH",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::NativeMip,
        values: "a `BranchConfig::apply` spec string",
        doc: concat!(
            "Branching-spec string applied to the native MIP's `BranchConfig` before its flags, for ",
            "A/B runs through scripts that pass no flags; a bad spec is a usage error."
        ),
    },
    EnvVar {
        name: "RUSOLVE_BT_GATE",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a ratio (default 4.0)",
        doc: concat!(
            "Cost gate on the BTRAN eta pattern walk: it switches to the pattern walk when that is ",
            "this factor cheaper; bit-exact either way, so it moves cost, not the pivot path."
        ),
    },
    EnvVar {
        name: "RUSOLVE_CUT_DEBUG",
        read_by: &["crates/rusolve-mip/src/cuts/mod.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::NativeMip,
        values: "set = on",
        doc: concat!(
            "Native MIP cut generation: prints per-separator cut counts, violations and stripping ",
            "decisions to stderr."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DEBUG",
        read_by: &[
            "crates/rusolve-lp/src/dual.rs",
            "crates/rusolve-lp/src/presolve.rs",
            "crates/rusolve-lp/src/simplex.rs",
        ],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy-engine trace output: presolve rule applications, simplex and dual audits, ",
            "SH-autopsy and terminal-state checks. Read only by the legacy engine, so setting it also ",
            "selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DEBUG_SOL_ABORT",
        read_by: &["crates/rusolve-ekk/src/mip/debug_sol.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "set to exactly 1",
        doc: "Port MIP: makes a violated HighsDebugSol check panic instead of logging to stderr.",
    },
    EnvVar {
        name: "RUSOLVE_DENSE_LU",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "A/B lever back to the dense LU factorization for the legacy simplex, whose default is ",
            "the sparse engine."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DPROG",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy dual simplex: prints a progress line every 500 iterations. Read only by the ",
            "legacy engine, so setting it also selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DTRACE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy dual simplex: prints `dtrace` lines for the z identity, the leaving row and the ",
            "entering column of each pivot. Read only by the legacy engine, so setting it also ",
            "selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ACC",
        read_by: &[
            "crates/rusolve-cli/src/main.rs",
            "crates/rusolve-lp/src/dual.rs",
        ],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Dual accuracy regime; the CLI default is on and `=0` is the lever back to the ",
            "earlier accuracy handling."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ACC_PARTS",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a bitmask integer",
        doc: concat!(
            "Selects which parts of the dual accuracy regime stay on (bit 0 skips the refinement ",
            "rounds), used to price each part separately."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_AUDIT_INC",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an audit-pricing spec string",
        doc: concat!(
            "Lever: prices the dual audit incrementally instead of recomputing every audited ",
            "reduced cost."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_BACKTRACK",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: "Enables backtracking when a dual ratio step turns out not to be a genuine decrease.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_CORRECT",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: "Lever: repairs accumulated dual infeasibility by an explicit correction step.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_CTIE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Breaks CHUZC ties (equal |a| in the pivot group) by a hashed permutation of the variable ",
            "index instead of the lowest index."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_DEBTSHIFT",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: concat!(
            "Lever: a flip-only/saturate step's dual advance is treated as a debt only a pivot on ",
            "the same row pays, instead of being kept as the edge offset."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_DINF",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Diagnostic: traces the sign of the leaving row's dual before the step and the row's ",
            "infeasibility after it. Read only by the legacy engine, so setting it also selects that ",
            "engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_DSE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Dual steepest-edge weights (Forrest & Goldfarb 1992); `=0` is the lever back to Dantzig ",
            "pricing."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_FLIPPIV",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = flip-only, 1 = pivot in warm too",
        doc: concat!(
            "Lever on the flip-prefix exits: `=0` keeps the flip-only offset commit, `=1` also ",
            "pivots in warm (MIP node) re-solves."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_HANDOVER",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Hands a stalling solve over to the fresh-pi path at the first ladder rung, as a safety ",
            "net on the incremental-d path."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_HARRIS",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a Harris-relaxation factor",
        doc: concat!(
            "Lever: applies Harris' relaxed pivot test to the dual entering column above ",
            "the given strength."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_HARRIS_GATE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a ratio",
        doc: "Lever: the best-pivot gate below which the Harris relaxation is applied.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_HCHUZR",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = never, 1 = large LPs, 2 = always",
        doc: concat!(
            "Lever: hyper CHUZR (choose the leaving row from a hyper-sparse candidate set) ",
            "instead of the full violation scan."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_HGROUP",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: "Lever: HiGHS' Harris-relaxed BFRT grouping for the dual ratio test.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_INCD",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Keeps reduced costs incrementally instead of recomputing them from pi every pass, on LPs ",
            "above the nnz cutoff."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_INCP",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: concat!(
            "Lever: carries the dual point forward incrementally instead of re-anchoring pi by a ",
            "full BTRAN of the basic costs every iteration."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_INCP_CAD",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a re-anchor period",
        doc: "Cadence of the periodic pi re-anchor on the incremental-dual path.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_INF_LIST",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a fraction of m (default 0.05, 0 disables)",
        doc: concat!(
            "Hyper-sparse CHUZR (Hall & McKinnon 2005): take the leaving row from the ",
            "infeasible-basics list while it is at most this fraction of the rows."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_KEEPSHIFT",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Keeps a cost shift on a column after it enters the basis (as HiGHS does); `=0` restores ",
            "release-on-entry."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_LARGE_ALPHA",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: concat!(
            "HiGHS `chooseFinalLargeAlpha`: when the BFRT walk's final group holds no usable pivot, ",
            "back off to the last earlier group with a large alpha."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_NOPERT",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (disables perturbation)",
        doc: concat!(
            "Skips the hash-signed sub-tolerance cost perturbation that otherwise engages ",
            "once a plateau/stall survives the tie-scoring overlay."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_PERT",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: "HiGHS' up-front cost perturbation of cold dual solves; `=0` is the lever off.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_PERT_SEED",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an integer seed (default 0)",
        doc: "Seed of the HiGHS-style cost perturbation hash, for A/B runs over perturbation seeds.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_PICHECK",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "At every pi re-anchor on the incremental path, measures the drift the anchor erases and ",
            "reports it (a measurement lever)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_RESCAN_SKIP",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: "Lever: skips the fresh-pi ratio-pass rescan for columns already priced this pass.",
    },
    EnvVar {
        name: "RUSOLVE_DUAL_RFCAD",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Applies the dual engine's refactor-clock scale (4.0) to its factorizations; `=0` keeps ",
            "the primal's calibrated clock."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ROWPRICE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Row-wise PRICE for the dual's duals update on large LPs; `=0` falls back to column-wise ",
            "PRICE."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_RTIE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Breaks CHUZR ties (equal DSE score) by a hashed permutation of the row index instead of ",
            "the lowest index."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_SHIFT0",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, 1 = no shiftBack, 2 = default",
        doc: concat!(
            "Zero-step cost shift mode: 2 shifts and shifts back, 1 shifts only, 0 is the lever ",
            "off."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_STATS",
        read_by: &[
            "crates/rusolve-lp/src/dual.rs",
            "crates/rusolve-lp/src/pfi.rs",
            "crates/rusolve-lp/src/lu.rs",
        ],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy engine: per-solve cost counters for the dual, PFI and factorization, dumped as ",
            "one line each to stderr. Read only by the legacy engine, so setting it also selects that ",
            "engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_TA",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: concat!(
            "HiGHS' update-count dependent CHUZC pivot tolerance (HEkkDualRow::choosePossible) in ",
            "place of the flat pivot tolerance."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_WATCH",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a column index",
        doc: concat!(
            "Diagnostic: prints that column's dual, reduced cost, state and candidate position at ",
            "every pivot. Read only by the legacy engine, so setting it also selects that engine ",
            "(`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ZCHECK",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy dual simplex: checks the dual objective z against a recomputed identity and ",
            "reports drift. Read only by the legacy engine, so setting it also selects that engine ",
            "(`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ZFLOOR",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a multiplier K (see ZGAIN_FLOOR_K)",
        doc: concat!(
            "Lever: the z-gain noise-floor multiplier of the plateau ladder; `0` disables the ",
            "floor."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ZLAZY",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = per-visit pass, anything else = lazy",
        doc: concat!(
            "Lever: on large LPs the plateau ladder skips the z pass while the stall count is ",
            "below half the Bland threshold."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ZNET",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a lag in visits",
        doc: concat!(
            "Lever: the ladder's z test compares against the z of the last gain this many visits ",
            "back, so one dip does not re-arm the plateau window."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ZTHRESH",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a tolerance (default 1e-6)",
        doc: concat!(
            "How far the recomputed dual objective may drift from the maintained one before zcheck ",
            "reports it."
        ),
    },
    EnvVar {
        name: "RUSOLVE_DUAL_ZTRACE",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy dual simplex: prints each ladder event's delta-z alongside the state identity. ",
            "Read only by the legacy engine, so setting it also selects that engine ",
            "(`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_EKK",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::Any,
        values: "0 = legacy simplex, any other value = ekk",
        doc: concat!(
            "Engine choice for `rusolve solve`: forces the HiGHS-port ekk engine, or with `=0` the ",
            "legacy simplex; it beats the legacy-lever census rule and loses only to `--ipm`."
        ),
    },
    EnvVar {
        name: "RUSOLVE_EKK_HIPO",
        read_by: &["crates/rusolve-cli/src/ekk.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::Ekk,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Makes the HiPO interior point code available to the ekk LP path (`hipo_available`), so ",
            "`solver=ipm` runs HiPO then IPX as HiGHS 1.15+ does; off routes `solver=ipm` to IPX."
        ),
    },
    EnvVar {
        name: "RUSOLVE_HIPO_PROF",
        read_by: &["crates/rusolve-ekk/src/ipm/hipo/factorhighs/call_and_time_blas.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "unset = off, anything else = on",
        doc: "Times every HiPO BLAS call and the dense-factor phases (summed thread seconds and call counts) and prints the table in the HiPO summary.",
    },
    EnvVar {
        name: "RUSOLVE_PRESOLVE_PROF",
        read_by: &["crates/rusolve-ekk/src/presolve/hpresolve/presolve_prof.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Any,
        values: "unset = off, anything else = on",
        doc: "Prints one `presolve_prof` line on stderr per presolve: wall seconds and calls of each phase of the HiGHS-port presolve main loop.",
    },
    EnvVar {
        name: "RUSOLVE_EKK_LOG",
        read_by: &["crates/rusolve-cli/src/ekk.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "0 = off, anything else = on",
        doc: "Keeps the HiGHS log of the ekk LP path on the console (its `output_flag`).",
    },
    EnvVar {
        name: "RUSOLVE_EKK_OPTS",
        read_by: &["crates/rusolve-cli/src/ekk.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::Ekk,
        values: "comma-separated `key=value` HiGHS options",
        doc: concat!(
            "Sets HiGHS options by name on the ekk LP path before the solve; an unknown option or ",
            "invalid value is an error."
        ),
    },
    EnvVar {
        name: "RUSOLVE_EKK_PRESOLVE",
        read_by: &["crates/rusolve-cli/src/ekk.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::Ekk,
        values: "a HiGHS presolve string (`off` disables)",
        doc: concat!(
            "Overrides the ekk path's `presolve` option; `off` keeps the bit-exact presolve-off path, ",
            "unset keeps the HiGHS default (`choose`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_FIXTURES_DIR",
        read_by: &["crates/rusolve-ekk/tests/common/fixture.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::Ekk,
        values: "a directory path",
        doc: concat!(
            "Replaces the whole oracle fixture set: the manifest and every fixture are read from that ",
            "directory only, for running the port against another HiGHS version's oracle output."
        ),
    },
    EnvVar {
        name: "RUSOLVE_FT",
        read_by: &["crates/rusolve-lp/src/lu_ft.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "`1` or `on` = enabled",
        doc: concat!(
            "Selects the Forrest-Tomlin updated factorization for the legacy simplex instead of the ",
            "product-form eta file."
        ),
    },
    EnvVar {
        name: "RUSOLVE_GIT",
        read_by: &[
            "crates/rusolve-cli/src/envelope.rs",
            "crates/rusolve-cli/tests/envelope.rs",
        ],
        kind: EnvKind::Harness,
        scope: EnvScope::Any,
        values: "a 12-hex commit, optionally with -dirty",
        doc: concat!(
            "Compile-time only: build.rs sets it from `git rev-parse` and the envelope reports it ",
            "as `solver.git`. Setting it at run time changes nothing."
        ),
    },
    EnvVar {
        name: "RUSOLVE_HIGHS_FIXTURES",
        read_by: &["crates/rusolve-ekk/tests/common/fixture.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::Ekk,
        values: "a directory path",
        doc: concat!(
            "Extra directory searched first for target-only oracle fixtures (generated by ",
            "tools/highs-oracle), before `target/highs-fixtures`."
        ),
    },
    EnvVar {
        name: "RUSOLVE_HIPO_TRACE",
        read_by: &["crates/rusolve-ekk/src/ipm/hipo/oracle_trace.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "a file path",
        doc: concat!(
            "Port HiPO interior-point code: truncates and appends its oracle-trace lines to that ",
            "file, for the C++/Rust trace diff."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Routes `rusolve solve` to the interior-point LP instead of the simplex (the same ",
            "thing `--ipm` does)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_ADAPT",
        read_by: &["crates/rusolve-lp/src/ipm/solver.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a factor, 0 = off (default 1)",
        doc: concat!(
            "Raises rho and delta when a factorization hits Theta-cancellation pivots in the legacy ",
            "IPM."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_DELTA",
        read_by: &["crates/rusolve-lp/src/ipm/solver.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a number (default DELTA_DEFAULT)",
        doc: "Dual regularization delta (M += delta*I) of the legacy IPM.",
    },
    EnvVar {
        name: "RUSOLVE_IPM_DENSE",
        read_by: &["crates/rusolve-lp/src/ipm/normal.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = never, 1 = cost model (default), 2 = always",
        doc: "Dense-column splitting policy of the legacy IPM's normal-equations factorization.",
    },
    EnvVar {
        name: "RUSOLVE_IPM_DENSE_MODREL",
        read_by: &["crates/rusolve-lp/src/ipm/normal.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a relative threshold (default 1e-8)",
        doc: concat!(
            "A sparse-part pivot at or below this fraction of M_jj is modified rather than skipped ",
            "(Andersen 1996)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_DEPREL",
        read_by: &["crates/rusolve-lp/src/ipm/solver.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a relative threshold (default DEP_REL, 0 disables)",
        doc: concat!(
            "Tiny-pivot threshold of the structural-dependence probe on the legacy IPM's starting ",
            "factor."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_LOG",
        read_by: &[
            "crates/rusolve-lp/src/ipm/solver.rs",
            "crates/rusolve-lp/src/ipm/normal.rs",
        ],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy IPM: prints per-iteration and factorization diagnostics to stderr. Read only by ",
            "the legacy engine, so setting it also selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_MAXIT",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an iteration count (default 200)",
        doc: "Maximum IPM iterations for `rusolve solve --ipm`.",
    },
    EnvVar {
        name: "RUSOLVE_IPM_PIVPERT",
        read_by: &["crates/rusolve-lp/src/ipm/chol.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a relative perturbation (default 0 = off)",
        doc: concat!(
            "When above 0, tiny pivots in the IPM's normal-equations factorization are perturbed to ",
            "this fraction of M_jj instead of skipped."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_PIVREL",
        read_by: &["crates/rusolve-lp/src/ipm/chol.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a relative threshold (default 1e-14)",
        doc: concat!(
            "Relative tiny-pivot threshold: a pivot at or below this fraction of M_jj counts as ",
            "cancellation noise (Wright 1999)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_REG_DECAY",
        read_by: &["crates/rusolve-lp/src/ipm/solver.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a decay factor (default 1 = keep the raise)",
        doc: "How fast a regularization raise on rho/delta decays back toward its starting value.",
    },
    EnvVar {
        name: "RUSOLVE_IPM_RELAX",
        read_by: &["crates/rusolve-lp/src/ipm/chol.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an artificial-flop ratio (default 0.02, 0 = off)",
        doc: concat!(
            "Upper target for the artificial-flop ratio of relaxed amalgamation in the IPM's symbolic ",
            "factorization (HiPO kUpperRatioRelax)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_RETRY_REG",
        read_by: &["crates/rusolve-lp/src/ipm/solver.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a regularization value (<= 0 disables)",
        doc: concat!(
            "When the first IPM run ends NoProgress or Imprecise, retry once with rho = delta = this ",
            "value."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPM_RHO",
        read_by: &["crates/rusolve-lp/src/ipm/solver.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a number (default RHO_DEFAULT)",
        doc: "Primal regularization rho of the legacy IPM (Theta inverse += rho).",
    },
    EnvVar {
        name: "RUSOLVE_IPM_THREADS",
        read_by: &["crates/rusolve-lp/src/ipm/chol.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a thread count, 0 = all cores (default 1)",
        doc: concat!(
            "Factorization threads of the legacy IPM; the default keeps the run single-threaded like ",
            "every other rusolve path."
        ),
    },
    EnvVar {
        name: "RUSOLVE_IPX_TRACE",
        read_by: &["crates/rusolve-ekk/src/ipm/ipx/oracle_trace.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "a file path",
        doc: concat!(
            "Port IPX interior-point code: truncates and appends its oracle-trace lines to that file, ",
            "for the C++/Rust trace diff."
        ),
    },
    EnvVar {
        name: "RUSOLVE_LOCALMIP_LOG",
        read_by: &["crates/rusolve-mip/src/bnb.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::NativeMip,
        values: "set = on",
        doc: "Native MIP: prints the local-MIP sub-solve's progress lines to stderr.",
    },
    EnvVar {
        name: "RUSOLVE_LP_LOG",
        read_by: &["crates/rusolve-lp/src/simplex.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy engine: prints the presolve reduction (rows, columns and nonzeros before and ",
            "after) for the LP-gap decomposition. Read only by the legacy engine, so setting it also ",
            "selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_LS_STATS",
        read_by: &["crates/rusolve-ekk/src/mip/fpr_local_mip/stats.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "set to exactly 1",
        doc: concat!(
            "Port MIP: prints the FJ, local-MIP and FPR heuristic counters (calls, commits, effort ",
            "per second) at the end of the solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_LU_HSKIP",
        read_by: &["crates/rusolve-lp/src/lu.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a mode integer 0-3 (default 3)",
        doc: concat!(
            "Selects among the Markowitz high-skip search modes of the sparse LU factorization; the ",
            "modes are bit-identical, so only the cost moves."
        ),
    },
    EnvVar {
        name: "RUSOLVE_MIP_AC_FIXTURE",
        read_by: &["crates/rusolve-ekk/src/mip/ipm_replay.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::Ekk,
        values: "a path to a `.ac` fixture",
        doc: concat!(
            "Port MIP: replays the oracle's analytic-centre/IPM records from that file instead of ",
            "computing them, so a trace diff isolates divergence before the call."
        ),
    },
    EnvVar {
        name: "RUSOLVE_MIP_EKK",
        read_by: &["crates/rusolve-mip/src/engine_select.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::NativeMip,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Native MIP: uses the ported HiGHS simplex (ekk) as the node LP engine instead of the ",
            "rusolve engine."
        ),
    },
    EnvVar {
        name: "RUSOLVE_MIP_HIGHS",
        read_by: &["crates/rusolve-cli/src/mip_highs.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::PortMip,
        values: "0 = native search, anything else or unset = HiGHS port",
        doc: concat!(
            "`rusolve mip` runs the rusolve-ekk HiGHS MIP port by default (the envelope, `--opt` and the documented ",
            "exit codes exist only on that path); `0` selects the native rusolve-mip search."
        ),
    },
    EnvVar {
        name: "RUSOLVE_STD_LIBSTDCXX",
        read_by: &["crates/rusolve-ekk/src/util/cpp_std.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::PortMip,
        values: "1 = libstdc++ tie order, anything else or unset = libc++",
        doc: concat!(
            "Port MIP: the C++ heap operations (priority_queue, make/pop/sort_heap, partial_sort) break ties ",
            "as libstdc++ does (the Linux highspy wheel) instead of libc++ (the macOS oracle the identity ",
            "references come from). Only heap operations switch (not HighsMultiBuckets order, the pdqsort heap fallback, ",
            "IPX RandomPermute, or FMA contraction), so it is not a general bit-exact oracle. Diagnostic, default off."
        ),
    },
    EnvVar {
        name: "RUSOLVE_MIP_LOG",
        read_by: &[
            "crates/rusolve-mip/src/bnb.rs",
            "crates/rusolve-mip/src/presolve.rs",
            "crates/rusolve-mip/src/branching.rs",
        ],
        kind: EnvKind::Debug,
        scope: EnvScope::NativeMip,
        values: "set = on",
        doc: concat!(
            "Native MIP: prints presolve, cut-round, restart, clique and branching-tree progress ",
            "lines to stderr."
        ),
    },
    EnvVar {
        name: "RUSOLVE_NO_P1REPERT",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: concat!(
            "Turns off phase-1 re-perturbation, which by default re-perturbs a frozen epsilon ",
            "across a stall window."
        ),
    },
    EnvVar {
        name: "RUSOLVE_NO_P2PERT",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: concat!(
            "Turns off phase-2 entering-tie perturbation, which by default engages after 4x the ",
            "Bland stall threshold."
        ),
    },
    EnvVar {
        name: "RUSOLVE_NO_RELAXED",
        read_by: &["crates/rusolve-lp/src/simplex.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (disables the escape)",
        doc: concat!(
            "Suppresses the phase-1 escape into relaxed row choice after a long frozen-and-stagnant ",
            "stall."
        ),
    },
    EnvVar {
        name: "RUSOLVE_NO_SH",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: concat!(
            "Turns off SH (strong-Harris) pricing, falling back to Dantzig pricing in the legacy ",
            "simplex."
        ),
    },
    EnvVar {
        name: "RUSOLVE_P1PROBE",
        read_by: &["crates/rusolve-lp/src/simplex.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy phase 1: every 20,000 stalled pivots at a frozen vertex, probes whether any ",
            "improver column has a strictly positive ratio. Read only by the legacy engine, so ",
            "setting it also selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_P1TRACE",
        read_by: &["crates/rusolve-lp/src/simplex.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy phase 1: collects and prints the phase-1 trace (w probe and stall census) for the ",
            "artificial phase. Read only by the legacy engine, so setting it also selects that engine ",
            "(`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_P1_ROUTE_DUAL",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Policy arm: with RUSOLVE_PRIMAL=1, cold primal-infeasible starts go to the dual ",
            "engine instead of paying an artificial-basis phase 1."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PARCOLDOM_FUZZ_DUMP",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a seed number",
        doc: concat!(
            "Seed whose generated parallel-column LP the dominance fuzz test dumps to stderr before ",
            "checking it. rusolve-lp unit test only; never read by a shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PARCOLDOM_FUZZ_SEEDS",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a seed count (default 400)",
        doc: concat!(
            "How many random parallel-column LPs the column-dominance differential fuzz test checks ",
            "against presolve off. rusolve-lp unit test only; never read by a shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PDLP_TRACE",
        read_by: &["crates/rusolve-ekk/src/pdlp/oracle_trace.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "a file path",
        doc: concat!(
            "Port PDLP code: truncates and appends its oracle-trace lines to that file, for the ",
            "C++/Rust trace diff."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PFI_DEBUG",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy engine: prints what the factorization did with a singular sparse LU, including ",
            "the no-dense-fallback case. Read only by the legacy engine, so setting it also selects ",
            "that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PFI_DENSE_MAX",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a basis dimension (default 4096)",
        doc: concat!(
            "Above this basis size a singular sparse LU is no longer retried with the dense engine, ",
            "and the singularity propagates."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PRESOLVE_AGG",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "0 = off, anything else = on",
        doc: concat!(
            "Legacy presolve's general-scale duplicate-column aggregator; `=0` is the lever ",
            "back."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PRESOLVE_HASH",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy presolve: prints an FNV-1a hash over the exact bits of the reduced LP, so two ",
            "builds can be proven to apply identical reductions. Read only by the legacy engine, so ",
            "setting it also selects that engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PRESOLVE_STATS",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy presolve: prints per-rule effectiveness counters (what each reduction removed). ",
            "Read only by the legacy engine, so setting it also selects that engine ",
            "(`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PRIMAL",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: concat!(
            "Runs `rusolve solve` with the primal simplex; the default engine is the dual simplex ",
            "by measurement on our benchmark set."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PROFILE",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Profile,
        scope: EnvScope::Any,
        values: "`recommended` (default) or `highs-1.15.1`",
        doc: concat!(
            "Base HiGHS option set every CLI solve starts from; `highs-1.15.1` turns every H16 option ",
            "off, and any other value is a usage error."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_IMPEQ",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: concat!(
            "Enables the implied-equality conversion in the legacy presolve, which measured ",
            "harmful and stayed behind the lever."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_MEMGUARD_BYTES",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a byte budget (0 disables the guard)",
        doc: concat!(
            "Overrides the presolve memory-guard budget above which the reduction loop is skipped for ",
            "the lpopt giant tier."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_FOLD",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "A/B lever: disables the legacy presolve's bound-fold mutation only, for per-instance ",
            "regression hunting."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_FORCECOL",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "A/B lever: disables the forcing-column mutation only; the detection probe keeps ",
            "counting under RUSOLVE_PRESOLVE_STATS."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_MEMGUARD",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Disables the presolve memory guard so giant LPs are reduced regardless of the fitted ",
            "working-set estimate."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_PARCOL",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: "Lever back for the general-scale duplicate-column merge in the legacy presolve.",
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_PARCOLDOM",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: "Lever back for the parallel-column dominance arms of the legacy presolve.",
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_PARROW",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: "A/B lever: disables the parallel-row merge mutation only.",
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_PARROWEQ",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "A/B lever: disables only the cross-sense interval-collapse arm of the parallel-row rule ",
            "(the same-sense merge stays under RUSOLVE_PS_NO_PARROW)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_PS_NO_VUBEQ",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (lever back)",
        doc: concat!(
            "Lever back for the singleton-equality (variable-bound equality) rule in the legacy ",
            "presolve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_QP_TRACE",
        read_by: &["crates/rusolve-ekk/src/qpsolver/oracle_trace.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "a file path",
        doc: concat!(
            "Port QP solver code: truncates and appends its oracle-trace lines to that file, for the ",
            "C++/Rust trace diff."
        ),
    },
    EnvVar {
        name: "RUSOLVE_REFINE_AGE",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "a fraction in (0, 1]",
        doc: concat!(
            "Skips the iterative refinement of a PFI solve whose basis is younger than this fraction ",
            "of the refactor age; path-visible, so the iteration gate adjudicates it."
        ),
    },
    EnvVar {
        name: "RUSOLVE_REFINE_GATE",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set to exactly 1",
        doc: concat!(
            "Makes the PFI round-0 correction conditional on the double-word residual clearing the fp ",
            "floor, saving a raw solve when it already sits there."
        ),
    },
    EnvVar {
        name: "RUSOLVE_REQUIRE_TARGET_FIXTURES",
        read_by: &[
            "crates/rusolve-ekk/tests/common/fixture.rs",
            "crates/rusolve-ekk/tests/port/w1.rs",
            "crates/rusolve-cli/tests/mip_trace_threads.rs",
            "crates/rusolve-cli/tests/conflict_resolution.rs",
            "crates/rusolve-cli/tests/cut_selection.rs",
            "crates/rusolve-cli/tests/fpr_local_mip.rs",
            "crates/rusolve-cli/tests/scip_heuristics.rs",
        ],
        kind: EnvKind::Harness,
        scope: EnvScope::Any,
        values: "set to exactly 1",
        doc: concat!(
            "Turns a missing target-only oracle fixture from a SKIP banner into a test failure, for ",
            "CI gates."
        ),
    },
    EnvVar {
        name: "RUSOLVE_RFCAD",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "`K[:M]` (K >= 1, M a minimum-update floor)",
        doc: concat!(
            "Retunes the refactor clock of every factorization: K raises the clock threshold to K * ",
            "build_work, M floors the update count."
        ),
    },
    EnvVar {
        name: "RUSOLVE_RFDROP",
        read_by: &["crates/rusolve-lp/src/pfi.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Enables the per-refactor drift probe, which needs RUSOLVE_DUAL_STATS=1 and exists for ",
            "measurement builds only. Read only by the legacy engine, so setting it also selects that ",
            "engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_S2_STATS",
        read_by: &["crates/rusolve-ekk/src/mip/conflict_resolution/output.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "set to exactly 1",
        doc: concat!(
            "Port MIP conflict resolution: prints the `s2-stats` counters line when the conflict separator ",
            "adds rows."
        ),
    },
    EnvVar {
        name: "RUSOLVE_S2_TRACE",
        read_by: &["crates/rusolve-ekk/src/mip/conflict_resolution/mod.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "set to exactly 1",
        doc: concat!(
            "Port MIP conflict resolution: prints the first invalid conflict row as a diagnostic, ",
            "without changing the solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_S2_TRACE_DUMP",
        read_by: &["crates/rusolve-ekk/src/mip/conflict_resolution/output.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "set to exactly 1",
        doc: "With RUSOLVE_S2_TRACE=1, also dumps the terms of the first invalid conflict row.",
    },
    EnvVar {
        name: "RUSOLVE_S3_DEBUG_SOL",
        read_by: &["crates/rusolve-ekk/src/mip/symmetry_dynamic/debugsol.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::Ekk,
        values: "a path to a solution file",
        doc: concat!(
            "Port MIP dynamic symmetry: reads a known-optimal solution and checks the dynamic ",
            "symmetry-breaking decisions against it."
        ),
    },
    EnvVar {
        name: "RUSOLVE_SCALE",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on (default off)",
        doc: concat!(
            "Enables power-of-2 equilibration scaling in the legacy simplex, which helped some ",
            "instances and stalled others."
        ),
    },
    EnvVar {
        name: "RUSOLVE_SCALE_DEBUG",
        read_by: &["crates/rusolve-lp/src/scaling.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy scaling: prints whether the improvement gate accepted the scaling and by how ",
            "much. Read only by the legacy engine, so setting it also selects that engine ",
            "(`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_SPARSE",
        read_by: &["crates/rusolve-cli/src/main.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Lever: runs the direct CSC-only `Standard` build for `rusolve solve` when presolve ",
            "and dual are off."
        ),
    },
    EnvVar {
        name: "RUSOLVE_STATS_EVERY",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an iteration period",
        doc: concat!(
            "With RUSOLVE_DUAL_STATS=1, dumps the cost counters every N iterations so a time-limited ",
            "run still reports them. Read only by the legacy engine, so setting it also selects that ",
            "engine (`use_ekk_engine`)."
        ),
    },
    EnvVar {
        name: "RUSOLVE_SUBMIP_LOG",
        read_by: &["crates/rusolve-mip/src/bnb.rs"],
        kind: EnvKind::Debug,
        scope: EnvScope::NativeMip,
        values: "set = on",
        doc: "Native MIP: prints the root sub-MIP's progress lines to stderr.",
    },
    EnvVar {
        name: "RUSOLVE_TEST_PROFILE",
        read_by: &["crates/rusolve-ekk/tests/port/mod.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::Ekk,
        values: "`default`, `recommended` or `all-h16`",
        doc: concat!(
            "Base HighsOptions the port's replay test builders start from, before the oracle records ",
            "override what they carry."
        ),
    },
    EnvVar {
        name: "RUSOLVE_VUBEQ_FUZZ_DUMP",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a seed number",
        doc: concat!(
            "Seed whose generated LP the singleton-equality fuzz test dumps to stderr before checking ",
            "it. rusolve-lp unit test only; never read by a shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_VUBEQ_FUZZ_SEEDS",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "a seed count (default 400)",
        doc: concat!(
            "How many random LPs with an appended singleton-cost row the singleton-equality fuzz test ",
            "checks against presolve off. rusolve-lp unit test only; never read by a shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_VUBEQ_FUZZ_VERBOSE",
        read_by: &["crates/rusolve-lp/src/presolve.rs"],
        kind: EnvKind::Harness,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Prints each seed as the singleton-equality fuzz test runs it, so a failure names its ",
            "seed. rusolve-lp unit test only; never read by a shipped solve."
        ),
    },
    EnvVar {
        name: "RUSOLVE_WARM_INCD_MIN_NNZ",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "an nnz count (default 10000)",
        doc: concat!(
            "Nonzero count from which the incremental reduced costs also engage on warm (MIP node) ",
            "re-solves."
        ),
    },
    EnvVar {
        name: "RUSOLVE_WARM_STATS",
        read_by: &["crates/rusolve-lp/src/dual.rs"],
        kind: EnvKind::Lever,
        scope: EnvScope::LegacyLp,
        values: "set = on",
        doc: concat!(
            "Legacy dual simplex: accumulates the warm-start census (calls, warm vs fallback, timing) ",
            "and prints it at the end. Read only by the legacy engine, so setting it also selects ",
            "that engine (`use_ekk_engine`)."
        ),
    },
];

/// The registry entry for `name`, if the solver honours that name.
pub fn find(name: &str) -> Option<&'static EnvVar> {
    ENV_VARS.iter().find(|v| v.name == name)
}

/// `true` when `name` is in the registry.
pub fn is_known(name: &str) -> bool {
    find(name).is_some()
}

/// `NAME  kind  scope  values  doc`, one variable per line, columns aligned.
pub fn table_text() -> String {
    let w = ENV_VARS.iter().map(|v| v.name.len()).max().unwrap_or(0);
    let mut out = String::new();
    for v in ENV_VARS {
        out.push_str(v.name);
        for _ in v.name.len()..w {
            out.push(' ');
        }
        out.push_str(&format!(
            "  {}  {}  {}  {}\n",
            v.kind.label(),
            v.scope.label(),
            v.values,
            v.doc
        ));
    }
    out
}

/// The whole registry as a JSON array of objects with the registry fields.
pub fn table_json() -> String {
    serde_json::to_string_pretty(ENV_VARS).unwrap_or_else(|e| panic!("env registry json: {e}"))
}

/// Levenshtein distance over chars (the registry names are ASCII, but the
/// caller's typo need not be). Private on purpose: the option-name suggester
/// is a separate concern.
fn levenshtein(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0usize; b.len() + 1];
    for (i, ca) in a.iter().enumerate() {
        cur[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let cost = usize::from(ca != cb);
            cur[j + 1] = (prev[j] + cost).min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

/// The closest registry name to `name`, when one is close enough to be worth
/// naming (distance at most [`MAX_SUGGEST_DISTANCE`]). Ties break on the
/// shorter name, then lexicographically, so the suggestion is deterministic.
pub fn suggest(name: &str) -> Option<&'static EnvVar> {
    let mut best: Option<(usize, usize, &'static EnvVar)> = None;
    for v in ENV_VARS {
        let d = levenshtein(name, v.name);
        if d > MAX_SUGGEST_DISTANCE {
            continue;
        }
        let better = match best {
            None => true,
            Some((bd, blen, bv)) => {
                d < bd
                    || (d == bd
                        && (v.name.len() < blen || (v.name.len() == blen && v.name < bv.name)))
            }
        };
        if better {
            best = Some((d, v.name.len(), v));
        }
    }
    best.map(|(_, _, v)| v)
}

/// The most distant name still worth a "did you mean".
const MAX_SUGGEST_DISTANCE: usize = 3;

/// The startup check: a `RUSOLVE_*` variable in the process environment
/// that no crate reads is a typo, not a lever. Returns the message and, when a
/// registered name is close, the "did you mean" clause (`main` prints them as
/// an RS0005 diagnostic and exits 2, the usage code of the exit-code table in
/// `mip_highs`). Names are checked in sorted order, so the result does not
/// depend on environ order.
pub fn unknown_env() -> Option<(String, Option<String>)> {
    let mut unknown: Vec<String> = std::env::vars_os()
        .filter_map(|(k, _)| k.into_string().ok())
        .filter(|k| k.starts_with("RUSOLVE_") && !is_known(k))
        .collect();
    if unknown.is_empty() {
        return None;
    }
    unknown.sort();
    let name = &unknown[0];
    Some((
        format!("rusolve: unknown environment variable {name} (rusolve env lists them all)"),
        suggest(name).map(|s| format!("did you mean {}?", s.name)),
    ))
}

/// The first set registry variable that has to move `rusolve solve` to the
/// legacy simplex, in sorted name order so the notice does not depend on
/// environ order.
pub fn legacy_engine_lever() -> Option<&'static EnvVar> {
    let mut hits: Vec<&'static EnvVar> = std::env::vars_os()
        .filter_map(|(k, _)| k.into_string().ok())
        .filter_map(|k| find(&k))
        .filter(|v| v.is_legacy_engine_lever())
        .collect();
    hits.sort_by_key(|v| v.name);
    hits.first().copied()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};

    /// `use_ekk_engine` used to switch `rusolve solve` to the legacy simplex
    /// for every `RUSOLVE_*` name outside this list. It is the reference for
    /// [`env_registry_engine_choice_unchanged`], kept here and nowhere else:
    /// the registry replaces it.
    /// Written as one whitespace-separated string rather than a `&[&str]` so
    /// this file's own scan does not read the copied names as reads of them.
    const OLD_EKK_NEUTRAL_ENV: &str = "RUSOLVE_EKK RUSOLVE_EKK_LOG RUSOLVE_EKK_PRESOLVE \
        RUSOLVE_EKK_OPTS RUSOLVE_PROFILE RUSOLVE_BIN RUSOLVE_BENCH_MPS \
        RUSOLVE_HIGHS_FIXTURES RUSOLVE_REQUIRE_TARGET_FIXTURES RUSOLVE_MIP_LOG \
        RUSOLVE_LOCALMIP_LOG RUSOLVE_SUBMIP_LOG RUSOLVE_CUT_DEBUG RUSOLVE_AUDIT RUSOLVE_BRANCH";

    /// [`OLD_EKK_NEUTRAL_ENV`] as a name list.
    fn old_ekk_neutral_env() -> Vec<&'static str> {
        OLD_EKK_NEUTRAL_ENV.split_whitespace().collect()
    }

    /// `crates/`, the tree the registry has to cover.
    fn crates_root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("rusolve-cli is a workspace member")
            .to_path_buf()
    }

    /// `src` with line and block comments removed, so a name mentioned in prose
    /// cannot pass for a read. String literals are kept.
    fn strip_comments(src: &str) -> String {
        #[derive(PartialEq)]
        enum St {
            Code,
            Str,
            Line,
            Block,
        }
        let mut st = St::Code;
        let mut out = String::with_capacity(src.len());
        let mut esc = false;
        let mut chars = src.chars().peekable();
        while let Some(c) = chars.next() {
            match st {
                St::Code => {
                    if c == '"' {
                        st = St::Str;
                        out.push(c);
                    } else if c == '/' && chars.peek() == Some(&'/') {
                        st = St::Line;
                        chars.next();
                    } else if c == '/' && chars.peek() == Some(&'*') {
                        st = St::Block;
                        chars.next();
                    } else {
                        out.push(c);
                    }
                }
                St::Str => {
                    out.push(c);
                    if esc {
                        esc = false;
                    } else if c == '\\' {
                        esc = true;
                    } else if c == '"' {
                        st = St::Code;
                    }
                }
                St::Line => {
                    if c == '\n' {
                        st = St::Code;
                        out.push(c);
                    }
                }
                St::Block => {
                    if c == '*' && chars.peek() == Some(&'/') {
                        st = St::Code;
                        chars.next();
                    } else if c == '\n' {
                        out.push(c);
                    }
                }
            }
        }
        out
    }

    /// Every `.rs` file under `crates/` with its comments stripped, as
    /// (repo-relative path, text). `target/` and vendored trees are not our code.
    fn sources() -> Vec<(String, String)> {
        fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
            let Ok(rd) = std::fs::read_dir(dir) else {
                return;
            };
            for e in rd.flatten() {
                let p = e.path();
                let name = e.file_name().to_string_lossy().into_owned();
                if p.is_dir() {
                    if !matches!(name.as_str(), "target" | "rea" | ".git" | "node_modules") {
                        walk(&p, out);
                    }
                } else if name.ends_with(".rs") {
                    out.push(p);
                }
            }
        }
        let root = crates_root();
        let mut files = Vec::new();
        walk(&root, &mut files);
        files.sort();
        files
            .iter()
            .map(|p| {
                let text =
                    std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
                let rel = p
                    .strip_prefix(
                        root.parent()
                            .expect("crates/ sits under the workspace root"),
                    )
                    .expect("file under crates/")
                    .to_string_lossy()
                    .into_owned();
                (rel, strip_comments(&text))
            })
            .collect()
    }

    /// Every `"NAME"` string literal whose NAME is `RUSOLVE_…`.
    fn quoted_env_names(text: &str) -> Vec<String> {
        let mut out = Vec::new();
        let b = text.as_bytes();
        let mut i = 0;
        while i < b.len() {
            if b[i] != b'"' {
                i += 1;
                continue;
            }
            let start = i + 1;
            let mut j = start;
            while j < b.len() && b[j] != b'"' && b[j] != b'\n' {
                j += 1;
            }
            if j < b.len() && b[j] == b'"' {
                let s = &text[start..j];
                if s.starts_with("RUSOLVE_")
                    && s.len() > "RUSOLVE_".len()
                    && s.bytes()
                        .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == b'_')
                {
                    out.push(s.to_string());
                }
                i = j + 1;
            } else {
                i += 1;
            }
        }
        out
    }

    /// Every `const NAME: &str = "VALUE"` in `text`.
    fn const_strings(text: &str) -> Vec<(String, String)> {
        let mut out = Vec::new();
        for line in text.lines() {
            let t = line.trim_start();
            let t = ["pub(crate) ", "pub(super) ", "pub "]
                .iter()
                .find_map(|p| t.strip_prefix(p))
                .unwrap_or(t);
            let Some(rest) = t.strip_prefix("const ") else {
                continue;
            };
            let ident: String = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                .collect();
            if ident.is_empty() {
                continue;
            }
            let Some(eq) = rest.find('=') else { continue };
            let after = rest[eq + 1..].trim_start();
            let Some(after) = after.strip_prefix('"') else {
                continue;
            };
            let Some(end) = after.find('"') else { continue };
            out.push((ident, after[..end].to_string()));
        }
        out
    }

    /// The keys of every `env::var(...)` / `env::var_os(...)` call in `text`,
    /// as (string-literal keys, identifier keys). The literals are what catch
    /// the handful of variables without the `RUSOLVE_` prefix (`PD_DETAIL`,
    /// `P42_REPS`); the identifiers are the `&str` constants a key is named
    /// through (`FIXTURE_DIR_ENV`, `kPdlpTraceEnv`), and lowercase ones are
    /// function parameters (e.g. `ipm/solver.rs`'s `env_f64(key, ..)`), whose
    /// call sites all pass a literal.
    ///
    /// The call name is spelled in two pieces on purpose: written whole, this
    /// file's own source would read as an `env::var` call whose key is the
    /// literal that follows it.
    fn env_var_keys(text: &str) -> (Vec<String>, Vec<String>) {
        const CALL: &str = concat!("env::va", "r");
        let mut literals = Vec::new();
        let mut idents = Vec::new();
        let mut from = 0;
        while let Some(k) = text[from..].find(CALL) {
            let after_name = from + k + CALL.len();
            from = after_name;
            let Some(paren) = text[after_name..].find('(') else {
                break;
            };
            let arg = text[after_name + paren + 1..].trim_start();
            if let Some(rest) = arg.strip_prefix('"') {
                if let Some(end) = rest.find('"') {
                    literals.push(rest[..end].to_string());
                }
            } else {
                let ident: String = arg
                    .chars()
                    .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                    .collect();
                if !ident.is_empty() {
                    idents.push(ident);
                }
            }
        }
        (literals, idents)
    }

    /// The registry's own sources, exempt from the literal scan: its tests name
    /// variables that must NOT be registered (`RUSOLVE_PROFLE`, the typo the CLI
    /// has to reject). The exemption is closed by the assertion below that these files
    /// never name an environment variable in an `env::var` call.
    const SELF_TESTS: &[&str] = &[
        "crates/rusolve-cli/src/env_registry.rs",
        "crates/rusolve-cli/tests/env_registry.rs",
    ];

    /// Every variable the sources read is in the registry, every registry entry
    /// is still read somewhere, and every `read_by` file is real.
    #[test]
    fn env_registry_covers_every_read() {
        let srcs = sources();
        assert!(
            srcs.len() > 100,
            "expected to walk crates/, got {}",
            srcs.len()
        );
        let mut consts: Vec<(String, String, String)> = Vec::new();
        for (path, text) in &srcs {
            for (k, v) in const_strings(text) {
                consts.push((k, v, path.clone()));
            }
        }
        let resolve = |ident: &str| -> Option<String> {
            consts
                .iter()
                .find(|(k, _, _)| k == ident)
                .map(|(_, v, _)| v.clone())
        };

        // (variable, where it is read)
        let mut read: Vec<(String, String)> = Vec::new();
        let mut unresolved: Vec<String> = Vec::new();
        for (path, text) in &srcs {
            let (literals, key_idents) = env_var_keys(text);
            // Lowercase identifiers are function parameters, not keys, and the
            // scanner can land on one when a source mentions the call inside a
            // string. Only constant-shaped names are evidence of a read.
            let const_keys: Vec<&String> = key_idents
                .iter()
                .filter(|i| i.chars().any(|c| c.is_ascii_uppercase()))
                .collect();
            if SELF_TESTS.contains(&path.as_str()) {
                assert!(
                    literals.is_empty() && const_keys.is_empty(),
                    "{path} is exempt from the literal scan, so it may not read an \
                     environment variable by name: {literals:?} {const_keys:?}"
                );
                continue;
            }
            for n in quoted_env_names(text) {
                read.push((n, path.clone()));
            }
            for n in literals {
                read.push((n, path.clone()));
            }
            for ident in const_keys {
                match resolve(ident) {
                    Some(v) => read.push((v, format!("{path} (via const {ident})"))),
                    None => unresolved.push(format!("{path}: env::var({ident})")),
                }
            }
        }
        // A key named through a constant the scan cannot resolve would be a
        // hole in the check; only lowercase function parameters may be left.
        for u in &unresolved {
            let ident = u.rsplit('(').next().unwrap_or("").trim_end_matches(')');
            assert!(
                !ident.chars().any(|c| c.is_ascii_uppercase()),
                "env::var argument {u} is neither a literal nor a known &str constant"
            );
        }

        let unregistered: Vec<String> = read
            .iter()
            .filter(|(n, _)| !is_known(n))
            .map(|(n, w)| format!("{w}: {n}"))
            .collect();
        assert!(
            unregistered.is_empty(),
            "{} read(s) missing from the registry:\n{}",
            unregistered.len(),
            unregistered.join("\n")
        );

        let mut distinct: Vec<&str> = read.iter().map(|(n, _)| n.as_str()).collect();
        distinct.sort_unstable();
        distinct.dedup();
        // Entries only the tools/ scripts read (RUSOLVE_BIN) are registered so their
        // inheriting children do not reject them; the crate scan cannot see them.
        let script_only = |v: &EnvVar| v.read_by.iter().all(|f| !f.starts_with("crates/"));
        let crate_entries = ENV_VARS.iter().filter(|v| !script_only(v)).count();
        assert_eq!(
            distinct.len(),
            crate_entries,
            "the scan found {} distinct variables, the registry lists {} read by crates",
            distinct.len(),
            crate_entries
        );

        let all_text: String = srcs
            .iter()
            .map(|(_, t)| t.as_str())
            .collect::<Vec<_>>()
            .join("\n");
        let dead: Vec<&str> = ENV_VARS
            .iter()
            .filter(|v| !script_only(v))
            .map(|v| v.name)
            .filter(|n| !all_text.contains(&format!("\"{n}\"")))
            .collect();
        assert!(
            dead.is_empty(),
            "registry names no source reads any more: {dead:?}"
        );

        let workspace_root = crates_root().parent().unwrap().to_path_buf();
        for v in ENV_VARS {
            for f in v.read_by {
                assert!(
                    workspace_root.join(f).exists(),
                    "{} lists a reader that does not exist: {f}",
                    v.name
                );
                if !f.starts_with("crates/") {
                    let text = std::fs::read_to_string(workspace_root.join(f)).unwrap();
                    assert!(
                        text.contains(v.name),
                        "{} lists {f}, which does not name it",
                        v.name
                    );
                    continue;
                }
                let reads = srcs.iter().any(|(p, t)| {
                    p == f
                        && (t.contains(v.name)
                            || consts
                                .iter()
                                .any(|(k, value, _)| value == v.name && t.contains(k.as_str())))
                });
                assert!(
                    reads,
                    "{} lists {f} as a reader, but it does not name it",
                    v.name
                );
            }
        }
    }

    /// The registry names are unique, non-empty and each has a reader listed.
    #[test]
    fn env_registry_entries_are_well_formed() {
        let mut seen: Vec<&str> = Vec::new();
        for v in ENV_VARS {
            assert!(!v.name.is_empty(), "entry with an empty name");
            assert!(!v.read_by.is_empty(), "{} lists no reader", v.name);
            assert!(!v.values.is_empty(), "{} lists no values", v.name);
            assert!(!v.doc.is_empty(), "{} lists no doc", v.name);
            assert!(!seen.contains(&v.name), "duplicate entry {}", v.name);
            seen.push(v.name);
        }
        assert!(ENV_VARS.len() > 100, "registry lost entries");
    }

    /// The engine choice is unchanged for every variable that exists on the
    /// legacy path: the fix is that a variable the legacy engine never reads
    /// can no longer switch `rusolve solve` off ekk, and nothing else moved.
    #[test]
    fn env_registry_engine_choice_unchanged() {
        let old_list = old_ekk_neutral_env();
        let old = |name: &str| !old_list.contains(&name);
        let mut changed: Vec<&str> = Vec::new();
        for v in ENV_VARS {
            let new = v.is_legacy_engine_lever();
            if old(v.name) == new {
                continue;
            }
            assert!(
                !new,
                "{}: newly switches engines, which the registry forbids",
                v.name
            );
            let legacy_read = v.read_by.iter().any(|f| f.contains("crates/rusolve-lp"));
            assert!(
                v.kind == EnvKind::Harness || (!legacy_read && v.scope != EnvScope::LegacyLp),
                "{}: {} {} changed the engine choice, but the legacy engine reads it",
                v.name,
                v.kind.label(),
                v.scope.label()
            );
            changed.push(v.name);
        }
        assert!(
            !changed.is_empty(),
            "the inverse allowlist was the bug: nothing changed means the registry is wrong"
        );

        // Every old allowlist name must come out as NOT a legacy-engine lever,
        // and every one is registered (RUSOLVE_BIN as a script-only harness entry).
        assert_eq!(old_list.len(), 15, "the copied allowlist lost a name");
        for name in old_list {
            match find(name) {
                Some(v) => assert!(
                    !v.is_legacy_engine_lever(),
                    "{} was engine-neutral and now switches engines",
                    name
                ),
                None => panic!("{name} vanished from the registry"),
            }
        }

        // The levers census scripts set to study the legacy engine
        // still select it, so those runs measure what they always measured.
        for name in [
            "RUSOLVE_PRIMAL",
            "RUSOLVE_SPARSE",
            "RUSOLVE_SCALE",
            "RUSOLVE_NO_SH",
            "RUSOLVE_NO_P2PERT",
            "RUSOLVE_BLAND_AFTER",
            "RUSOLVE_DUAL_STATS",
            "RUSOLVE_LP_LOG",
            "RUSOLVE_PRESOLVE_STATS",
            "RUSOLVE_STATS_EVERY",
            "RUSOLVE_RFDROP",
        ] {
            assert!(
                find(name).is_some_and(|v| v.is_legacy_engine_lever()),
                "{name} must stay a legacy-engine lever"
            );
        }

        // The variables whose silent switch was the complaint stay neutral.
        for name in [
            "RUSOLVE_PDLP_TRACE",
            "RUSOLVE_QP_TRACE",
            "RUSOLVE_MIP_LOG",
            "RUSOLVE_CUT_DEBUG",
            "RUSOLVE_EKK_LOG",
            "RUSOLVE_EKK_OPTS",
            "RUSOLVE_MIP_HIGHS",
            "RUSOLVE_AUDIT",
            "RUSOLVE_BRANCH",
            "RUSOLVE_PROFILE",
            "RUSOLVE_TEST_PROFILE",
        ] {
            assert!(
                find(name).is_some_and(|v| !v.is_legacy_engine_lever()),
                "{name} must not switch engines"
            );
        }
    }

    /// The suggestion is the nearest name, and only when it is close.
    #[test]
    fn env_registry_suggests_the_nearest_name() {
        assert_eq!(
            suggest("RUSOLVE_PROFLE").map(|v| v.name),
            Some("RUSOLVE_PROFILE")
        );
        assert_eq!(
            suggest("RUSOLVE_DUAL_ACC").map(|v| v.name),
            Some("RUSOLVE_DUAL_ACC")
        );
        assert_eq!(
            suggest("RUSOLVE_PRIMAl").map(|v| v.name),
            Some("RUSOLVE_PRIMAL")
        );
        assert_eq!(suggest("RUSOLVE_NOT_A_REAL_LEVER_AT_ALL"), None);
        assert_eq!(levenshtein("kitten", "sitting"), 3);
        assert_eq!(levenshtein("", "abc"), 3);
        assert_eq!(levenshtein("RUSOLVE_PROFILE", "RUSOLVE_PROFILE"), 0);
    }

    /// The startup check names the offender and the fix, and stays quiet on a
    /// registered name.
    #[test]
    fn env_registry_startup_check_text() {
        let msg = format!(
            "rusolve: unknown environment variable {}; did you mean {}? (rusolve env lists them all)",
            "RUSOLVE_PROFLE",
            suggest("RUSOLVE_PROFLE").expect("close match").name
        );
        assert_eq!(
            msg,
            "rusolve: unknown environment variable RUSOLVE_PROFLE; did you mean RUSOLVE_PROFILE? (rusolve env lists them all)"
        );
        assert!(is_known("RUSOLVE_PROFILE"));
        assert!(!is_known("RUSOLVE_PROFLE"));
    }

    /// `rusolve env` renders every entry, in both shapes.
    #[test]
    fn env_registry_renders_every_entry() {
        let text = table_text();
        for v in ENV_VARS {
            assert!(
                text.contains(v.name),
                "{} missing from the text table",
                v.name
            );
            assert!(text.contains(v.kind.label()), "{} kind missing", v.name);
        }
        let json: serde_json::Value =
            serde_json::from_str(&table_json()).expect("the env table is valid JSON");
        let arr = json.as_array().expect("the env table is a JSON array");
        assert_eq!(arr.len(), ENV_VARS.len());
        for (v, e) in ENV_VARS.iter().zip(arr) {
            assert_eq!(e["name"].as_str(), Some(v.name));
            assert_eq!(
                e["kind"].as_str(),
                Some(v.kind.label().to_lowercase().as_str())
            );
            assert_eq!(e["values"].as_str(), Some(v.values));
            assert_eq!(e["doc"].as_str(), Some(v.doc));
            assert_eq!(
                e["read_by"].as_array().map(|a| a.len()),
                Some(v.read_by.len())
            );
        }
    }
}
