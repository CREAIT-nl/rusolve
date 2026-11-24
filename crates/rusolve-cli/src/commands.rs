//! The command registry: every subcommand of `rusolve` with its
//! one-line summary, its flags (name, value kind, default, description, which solver
//! stack honours it) and the id of the JSON schema its `--format json` output follows.
//!
//! One table, four readers:
//!
//! * [`usage_text`] / [`command_usage`]: what `rusolve --help` and `rusolve <cmd> --help`
//!   print, and the usage a bad invocation falls back to;
//! * [`help_json`] / [`command_help_json_doc`]: `rusolve --help --json`, schema
//!   `rusolve.help/1` (`docs/schema/rusolve-help-1.schema.json`);
//! * [`flag_stack_error`]: the message a stack gives for a flag only the other stack has;
//! * [`parse_flags`]: the parser itself. It looks every flag up here, and checks presence,
//!   kind, choice list and repeatability from [`ValueKind`], so a flag added to the table is
//!   parsed, documented and typed everywhere without a line of parser code. The `mip` port
//!   path (`mip_highs::run_mip_highs`) already parses this way; the other commands keep their
//!   hand-written loops, which the tests below check against this table.
//!
//! The option registry stays the source of `--opt` names
//! (`rusolve_ekk::lp_data::highs_option_records`); the `--opt` and `--options-file` flags
//! point at it (`rusolve options --json`) instead of copying it, and
//! [`option_properties`] reads it for the MCP `solve` tool. The environment variables are
//! `env_registry`.
//!
//! ## Keeping it complete
//!
//! `every_dispatch_arm_has_a_registry_entry_and_vice_versa` scans `main.rs` for the
//! `cmd == "<name>"` dispatch arms and compares them with [`COMMANDS`];
//! `every_mip_parser_flag_is_registered` checks every `"--flag" =>` arm of the native
//! `run_mip` parser against this table, and asks [`parse_flags`] itself whether the port
//! accepts each flag of its stack (its parser is this table, by design).
//! `a0_every_flag_kind_matches_its_placeholder` and `a0_help_json_types_every_flag` keep the
//! value kinds and the `--help --json` types in step with the placeholders, and
//! `a0_every_port_mip_flag_parses_from_its_registry_entry_alone` parses every port flag of
//! `mip` from nothing but its entry here. A new subcommand adds one arm in `main.rs` and one
//! entry here; a new flag of the port `mip` adds one entry here and the code that acts on it.

use serde_json::{Value, json};

/// The `schema_version` of `rusolve --help --json`.
pub const HELP_SCHEMA_VERSION: &str = "rusolve.help/1";

/// Which solver stack honours a flag. `rusolve mip` runs on the native `rusolve-mip`
/// port by default and on the native search with `RUSOLVE_MIP_HIGHS=0`; `rusolve solve` runs
/// on the HiGHS-port simplex by default (`ekk`) and on the legacy simplex or IPM with
/// `--legacy-lp` / `--ipm`. `Both` flags mean the same on both stacks.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Stack {
    Both,
    /// The HiGHS port only (the default engine for `mip` and `solve`).
    Port,
    /// The native `rusolve-mip` / legacy `rusolve-lp` stack only.
    Native,
}

impl Stack {
    pub fn as_str(self) -> &'static str {
        match self {
            Stack::Both => "both",
            Stack::Port => "port",
            Stack::Native => "native",
        }
    }
}

/// What the value of a flag must look like: the kind is what [`parse_flags`] checks,
/// what `--help --json` reports as the flag's `type`, and what a bindings generator would
/// use. The placeholder in usage is the human form of the same fact, and [`kind_of`] derives
/// the kind from it, so the two cannot drift apart.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ValueKind {
    /// No value: the flag is on or off.
    Switch,
    /// A signed integer, `--node-limit 100`.
    Int,
    /// A non-negative integer, `--wasm-fuel 1000`.
    UInt,
    /// A number, `--time-limit 12.5`.
    Float,
    /// A file, or `-` for stdin: any text. The file itself is checked when it is read
    /// (that is an input error, exit 3, not a usage error).
    Path,
    /// Any text, `--branch most_fractional`.
    Text,
    /// One of a fixed list, `--format json`.
    Choice(&'static [&'static str]),
    /// `key=value` with a non-empty key, `--opt mip_rel_gap=0.01`.
    KeyValue,
}

impl ValueKind {
    /// The `type` of this kind in `rusolve.help/1` and in a generated schema.
    pub fn type_name(self) -> &'static str {
        match self {
            ValueKind::Switch => "switch",
            ValueKind::Int => "int",
            ValueKind::UInt => "uint",
            ValueKind::Float => "float",
            ValueKind::Path => "path",
            ValueKind::Text => "text",
            ValueKind::Choice(_) => "choice",
            ValueKind::KeyValue => "key=value",
        }
    }

    /// The fixed list of a [`ValueKind::Choice`], empty for every other kind.
    pub fn choices(self) -> &'static [&'static str] {
        match self {
            ValueKind::Choice(list) => list,
            _ => &[],
        }
    }

    /// What the flag needs, worded for the end of "needs ...".
    pub fn expectation(self) -> String {
        match self {
            ValueKind::Switch => "no value".to_string(),
            ValueKind::Int => "an integer".to_string(),
            ValueKind::UInt => "a non-negative integer".to_string(),
            ValueKind::Float => "a number".to_string(),
            ValueKind::Path => "a path".to_string(),
            ValueKind::Text => "a value".to_string(),
            ValueKind::KeyValue => "key=value".to_string(),
            ValueKind::Choice(list) => format!("one of {}", or_list(list)),
        }
    }

    /// `true` iff `raw` is a value of this kind. Syntactic only: a path is not opened here
    /// and a number is not range-checked (the option registry, `--threads >= 1` and the
    /// cross-flag checks do that, so the messages stay with the checks that own them).
    fn accepts(self, raw: &str) -> bool {
        match self {
            ValueKind::Switch => false,
            ValueKind::Int => raw.parse::<i64>().is_ok(),
            ValueKind::UInt => raw.parse::<u64>().is_ok(),
            ValueKind::Float => raw.parse::<f64>().is_ok(),
            ValueKind::Path | ValueKind::Text => true,
            ValueKind::Choice(list) => list.contains(&raw),
            ValueKind::KeyValue => raw.split_once('=').is_some_and(|(k, _)| !k.is_empty()),
        }
    }
}

/// `["text", "csv", "json"]` as "text, csv or json".
fn or_list(list: &[&str]) -> String {
    match list {
        [] => String::new(),
        [one] => (*one).to_string(),
        [rest @ .., last] => format!("{} or {}", rest.join(", "), last),
    }
}

/// One flag of a command. `value` is the placeholder shown in usage (`None` = a switch) and
/// `kind` the same fact in the form [`parse_flags`] checks (see [`kind_of`]).
pub struct Flag {
    pub name: &'static str,
    pub value: Option<&'static str>,
    pub kind: ValueKind,
    pub default: Option<&'static str>,
    pub description: &'static str,
    pub stack: Stack,
    /// May be given more than once.
    pub repeatable: bool,
}

/// A positional argument.
pub struct Positional {
    pub name: &'static str,
    pub description: &'static str,
    pub required: bool,
}

/// One subcommand.
pub struct Command {
    pub name: &'static str,
    pub summary: &'static str,
    pub positionals: &'static [Positional],
    pub flags: &'static [Flag],
    /// The `schema_version` of its JSON output (`--format json`), if it has one.
    pub output_schema: Option<&'static str>,
}

/// String equality in a `const fn` (a `match` on `&str` is not const).
const fn streq(a: &str, b: &str) -> bool {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    if a.len() != b.len() {
        return false;
    }
    let mut i = 0;
    while i < a.len() {
        if a[i] != b[i] {
            return false;
        }
        i += 1;
    }
    true
}

/// The value kind a usage placeholder stands for, so that the human-readable placeholder and
/// the machine-readable kind of a flag are one fact, not two. Anything unrecognised is
/// `Text`, which accepts any value: a new placeholder is then at worst under-checked, never
/// wrong, and `every_flag_kind_matches_its_placeholder` below names it.
///
/// A placeholder that lists alternatives (`text|json`) is a [`ValueKind::Choice`], which the
/// table builds with [`val_c`]/[`val_cd`] because a `const fn` cannot split a string into a
/// `&'static [&'static str]`; `PATH|-` is the one exception, an alternative that is not a
/// value set but stdin.
const fn kind_of(placeholder: &str) -> ValueKind {
    if streq(placeholder, "N") || streq(placeholder, "K") {
        ValueKind::Int
    } else if streq(placeholder, "BYTES") {
        ValueKind::UInt
    } else if streq(placeholder, "SECONDS")
        || streq(placeholder, "S")
        || streq(placeholder, "F")
        || streq(placeholder, "W")
        || streq(placeholder, "X")
        || streq(placeholder, "G")
        || streq(placeholder, "V")
        || streq(placeholder, "M")
    {
        ValueKind::Float
    } else if streq(placeholder, "PATH")
        || streq(placeholder, "PATH|-")
        || streq(placeholder, "FILE.json")
    {
        ValueKind::Path
    } else if streq(placeholder, "KEY=VALUE") {
        ValueKind::KeyValue
    } else {
        ValueKind::Text
    }
}

const fn sw(name: &'static str, description: &'static str) -> Flag {
    Flag {
        name,
        value: None,
        kind: ValueKind::Switch,
        default: None,
        description,
        stack: Stack::Both,
        repeatable: false,
    }
}
const fn val(name: &'static str, value: &'static str, description: &'static str) -> Flag {
    Flag {
        name,
        value: Some(value),
        kind: kind_of(value),
        default: None,
        description,
        stack: Stack::Both,
        repeatable: false,
    }
}
const fn val_d(
    name: &'static str,
    value: &'static str,
    default: &'static str,
    description: &'static str,
) -> Flag {
    Flag {
        name,
        value: Some(value),
        kind: kind_of(value),
        default: Some(default),
        description,
        stack: Stack::Both,
        repeatable: false,
    }
}
/// A flag whose value is one of `choices`; `value` is the same list joined by `|`, the
/// placeholder `--help` prints.
const fn val_c(
    name: &'static str,
    value: &'static str,
    choices: &'static [&'static str],
    description: &'static str,
) -> Flag {
    Flag {
        kind: ValueKind::Choice(choices),
        ..val(name, value, description)
    }
}
const fn val_cd(
    name: &'static str,
    value: &'static str,
    choices: &'static [&'static str],
    default: &'static str,
    description: &'static str,
) -> Flag {
    Flag {
        kind: ValueKind::Choice(choices),
        ..val_d(name, value, default, description)
    }
}
/// The kind of a flag whose placeholder does not say it, `--wasm-fuel N` being a non-negative
/// integer rather than any integer.
const fn kind(k: ValueKind, f: Flag) -> Flag {
    Flag { kind: k, ..f }
}
const fn on(stack: Stack, f: Flag) -> Flag {
    Flag { stack, ..f }
}
const fn many(f: Flag) -> Flag {
    Flag {
        repeatable: true,
        ..f
    }
}

const fn pos(name: &'static str, description: &'static str) -> Positional {
    Positional {
        name,
        description,
        required: true,
    }
}

const MODEL: Positional = pos(
    "model",
    "model file: MPS, LP, AMPL .nl or MathOptFormat .mof.json, optionally .gz/.bz2/.xz/.zst",
);
/// The oracle-replay commands: the C++ HiGHS side reads MPS.
const MPS_MODEL: Positional = pos("model", "model file: MPS format, optionally .gz");

const FORMAT_TEXT_JSON: Flag = val_cd("--format", "text|json", TEXT_JSON, "text", "output format");
const TIME_LIMIT: Flag = val("--time-limit", "SECONDS", "wall-clock limit in seconds");

use Stack::{Native as N, Port as P};

/// The value lists of the `Choice` flags, in the order `--help` prints them.
static TEXT_JSON: &[&str] = &["text", "json"];
static TEXT_CSV_JSON: &[&str] = &["text", "csv", "json"];
static TEXT_JSONL: &[&str] = &["text", "jsonl"];
static DEFAULT_RECOMMENDED: &[&str] = &["default", "recommended"];
static RECOMMENDED_DEFAULT: &[&str] = &["recommended", "default"];
static SEVERITY: &[&str] = &["info", "warning", "error"];
static RELAX_MODE: &[&str] = &["weighted", "min-count"];

static MIP_FLAGS: &[Flag] = &[
    // both stacks accept it: on the native stack it is a no-op (no lint pass there)
    sw(
        "--no-warnings",
        "do not print the lint warnings (big-M, coefficient range, huge or free columns) before the solve and leave them out of the json diagnostics",
    ),
    // The shared surface: both stacks accept these with the same meaning.
    sw(
        "--csv",
        "print the harness CSV line `status,objective,bound,nodes,lp_iters,time_s,x`",
    ),
    val_cd(
        "--format",
        "text|csv|json",
        TEXT_CSV_JSON,
        "text",
        "output format; json prints one rusolve.result/1 object on stdout and the log on stderr",
    ),
    sw(
        "--solution",
        "include the incumbent by column name in the json result",
    ),
    sw("--include-solution", "same as --solution"),
    val(
        "--solution-file",
        "PATH",
        "write the incumbent to PATH: the json solution object when PATH ends in .json, else (port only) the HiGHS raw solution file",
    ),
    val("--time-limit", "SECONDS", "wall-clock limit in seconds"),
    val("--node-limit", "N", "branch-and-bound node limit"),
    val_d("--seed", "N", "0", "random seed"),
    val(
        "--threads",
        "N",
        "worker threads of the MIP search (default 1: serial; N > 1 runs the multi-worker tree search; the count is recorded in run.threads of the json result)",
    ),
    sw("--no-presolve", "turn presolve off"),
    // Port only.
    on(
        P,
        val(
            "--mip-rel-gap",
            "X",
            "relative MIP gap tolerance (HiGHS option mip_rel_gap)",
        ),
    ),
    on(
        P,
        val(
            "--mip-abs-gap",
            "X",
            "absolute MIP gap tolerance (HiGHS option mip_abs_gap)",
        ),
    ),
    many(on(
        P,
        val(
            "--opt",
            "KEY=VALUE",
            "set a HiGHS option by name; the names, types and ranges are `rusolve options --json`",
        ),
    )),
    on(
        P,
        val(
            "--options-file",
            "PATH",
            "read HiGHS options from a file (applied before --opt)",
        ),
    ),
    // Profile files (profile.rs)
    on(
        P,
        val(
            "--profile",
            "PATH",
            "apply a rusolve.profile/1 file (before --options-file and --opt)",
        ),
    ),
    // Progress events and the watchdog (progress.rs)
    on(
        P,
        val_c(
            "--log-format",
            "text|jsonl",
            TEXT_JSONL,
            "jsonl streams rusolve.event/1 progress events as NDJSON on stderr instead of the text log",
        ),
    ),
    on(
        P,
        sw(
            "--watchdog",
            "enforce --time-limit from a wall-clock thread: interrupt at the limit, hard stop after the grace",
        ),
    ),
    on(
        P,
        val(
            "--watchdog-grace",
            "SECONDS",
            "seconds after the time limit before the watchdog's hard stop (default 5)",
        ),
    ),
    // Warm-start inputs read before the solve (warm.rs)
    on(
        P,
        val(
            "--start",
            "PATH",
            "MIP start by column name: HiGHS solution file or the JSON --solution-file writes",
        ),
    ),
    on(
        P,
        val(
            "--basis",
            "PATH",
            "HiGHS basis file (v1/v2) for an LP warm start; read and checked but unused on a MIP",
        ),
    ),
    on(
        P,
        val(
            "--add-rows",
            "PATH",
            "add rows (MPS fragment or JSON list) before the solve, each through the row validity gate",
        ),
    ),
    // A sandboxed Wasm hook component (wasm-hooks.md).
    on(
        P,
        val(
            "--wasm-hook",
            "PATH",
            "run the .wasm/.wat hook component at PATH for branching, cut scoring and node selection (docs/design/wasm-hooks.md); the json result gains stats.wasm",
        ),
    ),
    on(
        P,
        val_d(
            "--gc-bigm-max",
            "X",
            "1000000",
            "largest big-M a linearized general constraint (indicator, SOS, PWL, abs/min/max) may carry; a larger one is an error naming the object",
        ),
    ),
    on(
        P,
        kind(
            ValueKind::UInt,
            val_d(
                "--wasm-fuel",
                "N",
                "1000000",
                "fuel (interpreter steps) per hook call with --wasm-hook; out of fuel falls back to the HiGHS rule for that call",
            ),
        ),
    ),
    on(
        P,
        val_d(
            "--wasm-memory",
            "BYTES",
            "1048576",
            "linear-memory cap in bytes for the --wasm-hook component",
        ),
    ),
    // Native levers: search switches measured on our benchmark set. The port reaches the same
    // ideas through --opt (see `rusolve options`).
    on(
        N,
        val(
            "--heur-off",
            "LIST",
            "disable heuristics: trivial,fj,simple,randomized,zi,shifting,diving or all",
        ),
    ),
    on(
        N,
        val(
            "--heur-on",
            "LIST",
            "enable heuristics (same names as --heur-off)",
        ),
    ),
    on(
        N,
        val(
            "--fj-restart-effort",
            "F",
            "feasibility-jump restart effort, >= 0",
        ),
    ),
    on(N, sw("--no-cuts", "no cutting planes")),
    on(N, sw("--tree-cuts", "cuts inside the tree")),
    on(N, sw("--no-tree-cuts", "no cuts inside the tree")),
    on(N, val("--tree-cut-effort", "X", "tree cut effort")),
    on(
        N,
        val(
            "--lp-age-limit",
            "N",
            "cut row aging limit (HiGHS mip_lp_age_limit)",
        ),
    ),
    on(N, sw("--tree-age-ticks", "age tree cuts by ticks")),
    on(N, sw("--no-tree-age-ticks", "age tree cuts by rounds")),
    on(N, sw("--modk", "mod-k cuts on")),
    on(N, sw("--no-modk", "mod-k cuts off")),
    on(N, sw("--probe-cliques", "clique probing on")),
    on(N, sw("--no-probe-cliques", "clique probing off")),
    on(N, sw("--clique-extend", "clique extension on")),
    on(N, sw("--no-clique-extend", "clique extension off")),
    on(N, sw("--gomory", "Gomory cuts on")),
    on(N, sw("--no-gomory", "Gomory cuts off")),
    on(
        N,
        val("--gomory-density", "X", "Gomory cut density bar, >= 0"),
    ),
    on(N, sw("--gomory-park", "park Gomory cuts")),
    on(N, sw("--no-gomory-park", "do not park Gomory cuts")),
    on(
        N,
        sw(
            "--drop-all-cuts-at-end",
            "drop every cut after the root loop",
        ),
    ),
    on(
        N,
        sw(
            "--drop-gomory-at-end",
            "drop Gomory cuts after the root loop",
        ),
    ),
    on(N, sw("--efficacy-rounds", "efficacy filter in the rounds")),
    on(N, sw("--efficacy-end", "efficacy filter at the end")),
    on(N, sw("--no-efficacy-end", "no efficacy filter at the end")),
    on(
        N,
        sw("--root-cutoff-row", "objective cutoff row at the root"),
    ),
    on(
        N,
        sw(
            "--no-root-cutoff-row",
            "no objective cutoff row at the root",
        ),
    ),
    on(
        N,
        sw("--root-cutoff-tighten", "tighten with the root cutoff"),
    ),
    on(
        N,
        sw("--no-root-cutoff-tighten", "no root cutoff tightening"),
    ),
    on(
        N,
        val(
            "--weak-rounds",
            "N",
            "weak root rounds before the loop stops",
        ),
    ),
    on(
        N,
        val(
            "--flat-weak-rounds",
            "N",
            "weak rounds allowed before the root bound first moves",
        ),
    ),
    on(N, val("--max-rounds", "N", "root cut round cap")),
    on(N, sw("--no-warm", "cold-start every node LP")),
    on(
        N,
        sw(
            "--no-warm-dse",
            "no dual steepest-edge weights on warm starts",
        ),
    ),
    on(
        N,
        sw("--warm-dse", "dual steepest-edge weights on warm starts"),
    ),
    on(
        N,
        sw(
            "--warm-unfit-last",
            "warm-start unfit snapshots from the latest basis",
        ),
    ),
    on(N, sw("--dual-lp", "dual simplex for node LPs")),
    on(N, sw("--root-dual", "dual simplex at the root")),
    on(N, sw("--no-root-dual", "primal simplex at the root")),
    on(N, sw("--no-redcost", "no reduced-cost fixing")),
    on(N, sw("--submip", "RENS/RINS sub-MIPs on")),
    on(N, sw("--no-submip", "RENS/RINS sub-MIPs off")),
    on(
        N,
        sw(
            "--root-dive",
            "dive-shaped root heuristic chain (implies --submip)",
        ),
    ),
    on(N, sw("--no-root-dive", "no root dive")),
    on(N, sw("--alns", "ALNS on")),
    on(N, sw("--no-alns", "ALNS off")),
    on(N, sw("--localmip", "local MIP on")),
    on(N, sw("--no-localmip", "local MIP off")),
    on(N, sw("--conflict", "conflict analysis on")),
    on(N, sw("--no-conflict", "conflict analysis off")),
    on(N, sw("--conflict-cutoff", "conflict proofs from cutoffs")),
    on(
        N,
        sw("--no-conflict-cutoff", "no conflict proofs from cutoffs"),
    ),
    on(N, sw("--conflict-score", "conflict-based branching score")),
    on(
        N,
        sw("--no-conflict-score", "no conflict-based branching score"),
    ),
    on(N, sw("--restart", "restarts on")),
    on(N, sw("--no-restart", "restarts off")),
    on(
        N,
        val("--restart-min-gap", "G", "root gap bar for a restart"),
    ),
    on(N, sw("--clique-fixings", "clique fixings on")),
    on(N, sw("--no-clique-fixings", "clique fixings off")),
    on(N, sw("--root-probing", "root probing on")),
    on(N, sw("--no-root-probing", "root probing off")),
    on(N, sw("--sb-invest", "strong-branching investment")),
    on(N, sw("--no-sb-invest", "no strong-branching investment")),
    on(
        N,
        val(
            "--cutoff",
            "V",
            "known objective cutoff, in the caller's sense",
        ),
    ),
    on(
        N,
        val(
            "--branch",
            "SPEC",
            "branching rule spec (BranchConfig::apply)",
        ),
    ),
];

static SOLVE_FLAGS: &[Flag] = &[
    sw(
        "--no-warnings",
        "do not print the lint warnings (big-M, coefficient range, huge or free columns) before the solve and leave them out of the json diagnostics",
    ),
    sw(
        "--csv",
        "print the harness CSV line `status,objective,phase1_iters,phase2_iters,x`",
    ),
    val_cd(
        "--format",
        "text|csv|json",
        TEXT_CSV_JSON,
        "text",
        "output format; json prints one rusolve.result/1 object on stdout",
    ),
    sw(
        "--solution",
        "include the solution by column name in the json result",
    ),
    sw("--include-solution", "same as --solution"),
    sw("--duals", "also print row duals and reduced costs"),
    val(
        "--solution-file",
        "PATH",
        "write the HiGHS raw solution file (primal and dual values) to PATH, for `rusolve check` (HiGHS-port engine)",
    ),
    sw(
        "--no-presolve",
        "turn presolve off (legacy engine and --ipm; the default engine ignores it)",
    ),
    val("--time-limit", "SECONDS", "wall-clock limit in seconds"),
    val(
        "--iter-limit",
        "N",
        "simplex iteration limit (default: none, the time limit governs)",
    ),
    on(
        N,
        sw("--ipm", "interior-point method instead of the simplex"),
    ),
    on(
        N,
        sw(
            "--legacy-lp",
            "the legacy rusolve-lp simplex instead of the HiGHS-port simplex",
        ),
    ),
    sw(
        "--ekk",
        "kept as a no-op: the HiGHS-port simplex is the default",
    ),
];

static CLASSIFY_FLAGS: &[Flag] = &[
    sw(
        "--json",
        "print the classification as the legacy JSON object",
    ),
    val_cd(
        "--format",
        "text|json",
        TEXT_JSON,
        "text",
        "json prints a rusolve.classify/1 object",
    ),
    sw(
        "--miplib-strict",
        "use the strict MIPLIB 2017 row-type rules",
    ),
];

// Profile files (profile.rs)
static PROFILE_FLAGS: &[Flag] = &[
    val_c(
        "--format",
        "text|json",
        TEXT_JSON,
        "output format of diff and show",
    ),
    sw("--json", "same as --format json"),
    val_cd(
        "--base",
        "default|recommended",
        DEFAULT_RECOMMENDED,
        "recommended",
        "base profile of `profile new`",
    ),
    val("--name", "NAME", "name recorded in the new profile"),
    many(val(
        "--opt",
        "KEY=VALUE",
        "an option of the new profile, validated against `rusolve options`",
    )),
    val(
        "--output",
        "PATH",
        "write the new profile to PATH (also -o)",
    ),
];

// Tuning (tune.rs, tune_race.rs)
static TUNE_FLAGS: &[Flag] = &[
    val(
        "--output",
        "PATH",
        "write the tuned rusolve.profile/1 to PATH (also -o)",
    ),
    val_c("--format", "text|json", TEXT_JSON, "output format"),
    val_cd(
        "--base",
        "recommended|default",
        RECOMMENDED_DEFAULT,
        "recommended",
        "profile the candidates change",
    ),
    val_d("--seeds", "N", "3", "seeds per instance"),
    val_d(
        "--holdout",
        "F",
        "0.25",
        "fraction of instances held out (at most 0.5)",
    ),
    val(
        "--budget-work",
        "W",
        "total work budget of the race, in the cost metric's unit (default 2e12 with \
         c1_work at --jobs 1, 1e7 with iters_nodes at --jobs > 1)",
    ),
    val(
        "--run-work-cap",
        "W",
        "work cap of one run, in the cost metric's unit (default 4e11 with c1_work, 2e6 with \
         iters_nodes)",
    ),
    val_d("--time-limit", "S", "30", "time limit of one run"),
    val_d("--jobs", "N", "1", "parallel runs"),
    val("--known", "PATH", "known optima for the wrong-answer check"),
    val(
        "--candidates",
        "PATH|-",
        "extra candidate option sets (e.g. an agent's proposals), one per line or JSON",
    ),
    sw(
        "--candidates-only",
        "race only the --candidates list, not the built-in neighbourhood",
    ),
    val_d("--max-candidates", "N", "16", "built-in neighbourhood size"),
    val_d("--race-seed", "N", "1", "seed of the race's own randomness"),
    val_d(
        "--iterations",
        "K",
        "1",
        "racing rounds: every round after the first races the elite of the previous round plus \
         candidates with one more option changed (iterated racing, irace)",
    ),
    val_d(
        "--elite",
        "N",
        "3",
        "survivors of a round that are carried on and mutated into the next round",
    ),
    val("--name", "NAME", "name recorded in the tuned profile"),
    sw(
        "--space",
        "print the options a race mutates and their values as JSON, then exit (no instances)",
    ),
];

static OPTIONS_FLAGS: &[Flag] = &[
    sw("--json", "print the option records as JSON"),
    val_c(
        "--profile",
        "default|recommended",
        DEFAULT_RECOMMENDED,
        "option profile to list (default: the profile a solve starts from)",
    ),
    sw(
        "--non-default",
        "only options that differ from the registry default",
    ),
];

static ENV_FLAGS: &[Flag] = &[sw(
    "--json",
    "print the environment variable registry as JSON",
)];

static CHECK_FLAGS: &[Flag] = &[
    FORMAT_TEXT_JSON,
    val("--max-findings", "N", "list at most N violations per kind"),
    val(
        "--primal-feasibility-tolerance",
        "X",
        "override primal_feasibility_tolerance",
    ),
    val(
        "--mip-feasibility-tolerance",
        "X",
        "override mip_feasibility_tolerance",
    ),
];

static LINT_FLAGS: &[Flag] = &[
    FORMAT_TEXT_JSON,
    many(val("--rule", "ID", "only run this rule")),
    val_c(
        "--min-severity",
        "info|warning|error",
        SEVERITY,
        "drop findings below this severity",
    ),
    val("--max-items", "N", "list at most N items per finding"),
    sw(
        "--with-solve",
        "also solve the model for the rules that need a result",
    ),
];

static EXPLAIN_FLAGS: &[Flag] = &[
    FORMAT_TEXT_JSON,
    sw(
        "--no-verify",
        "skip the independent verification of the infeasibility certificate",
    ),
    val("--max-members", "N", "list at most N members of the IIS"),
    sw(
        "--mip-iis",
        "explain an integer-only infeasibility (IIS by branch-and-bound tests)",
    ),
    val(
        "--per-test-node-limit",
        "N",
        "node limit of each --mip-iis feasibility test",
    ),
    TIME_LIMIT,
];

static DIAGNOSE_FLAGS: &[Flag] = &[FORMAT_TEXT_JSON, TIME_LIMIT];

static PROBE_FLAGS: &[Flag] = &[
    FORMAT_TEXT_JSON,
    val_d("--node-limit", "N", "50", "node limit of the probe run"),
    val_d(
        "--time-limit",
        "SECONDS",
        "60",
        "time limit of the probe run",
    ),
    val_d("--seed", "N", "0", "random seed"),
];

static RANGING_FLAGS: &[Flag] = &[
    FORMAT_TEXT_JSON,
    many(val(
        "--rows",
        "NAMES",
        "comma-separated row names or patterns to report",
    )),
    many(val(
        "--cols",
        "NAMES",
        "comma-separated column names or patterns to report",
    )),
];

static FEASRELAX_FLAGS: &[Flag] = &[
    FORMAT_TEXT_JSON,
    val_cd(
        "--mode",
        "weighted|min-count",
        RELAX_MODE,
        "weighted",
        "minimize the weighted violation or the number of relaxed items",
    ),
    val(
        "--global-lower-penalty",
        "W",
        "penalty of every lower bound",
    ),
    val(
        "--global-upper-penalty",
        "W",
        "penalty of every upper bound",
    ),
    val("--global-rhs-penalty", "W", "penalty of every row side"),
    val("--penalties", "FILE.json", "per-name penalties"),
    val("--big-m", "M", "big-M for the min-count mode"),
    TIME_LIMIT,
];

/// Every subcommand, in the order `rusolve --help` lists them.
pub static COMMANDS: &[Command] = &[
    Command {
        name: "solve",
        summary: "solve an LP (HiGHS-port simplex by default; --legacy-lp / --ipm for the native engines)",
        positionals: &[MODEL],
        flags: SOLVE_FLAGS,
        output_schema: Some("rusolve.result/1"),
    },
    Command {
        name: "mip",
        summary: "solve a MIP (HiGHS port with rusolve's levers by default; RUSOLVE_MIP_HIGHS=0 for the native search)",
        positionals: &[MODEL],
        flags: MIP_FLAGS,
        output_schema: Some("rusolve.result/1"),
    },
    Command {
        name: "check",
        summary: "check a solution against a model: every violated bound, row and integrality by name",
        positionals: &[
            MODEL,
            pos(
                "solution",
                "HiGHS raw solution file or JSON solution by column name",
            ),
        ],
        flags: CHECK_FLAGS,
        output_schema: Some("rusolve.check/1"),
    },
    Command {
        name: "explain",
        summary: "explain an infeasible model: an irreducible infeasible subsystem, by name",
        positionals: &[MODEL],
        flags: EXPLAIN_FLAGS,
        output_schema: Some("rusolve.explain/1"),
    },
    Command {
        name: "convert",
        summary: "write a model as CPLEX LP or MPS, by the output extension (.lp, .mps, .gz)",
        positionals: &[
            MODEL,
            pos(
                "output",
                "output file: .lp, .mps, .lp.gz or .mps.gz; the extension picks the format",
            ),
        ],
        flags: &[FORMAT_TEXT_JSON],
        output_schema: Some("rusolve.convert/1"),
    },
    Command {
        name: "lint",
        summary: "static model lint: numerics, structure and modelling smells, by name",
        positionals: &[MODEL],
        flags: LINT_FLAGS,
        output_schema: Some("rusolve.lint/1"),
    },
    Command {
        name: "diagnose",
        summary: "formulation or data? classify why a model is infeasible or unbounded",
        positionals: &[MODEL],
        flags: DIAGNOSE_FLAGS,
        output_schema: Some("rusolve.diagnose/1"),
    },
    Command {
        name: "ranging",
        summary: "sensitivity ranging of an LP optimum, by row and column name",
        positionals: &[MODEL],
        flags: RANGING_FLAGS,
        output_schema: Some("rusolve.ranging/1"),
    },
    Command {
        name: "feasrelax",
        summary: "cheapest feasibility relaxation of an infeasible model",
        positionals: &[MODEL],
        flags: FEASRELAX_FLAGS,
        output_schema: Some("rusolve.feasrelax/1"),
    },
    Command {
        name: "classify",
        summary: "static instance features: MIPLIB row types, column and objective statistics",
        positionals: &[MODEL],
        flags: CLASSIFY_FLAGS,
        output_schema: Some("rusolve.classify/1"),
    },
    Command {
        name: "probe",
        summary: "short run on the port that reports dynamic instance features",
        positionals: &[MODEL],
        flags: PROBE_FLAGS,
        output_schema: Some("rusolve.probe/1"),
    },
    Command {
        name: "profile",
        summary: "show, validate, create and diff rusolve.profile/1 option profiles",
        positionals: &[
            pos(
                "subcommand",
                "diff <a> <b> | show <p> | validate <file> | new",
            ),
            Positional {
                name: "profiles",
                description: "a profile file, `default` or `recommended`",
                required: false,
            },
        ],
        flags: PROFILE_FLAGS,
        output_schema: Some("rusolve.profile-diff/1"),
    },
    Command {
        name: "tune",
        summary: "race option candidates on an instance family and write the winning profile",
        positionals: &[pos(
            "instances",
            "a directory of models or a text file listing them",
        )],
        flags: TUNE_FLAGS,
        output_schema: Some("rusolve.tune/1"),
    },
    Command {
        name: "options",
        summary: "list every solver option with type, range, default and what it affects",
        positionals: &[],
        flags: OPTIONS_FLAGS,
        output_schema: None,
    },
    Command {
        name: "env",
        summary: "list every environment variable the solver reads",
        positionals: &[],
        flags: ENV_FLAGS,
        output_schema: None,
    },
    Command {
        name: "version",
        summary: "print the solver version and whether the CP-SAT engine's OR-Tools library is installed (no solve)",
        positionals: &[],
        flags: &[FORMAT_TEXT_JSON],
        output_schema: None,
    },
    Command {
        name: "mcp",
        summary: "Model Context Protocol server over stdio",
        positionals: &[],
        flags: &[],
        output_schema: None,
    },
    Command {
        name: "mip-trace",
        summary: "replay the HiGHS MIP oracle trace on the port (acceptance harness)",
        positionals: &[
            MPS_MODEL,
            Positional {
                name: "time_limit_s",
                description: "time limit in seconds",
                required: false,
            },
        ],
        flags: &[],
        output_schema: None,
    },
    Command {
        name: "presolve-dump",
        summary: "print the presolve reduction sequence of the port (acceptance harness)",
        positionals: &[MPS_MODEL],
        flags: &[],
        output_schema: None,
    },
    Command {
        name: "ekk-seq",
        summary: "replay the HiGHS API oracle sequence on the port (acceptance harness)",
        positionals: &[MPS_MODEL],
        flags: &[],
        output_schema: None,
    },
    Command {
        name: "dump",
        summary: "machine-readable model dump for differential testing",
        positionals: &[MODEL],
        flags: &[],
        output_schema: None,
    },
];

/// The exit-code table, the same for every command.
pub const EXIT_CODES: &[(i32, &str)] = &[
    (
        0,
        "the command ran; the result status is in the output, not in the code",
    ),
    (
        2,
        "usage error: unknown command or flag, a missing or unparsable value",
    ),
    (
        3,
        "input error: the model or solution file is missing, unreadable or does not parse",
    ),
    (4, "internal error"),
];

/// The registry entry of a command.
pub fn find(name: &str) -> Option<&'static Command> {
    COMMANDS.iter().find(|c| c.name == name)
}

/// The registered flag `flag` of `command`, if any.
pub fn find_flag(command: &str, flag: &str) -> Option<&'static Flag> {
    find(command)?.flags.iter().find(|f| f.name == flag)
}

/// The message for a flag that exists on the other stack only, `None` for any other flag.
/// `stack` is the stack the caller is. The text names the stack and the way to select it, so
/// an agent does not have to guess why a registered flag was refused.
pub fn flag_stack_error(command: &str, flag: &str, stack: Stack) -> Option<String> {
    let f = find_flag(command, flag)?;
    if f.stack == Stack::Both || f.stack == stack {
        return None;
    }
    let (have, how) = match (command, f.stack) {
        ("mip", Stack::Port) => ("the HiGHS port", "unset RUSOLVE_MIP_HIGHS"),
        ("mip", _) => ("the native search", "set RUSOLVE_MIP_HIGHS=0"),
        ("solve", Stack::Port) => ("the HiGHS-port simplex", "drop --legacy-lp and --ipm"),
        _ => ("the legacy engines", "pass --legacy-lp or --ipm"),
    };
    Some(format!(
        "rusolve {command}: flag '{flag}' is only supported by {have} ({how}; see `rusolve {command} --help`)"
    ))
}

fn flag_usage(f: &Flag) -> String {
    match f.value {
        Some(v) => format!("{} <{v}>", f.name),
        None => f.name.to_string(),
    }
}

// ------------------------------------------------------------------ the parser

/// A command line the registry refuses (the registry is what parses, so every surface
/// rejects the same way). The message always names the flag, what was given, and what the
/// flag expects, in the wording the usage error (exit 2) uses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UsageError {
    pub message: String,
    /// The diagnostic code (`RS0002` unknown flag, `RS0003` bad value, ...).
    pub code: &'static str,
    /// `help:` lines (did-you-mean).
    pub help: Vec<String>,
}

impl UsageError {
    fn say(message: impl Into<String>) -> Self {
        UsageError {
            message: message.into(),
            code: "RS0003",
            help: Vec::new(),
        }
    }
    fn coded(code: &'static str, message: impl Into<String>) -> Self {
        UsageError {
            code,
            ..Self::say(message)
        }
    }
    /// The diagnostic of this error: the message, its code and its did-you-mean help.
    pub fn diagnostic(&self) -> rusolve_core::diagnostic::Diagnostic {
        let mut d = rusolve_core::diagnostic::Diagnostic::error(self.code, self.message.clone());
        d.help = self.help.clone();
        d
    }
}

impl std::fmt::Display for UsageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for UsageError {}

/// One flag as it appeared on the command line. `value` is the raw text (`None` for a
/// switch); [`ParsedArgs`]'s typed getters read it.
#[derive(Clone, Debug)]
pub struct Given {
    pub name: &'static str,
    pub value: Option<String>,
}

/// A command line read through the registry: the positionals in order, and the flags **in
/// command-line order**, because the order of `--opt`, `--mip-rel-gap` and `--mip-abs-gap`
/// is the order their settings are applied in.
#[derive(Clone, Debug, Default)]
pub struct ParsedArgs {
    positionals: Vec<String>,
    flags: Vec<Given>,
}

impl ParsedArgs {
    /// The arguments that are not flags, in order.
    pub fn positionals(&self) -> &[String] {
        &self.positionals
    }

    /// Every flag given, in command-line order, repeats of a repeatable flag included.
    pub fn flags(&self) -> &[Given] {
        &self.flags
    }

    /// Was `name` given? Use this for a switch; the typed getters say `None` for one, as it
    /// carries no value.
    pub fn has(&self, name: &str) -> bool {
        self.flags.iter().any(|g| g.name == name)
    }

    fn value_of(&self, name: &str) -> Option<&str> {
        self.flags
            .iter()
            .find(|g| g.name == name)
            .and_then(|g| g.value.as_deref())
    }

    /// The value of `name` as text: the first occurrence, `None` when it was not given or
    /// is a switch.
    pub fn get_str(&self, name: &str) -> Option<&str> {
        self.value_of(name)
    }

    /// The value of `name` as an integer (`Int` flags). `None` when it was not given; a
    /// value that is not an integer never reaches here, [`parse_flags`] refused it.
    pub fn get_i64(&self, name: &str) -> Option<i64> {
        self.value_of(name).and_then(|v| v.parse().ok())
    }

    /// The value of `name` as a non-negative integer (`UInt` flags).
    pub fn get_u64(&self, name: &str) -> Option<u64> {
        self.value_of(name).and_then(|v| v.parse().ok())
    }

    /// The value of `name` as a number (`Float` flags).
    pub fn get_f64(&self, name: &str) -> Option<f64> {
        self.value_of(name).and_then(|v| v.parse().ok())
    }

    /// Every value given for a repeatable flag `name`, in command-line order.
    // The port `mip` parser reads its repeatable `--opt` through `flags()` (the order of the
    // option settings runs across flag names); this is the per-name reader the registry tests use
    // and the one the next command that moves its parser here needs.
    #[allow(dead_code)]
    pub fn all(&self, name: &str) -> Vec<&str> {
        self.flags
            .iter()
            .filter(|g| g.name == name)
            .filter_map(|g| g.value.as_deref())
            .collect()
    }
}

/// The message for a flag this command does not have at all.
pub fn unknown_flag_message(command: &str, flag: &str) -> String {
    format!("rusolve {command}: unknown flag '{flag}' (see `rusolve {command} --help`)")
}

fn value_error(command: &str, flag: &Flag, given: Option<&str>) -> UsageError {
    UsageError::say(format!(
        "rusolve {command}: {} needs {}, got {} (see `rusolve {command} --help`)",
        flag.name,
        flag.kind.expectation(),
        match given {
            Some(v) => format!("'{v}'"),
            None => "nothing".to_string(),
        }
    ))
}

/// Read `args` (everything after the command name) against the registry: positionals and
/// flags separated in one pass, every flag looked up here, its presence, kind, choice list
/// and repeatability checked here. What a flag *means* stays with the command (`mip_highs`
/// reads the result), and so do the checks between flags.
///
/// A token is a flag when it starts with `-`; the value of a flag that takes one is the next
/// token whatever it looks like, so `--time-limit -1` and `--candidates -` work. Each
/// rejection is a [`UsageError`] (exit 2) naming the flag: an unknown flag, a flag of the
/// other stack (`stack`, see [`flag_stack_error`]), a missing value, a value that is not of
/// the flag's kind or not in its list, and a repeat of a value flag that is not repeatable (a repeated switch is accepted).
pub fn parse_flags(command: &str, stack: Stack, args: &[String]) -> Result<ParsedArgs, UsageError> {
    let cmd = find(command).ok_or_else(|| {
        UsageError::say(format!(
            "rusolve {command}: not a registered command (see `rusolve --help`)"
        ))
    })?;
    let mut out = ParsedArgs {
        positionals: Vec::new(),
        flags: Vec::new(),
    };
    let mut next = 0;
    while next < args.len() {
        let token = args[next].as_str();
        next += 1;
        if !token.starts_with('-') {
            out.positionals.push(token.to_string());
            continue;
        }
        let flag = match cmd.flags.iter().find(|f| f.name == token) {
            Some(flag) => flag,
            None => {
                let mut e = UsageError::coded("RS0002", unknown_flag_message(cmd.name, token));
                // The closest flag of this command
                e.help.extend(rusolve_core::diagnostic::did_you_mean(
                    token,
                    cmd.flags.iter().map(|f| f.name),
                ));
                return Err(e);
            }
        };
        if let Some(msg) = flag_stack_error(cmd.name, flag.name, stack) {
            return Err(UsageError::coded("RS0002", msg));
        }
        // A repeated switch says the same thing twice and stays accepted, as before the registry; a
        // repeated value flag is ambiguous (which value wins?) and is refused.
        if out.has(flag.name) && flag.kind == ValueKind::Switch {
            continue;
        }
        if out.has(flag.name) && !flag.repeatable {
            return Err(UsageError::say(format!(
                "rusolve {}: flag '{}' may be given only once{} (see `rusolve {} --help`)",
                cmd.name,
                token,
                match out.value_of(flag.name) {
                    Some(v) => format!(", already given as '{v}'"),
                    None => String::new(),
                },
                cmd.name
            )));
        }
        let value = match flag.kind {
            ValueKind::Switch => None,
            kind => {
                let raw = args.get(next).map(String::as_str);
                match raw {
                    Some(raw) => {
                        next += 1;
                        if !kind.accepts(raw) {
                            return Err(value_error(cmd.name, flag, Some(raw)));
                        }
                        Some(raw.to_string())
                    }
                    None => return Err(value_error(cmd.name, flag, None)),
                }
            }
        };
        out.flags.push(Given {
            name: flag.name,
            value,
        });
    }
    Ok(out)
}

/// `usage: rusolve <cmd> <model> [flags]` for one command, with the flag table.
pub fn command_usage(c: &Command) -> String {
    let mut s = format!("usage: rusolve {}", c.name);
    for p in c.positionals {
        if p.required {
            s.push_str(&format!(" <{}>", p.name));
        } else {
            s.push_str(&format!(" [{}]", p.name));
        }
    }
    if !c.flags.is_empty() {
        s.push_str(" [flags]");
    }
    s.push_str("\n\n");
    s.push_str(c.summary);
    s.push('\n');
    if !c.positionals.is_empty() {
        s.push_str("\narguments:\n");
        for p in c.positionals {
            s.push_str(&format!("  {:<32} {}\n", p.name, p.description));
        }
    }
    if !c.flags.is_empty() {
        let stacked = c.flags.iter().any(|f| f.stack != Stack::Both);
        s.push_str("\nflags");
        if stacked {
            s.push_str(
                " ([port] / [native]: only that stack honours the flag; the others work on both)",
            );
        }
        s.push_str(":\n");
        for f in c.flags {
            let tag = match f.stack {
                Stack::Both => "",
                Stack::Port => "[port] ",
                Stack::Native => "[native] ",
            };
            let default = f
                .default
                .map(|d| format!(" (default {d})"))
                .unwrap_or_default();
            let rep = if f.repeatable { " (repeatable)" } else { "" };
            s.push_str(&format!(
                "  {:<32} {tag}{}{default}{rep}\n",
                flag_usage(f),
                f.description
            ));
        }
    }
    if c.name == "mip" {
        s.push_str(
            "\n`rusolve options` lists the names --opt accepts. RUSOLVE_MIP_HIGHS=0 selects the native search.\n",
        );
    }
    if let Some(schema) = c.output_schema {
        s.push_str(&format!("\njson output: schema_version \"{schema}\"\n"));
    }
    s.push_str("\nexit codes: 0 ran, 2 usage error, 3 input error, 4 internal error\n");
    s
}

/// The text of `rusolve --help` (also printed, to stderr, for a bare `rusolve`).
pub fn usage_text() -> String {
    let mut s = String::from(
        "usage: rusolve <command> [args]\n       rusolve <command> --help\n       rusolve --help --json\n\ncommands:\n",
    );
    for c in COMMANDS {
        s.push_str(&format!("  {:<14} {}\n", c.name, c.summary));
    }
    s.push_str("\n`rusolve --help --json` prints every command and flag as rusolve.help/1.\n");
    s.push_str(
        "`rusolve options` lists the solver options, `rusolve env` the environment variables.\n",
    );
    s.push_str("exit codes: 0 ran, 2 usage error, 3 input error, 4 internal error\n");
    s
}

fn flag_json(f: &Flag) -> Value {
    let mut json = json!({
        "name": f.name,
        "value": f.value,
        "type": f.kind.type_name(),
        "default": f.default,
        "description": f.description,
        "stack": f.stack.as_str(),
        "repeatable": f.repeatable,
    });
    let choices = f.kind.choices();
    if !choices.is_empty() {
        json["choices"] = json!(choices);
    }
    json
}

/// One command as an entry of `rusolve.help/1`.
pub fn command_help_json(c: &Command) -> Value {
    let usage = command_usage(c);
    let usage_line = usage.lines().next().unwrap_or_default();
    json!({
        "name": c.name,
        "summary": c.summary,
        "usage": usage_line,
        "positionals": c.positionals.iter().map(|p| json!({
            "name": p.name, "description": p.description, "required": p.required,
        })).collect::<Vec<_>>(),
        "flags": c.flags.iter().map(flag_json).collect::<Vec<_>>(),
        "output_schema": c.output_schema,
    })
}

fn help_doc(commands: Vec<Value>) -> Value {
    let exit_codes: serde_json::Map<String, Value> = EXIT_CODES
        .iter()
        .map(|(code, text)| (code.to_string(), Value::String((*text).to_string())))
        .collect();
    json!({
        "schema_version": HELP_SCHEMA_VERSION,
        "name": "rusolve",
        "version": env!("CARGO_PKG_VERSION"),
        "commands": commands,
        "exit_codes": exit_codes,
        "options": "rusolve options --json",
        "environment": "rusolve env --json",
    })
}

/// `rusolve --help --json`: the whole command registry, schema `rusolve.help/1`.
pub fn help_json() -> Value {
    help_doc(COMMANDS.iter().map(command_help_json).collect())
}

/// `rusolve <cmd> --help --json`: the same object with the one command.
pub fn command_help_json_doc(c: &Command) -> Value {
    help_doc(vec![command_help_json(c)])
}

/// The help request in `args` (everything after the program name), if there is one:
/// `Some(exit code)` after printing it. `rusolve --help [--json]`, `rusolve -h`,
/// `rusolve <cmd> --help [--json]`. A command that is not registered falls through (the
/// dispatcher reports it).
pub fn handle_help(args: &[String]) -> Option<std::process::ExitCode> {
    let first = args.first()?;
    let is_help = |a: &String| a == "--help" || a == "-h";
    let json = args.iter().any(|a| a == "--json");
    let text = if is_help(first) {
        if json {
            help_json().to_string()
        } else {
            usage_text().trim_end().to_string()
        }
    } else {
        let c = find(first)?;
        if !args[1..].iter().any(is_help) {
            return None;
        }
        if json {
            command_help_json_doc(c).to_string()
        } else {
            command_usage(c).trim_end().to_string()
        }
    };
    println!("{text}");
    Some(std::process::ExitCode::SUCCESS)
}

/// Every HiGHS option name with its JSON-schema type, for the MCP `solve` tool's `options`
/// property: one source, `highs_option_records::OPTION_RECORDS`. Types are permissive
/// (a string is accepted for every type, as the HiGHS setters do).
pub fn option_properties() -> serde_json::Map<String, Value> {
    use rusolve_ekk::lp_data::highs_option_records::{OPTION_RECORDS, RecordKind};
    let mut m = serde_json::Map::new();
    for r in OPTION_RECORDS {
        let ty = match r.kind {
            RecordKind::Bool { .. } => json!(["boolean", "string"]),
            RecordKind::Int { .. } => json!(["integer", "string"]),
            RecordKind::Double { .. } => json!(["number", "string"]),
            RecordKind::Str { .. } => json!("string"),
        };
        m.insert(r.name.to_string(), json!({ "type": ty }));
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    /// The names after `needle` up to the next `"` in `src`, e.g. `cmd == "` -> `solve`.
    fn quoted_after(src: &str, needle: &str) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        let mut rest = src;
        while let Some(at) = rest.find(needle) {
            rest = &rest[at + needle.len()..];
            if let Some(end) = rest.find('"') {
                out.insert(rest[..end].to_string());
            }
        }
        out
    }

    #[test]
    fn every_dispatch_arm_has_a_registry_entry_and_vice_versa() {
        let main = include_str!("main.rs");
        // the scan must not see this test module's own needles: main.rs is a different file
        let mut arms = quoted_after(main, "cmd == \"");
        arms.extend(quoted_after(main, "cmd != \""));
        let registered: BTreeSet<String> = COMMANDS.iter().map(|c| c.name.to_string()).collect();
        assert_eq!(
            registered.len(),
            COMMANDS.len(),
            "a command is registered twice"
        );
        let missing: Vec<_> = arms.difference(&registered).collect();
        let dead: Vec<_> = registered.difference(&arms).collect();
        assert!(
            missing.is_empty(),
            "dispatch arms in main.rs with no entry in commands::COMMANDS: {missing:?}"
        );
        assert!(
            dead.is_empty(),
            "commands::COMMANDS entries with no dispatch arm in main.rs: {dead:?}"
        );
    }

    /// The `"--flag" =>` arms (and `"--a" | "--b" =>` groups) of a parser source.
    fn flag_arms(src: &str) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        for line in src.lines() {
            let t = line.trim_start().trim_start_matches("| ");
            if !t.starts_with("\"--") {
                continue;
            }
            let head = t.split("=>").next().unwrap();
            for part in head.split('|') {
                let part = part.trim();
                if part.starts_with("\"--") && part.ends_with('"') && !part.contains(' ') {
                    out.insert(part.trim_matches('"').to_string());
                }
            }
        }
        out
    }

    #[test]
    fn every_mip_parser_flag_is_registered() {
        let mip = find("mip").unwrap();
        let registered: BTreeSet<String> = mip.flags.iter().map(|f| f.name.to_string()).collect();
        let main = include_str!("main.rs");
        let start = main.find("fn run_mip(").unwrap();
        let end = main[start..].find("fn run_classify").unwrap() + start;
        let native = flag_arms(&main[start..end]);
        let unregistered: Vec<_> = native.difference(&registered).collect();
        assert!(
            unregistered.is_empty(),
            "native `mip` parser flags missing from commands::MIP_FLAGS: {unregistered:?}"
        );
        // The port parser *is* this table (`commands::parse_flags`), so "parsed by the
        // port" is "accepted on the Port stack"; the native parser is still a hand-written
        // `match`, so its arms are read from the source.
        for f in mip.flags {
            let by_port = parse_flags("mip", Stack::Port, &flag_on_a_command_line(f)).is_ok();
            match f.stack {
                Stack::Both => {
                    assert!(
                        native.contains(f.name),
                        "{} is `both` but not parsed natively",
                        f.name
                    );
                    assert!(by_port, "{} is `both` but the port refuses it", f.name);
                }
                Stack::Native => {
                    assert!(
                        native.contains(f.name),
                        "{} is native but not parsed natively",
                        f.name
                    );
                    assert!(!by_port, "{} is native but the port accepts it", f.name);
                }
                Stack::Port => assert!(by_port, "{} is port but the port refuses it", f.name),
            }
        }
    }

    #[test]
    fn flag_names_are_unique_per_command() {
        for c in COMMANDS {
            let names: BTreeSet<_> = c.flags.iter().map(|f| f.name).collect();
            assert_eq!(names.len(), c.flags.len(), "{} repeats a flag", c.name);
        }
    }

    // ------------------------------------------------------------------ registry

    fn argv(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| s.to_string()).collect()
    }

    /// A value of each kind, so that a flag that exists only in this table still parses.
    fn sample_value(kind: ValueKind) -> Option<String> {
        Some(match kind {
            ValueKind::Switch => return None,
            ValueKind::Int => "7".to_string(),
            ValueKind::UInt => "1000".to_string(),
            ValueKind::Float => "1.5".to_string(),
            ValueKind::Path => "x.mps".to_string(),
            ValueKind::Text => "most_fractional".to_string(),
            ValueKind::Choice(list) => list[0].to_string(),
            ValueKind::KeyValue => "mip_rel_gap=0.01".to_string(),
        })
    }

    /// The same flag as it would appear on a command line: `--flag` or `--flag value`.
    fn flag_on_a_command_line(f: &Flag) -> Vec<String> {
        let mut args = vec![f.name.to_string()];
        args.extend(sample_value(f.kind));
        args
    }

    /// The kind of every flag must be what its usage placeholder says: a switch is a flag
    /// with no placeholder and every other kind has one, a `Choice` says its list in the
    /// placeholder, and the rest are what [`kind_of`] derives from it. One exception is
    /// spelled out below, and `PATH|-` is not an exception but stdin as a path.
    #[test]
    fn a0_every_flag_kind_matches_its_placeholder() {
        for c in COMMANDS {
            for f in c.flags {
                let what = format!("{} {}", c.name, f.name);
                let Some(placeholder) = f.value else {
                    assert_eq!(f.kind, ValueKind::Switch, "{what} has no placeholder");
                    continue;
                };
                assert!(!placeholder.is_empty(), "{what} has an empty placeholder");
                assert_ne!(f.kind, ValueKind::Switch, "{what} takes no value");
                let choices = f.kind.choices();
                if !choices.is_empty() {
                    assert_eq!(placeholder, choices.join("|"), "{what}");
                    continue;
                }
                if placeholder.contains('|') {
                    // `--candidates PATH|-`: an alternative that is a path (stdin), not a list
                    assert_eq!(f.kind, ValueKind::Path, "{what}");
                    continue;
                }
                let expected = match f.name {
                    // the one flag whose placeholder is weaker than its kind: fuel is counted
                    "--wasm-fuel" => ValueKind::UInt,
                    _ => kind_of(placeholder),
                };
                assert_eq!(f.kind, expected, "{what} takes <{placeholder}>");
            }
        }
    }

    /// The same fact in the machine-readable surface: every flag of `rusolve --help --json`
    /// carries its `type`, and a `Choice` its `choices` (the schema a bindings
    /// generator reads is the table the parser reads).
    #[test]
    fn a0_help_json_types_every_flag() {
        for c in COMMANDS {
            let flags = command_help_json(c)["flags"].as_array().unwrap().clone();
            assert_eq!(flags.len(), c.flags.len(), "{} flags", c.name);
            for (f, json) in c.flags.iter().zip(&flags) {
                let what = format!("{} {}", c.name, f.name);
                assert_eq!(json["type"], f.kind.type_name(), "{what}");
                match f.kind {
                    ValueKind::Choice(list) => assert_eq!(json["choices"], json!(list), "{what}"),
                    _ => assert!(json.get("choices").is_none(), "{what}"),
                }
            }
        }
    }

    /// Every way a flag can be wrong is a refusal that names the flag, the value and what
    /// the flag expects (exit 2 wording, "no silent switches").
    #[test]
    fn a0_parse_flags_refuses_every_malformed_use_of_a_flag() {
        let cases: Vec<(&[&str], &str, &str)> = vec![
            (&["--no-such-flag"], "--no-such-flag", "unknown flag"),
            (&["--seed"], "--seed", "an integer"),
            (&["--seed", "abc"], "--seed", "an integer"),
            (&["--time-limit", "soon"], "--time-limit", "a number"),
            (
                &["--format", "yaml"],
                "--format",
                "one of text, csv or json",
            ),
            (
                &["--wasm-fuel", "-1"],
                "--wasm-fuel",
                "a non-negative integer",
            ),
            (&["--opt", "mip_rel_gap"], "--opt", "key=value"),
            (&["--no-cuts"], "--no-cuts", "native search"),
            (
                &["--seed", "1", "--seed", "2"],
                "--seed",
                "may be given only once",
            ),
        ];
        for (args, flag, expects) in cases {
            let err = match parse_flags("mip", Stack::Port, &argv(args)) {
                Ok(parsed) => panic!("`{args:?}` parsed as {parsed:?}"),
                Err(err) => err,
            };
            assert!(err.message.contains(flag), "`{args:?}` names not: {err}");
            assert!(
                err.message.contains(expects),
                "`{args:?}` says not {expects}: {err}"
            );
            assert!(err.message.starts_with("rusolve mip:"), "{err}");
        }
        // and the same on the native stack, where the port's flags are the other stack's
        let err = parse_flags("mip", Stack::Native, &argv(&["--opt", "presolve=off"])).unwrap_err();
        assert!(
            err.message.contains("--opt") && err.message.contains("HiGHS port"),
            "{err}"
        );
    }

    /// The registry contract, in one line: a flag that is only in this table is a flag the port
    /// `mip` parser accepts, with its value typed. Nothing in `mip_highs.rs` names them.
    #[test]
    fn a0_every_port_mip_flag_parses_from_its_registry_entry_alone() {
        let mip = find("mip").unwrap();
        let mut parsed_line = vec!["tiny.mps".to_string()];
        for f in mip.flags.iter().filter(|f| f.stack != Stack::Native) {
            let args = flag_on_a_command_line(f);
            let parsed = parse_flags("mip", Stack::Port, &args).unwrap_or_else(|e| {
                panic!("{} <{:?}> is a registered port flag: {e}", f.name, f.value)
            });
            assert!(parsed.has(f.name), "{} is not read", f.name);
            assert_eq!(
                parsed.get_str(f.name),
                sample_value(f.kind).as_deref(),
                "{}",
                f.name
            );
            parsed_line.extend(args);
        }
        // together, in one command line: order kept, repeats of a repeatable flag kept
        let all = parse_flags("mip", Stack::Port, &parsed_line).unwrap();
        assert_eq!(all.positionals(), ["tiny.mps".to_string()]);
        assert_eq!(
            all.flags().len(),
            mip.flags
                .iter()
                .filter(|f| f.stack != Stack::Native)
                .count()
        );
    }

    /// Typed getters, command-line order, and the repeatable flag's every value.
    #[test]
    fn a0_a_repeated_switch_is_accepted_once() {
        let parsed = parse_flags("mip", Stack::Port, &argv(&["--solution", "--solution"]))
            .expect("a repeated switch parses");
        assert!(parsed.has("--solution"));
    }

    #[test]
    fn a0_parsed_args_reads_values_in_argument_order() {
        let parsed = parse_flags(
            "mip",
            Stack::Port,
            &argv(&[
                "tiny.mps",
                "--opt",
                "mip_rel_gap=0.5",
                "--no-presolve",
                "--opt",
                "presolve=off",
                "--node-limit",
                "10",
                "--time-limit",
                "2.5",
                "--wasm-fuel",
                "4096",
                "--format",
                "json",
            ]),
        )
        .unwrap();
        assert_eq!(parsed.positionals(), ["tiny.mps".to_string()]);
        assert_eq!(
            parsed.all("--opt"),
            ["mip_rel_gap=0.5", "presolve=off"],
            "a repeatable flag keeps every value in order"
        );
        assert_eq!(parsed.get_i64("--node-limit"), Some(10));
        assert_eq!(parsed.get_f64("--time-limit"), Some(2.5));
        assert_eq!(parsed.get_u64("--wasm-fuel"), Some(4096));
        assert_eq!(parsed.get_str("--format"), Some("json"));
        assert!(parsed.has("--no-presolve"));
        assert!(!parsed.has("--csv"));
        assert_eq!(
            parsed.get_str("--no-presolve"),
            None,
            "a switch carries no value"
        );
        assert_eq!(
            parsed.flags().iter().map(|g| g.name).collect::<Vec<_>>(),
            [
                "--opt",
                "--no-presolve",
                "--opt",
                "--node-limit",
                "--time-limit",
                "--wasm-fuel",
                "--format"
            ],
            "the settings order is the argument order"
        );
        // a value that starts with a dash is still the value of the flag before it
        let negative = parse_flags("mip", Stack::Port, &argv(&["--time-limit", "-1"])).unwrap();
        assert_eq!(negative.get_f64("--time-limit"), Some(-1.0));
    }

    #[test]
    fn mcp_option_properties_come_from_the_option_registry() {
        let props = option_properties();
        assert_eq!(
            props.len(),
            rusolve_ekk::lp_data::highs_option_records::option_record_count()
        );
        assert!(props.contains_key("mip_rel_gap"));
        assert!(props.contains_key("presolve"));
    }
}
