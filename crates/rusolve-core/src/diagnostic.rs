//! The shared rustc-style diagnostic (`docs/design/diagnostics.md`).
//!
//! One [`Diagnostic`] value (code `RSnnnn`, severity, optional [`Location`], message, notes, help
//! lines) renders three ways: [`Diagnostic::render_text`] (rustc layout: source line, caret,
//! `= help:`), [`Diagnostic::to_json`] (`rusolve.diagnostics/1`, schema in
//! `docs/schema/rusolve-diagnostics-1.schema.json`) and [`Diagnostic::mcp_error_content`] (the text
//! of an MCP tool result with `isError: true`, carrying the JSON).
//!
//! Every code lives in [`REGISTRY`] with a long explanation (`rusolve --explain RSnnnn`,
//! `docs/agent/diagnostics.md`). Tests fail when a source file names an unregistered code or a
//! registered code has no explanation or is missing from the generated page.
//!
//! Provenance: the layout follows the rustc diagnostic guide
//! (<https://rustc-dev-guide.rust-lang.org/diagnostics.html>, error codes and `--explain`:
//! <https://doc.rust-lang.org/error_codes/>); the did-you-mean threshold is rustc's
//! `find_best_match_for_name` (distance at most `max(len, 3) / 3`, over the optimal string
//! alignment distance, which counts a transposition as one edit; Damerau 1964).
//! No dependencies: the JSON writer is hand-rolled so `rusolve-io` and `rusolve-model` can use this.

use std::fmt::Write as _;

/// How serious a diagnostic is. `Error` fails the command; the others never do.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
    Note,
    Help,
}

impl Severity {
    pub fn as_str(self) -> &'static str {
        match self {
            Severity::Error => "error",
            Severity::Warning => "warning",
            Severity::Note => "note",
            Severity::Help => "help",
        }
    }
}

/// Where a diagnostic points.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Location {
    /// A line (1-based) and column (1-based, in chars) of a text file, with the source line.
    /// `len` is the width of the caret underline (at least 1).
    File {
        path: String,
        line: usize,
        column: usize,
        len: usize,
        source_line: String,
    },
    /// A constraint by name, with the row printed out (`c1: x + y <= 4`).
    Row { name: String, rendered: String },
    /// A column (variable) by name, with its declaration printed out (`x in [0, 10], integer`).
    Column { name: String, rendered: String },
    /// A line of Python source (a model built with `rusolve.Model`).
    Python {
        path: String,
        line: usize,
        source_line: String,
    },
}

/// One diagnostic. Build with [`Diagnostic::error`] and the `with_*` methods.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Diagnostic {
    pub code: &'static str,
    pub severity: Severity,
    pub message: String,
    pub location: Option<Location>,
    pub notes: Vec<String>,
    pub help: Vec<String>,
}

impl Diagnostic {
    pub fn new(severity: Severity, code: &'static str, message: impl Into<String>) -> Self {
        Diagnostic {
            code,
            severity,
            message: message.into(),
            location: None,
            notes: Vec::new(),
            help: Vec::new(),
        }
    }
    pub fn error(code: &'static str, message: impl Into<String>) -> Self {
        Self::new(Severity::Error, code, message)
    }
    pub fn warning(code: &'static str, message: impl Into<String>) -> Self {
        Self::new(Severity::Warning, code, message)
    }
    pub fn with_location(mut self, location: Location) -> Self {
        self.location = Some(location);
        self
    }
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }
    pub fn with_help(mut self, help: impl Into<String>) -> Self {
        self.help.push(help.into());
        self
    }

    /// rustc-style text, ending in a newline:
    ///
    /// ```text
    /// error[RS0002]: unknown flag '--time-limt'
    ///   --> model.mps:3:5
    ///    |
    ///  3 | ROWS x
    ///    |      ^
    ///    = note: ...
    ///    = help: did you mean '--time-limit'?
    /// ```
    pub fn render_text(&self) -> String {
        let mut o = String::new();
        let _ = writeln!(
            o,
            "{}[{}]: {}",
            self.severity.as_str(),
            self.code,
            self.message
        );
        let mut gutter = 1usize;
        match &self.location {
            Some(Location::File {
                path,
                line,
                column,
                len,
                source_line,
            }) => {
                gutter = line.to_string().len();
                let pad = " ".repeat(gutter);
                let _ = writeln!(o, "{pad}--> {path}:{line}:{column}");
                let _ = writeln!(o, "{pad} |");
                let _ = writeln!(o, "{line} | {source_line}");
                let _ = writeln!(
                    o,
                    "{pad} | {}{}",
                    " ".repeat(column.saturating_sub(1)),
                    "^".repeat((*len).max(1))
                );
            }
            Some(Location::Python {
                path,
                line,
                source_line,
            }) => {
                gutter = line.to_string().len();
                let pad = " ".repeat(gutter);
                let _ = writeln!(o, "{pad}--> {path}:{line}");
                let _ = writeln!(o, "{pad} |");
                let _ = writeln!(o, "{line} | {source_line}");
            }
            Some(Location::Row { name, rendered }) => {
                let _ = writeln!(o, " --> row '{name}'");
                let _ = writeln!(o, "  = row: {rendered}");
            }
            Some(Location::Column { name, rendered }) => {
                let _ = writeln!(o, " --> column '{name}'");
                let _ = writeln!(o, "  = column: {rendered}");
            }
            None => {}
        }
        let pad = " ".repeat(gutter);
        for n in &self.notes {
            let _ = writeln!(o, "{pad} = note: {n}");
        }
        for h in &self.help {
            let _ = writeln!(o, "{pad} = help: {h}");
        }
        o
    }

    /// The `rusolve.diagnostics/1` object of this diagnostic (no envelope).
    pub fn to_json(&self) -> String {
        let mut o = String::new();
        o.push('{');
        let _ = write!(
            o,
            "\"code\":{},\"severity\":{},\"message\":{}",
            json_str(self.code),
            json_str(self.severity.as_str()),
            json_str(&self.message)
        );
        if let Some(l) = &self.location {
            o.push_str(",\"location\":");
            match l {
                Location::File {
                    path,
                    line,
                    column,
                    len,
                    source_line,
                } => {
                    let _ = write!(
                        o,
                        "{{\"kind\":\"file\",\"path\":{},\"line\":{line},\"column\":{column},\"length\":{},\"source_line\":{}}}",
                        json_str(path),
                        (*len).max(1),
                        json_str(source_line)
                    );
                }
                Location::Row { name, rendered } => {
                    let _ = write!(
                        o,
                        "{{\"kind\":\"row\",\"name\":{},\"rendered\":{}}}",
                        json_str(name),
                        json_str(rendered)
                    );
                }
                Location::Column { name, rendered } => {
                    let _ = write!(
                        o,
                        "{{\"kind\":\"column\",\"name\":{},\"rendered\":{}}}",
                        json_str(name),
                        json_str(rendered)
                    );
                }
                Location::Python {
                    path,
                    line,
                    source_line,
                } => {
                    let _ = write!(
                        o,
                        "{{\"kind\":\"python\",\"path\":{},\"line\":{line},\"source_line\":{}}}",
                        json_str(path),
                        json_str(source_line)
                    );
                }
            }
        }
        let _ = write!(
            o,
            ",\"notes\":{},\"help\":{}}}",
            json_arr(&self.notes),
            json_arr(&self.help)
        );
        o
    }

    /// The text of an MCP tool result with `isError: true`: the same JSON document as
    /// [`diagnostics_json`] for this one diagnostic, so an agent parses one shape on every surface.
    pub fn mcp_error_content(&self, exit_code: i32) -> String {
        diagnostics_json(std::slice::from_ref(self), exit_code)
    }
}

/// The `rusolve.diagnostics/1` document: what `--format json` prints on stdout for an error
/// (`exit_code` is the process exit code: 2 usage, 3 input, 4 internal).
pub fn diagnostics_json(diags: &[Diagnostic], exit_code: i32) -> String {
    let items: Vec<String> = diags.iter().map(Diagnostic::to_json).collect();
    format!(
        "{{\"schema_version\":\"rusolve.diagnostics/1\",\"exit_code\":{exit_code},\"diagnostics\":[{}]}}",
        items.join(",")
    )
}

fn json_arr(items: &[String]) -> String {
    let v: Vec<String> = items.iter().map(|s| json_str(s)).collect();
    format!("[{}]", v.join(","))
}

/// A JSON string literal with the escapes RFC 8259 requires.
pub fn json_str(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(o, "\\u{:04x}", c as u32);
            }
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

// ---------------------------------------------------------------------------------- did you mean

/// Optimal string alignment distance (insert, delete, substitute, swap of two neighbours).
pub fn edit_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let mut d = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i;
    }
    for (j, cell) in d[0].iter_mut().enumerate() {
        *cell = j;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            let mut v = (d[i - 1][j] + 1)
                .min(d[i][j - 1] + 1)
                .min(d[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                v = v.min(d[i - 2][j - 2] + 1);
            }
            d[i][j] = v;
        }
    }
    d[a.len()][b.len()]
}

/// Up to `max` candidates close to `name`, closest first (ties: shorter, then alphabetical). A
/// candidate qualifies at distance 1..=`max(len, 3) / 3` (rustc's rule), or when it differs from
/// `name` only in case. The name itself never qualifies.
pub fn suggestions<'a, I>(name: &str, candidates: I, max: usize) -> Vec<String>
where
    I: IntoIterator<Item = &'a str>,
{
    let limit = name.chars().count().max(3) / 3;
    let lower = name.to_lowercase();
    let mut scored: Vec<(usize, usize, String)> = candidates
        .into_iter()
        .filter(|c| *c != name)
        .filter_map(|c| {
            let d = if c.to_lowercase() == lower {
                0
            } else {
                edit_distance(name, c)
            };
            (d <= limit).then(|| (d, c.len(), c.to_string()))
        })
        .collect();
    scored.sort();
    scored.dedup();
    scored.into_iter().take(max).map(|(_, _, c)| c).collect()
}

/// The `help:` line for a did-you-mean over `candidates`, `None` when nothing is close.
/// One candidate reads `did you mean '--x'?`, several `closest matches are '--x', '--y'`.
pub fn did_you_mean<'a, I>(name: &str, candidates: I) -> Option<String>
where
    I: IntoIterator<Item = &'a str>,
{
    let s = suggestions(name, candidates, 3);
    match s.as_slice() {
        [] => None,
        [one] => Some(format!("did you mean '{one}'?")),
        many => Some(format!(
            "closest matches are {}",
            many.iter()
                .map(|m| format!("'{m}'"))
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
}

// ---------------------------------------------------------------------------------- the registry

/// One registered code.
#[derive(Debug)]
pub struct CodeEntry {
    pub code: &'static str,
    pub severity: Severity,
    pub title: &'static str,
    /// The long explanation `rusolve --explain` prints: what it means, why it happens, the fix.
    pub explanation: &'static str,
}

/// Every code the solver emits. RS0xxx: command line and files; RS1xxx: model readers.
/// A new code is added here and nowhere else; `tools/gen_agent_docs.py` regenerates
/// `docs/agent/diagnostics.md`.
pub static REGISTRY: &[CodeEntry] = &[
    CodeEntry {
        code: "RS0001",
        severity: Severity::Error,
        title: "unknown command",
        explanation: "The first argument is not a rusolve command. `rusolve --help` lists the commands (`rusolve --help --json` for the registry). The help line names the closest command when one is within a few edits, for example `slove` for `solve`.",
    },
    CodeEntry {
        code: "RS0002",
        severity: Severity::Error,
        title: "unknown flag",
        explanation: "The command does not have this flag. `rusolve <command> --help` lists its flags. The help line names the closest flag when one is within a few edits, for example `--time-limt` for `--time-limit`. A flag that exists only on the other stack (HiGHS port or native/legacy) gets its own message naming the switch that selects that stack.",
    },
    CodeEntry {
        code: "RS0003",
        severity: Severity::Error,
        title: "missing or invalid flag value",
        explanation: "A flag needs a value it did not get, got one of the wrong kind (a number expected, text given), one outside its list of choices, or was given twice although it takes one value. Other flag combinations the command refuses (`--csv` with `--format json`) also use this code. The message names the flag and what it expects; `rusolve <command> --help` shows the type and choices.",
    },
    CodeEntry {
        code: "RS0004",
        severity: Severity::Error,
        title: "unknown solver option",
        explanation: "`--opt name=value`, `--options-file` or the MCP `options` map named an option the registry does not have, or gave a value it rejects. `rusolve options` lists every name with type, range and default. The help line names the closest option, for example `mip_rel_gapp` for `mip_rel_gap`.",
    },
    CodeEntry {
        code: "RS0005",
        severity: Severity::Error,
        title: "unknown environment variable",
        explanation: "A `RUSOLVE_*` variable is set that no crate reads. A typo in such a name would silently change nothing (or the engine), so the CLI exits before running. `rusolve env` lists the registered variables; the help line names the closest. Prefer an option (`--opt`) over a variable.",
    },
    CodeEntry {
        code: "RS0006",
        severity: Severity::Error,
        title: "missing or unexpected argument",
        explanation: "A command needs a model path it did not get, or got a positional argument it does not take. `rusolve <command> --help` shows the usage line.",
    },
    CodeEntry {
        code: "RS0010",
        severity: Severity::Error,
        title: "model or input file cannot be read",
        explanation: "The model, solution, basis, options or hook file is missing, unreadable, or its format is not recognised. Check the path and the extension: MPS, LP, .nl and MathOptFormat JSON are read, optionally gzip, bzip2, xz or zstd compressed.",
    },
    CodeEntry {
        code: "RS0011",
        severity: Severity::Error,
        title: "model file does not parse",
        explanation: "The file was read but its content is not a valid model in its format. The message carries the reader's own text; readers that report a file line and column point at it, with the declared names as help when a name is undeclared (RS1004). `rusolve convert` or a re-export from the modelling tool usually fixes a generated file.",
    },
    CodeEntry {
        code: "RS0012",
        severity: Severity::Error,
        title: "model is not supported by this command",
        explanation: "The model parses but belongs to a class the command cannot solve (nonconvex or mixed-integer quadratic, mixed-integer conic, general constraints with quadratic terms, or a Hessian given to `mip`). `rusolve classify <model>` names the class and the engine that would take it.",
    },
    CodeEntry {
        code: "RS0090",
        severity: Severity::Error,
        title: "internal error",
        explanation: "The solver failed in a way that is not the caller's mistake: an unexpected solver status, a result file that cannot be written, or an internal invariant. Re-run with the same input and report it with the model and the command line; the exit code is 4.",
    },
    CodeEntry {
        code: "RS1004",
        severity: Severity::Error,
        title: "name used before it is declared",
        explanation: "A model file refers to a row or column name that no ROWS or COLUMNS section declared (MPS), for example a COLUMNS entry on row `R7` when only `R1`..`R6` exist. The diagnostic points at the line and column of the name and lists the declared names as help. Declare the row, or fix the spelling; the help line suggests the closest declared name.",
    },
    CodeEntry {
        code: "RS1005",
        severity: Severity::Error,
        title: "value is not a number",
        explanation: "A coefficient, right-hand side, range or bound in a model file is not a number the reader can parse (`abc`, `1,5`, `1e`). The diagnostic points at the token. Write decimals with a point and exponents as `1e6`; `inf` and `-inf` are accepted for bounds. The reader reports every such value in the file (up to 20), not only the first.",
    },
    CodeEntry {
        code: "RS1006",
        severity: Severity::Error,
        title: "model file syntax error",
        explanation: "A model file breaks its format's grammar in a way no more specific code covers: an unknown section or bound type, a malformed row or expression, an invalid JSON document. The diagnostic points at the line and column (MathOptFormat: the JSON path) and quotes the reader's message. Structural errors end the list, because what follows would only be a cascade of the first; fix it and read again.",
    },
    CodeEntry {
        code: "RS2001",
        severity: Severity::Warning,
        title: "weak big-M coefficient",
        explanation: "An integer column carries a coefficient far above what the row can ever use (the row's implied activity bound minus its side), so the LP relaxation is weaker than it needs to be and the coefficient is numerically risky. The warning names the row, the column and the tightest valid coefficient and bound. Source: coefficient tightening in Achterberg et al., Presolve reductions in mixed integer programming, INFORMS J. Computing 2020. `rusolve lint` lists every finding; `--no-warnings` silences the warning on solve and mip.",
    },
    CodeEntry {
        code: "RS2002",
        severity: Severity::Warning,
        title: "wide coefficient range",
        explanation: "The constraint matrix mixes very large and very small nonzeros, which costs accuracy in the basis factorisation and the tolerances. The warning names the extremes. Rescale rows and columns (change units) so the nonzeros fall in about 1e-3 to 1e6. `rusolve lint` lists every finding; `--no-warnings` silences the warning.",
    },
    CodeEntry {
        code: "RS2003",
        severity: Severity::Warning,
        title: "huge bound or free column",
        explanation: "A column has a very large finite bound, or is free and appears with no cost in a single row, so it does no work or invites unbounded values. Tighten huge bounds to the largest value the variable can take, or delete a free column singleton. `rusolve lint` lists every finding; `--no-warnings` silences the warning.",
    },
    CodeEntry {
        code: "RS2010",
        severity: Severity::Warning,
        title: "model is infeasible",
        explanation: "The solve ended Infeasible. When cheap (a bounded post-solve budget shared with the relaxation re-solve) the diagnostic names an irreducible infeasible set of rows: dropping any one member makes the rest feasible, so one of them has to change. When the budget ran out first it names an infeasible subset instead: the message says 'a subset of' and a note says 'not reduced'. Those rows conflict, but not every one is needed; `rusolve explain` computes the irreducible set. When the LP relaxation is feasible the conflict involves integrality. Source: Chinneck and Dravnieks 1991, Locating minimal infeasible constraint sets in linear programs. `rusolve explain` verifies the set; `rusolve feasrelax` finds the cheapest relaxation.",
    },
    CodeEntry {
        code: "RS2011",
        severity: Severity::Warning,
        title: "model is unbounded",
        explanation: "The solve ended Unbounded. When cheap the diagnostic names the columns of a primal ray: moving along it keeps every row satisfied and improves the objective without limit. Bound those columns or add the missing constraint. `rusolve explain` verifies the ray and prints the objective rate.",
    },
    CodeEntry {
        code: "RS2012",
        severity: Severity::Warning,
        title: "returned solution violates a constraint",
        explanation: "The solution the solver returned violates a row or a column bound by more than 1e-6. The warning names the worst row and the count. Usually the model is badly scaled; `rusolve lint` shows where, and `rusolve check` re-verifies the solution independently of the solver. `--no-warnings` silences this warning.",
    },
    CodeEntry {
        code: "RS3001",
        severity: Severity::Error,
        title: "MCP tool argument is missing or invalid",
        explanation: "An MCP tool call carried an argument that is absent, of the wrong type or out of range, or named no known tool. The message names the argument path (for example `variables[2].lb`) and what was expected; the tool's inputSchema in tools/list states the full shape. Fix the argument and call again; nothing was changed by the failed call.",
    },
    CodeEntry {
        code: "RS3002",
        severity: Severity::Error,
        title: "MCP model handle is unknown",
        explanation: "The model_id given to an MCP tool names no loaded model, because it was never loaded, was freed with free_model, or the server restarted. Handles live only as long as the server process. Call load_model or build_model again and use the model_id it returns.",
    },
    CodeEntry {
        code: "RS4001",
        severity: Severity::Error,
        title: "invalid argument to a Python model-building call",
        explanation: "A `rusolve.Model` method (add_var, add_constr, set_bounds, ...) got an argument it cannot use: bounds with lb above ub, a vector of the wrong length, an unknown sense, a NaN, a duplicate index. Python raises `rusolve.ModelError` (also a `ValueError`) whose `.diagnostics[0].location` is the line of your code that made the call. The message names the argument; the docstring of the method lists what it accepts. Fix the argument and call again; the model is unchanged when a call raises.",
    },
    CodeEntry {
        code: "RS4002",
        severity: Severity::Error,
        title: "row or column index out of range in a Python call",
        explanation: "An index passed to a `rusolve.Model`, `Vars` or `Solution` call is below zero or at least the number of rows or columns added so far. Python raises `rusolve.ModelIndexError` (also an `IndexError`). Indices count from 0 in the order the rows or columns were added; `m.num_rows` and `m.num_cols` give the limits. A row index from before a failed call may be stale because a failing `add_constrs` removes the rows it had added.",
    },
    CodeEntry {
        code: "RS4003",
        severity: Severity::Error,
        title: "variables of two different Models combined",
        explanation: "A Var, LinExpr or Constraint made from one `rusolve.Model` was used with another Model, or a Solution of one Model was indexed with a Var of another. Indices are only meaningful inside their own Model, so mixing them would silently read the wrong column. Create the expression from the Model it is added to, or index the Solution with variables of the Model that produced it. Python raises `rusolve.ModelError`.",
    },
    CodeEntry {
        code: "RS4004",
        severity: Severity::Error,
        title: "the model is infeasible (Python)",
        explanation: "`Solution.raise_for_status()` found status Infeasible: no assignment satisfies every row and bound. Python raises `rusolve.InfeasibleError`. `Model.explain()` returns an irreducible infeasible subsystem (a small set of rows and bounds that conflict), with the Python line that added each row when `Model.capture_source` is on. `Model.feasrelax()` finds the cheapest relaxation that makes the model feasible.",
    },
    CodeEntry {
        code: "RS4005",
        severity: Severity::Error,
        title: "the model is unbounded or infeasible-or-unbounded (Python)",
        explanation: "`Solution.raise_for_status()` found status Unbounded or InfeasibleOrUnbounded: the objective can improve without limit, or presolve could not tell that from infeasibility. Python raises `rusolve.UnboundedError`. `Model.explain()` returns the unbounded ray or the infeasible subsystem. Usually a bound or a constraint is missing on the variables the ray moves along; re-solving with `presolve=False` separates the two cases.",
    },
    CodeEntry {
        code: "RS4006",
        severity: Severity::Error,
        title: "the solve stopped without an optimal solution (Python)",
        explanation: "`Solution.raise_for_status()` found a status other than Optimal, Infeasible or Unbounded: a time, iteration, node or solution limit, an interrupt, or an internal error. Python raises `rusolve.LimitError`. When `Solution.ok` is true the best primal solution found so far is in `Solution.x` with its gap in `Solution.gap`. Raise `time_limit` (or the named limit) to continue, or relax the gap.",
    },
    CodeEntry {
        code: "RS4007",
        severity: Severity::Error,
        title: "no primal solution to read (Python)",
        explanation: "A value was read from a `Solution` that has no primal vector, for example `sol[var]` after an infeasible solve (`sol.x` and `sol.duals` are None when there is no solution, they do not raise). Python raises `rusolve.NoSolutionError` (also a `ValueError`). Check `Solution.ok` before reading values, or call `Solution.raise_for_status()` to get the typed error for the status.",
    },
    CodeEntry {
        code: "RS4008",
        severity: Severity::Error,
        title: "a Python callback hook raised an exception",
        explanation: "A hook of the object passed as `Model.solve(callback=...)` (on_progress, on_improving_solution, on_cut_pool, user_solution, separate, branching_score, select_node, lazy_constraints) raised, or returned a value of the wrong type or length. The solver never sees the exception: the first one is kept, no later hook is called (each falls back to the solver's own rule), the solve is interrupted at its next limit check and `Model.solve` raises `rusolve.CallbackError` whose `__cause__` is the original exception and whose `.solution` is the partial result. `Model.callback_stats` counts calls and errors per hook. Fix the hook: it must return what its docstring names (a list of cuts, one score per candidate, an index into `ctx.open_nodes`, a full-length solution or None).",
    },
];

/// The registry entry of `code` (case-insensitive, `rs1004` works).
pub fn lookup(code: &str) -> Option<&'static CodeEntry> {
    REGISTRY.iter().find(|e| e.code.eq_ignore_ascii_case(code))
}

/// The text `rusolve --explain RSnnnn` prints.
pub fn explain_text(entry: &CodeEntry) -> String {
    format!(
        "{}: {} ({})\n\n{}\n",
        entry.code,
        entry.title,
        entry.severity.as_str(),
        entry.explanation
    )
}

/// The registry as a JSON array (`rusolve --explain all --format json`), which
/// `tools/gen_agent_docs.py` renders into `docs/agent/diagnostics.md`.
pub fn registry_json() -> String {
    let items: Vec<String> = REGISTRY
        .iter()
        .map(|e| {
            format!(
                "{{\"code\":{},\"severity\":{},\"title\":{},\"explanation\":{}}}",
                json_str(e.code),
                json_str(e.severity.as_str()),
                json_str(e.title),
                json_str(e.explanation)
            )
        })
        .collect();
    format!("[{}]", items.join(","))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};

    fn rs_files(dir: &Path, out: &mut Vec<PathBuf>) {
        for e in std::fs::read_dir(dir).unwrap().flatten() {
            let p = e.path();
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            if p.is_dir() {
                if name != "target" && name != "fixtures" {
                    rs_files(&p, out);
                }
            } else if name.ends_with(".rs") {
                out.push(p);
            }
        }
    }

    /// Every `"RSdddd"` literal in the workspace sources.
    fn codes_in_sources() -> Vec<(String, PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let mut files = Vec::new();
        rs_files(&root, &mut files);
        let mut found = Vec::new();
        for f in files {
            let text = std::fs::read_to_string(&f).unwrap_or_default();
            let b = text.as_bytes();
            for i in 0..b.len().saturating_sub(7) {
                if b[i] == b'"'
                    && &b[i + 1..i + 3] == b"RS"
                    && b[i + 3..i + 7].iter().all(u8::is_ascii_digit)
                    && b[i + 7] == b'"'
                {
                    found.push((text[i + 1..i + 7].to_string(), f.clone()));
                }
            }
        }
        found
    }

    /// Codes named at an emit site: `"RSdddd"` literals in non-test source, outside the
    /// registry's own `code:` lines. A file's inline test module (from `#[cfg(test)] mod`) and
    /// `tests/` directories do not count: a code only a test emits is dead.
    fn emitted_codes(files: &[(PathBuf, String)]) -> std::collections::HashSet<String> {
        let mut out = std::collections::HashSet::new();
        for (path, text) in files {
            if path.components().any(|c| c.as_os_str() == "tests") {
                continue;
            }
            let body = text.split("#[cfg(test)]\nmod ").next().unwrap_or("");
            for line in body.lines() {
                if line.trim_start().starts_with("code: \"RS") {
                    continue;
                }
                let b = line.as_bytes();
                for i in 0..b.len().saturating_sub(7) {
                    if b[i] == b'"'
                        && &b[i + 1..i + 3] == b"RS"
                        && b[i + 3..i + 7].iter().all(u8::is_ascii_digit)
                        && b[i + 7] == b'"'
                    {
                        out.insert(line[i + 1..i + 7].to_string());
                    }
                }
            }
        }
        out
    }

    fn workspace_sources() -> Vec<(PathBuf, String)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let mut files = Vec::new();
        rs_files(&root, &mut files);
        files
            .into_iter()
            .map(|f| {
                let t = std::fs::read_to_string(&f).unwrap_or_default();
                (f, t)
            })
            .collect()
    }

    #[test]
    fn every_registered_code_is_emitted_somewhere() {
        let emitted = emitted_codes(&workspace_sources());
        for e in REGISTRY {
            assert!(
                emitted.contains(e.code),
                "{} is in diagnostic::REGISTRY but no emit site names it: emit it or remove it",
                e.code
            );
        }
    }

    /// The scan fails for a code nothing emits (what a dummy registry entry would be).
    #[test]
    fn the_emitted_scan_misses_a_dummy_code() {
        let emitted = emitted_codes(&workspace_sources());
        assert!(!emitted.contains(concat!("RS99", "99")));
        // a registry-style line and a test-module literal are not emit sites
        let fake = vec![(
            PathBuf::from("x.rs"),
            "code: \"RS9998\",\n#[cfg(test)]\nmod t { fn f() { let _ = \"RS9997\"; } }".to_string(),
        )];
        assert!(emitted_codes(&fake).is_empty());
    }

    #[test]
    fn every_emitted_code_is_registered() {
        for (code, file) in codes_in_sources() {
            assert!(
                REGISTRY.iter().any(|e| e.code == code),
                "{code} is named in {} but is not in diagnostic::REGISTRY",
                file.display()
            );
        }
    }

    #[test]
    fn every_registered_code_is_explained_and_documented() {
        let page = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/agent/diagnostics.md"),
        )
        .expect("docs/agent/diagnostics.md: run python3 tools/gen_agent_docs.py");
        let mut seen = std::collections::HashSet::new();
        for e in REGISTRY {
            assert!(seen.insert(e.code), "{} registered twice", e.code);
            assert!(
                e.code.len() == 6 && e.code.starts_with("RS"),
                "bad code {}",
                e.code
            );
            assert!(!e.title.is_empty(), "{} has no title", e.code);
            assert!(
                e.explanation.len() >= 60,
                "{} explanation too short",
                e.code
            );
            assert!(
                page.contains(&format!("### {}", e.code)),
                "{} is registered but missing from docs/agent/diagnostics.md",
                e.code
            );
        }
    }

    #[test]
    fn text_renderer_snapshot() {
        let d = Diagnostic::error("RS1004", "row 'R7' is not declared")
            .with_location(Location::File {
                path: "m.mps".into(),
                line: 12,
                column: 15,
                len: 2,
                source_line: "    X1        R7        1.0".into(),
            })
            .with_note("ROWS declared 3 rows")
            .with_help("did you mean 'R1'?");
        assert_eq!(
            d.render_text(),
            "error[RS1004]: row 'R7' is not declared\n\
             \x20 --> m.mps:12:15\n\
             \x20  |\n\
             12 |     X1        R7        1.0\n\
             \x20  |               ^^\n\
             \x20  = note: ROWS declared 3 rows\n\
             \x20  = help: did you mean 'R1'?\n"
        );
        let r = Diagnostic::warning("RS0090", "big-M")
            .with_location(Location::Row {
                name: "c1".into(),
                rendered: "c1: 1e6 x + y <= 4".into(),
            })
            .with_help("use an indicator");
        assert_eq!(
            r.render_text(),
            "warning[RS0090]: big-M\n --> row 'c1'\n  = row: c1: 1e6 x + y <= 4\n  = help: use an indicator\n"
        );
    }

    #[test]
    fn json_escapes_and_envelope() {
        let d = Diagnostic::error("RS0002", "bad \"flag\"\n").with_help("h");
        let j = diagnostics_json(&[d], 2);
        assert_eq!(
            j,
            "{\"schema_version\":\"rusolve.diagnostics/1\",\"exit_code\":2,\"diagnostics\":[{\"code\":\"RS0002\",\"severity\":\"error\",\"message\":\"bad \\\"flag\\\"\\n\",\"notes\":[],\"help\":[\"h\"]}]}"
        );
    }

    #[test]
    fn did_you_mean_finds_typos() {
        let flags = ["--time-limit", "--iter-limit", "--format", "--opt"];
        assert_eq!(
            did_you_mean("--time-limt", flags.iter().copied()).as_deref(),
            Some("did you mean '--time-limit'?")
        );
        assert_eq!(did_you_mean("--zzzzzz", flags.iter().copied()), None);
        let opts = ["mip_rel_gap", "mip_abs_gap", "time_limit"];
        assert_eq!(
            suggestions("mip_rel_gapp", opts.iter().copied(), 3)[0],
            "mip_rel_gap"
        );
        assert_eq!(edit_distance("slove", "solve"), 1);
    }

    #[test]
    fn explain_finds_codes_case_insensitively() {
        assert!(lookup("rs1004").is_some());
        assert!(lookup(concat!("RS99", "99")).is_none());
        assert!(explain_text(lookup("RS0002").unwrap()).starts_with("RS0002: unknown flag"));
    }
}
