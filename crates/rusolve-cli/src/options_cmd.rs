//! `rusolve options` — the option registry as an agent sees it.
//!
//! The same table that validates `--opt` and options files is printed here, so what an agent
//! reads is exactly what the solver accepts: name, current value, type, range, default,
//! whether the option can move the search, and where it comes from (HiGHS 1.16.0, or the
//! rusolve change that added a lever).
//!
//! ```text
//! rusolve options                     # every option, current values
//! rusolve options --non-default       # only what differs from the registry default
//! rusolve options --profile default   # the HiGHS 1.15.1 profile instead of the best-measured one
//! rusolve options --json              # machine-readable, one entry per option
//! ```

use rusolve_ekk::lp_data::highs_option_records as registry;
use rusolve_ekk::lp_data::highs_options::HighsOptions;
use std::io::Write;
use std::process::ExitCode;

const USAGE: &str =
    "usage: rusolve options [--json] [--profile default|recommended] [--non-default]";

pub fn run_options(args: &[String]) -> ExitCode {
    let wj = crate::diag::wants_json(args);
    let usage = |msg: &str| {
        crate::diag::fail_usage_text(
            wj || args.iter().any(|a| a == "--json"),
            "options",
            &format!("{msg}\n{USAGE}"),
        )
    };
    let mut json = false;
    let mut non_default = false;
    let mut profile: Option<String> = None;
    let mut rest = args.iter();
    while let Some(flag) = rest.next() {
        match flag.as_str() {
            "--json" => json = true,
            "--non-default" => non_default = true,
            "--profile" => {
                let Some(value) = rest.next() else {
                    return usage("--profile needs default or recommended");
                };
                if value != "default" && value != "recommended" {
                    return usage(&format!("--profile {value}: needs default or recommended"));
                }
                profile = Some(value.clone());
            }
            "--help" | "-h" => {
                println!("{USAGE}");
                return ExitCode::SUCCESS;
            }
            other => {
                if !other.starts_with('-') {
                    return usage(&format!("unexpected extra argument '{other}'"));
                }
                return usage(&format!("flag '{other}' is not supported"));
            }
        }
    }
    // Without --profile, the same profile the solver itself would start from, so the listing
    // describes the run an agent is about to do.
    let options: HighsOptions = match profile.as_deref() {
        Some("default") => HighsOptions::default(),
        Some("recommended") => HighsOptions::recommended(),
        _ => crate::base_options(),
    };
    let text = if json {
        match registry::option_records_json(&options, Some(non_default)) {
            Ok(text) => text + "\n",
            Err(err) => {
                return crate::diag::fail_msg(json, "RS0090", 4, format!("rusolve options: {err}"));
            }
        }
    } else {
        registry::report_options_table(&options, Some(non_default))
    };
    // Listing into `head` is the normal way to read this: a closed pipe is success.
    let stdout = std::io::stdout();
    let mut stdout = stdout.lock();
    if let Err(err) = stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.flush())
        && err.kind() != std::io::ErrorKind::BrokenPipe
    {
        return crate::diag::fail_msg(json, "RS0090", 4, format!("rusolve options: {err}"));
    }
    ExitCode::SUCCESS
}
