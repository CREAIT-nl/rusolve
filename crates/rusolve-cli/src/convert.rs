//! `rusolve convert <in> <out>`: read a model (MPS or LP, optionally `.gz`)
//! and write it in the format the output extension names: `.lp`, `.mps`,
//! `.lp.gz`, `.mps.gz`.
//!
//! The LP output is the CPLEX dialect of `rusolve_io::write_lp` (columns,
//! rows, names, values and types read back identically by LP-format readers;
//! ranged rows carry an `Rg<row>` column).
//! Text by default; `--format json` prints one `rusolve.convert/1` object
//! with what the writer renamed or could not keep exact.
//!
//! Exit codes: 2 for a usage error (missing positional, unknown output
//! extension), 3 when the input cannot be read, 4 when the output cannot be
//! written.

use std::process::ExitCode;

use serde_json::json;

use crate::commands::{Stack, parse_flags};

/// `Some(true)` for an LP output path, `Some(false)` for MPS, `None` otherwise.
fn output_is_lp(path: &str) -> Option<bool> {
    let lower = path.to_ascii_lowercase();
    let stem = lower.strip_suffix(".gz").unwrap_or(&lower);
    if stem.ends_with(".lp") {
        Some(true)
    } else if stem.ends_with(".mps") {
        Some(false)
    } else {
        None
    }
}

pub fn run(args: &[String]) -> ExitCode {
    let parsed = match parse_flags("convert", Stack::Both, args) {
        Ok(p) => p,
        Err(e) => return crate::diag::fail_usage(crate::diag::wants_json(args), &e),
    };
    let json = parsed.get_str("--format") == Some("json");
    let [input, output] = parsed.positionals() else {
        return crate::diag::fail_msg(
            json,
            "RS0006",
            2,
            "rusolve convert: expected <model> <output> (see `rusolve convert --help`)",
        );
    };
    let Some(lp) = output_is_lp(output) else {
        return crate::diag::fail_msg(
            json,
            "RS0003",
            2,
            format!(
                "rusolve convert: {output}: unknown output format (expected .lp, .mps, .lp.gz or .mps.gz)"
            ),
        );
    };
    let m = match rusolve_io::read_model_file_conic(input) {
        Ok(m) => m,
        Err(e) => {
            return crate::diag::fail_read(json, input, &e.to_string());
        }
    };
    let report = match rusolve_io::write_model_file_report(output, &m) {
        Ok(r) => r.unwrap_or_default(),
        Err(e) => {
            return crate::diag::fail_msg(
                json,
                "RS0090",
                4,
                format!("rusolve convert: write {output}: {e}"),
            );
        }
    };
    if json {
        let v = json!({
            "schema_version": "rusolve.convert/1",
            "input": input,
            "output": output,
            "format": if lp { "lp" } else { "mps" },
            "gzip": output.to_ascii_lowercase().ends_with(".gz"),
            "columns": m.num_variables,
            "rows": m.num_constraints,
            "integers": m.bounds.iter().filter(|b| b.integer).count(),
            "range_columns": report.ranged_rows,
            "renamed_columns": report.renamed_columns,
            "renamed_rows": report.renamed_rows,
            "renamed_objective": report.renamed_objective,
            "cplex_only_numbers": report.cplex_only_numbers,
            "unreadable_numbers": report.unreadable_numbers,
        });
        println!("{v}");
    } else {
        println!(
            "wrote {output} ({}): {} columns, {} rows",
            if lp { "LP" } else { "MPS" },
            m.num_variables,
            m.num_constraints
        );
        if lp {
            println!(
                "range columns {}, renamed columns {}, renamed rows {}{}, cplex-only numbers {}, unreadable numbers {}",
                report.ranged_rows,
                report.renamed_columns,
                report.renamed_rows,
                if report.renamed_objective {
                    " (objective renamed)"
                } else {
                    ""
                },
                report.cplex_only_numbers,
                report.unreadable_numbers
            );
        }
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_format_follows_the_extension() {
        assert_eq!(output_is_lp("a.lp"), Some(true));
        assert_eq!(output_is_lp("a.LP.gz"), Some(true));
        assert_eq!(output_is_lp("a.mps.gz"), Some(false));
        assert_eq!(output_is_lp("a.txt"), None);
        assert_eq!(output_is_lp("a.gz"), None);
    }
}
