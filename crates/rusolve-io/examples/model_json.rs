//! `cargo run --release -p rusolve-io --example model_json -- <model> [<lp-out>]`
//!
//! Prints rusolve's reading of `<model>` (MPS or LP, optionally `.gz`) as
//! one JSON object, with the names `write_lp` gives every column and row. With
//! `<lp-out>`, also writes the LP text there and adds the write time and the
//! writer's report, so an external script can compare another reader's
//! reading of the LP text with rusolve's in-memory model. Non-finite numbers
//! are `Infinity` / `-Infinity` / `NaN` (accepted by Python's `json`).

use std::fmt::Write as _;
use std::time::Instant;

fn num(out: &mut String, v: f64) {
    if v.is_nan() {
        out.push_str("NaN");
    } else if v == f64::INFINITY {
        out.push_str("Infinity");
    } else if v == f64::NEG_INFINITY {
        out.push_str("-Infinity");
    } else {
        let _ = write!(out, "{v:?}");
    }
}

fn string(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn strings<'a>(out: &mut String, it: impl Iterator<Item = &'a str>) {
    out.push('[');
    for (k, s) in it.enumerate() {
        if k > 0 {
            out.push(',');
        }
        string(out, s);
    }
    out.push(']');
}

fn nums(out: &mut String, it: impl Iterator<Item = f64>) {
    out.push('[');
    for (k, v) in it.enumerate() {
        if k > 0 {
            out.push(',');
        }
        num(out, v);
    }
    out.push(']');
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(path) = args.first() else {
        eprintln!("usage: model_json <model> [<lp-out>]");
        std::process::exit(2);
    };
    let m = match rusolve_io::read_model_file(path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("read {path}: {e}");
            std::process::exit(3);
        }
    };
    let names = rusolve_io::lp_names(&m);
    let mut o = String::with_capacity(1 << 20);
    o.push('{');
    o.push_str("\"name\":");
    string(&mut o, &m.name);
    o.push_str(",\"maximize\":");
    o.push_str(if m.obj_sense == rusolve_io::ObjSense::Maximize {
        "true"
    } else {
        "false"
    });
    o.push_str(",\"obj_name\":");
    match &m.obj_name {
        Some(s) => string(&mut o, s),
        None => o.push_str("null"),
    }
    o.push_str(",\"lp_obj_name\":");
    match &names.objective {
        Some(s) => string(&mut o, s),
        None => o.push_str("null"),
    }
    o.push_str(",\"obj_offset\":");
    num(&mut o, m.obj_offset);
    o.push_str(",\"col_names\":");
    strings(&mut o, m.col_names.iter().map(String::as_str));
    o.push_str(",\"lp_col_names\":");
    strings(&mut o, names.columns.iter().map(String::as_str));
    o.push_str(",\"row_names\":");
    strings(&mut o, m.row_names.iter().map(String::as_str));
    o.push_str(",\"lp_row_names\":");
    strings(&mut o, names.rows.iter().map(String::as_str));
    o.push_str(",\"lp_range_columns\":");
    strings(
        &mut o,
        names
            .range_columns
            .iter()
            .map(|r| r.as_deref().unwrap_or("")),
    );
    o.push_str(",\"objective\":");
    nums(&mut o, m.objective.iter().copied());
    o.push_str(",\"lower\":");
    nums(&mut o, m.bounds.iter().map(|b| b.lower));
    o.push_str(",\"upper\":");
    nums(&mut o, m.bounds.iter().map(|b| b.upper));
    o.push_str(",\"integer\":[");
    for (k, b) in m.bounds.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        o.push_str(if b.integer { "1" } else { "0" });
    }
    o.push_str("],\"senses\":\"");
    for s in &m.senses {
        o.push(match s {
            rusolve_io::RowSense::LessEqual => 'L',
            rusolve_io::RowSense::GreaterEqual => 'G',
            rusolve_io::RowSense::Equal => 'E',
        });
    }
    o.push_str("\",\"rhs\":");
    nums(&mut o, m.rhs.iter().copied());
    o.push_str(",\"ranges\":");
    nums(&mut o, m.ranges.iter().copied());
    let a = &m.matrix;
    o.push_str(",\"col_ptr\":[");
    for (k, p) in a.col_ptr().iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        let _ = write!(o, "{p}");
    }
    o.push_str("],\"row_idx\":[");
    for (k, r) in a.row_idx().iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        let _ = write!(o, "{r}");
    }
    o.push_str("],\"values\":");
    nums(&mut o, a.values().iter().copied());
    o.push_str(",\"hessian\":[");
    for (k, (r, c, v)) in m.hessian.iter().enumerate() {
        if k > 0 {
            o.push(',');
        }
        let _ = write!(o, "[{r},{c},");
        num(&mut o, *v);
        o.push(']');
    }
    o.push(']');
    if let Some(out) = args.get(1) {
        let t = Instant::now();
        let (text, r) = rusolve_io::write_lp_report(&m);
        let lp_secs = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let mps = rusolve_io::write_mps(&m);
        let mps_secs = t.elapsed().as_secs_f64();
        if let Err(e) = std::fs::write(out, &text) {
            eprintln!("write {out}: {e}");
            std::process::exit(4);
        }
        let _ = write!(
            o,
            ",\"write\":{{\"lp_bytes\":{},\"lp_secs\":{lp_secs},\"mps_bytes\":{},\"mps_secs\":{mps_secs},\
             \"renamed_columns\":{},\"renamed_rows\":{},\"renamed_objective\":{},\"ranged_rows\":{},\
             \"cplex_only_numbers\":{},\"unreadable_numbers\":{}}}",
            text.len(),
            mps.len(),
            r.renamed_columns,
            r.renamed_rows,
            r.renamed_objective,
            r.ranged_rows,
            r.cplex_only_numbers,
            r.unreadable_numbers
        );
    }
    o.push('}');
    println!("{o}");
}
