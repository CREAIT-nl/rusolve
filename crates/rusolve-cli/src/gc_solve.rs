//! Solve path for models with general constraints (indicator, SOS, PWL,
//! abs/min/max/and/or/not): linearize into the port's `HighsLp`, solve, cut
//! the solution back to the user's columns and rows, and apply the
//! trickle-flow repair.
//!
//! A model without general constraints never reaches this module
//! (`mip_highs::run_mip_highs` checks `GcModel::has_general`), so its path is
//! exactly the port's.

use rusolve_ekk::lp_data::hconst::{HighsVarType, MatrixFormat};
use rusolve_ekk::lp_data::highs_lp::HighsLp;
use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_io::gc_lin::{self, Extension, GcError, LinInput, LinOptions};
use rusolve_io::{GcModel, GcTol};
use serde_json::{Value, json};

/// What `prepare` keeps for the solve and the report.
pub struct GcPlan {
    /// The user's model, before any derived column, row or bound change.
    pub user_lp: HighsLp,
    pub gc: GcModel,
    pub ext: Extension,
}

/// Linearizes `gc` into `lp` (appending columns and rows, tightening the
/// domain of PWL operands). Returns the plan holding the user's model.
pub fn prepare(lp: &mut HighsLp, gc: GcModel, bigm_max: f64) -> Result<GcPlan, GcError> {
    let user_lp = lp.clone();
    let n = lp.num_col_ as usize;
    let m = lp.num_row_ as usize;
    let col_start: Vec<usize> = lp.a_matrix_.start_[..=n]
        .iter()
        .map(|&v| v as usize)
        .collect();
    let nnz = col_start[n];
    let col_index: Vec<usize> = lp.a_matrix_.index_[..nnz]
        .iter()
        .map(|&v| v as usize)
        .collect();
    let (row_start, row_index, row_value) =
        gc_lin::csr_from_csc(m, &col_start, &col_index, &lp.a_matrix_.value_[..nnz]);
    let names: Vec<String> = if lp.col_names_.len() == n {
        lp.col_names_.clone()
    } else {
        Vec::new()
    };
    let ext = gc_lin::linearize(
        &gc,
        &LinInput {
            col_lower: &lp.col_lower_,
            col_upper: &lp.col_upper_,
            col_cost: &lp.col_cost_,
            maximize: lp.sense_ == rusolve_ekk::lp_data::hconst::ObjSense::kMaximize,
            row_start: &row_start,
            row_index: &row_index,
            row_value: &row_value,
            row_lower: &lp.row_lower_,
            row_upper: &lp.row_upper_,
            col_names: &names,
        },
        LinOptions {
            bigm_max,
            ..LinOptions::default()
        },
    )?;
    extend_lp(lp, &ext);
    Ok(GcPlan { user_lp, gc, ext })
}

/// Appends the extension's columns and rows to `lp` (column-wise matrix).
fn extend_lp(lp: &mut HighsLp, ext: &Extension) {
    let n = lp.num_col_ as usize;
    let m = lp.num_row_ as usize;
    for &(j, lo, hi) in &ext.bound_changes {
        lp.col_lower_[j] = lp.col_lower_[j].max(lo);
        lp.col_upper_[j] = lp.col_upper_[j].min(hi);
    }
    if ext.cols.is_empty() && ext.rows.is_empty() {
        return;
    }
    debug_assert_eq!(lp.a_matrix_.format_, MatrixFormat::kColwise);
    let n2 = n + ext.cols.len();
    let mut extra: Vec<Vec<(i32, f64)>> = vec![Vec::new(); n2];
    for (r, row) in ext.rows.iter().enumerate() {
        for &(j, c) in &row.terms {
            extra[j].push(((m + r) as i32, c));
        }
    }
    let old_start = std::mem::take(&mut lp.a_matrix_.start_);
    let old_index = std::mem::take(&mut lp.a_matrix_.index_);
    let old_value = std::mem::take(&mut lp.a_matrix_.value_);
    let mut start = Vec::with_capacity(n2 + 1);
    let mut index = Vec::with_capacity(old_index.len() + ext.rows.len() * 3);
    let mut value = Vec::with_capacity(old_index.len() + ext.rows.len() * 3);
    start.push(0);
    for (j, e) in extra.iter().enumerate() {
        if j < n {
            let (a, b) = (old_start[j] as usize, old_start[j + 1] as usize);
            index.extend_from_slice(&old_index[a..b]);
            value.extend_from_slice(&old_value[a..b]);
        }
        for &(i, c) in e {
            index.push(i);
            value.push(c);
        }
        start.push(index.len() as i32);
    }
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = index;
    lp.a_matrix_.value_ = value;
    lp.a_matrix_.p_end_.clear();
    lp.a_matrix_.num_col_ = n2 as i32;
    lp.a_matrix_.num_row_ = (m + ext.rows.len()) as i32;
    lp.num_col_ = n2 as i32;
    lp.num_row_ = (m + ext.rows.len()) as i32;
    let had_names = lp.col_names_.len() == n;
    let had_rows = lp.row_names_.len() == m;
    if ext.cols.iter().any(|c| c.integer) && lp.integrality_.is_empty() {
        lp.integrality_ = vec![HighsVarType::kContinuous; n];
    }
    for c in &ext.cols {
        lp.col_cost_.push(0.0);
        lp.col_lower_.push(c.lower);
        lp.col_upper_.push(c.upper);
        if !lp.integrality_.is_empty() {
            lp.integrality_.push(if c.integer {
                HighsVarType::kInteger
            } else {
                HighsVarType::kContinuous
            });
        }
        if had_names {
            lp.col_names_.push(format!("gc_{}", c.name));
        }
    }
    for r in &ext.rows {
        lp.row_lower_.push(r.lower);
        lp.row_upper_.push(r.upper);
        if had_rows {
            lp.row_names_.push(format!("gc_{}", r.name));
        }
    }
}

/// `routes` of the envelope's `general_constraints` object.
pub fn routes_json(ext: &Extension) -> Vec<Value> {
    ext.routes
        .iter()
        .map(|r| {
            json!({
                "object": r.object,
                "route": r.route,
                "max_big_m": r.max_big_m,
            })
        })
        .collect()
}

/// The outcome of [`finish`].
pub struct GcOutcome {
    pub report: Value,
}

/// After the solve: check the incumbent against the user's objects, apply the
/// trickle-flow repair once, and cut the solution to the user's size.
///
/// Trickle flow: a binary that is `1e-7` off its value lets a big-M row hold a
/// little slack, and the user's indicator then fails on the rounded binary.
/// The repair fixes the binary of each violated indicator at its rounded value
/// and solves once more; the event is reported either way.
pub fn finish(plan: &GcPlan, highs: &mut Highs, tol: GcTol) -> GcOutcome {
    let n = plan.user_lp.num_col_ as usize;
    let m = plan.user_lp.num_row_ as usize;
    let mut trickle = None;
    let check = |h: &Highs| -> Vec<rusolve_io::GcViolation> {
        plan.gc.check(
            &h.solution_.col_value[..n],
            &plan.user_lp.col_lower_,
            &plan.user_lp.col_upper_,
            &plan.user_lp.col_names_,
            tol,
        )
    };
    let usable = |h: &Highs| h.solution_.value_valid && h.solution_.col_value.len() >= n;
    if usable(highs) {
        let before = check(highs);
        let fix: Vec<(usize, &str)> = before
            .iter()
            .filter(|v| v.kind == "indicator")
            .filter_map(|v| {
                plan.gc
                    .indicators
                    .iter()
                    .find(|o| o.name == v.name)
                    .map(|o| (o.bin, o.name.as_str()))
            })
            .collect();
        if !fix.is_empty() {
            let mut names = Vec::new();
            for &(bin, name) in &fix {
                let r = highs.solution_.col_value[bin].round().clamp(0.0, 1.0);
                highs.change_col_bounds(bin as i32, r, r);
                names.push(name.to_string());
            }
            // the repaired solve starts from scratch: the incumbent is the
            // very point being repaired
            highs.solution_ = Default::default();
            highs.run();
            let after = if usable(highs) { check(highs).len() } else { 0 };
            trickle = Some(json!({
                "fixed_indicators": names,
                "violations_before": before.len(),
                "violations_after": after,
                "status": format!("{:?}", highs.get_model_status()),
            }));
        }
    }
    let mut violations = Vec::new();
    if usable(highs) {
        for v in check(highs) {
            violations.push(
                json!({"kind": v.kind, "name": v.name, "amount": v.amount, "detail": v.detail}),
            );
        }
    }
    // cut the solution to the user's columns and rows
    let s = &mut highs.solution_;
    s.col_value.truncate(n);
    s.row_value.truncate(m);
    s.col_dual.truncate(n);
    s.row_dual.truncate(m);
    let mut report = json!({
        "objects": plan.gc.counts().iter().filter(|(_, c)| *c > 0).map(|(k, c)| (k.to_string(), json!(c))).collect::<serde_json::Map<_, _>>(),
        "routes": routes_json(&plan.ext),
        "added_columns": plan.ext.cols.len(),
        "added_rows": plan.ext.rows.len(),
        "violations": violations,
        "verified": violations.is_empty(),
    });
    if let Some(t) = trickle {
        report["trickle_flow"] = t;
    }
    GcOutcome { report }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusolve_ekk::lp_data::highs_options::HighsOptions;

    /// Trickle flow: a binary at 1e-7 rounds to 0 in the checker, so the
    /// `b = 0 -> x >= 5` indicator is active while the big-M row let `x = 1`
    /// through. `finish` must fix `b = 0` and solve once more.
    #[test]
    fn trickle_flow_fixes_the_binary_and_resolves() {
        let dir = std::env::temp_dir().join(format!("rusolve-gc-trickle-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.lp");
        std::fs::write(
            &path,
            "Minimize\n obj: x + 3 b\nSubject To\n c: x + b >= 0\n i: b = 0 -> x >= 5\nBounds\n 0 <= x <= 10\nBinaries\n b\nEnd\n",
        )
        .unwrap();
        let (mut lp, _, gc) = crate::mip_highs::read_mip_full(path.to_str().unwrap()).unwrap();
        let plan = prepare(&mut lp, gc, 1e6).unwrap();
        let mut highs = Highs::with_options(HighsOptions {
            output_flag: false,
            threads: 1,
            ..HighsOptions::default()
        });
        highs.pass_model(lp);
        highs.run();
        // tamper with the incumbent: x = 1, b = 1e-7
        let (xi, bi) = (0, 1);
        highs.solution_.col_value[xi] = 1.0;
        highs.solution_.col_value[bi] = 1e-7;
        let out = finish(&plan, &mut highs, GcTol::default());
        let t = &out.report["trickle_flow"];
        assert_eq!(t["fixed_indicators"][0], "i", "{}", out.report);
        assert_eq!(t["violations_before"], 1);
        assert_eq!(t["violations_after"], 0);
        // b fixed at 0 forces x >= 5: the re-solved incumbent is x = 5, b = 0
        assert!((highs.solution_.col_value[xi] - 5.0).abs() < 1e-6);
        assert_eq!(
            highs.solution_.col_value.len(),
            2,
            "cut to the user's columns"
        );
        assert_eq!(out.report["verified"], true);
    }
}
