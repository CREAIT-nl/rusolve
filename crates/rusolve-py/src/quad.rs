//! Quadratic objective input for `Model.set_quadratic_objective` (convex QP). Kept apart from `model.rs` so the LP/MIP front end is
//! untouched: a model that never calls it never reaches this file.

use std::collections::BTreeMap;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rusolve_core::qp_convexity::{Convexity, Method, Verdict, analyze, class_code};
use rusolve_ekk::HighsInt;
use rusolve_ekk::lp_data::hconst::HighsVarType;
use rusolve_ekk::prelude::{HessianFormat, HighsHessian, HighsLp};

use crate::ingest::{F64s, Ints};
use crate::model::{MatKind, read_matrix};

/// Read a quadratic objective into the full symmetric entry list `(row, col, value)` with both
/// triangles stored, zeros dropped, duplicates summed, and `(Q + Q')/2` for an asymmetric Q.
/// `None` when no input was given (the Hessian is cleared). `n` is the current variable count.
pub(crate) fn read_hessian(
    q: Option<&Bound<'_, PyAny>>,
    rows: Option<&Bound<'_, PyAny>>,
    cols: Option<&Bound<'_, PyAny>>,
    vals: Option<&Bound<'_, PyAny>>,
    model: u64,
    n: usize,
) -> PyResult<Option<Vec<(usize, usize, f64)>>> {
    let q = q.filter(|q| !q.is_none());
    let triplet = rows.is_some() || cols.is_some() || vals.is_some();
    if q.is_some() && triplet {
        return Err(crate::diag::model_err(
            "set_quadratic_objective: give Q or rows/cols/vals, not both",
        ));
    }
    // raw[(i, j)] = sum of the entries given at (i, j); the Hessian is (raw + raw') / 2.
    let mut raw: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    let mut put = |i: i64, j: i64, v: f64| -> PyResult<()> {
        if i < 0 || j < 0 || i as u64 >= n as u64 || j as u64 >= n as u64 {
            return Err(crate::diag::index_err(format!(
                "quadratic objective: entry ({i}, {j}) is out of range (there are {n} variables)"
            )));
        }
        if !v.is_finite() {
            return Err(crate::diag::model_err(format!(
                "quadratic objective: entry ({i}, {j}) is {v}; values must be finite"
            )));
        }
        if v != 0.0 {
            *raw.entry((i as usize, j as usize)).or_insert(0.0) += v;
        }
        Ok(())
    };
    if triplet {
        let (Some(r), Some(c), Some(v)) = (rows, cols, vals) else {
            return Err(crate::diag::model_err(
                "set_quadratic_objective: rows, cols and vals go together",
            ));
        };
        let r = Ints::read(r, "rows", Some(model))?;
        let c = Ints::read(c, "cols", Some(model))?;
        let v = F64s::read_array(v, "vals")?;
        let (rs, cs, vs) = (r.slice(), c.slice(), v.slice());
        if rs.len() != cs.len() || rs.len() != vs.len() {
            return Err(crate::diag::model_err(format!(
                "set_quadratic_objective: rows, cols, vals have {}, {} and {} entries",
                rs.len(),
                cs.len(),
                vs.len()
            )));
        }
        for (k, &x) in vs.iter().enumerate() {
            put(rs.get(k), cs.get(k), x)?;
        }
    } else if let Some(q) = q {
        if let Ok(d) = q.cast::<PyDict>() {
            for (k, v) in d.iter() {
                let (i, j): (i64, i64) = k.extract().map_err(|_| {
                    PyTypeError::new_err("Q dict keys must be (i, j) pairs of variable indices")
                })?;
                put(i, j, crate::ingest::coef(&v)?)?;
            }
        } else {
            let m = read_matrix(q)?;
            let (nr, nc) = m.shape;
            if nr != nc {
                return Err(crate::diag::model_err(format!(
                    "Q must be square, got {nr} x {nc}"
                )));
            }
            if nr > n {
                return Err(crate::diag::model_err(format!(
                    "Q is {nr} x {nr} but the model has {n} variables"
                )));
            }
            let data = m.data.slice();
            match &m.kind {
                MatKind::Csr(p, i) => {
                    let (p, i) = (p.slice(), i.slice());
                    for r in 0..nr {
                        let (a, b) = (p.get(r) as usize, p.get(r + 1) as usize);
                        for (k, &x) in data.iter().enumerate().take(b).skip(a) {
                            put(r as i64, i.get(k), x)?;
                        }
                    }
                }
                MatKind::Csc(p, i) => {
                    let (p, i) = (p.slice(), i.slice());
                    for c in 0..nc {
                        let (a, b) = (p.get(c) as usize, p.get(c + 1) as usize);
                        for (k, &x) in data.iter().enumerate().take(b).skip(a) {
                            put(i.get(k), c as i64, x)?;
                        }
                    }
                }
                MatKind::Coo(r, c) => {
                    let (r, c) = (r.slice(), c.slice());
                    for (k, &x) in data.iter().enumerate() {
                        put(r.get(k), c.get(k), x)?;
                    }
                }
                MatKind::Dense => {
                    for r in 0..nr {
                        for c in 0..nc {
                            put(r as i64, c as i64, data[r * nc + c])?;
                        }
                    }
                }
            }
        }
    } else {
        return Ok(None);
    }
    // symmetrise: H[i][j] = H[j][i] = (raw[i][j] + raw[j][i]) / 2
    let mut sym: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for (&(i, j), &v) in &raw {
        *sym.entry((i, j)).or_insert(0.0) += v / 2.0;
        *sym.entry((j, i)).or_insert(0.0) += v / 2.0;
    }
    let out: Vec<(usize, usize, f64)> = sym
        .into_iter()
        .filter(|&(_, v)| v != 0.0)
        .map(|((i, j), v)| (i, j, v))
        .collect();
    Ok(if out.is_empty() { None } else { Some(out) })
}

/// The square column-wise `HighsHessian` (dimension `n`) of `entries`; dimension 0 when empty.
pub(crate) fn highs_hessian(entries: &[(usize, usize, f64)], n: usize) -> HighsHessian {
    if entries.is_empty() {
        return HighsHessian::default();
    }
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(i, j, v) in entries {
        cols[j].push((i, v));
    }
    let mut h = HighsHessian {
        dim_: n as HighsInt,
        format_: HessianFormat::kSquare,
        start_: Vec::with_capacity(n + 1),
        index_: Vec::with_capacity(entries.len()),
        value_: Vec::with_capacity(entries.len()),
    };
    h.start_.push(0);
    for c in &mut cols {
        c.sort_by_key(|&(i, _)| i);
        for &(i, v) in c.iter() {
            h.index_.push(i as HighsInt);
            h.value_.push(v);
        }
        h.start_.push(h.index_.len() as HighsInt);
    }
    h
}

type Entries = Vec<(usize, usize, f64)>;

/// What the convexity check found, kept for `Model.quad_info()`.
#[derive(Clone)]
pub(crate) struct QuadInfo {
    pub class: String,
    pub verdict: &'static str,
    pub method: &'static str,
    pub shift: f64,
    pub min_eig: Option<f64>,
}

fn col_label(lp: &HighsLp, j: usize) -> String {
    match lp.col_names_.get(j) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("c{j}"),
    }
}

/// The convexity check for a quadratic objective on the solver's current model: refuse an
/// integer variable (MIQP is not supported) every time; run the convexity analysis only when `analyse`.
/// Returns the Hessian entries to pass (shifted if the check asked for it) and the report.
pub(crate) fn assess(
    lp: &HighsLp,
    entries: &[(usize, usize, f64)],
    maximize: bool,
    analyse: bool,
) -> PyResult<Option<(Entries, QuadInfo)>> {
    let n = lp.num_col_ as usize;
    let (mut cont, mut bin, mut other) = (0, 0, 0);
    for j in 0..n {
        let int = lp
            .integrality_
            .get(j)
            .is_some_and(|t| *t != HighsVarType::kContinuous);
        if !int {
            cont += 1;
        } else if lp.col_lower_[j] >= 0.0 && lp.col_upper_[j] <= 1.0 {
            bin += 1;
        } else {
            other += 1;
        }
    }
    if bin + other > 0 {
        return Err(crate::diag::model_err(format!(
            "quadratic objective with {} integer variable(s): MIQP is not supported yet",
            bin + other
        )));
    }
    if !analyse {
        return Ok(None);
    }
    let c: Convexity = analyze(n, entries, maximize);
    let finite = lp
        .col_lower_
        .iter()
        .chain(&lp.col_upper_)
        .any(|v| v.is_finite());
    let class = class_code(&c, cont, bin, other, lp.num_row_ as usize, finite);
    if c.verdict == Verdict::Nonconvex {
        let vars: Vec<String> = c.witness.iter().map(|&(j, _)| col_label(lp, j)).collect();
        let lam = c.min_eig_estimate.unwrap_or(f64::NAN);
        let (kind, shown, need) = if maximize {
            ("largest", -lam, "negative")
        } else {
            ("smallest", lam, "positive")
        };
        return Err(crate::diag::model_err(format!(
            "nonconvex quadratic objective: Q has {kind} eigenvalue {shown:.6e} on variables {}; \
             this {} needs Q {need} semidefinite (class {class})",
            vars.join(", "),
            if maximize {
                "maximisation"
            } else {
                "minimisation"
            },
        )));
    }
    let mut out = entries.to_vec();
    if c.shift != 0.0 {
        let sh = if maximize { -c.shift } else { c.shift };
        let mut has = vec![false; n];
        for e in out.iter_mut().filter(|e| e.0 == e.1) {
            e.2 += sh;
            has[e.0] = true;
        }
        for (j, &h) in has.iter().enumerate() {
            if !h {
                out.push((j, j, sh));
            }
        }
    }
    let info = QuadInfo {
        class,
        verdict: match c.verdict {
            Verdict::Convex => "convex",
            Verdict::ConvexWithShift => "convex_after_shift",
            Verdict::Unverified => "unverified",
            Verdict::Nonconvex => "nonconvex",
        },
        method: match c.method {
            Method::Empty => "empty",
            Method::Diagonal => "diagonal",
            Method::Gershgorin => "gershgorin",
            Method::Ldlt => "ldlt",
            Method::Lanczos => "lanczos",
        },
        shift: c.shift,
        min_eig: c.min_eig_estimate,
    };
    Ok(Some((out, info)))
}
