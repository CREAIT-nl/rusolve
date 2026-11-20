//! Shared by the format tests: a model as a problem, independent of how a
//! writer spelled it.

use rusolve_io::{Model, ObjSense, RowSense};

const INF: f64 = f64::INFINITY;

/// Single-sided rows `(terms by name, sense, rhs)` with ranged and
/// two-sided rows split in two, columns named `ONE_VAR_CONSTANT` folded
/// into the objective constant, and the bounds/costs/integrality by name.
#[derive(Debug, PartialEq)]
pub struct Canon {
    sense: ObjSense,
    offset: f64,
    rows: Vec<CanonRow>,
    cols: Vec<CanonCol>,
}

/// (terms by name, 'E' / 'G' / 'L', rhs)
pub type CanonRow = (Vec<(String, f64)>, char, f64);
/// (name, cost, lower, upper, integer)
pub type CanonCol = (String, f64, f64, f64, bool);

pub fn canon(m: &Model) -> Canon {
    const ONE: &str = "ONE_VAR_CONSTANT";
    let mut terms: Vec<Vec<(String, f64)>> = vec![Vec::new(); m.num_constraints];
    let (ptr, idx, val) = (m.matrix.col_ptr(), m.matrix.row_idx(), m.matrix.values());
    for j in 0..m.num_variables {
        for k in ptr[j]..ptr[j + 1] {
            terms[idx[k] as usize].push((m.col_names[j].clone(), val[k]));
        }
    }
    let mut rows = Vec::new();
    for (i, mut t) in terms.into_iter().enumerate() {
        if t.iter().all(|(n, _)| n == ONE) {
            continue; // Pyomo MPS: `ONE_VAR_CONSTANT = 1`
        }
        t.sort_by(|a, b| a.0.cmp(&b.0));
        let (b, r) = (m.rhs[i], m.ranges[i]);
        let (lo, hi) = match m.senses[i] {
            RowSense::LessEqual if r != 0.0 => (b - r.abs(), b),
            RowSense::LessEqual => (-INF, b),
            RowSense::GreaterEqual if r != 0.0 => (b, b + r.abs()),
            RowSense::GreaterEqual => (b, INF),
            RowSense::Equal if r > 0.0 => (b, b + r),
            RowSense::Equal if r < 0.0 => (b + r, b),
            RowSense::Equal => (b, b),
        };
        if lo == hi {
            rows.push((t, 'E', lo));
            continue;
        }
        if lo > -INF {
            rows.push((t.clone(), 'G', lo));
        }
        if hi < INF {
            rows.push((t, 'L', hi));
        }
    }
    rows.sort_by(|a, b| format!("{a:?}").cmp(&format!("{b:?}")));
    let mut offset = m.obj_offset;
    let mut cols = Vec::new();
    for j in 0..m.num_variables {
        let b = m.bounds[j];
        if m.col_names[j] == ONE {
            offset += m.objective[j];
            continue;
        }
        cols.push((
            m.col_names[j].clone(),
            m.objective[j],
            b.lower,
            b.upper,
            b.integer,
        ));
    }
    cols.sort_by(|a, b| a.0.cmp(&b.0));
    Canon {
        sense: m.obj_sense,
        offset,
        rows,
        cols,
    }
}
