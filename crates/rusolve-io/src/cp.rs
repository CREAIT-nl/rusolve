//! Scheduling and combinatorial constraints (the CP layer) of a model.
//!
//! Like [`crate::gc`], this stores the *user's* objects; nothing here is a
//! reformulation. A constraint-programming engine reads them as they are, a MIP
//! fallback linearizes them, and the checker reports violations on these
//! objects, not on derived rows.
//!
//! Every column referenced here must be integer with finite bounds: the objects
//! are defined over integers. Indices are column indices of the model; interval
//! references are indices into [`CpModel::intervals`].
//!
//! Object semantics (the standard constraint-programming definitions, as in the
//! Global Constraint Catalogue, Beldiceanu, Carlsson and Rampon, SICS T2005:08,
//! and the public OR-Tools CP-SAT documentation):
//! - interval `i`: when present, `start + size = end`; when absent (its presence
//!   literal is 0) it takes part in no constraint and its columns are free
//!   within their bounds;
//! - no-overlap: no two present intervals of the set overlap in time (an
//!   interval of size 0 overlaps nothing);
//! - cumulative: at every time `t`, the demands of the present intervals with
//!   `start <= t < end` sum to at most the capacity;
//! - all-different: the columns take pairwise different values;
//! - element: `target = values[index]`, with `index` 0-based;
//! - table: the tuple of the columns' values is one of the rows of `tuples`
//!   (or none of them when `negated`).

use crate::model::VariableBound;

/// A constant or a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CpTerm {
    Const(i64),
    Col(usize),
}

/// An interval: `start + size = end` when present.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub name: String,
    pub start: usize,
    pub size: CpTerm,
    pub end: usize,
    /// A binary column; `None` means always present.
    pub presence: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoOverlap {
    pub name: String,
    pub intervals: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cumulative {
    pub name: String,
    pub intervals: Vec<usize>,
    /// One demand per interval, same order as `intervals`.
    pub demands: Vec<CpTerm>,
    pub capacity: CpTerm,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllDifferent {
    pub name: String,
    pub cols: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Element {
    pub name: String,
    pub index: usize,
    pub values: Vec<CpTerm>,
    pub target: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    pub name: String,
    pub cols: Vec<usize>,
    /// Each row has one value per column of `cols`.
    pub tuples: Vec<Vec<i64>>,
    pub negated: bool,
}

/// The CP objects of a model. Empty for a plain LP/MIP.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CpModel {
    pub intervals: Vec<Interval>,
    pub no_overlap: Vec<NoOverlap>,
    pub cumulative: Vec<Cumulative>,
    pub all_different: Vec<AllDifferent>,
    pub element: Vec<Element>,
    pub table: Vec<Table>,
}

/// A typed refusal naming the object: an invalid object (`Invalid`, a model error) or a valid one
/// the MIP linearization cannot take (`Unsupported`, the CP-SAT engine can).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CpError {
    pub kind: CpErrorKind,
    /// `interval 'j1'`, `no_overlap 'machine1'`, ...
    pub object: String,
    pub message: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CpErrorKind {
    Invalid,
    Unsupported,
}

impl CpError {
    pub fn invalid(object: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            kind: CpErrorKind::Invalid,
            object: object.into(),
            message: message.into(),
        }
    }

    pub fn unsupported(object: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            kind: CpErrorKind::Unsupported,
            object: object.into(),
            message: message.into(),
        }
    }
}

impl std::fmt::Display for CpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.object, self.message)
    }
}

impl std::error::Error for CpError {}

/// One violated object found by [`CpModel::check`].
#[derive(Debug, Clone, PartialEq)]
pub struct CpViolation {
    /// `interval`, `no_overlap`, `cumulative`, `all_different`, `element` or `table`.
    pub kind: &'static str,
    /// The object's name.
    pub name: String,
    /// How far from satisfying it the point is (overlap length, excess demand, difference, 1 for
    /// a combinatorial violation).
    pub amount: f64,
    pub detail: String,
}

/// Tolerances of [`CpModel::check`].
#[derive(Debug, Clone, Copy)]
pub struct CpTol {
    /// A column of an object is rounded to the nearest integer when within this; farther is a
    /// violation of the object (the objects are defined over integers).
    pub integrality: f64,
}

impl Default for CpTol {
    fn default() -> Self {
        Self { integrality: 1e-6 }
    }
}

/// A cumulative's capacity and its present intervals as (start, end, demand).
type CumulativeLoad = (i64, Vec<(i64, i64, i64)>);

impl CpModel {
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
            && self.no_overlap.is_empty()
            && self.cumulative.is_empty()
            && self.all_different.is_empty()
            && self.element.is_empty()
            && self.table.is_empty()
    }

    /// Per-kind object counts, for messages and the envelope.
    pub fn counts(&self) -> [(&'static str, usize); 6] {
        [
            ("interval", self.intervals.len()),
            ("no_overlap", self.no_overlap.len()),
            ("cumulative", self.cumulative.len()),
            ("all_different", self.all_different.len()),
            ("element", self.element.len()),
            ("table", self.table.len()),
        ]
    }

    /// The first object of the model, named, for typed "not supported here" errors.
    pub fn first_object(&self) -> Option<String> {
        if let Some(o) = self.no_overlap.first() {
            return Some(format!("no_overlap '{}'", o.name));
        }
        if let Some(o) = self.cumulative.first() {
            return Some(format!("cumulative '{}'", o.name));
        }
        if let Some(o) = self.all_different.first() {
            return Some(format!("all_different '{}'", o.name));
        }
        if let Some(o) = self.element.first() {
            return Some(format!("element '{}'", o.name));
        }
        if let Some(o) = self.table.first() {
            return Some(format!("table '{}'", o.name));
        }
        self.intervals
            .first()
            .map(|o| format!("interval '{}'", o.name))
    }

    /// Typed refusal for commands that cannot honour CP objects (they would drop them and report
    /// a wrong answer). `who` is the command or call.
    pub fn reject(&self, who: &str) -> Result<(), String> {
        match self.first_object() {
            Some(what) => Err(format!(
                "{who}: the model has {what}; scheduling and combinatorial (CP) objects are only supported by solve and check"
            )),
            None => Ok(()),
        }
    }

    /// Structural checks on the model as it is now: every referenced column exists, is integer
    /// with finite bounds, presence columns are binary, lengths match, sizes and demands are
    /// nonnegative, interval references are in range. The error names the object.
    pub fn validate(&self, bounds: &[VariableBound], col_names: &[String]) -> Result<(), CpError> {
        let n = bounds.len();
        let nm = |j: usize| {
            col_names
                .get(j)
                .filter(|s| !s.is_empty())
                .cloned()
                .unwrap_or_else(|| format!("c{j}"))
        };
        // an integer column with finite bounds
        let dom = |obj: &str, what: &str, j: usize| -> Result<(), CpError> {
            if j >= n {
                return Err(CpError::invalid(
                    obj,
                    format!("{what}: column index {j} out of range ({n} columns)"),
                ));
            }
            let b = &bounds[j];
            if !b.integer {
                return Err(CpError::invalid(
                    obj,
                    format!(
                        "{what}: variable '{}' is not integer (CP objects are defined over integer variables)",
                        nm(j)
                    ),
                ));
            }
            if !(b.lower.is_finite() && b.upper.is_finite()) {
                return Err(CpError::invalid(
                    obj,
                    format!(
                        "{what}: variable '{}' needs finite bounds (it has [{}, {}])",
                        nm(j),
                        b.lower,
                        b.upper
                    ),
                ));
            }
            Ok(())
        };
        let term = |obj: &str, what: &str, t: &CpTerm, nonneg: bool| -> Result<(), CpError> {
            match *t {
                CpTerm::Const(c) => {
                    if nonneg && c < 0 {
                        return Err(CpError::invalid(
                            obj,
                            format!("{what} is {c}; it must be nonnegative"),
                        ));
                    }
                    Ok(())
                }
                CpTerm::Col(j) => {
                    dom(obj, what, j)?;
                    if nonneg && bounds[j].lower < 0.0 {
                        return Err(CpError::invalid(
                            obj,
                            format!(
                                "{what}: variable '{}' can be negative (lower bound {}); it must be nonnegative",
                                nm(j),
                                bounds[j].lower
                            ),
                        ));
                    }
                    Ok(())
                }
            }
        };
        for it in &self.intervals {
            let obj = format!("interval '{}'", it.name);
            dom(&obj, "start", it.start)?;
            dom(&obj, "end", it.end)?;
            term(&obj, "size", &it.size, true)?;
            if let Some(p) = it.presence {
                dom(&obj, "presence", p)?;
                let b = &bounds[p];
                if b.lower < 0.0 || b.upper > 1.0 {
                    return Err(CpError::invalid(
                        &obj,
                        format!(
                            "presence variable '{}' is not binary (needs bounds within [0, 1])",
                            nm(p)
                        ),
                    ));
                }
            }
        }
        let k = self.intervals.len();
        let ivs = |obj: &str, list: &[usize]| -> Result<(), CpError> {
            for &i in list {
                if i >= k {
                    return Err(CpError::invalid(
                        obj,
                        format!("interval index {i} out of range ({k} intervals)"),
                    ));
                }
            }
            let mut seen = list.to_vec();
            seen.sort_unstable();
            if let Some(w) = seen.windows(2).find(|w| w[0] == w[1]) {
                return Err(CpError::invalid(
                    obj,
                    format!("interval '{}' is listed twice", self.intervals[w[0]].name),
                ));
            }
            Ok(())
        };
        for o in &self.no_overlap {
            ivs(&format!("no_overlap '{}'", o.name), &o.intervals)?;
        }
        for o in &self.cumulative {
            let obj = format!("cumulative '{}'", o.name);
            ivs(&obj, &o.intervals)?;
            if o.demands.len() != o.intervals.len() {
                return Err(CpError::invalid(
                    &obj,
                    format!(
                        "{} demands for {} intervals (one demand per interval)",
                        o.demands.len(),
                        o.intervals.len()
                    ),
                ));
            }
            for d in &o.demands {
                term(&obj, "demand", d, true)?;
            }
            term(&obj, "capacity", &o.capacity, true)?;
        }
        for o in &self.all_different {
            let obj = format!("all_different '{}'", o.name);
            for &j in &o.cols {
                dom(&obj, "variable", j)?;
            }
        }
        for o in &self.element {
            let obj = format!("element '{}'", o.name);
            if o.values.is_empty() {
                return Err(CpError::invalid(&obj, "needs at least one value"));
            }
            dom(&obj, "index", o.index)?;
            dom(&obj, "target", o.target)?;
            for v in &o.values {
                term(&obj, "value", v, false)?;
            }
        }
        for o in &self.table {
            let obj = format!("table '{}'", o.name);
            if o.cols.is_empty() {
                return Err(CpError::invalid(&obj, "needs at least one column"));
            }
            for &j in &o.cols {
                dom(&obj, "variable", j)?;
            }
            for (r, t) in o.tuples.iter().enumerate() {
                if t.len() != o.cols.len() {
                    return Err(CpError::invalid(
                        &obj,
                        format!(
                            "tuple {r} has {} values for {} columns",
                            t.len(),
                            o.cols.len()
                        ),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Evaluates every object on the point `x` (all model columns), on the user's objects.
    ///
    /// Every column an object reads is rounded to the nearest integer when within
    /// `tol.integrality`; a column farther from an integer is itself a violation of that object.
    /// An absent interval (presence 0) takes part in nothing. No-overlap: two present intervals
    /// of positive size overlap when `max(start) < min(end)`. Cumulative: the demand in force at
    /// every interval start (the load only rises there) is compared with the capacity.
    pub fn check(&self, x: &[f64], col_names: &[String], tol: CpTol) -> Vec<CpViolation> {
        let name = |j: usize| -> String {
            match col_names.get(j) {
                Some(n) if !n.is_empty() => n.clone(),
                _ => format!("c{j}"),
            }
        };
        let mut out = Vec::new();
        // The rounded value of a column, or the violation text when it is not integral.
        let val = |j: usize| -> Result<i64, String> {
            let v = x.get(j).copied().unwrap_or(f64::NAN);
            let r = v.round();
            if (v - r).abs() <= tol.integrality && r.is_finite() {
                Ok(r as i64)
            } else {
                Err(format!("'{}' = {v} is not an integer", name(j)))
            }
        };
        let tval = |t: &CpTerm| -> Result<i64, String> {
            match *t {
                CpTerm::Const(c) => Ok(c),
                CpTerm::Col(j) => val(j),
            }
        };
        // Evaluated intervals: None = the interval's own columns are not integral (reported on it).
        struct Iv {
            present: bool,
            s: i64,
            e: i64,
            size: i64,
        }
        let mut ev: Vec<Option<Iv>> = Vec::new();
        for it in &self.intervals {
            let r = (|| -> Result<Iv, String> {
                let present = match it.presence {
                    None => true,
                    Some(p) => match val(p)? {
                        0 => false,
                        1 => true,
                        v => return Err(format!("presence '{}' = {v} is not 0 or 1", name(p))),
                    },
                };
                Ok(Iv {
                    present,
                    s: val(it.start)?,
                    e: val(it.end)?,
                    size: tval(&it.size)?,
                })
            })();
            match r {
                Ok(iv) => {
                    if iv.present && iv.s + iv.size != iv.e {
                        out.push(CpViolation {
                            kind: "interval",
                            name: it.name.clone(),
                            amount: (iv.s + iv.size - iv.e).abs() as f64,
                            detail: format!(
                                "start {} + size {} = {} but end is {}",
                                iv.s,
                                iv.size,
                                iv.s + iv.size,
                                iv.e
                            ),
                        });
                    }
                    ev.push(Some(iv));
                }
                Err(d) => {
                    out.push(CpViolation {
                        kind: "interval",
                        name: it.name.clone(),
                        amount: 1.0,
                        detail: d,
                    });
                    ev.push(None);
                }
            }
        }
        let iname = |i: usize| self.intervals[i].name.as_str();
        for o in &self.no_overlap {
            let mut bad = None;
            let mut live: Vec<usize> = Vec::new();
            for &i in &o.intervals {
                match &ev[i] {
                    None => {
                        bad = Some(format!("interval '{}' has a non-integral column", iname(i)));
                    }
                    Some(iv) if iv.present && iv.size > 0 => live.push(i),
                    _ => {}
                }
            }
            if let Some(d) = bad {
                out.push(CpViolation {
                    kind: "no_overlap",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: d,
                });
                continue;
            }
            for (a, &i) in live.iter().enumerate() {
                for &j in &live[a + 1..] {
                    let (p, q) = (ev[i].as_ref().unwrap(), ev[j].as_ref().unwrap());
                    let overlap = p.e.min(q.e) - p.s.max(q.s);
                    if overlap > 0 {
                        out.push(CpViolation {
                            kind: "no_overlap",
                            name: o.name.clone(),
                            amount: overlap as f64,
                            detail: format!(
                                "intervals '{}' [{}, {}) and '{}' [{}, {}) overlap by {}",
                                iname(i),
                                p.s,
                                p.e,
                                iname(j),
                                q.s,
                                q.e,
                                overlap
                            ),
                        });
                    }
                }
            }
        }
        for o in &self.cumulative {
            let r = (|| -> Result<CumulativeLoad, String> {
                let cap = tval(&o.capacity)?;
                let mut live = Vec::new();
                for (k, &i) in o.intervals.iter().enumerate() {
                    let Some(iv) = &ev[i] else {
                        return Err(format!("interval '{}' has a non-integral column", iname(i)));
                    };
                    let d = tval(&o.demands[k])?;
                    if iv.present && iv.size > 0 && d > 0 {
                        live.push((iv.s, iv.e, d));
                    }
                }
                Ok((cap, live))
            })();
            match r {
                Err(d) => out.push(CpViolation {
                    kind: "cumulative",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: d,
                }),
                Ok((cap, live)) => {
                    let mut worst: Option<(i64, i64)> = None; // (excess, time)
                    for &(t, _, _) in &live {
                        let load: i64 = live
                            .iter()
                            .filter(|&&(s, e, _)| s <= t && t < e)
                            .map(|&(_, _, d)| d)
                            .sum();
                        if load > cap && worst.is_none_or(|(w, _)| load - cap > w) {
                            worst = Some((load - cap, t));
                        }
                    }
                    if let Some((excess, t)) = worst {
                        out.push(CpViolation {
                            kind: "cumulative",
                            name: o.name.clone(),
                            amount: excess as f64,
                            detail: format!(
                                "demand exceeds the capacity {cap} by {excess} at time {t}"
                            ),
                        });
                    }
                }
            }
        }
        for o in &self.all_different {
            let mut vals = Vec::new();
            let mut bad = None;
            for &j in &o.cols {
                match val(j) {
                    Ok(v) => vals.push((j, v)),
                    Err(d) => {
                        bad = Some(d);
                        break;
                    }
                }
            }
            if let Some(d) = bad {
                out.push(CpViolation {
                    kind: "all_different",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: d,
                });
                continue;
            }
            for (a, &(i, vi)) in vals.iter().enumerate() {
                for &(j, vj) in &vals[a + 1..] {
                    if vi == vj {
                        out.push(CpViolation {
                            kind: "all_different",
                            name: o.name.clone(),
                            amount: 1.0,
                            detail: format!("'{}' and '{}' are both {vi}", name(i), name(j)),
                        });
                    }
                }
            }
        }
        for o in &self.element {
            let r = (|| -> Result<(i64, i64, Option<i64>), String> {
                let idx = val(o.index)?;
                let t = val(o.target)?;
                let v = if idx >= 0 && (idx as usize) < o.values.len() {
                    Some(tval(&o.values[idx as usize])?)
                } else {
                    None
                };
                Ok((idx, t, v))
            })();
            match r {
                Err(d) => out.push(CpViolation {
                    kind: "element",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: d,
                }),
                Ok((idx, _, None)) => out.push(CpViolation {
                    kind: "element",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: format!(
                        "index '{}' = {idx} is outside 0..{}",
                        name(o.index),
                        o.values.len() - 1
                    ),
                }),
                Ok((idx, t, Some(v))) => {
                    if t != v {
                        out.push(CpViolation {
                            kind: "element",
                            name: o.name.clone(),
                            amount: (t - v).abs() as f64,
                            detail: format!(
                                "target '{}' = {t} but values[{idx}] = {v}",
                                name(o.target)
                            ),
                        });
                    }
                }
            }
        }
        for o in &self.table {
            let mut row = Vec::new();
            let mut bad = None;
            for &j in &o.cols {
                match val(j) {
                    Ok(v) => row.push(v),
                    Err(d) => {
                        bad = Some(d);
                        break;
                    }
                }
            }
            if let Some(d) = bad {
                out.push(CpViolation {
                    kind: "table",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: d,
                });
                continue;
            }
            let hit = o.tuples.contains(&row);
            if hit == o.negated {
                out.push(CpViolation {
                    kind: "table",
                    name: o.name.clone(),
                    amount: 1.0,
                    detail: if o.negated {
                        format!("{row:?} is one of the forbidden tuples")
                    } else {
                        format!(
                            "{row:?} is not one of the {} allowed tuples",
                            o.tuples.len()
                        )
                    },
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|j| format!("x{j}")).collect()
    }

    fn iv(name: &str, s: usize, size: i64, e: usize) -> Interval {
        Interval {
            name: name.into(),
            start: s,
            size: CpTerm::Const(size),
            end: e,
            presence: None,
        }
    }

    fn int(lo: f64, hi: f64) -> VariableBound {
        VariableBound {
            lower: lo,
            upper: hi,
            integer: true,
        }
    }

    #[test]
    fn interval_and_no_overlap() {
        let cp = CpModel {
            intervals: vec![iv("a", 0, 3, 1), iv("b", 2, 2, 3)],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        let t = CpTol::default();
        // a = [0,3), b = [3,5): fine
        assert!(cp.check(&[0.0, 3.0, 3.0, 5.0], &names(4), t).is_empty());
        // b = [2,4): overlap 1
        let v = cp.check(&[0.0, 3.0, 2.0, 4.0], &names(4), t);
        assert_eq!(v.len(), 1);
        assert_eq!((v[0].kind, v[0].amount), ("no_overlap", 1.0));
        // a with a wrong end is an interval violation
        let v = cp.check(&[0.0, 4.0, 5.0, 7.0], &names(4), t);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].kind, "interval");
        // fractional column
        let v = cp.check(&[0.0, 3.0, 3.5, 5.5], &names(4), t);
        assert!(v.iter().any(|v| v.detail.contains("not an integer")));
    }

    #[test]
    fn absent_interval_takes_part_in_nothing() {
        let mut a = iv("a", 0, 3, 1);
        a.presence = Some(4);
        let cp = CpModel {
            intervals: vec![a, iv("b", 2, 2, 3)],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        let t = CpTol::default();
        // absent a: its columns are free (end 9 does not match), no overlap with b
        assert!(
            cp.check(&[0.0, 9.0, 0.0, 2.0, 0.0], &names(5), t)
                .is_empty()
        );
        // present a with the same columns breaks the interval and the overlap
        let v = cp.check(&[0.0, 9.0, 0.0, 2.0, 1.0], &names(5), t);
        assert!(v.iter().any(|v| v.kind == "interval"));
    }

    #[test]
    fn zero_size_overlaps_nothing() {
        let cp = CpModel {
            intervals: vec![iv("a", 0, 4, 1), iv("z", 2, 0, 3)],
            no_overlap: vec![NoOverlap {
                name: "m".into(),
                intervals: vec![0, 1],
            }],
            ..Default::default()
        };
        assert!(
            cp.check(&[0.0, 4.0, 2.0, 2.0], &names(4), CpTol::default())
                .is_empty()
        );
    }

    #[test]
    fn cumulative_capacity() {
        let cp = CpModel {
            intervals: vec![iv("a", 0, 2, 1), iv("b", 2, 2, 3), iv("c", 4, 1, 5)],
            cumulative: vec![Cumulative {
                name: "r".into(),
                intervals: vec![0, 1, 2],
                demands: vec![CpTerm::Const(2), CpTerm::Const(1), CpTerm::Const(1)],
                capacity: CpTerm::Const(3),
            }],
            ..Default::default()
        };
        let t = CpTol::default();
        let x = |s: [f64; 3]| vec![s[0], s[0] + 2.0, s[1], s[1] + 2.0, s[2], s[2] + 1.0];
        assert!(cp.check(&x([0.0, 2.0, 1.0]), &names(6), t).is_empty());
        let v = cp.check(&x([0.0, 1.0, 1.0]), &names(6), t);
        assert_eq!(v.len(), 1, "{v:?}");
        assert_eq!((v[0].kind, v[0].amount), ("cumulative", 1.0));
    }

    #[test]
    fn all_different_element_table() {
        let t = CpTol::default();
        let cp = CpModel {
            all_different: vec![AllDifferent {
                name: "d".into(),
                cols: vec![0, 1, 2],
            }],
            element: vec![Element {
                name: "e".into(),
                index: 0,
                values: vec![CpTerm::Const(5), CpTerm::Const(7), CpTerm::Col(1)],
                target: 3,
            }],
            table: vec![
                Table {
                    name: "t".into(),
                    cols: vec![0, 1],
                    tuples: vec![vec![1, 2], vec![0, 2]],
                    negated: false,
                },
                Table {
                    name: "n".into(),
                    cols: vec![0, 2],
                    tuples: vec![vec![1, 0]],
                    negated: true,
                },
            ],
            ..Default::default()
        };
        // x0 = 1 -> values[1] = 7 = x3; x1 = 2, x2 = 3 (the pair (1, 3) is not forbidden)
        assert!(cp.check(&[1.0, 2.0, 3.0, 7.0], &names(4), t).is_empty());
        // x0 = 2 -> values[2] = x1 = 2 = x3; (2, 2) not in table t, and not all-different
        let v = cp.check(&[2.0, 2.0, 0.0, 2.0], &names(4), t);
        let kinds: Vec<_> = v.iter().map(|v| v.kind).collect();
        assert!(
            kinds.contains(&"all_different") && kinds.contains(&"table"),
            "{v:?}"
        );
        assert!(!kinds.contains(&"element"), "{v:?}");
        // element target wrong (values[1] = 7, not 5), forbidden tuple (1, 0)
        let v = cp.check(&[1.0, 2.0, 0.0, 5.0], &names(4), t);
        let kinds: Vec<_> = v.iter().map(|v| (v.kind, v.name.as_str())).collect();
        assert_eq!(kinds, vec![("element", "e"), ("table", "n")], "{v:?}");
    }

    #[test]
    fn validation_names_the_object() {
        let bounds = vec![int(0.0, 10.0), int(0.0, 10.0), VariableBound::default()];
        let nm = names(3);
        let mut cp = CpModel {
            intervals: vec![iv("a", 0, 3, 1)],
            ..Default::default()
        };
        assert!(cp.validate(&bounds, &nm).is_ok());
        cp.intervals[0].end = 2; // continuous, infinite upper
        let e = cp.validate(&bounds, &nm).unwrap_err();
        assert!(
            e.object == "interval 'a'" && e.message.contains("not integer"),
            "{e}"
        );
        cp.intervals[0].end = 1;
        cp.intervals[0].size = CpTerm::Const(-1);
        assert!(
            cp.validate(&bounds, &nm)
                .unwrap_err()
                .message
                .contains("nonnegative")
        );
        cp.intervals[0].size = CpTerm::Const(1);
        cp.intervals[0].presence = Some(1); // [0, 10] is not binary
        assert!(
            cp.validate(&bounds, &nm)
                .unwrap_err()
                .message
                .contains("not binary")
        );
        cp.intervals[0].presence = None;
        cp.no_overlap.push(NoOverlap {
            name: "m".into(),
            intervals: vec![0, 3],
        });
        let e = cp.validate(&bounds, &nm).unwrap_err();
        assert!(
            e.object == "no_overlap 'm'" && e.message.contains("out of range"),
            "{e}"
        );
        cp.no_overlap.clear();
        cp.cumulative.push(Cumulative {
            name: "r".into(),
            intervals: vec![0],
            demands: vec![],
            capacity: CpTerm::Const(1),
        });
        assert!(
            cp.validate(&bounds, &nm)
                .unwrap_err()
                .message
                .contains("demands")
        );
        cp.cumulative.clear();
        cp.table.push(Table {
            name: "t".into(),
            cols: vec![0, 1],
            tuples: vec![vec![1]],
            negated: false,
        });
        assert!(
            cp.validate(&bounds, &nm)
                .unwrap_err()
                .message
                .contains("tuple 0")
        );
    }
}
