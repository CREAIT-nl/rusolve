//! `(HighsLp, Option<CpModel>)` to `CpModelProto`.
//!
//! Variables: column `j` of the model is variable `j` of the proto, so a CP-SAT solution is read
//! back by position. Variables appended after the columns (see the no-overlap rule) are
//! auxiliary and counted in [`Translation::aux_variables`].
//!
//! Rows: `lo <= a x <= hi` becomes one linear constraint with the domain `[lo, hi]` (an infinite
//! side is `i64::MIN` / `i64::MAX`); a free row is dropped. Coefficients, row bounds and (for the
//! integer objective) costs must be integers of magnitude at most 2^53, else a [`Refusal`] names
//! the first one. A fractional bound of an integer column is rounded inward (`ceil` of the lower,
//! `floor` of the upper), which is exact for an integer column.
//!
//! Objective: integral costs go in as `CpObjectiveProto` (`scaling_factor` -1 and negated costs
//! for a maximization, so the response reports the objective in the model's sense with its
//! offset); otherwise as `FloatObjectiveProto`, and CP-SAT scales it itself.
//!
//! CP objects (`rusolve_io::cp`, semantics in its module docs):
//! * interval: an `IntervalConstraintProto` enforced by the presence column, plus the two linear
//!   constraints CP-SAT leaves to the client (`present => start + size = end`,
//!   `present => size >= 0`);
//! * no-overlap: CP-SAT counts a zero-size interval as one that must not lie strictly inside
//!   another, ours overlaps nothing. An interval whose size can be 0 therefore enters the
//!   constraint through an auxiliary optional interval, present when the interval is present and
//!   its size is at least 1 (an auxiliary Boolean `b <=> size >= 1`, and `q <=> present and b`);
//! * cumulative, all-different, element, table: the matching CP-SAT constraint; CP-SAT ignores
//!   zero-size intervals in a cumulative, as `cp.rs` defines.

use crate::pb::{self, constraint_proto::Constraint};
use rusolve_ekk::prelude::{HighsLp, HighsVarType, MatrixFormat, ObjSense, kHighsInf};
use rusolve_io::cp::{CpModel, CpTerm};
use std::fmt;

/// 2^53: the largest magnitude below which every integer is an exact `f64`.
const EXACT: f64 = 9_007_199_254_740_992.0;

/// How the objective was passed to CP-SAT.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObjectiveKind {
    /// All costs integral: `CpObjectiveProto`, exact.
    Integer,
    /// Some cost is fractional: `FloatObjectiveProto` (CP-SAT scales it).
    FloatingPoint,
}

impl ObjectiveKind {
    pub fn name(self) -> &'static str {
        match self {
            ObjectiveKind::Integer => "integer",
            ObjectiveKind::FloatingPoint => "floating_point",
        }
    }
}

/// The first thing CP-SAT cannot take exactly.
#[derive(Clone, Debug, PartialEq)]
pub enum Refusal {
    ContinuousColumn {
        col: usize,
        name: String,
    },
    InfiniteColumnBound {
        col: usize,
        name: String,
        side: &'static str,
    },
    EmptyDomain {
        col: usize,
        name: String,
        lower: f64,
        upper: f64,
    },
    NonIntegralCoefficient {
        row: usize,
        name: String,
        col_name: String,
        value: f64,
    },
    NonIntegralRowBound {
        row: usize,
        name: String,
        side: &'static str,
        value: f64,
    },
    OutOfRange {
        what: String,
        value: f64,
    },
    /// A CP object that cannot be translated as given (bad index, size, shape).
    Cp {
        object: String,
        detail: String,
    },
}

impl fmt::Display for Refusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Refusal::ContinuousColumn { col, name } => write!(
                f,
                "column {name} (index {col}) is continuous: CP-SAT takes integer columns only"
            ),
            Refusal::InfiniteColumnBound { col, name, side } => write!(
                f,
                "integer column {name} (index {col}) has no finite {side} bound: CP-SAT needs finite domains"
            ),
            Refusal::EmptyDomain {
                col,
                name,
                lower,
                upper,
            } => write!(
                f,
                "integer column {name} (index {col}) has no integer in [{lower}, {upper}]"
            ),
            Refusal::NonIntegralCoefficient {
                row,
                name,
                col_name,
                value,
            } => write!(
                f,
                "row {name} (index {row}) has the non-integer coefficient {value} on column {col_name}: CP-SAT takes integer coefficients"
            ),
            Refusal::NonIntegralRowBound {
                row,
                name,
                side,
                value,
            } => write!(
                f,
                "row {name} (index {row}) has the non-integer {side} bound {value}: CP-SAT takes integer right-hand sides"
            ),
            Refusal::OutOfRange { what, value } => write!(
                f,
                "{what} = {value} is beyond 2^53, the largest magnitude CP-SAT receives exactly"
            ),
            Refusal::Cp { object, detail } => write!(f, "CP object {object}: {detail}"),
        }
    }
}

/// A translated model.
#[derive(Clone, Debug)]
pub struct Translation {
    pub model: pb::CpModelProto,
    pub objective_kind: ObjectiveKind,
    /// Variables beyond the model's columns.
    pub aux_variables: usize,
}

fn col_name(lp: &HighsLp, j: usize) -> String {
    match lp.col_names_.get(j) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("c{j}"),
    }
}

fn row_name(lp: &HighsLp, i: usize) -> String {
    match lp.row_names_.get(i) {
        Some(n) if !n.is_empty() => n.clone(),
        _ => format!("r{i}"),
    }
}

fn is_integral(v: f64) -> bool {
    v == v.round()
}

fn expr_col(j: usize) -> pb::LinearExpressionProto {
    pb::LinearExpressionProto {
        vars: vec![j as i32],
        coeffs: vec![1],
        offset: 0,
    }
}

fn expr_const(c: i64) -> pb::LinearExpressionProto {
    pb::LinearExpressionProto {
        vars: vec![],
        coeffs: vec![],
        offset: c,
    }
}

fn linear(vars: Vec<i32>, coeffs: Vec<i64>, lo: i64, hi: i64) -> pb::LinearConstraintProto {
    pb::LinearConstraintProto {
        vars,
        coeffs,
        domain: vec![lo, hi],
    }
}

fn constraint(enforcement: Vec<i32>, c: Constraint) -> pb::ConstraintProto {
    pb::ConstraintProto {
        name: String::new(),
        enforcement_literal: enforcement,
        constraint: Some(c),
    }
}

/// The negation of literal `b` in the proto's encoding.
fn not(b: i32) -> i32 {
    -b - 1
}

struct Ctx<'a> {
    lp: &'a HighsLp,
    n: usize,
    /// `[lower, upper]` of every column after rounding.
    dom: Vec<(i64, i64)>,
    vars: Vec<pb::IntegerVariableProto>,
    cons: Vec<pb::ConstraintProto>,
    /// Constraint index of each interval of the CP model.
    interval_ix: Vec<usize>,
    /// The auxiliary interval standing for an interval in no-overlap sets, made once.
    nonzero_ix: Vec<Option<usize>>,
}

impl Ctx<'_> {
    fn cp_err<T>(&self, object: &str, detail: String) -> Result<T, Refusal> {
        Err(Refusal::Cp {
            object: object.to_string(),
            detail,
        })
    }

    fn col(&self, object: &str, j: usize) -> Result<usize, Refusal> {
        if j < self.n {
            Ok(j)
        } else {
            self.cp_err(
                object,
                format!("column index {j} is outside the model's {} columns", self.n),
            )
        }
    }

    fn term(&self, object: &str, t: CpTerm) -> Result<pb::LinearExpressionProto, Refusal> {
        match t {
            CpTerm::Const(c) => {
                if (c as f64).abs() > EXACT {
                    return Err(Refusal::OutOfRange {
                        what: format!("a constant of {object}"),
                        value: c as f64,
                    });
                }
                Ok(expr_const(c))
            }
            CpTerm::Col(j) => Ok(expr_col(self.col(object, j)?)),
        }
    }

    /// Lower bound of a term over the column domains.
    fn term_lower(&self, t: CpTerm) -> i64 {
        match t {
            CpTerm::Const(c) => c,
            CpTerm::Col(j) => self.dom.get(j).map_or(i64::MIN, |d| d.0),
        }
    }

    fn new_bool(&mut self, name: &str) -> i32 {
        self.vars.push(pb::IntegerVariableProto {
            name: name.to_string(),
            domain: vec![0, 1],
        });
        (self.vars.len() - 1) as i32
    }

    fn presence(&self, object: &str, p: Option<usize>) -> Result<Vec<i32>, Refusal> {
        match p {
            None => Ok(vec![]),
            Some(j) => {
                let j = self.col(object, j)?;
                let (lo, hi) = self.dom[j];
                if lo < 0 || hi > 1 {
                    return self.cp_err(
                        object,
                        format!(
                            "presence column {} has domain [{lo}, {hi}], not a subset of [0, 1]",
                            col_name(self.lp, j)
                        ),
                    );
                }
                Ok(vec![j as i32])
            }
        }
    }
}

/// Translates the model; see the module docs for what is exact and what is refused.
pub fn translate(lp: &HighsLp, cp: Option<&CpModel>) -> Result<Translation, Refusal> {
    let n = lp.num_col_ as usize;
    let m = lp.num_row_ as usize;
    // columns
    let mut dom = Vec::with_capacity(n);
    let mut vars = Vec::with_capacity(n);
    for j in 0..n {
        let name = col_name(lp, j);
        if lp.integrality_.get(j) != Some(&HighsVarType::kInteger) {
            return Err(Refusal::ContinuousColumn { col: j, name });
        }
        let (lo, hi) = (lp.col_lower_[j], lp.col_upper_[j]);
        if !(lo.is_finite() && lo > -kHighsInf) {
            return Err(Refusal::InfiniteColumnBound {
                col: j,
                name,
                side: "lower",
            });
        }
        if !(hi.is_finite() && hi < kHighsInf) {
            return Err(Refusal::InfiniteColumnBound {
                col: j,
                name,
                side: "upper",
            });
        }
        for (what, v) in [("lower", lo), ("upper", hi)] {
            if v.abs() > EXACT {
                return Err(Refusal::OutOfRange {
                    what: format!("the {what} bound of column {name}"),
                    value: v,
                });
            }
        }
        let (l, u) = (lo.ceil(), hi.floor());
        if l > u {
            return Err(Refusal::EmptyDomain {
                col: j,
                name,
                lower: lo,
                upper: hi,
            });
        }
        dom.push((l as i64, u as i64));
        vars.push(pb::IntegerVariableProto {
            name,
            domain: vec![l as i64, u as i64],
        });
    }

    // rows, by row
    let mut rows: Vec<Vec<(i32, i64)>> = vec![Vec::new(); m];
    let a = &lp.a_matrix_;
    let outer = if a.format_ == MatrixFormat::kRowwise {
        m
    } else {
        n
    };
    for k in 0..outer {
        let (s, e) = (a.start_[k] as usize, a.start_[k + 1] as usize);
        for p in s..e {
            let (i, j) = if a.format_ == MatrixFormat::kRowwise {
                (k, a.index_[p] as usize)
            } else {
                (a.index_[p] as usize, k)
            };
            let v = a.value_[p];
            if v == 0.0 {
                continue;
            }
            if !is_integral(v) {
                return Err(Refusal::NonIntegralCoefficient {
                    row: i,
                    name: row_name(lp, i),
                    col_name: col_name(lp, j),
                    value: v,
                });
            }
            if v.abs() > EXACT {
                return Err(Refusal::OutOfRange {
                    what: format!(
                        "the coefficient of column {} in row {}",
                        col_name(lp, j),
                        row_name(lp, i)
                    ),
                    value: v,
                });
            }
            rows[i].push((j as i32, v as i64));
        }
    }
    let mut cons = Vec::with_capacity(m);
    for (i, entries) in rows.into_iter().enumerate() {
        let (lo, hi) = (lp.row_lower_[i], lp.row_upper_[i]);
        let side = |v: f64, name: &'static str, inf: i64| -> Result<i64, Refusal> {
            if v.abs() >= kHighsInf {
                return Ok(inf);
            }
            if !is_integral(v) {
                return Err(Refusal::NonIntegralRowBound {
                    row: i,
                    name: row_name(lp, i),
                    side: name,
                    value: v,
                });
            }
            if v.abs() > EXACT {
                return Err(Refusal::OutOfRange {
                    what: format!("the {name} bound of row {}", row_name(lp, i)),
                    value: v,
                });
            }
            Ok(v as i64)
        };
        let (l, u) = (side(lo, "lower", i64::MIN)?, side(hi, "upper", i64::MAX)?);
        if l == i64::MIN && u == i64::MAX {
            continue;
        }
        if entries.is_empty() && l <= 0 && 0 <= u {
            continue;
        }
        let (vs, cs): (Vec<i32>, Vec<i64>) = entries.into_iter().unzip();
        let mut c = constraint(vec![], Constraint::Linear(linear(vs, cs, l, u)));
        c.name = row_name(lp, i);
        cons.push(c);
    }

    // objective
    let maximize = lp.sense_ == ObjSense::kMaximize;
    let sign = if maximize { -1.0 } else { 1.0 };
    let mut obj_vars = Vec::new();
    let mut costs = Vec::new();
    for (j, &c) in lp.col_cost_.iter().enumerate() {
        if c != 0.0 {
            obj_vars.push(j as i32);
            costs.push(c);
        }
    }
    let integral_costs = costs.iter().all(|&c| is_integral(c) && c.abs() <= EXACT);
    let (objective, floating, kind) = if integral_costs {
        (
            Some(pb::CpObjectiveProto {
                vars: obj_vars,
                coeffs: costs.iter().map(|&c| (sign * c) as i64).collect(),
                // user objective = scaling_factor * (sum + offset)
                offset: sign * lp.offset_,
                scaling_factor: sign,
                ..Default::default()
            }),
            None,
            ObjectiveKind::Integer,
        )
    } else {
        (
            None,
            Some(pb::FloatObjectiveProto {
                vars: obj_vars,
                coeffs: costs,
                offset: lp.offset_,
                maximize,
            }),
            ObjectiveKind::FloatingPoint,
        )
    };

    let mut cx = Ctx {
        lp,
        n,
        dom,
        vars,
        cons,
        interval_ix: Vec::new(),
        nonzero_ix: Vec::new(),
    };
    if let Some(cp) = cp {
        translate_cp(&mut cx, cp)?;
    }
    let aux_variables = cx.vars.len() - n;
    Ok(Translation {
        model: pb::CpModelProto {
            name: lp.model_name_.clone(),
            variables: cx.vars,
            constraints: cx.cons,
            objective,
            floating_point_objective: floating,
            ..Default::default()
        },
        objective_kind: kind,
        aux_variables,
    })
}

fn translate_cp(cx: &mut Ctx, cp: &CpModel) -> Result<(), Refusal> {
    // intervals first: every later object refers to them by constraint index
    for iv in &cp.intervals {
        let obj = format!("interval {}", iv.name);
        let start = cx.col(&obj, iv.start)?;
        let end = cx.col(&obj, iv.end)?;
        let size = cx.term(&obj, iv.size)?;
        let enforcement = cx.presence(&obj, iv.presence)?;
        if let CpTerm::Const(c) = iv.size
            && c < 0
        {
            return cx.cp_err(&obj, format!("the constant size {c} is negative"));
        }
        cx.interval_ix.push(cx.cons.len());
        cx.nonzero_ix.push(None);
        cx.cons.push(constraint(
            enforcement.clone(),
            Constraint::Interval(pb::IntervalConstraintProto {
                start: Some(expr_col(start)),
                end: Some(expr_col(end)),
                size: Some(size),
            }),
        ));
        // present => start + size = end, and size >= 0
        match iv.size {
            CpTerm::Const(c) => cx.cons.push(constraint(
                enforcement.clone(),
                Constraint::Linear(linear(vec![start as i32, end as i32], vec![1, -1], -c, -c)),
            )),
            CpTerm::Col(s) => {
                cx.cons.push(constraint(
                    enforcement.clone(),
                    Constraint::Linear(linear(
                        vec![start as i32, s as i32, end as i32],
                        vec![1, 1, -1],
                        0,
                        0,
                    )),
                ));
                if cx.dom[s].0 < 0 {
                    cx.cons.push(constraint(
                        enforcement.clone(),
                        Constraint::Linear(linear(vec![s as i32], vec![1], 0, i64::MAX)),
                    ));
                }
            }
        }
    }
    let ix = |cx: &Ctx, object: &str, i: usize| -> Result<usize, Refusal> {
        cx.interval_ix.get(i).copied().map_or_else(
            || {
                cx.cp_err(
                    object,
                    format!(
                        "interval index {i} is outside the model's {} intervals",
                        cx.interval_ix.len()
                    ),
                )
            },
            Ok,
        )
    };
    for no in &cp.no_overlap {
        let obj = format!("no-overlap {}", no.name);
        let mut set = Vec::with_capacity(no.intervals.len());
        for &i in &no.intervals {
            ix(cx, &obj, i)?;
            if let Some(k) = no_overlap_member(cx, &cp.intervals[i], i)? {
                set.push(k as i32);
            }
        }
        cx.cons.push(constraint(
            vec![],
            Constraint::NoOverlap(pb::NoOverlapConstraintProto { intervals: set }),
        ));
    }
    for cu in &cp.cumulative {
        let obj = format!("cumulative {}", cu.name);
        if cu.demands.len() != cu.intervals.len() {
            return cx.cp_err(
                &obj,
                format!(
                    "{} demands for {} intervals",
                    cu.demands.len(),
                    cu.intervals.len()
                ),
            );
        }
        let mut set = Vec::new();
        for &i in &cu.intervals {
            set.push(ix(cx, &obj, i)? as i32);
        }
        let mut demands = Vec::new();
        for &d in &cu.demands {
            if cx.term_lower(d) < 0 {
                return cx.cp_err(
                    &obj,
                    "a demand can be negative (CP-SAT needs demands >= 0)".into(),
                );
            }
            demands.push(cx.term(&obj, d)?);
        }
        let capacity = cx.term(&obj, cu.capacity)?;
        cx.cons.push(constraint(
            vec![],
            Constraint::Cumulative(pb::CumulativeConstraintProto {
                capacity: Some(capacity),
                intervals: set,
                demands,
            }),
        ));
    }
    for ad in &cp.all_different {
        let obj = format!("all-different {}", ad.name);
        let mut exprs = Vec::new();
        for &j in &ad.cols {
            exprs.push(expr_col(cx.col(&obj, j)?));
        }
        cx.cons.push(constraint(
            vec![],
            Constraint::AllDiff(pb::AllDifferentConstraintProto { exprs }),
        ));
    }
    for el in &cp.element {
        let obj = format!("element {}", el.name);
        let index = expr_col(cx.col(&obj, el.index)?);
        let target = expr_col(cx.col(&obj, el.target)?);
        let mut exprs = Vec::new();
        for &v in &el.values {
            exprs.push(cx.term(&obj, v)?);
        }
        if exprs.is_empty() {
            return cx.cp_err(&obj, "no values".into());
        }
        cx.cons.push(constraint(
            vec![],
            Constraint::Element(pb::ElementConstraintProto {
                linear_index: Some(index),
                linear_target: Some(target),
                exprs,
                ..Default::default()
            }),
        ));
    }
    for ta in &cp.table {
        let obj = format!("table {}", ta.name);
        let mut exprs = Vec::new();
        for &j in &ta.cols {
            exprs.push(expr_col(cx.col(&obj, j)?));
        }
        let mut values = Vec::with_capacity(ta.tuples.len() * ta.cols.len());
        for (r, t) in ta.tuples.iter().enumerate() {
            if t.len() != ta.cols.len() {
                return cx.cp_err(
                    &obj,
                    format!(
                        "tuple {r} has {} values for {} columns",
                        t.len(),
                        ta.cols.len()
                    ),
                );
            }
            values.extend_from_slice(t);
        }
        cx.cons.push(constraint(
            vec![],
            Constraint::Table(pb::TableConstraintProto {
                values,
                exprs,
                negated: ta.negated,
                ..Default::default()
            }),
        ));
    }
    Ok(())
}

/// The interval constraint a no-overlap set holds for interval `i`: the interval itself when its
/// size is at least 1 whenever it is present, none for a constant size 0, else an auxiliary interval present only when the
/// interval is present and its size is at least 1.
fn no_overlap_member(
    cx: &mut Ctx,
    iv: &rusolve_io::cp::Interval,
    i: usize,
) -> Result<Option<usize>, Refusal> {
    if cx.term_lower(iv.size) >= 1 {
        return Ok(Some(cx.interval_ix[i]));
    }
    if iv.size == CpTerm::Const(0) {
        // a zero-size interval overlaps nothing: it is not in the set
        return Ok(None);
    }
    if let Some(k) = cx.nonzero_ix[i] {
        return Ok(Some(k));
    }
    let obj = format!("interval {}", iv.name);
    let enforcement = cx.presence(&obj, iv.presence)?;
    let size = cx.term(&obj, iv.size)?;
    // b <=> size >= 1
    let b = cx.new_bool(&format!("{}_nonzero", iv.name));
    let (svars, scoef, soff) = (size.vars.clone(), size.coeffs.clone(), size.offset);
    cx.cons.push(constraint(
        vec![b],
        Constraint::Linear(linear(svars.clone(), scoef.clone(), 1 - soff, i64::MAX)),
    ));
    cx.cons.push(constraint(
        vec![not(b)],
        Constraint::Linear(linear(svars, scoef, i64::MIN, -soff)),
    ));
    // q <=> present and b
    let q = match enforcement.first() {
        None => b,
        Some(&p) => {
            let q = cx.new_bool(&format!("{}_counted", iv.name));
            cx.cons.push(constraint(
                vec![q],
                Constraint::BoolAnd(pb::BoolArgumentProto {
                    literals: vec![p, b],
                }),
            ));
            cx.cons.push(constraint(
                vec![],
                Constraint::BoolOr(pb::BoolArgumentProto {
                    literals: vec![not(p), not(b), q],
                }),
            ));
            q
        }
    };
    let start = expr_col(cx.col(&obj, iv.start)?);
    let end = expr_col(cx.col(&obj, iv.end)?);
    let k = cx.cons.len();
    cx.cons.push(constraint(
        vec![q],
        Constraint::Interval(pb::IntervalConstraintProto {
            start: Some(start),
            end: Some(end),
            size: Some(size),
        }),
    ));
    cx.nonzero_ix[i] = Some(k);
    Ok(Some(k))
}
