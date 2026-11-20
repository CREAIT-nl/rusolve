//! `Var`, `Vars`, `LinExpr`, `Constraint` and `quicksum` (docs/design/python-api.md,
//! "Expression layer"; it avoids a Python object per term).
//!
//! - A `Var` is a column index plus the id of its model; it allocates nothing else.
//! - A `Vars` is a strided range of columns; no per-element object exists until it is indexed.
//! - A `LinExpr` is a term list plus a constant. Terms live in a shared, append-only buffer
//!   ([`Terms`]): an expression is a prefix of a buffer, so `e + x` appends to `e`'s buffer when
//!   nothing was appended after `e` yet, and copies only otherwise. Chains `a + b + c ...` and
//!   the builtin `sum()` (which never calls `+=`) are therefore linear, not quadratic, and no
//!   earlier expression ever sees a later append (it reads only its own prefix).
//! - A `Constraint` holds the body (`lhs - rhs`, one or two inline terms or a shared prefix)
//!   and the bounds `lo <= body <= hi` with the constants moved into the bounds.
//!
//! Duplicate terms are kept here and merged by `rusolve_model::ModelBuilder` when a row is
//! built.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};

use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyTypeError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyFloat, PyList, PyRange, PySlice, PyTuple};

use crate::ingest;

/// A term: (column index, coefficient).
pub(crate) type Term = (i32, f64);

static NEXT_MODEL_ID: AtomicU64 = AtomicU64::new(1);

/// A fresh, never reused model id (0 means "no variable yet" in an expression).
pub(crate) fn new_model_id() -> u64 {
    NEXT_MODEL_ID.fetch_add(1, Ordering::Relaxed)
}

pub(crate) const MIXED_MODELS: &str = "cannot combine variables of two different Models";

/// The model id of a combination of `a` and `b` (0 = no variable).
#[inline]
pub(crate) fn join(a: u64, b: u64) -> PyResult<u64> {
    if a == 0 || a == b {
        Ok(b)
    } else if b == 0 {
        Ok(a)
    } else {
        Err(crate::diag::model_err(MIXED_MODELS))
    }
}

// ---------------------------------------------------------------- Terms

/// A prefix (`len` terms) of a shared append-only term buffer. See the module docs.
#[derive(Clone)]
pub(crate) struct Terms {
    buf: Arc<Mutex<Vec<Term>>>,
    len: usize,
}

impl Terms {
    pub(crate) fn from_vec(v: Vec<Term>) -> Self {
        let len = v.len();
        Terms {
            buf: Arc::new(Mutex::new(v)),
            len,
        }
    }

    fn lock(&self) -> MutexGuard<'_, Vec<Term>> {
        // A poisoned lock only means a panic elsewhere; the terms are still valid.
        self.buf.lock().unwrap_or_else(|e| e.into_inner())
    }

    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// Run `f` on this expression's terms.
    pub(crate) fn with<R>(&self, f: impl FnOnce(&[Term]) -> R) -> R {
        let g = self.lock();
        f(&g[..self.len])
    }

    /// These terms followed by those `f` pushes. In place when this prefix is the whole
    /// buffer (amortised O(pushed)), otherwise the prefix is copied first.
    pub(crate) fn extended(&self, extra: usize, f: impl FnOnce(&mut Vec<Term>)) -> Terms {
        let mut g = self.lock();
        if g.len() == self.len {
            f(&mut g);
            let len = g.len();
            drop(g);
            Terms {
                buf: self.buf.clone(),
                len,
            }
        } else {
            let mut v = Vec::with_capacity(self.len + extra);
            v.extend_from_slice(&g[..self.len]);
            drop(g);
            f(&mut v);
            Terms::from_vec(v)
        }
    }

    /// `self + scale * other`.
    pub(crate) fn plus(&self, other: &Terms, scale: f64) -> Terms {
        if Arc::ptr_eq(&self.buf, &other.buf) {
            // same buffer: copy `other` out first (the lock is not reentrant)
            let tmp = other.with(|t| t.to_vec());
            return self.extended(tmp.len(), |v| push_scaled(v, &tmp, scale));
        }
        let og = other.lock();
        let src = &og[..other.len];
        self.extended(src.len(), |v| push_scaled(v, src, scale))
    }

    /// `scale * self`, always a fresh buffer.
    pub(crate) fn scaled(&self, scale: f64) -> Terms {
        self.with(|t| {
            let mut v = Vec::with_capacity(t.len());
            push_scaled(&mut v, t, scale);
            Terms::from_vec(v)
        })
    }
}

#[inline]
fn push_scaled(v: &mut Vec<Term>, src: &[Term], scale: f64) {
    if scale == 1.0 {
        v.extend_from_slice(src);
    } else {
        v.extend(src.iter().map(|&(j, c)| (j, c * scale)));
    }
}

// ---------------------------------------------------------------- operands

/// A right-hand operand of an arithmetic or comparison operator.
pub(crate) enum Operand {
    Num(f64),
    Var(i32, u64),
    Expr(Terms, f64, u64),
}

/// `None` when `o` is not a number, `Var` or `LinExpr` (the operator then returns
/// NotImplemented). Var is checked before numbers because it has `__index__`.
pub(crate) fn operand(o: &Bound<'_, PyAny>) -> PyResult<Option<Operand>> {
    if let Ok(f) = o.cast_exact::<PyFloat>() {
        return Ok(Some(Operand::Num(f.value())));
    }
    if let Ok(v) = o.cast_exact::<Var>() {
        let v = v.get();
        return Ok(Some(Operand::Var(v.index, v.model)));
    }
    if let Ok(e) = o.cast_exact::<LinExpr>() {
        let e = e.try_borrow()?;
        return Ok(Some(Operand::Expr(e.terms.clone(), e.constant, e.model)));
    }
    Ok(ingest::number(o)?.map(Operand::Num))
}

fn quadratic() -> PyErr {
    PyTypeError::new_err(
        "rusolve expressions are linear: a variable or expression can only be multiplied or \
         divided by a number",
    )
}

fn strict() -> PyErr {
    PyTypeError::new_err("strict inequalities are not supported: use <=, >= or ==")
}

fn not_equal() -> PyErr {
    PyTypeError::new_err("!= is not a linear constraint")
}

fn divisor(c: f64) -> PyResult<f64> {
    if c == 0.0 {
        Err(PyZeroDivisionError::new_err(
            "division of an expression by zero",
        ))
    } else {
        Ok(1.0 / c)
    }
}

fn expr_obj(py: Python<'_>, terms: Terms, constant: f64, model: u64) -> PyResult<Py<PyAny>> {
    Ok(Py::new(
        py,
        LinExpr {
            terms,
            constant,
            model,
        },
    )?
    .into_any())
}

fn not_implemented(py: Python<'_>) -> Py<PyAny> {
    py.NotImplemented()
}

// ---------------------------------------------------------------- Var

/// A model variable: a column index and the id of its model.
#[pyclass(frozen, module = "rusolve")]
pub struct Var {
    pub(crate) index: i32,
    pub(crate) model: u64,
}

impl Var {
    pub(crate) fn new(index: i32, model: u64) -> Self {
        Var { index, model }
    }

    fn one(&self, coef: f64, constant: f64) -> LinExpr {
        LinExpr {
            terms: Terms::from_vec(vec![(self.index, coef)]),
            constant,
            model: self.model,
        }
    }

    /// `self - other` as a constraint body with bounds `[lo, hi]` for `op`.
    fn compare(&self, other: Operand, op: CompareOp) -> PyResult<Constraint> {
        let me = (self.index, 1.0);
        let (body, shift, flip, model, truth) = match other {
            // x op c
            Operand::Num(c) => (Body::One(me), -c, false, self.model, None),
            // x - y op 0
            Operand::Var(j, m) => (
                Body::Two(me, (j, -1.0)),
                0.0,
                false,
                join(self.model, m)?,
                matches!(op, CompareOp::Eq).then_some(j == self.index),
            ),
            // x op e  <=>  e - x (flipped op) -e.c: append to e's buffer instead of negating it
            Operand::Expr(t, c, m) => (
                Body::Shared(t.extended(1, |v| v.push((self.index, -1.0)))),
                c,
                true,
                join(self.model, m)?,
                None,
            ),
        };
        Constraint::from_cmp(body, shift, op, flip, model, truth)
    }
}

#[pymethods]
impl Var {
    /// The column index.
    #[getter]
    fn index(&self) -> i32 {
        self.index
    }
    fn __index__(&self) -> i64 {
        self.index as i64
    }
    fn __int__(&self) -> i64 {
        self.index as i64
    }
    fn __hash__(&self) -> isize {
        self.index as isize
    }
    fn __repr__(&self) -> String {
        format!("Var({})", self.index)
    }

    fn __add__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        Ok(match operand(other)? {
            None => not_implemented(py),
            Some(Operand::Num(c)) => Py::new(py, self.one(1.0, c))?.into_any(),
            Some(Operand::Var(j, m)) => expr_obj(
                py,
                Terms::from_vec(vec![(self.index, 1.0), (j, 1.0)]),
                0.0,
                join(self.model, m)?,
            )?,
            Some(Operand::Expr(t, c, m)) => expr_obj(
                py,
                t.extended(1, |v| v.push((self.index, 1.0))),
                c,
                join(self.model, m)?,
            )?,
        })
    }
    fn __radd__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.__add__(py, other)
    }
    fn __sub__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        Ok(match operand(other)? {
            None => not_implemented(py),
            Some(Operand::Num(c)) => Py::new(py, self.one(1.0, -c))?.into_any(),
            Some(Operand::Var(j, m)) => expr_obj(
                py,
                Terms::from_vec(vec![(self.index, 1.0), (j, -1.0)]),
                0.0,
                join(self.model, m)?,
            )?,
            Some(Operand::Expr(t, c, m)) => {
                let mut v = Vec::with_capacity(t.len() + 1);
                v.push((self.index, 1.0));
                t.with(|s| push_scaled(&mut v, s, -1.0));
                expr_obj(py, Terms::from_vec(v), -c, join(self.model, m)?)?
            }
        })
    }
    fn __rsub__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        Ok(match operand(other)? {
            None => not_implemented(py),
            Some(Operand::Num(c)) => Py::new(py, self.one(-1.0, c))?.into_any(),
            Some(Operand::Var(j, m)) => expr_obj(
                py,
                Terms::from_vec(vec![(j, 1.0), (self.index, -1.0)]),
                0.0,
                join(self.model, m)?,
            )?,
            Some(Operand::Expr(t, c, m)) => expr_obj(
                py,
                t.extended(1, |v| v.push((self.index, -1.0))),
                c,
                join(self.model, m)?,
            )?,
        })
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(Operand::Num(c)) => Ok(Py::new(py, self.one(c, 0.0))?.into_any()),
            Some(_) => Err(quadratic()),
        }
    }
    fn __rmul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.__mul__(py, other)
    }
    fn __truediv__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(Operand::Num(c)) => Ok(Py::new(py, self.one(divisor(c)?, 0.0))?.into_any()),
            Some(_) => Err(quadratic()),
        }
    }
    fn __neg__(&self) -> LinExpr {
        self.one(-1.0, 0.0)
    }
    fn __pos__(&self) -> LinExpr {
        self.one(1.0, 0.0)
    }

    fn __richcmp__(
        &self,
        py: Python<'_>,
        other: &Bound<'_, PyAny>,
        op: CompareOp,
    ) -> PyResult<Py<PyAny>> {
        let Some(o) = operand(other)? else {
            return Ok(not_implemented(py));
        };
        if matches!(op, CompareOp::Ne) {
            // `x != y` keeps plain identity semantics so `in`/`remove` work on lists of Vars
            if let Operand::Var(j, m) = o {
                return Ok(
                    pyo3::types::PyBool::new(py, j != self.index || m != self.model)
                        .to_owned()
                        .into_any()
                        .unbind(),
                );
            }
            return Err(not_equal());
        }
        if let (CompareOp::Eq, Operand::Var(_, m)) = (op, &o)
            && *m != self.model
        {
            // Vars of two Models are never equal (and make no constraint): plain `False`
            // keeps dict lookups that mix models (equal hashes) working.
            return Ok(pyo3::types::PyBool::new(py, false)
                .to_owned()
                .into_any()
                .unbind());
        }
        Ok(Py::new(py, self.compare(o, op)?)?.into_any())
    }

    /// `lo <= self <= hi` as a Constraint.
    fn between(&self, lo: f64, hi: f64) -> PyResult<Constraint> {
        Constraint::ranged(Body::One((self.index, 1.0)), lo, hi, 0.0, self.model)
    }
}

// ---------------------------------------------------------------- Vars

/// A block of variables: columns `start, start + step, ...` (`len` of them).
#[pyclass(frozen, module = "rusolve")]
pub struct Vars {
    pub(crate) start: i64,
    pub(crate) step: i64,
    pub(crate) len: usize,
    pub(crate) model: u64,
}

impl Vars {
    #[inline]
    pub(crate) fn col(&self, k: usize) -> i32 {
        (self.start + self.step * k as i64) as i32
    }
    fn end(&self) -> i64 {
        self.start + self.step * self.len as i64
    }
    /// Contiguous ascending block `start..start+len`?
    pub(crate) fn as_range(&self) -> Option<std::ops::Range<usize>> {
        (self.step == 1 || self.len <= 1)
            .then(|| self.start as usize..self.start as usize + self.len)
    }
}

#[pymethods]
impl Vars {
    fn __len__(&self) -> usize {
        self.len
    }

    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(s) = key.cast::<PySlice>() {
            let ind = s.indices(self.len as isize)?;
            return Ok(Py::new(
                py,
                Vars {
                    start: self.start + self.step * ind.start as i64,
                    step: self.step * ind.step as i64,
                    len: ind.slicelength,
                    model: self.model,
                },
            )?
            .into_any());
        }
        if key.is_instance_of::<PyTuple>() {
            return Err(PyTypeError::new_err(
                "Vars is one-dimensional: index it with an int or a slice",
            ));
        }
        let k: isize = key.extract().map_err(|_| {
            PyTypeError::new_err(format!(
                "Vars indices must be int or slice, not {}",
                ingest::type_name(key)
            ))
        })?;
        let n = self.len as isize;
        let k = if k < 0 { k + n } else { k };
        if k < 0 || k >= n {
            return Err(crate::diag::index_err(format!(
                "Vars index out of range (len {})",
                self.len
            )));
        }
        Ok(Py::new(py, Var::new(self.col(k as usize), self.model))?.into_any())
    }

    fn __iter__(&self) -> VarsIter {
        VarsIter {
            next: self.start,
            step: self.step,
            left: self.len,
            model: self.model,
        }
    }

    fn __repr__(&self) -> String {
        if self.step == 1 {
            format!("Vars({}..{})", self.start, self.end())
        } else {
            format!("Vars({}..{} step {})", self.start, self.end(), self.step)
        }
    }

    /// The column indices, as a `range`.
    #[getter]
    fn indices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyRange>> {
        PyRange::new_with_step(
            py,
            self.start as isize,
            self.end() as isize,
            self.step as isize,
        )
    }
    #[getter]
    fn start(&self) -> i64 {
        self.start
    }
    #[getter]
    fn stop(&self) -> i64 {
        self.end()
    }

    /// The sum of the block.
    fn sum(&self) -> LinExpr {
        LinExpr {
            terms: Terms::from_vec((0..self.len).map(|k| (self.col(k), 1.0)).collect()),
            constant: 0.0,
            model: self.model,
        }
    }

    /// `sum(coefs[k] * self[k])`; `coefs` a sequence or buffer of length `len(self)`.
    fn dot(&self, coefs: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        let c = ingest::F64s::read(coefs, "coefs")?;
        let c = c.per_item(self.len, "coefs")?;
        Ok(LinExpr {
            terms: Terms::from_vec((0..self.len).map(|k| (self.col(k), c.get(k))).collect()),
            constant: 0.0,
            model: self.model,
        })
    }
}

#[pyclass(module = "rusolve")]
pub struct VarsIter {
    next: i64,
    step: i64,
    left: usize,
    model: u64,
}

#[pymethods]
impl VarsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self) -> Option<Var> {
        if self.left == 0 {
            return None;
        }
        let v = Var::new(self.next as i32, self.model);
        self.next += self.step;
        self.left -= 1;
        Some(v)
    }
}

// ---------------------------------------------------------------- LinExpr

/// A linear expression `sum(coef * var) + constant`.
#[pyclass(module = "rusolve")]
pub struct LinExpr {
    pub(crate) terms: Terms,
    pub(crate) constant: f64,
    pub(crate) model: u64,
}

impl LinExpr {
    /// `self + scale * other` (fresh object).
    fn combine(&self, other: Operand, scale: f64) -> PyResult<LinExpr> {
        Ok(match other {
            Operand::Num(c) => LinExpr {
                terms: self.terms.clone(),
                constant: self.constant + scale * c,
                model: self.model,
            },
            Operand::Var(j, m) => LinExpr {
                terms: self.terms.extended(1, |v| v.push((j, scale))),
                constant: self.constant,
                model: join(self.model, m)?,
            },
            Operand::Expr(t, c, m) => LinExpr {
                terms: self.terms.plus(&t, scale),
                constant: self.constant + scale * c,
                model: join(self.model, m)?,
            },
        })
    }

    fn scaled(&self, s: f64) -> LinExpr {
        LinExpr {
            terms: if s == 1.0 {
                self.terms.clone()
            } else {
                self.terms.scaled(s)
            },
            constant: self.constant * s,
            model: self.model,
        }
    }

    /// The same object's operand, read before `self` is borrowed mutably (`e += e`).
    fn self_operand(&self) -> Operand {
        Operand::Expr(self.terms.clone(), self.constant, self.model)
    }

    fn in_place(
        slf: &Bound<'_, Self>,
        other: &Bound<'_, PyAny>,
        scale: f64,
        opname: &str,
    ) -> PyResult<()> {
        let o = if other.is(slf) {
            slf.borrow().self_operand()
        } else {
            match operand(other)? {
                Some(o) => o,
                None => {
                    return Err(PyTypeError::new_err(format!(
                        "unsupported operand type for {opname}: 'LinExpr' and '{}'",
                        ingest::type_name(other)
                    )));
                }
            }
        };
        let mut me = slf.borrow_mut();
        let new = me.combine(o, scale)?;
        *me = new;
        Ok(())
    }

    fn compare(&self, other: Operand, op: CompareOp) -> PyResult<Constraint> {
        let (terms, shift, model) = match other {
            Operand::Num(c) => (self.terms.clone(), self.constant - c, self.model),
            Operand::Var(j, m) => (
                self.terms.extended(1, |v| v.push((j, -1.0))),
                self.constant,
                join(self.model, m)?,
            ),
            Operand::Expr(t, c, m) => (
                self.terms.plus(&t, -1.0),
                self.constant - c,
                join(self.model, m)?,
            ),
        };
        Constraint::from_cmp(Body::Shared(terms), shift, op, false, model, None)
    }
}

#[pymethods]
impl LinExpr {
    /// `LinExpr(terms=None, constant=0.0)`; `terms` as `quicksum` accepts, or a dict
    /// `{var: coef}`.
    #[new]
    #[pyo3(signature = (terms=None, constant=0.0))]
    fn new(terms: Option<&Bound<'_, PyAny>>, constant: f64) -> PyResult<Self> {
        let mut acc = Acc::default();
        if let Some(t) = terms {
            if let Ok(d) = t.cast::<PyDict>() {
                for (k, v) in d.iter() {
                    let (j, m) = ingest::var_key(&k)?;
                    acc.model = join(acc.model, m)?;
                    acc.terms.push((j, ingest::coef(&v)?));
                }
            } else {
                acc.add_iterable(t)?;
            }
        }
        Ok(acc.into_expr(constant))
    }

    #[getter]
    fn constant(&self) -> f64 {
        self.constant
    }
    #[setter]
    fn set_constant(&mut self, c: f64) {
        self.constant = c;
    }

    /// Number of terms (duplicates counted).
    fn __len__(&self) -> usize {
        self.terms.len()
    }

    /// The terms as a list of `(index, coef)` (duplicates not merged).
    fn terms(&self) -> Vec<(i32, f64)> {
        self.terms.with(|t| t.to_vec())
    }

    fn __repr__(&self) -> String {
        let mut s = String::from("LinExpr(");
        self.terms.with(|t| {
            for (k, &(j, c)) in t.iter().take(8).enumerate() {
                if k > 0 {
                    s += if c < 0.0 { " - " } else { " + " };
                    s += &format!("{} x{j}", c.abs());
                } else {
                    s += &format!("{c} x{j}");
                }
            }
            if t.len() > 8 {
                s += &format!(" + ... ({} terms)", t.len());
            }
            if t.is_empty() || self.constant != 0.0 {
                if !t.is_empty() {
                    s += if self.constant < 0.0 { " - " } else { " + " };
                    s += &format!("{}", self.constant.abs());
                } else {
                    s += &format!("{}", self.constant);
                }
            }
        });
        s + ")"
    }

    fn __add__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(o) => Ok(Py::new(py, self.combine(o, 1.0)?)?.into_any()),
        }
    }
    fn __radd__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.__add__(py, other)
    }
    fn __sub__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(o) => Ok(Py::new(py, self.combine(o, -1.0)?)?.into_any()),
        }
    }
    fn __rsub__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(o) => Ok(Py::new(py, self.scaled(-1.0).combine(o, 1.0)?)?.into_any()),
        }
    }
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(Operand::Num(c)) => Ok(Py::new(py, self.scaled(c))?.into_any()),
            Some(_) => Err(quadratic()),
        }
    }
    fn __rmul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.__mul__(py, other)
    }
    fn __truediv__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            None => Ok(not_implemented(py)),
            Some(Operand::Num(c)) => Ok(Py::new(py, self.scaled(divisor(c)?))?.into_any()),
            Some(_) => Err(quadratic()),
        }
    }
    fn __neg__(&self) -> LinExpr {
        self.scaled(-1.0)
    }
    fn __pos__(&self) -> LinExpr {
        self.scaled(1.0)
    }

    fn __iadd__(slf: &Bound<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        Self::in_place(slf, other, 1.0, "+=")
    }
    fn __isub__(slf: &Bound<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        Self::in_place(slf, other, -1.0, "-=")
    }
    fn __imul__(&mut self, other: &Bound<'_, PyAny>) -> PyResult<()> {
        match operand(other)? {
            Some(Operand::Num(c)) => {
                *self = self.scaled(c);
                Ok(())
            }
            Some(_) => Err(quadratic()),
            None => Err(PyTypeError::new_err(format!(
                "unsupported operand type for *=: 'LinExpr' and '{}'",
                ingest::type_name(other)
            ))),
        }
    }
    fn __itruediv__(&mut self, other: &Bound<'_, PyAny>) -> PyResult<()> {
        match operand(other)? {
            Some(Operand::Num(c)) => {
                *self = self.scaled(divisor(c)?);
                Ok(())
            }
            Some(_) => Err(quadratic()),
            None => Err(PyTypeError::new_err(format!(
                "unsupported operand type for /=: 'LinExpr' and '{}'",
                ingest::type_name(other)
            ))),
        }
    }

    fn __richcmp__(
        &self,
        py: Python<'_>,
        other: &Bound<'_, PyAny>,
        op: CompareOp,
    ) -> PyResult<Py<PyAny>> {
        let Some(o) = operand(other)? else {
            return Ok(not_implemented(py));
        };
        if matches!(op, CompareOp::Ne) {
            return Err(not_equal());
        }
        Ok(Py::new(py, self.compare(o, op)?)?.into_any())
    }

    /// `lo <= self <= hi` as a Constraint.
    fn between(&self, lo: f64, hi: f64) -> PyResult<Constraint> {
        Constraint::ranged(
            Body::Shared(self.terms.clone()),
            lo,
            hi,
            self.constant,
            self.model,
        )
    }
}

// ---------------------------------------------------------------- Constraint

/// The body of a constraint: up to two inline terms (no allocation), or a shared prefix.
pub(crate) enum Body {
    One(Term),
    Two(Term, Term),
    Shared(Terms),
}

impl Body {
    pub(crate) fn with<R>(&self, f: impl FnOnce(&[Term]) -> R) -> R {
        match self {
            Body::One(t) => f(std::slice::from_ref(t)),
            Body::Two(a, b) => f(&[*a, *b]),
            Body::Shared(t) => t.with(f),
        }
    }
}

/// `lo <= body <= hi`. Made by comparing expressions or by `between`.
#[pyclass(frozen, module = "rusolve")]
pub struct Constraint {
    pub(crate) body: Body,
    pub(crate) lo: f64,
    pub(crate) hi: f64,
    pub(crate) model: u64,
    /// `Var == Var`: the truth value of the comparison (index equality), so Vars work as
    /// dict keys and in `in` tests. `None` for every other constraint.
    truth: Option<bool>,
}

impl Constraint {
    /// `body + shift (op) 0`, or with `flip` the operator reversed (`body` is then `rhs - lhs`).
    fn from_cmp(
        body: Body,
        shift: f64,
        op: CompareOp,
        flip: bool,
        model: u64,
        truth: Option<bool>,
    ) -> PyResult<Self> {
        let inf = f64::INFINITY;
        let op = match (op, flip) {
            (CompareOp::Le, true) => CompareOp::Ge,
            (CompareOp::Ge, true) => CompareOp::Le,
            (o, _) => o,
        };
        let (lo, hi) = match op {
            CompareOp::Le => (-inf, -shift),
            CompareOp::Ge => (-shift, inf),
            CompareOp::Eq => (-shift, -shift),
            CompareOp::Lt | CompareOp::Gt => return Err(strict()),
            CompareOp::Ne => return Err(not_equal()),
        };
        Ok(Constraint {
            body,
            lo: lo + 0.0,
            hi: hi + 0.0,
            model,
            truth,
        })
    }

    fn ranged(body: Body, lo: f64, hi: f64, constant: f64, model: u64) -> PyResult<Self> {
        if lo.is_nan() || hi.is_nan() {
            return Err(crate::diag::model_err("between: bounds must not be NaN"));
        }
        Ok(Constraint {
            body,
            lo: lo - constant,
            hi: hi - constant,
            model,
            truth: None,
        })
    }
}

#[pymethods]
impl Constraint {
    #[getter]
    fn lo(&self) -> f64 {
        self.lo
    }
    #[getter]
    fn hi(&self) -> f64 {
        self.hi
    }
    /// The body as a LinExpr (constant 0; the constants are in `lo`/`hi`).
    #[getter]
    fn expr(&self) -> LinExpr {
        LinExpr {
            terms: self.body.with(|t| Terms::from_vec(t.to_vec())),
            constant: 0.0,
            model: self.model,
        }
    }

    fn __bool__(&self) -> PyResult<bool> {
        self.truth.ok_or_else(|| {
            PyTypeError::new_err(
                "a Constraint has no truth value. For a ranged row write \
                 expr.between(lo, hi) instead of lo <= expr <= hi; to add it, pass it to \
                 Model.add_constr",
            )
        })
    }

    fn __repr__(&self) -> String {
        let e = self.expr().__repr__();
        if self.lo == self.hi {
            format!("Constraint({e} == {})", self.hi)
        } else if self.lo == f64::NEG_INFINITY {
            format!("Constraint({e} <= {})", self.hi)
        } else if self.hi == f64::INFINITY {
            format!("Constraint({e} >= {})", self.lo)
        } else {
            format!("Constraint({} <= {e} <= {})", self.lo, self.hi)
        }
    }
}

// ---------------------------------------------------------------- quicksum

/// An accumulating sum (quicksum, LinExpr(...)).
#[derive(Default)]
pub(crate) struct Acc {
    pub(crate) terms: Vec<Term>,
    pub(crate) constant: f64,
    pub(crate) model: u64,
}

impl Acc {
    fn into_expr(self, extra_constant: f64) -> LinExpr {
        LinExpr {
            terms: Terms::from_vec(self.terms),
            constant: self.constant + extra_constant,
            model: self.model,
        }
    }

    fn add_scaled(&mut self, o: Operand, s: f64) -> PyResult<()> {
        match o {
            Operand::Num(c) => self.constant += s * c,
            Operand::Var(j, m) => {
                self.model = join(self.model, m)?;
                self.terms.push((j, s));
            }
            Operand::Expr(t, c, m) => {
                self.model = join(self.model, m)?;
                t.with(|src| push_scaled(&mut self.terms, src, s));
                self.constant += s * c;
            }
        }
        Ok(())
    }

    fn add_item(&mut self, item: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Ok(v) = item.cast_exact::<Var>() {
            let v = v.get();
            self.model = join(self.model, v.model)?;
            self.terms.push((v.index, 1.0));
            return Ok(());
        }
        if let Some(o) = operand(item)? {
            return self.add_scaled(o, 1.0);
        }
        if let Ok(t) = item.cast_exact::<PyTuple>()
            && t.len() == 2
        {
            let (a, b) = (t.get_item(0)?, t.get_item(1)?);
            // (coef, var) or (var, coef)
            let (c, e) = match (ingest::number(&a)?, ingest::number(&b)?) {
                (Some(c), _) if !a.is_instance_of::<Var>() => (c, b),
                (_, Some(c)) if !b.is_instance_of::<Var>() => (c, a),
                _ => return Err(bad_item(item)),
            };
            return match operand(&e)? {
                Some(o @ (Operand::Var(..) | Operand::Expr(..))) => self.add_scaled(o, c),
                _ => Err(bad_item(item)),
            };
        }
        if let Ok(vs) = item.cast_exact::<Vars>() {
            let vs = vs.get();
            self.model = join(self.model, vs.model)?;
            self.terms.extend((0..vs.len).map(|k| (vs.col(k), 1.0)));
            return Ok(());
        }
        Err(bad_item(item))
    }

    pub(crate) fn add_iterable(&mut self, it: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Ok(l) = it.cast_exact::<PyList>() {
            self.terms.reserve(l.len());
            for k in 0..l.len() {
                self.add_item(&l.get_item(k)?)?;
            }
        } else if let Ok(t) = it.cast_exact::<PyTuple>() {
            for item in t.iter() {
                self.add_item(&item)?;
            }
        } else if let Ok(vs) = it.cast_exact::<Vars>() {
            self.add_item(vs.as_any())?;
        } else {
            for item in it.try_iter()? {
                self.add_item(&item?)?;
            }
        }
        Ok(())
    }
}

fn bad_item(item: &Bound<'_, PyAny>) -> PyErr {
    PyTypeError::new_err(format!(
        "quicksum: items must be Var, LinExpr, numbers or (coef, var) pairs, got {}",
        ingest::type_name(item)
    ))
}

/// Sum Vars, LinExprs, numbers and `(coef, var)` pairs into one LinExpr, without
/// intermediate objects.
#[pyfunction]
fn quicksum(iterable: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
    let mut acc = Acc::default();
    acc.add_iterable(iterable)?;
    Ok(acc.into_expr(0.0))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Var>()?;
    m.add_class::<Vars>()?;
    m.add_class::<LinExpr>()?;
    m.add_class::<Constraint>()?;
    m.add_function(wrap_pyfunction!(quicksum, m)?)?;
    Ok(())
}
