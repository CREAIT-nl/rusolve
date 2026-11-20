//! Reading Python values fast (docs/design/python-api.md): scalars, sequences and buffers to
//! `f64`/`i32` (numpy arrays through `PyBuffer`, without a copy), and rows (dict,
//! `(indices, values)` pair, `LinExpr`, `Var`, `Constraint`) straight into
//! `rusolve_model::ModelBuilder::push`.
//!
//! The rule (measured with `bench/compare_libs.py`): never go
//! through the generic sequence protocol when a faster form exists. Order of attempts for a
//! value list: exact `list`/`tuple` (indexed access, exact `float`/`int` fast paths), then
//! the buffer protocol (float64/int32/int64 read in place; other numeric dtypes converted
//! once), then any iterable.

use std::ops::Range;

use pyo3::buffer::{Element, ElementType, PyBuffer, PyUntypedBuffer};
use pyo3::exceptions::PyTypeError;
use pyo3::ffi;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyFloat, PyInt, PyList, PyRange, PyString, PyTuple};
use rusolve_model::{BuildError, Fill, ModelBuilder};

use crate::expr::{Constraint, LinExpr, MIXED_MODELS, Var, Vars};

pub(crate) fn type_name(o: &Bound<'_, PyAny>) -> String {
    o.get_type()
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|_| "?".into())
}

// ---------------------------------------------------------------- errors

/// A builder error as a Python exception: IndexError for an index out of range, ValueError
/// otherwise. `row_name` (if the row has one) is appended.
pub(crate) fn build_err(e: BuildError, row_name: Option<&str>) -> PyErr {
    let mut msg = e.to_string();
    if let Some(n) = row_name {
        msg += &format!(" (constraint {n:?})");
    }
    match e {
        BuildError::ColumnOutOfRange { .. } | BuildError::RowOutOfRange { .. } => {
            crate::diag::index_err(msg)
        }
        _ => crate::diag::model_err(msg),
    }
}

fn wrong_model() -> PyErr {
    crate::diag::model_err(format!("{MIXED_MODELS}: this Var belongs to another Model"))
}

// ---------------------------------------------------------------- scalars

/// `o` as a number, or `None` if it is not one. Exact float and int first; Var and the
/// other expression types are never numbers (Var has `__index__`); then anything with
/// `__float__`/`__index__` (numpy scalars, bool, Fraction, ...).
#[inline]
pub(crate) fn number(o: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    if let Ok(f) = o.cast_exact::<PyFloat>() {
        return Ok(Some(f.value()));
    }
    if o.is_exact_instance_of::<PyInt>() {
        return Ok(Some(o.extract::<f64>()?));
    }
    if o.is_instance_of::<Var>()
        || o.is_instance_of::<LinExpr>()
        || o.is_instance_of::<Constraint>()
        || o.is_instance_of::<Vars>()
        || o.is_instance_of::<PyString>()
    {
        return Ok(None);
    }
    Ok(o.extract::<f64>().ok())
}

/// A coefficient or bound: a number, else TypeError.
#[inline]
pub(crate) fn coef(o: &Bound<'_, PyAny>) -> PyResult<f64> {
    if let Ok(f) = o.cast_exact::<PyFloat>() {
        return Ok(f.value());
    }
    number(o)?
        .ok_or_else(|| PyTypeError::new_err(format!("expected a number, got {}", type_name(o))))
}

/// A dict key / index: a Var (index, model) or an int (index, 0).
#[inline]
pub(crate) fn var_key(o: &Bound<'_, PyAny>) -> PyResult<(i32, u64)> {
    if let Ok(v) = o.cast_exact::<Var>() {
        let v = v.get();
        return Ok((v.index, v.model));
    }
    if (o.is_instance_of::<PyInt>() || !o.is_instance_of::<PyFloat>())
        && let Ok(j) = o.extract::<i64>()
    {
        return Ok((clamp_col(j), 0));
    }
    Err(PyTypeError::new_err(format!(
        "a variable must be a Var or an int column index, got {}",
        type_name(o)
    )))
}

/// A column index as `i32`, keeping out-of-range values out of range (the builder then
/// reports them) instead of wrapping.
#[inline]
fn clamp_col(j: i64) -> i32 {
    j.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

/// A Var of model `model`, or an int index.
#[inline]
pub(crate) fn col_of(o: &Bound<'_, PyAny>, model: u64) -> PyResult<i64> {
    if let Ok(v) = o.cast_exact::<Var>() {
        let v = v.get();
        if v.model != model {
            return Err(wrong_model());
        }
        return Ok(v.index as i64);
    }
    if (o.is_instance_of::<PyInt>() || !o.is_instance_of::<PyFloat>())
        && let Ok(j) = o.extract::<i64>()
    {
        return Ok(j);
    }
    Err(PyTypeError::new_err(format!(
        "a variable must be a Var or an int column index, got {}",
        type_name(o)
    )))
}

// ---------------------------------------------------------------- buffers

/// A borrowed view of a typed, C-contiguous buffer.
fn buf_slice<T: Element>(b: &PyBuffer<T>) -> &[T] {
    let n = b.item_count();
    if n == 0 {
        return &[];
    }
    // SAFETY: `b` was checked (by `into_typed`) to hold items of exactly `T`'s size and format,
    // aligned for `T`, and (by the callers) to be C-contiguous, so `buf_ptr` points to `n`
    // initialised `T`s. The exporter keeps that memory alive and in place while the
    // `Py_buffer` is held, and the slice borrows `b`. The slice is only read while this thread
    // holds the GIL and runs no Python code, so no Python code can write to it meanwhile; the
    // remaining case, native code of another thread writing the same array with the GIL
    // released, is the same contract as numpy's own C API and rust-numpy's readonly arrays.
    unsafe { std::slice::from_raw_parts(b.buf_ptr() as *const T, n) }
}

/// A 1-D numeric buffer of `o`, if it exports one (not str/bytes).
fn numeric_buffer(o: &Bound<'_, PyAny>) -> Option<PyUntypedBuffer> {
    let b = PyUntypedBuffer::get(o).ok()?;
    match ElementType::from_format(b.format()) {
        ElementType::Float { .. } | ElementType::SignedInteger { .. } => Some(b),
        ElementType::UnsignedInteger { bytes } if bytes > 1 => Some(b),
        _ => None,
    }
}

macro_rules! convert_buffer {
    ($b:expr, $py:expr, $t:ty, $out:ty) => {
        $b.as_typed::<$t>()?
            .to_vec($py)?
            .into_iter()
            .map(|x| x as $out)
            .collect::<Vec<$out>>()
    };
}

/// Numeric buffer contents converted to `f64` (any float or int dtype).
fn buffer_to_f64(py: Python<'_>, b: &PyUntypedBuffer) -> PyResult<Option<Vec<f64>>> {
    Ok(Some(match ElementType::from_format(b.format()) {
        ElementType::Float { bytes: 8 } => b.as_typed::<f64>()?.to_vec(py)?,
        ElementType::Float { bytes: 4 } => convert_buffer!(b, py, f32, f64),
        ElementType::SignedInteger { bytes: 8 } => convert_buffer!(b, py, i64, f64),
        ElementType::SignedInteger { bytes: 4 } => convert_buffer!(b, py, i32, f64),
        ElementType::SignedInteger { bytes: 2 } => convert_buffer!(b, py, i16, f64),
        ElementType::SignedInteger { bytes: 1 } => convert_buffer!(b, py, i8, f64),
        ElementType::UnsignedInteger { bytes: 8 } => convert_buffer!(b, py, u64, f64),
        ElementType::UnsignedInteger { bytes: 4 } => convert_buffer!(b, py, u32, f64),
        ElementType::UnsignedInteger { bytes: 2 } => convert_buffer!(b, py, u16, f64),
        _ => return Ok(None),
    }))
}

/// Numeric (integer) buffer contents converted to `i64`.
fn buffer_to_i64(py: Python<'_>, b: &PyUntypedBuffer) -> PyResult<Option<Vec<i64>>> {
    Ok(Some(match ElementType::from_format(b.format()) {
        ElementType::SignedInteger { bytes: 8 } => b.as_typed::<i64>()?.to_vec(py)?,
        ElementType::SignedInteger { bytes: 4 } => convert_buffer!(b, py, i32, i64),
        ElementType::SignedInteger { bytes: 2 } => convert_buffer!(b, py, i16, i64),
        ElementType::SignedInteger { bytes: 1 } => convert_buffer!(b, py, i8, i64),
        ElementType::UnsignedInteger { bytes: 4 } => convert_buffer!(b, py, u32, i64),
        ElementType::UnsignedInteger { bytes: 2 } => convert_buffer!(b, py, u16, i64),
        ElementType::UnsignedInteger { bytes: 8 } => b
            .as_typed::<u64>()?
            .to_vec(py)?
            .into_iter()
            .map(|x| x.min(i64::MAX as u64) as i64)
            .collect(),
        _ => return Ok(None),
    }))
}

fn is_typed<T: Element>(b: &PyUntypedBuffer) -> bool {
    b.as_typed::<T>().is_ok() && b.is_c_contiguous()
}

// ---------------------------------------------------------------- f64 arrays

/// A scalar or a 1-D list of numbers: borrowed from a C-contiguous float64 buffer, or owned.
pub(crate) enum F64s {
    Scalar(f64),
    Buf(PyBuffer<f64>),
    Owned(Vec<f64>),
}

impl F64s {
    /// A scalar, a float64 buffer (read in place), any numeric buffer, list, tuple or
    /// iterable of numbers.
    pub(crate) fn read(o: &Bound<'_, PyAny>, what: &str) -> PyResult<F64s> {
        if let Ok(f) = o.cast_exact::<PyFloat>() {
            return Ok(F64s::Scalar(f.value()));
        }
        if o.is_exact_instance_of::<PyInt>() {
            return Ok(F64s::Scalar(o.extract()?));
        }
        Self::read_seq(o, what, true)
    }

    /// As [`Self::read`], but a scalar is an error.
    pub(crate) fn read_array(o: &Bound<'_, PyAny>, what: &str) -> PyResult<F64s> {
        Self::read_seq(o, what, false)
    }

    fn read_seq(o: &Bound<'_, PyAny>, what: &str, scalar_ok: bool) -> PyResult<F64s> {
        let py = o.py();
        if let Ok(l) = o.cast_exact::<PyList>() {
            let mut v = Vec::with_capacity(l.len());
            for k in 0..l.len() {
                v.push(item_f64(&l.get_item(k)?, what, k)?);
            }
            return Ok(F64s::Owned(v));
        }
        if let Ok(t) = o.cast_exact::<PyTuple>() {
            let mut v = Vec::with_capacity(t.len());
            for (k, x) in t.iter().enumerate() {
                v.push(item_f64(&x, what, k)?);
            }
            return Ok(F64s::Owned(v));
        }
        if let Some(b) = numeric_buffer(o) {
            if b.dimensions() > 1 {
                return Err(crate::diag::model_err(format!(
                    "{what}: expected a 1-D array, got {} dimensions",
                    b.dimensions()
                )));
            }
            if b.dimensions() == 0 {
                if scalar_ok && let Some(v) = buffer_to_f64(py, &b)? {
                    return Ok(F64s::Scalar(v.first().copied().unwrap_or(0.0)));
                }
            } else if is_typed::<f64>(&b) {
                return Ok(F64s::Buf(b.into_typed::<f64>()?));
            } else if let Some(v) = buffer_to_f64(py, &b)? {
                return Ok(F64s::Owned(v));
            }
        }
        if scalar_ok && let Some(c) = number(o)? {
            return Ok(F64s::Scalar(c));
        }
        if o.is_instance_of::<PyString>() {
            return Err(PyTypeError::new_err(format!(
                "{what}: expected a number or a sequence of numbers, got str"
            )));
        }
        let it = o.try_iter().map_err(|_| {
            PyTypeError::new_err(format!(
                "{what}: expected a number or a sequence of numbers, got {}",
                type_name(o)
            ))
        })?;
        let mut v = Vec::new();
        for (k, x) in it.enumerate() {
            v.push(item_f64(&x?, what, k)?);
        }
        Ok(F64s::Owned(v))
    }

    pub(crate) fn len(&self) -> Option<usize> {
        match self {
            F64s::Scalar(_) => None,
            F64s::Buf(b) => Some(b.item_count()),
            F64s::Owned(v) => Some(v.len()),
        }
    }

    /// The values (empty for a scalar).
    pub(crate) fn slice(&self) -> &[f64] {
        match self {
            F64s::Scalar(_) => &[],
            F64s::Buf(b) => buf_slice(b),
            F64s::Owned(v) => v,
        }
    }

    pub(crate) fn fill(&self) -> Fill<'_> {
        match self {
            F64s::Scalar(c) => Fill::Scalar(*c),
            _ => Fill::Slice(self.slice()),
        }
    }

    /// Check that this is a scalar or has exactly `n` values; a view for per-item access.
    pub(crate) fn per_item(&self, n: usize, what: &str) -> PyResult<PerItem<'_>> {
        match self {
            F64s::Scalar(c) => Ok(PerItem::Scalar(*c)),
            _ => {
                let s = self.slice();
                if s.len() != n {
                    return Err(crate::diag::model_err(format!(
                        "{what}: expected {n} values, got {}",
                        s.len()
                    )));
                }
                Ok(PerItem::Slice(s))
            }
        }
    }
}

/// A scalar or one value per item.
#[derive(Clone, Copy)]
pub(crate) enum PerItem<'a> {
    Scalar(f64),
    Slice(&'a [f64]),
}

impl PerItem<'_> {
    #[inline]
    pub(crate) fn get(&self, k: usize) -> f64 {
        match self {
            PerItem::Scalar(c) => *c,
            PerItem::Slice(s) => s[k],
        }
    }
}

#[inline]
fn item_f64(x: &Bound<'_, PyAny>, what: &str, k: usize) -> PyResult<f64> {
    if let Ok(f) = x.cast_exact::<PyFloat>() {
        return Ok(f.value());
    }
    number(x)?.ok_or_else(|| {
        PyTypeError::new_err(format!(
            "{what}[{k}]: expected a number, got {}",
            type_name(x)
        ))
    })
}

/// Optional f64 input (None → None).
pub(crate) fn opt_f64s(o: Option<&Bound<'_, PyAny>>, what: &str) -> PyResult<Option<F64s>> {
    o.filter(|o| !o.is_none())
        .map(|o| F64s::read(o, what))
        .transpose()
}

// ---------------------------------------------------------------- bool arrays

/// A flag for all, or one per item.
pub(crate) enum Flags {
    Scalar(bool),
    Each(Vec<bool>),
}

impl Flags {
    pub(crate) fn read(o: &Bound<'_, PyAny>, what: &str) -> PyResult<Flags> {
        if let Ok(b) = o.cast_exact::<pyo3::types::PyBool>() {
            return Ok(Flags::Scalar(b.is_true()));
        }
        if o.is_instance_of::<PyInt>() || o.is_instance_of::<PyFloat>() {
            return Ok(Flags::Scalar(o.is_truthy()?));
        }
        if let Some(b) = numeric_buffer(o)
            && b.dimensions() == 1
            && let Some(v) = buffer_to_f64(o.py(), &b)?
        {
            return Ok(Flags::Each(v.into_iter().map(|x| x != 0.0).collect()));
        }
        if o.is_instance_of::<PyString>() {
            return Err(PyTypeError::new_err(format!(
                "{what}: expected a bool or a sequence of bool"
            )));
        }
        match o.try_iter() {
            Ok(it) => {
                let mut v = Vec::new();
                for x in it {
                    v.push(x?.is_truthy()?);
                }
                Ok(Flags::Each(v))
            }
            // numpy bool scalars and other truthy objects
            Err(_) => Ok(Flags::Scalar(o.is_truthy()?)),
        }
    }
}

// ---------------------------------------------------------------- integer arrays

/// A 1-D list of integers: borrowed from a C-contiguous int32/int64 buffer, or owned.
pub(crate) enum Ints {
    B32(PyBuffer<i32>),
    B64(PyBuffer<i64>),
    V32(Vec<i32>),
    V64(Vec<i64>),
}

/// The integers of an [`Ints`].
#[derive(Clone, Copy)]
pub(crate) enum IntSlice<'a> {
    I32(&'a [i32]),
    I64(&'a [i64]),
}

impl IntSlice<'_> {
    pub(crate) fn len(&self) -> usize {
        match self {
            IntSlice::I32(s) => s.len(),
            IntSlice::I64(s) => s.len(),
        }
    }
    #[inline]
    pub(crate) fn get(&self, k: usize) -> i64 {
        match self {
            IntSlice::I32(s) => s[k] as i64,
            IntSlice::I64(s) => s[k],
        }
    }
}

impl Ints {
    /// Integers from a buffer, a list/tuple/iterable of ints, Vars (of `model`, if given), a
    /// `range` or a `Vars` block.
    pub(crate) fn read(o: &Bound<'_, PyAny>, what: &str, model: Option<u64>) -> PyResult<Ints> {
        let item = |x: &Bound<'_, PyAny>, k: usize| -> PyResult<i64> {
            if let Some(m) = model {
                return col_of(x, m);
            }
            if (x.is_instance_of::<PyInt>() || !x.is_instance_of::<PyFloat>())
                && let Ok(j) = x.extract::<i64>()
            {
                return Ok(j);
            }
            Err(PyTypeError::new_err(format!(
                "{what}[{k}]: expected an int, got {}",
                type_name(x)
            )))
        };
        if let Ok(l) = o.cast_exact::<PyList>() {
            let mut v = Vec::with_capacity(l.len());
            for k in 0..l.len() {
                v.push(item(&l.get_item(k)?, k)?);
            }
            return Ok(Ints::V64(v));
        }
        if let Ok(t) = o.cast_exact::<PyTuple>() {
            let mut v = Vec::with_capacity(t.len());
            for (k, x) in t.iter().enumerate() {
                v.push(item(&x, k)?);
            }
            return Ok(Ints::V64(v));
        }
        if let Ok(vs) = o.cast_exact::<Vars>() {
            let vs = vs.get();
            if let Some(m) = model
                && vs.model != m
            {
                return Err(wrong_model());
            }
            return Ok(Ints::V32((0..vs.len).map(|k| vs.col(k)).collect()));
        }
        if let Ok(r) = o.cast_exact::<PyRange>() {
            let (start, stop, step): (i64, i64, i64) = (
                r.getattr("start")?.extract()?,
                r.getattr("stop")?.extract()?,
                r.getattr("step")?.extract()?,
            );
            let n = r.len()?;
            let _ = stop;
            return Ok(Ints::V64((0..n as i64).map(|k| start + k * step).collect()));
        }
        if let Some(b) = numeric_buffer(o) {
            if b.dimensions() != 1 {
                return Err(crate::diag::model_err(format!(
                    "{what}: expected a 1-D integer array"
                )));
            }
            if is_typed::<i32>(&b) {
                return Ok(Ints::B32(b.into_typed()?));
            }
            if is_typed::<i64>(&b) {
                return Ok(Ints::B64(b.into_typed()?));
            }
            if let Some(v) = buffer_to_i64(o.py(), &b)? {
                return Ok(Ints::V64(v));
            }
            return Err(PyTypeError::new_err(format!(
                "{what}: expected an integer array, got a float array"
            )));
        }
        if o.is_instance_of::<PyString>() {
            return Err(PyTypeError::new_err(format!(
                "{what}: expected a sequence of ints"
            )));
        }
        let it = o.try_iter().map_err(|_| {
            PyTypeError::new_err(format!(
                "{what}: expected a sequence of ints, got {}",
                type_name(o)
            ))
        })?;
        let mut v = Vec::new();
        for (k, x) in it.enumerate() {
            v.push(item(&x?, k)?);
        }
        Ok(Ints::V64(v))
    }

    pub(crate) fn slice(&self) -> IntSlice<'_> {
        match self {
            Ints::B32(b) => IntSlice::I32(buf_slice(b)),
            Ints::B64(b) => IntSlice::I64(buf_slice(b)),
            Ints::V32(v) => IntSlice::I32(v),
            Ints::V64(v) => IntSlice::I64(v),
        }
    }

    /// As `i32` column indices (a copy unless already int32), each checked against
    /// `0..limit`.
    pub(crate) fn to_cols(&self, limit: usize, what: &str) -> PyResult<Vec<i32>> {
        let s = self.slice();
        let mut v = Vec::with_capacity(s.len());
        for k in 0..s.len() {
            let j = s.get(k);
            if j < 0 || j as u64 >= limit as u64 {
                return Err(crate::diag::index_err(format!(
                    "{what}[{k}] = {j} is out of range (there are {limit})"
                )));
            }
            v.push(j as i32);
        }
        Ok(v)
    }
}

/// Indices of variables: a Var, an int, a Vars, or a sequence/buffer of Vars or ints, each
/// checked against `0..limit`. A contiguous `Vars` is returned as a range without a list.
pub(crate) enum Cols {
    Range(Range<usize>),
    List(Vec<i32>),
}

impl Cols {
    pub(crate) fn read(
        o: &Bound<'_, PyAny>,
        model: u64,
        limit: usize,
        what: &str,
    ) -> PyResult<Cols> {
        if let Ok(vs) = o.cast_exact::<Vars>() {
            let vs = vs.get();
            if vs.model != model {
                return Err(wrong_model());
            }
            if let Some(r) = vs.as_range() {
                return Ok(Cols::Range(r));
            }
        }
        if o.is_instance_of::<Var>() || o.is_instance_of::<PyInt>() {
            let j = col_of(o, model)?;
            if j < 0 || j as u64 >= limit as u64 {
                return Err(crate::diag::index_err(format!(
                    "{what} {j} is out of range (there are {limit})"
                )));
            }
            return Ok(Cols::List(vec![j as i32]));
        }
        Ok(Cols::List(
            Ints::read(o, what, Some(model))?.to_cols(limit, what)?,
        ))
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Cols::Range(r) => r.len(),
            Cols::List(v) => v.len(),
        }
    }

    #[inline]
    pub(crate) fn get(&self, k: usize) -> usize {
        match self {
            Cols::Range(r) => r.start + k,
            Cols::List(v) => v[k] as usize,
        }
    }
}

// ---------------------------------------------------------------- rows

/// What a row object said about its bounds.
pub(crate) enum RowKind {
    /// A plain row (dict, pair, LinExpr, Var): its constant (moved to the bounds).
    Plain(f64),
    /// A Constraint: its own bounds.
    Bounded(f64, f64),
}

/// Push the entries of `row` into the builder's open row. On error the open row is aborted.
pub(crate) fn push_row(
    b: &mut ModelBuilder,
    row: &Bound<'_, PyAny>,
    model: u64,
) -> PyResult<RowKind> {
    let r = push_row_inner(b, row, model);
    if r.is_err() {
        b.abort_row();
    }
    r
}

fn push_row_inner(b: &mut ModelBuilder, row: &Bound<'_, PyAny>, model: u64) -> PyResult<RowKind> {
    if let Ok(d) = row.cast_exact::<PyDict>() {
        push_dict(b, d, model)?;
        return Ok(RowKind::Plain(0.0));
    }
    if let Ok(c) = row.cast_exact::<Constraint>() {
        let c = c.get();
        check_model(c.model, model)?;
        c.body.with(|t| push_terms(b, t))?;
        return Ok(RowKind::Bounded(c.lo, c.hi));
    }
    if let Ok(e) = row.cast_exact::<LinExpr>() {
        let e = e.try_borrow()?;
        check_model(e.model, model)?;
        e.terms.with(|t| push_terms(b, t))?;
        return Ok(RowKind::Plain(e.constant));
    }
    if let Ok(v) = row.cast_exact::<Var>() {
        let v = v.get();
        check_model(v.model, model)?;
        b.push(v.index as i64, 1.0)
            .map_err(|e| build_err(e, None))?;
        return Ok(RowKind::Plain(0.0));
    }
    if let Ok(t) = row.cast_exact::<PyTuple>()
        && t.len() == 2
    {
        push_pair(b, &t.get_item(0)?, &t.get_item(1)?, model)?;
        return Ok(RowKind::Plain(0.0));
    }
    if let Ok(d) = row.cast::<PyDict>() {
        push_dict(b, d, model)?;
        return Ok(RowKind::Plain(0.0));
    }
    Err(PyTypeError::new_err(format!(
        "a row must be a dict {{var: coef}}, an (indices, values) pair, a LinExpr, a Var or a \
         Constraint, got {}",
        type_name(row)
    )))
}

#[inline]
fn check_model(got: u64, model: u64) -> PyResult<()> {
    if got != 0 && got != model {
        Err(wrong_model())
    } else {
        Ok(())
    }
}

#[inline]
fn push_terms(b: &mut ModelBuilder, t: &[(i32, f64)]) -> PyResult<()> {
    for &(j, c) in t {
        b.push(j as i64, c).map_err(|e| build_err(e, None))?;
    }
    Ok(())
}

/// The hot path of dict rows: one raw `PyDict_Next` per entry with borrowed references, and
/// exact-type checks by type pointer for the common `int -> float|int` entry.
///
/// Why raw FFI: profiling the pyo3 iterator on the 19k-nonzero production row set
/// showed ~45% of build time in
/// `Py_IncRef`/`Py_DecRef` (out-of-line under abi3) plus a `PyDict_Size` call per entry, and int
/// coefficients going through a temporary float object. `PyDict_Next` hands out borrowed
/// references (CPython docs, "Dictionary Objects"), so the fast path touches no refcounts.
#[inline]
fn push_dict(b: &mut ModelBuilder, d: &Bound<'_, PyDict>, model: u64) -> PyResult<()> {
    let py = d.py();
    let dp = d.as_ptr();
    let mut pos: ffi::Py_ssize_t = 0;
    let mut k: *mut ffi::PyObject = std::ptr::null_mut();
    let mut v: *mut ffi::PyObject = std::ptr::null_mut();
    // SAFETY: `dp` is a live dict (we hold `d` and the GIL). PyDict_Next stores borrowed
    // references to the entry's key and value, valid while the dict is not mutated. The fast
    // path calls only PyLong_AsLongLong/PyFloat_AsDouble/PyLong_AsDouble on exact int/float
    // objects, which run no Python code, so the dict cannot change under us. Every path that
    // can run Python code (`col_of`, `coef` on other types) first takes owned references
    // (`to_owned`), so a mutation there cannot free the objects it is using; PyDict_Next itself
    // bounds-checks `pos` against the current table, so iteration stays memory-safe.
    while unsafe { ffi::PyDict_Next(dp, &mut pos, &mut k, &mut v) } != 0 {
        let j = if unsafe { ffi::PyLong_CheckExact(k) } != 0 {
            let j = unsafe { ffi::PyLong_AsLongLong(k) };
            if j == -1 && PyErr::take(py).is_some() {
                // out of i64 range: the generic path raises the usual error
                col_of(&unsafe { Borrowed::from_ptr(py, k) }.to_owned(), model)?
            } else {
                j
            }
        } else {
            let key = unsafe { Borrowed::from_ptr(py, k) }.to_owned();
            if let Ok(var) = key.cast_exact::<Var>() {
                let var = var.get();
                if var.model != model {
                    return Err(wrong_model());
                }
                var.index as i64
            } else {
                col_of(&key, model)?
            }
        };
        let c = if unsafe { ffi::PyFloat_CheckExact(v) } != 0 {
            unsafe { ffi::PyFloat_AsDouble(v) }
        } else if unsafe { ffi::PyLong_CheckExact(v) } != 0 {
            let c = unsafe { ffi::PyLong_AsDouble(v) };
            if c == -1.0
                && let Some(e) = PyErr::take(py)
            {
                return Err(e);
            }
            c
        } else {
            coef(&unsafe { Borrowed::from_ptr(py, v) }.to_owned())?
        };
        b.push(j, c).map_err(|e| build_err(e, None))?;
    }
    Ok(())
}

fn push_pair(
    b: &mut ModelBuilder,
    idx: &Bound<'_, PyAny>,
    val: &Bound<'_, PyAny>,
    model: u64,
) -> PyResult<()> {
    let idx = Ints::read(idx, "indices", Some(model))?;
    let val = F64s::read_array(val, "values")?;
    let (is, vs) = (idx.slice(), val.slice());
    if is.len() != vs.len() {
        return Err(crate::diag::model_err(format!(
            "row (indices, values): {} indices but {} values",
            is.len(),
            vs.len()
        )));
    }
    for (k, &v) in vs.iter().enumerate() {
        b.push(is.get(k), v).map_err(|e| build_err(e, None))?;
    }
    Ok(())
}
