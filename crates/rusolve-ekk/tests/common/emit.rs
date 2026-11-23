//! A JSON builder that produces records in exactly the oracle's encoding
//! (`tools/highs-oracle/src/emit.{h,cpp}`), so actual values from the Rust
//! port can be diffed against fixture records with [`super::check::diff_values`].
//!
//! Method names follow emit.h: `i` (integer), `b` (bool as 0/1), `f` (f64 as
//! 16 hex digits of the bit pattern), `h` (hash as 16 hex digits), `s`
//! (string), `vi` / `vi8` / `vf` (vectors), `o` (sub-object), `hv` (HVector),
//! `mat` (HighsSparseMatrix).

use serde_json::{Map, Value};

use super::records::{HVectorRec, MatrixRec};
use super::{f64_hex, hash::hash_hex};

/// Borrowed view of the `pack*` part of an HVector.
#[derive(Clone, Copy, Debug)]
pub struct PackView<'a> {
    pub pack_flag: bool,
    pub pack_count: i32,
    /// At least `pack_count` long; only `[0..pack_count)` is emitted.
    pub pack_index: &'a [i32],
    /// At least `pack_count` long; only `[0..pack_count)` is emitted.
    pub pack_value: &'a [f64],
}

/// Borrowed view of an HVector, enough to emit it as the oracle does.
#[derive(Clone, Copy, Debug)]
pub struct HVectorView<'a> {
    pub size: i32,
    pub count: i32,
    /// The whole `index` vector (HiGHS sizes it to `size`); the oracle emits
    /// `index[0..count)` only when `0 <= count <= index.len()`.
    pub index: &'a [i32],
    /// The whole `array` (length `size`).
    pub array: &'a [f64],
    pub synthetic_tick: f64,
    pub pack: Option<PackView<'a>>,
}

/// Borrowed view of a HighsSparseMatrix.
#[derive(Clone, Copy, Debug)]
pub struct MatrixView<'a> {
    /// 1 colwise, 2 rowwise, 3 rowwise partitioned.
    pub format: i32,
    pub num_col: i32,
    pub num_row: i32,
    pub start: &'a [i32],
    /// Emitted only when `format == 3`.
    pub p_end: Option<&'a [i32]>,
    pub index: &'a [i32],
    pub value: &'a [f64],
}

impl MatrixRec {
    pub fn view(&self) -> MatrixView<'_> {
        MatrixView {
            format: self.format,
            num_col: self.num_col,
            num_row: self.num_row,
            start: &self.start,
            p_end: self.p_end.as_deref(),
            index: &self.index,
            value: &self.value,
        }
    }
}

/// JSON object builder (by value, chainable). Key order is irrelevant for
/// comparison; values are encoded exactly like the oracle.
#[derive(Clone, Debug, Default)]
pub struct Obj {
    m: Map<String, Value>,
}

impl Obj {
    pub fn new() -> Self {
        Obj { m: Map::new() }
    }

    /// A record with field `t` set (`Obj(const char* type)` in emit.h).
    pub fn typed(t: &str) -> Self {
        Obj::new().s("t", t)
    }

    pub fn i(mut self, k: &str, v: impl Into<i64>) -> Self {
        self.m.insert(k.into(), Value::from(v.into()));
        self
    }

    pub fn b(self, k: &str, v: bool) -> Self {
        self.i(k, i64::from(v))
    }

    pub fn f(mut self, k: &str, v: f64) -> Self {
        self.m.insert(k.into(), Value::String(f64_hex(v)));
        self
    }

    /// Emit `k` only when `v` is `Some` (fields the oracle writes
    /// conditionally, e.g. `edge_weight_row_out_after`, `edge_weight_in`).
    pub fn opt_f(self, k: &str, v: Option<f64>) -> Self {
        match v {
            Some(x) => self.f(k, x),
            None => self,
        }
    }

    /// Emit `k` only when `v` is `Some` (primal `iter_none` pivot fields).
    pub fn opt_i(self, k: &str, v: Option<i32>) -> Self {
        match v {
            Some(x) => self.i(k, x),
            None => self,
        }
    }

    /// Emit `k` only when `v` is `Some` (`dual_start`'s
    /// `allow_dual_steepest_edge_to_devex_switch`).
    pub fn opt_b(self, k: &str, v: Option<bool>) -> Self {
        match v {
            Some(x) => self.b(k, x),
            None => self,
        }
    }

    pub fn h(mut self, k: &str, v: u64) -> Self {
        self.m.insert(k.into(), Value::String(hash_hex(v)));
        self
    }

    pub fn s(mut self, k: &str, v: &str) -> Self {
        self.m.insert(k.into(), Value::String(v.into()));
        self
    }

    pub fn vi(mut self, k: &str, v: &[i32]) -> Self {
        self.m.insert(
            k.into(),
            Value::Array(v.iter().map(|&x| Value::from(x)).collect()),
        );
        self
    }

    pub fn vi8(mut self, k: &str, v: &[i8]) -> Self {
        self.m.insert(
            k.into(),
            Value::Array(v.iter().map(|&x| Value::from(x)).collect()),
        );
        self
    }

    pub fn vf(mut self, k: &str, v: &[f64]) -> Self {
        self.m.insert(
            k.into(),
            Value::Array(v.iter().map(|&x| Value::String(f64_hex(x))).collect()),
        );
        self
    }

    pub fn o(mut self, k: &str, sub: Obj) -> Self {
        self.m.insert(k.into(), sub.build());
        self
    }

    /// Insert an already-built JSON value (e.g. the raw `options` object).
    pub fn raw(mut self, k: &str, v: Value) -> Self {
        self.m.insert(k.into(), v);
        self
    }

    pub fn list(mut self, k: &str, items: Vec<Obj>) -> Self {
        self.m.insert(
            k.into(),
            Value::Array(items.into_iter().map(Obj::build).collect()),
        );
        self
    }

    /// An HVector exactly as `Obj::hv` in emit.cpp: `size`, `count`,
    /// `index[0..count)` (only if `0 <= count <= index.len()`),
    /// `synthetic_tick`, `nz_i`/`nz_v` (every position whose **bit pattern**
    /// is non-zero, ascending), and with `with_pack` also
    /// `packFlag`/`packCount`/`packIndex`/`packValue`.
    pub fn hv(self, k: &str, v: &HVectorView, with_pack: bool) -> Self {
        self.o(k, hv_obj(v, with_pack))
    }

    /// Re-emit a decoded fixture HVector (used for round-trip checks and to
    /// copy expected inputs into actual records).
    pub fn hv_rec(self, k: &str, v: &HVectorRec) -> Self {
        let array = v.array();
        let index = v.index.clone().unwrap_or_default();
        let pack = v.pack.as_ref().map(|p| PackView {
            pack_flag: p.pack_flag,
            pack_count: p.pack_count,
            pack_index: &p.pack_index,
            pack_value: &p.pack_value,
        });
        let view = HVectorView {
            size: v.size,
            count: if v.index.is_some() { v.count } else { -1 },
            index: &index,
            array: &array,
            synthetic_tick: v.synthetic_tick,
            pack,
        };
        let mut o = hv_obj(&view, v.pack.is_some());
        if v.index.is_none() {
            // Dense mode: restore the original negative count.
            o = o.i("count", v.count);
        }
        self.o(k, o)
    }

    /// A HighsSparseMatrix exactly as `Obj::mat` in emit.cpp.
    pub fn mat(self, k: &str, m: &MatrixView) -> Self {
        let mut o = Obj::new()
            .i("format", m.format)
            .i("num_col", m.num_col)
            .i("num_row", m.num_row)
            .vi("start", m.start);
        if m.format == 3 {
            let p_end = m
                .p_end
                .expect("MatrixView: format 3 (rowwise partitioned) needs p_end");
            o = o.vi("p_end", p_end);
        }
        self.o(k, o.vi("index", m.index).vf("value", m.value))
    }

    pub fn build(self) -> Value {
        Value::Object(self.m)
    }
}

fn hv_obj(v: &HVectorView, with_pack: bool) -> Obj {
    let mut o = Obj::new().i("size", v.size).i("count", v.count);
    if v.count >= 0 && (v.count as usize) <= v.index.len() {
        o = o.vi("index", &v.index[..v.count as usize]);
    }
    o = o.f("synthetic_tick", v.synthetic_tick);
    let mut nz_i = Vec::new();
    let mut nz_v = Vec::new();
    for (j, &x) in v.array.iter().enumerate() {
        if x.to_bits() != 0 {
            nz_i.push(j as i32);
            nz_v.push(x);
        }
    }
    o = o.vi("nz_i", &nz_i).vf("nz_v", &nz_v);
    if with_pack {
        let p = v
            .pack
            .expect("HVectorView: with_pack requested but pack is None");
        let n = p.pack_count.max(0) as usize;
        o = o
            .b("packFlag", p.pack_flag)
            .i("packCount", p.pack_count)
            .vi("packIndex", &p.pack_index[..n])
            .vf("packValue", &p.pack_value[..n]);
    }
    o
}
