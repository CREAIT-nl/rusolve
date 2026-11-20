//! Golden test for manifest row 1, `src/util/highs_cdouble.rs`
//! (HighsCDouble.h).
//!
//! There is no oracle fixture for HighsCDouble (signature/fixture gap G-CD:
//! the oracle dumps its effect only through the `quad = 1` PRICE records of
//! `t_sparse_matrix`). This test therefore compares the port, bit for bit
//! (`hi` and `lo` separately), with `Ref`, a line-by-line transcription of
//! HighsCDouble.h (HiGHS 1.15.1, MIT License, Copyright (cc) HiGHS; see
//! THIRD_PARTY_NOTICES.md). Rust f64 arithmetic is IEEE-754 binary64 with
//! round-to-nearest and no contraction, which is what the oracle's
//! `-ffp-contract=off` build gives the C++ code, so `Ref` reproduces HiGHS.
//!
//! Inputs: a fixed set of doubles (signed zeros, units, non-representable
//! decimals, huge/tiny/subnormal magnitudes, near-ties) and compensated
//! values built from them, so `lo != 0` paths are exercised. Every operator
//! of the header is checked over all (lhs, rhs) pairs.

mod common;
mod port;

use rusolve_ekk::prelude::*;

/// Transcription of `class HighsCDouble` (HighsCDouble.h). Each method
/// names its C++ counterpart.
#[derive(Clone, Copy, Debug)]
struct Ref {
    hi: f64,
    lo: f64,
}

impl Ref {
    fn two_sum(x: &mut f64, y: &mut f64, a: f64, b: f64) {
        *x = a + b;
        let z = *x - a;
        *y = (a - (*x - z)) + (b - z);
    }
    fn split(x: &mut f64, y: &mut f64, a: f64) {
        let factor = ((1i64 << 27) + 1) as f64;
        let cc = factor * a;
        *x = cc - (cc - a);
        *y = a - *x;
    }
    fn two_product(x: &mut f64, y: &mut f64, a: f64, b: f64) {
        *x = a * b;
        let (mut a1, mut a2, mut b1, mut b2) = (0.0, 0.0, 0.0, 0.0);
        Self::split(&mut a1, &mut a2, a);
        Self::split(&mut b1, &mut b2, b);
        *y = a2 * b2 - (((*x - a1 * b1) - a2 * b1) - a1 * b2);
    }
    fn new(hi: f64, lo: f64) -> Ref {
        Ref { hi, lo }
    }
    /// `HighsCDouble(double val)`
    fn from(v: f64) -> Ref {
        Ref { hi: v, lo: 0.0 }
    }
    /// `explicit operator double()`
    fn dbl(self) -> f64 {
        self.hi + self.lo
    }
    /// `operator+=(double)`
    fn add_assign_d(&mut self, v: f64) {
        let mut cc = 0.0;
        let hi = self.hi;
        Self::two_sum(&mut self.hi, &mut cc, v, hi);
        self.lo += cc;
    }
    /// `operator+=(const HighsCDouble&)`
    fn add_assign(&mut self, v: Ref) {
        self.add_assign_d(v.hi);
        self.lo += v.lo;
    }
    /// `operator-=(double)`
    fn sub_assign_d(&mut self, v: f64) {
        self.add_assign_d(-v);
    }
    /// `operator-=(const HighsCDouble&)`
    fn sub_assign(&mut self, v: Ref) {
        self.sub_assign_d(v.hi);
        self.lo -= v.lo;
    }
    /// `operator*=(double)`
    fn mul_assign_d(&mut self, v: f64) {
        let cc = self.lo * v;
        let hi = self.hi;
        Self::two_product(&mut self.hi, &mut self.lo, hi, v);
        self.add_assign_d(cc);
    }
    /// `operator*=(const HighsCDouble&)`
    fn mul_assign(&mut self, v: Ref) {
        let c1 = self.hi * v.lo;
        let c2 = self.lo * v.hi;
        let hi = self.hi;
        Self::two_product(&mut self.hi, &mut self.lo, hi, v.hi);
        self.add_assign_d(c1);
        self.add_assign_d(c2);
    }
    /// `operator/=(double)`
    fn div_assign_d(&mut self, v: f64) {
        let d = Ref::new(self.hi / v, self.lo / v);
        let mut cc = d.mul_d(v).sub(*self);
        cc.hi /= v;
        cc.lo /= v;
        *self = d.sub(cc);
    }
    /// `operator/=(const HighsCDouble&)`
    fn div_assign(&mut self, v: Ref) {
        let vdbl = v.hi + v.lo;
        let d = Ref::new(self.hi / vdbl, self.lo / vdbl);
        let mut cc = d.mul(v).sub(*self);
        cc.hi /= vdbl;
        cc.lo /= vdbl;
        *self = d.sub(cc);
    }
    /// unary `operator-`
    fn neg(self) -> Ref {
        Ref::new(-self.hi, -self.lo)
    }
    /// `operator+(double) const`
    fn add_d(self, v: f64) -> Ref {
        let mut res = Ref::new(0.0, 0.0);
        Self::two_sum(&mut res.hi, &mut res.lo, self.hi, v);
        res.lo += self.lo;
        res
    }
    /// `operator+(const HighsCDouble&) const`
    fn add(self, v: Ref) -> Ref {
        let mut res = self.add_d(v.hi);
        res.lo += v.lo;
        res
    }
    /// `operator-(double) const`
    fn sub_d(self, v: f64) -> Ref {
        let mut res = Ref::new(0.0, 0.0);
        Self::two_sum(&mut res.hi, &mut res.lo, self.hi, -v);
        res.lo += self.lo;
        res
    }
    /// `operator-(const HighsCDouble&) const`
    fn sub(self, v: Ref) -> Ref {
        let mut res = self.sub_d(v.hi);
        res.lo -= v.lo;
        res
    }
    /// `operator*(double) const`
    fn mul_d(self, v: f64) -> Ref {
        let mut res = Ref::new(0.0, 0.0);
        Self::two_product(&mut res.hi, &mut res.lo, self.hi, v);
        res.add_assign_d(self.lo * v);
        res
    }
    /// `operator*(const HighsCDouble&) const`
    fn mul(self, v: Ref) -> Ref {
        let mut res = self.mul_d(v.hi);
        res.add_assign_d(self.hi * v.lo);
        res
    }
    /// `operator/(double) const`
    fn div_d(self, v: f64) -> Ref {
        let mut res = self;
        res.div_assign_d(v);
        res
    }
    /// `operator/(const HighsCDouble&) const`
    fn div(self, v: Ref) -> Ref {
        let mut res = self;
        res.div_assign(v);
        res
    }
    /// `renormalize()`
    fn renormalize(&mut self) {
        let (hi, lo) = (self.hi, self.lo);
        Self::two_sum(&mut self.hi, &mut self.lo, hi, lo);
    }
    /// `abs(const HighsCDouble&)`: `v < 0 ? -v : v`
    fn abs(self) -> Ref {
        if self.dbl() < 0.0 { self.neg() } else { self }
    }
    /// `sqrt(const HighsCDouble&)`
    fn sqrt(self) -> Ref {
        let cc = (self.hi + self.lo).sqrt();
        if cc == 0.0 {
            return Ref::from(0.0);
        }
        let mut res = self.div_d(cc);
        res.add_assign_d(cc);
        res.hi *= 0.5;
        res.lo *= 0.5;
        res
    }
    /// `floor(const HighsCDouble&)`
    fn floor(self) -> Ref {
        if self.abs().dbl() < 1.0 {
            if self.dbl() == 0.0 || self.dbl() > 0.0 {
                return Ref::from(0.0);
            }
            return Ref::from(-1.0);
        }
        let floor_x = self.dbl().floor();
        let mut res = Ref::new(0.0, 0.0);
        Self::two_sum(
            &mut res.hi,
            &mut res.lo,
            floor_x,
            self.sub_d(floor_x).dbl().floor(),
        );
        res
    }
    /// `ceil(const HighsCDouble&)`
    fn ceil(self) -> Ref {
        if self.abs().dbl() < 1.0 {
            if self.dbl() == 0.0 || self.dbl() < 0.0 {
                return Ref::from(0.0);
            }
            return Ref::from(1.0);
        }
        let ceil_x = self.dbl().ceil();
        let mut res = Ref::new(0.0, 0.0);
        Self::two_sum(
            &mut res.hi,
            &mut res.lo,
            ceil_x,
            self.sub_d(ceil_x).dbl().ceil(),
        );
        res
    }
    /// `round(const HighsCDouble&)`: `floor(x + 0.5)`
    fn round(self) -> Ref {
        self.add_d(0.5).floor()
    }
    /// `ldexp(const HighsCDouble&, int)`
    fn ldexp(self, exp: i32) -> Ref {
        Ref::new(ldexp(self.hi, exp), ldexp(self.lo, exp))
    }
}

/// `std::ldexp` for the exponents used here (|exp| small, results normal
/// or exactly representable): exact scaling by a power of two.
fn ldexp(x: f64, exp: i32) -> f64 {
    let mut r = x;
    let mut e = exp;
    while e > 0 {
        let s = e.min(1000);
        r *= 2f64.powi(s);
        e -= s;
    }
    while e < 0 {
        let s = (-e).min(1000);
        r *= 2f64.powi(-s);
        e += s;
    }
    r
}

fn port_of(r: Ref) -> HighsCDouble {
    let mut p = HighsCDouble::from(0.0);
    p.hi = r.hi;
    p.lo = r.lo;
    p
}

fn bits(x: f64) -> String {
    format!("{:016x} ({x:e})", x.to_bits())
}

struct Check {
    n: usize,
    fails: Vec<String>,
}

impl Check {
    fn same(&mut self, what: &str, e: Ref, a: HighsCDouble) {
        self.n += 1;
        if e.hi.to_bits() != a.hi.to_bits() || e.lo.to_bits() != a.lo.to_bits() {
            self.fails.push(format!(
                "{what}\n    expected hi={} lo={}\n    actual   hi={} lo={}",
                bits(e.hi),
                bits(e.lo),
                bits(a.hi),
                bits(a.lo)
            ));
        }
    }
    fn same_f64(&mut self, what: &str, e: f64, a: f64) {
        self.n += 1;
        if e.to_bits() != a.to_bits() {
            self.fails.push(format!(
                "{what}\n    expected {}\n    actual   {}",
                bits(e),
                bits(a)
            ));
        }
    }
    fn same_bool(&mut self, what: &str, e: bool, a: bool) {
        self.n += 1;
        if e != a {
            self.fails
                .push(format!("{what}\n    expected {e}\n    actual   {a}"));
        }
    }
    #[track_caller]
    fn finish(&self, what: &str) {
        assert!(self.n > 0);
        if self.fails.is_empty() {
            return;
        }
        let shown: Vec<_> = self.fails.iter().take(20).cloned().collect();
        panic!(
            "{what}: {} of {} comparisons differ from HighsCDouble.h (first {} shown)\n{}",
            self.fails.len(),
            self.n,
            shown.len(),
            shown.join("\n")
        );
    }
}

// Some literals carry more digits than f64 holds on purpose: the compiler
// rounds them, giving inputs with long mantissas.
#[allow(clippy::excessive_precision)]
fn doubles() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        1.0,
        -1.0,
        0.5,
        -0.5,
        0.1,
        -0.3,
        1.0 / 3.0,
        2.0 / 3.0,
        std::f64::consts::PI,
        -std::f64::consts::E,
        3.0,
        7.25,
        -2.5,
        1e16,
        -1e16 + 3.0,
        9007199254740993.0, // 2^53 + 1 (rounds to 2^53)
        1e-16,
        -7e-17,
        1e300,
        -1e-300,
        f64::MIN_POSITIVE,
        5e-324,
        0.999_999_999_999_999_9,
        -1.000_000_000_000_000_2,
        123456789.123456789,
        -0.75,
        2.0f64.powi(27) + 1.0,
        1e-8,
    ]
}

/// Compensated values: `HighsCDouble(a) + b` for a spread of pairs, plus
/// plain conversions.
fn cdoubles() -> Vec<Ref> {
    let d = doubles();
    let mut out: Vec<Ref> = d.iter().map(|&x| Ref::from(x)).collect();
    for (i, &a) in d.iter().enumerate() {
        let b = d[(i * 7 + 3) % d.len()];
        out.push(Ref::from(a).add_d(b * 1e-17));
        out.push(Ref::from(a).add_d(b));
    }
    out.push(Ref::from(1.0).add_d(1e-20).mul_d(3.0));
    out.push(Ref::from(0.1).mul_d(0.1));
    out
}

fn port_ops_match_reference(cc: &mut Check) {
    let ds = doubles();
    let cs = cdoubles();
    for (i, &x) in cs.iter().enumerate() {
        let px = port_of(x);
        let tag = format!("x=cs[{i}](hi={}, lo={})", bits(x.hi), bits(x.lo));
        cc.same_f64(&format!("double({tag})"), x.dbl(), f64::from(px));
        cc.same(&format!("-{tag}"), x.neg(), -px);
        cc.same(&format!("abs({tag})"), x.abs(), px.abs());
        if x.dbl() >= 0.0 {
            cc.same(&format!("sqrt({tag})"), x.sqrt(), px.sqrt());
        }
        if x.dbl().abs() < 1e290 {
            cc.same(&format!("floor({tag})"), x.floor(), px.floor());
            cc.same(&format!("ceil({tag})"), x.ceil(), px.ceil());
            cc.same(&format!("round({tag})"), x.round(), px.round());
        }
        for e in [-60, -1, 0, 1, 3, 52] {
            cc.same(&format!("ldexp({tag}, {e})"), x.ldexp(e), px.ldexp(e));
        }
        let mut r = x;
        r.renormalize();
        let mut pr = px;
        pr.renormalize();
        cc.same(&format!("renormalize({tag})"), r, pr);

        for (j, &v) in ds.iter().enumerate() {
            let t = format!("{tag}, v=ds[{j}]={}", bits(v));
            cc.same(&format!("x + v [{t}]"), x.add_d(v), px + v);
            cc.same(&format!("v + x [{t}]"), x.add_d(v), v + px);
            cc.same(&format!("x - v [{t}]"), x.sub_d(v), px - v);
            cc.same(&format!("v - x [{t}]"), x.neg().add_d(v), v - px);
            cc.same(&format!("x * v [{t}]"), x.mul_d(v), px * v);
            cc.same(&format!("v * x [{t}]"), x.mul_d(v), v * px);
            let mut e = x;
            e.add_assign_d(v);
            let mut a = px;
            a += v;
            cc.same(&format!("x += v [{t}]"), e, a);
            let mut e = x;
            e.sub_assign_d(v);
            let mut a = px;
            a -= v;
            cc.same(&format!("x -= v [{t}]"), e, a);
            let mut e = x;
            e.mul_assign_d(v);
            let mut a = px;
            a *= v;
            cc.same(&format!("x *= v [{t}]"), e, a);
            if v != 0.0 && v.abs() > 1e-290 {
                cc.same(&format!("x / v [{t}]"), x.div_d(v), px / v);
                let mut e = x;
                e.div_assign_d(v);
                let mut a = px;
                a /= v;
                cc.same(&format!("x /= v [{t}]"), e, a);
            }
            if x.dbl() != 0.0 && x.dbl().abs() > 1e-290 {
                cc.same(&format!("v / x [{t}]"), Ref::from(v).div(x), v / px);
            }
            cc.same_bool(&format!("x < v [{t}]"), x.dbl() < v, px < v);
            cc.same_bool(&format!("x > v [{t}]"), x.dbl() > v, px > v);
            cc.same_bool(&format!("x == v [{t}]"), x.dbl() == v, px == v);
            cc.same_bool(&format!("v < x [{t}]"), v < x.dbl(), v < px);
        }
        for (j, &y) in cs.iter().enumerate().step_by(3) {
            let py = port_of(y);
            let t = format!("{tag}, y=cs[{j}](hi={}, lo={})", bits(y.hi), bits(y.lo));
            cc.same(&format!("x + y [{t}]"), x.add(y), px + py);
            cc.same(&format!("x - y [{t}]"), x.sub(y), px - py);
            cc.same(&format!("x * y [{t}]"), x.mul(y), px * py);
            let mut e = x;
            e.add_assign(y);
            let mut a = px;
            a += py;
            cc.same(&format!("x += y [{t}]"), e, a);
            let mut e = x;
            e.sub_assign(y);
            let mut a = px;
            a -= py;
            cc.same(&format!("x -= y [{t}]"), e, a);
            let mut e = x;
            e.mul_assign(y);
            let mut a = px;
            a *= py;
            cc.same(&format!("x *= y [{t}]"), e, a);
            if y.dbl() != 0.0 && y.dbl().abs() > 1e-290 {
                cc.same(&format!("x / y [{t}]"), x.div(y), px / py);
                let mut e = x;
                e.div_assign(y);
                let mut a = px;
                a /= py;
                cc.same(&format!("x /= y [{t}]"), e, a);
            }
            cc.same_bool(&format!("x < y [{t}]"), x.dbl() < y.dbl(), px < py);
            cc.same_bool(&format!("x > y [{t}]"), x.dbl() > y.dbl(), px > py);
            cc.same_bool(&format!("x == y [{t}]"), x.dbl() == y.dbl(), px == py);
        }
    }
}

#[test]
fn highs_cdouble_matches_header_bit_for_bit() {
    let mut cc = Check {
        n: 0,
        fails: Vec::new(),
    };
    port_ops_match_reference(&mut cc);
    cc.finish("HighsCDouble (HighsCDouble.h)");
}

/// The error-free transforms, which every compensated sum relies on.
#[test]
fn error_free_transforms_match_header() {
    let mut cc = Check {
        n: 0,
        fails: Vec::new(),
    };
    let ds = doubles();
    for &a in &ds {
        for &b in &ds {
            let (mut ex, mut ey, mut ax, mut ay) = (0.0, 0.0, 0.0, 0.0);
            Ref::two_sum(&mut ex, &mut ey, a, b);
            HighsCDouble::two_sum(&mut ax, &mut ay, a, b);
            let t = format!("a={}, b={}", bits(a), bits(b));
            cc.same_f64(&format!("two_sum.x [{t}]"), ex, ax);
            cc.same_f64(&format!("two_sum.y [{t}]"), ey, ay);
            if a.abs() < 1e290 && b.abs() < 1e290 {
                Ref::two_product(&mut ex, &mut ey, a, b);
                HighsCDouble::two_product(&mut ax, &mut ay, a, b);
                cc.same_f64(&format!("two_product.x [{t}]"), ex, ax);
                cc.same_f64(&format!("two_product.y [{t}]"), ey, ay);
            }
        }
        if a.abs() < 1e290 {
            let (mut ex, mut ey, mut ax, mut ay) = (0.0, 0.0, 0.0, 0.0);
            Ref::split(&mut ex, &mut ey, a);
            HighsCDouble::split(&mut ax, &mut ay, a);
            cc.same_f64(&format!("split.x [a={}]", bits(a)), ex, ax);
            cc.same_f64(&format!("split.y [a={}]", bits(a)), ey, ay);
        }
    }
    cc.finish("HighsCDouble::two_sum/split/two_product (HighsCDouble.h)");
}
