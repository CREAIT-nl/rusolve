//! HiGHS source: util/HighsCDouble.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Compensated ("double-double") arithmetic. Every C++ operator overload maps
//! to the matching `std::ops` trait; each trait body translates the C++
//! operator named in its doc line. Comparison operators in C++ compare
//! `double(*this)` with `double(other)`; `PartialEq`/`PartialOrd` reproduce
//! that (the C++ class has no total order).
//!
//! Translate `double(x)` / `static_cast<double>(x)` as `f64::from(x)`, and
//! `HighsCDouble(v)` (implicit or explicit) as `HighsCDouble::from(v)`.
//! Rust never fuses `a * b + c` into an FMA, which the error-free
//! transformations below rely on (see crate docs).

use core::cmp::Ordering;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// HiGHS: HighsCDouble.h::HighsCDouble
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsCDouble {
    pub hi: f64,
    pub lo: f64,
}

impl HighsCDouble {
    /// HiGHS: HighsCDouble.h::HighsCDouble::two_sum (lines 33-39); out-params `x`, `y`.
    pub fn two_sum(x: &mut f64, y: &mut f64, a: f64, b: f64) {
        *x = a + b;
        let z = *x - a;
        *y = (a - (*x - z)) + (b - z);
    }

    /// HiGHS: HighsCDouble.h::HighsCDouble::split (lines 41-49); out-params `x`, `y`.
    pub fn split(x: &mut f64, y: &mut f64, a: f64) {
        let factor: f64 = ((1i32 << 27) + 1) as f64;
        let c = factor * a;
        *x = c - (c - a);
        *y = a - *x;
    }

    /// HiGHS: HighsCDouble.h::HighsCDouble::two_product (lines 51-57); out-params `x`, `y`.
    pub fn two_product(x: &mut f64, y: &mut f64, a: f64, b: f64) {
        *x = a * b;
        let (mut a1, mut a2, mut b1, mut b2) = (0.0, 0.0, 0.0, 0.0);
        Self::split(&mut a1, &mut a2, a);
        Self::split(&mut b1, &mut b2, b);
        *y = a2 * b2 - (((*x - a1 * b1) - a2 * b1) - a1 * b2);
    }

    /// HiGHS: HighsCDouble.h::HighsCDouble::HighsCDouble(double hi_, double lo_)
    pub const fn new(hi: f64, lo: f64) -> Self {
        Self { hi, lo }
    }

    /// HiGHS: HighsCDouble.h::HighsCDouble::renormalize (line 269)
    pub fn renormalize(&mut self) {
        let (a, b) = (self.hi, self.lo);
        Self::two_sum(&mut self.hi, &mut self.lo, a, b);
    }

    /// HiGHS: HighsCDouble.h friend abs (line 271)
    pub fn abs(self) -> HighsCDouble {
        if self < 0.0 { -self } else { self }
    }

    /// HiGHS: HighsCDouble.h friend sqrt (lines 273-286)
    pub fn sqrt(self) -> HighsCDouble {
        let c = f64::from(self).sqrt();

        // guard against division by zero
        if c == 0.0 {
            return HighsCDouble::from(0.0);
        }

        // calculate precise square root by newton step
        let mut res = self / c;
        res += c;
        // multiplication by 0.5 is exact
        res.hi *= 0.5;
        res.lo *= 0.5;
        res
    }

    /// HiGHS: HighsCDouble.h friend floor (lines 288-300)
    #[allow(clippy::double_comparisons)] // kept parallel to the C++
    pub fn floor(self) -> HighsCDouble {
        if self.abs() < 1.0 {
            if self == 0.0 || self > 0.0 {
                return HighsCDouble::from(0.0);
            }
            return HighsCDouble::from(-1.0);
        }
        let floor_x = f64::from(self).floor();
        let mut res = HighsCDouble::default();

        Self::two_sum(
            &mut res.hi,
            &mut res.lo,
            floor_x,
            f64::from(self - floor_x).floor(),
        );
        res
    }

    /// HiGHS: HighsCDouble.h friend ceil (lines 302-314)
    #[allow(clippy::double_comparisons)] // kept parallel to the C++
    pub fn ceil(self) -> HighsCDouble {
        if self.abs() < 1.0 {
            if self == 0.0 || self < 0.0 {
                return HighsCDouble::from(0.0);
            }
            return HighsCDouble::from(1.0);
        }
        let ceil_x = f64::from(self).ceil();
        let mut res = HighsCDouble::default();

        Self::two_sum(
            &mut res.hi,
            &mut res.lo,
            ceil_x,
            f64::from(self - ceil_x).ceil(),
        );
        res
    }

    /// HiGHS: HighsCDouble.h friend round (line 316)
    pub fn round(self) -> HighsCDouble {
        (self + 0.5).floor()
    }

    /// HiGHS: HighsCDouble.h friend ldexp (lines 318-320)
    pub fn ldexp(self, exp: i32) -> HighsCDouble {
        HighsCDouble::new(ldexp_f64(self.hi, exp), ldexp_f64(self.lo, exp))
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::HighsCDouble(double val) (line 64)
impl From<f64> for HighsCDouble {
    fn from(val: f64) -> Self {
        Self { hi: val, lo: 0.0 }
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator double (line 66): `hi + lo`
impl From<HighsCDouble> for f64 {
    fn from(v: HighsCDouble) -> f64 {
        v.hi + v.lo
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator+=(double) (lines 68-73)
impl AddAssign<f64> for HighsCDouble {
    fn add_assign(&mut self, v: f64) {
        let mut c = 0.0;
        let hi = self.hi;
        HighsCDouble::two_sum(&mut self.hi, &mut c, v, hi);
        self.lo += c;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator+=(const HighsCDouble&) (lines 75-79)
impl AddAssign<HighsCDouble> for HighsCDouble {
    fn add_assign(&mut self, v: HighsCDouble) {
        *self += v.hi;
        self.lo += v.lo;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator-=(double) (lines 81-84)
impl SubAssign<f64> for HighsCDouble {
    fn sub_assign(&mut self, v: f64) {
        *self += -v;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator-=(const HighsCDouble&) (lines 86-90)
impl SubAssign<HighsCDouble> for HighsCDouble {
    fn sub_assign(&mut self, v: HighsCDouble) {
        *self -= v.hi;
        self.lo -= v.lo;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator*=(double) (lines 92-97)
impl MulAssign<f64> for HighsCDouble {
    fn mul_assign(&mut self, v: f64) {
        let c = self.lo * v;
        let hi = self.hi;
        HighsCDouble::two_product(&mut self.hi, &mut self.lo, hi, v);
        *self += c;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator*=(const HighsCDouble&) (lines 99-106)
impl MulAssign<HighsCDouble> for HighsCDouble {
    fn mul_assign(&mut self, v: HighsCDouble) {
        let c1 = self.hi * v.lo;
        let c2 = self.lo * v.hi;
        let hi = self.hi;
        HighsCDouble::two_product(&mut self.hi, &mut self.lo, hi, v.hi);
        *self += c1;
        *self += c2;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator/=(double) (lines 108-115)
impl DivAssign<f64> for HighsCDouble {
    fn div_assign(&mut self, v: f64) {
        let d = HighsCDouble::new(self.hi / v, self.lo / v);
        let mut c = d * v - *self;
        c.hi /= v;
        c.lo /= v;
        *self = d - c;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator/=(const HighsCDouble&) (lines 117-125)
impl DivAssign<HighsCDouble> for HighsCDouble {
    fn div_assign(&mut self, v: HighsCDouble) {
        let vdbl = v.hi + v.lo;
        let d = HighsCDouble::new(self.hi / vdbl, self.lo / vdbl);
        let mut c = d * v - *self;
        c.hi /= vdbl;
        c.lo /= vdbl;
        *self = d - c;
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator-() (line 127)
impl Neg for HighsCDouble {
    type Output = HighsCDouble;
    fn neg(self) -> HighsCDouble {
        HighsCDouble::new(-self.hi, -self.lo)
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator+(double) (lines 129-136)
impl Add<f64> for HighsCDouble {
    type Output = HighsCDouble;
    fn add(self, v: f64) -> HighsCDouble {
        let mut res = HighsCDouble::default();

        HighsCDouble::two_sum(&mut res.hi, &mut res.lo, self.hi, v);
        res.lo += self.lo;

        res
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator+(const HighsCDouble&) (lines 138-143)
impl Add<HighsCDouble> for HighsCDouble {
    type Output = HighsCDouble;
    fn add(self, v: HighsCDouble) -> HighsCDouble {
        let mut res = self + v.hi;
        res.lo += v.lo;

        res
    }
}

/// HiGHS: HighsCDouble.h friend operator+(double, const HighsCDouble&) (lines 145-147)
impl Add<HighsCDouble> for f64 {
    type Output = HighsCDouble;
    fn add(self, b: HighsCDouble) -> HighsCDouble {
        b + self
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator-(double) (lines 149-156)
impl Sub<f64> for HighsCDouble {
    type Output = HighsCDouble;
    fn sub(self, v: f64) -> HighsCDouble {
        let mut res = HighsCDouble::default();

        HighsCDouble::two_sum(&mut res.hi, &mut res.lo, self.hi, -v);
        res.lo += self.lo;

        res
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator-(const HighsCDouble&) (lines 158-163)
impl Sub<HighsCDouble> for HighsCDouble {
    type Output = HighsCDouble;
    fn sub(self, v: HighsCDouble) -> HighsCDouble {
        let mut res = self - v.hi;
        res.lo -= v.lo;

        res
    }
}

/// HiGHS: HighsCDouble.h friend operator-(double, const HighsCDouble&) (lines 165-167)
impl Sub<HighsCDouble> for f64 {
    type Output = HighsCDouble;
    fn sub(self, b: HighsCDouble) -> HighsCDouble {
        -b + self
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator*(double) (lines 169-176)
impl Mul<f64> for HighsCDouble {
    type Output = HighsCDouble;
    fn mul(self, v: f64) -> HighsCDouble {
        let mut res = HighsCDouble::default();

        HighsCDouble::two_product(&mut res.hi, &mut res.lo, self.hi, v);
        res += self.lo * v;

        res
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator*(const HighsCDouble&) (lines 178-183)
impl Mul<HighsCDouble> for HighsCDouble {
    type Output = HighsCDouble;
    fn mul(self, v: HighsCDouble) -> HighsCDouble {
        let mut res = self * v.hi;
        res += self.hi * v.lo;

        res
    }
}

/// HiGHS: HighsCDouble.h friend operator*(double, const HighsCDouble&) (lines 185-187)
impl Mul<HighsCDouble> for f64 {
    type Output = HighsCDouble;
    fn mul(self, b: HighsCDouble) -> HighsCDouble {
        b * self
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator/(double) (lines 189-195)
impl Div<f64> for HighsCDouble {
    type Output = HighsCDouble;
    fn div(self, v: f64) -> HighsCDouble {
        let mut res = self;

        res /= v;

        res
    }
}

/// HiGHS: HighsCDouble.h::HighsCDouble::operator/(const HighsCDouble&) (lines 197-203)
impl Div<HighsCDouble> for HighsCDouble {
    type Output = HighsCDouble;
    fn div(self, v: HighsCDouble) -> HighsCDouble {
        let mut res = self;

        res /= v;

        res
    }
}

/// HiGHS: HighsCDouble.h friend operator/(double, const HighsCDouble&) (lines 205-207)
impl Div<HighsCDouble> for f64 {
    type Output = HighsCDouble;
    fn div(self, b: HighsCDouble) -> HighsCDouble {
        HighsCDouble::from(self) / b
    }
}

/// HiGHS: HighsCDouble.h operator==/!= (HighsCDouble, HighsCDouble) (lines 249-262)
impl PartialEq for HighsCDouble {
    fn eq(&self, other: &HighsCDouble) -> bool {
        f64::from(*self) == f64::from(*other)
    }
}

/// HiGHS: HighsCDouble.h operator==/!= (HighsCDouble, double) (lines 253, 263)
impl PartialEq<f64> for HighsCDouble {
    fn eq(&self, other: &f64) -> bool {
        f64::from(*self) == *other
    }
}

/// HiGHS: HighsCDouble.h friend operator==/!= (double, HighsCDouble) (lines 255-267)
impl PartialEq<HighsCDouble> for f64 {
    fn eq(&self, other: &HighsCDouble) -> bool {
        *self == f64::from(*other)
    }
}

/// HiGHS: HighsCDouble.h operator</>/<=/>= (HighsCDouble, HighsCDouble) (lines 209-247)
impl PartialOrd for HighsCDouble {
    fn partial_cmp(&self, other: &HighsCDouble) -> Option<Ordering> {
        f64::from(*self).partial_cmp(&f64::from(*other))
    }
}

/// HiGHS: HighsCDouble.h operator</>/<=/>= (HighsCDouble, double) (lines 213, 223, 233, 243)
impl PartialOrd<f64> for HighsCDouble {
    fn partial_cmp(&self, other: &f64) -> Option<Ordering> {
        f64::from(*self).partial_cmp(other)
    }
}

/// HiGHS: HighsCDouble.h friend operator</>/<=/>= (double, HighsCDouble) (lines 215, 225, 235, 245)
impl PartialOrd<HighsCDouble> for f64 {
    fn partial_cmp(&self, other: &HighsCDouble) -> Option<Ordering> {
        self.partial_cmp(&f64::from(*other))
    }
}

/// C `ldexp` for doubles (`std::ldexp`): exact scaling by 2^exp, including
/// subnormal results (same algorithm as musl `scalbn`). Pinned against Apple libm on arm64
/// (the oracle platform) by the unit test `ldexp_matches_apple_libm`.
pub fn ldexp_f64(x: f64, exp: i32) -> f64 {
    let mut y = x;
    let mut n = exp;
    if n > 1023 {
        y *= f64::from_bits(0x7fe0_0000_0000_0000); // 2^1023
        n -= 1023;
        if n > 1023 {
            y *= f64::from_bits(0x7fe0_0000_0000_0000);
            n -= 1023;
            if n > 1023 {
                n = 1023;
            }
        }
    } else if n < -1022 {
        // 2^-1022 * 2^53 keeps the intermediate rounding single
        let k = f64::from_bits(0x0360_0000_0000_0000); // 2^(-1022+53)
        y *= k;
        n += 1022 - 53;
        if n < -1022 {
            y *= k;
            n += 1022 - 53;
            if n < -1022 {
                n = -1022;
            }
        }
    }
    y * f64::from_bits(((0x3ff + n) as u64) << 52)
}

#[cfg(test)]
mod ldexp_tests {
    use super::ldexp_f64;

    /// `(bits(x), e, bits(std::ldexp(x, e)))` from a clang -O2 probe on arm64 macOS.
    const CASES: &[(u64, i32, u64)] = &[
        (0x3ff0000000000000, 0, 0x3ff0000000000000),
        (0x3ff0000000000000, 5, 0x4040000000000000),
        (0x3ff0000000000000, -1074, 0x0000000000000001),
        (0x3ff0000000000000, -1075, 0x0000000000000000),
        (0x3ff0000000000000, -1060, 0x0000000000004000),
        (0x3ff0000000000000, -1022, 0x0010000000000000),
        (0x3ff0000000000000, 1023, 0x7fe0000000000000),
        (0x3ff0000000000000, 1024, 0x7ff0000000000000),
        (0x3ff0000000000000, -2000, 0x0000000000000000),
        (0x3ff0000000000000, 2000, 0x7ff0000000000000),
        (0x3ff0000000000000, -1073, 0x0000000000000002),
        (0x4008000000000000, 0, 0x4008000000000000),
        (0x4008000000000000, 5, 0x4058000000000000),
        (0x4008000000000000, -1074, 0x0000000000000003),
        (0x4008000000000000, -1075, 0x0000000000000002),
        (0x4008000000000000, -1060, 0x000000000000c000),
        (0x4008000000000000, -1022, 0x0028000000000000),
        (0x4008000000000000, 1023, 0x7ff0000000000000),
        (0x4008000000000000, 1024, 0x7ff0000000000000),
        (0x4008000000000000, -2000, 0x0000000000000000),
        (0x4008000000000000, 2000, 0x7ff0000000000000),
        (0x4008000000000000, -1073, 0x0000000000000006),
        (0x3fe6666666666666, 0, 0x3fe6666666666666),
        (0x3fe6666666666666, 5, 0x4036666666666666),
        (0x3fe6666666666666, -1074, 0x0000000000000001),
        (0x3fe6666666666666, -1075, 0x0000000000000000),
        (0x3fe6666666666666, -1060, 0x0000000000002ccd),
        (0x3fe6666666666666, -1022, 0x000b333333333333),
        (0x3fe6666666666666, 1023, 0x7fd6666666666666),
        (0x3fe6666666666666, 1024, 0x7fe6666666666666),
        (0x3fe6666666666666, -2000, 0x0000000000000000),
        (0x3fe6666666666666, 2000, 0x7ff0000000000000),
        (0x3fe6666666666666, -1073, 0x0000000000000001),
        (0xbff8000000000000, 0, 0xbff8000000000000),
        (0xbff8000000000000, 5, 0xc048000000000000),
        (0xbff8000000000000, -1074, 0x8000000000000002),
        (0xbff8000000000000, -1075, 0x8000000000000001),
        (0xbff8000000000000, -1060, 0x8000000000006000),
        (0xbff8000000000000, -1022, 0x8018000000000000),
        (0xbff8000000000000, 1023, 0xffe8000000000000),
        (0xbff8000000000000, 1024, 0xfff0000000000000),
        (0xbff8000000000000, -2000, 0x8000000000000000),
        (0xbff8000000000000, 2000, 0xfff0000000000000),
        (0xbff8000000000000, -1073, 0x8000000000000003),
        (0x3ff0000000000001, 0, 0x3ff0000000000001),
        (0x3ff0000000000001, 5, 0x4040000000000001),
        (0x3ff0000000000001, -1074, 0x0000000000000001),
        (0x3ff0000000000001, -1075, 0x0000000000000001),
        (0x3ff0000000000001, -1060, 0x0000000000004000),
        (0x3ff0000000000001, -1022, 0x0010000000000001),
        (0x3ff0000000000001, 1023, 0x7fe0000000000001),
        (0x3ff0000000000001, 1024, 0x7ff0000000000000),
        (0x3ff0000000000001, -2000, 0x0000000000000000),
        (0x3ff0000000000001, 2000, 0x7ff0000000000000),
        (0x3ff0000000000001, -1073, 0x0000000000000002),
        (0x0000000000000001, 0, 0x0000000000000001),
        (0x0000000000000001, 5, 0x0000000000000020),
        (0x0000000000000001, -1074, 0x0000000000000000),
        (0x0000000000000001, -1075, 0x0000000000000000),
        (0x0000000000000001, -1060, 0x0000000000000000),
        (0x0000000000000001, -1022, 0x0000000000000000),
        (0x0000000000000001, 1023, 0x3cc0000000000000),
        (0x0000000000000001, 1024, 0x3cd0000000000000),
        (0x0000000000000001, -2000, 0x0000000000000000),
        (0x0000000000000001, 2000, 0x79d0000000000000),
        (0x0000000000000001, -1073, 0x0000000000000000),
        (0x0010000000000000, 0, 0x0010000000000000),
        (0x0010000000000000, 5, 0x0060000000000000),
        (0x0010000000000000, -1074, 0x0000000000000000),
        (0x0010000000000000, -1075, 0x0000000000000000),
        (0x0010000000000000, -1060, 0x0000000000000000),
        (0x0010000000000000, -1022, 0x0000000000000000),
        (0x0010000000000000, 1023, 0x4000000000000000),
        (0x0010000000000000, 1024, 0x4010000000000000),
        (0x0010000000000000, -2000, 0x0000000000000000),
        (0x0010000000000000, 2000, 0x7d10000000000000),
        (0x0010000000000000, -1073, 0x0000000000000000),
        (0x7fefffffffffffff, 0, 0x7fefffffffffffff),
        (0x7fefffffffffffff, 5, 0x7ff0000000000000),
        (0x7fefffffffffffff, -1074, 0x3ccfffffffffffff),
        (0x7fefffffffffffff, -1075, 0x3cbfffffffffffff),
        (0x7fefffffffffffff, -1060, 0x3dafffffffffffff),
        (0x7fefffffffffffff, -1022, 0x400fffffffffffff),
        (0x7fefffffffffffff, 1023, 0x7ff0000000000000),
        (0x7fefffffffffffff, 1024, 0x7ff0000000000000),
        (0x7fefffffffffffff, -2000, 0x02efffffffffffff),
        (0x7fefffffffffffff, 2000, 0x7ff0000000000000),
        (0x7fefffffffffffff, -1073, 0x3cdfffffffffffff),
        (0x3fe8000000000000, 0, 0x3fe8000000000000),
        (0x3fe8000000000000, 5, 0x4038000000000000),
        (0x3fe8000000000000, -1074, 0x0000000000000001),
        (0x3fe8000000000000, -1075, 0x0000000000000000),
        (0x3fe8000000000000, -1060, 0x0000000000003000),
        (0x3fe8000000000000, -1022, 0x000c000000000000),
        (0x3fe8000000000000, 1023, 0x7fd8000000000000),
        (0x3fe8000000000000, 1024, 0x7fe8000000000000),
        (0x3fe8000000000000, -2000, 0x0000000000000000),
        (0x3fe8000000000000, 2000, 0x7ff0000000000000),
        (0x3fe8000000000000, -1073, 0x0000000000000002),
    ];

    #[test]
    fn ldexp_matches_apple_libm() {
        for &(x, e, r) in CASES {
            let got = ldexp_f64(f64::from_bits(x), e).to_bits();
            assert_eq!(got, r, "ldexp({:e}, {e})", f64::from_bits(x));
        }
    }

    /// The node-queue / search `treeweight` increments `std::ldexp(1.0, 1 - depth)` used to go
    /// through `HighsCDouble::ldexp` on 1.0 (and a bit-construction copy); both must agree with
    /// `ldexp_f64` for every depth, including the subnormal and underflow range.
    #[test]
    fn ldexp_one_matches_cdouble_ldexp() {
        use super::HighsCDouble;
        for e in -1100..=1100 {
            let via_cdouble = f64::from(HighsCDouble::from(1.0).ldexp(e));
            assert_eq!(
                ldexp_f64(1.0, e).to_bits(),
                via_cdouble.to_bits(),
                "e = {e}"
            );
        }
    }
}
