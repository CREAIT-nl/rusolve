//! HiGHS source: util/HighsIntegers.h (lines 1-212)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//! All members are `static`: `HighsIntegers::gcd(a, b)` stays a path call.

use crate::prelude::*;

/// HiGHS: HighsIntegers.h::HighsIntegers (lines 21-210; static helpers only)
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsIntegers;

impl HighsIntegers {
    /// HiGHS: HighsIntegers.h::HighsIntegers::mod (lines 23-26)
    /// C++: `static int64_t mod(int64_t a, int64_t m)`
    pub fn mod_i64(a: i64, m: i64) -> i64 {
        let r = a % m;
        r + i64::from(r < 0) * m
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::mod (lines 28-30)
    /// C++: `static double mod(double a, double m)`
    pub fn mod_f64(a: f64, m: f64) -> f64 {
        // `(a < 0) * m`: the bool promotes to 0 or 1 before the multiply
        (a % m).trunc() + (if a < 0.0 { 1.0 } else { 0.0 }) * m
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::nearestInteger (lines 32-34)
    /// C++: `static int64_t nearestInteger(double x)`
    pub fn nearest_integer(x: f64) -> i64 {
        (x + 0.5f64.copysign(x)) as i64
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::isIntegral (lines 36-39)
    /// C++: `static bool isIntegral(double x, double eps)`
    pub fn is_integral(x: f64, eps: f64) -> bool {
        let y = (x - x.trunc()).abs();
        cpp_min(y, 1.0 - y) <= eps
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::modularInverse (lines 41-65)
    /// C++: `static int64_t modularInverse(int64_t a, int64_t m)`
    pub fn modular_inverse(a: i64, m: i64) -> i64 {
        let mut a = a;
        let mut m = m;
        let mut y: i64 = 0;
        let mut x: i64 = 1;

        if m == 1 {
            return 0;
        }

        a = Self::mod_i64(a, m);

        while a > 1 {
            // compute quotient q = a / m and remainder r = a % m
            let q = a / m;
            let mut r = a - q * m;

            // update (a,m) = (m,r)
            a = m;
            m = r;

            // update x and y of extended euclidean algorithm
            r = x - q * y;
            x = y;
            y = r;
        }

        x
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::gcd (lines 67-86)
    /// C++: `static int64_t gcd(int64_t a, int64_t b)`
    pub fn gcd(a: i64, b: i64) -> i64 {
        debug_assert!(a != i64::MIN);
        debug_assert!(b != i64::MIN);

        let mut a = a;
        let mut b = b;
        if a < 0 {
            a = -a;
        }
        if b < 0 {
            b = -b;
        }

        if a == 0 {
            return b;
        }
        if b == 0 {
            return a;
        }

        // C++ do { h = a % b; a = b; b = h; } while (b != 0);
        loop {
            let h = a % b;
            a = b;
            b = h;
            if b == 0 {
                break;
            }
        }

        a
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::denominator (lines 88-125)
    /// C++: `static int64_t denominator(double x, double eps, int64_t maxdenom)`
    pub fn denominator(x: f64, eps: f64, maxdenom: i64) -> i64 {
        // computes a rational approximation with given maximal denominator
        let mut x = x;
        let mut ai = x as i64;
        let mut m: [i64; 4] = [ai, 1, 1, 0];

        let mut xi = HighsCDouble::from(x);
        let mut fraction = xi - (ai as f64);

        while fraction > eps {
            xi = 1.0 / fraction;
            if f64::from(xi) > (1i64 << 53) as f64 {
                break;
            }

            ai = f64::from(xi) as i64;
            let mut t = m[2] * ai + m[3];
            if t > maxdenom {
                break;
            }

            m[3] = m[2];
            m[2] = t;

            t = m[0] * ai + m[1];
            m[1] = m[0];
            m[0] = t;

            fraction = xi - (ai as f64);
        }

        ai = (maxdenom - m[3]) / m[2];
        m[1] += m[0] * ai;
        m[3] += m[2] * ai;

        let x0 = m[0] as f64 / m[2] as f64;
        let x1 = m[1] as f64 / m[3] as f64;
        x = x.abs();
        let err0 = (x - x0).abs();
        let err1 = (x - x1).abs();

        if err0 < err1 {
            return m[2];
        }
        m[3]
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::integralScale (lines 127-204)
    /// C++: `static double integralScale(const double* vals, HighsInt numVals, double deltadown, double deltaup)`
    /// (`HighsCDouble` arithmetic inside: keep it.)
    pub fn integral_scale(vals: &[f64], num_vals: HighsInt, deltadown: f64, deltaup: f64) -> f64 {
        if num_vals == 0 {
            return 0.0;
        }
        let vals = &vals[..num_vals as usize];

        // std::minmax_element with comparator |a| < |b|: the first smallest
        // element and the last largest element (standard guarantee)
        let mut min_idx = 0usize;
        let mut max_idx = 0usize;
        for i in 1..vals.len() {
            if vals[i].abs() < vals[min_idx].abs() {
                min_idx = i;
            }
            if !(vals[i].abs() < vals[max_idx].abs()) {
                max_idx = i;
            }
        }
        let minval = vals[min_idx];
        let maxval = vals[max_idx];

        let mut expshift: i32 = 0;

        // to cover many small denominators at once use a denominator of 75 * 2^n
        // with n-3 being large enough so that the smallest value is not below 0.5
        // but ignore tiny values bew deltadown/deltaup.
        if minval < -deltadown || minval > deltaup {
            expshift = crate::util::highs_hash::frexp(minval).1;
        }
        expshift = cpp_max(-expshift, 0) + 3;

        // guard against making the largest value too big which may cause overflows
        // with intermediate gcd values
        let mut exp_max_val = crate::util::highs_hash::frexp(maxval).1;
        exp_max_val = cpp_min(exp_max_val, 32);
        if exp_max_val + expshift > 32 {
            expshift = 32 - exp_max_val;
        }

        // C++ `uint64_t{75} << expshift`: a shift count >= 64 is undefined in C++;
        // wrapping_shl matches the x86-64 hardware masking of the count
        let mut denom: u64 = 75u64.wrapping_shl(expshift as u32);
        let mut startdenom: i64 = denom as i64;

        // now check if the values are integral and if not compute a common
        // denominator for their remaining fraction
        let mut val = (startdenom as f64) * HighsCDouble::from(vals[0]);
        let mut downval = (val + deltaup).floor();
        let mut fraction = val - downval;

        if fraction > deltadown {
            // use a continued fraction algorithm to compute small missing
            // denominators for the remaining fraction
            denom =
                denom.wrapping_mul(Self::denominator(f64::from(fraction), deltaup, 1000) as u64);
            val = (denom as f64) * HighsCDouble::from(vals[0]);
            downval = (val + deltaup).floor();
            fraction = val - downval;

            // if this is not sufficient for reaching integrality, we stop here
            if fraction > deltadown {
                return 0.0;
            }
        }

        let mut currgcd: u64 = f64::from(downval).abs() as u64;

        for i in 1..num_vals as usize {
            val = (denom as f64) * HighsCDouble::from(vals[i]);
            downval = (val + deltaup).floor();
            fraction = val - downval;

            if fraction > deltadown {
                val = (startdenom as f64) * HighsCDouble::from(vals[i]);
                fraction = val - val.floor();
                denom = denom
                    .wrapping_mul(Self::denominator(f64::from(fraction), deltaup, 1000) as u64);
                val = (denom as f64) * HighsCDouble::from(vals[i]);
                downval = (val + deltaup).floor();
                fraction = val - downval;

                if fraction > deltadown {
                    return 0.0;
                }
            }

            if currgcd != 1 {
                currgcd = Self::gcd(currgcd as i64, f64::from(downval) as i64) as u64;

                // if the denominator is large, divide by the current gcd to prevent
                // unnecessary overflows
                if denom > u64::from(u32::MAX) {
                    denom /= currgcd;
                    if startdenom != 1 {
                        startdenom /= Self::gcd(currgcd as i64, startdenom);
                    }
                    currgcd = 1;
                }
            }
        }

        (denom as f64) / (currgcd as f64)
    }

    /// HiGHS: HighsIntegers.h::HighsIntegers::integralScale (lines 206-209)
    /// C++: `static double integralScale(const std::vector<double>& vals, double deltadown, double deltaup)`
    pub fn integral_scale_vec(vals: &[f64], deltadown: f64, deltaup: f64) -> f64 {
        Self::integral_scale(vals, vals.len() as HighsInt, deltadown, deltaup)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mod_and_rounding() {
        assert_eq!(HighsIntegers::mod_i64(-7, 3), 2);
        assert_eq!(HighsIntegers::mod_i64(7, 3), 1);
        assert_eq!(HighsIntegers::mod_f64(-7.5, 3.0), 2.0);
        assert_eq!(HighsIntegers::nearest_integer(2.5), 3);
        assert_eq!(HighsIntegers::nearest_integer(-2.5), -3);
        assert_eq!(HighsIntegers::nearest_integer(2.4), 2);
        assert!(HighsIntegers::is_integral(2.0000001, 1e-6));
        assert!(!HighsIntegers::is_integral(2.5, 1e-6));
    }

    #[test]
    fn gcd_and_inverse() {
        assert_eq!(HighsIntegers::gcd(12, -18), 6);
        assert_eq!(HighsIntegers::gcd(0, 5), 5);
        assert_eq!(HighsIntegers::gcd(7, 0), 7);
        let inv = HighsIntegers::modular_inverse(3, 7);
        assert_eq!((3 * inv).rem_euclid(7), 1);
        assert_eq!(HighsIntegers::modular_inverse(5, 1), 0);
    }

    #[test]
    fn rational_denominator() {
        assert_eq!(HighsIntegers::denominator(0.5, 1e-9, 1000), 2);
        assert_eq!(HighsIntegers::denominator(0.25, 1e-9, 1000), 4);
    }

    #[test]
    fn integral_scale_of_quarters() {
        // 0.5, 0.25, 1.0 scaled by 4 are integral
        assert_eq!(
            HighsIntegers::integral_scale_vec(&[0.5, 0.25, 1.0], 1e-9, 1e-9),
            4.0
        );
        assert_eq!(HighsIntegers::integral_scale(&[], 0, 1e-9, 1e-9), 0.0);
    }
}
