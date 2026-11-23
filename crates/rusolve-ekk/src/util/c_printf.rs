//! C `printf` floating-point conversions (`%f`, `%e`, `%g` with precision and width), as
//! printed by the platform libc the HiGHS oracle is built against (Apple libc, arm64).
//!
//! Not a HiGHS source file: HiGHS logs through `snprintf`, and Rust's `{}` float form is
//! the shortest round-trip string, which is not what `printf` prints. Every log line that
//! the port reproduces formats its doubles through these helpers.
//!
//! Verified value-for-value against Apple libc `snprintf` with a C probe over special
//! values (signed zeros, ties, subnormals, extremes, inf/nan) and 20000 pseudo-random
//! doubles, for `%g %.Ng %.Nf %e %.Ne %W.Pg %-W.Pg`. Observed libc facts the helpers
//! encode:
//! - rounding is exact decimal round-half-even (`%.1f` of 0.25 is `0.2`, of 0.35 `0.3`);
//! - every NaN prints as `nan`, whatever its sign bit (no `-nan` on Apple libc);
//! - infinities print as `inf` / `-inf`; `-0.0` keeps its sign (`-0`);
//! - the `%e` exponent has a sign and at least two digits.

/// The libc spelling of a non-finite value, or `None` for a finite one.
fn non_finite(val: f64) -> Option<&'static str> {
    if val.is_nan() {
        Some("nan")
    } else if val.is_infinite() {
        Some(if val < 0.0 { "-inf" } else { "inf" })
    } else {
        None
    }
}

/// Remove trailing zeros after a decimal point, and the point itself if nothing is left
/// (the `%g` rule without the `#` flag).
fn strip_fraction_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// Rust's `{:.*e}` mantissa and exponent, the exponent rendered C-style (`e+05`).
fn exp_parts(val: f64, prec: usize) -> (String, i32) {
    let sci = format!("{:.*e}", prec, val);
    let (mantissa, exponent) = sci.split_once('e').expect("`{:e}` always has an exponent");
    (
        mantissa.to_string(),
        exponent.parse().expect("`{:e}` exponent is an integer"),
    )
}

fn c_exponent(x: i32) -> String {
    format!("e{}{:02}", if x < 0 { '-' } else { '+' }, x.unsigned_abs())
}

/// C `%.{prec}f`.
pub fn c_fmt_f(val: f64, prec: usize) -> String {
    match non_finite(val) {
        Some(s) => s.to_string(),
        None => format!("{:.*}", prec, val),
    }
}

/// C `%.{prec}e`.
pub fn c_fmt_e(val: f64, prec: usize) -> String {
    if let Some(s) = non_finite(val) {
        return s.to_string();
    }
    let (mantissa, x) = exp_parts(val, prec);
    format!("{}{}", mantissa, c_exponent(x))
}

/// C `%.{prec}g` (no flags): `prec` significant digits (0 counts as 1), the `%e` form
/// when the exponent X of the value rounded to `prec` digits is `< -4` or `>= prec`,
/// trailing fractional zeros removed (C11 7.21.6.1).
pub fn c_fmt_g(val: f64, prec: usize) -> String {
    if let Some(s) = non_finite(val) {
        return s.to_string();
    }
    let p = prec.max(1);
    let (mantissa, x) = exp_parts(val, p - 1);
    if x < -4 || x >= p as i32 {
        format!("{}{}", strip_fraction_zeros(&mantissa), c_exponent(x))
    } else {
        let decimals = (p as i32 - 1 - x) as usize;
        strip_fraction_zeros(&format!("{:.*}", decimals, val)).to_string()
    }
}

/// C `%{width}s` of an already converted field: right-justified with spaces.
pub fn c_pad_right(s: &str, width: usize) -> String {
    format!("{:>width$}", s, width = width)
}

/// C `%-{width}s` of an already converted field: left-justified with spaces.
pub fn c_pad_left(s: &str, width: usize) -> String {
    format!("{:<width$}", s, width = width)
}

/// C `%{width}.{prec}f`.
pub fn c_fmt_f_width(val: f64, prec: usize, width: usize) -> String {
    c_pad_right(&c_fmt_f(val, prec), width)
}

/// C `%{width}.{prec}g`.
pub fn c_fmt_g_width(val: f64, prec: usize, width: usize) -> String {
    c_pad_right(&c_fmt_g(val, prec), width)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values copied from the Apple libc probe (see the module doc).
    #[test]
    fn matches_apple_libc() {
        assert_eq!(c_fmt_g(0.0, 6), "0");
        assert_eq!(c_fmt_g(-0.0, 6), "-0");
        assert_eq!(c_fmt_g(f64::NAN, 6), "nan");
        assert_eq!(c_fmt_g(-f64::NAN, 17), "nan");
        assert_eq!(c_fmt_g(f64::NEG_INFINITY, 6), "-inf");
        assert_eq!(c_fmt_g(0.35, 6), "0.35");
        assert_eq!(c_fmt_f(0.35, 1), "0.3");
        assert_eq!(c_fmt_f(0.25, 1), "0.2");
        assert_eq!(c_fmt_f(2.5, 0), "2");
        assert_eq!(c_fmt_g(1e-5, 6), "1e-05");
        assert_eq!(c_fmt_g(1e-4, 6), "0.0001");
        assert_eq!(c_fmt_g(123456.0, 6), "123456");
        assert_eq!(c_fmt_g(1234567.0, 6), "1.23457e+06");
        assert_eq!(c_fmt_g(999999.5, 6), "1e+06");
        assert_eq!(c_fmt_g(67.699168357953639, 17), "67.699168357953639");
        assert_eq!(c_fmt_g(1e100, 17), "1e+100");
        assert_eq!(c_fmt_e(1e100, 2), "1.00e+100");
        assert_eq!(c_fmt_e(0.0, 6), "0.000000e+00");
        assert_eq!(c_fmt_g_width(0.5, 3, 8), "     0.5");
        assert_eq!(c_pad_left(&c_fmt_g(0.5, 4), 12), "0.5         ");
    }
}
