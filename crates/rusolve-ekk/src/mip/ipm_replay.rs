//! Oracle replay of the IPM (IPX) results that the MIP solver uses (ours, not a HiGHS file).
//!
//! **A cross-check only**: the MIP runs the ported IPX at
//! both sites and no longer uses the records. When [`kMipAcFixtureEnv`] is set, each site
//! calls [`ipm_replay_check`] after its IPM solve, which panics at the first bit that
//! differs from the oracle's IPM result: the first place to look when a MIP trace
//! diverges. Unset (the default), nothing is read.
//!
//! The MIP solver runs IPX in two places:
//!
//! * `HighsLpRelaxation::run`, the IPM fallback after `kIterationLimit`
//!   (HighsLpRelaxation.cpp:1358-1421). The oracle record is `t = "lp_ipm"`.
//! * `HighsMipSolverData::startAnalyticCenterComputation` (HighsMipSolverData.cpp:409-539). The oracle record is `t = "ac"`.
//!
//! The oracle (`tools/highs-oracle`) writes one JSON line per IPM call when
//! `HIGHS_ORACLE_AC_DUMP=<path>` is set. The format is in `tools/highs-oracle/FORMAT.md`,
//! section "MIP IPM records". [`ipm_replay_next`] reads that file from the path in the
//! environment variable [`kMipAcFixtureEnv`] (`RUSOLVE_MIP_AC_FIXTURE`) and returns the
//! records in order.
//!
//! # Matching
//!
//! A record's `call` is its 0-based index over all IPM records of the oracle process, in
//! emission order (sub-MIPs included). [`ipm_replay_next`] keeps the same counter per
//! thread. The record found must have the expected `t` and `submip_level`; anything else
//! means the caller diverged from the oracle before this call, and it panics.
//!
//! # Record use (as the earlier record-driven replay did)
//!
//! * `lp_ipm`: [`ipm_basis_from_record`] is the basis the C++ passes to
//!   `lpsolver.setBasis(ipm.getBasis(), ...)` (HighsLpRelaxation.cpp:1419).
//! * `ac`: `stored == 1` means the centre (`centre`, every column) was assigned, with
//!   [`model_status_from_i64`]`(model_status)`.
//!
//! Vectors of doubles are C `%a` hexfloat strings ([`parse_hexfloat`]); `f64::from_str` does
//! not parse them.

use crate::prelude::*;
use std::cell::RefCell;

/// Environment variable naming the `.ac` fixture to replay (FORMAT.md "MIP IPM records").
pub const kMipAcFixtureEnv: &str = "RUSOLVE_MIP_AC_FIXTURE";

/// `t` of the `HighsLpRelaxation::run` IPM fallback record.
pub const kMipIpmRecordLpIpm: &str = "lp_ipm";
/// `t` of the `startAnalyticCenterComputation` record.
pub const kMipIpmRecordAc: &str = "ac";

/// One line of a `.ac` fixture. Fields that the record kind does not carry are `None` or
/// empty.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MipIpmRecord {
    /// `t`: `"ac"` or `"lp_ipm"`.
    pub t: String,
    /// `call`: process-wide 0-based index over all records.
    pub call: i64,
    /// `submip` (0/1).
    pub submip: i64,
    /// `submip_level` (`HighsMipSolver::submip_level`, 0 at the root).
    pub submip_level: i64,
    /// `model_status` (`HighsModelStatus` as an integer; 7 = optimal).
    pub model_status: i64,
    /// `ac`: `stored` (1 when the centre was assigned).
    pub stored: Option<i64>,
    /// `ac`: `centre` (the IPM `col_value`).
    pub centre: Vec<f64>,
    /// `lp_ipm`: `ipm_solver` (`"ipx"` or `"hipo"`).
    pub ipm_solver: Option<String>,
    /// `lp_ipm`: `basis_valid` (0/1).
    pub basis_valid: Option<i64>,
    /// `lp_ipm`: `col_status` (`HighsBasisStatus` as integers, after crossover).
    pub col_status: Vec<i64>,
    /// `lp_ipm`: `row_status`.
    pub row_status: Vec<i64>,
    /// `lp_ipm`: `col_value` (unscaled IPM solution).
    pub col_value: Vec<f64>,
    /// `lp_ipm`: `row_value`.
    pub row_value: Vec<f64>,
}

impl MipIpmRecord {
    /// `model_status == HighsModelStatus::kOptimal` (7).
    pub fn is_optimal(&self) -> bool {
        self.model_status == 7
    }

    /// Parse one JSON line of a `.ac` fixture.
    pub fn parse(line: &str) -> Result<Self, String> {
        let fields = parse_flat_json_object(line)?;
        let mut rec = MipIpmRecord::default();
        for (key, val) in fields {
            match key.as_str() {
                "t" => rec.t = val.as_str(&key)?.to_string(),
                "call" => rec.call = val.as_int(&key)?,
                "submip" => rec.submip = val.as_int(&key)?,
                "submip_level" => rec.submip_level = val.as_int(&key)?,
                "model_status" => rec.model_status = val.as_int(&key)?,
                "stored" => rec.stored = Some(val.as_int(&key)?),
                "centre" => rec.centre = val.as_hexfloat_vec(&key)?,
                "ipm_solver" => rec.ipm_solver = Some(val.as_str(&key)?.to_string()),
                "basis_valid" => rec.basis_valid = Some(val.as_int(&key)?),
                "col_status" => rec.col_status = val.as_int_vec(&key)?,
                "row_status" => rec.row_status = val.as_int_vec(&key)?,
                "col_value" => rec.col_value = val.as_hexfloat_vec(&key)?,
                "row_value" => rec.row_value = val.as_hexfloat_vec(&key)?,
                // Unknown keys are ignored so that the oracle may add fields.
                _ => {}
            }
        }
        if rec.t.is_empty() {
            return Err("record has no \"t\"".to_string());
        }
        Ok(rec)
    }
}

/// All records of a `.ac` fixture, in file order.
pub fn parse_mip_ipm_fixture(text: &str) -> Result<Vec<MipIpmRecord>, String> {
    text.lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty())
        .map(|(i, l)| MipIpmRecord::parse(l).map_err(|e| format!("line {}: {}", i + 1, e)))
        .collect()
}

struct ReplayState {
    /// `None` until the first call; `Some(None)` when no fixture is configured.
    records: Option<Option<Vec<MipIpmRecord>>>,
    /// Next process-wide call index.
    next_call: i64,
}

// Kept per thread (it is not part of any MIP type, so `HighsMipSolver` stays `Send`). The
// state lives on the thread that owns the `HighsMipSolver` and runs its IPM calls; record
// checking at threads > 1 is not implemented.
thread_local! {
    static REPLAY: RefCell<ReplayState> = const {
        RefCell::new(ReplayState { records: None, next_call: 0 })
    };
}

/// Install `records` as this thread's fixture (instead of reading [`kMipAcFixtureEnv`]) and
/// reset the call counter. `None` means "no fixture". For tests and for a driver that solves
/// several models on one thread (call it before each solve).
pub fn ipm_replay_reset(records: Option<Vec<MipIpmRecord>>) {
    REPLAY.with(|r| {
        let mut r = r.borrow_mut();
        r.records = Some(records);
        r.next_call = 0;
    });
}

/// Count one IPM call of kind `t` (`"ac"` / `"lp_ipm"`) made at `submip_level` and return the
/// oracle record for it, or `None` when no fixture is configured (or the fixture has no
/// record with this call index: the oracle run stopped earlier).
///
/// The fixture is read on the first call of the thread from the path in
/// [`kMipAcFixtureEnv`] (unless [`ipm_replay_reset`] installed one). Panics when the
/// variable names a file that cannot be read or parsed, or when the record's `t` /
/// `submip_level` differ from the call (the port diverged from the oracle).
pub fn ipm_replay_next(t: &str, submip_level: HighsInt) -> Option<MipIpmRecord> {
    REPLAY.with(|r| {
        let mut r = r.borrow_mut();
        if r.records.is_none() {
            let loaded = std::env::var_os(kMipAcFixtureEnv).map(|path| {
                let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
                    panic!(
                        "{kMipAcFixtureEnv}={}: cannot read: {e}",
                        path.to_string_lossy()
                    )
                });
                parse_mip_ipm_fixture(&text).unwrap_or_else(|e| {
                    panic!("{kMipAcFixtureEnv}={}: {e}", path.to_string_lossy())
                })
            });
            r.records = Some(loaded);
        }
        let call = r.next_call;
        r.next_call += 1;
        let records = r.records.as_ref().and_then(|x| x.as_ref())?;
        let rec = records.iter().find(|x| x.call == call)?;
        assert!(
            rec.t == t && rec.submip_level == submip_level as i64,
            "IPM replay: call {call} is ({t}, submip_level {submip_level}) in the port but \
             ({}, submip_level {}) in the oracle: the port diverged before this IPM call",
            rec.t,
            rec.submip_level
        );
        Some(rec.clone())
    })
}

/// `HighsModelStatus` from its C++ integer value (oracle records carry the int).
pub fn model_status_from_i64(v: i64) -> HighsModelStatus {
    use HighsModelStatus::*;
    const ALL: [HighsModelStatus; 20] = [
        kNotset,
        kLoadError,
        kModelError,
        kPresolveError,
        kSolveError,
        kPostsolveError,
        kModelEmpty,
        kOptimal,
        kInfeasible,
        kUnboundedOrInfeasible,
        kUnbounded,
        kObjectiveBound,
        kObjectiveTarget,
        kTimeLimit,
        kIterationLimit,
        kUnknown,
        kSolutionLimit,
        kInterrupt,
        kMemoryLimit,
        kHighsInterrupt,
    ];
    match ALL.iter().find(|s| **s as i64 == v) {
        Some(s) => *s,
        None => panic!("model status {v} out of range"),
    }
}

/// The `HighsBasis` of an `lp_ipm` record (what `lpsolver.setBasis(ipm.getBasis(), ..)`
/// received in the oracle, HighsLpRelaxation.cpp:1419). `HighsBasisStatus`
/// discriminants: kLower=0, kBasic=1, kUpper=2, kZero=3, kNonbasic=4 (HConst.h).
pub fn ipm_basis_from_record(rec: &MipIpmRecord) -> HighsBasis {
    fn status(v: i64) -> HighsBasisStatus {
        match v {
            0 => HighsBasisStatus::kLower,
            1 => HighsBasisStatus::kBasic,
            2 => HighsBasisStatus::kUpper,
            3 => HighsBasisStatus::kZero,
            4 => HighsBasisStatus::kNonbasic,
            _ => panic!("ipm record: HighsBasisStatus {v} out of range"),
        }
    }
    let mut basis = HighsBasis::default();
    basis.valid = rec.basis_valid.unwrap_or(0) != 0;
    basis.col_status = rec.col_status.iter().map(|&v| status(v)).collect();
    basis.row_status = rec.row_status.iter().map(|&v| status(v)).collect();
    basis
}

/// Cross-check one real IPM result against the oracle record of the same call (a no-op
/// without [`kMipAcFixtureEnv`]). Both sites call it exactly once per IPM call, after the
/// solve, so the call counter matches the oracle's. Panics on the first difference.
///
/// * `ac`: `stored` is the `sol.size() == numCol()` check; when stored, `model_status` and
///   every bit of `col_value` must match.
/// * `lp_ipm`: `model_status`, `basis` (validity and statuses) and every bit of `col_value`
///   / `row_value` must match.
pub fn ipm_replay_check(
    t: &str,
    submip_level: HighsInt,
    stored: bool,
    model_status: HighsModelStatus,
    solution: &HighsSolution,
    basis: &HighsBasis,
) {
    let Some(rec) = ipm_replay_next(t, submip_level) else {
        return;
    };
    let call = rec.call;
    let same = |what: &str, want: &[f64], got: &[f64]| {
        assert_eq!(
            want.len(),
            got.len(),
            "IPM replay check: call {call} {t}: {what} length"
        );
        for (i, (w, g)) in want.iter().zip(got).enumerate() {
            assert_eq!(
                w.to_bits(),
                g.to_bits(),
                "IPM replay check: call {call} {t}: {what}[{i}] oracle {w:e} port {g:e}"
            );
        }
    };
    if t == kMipIpmRecordAc {
        assert_eq!(
            rec.stored,
            Some(stored as i64),
            "IPM replay check: call {call} ac: stored"
        );
        if !stored {
            return;
        }
        same("centre", &rec.centre, &solution.col_value);
    } else {
        let want = ipm_basis_from_record(&rec);
        assert_eq!(
            want.valid, basis.valid,
            "IPM replay check: call {call} lp_ipm: basis_valid"
        );
        assert!(
            want.col_status == basis.col_status && want.row_status == basis.row_status,
            "IPM replay check: call {call} lp_ipm: basis statuses differ"
        );
        same("col_value", &rec.col_value, &solution.col_value);
        same("row_value", &rec.row_value, &solution.row_value);
    }
    assert_eq!(
        model_status_from_i64(rec.model_status),
        model_status,
        "IPM replay check: call {call} {t}: model_status"
    );
}

// ------------------------------------------------------------------ hexfloat

/// Parse a C `printf("%a")` double: `[-]0x<h>[.<hhh>]p<[+-]d>`, `inf`, `-inf`, `nan`
/// (case-insensitive; a leading `+` is accepted). Exact for every finite double that `%a`
/// prints (normal and subnormal), and keeps the sign of `-0x0p+0`.
pub fn parse_hexfloat(s: &str) -> Option<f64> {
    let s = s.trim();
    let (neg, body) = match s.as_bytes().first()? {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    let lower = body.to_ascii_lowercase();
    let mag = match lower.as_str() {
        "inf" | "infinity" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let rest = lower.strip_prefix("0x")?;
            let (digits, exp) = match rest.split_once('p') {
                Some((d, e)) => (d, e.parse::<i64>().ok()?),
                None => (rest, 0),
            };
            let (int_part, frac_part) = match digits.split_once('.') {
                Some((i, f)) => (i, f),
                None => (digits, ""),
            };
            if int_part.is_empty() && frac_part.is_empty() {
                return None;
            }
            // Mantissa as an integer; drop digits past 64 bits into a sticky shift.
            let mut mant: u64 = 0;
            let mut shift: i64 = 0;
            for (k, ch) in int_part.chars().chain(frac_part.chars()).enumerate() {
                let d = ch.to_digit(16)? as u64;
                if mant >> 60 == 0 {
                    mant = (mant << 4) | d;
                    if k >= int_part.len() {
                        shift -= 4;
                    }
                } else if k < int_part.len() {
                    shift += 4;
                }
            }
            ldexp(mant as f64, exp + shift)
        }
    };
    Some(if neg { -mag } else { mag })
}

/// `x * 2^e`, rounding once (the steps of 2^±1000 are exact while the value stays normal).
fn ldexp(mut x: f64, mut e: i64) -> f64 {
    if x == 0.0 {
        return x;
    }
    while e > 1000 {
        x *= f64::from_bits(((1000 + 1023) as u64) << 52);
        e -= 1000;
        if x.is_infinite() {
            return x;
        }
    }
    while e < -1000 {
        x *= f64::from_bits(((-1000 + 1023) as u64) << 52);
        e += 1000;
        if x == 0.0 {
            return x;
        }
    }
    x * f64::from_bits(((e + 1023) as u64) << 52)
}

// ------------------------------------------------------------------ minimal JSON

/// A JSON value of the flat records: numbers, strings and arrays of them.
#[derive(Clone, Debug, PartialEq)]
enum JsonVal {
    Int(i64),
    Num(f64),
    Str(String),
    Arr(Vec<JsonVal>),
}

impl JsonVal {
    fn as_int(&self, key: &str) -> Result<i64, String> {
        match self {
            JsonVal::Int(v) => Ok(*v),
            _ => Err(format!("\"{key}\": expected an integer")),
        }
    }
    fn as_str(&self, key: &str) -> Result<&str, String> {
        match self {
            JsonVal::Str(v) => Ok(v),
            _ => Err(format!("\"{key}\": expected a string")),
        }
    }
    fn as_int_vec(&self, key: &str) -> Result<Vec<i64>, String> {
        match self {
            JsonVal::Arr(v) => v.iter().map(|x| x.as_int(key)).collect(),
            _ => Err(format!("\"{key}\": expected an array")),
        }
    }
    fn as_hexfloat_vec(&self, key: &str) -> Result<Vec<f64>, String> {
        match self {
            JsonVal::Arr(v) => {
                v.iter()
                    .map(|x| match x {
                        JsonVal::Str(s) => parse_hexfloat(s)
                            .ok_or_else(|| format!("\"{key}\": bad hexfloat {s:?}")),
                        JsonVal::Int(i) => Ok(*i as f64),
                        JsonVal::Num(f) => Ok(*f),
                        JsonVal::Arr(_) => Err(format!("\"{key}\": nested array")),
                    })
                    .collect()
            }
            _ => Err(format!("\"{key}\": expected an array")),
        }
    }
}

struct JsonParser<'a> {
    b: &'a [u8],
    i: usize,
}

impl JsonParser<'_> {
    fn ws(&mut self) {
        while self.i < self.b.len() && self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }
    fn eat(&mut self, ch: u8) -> Result<(), String> {
        self.ws();
        if self.b.get(self.i) == Some(&ch) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected '{}' at byte {}", ch as char, self.i))
        }
    }
    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.b.get(self.i).copied()
    }
    fn string(&mut self) -> Result<String, String> {
        self.eat(b'"')?;
        let mut out = String::new();
        loop {
            let ch = *self.b.get(self.i).ok_or("unterminated string")?;
            self.i += 1;
            match ch {
                b'"' => return Ok(out),
                b'\\' => {
                    let e = *self.b.get(self.i).ok_or("unterminated escape")?;
                    self.i += 1;
                    out.push(match e {
                        b'"' => '"',
                        b'\\' => '\\',
                        b'/' => '/',
                        b'n' => '\n',
                        b't' => '\t',
                        b'r' => '\r',
                        _ => return Err(format!("unsupported escape '\\{}'", e as char)),
                    });
                }
                _ => {
                    // The records are ASCII; keep any UTF-8 byte sequence intact.
                    let start = self.i - 1;
                    let mut end = self.i;
                    while end < self.b.len() && (self.b[end] & 0xC0) == 0x80 {
                        end += 1;
                    }
                    out.push_str(
                        std::str::from_utf8(&self.b[start..end]).map_err(|e| e.to_string())?,
                    );
                    self.i = end;
                }
            }
        }
    }
    fn number(&mut self) -> Result<JsonVal, String> {
        self.ws();
        let start = self.i;
        while self.i < self.b.len()
            && matches!(
                self.b[self.i],
                b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9'
            )
        {
            self.i += 1;
        }
        let s = std::str::from_utf8(&self.b[start..self.i]).map_err(|e| e.to_string())?;
        if let Ok(v) = s.parse::<i64>() {
            return Ok(JsonVal::Int(v));
        }
        s.parse::<f64>()
            .map(JsonVal::Num)
            .map_err(|_| format!("bad number {s:?} at byte {start}"))
    }
    fn value(&mut self) -> Result<JsonVal, String> {
        match self.peek() {
            Some(b'"') => self.string().map(JsonVal::Str),
            Some(b'[') => {
                self.i += 1;
                let mut v = Vec::new();
                if self.peek() == Some(b']') {
                    self.i += 1;
                    return Ok(JsonVal::Arr(v));
                }
                loop {
                    v.push(self.value()?);
                    match self.peek() {
                        Some(b',') => self.i += 1,
                        Some(b']') => {
                            self.i += 1;
                            return Ok(JsonVal::Arr(v));
                        }
                        _ => return Err(format!("expected ',' or ']' at byte {}", self.i)),
                    }
                }
            }
            Some(_) => self.number(),
            None => Err("unexpected end of input".to_string()),
        }
    }
}

/// Parse one flat JSON object (`{"k": v, ...}` with number / string / array values).
fn parse_flat_json_object(line: &str) -> Result<Vec<(String, JsonVal)>, String> {
    let mut p = JsonParser {
        b: line.as_bytes(),
        i: 0,
    };
    p.eat(b'{')?;
    let mut out = Vec::new();
    if p.peek() == Some(b'}') {
        return Ok(out);
    }
    loop {
        let key = p.string()?;
        p.eat(b':')?;
        out.push((key, p.value()?));
        match p.peek() {
            Some(b',') => p.i += 1,
            Some(b'}') => {
                p.i += 1;
                break;
            }
            _ => return Err(format!("expected ',' or '}}' at byte {}", p.i)),
        }
    }
    p.ws();
    if p.i != p.b.len() {
        return Err(format!("trailing data at byte {}", p.i));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hexfloat_exact() {
        assert_eq!(parse_hexfloat("0x1p+0"), Some(1.0));
        assert_eq!(parse_hexfloat("0x1.8p+1"), Some(3.0));
        assert_eq!(parse_hexfloat("-0x1.4p-2"), Some(-0.3125));
        assert_eq!(
            parse_hexfloat("0x1.833b177af0ffep-1").map(f64::to_bits),
            Some(0x3FE8_33B1_77AF_0FFE)
        );
        let z = parse_hexfloat("-0x0p+0").unwrap();
        assert_eq!(z, 0.0);
        assert!(z.is_sign_negative());
        assert!(parse_hexfloat("0x0p+0").unwrap().is_sign_positive());
        assert_eq!(parse_hexfloat("0x1.fffffffffffffp+1023"), Some(f64::MAX));
        assert_eq!(parse_hexfloat("0x1p-1022"), Some(f64::MIN_POSITIVE));
        // glibc prints subnormals as 0x0.<frac>p-1022, macOS as 0x1.<frac>p-10xx.
        assert_eq!(
            parse_hexfloat("0x0.0000000000001p-1022").map(f64::to_bits),
            Some(1)
        );
        assert_eq!(parse_hexfloat("0x1p-1074").map(f64::to_bits), Some(1));
        assert_eq!(
            parse_hexfloat("0x1.ffffffffffffep-1023").map(f64::to_bits),
            Some(0x000F_FFFF_FFFF_FFFF)
        );
        assert_eq!(parse_hexfloat("0X1P+3"), Some(8.0));
        assert_eq!(parse_hexfloat("inf"), Some(f64::INFINITY));
        assert_eq!(parse_hexfloat("-inf"), Some(f64::NEG_INFINITY));
        assert!(parse_hexfloat("nan").unwrap().is_nan());
        assert_eq!(parse_hexfloat("1.5"), None);
        assert_eq!(parse_hexfloat("0xp+0"), None);
        assert_eq!(parse_hexfloat(""), None);
        // Round trip over a spread of bit patterns (%a prints 13 fraction digits).
        for bits in [
            0x3FF0_0000_0000_0001u64,
            0x7FEF_0123_4567_89AB,
            0x0010_0000_0000_0000,
            0x4123_4567_89AB_CDEF,
        ] {
            let v = f64::from_bits(bits);
            let exp = ((bits >> 52) & 0x7FF) as i64 - 1023;
            let s = format!("0x1.{:013x}p{:+}", bits & 0x000F_FFFF_FFFF_FFFF, exp);
            assert_eq!(
                parse_hexfloat(&s).map(f64::to_bits),
                Some(v.to_bits()),
                "{s}"
            );
        }
    }

    #[test]
    fn record_parse() {
        let line = r#"{"t":"lp_ipm","call":1,"submip":0,"submip_level":0,"model_status":7,"ipm_solver":"ipx","basis_valid":1,"col_status":[2,0],"row_status":[1],"col_value":["0x1p+0","-0x0p+0"],"row_value":["0x1.8p+1"]}"#;
        let r = MipIpmRecord::parse(line).unwrap();
        assert_eq!(r.t, kMipIpmRecordLpIpm);
        assert_eq!(r.call, 1);
        assert!(r.is_optimal());
        assert_eq!(r.ipm_solver.as_deref(), Some("ipx"));
        assert_eq!(r.basis_valid, Some(1));
        assert_eq!(r.col_status, vec![2, 0]);
        assert_eq!(r.row_status, vec![1]);
        assert_eq!(r.col_value[0], 1.0);
        assert!(r.col_value[1].is_sign_negative());
        assert_eq!(r.row_value, vec![3.0]);
        assert!(MipIpmRecord::parse(r#"{"call":1}"#).is_err());
        assert!(MipIpmRecord::parse(r#"{"t":"ac","call":"x"}"#).is_err());
    }

    #[test]
    fn replay_matches_by_call_index() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/mip/enlight_hard.ac"
        );
        let recs = parse_mip_ipm_fixture(&std::fs::read_to_string(path).unwrap()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!((recs[0].t.as_str(), recs[0].call), (kMipIpmRecordAc, 0));
        assert_eq!(recs[0].stored, Some(1));
        assert_eq!(recs[0].centre.len(), 200);
        let lp = &recs[1];
        assert_eq!(
            (lp.t.as_str(), lp.call, lp.model_status),
            (kMipIpmRecordLpIpm, 1, 7)
        );
        assert_eq!((lp.col_status.len(), lp.row_status.len()), (200, 285));
        assert_eq!((lp.col_value.len(), lp.row_value.len()), (200, 285));

        ipm_replay_reset(Some(recs));
        let ac = ipm_replay_next(kMipIpmRecordAc, 0).unwrap();
        assert_eq!(ac.call, 0);
        let ipm = ipm_replay_next(kMipIpmRecordLpIpm, 0).unwrap();
        assert_eq!(ipm.call, 1);
        assert!(ipm_replay_next(kMipIpmRecordLpIpm, 0).is_none());
        ipm_replay_reset(None);
        assert!(ipm_replay_next(kMipIpmRecordAc, 0).is_none());
    }

    #[test]
    #[should_panic(expected = "diverged")]
    fn replay_kind_mismatch_panics() {
        let rec = MipIpmRecord {
            t: kMipIpmRecordAc.to_string(),
            ..Default::default()
        };
        ipm_replay_reset(Some(vec![rec]));
        ipm_replay_next(kMipIpmRecordLpIpm, 0);
    }
}
