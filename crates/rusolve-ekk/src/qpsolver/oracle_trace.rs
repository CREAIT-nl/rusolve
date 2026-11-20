//! Rust mirror of the oracle QUASS trace:
//! `tools/highs-oracle/highs-oracle-w1-hooks.patch`
//! (`highs/qpsolver/HighsOracleQpDump.h`). Not part of upstream HiGHS.
//!
//! The patched C++ appends one `qp_iter` JSON line at the end of every pass of
//! the `while (true)` loop of `Quass::solve` (quass.cpp 351-524: after
//! `redcosts.update()` on a step, after `atfsep = true` on a null step), and
//! one `qp_end` line after the post-loop `loginformation` (quass.cpp 526).
//! The Rust writes the identical records at the same points, so a
//! fixture comparison is a byte diff. Doubles are the 16 hex digits of their
//! bit pattern; vectors are FNV-1a 64 hashes of the little-endian bit
//! patterns of the dense `value` array (FORMAT.md "Hashes").
//!
//! Sinks: [`capture_begin`] / [`capture_end`] install an in-memory buffer on
//! the current thread (tests). Otherwise, when the environment variable
//! [`kQpTraceEnv`] (`RUSOLVE_QP_TRACE`) names a file, it is truncated by the
//! first record of the thread and appended to afterwards (the C++
//! `HIGHS_ORACLE_QP_DUMP` file is only appended to, since the oracle writes
//! the `w1_header` line into it before the run). Without either nothing is written and nothing is
//! evaluated. Same design as `ipm/ipx/oracle_trace.rs`.

use super::qpvector::QpVector;
use std::cell::RefCell;
use std::fmt::Write as _;
use std::io::Write as _;

/// Environment variable naming the trace file (C++ `HIGHS_ORACLE_QP_DUMP`).
#[allow(non_upper_case_globals)]
pub const kQpTraceEnv: &str = "RUSOLVE_QP_TRACE";

enum Sink {
    Unopened,
    Off,
    File(std::fs::File),
    Buffer(String),
}

struct Trace {
    sink: Sink,
    call: i64,
}

thread_local! {
    static TRACE: RefCell<Trace> = const { RefCell::new(Trace { sink: Sink::Unopened, call: 0 }) };
}

/// Starts capturing the trace of this thread in memory (call counter reset).
pub fn capture_begin() {
    TRACE.with(|t| {
        *t.borrow_mut() = Trace {
            sink: Sink::Buffer(String::new()),
            call: 0,
        }
    });
}

/// Ends an in-memory capture and returns the records written since
/// [`capture_begin`]; the thread falls back to the environment variable.
pub fn capture_end() -> String {
    TRACE.with(|t| {
        let old = std::mem::replace(
            &mut *t.borrow_mut(),
            Trace {
                sink: Sink::Unopened,
                call: 0,
            },
        );
        match old.sink {
            Sink::Buffer(s) => s,
            _ => String::new(),
        }
    })
}

/// True when a record would be written (opens the env file on first use).
pub fn enabled() -> bool {
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        if let Sink::Unopened = t.sink {
            t.sink = match std::env::var_os(kQpTraceEnv) {
                Some(p) if !p.is_empty() => match std::fs::File::create(&p) {
                    Ok(f) => Sink::File(f),
                    Err(_) => Sink::Off,
                },
                _ => Sink::Off,
            };
        }
        !matches!(t.sink, Sink::Off)
    })
}

/// Writes one record; `body` gets the call number.
fn emit(head: &str, body: impl FnOnce(&mut String)) {
    if !enabled() {
        return;
    }
    TRACE.with(|t| {
        let mut t = t.borrow_mut();
        let call = t.call;
        t.call += 1;
        let mut s = String::new();
        let _ = write!(s, "{{\"t\":\"{head}\",\"call\":{call}");
        body(&mut s);
        s.push_str("}\n");
        match &mut t.sink {
            Sink::Buffer(b) => b.push_str(&s),
            Sink::File(f) => {
                let _ = f.write_all(s.as_bytes());
                let _ = f.flush();
            }
            _ => {}
        }
    });
}

fn put_f(s: &mut String, key: &str, v: f64) {
    let _ = write!(s, ",\"{key}\":\"{:016x}\"", v.to_bits());
}

fn put_i(s: &mut String, key: &str, v: i64) {
    let _ = write!(s, ",\"{key}\":{v}");
}

/// FNV-1a 64 over the little-endian bit patterns of `v.value[0..v.dim]`.
fn put_h(s: &mut String, key: &str, v: &QpVector) {
    let mut h: u64 = 0xcbf29ce484222325;
    let n = (v.dim.max(0) as usize).min(v.value.len());
    for x in &v.value[..n] {
        let bits = x.to_bits();
        for b in 0..8 {
            h ^= (bits >> (8 * b)) & 0xff;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    let _ = write!(s, ",\"{key}\":\"{h:016x}\"");
}

/// One pass of the QUASS main loop (record `qp_iter`, FORMAT.md "W1 records").
///
/// `minidx` is the pricing result of a major pass (`-1` on a minor pass);
/// `limiting` is `stepres.limitingconstraint` (`-2` when no ratio test ran,
/// i.e. a null step); `alpha` is `stepres.alpha` (`0.0` on a null step).
/// `gradient` is the stored gradient vector as it stands (C++
/// `Gradient::oracleGradient()`, added by the oracle hook patch): the caller must not
/// go through `get_gradient`, whose recompute would perturb the run.
#[allow(clippy::too_many_arguments)]
pub fn qp_iter(
    num_iterations: i64,
    atfsep: bool,
    minidx: i64,
    limiting: i64,
    alpha: f64,
    num_active: i64,
    num_inactive: i64,
    primal: &QpVector,
    rowactivity: &QpVector,
    gradient: &QpVector,
    p: &QpVector,
) {
    emit("qp_iter", |s| {
        put_i(s, "iter", num_iterations);
        put_i(s, "atfsep", atfsep as i64);
        put_i(s, "minidx", minidx);
        put_i(s, "limiting", limiting);
        put_f(s, "alpha", alpha);
        put_i(s, "num_active", num_active);
        put_i(s, "num_inactive", num_inactive);
        put_h(s, "h_primal", primal);
        put_h(s, "h_rowactivity", rowactivity);
        put_h(s, "h_gradient", gradient);
        put_h(s, "h_p", p);
    });
}

/// End of `Quass::solve` main loop (record `qp_end`): `runtime.status` as
/// its `QpModelStatus` integer value.
pub fn qp_end(num_iterations: i64, status: i64, primal: &QpVector, rowactivity: &QpVector) {
    emit("qp_end", |s| {
        put_i(s, "iter", num_iterations);
        put_i(s, "status", status);
        put_h(s, "h_primal", primal);
        put_h(s, "h_rowactivity", rowactivity);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capture_round_trip() {
        capture_begin();
        let v = QpVector {
            num_nz: 1,
            dim: 2,
            index: vec![0, 0],
            value: vec![1.0, 0.0],
        };
        qp_end(3, 1, &v, &v);
        let out = capture_end();
        assert!(out.starts_with("{\"t\":\"qp_end\",\"call\":0,\"iter\":3,\"status\":1,"));
        assert!(out.ends_with("}\n"));
        assert!(capture_end().is_empty());
    }
}
