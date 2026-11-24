//! Progress events (`--log-format jsonl`) and the watchdog (`--watchdog`) for the port MIP path
//! (`RUSOLVE_MIP_HIGHS=1 rusolve mip`, `crate::mip_highs::run_mip_highs`).
//!
//! Both are built on the ported HiGHS callbacks (`rusolve_ekk::lp_data::highs_callback`,
//! HiGHS `lp_data/HighsCallback.{h,cpp}`), so the search code is untouched:
//!
//! * `kCallbackMipImprovingSolution` (fired from `HighsMipSolverData::saveReportMipSolution`,
//!   HighsMipSolverData.cpp:2666) -> an `incumbent` event;
//! * `kCallbackMipInterrupt` (fired from `HighsMipSolverData::checkLimits`,
//!   HighsMipSolverData.cpp:2534, on every limit check) -> `bound` and `phase` events, and the
//!   interrupt channel: `HighsCallbackInput::user_interrupt` makes `checkLimits` return true
//!   with `kInterrupt` (HighsMipSolverData.cpp:2560-2575).
//!
//! No callback is registered unless one of the two flags is given, so a default run is the
//! run it always was (`tools/bench/identity.sh`, 20 of 20). The callbacks only read the
//! `HighsCallbackOutput` and never set `user_has_solution`.
//!
//! ## Events
//!
//! NDJSON on stderr, one object per line, schema `docs/schema/rusolve-event-1.schema.json`
//! (`schema_version` `rusolve.event/1`). Every event has `event`, `seq` and `t` (seconds since
//! the solve started). Events: `phase` (`presolve`, `root`, `tree`), `incumbent`, `bound`,
//! `watchdog`, `final`. HiGHS has no phase callback, so the phases are read off the callback
//! sequence: `presolve` is announced before `Highs::run`, `root` at the first MIP interrupt
//! callback (the solver has finished presolve and is evaluating the root node) and `tree`
//! when a callback first reports a positive node count. The HiGHS text log is not written in
//! this mode (stderr carries events only).
//!
//! The `final` event is built from the same `Envelope` as `--format json`, so its status,
//! objective, bound and gap are the envelope's.
//!
//! ## Watchdog
//!
//! With a finite time limit, a watcher thread sets an atomic flag at the limit (measured from
//! the solve start); the interrupt callback copies it into `user_interrupt`, and HiGHS ends the
//! solve at its next limit check. HiGHS's own `time_limit` is the same limit, so that
//! interrupt is a second, independent way out. HiGHS 1.16.0 has loops that never reach a
//! limit check (a `strengthenInequalities` cover loop that never ends on fiball
//! runs in a sub-MIP presolve, where callbacks are never active, `HighsMipSolverData.h`
//! `submip`), so the watcher also has a hard deadline: when the solve has not returned at
//! limit plus grace, it prints the result from the last state the callbacks saw
//! (`limit.enforced_by = "watchdog_kill"`, no `solution`) and ends the process.

use crate::envelope::{Envelope, Improving, Limit, limit_kind};
use rusolve_ekk::mip::mip_solver_data::*;
use rusolve_ekk::prelude::*;
use serde_json::{Map, Value, json};
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// `schema_version` of every event line.
pub const EVENT_SCHEMA: &str = "rusolve.event/1";

/// Minimum seconds between two `bound` events (the bound can move at every node).
const BOUND_EVENT_INTERVAL_S: f64 = 0.05;

/// The text of a `kSolutionSource*` code and its one-letter code in the HiGHS log,
/// HighsMipSolverData.cpp::solutionSourceToString (lines 87-162, ported at
/// `rusolve_ekk::mip::mip_solver_data::part1`).
#[allow(non_upper_case_globals)]
fn solution_source(code: HighsInt) -> Option<(&'static str, &'static str)> {
    Some(match code {
        kSolutionSourceBranching => ("B", "Branching"),
        kSolutionSourceCentralRounding => ("C", "Central rounding"),
        kSolutionSourceFeasibilityPump => ("F", "Feasibility pump"),
        kSolutionSourceHeuristic => ("H", "Heuristic"),
        kSolutionSourceShifting => ("I", "Shifting"),
        kSolutionSourceFeasibilityJump => ("J", "Feasibility jump"),
        kSolutionSourceSubMip => ("L", "Sub-MIP"),
        kSolutionSourceEmptyMip => ("P", "Empty MIP"),
        kSolutionSourceRandomizedRounding => ("R", "Randomized rounding"),
        kSolutionSourceSolveLp => ("S", "Solve LP"),
        kSolutionSourceEvaluateNode => ("T", "Evaluate node"),
        kSolutionSourceUnbounded => ("U", "Unbounded"),
        kSolutionSourceUserSolution => ("X", "User solution"),
        kSolutionSourceHighsSolution => ("Y", "HiGHS solution"),
        kSolutionSourceZiRound => ("Z", "ZI Round"),
        kSolutionSourceTrivialZ => ("z", "Trivial zero"),
        kSolutionSourceTrivialL => ("l", "Trivial lower"),
        kSolutionSourceTrivialU => ("u", "Trivial upper"),
        kSolutionSourceTrivialP => ("p", "Trivial point"),
        kSolutionSourceFPR => ("A", "FPR"),
        kSolutionSourceLocalMIP => ("M", "Local MIP"),
        kSolutionSourceShiftProp => ("Q", "Shift and propagate"),
        kSolutionSourceNodeRounding => ("D", "Node rounding"),
        _ => return None,
    })
}

fn finite(v: f64) -> Option<f64> {
    v.is_finite().then_some(v)
}

/// What the callbacks last saw.
struct State {
    seq: u64,
    phase: &'static str,
    improvements: u64,
    last_improvement_s: Option<f64>,
    objective: Option<f64>,
    nodes: i64,
    lp_iterations: i64,
    dual_bound: f64,
    gap: f64,
    last_bound_dual: f64,
    last_bound_event_s: f64,
}

/// The state the watchdog reports on a hard stop.
pub struct Snapshot {
    pub objective: Option<f64>,
    pub dual_bound: Option<f64>,
    pub gap: Option<f64>,
    pub nodes: i64,
    pub lp_iterations: i64,
}

/// Progress events and the watchdog flag, shared between the solver thread (callbacks) and
/// the watcher thread.
pub struct Tracker {
    start: Instant,
    jsonl: bool,
    state: Mutex<State>,
    interrupt: AtomicBool,
    interrupt_requested: AtomicBool,
    killed: AtomicBool,
    final_claimed: AtomicBool,
}

impl Tracker {
    pub fn new(start: Instant, jsonl: bool) -> Arc<Self> {
        Arc::new(Self {
            start,
            jsonl,
            state: Mutex::new(State {
                seq: 0,
                phase: "",
                improvements: 0,
                last_improvement_s: None,
                objective: None,
                nodes: 0,
                lp_iterations: 0,
                dual_bound: f64::NEG_INFINITY,
                gap: f64::INFINITY,
                last_bound_dual: f64::NAN,
                last_bound_event_s: f64::NEG_INFINITY,
            }),
            interrupt: AtomicBool::new(false),
            interrupt_requested: AtomicBool::new(false),
            killed: AtomicBool::new(false),
            final_claimed: AtomicBool::new(false),
        })
    }

    fn now(&self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }

    /// Writes one event line (a no-op without `--log-format jsonl`). One `write` under the
    /// stderr lock, so lines from the two threads never interleave.
    fn emit(&self, st: &mut State, event: &str, t: f64, fields: Value) {
        if !self.jsonl {
            return;
        }
        st.seq += 1;
        let mut map = Map::new();
        map.insert("schema_version".into(), json!(EVENT_SCHEMA));
        map.insert("event".into(), json!(event));
        map.insert("seq".into(), json!(st.seq));
        map.insert("t".into(), json!(t));
        if let Value::Object(extra) = fields {
            map.extend(extra);
        }
        let mut line = Value::Object(map).to_string();
        line.push('\n');
        let _ = std::io::stderr().lock().write_all(line.as_bytes());
    }

    /// `phase` event, once per phase.
    pub fn phase(&self, name: &'static str) {
        let mut st = self.state.lock().unwrap();
        self.set_phase(&mut st, name, self.now());
    }

    fn set_phase(&self, st: &mut State, name: &'static str, t: f64) {
        if st.phase != name {
            st.phase = name;
            self.emit(st, "phase", t, json!({ "phase": name }));
        }
    }

    fn on_progress(&self, out: &HighsCallbackOutput) {
        let t = self.now();
        let mut st = self.state.lock().unwrap();
        st.nodes = out.mip_node_count;
        st.lp_iterations = out.mip_total_lp_iterations;
        st.dual_bound = out.mip_dual_bound;
        st.gap = out.mip_gap;
        if st.phase == "presolve" {
            self.set_phase(&mut st, "root", t);
        }
        if st.phase == "root" && out.mip_node_count > 0 {
            self.set_phase(&mut st, "tree", t);
        }
        let moved = out.mip_dual_bound.is_finite()
            && out.mip_dual_bound.to_bits() != st.last_bound_dual.to_bits();
        if moved && t - st.last_bound_event_s >= BOUND_EVENT_INTERVAL_S {
            st.last_bound_dual = out.mip_dual_bound;
            st.last_bound_event_s = t;
            let fields = json!({
                "dual_bound": finite(out.mip_dual_bound),
                "primal_bound": finite(out.mip_primal_bound),
                "gap": finite(out.mip_gap),
                "nodes": out.mip_node_count,
                "lp_iterations": out.mip_total_lp_iterations,
            });
            self.emit(&mut st, "bound", t, fields);
        }
    }

    fn on_incumbent(&self, out: &HighsCallbackOutput) {
        let t = self.now();
        let mut st = self.state.lock().unwrap();
        st.improvements += 1;
        st.last_improvement_s = Some(t);
        st.objective = finite(out.objective_function_value);
        st.nodes = out.mip_node_count;
        st.lp_iterations = out.mip_total_lp_iterations;
        st.dual_bound = out.mip_dual_bound;
        st.gap = out.mip_gap;
        let source = solution_source(out.mip_solution_source);
        let fields = json!({
            "objective": st.objective,
            "source": source.map(|s| s.1),
            "source_code": source.map(|s| s.0),
            "dual_bound": finite(out.mip_dual_bound),
            "gap": finite(out.mip_gap),
            "nodes": out.mip_node_count,
            "lp_iterations": out.mip_total_lp_iterations,
            "improvement": st.improvements,
        });
        self.emit(&mut st, "incumbent", t, fields);
    }

    /// The callback to install with `Highs::set_callback`; start
    /// `kCallbackMipImprovingSolution` and `kCallbackMipInterrupt` after it.
    pub fn callback(self: &Arc<Self>) -> HighsCallbackFunctionType {
        let tracker = Arc::clone(self);
        HighsCallbackFunctionType::new(
            move |kind: HighsInt,
                  _message: &str,
                  out: &HighsCallbackOutput,
                  cb_in: &mut HighsCallbackInput| {
                if kind == kCallbackMipImprovingSolution {
                    tracker.on_incumbent(out);
                } else if kind == kCallbackMipInterrupt {
                    tracker.on_progress(out);
                    if tracker.interrupt.load(Ordering::Relaxed) {
                        cb_in.user_interrupt = true;
                    }
                }
            },
        )
    }

    /// Starts the watcher thread. At `limit` seconds after the solve start it requests the
    /// interrupt; at `limit + grace` it calls `on_kill` (which prints the result), then ends
    /// the process with exit code 0 unless the solve has already finished.
    pub fn spawn_watchdog(
        self: &Arc<Self>,
        limit: f64,
        grace: f64,
        on_kill: Box<dyn FnOnce(&Tracker) + Send>,
    ) {
        let tracker = Arc::clone(self);
        std::thread::spawn(move || {
            let mut on_kill = Some(on_kill);
            loop {
                if tracker.final_claimed.load(Ordering::SeqCst) {
                    return;
                }
                let t = tracker.now();
                if t >= limit && !tracker.interrupt_requested.swap(true, Ordering::SeqCst) {
                    tracker.interrupt.store(true, Ordering::Relaxed);
                    let mut st = tracker.state.lock().unwrap();
                    tracker.emit(
                        &mut st,
                        "watchdog",
                        t,
                        json!({ "action": "interrupt_requested", "limit_s": limit }),
                    );
                }
                if t >= limit + grace {
                    if tracker.final_claimed.swap(true, Ordering::SeqCst) {
                        return;
                    }
                    tracker.killed.store(true, Ordering::SeqCst);
                    {
                        let mut st = tracker.state.lock().unwrap();
                        tracker.emit(
                            &mut st,
                            "watchdog",
                            t,
                            json!({ "action": "hard_stop", "grace_s": grace }),
                        );
                    }
                    if let Some(f) = on_kill.take() {
                        f(&tracker);
                    }
                    let _ = std::io::stdout().flush();
                    std::process::exit(0);
                }
                std::thread::sleep(Duration::from_millis(5));
            }
        });
    }

    /// The main thread calls this when `Highs::run` returned. `false` means the watchdog is
    /// already printing the result and ending the process: the caller must not print.
    pub fn claim_final(&self) -> bool {
        !self.final_claimed.swap(true, Ordering::SeqCst)
    }

    /// What the callbacks last saw.
    pub fn snapshot(&self) -> Snapshot {
        let st = self.state.lock().unwrap();
        Snapshot {
            objective: st.objective,
            dual_bound: finite(st.dual_bound),
            gap: finite(st.gap),
            nodes: st.nodes,
            lp_iterations: st.lp_iterations,
        }
    }

    fn improving(&self, time_s: f64, has_solution: bool) -> Improving {
        let st = self.state.lock().unwrap();
        let since = st.last_improvement_s.map(|l| (time_s - l).max(0.0));
        Improving {
            tracked: true,
            count: Some(st.improvements),
            last_improvement_s: st.last_improvement_s,
            since_last_improvement_s: since,
            still_improving: Some(has_solution && since.is_some_and(|s| s <= 0.1 * time_s)),
        }
    }

    /// Fills the tracker-only parts of a finished solve's envelope: `improving`, and a
    /// `kInterrupt` the watchdog's own request caused is reported as the time limit it
    /// enforced (`Interrupted` is kept for an interrupt nobody here asked for).
    pub fn annotate(&self, env: &mut Envelope) {
        if env.status == "Interrupted" && self.interrupt_requested.load(Ordering::SeqCst) {
            env.status = "TimeLimit";
            env.limit = limit_kind(env.status).map(|kind| Limit {
                kind,
                enforced_by: "watchdog_interrupt",
                improving: self.improving(env.stats.time_s, env.has_solution),
            });
        }
        if let Some(limit) = &mut env.limit {
            limit.improving = self.improving(env.stats.time_s, env.has_solution);
        }
    }

    /// The envelope of a hard stop: `template` is the envelope of the model with no solve
    /// result, completed from the last state the callbacks saw.
    pub fn hard_stop_envelope(&self, mut template: Envelope) -> Envelope {
        let snap = self.snapshot();
        let time_s = self.now();
        template.status = "TimeLimit";
        template.objective = snap.objective;
        template.has_solution = snap.objective.is_some();
        template.dual_bound = snap.dual_bound;
        template.gap = if template.has_solution {
            snap.gap
        } else {
            None
        };
        template.stats.nodes = snap.nodes;
        template.stats.lp_iterations = snap.lp_iterations;
        template.stats.time_s = time_s;
        template.solution = None;
        template.limit = Some(Limit {
            kind: "time",
            enforced_by: "watchdog_kill",
            improving: self.improving(time_s, template.has_solution),
        });
        template
    }

    /// The `final` event, built from the envelope of the same run.
    pub fn emit_final(&self, env: &Envelope) {
        let mut st = self.state.lock().unwrap();
        let fields = json!({
            "status": env.status,
            "objective": env.objective,
            "dual_bound": env.dual_bound,
            "gap": env.gap,
            "has_solution": env.has_solution,
            "nodes": env.stats.nodes,
            "lp_iterations": env.stats.lp_iterations,
            "time_s": env.stats.time_s,
            "limit": env.limit,
        });
        let t = self.now();
        self.emit(&mut st, "final", t, fields);
    }
}
