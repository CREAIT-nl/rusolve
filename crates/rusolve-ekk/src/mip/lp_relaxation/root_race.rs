//! rusolve, not a HiGHS file: concurrent root LP.
//!
//! The root LP of a large MIP is where it stalls, and the winner there is model-dependent:
//! dual simplex or interior point plus crossover (gfd-schedulen180f7d50m30k18: 69 s against
//! 10.6 s, measured on our benchmark set). Running both and keeping the first to finish is the
//! concurrent LP approach of Bixby, "Solving real-world linear programs: a decade and more of
//! progress", Oper. Res. 50(1), 2002.
//!
//! Two fresh `Highs` instances on two threads solve a copy of the root LP, one with dual
//! simplex, one with IPX (or HiPO per `mip_ipm_solver`) and crossover. The first to return an
//! optimal valid basis sets a stop flag that the other observes through the simplex / IPM
//! interrupt callbacks. The winner's basis is installed in the relaxation's own `Highs`, so the
//! normal path then re-solves warm from an optimal basis (zero or a handful of iterations) and
//! every later re-solve is warm as usual.

use crate::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicI32, Ordering};

impl HighsLpRelaxation {
    /// Race the root LP; on success the relaxation's `lpsolver` holds an optimal basis.
    /// Returns false (nothing changed) when neither solver produced one.
    pub fn race_root_lp(&mut self, ms: &HighsMipSolver, time_limit: f64) -> bool {
        let lp = self.lpsolver.get_lp().clone();
        let base_options = self.lpsolver.options_.clone();
        let ipm = if ms.options_mip_.mip_ipm_solver == kHipoString {
            kHipoString
        } else {
            kIpxString
        };
        let stop = std::sync::Arc::new(AtomicBool::new(false));
        let winner = AtomicI32::new(-1);
        let start = std::time::Instant::now();
        let results: Vec<(Option<HighsBasis>, f64)> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..2)
                .map(|k| {
                    let (lp, base_options, stop, winner) = (&lp, &base_options, &stop, &winner);
                    std::thread::Builder::new()
                        .stack_size(256 << 20)
                        .spawn_scoped(scope, move || {
                            crate::parallel::highs_parallel::initialize_scheduler(1);
                            let t0 = std::time::Instant::now();
                            let mut h = Highs::default();
                            h.options_ = base_options.clone();
                            h.set_option_value_bool("output_flag", false);
                            h.set_option_value_int("threads", 1);
                            h.set_option_value_double("time_limit", time_limit);
                            h.set_option_value_str("presolve", kHighsOffString);
                            h.set_option_value_str(
                                "solver",
                                if k == 0 { kSimplexString } else { ipm },
                            );
                            if k == 0 {
                                h.set_option_value_int("simplex_strategy", kSimplexStrategyDual);
                            }
                            h.pass_model(lp.clone());
                            let stop_cb = std::sync::Arc::clone(stop);
                            h.set_callback(HighsCallbackFunctionType::new(
                                move |_ty, _msg, _out, input| {
                                    if stop_cb.load(Ordering::Relaxed) {
                                        input.user_interrupt = true;
                                    }
                                },
                            ));
                            h.start_callback(if k == 0 {
                                kCallbackSimplexInterrupt
                            } else {
                                kCallbackIpmInterrupt
                            });
                            let status = h.optimize_lp();
                            let ok = status != HighsStatus::kError
                                && h.get_model_status() == HighsModelStatus::kOptimal
                                && h.get_basis().valid;
                            let basis = if ok {
                                if winner
                                    .compare_exchange(-1, k, Ordering::AcqRel, Ordering::Acquire)
                                    .is_ok()
                                {
                                    stop.store(true, Ordering::Release);
                                }
                                Some(h.get_basis().clone())
                            } else {
                                None
                            };
                            (basis, t0.elapsed().as_secs_f64())
                        })
                        .expect("spawn root LP racer")
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("root LP racer panicked"))
                .collect()
        });
        let w = winner.load(Ordering::Acquire);
        crate::highs_log_user!(
            &ms.options_mip_.log_options,
            HighsLogType::kInfo,
            "Concurrent root LP: winner {} after {:.2}s (simplex {}, {} {})\n",
            match w {
                0 => "dual simplex",
                1 => ipm,
                _ => "none",
            },
            start.elapsed().as_secs_f64(),
            if results[0].0.is_some() {
                "finished"
            } else {
                "stopped"
            },
            ipm,
            if results[1].0.is_some() {
                "finished"
            } else {
                "stopped"
            },
        );
        if w < 0 {
            return false;
        }
        let basis = results[w as usize].0.clone().expect("winner has a basis");
        self.lpsolver
            .set_basis(&basis, "HighsLpRelaxation::race_root_lp");
        true
    }
}
