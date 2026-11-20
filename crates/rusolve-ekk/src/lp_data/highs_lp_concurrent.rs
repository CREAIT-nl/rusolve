//! rusolve lever `lp_concurrent`, not a HiGHS file: concurrent LP
//! optimizer. Two fresh `Highs` instances on two threads solve a copy of the LP with their
//! own presolve and clean-up, one with dual simplex and one with the HiPO interior point plus
//! crossover; the first to return Optimal sets a stop flag that the other observes through the
//! simplex / IPM interrupt callbacks, and its solution and basis become this instance's.
//!
//! Source: Bixby, "Solving real-world linear programs: a decade and more of
//! progress", Oper. Res. 50(1), 2002; the same scheme at the MIP root is
//! `mip/lp_relaxation/root_race.rs`.

use crate::prelude::*;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

impl Highs {
    /// Structural rule of `lp_solver_auto`, from the LP before presolve. HiPO wins on sparse,
    /// not tall LPs whose predicted Cholesky factor is cheap: A has 3 to 10 nonzeros per
    /// column (more is allowed with up to 6e6 nonzeros; fewer than 3 is
    /// network or multicommodity-flow structure, where simplex wins), rows < 3x
    /// columns (tall LPs are dualised by simplex, `lp_dualize_auto`), and predicted factor flops
    /// at most 1e11 (above that one IPM iteration costs seconds), and at most 3e6 nonzeros
    /// (ordering plus analysis took 60-200 s and 7+ GB on L1_sixm250obs / fhnw-binschedule1, at
    /// 4.3M / 8.6M nonzeros; every measured win has at most 1.8M).
    /// Thresholds fitted on our benchmark set; only detected structure is
    /// used.
    pub(crate) fn lp_auto_picks_ipm(&mut self) -> bool {
        let lp = &self.model_lp_;
        let (nr, nc, nz) = (
            lp.num_row_ as f64,
            lp.num_col_ as f64,
            lp.a_matrix_.num_nz() as f64,
        );
        // More than 10 nonzeros per column makes the normal equations dense; that stays with
        // simplex unless the predicted factor is cheap (neos-5052403-cygnet, 149 per column,
        // 8.5e10 predicted flops before presolve: HiPO 62 s, simplex unsolved at 300 s).
        // The dense case allows 6e6 nonzeros.
        let dense = nz > 10.0 * nc;
        let (mut nz_cap, flops_cap) = if dense {
            (6.0e6, 1.0e11)
        } else {
            (3.0e6, 1.0e11)
        };
        // With ipm_hipo_amd_first the analysis orders the normal equations first and skips AS,
        // so it stays cheap when A A^T is small (rail4284: 11.3M nonzeros, 4284 rows, analysis
        // 0.6 s instead of 95 s)
        if self.options_.ipm_hipo_amd_first && nr * (nr + 1.0) / 2.0 <= 1.5e7 {
            nz_cap = 3.0e7;
        }
        if nc < 1.0 || nz < 1.0 || nr >= 3.0 * nc || nz < 3.0 * nc || nz > nz_cap {
            return false;
        }
        // The ordering and symbolic analysis poll nothing and cost 20-200 s on the largest
        // models, so with a finite limit they run on a worker the
        // caller stops waiting for after a quarter of the limit; a late answer counts as
        // "too expensive" and the worker ends on its own.
        let opts = self.options_.clone();
        let budget = opts.time_limit * 0.25;
        let f = if budget.is_finite() {
            let (tx, rx) = std::sync::mpsc::channel();
            let (lp, mut timer) = (self.model_lp_.clone(), self.timer_.clone());
            let spawned = std::thread::Builder::new()
                .stack_size(256 << 20)
                .spawn(move || {
                    crate::parallel::highs_parallel::initialize_scheduler(1);
                    let r = crate::ipm::hipo_wrapper::hipo_analysed_flops(&opts, &mut timer, &lp);
                    let _ = tx.send(r);
                });
            match spawned {
                Ok(_) => rx
                    .recv_timeout(std::time::Duration::from_secs_f64(budget.max(0.0)))
                    .ok()
                    .flatten()
                    .map_or(f64::INFINITY, |(f, _)| f),
                Err(_) => f64::INFINITY,
            }
        } else {
            crate::ipm::hipo_wrapper::hipo_analysed_flops(&opts, &mut self.timer_, &self.model_lp_)
                .map_or(f64::INFINITY, |(f, _)| f)
        };
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "LP auto: rows {} cols {} nnz {} flops {:.3e}\n",
            nr,
            nc,
            nz,
            f
        );
        f <= flops_cap
    }

    /// Race simplex against IPM on the incumbent LP. Returns `None` (nothing changed) when the
    /// lever is off, fewer than two threads are allowed, the model is not a plain LP, or neither
    /// racer reached optimal.
    pub(crate) fn try_lp_concurrent(&mut self) -> Option<HighsStatus> {
        if !self.options_.lp_concurrent
            || self.options_.threads < 2
            || self.model_lp_.is_mip()
            || self.model_hessian_.dim_ != 0
            || self.model_lp_.has_infinite_cost_
            || self.basis_.valid
            || self.solution_.value_valid
            || self.options_.solver != kHighsChooseString
            // two model copies and two solves: keep the race away from LPs whose copies are large
            || self.model_lp_.a_matrix_.num_nz() > 30_000_000
        {
            return None;
        }
        let lp = Arc::new(self.model_lp_.clone());
        let mut base = self.options_.clone();
        base.lp_concurrent = false;
        let user_limit = self.options_.time_limit;
        let stop = Arc::new(AtomicBool::new(false));
        let start = std::time::Instant::now();
        let (tx, rx) = std::sync::mpsc::channel::<(i32, Option<Highs>)>();
        // The racers are not joined: the loser may sit in a phase that does not poll the
        // interrupt (HiPO analyse / factorisation took 68 s to notice on neos-5251015), and
        // waiting for it would hand back the time the race won. It ends at its next poll.
        for k in 0..2i32 {
            let (lp, base, stop, tx) =
                (Arc::clone(&lp), base.clone(), Arc::clone(&stop), tx.clone());
            std::thread::Builder::new()
                .stack_size(256 << 20)
                .spawn(move || {
                    crate::parallel::highs_parallel::initialize_scheduler(1);
                    let mut h = Highs::with_options(base);
                    h.options_.hipo_available = true;
                    h.set_option_value_bool("output_flag", false);
                    h.set_option_value_int("threads", 1);
                    h.set_option_value_double("time_limit", user_limit);
                    h.set_option_value_str(
                        "solver",
                        if k == 0 { kSimplexString } else { kIpmString },
                    );
                    h.pass_model((*lp).clone());
                    let stop_cb = Arc::clone(&stop);
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
                    let status = h.run();
                    let ok = status != HighsStatus::kError
                        && h.get_model_status() == HighsModelStatus::kOptimal;
                    let _ = tx.send((k, if ok { Some(h) } else { None }));
                })
                .expect("spawn LP racer");
        }
        drop(tx);
        let mut w: i32 = -1;
        let mut won: Option<Highs> = None;
        for (k, h) in rx.iter() {
            if let Some(h) = h {
                w = k;
                won = Some(h);
                stop.store(true, Ordering::Release);
                break;
            }
        }
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "Concurrent LP: winner {} after {:.2}s\n",
            match w {
                0 => "dual simplex",
                1 => "ipm",
                _ => "none",
            },
            start.elapsed().as_secs_f64()
        );
        if w < 0 {
            return None;
        }
        let sub = won?;
        self.clear_solver();
        self.solution_ = sub.solution_.clone();
        self.basis_ = sub.basis_.clone();
        self.info_ = sub.info_.clone();
        self.model_status_ = sub.model_status_;
        Some(self.return_from_highs(HighsStatus::kOk))
    }
}
