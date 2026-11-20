//! HiGHS source: ipm/hipo/ipm/FactorHiGHSSolver.cpp (lines 357-524)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Choice of the fill-reducing ordering (`chooseOrdering`; the choice of
//! the linear system, `chooseNla`, is in `parallel.rs`).
//! With `options.parallel == "off"` the orderings run in sequence, in
//! `orderings_to_try` order; otherwise in `highs::parallel::for_each` with
//! grain 1 (C++ lines 455-461), which at one worker visits the
//! indices in the same order. The `hipo_ord` trace records are buffered per
//! index and written after the loop, in index order, from the calling thread
//! (see `tools/highs-oracle/highs-oracle-w5-hooks.patch`).

use super::FactorHiGHSSolver;
use crate::HighsInt;
use crate::ipm::hipo::auxiliary::fullFromLower;
use crate::ipm::hipo::auxiliary::logger::Logger;
use crate::ipm::hipo::factorhighs::factor_highs::FHsolver;
use crate::ipm::hipo::factorhighs::symbolic::Symbolic;
use crate::ipm::hipo::ipm::options::Options;
use crate::ipm::hipo::ipm::parameters::{kFlopsOrderingThresh, kLargeStorageGB};
use crate::ipm::hipo::ipm::status::{kStatusErrorAnalyse, kStatusOk};
use crate::ipm::hipo::oracle_trace;
use crate::ipm::hipo::orderings::amd::{AMD_CONTROL, AMD_INFO, AMD_OK};
use crate::ipm::hipo::orderings::extras::{amd, metis, rcm};
use crate::ipm::hipo::orderings::metis::{
    METIS_NOPTIONS, METIS_OK, METIS_OPTION_DBGLVL, METIS_OPTION_NO2HOP, METIS_OPTION_SEED, idx_t,
};
use crate::ipm::hipo::{Int, Int64};
use crate::lp_data::hconst::{kHighsChooseString, kHighsInf, kHighsOffString};
use crate::lp_data::highs_options::{kHipoAmdString, kHipoMetisString, kHipoRcmString};
use crate::parallel::highs_parallel;
use crate::util::c_printf::{c_fmt_e, c_fmt_f};

/// Ours (no C++ counterpart): a raw pointer the `run_ordering_and_analyse`
/// tasks share; the C++ lambda captures the vectors and `this`
/// by reference. Every dereference names its protocol in a `SAFETY` comment.
struct SlotPtr<T>(*mut T);

// SAFETY: the pointee is reached only under the protocol named at each use
// (disjoint index i, or the read-only `FHsolver::analyse`).
unsafe impl<T> Send for SlotPtr<T> {}
// SAFETY: as above.
unsafe impl<T> Sync for SlotPtr<T> {}

impl<T> SlotPtr<T> {
    /// C++ `&v[i]` (the pointer offset by `i` elements).
    ///
    /// # Safety
    /// `i` is in bounds of the allocation the pointer was taken from.
    unsafe fn at(&self, i: usize) -> *mut T {
        // SAFETY: forwarded to the caller
        unsafe { self.0.add(i) }
    }
}

impl FactorHiGHSSolver {
    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::chooseOrdering (lines 357-524)
    /// `FH_` and `ordering` (`ordering_AS_` / `ordering_NE_`) are passed
    /// explicitly; `NE_nz_limit` is `kkt_.NE_nz_limit.load()`. Orderings are
    /// called through `crate::ipm::hipo::orderings` (amd/metis/rcm, the
    /// `HighsExtrasExternalDeps.h` API). Oracle trace hook: one `oracle_trace::hipo_ord`
    /// record per ordering, right after the "Finished %s for %s" log line
    /// (before the `failure[i]` test and `FH_.analyse`).
    pub fn chooseOrdering(
        FH_: &mut FHsolver,
        options: &Options,
        logger: &Logger,
        NE_nz_limit: Int64,
        rows: &[Int],
        ptr: &[Int],
        signs: &[Int],
        S: &mut Symbolic,
        ordering: &mut String,
        nla: &str,
    ) -> Int {
        // Run analyse phase.
        // - If ordering is "amd", "metis", "rcm" run only the ordering requested.
        // - If ordering is "choose", run "amd", "metis", and choose the best.

        // select which fill-reducing orderings should be tried
        let mut orderings_to_try: Vec<String> = Vec::new();
        if options.ordering != kHighsChooseString {
            orderings_to_try.push(options.ordering.clone());
        } else {
            orderings_to_try.push(kHipoAmdString.to_string());
            orderings_to_try.push(kHipoMetisString.to_string());

            // rcm is much worse in general, so no point in trying for now
        }

        // vector<bool> is not thread-safe
        // (std::vector<char>: stores are truncated to 8 bits as in C++)
        let mut failure: Vec<i8> = vec![0; orderings_to_try.len()];

        if nla == "NE" {
            if (ptr[ptr.len() - 1] as Int64) >= NE_nz_limit {
                logger.printInfo("NE interrupted\n");
                return kStatusErrorAnalyse;
            }
        }

        // compute full-format matrix without diagonal entries
        let mut full_ptr: Vec<Int> = Vec::new();
        let mut full_rows: Vec<Int> = Vec::new();
        fullFromLower(ptr, rows, &mut full_ptr, &mut full_rows);
        let n: Int = full_ptr.len() as Int - 1;

        let mut permutations: Vec<Vec<Int>> = vec![vec![0; n as usize]; orderings_to_try.len()];

        let mut symbolics: Vec<Symbolic> = vec![S.clone(); orderings_to_try.len()];

        // rusolve oracle (highs-oracle-w5-hooks.patch): the hipo_ord record of
        // each ordering, buffered per index and written after the loop below, in
        // index order, from this thread (at threads > 1 the orderings run
        // concurrently in the for_each; the HiPO trace sink is thread_local).
        let mut oracle_ord_perm: Vec<Vec<Int>> = vec![Vec::new(); orderings_to_try.len()];
        let mut oracle_ord_failure: Vec<i8> = vec![0; orderings_to_try.len()];

        // Per-slot raw pointers (disjoint index): the C++ lambda
        // captures these vectors by reference and each call writes only its
        // own index i. None of the vectors is resized while the pointers live.
        let failure_p = SlotPtr(failure.as_mut_ptr());
        let permutations_p = SlotPtr(permutations.as_mut_ptr());
        let symbolics_p = SlotPtr(symbolics.as_mut_ptr());
        let oracle_ord_perm_p = SlotPtr(oracle_ord_perm.as_mut_ptr());
        let oracle_ord_failure_p = SlotPtr(oracle_ord_failure.as_mut_ptr());
        // C++ `FH_` (the member the lambda reaches through `this`), shared by
        // every call.
        let fh_p = SlotPtr(FH_ as *mut FHsolver);
        let skip_metis = options.amd_first && options.ordering == kHighsChooseString;
        let orderings_to_try_r: &Vec<String> = &orderings_to_try;
        let full_ptr_r: &Vec<Int> = &full_ptr;
        let full_rows_r: &Vec<Int> = &full_rows;

        // HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::chooseOrdering (lines 398-453)
        // lambda run_ordering_and_analyse
        let run_ordering_and_analyse = |i: Int| {
            let iu = i as usize;
            // rusolve: with AMD first (forced serial, the "choose" list
            // amd, metis), Metis is skipped when AMD succeeded with predicted flops below 1e9.
            // SAFETY: serial loop forced when skip_metis (below): nothing races with slot 0.
            if iu == 1 && skip_metis {
                let f0: i8 = unsafe { *failure_p.at(0) };
                let fl: f64 = unsafe { (*symbolics_p.at(0)).flops() };
                if f0 == 0 && fl < 1.0e9 {
                    unsafe { *failure_p.at(1) = 2 };
                    logger.printInfo(&format!(
                        "Skipped metis for {} (AMD flops {:.2e})\n",
                        nla, fl
                    ));
                    return;
                }
            }
            // SAFETY: disjoint index (C++ for_each grain 1, or the serial loop):
            // only the call for index i touches failure[i], permutations[i],
            // symbolics[i], oracle_ord_perm[i] and oracle_ord_failure[i]; the
            // vectors outlive the loop / for_each and are not resized.
            let failure_i: &mut i8 = unsafe { &mut *failure_p.at(iu) };
            // SAFETY: as above (disjoint index i).
            let permutation_i: &mut Vec<Int> = unsafe { &mut *permutations_p.at(iu) };
            // SAFETY: as above (disjoint index i).
            let symbolic_i: &mut Symbolic = unsafe { &mut *symbolics_p.at(iu) };

            logger.printInfo(&format!("Running {} for {}\n", orderings_to_try_r[iu], nla));

            if orderings_to_try_r[iu] == kHipoMetisString {
                let mut options_metis: Vec<idx_t> = vec![0; METIS_NOPTIONS];
                metis::set_default_options(&mut options_metis);
                options_metis[METIS_OPTION_SEED] = crate::ipm::hipo::factorhighs::kMetisSeed;

                options_metis[METIS_OPTION_DBGLVL] = 0;

                // no2hop improves the quality of ordering in general
                options_metis[METIS_OPTION_NO2HOP] = 1;

                let mut iperm: Vec<Int> = vec![0; n as usize];
                // C++ passes `&n` (the shared local); METIS only reads `*nvtxs`
                // (ometis.c METIS_NodeND), so a copy is the same call.
                let mut n_metis: Int = n;
                let status = metis::nodeND(
                    &mut n_metis,
                    full_ptr_r,
                    full_rows_r,
                    None,
                    Some(&options_metis),
                    permutation_i,
                    &mut iperm,
                );
                if status != METIS_OK {
                    *failure_i = 1;
                }
            } else if orderings_to_try_r[iu] == kHipoAmdString {
                let mut control: Vec<f64> = vec![0.0; AMD_CONTROL];
                amd::set_defaults(&mut control);
                let mut info: Vec<f64> = vec![0.0; AMD_INFO];

                let status = amd::order(
                    n,
                    full_ptr_r,
                    full_rows_r,
                    permutation_i,
                    Some(&control),
                    Some(&mut info),
                );
                if status != AMD_OK {
                    *failure_i = 1;
                }
            } else if orderings_to_try_r[iu] == kHipoRcmString {
                let status = rcm::genrcm(
                    n,
                    full_ptr_r[full_ptr_r.len() - 1],
                    full_ptr_r,
                    full_rows_r,
                    permutation_i,
                );
                if status != 0 {
                    *failure_i = 1;
                }
            } else {
                debug_assert!(false);
            }

            logger.printInfo(&format!(
                "Finished {} for {}\n",
                orderings_to_try_r[iu], nla
            ));
            // rusolve oracle: buffer the hipo_ord record of index i.
            // SAFETY: disjoint index i (as above).
            unsafe {
                *oracle_ord_perm_p.at(iu) = permutation_i.clone();
                *oracle_ord_failure_p.at(iu) = *failure_i;
            }

            if *failure_i != 0 {
                logger.printInfo(&format!(
                    "Error with {} for {}\n",
                    orderings_to_try_r[iu], nla
                ));
                return;
            }

            // SAFETY: C++ calls `FH_.analyse` on the shared member from every
            // ordering task; `FHsolver::analyse` only reads `nb_`, `logger_` and
            // `data_` (FactorHiGHS.cpp 37-43), so the concurrent calls do not
            // race. The `&mut` that the Rust signature asks for is aliased
            // across tasks; FH_ outlives the loop.
            *failure_i =
                unsafe { (*fh_p.at(0)).analyse(symbolic_i, rows, ptr, signs, permutation_i) } as i8;

            if *failure_i != 0 && logger.debug(2) {
                logger.print("Failed symbolic:");
                symbolic_i.print(logger, true);
            }
        };

        if options.parallel == kHighsOffString || skip_metis {
            for i in 0..orderings_to_try.len() as Int {
                run_ordering_and_analyse(i);
            }
        } else {
            highs_parallel::for_each(
                0,
                orderings_to_try.len() as HighsInt,
                &|start: HighsInt, _end: HighsInt| run_ordering_and_analyse(start),
                1,
            );
        }
        for i in 0..orderings_to_try.len() as Int {
            // rusolve oracle
            if failure[i as usize] == 2 {
                continue;
            }
            oracle_trace::hipo_ord(
                nla,
                &orderings_to_try[i as usize],
                n,
                &full_ptr,
                &full_rows,
                &oracle_ord_perm[i as usize],
                oracle_ord_failure[i as usize] != 0,
            );
        }

        let mut num_success: Int = 0;
        for b in failure.iter() {
            if *b == 0 {
                num_success += 1;
            }
        }

        if num_success > 0 {
            for i in 0..orderings_to_try.len() {
                if failure[i] == 0 {
                    logger.printInfo(&format!(
                        "{:>20} for {}: {} {}\n",
                        orderings_to_try[i],
                        nla,
                        c_fmt_e(symbolics[i].flops(), 2),
                        c_fmt_f(symbolics[i].size() as f64 / symbolics[i].sn() as f64, 2)
                    ));
                }
            }

            // find the ordering with best flops
            let mut best_flops: f64 = kHighsInf;
            for i in 0..orderings_to_try.len() {
                if failure[i] == 0 && symbolics[i].flops() < best_flops {
                    best_flops = symbolics[i].flops();
                }
            }

            // find orderings with flops within kFlopsOrderingThresh of the best
            let mut consider: Vec<Int> = Vec::new();
            for i in 0..orderings_to_try.len() {
                if failure[i] == 0 && symbolics[i].flops() <= kFlopsOrderingThresh * best_flops {
                    consider.push(i as Int);
                }
            }

            // among them, find the one with largest supernodes
            let mut best_sn_size: f64 = 0.0;
            let mut chosen: Int = -1;
            for &i in consider.iter() {
                let sn_avg_size: f64 =
                    symbolics[i as usize].size() as f64 / symbolics[i as usize].sn() as f64;
                if sn_avg_size > best_sn_size {
                    best_sn_size = sn_avg_size;
                    chosen = i;
                }
            }

            // fix selection if one or more require too much memory
            let _bytes_thresh: f64 = kLargeStorageGB * 1024.0 * 1024.0 * 1024.0;
            let mut best_memory: f64 = kHighsInf;
            let mut _ind_best_memory: Int = -1;
            for i in 0..orderings_to_try.len() {
                if symbolics[i].storage() < best_memory {
                    best_memory = symbolics[i].storage();
                    _ind_best_memory = i as Int;
                }
            }

            debug_assert!(chosen >= 0 && chosen < orderings_to_try.len() as Int);

            *S = std::mem::take(&mut symbolics[chosen as usize]);
            *ordering = orderings_to_try[chosen as usize].clone();
        }

        if num_success > 0 {
            kStatusOk
        } else {
            kStatusErrorAnalyse
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.

    use super::*;
    use crate::ipm::hipo::factorhighs::factor_highs::FHsolver;
    use crate::ipm::hipo::factorhighs::symbolic::Symbolic;
    use crate::ipm::hipo::ipm::options::Options;
    use crate::ipm::hipo::ipm::status::kStatusOk;
    use crate::lp_data::hconst::kHighsOnString;
    use crate::lp_data::highs_options::kHipoAmdString;
    use crate::lp_data::highs_options::kHipoMetisString;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    /// g x g 5-point grid (n = g * g), lower triangle with diagonal in CSC:
    /// column j lists rows i >= j (diagonal, right neighbour, lower neighbour).
    fn grid_lower(g: Int) -> (Vec<Int>, Vec<Int>) {
        let n = g * g;
        let (mut rows, mut ptr): (Vec<Int>, Vec<Int>) = (Vec::new(), vec![0 as Int]);
        for j in 0..n {
            let (r, c) = (j / g, j % g);
            rows.push(j);
            if c + 1 < g {
                rows.push(j + 1);
            }
            if r + 1 < g {
                rows.push(j + g);
            }
            ptr.push(rows.len() as Int);
        }
        (rows, ptr)
    }

    /// Everything `chooseOrdering` hands back, plus the hipo_ord trace of
    /// this thread: (status, ordering, flops bits, size, sn, storage bits, trace).
    type Outcome = (Int, String, u64, Int, Int, u64, String);

    fn run_choose(parallel: &str, ordering_opt: &str, nla: &str, ne_nz_limit: Int64) -> Outcome {
        let (rows, ptr) = grid_lower(10);
        let signs: Vec<Int> = vec![1; ptr.len() - 1];

        let mut S = Symbolic::new();
        let mut FH_ = FHsolver::new(None, 8);
        let mut options = Options::default();
        options.ordering = ordering_opt.to_string();
        options.parallel = parallel.to_string();
        let logger = Logger { log_options_: None };
        let mut ordering = String::new();

        oracle_trace::capture_begin();
        let status = FactorHiGHSSolver::chooseOrdering(
            &mut FH_,
            &options,
            &logger,
            ne_nz_limit,
            &rows,
            &ptr,
            &signs,
            &mut S,
            &mut ordering,
            nla,
        );
        let trace = oracle_trace::capture_end();
        (
            status,
            ordering,
            S.flops().to_bits(),
            S.size(),
            S.sn(),
            S.storage().to_bits(),
            trace,
        )
    }

    /// Runs `run_choose` on a fresh thread with its own `threads`-worker executor.
    fn run_on_threads(threads: i32, parallel: &'static str) -> Outcome {
        std::thread::spawn(move || {
            crate::parallel::highs_parallel::initialize_scheduler(threads);
            let out = run_choose(parallel, kHighsChooseString, "AS", Int64::MAX);
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap()
    }

    #[test]
    fn w5_hipo_ord_choose_ordering_for_each() {
        // reference: C++ lines 455-457 (parallel == "off"), the serial loop
        let reference = run_on_threads(1, kHighsOffString);
        assert_eq!(reference.0, kStatusOk);
        assert!(
            reference.1 == kHipoAmdString || reference.1 == kHipoMetisString,
            "ordering was: {}",
            reference.1
        );
        assert!(f64::from_bits(reference.2).is_finite());
        assert!(f64::from_bits(reference.2) > 0.0);

        // one hipo_ord record per ordering, in index order (amd, metis), each
        // before the failure test and analyse (failure = 0 here)
        let lines: Vec<&str> = reference.6.lines().collect();
        assert_eq!(lines.len(), 2, "trace: {}", reference.6);
        assert!(lines[0].contains("\"ordering\":\"amd\""), "{}", lines[0]);
        assert!(lines[1].contains("\"ordering\":\"metis\""), "{}", lines[1]);

        // C++ lines 458-461: for_each grain 1. Same result, and the same
        // trace (buffered per index, written after the for_each from the
        // calling thread), at every thread count and on every run.
        for threads in [1, 2, 4, 8] {
            for _run in 0..2 {
                let out = run_on_threads(threads, kHighsOnString);
                assert_eq!(out, reference, "threads = {threads}");
            }
        }
    }

    #[test]
    fn w5_hipo_ord_single_ordering_and_ne_interrupt() {
        // C++ lines 367-368: a named ordering is the only one tried
        let out = std::thread::spawn(|| {
            let out = run_choose(kHighsOnString, kHipoMetisString, "AS", Int64::MAX);
            HighsTaskExecutor::shutdown(true);
            out
        })
        .join()
        .unwrap();
        assert_eq!(out.0, kStatusOk);
        assert_eq!(out.1, kHipoMetisString);
        assert_eq!(out.6.lines().count(), 1);

        // C++ lines 380-385: NE with ptr.back() >= NE_nz_limit returns before
        // any ordering runs (no trace record, ordering untouched)
        let out = run_choose(kHighsOffString, kHighsChooseString, "NE", 0);
        assert_eq!(out.0, kStatusErrorAnalyse);
        assert_eq!(out.1, "");
        assert_eq!(out.6, "");
    }
}
