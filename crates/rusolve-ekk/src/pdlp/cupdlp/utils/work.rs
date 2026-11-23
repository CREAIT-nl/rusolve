//! HiGHS source: pdlp/cupdlp/cupdlp_utils.c (lines 1032-1158)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `vec_Alloc`, `PDHG_Alloc`, `PDHG_Create`, `PDHG_Destroy`,
//! `PDHG_Init_Data` and the clock. `CUPDLP_TIMER` is defined
//! (`cupdlp_defs.h` line 7), so `my_clock` is `time(NULL)`: whole seconds
//! of wall time. The clock feeds only `dSolvingTime > dTimeLim` tests,
//! and `dTimeLim` is `HighsOptions::time_limit` (infinite by default), so
//! traces are deterministic unless a finite time limit is set.

use super::super::defs::*;
use super::*;

/// HiGHS: cupdlp_utils.c::vec_Alloc (lines 1032-1043)
///
/// `cupdlp_init_zero_vec_double` is `calloc`.
pub fn vec_Alloc(vec: &mut CUPDLPvec, n: cupdlp_int) -> cupdlp_retcode {
    let retcode = RETCODE_OK;
    vec.data = vec![0.0; n.max(0) as usize];
    vec.len = n;
    retcode
}

/// HiGHS: cupdlp_utils.c::PDHG_Alloc (lines 1045-1096)
///
/// CPU build: `buffer3` is allocated (`defined(CUPDLP_CPU)`), both
/// `buffer2` and `buffer3` with `MAX(2048, MAX(nCols, nRows))` zeros; the
/// `begin = getTimeStamp()` it takes feeds only a GPU timer and is dropped.
pub fn PDHG_Alloc(w: &mut CUPDLPwork) -> cupdlp_retcode {
    let mut retcode: cupdlp_retcode;

    'exit_cleanup: {
        // CUPDLP_INIT_SETTINGS(w->settings, 1);
        // CUPDLP_INIT_RESOBJ(w->resobj, 1);
        // CUPDLP_INIT_ITERATES(w->iterates, 1);
        // CUPDLP_INIT_STEPSIZE(w->stepsize, 1);
        // CUPDLP_INIT_TIMERS(w->timers, 1);
        // Allocation that cannot fail in Rust becomes Default; NULL tests disappear.
        w.settings = CUPDLPsettings::default();
        w.resobj = CUPDLPresobj::default();
        w.iterates = CUPDLPiterates::default();
        w.stepsize = CUPDLPstepsize::default();
        w.timers = CUPDLPtimers::default();

        // CUPDLP_CALL(timers_Alloc(w->timers));
        retcode = timers_Alloc(&mut w.timers);
        if retcode != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // cupdlp_float begin = getTimeStamp();
        // Dropped: feeds only GPU timer (line 1091, inside #ifndef CUPDLP_CPU).

        // buffer
        // CUPDLP_INIT_CUPDLP_VEC(w->buffer, 1);
        w.buffer = CUPDLPvec::default();
        // CUPDLP_CALL(vec_Alloc(w->buffer, w->problem->data->nRows));
        retcode = vec_Alloc(&mut w.buffer, w.problem.data.nRows);
        if retcode != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // cupdlp_init_zero_vec_double(w->buffer2,
        //                      MAX(2048, MAX(w->problem->data->nCols, w->problem->data->nRows)));
        let n_cols = w.problem.data.nCols;
        let n_rows = w.problem.data.nRows;
        let size_buf = MAX(2048, MAX(n_cols, n_rows));
        w.buffer2 = vec![0.0; size_buf as usize];

        // #if defined(CUPDLP_CPU) || !(USE_KERNELS) || ...
        // CPU build: buffer3 is allocated.
        // cupdlp_init_zero_vec_double(w->buffer3,
        //                      MAX(2048, MAX(w->problem->data->nCols, w->problem->data->nRows)));
        w.buffer3 = vec![0.0; size_buf as usize];

        // for scaling
        // cupdlp_init_zero_vec_double(w->colScale, w->problem->data->nCols);
        w.colScale = vec![0.0; n_cols as usize];
        // cupdlp_init_zero_vec_double(w->rowScale, w->problem->data->nRows);
        w.rowScale = vec![0.0; n_rows as usize];

        // CUPDLP_CALL(settings_Alloc(w->settings));
        retcode = settings_Alloc(&mut w.settings);
        if retcode != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(resobj_Alloc(w->resobj, w->problem, w->problem->data->nCols,
        //                          w->problem->data->nRows));
        retcode = resobj_Alloc(&mut w.resobj, &w.problem, n_cols, n_rows);
        if retcode != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(iterates_Alloc(w->iterates, w->problem->data->nCols,
        //                            w->problem->data->nRows));
        retcode = iterates_Alloc(&mut w.iterates, n_cols, n_rows);
        if retcode != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // CUPDLP_CALL(stepsize_Alloc(w->stepsize));
        retcode = stepsize_Alloc(&mut w.stepsize);
        if retcode != RETCODE_OK {
            retcode = RETCODE_FAILED;
            break 'exit_cleanup;
        }

        // #ifndef CUPDLP_CPU ... #endif
        // GPU code dropped.
    }

    retcode
}

/// HiGHS: cupdlp_utils.c::PDHG_Create (lines 1098-1109)
///
/// `CUPDLPwork **ww` becomes the returned box; `lp` and `scaling` are
/// moved in.
pub fn PDHG_Create(
    ww: &mut Option<Box<CUPDLPwork>>,
    lp: CUPDLPproblem,
    scaling: CUPDLPscaling,
) -> cupdlp_retcode {
    let retcode = RETCODE_OK;

    // CUPDLP_INIT_ZERO_CUPDLP_WORK(*ww, 1);
    // Allocation that cannot fail in Rust becomes Default.
    *ww = Some(Box::new(CUPDLPwork::default()));

    // CUPDLPwork *w = *ww;
    let w = ww.as_mut().unwrap();

    // w->problem = lp;
    w.problem = lp;
    // w->scaling = scaling;
    w.scaling = scaling;

    // exit_cleanup:
    retcode
}

/// HiGHS: cupdlp_utils.c::PDHG_Destroy (lines 1111-1115)
///
/// `PDHG_Clear` is `Drop`: this takes the box and drops it.
pub fn PDHG_Destroy(w: &mut Option<Box<CUPDLPwork>>) {
    // if (w && *w) {
    //   PDHG_Clear(*w);
    // }
    // PDHG_Clear is Drop in Rust. Taking the box drops it.
    let _ = w.take();
}

/// HiGHS: cupdlp_utils.c::PDHG_Init_Data (line 1117) (empty upstream)
pub fn PDHG_Init_Data(work: &mut CUPDLPwork) {
    // Empty upstream.
}

/// HiGHS: cupdlp_utils.c::my_clock (lines 1138-1156)
///
/// `CUPDLP_TIMER` defined, not `WIN32`: `(double)time(NULL)`.
///
/// `time(NULL)` is the whole seconds
/// since the epoch.
pub fn my_clock() -> f64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as f64)
        .unwrap_or(0.0)
}

/// HiGHS: cupdlp_utils.c::getTimeStamp (line 1158)
///
pub fn getTimeStamp() -> f64 {
    my_clock()
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn skeleton_vec_alloc_and_clock() {
        let mut v = CUPDLPvec::default();
        assert_eq!(vec_Alloc(&mut v, 3), RETCODE_OK);
        assert_eq!((v.len, v.data.clone()), (3, vec![0.0; 3]));
        assert!(getTimeStamp() > 1.0e9);
        assert_eq!(getTimeStamp().fract(), 0.0);
    }

    #[test]
    fn pdhg_alloc_sizes_buffers() {
        let mut p = CUPDLPproblem::default();
        p.nRows = 3;
        p.nCols = 2;
        p.data.nRows = 3;
        p.data.nCols = 2;
        p.lower = vec![0.0; 2];
        p.upper = vec![1.0; 2];
        let mut ww = None;
        assert_eq!(
            PDHG_Create(&mut ww, p, CUPDLPscaling::default()),
            RETCODE_OK
        );
        let w = ww.as_mut().unwrap();
        assert_eq!(PDHG_Alloc(w), RETCODE_OK);
        assert_eq!(w.buffer.len, 3);
        assert_eq!(w.buffer2.len(), 2048);
        assert_eq!(w.buffer3.len(), 2048);
        assert_eq!(w.colScale.len(), 2);
        assert_eq!(w.rowScale.len(), 3);
        assert_eq!(w.stepsize.eLineSearchMethod, PDHG_ADAPTIVE_LINESEARCH);
        assert_eq!(w.settings.nIterLim, i32::MAX);
        assert!(getTimeStamp() >= 0.0);
        PDHG_Destroy(&mut ww);
        assert!(ww.is_none());
    }
}
