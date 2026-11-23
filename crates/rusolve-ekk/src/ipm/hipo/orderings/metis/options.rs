//! HiGHS source: extern/metis/libmetis/options.c, extern/metis/libmetis/auxapi.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Control structure setup. `SetupCtrl` mixes `double` literals into `real_t`
//! fields (`pfactor = 0.1*GETOPTION(..)`, `ubfactors[i] += 0.0000499`): the
//! expression is evaluated in `f64` and rounded once to `f32` on store.

use super::{ctrl_t, graph_t, idx_t, real_t};

/// HiGHS: options.c::SetupCtrl (lines 17-82)
/// `options` is `None` for the C `NULL`; returns `None` when `CheckParams`
/// fails (the C `FreeCtrl` + `NULL`).
pub fn SetupCtrl(
    optype: i32,
    options: Option<&[idx_t]>,
    ncon: idx_t,
    nparts: idx_t,
    tpwgts: Option<&[real_t]>,
    ubvec: Option<&[real_t]>,
) -> Option<Box<ctrl_t>> {
    // HiGHS: options.c lines 20-21: idx_t i, j; ctrl_t *ctrl;
    // C line 23: ctrl = (ctrl_t *)malloc(sizeof(ctrl_t));
    // C line 25: memset((void *)ctrl, 0, sizeof(ctrl_t));
    let mut ctrl = Box::new(ctrl_t::default());

    // HiGHS: options.c lines 27-43: switch (optype)
    match optype {
        super::METIS_OP_OMETIS => {
            ctrl.rtype = super::GETOPTION(
                options,
                super::METIS_OPTION_RTYPE,
                super::METIS_RTYPE_SEP1SIDED,
            );
            ctrl.iptype = super::GETOPTION(
                options,
                super::METIS_OPTION_IPTYPE,
                super::METIS_IPTYPE_EDGE,
            );
            ctrl.nseps = super::GETOPTION(options, super::METIS_OPTION_NSEPS, 1);
            ctrl.niter = super::GETOPTION(options, super::METIS_OPTION_NITER, 10);
            ctrl.ufactor = super::GETOPTION(
                options,
                super::METIS_OPTION_UFACTOR,
                super::OMETIS_DEFAULT_UFACTOR,
            );
            ctrl.compress = super::GETOPTION(options, super::METIS_OPTION_COMPRESS, 1);
            ctrl.ccorder = super::GETOPTION(options, super::METIS_OPTION_CCORDER, 0);
            // C line 36: ctrl->pfactor = 0.1*GETOPTION(options, METIS_OPTION_PFACTOR, 0);
            // 0.1 is a double literal; the product is f64, rounded to f32 on store.
            ctrl.pfactor = (0.1
                * (super::GETOPTION(options, super::METIS_OPTION_PFACTOR, 0) as f64))
                as real_t;

            ctrl.CoarsenTo = 100;
        }

        _ => {
            // C line 42: GK_ERREXIT("Unknown optype of %d\n", optype);
            super::gklib::gk_errexit(&format!("Unknown optype of {}\n", optype));
        }
    }

    // HiGHS: options.c lines 45-50: common options
    ctrl.ctype = super::GETOPTION(options, super::METIS_OPTION_CTYPE, super::METIS_CTYPE_SHEM);
    ctrl.no2hop = super::GETOPTION(options, super::METIS_OPTION_NO2HOP, 0);
    ctrl.seed = super::GETOPTION(options, super::METIS_OPTION_SEED, -1);
    ctrl.dbglvl = super::GETOPTION(options, super::METIS_OPTION_DBGLVL, 0);
    ctrl.dropedges = super::GETOPTION(options, super::METIS_OPTION_DROPEDGES, 0);

    // HiGHS: options.c lines 52-56: set non-option information
    ctrl.optype = optype;
    ctrl.ncon = ncon;
    ctrl.nparts = nparts;
    ctrl.maxvwgt = super::gklib::ismalloc(ncon as usize, 0);

    // HiGHS: options.c line 58: ctrl->rng_state = (ctrl->seed == -1 ? 4321 : ctrl->seed);
    ctrl.rng_state = if ctrl.seed == -1 {
        4321
    } else {
        ctrl.seed as u32
    };

    // HiGHS: options.c lines 61-66: setup the ubfactors
    // C line 62: ctrl->ubfactors = rsmalloc(ctrl->ncon, I2RUBFACTOR(ctrl->ufactor));
    // I2RUBFACTOR returns f64; rsmalloc takes real_t, so the C implicit conversion rounds once.
    ctrl.ubfactors = super::gklib::rsmalloc(
        ctrl.ncon as usize,
        super::I2RUBFACTOR(ctrl.ufactor) as real_t,
    );
    // C lines 63-64: if (ubvec) rcopy(ctrl->ncon, ubvec, ctrl->ubfactors);
    if let Some(ubvec) = ubvec {
        super::gklib::rcopy(ctrl.ncon as usize, ubvec, &mut ctrl.ubfactors);
    }
    // C lines 65-66: for (i=0; i<ctrl->ncon; i++) ctrl->ubfactors[i] += 0.0000499;
    // 0.0000499 is a double literal; the += is f64 arithmetic, result rounded to f32 on store.
    for i in 0..ctrl.ncon {
        ctrl.ubfactors[i as usize] = (ctrl.ubfactors[i as usize] as f64 + 0.0000499) as real_t;
    }

    // HiGHS: options.c lines 68-71: Allocate memory for balance multipliers
    // C line 71: ctrl->pijbm = rmalloc(nparts*ncon);
    ctrl.pijbm = super::gklib::rmalloc((nparts * ncon) as usize);

    // HiGHS: options.c line 73: IFSET(ctrl->dbglvl, METIS_DBG_INFO, PrintCtrl(ctrl));
    if ctrl.dbglvl & super::METIS_DBG_INFO != 0 {
        PrintCtrl(&ctrl);
    }

    // HiGHS: options.c lines 75-81: if (!CheckParams(ctrl)) { FreeCtrl(&ctrl); return NULL; } else { return ctrl; }
    if CheckParams(&ctrl) == 0 {
        FreeCtrl(&mut Some(ctrl));
        None
    } else {
        Some(ctrl)
    }
}

/// HiGHS: options.c::Setup2WayBalMultipliers (lines 87-95)
pub fn Setup2WayBalMultipliers(ctrl: &mut ctrl_t, graph: &graph_t, tpwgts: &[real_t]) {
    // HiGHS: options.c lines 89-90: idx_t i, j;
    // C lines 91-94: for (i=0; i<2; i++) { for (j=0; j<graph->ncon; j++) ctrl->pijbm[i*graph->ncon+j] = graph->invtvwgt[j]/tpwgts[i*graph->ncon+j]; }
    for i in 0..2 {
        for j in 0..graph.ncon {
            ctrl.pijbm[(i * graph.ncon + j) as usize] =
                graph.invtvwgt[j as usize] / tpwgts[(i * graph.ncon + j) as usize];
        }
    }
}

/// HiGHS: options.c::PrintCtrl (lines 101-164)
/// Debug printer (`dbglvl & METIS_DBG_INFO` only); never reached from HiPO.
pub fn PrintCtrl(ctrl: &ctrl_t) {
    // HIGHS_ORDERING_PRINT expands to nothing under NDEBUG (oracle build).
    // C lines 101-164: all output is debug-only.
    let _ = ctrl;
}

/// HiGHS: options.c::CheckParams (lines 170-239)
pub fn CheckParams(ctrl: &ctrl_t) -> i32 {
    // HiGHS: options.c lines 172-174: idx_t i, j; real_t sum; mdbglvl_et dbglvl=METIS_DBG_INFO;
    // C line 174: mdbglvl_et dbglvl=METIS_DBG_INFO;
    let _dbglvl: idx_t = super::METIS_DBG_INFO;

    // HiGHS: options.c lines 176-236: switch (ctrl->optype)
    match ctrl.optype {
        super::METIS_OP_OMETIS => {
            // C lines 178-181
            if ctrl.ctype != super::METIS_CTYPE_RM && ctrl.ctype != super::METIS_CTYPE_SHEM {
                // IFSET(dbglvl, METIS_DBG_INFO, ...) is off under NDEBUG.
                return 0;
            }
            // C lines 182-185
            if ctrl.iptype != super::METIS_IPTYPE_EDGE && ctrl.iptype != super::METIS_IPTYPE_NODE {
                return 0;
            }
            // C lines 186-189
            if ctrl.rtype != super::METIS_RTYPE_SEP1SIDED
                && ctrl.rtype != super::METIS_RTYPE_SEP2SIDED
            {
                return 0;
            }
            // C lines 190-193
            if ctrl.nseps <= 0 {
                return 0;
            }
            // C lines 194-197
            if ctrl.niter <= 0 {
                return 0;
            }
            // C lines 198-201
            if ctrl.ufactor <= 0 {
                return 0;
            }
            // C lines 202-205
            if ctrl.nparts != 3 {
                return 0;
            }
            // C lines 206-209
            if ctrl.ncon != 1 {
                return 0;
            }
            // C lines 210-213
            if ctrl.compress != 0 && ctrl.compress != 1 {
                return 0;
            }
            // C lines 214-217
            if ctrl.ccorder != 0 && ctrl.ccorder != 1 {
                return 0;
            }
            // C lines 218-221: if (ctrl->pfactor < 0.0)
            if ctrl.pfactor < 0.0 {
                return 0;
            }

            // C lines 223-229: for (i=0; i<ctrl->ncon; i++) { if (ctrl->ubfactors[i] <= 1.0) { ... return 0; } }
            for i in 0..ctrl.ncon {
                if ctrl.ubfactors[i as usize] <= 1.0 {
                    return 0;
                }
            }
        }

        _ => {
            // C lines 233-235: default: IFSET(...); return 0;
            return 0;
        }
    }

    // C line 238: return 1;
    1
}

/// HiGHS: options.c::FreeCtrl (lines 245-258)
pub fn FreeCtrl(r_ctrl: &mut Option<Box<ctrl_t>>) {
    // HiGHS: options.c line 247: ctrl_t *ctrl = *r_ctrl;
    // C line 249: FreeWorkSpace(ctrl);
    if let Some(ctrl) = r_ctrl.as_mut() {
        super::wspace::FreeWorkSpace(ctrl);
    }
    // C lines 251-255: gk_free of tpwgts, pijbm, ubfactors, maxvwgt, ctrl (Drop).
    // C line 257: *r_ctrl = NULL;
    *r_ctrl = None;
}

/// HiGHS: auxapi.c::Highs_METIS_SetDefaultOptions (lines 20-25)
pub fn Highs_METIS_SetDefaultOptions(options: &mut [idx_t]) -> i32 {
    // HiGHS: auxapi.c line 22: iset(METIS_NOPTIONS, -1, options);
    super::gklib::iset(super::METIS_NOPTIONS, -1, options);

    // C line 24: return METIS_OK;
    super::METIS_OK
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn default_options_are_minus_one() {
        let mut o = vec![0 as idx_t; crate::ipm::hipo::orderings::metis::METIS_NOPTIONS];
        assert_eq!(
            Highs_METIS_SetDefaultOptions(&mut o),
            crate::ipm::hipo::orderings::metis::METIS_OK
        );
        assert!(o.iter().all(|&v| v == -1));
    }

    #[test]
    fn setup_ctrl_ometis_defaults() {
        let ctrl = SetupCtrl(
            crate::ipm::hipo::orderings::metis::METIS_OP_OMETIS,
            None,
            1,
            3,
            None,
            None,
        )
        .expect("valid");
        assert_eq!(
            (
                ctrl.rtype,
                ctrl.nseps,
                ctrl.niter,
                ctrl.compress,
                ctrl.CoarsenTo
            ),
            (1, 1, 10, 1, 100)
        );
        assert_eq!(ctrl.rng_state, 4321);
    }
}
